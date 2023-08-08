#include "segmentation.h"
#include "cmath"
#include "floyd.h"
#include "dijkstra.h"
#include "utils.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include<Eigen/Eigenvalues>
#include "queue"



namespace topomesh {
	SpectralClusteringCuts::SpectralClusteringCuts(MMeshT* mesh, float delte, float eta):_delta(delte),_eta(eta)
	{		
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		if (!mesh->is_HalfEdge()) mesh->init_halfedge();
		float ave_geod_dist;
		float ave_ang_dist;
		int edge = 0;
		for (MMeshFace& f : mesh->faces)
		{
			trimesh::point c = (f.V0(0)->p + f.V1(0)->p + f.V2(0)->p)/3.f;		
			MMeshHalfEdge* halfedge_ptr = f.f_mhe;
			do
			{
				if (halfedge_ptr->IsS())
				{
					halfedge_ptr = halfedge_ptr->next;
					continue;
				}
				if (halfedge_ptr->opposite == nullptr || halfedge_ptr->opposite->IsD())
				{
					halfedge_ptr->SetS();
					continue;
				}
				trimesh::point m = (halfedge_ptr->edge_vertex.first->p + halfedge_ptr->edge_vertex.second->p) / 2.f;
				trimesh::point cf = (halfedge_ptr->opposite->indication_face->V0(0)->p + halfedge_ptr->opposite->indication_face->V1(0)->p +
					halfedge_ptr->opposite->indication_face->V2(0)->p) / 3.f;
				float dist = trimesh::dist(cf, m) + trimesh::dist(c, m);
				float angle = f.dihedral(halfedge_ptr->opposite->indication_face);
				float alph = angle > M_PI/2.0 ? _eta : 1.0f;
				float ang_dist = alph * (1.f - std::cos(angle));
				edge++;
				ave_geod_dist += dist;
				ave_ang_dist += ang_dist;
				halfedge_ptr->attritube_vec.push_back(dist);
				halfedge_ptr->attritube_vec.push_back(ang_dist);
				halfedge_ptr->opposite->attritube_vec.push_back(dist);
				halfedge_ptr->opposite->attritube_vec.push_back(ang_dist);
				halfedge_ptr->SetS();
				halfedge_ptr->opposite->SetS();
				halfedge_ptr = halfedge_ptr->next;
			} while (halfedge_ptr!=f.f_mhe);			
		}

		ave_geod_dist = ave_geod_dist*1.f / edge * 1.f;
		ave_ang_dist = ave_ang_dist * 1.f / edge * 1.f;
		for (MMeshHalfEdge& he : mesh->half_edge)
			he.ClearS();		
		Eigen::SparseMatrix<float> D(mesh->faces.size(),mesh->faces.size());
		typedef Eigen::Triplet<float> Tr;
		std::vector<Tr> D_triplet;
		for (MMeshHalfEdge& he : mesh->half_edge)
		{
			if (!he.IsS())
			{
				float value = _delta * (he.attritube_vec[0] * 1.f / (ave_geod_dist * 1.f))
					+ (1.f - _delta) * (he.attritube_vec[1] * 1.f / (ave_ang_dist * 1.f));
				D_triplet.push_back(Tr(he.indication_face->index, he.opposite->indication_face->index, value));
				he.SetS();
			}
		}
		D.setFromTriplets(D_triplet.begin(),D_triplet.end());

		/*topomesh::Dijkstra<Eigen::SparseMatrix<float>> dijk(D);
		Eigen::SparseMatrix<float>  r = dijk.get_result();*/
		/*topomesh::floyd<Eigen::SparseMatrix<float>> floyd(D);
		Eigen::MatrixXf result = floyd.get_result();*/
		//------- 矩阵太大，计算时间长
		float sigma = D.sum() * 1.f / (1.f*std::pow(edge,2));
		Eigen::SparseMatrix<float> W;
		std::vector<Tr> W_triplet;
		std::vector<Tr> diag_triplet;
		for (int k = 0; k < D.outerSize(); ++k)
			for (Eigen::SparseMatrix<float>::InnerIterator it(D, k); it; ++it)
			{
				float value = it.value();
				value = std::exp((-1.f * value) / (2.f * std::pow(sigma, 2)));
				W_triplet.push_back(Tr(it.row(), it.col(), value));
			}

		W.setFromTriplets(W_triplet.begin(), W_triplet.end());
		for (int i = 0; i < D.rows(); i++)
			diag_triplet.push_back(Tr(i, i, 1));
		W.setFromTriplets(diag_triplet.begin(), diag_triplet.end());
		Eigen::VectorXf rowSums = W * Eigen::VectorXf::Ones(W.cols());
		Eigen::SparseMatrix<float> LD;
		std::vector<Tr> LD_triplet;
		for (int i = 0; i < rowSums.size(); i++)
			LD_triplet.push_back(Tr(i,i,rowSums(i)));
		LD.setFromTriplets(LD_triplet.begin(), LD_triplet.end());
		Eigen::SparseMatrix<float> L;
		L = ((((W*LD).pruned()).transpose()*LD).pruned()).transpose();
		
		Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<float>> es(L);
		std::cout<<"vec : "<<es.eigenvectors()<<"\n";
	}


	void SpectralClusteringCuts::BlockSpectralClusteringCuts(MMeshT* mesh)
	{
		//还没所有清除标记位
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		//if (!mesh->is_HalfEdge()) mesh->init_halfedge();
#if 1
		int make_user = 1;//user <=65536
		for (MMeshFace& f : mesh->faces)
		{
			if (f.IsV()) continue;
			topomesh::FacePatch path;
			findNeighborFacesOfSameAsNormal(mesh,f.index,path,10.f,true);
			for (int i : path)
			{
				mesh->faces[i].SetV();
				mesh->faces[i].SetU(make_user);
			}
			make_user++;
			result.push_back(path);
		}

		/*for (MMeshHalfEdge& he : mesh->half_edge)
		{

		}*/
#else
		for (topomesh::MMeshHalfEdge& he : mesh->half_edge)
		{
			if (he.IsS() || he.opposite==nullptr) continue;
			he.SetS(); he.opposite->SetS();
			float angle = he.indication_face->dihedral(he.opposite->indication_face);
			if (angle < 170.f)
			{
				he.indication_face->SetA();
				he.opposite->indication_face->SetA();
			}
		}
		std::vector<int> faceindex;
		for (topomesh::MMeshFace& f : mesh->faces)
		{
			int a = f.getA();
			if (a>0)
				faceindex.push_back(f.index);
		}
		result.push_back(faceindex);
#endif // 0
	}

	void SpectralClusteringCuts::test()
	{

	}

	Segmentation::Segmentation(trimesh::TriMesh* mesh)
	{

	}

	Segmentation::~Segmentation()
	{

	}

	void Segmentation::autoSegment(int num)   // num  < 0 auto
	{

	}
	//
	int Segmentation::createGroup()
	{
		return -1;
	}

	void Segmentation::removeGroup(int index)
	{

	}

	void Segmentation::addSeed2Group(int groupIndex, int index)
	{

	}

	void Segmentation::addSeeds2Group(int groupIndex, const std::vector<int>& indices)
	{

	}

	void Segmentation::removeGroupSeed(int groupIndex, int index)
	{

	}

	const FacePatchs& Segmentation::currentPatches()
	{
		return m_patches;
	}
}