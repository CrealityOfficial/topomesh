#include "segmentation.h"
#include "cmath"
#include "floyd.h"
#include "dijkstra.h"
#include "utils.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
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
				if (halfedge_ptr->IsS()) continue;
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
		Eigen::SparseMatrix<float> D;
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
		
	}


	void SpectralClusteringCuts::BlockSpectralClusteringCuts(MMeshT* mesh)
	{
#if 0
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		std::vector<float> face_smooth;
		//----is not isD------------
		for (MMeshFace& f : mesh->faces)
		{
			float angle;
			for (MMeshFace* ff : f.connect_face)
			{
				float arc = trimesh::normalized(f.normal) ^ trimesh::normalized(ff->normal);
				arc = arc >= 1.f ? 1.f : arc;
				arc = arc <= -1.f ? -1.f : arc;
				angle += (std::acos(arc) * 180 / M_PI);
			}
			angle = angle * 1.f / (f.connect_face.size() * 1.f);
			face_smooth.push_back(angle);
		}

		for (MMeshFace& f : mesh->faces)
		{
			if (f.IsS()) continue;
			f.SetS();
			std::vector<int> block;
			std::queue<int> queue;
			queue.push(f.index);
			while (!queue.empty())
			{
				block.push_back(queue.front());
				for (int i = 0; i < mesh->faces[queue.front()].connect_face.size(); i++)
				{
					if (!mesh->faces[queue.front()].connect_face[i]->IsS())
					{
						float agdiff = face_smooth[f.index] - face_smooth[mesh->faces[queue.front()].connect_face[i]->index];
						agdiff = std::abs(agdiff);
						if (agdiff < 5.f)
						{
							//block.push_back(mesh->faces[queue.front()].connect_face[i]->index);
							queue.push(mesh->faces[queue.front()].connect_face[i]->index);
							mesh->faces[queue.front()].connect_face[i]->SetS();
						}
					}
				}
				queue.pop();
			}

			result.push_back(block);
		}
#else
		for (MMeshFace& f : mesh->faces)
		{
			if (f.IsS()) continue;
			std::vector<int> block_face;
			for (MMeshFace* ff : f.connect_face)
			{
				float angle = f.dihedral(ff);
				if (angle < 180)
				{
					f.SetS(); ff->SetS();
				}
			}
		}

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