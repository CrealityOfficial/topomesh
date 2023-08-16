#include "segmentation.h"
#include "cmath"
#include "floyd.h"
#include "dijkstra.h"
#include "topomesh/clustering/k-means.h"
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
					halfedge_ptr = halfedge_ptr->next;
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
			if (!he.IsS()&&he.attritube_vec.size()==2)
			{
				/*std::cout << "face : " << he.indication_face->index << "to face : " << he.opposite->indication_face->index 
					<< " distance : " << he.attritube_vec[0] << " angle : " << he.attritube_vec[1] << "\n";*/
				/*float value = _delta * (he.attritube_vec[0] * 1.f / (ave_geod_dist * 1.f))
					+ (1.f - _delta) * (he.attritube_vec[1] * 1.f / (ave_ang_dist * 1.f));*/
				float value = he.attritube_vec[1]*50.f;
				D_triplet.push_back(Tr(he.indication_face->index, he.opposite->indication_face->index, value));
				he.SetS();
			}
		}
		D.setFromTriplets(D_triplet.begin(),D_triplet.end());
		//std::cout << "D : \n" << D << "\n";
		/*topomesh::Dijkstra<Eigen::SparseMatrix<float>> dijk(D);
		Eigen::SparseMatrix<float>  r = dijk.get_result();*/
		/*topomesh::floyd<Eigen::SparseMatrix<float>> floyd(D);
		Eigen::MatrixXf result = floyd.get_result();*/
		//------- ����̫�󣬼���ʱ�䳤
	/*	float sigma = D.sum() * 1.f / (1.f*std::pow(edge,2));
		Eigen::SparseMatrix<float> W(D.rows(), D.cols());
		std::vector<Tr> W_triplet;
		std::vector<Tr> diag_triplet;
		for (int k = 0; k < D.outerSize(); ++k)
			for (Eigen::SparseMatrix<float>::InnerIterator it(D, k); it; ++it)
			{
				float value = it.value();
				std::cout << "before value : " << value << "\n";
				value = std::exp((-1.f * value) / (2.f * std::pow(sigma, 2)));
				std::cout << "value : " << value <<" row :"<<it.row()<<" col :"<< it.col() << "\n";
				W_triplet.push_back(Tr(it.row(), it.col(), value));
			}

		W.setFromTriplets(W_triplet.begin(), W_triplet.end());		
		std::cout << "W : \n" << W << "\n";*/
		Eigen::SparseMatrix<float> LD(D.rows(), D.cols());
		Eigen::VectorXf rowsum(D.rows());
		rowsum.setZero();
		std::vector<Tr> LD_triplet;
		for (int i = 0; i < D.outerSize(); i++)
		{
			for (Eigen::SparseMatrix<float>::InnerIterator it(D, i); it; ++it) {
				rowsum(it.row()) += it.value();
				/*std::cout << "value : " << it.value() << " row :" << it.row() << " col :" << it.col() << "\n";
				std::cout << "rowsum : " << rowsum(it.row()) << "\n";*/
			}
		}
		for (int i = 0; i < rowsum.size(); i++)
			LD_triplet.push_back(Tr(i,i,rowsum(i)));
		LD.setFromTriplets(LD_triplet.begin(), LD_triplet.end());
		//std::cout << "LD : \n" << LD << "\n";
		Eigen::SparseMatrix<float> L(D.rows(), D.cols());
		L = LD - D;
		//L = ((((W*LD).pruned()).transpose()*LD).pruned()).transpose();
		
		Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<float>> es(L);
		//std::cout << "L : \n" << L << "\n";
		//std::cout << "value : \n" << es.eigenvalues() << "\n";
		//std::cout<<"vec : \n"<<es.eigenvectors()<<"\n";
		//std::cout << "block : \n" << es.eigenvectors().block(0,1,L.rows(),3)<<"\n";
		Eigen::MatrixXf block = es.eigenvectors().block(0, 1, L.rows(), 4);
		result.resize(16);
		for (int i = 0; i < block.rows(); i++)
		{
			int index = 0;
			for (int j = 0; j < block.cols(); j++)
			{
				if (block(i, j) > 0)
					index += std::pow(2, j);
			}
			result[index].push_back(i);
		}
		for (int i = 0; i < result.size(); i++)
		{
			if (result[i].empty())
			{
				result.erase(result.begin() + i); i--;
			}
		}
		/*for (int i = 0; i < result.size(); i++)
		{
			for (int j = 0; j < result[i].size(); j++)
				std::cout << "i : " << i << " j : " << result[i][j] << "\n";
		}*/

		/*topomesh::kmeansClustering<Eigen::MatrixXf> kmeans(block, 9);
		this->result.resize(kmeans.get_result()->size());
		for (int i = 0; i < kmeans.get_result()->size(); i++)
		{
			for (int j = 0; j < kmeans.get_result()->at(i).size(); j++)
			{
				std::cout << "i : " << i << " j : " << kmeans.get_result()->at(i).at(j) << "\n";
				this->result[i].push_back(kmeans.get_result()->at(i).at(j));
			}
		}*/
		
	}
	
	void SpectralClusteringCuts::BlockSpectralClusteringCuts(MMeshT* mesh)
	{				
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		//if (!mesh->is_HalfEdge()) mesh->init_halfedge();
#if 1				
		std::vector<topomesh::FacePatch> container;
		int make_user = 1;//user <=65536
		for (MMeshFace& f : mesh->faces)
		{
			if (f.IsV()) continue;
			topomesh::FacePatch path;
			findNeighborFacesOfConsecutive(mesh,f.index,path,7.0f,true);
			for (int i : path)
			{
				mesh->faces[i].SetV();
				mesh->faces[i].SetU(make_user);
			}
			if (make_user < 65535)
				make_user++;
			container.push_back(path);
			//result.push_back(path);
		}	
		for (int i = 0; i < container.size(); i++)
		{
			if (container[i].size() <28)//�ñ�����ϵ
			{				
				for (int j = 0; j < container[i].size(); j++)
				{
					mesh->faces[container[i][j]].SetM();
				}				
			}
		}
	
		for (MMeshFace& f : mesh->faces)
		{
			if (f.IsM()&&!f.IsS())
			{
				int user=-1;
				std::vector<int> changeuser;
				std::queue<int> que;
				que.push(f.index);
				while (!que.empty())
				{
					mesh->faces[que.front()].SetS();
					for (MMeshFace* ff : mesh->faces[que.front()].connect_face)
					{
						if (!ff->IsS())
						{
							ff->SetS();
							if (ff->IsM())
								que.push(ff->index);
							else
								user = ff->GetU();
						}
					}
					changeuser.push_back(que.front());
					que.pop();
				}
				if(user!=-1)
					for (int i = 0; i < changeuser.size(); i++)
					{
						mesh->faces[changeuser[i]].SetU(user);
					}

			}
		}
		result.resize(make_user );
		for (MMeshFace& f : mesh->faces)
		{			
			result[f.GetU() - 1].push_back(f.index);
		}
		for (int i = 0; i < result.size(); i++)
		{
			if (result[i].empty())
			{
				result.erase(result.begin() + i);
				i--;
			}
		}

		for (MMeshFace& f : mesh->faces)
			f.ClearALL();
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
		Eigen::MatrixXf L(12,3);
		L << -0.291691, -0.173686, 0.367083,
			-0.286053, -0.175397, 0.370689,
			-0.291229, 0.40508, -0.0331007,
			-0.285586, 0.409053, -0.0334231,
			-0.291525, -0.230866, -0.334235,
			-0.285886, -0.233138, -0.337515,
			0.291694, 0.173695, -0.367075,
			0.286056, 0.175406, -0.370682,
			0.291523, 0.230861, 0.334242,
			0.285884, 0.233132, 0.337521,
			0.291226, -0.405083, 0.0330848,
			0.285583, -0.409055, 0.0334072;
		topomesh::kmeansClustering<Eigen::MatrixXf> kmean(L, 6);
		for (int i = 0; i < kmean.get_result()->size(); i++)
		{
			for (int j = 0; j < kmean.get_result()->at(i).size(); j++)
			{
				std::cout << "i : " << i << " j : " << kmean.get_result()->at(i).at(j) << "\n";
			}
		}
		/*Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(L);
		std::cout << "v : \n" << es.eigenvalues() << "\n";
		std::cout << "vec : \n" << es.eigenvectors() << "\n";*/
	}

	
	InteractionCuts::InteractionCuts(MMeshT* mesh, std::vector<int> indication, float alpth, float delta) :_alpth(alpth), _delta(delta)
	{
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		for (int i = 0; i < indication.size(); i++)
		{
			std::vector<int> container = { indication[i] };
			std::vector<int> next_container;
			int score = 18000;
			mesh->faces[container[0]].SetU(score);
			while (!container.empty())
			{
				for (int j = 0; j < container.size(); j++)
				{
					MMeshFace& f = mesh->faces[container[j]];
					f.SetS();
					for (MMeshFace* ff : f.connect_face)if(!ff->IsS())
					{
						ff->SetS();
						float angle = f.dihedral(ff);
						int new_score = f.GetU()-10;
						if (angle >= 180)
							new_score = new_score - _alpth*(angle - 180.f)*1.f/180.f ;
						else
							new_score = new_score - _delta * (angle - 180.f) * 1.f / 180.f;
						ff->SetU(new_score);
						next_container.push_back(ff->index);
					}
				}
				container.clear();
				container.swap(next_container);
			}
			for (MMeshFace& f : mesh->faces)
				f.ClearS();
		}
		for (MMeshFace& f : mesh->faces)
		{
			if (f.GetU() > 8000)
				result.push_back(f.index);
		}
	}

	SegmentationGroup::SegmentationGroup(int id, ManualSegmentation* segmentation)
			: m_segmentation(segmentation)
			, m_id(id)
	{
		
	}	


	Segmentation::Segmentation(trimesh::TriMesh* mesh)
	{

	}

	SegmentationGroup::~SegmentationGroup()
	{

	}

	bool SegmentationGroup::addSeed(int index)
	{
		std::vector<int> indices;
		indices.push_back(index);

		return addSeeds(indices);
	}

	bool SegmentationGroup::addSeeds(const std::vector<int>& indices)
	{
		for (int index : indices)
		{
			assert(m_segmentation->checkValid(index)
				&& m_segmentation->checkEmpty(index));
		}

		m_segmentation->replaceSeeds(indices, m_id);
		return true;
	}

	bool SegmentationGroup::removeSeed(int index)
	{
		assert(m_segmentation->checkValid(index)
			&& m_segmentation->checkIn(index, m_id));

		m_segmentation->removeSeed(index, m_id);
		return true;
	}

	ManualSegmentation::ManualSegmentation(TopoTriMeshPtr mesh)
		: m_mesh(mesh)
		, m_debugger(nullptr)
		, m_usedGroupId(0)
	{
		m_faceSize = (int)m_mesh->faces.size();
		m_groupMap.resize(m_faceSize, -1);
	}

	ManualSegmentation::~ManualSegmentation()
	{
		for (SegmentationGroup* group : m_groups)
			delete group;

		m_groups.clear();
	}

	void ManualSegmentation::setDebugger(ManualSegmentationDebugger* debugger)
	{
		m_debugger = debugger;
	}

	//
	SegmentationGroup* ManualSegmentation::addGroup()
	{
		SegmentationGroup* group = new SegmentationGroup(m_usedGroupId, this);
		m_groups.push_back(group);

		++m_usedGroupId;
		return group;
	}

	void ManualSegmentation::removeGroup(SegmentationGroup* group)
	{
		if (!group)
			return;

		std::vector<SegmentationGroup*>::iterator it = std::find(m_groups.begin(), m_groups.end(), group);
		if (it != m_groups.end())
		{
			int id = (*it)->m_id;
			for (int& index : m_groupMap)
			{
				if (index == id)
					index = -1;
			}
			m_groups.erase(it);
		}

		delete group;
	}

	SegmentationGroup* ManualSegmentation::checkIndex(int index)
	{
		if (!checkValid(index))
			return nullptr;

		int id = m_groupMap.at(index);
		for (SegmentationGroup* group : m_groups)
		{
			if (group->m_id == id)
				return group;
		}
		return nullptr;
	}

	void ManualSegmentation::segment()
	{

	}

	bool ManualSegmentation::checkOther(int index, int groupID)
	{
		return m_groupMap.at(index) == -1 || m_groupMap.at(index) == groupID;
	}

	bool ManualSegmentation::checkEmpty(int index)
	{
		return m_groupMap.at(index) == -1;
	}

	bool ManualSegmentation::checkIn(int index, int groupID)
	{
		return m_groupMap.at(index) == groupID;
	}

	bool ManualSegmentation::checkValid(int index)
	{
		return index >= 0 && index < m_faceSize;
	}

	void ManualSegmentation::replaceSeeds(const std::vector<int>& indices, int groupID)
	{
		for (int index : indices)
		{
			m_groupMap.at(index) = groupID;
		}
	}

	void ManualSegmentation::removeSeed(int index, int groupID)
	{
		m_groupMap.at(index) = -1;
	}
}