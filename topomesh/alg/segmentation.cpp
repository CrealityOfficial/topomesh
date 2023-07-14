#include "segmentation.h"

namespace topomesh {
	void SpectralClusteringCuts(MMeshT* mesh,std::vector<std::vector<int>>& bloack)
	{		
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		if (!mesh->is_HalfEdge()) mesh->init_halfedge();
		float ave_geod_dist;
		float ave_ang_dist;
		for (MMeshFace& f : mesh->faces)
		{
			trimesh::point c = (f.V0(0)->p + f.V1(0)->p + f.V2(0)->p)/3.f;		
			MMeshHalfEdge* halfedge_ptr = f.f_mhe;
			do
			{
				if (halfedge_ptr->IsS()) continue;
				trimesh::point m = (halfedge_ptr->edge_vertex.first->p + halfedge_ptr->edge_vertex.second->p) / 2.f;
				trimesh::point cf = (halfedge_ptr->opposite->indication_face->V0(0)->p + halfedge_ptr->opposite->indication_face->V1(0)->p +
					halfedge_ptr->opposite->indication_face->V2(0)->p) / 3.f;
				float dist = trimesh::dist(cf, m) + trimesh::dist(c, m);
				float angle = f.dihedral(halfedge_ptr->opposite->indication_face);
				float alph = angle > M_PI ? 0.2f : 1.0f;
				float ang_dist = alph * (1.f - std::cos(angle));
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
	}
}