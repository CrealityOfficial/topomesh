#include "remesh.h"

namespace topomesh {
	void  SimpleRemeshing(MMeshT* mesh,const std::vector<int>& faceindexs,float threshlod)
	{
		if (faceindexs.empty()) return;
		if (!mesh->is_HalfEdge()) mesh->init_halfedge();
		for (int i : faceindexs)
			mesh->faces[i].SetS();		
		for (int i : faceindexs)
		{
			MMeshHalfEdge* halfedge = mesh->faces[i].f_mhe;
			int splitEdge = 0;
			do
			{			
				if (halfedge->IsS())
				{
					halfedge = halfedge->next; splitEdge++; continue;
				}
				float len = trimesh::dist2(halfedge->edge_vertex.first->p, halfedge->edge_vertex.second->p);
				if (len > threshlod * 4.f / 3.f)
				{
					halfedge->SetS();
					trimesh::point mid = (halfedge->edge_vertex.first->p + halfedge->edge_vertex.second->p) / 2.f;
					mesh->appendVertex(mid);
					halfedge->attritube = mesh->vertices.size() - 1;
					if (halfedge->opposite != nullptr)
					{
						halfedge->opposite->SetS();
						halfedge->opposite->attritube = mesh->vertices.size() - 1;
					}
					splitEdge++;
				}				
				halfedge = halfedge->next;
			} while (halfedge!=mesh->faces[i].f_mhe);
			if (splitEdge > 0) mesh->deleteFace(i);
			halfedge = mesh->faces[i].f_mhe;
			switch (splitEdge)
			{
			case 1:
				int appendVertex;
				do
				{
					if (halfedge->IsS())
						appendVertex = halfedge->attritube;
					halfedge = halfedge->next;
				} while (halfedge != mesh->faces[i].f_mhe);
				halfedge = mesh->faces[i].f_mhe;
				do
				{					
					if(!halfedge->IsS())
					{
						mesh->appendFace(halfedge->edge_vertex.first->index, halfedge->edge_vertex.second->index, appendVertex);
					}
					halfedge = halfedge->next;
				} while (halfedge!=mesh->faces[i].f_mhe);
				break;
			case  2:
				do
				{					
					if (halfedge->IsS())
					{
						int connect_vexter;
						if (halfedge->next->IsS()) connect_vexter = halfedge->next->attritube;
						if (halfedge->next->next->IsS()) connect_vexter = halfedge->next->next->attritube;
						mesh->appendFace(halfedge->edge_vertex.first->index, halfedge->attritube, connect_vexter);						
					}
					else
					{
						mesh->appendFace(halfedge->edge_vertex.first->index, halfedge->edge_vertex.second->index, halfedge->next->next->attritube);
					}
					halfedge = halfedge->next;
				} while (halfedge != mesh->faces[i].f_mhe);
				break;
			case 3:
				do
				{
					mesh->appendFace(halfedge->attritube, halfedge->edge_vertex.second->index, halfedge->next->attritube);					
					halfedge = halfedge->next;
				} while (halfedge!=mesh->faces[i].f_mhe);
				mesh->appendFace(halfedge->attritube, halfedge->next->attritube, halfedge->next->next->attritube);
				break;
			}

			halfedge = mesh->faces[i].f_mhe;
			do
			{
				if (halfedge->opposite == nullptr) {
					halfedge = halfedge->next;
					continue;
				}
				MMeshHalfEdge* heo = halfedge->opposite;
				if (halfedge->IsS() && !heo->indication_face->IsS())
				{
					mesh->deleteFace(heo->indication_face->index);
					mesh->appendFace(heo->next->edge_vertex.first->index, heo->next->edge_vertex.second->index, heo->attritube);
					mesh->appendFace(heo->next->next->edge_vertex.first->index, heo->next->next->edge_vertex.second->index, heo->attritube);
				}
				halfedge = halfedge->next;
			} while (halfedge != mesh->faces[i].f_mhe);

		}
		for (int i : faceindexs)
		{
			mesh->faces[i].ClearS();
			MMeshHalfEdge* halfedge = mesh->faces[i].f_mhe;
			halfedge->ClearS();
			halfedge->next->ClearS();
			halfedge->next->next->ClearS();
		}

	}
}