#include "letter.h"

namespace topomesh
{
	void concaveOrConvexOfFaces(MMeshT* mt, std::vector<int>& faces, bool concave)
	{
		trimesh::point ave_normal;
		for (int i = 0; i < faces.size(); i++)
		{
			ave_normal += trimesh::trinorm(mt->faces[faces[i]].connect_vertex[0]->p, mt->faces[faces[i]].connect_vertex[1]->p, mt->faces[faces[i]].connect_vertex[2]->p);
			mt->faces[faces[i]].SetS();
			mt->faces[faces[i]].V0(0)->SetS();
			mt->faces[faces[i]].V0(1)->SetS();
			mt->faces[faces[i]].V0(2)->SetS();
		}
		ave_normal /= faces.size();
		for (MMeshVertex& v : mt->vertices)
		{
			if (v.IsS())
				for (MMeshVertex* vv : v.connected_vertex)
					if (!vv->IsS())
					{
						v.SetB(); break;
					}
		}
		for (MMeshVertex& v : mt->vertices)
		{
			if (v.IsS())
			{
				if (!v.IsB())
					v.p -= ave_normal;
				else
					splitPoint(mt, &v, ave_normal);
			}
		}

	}

	void splitPoint(MMeshT* mt, MMeshVertex* v, trimesh::point ori)
	{
		mt->appendVertex(trimesh::point(v->p - ori));
		for (MMeshFace* f : v->connected_face)
		{
			f->SetV();
			if (f->IsS())
			{
				f->V1(v)->SetA(); f->V2(v)->SetA();
				int vin = f->getVFindex(v);
				f->connect_vertex[vin] = &mt->vertices.back();
				mt->vertices.back().connected_face.push_back(f);
			}
		}

		for (MMeshFace* f : v->connected_face)
		{
			std::vector<MMeshFace*>::iterator it;
			if (f->IsS())
				for (it = f->connect_face.begin(); it != f->connect_face.end();)
				{
					if ((*it)->IsV() && !(*it)->IsS())
						it = f->connect_face.erase(it);
					else
						it++;
				}
			if (!f->IsS())
				for (it = f->connect_face.begin(); it != f->connect_face.end(); )
				{
					if ((*it)->IsV() && (*it)->IsS())
						it = f->connect_face.erase(it);
					else
						it++;
				}
		}
		
		for (MMeshVertex* vc : v->connected_vertex)
		{
			if (vc->IsA(1))
				mt->appendFace(vc->index, v->index, mt->vertices.size() - 1);
			if (vc->IsA(1) || vc->IsA(2))
				mt->vertices.back().connected_vertex.push_back(vc);			
		}

		for (MMeshFace* f : v->connected_face)
		{
			f->ClearV();
			f->connect_vertex[0]->ClearA();
			f->connect_vertex[1]->ClearA();
			f->connect_vertex[2]->ClearA();
		}
	}
}