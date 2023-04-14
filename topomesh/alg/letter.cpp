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

	void lettering(MMeshT* mesh, const std::vector<ClipperLibXYZ::Paths>& paths, const CameraParam& camera, const LetterParam& param, std::vector<int>* faceindex)
	{
		int n_path = paths.size();
		trimesh::ivec2 k_wh = camera.p2 - camera.p1;
		trimesh::vec2 s_wh = trimesh::vec2((n_path * param.height)/(k_wh.x*1.0f), param.height/ k_wh.y * 1.0f );
		std::vector<trimesh::ivec2> screenPoint;
		for (unsigned i=0;i<paths.size();i++)
			for (unsigned j = 0; j < paths[i].size(); j++)
				for (unsigned k = 0; k < paths[i][j].size(); j++)
					screenPoint.push_back(camera.p1 +trimesh::ivec2((paths[i][j][k].X + i * param.height) / s_wh.x, paths[i][j][k].Y / s_wh.y));
		
		trimesh::point screen_begin = camera.pos + camera.n * camera.look + (camera.t - camera.b) * camera.up + (camera.r - camera.l) * -camera.right;
		for (trimesh::ivec2 p : screenPoint)
		{
			trimesh::point world_point = screen_begin + p.x * camera.right * camera.dx + p.y * -camera.up * camera.dy;
			//if(intersectionTriangle(mesh, world_point, camera.look))
		}

	}


	void screenToWorldPoint(MMeshT* mesh, const CameraParam& camera)
	{
		Eigen::Matrix4f c;
		c << camera.right.x, camera.right.y, camera.right.z, -camera.pos DOT camera.right,
			camera.up.x, camera.up.y, camera.up.z, -camera.pos DOT camera.up,
			camera.look.x, camera.look.y, camera.look.z, -camera.pos DOT camera.look,
			0, 0, 0, 1;
		Eigen::Matrix4f m;
		m << 2.0 * camera.n / (camera.r - camera.l), 0, -(camera.l + camera.r) * 1.0 / (camera.r - camera.l), 0,
			0, 2.0 * camera.n / (camera.t - camera.b), -(camera.t+ camera.b) * 1.0 / (camera.t - camera.b), 0,
			0, 0, 2.0 * (camera.n + camera.f) / (camera.f - camera.n), -2.0 * camera.n * camera.f / (camera.f - camera.n),
			0, 0, 1, 0;
		Eigen::Matrix4f s;

	}

	bool intersectionTriangle(MMeshT* mt, trimesh::point p, trimesh::point normal)
	{
		if (!mt->is_FaceNormals()) mt->getFacesNormals();
		for (MMeshFace& f : mt->faces)
		{
			float a = f.normal ^ normal;
			if (a > 0) continue;
			trimesh::point v01 = f.V0(1)->p - f.V0(0)->p;
			trimesh::point v02 = f.V0(2)->p - f.V0(0)->p;
			trimesh::point v0p = p - f.V0(0)->p;

			trimesh::point pe = normal % v02;
			trimesh::point pd = v0p % v01;
			float pq = 1.0/std::abs(pe ^ v01);
			float t = pq * (pd ^ v02);
			float u = pq * (pe ^ v0p);
			float v = pq * (pd ^ normal);
			if (u > 0 && v > 0 && (u + v) < 1)
			{
				mt->appendVertex(trimesh::point(f.V0(0)->p + u * v01 + v * v02));
				return true;
			}
		}
		return false;
	}
}