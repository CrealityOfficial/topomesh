#include "letter.h"

namespace topomesh
{
	void concaveOrConvexOfFaces(MMeshT* mt, std::vector<int>& faces, bool concave)
	{
		trimesh::point ave_normal;
		for (int i = 0; i < faces.size(); i++)if(!mt->faces[faces[i]].IsD())
		{
			ave_normal += trimesh::trinorm(mt->faces[faces[i]].connect_vertex[0]->p, mt->faces[faces[i]].connect_vertex[1]->p, mt->faces[faces[i]].connect_vertex[2]->p);
			mt->faces[faces[i]].SetS();
			mt->faces[faces[i]].V0(0)->SetS();
			mt->faces[faces[i]].V0(1)->SetS();
			mt->faces[faces[i]].V0(2)->SetS();
		}
		ave_normal /= faces.size();
		for (MMeshVertex& v : mt->vertices)if(!v.IsD())
		{
			if (v.IsS())
				for (MMeshVertex* vv : v.connected_vertex)
					if (!vv->IsS())
					{
						v.SetB(); break;
					}
		}
		for (MMeshVertex& v : mt->vertices)if(!v.IsD())
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
		std::vector<trimesh::point> worldPointOfScreen;
		std::vector<int> embeding_vertex;
		wordToWorldPoint(camera,param,paths,worldPointOfScreen);
		for (trimesh::point& p : worldPointOfScreen)
		{
			trimesh::point orient = p - camera.pos;
			if (intersectionTriangle(mesh, p, orient))
			{
				//得到映射到mesh的字体点
				embeding_vertex.push_back(mesh->vertices.back().index);
			}
		}
		std::vector<trimesh::ivec2> e;
		mesh->getEdge(e);
		trimesh::quaternion q = q.rotationTo(trimesh::vec3(0, 0, 1.0f), camera.look);
		trimesh::fxform fx = mmesh::fromQuaterian(q);
		for (MMeshVertex& v : mesh->vertices)
		{
			v.p = fx * v.p;
		}

	}


	void WorldPointToscreen(MMeshT* mesh, const CameraParam& camera, std::vector<trimesh::point>& lines)
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
		/*Eigen::Matrix4f s;
		s << camera.w / 2.0f, 0, 0, camera.w / 2.0f,
			-camera.h / 2.0f, 0, 0, camera.h,
			0, 0, 1, 0,
			0, 0, 0, 0;*/
		std::vector<trimesh::point > line;
		for (trimesh::point& l : lines)
		{
			Eigen::Vector4f l4(l.x, l.y, l.z, 1);
			Eigen::Vector4f sl4 = m * c * l4;
			line.push_back(trimesh::point(sl4.x(), sl4.y(), 0));
		}
		for (MMeshFace& f : mesh->faces)if (!f.IsD())
		{
			float a = f.normal ^ camera.look;
			if (a < 0)
			{
				f.V0(0)->SetS(); f.V1(0)->SetS(); f.V2(0)->SetS();
			}
		}
		std::vector<trimesh::ivec2> edge;
		mesh->getEdge(edge,true);
		for (MMeshVertex& v : mesh->vertices)if (!v.IsD())
		{
			Eigen::Vector4f v4(v.p.x, v.p.y, v.p.z, 1);
			Eigen::Vector4f s4 = m * c * v4;
			v.p = trimesh::point(s4.x(), s4.y(), s4.z());
			v.ClearS();
		}
		std::vector<trimesh::vec3> corsspoints;
		mesh->calculateCrossPoint(edge, line, corsspoints);
		
		/*if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		for (MMeshFace& f : mesh->faces)if(!f.IsD())
		{
			float a = f.normal ^ camera.look;
			if (a > 0) continue;
			for (int i = 0; i < 3; i++)
			{
				Eigen::Vector4f v4(f.V0(0)->p.x, f.V0(0)->p.y, f.V0(0)->p.z, 1);
				Eigen::Vector4f s4 = m*c * v4;
				if (s4.x() < camera.p2.x && s4.x() > camera.p1.x && s4.y() < camera.p1.y && s4.y() > camera.p1.y)
				{
					f.SetS();
					std::cout << " f index: " << f.index << "\n";
					break;
				}
			}
		}*/
	}

	bool intersectionTriangle(MMeshT* mt, trimesh::point p, trimesh::point normal)
	{
		if (!mt->is_FaceNormals()) mt->getFacesNormals();
		for (MMeshFace& f : mt->faces)if(!f.IsD())
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
				//mt->appendVertex(trimesh::point(f.V0(0)->p + u * v01 + v * v02));
				mt->deleteFace(f.index);
				return true;
			}
		}
		return false;
	}

	void wordToWorldPoint(const CameraParam& camera, const LetterParam& letter, const std::vector<ClipperLibXYZ::Paths>& paths, std::vector<trimesh::point>& points)
	{
		int len = paths.size();
		trimesh::point m1 = getWorldPoint(camera, camera.p1);
		trimesh::point m2 = getWorldPoint(camera, trimesh::ivec2(camera.p2.x,camera.p1.y));
		trimesh::point m3 = getWorldPoint(camera, trimesh::ivec2(camera.p1.x, camera.p2.y));
		trimesh::point m4 = getWorldPoint(camera, camera.p2);
		trimesh::point ori_x = trimesh::normalized(m2 - m1);
		trimesh::point ori_y = trimesh::normalized(m3 - m1);
		float w = trimesh::distance(m1, m2);
		float h = trimesh::distance(m1, m3);
		float wordSize = w / (len * 1.0f);
		float pointLen = wordSize / (letter.height * 1.0f);
		float span = (h - wordSize) / 2.0f;
		for (unsigned i = 0; i < len; i++)
		{
			trimesh::point begin = m1 + span * ori_y + i * wordSize * ori_x;
			for (unsigned j = 0; j < paths[i].size(); j++)
			{
				for (unsigned k = 0; k < paths[i][j].size(); k++)
				{
					trimesh::point point = begin + paths[i][j][k].X * pointLen * ori_x + (letter.height - paths[i][j][k].Y) * pointLen * ori_y;
					points.push_back(point);
				}
			}

		}
	}

	trimesh::point getWorldPoint(const CameraParam& camera,trimesh::ivec2 p)
	{
		trimesh::point dir = trimesh::normalized(camera.look);
		trimesh::point screenCenter = camera.pos + camera.n * dir;
		trimesh::point left = dir % camera.up;
		trimesh::normalize(left);
		float h = 2.0f * camera.n * std::tanf(camera.fov * M_PI / 2.0f / 180.0f);
		float w = h * camera.aspect;
		float x_ratio = (float)p.x / (float)camera.w - 0.5f;
		float y_ratio = (float)p.y / (float)camera.h - 0.5f;
		return trimesh::point(screenCenter - camera.up * y_ratio * h + left * x_ratio * w);
	}
}