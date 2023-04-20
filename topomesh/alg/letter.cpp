#include "letter.h"


#define FLOATERR 0.00001f

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
						v.SetB();break;
					}
		}
		for (MMeshVertex& v : mt->vertices)if(!v.IsD())
		{
			if (v.IsS())
			{
				mt->appendVertex(trimesh::point(v.p + 20 * ave_normal));
				if (!v.IsB())
					v.p -= 50 * ave_normal;					
				else
					splitPoint(mt, &v, ave_normal);					
			}
		}

	}

	void splitPoint(MMeshT* mt, MMeshVertex* v, trimesh::point ori)
	{
		mt->appendVertex(trimesh::point(v->p - 50*ori));		
		for (MMeshFace* f : v->connected_face)if(!f->IsD())
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

		for (MMeshFace* f : v->connected_face)if (!f->IsD())
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
		for (MMeshVertex* vc : v->connected_vertex)if (!vc->IsD())
		{			
			if (vc->IsA(1))
				mt->appendFace(vc->index, v->index, mt->vertices.size() - 1);
			if (vc->IsA(1) || vc->IsA(2))
				mt->vertices.back().connected_vertex.push_back(vc);			
		}

		for (MMeshFace* f : v->connected_face)if (!f->IsD())
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
	//void embedingAndCutting(MMeshT* mesh, const CameraParam& camera, std::vector<trimesh::point>& lines)
	void embedingAndCutting(MMeshT* mesh, std::vector<trimesh::point>& lines)
	{
		/*Eigen::Matrix4f c;
		c << camera.right.x, camera.right.y, camera.right.z, -camera.pos DOT camera.right,
			camera.up.x, camera.up.y, camera.up.z, -camera.pos DOT camera.up,
			camera.look.x, camera.look.y, camera.look.z, -camera.pos DOT camera.look,
			0, 0, 0, 1;
		Eigen::Matrix4f m;
		m << 2.0 * camera.n / (camera.r - camera.l), 0, -(camera.l + camera.r) * 1.0 / (camera.r - camera.l), 0,
			0, 2.0 * camera.n / (camera.t - camera.b), -(camera.t+ camera.b) * 1.0 / (camera.t - camera.b), 0,
			0, 0, 2.0 * (camera.n + camera.f) / (camera.f - camera.n), -2.0 * camera.n * camera.f / (camera.f - camera.n),
			0, 0, 1, 0;*/
		/*Eigen::Matrix4f s;
		s << camera.w / 2.0f, 0, 0, camera.w / 2.0f,
			-camera.h / 2.0f, 0, 0, camera.h,
			0, 0, 1, 0,
			0, 0, 0, 0;*/
		std::vector<std::pair<int,trimesh::point>> line;
		/*for (trimesh::point& l : lines)
		{
			Eigen::Vector4f l4(l.x, l.y, l.z, 1);
			Eigen::Vector4f sl4 = m * c * l4;
			trimesh::point normal = l - camera.pos;
			if (intersectionTriangle(mesh, l, normal))
				line.push_back(std::make_pair(mesh->VN()-1, trimesh::point(sl4.x(), sl4.y(), 0)));
			else
				line.push_back(std::make_pair(-1, trimesh::point(sl4.x(), sl4.y(), 0)));
		}

		for (MMeshFace& f : mesh->faces)if (!f.IsD())
		{
			float a = f.normal ^ camera.look;
			if (a < 0)
			{
				f.V0(0)->SetS(); f.V1(0)->SetS(); f.V2(0)->SetS();
			}
		}*/
		std::vector<trimesh::ivec2> edge;
		/*mesh->getEdge(edge,true);
		for (MMeshVertex& v : mesh->vertices)if (!v.IsD())
		{
			Eigen::Vector4f v4(v.p.x, v.p.y, v.p.z, 1);
			Eigen::Vector4f s4 = m * c * v4;
			v.p = trimesh::point(s4.x(), s4.y(), s4.z());
			v.ClearS();
		}*/
		mesh->getEdge(edge);
		line.push_back(std::make_pair(1, trimesh::point(0.15f, 0.05f, 0)));
		line.push_back(std::make_pair(1, trimesh::point(-0.15f, 0.05f, 0)));
		line.push_back(std::make_pair(1, trimesh::point(-0.1f, -0.15f, 0)));
		line.push_back(std::make_pair(1, trimesh::point(0.2f, -0.15f, 0)));
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		trimesh::point ave_normal;
		for (MMeshFace& f : mesh->faces)
		{
			ave_normal += trimesh::normalized(f.normal);
		}
		ave_normal /= mesh->FN();
		ave_normal = -ave_normal;
		trimesh::quaternion q = q.rotationTo(trimesh::point(0, 0, 1), ave_normal);
		trimesh::fxform xf = trimesh::fromQuaterian(q);
		for (MMeshVertex& v : mesh->vertices)
		{
			v.p = xf * v.p;
		}
		std::vector<trimesh::point> l = { trimesh::point(0.15f,0.05f,0),trimesh::point(-0.15f,0.05f,0), trimesh::point(-0.1f, -0.15f, 0) ,trimesh::point(0.2f, -0.15f, 0) };
		auto crossProduct = [=](trimesh::vec2 p1, trimesh::vec2 p2) ->float {
			return p1.x * p2.y - p1.y * p2.x;
		};		
		for (MMeshFace& f : mesh->faces)
		{
			for (int i = 0; i < l.size(); i++)
			{
				trimesh::vec2 v10 = trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y) - trimesh::vec2(l[i].x, l[i].y);
				trimesh::vec2 v20 = trimesh::vec2(f.V0(1)->p.x, f.V0(1)->p.y) - trimesh::vec2(l[i].x, l[i].y);
				trimesh::vec2 v30 = trimesh::vec2(f.V0(2)->p.x, f.V0(2)->p.y) - trimesh::vec2(l[i].x, l[i].y);
				trimesh::vec2 v12 = trimesh::vec2(f.V0(1)->p.x, f.V0(1)->p.y) - trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y);//0->1
				trimesh::vec2 v13 = trimesh::vec2(f.V0(2)->p.x, f.V0(2)->p.y) - trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y);//0->2
				trimesh::point vv01 = f.V1(0)->p - f.V0(0)->p;
				trimesh::point vv02 = f.V2(0)->p - f.V0(0)->p;
				if ((crossProduct(v10, v20) >= 0 && crossProduct(v20, v30) >= 0 && crossProduct(v30, v10) >= 0) ||
					(crossProduct(v10, v20) < 0 && crossProduct(v20, v30) < 0 && crossProduct(v30, v10) < 0))
				{
					f.SetV();
					Eigen::Matrix2f e;
					e << v12.x, v13.x, v12.y, v13.y;
					Eigen::Vector2f b = { l[i].x - f.V0(0)->p.x ,l[i].y - f.V0(0)->p.y };				
					Eigen::Vector2f x = e.fullPivLu().solve(b);
					mesh->appendVertex(trimesh::point(f.V0(0)->p + x.x() * vv01 + x.y() * vv02));
					f.inner_vertex.push_back(mesh->VN() - 1);
				}
			}
		}
		for (int i = 0; i < line.size(); i++)
		{
			std::vector<std::pair<float, trimesh::ivec2>> corsspoint;
			mesh->calculateCrossPoint(edge, std::make_pair(line[i].second,line[(i+1)%line.size()].second), corsspoint);
			for(std::pair<float, trimesh::ivec2>& cp:corsspoint)
			{	
				bool cover = false;
				for (MMeshFace* f : mesh->vertices[cp.second.x].connected_face)if(!f->IsD())
					f->SetA();
				for (MMeshFace* f : mesh->vertices[cp.second.y].connected_face)if (!f->IsD())
					f->SetA();
				for (MMeshFace* f : mesh->vertices[cp.second.x].connected_face)if (f->IsA(2)&&!f->IsD())
				{									
					f->SetS();
					trimesh::point d = mesh->vertices[cp.second.y].p - mesh->vertices[cp.second.x].p;
					int index1 = f->getVFindex(&mesh->vertices[cp.second.x]);
					int index2 = f->getVFindex(&mesh->vertices[cp.second.y]);
					if (!cover)
					{
						mesh->appendVertex(trimesh::point(mesh->vertices[cp.second.x].p + cp.first * d)); cover = true;
					}
					if (f->V1(index1) == &mesh->vertices[cp.second.y])
						f->uv_coord.push_back(trimesh::vec3(index1, mesh->vertices.size() - 1,i));
					else
						f->uv_coord.push_back(trimesh::vec3(index2, mesh->vertices.size() - 1,i));
				}
				for (MMeshFace* f : mesh->vertices[cp.second.x].connected_face)if (!f->IsD())
					f->ClearA();
				for (MMeshFace* f : mesh->vertices[cp.second.y].connected_face)if (!f->IsD())
					f->ClearA();

			}
		}	
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		for (MMeshFace& f : mesh->faces)if(!f.IsD()&&f.IsS())
		{
			/*float a = f.normal ^ camera.look;
			if (a > 0) continue;*/			
			mesh->deleteFace(f);
			f.ClearS();
			if(!f.IsV())
				if (f.uv_coord.size() == 2)
				{
					if (f.V1(f.uv_coord[0].x)->index == f.V0(f.uv_coord[1].x)->index)
					{						
						mesh->appendFace(f.V0(f.uv_coord[0].x)->index, f.uv_coord[0].y, f.V2(f.uv_coord[0].x)->index);
						mesh->appendFace(f.uv_coord[0].y, f.V0(f.uv_coord[1].x)->index, f.uv_coord[1].y);
						mesh->appendFace(f.uv_coord[0].y, f.uv_coord[1].y, f.V2(f.uv_coord[0].x)->index);
					}
					else
					{					
						mesh->appendFace(f.V0(f.uv_coord[1].x)->index, f.uv_coord[1].y, f.V1(f.uv_coord[0].x)->index);
						mesh->appendFace(f.V0(f.uv_coord[0].x)->index, f.uv_coord[0].y, f.uv_coord[1].y);
						mesh->appendFace(f.V1(f.uv_coord[0].x)->index, f.uv_coord[0].y, f.uv_coord[1].y);
					}
				}
				else//----一个三角面被多条线切割			
				{										
					std::vector<std::pair<int,float>> faceVertexSque;					
					for (int j = 0; j < 3; j++)
					{					
						trimesh::point v1 = f.V1(j)->p - f.V0(j)->p;
						trimesh::point v2 = f.V2(j)->p - f.V0(j)->p;
						float a = std::acosf(trimesh::normalized(v1) ^ trimesh::normalized(v2));					
						faceVertexSque.push_back(std::make_pair(f.V0(j)->index,a));
						std::vector<std::pair<int, float>> verVertexSque;
						for (int i = 0; i < f.uv_coord.size(); i++)
						{						
							if (j == f.uv_coord[i].x)
							{								
								verVertexSque.push_back(std::make_pair((int)f.uv_coord[i].y,M_PI));
							}
						}
						std::sort(verVertexSque.begin(), verVertexSque.end(), [&](std::pair<int, float> a, std::pair<int, float> b)->bool {
							float ad = trimesh::distance2(mesh->vertices[a.first].p, f.V0(j)->p);
							float bd = trimesh::distance2(mesh->vertices[b.first].p, f.V0(j)->p);
							return ad < bd;						
							});
						faceVertexSque.insert(faceVertexSque.end(), verVertexSque.begin(), verVertexSque.end());						
					}	
					
					std::vector<std::pair<int, int>> lines;
					for(int i=0;i<f.uv_coord.size();i++)
						for (int j = i + 1; j < f.uv_coord.size(); j++)
						{
							if (f.uv_coord[i].z == f.uv_coord[j].z)
								lines.push_back(std::make_pair(f.uv_coord[i].y, f.uv_coord[j].y));
						}					
					for (int i = 0; i < lines.size(); i++)
					{
						for (int j = 0; j < faceVertexSque.size(); j++)
						{
							if (lines[i].first == faceVertexSque[j].first)
							{
								mesh->appendFace(lines[i].first, faceVertexSque[(j + 1) % faceVertexSque.size()].first, lines[i].second);								
								mesh->vertices[lines[i].first].inner.push_back(mesh->faces.size() - 1);  //是下标不是index face index  不一定等于下标
								mesh->vertices[faceVertexSque[(j + 1) % faceVertexSque.size()].first].inner.push_back(mesh->faces.size() - 1);
								mesh->vertices[lines[i].second].inner.push_back(mesh->faces.size() - 1);
							}
							if (lines[i].second == faceVertexSque[j].first)
							{
								mesh->appendFace(lines[i].first, lines[i].second, faceVertexSque[(j + 1) % faceVertexSque.size()].first);
								mesh->vertices[lines[i].first].inner.push_back(mesh->faces.size() - 1);
								mesh->vertices[faceVertexSque[(j + 1) % faceVertexSque.size()].first].inner.push_back(mesh->faces.size() - 1);
								mesh->vertices[lines[i].second].inner.push_back(mesh->faces.size() - 1);
							}
						}
					}				
					std::vector<int> last_face;
					for (int i = 0; i < faceVertexSque.size(); i++)
					{			
						float a = 0;
						MMeshVertex* v = &mesh->vertices[faceVertexSque[i].first];
						for (int j = 0; j < v->inner.size(); j++)
						{							
							trimesh::point v01= trimesh::normalized(mesh->faces[v->inner[j]].V1(v)->p-v->p);
							trimesh::point v02 = trimesh::normalized(mesh->faces[v->inner[j]].V2(v)->p-v->p);												
							a += std::acosf(v01 ^ v02);
						}
						if (a < faceVertexSque[i].second - FLOATERR || a > faceVertexSque[i].second + FLOATERR)
							last_face.push_back(faceVertexSque[i].first);
					}
					if(last_face.size()==3)
						mesh->appendFace(last_face[0], last_face[1], last_face[2]);
				}
			else//顶点			
			{
				if (f.inner_vertex.size() == 1)//只有一个点 直接剖分
				{
					std::vector<int> faceVertexSque ;
					for (int i = 0; i < f.connect_vertex.size(); i++)
					{
						faceVertexSque.push_back(f.V0(i)->index);
						for (int j = 0; j < f.uv_coord.size(); j++)
						{
							if (f.uv_coord[j].x == i)
								faceVertexSque.push_back(f.uv_coord[j].y);
						}
					}

					for (int i = 0; i < faceVertexSque.size(); i++)
					{
						mesh->appendFace(faceVertexSque[i], faceVertexSque[(i + 1) % faceVertexSque.size()], f.inner_vertex[0]);
					}
				}
				else//内部多个点
				{

				}
			}
			f.ClearV();
		}		
		/*for (MMeshVertex& v : mesh->vertices)if (!v.IsD())
		{
			Eigen::Vector4f v4(v.p.x, v.p.y, v.p.z, 1);
			Eigen::Vector4f s4 = c.inverse() * m.inverse() * v4;
			v.p = trimesh::point(s4.x(), s4.y(), s4.z());
			v.ClearS();
		}*/
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
				f.SetV();
				f.inner_vertex.push_back(mt->VN() - 1);
				//先不删除顶点face
				/*if(!f.IsD())
					mt->deleteFace(f.index);*/
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

	void polygonInnerFaces(MMeshT* mt, std::vector<trimesh::vec2>& poly, std::vector<int>& faceIndex)
	{
		for (MMeshFace& f : mt->faces)if(!f.IsD())
		{
			trimesh::point c= (f.V0(0)->p + f.V0(1)->p + f.V0(2)->p) / 3.0;			
			int rayCorssPoint = 0;
			for (int i = 0; i < poly.size(); i++)//左射线
			{
				if (std::abs(poly[(i + 1) % poly.size()].y - poly[i].y) < FLOATERR) continue;	
				if (std::min(poly[i].y, poly[(i + 1) % poly.size()].y) < c.y && std::max(poly[i].y, poly[(i + 1) % poly.size()].y) > c.y)
					if (std::min(poly[i].x, poly[(i + 1) % poly.size()].x) < c.x)					
					{
						float k = (poly[(i + 1) % poly.size()].y - poly[i].y) / (poly[(i + 1) % poly.size()].x - poly[i].x);
						if (c.y > k * (c.x - poly[i].x) + poly[i].y)
							rayCorssPoint++;
					}
			}
			if ((rayCorssPoint % 2) != 0)
			{
				faceIndex.push_back(f.index);
				//mt->deleteFace(f);
			}
			
		}
	}
}