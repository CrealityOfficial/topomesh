#include "letter.h"


#define FLOATERR 0.0001f

namespace topomesh
{
	void concaveOrConvexOfFaces(MMeshT* mt, std::vector<int>& faces, bool concave)
	{

		trimesh::point ave_normal;
		for (int i = 0; i < faces.size(); i++)if (!mt->faces[faces[i]].IsD())
		{
			ave_normal += trimesh::trinorm(mt->faces[faces[i]].connect_vertex[0]->p, mt->faces[faces[i]].connect_vertex[1]->p, mt->faces[faces[i]].connect_vertex[2]->p);
			mt->faces[faces[i]].SetS();
			mt->faces[faces[i]].V0(0)->SetS();
			mt->faces[faces[i]].V0(1)->SetS();
			mt->faces[faces[i]].V0(2)->SetS();			
		}		
		ave_normal /= faces.size();

		for (int i = 0; i < faces.size(); i++)if (!mt->faces[faces[i]].IsD())
		{
			for (int j = 0; j < 3; j++)
			{
				for (MMeshFace* vv : mt->faces[faces[i]].V0(j)->connected_face)
					if (!vv->IsS())
					{
						mt->faces[faces[i]].V0(j)->SetV(); break;
					}
			}
		}
		for (MMeshVertex& v : mt->vertices)if (!v.IsD())
		{
			if (v.IsS())
			{
				if (!v.IsV())
					v.p -= 120*ave_normal;
					//continue;
				else
					splitPoint(mt, &v, 120*ave_normal);
					//continue;
					//v.p -= 120*ave_normal;
			}
		}

	}

	void splitPoint(MMeshT* mt, MMeshVertex* v, trimesh::point ori)
	{
		mt->appendVertex(trimesh::point(v->p - ori));
		for (MMeshFace* f : v->connected_face)if (!f->IsD())
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
		//for (MMeshVertex* vc : v->connected_vertex)if (!vc->IsD())
		for(unsigned i=0;i<v->connected_vertex.size();i++)if(!v->connected_vertex[i]->IsD())
		{
			if (v->connected_vertex[i]->IsA(1))
			{
				mt->appendFace(v->connected_vertex[i]->index, v->index, mt->vertices.size() - 1);
			}
			if (v->connected_vertex[i]->IsA(1) || v->connected_vertex[i]->IsA(2))
				mt->vertices.back().connected_vertex.push_back(v->connected_vertex[i]);
		}

		for (MMeshFace* f : v->connected_face)if (!f->IsD())
		{
			f->ClearV();
			f->connect_vertex[0]->ClearA();
			f->connect_vertex[1]->ClearA();
			f->connect_vertex[2]->ClearA();
		}
	}

	void lettering(MMeshT* mesh, const std::vector<ClipperLibXYZ::Paths>& paths, CameraParam& camera, const LetterParam& Letter, std::vector<int>* faceindex)
	{
		//std::vector<trimesh::point> worldPointOfScreen;
		//std::vector<int> embeding_vertex;
		//wordToWorldPoint(camera,param,paths,worldPointOfScreen);
		//for (trimesh::point& p : worldPointOfScreen)
		//{
		//	trimesh::point orient = p - camera.pos;
		//	if (intersectionTriangle(mesh, p, orient))
		//	{
		//		//得到映射到mesh的字体点
		//		embeding_vertex.push_back(mesh->vertices.back().index);
		//	}
		//}
		//std::vector<trimesh::ivec2> e;
		//mesh->getEdge(e);
		//trimesh::quaternion q = q.rotationTo(trimesh::vec3(0, 0, 1.0f), camera.lookAt);
		//trimesh::fxform fx = mmesh::fromQuaterian(q);
		//for (MMeshVertex& v : mesh->vertices)
		//{
		//	v.p = fx * v.p;
		//}

		loadCameraParam(camera);
		std::vector<std::vector<trimesh::point>> wordPos;
		wordToWorldPoint(camera, Letter, paths, wordPos);
		Eigen::Matrix4f viewMatrix;
		Eigen::Matrix4f projectionMatrix;
		getViewMatrixAndProjectionMatrix(camera, viewMatrix, projectionMatrix);
		std::vector<std::vector<trimesh::vec2>> wordScrennPos;
		getEmbedingPoint(wordPos, viewMatrix, projectionMatrix, wordScrennPos);
		embedingAndCutting(mesh, wordScrennPos, viewMatrix, projectionMatrix,camera);
		std::vector<int> facesIndex;
		polygonInnerFaces(mesh, wordScrennPos, facesIndex, camera);
		concaveOrConvexOfFaces(mesh, facesIndex);
		unTransformationMesh(mesh, viewMatrix, projectionMatrix);
	}
	void embedingAndCutting(MMeshT* mesh, std::vector<std::vector<trimesh::vec2>>& lines, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix, const CameraParam& camera)
	{
		
		for (MMeshVertex& v : mesh->vertices)if (!v.IsD())
		{
			Eigen::Vector4f vPoint = { v.p.x,v.p.y,v.p.z,1.0 };
			Eigen::Vector4f point = ProjectionMatrix * ViewMatrix * vPoint;
			v.p = trimesh::point(point.x() * (1.0f / point.w()), point.y() * (1.0f / point.w()), point.z() * (1.0f / point.w()));
		}

		auto crossProduct = [=](trimesh::vec2 p1, trimesh::vec2 p2) ->float {
			return p1.x * p2.y - p1.y * p2.x;
		};
#if 0
		for (MMeshFace& f : mesh->faces)
		{
			/*float a = f.normal ^ camera.dir;
			if (a > 0) continue;*/
			for (int i = 0; i < lines.size(); i++)
			{
				for (int j = 0; j < lines[i].size(); j++)
				{
					trimesh::vec2 v10 = trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y) - trimesh::vec2(lines[i][j].x, lines[i][j].y);
					trimesh::vec2 v20 = trimesh::vec2(f.V0(1)->p.x, f.V0(1)->p.y) - trimesh::vec2(lines[i][j].x, lines[i][j].y);
					trimesh::vec2 v30 = trimesh::vec2(f.V0(2)->p.x, f.V0(2)->p.y) - trimesh::vec2(lines[i][j].x, lines[i][j].y);
					trimesh::vec2 v12 = trimesh::vec2(f.V0(1)->p.x, f.V0(1)->p.y) - trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y);//0->1
					trimesh::vec2 v13 = trimesh::vec2(f.V0(2)->p.x, f.V0(2)->p.y) - trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y);//0->2
					trimesh::point vv01 = f.V1(0)->p - f.V0(0)->p;
					trimesh::point vv02 = f.V2(0)->p - f.V0(0)->p;
					if ((crossProduct(v10, v20) >= 0 && crossProduct(v20, v30) >= 0 && crossProduct(v30, v10) >= 0) ||
						(crossProduct(v10, v20) < 0 && crossProduct(v20, v30) < 0 && crossProduct(v30, v10) < 0))
					{
						if(j==0)
							f.SetV();
						Eigen::Matrix2f e;
						e << v12.x, v13.x, v12.y, v13.y;
						Eigen::Vector2f b = { lines[i][j].x - f.V0(0)->p.x ,lines[i][j].y - f.V0(0)->p.y };
						Eigen::Vector2f x = e.fullPivLu().solve(b);
						mesh->appendVertex(trimesh::point(f.V0(0)->p + x.x() * vv01 + x.y() * vv02));						
						f.uv_coord.push_back(trimesh::vec4(-1, mesh->vertices.back().index, j, i));
					}
				}
			}
		}
		std::vector<trimesh::ivec2> edge;
		mesh->getEdge(edge);
		for (int i = 0; i < lines.size(); i++)
		{
			for (int j = 0; j < lines[i].size(); j++)
			{
				std::vector<std::pair<float, trimesh::ivec2>> corsspoint;
				mesh->calculateCrossPoint(edge, std::make_pair(trimesh::point(lines[i][j].x, lines[i][j].y, 0), trimesh::point(lines[i][(j + 1) % lines[i].size()].x, lines[i][(j + 1) % lines[i].size()].y, 0)), corsspoint);
				for (std::pair<float, trimesh::ivec2>& cp : corsspoint)
				{
					bool cover = false;
					for (MMeshFace* f : mesh->vertices[cp.second.x].connected_face)if (!f->IsD())
						f->SetA();
					for (MMeshFace* f : mesh->vertices[cp.second.y].connected_face)if (!f->IsD())
						f->SetA();
					for (MMeshFace* f : mesh->vertices[cp.second.x].connected_face)if (f->IsA(2) && !f->IsD())
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
							f->uv_coord.push_back(trimesh::vec4(index1, mesh->vertices.size() - 1, j, i));
						else
							f->uv_coord.push_back(trimesh::vec4(index2, mesh->vertices.size() - 1, j, i));
					}
					for (MMeshFace* f : mesh->vertices[cp.second.x].connected_face)if (!f->IsD())
						f->ClearA();
					for (MMeshFace* f : mesh->vertices[cp.second.y].connected_face)if (!f->IsD())
						f->ClearA();

				}
			}
		}
		
#else
		//---------new-----
		for(int i = 0; i < lines.size(); i++)
		{
			for (int j = 0; j < lines[i].size(); j++)
			{
				std::vector<trimesh::ivec3> push_lines;
				for (MMeshFace& f : mesh->faces)
				{
					trimesh::vec2 v10 = trimesh::vec2(lines[i][j].x, lines[i][j].y)-trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y);
					trimesh::vec2 v20 = trimesh::vec2(lines[i][j].x, lines[i][j].y)-trimesh::vec2(f.V0(1)->p.x, f.V0(1)->p.y);
					trimesh::vec2 v30 = trimesh::vec2(lines[i][j].x, lines[i][j].y)-trimesh::vec2(f.V0(2)->p.x, f.V0(2)->p.y);
					trimesh::vec2 v12 = trimesh::vec2(f.V0(1)->p.x, f.V0(1)->p.y) - trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y);//1->2
					trimesh::vec2 v13 = trimesh::vec2(f.V0(2)->p.x, f.V0(2)->p.y) - trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y);//0->2
					trimesh::vec2 v23 = trimesh::vec2(f.V0(2)->p.x, f.V0(2)->p.y) - trimesh::vec2(f.V0(1)->p.x, f.V0(1)->p.y);//2->3
					trimesh::vec2 v31 = trimesh::vec2(f.V0(0)->p.x, f.V0(0)->p.y) - trimesh::vec2(f.V0(2)->p.x, f.V0(2)->p.y);//3->1
					trimesh::point vv01 = f.V1(0)->p - f.V0(0)->p;
					trimesh::point vv02 = f.V2(0)->p - f.V0(0)->p;
					if ((crossProduct(v12, v10) >= 0 && crossProduct(v23, v20) >= 0 && crossProduct(v31, v30) >= 0) ||
						(crossProduct(v12, v10) < 0 && crossProduct(v23, v20) < 0 && crossProduct(v31, v30) < 0))
					{
						if(j==0)
							f.SetV();
						Eigen::Matrix2f e;
						e << v12.x, v13.x, v12.y, v13.y;
						Eigen::Vector2f b = { lines[i][j].x - f.V0(0)->p.x ,lines[i][j].y - f.V0(0)->p.y };
						Eigen::Vector2f x = e.fullPivLu().solve(b);
						mesh->appendVertex(trimesh::point(f.V0(0)->p + x.x() * vv01 + x.y() * vv02));
						f.uv_coord.push_back(trimesh::vec4(-1, mesh->vertices.back().index, j, i));
					}
					std::vector<trimesh::ivec2> edge = {trimesh::ivec2(f.V0(0)->index,f.V1(0)->index),trimesh::ivec2(f.V1(0)->index,f.V2(0)->index) ,trimesh::ivec2(f.V2(0)->index,f.V0(0)->index) };
					std::vector<std::pair<float, trimesh::ivec2>> corsspoint;
					mesh->calculateCrossPoint(edge, std::make_pair(trimesh::point(lines[i][j].x, lines[i][j].y, 0), trimesh::point(lines[i][(j + 1) % lines[i].size()].x, lines[i][(j + 1) % lines[i].size()].y, 0)), corsspoint);
					std::vector<trimesh::vec4> dis;
					for (std::pair<float, trimesh::ivec2>& cp : corsspoint)
					{
						f.SetS();
						int index = f.getVFindex(&mesh->vertices[cp.second.x]);
						bool pass = false;
						for(trimesh::ivec3& lp: push_lines)
							if (std::min(cp.second.x, cp.second.y) == lp.x && std::max(cp.second.x, cp.second.y) == lp.y)
							{
								//f.uv_coord.push_back(trimesh::vec4(index, lp.z, j, i));
								dis.push_back(trimesh::vec4(index, lp.z, j, i));
								pass = true; break;
							}
						if (pass)
							continue;
						trimesh::point d = mesh->vertices[cp.second.y].p - mesh->vertices[cp.second.x].p;
						mesh->appendVertex(trimesh::point(mesh->vertices[cp.second.x].p + cp.first * d));
						//f.uv_coord.push_back(trimesh::vec4(index, mesh->vertices.size() - 1, j, i));
						dis.push_back(trimesh::vec4(index, mesh->vertices.size() - 1, j, i));
						push_lines.push_back(trimesh::ivec3(std::min(cp.second.x, cp.second.y), std::max(cp.second.x, cp.second.y), mesh->vertices.size() - 1));
					}
					if (dis.empty()) continue;
					else if (dis.size() == 1)
						f.uv_coord.push_back(dis[0]);
					else
					{
						float a = trimesh::distance2(lines[i][j],trimesh::vec2(mesh->vertices[dis[0].y].p.x, mesh->vertices[dis[0].y].p.y));
						float b = trimesh::distance2(lines[i][j], trimesh::vec2(mesh->vertices[dis[1].y].p.x, mesh->vertices[dis[1].y].p.y));
						if (a < b)
						{
							f.uv_coord.push_back(dis[0]); f.uv_coord.push_back(dis[1]);
						}
						else
						{
							f.uv_coord.push_back(dis[1]); f.uv_coord.push_back(dis[0]);
						}
					}
				}
			}
		}
#endif
		for (MMeshFace& f : mesh->faces)if (!f.IsD() && f.IsS())
		{				
#if 0
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
							if ((f.uv_coord[i].z == f.uv_coord[j].z)&& (f.uv_coord[i].w == f.uv_coord[j].w))
								lines.push_back(std::make_pair(f.uv_coord[i].y, f.uv_coord[j].y));
						}
					for (int i = 0; i < faceVertexSque.size(); i++)
						mesh->vertices[faceVertexSque[i].first].inner.clear();
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
						std::vector<int> verVertexSque;
						for (int j = 0; j < f.uv_coord.size(); j++)
						{
							if (f.uv_coord[j].x == i)
								verVertexSque.push_back(f.uv_coord[j].y);
						}
						std::sort(verVertexSque.begin(), verVertexSque.end(), [&](int a, int b)->bool {
							float ad = trimesh::distance2(mesh->vertices[a].p, f.V0(i)->p);
							float bd = trimesh::distance2(mesh->vertices[b].p, f.V0(i)->p);
							return ad < bd;
							});
						faceVertexSque.insert(faceVertexSque.end(), verVertexSque.begin(), verVertexSque.end());
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

#else
			mesh->deleteFace(f);
			if (f.V0(0)->index == 277)
				std::cout << "277" << std::endl;
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
				continue;
			}
						
			std::vector<int> faceVertexSque;
			for (int i = 0; i < f.connect_vertex.size(); i++)
			{
				f.connect_vertex[i]->SetV();
				faceVertexSque.push_back(f.connect_vertex[i]->index);
				std::vector<int> verVertexSque;
				for (int j = 0; j < f.uv_coord.size(); j++)
				{
					if (f.uv_coord[j].x == i)
						verVertexSque.push_back(f.uv_coord[j].y);
				}
				std::sort(verVertexSque.begin(), verVertexSque.end(), [&](int a, int b)->bool {
					float ad = trimesh::distance2(mesh->vertices[a].p, f.V0(i)->p);
					float bd = trimesh::distance2(mesh->vertices[b].p, f.V0(i)->p);
					return ad < bd;
					});
				faceVertexSque.insert(faceVertexSque.end(), verVertexSque.begin(), verVertexSque.end());
			}
			int l = 1;
			std::vector<std::vector<int>> innerPoint(f.uv_coord.size());
			std::vector<bool> curve(f.uv_coord.size(), false);
			std::vector<int>  lastIndex(f.uv_coord.size(), -1);
			if (!f.IsV())
			{
				int n=0;
				for (int i = 0; i < f.uv_coord.size(); i++)
				{
					mesh->vertices[f.uv_coord[i].y].SetU(l);
					lastIndex[l] = f.uv_coord[i].y;
					if (f.uv_coord[i].x != -1)
						n++;
					else
					{
						curve[l] = true;
						innerPoint[l].push_back(f.uv_coord[i].y);
					}
					if (n == 2)
					{
						n = 0; l++;
					}
				}			
			}
			else
			{				
				bool pass = false;
				int n=0;
				for (int i = 0; i < f.uv_coord.size(); i++)
				{
					if (f.uv_coord[i].x != -1&&pass==false)
					{
						pass = true;	
						continue;
					}
					if (pass)
					{
						mesh->vertices[f.uv_coord[i].y].SetU(l);
						lastIndex[l] = f.uv_coord[i].y;
						if (f.uv_coord[i].x != -1)						
							n++;						
						else
						{
							curve[l] = true;
							innerPoint[l].push_back(f.uv_coord[i].y);
						}
						if (n == 2)
						{
							n = 0; l++;
						}
					}
				}
				for (int i = 0; i < f.uv_coord.size(); i++)
				{
					mesh->vertices[f.uv_coord[i].y].SetU(l);
					lastIndex[l] = f.uv_coord[i].y;
					curve[l] = true;
					if (f.uv_coord[i].x != -1)
						break;
					innerPoint[l].push_back(f.uv_coord[i].y);
				}
				l++;
			}			
			std::vector<std::vector<int>> polygon;
			polygon.push_back(faceVertexSque);

			for (int u = 1; u < l; u++)
			{
				int polygon_size = polygon.size();
				for (int i = 0; i < polygon_size; i++)
				{
					int begin = -1, end = -1;
					for (int j = 0; j < polygon[i].size(); j++)
					{
						if (mesh->vertices[polygon[i][j]].IsU(u))
						{
							if (begin == -1)
								begin = j;
							else
							{
								end = j; break;
							}
						}
					}
					if (begin == -1 || end == -1) continue;
					std::vector<int> subpoly1;
					bool push1 = false;
					for (int j = begin; j <= end; j++)
					{
						subpoly1.push_back(polygon[i][j]);
						int getu = mesh->vertices[polygon[i][j]].GetU();
						if (!mesh->vertices[polygon[i][j]].IsV() && getu > u)
							push1 = true;
					}
					if (curve[u])
					{
						if (subpoly1[subpoly1.size() - 1] == lastIndex[u])
							for (int j = innerPoint[u].size() - 1; j > -1; j--)
								subpoly1.push_back(innerPoint[u][j]);
						else
							subpoly1.insert(subpoly1.end(), innerPoint[u].begin(), innerPoint[u].end());
					}
					if (push1)
					{
						polygon.push_back(subpoly1);
					}
					else {					
						if (subpoly1.size() == 3)
							mesh->appendFace(subpoly1[0], subpoly1[1], subpoly1[2]);						
						else
						{							
							fillTriangle(mesh, subpoly1);							
						}
					}

					std::vector<int> subpoly2;
					bool push2 = false;
					for (int j = end; j <= begin + polygon[i].size(); j++)
					{
						subpoly2.push_back(polygon[i][j % polygon[i].size()]);
						int getu = mesh->vertices[polygon[i][j % polygon[i].size()]].GetU();
						if (!mesh->vertices[polygon[i][j % polygon[i].size()]].IsV() && getu > u)
							push2 = true;
					}
					if (curve[u])
					{
						if (subpoly2[subpoly2.size() - 1] == lastIndex[u])
							for (int j = innerPoint[u].size() - 1; j > -1; j--)
								subpoly2.push_back(innerPoint[u][j]);
						else
							subpoly2.insert(subpoly2.end(), innerPoint[u].begin(), innerPoint[u].end());
					}
					if (push2)
					{
						polygon.push_back(subpoly2);
					}
					else {
						if (subpoly2.size() == 3)
							mesh->appendFace(subpoly2[0], subpoly2[1], subpoly2[2]);						
						else
						{							
							fillTriangle(mesh, subpoly2);
						}
					}
					if (push1 || push2)
					{
						polygon.erase(polygon.begin() + i);
						break;
					}
				}
			}
			for (int i = 0; i < faceVertexSque.size(); i++)
			{
				mesh->vertices[faceVertexSque[i]].ClearU();
				mesh->vertices[faceVertexSque[i]].ClearV();
			}
#endif
		}


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
			float pq = 1.0 / std::abs(pe ^ v01);
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

	void wordToWorldPoint(const CameraParam& camera, const LetterParam& letter, const std::vector<ClipperLibXYZ::Paths>& paths, std::vector<std::vector<trimesh::point>>& points)
	{
		int len = paths.size();
		points.resize(paths[0].size());//暂时单字
		trimesh::point m1 = getWorldPoint(camera, camera.p1);
		trimesh::point m2 = getWorldPoint(camera, trimesh::ivec2(camera.p2.x, camera.p1.y));
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
					points[j].push_back(point);
				}
			}

		}
	}

	trimesh::point getWorldPoint(const CameraParam& camera, trimesh::ivec2 p)
	{
		trimesh::point screenCenter = camera.pos + camera.n * camera.dir;
		trimesh::point left = camera.dir % camera.up;
		trimesh::normalize(left);
		std::pair<float, float> screenhw;
		getScreenWidthAndHeight(camera, screenhw);
		float x_ratio = (float)p.x / (float)camera.ScreenSize.y - 0.5f;
		float y_ratio = (float)p.y / (float)camera.ScreenSize.x - 0.5f;
		return trimesh::point(screenCenter - camera.up * y_ratio * screenhw.first + left * x_ratio * screenhw.second);
	}

	void polygonInnerFaces(MMeshT* mt, std::vector<std::vector<trimesh::vec2>>& poly, std::vector<int>& faceIndex, const CameraParam& camera)
	{
		for (MMeshFace& f : mt->faces)if (!f.IsD())
		{
			/*float a = f.normal ^ camera.dir;
			if (a > 0) continue;*/
			trimesh::point c = (f.V0(0)->p + f.V0(1)->p + f.V0(2)->p) / 3.0;
			int rayCorssPoint = 0;
			for (int i = 0; i < poly.size(); i++)//左射线
			{
				for (int j = 0; j < poly[i].size(); j++)
				{
					if (std::abs(poly[i][(j + 1) % poly[i].size()].y - poly[i][j].y) < FLOATERR) continue;
					if (c.y < std::min(poly[i][j].y, poly[i][(j + 1) % poly[i].size()].y)) continue;
					if (c.y > std::max(poly[i][j].y, poly[i][(j + 1) % poly[i].size()].y)) continue;
					double x = (c.y - poly[i][j].y) * (poly[i][(j + 1) % poly[i].size()].x - poly[i][j].x) / (poly[i][(j + 1) % poly[i].size()].y - poly[i][j].y) + poly[i][j].x;
					if (x > c.x)
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

	void getScreenWidthAndHeight(const CameraParam& camera, std::pair<float, float>& wh)
	{
		wh.first = 2.0f * camera.n * std::tanf(camera.fov * M_PIf / 2.0f / 180.0f);
		wh.second = wh.first * camera.aspect;
	}

	void getViewMatrixAndProjectionMatrix(const CameraParam& camera, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix)
	{
		trimesh::point undir = -1.0f * camera.dir;
		ViewMatrix << camera.right.x, camera.right.y, camera.right.z, -1.0f * camera.pos DOT camera.right,
			camera.up.x, camera.up.y, camera.up.z, -1.0f * camera.pos DOT camera.up,
			undir.x, undir.y, undir.z, -1.0f * camera.pos DOT undir,
			0, 0, 0, 1;
		std::pair<float, float> ScreenSize;// h w
		getScreenWidthAndHeight(camera, ScreenSize);
		ProjectionMatrix << 2.0f * camera.n / (1.0f * ScreenSize.second), 0, 0, 0,
			0, 2.0f * camera.n / (1.0f * ScreenSize.first), 0, 0,
			0, 0, -1.0f * (camera.n + camera.f) / (1.0f * (camera.f - camera.n)), -2.0f * camera.n * camera.f / (1.0f * (camera.f - camera.n)),
			0, 0, -1, 0;

	}

	void loadCameraParam(CameraParam& camera)
	{
		camera.dir = trimesh::normalized(camera.lookAt - camera.pos);
		trimesh::point undir = -1.0f * camera.dir;
		camera.right = trimesh::normalized(camera.dir % camera.up);
		trimesh::normalize(camera.up);
	}

	void getEmbedingPoint(std::vector<std::vector<trimesh::point>>& lines, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix, std::vector<std::vector<trimesh::vec2>>& poly)
	{
		poly.resize(lines.size());
		for (unsigned i = 0; i < lines.size(); i++)
		{
			for (trimesh::point& p : lines[i])
			{
				Eigen::Vector4f linePoint = { p.x,p.y,p.z,1.0 };
				Eigen::Vector4f point = ProjectionMatrix * ViewMatrix * linePoint;
				poly[i].push_back(trimesh::vec2(point.x() * (1.0f / point.w()), point.y() * (1.0f / point.w())));
			}
		}
	}

	void unTransformationMesh(MMeshT* mesh, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix)
	{
		for (MMeshVertex& v : mesh->vertices)if (!v.IsD())
		{
			Eigen::Vector4f vPoint = { v.p.x,v.p.y,v.p.z,1.0 };
			Eigen::Vector4f point = ViewMatrix.inverse() * ProjectionMatrix.inverse() * vPoint;
			v.p = trimesh::point(point.x() * (1.0f / point.w()), point.y() * (1.0f / point.w()), point.z() * (1.0f / point.w()));
		}
	}

	trimesh::TriMesh* letter(trimesh::TriMesh* mesh, const SimpleCamera& camera, const TriPolygons& polygons,
		LetterDebugger* debugger, ccglobal::Tracer* tracer)
	{
		mesh->need_adjacentfaces();
		mesh->need_neighbors();
		mesh->need_normals();
		MMeshT mt(mesh);
		CameraParam cp;
		cp.lookAt = camera.center;
		cp.pos = camera.pos;
		cp.up = camera.up;
		cp.n = camera.n; cp.f = camera.f;
		cp.fov = camera.fov; cp.aspect = camera.aspect;
		loadCameraParam(cp);
		Eigen::Matrix4f viewMatrix;
		Eigen::Matrix4f projectionMatrix;
		getViewMatrixAndProjectionMatrix(cp, viewMatrix, projectionMatrix);	
		/*viewMatrix << 0.998806, -0.04885, 0, 4.65661e-10,
			-0.0487726, -0.997224, 0.05627, -7.45058e-09,
			-0.00274879, -0.0562028, -0.998416, -1.87256,
			0, 0, 0, 1;
		projectionMatrix << 1.93406f, 0, 0, 0,
			0, 3.73205, 0, 0,
			0, 0, -1.00076, -2.29203,
			0, 0, -1, 0;*/
		std::cout << "ViewMatrix : " << std::endl;
		std::cout << viewMatrix << std::endl;
		std::cout << "ProjectionMatrix : " << std::endl;
		std::cout << projectionMatrix << std::endl;
		std::vector<std::vector<trimesh::vec2>> poly;
		poly.resize(polygons.size());		
		for (int i = 0; i < polygons.size(); i++) {
			
			for (int j = polygons[i].size() - 1; j > 0; j--)
			{
				poly[i].push_back(trimesh::vec2(polygons[i][j].x, polygons[i][j].y));
			}
		}
		embedingAndCutting(&mt, poly, viewMatrix, projectionMatrix,cp);
		std::vector<int> facesIndex;
		polygonInnerFaces(&mt, poly, facesIndex,cp);
		concaveOrConvexOfFaces(&mt, facesIndex);
		unTransformationMesh(&mt, viewMatrix, projectionMatrix);
		trimesh::TriMesh* newmesh = new trimesh::TriMesh();
		mt.mmesh2trimesh(newmesh);
		newmesh->write("visualizationmesh.ply");
		return newmesh;
	}

	void fillTriangle(MMeshT* mesh, std::vector<int>& vindex)
	{
		int size = vindex.size();
		if (size == 0) return;
		if (size == 3)
		{
			mesh->appendFace(vindex[0], vindex[1], vindex[2]); return;
		}
		/*auto crossProduct = [=](trimesh::vec2 p1, trimesh::vec2 p2) ->float {
			return p1.x * p2.y - p1.y * p2.x;
		};*/
		int index = -1;
		for (int i = 0; i < size; i++)
		{
			trimesh::vec2 v1 = trimesh::vec2(mesh->vertices[vindex[(i + 1) % size]].p.x, mesh->vertices[vindex[(i + 1) % size]].p.y) - trimesh::vec2(mesh->vertices[vindex[i]].p.x, mesh->vertices[vindex[i]].p.y);
			trimesh::vec2 v2 = trimesh::vec2(mesh->vertices[vindex[(i + size - 1) % size]].p.x, mesh->vertices[vindex[(i + size - 1) % size]].p.y) - trimesh::vec2(mesh->vertices[vindex[i]].p.x, mesh->vertices[vindex[i]].p.y);			
			float a = std::acosf(trimesh::normalized(v1) ^ trimesh::normalized(v2));
			float b = trimesh::point((mesh->vertices[vindex[(i + 1) % size]].p - mesh->vertices[vindex[i]].p)%(mesh->vertices[vindex[(i + size - 1) % size]].p- mesh->vertices[vindex[i]].p)).z;
			//if (a >= M_PIf-FLOATERR)
			if(b<=0)
				continue;
			bool pass = false;
			std::vector<trimesh::vec2> triangle = { trimesh::vec2(mesh->vertices[vindex[i]].p.x,mesh->vertices[vindex[i]].p.y),trimesh::vec2(mesh->vertices[vindex[(i + 1) % size]].p.x, mesh->vertices[vindex[(i + 1) % size]].p.y),
			trimesh::vec2(mesh->vertices[vindex[(i + size - 1) % size]].p.x, mesh->vertices[vindex[(i + size - 1) % size]].p.y) };
			for (int j = 0; j < size; j++)
			{
				if (j != i && j != ((i + 1) % size) && j != ((i + size - 1) % size))
				{										
					int c = 0;
					for (int k = 0; k < triangle.size(); k++)
					{
						if (std::abs(triangle[k].y - triangle[(k + 1) % 3].y) < FLOATERR) continue;
						if (mesh->vertices[vindex[j]].p.y < std::min(triangle[k].y, triangle[(k + 1) % 3].y))continue;
						if (mesh->vertices[vindex[j]].p.y > std::max(triangle[k].y, triangle[(k + 1) % 3].y)) continue;
						double x = (mesh->vertices[vindex[j]].p.y - triangle[k].y) * (triangle[(k + 1) % 3].x - triangle[k].x) / (triangle[(k + 1) % 3].y - triangle[k].y) + triangle[k].x;
						if (x > mesh->vertices[vindex[j]].p.x)
							c++;
					}
					if (c % 2 != 0)
					{
						pass = true;
						break;
					}
				}
			}
			if (pass)
				continue;
			else
			{
				mesh->appendFace(vindex[i], vindex[(i + 1) % size], vindex[(i + size - 1) % size]);
				index = i;
				break;
			}
		}
		if (index != -1)
		{
			vindex.erase(vindex.begin() + index);
			fillTriangle(mesh, vindex);
		}
		
	}
}