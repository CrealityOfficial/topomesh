#include "new_letter.h"
#include "letter.h"
#include "earclipping.h"
#include "trimesh2/TriMesh_algo.h"

namespace topomesh {
	trimesh::TriMesh* CreateFontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter, float height,
		std::vector<float>& word_location , std::vector<int>& mesh_vertex_sizes, std::vector<trimesh::vec3>& word_mesh_center)
	{
		trimesh::TriMesh* _return_mesh = new trimesh::TriMesh();					
		mesh_vertex_sizes.push_back(0);

		//trimesh::TriMesh* bbx_center_mesh = new trimesh::TriMesh();

		for (int li = 0; li < letter.size(); li++)
		{
			MMeshT mt(5000, 10000);
			mt.set_VFadjacent(true);
			std::vector<std::vector<trimesh::vec2>> totalpoly = letter[li];
			trimesh::vec2 wordbbx_min(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
			trimesh::vec2 wordbbx_max(std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
			for(int pi=0;pi<totalpoly.size();pi++)
				for (int ppi = 0; ppi < totalpoly[pi].size(); ppi++)
				{
					if (totalpoly[pi][ppi].x < wordbbx_min.x)
						wordbbx_min.x = totalpoly[pi][ppi].x;
					if (totalpoly[pi][ppi].x > wordbbx_max.x)
						wordbbx_max.x = totalpoly[pi][ppi].x;

					if (totalpoly[pi][ppi].y < wordbbx_min.y)
						wordbbx_min.y = totalpoly[pi][ppi].y;
					if (totalpoly[pi][ppi].x > wordbbx_max.x)
						wordbbx_max.y = totalpoly[pi][ppi].y;
				}
			trimesh::TriMesh* _word_mesh = new trimesh::TriMesh();
			mt.appendVertex(trimesh::point(wordbbx_min.x-0.2f, wordbbx_max.y+0.2f, 0));
			mt.appendVertex(trimesh::point(wordbbx_min.x-0.2f, wordbbx_min.y-0.2f, 0));
			mt.appendVertex(trimesh::point(wordbbx_max.x+0.2f, wordbbx_min.y-0.2f, 0));
			mt.appendVertex(trimesh::point(wordbbx_max.x+0.2f, wordbbx_max.y+0.2f, 0));
			mt.appendFace(0, 1, 2);
			mt.appendFace(0, 2, 3);
			std::vector<int> faceindex;
			for (int i = 0; i < mt.faces.size(); i++)
				faceindex.push_back(i);
					
			setMark(totalpoly);
			embedingAndCutting(&mt, totalpoly, faceindex);
			faceindex.clear();
			for (int i = 0; i < mt.faces.size(); i++)
				faceindex.push_back(i);
			std::vector<int> outfacesIndex;
			std::vector<std::vector<std::vector<trimesh::vec2>>> word_letter = { totalpoly };
			polygonInnerFaces(&mt, word_letter, faceindex, outfacesIndex);

			std::vector<int> least;
			std::set_difference(faceindex.begin(), faceindex.end(), outfacesIndex.begin(), outfacesIndex.end(), std::inserter(least, least.begin()));
			for (int fi : least)
				mt.deleteFace(fi);

			for (MMeshVertex& v : mt.vertices)
				if (v.connected_face.empty())
					mt.deleteVertex(v);
			mt.shrinkMesh();
			mt.init_halfedge();
			mt.set_HalfEdge(false);
			for (MMeshVertex& v : mt.vertices)
				v.p.z += height;
			int face_size = mt.faces.size();
			std::vector<int> compara_table(mt.vertices.size());
			for (int fi = 0; fi < face_size; fi++)
			{
				MMeshFace& f = mt.faces[fi];
				if (f.IsD()) continue;
				int v0, v1, v2;
				if (!f.V0(0)->IsV())
				{
					mt.appendVertex(trimesh::point(f.V0(0)->p.x, f.V0(0)->p.y, 0));
					f.V0(0)->SetV();
					compara_table[f.V0(0)->index] = mt.vertices.size() - 1;
					v0 = mt.vertices.size() - 1;
				}
				else
				{
					v0 = compara_table[f.V0(0)->index];
				}

				if (!f.V1(0)->IsV())
				{
					mt.appendVertex(trimesh::point(f.V1(0)->p.x, f.V1(0)->p.y, 0));
					f.V1(0)->SetV();
					compara_table[f.V1(0)->index] = mt.vertices.size() - 1;
					v1 = mt.vertices.size() - 1;
				}
				else
				{
					v1 = compara_table[f.V1(0)->index];
				}

				if (!f.V2(0)->IsV())
				{
					mt.appendVertex(trimesh::point(f.V2(0)->p.x, f.V2(0)->p.y, 0));
					f.V2(0)->SetV();
					compara_table[f.V2(0)->index] = mt.vertices.size() - 1;
					v2 = mt.vertices.size() - 1;
				}
				else
				{
					v2 = compara_table[f.V2(0)->index];
				}
				mt.appendFace(v2, v1, v0);
				MMeshHalfEdge* p = f.f_mhe;
				int n = 3;
				while (1)
				{
					if (p->opposite == nullptr)
					{
						int i1 = p->edge_vertex.second->index;
						int i2 = p->edge_vertex.first->index;
						int i3 = 0;
						int i4 = 0;
						if (n == 3)
						{
							i3 = v0;
							i4 = v1;
						}
						else if (n == 2)
						{
							i3 = v1;
							i4 = v2;
						}
						else if (n == 1)
						{
							i3 = v2;
							i4 = v0;
						}
						
						mt.appendFace(i1, i2, i3);
						mt.appendFace(i3, i4, i1);
					}
					n--;
					if (p->next == f.f_mhe)
						break;
					p = p->next;
				}

			}
			mt.mmesh2trimesh(_word_mesh);
			
			_word_mesh->need_bbox();
			word_mesh_center.push_back(_word_mesh->bbox.center());
			//bbx_center_mesh->vertices.push_back(_word_mesh->bbox.center());
			int vsize = _return_mesh->vertices.size();	
			for (trimesh::point v : _word_mesh->vertices)
				_return_mesh->vertices.push_back(v);
			for (trimesh::TriMesh::Face f : _word_mesh->faces)
			{
				_return_mesh->faces.push_back(trimesh::TriMesh::Face(vsize+f[0], vsize + f[1], vsize + f[2]));
			}
			mesh_vertex_sizes.push_back(_return_mesh->vertices.size());
		}
		
		_return_mesh->need_bbox();
		trimesh::vec3 bbx_center = _return_mesh->bbox.center();
		trimesh::trans(_return_mesh, -bbx_center);
		for (int i = 0; i < word_mesh_center.size(); i++)
		{
			word_location.push_back((word_mesh_center[i].x - bbx_center.x));
			word_mesh_center[i]-=bbx_center;			
			//bbx_center_mesh->vertices.push_back(word_mesh_center[i]);
		}
		//bbx_center_mesh->write("bbx_center.ply");
		//_return_mesh->write("returnmesh.ply");
		return _return_mesh;
	}


	void MeshTransform(trimesh::TriMesh* traget_meshes,  trimesh::TriMesh* font_mesh, int face_id,
		trimesh::vec3 location, trimesh::vec3 dir, std::vector<float>& word_location, std::vector<int>& mesh_vertex_sizes, std::vector<trimesh::vec3>& word_mesh_center,
		trimesh::vec3 up,bool is_surround, float angle)
	{		
		if (!is_surround)
		{
			font_mesh->need_bbox();
			float hz = font_mesh->bbox.max.z - font_mesh->bbox.min.z;
			trimesh::vec3 dirTo(0, 1, 0);
			trimesh::vec3 faceTo(0, 0, -1);
			float cos = trimesh::vec3(0, 0, 1).dot(dir);
			if (cos < 0.f)
			{
				cos = trimesh::vec3(0, 0, -1).dot(dir);
				trimesh::vec3 scale_dir = dir * cos;
				dirTo = trimesh::vec3(0, 0, -1) - scale_dir;
				if (std::abs(-cos -1.f)< 1e-3)
					dirTo = trimesh::vec3(0, 1, 0);
			}
			else
			{
				trimesh::vec3 scale_dir = dir * cos;
				dirTo = trimesh::vec3(0, 0, 1) - scale_dir;
				if (std::abs(cos-1.f) < 1e-3)
					dirTo = trimesh::vec3(0, 1, 0);
			}
			
			trimesh::xform xf = trimesh::xform::rot_into(up, dirTo);
			faceTo = xf * faceTo;
			trimesh::vec3 fn=trimesh::normalized(traget_meshes->trinorm(face_id));
			trimesh::xform rot = trimesh::xform::rot_into(faceTo,fn);

			trimesh::apply_xform(font_mesh, xf);
			trimesh::apply_xform(font_mesh, rot);
			//trimesh::apply_xform(font_mesh, trimesh::xform::rot_into(trimesh::vec3(0, 0, 1), dir));
			trimesh::vec3 trans = location + ((2.0f * hz) / 5.0f) * dir;
			trimesh::trans(font_mesh, trans);
		}
		else
		{
			trimesh::TriMesh* _copy_mesh = new trimesh::TriMesh();
			_copy_mesh = traget_meshes;			

			//trimesh::xform xf = trimesh::xform::rot_into(trimesh::vec3(0,0,1), trimesh::vec3(0, 0, 1));
			trimesh::xform rot = trimesh::xform::rot_into(dir,trimesh::vec3(0,-1,0));

			trimesh::apply_xform(_copy_mesh, rot);
			//up = trimesh::normalized(rot * up);
			dir = trimesh::normalized(rot * dir);
			location = rot * location;
			trimesh::xform xf = trimesh::xform::rot(angle, dir.x, dir.y, dir.z);
			float height = (xf * location).z;
			
			_copy_mesh->need_across_edge();
			std::vector<std::pair<int, std::pair<float,float>>> left_clip_faces;
			std::vector<std::pair<int, std::pair<float,float>>> right_clip_faces;
			std::vector<std::pair<trimesh::vec3, trimesh::vec3>> left_face_corss_point;
			std::vector<std::pair<trimesh::vec3, trimesh::vec3>> right_face_corss_point;
			float left_len = 0;
			float right_len = 0;
			std::vector<int> face_marks(_copy_mesh->faces.size(),false);
			std::vector<trimesh::vec3> corss_points(2);
			int nn = 0;
			std::queue<int> que;
			que.push(face_id);
			while (!que.empty())
			{
				int f = que.front();
				que.pop();
				face_marks[f] = true;				
				
				for (int vi = 0; vi < 3; vi++)
				{
					int v = _copy_mesh->faces[f].at(vi);
					int v_n = _copy_mesh->faces[f].at((vi + 1) % 3);
					trimesh::vec3 dir = _copy_mesh->vertices[v] - _copy_mesh->vertices[v_n];
					float h = std::abs(_copy_mesh->vertices[v].z - _copy_mesh->vertices[v_n].z);
					if (h < 1e-6)
						continue;
					float scale = (std::abs(height - _copy_mesh->vertices[v_n].z)) / h * 1.0f;
					trimesh::vec3 new_position = _copy_mesh->vertices[v_n] + dir * scale;
					if (nn<2)
					{
						//0 --left    1 --right
						if (new_position.x < location.x)
						{
							corss_points[0] = new_position;
							float d= trimesh::distance(new_position, location);
							left_clip_faces.push_back(std::make_pair(f, std::make_pair(left_len, left_len + d)));
							left_face_corss_point.push_back(std::make_pair(location, new_position));
							left_len += d;
						}
						else
						{
							corss_points[1] = new_position;
							float d = trimesh::distance(new_position, location);
							right_clip_faces.push_back(std::make_pair(f, std::make_pair(right_len, right_len + d)));
							right_face_corss_point.push_back(std::make_pair(location, new_position));
							right_len += d;
						}
						nn++;
					}
					else
					{					
						float d1 = trimesh::distance(corss_points[0], new_position);
						float d2 = trimesh::distance(corss_points[1], new_position);
						if (d1 < d2)
						{
							if (d1 < 1e-6)
								continue;
							left_clip_faces.push_back(std::make_pair(f, std::make_pair(left_len,left_len+d1)));
							left_face_corss_point.push_back(std::make_pair(corss_points[0], new_position));
							corss_points[0] = new_position;		
							left_len += d1;
						}
						else
						{
							if (d2 < 1e-6)
								continue;
							right_clip_faces.push_back(std::make_pair(f, std::make_pair(right_len, right_len + d2)));
							right_face_corss_point.push_back(std::make_pair(corss_points[1], new_position));
							corss_points[1] = new_position;
							right_len += d2;
						}
					}
				}
				for (int fi = 0; fi < 3; fi++)
				{
					int ff = _copy_mesh->across_edge[f][fi];
					
					if (ff==-1||face_marks[ff])
						continue;
					//bool pass = false;				
					
					float z0 = _copy_mesh->vertices[_copy_mesh->faces[ff][0]].z;
					float z1 = _copy_mesh->vertices[_copy_mesh->faces[ff][1]].z;
					float z2 = _copy_mesh->vertices[_copy_mesh->faces[ff][2]].z;
					float max_z = std::max({ z0,z1,z2 });
					float min_z = std::min({ z0,z1,z2 });
					if (height >= min_z && height <= max_z)
						que.push(ff);
				}					
			}

			//trimesh::TriMesh* pointmesh = new trimesh::TriMesh();

			//for (int mi = 0; mi < 1; mi++)
			for (int mi = 0; mi < word_location.size(); mi++)
			{
				float wl=word_location[mi];
				float abs_wl = std::abs(wl);
				int f = -1;
				trimesh::vec3 new_wordpoint(0,0,0);
				if (wl <= 0)
				{
					for (int li = 0; li < left_clip_faces.size(); li++)
					{
						if (abs_wl<left_clip_faces[li].second.second && abs_wl>left_clip_faces[li].second.first)
						{
							f = left_clip_faces[li].first;
							float vv_sacle = (abs_wl - left_clip_faces[li].second.first) / (left_clip_faces[li].second.second - left_clip_faces[li].second.first);
							trimesh::vec3 vv_dir = left_face_corss_point[li].second-left_face_corss_point[li].first;
							new_wordpoint = left_face_corss_point[li].first + vv_sacle * vv_dir;
							break;
						}
					}
				}
				else
				{
					for (int ri = 0; ri < right_clip_faces.size(); ri++)
					{
						if (abs_wl<right_clip_faces[ri].second.second && abs_wl>right_clip_faces[ri].second.first)
						{
							f = right_clip_faces[ri].first;
							float vv_sacle = (abs_wl - right_clip_faces[ri].second.first) / (right_clip_faces[ri].second.second - right_clip_faces[ri].second.first);
							trimesh::vec3 vv_dir = right_face_corss_point[ri].second - right_face_corss_point[ri].first;
							new_wordpoint = right_face_corss_point[ri].first + vv_sacle * vv_dir;
							break;
						}
					}
				}
				

				if (f != -1)
				{
					//back zero vec
					int b = mi + 1;
					for (int vi = mesh_vertex_sizes[b - 1]; vi < mesh_vertex_sizes[b]; vi++)
					{
						font_mesh->vertices[vi] -= word_mesh_center[mi];
					}
					
					//font_mesh->write("fontmesh.ply");

					trimesh::vec3 f_n = trimesh::normalized(traget_meshes->trinorm(f));
					trimesh::xform xf_f = trimesh::xform::rot_into(up,trimesh::vec3(0, 0, 1));
					trimesh::vec3 new_dir = xf_f*trimesh::vec3(0, 0, -1);
					trimesh::xform face_new_dir = trimesh::xform::rot_into(new_dir, f_n);

					
					for (int vi = mesh_vertex_sizes[b-1]; vi < mesh_vertex_sizes[b]; vi++)
					{
						trimesh::vec3 new_point = face_new_dir*xf_f * font_mesh->vertices[vi];
						//font_mesh->vertices[vi] = new_point;
					}
					//trimesh::vec3 new_word_bbx_center = xf_f * word_mesh_center[mi];
					//trimesh::vec3 word_trans = new_wordpoint - word_mesh_center[mi];

					//pointmesh->vertices.push_back(new_wordpoint);
					for (int vi = mesh_vertex_sizes[b-1]; vi < mesh_vertex_sizes[b]; vi++)
					{
						font_mesh->vertices[vi] += new_wordpoint;
					}
					word_mesh_center[mi] = new_wordpoint;
					//pointmesh->vertices.push_back(new_wordpoint);
				}
			}

			trimesh::xform rott=trimesh::inv(rot);
			trimesh::apply_xform(font_mesh, rott);

			//pointmesh->write("pointmesh.ply");
			//_copy_mesh->write("copymesh.ply");
		}


		//font_mesh->write("fontmesh.ply");
		//traget_meshes->write("targetmesh.ply");
	}




	FontMesh::FontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter, float height,
		trimesh::vec3 face_to, trimesh::vec3 up) :Height(height)
	{
		_return_mesh = new trimesh::TriMesh();
		for (int li = 0; li < letter.size(); li++)
		{
			MMeshT mt(5000, 10000);
			mt.set_VFadjacent(true);
			std::vector<std::vector<trimesh::vec2>> totalpoly = letter[li];
			trimesh::vec2 wordbbx_min(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
			trimesh::vec2 wordbbx_max(std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
			for (int pi = 0; pi < totalpoly.size(); pi++)
				for (int ppi = 0; ppi < totalpoly[pi].size(); ppi++)
				{
					if (totalpoly[pi][ppi].x < wordbbx_min.x)
						wordbbx_min.x = totalpoly[pi][ppi].x;
					if (totalpoly[pi][ppi].x > wordbbx_max.x)
						wordbbx_max.x = totalpoly[pi][ppi].x;

					if (totalpoly[pi][ppi].y < wordbbx_min.y)
						wordbbx_min.y = totalpoly[pi][ppi].y;
					if (totalpoly[pi][ppi].x > wordbbx_max.x)
						wordbbx_max.y = totalpoly[pi][ppi].y;
				}
			mt.appendVertex(trimesh::point(wordbbx_min.x - 0.2f, wordbbx_max.y + 0.2f, 0));
			mt.appendVertex(trimesh::point(wordbbx_min.x - 0.2f, wordbbx_min.y - 0.2f, 0));
			mt.appendVertex(trimesh::point(wordbbx_max.x + 0.2f, wordbbx_min.y - 0.2f, 0));
			mt.appendVertex(trimesh::point(wordbbx_max.x + 0.2f, wordbbx_max.y + 0.2f, 0));
			mt.appendFace(0, 1, 2);
			mt.appendFace(0, 2, 3);
			std::vector<int> faceindex;
			for (int i = 0; i < mt.faces.size(); i++)
				faceindex.push_back(i);

			setMark(totalpoly);
			embedingAndCutting(&mt, totalpoly, faceindex);
			faceindex.clear();
			for (int i = 0; i < mt.faces.size(); i++)
				faceindex.push_back(i);
			std::vector<int> outfacesIndex;
			std::vector<std::vector<std::vector<trimesh::vec2>>> word_letter = { totalpoly };
			polygonInnerFaces(&mt, word_letter, faceindex, outfacesIndex);

			std::vector<int> least;
			std::set_difference(faceindex.begin(), faceindex.end(), outfacesIndex.begin(), outfacesIndex.end(), std::inserter(least, least.begin()));
			for (int fi : least)
				mt.deleteFace(fi);

			for (MMeshVertex& v : mt.vertices)
				if (v.connected_face.empty())
					mt.deleteVertex(v);
			mt.shrinkMesh();
			mt.init_halfedge();
			mt.set_HalfEdge(false);
			for (MMeshVertex& v : mt.vertices)
				v.p.z += Height;
			int face_size = mt.faces.size();
			std::vector<int> compara_table(mt.vertices.size());
			for (int fi = 0; fi < face_size; fi++)
			{
				MMeshFace& f = mt.faces[fi];
				if (f.IsD()) continue;
				int v0, v1, v2;
				if (!f.V0(0)->IsV())
				{
					mt.appendVertex(trimesh::point(f.V0(0)->p.x, f.V0(0)->p.y, 0));
					f.V0(0)->SetV();
					compara_table[f.V0(0)->index] = mt.vertices.size() - 1;
					v0 = mt.vertices.size() - 1;
				}
				else
				{
					v0 = compara_table[f.V0(0)->index];
				}

				if (!f.V1(0)->IsV())
				{
					mt.appendVertex(trimesh::point(f.V1(0)->p.x, f.V1(0)->p.y, 0));
					f.V1(0)->SetV();
					compara_table[f.V1(0)->index] = mt.vertices.size() - 1;
					v1 = mt.vertices.size() - 1;
				}
				else
				{
					v1 = compara_table[f.V1(0)->index];
				}

				if (!f.V2(0)->IsV())
				{
					mt.appendVertex(trimesh::point(f.V2(0)->p.x, f.V2(0)->p.y, 0));
					f.V2(0)->SetV();
					compara_table[f.V2(0)->index] = mt.vertices.size() - 1;
					v2 = mt.vertices.size() - 1;
				}
				else
				{
					v2 = compara_table[f.V2(0)->index];
				}
				mt.appendFace(v2, v1, v0);
				MMeshHalfEdge* p = f.f_mhe;
				int n = 3;
				while (1)
				{
					if (p->opposite == nullptr)
					{
						int i1 = p->edge_vertex.second->index;
						int i2 = p->edge_vertex.first->index;
						int i3 = 0;
						int i4 = 0;
						if (n == 3)
						{
							i3 = v0;
							i4 = v1;
						}
						else if (n == 2)
						{
							i3 = v1;
							i4 = v2;
						}
						else if (n == 1)
						{
							i3 = v2;
							i4 = v0;
						}

						mt.appendFace(i1, i2, i3);
						mt.appendFace(i3, i4, i1);
					}
					n--;
					if (p->next == f.f_mhe)
						break;
					p = p->next;
				}
			}
			trimesh::TriMesh* _word_mesh = new trimesh::TriMesh();
			mt.mmesh2trimesh(_word_mesh);
			_word_mesh->need_bbox();
			font_meshs.push_back(_word_mesh);
			FaceTo.push_back(face_to);
			Up.push_back(up);			

			word_absolute_location.push_back(_word_mesh->bbox.center());
			trimesh::TriMesh* _word_mesh_t=new trimesh::TriMesh();
			mt.mmesh2trimesh(_word_mesh_t);
			init_font_meshs.push_back(_word_mesh_t);
			bbx += _word_mesh->bbox;
			
		}
		for (int li = 0; li < letter.size(); li++)
		{			
			word_relative_location.push_back(word_absolute_location[li]- bbx.center());
		}
	}

	FontMesh::FontMesh(const FontMesh& other)
	{
		font_meshs = other.font_meshs;
		init_font_meshs =other.init_font_meshs;
		word_relative_location =other.word_relative_location;
		word_absolute_location = other.word_absolute_location;
		FaceTo=other.FaceTo;
		Up=other.Up;
		_return_mesh = other._return_mesh;
	}


	FontMesh::~FontMesh()
	{
		for (auto m : font_meshs)
		{
			m->clear(); m = nullptr;
		}
		font_meshs.clear();
		for (auto im : init_font_meshs)
		{
			im->clear(); im = nullptr;
		}
		init_font_meshs.clear();
		word_relative_location.clear();
		FaceTo.clear();
		Up.clear();
		_return_mesh->clear();
		_return_mesh = nullptr;
	}


	trimesh::TriMesh* FontMesh::getFontMesh()
	{
		if (is_change)
		{
			_return_mesh->clear();
			for (int wi = 0; wi < font_meshs.size(); wi++)
			{
				int vsize = _return_mesh->vertices.size();
				for (trimesh::point v : font_meshs[wi]->vertices)
					_return_mesh->vertices.push_back(v);
				for (trimesh::TriMesh::Face f : font_meshs[wi]->faces)
				{
					_return_mesh->faces.push_back(trimesh::TriMesh::Face(vsize + f[0], vsize + f[1], vsize + f[2]));
				}
			}
		}
		return _return_mesh;			
	}


	void FontMesh::FontTransform(trimesh::TriMesh* traget_meshes, int face_id, trimesh::vec3 location, bool is_surround)
	{
		trimesh::vec3 fn = trimesh::normalized(traget_meshes->trinorm(face_id));
		is_change = true;
		if (!is_surround)
		{
			_return_mesh->clear();			
			for (int wi = 0; wi < init_font_meshs.size(); wi++)
			{
				int vsize = _return_mesh->vertices.size();
				for (trimesh::point v : init_font_meshs[wi]->vertices)
					_return_mesh->vertices.push_back(v);
				for (trimesh::TriMesh::Face f : init_font_meshs[wi]->faces)
				{
					_return_mesh->faces.push_back(trimesh::TriMesh::Face(vsize + f[0], vsize + f[1], vsize + f[2]));
				}
			}
					
			trimesh::trans(_return_mesh, -bbx.center());
			//---rot---	
			trimesh::xform rot=trimesh::xform::rot_into(trimesh::vec3(0,0,-1),fn);
			trimesh::apply_xform(_return_mesh, rot);		
			trimesh::vec3 new_up = rot * trimesh::vec3(0, -1, 0);
			trimesh::normalize(new_up);
			float zcos = new_up.dot(trimesh::vec3(0,0,1));
			if (std::abs(zcos) < 1e-2)
			{
				float y_angle = std::acos(new_up.dot(trimesh::vec3(0, 1, 0)));
				trimesh::rot(_return_mesh,y_angle,fn);
			}
			else
			{
				float z_angle = std::acos(zcos);
				trimesh::rot(_return_mesh, z_angle, fn);
			}

			trimesh::trans(_return_mesh, location);
			is_change = false;						
		}
		else {
			trimesh::TriMesh* _copy_mesh = new trimesh::TriMesh();
			_copy_mesh = traget_meshes;
			_copy_mesh->need_bbox();
			trimesh::trans(_copy_mesh, -_copy_mesh->bbox.center());
			trimesh::xform xf = trimesh::xform::rot_into(fn,trimesh::vec3(0,-1,0));
			trimesh::apply_xform(_copy_mesh, xf);
			float height = (xf * location).z;


			_copy_mesh->need_across_edge();
			std::vector<int> face_marks(_copy_mesh->faces.size(), false);
			std::vector<std::pair<int, std::pair<float, float>>> face_line_len;
			std::vector<std::pair<trimesh::vec3, trimesh::vec3>> face_corss_point;
			std::queue<int> que;
			trimesh::vec3 front_cross=location;			
			float len = 0;
			bool is_frist=true;
			bool is_frist2 = true;
			trimesh::ivec2 frist_mark_v;
			que.push(face_id);
			while (!que.empty())
			{
				int f = que.front();
				que.pop();
				face_marks[f] = true;

				std::vector<trimesh::vec3> double_cross;
				
				for (int vi = 0; vi < 3; vi++)
				{
					int v = _copy_mesh->faces[f].at(vi);
					int v_n = _copy_mesh->faces[f].at((vi + 1) % 3);
					trimesh::vec3 dir = _copy_mesh->vertices[v] - _copy_mesh->vertices[v_n];
					float h = std::abs(_copy_mesh->vertices[v].z - _copy_mesh->vertices[v_n].z);
					if (h < 1e-4)
						continue;
					float scale = (std::abs(height - _copy_mesh->vertices[v_n].z)) / h * 1.0f;
					trimesh::vec3 new_position = _copy_mesh->vertices[v_n] + dir * scale;
					if (is_frist)
					{
						if (new_position.x > location.x)
						{
							face_corss_point.push_back(std::make_pair(front_cross,new_position));
							float d = trimesh::distance(front_cross, new_position);
							len += d;
							face_line_len.push_back(std::make_pair(f,std::make_pair(0,d)));
							front_cross = new_position;
							is_frist = false;		
							frist_mark_v = trimesh::ivec2(v,v_n);
						}
						/*else
						{
							back_cross = new_position;
						}*/
					}
					else
					{
						double_cross.push_back(new_position);
					}
				}
				if (double_cross.size()==2)
				{
					float d1 = trimesh::distance(double_cross[0], front_cross);
					float d2 = trimesh::distance(double_cross[1], front_cross);
					if (d1 > d2)
					{
						face_corss_point.push_back(std::make_pair(double_cross[1], double_cross[0]));
						front_cross = double_cross[0];
						face_line_len.push_back(std::make_pair(f, std::make_pair(len,len + d1)));
						len += d1;
					}
					else
					{
						face_corss_point.push_back(std::make_pair(double_cross[0], double_cross[1]));
						front_cross = double_cross[1];
						face_line_len.push_back(std::make_pair(f, std::make_pair(len, len + d2)));
						len += d2;
					}
				}

				for (int fi = 0; fi < 3; fi++)
				{
					int ff = _copy_mesh->across_edge[f][fi];
					if (ff == -1 || face_marks[ff])
						continue;
					if (is_frist2)
					{
						int pass_n = 0;
						for (int fvi = 0; fvi < 3; fvi++)
						{
							int fv = _copy_mesh->faces[ff][fvi];
							if (fv == frist_mark_v[0] || fv == frist_mark_v[1])
								pass_n++;
						}
						if (pass_n!=2)
							continue;
						else
						{
							is_frist2 = false;
						}
					}
					float z0 = _copy_mesh->vertices[_copy_mesh->faces[ff][0]].z;
					float z1 = _copy_mesh->vertices[_copy_mesh->faces[ff][1]].z;
					float z2 = _copy_mesh->vertices[_copy_mesh->faces[ff][2]].z;
					float max_z = std::max({ z0,z1,z2 });
					float min_z = std::min({ z0,z1,z2 });
					if (height >= min_z && height <= max_z)
						que.push(ff);
				}

			}
			float last_d = trimesh::distance(front_cross, location);
			face_corss_point.push_back(std::make_pair(front_cross, location));
			face_line_len.push_back(std::make_pair(face_id,std::make_pair(len, len+last_d)));
			len += last_d;
								
			font_meshs.clear();
			//trimesh::TriMesh* points = new trimesh::TriMesh();
			for (int wi = 0; wi < init_font_meshs.size(); wi++)
			{
				float loc = word_relative_location[wi].x;
				int sel_f = -1;
				trimesh::vec3 word_new_location(0,0,0);
				if (loc < 0)
				{
					loc = loc + len;
				}				
				for (int fl_i = 0; fl_i < face_line_len.size(); fl_i++)
				{
					if (loc<face_line_len[fl_i].second.second && loc>face_line_len[fl_i].second.first)
					{
						sel_f = face_line_len[fl_i].first;
						float vv_scale = (loc - face_line_len[fl_i].second.first) / (face_line_len[fl_i].second.second - face_line_len[fl_i].second.first);
						trimesh::vec3 dirTo = face_corss_point[fl_i].second - face_corss_point[fl_i].first;
						word_new_location = face_corss_point[fl_i].first + vv_scale * dirTo;
						break;
					}
				}
				//points->vertices.push_back(word_new_location);
				//trimesh::vec3 sel_fn= trimesh::normalized(traget_meshes->trinorm(sel_f));
				trimesh::TriMesh* new_font_mesh=init_font_meshs[wi];
				new_font_mesh->need_bbox();
				word_new_location-=new_font_mesh->bbox.center();
				trimesh::trans(new_font_mesh, word_new_location);
				//trimesh::xform xxf = trimesh::inv(xf);
				//trimesh::apply_xform(new_font_mesh, xxf);
				font_meshs.push_back(new_font_mesh);
			}
			//points->write("points.ply");
			//_copy_mesh->write("copymesh.ply");
		}
	}

}