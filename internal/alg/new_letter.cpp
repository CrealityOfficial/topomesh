#include "new_letter.h"
#include "letter.h"
#include "earclipping.h"
#include "trimesh2/TriMesh_algo.h"
#include <Eigen/Dense>
#include <unordered_map>
#include "msbase/mesh/dumplicate.h"

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

	

	FontMesh::FontMesh(float height, float depth, float angle)
	{
		_return_mesh = new trimesh::TriMesh();
		_return_surround_mesh = new trimesh::TriMesh();
		m_config.height = height;
		m_config.distance = depth;
		m_config.angle = angle; 
	}

	FontMesh::FontMesh(const FontMesh& other)
	{		
		init_font_meshs =other.init_font_meshs;
		word_init_location =other.word_init_location;
		word_absolute_location = other.word_absolute_location;
		FaceTo=other.FaceTo;
		word_FaceTo = other.word_FaceTo;
		word_Up = other.word_Up;
		Up = other.Up;
		click_location = other.click_location;
		_return_surround_mesh = new trimesh::TriMesh;
		_return_surround_mesh = other._return_surround_mesh;
		_return_mesh = new trimesh::TriMesh;
		_return_mesh = other._return_mesh;
		_m_xform = other._m_xform;
		m_config = other.m_config;
		
	}


	FontMesh::~FontMesh()
	{		
		for (auto im : init_font_meshs)
		{
			im->clear(); im = nullptr;
		}
		init_font_meshs.clear();
		word_init_location.clear();
		word_FaceTo.clear();
		word_Up.clear();
		_return_mesh->clear();
		_return_mesh = nullptr;
		_return_surround_mesh->clear();
		_return_surround_mesh = nullptr;
	}


	void FontMesh::setText(const std::string& text)
	{
		m_config.text = text;
	}

	std::string FontMesh::text() const
	{
		return m_config.text;
	}

	float FontMesh::angle()
	{
		return m_config.angle;
	}

	float FontMesh::height()
	{
		return m_config.height;
	}

	trimesh::vec3 FontMesh::currentFaceTo()
	{
		return current_faceto;
	}

	void FontMesh::setState(int state)
	{
		if (state != m_config.state)
		{
			m_config.state = state;
			is_change_state = true;
		}
	}

	float FontMesh::depth()
	{
		return m_config.distance;
	}

	trimesh::TriMesh* FontMesh::getFontMesh()
	{
		if (!m_config.state)
		{				
			trimesh::TriMesh* result = new trimesh::TriMesh;
			*result = *_return_mesh;		
			
			trimesh::xform xf = trimesh::xform::rot_into(trimesh::vec3(0,0,1), FaceTo.second);
			trimesh::vec3 dir_up = xf * trimesh::vec3(0, 1, 0);
			trimesh::xform xxf = trimesh::xform::rot_into(dir_up, Up.second);
			trimesh::apply_xform(result, xxf*xf);
			
			trimesh::xform rot_xf = trimesh::xform::rot(m_config.angle *M_PI*1.0f, FaceTo.second);
			trimesh::apply_xform(result, rot_xf);
			float half = m_config.height / 2.0f;
			trimesh::vec3 dirto = (m_config.distance + half) * FaceTo.second;
			trimesh::trans(result, click_location + dirto);
			//result->write("result.ply");

			return result;
		}
		else {
			_return_surround_mesh->clear();
			int v_size = 0;
			//trimesh::TriMesh* points = new trimesh::TriMesh();
			for (int mi = 0; mi < init_font_meshs.size(); mi++)
			{
				float half = m_config.height / 2.0f;
				trimesh::vec3 dirto = (m_config.distance + half) * word_FaceTo[mi];
				trimesh::vec3 transTo = word_absolute_location[mi]+ dirto;
				//points->vertices.push_back(word_absolute_location[mi]);
				trimesh::xform face_xf = trimesh::xform::rot_into(trimesh::vec3(0,0,1), word_FaceTo[mi]);
				trimesh::vec3 new_up = face_xf * trimesh::vec3(0, 1, 0);	

				//trimesh::xform rot_xf = trimesh::xform::rot(m_config.angle * M_PI * 1.0f, word_FaceTo[mi]);
				trimesh::vec3 new_word_up = /*rot_xf * */word_Up[mi];

				trimesh::xform angle_xf = trimesh::xform::rot_into(new_up, new_word_up);
				for (int vi = 0; vi < init_font_meshs[mi]->vertices.size(); vi++)
				{
					trimesh::vec3 v = init_font_meshs[mi]->vertices[vi];
					v = angle_xf * face_xf * v;
					 _return_surround_mesh->vertices.push_back(v + transTo);
				}
				for (int fi = 0; fi < init_font_meshs[mi]->faces.size(); fi++)
				{
					_return_surround_mesh->faces.push_back(trimesh::ivec3(init_font_meshs[mi]->faces[fi][0] + v_size, init_font_meshs[mi]->faces[fi][1] + v_size,
						init_font_meshs[mi]->faces[fi][2] + v_size));
				}
				v_size += init_font_meshs[mi]->vertices.size();				
			}
			trimesh::TriMesh* result = new trimesh::TriMesh;
			*result = *_return_surround_mesh;
			//result->write("surround.ply");
			//points->write("surrend_point.ply");
			return result;
		}
		
	}

	bool FontMesh::calRelativeCoord(trimesh::TriMesh* traget_meshes, int face_id, trimesh::vec3 location)
	{
		
		int v0 = traget_meshes->faces[face_id][0];
		int v1 = traget_meshes->faces[face_id][1];
		int v2 = traget_meshes->faces[face_id][2];
		trimesh::vec3 p1 = traget_meshes->vertices[v1] - traget_meshes->vertices[v0];
		trimesh::vec3 p2 = traget_meshes->vertices[v2] - traget_meshes->vertices[v0];

		Eigen::Matrix2f e;
		e << p1.x, p2.x, p1.y, p2.y;
		Eigen::Vector2f b = { location.x ,location.y };
		Eigen::Vector2f x = e.fullPivLu().solve(b);

		relative_coord = trimesh::vec2(x.x(),x.y());
		if ((x.x() + x.y()) > 1.0f || x.x() < 0 || x.y() < 0)
			return true;
		return false;
	}


	trimesh::vec3 FontMesh::getRelativeCoord(trimesh::TriMesh* traget_meshes)
	{
		int v0 = traget_meshes->faces[sel_faceid][0];
		int v1 = traget_meshes->faces[sel_faceid][1];
		int v2 = traget_meshes->faces[sel_faceid][2];
		trimesh::vec3 p1 = traget_meshes->vertices[v1] - traget_meshes->vertices[v0];
		trimesh::vec3 p2 = traget_meshes->vertices[v2] - traget_meshes->vertices[v0];
		click_location = traget_meshes->vertices[v0]+p1 * relative_coord.x + p2 * relative_coord.y;
		return click_location;
	}

	bool FontMesh::checkMistakes(trimesh::TriMesh* traget_meshes, int face_id, trimesh::vec3 location)
	{
		return calRelativeCoord(traget_meshes, face_id, location);
	}


	bool FontMesh::FontTransform(trimesh::TriMesh* traget_meshes, int face_id, trimesh::vec3 location, bool is_surround)
	{
		bool is_mistake = false;
		if (calRelativeCoord(traget_meshes, face_id, location))
		{
			/*traget_meshes->write("targetmesh.ply");
			trimesh::TriMesh* lines = new trimesh::TriMesh();
			lines->vertices.push_back(location);
			lines->write("lines.ply");*/
			is_mistake = true;
		}
		//if(1)
		if (!is_change_state || sel_faceid == -1)
		{
			click_location = location;
			sel_faceid = face_id;
		}	
		m_config.state = is_surround;
		is_change_state = false;
		trimesh::vec3 fn = trimesh::normalized(traget_meshes->trinorm(sel_faceid));
		current_faceto = fn;

		/*sel_faceid = face_id;
		int va = traget_meshes->faces[sel_faceid].at(0);
		int vb = traget_meshes->faces[sel_faceid].at(1);
		int vc = traget_meshes->faces[sel_faceid].at(2);
		click_location = (traget_meshes->vertices[va] + traget_meshes->vertices[vb] + traget_meshes->vertices[vc]) / 3.0f;
		trimesh::vec3 fn = trimesh::normalized(traget_meshes->trinorm(sel_faceid));
		current_faceto = fn;*/

		if (!m_config.state)
		{								
			trimesh::vec3 ori_faceTo = FaceTo.second;
			trimesh::xform faceto_xf = trimesh::xform::rot_into(ori_faceTo,fn);
			FaceTo.first = FaceTo.second;
			FaceTo.second = fn;					
			trimesh::vec3 new_up = faceto_xf * Up.second;
			trimesh::normalize(new_up);
			Up.first = new_up;
			float zcos = trimesh::vec3(0,0,1).dot(fn);
			trimesh::xform angle_xf;
			trimesh::vec3 up_dirto;
			if (std::abs(zcos-1.0f) < 1e-2)
			{
				up_dirto = trimesh::vec3(0, 1, 0);
			}
			else
			{	
				up_dirto = trimesh::vec3(0, 0, 1) + zcos * -fn;
			}			
			Up.second = trimesh::normalized(up_dirto);		
			return false;
		}
		else {
#if 0
			trimesh::TriMesh* _copy_mesh = new trimesh::TriMesh;
			*_copy_mesh = *traget_meshes;
			_copy_mesh->clear_bbox();
			_copy_mesh->need_bbox();
			//_copy_mesh->write("before_copymesh.ply");
			trimesh::vec3 trans_bbx_center = _copy_mesh->bbox.center();
			//msbase::mergeNearPoints(_copy_mesh, nullptr, 1e-4f);
			trimesh::trans(_copy_mesh, -_copy_mesh->bbox.center());
			trimesh::xform xf = trimesh::xform::rot_into(fn, trimesh::vec3(0, -1, 0));

			trimesh::vec3 new_face_to = xf * FaceTo.second;
			trimesh::xform rot_xf = trimesh::xform::rot(-m_config.angle * M_PI * 1.0f, new_face_to);
			trimesh::xform r_xxf = rot_xf * xf;
			trimesh::apply_xform(_copy_mesh, r_xxf);
			trimesh::vec3 _copy_location = r_xxf * (click_location - trans_bbx_center);
			float height = (r_xxf * click_location).z;

			/*trimesh::TriMesh* locationmesh = new trimesh::TriMesh();
			locationmesh->vertices.push_back(_copy_location);
			locationmesh->write("locationmesh.ply");
			_copy_mesh->write("_copymesh.ply");*/

			_copy_mesh->need_across_edge();
			std::vector<int> face_marks(_copy_mesh->faces.size(), false);
			std::vector<std::pair<int, std::pair<float, float>>> face_line_len;
			std::vector<std::pair<trimesh::vec3, trimesh::vec3>> face_corss_point;
			std::queue<int> que;
			trimesh::vec3 front_cross = _copy_location;

			//std::vector<int> inputfaces;

			float len = 0;
			bool is_frist = true;
			bool is_frist2 = true;
			trimesh::ivec2 frist_mark_v;
			trimesh::ivec2 next_mark_v;
			que.push(sel_faceid);
			while (!que.empty())
			{
				int f = que.front();
				que.pop();
				face_marks[f] = true;
				std::vector<trimesh::vec3> double_cross;
				std::vector<trimesh::ivec2> v_mark;
				for (int vi = 0; vi < 3; vi++)
				{
					int v = _copy_mesh->faces[f].at(vi);
					int v_n = _copy_mesh->faces[f].at((vi + 1) % 3);
					trimesh::vec3 dir = _copy_mesh->vertices[v] - _copy_mesh->vertices[v_n];
					float h = _copy_mesh->vertices[v].z - _copy_mesh->vertices[v_n].z;
					if (std::abs(h) < 1e-4)
						continue;
					float scale = (height - _copy_mesh->vertices[v_n].z) / h * 1.0f;
					if (scale > 1 || scale < 0)
						continue;
					trimesh::vec3 new_position = _copy_mesh->vertices[v_n] + dir * scale;
					if (is_frist)
					{
						if (new_position.x > _copy_location.x)
						{
							face_corss_point.push_back(std::make_pair(front_cross, new_position));
							float d = trimesh::distance(front_cross, new_position);
							len += d;
							face_line_len.push_back(std::make_pair(f, std::make_pair(0, d)));
							front_cross = new_position;
							is_frist = false;
							frist_mark_v = trimesh::ivec2(v, v_n);
						}
					}
					else
					{
						double_cross.push_back(new_position);
						v_mark.push_back(trimesh::ivec2(v, v_n));
					}
				}
				if (double_cross.size() == 2)
				{
					float d1 = trimesh::distance(double_cross[0], front_cross);
					float d2 = trimesh::distance(double_cross[1], front_cross);
					if (d1 > d2)
					{
						face_corss_point.push_back(std::make_pair(double_cross[1], double_cross[0]));
						front_cross = double_cross[0];
						face_line_len.push_back(std::make_pair(f, std::make_pair(len, len + d1)));
						len += d1;
						next_mark_v = v_mark[0];
					}
					else
					{
						face_corss_point.push_back(std::make_pair(double_cross[0], double_cross[1]));
						front_cross = double_cross[1];
						face_line_len.push_back(std::make_pair(f, std::make_pair(len, len + d2)));
						len += d2;
						next_mark_v = v_mark[1];
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
						if (pass_n != 2)
							continue;
						else
						{
							is_frist2 = false;
						}
					}
					else
					{
						int pass_c = 0;
						for (int fvi = 0; fvi < 3; fvi++)
						{
							int fv = _copy_mesh->faces[ff][fvi];
							if (fv == next_mark_v[0] || fv == next_mark_v[1])
								pass_c++;
						}
						if (pass_c != 2)
							continue;
					}

					float z0 = _copy_mesh->vertices[_copy_mesh->faces[ff][0]].z;
					float z1 = _copy_mesh->vertices[_copy_mesh->faces[ff][1]].z;
					float z2 = _copy_mesh->vertices[_copy_mesh->faces[ff][2]].z;
					float max_z = std::max({ z0,z1,z2 });
					float min_z = std::min({ z0,z1,z2 });
					if (height >= min_z && height <= max_z)
					{
						que.push(ff);
						//inputfaces.push_back(ff);
						break;
					}
				}
				v_mark.clear();
			}
			float last_d = trimesh::distance(front_cross, _copy_location);
			face_corss_point.push_back(std::make_pair(front_cross, _copy_location));
			face_line_len.push_back(std::make_pair(sel_faceid, std::make_pair(len, len + last_d)));
			len += last_d;

			/*trimesh::TriMesh* flines = new trimesh::TriMesh();
			for (int fci = 0; fci < face_corss_point.size(); fci++)
			{
				flines->vertices.push_back(face_corss_point[fci].first);
				flines->vertices.push_back(face_corss_point[fci].second);
			}
			flines->write("flines.ply");*/

			//trimesh::TriMesh* locationpoint = new trimesh::TriMesh();
			//locationpoint->vertices.push_back(location);
			trimesh::xform xxf = trimesh::inv(r_xxf);
			trimesh::apply_xform(_copy_mesh, xxf);
			float zaxis_cos = trimesh::vec3(0, 0, 1).dot(fn);
			trimesh::vec3 axis_to;
			trimesh::vec3 axis_up;
			if (std::abs(zaxis_cos - 1.f) < 1e-2)
			{
				axis_to = trimesh::vec3(0, 1, 0);
				axis_up = trimesh::vec3(0, 0, 1);
			}
			else
			{
				axis_to = trimesh::vec3(0, 0, 1);
				axis_up = trimesh::vec3(0, 1, 0);
			}

			for (int wi = 0; wi < init_font_meshs.size(); wi++)
			{
				float loc = word_init_location[wi].x - bbx.center().x;
				int sel_f = -1;
				trimesh::vec3 word_new_location(0, 0, 0);
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
				if (sel_f == -1)
					continue;
				trimesh::vec3 sel_fn = trimesh::normalized(traget_meshes->trinorm(sel_f));
				word_FaceTo[wi] = sel_fn;

				float axis_cos = axis_to.dot(sel_fn);
				trimesh::vec3 up_dirto;
				if (std::abs(axis_cos - 1.0f) < 1e-2)
				{
					up_dirto = axis_up;
				}
				else
				{
					up_dirto = axis_to + axis_cos * -sel_fn;
				}
				word_Up[wi] = /*xxf**/up_dirto;
				//locationpoint->vertices.push_back(word_new_location);
				word_new_location = xxf * word_new_location;
				word_absolute_location[wi] = word_new_location + trans_bbx_center;
				//locationpoint->vertices.push_back(word_new_location);

			}
			//locationpoint->write("locationpoint.ply");
			//points->write("points.ply");	
#elif 0
			trimesh::TriMesh* _copy_mesh = new trimesh::TriMesh;
			//traget_meshes->need_adjacentfaces();
			//traget_meshes->need_across_edge();
			*_copy_mesh = *traget_meshes;
			_copy_mesh->need_bbox();
			trimesh::vec3 trans_bbx_center = _copy_mesh->bbox.center();
			trimesh::trans(_copy_mesh, -trans_bbx_center);
			trimesh::xform xf = trimesh::xform::rot_into(fn, trimesh::vec3(0, -1, 0));

			trimesh::vec3 new_face_to = xf * current_faceto;
			trimesh::xform rot_xf = trimesh::xform::rot(-m_config.angle * M_PI * 1.0f, new_face_to);
			trimesh::xform r_xxf = rot_xf * xf;
			trimesh::apply_xform(_copy_mesh, r_xxf);
			trimesh::vec3 _copy_location = r_xxf * (click_location - trans_bbx_center);
			float height = _copy_location.z;

			/*_copy_mesh->write("_copymesh.ply");
			trimesh::TriMesh* point = new trimesh::TriMesh();
			point->vertices.push_back(_copy_location);
			point->write("points.ply");*/

			_copy_mesh->need_across_edge();
			msbase::mergeNearPoints(_copy_mesh, nullptr, 1e-4f);


			struct hash_func {
				size_t operator()(const trimesh::ivec2& v)const
				{
					return int(v.x * 100000000 + v.y);
				}
			};

			struct equal_ivec2 {
				bool operator()(const trimesh::ivec2& v1, const trimesh::ivec2& v2) const
				{
					return v1.x == v2.x && v1.y == v2.y;
				}
			};

			std::unordered_map<trimesh::ivec2, int, hash_func, equal_ivec2> edges;
			typedef typename std::unordered_map<trimesh::ivec2, int, hash_func, equal_ivec2>::value_type edge_value;
			std::vector<bool> del_faces(_copy_mesh->faces.size(), false);
			std::vector<int> center_faces;
			for (int fi = 0; fi < _copy_mesh->faces.size(); fi++)
			{
				int v0 = _copy_mesh->faces[fi].at(0);
				int v1 = _copy_mesh->faces[fi].at(1);
				int v2 = _copy_mesh->faces[fi].at(2);

				if ((_copy_mesh->vertices[v0].z < height && _copy_mesh->vertices[v1].z < height && _copy_mesh->vertices[v2].z < height) ||
					(_copy_mesh->vertices[v0].z > height && _copy_mesh->vertices[v1].z > height && _copy_mesh->vertices[v2].z > height))
					continue;

				center_faces.push_back(fi);
				for (int vi = 0; vi < 3; vi++)
				{
					int v = _copy_mesh->faces[fi].at(vi);
					int v_n = _copy_mesh->faces[fi].at((vi + 1) % 3);
					int min = std::min(v, v_n);
					int max = std::max(v, v_n);

					float h = _copy_mesh->vertices[v].z - _copy_mesh->vertices[v_n].z;
					if (std::abs(h) < 1e-5)
						continue;
					float scale = (height - _copy_mesh->vertices[v_n].z) / h * 1.0f;
					if (scale > 1 || scale < 0)
						continue;
					trimesh::vec3 dir = _copy_mesh->vertices[v] - _copy_mesh->vertices[v_n];
					trimesh::vec3 new_position = _copy_mesh->vertices[v_n] + dir * scale;

					auto it_e = edges.find(trimesh::ivec2(min, max));
					if (it_e != edges.end())
					{

					}
					else
					{
						_copy_mesh->vertices.push_back(new_position);
						int vsize = _copy_mesh->vertices.size() - 1;
						edges.insert(edge_value(trimesh::ivec2(min, max), vsize));
					}
				}
			}

			std::vector<int> begin_vertex(2);
			for (int i = 0; i < center_faces.size(); i++)
			{
				int f = center_faces[i];
				del_faces[f] = true;

				std::vector<int> n_v = { _copy_mesh->faces[f][0],_copy_mesh->faces[f][1],_copy_mesh->faces[f][2] };
				int n = 0;
				std::vector<bool> index(5, false);
				int sum_index = 0;
				for (int vi = 0; vi < 3; vi++)
				{
					int v = _copy_mesh->faces[f].at(vi);
					int v_n = _copy_mesh->faces[f].at((vi + 1) % 3);
					int min = std::min(v, v_n);
					int max = std::max(v, v_n);
					float h = _copy_mesh->vertices[v].z - _copy_mesh->vertices[v_n].z;
					if (std::abs(h) < 1e-5)
						continue;
					auto ite = edges.find(trimesh::ivec2(min, max));
					if (ite != edges.end())
					{
						n_v.insert(n_v.begin() + vi + n + 1, ite->second);
						index[vi + n + 1] = true;
						sum_index += (vi + n + 1);
						n++;
					}
				}


				for (int ni = 0; ni < n_v.size(); ni++)
				{
					if (index[ni])
					{
						int before_v = n_v[(ni - 1 + n_v.size()) % n_v.size()];
						int after_v = n_v[sum_index - ni];
						_copy_mesh->faces.push_back(trimesh::ivec3(before_v, n_v[ni], after_v));
						if (f == sel_faceid)
						{
							begin_vertex[0] = (n_v[ni]);
							begin_vertex[1] = (after_v);
						}
					}
					if (!index[ni] && !index[(ni + 1) % n_v.size()])
					{
						if (index[(ni - 1 + n_v.size()) % n_v.size()])
						{
							_copy_mesh->faces.push_back(trimesh::ivec3(n_v[(ni - 1 + n_v.size()) % n_v.size()], n_v[ni], n_v[(ni + 1) % n_v.size()]));
						}
					}

				}
			}

			int n = 0;
			int dsize = del_faces.size();
			for (int i = 0; i < _copy_mesh->faces.size(); i++)
			{
				if (i >= dsize)
					del_faces.push_back(false);
				else
				{
					if (del_faces[i])
						n++;
				}
			}
			trimesh::remove_faces(_copy_mesh, del_faces);



			std::vector<bool> delfaces(_copy_mesh->faces.size(), false);
			for (int fi = 0; fi < _copy_mesh->faces.size(); fi++)
			{
				int v0 = _copy_mesh->faces[fi].at(0);
				int v1 = _copy_mesh->faces[fi].at(1);
				int v2 = _copy_mesh->faces[fi].at(2);
				trimesh::vec3 c = (_copy_mesh->vertices[v0] + _copy_mesh->vertices[v1] + _copy_mesh->vertices[v2]) / 3.0f;
				if (c.z > height)
				{
					delfaces[fi] = true;
				}
			}

			trimesh::remove_faces(_copy_mesh, delfaces);


			//_copy_mesh->write("crossmesh.ply");

			_copy_mesh->need_neighbors();
			_copy_mesh->need_adjacentfaces();

			std::vector<bool> mark_vertex(_copy_mesh->vertices.size(), false);
			std::vector<bool> vis_vertex(_copy_mesh->vertices.size(), false);

			//trimesh::TriMesh* lines1 = new trimesh::TriMesh();
			for (int vi = 0; vi < _copy_mesh->vertices.size(); vi++)
			{
				if (std::abs(_copy_mesh->vertices[vi].z - height) < 1e-4 && _copy_mesh->adjacentfaces[vi].size() != _copy_mesh->neighbors[vi].size())
				{
					mark_vertex[vi] = true;
					//lines1->vertices.push_back(_copy_mesh->vertices[vi]);
				}
			}
			//lines1->write("lines1.ply");
			std::vector<std::pair<int, std::pair<float, float>>> face_line_len;
			std::vector<std::pair<trimesh::vec3, trimesh::vec3>> face_corss_point;

			int frist_vertex = -1;
			float len = 0;

			vis_vertex[begin_vertex[0]] = true;
			vis_vertex[begin_vertex[1]] = true;
			if (_copy_mesh->vertices[begin_vertex[0]].x > _copy_mesh->vertices[begin_vertex[1]].x)
			{
				frist_vertex = begin_vertex[0];
				float dist = trimesh::distance(_copy_mesh->vertices[begin_vertex[0]], _copy_location);
				len += dist;
				face_line_len.push_back(std::make_pair(_copy_mesh->adjacentfaces[frist_vertex][0], std::make_pair(0, len)));
				face_corss_point.push_back(std::make_pair(_copy_location, _copy_mesh->vertices[begin_vertex[0]]));
			}
			else
			{
				frist_vertex = begin_vertex[1];
				float dist = trimesh::distance(_copy_mesh->vertices[begin_vertex[1]], _copy_location);
				len += dist;
				face_line_len.push_back(std::make_pair(_copy_mesh->adjacentfaces[frist_vertex][0], std::make_pair(0, len)));
				face_corss_point.push_back(std::make_pair(_copy_location, _copy_mesh->vertices[begin_vertex[1]]));
			}



			trimesh::TriMesh* lines = new trimesh::TriMesh();

			while (true)
			{
				int before_vertex = frist_vertex;
				vis_vertex[frist_vertex] = true;

				for (int vi = 0; vi < _copy_mesh->neighbors[frist_vertex].size(); vi++)
				{
					int v = _copy_mesh->neighbors[frist_vertex][vi];
					if (mark_vertex[v] && !vis_vertex[v])
					{

						float d = trimesh::distance(_copy_mesh->vertices[v], _copy_mesh->vertices[frist_vertex]);
						face_line_len.push_back(std::make_pair(_copy_mesh->adjacentfaces[v][0], std::make_pair(len, len + d)));
						len += d;
						face_corss_point.push_back(std::make_pair(_copy_mesh->vertices[frist_vertex], _copy_mesh->vertices[v]));

						lines->vertices.push_back(_copy_mesh->vertices[frist_vertex] + trans_bbx_center);
						lines->vertices.push_back(_copy_mesh->vertices[v] + trans_bbx_center);

						frist_vertex = v;
						break;
					}
				}
				if (before_vertex == frist_vertex)
					break;
			}

			lines->write("lines.ply");

			trimesh::xform xxf = trimesh::inv(r_xxf);
			trimesh::apply_xform(_copy_mesh, xxf);
			//_copy_mesh->write("ratatemesh.ply");

			trimesh::TriMesh* resultp = new trimesh::TriMesh();

			float zaxis_cos = trimesh::vec3(0, 0, 1).dot(fn);
			trimesh::vec3 axis_to;
			trimesh::vec3 axis_up;
			if (std::abs(zaxis_cos - 1.f) < 1e-2)
			{
				axis_to = trimesh::vec3(0, 1, 0);
				axis_up = trimesh::vec3(0, 0, 1);
			}
			else
			{
				axis_to = trimesh::vec3(0, 0, 1);
				axis_up = trimesh::vec3(0, 1, 0);
			}
			for (int wi = 0; wi < init_font_meshs.size(); wi++)
			{
				float loc = word_init_location[wi].x - bbx.center().x;
				int sel_f = -1;
				trimesh::vec3 word_new_location(0, 0, 0);
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
				if (sel_f == -1)
					continue;
				trimesh::vec3 sel_fn = trimesh::normalized(_copy_mesh->trinorm(sel_f));
				word_FaceTo[wi] = sel_fn;

				float axis_cos = axis_to.dot(sel_fn);
				trimesh::vec3 up_dirto;
				if (std::abs(axis_cos - 1.0f) < 1e-2)
				{
					up_dirto = axis_up;
				}
				else
				{
					up_dirto = axis_to + axis_cos * -sel_fn;
				}
				word_Up[wi] = up_dirto;
				word_new_location = xxf * word_new_location;
				word_absolute_location[wi] = word_new_location + trans_bbx_center;
				resultp->vertices.push_back(word_new_location + trans_bbx_center);
			}
			resultp->write("resultp.ply");
			return false;
#elif 1
			
			trimesh::TriMesh* _copy_mesh = new trimesh::TriMesh;
			*_copy_mesh = *traget_meshes;
			_copy_mesh->clear_bbox();
			_copy_mesh->need_bbox();
			_copy_mesh->clear_across_edge();
			_copy_mesh->need_across_edge();	
			
			trimesh::vec3 trans_bbx_center = _copy_mesh->bbox.center();
			trimesh::trans(_copy_mesh, -trans_bbx_center);	
			//_copy_mesh->write("copymesh.ply");
			bool is_right = true;
			trimesh::xform xf = trimesh::xform::rot_into(fn, trimesh::vec3(0, -1, 0));
			if (trimesh::angle(trimesh::vec3(0, 1, 0), fn) < (M_PI / 2.0f))
			{
				xf = trimesh::xform::rot_into(fn, trimesh::vec3(0, 1, 0));
				is_right = false;
			}

			trimesh::vec3 new_face_to = xf * current_faceto;
			trimesh::xform rot_xf = trimesh::xform::rot(-m_config.angle * M_PI * 1.0f, new_face_to);
			trimesh::xform r_xxf = rot_xf * xf;
			trimesh::apply_xform(_copy_mesh, r_xxf);
			trimesh::vec3 _copy_location = r_xxf * (click_location - trans_bbx_center);
			float height = _copy_location.z;	
			
			trimesh::xform inve_xxf = trimesh::inv(r_xxf);
			//_copy_mesh->write("copymesh.ply");

			//trimesh::TriMesh* lines = new trimesh::TriMesh();					
			//---find height area and cal length
			std::unordered_map<int, std::pair<trimesh::vec3, trimesh::vec3>> faces_CrossPoints;		
			typedef typename std::unordered_map<int, std::pair<trimesh::vec3, trimesh::vec3>>::value_type faces_value;
			std::vector<bool> is_height_faces(_copy_mesh->faces.size(),false);
			int mistake_face = -1;
			float mistake_eps = std::numeric_limits<float>::max();
			for (int fi = 0; fi < _copy_mesh->faces.size(); fi++)
			{
				int v0 = _copy_mesh->faces[fi].at(0);
				int v1 = _copy_mesh->faces[fi].at(1);
				int v2 = _copy_mesh->faces[fi].at(2);

				if ((_copy_mesh->vertices[v0].z < height && _copy_mesh->vertices[v1].z < height && _copy_mesh->vertices[v2].z < height) ||
					(_copy_mesh->vertices[v0].z > height && _copy_mesh->vertices[v1].z > height && _copy_mesh->vertices[v2].z > height))
					continue;
				if (is_mistake)
				{
					trimesh::vec3 c= (_copy_mesh->vertices[v0] + _copy_mesh->vertices[v1] + _copy_mesh->vertices[v2]) / 3.0f;
					float d=trimesh::distance(c, _copy_location);
					if (d < mistake_eps)
					{
						mistake_eps = d;
						mistake_face = fi;
					}
				}
				is_height_faces[fi] = true;
				std::vector<trimesh::vec3> both_vertex;
				for (int vi = 0; vi < 3; vi++)
				{
					int v = _copy_mesh->faces[fi].at(vi);
					int v_n = _copy_mesh->faces[fi].at((vi + 1) % 3);
					trimesh::vec3 dir = _copy_mesh->vertices[v] - _copy_mesh->vertices[v_n];
					float h = _copy_mesh->vertices[v].z - _copy_mesh->vertices[v_n].z;
					if (std::abs(h) < 1e-4)
						continue;
					float scale = (height - _copy_mesh->vertices[v_n].z) / h * 1.0f;
					if (scale > 1 || scale < 0)
						continue;
					trimesh::vec3 new_position = _copy_mesh->vertices[v_n] + dir * scale;
					both_vertex.push_back(new_position);
				}
				if (both_vertex.size() == 2)
				{					
					faces_CrossPoints.insert(faces_value(fi, std::make_pair(both_vertex[0], both_vertex[1])));
					//lines->vertices.push_back(inve_xxf*both_vertex[0]);
					//lines->vertices.push_back(inve_xxf*both_vertex[1]);
				}
				
			}
			if (is_mistake)
				sel_faceid = mistake_face;
			if (sel_faceid == -1)
				return true;
			//lines->write("lines.ply");
			
			//---find orient and frist length
			float right_x = -std::numeric_limits<float>::max();
			float left_x = std::numeric_limits<float>::max();
			int right_face=-1;
			for (int ff = 0; ff < _copy_mesh->across_edge[sel_faceid].size(); ff++)
			{
				int ffi = _copy_mesh->across_edge[sel_faceid][ff];
				if (ffi == -1)
					continue;
				if (is_height_faces[ffi])
				{
					int v0 = _copy_mesh->faces[ffi].at(0);
					int v1 = _copy_mesh->faces[ffi].at(1);
					int v2 = _copy_mesh->faces[ffi].at(2);
					trimesh::vec3 c = (_copy_mesh->vertices[v0] + _copy_mesh->vertices[v1] + _copy_mesh->vertices[v2]) / 3.0f;
					if (is_right)
					{
						if (c.x > right_x)
						{
							right_x = c.x;
							right_face = ffi;
						}
					}
					else
					{
						if (c.x < left_x)
						{
							left_x = c.x;
							right_face = ffi;
						}
					}
				}
			}
			if (right_face < 0)
				return true;
		

			std::vector<std::pair<float,int>> faces_length;			
			auto pair_vertex = faces_CrossPoints[sel_faceid];
			trimesh::vec3 iter_vertex;
			trimesh::vec3 end_vertex;
			float len = 0;
			if (pair_vertex.first.x > _copy_location.x)
			{
				iter_vertex = pair_vertex.first;	
				end_vertex = pair_vertex.second;

				trimesh::vec3 temp_vertex = faces_CrossPoints[sel_faceid].first;
				faces_CrossPoints[sel_faceid].first = faces_CrossPoints[sel_faceid].second;
				faces_CrossPoints[sel_faceid].second = temp_vertex;
			}
			else
			{
				iter_vertex = pair_vertex.second;
				end_vertex = pair_vertex.first;
				
			}
			float dist = trimesh::distance(_copy_location, iter_vertex);
			len += dist;
			faces_length.push_back(std::make_pair(len, sel_faceid));

			//trimesh::TriMesh* lines1 = new trimesh::TriMesh();
			//---set step 
			std::vector<int> path_face_index;
			path_face_index.push_back(sel_faceid);
			path_face_index.push_back(right_face);
			std::vector<bool> vis_face(_copy_mesh->faces.size(),false);
			std::queue<int> que;
			que.push(right_face);
			vis_face[sel_faceid] = true;
			vis_face[right_face] = true;
						
			while (!que.empty())
			{			
				int f = que.front();
				auto pair_vertex = faces_CrossPoints[f];
				float d1 = trimesh::distance(pair_vertex.first, iter_vertex);
				float d2 = trimesh::distance(pair_vertex.second, iter_vertex);

				//lines1->vertices.push_back(pair_vertex.first);
				//lines1->vertices.push_back(pair_vertex.second);

				float dist = trimesh::distance(pair_vertex.first, pair_vertex.second);
				len += dist;
				faces_length.push_back(std::make_pair(len, f));
				if (d1 < d2)
				{
					iter_vertex = pair_vertex.second;
				}
				else
				{
					iter_vertex = pair_vertex.first;

					trimesh::vec3 temp_vertex = faces_CrossPoints[f].first;
					faces_CrossPoints[f].first = faces_CrossPoints[f].second;
					faces_CrossPoints[f].second = temp_vertex;
				}
				que.pop();				
				for (int fi = 0; fi < _copy_mesh->across_edge[f].size(); fi++)
				{
					int ff = _copy_mesh->across_edge[f][fi];
					if (ff == -1||vis_face[ff]|| !is_height_faces[ff])
						continue;
					que.push(ff);
					vis_face[ff] = true;			
					path_face_index.push_back(ff);
				}						
			}
			float last_dist = trimesh::distance(_copy_location, end_vertex);
			len += last_dist;
			faces_length.push_back(std::make_pair(len, sel_faceid));		
			//lines1->write("lines1.ply");
						

			//trimesh::TriMesh* lines1 = new trimesh::TriMesh();
			//set words location			
			trimesh::xform inv_xxf = trimesh::inv(r_xxf);
			for (int wi = 0; wi < init_font_meshs.size(); wi++)
			{
				float loc = word_init_location[wi].x - bbx.center().x;				
				trimesh::vec3 word_new_location(0, 0, 0);
				while(loc < 0)
				{
					loc = loc + len;
				}
				while(loc > len)
				{
					loc = loc - len;
				}
				for (int fl_i = 0; fl_i < faces_length.size(); fl_i++)
				{
					float before_length;
					if (fl_i == 0)
					{
						before_length = 0.f;
					}
					else
					{
						before_length = faces_length[fl_i - 1].first;
					}
					if (loc<faces_length[fl_i].first && loc>before_length)
					{
						float weight = (loc - before_length)/(faces_length[fl_i].first- before_length);
						int sf = faces_length[fl_i].second;
						auto pair_vertex = faces_CrossPoints[sf];
						trimesh::vec3 dirto = pair_vertex.second - pair_vertex.first;						
						trimesh::vec3 new_location = (pair_vertex.first + weight* dirto);												
						//lines1->vertices.push_back(new_location);
						trimesh::vec3 world_location = (inv_xxf * new_location) + trans_bbx_center;
					
						word_absolute_location[wi] = world_location;
						trimesh::vec3 sel_fn = trimesh::normalized(traget_meshes->trinorm(sf));
						word_FaceTo[wi] = sel_fn;

						trimesh::vec3 rot_sel_fn = trimesh::normalized(_copy_mesh->trinorm(sf));
						trimesh::vec3 up_dirto = trimesh::normalized(rot_sel_fn.cross(dirto));					
						trimesh::vec3 ori_up_dirto = inv_xxf * up_dirto;
						if (ori_up_dirto.z < 0)
						{
							trimesh::xform rota = trimesh::xform::rot( M_PI * 1.0f, sel_fn);
							ori_up_dirto = rota * ori_up_dirto;
						}
						word_Up[wi] = ori_up_dirto;

						/*float axis_cos = trimesh::vec3(0,0,1).dot(sel_fn);
						trimesh::vec3 up_dirto;
						if (std::abs(axis_cos - 1.0f) < 1e-2)
						{
							up_dirto = trimesh::vec3(0, 1, 0);
						}
						else
						{
							up_dirto = trimesh::vec3(0, 0, 1) + axis_cos * -sel_fn;
						}
						word_Up[wi] = up_dirto;	*/				
						break;
					}
				}
			}
			
			//lines1->write("lines1.ply");
			return false;
#endif
		}		
	}


	void FontMesh::updateXform(trimesh::xform xform)
	{
		_m_xform = xform;
		m_config.height = _m_xform[10] * m_config.height;
		click_location = _m_xform * click_location;
		FaceTo.second = trimesh::normalized(_m_xform * FaceTo.second);
		Up.second = trimesh::normalized(_m_xform * Up.second);
		for (int i = 0; i < word_FaceTo.size(); i++)
		{
			word_FaceTo[i] = trimesh::normalized(_m_xform* word_FaceTo[i]);
			word_Up[i] = trimesh::normalized(_m_xform * word_Up[i]);
		}
	}

	trimesh::xform FontMesh::xform()
	{
		return _m_xform;
	}

	void FontMesh::rotateFontMesh(trimesh::TriMesh* traget_mesh, float angle)
	{				
		m_config.angle = angle;
		if (!m_config.state)
		{		
			
		}
		else
		{
			if(traget_mesh)
				FontTransform(traget_mesh, sel_faceid, click_location, m_config.state);
		}
	}



	void FontMesh::updateFontPoly(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter)
	{		
		for (auto& mesh : init_font_meshs)
		{
			mesh->clear(); mesh = nullptr;
		}
		init_font_meshs.clear();		
		CreateFontMesh(letter, trimesh::vec3(0, 0, -1), trimesh::vec3(0, -1, 0), true,false);
	}

	void FontMesh::updateFontHeight(float height)
	{
		float cha = (height - m_config.height)/2.0f;
		for (int mi = 0; mi < init_font_meshs.size(); mi++)
		{
			int half_v = init_font_meshs[mi]->vertices.size() / 2;
			for (int vi = 0; vi < init_font_meshs[mi]->vertices.size(); vi++)
			{
				if (vi < half_v)
				{
					init_font_meshs[mi]->vertices[vi].z -= cha;
				}
				else
				{
					init_font_meshs[mi]->vertices[vi].z += cha;
				}			
			}	
		}		
		is_init_location = false;
		if(!m_config.state)
			InitFontMesh();
		m_config.height = height;
	}


	void FontMesh::updateFontDepth(float depth)
	{
		m_config.distance = depth;
	}


	void FontMesh::CreateFontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter,
		trimesh::vec3 face_to , trimesh::vec3 up, bool is_adjust, bool is_init)
	{
		word_FaceTo.clear();
		word_Up.clear();
		word_init_location.clear();
		word_absolute_location.clear();		
		init_font_meshs.clear();
		bbx.clear();		
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
					if (totalpoly[pi][ppi].y > wordbbx_max.y)
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
				v.p.z += m_config.height;
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
			//font_meshs.push_back(_word_mesh);

			word_init_location.push_back(_word_mesh->bbox.center());
			bbx += _word_mesh->bbox;
			trimesh::trans(_word_mesh,-_word_mesh->bbox.center());
			
			if (is_adjust)
			{
				trimesh::xform xf = trimesh::xform::rot_into(face_to,trimesh::vec3(0,0,1));
				trimesh::vec3 new_up = xf * up;
				trimesh::xform up_xf= trimesh::xform::rot_into(new_up, trimesh::vec3(0, 1, 0));
				trimesh::apply_xform(_word_mesh, up_xf * xf);
				word_FaceTo.push_back(trimesh::vec3(0,0,1));
				word_Up.push_back(trimesh::vec3(0,1,0));
			}
			else
			{
				word_FaceTo.push_back(face_to);
				word_Up.push_back(up);
			}

			word_absolute_location.push_back(trimesh::vec3(0, 0, 0));
			init_font_meshs.push_back(_word_mesh);				
		}		
		is_init_location = is_init;
		is_init_adjust = is_adjust;
		InitFontMesh();		
	}

	void FontMesh::InitFontMesh()
	{
		//is_change = false;
		_return_mesh->clear();
		for (int wi = 0; wi < init_font_meshs.size(); wi++)
		{
			int vsize = _return_mesh->vertices.size();
			trimesh::vec3 transTo = word_init_location[wi];
			for (trimesh::point v : init_font_meshs[wi]->vertices)
				_return_mesh->vertices.push_back(v + transTo);
			for (trimesh::TriMesh::Face f : init_font_meshs[wi]->faces)
			{
				_return_mesh->faces.push_back(trimesh::TriMesh::Face(vsize + f[0], vsize + f[1], vsize + f[2]));
			}
			word_absolute_location[wi] = word_init_location[wi];				
		}
		
		if (is_init_location)
		{
			sel_faceid = -1;
			click_location = trimesh::vec3(0, 0, 0);
			if (is_init_adjust)
			{
				FaceTo.first = trimesh::vec3(0, 0, 1);
				FaceTo.second = trimesh::vec3(0, 0, 1);
				Up.first = trimesh::vec3(0, 1, 0);
				Up.second = trimesh::vec3(0, 1, 0);
			}
			else
			{
				FaceTo.first = trimesh::vec3(0, 0, -1);
				FaceTo.second = trimesh::vec3(0, 0, -1);
				Up.first = trimesh::vec3(0, -1, 0);
				Up.second = trimesh::vec3(0,-1, 0);
			}
			
		}
		_return_mesh->need_bbox();
		trimesh::trans(_return_mesh, -_return_mesh->bbox.center());
	}

}