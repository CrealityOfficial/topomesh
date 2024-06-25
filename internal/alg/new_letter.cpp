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
			//_word_mesh->write("_word_mesh.ply");		
			_word_mesh->need_bbox();
			word_mesh_center.push_back(_word_mesh->bbox.center());
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
		}
		//_return_mesh->write("returnmesh.ply");
		return _return_mesh;
	}


	void MeshTransform(trimesh::TriMesh* traget_meshes,  trimesh::TriMesh* font_mesh, int face_id,
		trimesh::vec3 location, trimesh::vec3 dir, std::vector<float>& word_location, std::vector<int>& mesh_vertex_sizes, std::vector<trimesh::vec3>& word_mesh_center,
		trimesh::vec3 up,bool is_surround)
	{		
		if (!is_surround)
		{
			font_mesh->need_bbox();
			float hz = font_mesh->bbox.max.z - font_mesh->bbox.min.z;
			trimesh::vec3 dirTo(0, 1, 0);
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
			
			trimesh::apply_xform(font_mesh, trimesh::xform::rot_into(up, dirTo));
			//trimesh::apply_xform(font_mesh, trimesh::xform::rot_into(trimesh::vec3(0, 0, 1), dir));
			trimesh::vec3 trans = location + ((2.0f * hz) / 5.0f) * dir;
			trimesh::trans(font_mesh, trans);
		}
		else
		{
			trimesh::TriMesh* _copy_mesh = new trimesh::TriMesh();
			_copy_mesh = traget_meshes;			
			trimesh::xform xf = trimesh::xform::rot_into(trimesh::vec3(0,0,1), trimesh::vec3(0, 0, 1));
			trimesh::apply_xform(_copy_mesh, xf);
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
					bool pass = false;				
					
					float z0 = _copy_mesh->vertices[_copy_mesh->faces[ff][0]].z;
					float z1 = _copy_mesh->vertices[_copy_mesh->faces[ff][1]].z;
					float z2 = _copy_mesh->vertices[_copy_mesh->faces[ff][2]].z;
					float max_z = std::max({ z0,z1,z2 });
					float min_z = std::min({ z0,z1,z2 });
					if (height >= min_z && height <= max_z)
						que.push(ff);
				}					
			}

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
							float vv_sacle = (abs_wl - left_clip_faces[li].second.second) / (left_clip_faces[li].second.first - left_clip_faces[li].second.second);
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
							float vv_sacle = (abs_wl - right_clip_faces[ri].second.second) / (right_clip_faces[ri].second.first - right_clip_faces[ri].second.second);
							trimesh::vec3 vv_dir = right_face_corss_point[ri].second - right_face_corss_point[ri].first;
							new_wordpoint = right_face_corss_point[ri].first + vv_sacle * vv_dir;
							break;
						}
					}
				}

				if (f != -1)
				{
					trimesh::vec3 f_n = trimesh::normalized(traget_meshes->trinorm(f));
					trimesh::xform xf_f = trimesh::xform::rot_into(trimesh::vec3(0, 0, 1), f_n);
					int b = mi + 1;
					for (int vi = mesh_vertex_sizes[b-1]; vi < mesh_vertex_sizes[b]; vi++)
					{
						trimesh::vec3 new_point = xf_f * font_mesh->vertices[vi];
						font_mesh->vertices[vi] = new_point;
					}
					trimesh::vec3 new_word_bbx_center = xf_f * word_mesh_center[mi];
					trimesh::vec3 word_trans = new_wordpoint - new_word_bbx_center;
					for (int vi = mesh_vertex_sizes[b-1]; vi < mesh_vertex_sizes[b]; vi++)
					{
						font_mesh->vertices[vi] += word_trans;
					}
				}
			}

		}
		//trimesh::trans(font_mesh, location);
	}
}