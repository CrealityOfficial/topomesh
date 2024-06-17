#include "new_letter.h"
#include "letter.h"
#include "earclipping.h"
#include "trimesh2/TriMesh_algo.h"

namespace topomesh {
	trimesh::TriMesh* CreateFontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter, float height)
	{
		trimesh::TriMesh* _return_mesh = new trimesh::TriMesh();		
		for (int li = 0; li < letter.size(); li++)
		{
			MMeshT mt(5000, 10000);
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

			mt.shrinkMesh();
			mt.init_halfedge();
			mt.set_HalfEdge(false);
			for (MMeshVertex& v : mt.vertices)
				v.p.z += height;
			int face_size = mt.faces.size();
			for (int fi = 0; fi < face_size; fi++)
			{
				MMeshFace& f = mt.faces[fi];
				if (f.IsD()) continue;
				mt.appendVertex(trimesh::point(f.V0(0)->p.x, f.V0(0)->p.y, 0));
				mt.appendVertex(trimesh::point(f.V0(1)->p.x, f.V0(1)->p.y, 0));
				mt.appendVertex(trimesh::point(f.V0(2)->p.x, f.V0(2)->p.y, 0));
				mt.appendFace(mt.VN() - 1, mt.VN() - 2, mt.VN() - 3);
				MMeshHalfEdge* p = f.f_mhe;
				int n = 3;
				while (1)
				{
					if (p->opposite == nullptr)
					{
						int i1 = p->edge_vertex.second->index;
						int i2 = p->edge_vertex.first->index;
						int i3 = mt.VN() - n;
						int i4 = mt.VN() - (n - 1);
						if (n == 1)
							i4 = mt.VN() - 3;
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
			int vsize = _return_mesh->vertices.size();
			for (trimesh::point v : _word_mesh->vertices)
				_return_mesh->vertices.push_back(v);
			for (trimesh::TriMesh::Face f : _word_mesh->faces)
			{
				_return_mesh->faces.push_back(trimesh::TriMesh::Face(vsize+f[0], vsize + f[1], vsize + f[2]));
			}
		}
		
		
		_return_mesh->write("returnmesh.ply");
		return _return_mesh;
	}
}