#include "utils.h"
#include "topomesh/data/mmesht.h"
#include "topomesh/alg/letter.h"

namespace topomesh
{
	trimesh::TriMesh* generateColumnar(const TriPolygons& polys, const ColumnarParam& param)
	{
		/*TriPolygon poly = { trimesh::vec3(0.5f,0.0f,0.0f),trimesh::vec3(1.0f,0.0f,0),trimesh::vec3(1.5f,0.5f,0.0f),trimesh::vec3(1.0f,1.0f,0.0),
			trimesh::vec3(0.5f,1.0f,0.0),trimesh::vec3(0.0f,0.5f,0.0) };
		TriPolygons polys1 = {poly};*/
		MMeshT mt;
		for (int i = 0; i < polys.size(); i++)
		{			
			std::vector<int> vertex_index;			
			for (int j = 0; j < polys[i].size(); j++)
			{
				mt.appendVertex(polys[i][j]);
				vertex_index.push_back(mt.vertices.back().index);
			}			
			std::vector<int> copy_vertex = vertex_index;
			int before_face_index = mt.faces.size();			
			if (vertex_index.size() < 3) continue;		
			fillTriangleForTraverse(&mt, vertex_index);
			int after_face_index = mt.faces.size();
			for (int fi = before_face_index; fi < after_face_index; fi++)
			{
				if(((mt.faces[fi].V0(1)->p - mt.faces[fi].V0(0)->p)%(mt.faces[fi].V0(2)->p - mt.faces[fi].V0(0)->p)).z>0)
					std::swap(mt.faces[fi].connect_vertex[1], mt.faces[fi].connect_vertex[2]);
			}
			int before_vertex_index = mt.vertices.size();
			extendPoint(&mt, copy_vertex, param);
			int after_vertex_index = mt.vertices.size();
			vertex_index.clear();
			for (int vi = before_vertex_index; vi < after_vertex_index; vi++)
				vertex_index.push_back(vi);
			fillTriangleForTraverse(&mt, vertex_index);
		}
		trimesh::TriMesh* mesh = new trimesh::TriMesh();
		mt.quickTransform(mesh);
		//mesh->write("columnar.ply");
		return mesh;
	}

	void extendPoint(MMeshT* mesh, std::vector<int>& vertex_index, const ColumnarParam& param)
	{		
		std::vector<int> copy_vertex_index = vertex_index;
		int size = copy_vertex_index.size();
		for (int vi = 0; vi < vertex_index.size(); vi++)
		{
			mesh->appendVertex(mesh->vertices[vertex_index[vi]].p + trimesh::point(0, 0, param.zEnd));
			mesh->appendFace(vertex_index[vi], copy_vertex_index[(vi + 1) % size], mesh->vertices.back().index);
			mesh->appendFace(vertex_index[vi], copy_vertex_index[(vi + size - 1) % size], mesh->vertices.back().index);
			copy_vertex_index[vi] = mesh->vertices.back().index;
		}
	}


	void findNeighborFacesOfSameAsNormal(MMeshT* mesh, int indicate, std::vector<int>& faceIndexs, float angle_threshold)
	{
		if (!mesh->is_FaceNormals()) mesh->getFacesNormals();
		for (MMeshFace& f : mesh->faces)
			f.ClearS();
		trimesh::point normal=mesh->faces[indicate].normal;	
		std::queue<int> queue;
		queue.push(indicate);
		while (!queue.empty())
		{
			faceIndexs.push_back(queue.front());
			mesh->faces[queue.front()].SetS();
			for (int i = 0; i < mesh->faces[queue.front()].connect_face.size(); i++)
			{
				if (!mesh->faces[queue.front()].connect_face[i]->IsS())
				{
					mesh->faces[queue.front()].connect_face[i]->SetS();
					float arc = trimesh::normalized(mesh->faces[indicate].normal) ^ trimesh::normalized(mesh->faces[queue.front()].connect_face[i]->normal);
					if ((std::acos(arc) * 180 / M_PI) < angle_threshold)
						queue.push(mesh->faces[queue.front()].connect_face[i]->index);
				}
			}
			queue.pop();
		}	
	}
}