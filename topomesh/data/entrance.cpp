#include "entrance.h"
#include "topomesh/alg/subdivision.h"
#include "topomesh/alg/segmentation.h"
#include "topomesh/alg/remesh.h"

namespace topomesh {
	trimesh::TriMesh* InternelData::mmesht2trimesh(bool save)
	{
		trimesh::TriMesh* newmesh = new trimesh::TriMesh();
		if (!save)
			this->_mesh.quickTransform(newmesh);
		else
			this->_mesh.mmesh2trimesh(newmesh);
		return newmesh;
	}

	void InternelData::loopSubdivsion(std::vector<int>& faceindexs, std::vector<std::tuple<int, trimesh::point>>& vertex,
		std::vector<std::tuple<int, trimesh::ivec3>>& face_vertex, bool is_move, int iteration)
	{
		topomesh::loopSubdivision(&this->_mesh, faceindexs, vertex, face_vertex,is_move, iteration);
		return;
	}

	void InternelData::SimpleSubdivsion(std::vector<int>& faceindexs)
	{
		topomesh::SimpleMidSubdivision( &this->_mesh,faceindexs);
		return;
	}

	void InternelData::chunkedMesh(int n)
	{
		if (!this->_mesh.is_BoundingBox()) this->_mesh.getBoundingBox();
		float x = (this->_mesh.boundingbox.max_x - this->_mesh.boundingbox.min_x) / n * 1.f;
		float y = (this->_mesh.boundingbox.max_y - this->_mesh.boundingbox.min_y) / n * 1.f;
		float z = (this->_mesh.boundingbox.max_z - this->_mesh.boundingbox.min_z) / n * 1.f;

		for (MMeshFace& f : this->_mesh.faces)
		{
			if (f.GetU() > 0) continue;
			trimesh::point c = (f.V0(0)->p + f.V1(0)->p + f.V2(0)->p) / 3.f;
			int xi = c.x / x;
			int yi = c.y / y;
			int zi = c.z / z;
			int chunked = zi * std::pow(n, 2) + yi * n + xi + 1;
			f.SetU(chunked);
		}
	}

	void InternelData::getChunkFace(int ni ,std::vector<int>& faceindexs)
	{
		int deleteFnum = 0;
		for (MMeshFace& f : this->_mesh.faces)
		{
			if (f.IsD())
			{
				deleteFnum++; continue;
			}
			if (f.GetU() == ni)
				faceindexs.push_back(f.index-deleteFnum);
		}
	}

	void InternelData::getVertex(const std::vector<int>& faceindexs, std::vector<std::tuple<trimesh::ivec3, trimesh::point, trimesh::point, trimesh::point>>& vertexindexs)
	{
		for (int i : faceindexs)
		{
			MMeshFace& f = this->_mesh.faces[i];
			vertexindexs.push_back(std::make_tuple(trimesh::ivec3(f.V0(0)->index, f.V1(0)->index, f.V2(0)->index),
				f.V0(0)->p, f.V1(0)->p, f.V2(0)->p));
		}
	}

	void InternelData::SimpleRemeshing(const std::vector<int>& faceindexs, float thershold)
	{
		return topomesh::SimpleRemeshing(&this->_mesh, faceindexs, thershold);
	}
}