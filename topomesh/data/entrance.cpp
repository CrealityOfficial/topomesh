#include "entrance.h"
#include "topomesh/alg/subdivision.h"
#include "topomesh/alg/segmentation.h"

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
		this->_mesh.boundingbox;
	}
}