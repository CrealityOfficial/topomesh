#pragma once
#include "mmesht.h"

namespace topomesh {

	class BaseParameterClass {};

	class InternelData
	{
	public:
		InternelData() {};
		InternelData(trimesh::TriMesh* trimesh) :_mesh(trimesh){};
		~InternelData() { _mesh.clear(); };
	private:
		MMeshT _mesh;

	public:
		void chunkedMesh(int n);
		trimesh::TriMesh* mmesht2trimesh(bool save = false);
		void loopSubdivsion(std::vector<int>& faceindexs, std::vector<std::tuple<int, trimesh::point>>& vertex,
			std::vector<std::tuple<int, trimesh::ivec3>>& face_vertex, bool is_move=false,int iteration=1);
		void SimpleSubdivsion(std::vector<int>& faceindexs);
	};
}