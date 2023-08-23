#pragma once
#include<vector>
#include "trimesh2/TriMesh.h"


namespace topomesh {
	class EarClipping {
	public:
		EarClipping() {};
		EarClipping(std::vector<std::pair<trimesh::point, int>>& data);
		~EarClipping() {};

		std::vector<trimesh::ivec3>& const getResult() { return this->_result; }
	private:
		std::vector<trimesh::ivec3> _result;
	};
}