#pragma once
#include "trimesh2/TriMesh.h"

namespace topomesh {
	float getMeshVolume(trimesh::TriMesh* mesh,std::vector<int>& faces);
}