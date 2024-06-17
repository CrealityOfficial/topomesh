#pragma once
#include "trimesh2/TriMesh.h"


namespace topomesh {
	trimesh::TriMesh* CreateFontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter,float height);
}