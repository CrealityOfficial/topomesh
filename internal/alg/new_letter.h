#pragma once
#include "trimesh2/TriMesh.h"
#include "topomesh/interface.h"

namespace topomesh {
	TOPOMESH_API trimesh::TriMesh* CreateFontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter,float height);
}