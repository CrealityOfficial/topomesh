#pragma once
#include <vector>
#include "trimesh2/TriMesh.h"

typedef trimesh::point point;

namespace topomesh
{
	bool rayTriangleIntersect(const point&,const point&,const point&,const point&,const point&,point&);
}