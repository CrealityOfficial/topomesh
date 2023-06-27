#ifndef TOPOMESH_UTIL_1687839380331_H
#define TOPOMESH_UTIL_1687839380331_H
#include "topomesh/data/mmesht.h"
#include "topomesh/data/convert.h"

namespace topomesh
{
	trimesh::TriMesh* toTrimesh(const MMeshT& mesh);

	void merge(const TriPolygons& poly1, const TriPolygons& poly2, TriPolygons& out);
}

#endif // TOPOMESH_UTIL_1687839380331_H