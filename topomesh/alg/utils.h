#ifndef TOPOMESH_UTILS_1687934636754_H
#define TOPOMESH_UTILS_1687934636754_H
#include "topomesh/data/convert.h"

namespace topomesh
{
	struct ColumnarParam
	{
		float zStart = 0.0f;
		float zEnd = 10.0f;
	};

	trimesh::TriMesh* generateColumnar(const TriPolygons& polys, const ColumnarParam& param);
}

#endif // TOPOMESH_UTILS_1687934636754_H