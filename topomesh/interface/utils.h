#ifndef TOPOMESH_UTILS_1695108550958_H
#define TOPOMESH_UTILS_1695108550958_H
#include "topomesh/interface/idata.h"

namespace topomesh
{
	TOPOMESH_API void findNeignborFacesOfSameAsNormal(trimesh::TriMesh* trimesh, int indicate, float angle_threshold,
									/*out*/ std::vector<int>& faceIndexs);
}

#endif // TOPOMESH_UTILS_1695108550958_H