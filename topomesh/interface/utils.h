#ifndef TOPOMESH_UTILS_1695108550958_H
#define TOPOMESH_UTILS_1695108550958_H
#include "topomesh/interface/idata.h"

namespace topomesh
{
	TOPOMESH_API void findNeignborFacesOfSameAsNormal(trimesh::TriMesh* trimesh, int indicate, float angle_threshold,
									/*out*/ std::vector<int>& faceIndexs);

	TOPOMESH_API void findBoundary(trimesh::TriMesh* trimesh);
	TOPOMESH_API void triangulate(trimesh::TriMesh* trimesh,std::vector<int>& sequentials);
}

#endif // TOPOMESH_UTILS_1695108550958_H