#ifndef TOPOMESH_UTILS_1687934636754_H
#define TOPOMESH_UTILS_1687934636754_H
#include "topomesh/data/convert.h"
#include "topomesh/data/mmesht.h"

namespace topomesh
{
	struct ColumnarParam
	{
		float zStart = 0.0f;
		float zEnd = 10.0f;
	};

	trimesh::TriMesh* generateColumnar(const TriPolygons& polys, const ColumnarParam& param,const std::vector<std::vector<float>>* height=nullptr);
	void extendPoint(MMeshT* mesh,std::vector<int>& vertex_index, const ColumnarParam& param, const std::vector<float>* height = nullptr);
	void findNeighborFacesOfSameAsNormal(MMeshT* mesh, int indicate, std::vector<int>& faceIndexs,float angle_threshold,bool vis = false);
	void findNeighborFacesOfConsecutive(MMeshT* mesh, int indicate, std::vector<int>& faceIndexs, float angle_threshold, bool vis = false);
}

#endif // TOPOMESH_UTILS_1687934636754_H