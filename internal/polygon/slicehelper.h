#ifndef SLICE_SLICEHELPER_1598712992630_H
#define SLICE_SLICEHELPER_1598712992630_H
#include "trimesh2/TriMesh.h"
#include "trimesh2/quaternion.h"
#include "point2.h"
#include "vertex.h"
#include <unordered_map>

namespace topomesh
{
	class SliceHelper
	{
	public:
		SliceHelper();
		~SliceHelper();

		void prepare(trimesh::TriMesh* _mesh);
		void getMeshFace();
		std::vector<std::vector<std::pair<uint32_t, int>>> generateVertexConnectVertexData();
		void generateConcave(std::vector<trimesh::vec3>& concave, const trimesh::quaternion* rotation, const trimesh::vec3 scale);
		void buildMeshFaceHeightsRange(const trimesh::TriMesh* meshSrc, std::vector<Point2>& heightRanges);
	public:
		std::vector<MeshFace> faces;
		std::vector<std::vector<uint32_t>> vertexConnectFaceData;

		trimesh::TriMesh* meshSrc;
		std::vector<Point2> faceRanges;
	};
}

#endif // SLICE_SLICEHELPER_1598712992630_H