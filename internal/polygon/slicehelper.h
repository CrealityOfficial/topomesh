#ifndef SLICE_SLICEHELPER_1598712992630_H
#define SLICE_SLICEHELPER_1598712992630_H
#include "trimesh2/TriMesh.h"
#include "trimesh2/quaternion.h"
#include <unordered_map>

namespace topomesh
{
	class MeshFace
	{
	public:
		int vertex_index[3] = { -1 };
		int connected_face_index[3];
	};

	class SliceHelper
	{
	public:
		SliceHelper();
		~SliceHelper();

		void prepare(trimesh::TriMesh* _mesh);
		void getMeshFace();
		std::vector<std::vector<std::pair<uint32_t, int>>> generateVertexConnectVertexData();
		void generateConcave(std::vector<trimesh::vec3>& concave, const trimesh::quaternion* rotation, const trimesh::vec3 scale);
	protected:
		std::vector<MeshFace> faces;
		std::vector<std::vector<uint32_t>> vertexConnectFaceData;

		trimesh::TriMesh* meshSrc;
	};
}

#endif // SLICE_SLICEHELPER_1598712992630_H