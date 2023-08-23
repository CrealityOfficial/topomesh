#ifndef TOPOMESH_HONEYCOMB_1692753685514_H
#define TOPOMESH_HONEYCOMB_1692753685514_H
#include "topomesh/interface/idata.h"
#include "topomesh/interface/mmproxy.h"

namespace topomesh
{
	struct CombParam
	{
		float width = 1.0f;
		float radius = 1.0f;

		bool holeConnect = false;
		float holeHeight = 1.0f;
		float holeRadius = 1.0f;
		float holeGap = 1.0f;
	};

	class HoneyComb
	{
	public:
		HoneyComb(TopoTriMeshPtr mesh);
		~HoneyComb();

		void checkPlane(int indicate, std::vector<int>& faceIndexs);
	
		trimesh::TriMesh* generateModelFitComb(const CombParam& param, const trimesh::vec3& normal);

		trimesh::TriMesh* createInnerCombModel(const CombParam& param, const trimesh::vec3& normal);
		trimesh::TriMesh* createExtCombModel(const CombParam& param, const trimesh::vec3& normal, const std::vector<int>& indices);
	protected:
		std::unique_ptr<MMProxy> m_proxy;
	};
}

#endif // TOPOMESH_HONEYCOMB_1692753685514_H