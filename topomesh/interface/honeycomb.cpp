#include "honeycomb.h"

namespace topomesh
{
	HoneyComb::HoneyComb(TopoTriMeshPtr mesh)
	{
		//mesh->need_adjacentfaces();
		//mesh->need_neighbors();
		//mesh->need_normals();

		m_proxy.reset(new MMProxy(mesh));
	}

	HoneyComb::~HoneyComb()
	{

	}

	void HoneyComb::checkPlane(int indicate, std::vector<int>& faceIndexs)
	{
		return m_proxy->checkNeigbour(indicate, faceIndexs, 89.0f);
	}

	trimesh::TriMesh* HoneyComb::generateModelFitComb(const CombParam& param, const trimesh::vec3& normal)
	{
		return nullptr;
	}

	trimesh::TriMesh* HoneyComb::createInnerCombModel(const CombParam& param, const trimesh::vec3& normal)
	{
		return nullptr;
	}

	trimesh::TriMesh* HoneyComb::createExtCombModel(const CombParam& param, const trimesh::vec3& normal, const std::vector<int>& indices)
	{
		return nullptr;
	}
}