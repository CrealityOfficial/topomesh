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

	void HoneyComb::checkNeigbour(int indicate, std::vector<int>& faceIndexs, float angle_threshold)
	{
		return m_proxy->checkNeigbour(indicate, faceIndexs, angle_threshold);
	}
}