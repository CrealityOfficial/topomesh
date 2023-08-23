#ifndef TOPOMESH_HONEYCOMB_1692753685514_H
#define TOPOMESH_HONEYCOMB_1692753685514_H
#include "topomesh/interface/idata.h"
#include "topomesh/interface/mmproxy.h"

namespace topomesh
{
	class HoneyComb
	{
	public:
		HoneyComb(TopoTriMeshPtr mesh);
		~HoneyComb();

		void checkNeigbour(int indicate, std::vector<int>& faceIndexs, float angle_threshold);
	protected:
		std::unique_ptr<MMProxy> m_proxy;
	};
}

#endif // TOPOMESH_HONEYCOMB_1692753685514_H