#ifndef TOPOMESH_MMPROXY_1692753685513_H
#define TOPOMESH_MMPROXY_1692753685513_H
#include "topomesh/interface/idata.h"

namespace topomesh
{
	class MMeshT;
	class MMProxy
	{
	public:
		MMProxy(TopoTriMeshPtr m);
		~MMProxy();

		void checkNeigbour(int indicate, std::vector<int>& faceIndexs, float angle_threshold);
	protected:
		TopoTriMeshPtr mesh;
		std::shared_ptr<MMeshT> innerMesh;
	};
}

#endif // TOPOMESH_MMPROXY_1692753685513_H