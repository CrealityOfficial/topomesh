#include "mmproxy.h"
#include "topomesh/data/mmesht.h"

#include "topomesh/alg/utils.h"

namespace topomesh
{
	MMProxy::MMProxy(TopoTriMeshPtr m)
		: mesh(m)
	{
		innerMesh.reset(new MMeshT(mesh.get()));
	}

	MMProxy::~MMProxy()
	{

	}

	void MMProxy::checkNeigbour(int indicate, std::vector<int>& faceIndexs, float angle_threshold)
	{
		findNeighborFacesOfSameAsNormal(innerMesh.get(), indicate, faceIndexs, angle_threshold);
	}
}