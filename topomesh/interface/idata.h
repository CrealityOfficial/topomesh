#ifndef TOPOMESH_IDATA_1692613164081_H
#define TOPOMESH_IDATA_1692613164081_H
#include "topomesh/interface.h"
#include "trimesh2/TriMesh.h"
#include <memory>

namespace topomesh
{
	typedef std::vector<trimesh::vec3> TriPolygon;
	typedef std::vector<TriPolygon> TriPolygons;

	struct HexaPolygon {
		TriPolygon poly;
		trimesh::ivec3 coord;
		std::vector<int> neighbors;
		HexaPolygon() : neighbors(6, -1) {}
	};
	typedef std::vector<HexaPolygon> HexaPolygons;

	struct SimpleCamera
	{
		float f;
		float n;
		float fov;
		float aspect;

		trimesh::point pos;
		trimesh::point center;
		trimesh::point up;
	};

	typedef std::vector<int> FacePatch;
	typedef std::vector<FacePatch> FacePatchs;

	typedef std::shared_ptr<trimesh::TriMesh> TopoTriMeshPtr;
}

#endif // TOPOMESH_IDATA_1692613164081_H