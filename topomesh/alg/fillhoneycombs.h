#pragma once
#include "topomesh/data/mmesht.h"
#include "topomesh/alg/convert.h"
#include "topomesh/alg/letter.h"

#include "ccglobal/tracer.h"


typedef trimesh::point vec3;
typedef trimesh::vec2 vec2;

namespace topomesh {
	void GenerateHoneyCombs(const trimesh::TriMesh* mesh, trimesh::TriMesh& resultmesh, const TriPolygon& poly, vec3 axisDir = vec3(0, 0, 1),
		vec2 arrayDir = vec2(1, 0), double honeyCombRadius = 3, double nestWidth = 1, double shellThickness = 2);
	void findNeighVertex(trimesh::TriMesh* mesh, const std::vector<int>& upfaceid,const std::vector<int>& botfaceid);
	void innerHex(MMeshT* mesh, std::vector<std::vector<trimesh::vec2>>& poly, std::vector<int>& inFace, std::vector<int>& outFace,float len);
}