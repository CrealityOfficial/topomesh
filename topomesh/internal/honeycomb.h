#ifndef TOPOMESH_HONEYCOMB_1687839380318_H
#define TOPOMESH_HONEYCOMB_1687839380318_H
#include "topomesh/alg/fillhoneycombs.h"
#include "topomesh/data/mmesht.h"

namespace topomesh
{
	trimesh::TriMesh* innerGenerateHoneyCombs(trimesh::TriMesh* trimesh, const HoneyCombParam& honeyparams = HoneyCombParam(),
		ccglobal::Tracer* tracer = nullptr, HoneyCombDebugger* debugger = nullptr);

	void honeyCombLetter(MMeshT& mesh, const std::vector<int>& upper, const std::vector<int>& lower, 
		const TriPolygons& polygons, HoneyCombDebugger* debugger = nullptr);
}

#endif // TOPOMESH_HONEYCOMB_1687839380318_H