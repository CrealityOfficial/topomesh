#include "honeycomb.h"
#include "topomesh/internal/util.h"

namespace topomesh
{
	trimesh::TriMesh* innerGenerateHoneyCombs(trimesh::TriMesh* trimesh, const HoneyCombParam& honeyparams,
		ccglobal::Tracer* tracer, HoneyCombDebugger* debugger)
	{
		if (!trimesh)
			return nullptr;

		MMeshT mesh(trimesh);

		std::vector<int> upperIndex;
		std::vector<int> lowerIndex;
		TriPolygons meshPoly;

		if (debugger)
			debugger->onGenerateBottomPolygons(meshPoly);

		honeyCombLetter(mesh, upperIndex, lowerIndex, meshPoly);
		return toTrimesh(mesh);
	}

	void honeyCombLetter(MMeshT& mesh, const std::vector<int>& upper, const std::vector<int>& lower,
		const TriPolygons& meshPoly, HoneyCombDebugger* debugger)
	{
		TriPolygons combPoly;
		TriPolygons result;
		merge(meshPoly, combPoly, result);

		if (debugger)
			debugger->onGenerateInfillPolygons(meshPoly);
		//
	}
}