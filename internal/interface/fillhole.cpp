#include "topomesh/interface/fillhole.h"

#include "internal/data/CMesh.h"
namespace topomesh
{
	void findHoles(trimesh::TriMesh* tmesh, /*out*/ std::vector<std::vector<int>>& holes)
	{
		if (!tmesh)
			return;

		CMesh mesh(tmesh);
		std::vector<int> edges;
		mesh.SelectIndividualEdges(edges);
		mesh.GetSequentialPoints(edges, holes);
	}
}