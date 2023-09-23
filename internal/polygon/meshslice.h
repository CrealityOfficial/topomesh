#ifndef CX_MESHSLICE_1599726705111_H
#define CX_MESHSLICE_1599726705111_H
#include <vector>
#include "polygon.h"
#include "trimesh2/TriMesh.h"

namespace topomesh
{
	class SlicedMeshLayer
	{
	public:
		SlicedMeshLayer();
		~SlicedMeshLayer();

		Polygons polygons;
		Polygons openPolylines;
		ClipperLib::cInt z;
	};

	class SlicedMesh
	{
	public:
		SlicedMesh();
		~SlicedMesh();

		void save(int index, const std::string& prefix);

		std::vector<SlicedMeshLayer> m_layers;
	};

	class SliceHelper;
	void sliceMeshes_src(const std::vector<trimesh::TriMesh*>& meshes, std::vector<SlicedMesh>& slicedMeshes, std::vector<int>& z);
	void sliceMesh_src(trimesh::TriMesh* mesh, SlicedMesh& slicedMesh, std::vector<int>& z);
	void sliceMesh_src(SlicedMeshLayer& slicedMeshLayer, int z, SliceHelper* helper);
}

#endif // CX_MESHSLICE_1599726705111_H