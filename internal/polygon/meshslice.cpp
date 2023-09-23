#include "meshslice.h"
#include "slicepolygonbuilder.h"
#include "slicehelper.h"

#include "ccglobal/tracer.h"

namespace topomesh
{
	SlicedMeshLayer::SlicedMeshLayer()
		:z(0)
	{

	}

	SlicedMeshLayer::~SlicedMeshLayer()
	{

	}

	SlicedMesh::SlicedMesh()
	{

	}

	SlicedMesh::~SlicedMesh()
	{

	}

	void sliceMeshes_src(const std::vector<trimesh::TriMesh*>& meshes, std::vector<SlicedMesh>& slicedMeshes, std::vector<int>& z)
	{
		size_t meshCount = meshes.size();
		if (meshCount > 0)
		{
			if (slicedMeshes.size() != meshCount)
				slicedMeshes.resize(meshCount);

			for (size_t i = 0; i < meshCount; ++i)
				sliceMesh_src(meshes.at(i), slicedMeshes.at(i), z);
		}
	}


	void sliceMesh_src(trimesh::TriMesh* mesh, SlicedMesh& slicedMesh, std::vector<int>& z)
	{
		size_t size = z.size();
		if (size > 0)
		{
			slicedMesh.m_layers.resize(size);

			SliceHelper helper;
			helper.prepare(mesh);

#pragma omp parallel for
			for (int i = 0; i < size; ++i)
				sliceMesh_src(slicedMesh.m_layers.at(i), z.at(i), &helper);
		}
	}


	void sliceMesh_src(SlicedMeshLayer& slicedMeshLayer, int z, SliceHelper* helper)
	{
		SlicePolygonBuilder builder;
		builder.sliceOneLayer_dst(helper, z, &slicedMeshLayer.polygons, &slicedMeshLayer.openPolylines);
	}
}

