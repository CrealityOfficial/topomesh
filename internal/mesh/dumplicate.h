#ifndef MMESH_MNODE_DUMPLICATE_1622032440408_H
#define MMESH_MNODE_DUMPLICATE_1622032440408_H
#include "trimesh2/TriMesh.h"
#include <unordered_map>
#include "ccglobal/tracer.h"

namespace topomesh
{
    bool dumplicateMesh(trimesh::TriMesh* mesh, ccglobal::Tracer* tracer = nullptr, float ratio = 0.3f);
    bool mergeNearPoints(trimesh::TriMesh* mesh, ccglobal::Tracer* tracer = nullptr, float eps = 1E-8F);
}

#endif // MMESH_MNODE_DUMPLICATE_1622032440408_H