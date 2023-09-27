#ifndef MMESH_MNODE_DUMPLICATE_1622032440408_H
#define MMESH_MNODE_DUMPLICATE_1622032440408_H
#include "trimesh2/TriMesh.h"
#include <unordered_map>
#include "ccglobal/tracer.h"

namespace topomesh
{
    bool dumplicateMesh(trimesh::TriMesh* mesh, ccglobal::Tracer* tracer = nullptr, const float& ratio = 0.3f);
}

#endif // MMESH_MNODE_DUMPLICATE_1622032440408_H