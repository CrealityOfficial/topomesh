#ifndef TOPOMESH_LETTER_1680853426716_H
#define TOPOMESH_LETTER_1680853426716_H
#include "clipperxyz/clipper.hpp"
#include "mmesh/trimesh/polygon.h"
#include "topomesh/data/mmesht.h"

namespace topomesh
{
	void lettering(trimesh::TriMesh* mesh, ClipperLibXYZ::Paths* paths, std::vector<int>& faceIndexs, std::vector<trimesh::point>& mouseClick, int deep);
	void concaveOrConvexOfFaces(MMeshT* mt,std::vector<int>& faces, bool concave=false);
	void splitPoint(MMeshT* mt, MMeshVertex* v, trimesh::point ori);
}

#endif // TOPOMESH_LETTER_1680853426716_H