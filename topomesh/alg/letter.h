#ifndef TOPOMESH_LETTER_1680853426716_H
#define TOPOMESH_LETTER_1680853426716_H
#include "mmesh/trimesh/polygon.h"
#include "topomesh/data/mmesht.h"
#include "topomesh/data/convert.h"

namespace topomesh
{
	void letter(MMeshT* mesh, const ClipperLibXYZ::Paths& paths, const Camera& camera, const LetterParam& param);

	void lettering(trimesh::TriMesh* mesh, ClipperLibXYZ::Paths* paths, std::vector<int>& faceIndexs, std::vector<trimesh::point>& mouseClick, const LetterParam& param);
	void concaveOrConvexOfFaces(MMeshT* mt,std::vector<int>& faces, bool concave=false);
	void splitPoint(MMeshT* mt, MMeshVertex* v, trimesh::point ori);
}

#endif // TOPOMESH_LETTER_1680853426716_H