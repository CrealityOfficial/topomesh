#include "subdivision.h"

namespace topomesh {
	void SimpleMidSubdiv(MMeshT* mesh, std::vector<int>& faceindexs)
	{
		for (int fi : faceindexs)
		{
			trimesh::point c = (mesh->faces[fi].V0(0)->p + mesh->faces[fi].V1(0)->p + mesh->faces[fi].V2(0)->p)/3.0;
			mesh->appendVertex(c);
			mesh->deleteFace(fi);
			mesh->appendFace(mesh->faces[fi].V0(0)->index, mesh->faces[fi].V0(1)->index, mesh->vertices.back().index);
			mesh->appendFace(mesh->faces[fi].V0(1)->index, mesh->faces[fi].V0(2)->index, mesh->vertices.back().index);
			mesh->appendFace(mesh->faces[fi].V0(2)->index, mesh->faces[fi].V0(0)->index, mesh->vertices.back().index);
		}
	}

	void loopSubdiv(MMeshT* mesh, std::vector<int>& faceindexs, int iteration)
	{
		auto bate = [&](int n)->float {
			float alpth = 3.0f / 8.0f + 1.f / 4.f * std::cos(2.0f * M_PI / (n * 1.0f));
			return (5.0f / 8.0f - alpth * alpth)/(n*1.0f);
		};

		for (int fi : faceindexs)
		{
			mesh->faces[fi].SetS();

		}
	}

}