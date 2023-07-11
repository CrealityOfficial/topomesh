#pragma once
#include "topomesh/data/convert.h"
#include "topomesh/data/mmesht.h"


namespace topomesh {
	void SimpleMidSubdivision(MMeshT* mesh, std::vector<int>& faceindexs);
	void loopSubdivision(MMeshT* mesh, std::vector<int>& faceindexs, int iteration = 1);
	void sqrt3Subdivision(MMeshT* mesh, std::vector<int>& faceindexs);
}