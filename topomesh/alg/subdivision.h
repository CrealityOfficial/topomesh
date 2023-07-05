#pragma once
#include "topomesh/data/convert.h"
#include "topomesh/data/mmesht.h"


namespace topomesh {
	void SimpleMidSubdiv(MMeshT* mesh, std::vector<int>& faceindexs);
	void loopSubdiv(MMeshT* mesh, std::vector<int>& faceindexs, int iteration = 1);
}