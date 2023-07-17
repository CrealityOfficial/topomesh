#pragma once
#include "topomesh/data/convert.h"
#include "topomesh/data/mmesht.h"

namespace topomesh {
	void  SimpleRemeshing(MMeshT* mesh,const std::vector<int>& faceindexs,float thershold);
	void  IsotropicRemeshing(MMeshT* mesh, std::vector<int>& faceindexs);
}