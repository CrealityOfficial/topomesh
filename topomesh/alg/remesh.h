#pragma once
#include "topomesh/data/convert.h"
#include "topomesh/data/mmesht.h"
#include "set"

namespace topomesh {
	void  SimpleRemeshing(MMeshT* mesh,const std::vector<int>& faceindexs,float thershold,std::vector<int>& chunks);
	void  IsotropicRemeshing(MMeshT* mesh, std::vector<int>& faceindexs);
}