#pragma once
#include "topomesh/data/convert.h"
#include "topomesh/data/mmesht.h"

namespace topomesh {
	void SpectralClusteringCuts(MMeshT* mesh,std::vector<std::vector<int>>& bloack);
}