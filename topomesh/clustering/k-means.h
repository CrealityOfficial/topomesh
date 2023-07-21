#pragma once
#include "topomesh/data/mmesht.h"


namespace topomesh {
	class kmeansClustering {
	public:
		kmeansClustering() {};
		kmeansClustering(MMeshT* mesh,int centre_num = 2);
		~kmeansClustering() {};
	private:
		std::vector<std::vector<float>> result;
	};
}