#pragma once
#include "topomesh/alg/laplacian.h"


namespace topomesh {
	class kmeansClustering {
	public:
		kmeansClustering() {};
		kmeansClustering(LaplacianMatrix* laplacian,int centre_num = 2);
		~kmeansClustering() {};
	private:
		std::vector<std::vector<float>> result;
	};
}