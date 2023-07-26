#pragma once
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "topomesh/data/mmesht.h"


namespace topomesh {
	class kmeansClustering {
	public:
		kmeansClustering() {};
		kmeansClustering(Eigen::SparseMatrix<float>* data,int centre_num = 2);
		~kmeansClustering() {};
	private:
		std::vector<std::vector<float>> result;
	};
}