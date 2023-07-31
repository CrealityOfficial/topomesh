#pragma once
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "topomesh/data/mmesht.h"


namespace topomesh {
	//template <typename T>
	class kmeansClustering {
	public:
		kmeansClustering() {};
		//template <typename T>
		kmeansClustering(Eigen::MatrixXf* data,int centre_num = 2);//Eigen::SparseMatrix<float>
		~kmeansClustering() {};
		std::vector<std::vector<float>>* get_result() { return &result; }
	private:
		std::vector<std::vector<float>> result;
	};
}