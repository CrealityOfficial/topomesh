#pragma once
#include "topomesh/data/convert.h"
#include "topomesh/data/mmesht.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
//#include "Eigen/src/Cholesky/LLT.h"

namespace topomesh {

	class LaplacianMatrix
	{
	public:
		LaplacianMatrix() {};
		LaplacianMatrix(MMeshT* mesh);
		~LaplacianMatrix() {};

		void userDefinedMatrix(const int row, const int col);
		void normalzationLaplacian();
	private:
		int row;
		int col;
		Eigen::SparseMatrix<int>  D;
		Eigen::SparseMatrix<int>  W;
		Eigen::SparseMatrix<int>  L;
		Eigen::SparseMatrix<int>  NL;
	};
}