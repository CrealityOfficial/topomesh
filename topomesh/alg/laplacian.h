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
		LaplacianMatrix(MMeshT* mesh,bool weight=false);
		~LaplacianMatrix() {};

		void userDefinedMatrix(const int row, const int col);
		void normalzationLaplacian();

		inline const Eigen::SparseMatrix<float, Eigen::RowMajor>* getLaplacianPtr() { return &L; };
		inline const Eigen::SparseMatrix<float, Eigen::RowMajor>* getNormalizationLaplacianPtr() { return &NL; };

		inline const int getRow() { return row; }
		inline const int getCol() { return col; }
	private:
		int row;
		int col;
		Eigen::SparseMatrix<float, Eigen::RowMajor>  L;
		Eigen::SparseMatrix<float, Eigen::RowMajor>  NL;
		Eigen::SparseMatrix<float, Eigen::RowMajor>  Lsym;
		Eigen::SparseMatrix<float, Eigen::RowMajor>  Lrw;
	};
}