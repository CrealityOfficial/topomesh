#include "k-means.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"


namespace topomesh {
	kmeansClustering::kmeansClustering(LaplacianMatrix* laplacian, int centre_num)
	{
		Eigen::SparseMatrix<float, Eigen::RowMajor> centre(laplacian->getRow(), laplacian->getCol());
		std::vector<int> centredata;
		for (int i = 0; i < centre_num; i++)
		{

		}
	}
}