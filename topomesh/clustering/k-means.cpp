#include "k-means.h"
#include "topomesh/alg/laplacian.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "set"


namespace topomesh {
	kmeansClustering::kmeansClustering(MMeshT* mesh, int centre_num)
	{
		LaplacianMatrix laplacian(mesh);

		Eigen::SparseMatrix<float, Eigen::RowMajor> centre(laplacian.getRow(), laplacian.getCol());
		std::vector<int> centre_data;
		for (int i = 0; i < centre_num; i++)
		{
			int rand = std::rand()%centre_num;
		}
	}
}