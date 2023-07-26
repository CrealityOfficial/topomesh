#include "k-means.h"
#include "topomesh/alg/laplacian.h"
#include "set"


namespace topomesh {
	kmeansClustering::kmeansClustering(Eigen::SparseMatrix<float>* data, int centre_num)
	{		
		int data_num = data->rows();
		int data_dim = data->cols();
		//Eigen::SparseMatrix<float> centre(data_num, data_dim);
		std::vector<int> centre_data(centre_num);
		for (int i = 0; i < centre_num; i++)
		{
			int rand = std::rand()%centre_num;

		}
	}
}