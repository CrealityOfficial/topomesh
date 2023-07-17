#include "laplacian.h"

namespace topomesh {
	LaplacianMatrix::LaplacianMatrix(MMeshT* mesh)
	{
		if (mesh->VN() != mesh->vertices.size()) mesh->shrinkMesh();
		this->row = mesh->VN();
		this->col = mesh->VN();
		D.resize(this->row, this->col);
		W.resize(this->row, this->col);
		L.resize(this->row, this->col);
		typedef Eigen::Triplet<int> Tr;
		std::vector<Tr> D_tripletList;
		std::vector<Tr> W_tripletList;
		for (MMeshVertex& v : mesh->vertices)if(!v.IsD())
		{
			D_tripletList.push_back(Tr(v.index,v.index,v.connected_vertex.size()));
			for (MMeshVertex* vv : v.connected_vertex)
				W_tripletList.push_back(Tr(v.index, vv->index, 1));
		}
		D.setFromTriplets(D_tripletList.begin(), D_tripletList.end());
		W.setFromTriplets(W_tripletList.begin(), W_tripletList.end());

		L = D - W;
	};

	void LaplacianMatrix::normalzationLaplacian()
	{
		
	}

	void LaplacianMatrix::userDefinedMatrix(const int row, const int col)
	{

	}
}