#include "laplacian.h"

namespace topomesh {
	LaplacianMatrix::LaplacianMatrix(MMeshT* mesh, bool weight)
	{
		if (mesh->VN() != mesh->vertices.size()) mesh->shrinkMesh();
		if (weight && !mesh->is_HalfEdge()) return;
		this->row = mesh->VN();
		this->col = mesh->VN();
		Eigen::SparseMatrix<float>  D;
		Eigen::SparseMatrix<float>  W;
		Eigen::SparseMatrix<float>  D_sqrt;
		D.resize(this->row, this->col);
		W.resize(this->row, this->col);
		L.resize(this->row, this->col);
		typedef Eigen::Triplet<float> Tr;
		std::vector<Tr> D_tripletList;
		std::vector<Tr> Dsqrt_tripletList;
		std::vector<Tr> W_tripletList;
		for (MMeshVertex& v : mesh->vertices)if(!v.IsD())
		{
			float dre = v.connected_vertex.size();
			D_tripletList.push_back(Tr(v.index,v.index, dre));
			Dsqrt_tripletList.push_back(Tr(v.index, v.index, 1.f / std::sqrt(dre) * 1.f));
			if(weight)
				for (MMeshHalfEdge* he : v.v_mhe)
					W_tripletList.push_back(Tr(v.index, he->edge_vertex.second->index, he->attritube_f));
			else
				for (MMeshVertex* vv : v.connected_vertex)
					W_tripletList.push_back(Tr(v.index, vv->index, 1));
		}
		D.setFromTriplets(D_tripletList.begin(), D_tripletList.end());
		W.setFromTriplets(W_tripletList.begin(), W_tripletList.end());
		D_sqrt.setFromTriplets(Dsqrt_tripletList.begin(), Dsqrt_tripletList.end());
		L = D - W;
		NL = D_sqrt * L * D_sqrt;
	};

	void LaplacianMatrix::normalzationLaplacian()
	{
		
	}

	void LaplacianMatrix::userDefinedMatrix(const int row, const int col)
	{

	}
}