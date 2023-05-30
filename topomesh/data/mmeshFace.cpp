#include"mmeshFace.h"

#define FLOATERR 1e-8f

namespace topomesh {

	bool MMeshFace::innerFace(trimesh::point v)
	{
		int RightrayCorssPoint = 0,LeftrayCrossPoint = 0;
		for (int k = 0; k < this->connect_vertex.size(); k++)
		{
			if (std::abs(this->connect_vertex[k]->p.y - this->connect_vertex[k+1]->p.y) < FLOATERR) continue;
			if (v.y < std::min(this->connect_vertex[k]->p.y, this->connect_vertex[k + 1]->p.y))continue;
			if (v.y > std::max(this->connect_vertex[k]->p.y, this->connect_vertex[k + 1]->p.y)) continue;
			double x = (v.y - this->connect_vertex[k]->p.y) * (this->connect_vertex[k + 1]->p.x - this->connect_vertex[k]->p.x) / (this->connect_vertex[k + 1]->p.y - this->connect_vertex[k]->p.y) + this->connect_vertex[k]->p.x;
			if (x - v.x <= 0)
			{
				RightrayCorssPoint++;
			}
			else if (x -v.x >= 0)
			{
				LeftrayCrossPoint++;
			}
		}
		if (RightrayCorssPoint > 0 && LeftrayCrossPoint > 0)
		{
			return true;
		}
		return false;
	}


};