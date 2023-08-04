#include"mmeshFace.h"

#define FLOATERR 1e-8f

namespace topomesh {

	bool MMeshFace::innerFace(trimesh::point v)
	{
		int RightrayCorssPoint = 0,LeftrayCrossPoint = 0;
		for (size_t k = 0; k < this->connect_vertex.size(); k++)
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

	std::pair<int, int> MMeshFace::commonEdge(MMeshFace* f)
	{
		std::vector<int> id;
		for (int i = 0; i < this->connect_vertex.size(); i++)
		{
			for (int j = 0; j < f->connect_vertex.size(); j++)
			{
				if (this->connect_vertex[i]->index == f->connect_vertex[j]->index)
				{
					id.push_back(this->connect_vertex[i]->index);
				}
			}
		}
		return std::make_pair(id[0], id[1]);
	}

	float MMeshFace::dihedral(MMeshFace* f)
	{
		float arc = this->normal ^ f->normal;
		arc = arc > 1.f ? 1 : arc;
		arc = arc < -1.f ? -1.f : arc;
		float ang = std::acos(arc) * 180 / M_PI;
		return ang;
		/*if (arc < 0)
			return 180 + ang;
		else
			return 180 - ang;*/
	}
};