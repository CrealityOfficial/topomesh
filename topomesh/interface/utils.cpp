#include "utils.h"
#include "topomesh/alg/earclipping.h"
#include <queue>

namespace topomesh
{
	void findNeignborFacesOfSameAsNormal(trimesh::TriMesh* trimesh, int indicate, float angle_threshold, std::vector<int>& faceIndexs)
	{
		if (!trimesh)
			return;

		trimesh->need_across_edge();
		faceIndexs.push_back(indicate);
		trimesh::point normal = (trimesh->vertices[trimesh->faces[indicate][1]] - trimesh->vertices[trimesh->faces[indicate][0]]) %
			(trimesh->vertices[trimesh->faces[indicate][2]] - trimesh->vertices[trimesh->faces[indicate][0]]);
		trimesh::normalize(normal);
		std::vector<int> vis(trimesh->faces.size(), 0);
		std::queue<int> facequeue;
		facequeue.push(indicate);
		while (!facequeue.empty())
		{
			vis[facequeue.front()] = 1;
			for (int i = 0; i < trimesh->across_edge[facequeue.front()].size(); i++)
			{
				int face = trimesh->across_edge[facequeue.front()][i];
				if (vis[face]) continue;
				trimesh::point v1 = trimesh->vertices[trimesh->faces[face][1]] - trimesh->vertices[trimesh->faces[face][0]];
				trimesh::point v2 = trimesh->vertices[trimesh->faces[face][2]] - trimesh->vertices[trimesh->faces[face][0]];
				trimesh::point nor = v1 % v2;
				float arc = trimesh::normalized(nor) ^ normal;
				arc = arc >= 1.f ? 1.f : arc;
				arc = arc <= -1.f ? -1.f : arc;
				float ang = std::acos(arc) * 180 / M_PI;
				if (ang < angle_threshold)
				{
					vis[face] = 1;
					facequeue.push(face);
					faceIndexs.push_back(face);
				}
			}
			facequeue.pop();
		}
	}

	void findBoundary(trimesh::TriMesh* trimesh)
	{

	}

	void triangulate(trimesh::TriMesh* trimesh, std::vector<int>& sequentials)
	{
		std::vector<std::pair<trimesh::point, int>> lines;
		for (int i = 0; i < sequentials.size(); i++)
		{
			lines.push_back(std::make_pair(trimesh->vertices[sequentials[i]], sequentials[i]));
		}
		topomesh::EarClipping earclip(lines);
		std::vector<trimesh::ivec3> result = earclip.getResult();
		for (int fi = 0; fi < result.size(); fi++)
		{
			trimesh->faces.push_back(result[fi]);
		}
	}
}