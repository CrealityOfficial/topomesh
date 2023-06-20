#include "fillhoneycombs.h"

#define EPS 1e-8f

namespace topomesh {
	void findNeighVertex(trimesh::TriMesh* mesh, const std::vector<int>& upfaceid,const std::vector<int>& botfaceid)
	{
		std::vector<int> vertexidx;
		std::vector<int> compara_vertex;
		for (int i = 0; i < upfaceid.size(); i++)
		{
			vertexidx.push_back(mesh->faces[upfaceid[i]][0]);
			vertexidx.push_back(mesh->faces[upfaceid[i]][1]);
			vertexidx.push_back(mesh->faces[upfaceid[i]][2]);
		}
		std::sort(vertexidx.begin(), vertexidx.end());
		std::vector<int>::iterator itr = std::unique(vertexidx.begin(), vertexidx.end());
		vertexidx.resize(std::distance(vertexidx.begin(), itr));

		for (int i = 0; i < botfaceid.size(); i++)
		{
			compara_vertex.push_back(mesh->faces[botfaceid[i]][0]);
			compara_vertex.push_back(mesh->faces[botfaceid[i]][1]);
			compara_vertex.push_back(mesh->faces[botfaceid[i]][2]);
		}
		std::sort(compara_vertex.begin(), compara_vertex.end());
		std::vector<int>::iterator it = std::unique(compara_vertex.begin(), compara_vertex.end());
		compara_vertex.resize(std::distance(compara_vertex.begin(), it));
#if 0
		const int K = 10;
		std::list<Point> points;
		for (int i = 0; i < vertexidx.size(); i++)
		{
			points.push_back(Point(mesh->vertices[vertexidx[i]].x, mesh->vertices[vertexidx[i]].y));
		}
		Tree kdtree(points.begin(), points.end());
		Point query(0, 0);
		Neighbor_search search(kdtree, query, K);
		for (Neighbor_search::iterator it = search.begin(); it != search.end(); it++)
			std::cout << "point :" << it->first << " distance^2 : " << it->second << "\n";
#else
		const int width = 100, height = 100;
		std::vector<std::vector<float>> mapp(width, std::vector<float>(height, std::numeric_limits<float>::min()));
		std::vector<std::vector<int>> mapind(width, std::vector<int>(height, -1));
		mesh->need_bbox();
		float length_x = (mesh->bbox.max.x - mesh->bbox.min.x) / (width * 1.0f);
		float length_y = (mesh->bbox.max.y - mesh->bbox.min.y) / (height * 1.0f);

		float begin_x = mesh->bbox.min.x;
		float begin_y = mesh->bbox.min.y;

		for (int i = 0; i < vertexidx.size(); i++)
		{
			float x = mesh->vertices[vertexidx[i]].x - begin_x;
			float y = mesh->vertices[vertexidx[i]].y - begin_y;
			int x_ind = x / length_x;
			int y_ind = y / length_y;
			if (x_ind == width)
				x_ind--;
			if (y_ind == height)
				y_ind--;
			if (mesh->vertices[vertexidx[i]].z > mapp[x_ind][y_ind])
			{
				mapp[x_ind][y_ind] = mesh->vertices[vertexidx[i]].z;
				//std::cout << " x_ind :" << x_ind << " y_ind :" << y_ind << " z :" << mesh->vertices[vertexidx[i]].z << "\n";
				mapind[x_ind][y_ind] = vertexidx[i];
			}
		}

		/*for (int vi = 0; vi < compara_vertex.size(); vi++)
		{
			int xi = (mesh->vertices[compara_vertex[vi]].x - begin_x) / length_x;
			int yi = (mesh->vertices[compara_vertex[vi]].y - begin_y) / length_y;
			if (mapp[xi][yi] != std::numeric_limits<float>::max())
			{

			}
		}*/
		trimesh::vec2 v(-15.458, -16.137);
		int xi = (v.x - begin_x) / length_x;
		int yi = (v.y - begin_y) / length_y;
		if (mapp[xi][yi] != std::numeric_limits<float>::min())
		{
			std::cout << "mapind :" << mapind[xi][yi] << " z : " << mapp[xi][yi] << "\n";
		}
		else {//----------------遇到空白遍历周围单元格
			float max = std::numeric_limits<float>::min();
			int xi_x, yi_y;
			for (int i = 1; i < width; i++)
			{
				for(int ii=xi-i;ii<=xi+i;ii++)
					for (int jj = yi - i; jj <= yi + i; jj++)
					{
						if (std::abs(ii - xi) != i && std::abs(jj - yi) != i) continue;
						if (mapp[ii][jj] > max)
						{
							max = mapp[ii][jj];
							xi_x = ii; yi_y = jj;
						}
					}
				if (max != std::numeric_limits<float>::min())
				{
					std::cout << "mapind :" << mapind[xi_x][yi_y] << " z : " << max << "\n";
					break;
				}
			}
		}

		

#endif
	}

	//只支持六边形快速查找
	void innerHex(MMeshT* mesh, std::vector<std::vector<trimesh::vec2>>& poly, std::vector<int>& inFace, std::vector<int>& outFace, float len)
	{
		std::vector<trimesh::vec2> center(poly.size());
		trimesh::vec2 start(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
		trimesh::vec2 top(std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
		for (int i = 0; i < poly.size(); i++)
		{
			trimesh::vec2 v(0,0);
			for (int j = 0; j < poly[i].size(); j++)
			{
				v += poly[i][j];
				if (poly[i][j].x < start.x)
					start.x = poly[i][j].x;
				if (poly[i][j].y < start.y)
					start.y = poly[i][j].y;
				if (poly[i][j].x > top.x)
					top.x = poly[i][j].x;
				if (poly[i][j].y > top.y)
					top.y = poly[i][j].y;
			}
			v = v / poly[i].size() * 1.0f;
			center[i] = v;
		}
		int width = (top.x - start.x) / len + 1;
		int height = (top.y - start.y) / (std::sqrt(3) * len) + 1;

		std::vector<std::vector<std::pair<int,int>>> block(center.size());
		for (int i = 0; i < center.size(); i++)
		{
			trimesh::vec2 p = trimesh::vec2(center[i].x + 0.1 * len, center[i].y + 0.1*len) - start;
			int x = p.x / len;
			int y = p.y / (std::sqrt(3) * len);
			block[i].push_back(std::make_pair(x,y));
			block[i].push_back(std::make_pair(x + 1, y)); 
			block[i].push_back(std::make_pair(x - 1, y)); 
			block[i].push_back(std::make_pair(x - 2, y)); 
			block[i].push_back(std::make_pair(x, y - 1)); 
			block[i].push_back(std::make_pair(x + 1, y - 1)); 
			block[i].push_back(std::make_pair(x - 1, y - 1)); 
			block[i].push_back(std::make_pair(x - 2, y - 1)); 
		}
		std::vector<std::vector<std::vector<int>>> faceBlock(width,std::vector<std::vector<int>>(height));
		for (int fi = 0; fi < inFace.size(); fi++)
		{
			MMeshFace& f = mesh->faces[inFace[fi]];
			trimesh::point c = (f.V0(0)->p + f.V0(1)->p + f.V0(2)->p) / 3.0;
			int xi = (c.x - start.x) / len;
			int yi = (c.y - start.y) / (std::sqrt(3) * len);
			faceBlock[xi][yi].push_back(fi);
		}
		std::vector<std::vector<int>> polyInnerFaces(block.size());
		for (int i = 0; i < block.size(); i++)
		{
			for (int j = 0; j < block[i].size(); j++)
			{
				polyInnerFaces[i].insert(polyInnerFaces[i].end(), faceBlock[block[i][j].first][block[i][j].second].begin(), faceBlock[block[i][j].first][block[i][j].second].end());
			}
		}

		for (int i = 0; i < polyInnerFaces.size(); i++)
		{
			for (int j = 0; j < polyInnerFaces[i].size(); j++)
			{
				int leftCross = 0;int rightCross = 0;
				MMeshFace& f = mesh->faces[polyInnerFaces[i][j]];
				trimesh::point c= (f.V0(0)->p + f.V0(1)->p + f.V0(2)->p) / 3.0;
				for (int k = 0; k < poly[i].size(); k++)
				{
					if (std::abs(poly[i][k].y - poly[i][(k + 1) % poly[i].size()].y) < EPS) continue;
					if (c.y < std::min(poly[i][k].y, poly[i][(k + 1) % poly[i].size()].y)) continue;
					if (c.y > std::max(poly[i][k].y, poly[i][(k + 1) % poly[i].size()].y)) continue;
					double x = (c.y - poly[i][k].y) * (poly[i][(k + 1) % poly[i].size()].x - poly[i][k].x) / (poly[i][(k + 1) % poly[i].size()].y - poly[i][k].y) + poly[i][k].x;
					if (x - c.x <= 0)
						rightCross++;
					else if (x - c.x >= 0)
						leftCross++;
				}
				if (leftCross > 0 && rightCross > 0)
					outFace.push_back(polyInnerFaces[i][j]);
			}
		}

	}
}