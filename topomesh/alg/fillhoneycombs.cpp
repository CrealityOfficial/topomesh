#include "fillhoneycombs.h"
#include "topomesh/alg/letter.h"
#include "topomesh/data/mmesht.h"
#include "topomesh/honeycomb/Matrix.h"
#include "topomesh/honeycomb/Polyline.h"
#include "topomesh/honeycomb/HoneyComb.h"
#ifndef EPS
#define EPS 1e-8f
#endif // !EPS

namespace topomesh {
    struct honeyLetterOpt {
        std::vector<int>bottom;
        std::vector<int>others;
        struct hexagon {
            double radius = 1.0;
            std::vector<trimesh::vec3>borders;
        };
        std::vector<hexagon>hexgons;
    };

    honeycomb::Mesh ConstructFromTriMesh(const trimesh::TriMesh* trimesh)
    {
        honeycomb::Mesh mesh;
        auto& faces = mesh.GetFaces();
        auto& points = mesh.GetPoints();
        auto& normals = mesh.GetFaceNormals();
        auto& indexs = mesh.GetFaceVertexAdjacency();
        const auto& vertices = trimesh->vertices;
        const auto& faceVertexs = trimesh->faces;
        const auto& norms = trimesh->normals;
        faces.resize(faceVertexs.size());
        std::iota(faces.begin(), faces.end(), 0);
        points.reserve(vertices.size());
        for (const auto& v : vertices) {
            honeycomb::Point p(v[0], v[1], v[2]);
            points.emplace_back(std::move(p));
        }
        indexs.reserve(faceVertexs.size());
        for (const auto& fa : faceVertexs) {
            std::vector<int> vertexs;
            vertexs.reserve(3);
            vertexs.emplace_back(fa[0]);
            vertexs.emplace_back(fa[1]);
            vertexs.emplace_back(fa[2]);
            indexs.emplace_back(vertexs);
        }
        mesh.GenerateFaceNormals();
        //mesh.WriteSTLFile("trimesh转mesh");
        return mesh;
    }

    trimesh::TriMesh ConstructFromHoneyMesh(const honeycomb::Mesh& honeyMesh)
    {
        trimesh::TriMesh triMesh;
        const auto& bound = honeyMesh.Bound();
        const auto& points = honeyMesh.GetPoints();
        const auto& indexs = honeyMesh.GetFaceVertexAdjacency();
        triMesh.vertices.reserve(points.size());
        triMesh.faces.reserve(indexs.size());
        for (const auto& p : points) {
            const auto& pt = trimesh::vec3((float)p.x, (float)p.y, (float)p.z);
            triMesh.vertices.emplace_back(pt);
        }
        for (const auto& f : indexs) {
            const auto& fa = trimesh::vec3(f[0], f[1], f[2]);
            triMesh.faces.emplace_back(fa);
        }
        const auto& min = bound.Min();
        const auto& max = bound.Max();
        const auto& p1 = trimesh::vec3((float)min.x, (float)min.y, (float)min.z);
        const auto& p2 = trimesh::vec3((float)max.x, (float)max.y, (float)max.z);
        triMesh.bbox = trimesh::box3(p1, p2);
        return triMesh;
    }

    void GenerateHexagons(honeycomb::Mesh& honeyMesh, const HoneyCombParam& honeyparams, honeyLetterOpt& letteropts)
    {
        //第1步，寻找底面（最大平面）朝向
        std::vector<int>bottomFaces;
        honeycomb::Point dir = honeyMesh.FindBottomDirection(&bottomFaces);
        honeyMesh.Rotate(dir, honeycomb::Point(0, 0, -1));
        honeycomb::BoundBox bound = honeyMesh.Bound();
        honeyMesh.Translate(-bound.Min());
        letteropts.bottom.resize(bottomFaces.size());
        letteropts.bottom.assign(bottomFaces.begin(), bottomFaces.end());
        const auto& honeyFaces = honeyMesh.GetFaces();
        std::sort(bottomFaces.begin(), bottomFaces.end());
        std::vector<int> otherFaces(honeyFaces.size() - bottomFaces.size());
        std::set_difference(honeyFaces.begin(), honeyFaces.end(), bottomFaces.begin(), bottomFaces.end(), otherFaces.begin());
        letteropts.others = std::move(otherFaces);
        //第2步，平移至xoy平面后底面整平
        honeyMesh.FlatBottomSurface(&bottomFaces);
        honeycomb::Mesh cutMesh;
        cutMesh.MiniCopy(honeyMesh);
        //第3步，剪切掉底面得边界轮廓
        cutMesh.DeleteFaces(bottomFaces, true);
        cutMesh.WriteSTLFile("删除底面后剩余面片");
        std::vector<int> edges;
        cutMesh.SelectIndividualEdges(edges);
        //第4步，底面边界轮廓点索引排序
        std::vector<std::vector<int>>sequentials;
        cutMesh.GetSequentialPoints(edges, sequentials);
        //第5步，构建底面平面多边形轮廓
        std::vector<honeycomb::Poly2d> polys;
        polys.reserve(sequentials.size());
        const auto & points = cutMesh.GetPoints();
        for (const auto& seq : sequentials) {
            std::vector<honeycomb::Point2d> border;
            border.reserve(seq.size());
            for (const auto& v : seq) {
                const auto& p = points[v];
                border.emplace_back(p.x, p.y);
            }
            honeycomb::Poly2d poly(border);
            polys.emplace_back(poly);
        }
        honeycomb::ExPoly2d boundarys(polys);
        //第6步，生成六边形网格阵列
        honeycomb::HexagonOpt opt;
        opt.nestWidth = honeyparams.nestWidth;
        opt.hexagonRadius = honeyparams.honeyCombRadius;
        opt.shellThickness = honeyparams.shellThickness;
        const auto& edgelist = boundarys.Edges();
        std::vector<std::vector<honeycomb::Point2d>> hexborders;
        std::vector<honeycomb::Hexagon>& hexagons = boundarys.GenerateHexagons(opt, hexborders);
        letteropts.hexgons.reserve(hexagons.size());
        for (const auto& hex : hexagons) {
            const auto& border = hex.Border();
            std::vector<trimesh::vec3> triborders;
            triborders.reserve(border.size());
            for (const auto& p : border) {
                const auto& point = trimesh::vec3(p.x, p.y, 0);
                triborders.emplace_back(std::move(point));
            }
            letteropts.hexgons.emplace_back(std::move(honeyLetterOpt::hexagon{ opt.hexagonRadius, triborders }));
        }
        hexborders.insert(hexborders.end(), edgelist.begin(), edgelist.end());
        const std::vector<std::vector<honeycomb::Point>>& result = WrapToPoints(hexborders);
        honeycomb::Mesh hexagonMesh;
        honeycomb::SaveEdgesToMesh(result, hexagonMesh, 0.04, 5);
        hexagonMesh.WriteSTLFile("六边形网格");
    }

    trimesh::TriMesh* generateHoneyCombs(trimesh::TriMesh* trimesh, const HoneyCombParam& honeyparams,
        ccglobal::Tracer* tracer, HoneyCombDebugger* debugger)
	{ 
        honeyLetterOpt letterOpts;
        honeycomb::Mesh& inputMesh = ConstructFromTriMesh(trimesh);
        GenerateHexagons(inputMesh, honeyparams, letterOpts);
        trimesh::TriMesh& mesh = ConstructFromHoneyMesh(inputMesh);
		return nullptr;
	}

	void GenerateHoneyCombs(const trimesh::TriMesh* mesh, trimesh::TriMesh& resultmesh, const TriPolygon& poly, trimesh::vec3 axisDir ,
		trimesh::vec2 arrayDir, double honeyCombRadius , double nestWidth , double shellThickness)
	{
		trimesh::TriMesh* newMesh=new trimesh::TriMesh();
		*newMesh = *mesh;
		newMesh->need_normals();
		newMesh->need_adjacentfaces();
		newMesh->need_neighbors();
		std::vector<int> upfaceindex;
		std::vector<int> botfaceindex;
		for (int fi = 0; fi < newMesh->faces.size(); fi++)
		{
			trimesh::point n = trimesh::normalized(newMesh->trinorm(fi));
			if (std::abs(n.z+1)<10*EPS)
				botfaceindex.push_back(fi);
			else
				upfaceindex.push_back(fi);
		}
		MMeshT mt(newMesh);		
		//std::vector<std::vector<trimesh::vec2>> lines;
		//embedingAndCutting(&mt, lines, botfaceindex);
		std::vector<int> honeycombs = {21463,20775,21464,22129,22792,22128};
		
		for (int i = 0; i < honeycombs.size(); i++)
		{
			mt.faces[honeycombs[i]].V0(0)->SetS();
			mt.faces[honeycombs[i]].V1(0)->SetS();
			mt.faces[honeycombs[i]].V2(0)->SetS();
			mt.faces[honeycombs[i]].SetS();
		}		
		for (int vi = 0; vi < mt.vertices.size(); vi++)
		{
			if (mt.vertices[vi].IsS())
			{				
				for (MMeshVertex* v : mt.vertices[vi].connected_vertex)
				{
					if (!v->IsS())
					{
						mt.vertices[vi].SetV(); break;
					}
				}
			}
		}
		std::vector<std::pair<int, int>> vd;
		findNeighVertex(newMesh, upfaceindex, honeycombs,vd);
		for (int i = 0; i < vd.size(); i++)
		{
			float z = mt.vertices[vd[i].second].p.z-2.0f;
			if(mt.vertices[vd[i].first].IsV())
				splitPoint(&mt, &mt.vertices[vd[i].first],trimesh::point(0,0,z));
			else
				mt.vertices[vd[i].first].p += trimesh::point(0, 0, z);
		}
		//concaveOrConvexOfFaces(&mt, honeycombs, true, 20.f);
		mt.mmesh2trimesh(newMesh);
		newMesh->write("botmesh.ply");	
	}

	void findNeighVertex(trimesh::TriMesh* mesh, const std::vector<int>& upfaceid,const std::vector<int>& botfaceid, std::vector<std::pair<int, int>>& vertex_distance)
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

		for (int vi = 0; vi < compara_vertex.size(); vi++)
		{
			int xi = (mesh->vertices[compara_vertex[vi]].x - begin_x) / length_x;
			int yi = (mesh->vertices[compara_vertex[vi]].y - begin_y) / length_y;
			if (mapp[xi][yi] != std::numeric_limits<float>::min())
			{
				std::cout << "mapind :" << mapind[xi][yi] << " z : " << mapp[xi][yi] << "\n";
				vertex_distance.push_back(std::make_pair(compara_vertex[vi], mapind[xi][yi]));
			}
			else {//----------------锟斤拷锟斤拷锟秸白憋拷锟斤拷锟斤拷围锟斤拷元锟斤拷
				float max = std::numeric_limits<float>::min();
				int xi_x, yi_y;
				for (int i = 1; i < width; i++)
				{
					for (int ii = xi - i; ii <= xi + i; ii++)
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
						vertex_distance.push_back(std::make_pair(compara_vertex[vi], mapind[xi_x][yi_y]));
						break;
					}
				}
			}
		}		
		
#endif
	}

	//只支锟斤拷锟斤拷锟斤拷锟轿匡拷锟劫诧拷锟斤拷
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