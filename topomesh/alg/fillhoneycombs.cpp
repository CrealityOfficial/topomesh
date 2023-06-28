#include "fillhoneycombs.h"
#include "topomesh/math/polygon.h"
#include "topomesh/math/AABB.h"
#include "topomesh/math/SVG.h"
#include "topomesh/data/mmesht.h"
#include "topomesh/alg/letter.h"
#include "topomesh/honeycomb/Matrix.h"
#include "topomesh/honeycomb/Polyline.h"
#include "topomesh/honeycomb/HoneyComb.h"
#ifndef EPS
#define EPS 1e-8f
#endif // !EPS

namespace topomesh {
    struct honeyLetterOpt {
        std::vector<int>bottom; ///<底面大块平面的面片索引
        std::vector<int>others; ///<去掉底面后其余面片索引（已保留原模型对应的索引）
        /*每个六角网格结构体，
        包含收缩前的网格边长以及
        收缩后网格边界点坐标
        （不一定6个点，默认xoy坐标系逆时针排序）
        */
        struct hexagon {
            double radius = 1.0;
            std::vector<trimesh::vec3>borders;
        };
        std::vector<hexagon>hexgons; ///
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
        //mesh.WriteSTLFile("trimeshתmesh");
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

    void GenerateExHexagons(honeycomb::Mesh& honeyMesh, const HoneyCombParam& honeyparams, honeyLetterOpt& letterOpts, HoneyCombDebugger* debugger = nullptr)
    {
        //拷贝一份数据
        honeycomb::Mesh cutMesh;
        cutMesh.MiniCopy(honeyMesh);
        //第3步，剪切掉底面得边界轮廓
        cutMesh.DeleteFaces(letterOpts.bottom, true);
        //cutMesh.WriteSTLFile("删除底面后剩余面片");
        std::vector<int> edges;
        cutMesh.SelectIndividualEdges(edges);
        //第4步，底面轮廓边界所有边排序
        std::vector<std::vector<int>>sequentials;
        cutMesh.GetSequentialPoints(edges, sequentials);
        //第5步，构建底面边界轮廓多边形
        std::vector<honeycomb::Poly2d> polys;
        polys.reserve(sequentials.size());
        const auto& points = cutMesh.GetPoints();
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
        if (debugger) {
            //显示底面边界轮廓多边形
            TriPolygons polygons;
            polygons.reserve(polys.size());
            for (const auto& poly : polys) {
                const auto& points = poly.Points();
                TriPolygon polygon;
                polygon.reserve(points.size());
                for (const auto& p : points) {
                    polygon.emplace_back(trimesh::vec3((float)p.x, (float)p.y, 0));
                }
                polygons.emplace_back(std::move(polygon));
            }
            debugger->onGenerateInfillPolygons(polygons);
        }
        //第6步，底面边界轮廓抽壳
        const double resolution = honeyparams.resolution;
        const double radius = honeyparams.honeyCombRadius;
        const double thickness = honeyparams.shellThickness;
        const double side = radius - honeyparams.nestWidth / SQRT3;
        cxutil::Polygons bpolygons; ///<底面边界轮廓多边形
        cxutil::Polygons mpolygons; ///<地面边界抽壳多边形
        {
            for (const auto& poly : polys) {
                const auto& points = poly.Points();
                ClipperLib::Path path;
                for (const auto& p : points) {
                    const auto& x = int(p.x / resolution);
                    const auto & y = int(p.y / resolution);
                    path.emplace_back(ClipperLib::IntPoint(x, y));
                }
                bpolygons.paths.emplace_back(std::move(path));
            }
            cxutil::Polygons cpolygons(bpolygons);
            mpolygons = cpolygons.offset(-int(thickness / resolution));
        }
        if (debugger) {
            //显示底面轮廓抽壳多边形
            TriPolygons polygons;
            ClipperLib::Paths paths = mpolygons.paths;
            polygons.reserve(paths.size());
            for (const auto& path : paths) {
                TriPolygon poly;
                poly.reserve(path.size());
                for (const auto& p : path) {
                    poly.emplace_back(trimesh::vec3(p.X * resolution, p.Y * resolution, 0));
                }
                polygons.emplace_back(std::move(poly));
            }
            debugger->onGenerateBottomPolygons(polygons);
        }
        //第7步，底面边界抽壳同时，生成六角网格阵列
        auto& min = boundarys.Bound().Min();
        auto& max = boundarys.Bound().Max();
        if (honeyparams.polyline) {
            min.x = std::numeric_limits<float>::max();
            min.y = std::numeric_limits<float>::max();
            max.x = std::numeric_limits<float>::lowest();
            max.y = std::numeric_limits<float>::lowest();
            const auto& poly = *honeyparams.polyline;
            for (const auto& p : poly) {
                if (p.x < min.x) min.x = p.x;
                if (p.y < min.y) min.y = p.y;
                if (p.x > max.x) max.x = p.x;
                if (p.y > max.y) max.y = p.y;
            }
        }
        const auto& xdist = 3.0 / 2.0 * radius;
        const auto& ydist = SQRT3 / 2.0 * radius;
        trimesh::vec3& origin = trimesh::vec3(min.x, max.y, 0);

        HexagonArrayParam hexagonparams;
        hexagonparams.pos = origin;
        hexagonparams.radius = radius;
        hexagonparams.nestWidth = honeyparams.nestWidth;
        hexagonparams.ncols = std::ceil((max.x - min.x) / xdist);
        hexagonparams.nrows = (std::ceil((max.y - min.y) / ydist) + 1) / 2;
        TriPolygons polygons = GenerateHexagons(hexagonparams);
        if (debugger) {
            //显示六角网格阵列
            debugger->onGenerateBottomPolygons(polygons);
        }
        //第8步，计算网格阵列与抽壳轮廓交集
        cxutil::Polygons hpolygons; ///<整六角网格多边形
        cxutil::Polygons opolygons; ///<与抽壳边界求交的六角网格多边形
        {
            hpolygons.paths.reserve(polygons.size());
            for (const auto& poly : polygons) {
                ClipperLib::Path path;
                path.reserve(poly.size());
                for (const auto& p : poly) {
                    const auto& x = int(p.x / resolution);
                    const auto & y = int(p.y / resolution);
                    path.emplace_back(ClipperLib::IntPoint(x, y));
                }
                hpolygons.paths.emplace_back(std::move(path));
            }
        }
        cxutil::Polygons ipolygons = mpolygons.intersectionPolyLines(hpolygons);
        if (debugger) {
            //显示初始求交结果
            TriPolygons polygons;
            ClipperLib::Paths paths = ipolygons.paths;
            polygons.reserve(paths.size());
            for (const auto& path : paths) {
                TriPolygon poly;
                poly.reserve(path.size());
                for (const auto& p : path) {
                    poly.emplace_back(trimesh::vec3(p.X * resolution, p.Y * resolution, 0));
                }
                polygons.emplace_back(std::move(poly));
            }
            debugger->onGenerateBottomPolygons(polygons);
            /*std::string svgFile = "bottom_fill_in_hexagon.svg";
            cxutil::AABB box(bpolygons.min(), bpolygons.max());
            cxutil::SVG svg(svgFile, box, 0.01);
            svg.writePolygons(bpolygons, cxutil::SVG::Color::RED, 3.0);
            svg.writePolygons(mpolygons, cxutil::SVG::Color::GREEN, 3.0);
            svg.writePolygons(ipolygons, cxutil::SVG::Color::BLUE, 3.0);*/
        }
        //第9步，筛选符合需求的网格多边形
        TriPolygons outpolygons; ///<求交后网格多边形序列
        const double hexarea = 3.0 / 2.0 * SQRT3 * std::pow(side / resolution, 2);
        for (const auto& path : ipolygons.paths) {
            cxutil::Polygons polys;
            TriPolygon tripolygon;
            polys.paths.emplace_back(std::move(path));
            if (polys.area() < 0) {
                polys.clear();
                ClipperLib::Path tmp = std::move(path);
                std::reverse(tmp.begin(), tmp.end());
                polys.paths.emplace_back(std::move(path));
            }
            if (polys.area() >= 1.0 / 3.0 * hexarea) {
                opolygons.paths.emplace_back(std::move(path));
                for (const auto& p : path) {
                    const auto& point = trimesh::vec3(p.X * resolution, p.Y * resolution, 0);
                    tripolygon.emplace_back(std::move(point));
                }
                outpolygons.emplace_back(std::move(tripolygon));
            }
        }
        if (debugger) {
            //显示最终符合需求的网格多边形
            debugger->onGenerateBottomPolygons(outpolygons);
            /*std::string svgFile = "bottom_in_result_hexagon.svg";
            cxutil::AABB box(bpolygons.min(), bpolygons.max());
            cxutil::SVG svg(svgFile, box, 0.01);
            svg.writePolygons(bpolygons, cxutil::SVG::Color::RED, 3.0);
            svg.writePolygons(mpolygons, cxutil::SVG::Color::GREEN, 3.0);
            svg.writePolygons(opolygons, cxutil::SVG::Color::BLUE, 3.0);*/
        }
        letterOpts.hexgons.reserve(opolygons.size());
        for (const auto& poly : outpolygons) {
            letterOpts.hexgons.emplace_back(std::move(honeyLetterOpt::hexagon{ radius, poly }));
        }

    }
    trimesh::TriMesh* generateHoneyCombs(trimesh::TriMesh* trimesh, const HoneyCombParam& honeyparams,
        ccglobal::Tracer* tracer, HoneyCombDebugger* debugger)
    {
        honeyLetterOpt letterOpts;
        honeycomb::Mesh& inputMesh = ConstructFromTriMesh(trimesh);
        //第1步，寻找底面（最大平面）朝向
        std::vector<int>bottomFaces;
        honeycomb::Point dir = inputMesh.FindBottomDirection(&bottomFaces);
        inputMesh.Rotate(dir, honeycomb::Point(0, 0, -1));
        honeycomb::BoundBox bound = inputMesh.Bound();
        const auto & minPt = bound.Min();
        const auto & maxPt = bound.Max();
        inputMesh.Translate(-minPt);
        letterOpts.bottom.resize(bottomFaces.size());
        letterOpts.bottom.assign(bottomFaces.begin(), bottomFaces.end());
        const auto & honeyFaces = inputMesh.GetFaces();
        std::sort(bottomFaces.begin(), bottomFaces.end());
        std::vector<int> otherFaces(honeyFaces.size() - bottomFaces.size());
        std::set_difference(honeyFaces.begin(), honeyFaces.end(), bottomFaces.begin(), bottomFaces.end(), otherFaces.begin());
        letterOpts.others = std::move(otherFaces);
        //第2步，平移至xoy平面后底面整平
        inputMesh.FlatBottomSurface(&bottomFaces);
        if (debugger) {
            //显示底面区域
            TriPolygons polygons;
            const auto& inPoints = inputMesh.GetPoints();
            const auto& inIndexs = inputMesh.GetFaceVertexAdjacency();
            polygons.reserve(bottomFaces.size());
            for (const auto& f : bottomFaces) {
                TriPolygon poly;
                poly.reserve(3);
                for (int i = 0; i < 3; ++i) {
                    const auto& p = inPoints[inIndexs[f][i]];
                    poly.emplace_back(trimesh::vec3(p.x, p.y, 0));
                }
                polygons.emplace_back(std::move(poly));
            }
            debugger->onGenerateBottomPolygons(polygons);
        }
        //第3步，生成六角网格
        GenerateExHexagons(inputMesh, honeyparams, letterOpts, debugger);
        trimesh::TriMesh& mesh = ConstructFromHoneyMesh(inputMesh);

		mesh.need_adjacentfaces();
		mesh.need_neighbors();
		MMeshT mt(&mesh);
		std::vector<std::vector<trimesh::vec2>> honeycombs;
		for (int i = 0; i < letterOpts.hexgons.size(); i++)
		{
			std::vector<trimesh::vec2> tmp;
			for (int j = 0; j < letterOpts.hexgons[i].borders.size(); j++)
			{
				tmp.push_back(trimesh::vec2(letterOpts.hexgons[i].borders[j].x, letterOpts.hexgons[i].borders[j].y));
			}
			honeycombs.emplace_back(tmp);
		}	
		int old_face_size = mt.faces.size();
		embedingAndCutting(&mt,honeycombs,letterOpts.bottom);
		std::vector<std::vector<std::vector<trimesh::vec2>>> container;
		container.push_back(honeycombs);
		std::vector<int> outfacesIndex;
		int face_size = mt.faces.size();
		std::vector<int> faceindex;
		for (int i = 0; i < letterOpts.bottom.size(); i++)
			faceindex.push_back(letterOpts.bottom[i]);
		for (int i = old_face_size; i < face_size; i++)
			faceindex.push_back(i);
		polygonInnerFaces(&mt, container, faceindex, outfacesIndex);
		//concaveOrConvexOfFaces(&mt, outfacesIndex, true, 10.0f);
		for (int i = 0; i < outfacesIndex.size(); i++)
		{
			mt.faces[outfacesIndex[i]].V0(0)->SetS();
			mt.faces[outfacesIndex[i]].V1(0)->SetS();
			mt.faces[outfacesIndex[i]].V2(0)->SetS();
			mt.faces[outfacesIndex[i]].SetS();
		}		
		for (int vi = 0; vi < mt.vertices.size(); vi++)if(!mt.vertices[vi].IsD())
		{
			if (mt.vertices[vi].IsS())
			{
				for (MMeshFace* f : mt.vertices[vi].connected_face)
				{
					if (!f->IsS())
					{
						mt.vertices[vi].SetV(); break;
					}
				}
			}
		}
		std::vector<std::pair<int, int>> vd;
		findNeighVertex(&mt, letterOpts.others, outfacesIndex, vd);
		for (int i = 0; i < vd.size(); i++)
		{
			float z = mt.vertices[vd[i].second].p.z - 2.0f;
			if (mt.vertices[vd[i].first].IsV())
				splitPoint(&mt, &mt.vertices[vd[i].first], trimesh::point(0, 0, z));				
			else
				mt.vertices[vd[i].first].p += trimesh::point(0, 0, z);
		}		
		trimesh::TriMesh* newMesh = new trimesh::TriMesh();
		//mt.mmesh2trimesh(newMesh);
		mt.quickTransform(newMesh);
		newMesh->write("honeycombs.ply");
		return newMesh;
	}

    TriPolygons GenerateHexagons(const HexagonArrayParam& hexagonparams)
    {
        const float radius = hexagonparams.radius;
        const float nestWidth = hexagonparams.nestWidth;
        const float xdelta = 3.0 / 2.0 * radius;
        const float ydelta = SQRT3 / 2.0 * radius;
        const float side = radius - nestWidth / SQRT3;
        const size_t nrows = hexagonparams.nrows;
        const size_t ncols = hexagonparams.ncols;
        const auto& p0 = hexagonparams.pos;
        std::vector<std::vector<trimesh::point>> gridPoints;
        for (size_t i = 0; i < 2 * nrows - 1; ++i) {
            std::vector<trimesh::point> crossPts;
            for (size_t j = 0; j < ncols; ++j) {
                const auto& pt = trimesh::vec3(p0.x + xdelta * j, p0.y - ydelta * i, p0.z);
                crossPts.emplace_back(pt);
            }
            gridPoints.emplace_back(crossPts);
        }
        const size_t nums = (2 * nrows - 1) * ncols;
        TriPolygons polygons;
        polygons.reserve(nums);
        //计算六角网格对应边界包围盒坐标奇偶性质
        for (size_t i = 0; i < 2 * nrows - 1; i += 2) {
            const auto& rowPoints = gridPoints[i];
            for (size_t j = 0; j < ncols; ++j) {
                const auto& pt = rowPoints[j];
                if (j % 2 == 0) {
                    const auto& center = honeycomb::Point2d(pt.x, pt.y);
                    honeycomb::Hexagon hexagon(center, side);
                    const auto& border = hexagon.Border();
                    TriPolygon polygon;
                    polygon.reserve(border.size());
                    for (const auto& p : border) {
                        polygon.emplace_back(trimesh::vec3((float)p.x, (float)p.y, p0.z));
                    }
                    polygons.emplace_back(std::move(polygon));
                } else {
                    const auto& center = honeycomb::Point2d(pt.x, (double)(pt.y - ydelta));
                    honeycomb::Hexagon hexagon(center, side);
                    const auto & border = hexagon.Border();
                    TriPolygon polygon;
                    polygon.reserve(border.size());
                    for (const auto& p : border) {
                        polygon.emplace_back(trimesh::vec3((float)p.x, (float)p.y, p0.z));
                    }
                    polygons.emplace_back(std::move(polygon));
                }
            }
        }
        return polygons;
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
			//float z = mt.vertices[vd[i].second].p.z-2.0f;
			//if(mt.vertices[vd[i].first].IsV())
			//	splitPoint(&mt, &mt.vertices[vd[i].first],trimesh::point(0,0,z));
			//else
			//	mt.vertices[vd[i].first].p += trimesh::point(0, 0, z);
		}
		//concaveOrConvexOfFaces(&mt, honeycombs, true, 20.f);
		mt.mmesh2trimesh(newMesh);
		newMesh->write("botmesh.ply");	
	}


	void findNeighVertex(MMeshT* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, int>>& vertex_distance)
	{
		std::vector<int> vertexidx;
		std::vector<int> compara_vertex;
		for (int i = 0; i < upfaceid.size(); i++)
		{
			vertexidx.push_back(mesh->faces[upfaceid[i]].V0(0)->index);
			vertexidx.push_back(mesh->faces[upfaceid[i]].V1(0)->index);
			vertexidx.push_back(mesh->faces[upfaceid[i]].V2(0)->index);
		}
		std::sort(vertexidx.begin(), vertexidx.end());
		std::vector<int>::iterator itr = std::unique(vertexidx.begin(), vertexidx.end());
		vertexidx.resize(std::distance(vertexidx.begin(), itr));

		for (int i = 0; i < botfaceid.size(); i++)
		{
			compara_vertex.push_back(mesh->faces[botfaceid[i]].V0(0)->index);
			compara_vertex.push_back(mesh->faces[botfaceid[i]].V1(0)->index);
			compara_vertex.push_back(mesh->faces[botfaceid[i]].V2(0)->index);
		}
		std::sort(compara_vertex.begin(), compara_vertex.end());
		std::vector<int>::iterator it = std::unique(compara_vertex.begin(), compara_vertex.end());
		compara_vertex.resize(std::distance(compara_vertex.begin(), it));

		const int width = 100, height = 100;
		std::vector<std::vector<float>> mapp(width, std::vector<float>(height, std::numeric_limits<float>::max()));
		std::vector<std::vector<int>> mapind(width, std::vector<int>(height, -1));
		
		mesh->getBoundingBox();
		float length_x = (mesh->boundingbox.max_x - mesh->boundingbox.min_x) / (width * 1.0f);
		float length_y = (mesh->boundingbox.max_y - mesh->boundingbox.min_y) / (height * 1.0f);

		float begin_x = mesh->boundingbox.min_x;
		float begin_y = mesh->boundingbox.min_y;

		for (int i = 0; i < vertexidx.size(); i++)
		{
			float x = mesh->vertices[vertexidx[i]].p.x - begin_x;
			float y = mesh->vertices[vertexidx[i]].p.y - begin_y;
			int x_ind = x / length_x;
			int y_ind = y / length_y;
			if (x_ind == width)
				x_ind--;
			if (y_ind == height)
				y_ind--;
			if (mesh->vertices[vertexidx[i]].p.z < mapp[x_ind][y_ind])
			{
				mapp[x_ind][y_ind] = mesh->vertices[vertexidx[i]].p.z;				
				mapind[x_ind][y_ind] = vertexidx[i];
			}
		}

		for (int vi = 0; vi < compara_vertex.size(); vi++)
		{
			int xi = (mesh->vertices[compara_vertex[vi]].p.x - begin_x) / length_x;
			int yi = (mesh->vertices[compara_vertex[vi]].p.y - begin_y) / length_y;
			if (mapp[xi][yi] != std::numeric_limits<float>::max())
			{				
				vertex_distance.push_back(std::make_pair(compara_vertex[vi], mapind[xi][yi]));
			}
			else {
				float max = std::numeric_limits<float>::max();
				int xi_x, yi_y;
				for (int i = 1; i < width; i++)
				{
					for (int ii = xi - i; ii <= xi + i; ii++)
						for (int jj = yi - i; jj <= yi + i; jj++)
						{
							if (std::abs(ii - xi) != i && std::abs(jj - yi) != i) continue;
							if (mapp[ii][jj] < max)
							{
								max = mapp[ii][jj];
								xi_x = ii; yi_y = jj;
							}
						}
					if (max != std::numeric_limits<float>::max())
					{						
						vertex_distance.push_back(std::make_pair(compara_vertex[vi], mapind[xi_x][yi_y]));
						break;
					}
				}
			}
		}
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
			else {
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