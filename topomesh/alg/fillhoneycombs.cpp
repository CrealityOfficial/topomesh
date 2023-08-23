﻿#include "fillhoneycombs.h"
#include "topomesh/math/polygon.h"
#include "topomesh/math/AABB.h"
#include "topomesh/math/SVG.h"
#include "topomesh/data/mmesht.h"
#include "topomesh/alg/letter.h"
#include "topomesh/honeycomb/Matrix.h"
#include "topomesh/honeycomb/Polyline.h"
#include "topomesh/honeycomb/HoneyComb.h"

#include "topomesh/alg/utils.h"
#include "topomesh/data/CMesh.h"
//#include "topomesh/data/entrance.h"
//#include "topomesh/alg/remesh.h"

#ifndef EPS
#define EPS 1e-8f
#endif // !EPS

namespace topomesh {
    enum class HexDirection :int {
        NO_NEIGHBOR = -1,
        X_NEGATIVEZ = 0, ///<(1,0,-1)
        Y_NEGATIVEZ = 1, ///<(0,1,-1)
        Y_NEGATIVEX = 2, ///<(-1,1,0)
        Z_NEGATIVEX = 3, ///<(-1,0,1)
        Z_NEGATIVEY = 4, ///<(0,-1,1)
        X_NEGATIVEY = 5, ///<(1,-1,0)
    };
    HexDirection hex_neighbor(const trimesh::ivec3& va, const trimesh::ivec3& vb)
    {
        const trimesh::ivec3& hex = vb - va;
        const int dist = int(abs(hex.x) + abs(hex.y) + abs(hex.z));
        if (dist == 2) {
            if (hex.x == 1) {
                if (hex.y == 0) return HexDirection::X_NEGATIVEZ;
                else return HexDirection::X_NEGATIVEY;
            } else if (hex.y == 1) {
                if (hex.x == 0) return HexDirection::Y_NEGATIVEZ;
                else return HexDirection::Y_NEGATIVEX;
            } else if (hex.z == 1) {
                if (hex.x == 0) return HexDirection::Z_NEGATIVEY;
                else return HexDirection::Z_NEGATIVEX;
            }
        }
        return HexDirection::NO_NEIGHBOR;
    };
    struct honeyLetterOpt {
        std::vector<int>bottom; ///<底面大块平面的面片索引
        std::vector<int>others; ///<去掉底面后其余面片索引（已保留原模型对应的索引）
        /*
        hexagon: 每个六角网格结构体;
        radius: 包含收缩前的网格边长;
        borders: 收缩后轮廓顶点(不一定6个点，默认xoy坐标系逆时针排序);
        neighbors: 六个方向相邻的六角网格索引(6个方向不一定都有相邻，没有记-1);
        */
        struct hexagon {
            double radius = 1.0;
            std::vector<trimesh::vec3>borders;
            std::vector<int> neighbors;
            hexagon() : neighbors(6, -1) {}
        };
        //所有六角网格及邻居关系
        std::vector<hexagon>hexgons;
    };

    honeycomb::Mesh ConstructFromTriMesh(const trimesh::TriMesh* trimesh)
    {
        honeycomb::Mesh mesh;
        auto& faces = mesh.GetFaces();
        auto& points = mesh.GetPoints();
        auto& indexs = mesh.GetFaceVertexAdjacency();
        const auto& vertices = trimesh->vertices;
        const auto& faceVertexs = trimesh->faces;
        const int facenums = faceVertexs.size();
        faces.reserve(facenums);
        for (int i = 0; i < facenums; ++i) {
            faces.emplace_back(i);
        }
        points.reserve(vertices.size());
        for (const auto& v : vertices) {
            points.emplace_back(honeycomb::Point{ v[0], v[1], v[2] });
        }
        indexs.reserve(faceVertexs.size());
        for (const auto& fa : faceVertexs) {
            indexs.emplace_back(std::vector<int>{fa[0], fa[1], fa[2] });
        }
        mesh.GenerateFaceNormals();
        //mesh.WriteSTLFile("trimesh2mesh");
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
            triMesh.faces.emplace_back(trimesh::vec3(f[0], f[1], f[2]));
        }
        const auto& min = bound.Min();
        const auto& max = bound.Max();
        const auto& p1 = trimesh::vec3((float)min.x, (float)min.y, (float)min.z);
        const auto& p2 = trimesh::vec3((float)max.x, (float)max.y, (float)max.z);
        triMesh.bbox = std::move(trimesh::box3(p1, p2));
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
                const auto& pts = poly.Points();
                TriPolygon polygon;
                polygon.reserve(pts.size());
                for (const auto& p : pts) {
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
        cxutil::Polygons mpolygons; ///<底面边界抽壳多边形
        {
            for (const auto& poly : polys) {
                const auto& pts = poly.Points();
                ClipperLib::Path path;
                path.reserve(pts.size());
                for (const auto& p : pts) {
                    const auto& x = int(p.x / resolution);
                    const auto& y = int(p.y / resolution);
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
        auto min = boundarys.Bound().Min();
        auto max = boundarys.Bound().Max();
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
        trimesh::vec3 origin(min.x, max.y, 0);

        HexagonArrayParam hexagonparams;
        hexagonparams.pos = origin;
        hexagonparams.radius = radius;
        hexagonparams.nestWidth = honeyparams.nestWidth;
        hexagonparams.ncols = std::ceil((max.x - min.x) / xdist);
        hexagonparams.nrows = (std::ceil((max.y - min.y) / ydist) + 1) / 2;
        HexaPolygons hexagons = GenerateHexagons(hexagonparams);
        if (debugger) {
            //显示六角网格阵列
            debugger->onGenerateBottomPolygons(convertFromHexaPolygons(hexagons));
        }
        //第8步，计算网格阵列与抽壳轮廓交集
        std::vector<cxutil::Polygons> hpolygons; ///<完整六角网格多边形序列
        std::vector<cxutil::Polygons> ipolygons; ///<初次求交网格多边形序列
        std::vector<trimesh::ivec3> icoords;
        {
            hpolygons.reserve(hexagons.size());
            for (const auto& hexa : hexagons) {
                cxutil::Polygons polygons;
                ClipperLib::Path path;
                path.reserve(hexa.poly.size());
                for (const auto& p : hexa.poly) {
                    const auto& x = int(p.x / resolution);
                    const auto& y = int(p.y / resolution);
                    path.emplace_back(ClipperLib::IntPoint(x, y));
                }
                polygons.paths.emplace_back(std::move(path));
                const cxutil::Polygons & ipolys = mpolygons.intersection(polygons);
                if (!ipolys.empty()) {
                    icoords.emplace_back(hexa.coord);
                    ipolygons.emplace_back(std::move(ipolys));
                }
                hpolygons.emplace_back(std::move(polygons));
            }
        }
        if (debugger) {
            //显示初始求交结果
            TriPolygons tripolys;
            tripolys.reserve(ipolygons.size());
            for (const auto& ipolys : ipolygons) {
                TriPolygon poly;
                ClipperLib::Path path = ipolys.paths.front();
                poly.reserve(path.size());
                for (const auto& p : path) {
                    poly.emplace_back(trimesh::vec3(p.X * resolution, p.Y * resolution, 0));
                }
                tripolys.emplace_back(std::move(poly));
            }
            debugger->onGenerateBottomPolygons(tripolys);
        }
        //第9步，筛选符合需求的网格多边形
        TriPolygons outtripolys; ///<最终保留求交网格多边形序列
        const double minrate = honeyparams.keepHexagonRate;
        const double hexarea = 3.0 / 2.0 * SQRT3 * std::pow(side / resolution, 2);
        std::vector<trimesh::ivec3> ocoords;
        const int isize = ipolygons.size();
        ocoords.reserve(isize);
        for (int i = 0; i < isize; i++) {
            TriPolygon tripolys;
            cxutil::Polygons ipolys = ipolygons[i];
            ClipperLib::Path& path = ipolys.paths.front();
            if (ipolys.area() < 0) {
                ClipperLib::Path tmp;
                tmp.swap(path);
                std::reverse(tmp.begin(), tmp.end());
                ipolys.paths.emplace_back(std::move(tmp));
            }
            if (ipolys.area() >= minrate * hexarea) {
                for (const auto& p : path) {
                    const auto& point = trimesh::vec3(p.X * resolution, p.Y * resolution, 0);
                    tripolys.emplace_back(std::move(point));
                }
                ocoords.emplace_back(std::move(icoords[i]));
                outtripolys.emplace_back(std::move(tripolys));
            }
        }
        if (debugger) {
            //显示最终符合需求的网格多边形
            debugger->onGenerateBottomPolygons(outtripolys);
        }
        const int osize = outtripolys.size();
        letterOpts.hexgons.reserve(osize);
        for (const auto& poly : outtripolys) {
            honeyLetterOpt::hexagon hexa;
            hexa.borders = poly;
            hexa.radius = radius;
            letterOpts.hexgons.emplace_back(std::move(hexa));
        }
        std::vector<std::vector<int>> hoods(osize, std::vector<int>(osize, 1));
        for (int i = 0; i < osize; ++i) {
            for (int j = 0; j < osize && (j != i) && hoods[i][j]; ++j) {
                const int res = int(hex_neighbor(ocoords[i], ocoords[j]));
                if (res >= 0) {
                    letterOpts.hexgons[i].neighbors[res] = j;
                    letterOpts.hexgons[j].neighbors[(res + 3) % 6] = i;
                    hoods[i][j] = 0;
                    hoods[i][j] = 0;
                }
            }
        }
        return;
    }
    trimesh::TriMesh* generateHoneyCombs(trimesh::TriMesh* trimesh, const HoneyCombParam& honeyparams,
        ccglobal::Tracer* tracer, HoneyCombDebugger* debugger)
    {	        
        honeyLetterOpt letterOpts;
        honeycomb::Mesh&& inputMesh = ConstructFromTriMesh(trimesh);
        //inputMesh.WriteSTLFile("inputmesh");
        //第1步，寻找底面（最大平面）朝向
        std::vector<int>bottomFaces;
        honeycomb::Point dir = inputMesh.FindBottomDirection(&bottomFaces);
        inputMesh.Rotate(dir, honeycomb::Point(0, 0, -1));
        honeycomb::BoundBox bound = inputMesh.Bound();
        const auto& minPt = bound.Min();
        inputMesh.Translate(-minPt);
        letterOpts.bottom.resize(bottomFaces.size());
        letterOpts.bottom.assign(bottomFaces.begin(), bottomFaces.end());
        const auto& honeyFaces = inputMesh.GetFaces();
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
        trimesh::TriMesh&& mesh = ConstructFromHoneyMesh(inputMesh);
		
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
		std::vector<std::pair<int, float>> vd;
		findNeighVertex(&mt, letterOpts.others, outfacesIndex, vd);
		for (int i = 0; i < vd.size(); i++)
		{
			float z = vd[i].second - honeyparams.shellThickness;
			if (z < mt.vertices[vd[i].first].p.z) continue;;
			if (mt.vertices[vd[i].first].IsV())
				splitPoint(&mt, &mt.vertices[vd[i].first], trimesh::point(0, 0, z));	
				//mt.appendVertex(trimesh::point(mt.vertices[vd[i].first].p + trimesh::point(0, 0, z)));
			else
				mt.vertices[vd[i].first].p += trimesh::point(0, 0, z);
		}		
		trimesh::TriMesh* newMesh = new trimesh::TriMesh();
		//mt.mmesh2trimesh(newMesh);
		mt.quickTransform(newMesh);
        //restore original pose
        const auto& dir1 = trimesh::vec3(0, 0, -1);
        const auto& dir2 = trimesh::vec3(dir.x, dir.y, dir.z);
        trimesh::fxform && mat = trimesh::fxform::rot_into(dir1, dir2);
        for (auto& p : newMesh->vertices) {
            p += trimesh::vec3(minPt.x, minPt.y, minPt.z);
            p = mat * p;
        }
		newMesh->write("honeycombs.ply");
		return newMesh;
	}

    HexaPolygons GenerateHexagons(const HexagonArrayParam& hexagonparams)
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
        HexaPolygons polygons;
        polygons.reserve(nums);
        //计算六角网格对应边界包围盒坐标奇偶性质
        for (int i = 0; i < 2 * nrows - 1; i += 2) {
            const auto& rowPoints = gridPoints[i];
            for (int j = 0; j < ncols; ++j) {
                const auto& pt = rowPoints[j];
                if (j % 2 == 0) {
                    const auto& center = honeycomb::Point2d(pt.x, pt.y);
                    honeycomb::Hexagon hexagon(center, side);
                    const auto& border = hexagon.Border();
                    HexaPolygon hexa;
                    hexa.side = side;
                    hexa.poly.reserve(border.size());
                    for (const auto& p : border) {
                        hexa.poly.emplace_back(trimesh::vec3((float)p.x, (float)p.y, p0.z));
                    }
                    hexa.coord = trimesh::ivec3(j, -(i + j) / 2, (i - j) / 2);
                    polygons.emplace_back(std::move(hexa));
                } else {
                    const auto& center = honeycomb::Point2d(pt.x, (double)pt.y - ydelta);
                    honeycomb::Hexagon hexagon(center, side);
                    const auto & border = hexagon.Border();
                    HexaPolygon hexa;
                    hexa.side = side;
                    hexa.poly.reserve(border.size());
                    for (const auto& p : border) {
                        hexa.poly.emplace_back(trimesh::vec3((float)p.x, (float)p.y, p0.z));
                    }
                    hexa.coord = trimesh::ivec3(j, -(i + 1 + j) / 2, (i + 1 - j) / 2);
                    polygons.emplace_back(std::move(hexa));
                }
            }
        }
        return polygons;
    }
    void GenerateHexagonNeighbors(HexaPolygons& hexas, float cheight)
    {
        const int size = hexas.size();
        std::vector<std::vector<int>> hoods(size, std::vector<int>(size, 1));
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size && (j != i); ++j) {
                if (hoods[i][j]) {
                    const int res = int(hex_neighbor(hexas[i].coord, hexas[j].coord));
                    if (res >= 0) {
                        hexas[i].neighbors[res] = j;
                        hexas[j].neighbors[(res + 3) % 6] = i;
                        hoods[i][j] = 0;
                        hoods[i][j] = 0;
                    }
                }
            }
        }
        for (auto& hexa : hexas) {
            for (int j = 0; j < 6; j++) {
                if ((!hexa.canAdds[j]) && (hexa.neighbors[j] >= 0)) {
                    auto& hn = hexas[hexa.neighbors[j]];
                    hexa.ctop = hexa.side * hexa.ratio * 0.5f + cheight;
                    hn.ctop = hn.side * hn.ratio * 0.5f + cheight;
                    if ((hexa.depth > hexa.ctop) && (hn.depth > hn.ctop)) {
                        hexa.canAdds[j] = true;
                        hn.canAdds[(j + 3) % 6] = true;
                    }
                }
            }
        }
        return;
    }

    TriPolygons traitCurrentPolygons(const HexaPolygons& hexas, int index)
    {
        topomesh::TriPolygons polys;
        if (index >= 0 && index < (hexas.size())) {
            polys.emplace_back(hexas[index].poly);
        }
        return polys;
    }

    TriPolygons traitNeighborPolygons(const HexaPolygons& hexas, int index)
    {
        topomesh::TriPolygons polys;
        if (index >= 0 && index < (hexas.size()))             {
            const auto& neighbors = hexas.at(index).neighbors;
            int nums = neighbors.size();
            polys.reserve(nums);
            for (int i = 0; i < nums; ++i) {
                if (neighbors[i] >= 0)
                    polys.emplace_back(hexas.at(neighbors[i]).poly);
            }
        }
        return polys;
    }

    TriPolygons traitDirctionPolygon(const HexaPolygons& hexas, int index, int dir)
    {
        topomesh::TriPolygons polys;
        if (index >= 0 && index < (hexas.size())) {
            if (dir >= 0 && dir <= 5) {
                const auto& neighbors = hexas.at(index).neighbors;
                int val = neighbors[dir];
                if (val >= 0) {
                    polys.emplace_back(hexas[val].poly);
                }
            }
        }
        return polys;
    }

    TriPolygon traitPlanarCircle(const trimesh::vec3& center, float r, std::vector<int>& indexs, const trimesh::vec3& dir, int nums)
    {
        TriPolygon points;
        points.reserve(nums);
        trimesh::vec3 ydir = trimesh::normalized(dir);
        trimesh::vec3 zdir(0, 0, 1);
        const float phi = 2.0f * M_PI / float(nums);
        for (int i = 0; i < nums; ++i) {
            const auto& theta = phi * (float)i + (float)M_PI_2;
            const auto& y = ydir * r * std::cos(theta);
            const auto& z = zdir * r * std::sin(theta);
            points.emplace_back(std::move(center + y + z));
        }
        float miny = points.front() DOT ydir;
        float minz = points.front() DOT zdir;
        float maxy = points.front() DOT ydir;
        float maxz = points.front() DOT zdir;
        for (int i = 1; i < nums; ++i) {
            float yn = points[i] DOT ydir;
            float zn = points[i] DOT zdir;
            if (yn < miny) miny = yn;
            if (zn < minz) minz = zn;
            if (yn > maxy) maxy = yn;
            if (zn > maxz) maxz = zn;
        }
        if (nums % 2 == 0) indexs.reserve(4);
        else indexs.reserve(5);
        for (int i = 0; i < nums; ++i) {
            float zn = points[i] DOT zdir;
            if (std::abs(zn - maxz) < EPS) {
                indexs.emplace_back(i);
            }
        }
        for (int i = 0; i < nums; ++i) {
            float yn = points[i] DOT ydir;
            if (std::abs(yn - miny) < EPS) {
                indexs.emplace_back(i);
            }
        }
        for (int i = 0; i < nums; ++i) {
            float zn = points[i] DOT zdir;
            if (std::abs(zn - minz) < EPS) {
                indexs.emplace_back(i);
            }
        }
        for (int i = 0; i < nums; ++i) {
            float yn = points[i] DOT ydir;
            if (std::abs(yn - maxy) < EPS) {
                indexs.emplace_back(i);
            }
        }
        return points;
    }

    std::shared_ptr<trimesh::TriMesh> generateHolesColumnar(HexaPolygons& hexas, const ColumnarHoleParam& param)
    {
        const float cheight = param.cheight;
        GenerateHexagonNeighbors(hexas, cheight);
        std::vector<trimesh::vec3> points;
        std::vector<trimesh::ivec3> faces;
        int hexaspointnum = 0, num2 = 0;
        for (auto& hexa : hexas) {
            hexa.startIndex = num2;
            for (const auto& p : hexa.poly) {
                points.emplace_back(std::move(p));
                ++num2, ++hexaspointnum;
            }
        }

        for (const auto& hexa : hexas) {
            for (const auto& p : hexa.poly) {
                const auto& pt = trimesh::vec3(p.x, p.y, p.z + hexa.depth);
                points.emplace_back(std::move(pt));
                ++num2;
            }
        }
        int holesnum = 0, nslices = param.nslices;
        for (auto& hexa : hexas) {
            for (int i = 0; i < 6; ++i) {
                if (hexa.canAdds[i] && (!hexa.hasAdds[i]) && (hexa.neighbors[i] >= 0)) {
                    auto& hn = hexas[hexa.neighbors[i]];
                    const auto& a = hexa.poly[i];
                    const auto& b = hexa.poly[(i + 1) % 6];
                    const auto& r1 = hexa.side * hexa.ratio * 0.5f;
                    const auto& c1 = (a + b) / 2.0 + trimesh::vec3(0, 0, cheight);
                    std::vector<int> corner1;
                    const auto& poly1 = traitPlanarCircle(c1, r1, corner1, a - b, nslices);
                    hexa.corners[i].swap(corner1);
                    const auto& c = hn.poly[(i + 3) % 6];
                    const auto& d = hn.poly[(i + 4) % 6];
                    const auto& r2 = hn.side * hn.ratio * 0.5f;
                    const auto& c2 = (c + d) / 2.0 + trimesh::vec3(0, 0, cheight);
                    std::vector<int> corner2;
                    const auto& poly2 = traitPlanarCircle(c2, r2, corner2, c - d, nslices);
                    hn.corners[(i + 3) % 6].swap(corner2);
                    points.insert(points.end(), poly1.begin(), poly1.end());
                    points.insert(points.end(), poly2.begin(), poly2.end());
                    hexa.hasAdds[i] = true;
                    hn.hasAdds[(i + 3) % 6] = true;
                    hexa.holeIndexs[i] = holesnum;
                    hn.holeIndexs[(i + 3) % 6] = holesnum + 1;
                    holesnum += 2;
                }
            }
        }
        int holefacenums = holesnum * nslices;
        int rectfacenums = (hexaspointnum - holesnum) * 2;
        int holerectfacenums = holefacenums * (nslices + 4);
        int upperfacenums = hexaspointnum;
        int bottomfacenums = 0;
        int allfacenums = holefacenums + rectfacenums + holerectfacenums + upperfacenums + bottomfacenums;
        faces.reserve(allfacenums);
        for (int i = 0; i < hexas.size(); ++i) {
            const auto& hexa = hexas[i];
            const auto& poly = hexa.poly;
            int polynums = poly.size();
            const int& start = hexa.startIndex;
            //六角网格顶部
            const int& upstart = start + hexaspointnum;
            for (int j = 1; j < polynums - 1; ++j) {
                const int& cur = upstart + j;
                const int& next = upstart + j + 1;
                faces.emplace_back(trimesh::ivec3(upstart, next, cur));
            }
            //六角网格侧面
            for (int j = 0; j < polynums; ++j) {
                const int& a = start + j;
                const int& b = start + (j + 1) % polynums;
                const int& c = a + hexaspointnum;
                const int& d = b + hexaspointnum;
                if (hexa.hasAdds[j]) {
                    const auto& corner = hexa.corners[j];
                    const int& cstart = num2 + hexa.holeIndexs[j] * nslices;
                    if (nslices % 2 == 0) {
                        int r = corner[0];
                        int s = corner[1];
                        int t = corner[2];
                        int u = corner[3];
                        faces.emplace_back(trimesh::ivec3(cstart + r, c, d));
                        for (int k = 0; k < s; ++k) {
                            faces.emplace_back(trimesh::ivec3(cstart + k, d, cstart + k + 1));
                        }
                        faces.emplace_back(trimesh::ivec3(b, cstart + s, d));
                        for (int k = s; k < t; ++k) {
                            faces.emplace_back(trimesh::ivec3(cstart + k, b, cstart + k + 1));
                        }
                        faces.emplace_back(trimesh::ivec3(b, a, cstart + t));
                        for (int k = t; k < u; ++k) {
                            faces.emplace_back(trimesh::ivec3(a, cstart + k + 1, cstart + k));
                        }
                        faces.emplace_back(trimesh::ivec3(a, c, cstart + u));
                        for (int k = u; k < nslices; ++k) {
                            faces.emplace_back(trimesh::ivec3(c, cstart + (k + 1) % nslices, cstart + k));
                        }
                    } else {
                        int r = corner[0];
                        int s = corner[1];
                        int t = corner[2];
                        int u = corner[4];
                        faces.emplace_back(trimesh::ivec3(cstart + r, c, d));
                        for (int k = 0; k < s; ++k) {
                            faces.emplace_back(trimesh::ivec3(cstart + k, d, cstart + k + 1));
                        }
                        faces.emplace_back(trimesh::ivec3(b, cstart + s, d));
                        for (int k = s; k < t; ++k) {
                            faces.emplace_back(trimesh::ivec3(cstart + k, b, cstart + k + 1));
                        }
                        faces.emplace_back(trimesh::ivec3(b, a, cstart + t));
                        for (int k = t; k < u; ++k) {
                            faces.emplace_back(trimesh::ivec3(a, cstart + k + 1, cstart + k));
                        }
                        faces.emplace_back(trimesh::ivec3(a, c, cstart + u));
                        for (int k = u; k < nslices; ++k) {
                            faces.emplace_back(trimesh::ivec3(c, cstart + (k + 1) % nslices, cstart + k));
                        }
                    }
                } else {
                    faces.emplace_back(trimesh::ivec3(b, a, c));
                    faces.emplace_back(trimesh::ivec3(b, c, d));
                }
            }
        }
        //六角网格孔部
        for (int i = 0; i < holesnum; i += 2) {
            const int& start = num2 + i * nslices;
            for (int j = 0; j < nslices; ++j) {
                const int& lcur = start + j;
                const int& lnext = start + (j + 1) % nslices;
                const int& rcur = start + (nslices - j) % nslices + nslices;
                const int& rnext = start + (nslices - j - 1) % nslices + nslices;
                faces.emplace_back(trimesh::ivec3(lnext, rnext, rcur));
                faces.emplace_back(trimesh::ivec3(lnext, rcur, lcur));
            }
        }
        //六角网格底部

        std::shared_ptr<trimesh::TriMesh> triMesh(new trimesh::TriMesh());
        triMesh->vertices.swap(points);
        triMesh->faces.swap(faces);
        //triMesh->write("holes.stl");
        return triMesh;
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


	void findNeighVertex(MMeshT* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, float>>& vertex_distance)
	{		
		
		std::vector<int> compara_vertex;		
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
		//std::vector<std::vector<float>> mapp(width, std::vector<float>(height, std::numeric_limits<float>::max()));
		std::vector<std::vector<std::vector<int>>> mapind(width, std::vector<std::vector<int>>(height));
		
		mesh->getBoundingBox();
		float length_x = (mesh->boundingbox.max_x - mesh->boundingbox.min_x) / (width * 1.0f);
		float length_y = (mesh->boundingbox.max_y - mesh->boundingbox.min_y) / (height * 1.0f);

		float begin_x = mesh->boundingbox.min_x;
		float begin_y = mesh->boundingbox.min_y;


		for (int i = 0; i < upfaceid.size(); i++)
		{
			float max_x=std::numeric_limits<float>::min(), max_y = std::numeric_limits<float>::min();
			float min_x = std::numeric_limits<float>::max(), min_y = std::numeric_limits<float>::max();
			for (int j = 0; j < 3; j++)
			{
				if (mesh->faces[upfaceid[i]].V0(j)->p.x > max_x)
					max_x = mesh->faces[upfaceid[i]].V0(j)->p.x;
				if (mesh->faces[upfaceid[i]].V0(j)->p.x < min_x)
					min_x = mesh->faces[upfaceid[i]].V0(j)->p.x;
				if (mesh->faces[upfaceid[i]].V0(j)->p.y > max_y)
					max_y = mesh->faces[upfaceid[i]].V0(j)->p.y;
				if (mesh->faces[upfaceid[i]].V0(j)->p.y < min_y)
					min_y = mesh->faces[upfaceid[i]].V0(j)->p.y;
			}
			int x_min_id =(min_x - begin_x) / length_x;
			int y_min_id = (min_y - begin_y) / length_y;
			int x_max_id = (max_x - begin_x) / length_x;
			int y_max_id = (max_y - begin_y) / length_y;

			if (x_max_id == width)
				x_max_id--;
			if (y_max_id == height)
				y_max_id--;

			for (int y = y_min_id; y <= y_max_id; y++)
				for (int x = x_min_id; x <= x_max_id; x++)
					mapind[x][y].push_back(upfaceid[i]);
		}

		for (int vi = 0; vi < compara_vertex.size(); vi++)
		{
			int xi = (mesh->vertices[compara_vertex[vi]].p.x - begin_x) / length_x;
			int yi = (mesh->vertices[compara_vertex[vi]].p.y - begin_y) / length_y;
			if (xi == width)
				xi--;
			if (yi == height)
				yi--;
			float min_z = std::numeric_limits<float>::max();
			for (int ii = 0; ii < mapind[xi][yi].size(); ii++)
			{
				for (int k = 0; k < 3; k++)
				{
					if (mesh->faces[mapind[xi][yi][ii]].V0(k)->p.z < min_z)
						min_z = mesh->faces[mapind[xi][yi][ii]].V0(k)->p.z;
				}
			}
			vertex_distance.push_back(std::make_pair(compara_vertex[vi], min_z));
		}


		/*for (int i = 0; i < vertexidx.size(); i++)
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
		}*/
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
				//std::cout << "mapind :" << mapind[xi][yi] << " z : " << mapp[xi][yi] << "\n";
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
						//std::cout << "mapind :" << mapind[xi_x][yi_y] << " z : " << max << "\n";
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


	HoneyCombContext::HoneyCombContext(std::shared_ptr<trimesh::TriMesh> mesh)
	{
		m_mesh = mesh;
		mesh->need_adjacentfaces();
		mesh->need_neighbors();
		mesh->need_normals();
		innerMesh.reset(new MMeshT(mesh.get()));
	}

	HoneyCombContext::~HoneyCombContext()
	{

	}

	void HoneyCombContext::checkNeigbour(int indicate, std::vector<int>& faceIndexs, float angle_threshold)
	{
		findNeighborFacesOfSameAsNormal(innerMesh.get(), indicate, faceIndexs, angle_threshold);
	}

	trimesh::TriMesh* HoneyCombContext::data()
	{
		return m_mesh.get();
	}
}