#include "fillhoneycombs.h"
#include "topomesh/data/mmesht.h"
#include "topomesh/alg/letter.h"
#include "topomesh/alg/earclipping.h"

#include "topomesh/alg/utils.h"
#include "topomesh/data/CMesh.h"
#include "trimesh2/TriMesh_algo.h"
#include "topomesh/alg/solidtriangle.h"
#include "clipper/clipper.hpp"
#include "cxutil/math/polygon.h"

//#include "topomesh/data/entrance.h"
//#include "topomesh/alg/remesh.h"

#ifndef EPS
#define EPS 1e-8f
#endif // !EPS
#ifndef SQRT3
#define SQRT3 1.732050807568877
#endif // !SQRT3

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

    void GenerateExHexagons(CMesh& honeyMesh, const HoneyCombParam& honeyparams, honeyLetterOpt& letterOpts, HoneyCombDebugger* debugger = nullptr)
    {
        //拷贝一份数据
        CMesh cutMesh;
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
        TriPolygons polys;
        polys.reserve(sequentials.size());
        const auto& points = cutMesh.mpoints;
        for (const auto& seq : sequentials) {
            TriPolygon border;
            border.reserve(seq.size());
            for (const auto& v : seq) {
                const auto& p = points[v];
                border.emplace_back(p.x, p.y, 0);
            }
            polys.emplace_back(border);
        }
        if (debugger) {
            //显示底面边界轮廓多边形
            debugger->onGenerateInfillPolygons(polys);
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
                ClipperLib::Path path;
                path.reserve(poly.size());
                for (const auto& p : poly) {
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
        trimesh::dvec3 min, max;
        min.x = std::numeric_limits<double>::max();
        min.y = std::numeric_limits<double>::max();
        max.x = std::numeric_limits<double>::lowest();
        max.y = std::numeric_limits<double>::lowest();
        for (const auto& poly : polys) {
            for (const auto& p : poly) {
                if (p.x < min.x) min.x = p.x;
                if (p.y < min.y) min.y = p.y;
                if (p.x > max.x) max.x = p.x;
                if (p.y > max.y) max.y = p.y;
            }
        }
        if (honeyparams.polyline) {
            min.x = std::numeric_limits<double>::max();
            min.y = std::numeric_limits<double>::max();
            max.x = std::numeric_limits<double>::lowest();
            max.y = std::numeric_limits<double>::lowest();
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
            hpolygons.reserve(hexagons.polys.size());
            for (const auto& hexa : hexagons.polys) {
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

    trimesh::vec3 adjustHoneyCombParam(trimesh::TriMesh* trimesh,const HoneyCombParam& honeyparams)
    {
       // if (trimesh == nullptr) return ;
        if (!honeyparams.faces.empty())
        {
            trimesh::vec3 ave_normal(0, 0, 0);
            for (int fi : honeyparams.faces)
            {
                trimesh::vec3 v1 = trimesh->vertices[trimesh->faces[fi][1]] - trimesh->vertices[trimesh->faces[fi][0]];
                trimesh::vec3 v2 = trimesh->vertices[trimesh->faces[fi][2]] - trimesh->vertices[trimesh->faces[fi][0]];
                ave_normal += trimesh::normalized(v1 % v2);
            }
            ave_normal /= (ave_normal / (honeyparams.faces.size() * 1.f));
            trimesh::normalize(ave_normal);
           /* trimesh::vec3* dir = const_cast<trimesh::vec3*>(&honeyparams.axisDir);
            *dir = ave_normal;*/
            return ave_normal;
        }else
             return honeyparams.axisDir;
    }

    trimesh::TriMesh* GenerateHoneyCombs(trimesh::TriMesh* trimesh , const HoneyCombParam& honeyparams ,
        ccglobal::Tracer* tracer , HoneyCombDebugger* debugger )
    {
        //0.初始化Cmesh,并旋转
        trimesh::TriMesh* result_trimesh = new trimesh::TriMesh();
        if(trimesh == nullptr) return result_trimesh;
        CMesh cmesh(trimesh);
        //1.重新整理输入参数
        trimesh::vec3 dir=adjustHoneyCombParam(trimesh, honeyparams);
        cmesh.Rotate(dir, trimesh::vec3(0, 0, 1));
        //2.找到生成蜂窝指定的区域（自定义或者是用户自己指定）          
        cmesh.GenerateBoundBox();
        HexagonArrayParam hexagonparams;
        hexagonparams.dir = honeyparams.axisDir;
        hexagonparams.pos = trimesh::vec3(cmesh.mbox.min.x, cmesh.mbox.max.y, 0);
        hexagonparams.radius = honeyparams.honeyCombRadius;
        hexagonparams.nestWidth = honeyparams.nestWidth;
        const auto& xdist = 3.0 / 2.0 * hexagonparams.radius;
        const auto& ydist = SQRT3 / 2.0 * hexagonparams.radius;
        hexagonparams.ncols = std::ceil((cmesh.mbox.max.x - cmesh.mbox.min.x) / xdist);
        hexagonparams.nrows = (std::ceil((cmesh.mbox.max.y - cmesh.mbox.min.y) / ydist) + 1) / 2;
        HexaPolygons hexagons = GenerateHexagons(hexagonparams);

        
        return result_trimesh;
    }

    trimesh::TriMesh* generateHoneyCombs(trimesh::TriMesh* trimesh, const HoneyCombParam& honeyparams,
        ccglobal::Tracer* tracer, HoneyCombDebugger* debugger)
    {	        
        honeyLetterOpt letterOpts;
        CMesh inputMesh(trimesh);
        //inputMesh.WriteSTLFile("inputmesh");
        //第1步，寻找底面（最大平面）朝向
        std::vector<int>bottomFaces;
        trimesh::vec3 dir = inputMesh.FindBottomDirection(&bottomFaces);
        inputMesh.Rotate(dir, trimesh::vec3(0, 0, -1));
        const auto& minPt = inputMesh.mbox.min;
        inputMesh.Translate(-minPt);
        letterOpts.bottom.resize(bottomFaces.size());
        letterOpts.bottom.assign(bottomFaces.begin(), bottomFaces.end());
        std::vector<int> honeyFaces;
        honeyFaces.reserve(inputMesh.mfaces.size());
        for (int i = 0; i < inputMesh.mfaces.size(); ++i) {
            honeyFaces.emplace_back(i);
        }
        std::sort(bottomFaces.begin(), bottomFaces.end());
        std::vector<int> otherFaces(honeyFaces.size() - bottomFaces.size());
        std::set_difference(honeyFaces.begin(), honeyFaces.end(), bottomFaces.begin(), bottomFaces.end(), otherFaces.begin());
        letterOpts.others = std::move(otherFaces);
        //第2步，平移至xoy平面后底面整平
        inputMesh.FlatBottomSurface(&bottomFaces);
        if (debugger) {
            //显示底面区域
            TriPolygons polygons;
            const auto& inPoints = inputMesh.mpoints;
            const auto& inIndexs = inputMesh.mfaces;
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
        trimesh::TriMesh&& mesh = inputMesh.GetTriMesh();

        std::vector<bool> delectface(mesh.faces.size(), false);
        for (int i : bottomFaces)
            delectface[i] = true;
        trimesh::remove_faces(&mesh, delectface); 

        mesh.need_bbox();
        std::vector<std::tuple<trimesh::point, trimesh::point, trimesh::point>> data;
        for (trimesh::TriMesh::Face& f : mesh.faces)
            data.push_back(std::make_tuple(mesh.vertices[f[0]], mesh.vertices[f[1]], mesh.vertices[f[2]]));
        topomesh::SolidTriangle sd(&data, 100, 100, mesh.bbox.max.x, mesh.bbox.min.x, mesh.bbox.max.y, mesh.bbox.min.y);
        sd.work();
        std::vector<std::vector<float>> result = sd.getResult();
        trimesh::TriMesh* newmesh = new trimesh::TriMesh();
        

        for(int i=0;i<result.size();i++)
            for (int j = 0; j < result[i].size(); j++)
            {
                if (result[i][j] == std::numeric_limits<float>::max()) { result[i][j] = 0.f;/* std::cout << "i :" << i << " j : " << j << "\n";*/ }
                newmesh->vertices.push_back(trimesh::point(i, j, result[i][j]));
            }
    
       /* for (int i = 0; i < letterOpts.hexgons.size(); i++)
        {        
            for (int j = 0; j < letterOpts.hexgons[i].borders.size(); j++)
            {
                float h = sd.getCoordData(letterOpts.hexgons[i].borders[j].x, letterOpts.hexgons[i].borders[j].y);
                if (h == std::numeric_limits<float>::max()) h = 0.f;
                newmesh->vertices.push_back(trimesh::point(letterOpts.hexgons[i].borders[j].x, letterOpts.hexgons[i].borders[j].y, h));
            }          
        }*/
        newmesh->write("newhoneymesh.ply");   
        mesh.write("honeyhole.ply");
        return &mesh;
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
			if (z < mt.vertices[vd[i].first].p.z) continue;
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
        polygons.polys.reserve(nums);
        polygons.side = side;
        //计算六角网格对应边界包围盒坐标奇偶性质
        for (int i = 0; i < 2 * nrows - 1; i += 2) {
            const auto& rowPoints = gridPoints[i];
            for (int j = 0; j < ncols; ++j) {
                const auto& pt = rowPoints[j];
                if (j % 2 == 0) {
                    const auto& center = trimesh::vec2(pt.x, pt.y);
                    topomesh::Hexagon hexagon(center, side);
                    const auto& border = hexagon.border;
                    HexaPolygon hexa;
                    hexa.poly.reserve(border.size());
                    for (const auto& p : border) {
                        hexa.poly.emplace_back(trimesh::vec3((float)p.x, (float)p.y, p0.z));
                    }
                    hexa.coord = trimesh::ivec3(j, -(i + j) / 2, (i - j) / 2);
                    polygons.polys.emplace_back(std::move(hexa));
                } else {
                    const auto& center = trimesh::vec2(pt.x, (double)pt.y - ydelta);
                    topomesh::Hexagon hexagon(center, side);
                    const auto & border = hexagon.border;
                    HexaPolygon hexa;
                    hexa.poly.reserve(border.size());
                    for (const auto& p : border) {
                        hexa.poly.emplace_back(trimesh::vec3((float)p.x, (float)p.y, p0.z));
                    }
                    hexa.coord = trimesh::ivec3(j, -(i + 1 + j) / 2, (i + 1 - j) / 2);
                    polygons.polys.emplace_back(std::move(hexa));
                }
            }
        }
        return polygons;
    }
    
    void GenerateHexagonNeighbors(HexaPolygons& hexas, const ColumnarHoleParam& param)
    {
        const int size = hexas.polys.size();
        std::vector<std::vector<int>> hoods(size, std::vector<int>(size, 1));
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (hoods[i][j] && (j != i)) {
                    const int res = int(hex_neighbor(hexas.polys[i].coord, hexas.polys[j].coord));
                    if (res >= 0) {
                        hexas.polys[i].edges[res].neighbor = j;
                        hexas.polys[j].edges[(res + 3) % 6].neighbor = i;
                        hoods[i][j] = 0;
                        hoods[i][j] = 0;
                    }
                }
            }
        }
        const float cradius = hexas.side * param.ratio * 0.5f;
        const float cdelta = param.delta + 2 * cradius;
        const float cheight = param.height + cradius;
        float cmax = param.height + 2 * cradius;
        //更新六棱柱每个侧面最低点最高点坐标值
        for (auto& hexa : hexas.polys) {
            for (int j = 0; j < 6; ++j) {
                if (hexa.edges[j].neighbor >= 0 && (!hexa.edges[j].canAdd)) {
                    auto& edge = hexa.edges[j];
                    auto& next = hexa.edges[(j + 1) % 6];
                    auto& oh = hexas.polys[hexa.edges[j].neighbor];
                    auto& oe = oh.edges[(j + 3) % 6];
                    auto& onext = oh.edges[(j + 4) % 6];
                    edge.lowHeight = std::max(edge.lowHeight, next.lowHeight);
                    edge.topHeight = std::min(edge.topHeight, next.topHeight);
                    oe.lowHeight = std::max(oe.lowHeight, onext.lowHeight);
                    oe.topHeight = std::min(oe.topHeight, onext.topHeight);
                    float lowerlimit = std::max(edge.lowHeight, oe.lowHeight);
                    float upperlimit = std::min(edge.topHeight, oe.topHeight);
                    float medium = (lowerlimit + upperlimit) * 0.5;
                    float scope = (upperlimit - lowerlimit) * 0.5;
                    if (medium < cheight && (cmax - medium) < scope) {
                        hexa.edges[j].canAdd = true;
                        oh.edges[(j + 3) % 6].canAdd = true;
                    } else if (medium >= cheight) {
                        float appro = std::round((medium - cheight) / cdelta);
                        float dist = std::abs((medium - cheight) - cdelta * appro);
                        if (dist + cradius < scope) {
                            hexa.edges[j].blower = false;
                            hexa.edges[j].canAdd = true;
                            oh.edges[(j + 3) % 6].canAdd = true;
                        }
                    }
                }
            }
        }
        return;
    }

    TriPolygons traitCurrentPolygons(const HexaPolygons& hexas, int index)
    {
        topomesh::TriPolygons polys;
        if (index >= 0 && index < (hexas.polys.size())) {
            polys.emplace_back(hexas.polys[index].poly);
        }
        return polys;
    }

    TriPolygons traitNeighborPolygons(const HexaPolygons& hexas, int index)
    {
        topomesh::TriPolygons polys;
        if (index >= 0 && index < (hexas.polys.size()))             {
            const auto& edges = hexas.polys[index].edges;
            int nums = edges.size();
            polys.reserve(nums);
            for (int i = 0; i < nums; ++i) {
                if (edges[i].neighbor >= 0)
                    polys.emplace_back(hexas.polys[edges[i].neighbor].poly);
            }
        }
        return polys;
    }

    TriPolygons traitDirctionPolygon(const HexaPolygons& hexas, int index, int dir)
    {
        topomesh::TriPolygons polys;
        if (index >= 0 && index < (hexas.polys.size())) {
            if (dir >= 0 && dir <= 5) {
                const auto& edges = hexas.polys[index].edges;
                int val = edges[dir].neighbor;
                if (val >= 0) {
                    polys.emplace_back(hexas.polys[val].poly);
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
        std::vector<trimesh::vec3> points;
        std::vector<trimesh::ivec3> faces;
        int hexaEdges = 0, num2 = 0;
        for (auto& hexa : hexas.polys) {
            hexa.startIndex = num2;
            const auto& poly = hexa.poly;
            for (int i = 0; i < poly.size(); ++i) {
                const auto& p = poly[i];
                points.emplace_back(p.x, p.y, hexa.edges[i].lowHeight);
                ++num2, ++hexaEdges;
            }
        }

        for (const auto& hexa : hexas.polys) {
            const auto& poly = hexa.poly;
            for (int i = 0; i < poly.size(); ++i) {
                const auto& p = poly[i];
                points.emplace_back(p.x, p.y, hexa.edges[i].topHeight);
                ++num2;
            }
        }
        //每条边的最低点最高点可能被更新
        GenerateHexagonNeighbors(hexas, param);
        const float cradius = hexas.side * param.ratio * 0.5f;
        const float cdelta = param.delta + 2 * cradius; ///<相邻两层圆心高度差
        const float cheight = param.height + cradius; ///<第一层圆心高度
        float cmax = param.height + 2.0 * cradius; ///<第一层圆最高点
        float addz = cdelta / 2.0; ///<相邻圆心高度差的一半
        trimesh::vec3 trans(0, 0, cdelta); ///相邻两圆平移量
        const int nslices = param.nslices;
        const int singlenums = nslices + 4;
        std::vector<int> holeStarts; ///<每个孔的第一个顶点索引
        int singleHoleEdges = 0;
        int mutiHoleEdges = 0;
        int holesnum = 0;
        int bottomfacenums = 0;
        ///六角网格底部
        if (hexas.bSewBottom) {
            ///两个网格中间矩形区域
            for (auto& hexa : hexas.polys) {
                for (int i = 0; i < 6; ++i) {
                    if ((!hexa.edges[i].addRect) && (hexa.edges[i].neighbor >= 0)) {
                        auto& oh = hexas.polys[hexa.edges[i].neighbor];
                        const int& start = hexa.startIndex;
                        const int& a = start + i;
                        const int& b = start + (i + 1) % 6;
                        const int& ostart = oh.startIndex;
                        const int& c = ostart + (i + 3) % 6;
                        const int& d = ostart + (i + 4) % 6;
                        faces.emplace_back(trimesh::ivec3(b, d, a));
                        faces.emplace_back(trimesh::ivec3(b, c, d));
                        hexa.edges[i].addRect = true;
                        oh.edges[(i + 3) % 6].addRect = true;
                        bottomfacenums += 2;
                    }
                }
            }
            ///三个网格中间部分
            for (auto& hexa : hexas.polys) {
                for (int i = 0; i < 6; ++i) {
                    if ((!hexa.edges[i].addTriangle) && (hexa.edges[i].neighbor >= 0) && (hexa.edges[(i + 1) % 6].neighbor >= 0)) {
                        auto& ohb = hexas.polys[hexa.edges[i].neighbor];
                        auto& ohc = hexas.polys[hexa.edges[(i + 1) % 6].neighbor];
                        const int& a = hexa.startIndex + (i + 1) % 6; ///<线段终点
                        const int& b = ohb.startIndex + (i + 3) % 6; ///<线段始点
                        const int& c = ohc.startIndex + (i + 4 + 1) % 6; ///<线段终点
                        faces.emplace_back(trimesh::ivec3(a, c, b));
                        hexa.edges[i].addTriangle = true;
                        ohb.edges[(i + 3) % 6].addTriangle = true;
                        ohc.edges[(i + 4) % 6].addTriangle = true;
                        ++bottomfacenums;
                    }
                }
            }
        }
        ///添加圆孔的顶点
        for (auto& hexa : hexas.polys) {
            for (int i = 0; i < 6; ++i) {
                if (hexa.edges[i].canAdd && (!hexa.edges[i].hasAdd) && (hexa.edges[i].neighbor >= 0)) {
                    auto& edge = hexa.edges[i];
                    auto& next = hexa.edges[(i + 1) % 6];
                    auto& oh = hexas.polys[hexa.edges[i].neighbor];
                    auto& oe = oh.edges[(i + 3) % 6];
                    auto& onext = oh.edges[(i + 4) % 6];
                    float lowerlimit = std::max(edge.lowHeight, oe.lowHeight);
                    float upperlimit = std::min(edge.topHeight, oe.topHeight);
                    float medium = (lowerlimit + upperlimit) * 0.5;
                    float scope = (upperlimit - lowerlimit) * 0.5;
                    const auto & a = hexa.poly[i];
                    const auto & b = hexa.poly[(i + 1) % 6];
                    const auto & c = oh.poly[(i + 3) % 6];
                    const auto & d = oh.poly[(i + 4) % 6];
                    if (medium < cheight && (cmax - medium) < scope) {
                        ///此时最高点下方有几个圆就可增加几个圆
                        int n = std::floor((upperlimit - cmax) / cdelta);
                        edge.holenums = n + 1, oe.holenums = n + 1;
                        edge.starts.reserve((size_t)n + 1);
                        oe.starts.reserve((size_t)n + 1);
                        const auto & ab = (a + b) * 0.5;
                        const auto & c1 = trimesh::vec3(ab.x, ab.y, cheight);
                        std::vector<int> corner1;
                        auto & poly1 = traitPlanarCircle(c1, cradius, corner1, a - b, nslices);
                        edge.corners.swap(corner1);
                        edge.starts.emplace_back(points.size());
                        holeStarts.emplace_back(points.size());
                        points.insert(points.end(), poly1.begin(), poly1.end());

                        const auto & cd = (c + d) * 0.5;
                        const auto & c2 = trimesh::vec3(cd.x, cd.y, cheight);
                        std::vector<int> corner2;
                        auto & poly2 = traitPlanarCircle(c2, cradius, corner2, c - d, nslices);
                        oe.corners.swap(corner2);
                        oe.starts.emplace_back(points.size());
                        holeStarts.emplace_back(points.size());
                        points.insert(points.end(), poly2.begin(), poly2.end());
                        holesnum += 2;
                        if (n < 1) singleHoleEdges += 2;
                        else {
                            mutiHoleEdges += 2;
                            edge.bmutihole = true;
                            oe.bmutihole = true;
                            edge.addPoints.reserve(n);
                            oe.addPoints.reserve(n);
                        }
                        for (int j = 1; j <= n; ++j) {
                            const auto& z = cheight + float(2 * j - 1) * addz;
                            edge.addPoints.emplace_back(points.size());
                            points.emplace_back(trimesh::vec3(a.x, a.y, z));
                            points.emplace_back(trimesh::vec3(b.x, b.y, z));
                            oe.addPoints.emplace_back(points.size());
                            points.emplace_back(trimesh::vec3(c.x, c.y, z));
                            points.emplace_back(trimesh::vec3(d.x, d.y, z));

                            edge.starts.emplace_back(points.size());
                            holeStarts.emplace_back(points.size());
                            translateTriPolygon(poly1, trans);
                            points.insert(points.end(), poly1.begin(), poly1.end());
                            oe.starts.emplace_back(points.size());
                            holeStarts.emplace_back(points.size());
                            translateTriPolygon(poly2, trans);
                            points.insert(points.end(), poly2.begin(), poly2.end());
                            holesnum += 2;
                        }
                        hexa.edges[i].hasAdd = true;
                        oh.edges[(i + 3) % 6].hasAdd = true;
                    } else if (medium >= cheight) {
                        float ratio = (medium - cheight) / cdelta;
                        float appro = std::round(ratio);
                        float dist = std::abs((medium - cheight) - cdelta * appro);
                        if (dist + cradius < scope) {
                            ///最低点下方有几个圆孔
                            const int& n1 = lowerlimit < cmax ? 0 : ((lowerlimit - cmax) / cdelta + 1);
                            ///最高点下方有几个圆孔
                            const int& n2 = (upperlimit - cmax) / cdelta + 1;
                            const int& n = n2 - n1;
                            edge.holenums = n, oe.holenums = n;
                            edge.starts.reserve((size_t)n);
                            oe.starts.reserve((size_t)n);
                            const auto & ab = (a + b) * 0.5;
                            const auto & c1 = trimesh::vec3(ab.x, ab.y, cheight + n1 * cdelta);
                            std::vector<int> corner1;
                            auto & poly1 = traitPlanarCircle(c1, cradius, corner1, a - b, nslices);
                            edge.corners.swap(corner1);
                            edge.starts.emplace_back(points.size());
                            holeStarts.emplace_back(points.size());
                            points.insert(points.end(), poly1.begin(), poly1.end());

                            const auto & cd = (c + d) * 0.5;
                            const auto & c2 = trimesh::vec3(cd.x, cd.y, cheight + n1 * cdelta);
                            std::vector<int> corner2;
                            auto & poly2 = traitPlanarCircle(c2, cradius, corner2, c - d, nslices);
                            oe.corners.swap(corner2);
                            oe.starts.emplace_back(points.size());
                            holeStarts.emplace_back(points.size());
                            points.insert(points.end(), poly2.begin(), poly2.end());
                            holesnum += 2;
                            if (n <= 1) singleHoleEdges += 2;
                            else {
                                mutiHoleEdges += 2;
                                edge.bmutihole = true;
                                oe.bmutihole = true;
                                edge.addPoints.reserve(size_t(n) - 1);
                                oe.addPoints.reserve(size_t(n) - 1);
                            }
                            for (int j = n1; j < n2 - 1; ++j) {
                                const auto& z = cheight + float(2 * j + 1) * addz;
                                edge.addPoints.emplace_back(points.size());
                                points.emplace_back(trimesh::vec3(a.x, a.y, z));
                                points.emplace_back(trimesh::vec3(b.x, b.y, z));
                                oe.addPoints.emplace_back(points.size());
                                points.emplace_back(trimesh::vec3(c.x, c.y, z));
                                points.emplace_back(trimesh::vec3(d.x, d.y, z));

                                edge.starts.emplace_back(points.size());
                                holeStarts.emplace_back(points.size());
                                translateTriPolygon(poly1, trans);
                                points.insert(points.end(), poly1.begin(), poly1.end());
                                oe.starts.emplace_back(points.size());
                                holeStarts.emplace_back(points.size());
                                translateTriPolygon(poly2, trans);
                                points.insert(points.end(), poly2.begin(), poly2.end());
                                holesnum += 2;
                            }
                            hexa.edges[i].hasAdd = true;
                            oh.edges[(i + 3) % 6].hasAdd = true;

                        }
                    }
                }
            }
        }
        int noHoleEdges = hexaEdges - singleHoleEdges - mutiHoleEdges;
        int holeFaces = holesnum * singlenums;
        int rectfacenums = noHoleEdges * 2;
        int upperfacenums = hexaEdges;
        int allfacenums = holeFaces + rectfacenums + upperfacenums + bottomfacenums;
        faces.reserve(allfacenums);
        if (hexas.bSewTop) {
            for (int i = 0; i < hexas.polys.size(); ++i) {
                const auto& hexa = hexas.polys[i];
                const auto& poly = hexa.poly;
                int polynums = poly.size();
                const int& start = hexa.startIndex;
                ///六角网格顶部
                const int& upstart = start + hexaEdges;
                for (int j = 1; j < polynums - 1; ++j) {
                    const int& cur = upstart + j;
                    const int& next = upstart + j + 1;
                    faces.emplace_back(trimesh::ivec3(upstart, next, cur));
                }
            }
        }
        for (int i = 0; i < hexas.polys.size(); ++i) {
            const auto& hexa = hexas.polys[i];
            const auto& poly = hexa.poly;
            int polynums = poly.size();
            const int& start = hexa.startIndex;
            //六角网格侧面
            for (int j = 0; j < polynums; ++j) {
                ///侧面下方两个顶点索引
                int a = start + j;
                int b = start + (j + 1) % polynums;
                ///侧面上方两个顶点索引
                int c = a + hexaEdges;
                int d = b + hexaEdges;
                if (hexa.edges[j].canAdd) {
                    const auto& edge = hexa.edges[j];
                    const auto& corner = edge.corners;
                    int r = corner[0];
                    int s = corner[1];
                    int t = corner[2];
                    int u = corner[4];
                    if (nslices % 2 == 0) {
                        u = corner[3];
                    }
                    if (edge.bmutihole) {
                        //第一个圆孔的第一个顶点的索引
                        int cstart = edge.starts.front();
                        c = edge.addPoints.front(), d = c + 1;
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
                        //中间的圆孔构建面片
                        for (int num = 1; num < edge.holenums - 1; ++num) {
                            cstart = edge.starts[num];
                            a = c; b = d;
                            c = edge.addPoints[num]; d = c + 1;

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
                        //最上方圆孔构建面
                        cstart = edge.starts.back();
                        a = c; b = d;
                        c = start + j + hexaEdges;
                        d = start + (j + 1) % polynums + hexaEdges;
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
                        ///只有单孔的侧面
                        int cstart = edge.starts.front();
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
                    ///无孔的侧面
                    faces.emplace_back(trimesh::ivec3(b, a, c));
                    faces.emplace_back(trimesh::ivec3(b, c, d));
                }
            }
        }
        ///六角网格桥接的孔洞部
        for (int i = 0; i < holesnum; i += 2) {
            const int& start = holeStarts[i];
            for (int j = 0; j < nslices; ++j) {
                const int& lcur = start + j;
                const int& lnext = start + (j + 1) % nslices;
                const int& rcur = start + (nslices - j) % nslices + nslices;
                const int& rnext = start + (nslices - j - 1) % nslices + nslices;
                faces.emplace_back(trimesh::ivec3(lnext, rnext, rcur));
                faces.emplace_back(trimesh::ivec3(lnext, rcur, lcur));
            }
        }
        std::shared_ptr<trimesh::TriMesh> triMesh(new trimesh::TriMesh());
        triMesh->vertices.swap(points);
        triMesh->faces.swap(faces);
        //triMesh->write("holes.stl");
        return triMesh;
    }

    void findHoneyCombsCoord(trimesh::TriMesh* mesh, const honeyLetterOpt& honeycombs, std::vector<std::vector<std::pair<float, float>>>& coord)
    {
        std::vector<int> UpVertex;
        std::vector<int> BotVertex;
       /* mesh->need_normals();
        for (int vi = 0; vi < mesh->normals.size(); vi++)
        {
            trimesh::point VertexNormal = trimesh::normalized(mesh->normals[vi]);
            float dir = VertexNormal ^ trimesh::point(0, 0, 1);
            if (dir >= 0)
                UpVertex.push_back(vi);
            else
                BotVertex.push_back(vi);
        }*/

        for (int fi = 0; fi < mesh->faces.size(); fi++)
        {
            trimesh::point v1 = mesh->vertices[mesh->faces[fi][1]] - mesh->vertices[mesh->faces[fi][0]];
            trimesh::point v2 = mesh->vertices[mesh->faces[fi][2]] - mesh->vertices[mesh->faces[fi][0]];
            trimesh::point normal = v1 % v2;
            float dir = normal ^ trimesh::point(0, 0, 1);
            if (dir >= 0)
            {
                UpVertex.push_back(mesh->faces[fi][0]);
                UpVertex.push_back(mesh->faces[fi][1]);
                UpVertex.push_back(mesh->faces[fi][2]);
            }
            else
            {
                BotVertex.push_back(mesh->faces[fi][0]);
                BotVertex.push_back(mesh->faces[fi][1]);
                BotVertex.push_back(mesh->faces[fi][2]);
            }
        }
        std::sort(UpVertex.begin(), UpVertex.end());
        std::vector<int>::iterator itu = std::unique(UpVertex.begin(), UpVertex.end());
        UpVertex.resize(std::distance(UpVertex.begin(), itu));
        std::sort(BotVertex.begin(), BotVertex.end());
        std::vector<int>::iterator itb = std::unique(BotVertex.begin(), BotVertex.end());
        BotVertex.resize(std::distance(BotVertex.begin(), itb));
        
        const int width = 100, height = 100;
        mesh->need_bbox();
        float start_x = mesh->bbox.min.x;
        float start_y = mesh->bbox.min.y;
        float length_x = (mesh->bbox.max.x - mesh->bbox.min.x) / (width * 1.0f);
        float length_y = (mesh->bbox.max.y - mesh->bbox.min.y) / (height * 1.0f);

        std::vector<std::vector<float>> UpVertexMap(width,std::vector<float>(height,std::numeric_limits<float>::max()));
        std::vector<std::vector<float>> BotVertexMap(width, std::vector<float>(height,std::numeric_limits<float>::max()));
        for (int i = 0; i < UpVertex.size(); i++)
        {
            float x = mesh->vertices[UpVertex[i]].x - start_x;
            float y = mesh->vertices[UpVertex[i]].y - start_y;
            int xi = x / length_x;
            int yi = y / length_y;
            if (xi == 100) xi--;
            if (yi == 100) yi--;
            if(mesh->vertices[UpVertex[i]].z< UpVertexMap[xi][yi])
                UpVertexMap[xi][yi]=(mesh->vertices[UpVertex[i]].z);
        }
        for (int i = 0; i < BotVertex.size(); i++)
        {
            float x = mesh->vertices[BotVertex[i]].x - start_x;
            float y = mesh->vertices[BotVertex[i]].y - start_y;
            int xi = x / length_x;
            int yi = y / length_y;
            if (xi == 100) xi--;
            if (yi == 100) yi--;
            if (mesh->vertices[BotVertex[i]].z < BotVertexMap[xi][yi])
                BotVertexMap[xi][yi]=mesh->vertices[BotVertex[i]].z;
        }
        std::vector<std::vector<trimesh::vec2>> honeycombscoord;
        for (int i = 0; i < honeycombs.hexgons.size(); i++)
        {
            std::vector<trimesh::vec2> dim2;
            for (int j = 0; j < honeycombs.hexgons[i].borders.size(); j++)
                dim2.push_back(trimesh::vec2(honeycombs.hexgons[i].borders[j].x,honeycombs.hexgons[i].borders[j].y));
            honeycombscoord.emplace_back(dim2);
        }
        for (int i = 0; i < honeycombscoord.size(); i++)
        {
            std::vector<std::pair<float, float>> BeginAndEnd;
            for (int j = 0; j < honeycombscoord[i].size(); j++)
            {
                float x = honeycombscoord[i][j].x - start_x;
                float y = honeycombscoord[i][j].y - start_y;
                int xi = x / length_x;
                int yi = y / length_y;
                float top_z= UpVertexMap[xi][yi];
                float bot_z = BotVertexMap[xi][yi];

                if(top_z == std::numeric_limits<float>::max())
                {
                    for (int n = 1; n < width; n++)
                    {
                        for (int ii = xi - n; ii <= xi + n; ii++)
                        {
                            for (int jj = yi - n; jj <= yi + n; jj++)
                            {
                                if (std::abs(ii - xi) != n && std::abs(jj - yi) != n) continue;
                                if (ii < 0 || ii >= 100 || jj < 0 || jj <= 100)continue;
                                if (UpVertexMap[ii][jj] < top_z)
                                {
                                    top_z = UpVertexMap[ii][jj]; break;
                                }
                            }
                            if (top_z != std::numeric_limits<float>::max()) break;
                        }
                        if (top_z != std::numeric_limits<float>::max()) break;
                    }
                }
                if (top_z == std::numeric_limits<float>::max())
                    top_z = 50;
                bot_z = 0;
                if(bot_z == std::numeric_limits<float>::max())
                {
                    for (int n = 1; n < width; n++)
                    {
                        for (int ii = xi - n; ii <= xi + n; ii++)
                        {
                            for (int jj = yi - n; jj <= yi + n; jj++)
                            {
                                if (std::abs(ii - xi) != n && std::abs(jj - yi) != n) continue;
                                if (ii < 0 || ii >= 100 || jj < 0 || jj <= 100)continue;
                                if (BotVertexMap[ii][jj] < bot_z)
                                {
                                    bot_z = BotVertexMap[ii][jj]; break;
                                }
                            }
                            if (bot_z != std::numeric_limits<float>::max()) break;
                        }
                        if (bot_z != std::numeric_limits<float>::max()) break;
                    }
                }
                BeginAndEnd.push_back(std::make_pair(top_z, bot_z));
            }
            coord.push_back(BeginAndEnd);
        }
        //std::vector<std::vector<std::pair<trimesh::point, trimesh::point>>>& coord
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
}