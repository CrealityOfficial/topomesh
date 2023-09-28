#include "fillhoneycombs.h"
#include "topomesh/data/mmesht.h"
#include "topomesh/alg/letter.h"
#include "topomesh/alg/earclipping.h"
#include "topomesh/alg/subdivision.h"
#include "topomesh/alg/utils.h"
#include "trimesh2/TriMesh_algo.h"
#include "topomesh/alg/solidtriangle.h"
#include "internal/polygon/comb.h"
#include "internal/mesh/dumplicate.h"
#include "random"

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

    void GenerateTriPolygonsHexagons(const TriPolygons&polys, const HoneyCombParam& honeyparams, honeyLetterOpt& letterOpts, HoneyCombDebugger* debugger) {
        if (debugger) {
            //显示底面边界轮廓多边形
            debugger->onGenerateInfillPolygons(polys);
        }
        //第0步，底面边界轮廓抽壳
        const double resolution = honeyparams.resolution;
        const double side = honeyparams.honeyCombRadius;
        const double thickness = honeyparams.shellThickness;
        const double radius = side + honeyparams.nestWidth / SQRT3;
        Polygons bpolygons; ///<底面边界轮廓多边形
        Polygons mpolygons; ///<底面边界抽壳多边形
        {
            for (const auto& poly : polys) {
                ClipperLib::Path path;
                path.reserve(poly.size());
                for (const auto& p : poly) {
                    const auto& x = std::round(p.x / resolution);
                    const auto& y = std::round(p.y / resolution);
                    path.emplace_back(ClipperLib::IntPoint(x, y));
                }
                bpolygons.paths.emplace_back(std::move(path));
            }
            Polygons cpolygons(bpolygons);
            mpolygons = cpolygons.offset(-std::round(thickness / resolution));
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
        //第1步，底面边界抽壳同时，生成六角网格阵列
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
        const auto& xdist = 3.0 / 2.0 * radius;
        const auto& ydist = SQRT3 / 2.0 * radius;
        trimesh::vec3 origin(min.x, max.y, 0);

        HexagonArrayParam hexagonparams;
        hexagonparams.pos = origin;
        hexagonparams.radius = radius;
        hexagonparams.nestWidth = honeyparams.nestWidth;
        hexagonparams.ncols = std::ceil((max.x - min.x) / xdist);
        hexagonparams.nrows = (std::ceil((max.y - min.y) / ydist) + 1) / 2;
        HexaPolygons hexagons = GenerateHexagonsGridArray(hexagonparams);
        if (debugger) {
            //显示六角网格阵列
            debugger->onGenerateBottomPolygons(convertFromHexaPolygons(hexagons));
        }
        //第2步，计算网格阵列与抽壳轮廓交集
        std::vector<Polygons> hpolygons; ///<完整六角网格多边形序列
        std::vector<Polygons> ipolygons; ///<初次求交网格多边形序列
        std::vector<HexaPolygon> ihexagons;
        ihexagons.reserve(hexagons.polys.size());
        {
            hpolygons.reserve(hexagons.polys.size());
            for (auto& hexa : hexagons.polys) {
                Polygons polygons;
                ClipperLib::Path path;
                path.reserve(hexa.poly.size());
                for (const auto& p : hexa.poly) {
                    const auto& x = std::round(p.x / resolution);
                    const auto& y = std::round(p.y / resolution);
                    path.emplace_back(ClipperLib::IntPoint(x, y));
                }
                bool belong = true;
                for (const auto& p : path) {
                    if (!mpolygons.inside(p)) {
                        belong = false;
                        break;
                    }
                }
                polygons.paths.emplace_back(std::move(path));
                if (belong) {
                    hexa.standard = true;
                    ipolygons.emplace_back(std::move(polygons));
                    ihexagons.emplace_back(std::move(hexa));
                } else {
                    Polygons&& ipolys = mpolygons.intersection(polygons);
                    if (!ipolys.empty()) {
                        if (ipolys.area() < 0) {
                            ClipperLib::Path& path = ipolys.paths.front();
                            ClipperLib::Path tmp;
                            tmp.swap(path);
                            std::reverse(tmp.begin(), tmp.end());
                            ipolys.paths.emplace_back(std::move(tmp));
                        }

                        if (debugger) {
                            AABB box(polygons);
                            box.expand(50000);
                            std::string filename = "hexagon.svg";
                            SVG svg(filename, box, 0.01);
                            svg.writePolygons(mpolygons, SVG::Color::BLACK, 2);
                            svg.writePolygons(polygons, SVG::Color::RED, 2);
                            svg.writePolygons(ipolys, SVG::Color::BLUE, 2);
                            svg.writePoints(polygons, false, 2, SVG::Color::RAINBOW);
                            svg.writePoints(ipolys, true, 3, SVG::Color::GREEN);
                        }
                        auto&& edgemaps = GetHexagonEdgeMap(ipolys, polygons.paths[0], hexagons.side, resolution);
                        hexa.standard = false;
                        hexa.p2hPointMap.swap(edgemaps[0]);
                        hexa.h2pPointMap.swap(edgemaps[1]);
                        hexa.p2hEdgeMap.swap(edgemaps[2]);
                        hexa.h2pEdgeMap.swap(edgemaps[3]);
                        ipolygons.emplace_back(std::move(ipolys));
                        ihexagons.emplace_back(std::move(hexa));
                    }
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
        //第3步，筛选符合需求的网格多边形
        TriPolygons outtripolys; ///<最终保留求交网格多边形序列
        const double minrate = honeyparams.keepHexagonRate;
        const double keeprate = 3.0 / 2.0 * SQRT3 * std::pow(side / resolution, 2) * minrate;
        const double keeparea = honeyparams.keepHexagonArea / std::pow(resolution, 2);
        const double minhexarea = std::min(keeprate, keeparea);
        std::vector<HexaPolygon> ohexagons;
        const int isize = ipolygons.size();
        ohexagons.reserve(isize);
        for (int i = 0; i < isize; i++) {
            TriPolygon tripolys;
            Polygons ipolys = ipolygons[i];
            ClipperLib::Path& path = ipolys.paths.front();
            if (ipolys.area() >= minhexarea) {
                for (const auto& p : path) {
                    const auto& point = trimesh::vec3(p.X * resolution, p.Y * resolution, 0);
                    tripolys.emplace_back(std::move(point));
                }
                outtripolys.emplace_back(tripolys);
                ihexagons[i].poly.swap(tripolys);
                ohexagons.emplace_back(ihexagons[i]);
            }
        }
        if (debugger) {
            //显示最终符合需求的网格多边形
            debugger->onGenerateBottomPolygons(outtripolys);
        }
        const int osize = ohexagons.size();
        letterOpts.hexgons.reserve(osize);
        letterOpts.side = side;
        if (honeyparams.bKeepHexagon) {
            for (int i = 0; i < osize; ++i) {
                HexaPolygon& h = ohexagons[i];
                Hexagon hexagon(h.center, hexagons.side);
                h.hexagon.swap(hexagon.border);
                Polygons polygons = SaveTriPolygonToPolygons(h.hexagon);
                Polygons ipolys = SaveTriPolygonToPolygons(h.poly);
                AABB box(polygons);
                box.expand(50000);
                std::string filename = "test/hexagon" + std::to_string(i) + ".svg";
                SVG svg(filename, box, 0.01);
                svg.writePolygons(mpolygons, SVG::Color::BLACK, 2);
                svg.writePolygons(polygons, SVG::Color::RED, 2);
                svg.writePolygons(ipolys, SVG::Color::BLUE, 2);
                svg.writePoints(polygons, false, 2, SVG::Color::RAINBOW);
                svg.writePoints(ipolys, true, 3, SVG::Color::GREEN);
            }
        }
        for (const auto& hexa : ohexagons) {
            letterOpts.hexgons.emplace_back(std::move(hexa));
        }
        return;
    }

    void GenerateBottomHexagons(CMesh& honeyMesh, const HoneyCombParam& honeyparams, honeyLetterOpt& letterOpts, HoneyCombDebugger* debugger)
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
                border.emplace_back(p);
            }
            polys.emplace_back(border);
        }
        //第6步，在底面边界轮廓多边形内生成蜂窝六边形
        GenerateTriPolygonsHexagons(polys, honeyparams, letterOpts, debugger);
        honeyMesh.mbox = cutMesh.mbox;
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
             return trimesh::point(0,0,1);
    }

    trimesh::TriMesh* findOutlineOfDir(trimesh::TriMesh* mesh,std::vector<int>& botfaces)
    {
        mesh->need_adjacentfaces();
        mesh->need_across_edge();
        int index;
        float min_z = std::numeric_limits<float>::max();
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            if (mesh->vertices[i].z < min_z)
            {
                min_z = mesh->vertices[i].z;
                index = i;
            }
        }
        int faceindex;
        for (int i = 0; i < mesh->adjacentfaces[index].size(); i++)
        {
            int face = mesh->adjacentfaces[index][i];
            trimesh::point v1 = mesh->vertices[mesh->faces[face][1]] - mesh->vertices[mesh->faces[face][0]];
            trimesh::point v2 = mesh->vertices[mesh->faces[face][2]] - mesh->vertices[mesh->faces[face][0]];
            trimesh::point nor = v1 % v2;
            float arc = nor ^ trimesh::point(0, 0, 1);
            if (arc < 0)
            {
                faceindex = face;
                break;
            }
        }
        std::vector<int> botface = { faceindex };
        std::vector<int> vis(mesh->faces.size(), 0);      
        std::queue<int> facequeue;
        facequeue.push(faceindex);
        while (!facequeue.empty())
        {
            vis[facequeue.front()] = 1;
            for (int i = 0; i < mesh->across_edge[facequeue.front()].size(); i++)
            {
                int face = mesh->across_edge[facequeue.front()][i];
                if (vis[face]) continue;
                trimesh::point v1 = mesh->vertices[mesh->faces[face][1]] - mesh->vertices[mesh->faces[face][0]];
                trimesh::point v2 = mesh->vertices[mesh->faces[face][2]] - mesh->vertices[mesh->faces[face][0]];
                trimesh::point nor = v1 % v2;
                float arc = trimesh::normalized(nor) ^ trimesh::point(0, 0, 1);
                arc = arc >= 1.f ? 1.f : arc;
                arc = arc <= -1.f ? -1.f : arc;
                float ang = std::acos(arc) * 180 / M_PI;
                if (arc<0/*&&ang >120*/)                                
                {
                    vis[face] = 1;
                    facequeue.push(face);
                    botface.push_back(face);
                }
            }
            facequeue.pop();
        }

       
        trimesh::TriMesh* newmesh=new trimesh::TriMesh(*mesh);
        std::vector<bool> deleteFace(mesh->faces.size(), true);
        for (int i = 0; i < botface.size(); i++)
            deleteFace[botface[i]] = false; 
        trimesh::remove_faces(newmesh, deleteFace);
        trimesh::remove_unused_vertices(newmesh);
        for (trimesh::point& v : newmesh->vertices)
            v = trimesh::point(v.x, v.y, 0);
        botfaces.swap(botface);
        /*newmesh->write("removeface.ply");
        mesh->write("mesh.ply");*/
        return newmesh;
    }

    std::shared_ptr<trimesh::TriMesh> GenerateHoneyCombs(trimesh::TriMesh* trimesh , const HoneyCombParam& honeyparams ,
        ccglobal::Tracer* tracer , HoneyCombDebugger* debugger )
    {
        //0.初始化Cmesh,并旋转
        if(trimesh == nullptr) return nullptr;

        CMesh cmesh(trimesh);
        //1.重新整理输入参数
        /*trimesh::vec3 dir=adjustHoneyCombParam(trimesh, honeyparams);
        cmesh.Rotate(dir, trimesh::vec3(0, 0, 1));*/
        //2.找到生成蜂窝指定的区域（自定义或者是用户自己指定）          
        if (honeyparams.faces.empty()) {
            //自定义最大底面朝向
            if (honeyparams.axisDir == trimesh::vec3(0,0,0))
            {
                honeyLetterOpt letterOpts;
                //inputMesh.WriteSTLFile("inputmesh");
                //第1步，寻找底面（最大平面）朝向
                std::vector<int>bottomFaces;
                trimesh::vec3 dir = cmesh.FindBottomDirection(&bottomFaces);
                cmesh.Rotate(dir, trimesh::vec3(0, 0, -1));
                cmesh.GenerateBoundBox();
                const auto minPt = cmesh.mbox.min;
                cmesh.Translate(-minPt);
                letterOpts.bottom.resize(bottomFaces.size());
                letterOpts.bottom.assign(bottomFaces.begin(), bottomFaces.end());
                std::vector<int> honeyFaces;
                honeyFaces.reserve(cmesh.mfaces.size());
                for (int i = 0; i < cmesh.mfaces.size(); ++i) {
                    honeyFaces.emplace_back(i);
                }
                std::sort(bottomFaces.begin(), bottomFaces.end());
                std::vector<int> otherFaces(honeyFaces.size() - bottomFaces.size());
                std::set_difference(honeyFaces.begin(), honeyFaces.end(), bottomFaces.begin(), bottomFaces.end(), otherFaces.begin());
                letterOpts.others = std::move(otherFaces);
                //第2步，平移至xoy平面后底面整平
                cmesh.FlatBottomSurface(&bottomFaces);
                if (debugger) {
                    //显示底面区域
                    TriPolygons polygons;
                    const auto& inPoints = cmesh.mpoints;
                    const auto& inIndexs = cmesh.mfaces;
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
                //第3步，生成底面六角网格
                GenerateBottomHexagons(cmesh, honeyparams, letterOpts, debugger);
                trimesh::TriMesh&& mesh = cmesh.GetTriMesh();            
                trimesh = &mesh; 
                trimesh->bbox.valid = false;
                trimesh->need_bbox();             
                int row = 400;
                int col = 400;
                std::vector<std::tuple<trimesh::point, trimesh::point, trimesh::point>> Upfaces;
                for (int i : letterOpts.others)
                    Upfaces.push_back(std::make_tuple(trimesh->vertices[trimesh->faces[i][0]], trimesh->vertices[trimesh->faces[i][1]],
                        trimesh->vertices[trimesh->faces[i][2]]));
                topomesh::SolidTriangle upST(&Upfaces, row, col, trimesh->bbox.max.x, trimesh->bbox.min.x, trimesh->bbox.max.y, trimesh->bbox.min.y);
                upST.work();    
            
                HexaPolygons hexpolys;
                hexpolys.side = letterOpts.side;
                trimesh::point max_xy = trimesh->bbox.max;
                trimesh::point min_xy = trimesh->bbox.min;
                float lengthx = (trimesh->bbox.max.x - trimesh->bbox.min.x)/ (col*1.f);
                float lengthy = (trimesh->bbox.max.y - trimesh->bbox.min.y)/ (row*1.f);
               
                for (auto& hg : letterOpts.hexgons)
                {
                    std::vector<float> height;                                   
                    for (int i = 0; i < hg.poly.size(); i++)
                    {
                        trimesh::point p = hg.poly[i] - min_xy;
                        int xi = p.x / lengthx;
                        int yi = p.y / lengthy;
                        xi = xi == col ? --xi : xi;
                        yi = yi == row ? --yi : yi;                      
                        float min_z = upST.getDataMinZInterpolation(xi, yi);
                        if (min_z != std::numeric_limits<float>::max())
                        {
                            min_z -= honeyparams.shellThickness;
                            height.push_back(min_z);                                                  
                        }
                        else
                        {
                            height.push_back(0.f);                          
                        }
                       
                    }
#if false
                    for (int i = 0; i < height.size(); i++)
                    {
                        if (height[i] == 0.f) continue;
                        int next = (i + 1) % height.size();
                        trimesh::point c = hg.poly[next] - hg.poly[i];
                        float ch = height[next] - height[i];
                        int cx = coord[next].x - coord[i].x;
                        int cy = coord[next].y - coord[i].y;
                        if (cx == 0 && cy == 0) continue;
                        cx = std::abs(cx);
                        cy = std::abs(cy);
                        if (cx > cy)
                        {
                            cx += 1;
                            trimesh::point t = c / (cx*1.f);
                            float th = ch / (cx * 1.f);
                            float max_c = std::numeric_limits<float>::max();
                            for (int j = 1; j < cx;j++)
                            {
                                trimesh::point tt = hg.poly[i] + j*1.f * t;
                                float min_z = upST.getDataMinZ(tt.x,tt.y)- honeyparams.shellThickness;
                                float temp_z = height[i] + j * 1.f * th;
                                if (temp_z > min_z)
                                {
                                    if (max_c > (temp_z - min_z))
                                        max_c = (temp_z - min_z);
                                }
                            }
                            if (max_c != std::numeric_limits<float>::max())
                            {
                                if(height[next]> max_c)
                                    height[next] -= max_c;
                                if(height[i]>max_c)
                                    height[i] -= max_c;
                            }
                        }
                        else
                        {
                            cy += 1;
                            trimesh::point t = c / (cy * 1.f);
                            float th = ch / (cy * 1.f);
                            float max_c = std::numeric_limits<float>::max();
                            for (int j = 1; j < cy; j++)
                            {
                                trimesh::point tt = hg.poly[i] + j * 1.f * t;
                                float min_z = upST.getDataMinZ(tt.x, tt.y) - honeyparams.shellThickness;
                                float temp_z = height[i] + j * 1.f * th;
                                if (temp_z > min_z)
                                {
                                    if (max_c > (temp_z - min_z))
                                        max_c = (temp_z - min_z);
                                }
                            }
                            if (max_c != std::numeric_limits<float>::max())
                            {
                                if (height[next] > max_c)
                                    height[next] -= max_c;
                                if (height[i] > max_c)
                                    height[i] -= max_c;
                            }
                        }
                    }
#elif false
                    float last_z = std::numeric_limits<float>::max();
                    for (int i = min_yi; i <= max_yi; i++)
                    {
                        for (int j = min_xi; j <= max_xi; j++)
                        {
                            float min_z = upST.getDataMinZCoord(j, i);
                           // if (min_z == std::numeric_limits<float>::max()) min_z = 0.0f;
                            if (min_z < last_z)
                                last_z = min_z;                           
                        }
                    }
                    if(last_z>= (honeyparams.shellThickness+1.2f))
                       

                    last_z -= (honeyparams.shellThickness+1.2f);
                    last_z = last_z < 0.f ? 0.2f : last_z;
#endif
                    hg.edges.resize(hg.poly.size());
                    for (int i = 0; i < hg.edges.size(); i++)
                    {                       
                         hg.edges[i].topHeight = height[i];
                       // hg.edges[i].topHeight = last_z;
                    }
                    hexpolys.polys.push_back(hg);                   
                }
               

                topomesh::ColumnarHoleParam columnParam;
                columnParam.nslices = honeyparams.nslices;
                columnParam.ratio = honeyparams.ratio;
                columnParam.height = honeyparams.cheight;
                columnParam.delta = honeyparams.delta;
                columnParam.holeConnect = honeyparams.holeConnect;
                std::shared_ptr<trimesh::TriMesh> newmesh(topomesh::generateHolesColumnar(hexpolys, columnParam));
               

                std::vector<int> topfaces;
                topfaces.swap(hexpolys.topfaces);
                /*for (int fi = 0; fi < topfaces.size(); fi++)
                {
                    trimesh::point v0 = newmesh->vertices[newmesh->faces[topfaces[fi]][0]];
                    trimesh::point v1 = newmesh->vertices[newmesh->faces[topfaces[fi]][1]];
                    trimesh::point v2 = newmesh->vertices[newmesh->faces[topfaces[fi]][2]];
                    trimesh::vec3 a = v0 - v1;
                    trimesh::vec3 b = v1 - v2;
                    float are= sqrt(pow((a.y * b.z - a.z * b.y), 2) + pow((a.z * b.x - a.x * b.z), 2)
                        + pow((a.x * b.y - a.y * b.x), 2)) / 2.0f;
                    if (are < 0.1f)
                    {
                        topfaces.erase(topfaces.begin() + fi);
                        fi--;
                    }
                }
                */
                                        
                std::vector<int> outfaces;
               // topomesh::SimpleMidSubdivision(newmesh.get(), topfaces);
                for (int l = 0; l < 2; l++)
                {
                    topomesh::loopSubdivision(newmesh.get(), topfaces, outfaces);
                    outfaces.swap(topfaces);
                    outfaces.clear();
                }                         
                for (int fi : topfaces)
                {
                    for (int vi = 0; vi < 3; vi++)
                    {
                        float min_x = newmesh->vertices[newmesh->faces[fi][vi]].x - honeyparams.shellThickness * 3.f /5.f;
                        float min_y = newmesh->vertices[newmesh->faces[fi][vi]].y - honeyparams.shellThickness * 3.f /5.f;
                        float max_x = newmesh->vertices[newmesh->faces[fi][vi]].x + honeyparams.shellThickness * 3.f /5.f;
                        float max_y = newmesh->vertices[newmesh->faces[fi][vi]].y + honeyparams.shellThickness * 3.f /5.f;
                        int min_xi = (min_x - min_xy.x) / lengthx;
                        int min_yi = (min_y - min_xy.y) / lengthy;
                        int max_xi = (max_x - min_xy.x) / lengthx;
                        int max_yi = (max_y - min_xy.y) / lengthy;
                        float min_z = std::numeric_limits<float>::max();
                        for(int xii=min_xi;xii<=max_xi;xii++)
                            for (int yii = min_yi; yii <= max_yi; yii++)
                            {
                                float zz = upST.getDataMinZInterpolation(xii,yii);                              
                                if (zz < min_z)
                                    min_z = zz;
                            }
                       // float min_z = upST.getDataMinZ(newmesh->vertices[newmesh->faces[fi][vi]].x, newmesh->vertices[newmesh->faces[fi][vi]].y)/*- honeyparams.shellThickness*/;
                       // float min_z = upST.getDataMinZInterpolation(newmesh->vertices[newmesh->faces[fi][vi]].x, newmesh->vertices[newmesh->faces[fi][vi]].y) - 1.2*honeyparams.shellThickness;
                        if (min_z > 1.2*honeyparams.shellThickness)
                            min_z -= 1.2*honeyparams.shellThickness;
                        newmesh->vertices[newmesh->faces[fi][vi]] = trimesh::point(newmesh->vertices[newmesh->faces[fi][vi]].x, newmesh->vertices[newmesh->faces[fi][vi]].y, min_z);
                       // pointmesh->vertices.push_back(trimesh::point(newmesh->vertices[newmesh->faces[fi][vi]].x, newmesh->vertices[newmesh->faces[fi][vi]].y, min_z));
                                             
                    }
                }
               
                
#if 0
                std::vector<std::vector<int>> sequentials;
                getMeshBoundarys(*trimesh, sequentials);
                std::vector<std::vector<int>> sequentials2;
                getMeshBoundarys(*newmesh, sequentials2);

                if (sequentials.size() != sequentials2.size())
                    return nullptr;

                for (int m = 0; m < sequentials.size(); m++)
                {
                    trimesh::TriMesh* pointmesh1 = new trimesh::TriMesh();
                    for (int vii = 0; vii < sequentials[m].size(); vii++)
                    {
                        pointmesh1->vertices.push_back(trimesh->vertices[sequentials[m][vii]]);
                    }
                    int linesize = pointmesh1->vertices.size();
                    for (int vii = 0; vii < sequentials2[m].size(); vii++)
                    {
                        pointmesh1->vertices.push_back(newmesh->vertices[sequentials2[m][vii]]);
                    }
                    const int lineNum = 8;
                    std::vector<float> min_dist(lineNum, std::numeric_limits<float>::max());
                    std::vector<int> f_index;
                    for (int nn = 0; nn < lineNum; nn++)
                    {
                        int ln = (int)(linesize * nn / lineNum);
                        f_index.push_back(ln);
                    }
                    std::vector<int> near_index(lineNum, 0);
                    for (int vi = linesize; vi < pointmesh1->vertices.size(); vi++)
                    {
                        for (int nn = 0; nn < lineNum; nn++)
                        {
                            float dist = trimesh::distance(pointmesh1->vertices[f_index[nn]], pointmesh1->vertices[vi]);
                            if (dist < min_dist[nn])
                            {
                                near_index[nn] = vi;
                                min_dist[nn] = dist;
                            }
                        }
                    }
                    std::vector<std::vector<int>> sequent(lineNum, std::vector<int>());
                    //std::cout << "linesize : " << linesize << "\n";
                    for (int nn = 0; nn < lineNum; nn++)
                    {
                        //std::cout << "f_index : " << f_index[nn]<<" near_index : "<< near_index[nn]<< "\n";
                        if (nn != lineNum - 1)
                        {
                            for (int vi = f_index[nn]; vi <= f_index[nn + 1]; vi++)
                                sequent[nn].push_back(vi);
                            if (near_index[nn + 1] < near_index[nn])
                                for (int vi = near_index[nn + 1]; vi <= near_index[nn]; vi++)
                                    sequent[nn].push_back(vi);
                            else
                            {
                                for (int vi = near_index[nn + 1]; vi < pointmesh1->vertices.size(); vi++)
                                    sequent[nn].push_back(vi);
                                for (int vi = linesize; vi <= near_index[nn]; vi++)
                                    sequent[nn].push_back(vi);
                            }
                        }
                        else
                        {
                            for (int vi = f_index[nn]; vi < linesize; vi++)
                                sequent[nn].push_back(vi);
                            sequent[nn].push_back(f_index[0]);
                            if (near_index[0] < near_index[nn])
                                for (int vi = near_index[0]; vi <= near_index[nn]; vi++)
                                    sequent[nn].push_back(vi);
                            else
                            {
                                for (int vi = near_index[0]; vi < pointmesh1->vertices.size(); vi++)
                                    sequent[nn].push_back(vi);
                                for (int vi = linesize; vi <= near_index[nn]; vi++)
                                    sequent[nn].push_back(vi);
                            }
                        }
                    }
                    for (int nn = 0; nn < lineNum; nn++)
                    {
                        /* for (int ni = 0; ni < sequent[nn].size(); ni++)
                             std::cout << " " << sequent[nn][ni] << " ";
                         std::cout << "\n";*/
                        std::vector<std::pair<trimesh::point, int>> lines;
                        for (int vi = 0; vi < sequent[nn].size(); vi++)
                        {
                            lines.push_back(std::make_pair(pointmesh1->vertices[sequent[nn][vi]], sequent[nn][vi]));
                        }
                        topomesh::EarClipping earclip(lines);
                        std::vector<trimesh::ivec3> result = earclip.getResult();
                        for (int fi = 0; fi < result.size(); fi++)
                        {
                            result[fi] = trimesh::ivec3(result[fi].x, result[fi].z, result[fi].y);
                            pointmesh1->faces.push_back(result[fi]);
                        }
                    }

                   // std::string name = "pointmesh" + std::to_string(m) + ".ply";
                   // pointmesh1->write(name);

                    int vertexsize = newmesh->vertices.size();
                    for (int vi = 0; vi < pointmesh1->vertices.size(); vi++)
                        newmesh->vertices.push_back(pointmesh1->vertices[vi]);
                    for (int fi = 0; fi < pointmesh1->faces.size(); fi++)
                        newmesh->faces.push_back(trimesh::TriMesh::Face(pointmesh1->faces[fi][0] + vertexsize, pointmesh1->faces[fi][1] + vertexsize, pointmesh1->faces[fi][2] + vertexsize));

                }
#endif               
                JointBotMesh(trimesh,newmesh.get(), bottomFaces);
                trimesh::trans(newmesh.get(), minPt);
                trimesh::apply_xform(newmesh.get(), trimesh::xform::rot_into(trimesh::vec3(0, 0, -1), dir));
                //trimesh->write("trimesh.ply");
                //newmesh->write("holesColumnar.stl");
                return newmesh;
            }           
        } 
        else {
            //user indication faceindex
            trimesh::vec3 dir = honeyparams.axisDir;
            trimesh::apply_xform(trimesh, trimesh::xform::rot_into(dir, trimesh::vec3(0, 0, 1)));


            
        }
              
        return nullptr;
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
        //第3步，底面生成六角网格
        GenerateBottomHexagons(inputMesh, honeyparams, letterOpts, debugger);
        trimesh::TriMesh&& mesh = inputMesh.GetTriMesh();

        //std::vector<bool> delectface(mesh.faces.size(), false);
        //for (int i : bottomFaces)
        //    delectface[i] = true;
        //trimesh::remove_faces(&mesh, delectface); 

        //mesh.need_bbox();
        //std::vector<std::tuple<trimesh::point, trimesh::point, trimesh::point>> data;
        //for (trimesh::TriMesh::Face& f : mesh.faces)
        //    data.push_back(std::make_tuple(mesh.vertices[f[0]], mesh.vertices[f[1]], mesh.vertices[f[2]]));
        //topomesh::SolidTriangle sd(&data, 100, 100, mesh.bbox.max.x, mesh.bbox.min.x, mesh.bbox.max.y, mesh.bbox.min.y);
        //sd.work();
        //std::vector<std::vector<float>> result = sd.getResult();
        //trimesh::TriMesh* newmesh = new trimesh::TriMesh();
        //

        //for(int i=0;i<result.size();i++)
        //    for (int j = 0; j < result[i].size(); j++)
        //    {
        //        if (result[i][j] == std::numeric_limits<float>::max()) { result[i][j] = 0.f;/* std::cout << "i :" << i << " j : " << j << "\n";*/ }
        //        newmesh->vertices.push_back(trimesh::point(i, j, result[i][j]));
        //    }
    
        //newmesh->write("newhoneymesh.ply");   
        //mesh.write("honeyhole.ply");
        //return &mesh;
		mesh.need_adjacentfaces();
		mesh.need_neighbors();
		MMeshT mt(&mesh);
		std::vector<std::vector<trimesh::vec2>> honeycombs;
		for (int i = 0; i < letterOpts.hexgons.size(); i++)		
		{
			std::vector<trimesh::vec2> tmp;
			for (int j = 0; j < letterOpts.hexgons[i].poly.size(); j++)
			{
				tmp.push_back(trimesh::vec2(letterOpts.hexgons[i].poly[j].x, letterOpts.hexgons[i].poly[j].y));
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


    void JointBotMesh(trimesh::TriMesh* mesh, trimesh::TriMesh* newmesh,  std::vector<int>& botfaces)
    {
#if 0
        std::vector<bool> deletefaces(mesh->faces.size(), false);
        for (int f = 0; f < botfaces.size(); f++)
        {
            deletefaces[botfaces[f]] = true;
        }
        trimesh::remove_faces(mesh, deletefaces);
        trimesh::remove_unused_vertices(mesh);
        int facesize = newmesh->faces.size();
        int vertexsize = newmesh->vertices.size();
        for (int vi = 0; vi < mesh->vertices.size(); vi++)
            newmesh->vertices.push_back(mesh->vertices[vi]);
        for (int fi = 0; fi < mesh->faces.size(); fi++)
            newmesh->faces.push_back(trimesh::TriMesh::Face(mesh->faces[fi][0] + vertexsize, mesh->faces[fi][1] + vertexsize, mesh->faces[fi][2] + vertexsize));


        topomesh::MMeshT joinmesh(81960,81960);
        std::vector<std::vector<trimesh::point>> sequentials=GetOpenMeshBoundarys(*newmesh);
        std::vector<std::pair<float,int>> arc_array;
        std::vector<std::pair<int, int>> begin_and_end;
        std::vector<bool> is_reverse(sequentials.size(),false);
        for (int li = 0; li < sequentials.size(); li++)
        {
            float arc = 0.f;
            begin_and_end.push_back(std::make_pair(joinmesh.vertices.size(), joinmesh.vertices.size()+ sequentials[li].size()));
            for (int vi = 0; vi < sequentials[li].size(); vi++)
            {
                joinmesh.appendVertex(sequentials[li][vi]);
                int next_vi = (vi + 1) % sequentials[li].size();
                arc += (sequentials[li][vi].x * sequentials[li][next_vi].y - sequentials[li][vi].y * sequentials[li][next_vi].x);
            }
            if (arc < 0) is_reverse[li] = true;
            arc_array.push_back(std::make_pair(std::abs(arc),li));           
        }

        auto arc_sort = [](std::pair<float, int>& a, std::pair<float, int>& b) {
            return a.first > b.first;
        };

        std::sort(arc_array.begin(), arc_array.end(), arc_sort);

        for (int begin = 0; begin < arc_array.size(); begin++)
        {
            std::vector<int> outfaceIndex;
            if (begin != 0)
            {
                std::vector<std::vector<std::vector<trimesh::vec2>>> polygon(1, std::vector<std::vector<trimesh::vec2>>(1, std::vector<trimesh::vec2>()));
                for (int vi = begin_and_end[arc_array[begin].second].first; vi < begin_and_end[arc_array[begin].second].second; vi++)
                {
                    polygon[0][0].push_back(trimesh::vec2(joinmesh.vertices[vi].p.x, joinmesh.vertices[vi].p.y));
                }
                std::vector<int> infaceIndex(joinmesh.faces.size());
                std::iota(infaceIndex.begin(), infaceIndex.end(), 0);
                topomesh::embedingAndCutting(&joinmesh, polygon[0], infaceIndex);
                std::vector<int> newinfaceIndex(joinmesh.faces.size());
                std::iota(newinfaceIndex.begin(), newinfaceIndex.end(), 0);
                topomesh::polygonInnerFaces(&joinmesh, polygon, newinfaceIndex, outfaceIndex);
            }
            
            if (outfaceIndex.empty())
            {
                std::vector<std::pair<trimesh::point, int>> lines;
                if (!is_reverse[arc_array[begin].second])
                    for (int vi = begin_and_end[arc_array[begin].second].first; vi < begin_and_end[arc_array[begin].second].second; vi++)
                    {
                        lines.push_back(std::make_pair(joinmesh.vertices[vi].p, vi));
                    }
                else
                    for (int vi = begin_and_end[arc_array[begin].second].second - 1; vi >= begin_and_end[arc_array[begin].second].first; vi--)
                    {
                        lines.push_back(std::make_pair(joinmesh.vertices[vi].p, vi));
                    }
                topomesh::EarClipping earclip(lines);
                std::vector<trimesh::ivec3> result = earclip.getResult();
                for (int fi = 0; fi < result.size(); fi++)
                {
                    joinmesh.appendFace(result[fi].x, result[fi].y, result[fi].z);
                }
            }
            else
            {
                for (int fi = 0; fi < outfaceIndex.size(); fi++)
                    joinmesh.deleteFace(outfaceIndex[fi]);
            }
        }

        trimesh::TriMesh* jointmesh = new trimesh::TriMesh();
        joinmesh.quickTransform(jointmesh);     
        vertexsize = newmesh->vertices.size();
        for (int vi = 0; vi < jointmesh->vertices.size(); vi++)
            newmesh->vertices.push_back(jointmesh->vertices[vi]);
        for (int fi = 0; fi < jointmesh->faces.size(); fi++)
            newmesh->faces.push_back(trimesh::TriMesh::Face(jointmesh->faces[fi][0] + vertexsize, jointmesh->faces[fi][2] + vertexsize, jointmesh->faces[fi][1] + vertexsize));
        dumplicateMesh(newmesh);
        
#else
        std::map<int, int> vmap;
        std::map<int, int> fmap;        
        topomesh::MMeshT mt(mesh, botfaces, vmap, fmap);
        std::vector<bool> deletefaces(mesh->faces.size(), false);
        for (int f = 0; f < botfaces.size(); f++)
        {
            deletefaces[botfaces[f]] = true;
        }
        trimesh::remove_faces(mesh, deletefaces);
        trimesh::remove_unused_vertices(mesh);
             
       
        std::vector<std::vector<trimesh::point>> sequentials = GetOpenMeshBoundarys(*newmesh);
        std::vector<std::vector<std::vector<trimesh::vec2>>> polygon(1, std::vector<std::vector<trimesh::vec2>>(sequentials.size(), std::vector<trimesh::vec2>()));
        for (int i = 0; i < sequentials.size(); i++)
        {
            for(int j=0;j<sequentials[i].size();j++)
                polygon[0][i].push_back(trimesh::vec2(sequentials[i][j].x, sequentials[i][j].y));
        }
        std::vector<int> faces(mt.faces.size());
        std::iota(faces.begin(), faces.end(), 0);
        topomesh::embedingAndCutting(&mt, polygon[0], faces);
        std::vector<int> newinfaceIndex(mt.faces.size());
        std::iota(newinfaceIndex.begin(), newinfaceIndex.end(), 0);
        std::vector<int> outfaceIndex;
        topomesh::polygonInnerFaces(&mt, polygon, newinfaceIndex, outfaceIndex);
        for (int i = 0; i < outfaceIndex.size(); i++)
            mt.deleteFace(outfaceIndex[i]);
        trimesh::TriMesh* resultmesh=new trimesh::TriMesh();
        mt.quickTransform(resultmesh);

        for (int fi = 0; fi < resultmesh->faces.size(); fi++)
        {
            trimesh::point v1 = resultmesh->vertices[resultmesh->faces[fi][1]] - resultmesh->vertices[resultmesh->faces[fi][0]];
            trimesh::point v2 = resultmesh->vertices[resultmesh->faces[fi][2]] - resultmesh->vertices[resultmesh->faces[fi][0]];
            float z = (v1 % v2).z;
            if (z > 0)
                resultmesh->faces[fi] = trimesh::TriMesh::Face(resultmesh->faces[fi][0], resultmesh->faces[fi][2], resultmesh->faces[fi][1]);
        }

        int facesize = newmesh->faces.size();
        int vertexsize = newmesh->vertices.size();
        for (int vi = 0; vi < resultmesh->vertices.size(); vi++)
            newmesh->vertices.push_back(resultmesh->vertices[vi]);
        for (int fi = 0; fi < resultmesh->faces.size(); fi++)
            newmesh->faces.push_back(trimesh::TriMesh::Face(resultmesh->faces[fi][0] + vertexsize, resultmesh->faces[fi][1] + vertexsize, resultmesh->faces[fi][2] + vertexsize));
        int resulrfacesize = newmesh->faces.size();

        vertexsize = newmesh->vertices.size();
        for (int vi = 0; vi < mesh->vertices.size(); vi++)
            newmesh->vertices.push_back(mesh->vertices[vi]);
        for (int fi = 0; fi < mesh->faces.size(); fi++)
            newmesh->faces.push_back(trimesh::TriMesh::Face(mesh->faces[fi][0] + vertexsize, mesh->faces[fi][1] + vertexsize, mesh->faces[fi][2] + vertexsize));
        //newmesh->write("step1.stl");
        //newmesh->vertices[9421] = trimesh::vec3(newmesh->vertices[9421].x, newmesh->vertices[9421].y+0.000001, newmesh->vertices[9421].z);
        dumplicateMesh2(newmesh,nullptr,0.3f,2e-3f);
        //dumplicateMesh(newmesh, nullptr, 0.3f, 1e-4f);
#endif
        newmesh->write("step2.stl");
      
        newmesh->need_across_edge();
        newmesh->clear_neighbors();
        newmesh->need_neighbors();
        newmesh->clear_adjacentfaces();
        newmesh->need_adjacentfaces();

        std::vector<bool> is_boundarys(newmesh->vertices.size(),false);      
        for (int vi = 0; vi < newmesh->vertices.size(); vi++)
        {
            if (newmesh->adjacentfaces[vi].size() == newmesh->neighbors[vi].size()-1)
                is_boundarys[vi] = true;          
        }

        std::vector<bool> is_vis(newmesh->vertices.size(), false);      

        std::vector<bool> deleteface1(newmesh->faces.size(), false);
        for (int fi = 0; fi < facesize; fi++)
        {
            bool is_boundary = false;
            int fii = 0;
            for (; fii < newmesh->across_edge[fi].size(); fii++)
            {
                if (newmesh->across_edge[fi][fii] == -1)
                {
                    is_boundary = true;
                    break;
                }
            }
            if (is_boundary)
            {
                deleteface1[fi] = true;                              
                int oppovertex = fii;
                trimesh::point fn = (newmesh->vertices[newmesh->faces[fi][1]] - newmesh->vertices[newmesh->faces[fi][0]]) %
                    (newmesh->vertices[newmesh->faces[fi][2]] - newmesh->vertices[newmesh->faces[fi][0]]);
                trimesh::normalize(fn);
                int beginindex = newmesh->faces[fi][(oppovertex+1)%3];
                if (fi == 370)
                    std::cout << "\n";
                std::vector<int> vertexlines;
                int vsize = vertexlines.size();
                for (int vi = 0; vi < newmesh->neighbors[beginindex].size(); vi++)
                {
                    int index = newmesh->neighbors[beginindex][vi];
                    if (index == newmesh->faces[fi][(oppovertex + 2) % 3]) continue;
                    if (is_boundarys[index] && !is_vis[index])
                    {
                        trimesh::point tn = (newmesh->vertices[index] - newmesh->vertices[beginindex]) %
                            (newmesh->vertices[newmesh->faces[fi][oppovertex]]-newmesh->vertices[beginindex]);
                        trimesh::normalize(tn);
                        if (std::fabs(tn.x - fn.x) <= 1e-4f && std::fabs(tn.y - fn.y) <= 1e-4f && std::fabs(tn.z - fn.z) <= 1e-4f)
                        {
                            /*if(newmesh->neighbors[beginindex][vi]==9618)
                                std::cout << "\n";*/
                            vertexlines.push_back(newmesh->neighbors[beginindex][vi]); 
                            is_vis[newmesh->neighbors[beginindex][vi]] = true;
                            break;
                        }
                    }
                }
                while (vsize != vertexlines.size())
                {
                    vsize = vertexlines.size();
                    int index = vertexlines.back();
                    bool is_break = false;
                    if (index == newmesh->faces[fi][(oppovertex + 2) % 3]) break;
                    for (int vi = 0; vi < newmesh->neighbors[index].size(); vi++)
                    {
                        if (newmesh->neighbors[index][vi] == newmesh->faces[fi][(oppovertex + 2) % 3])
                        {
                            is_break = true;
                            break;
                        }
                    }
                    if (is_break)
                        break;
                    for (int vi = 0; vi < newmesh->neighbors[index].size(); vi++)
                    {                      
                        if (is_boundarys[newmesh->neighbors[index][vi]] && !is_vis[newmesh->neighbors[index][vi]])
                        {
                           /* if (newmesh->neighbors[index][vi] == 9618)
                                std::cout << "\n";*/
                            trimesh::point tn = (newmesh->vertices[newmesh->neighbors[index][vi]] - newmesh->vertices[beginindex]) %
                                (newmesh->vertices[newmesh->faces[fi][oppovertex]] - newmesh->vertices[beginindex]);
                            trimesh::normalize(tn);
                            if (std::fabs(tn.x - fn.x) <= 1e-4f && std::fabs(tn.y - fn.y) <= 1e-4f && std::fabs(tn.z - fn.z) <= 1e-4f)
                            {
                                vertexlines.push_back(newmesh->neighbors[index][vi]);
                                is_vis[newmesh->neighbors[index][vi]] = true;
                                break;
                            }
                        }
                    }
                }
                if (vertexlines.empty()) 
                    continue;

                newmesh->faces.push_back(trimesh::ivec3(newmesh->faces[fi][oppovertex], newmesh->faces[fi][(oppovertex+1)%3], vertexlines[0]));
                for (int ii = 0; ii < vertexlines.size() - 1; ii++)
                {
                    newmesh->faces.push_back(trimesh::ivec3(newmesh->faces[fi][oppovertex], vertexlines[ii], vertexlines[ii+1]));
                }
                newmesh->faces.push_back(trimesh::ivec3(newmesh->faces[fi][oppovertex], vertexlines.back(),newmesh->faces[fi][(oppovertex + 2) % 3]));
            }
        }
        int beginsize = deleteface1.size();
        for (int fi = beginsize; fi < newmesh->faces.size(); fi++)
            deleteface1.push_back(false);

        trimesh::remove_faces(newmesh, deleteface1);
        trimesh::remove_unused_vertices(newmesh);
        
       
       /* std::vector<std::vector<int>> sequentialsindex;
        getTriMeshBoundarys(*newmesh,sequentialsindex);
        for (int i = 0; i < sequentialsindex.size(); i++)
        {
            std::vector<std::pair<trimesh::point, int>> lines;
            for (int j = 0; j < sequentialsindex[i].size(); j++)
                lines.push_back(std::make_pair(newmesh->vertices[sequentialsindex[i][j]], sequentialsindex[i][j]));
            topomesh::EarClipping earclip(lines);
            std::vector<trimesh::ivec3> result = earclip.getResult();
            for (int fi = 0; fi < result.size(); fi++)
            {         
                newmesh->faces.push_back(result[fi]);
            }
        }*/           

        newmesh->write("newmesh.stl");
    }


    void getTriMeshBoundarys(trimesh::TriMesh& trimesh, std::vector<std::vector<int>>& sequentials)
    {
        CMesh mesh(&trimesh);
        std::vector<int> edges;
        mesh.SelectIndividualEdges(edges);             
        mesh.GetSequentialPoints(edges, sequentials);
    }

    TriPolygons GetOpenMeshBoundarys(const trimesh::TriMesh& triMesh, HoneyCombDebugger* debugger)
    {
        CMesh mesh(&triMesh);
        TriPolygons polys;
        std::vector<int> edges;
        mesh.SelectIndividualEdges(edges);
        //CMesh edgeMesh = mesh.SaveEdgesToMesh(edges);
        //edgeMesh.WriteSTLFile("edges.stl");
        //第0步，底面轮廓边界所有边排序
        std::vector<std::vector<int>>sequentials;
        mesh.GetSequentialPoints(edges, sequentials);
        //第1步，构建底面边界轮廓多边形
        polys.reserve(sequentials.size());
        const auto& points = mesh.mpoints;
        for (const auto& seq : sequentials) {
            TriPolygon border;
            border.reserve(seq.size());
            for (const auto& v : seq) {
                const auto& p = points[v];
                border.emplace_back(p);
            }
            polys.emplace_back(border);
        }
        if (debugger) {
            //显示底面边界轮廓多边形
            debugger->onGenerateInfillPolygons(polys);
        }
        return polys;
    }

    HexaPolygons GenerateTriPolygonsHexagons(const TriPolygons& polys, const HexagonArrayParam& hexagonparams)
    {

        return HexaPolygons();
    }

    HexaPolygons GenerateHexagonsGridArray(const HexagonArrayParam& hexagonparams)
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
                    const auto& center = trimesh::vec3(pt.x, pt.y, 0);
                    topomesh::Hexagon hexagon(center, side);
                    const auto& border = hexagon.border;
                    HexaPolygon hexa;
                    hexa.center = hexagon.centroid;
                    hexa.poly.reserve(border.size());
                    for (const auto& p : border) {
                        hexa.poly.emplace_back(trimesh::vec3((float)p.x, (float)p.y, p0.z));
                    }
                    hexa.coord = trimesh::ivec3(j, -(i + j) / 2, (i - j) / 2);
                    polygons.polys.emplace_back(std::move(hexa));
                } else {
                    const auto& center = trimesh::vec3(pt.x, (double)pt.y - ydelta, 0);
                    topomesh::Hexagon hexagon(center, side);
                    const auto& border = hexagon.border;
                    HexaPolygon hexa;
                    hexa.center = hexagon.centroid;
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
        const int nums = hexas.polys.size();
        for (int i = 0; i < nums; ++i) {
            if (hexas.polys[i].standard) {
                for (int j = 0; j < 6; ++j) {
                    hexas.polys[i].h2pPointMap.emplace(j, j);
                    hexas.polys[i].p2hPointMap.emplace(j, j);
                    hexas.polys[i].h2pEdgeMap.emplace(j, j);
                    hexas.polys[i].p2hEdgeMap.emplace(j, j);
                }
            }
            if (hexas.polys[i].edges.size() == 0) hexas.polys[i].edges.resize(6);
        }
        std::vector<std::vector<int>> neighborhoods(nums, std::vector<int>(nums, 1));
        std::vector<std::vector<int>> associatehoods(nums, std::vector<int>(nums, 1));
        for (int i = 0; i < nums; ++i) {
            for (int j = 0; j < nums; ++j) {
                //更新完整边对应关系
                if (neighborhoods[i][j] && (j != i)) {
                    const auto& hexa = hexas.polys[i];
                    const auto& hexb = hexas.polys[j];
                    int res = int(hex_neighbor(hexa.coord, hexb.coord));
                    if (res >= 0) {
                        const auto& h2pmap1 = hexa.h2pEdgeMap;
                        auto itr1 = h2pmap1.find(res); ///<当前线段对应六角网格边的索引
                        if (itr1 != h2pmap1.end()) {
                            const auto& h2pmap2 = hexb.h2pEdgeMap;
                            const int inx = (res + 3) % 6; ///<临近线段对应六角网格边的索引
                            const auto& itr2 = h2pmap2.find(inx);
                            if (itr2 != h2pmap2.end()) {
                                hexas.polys[i].edges[itr1->second].neighbor = j;
                                hexas.polys[j].edges[itr2->second].neighbor = i;
                                neighborhoods[i][j] = 0;
                                neighborhoods[j][i] = 0;
                            }
                        }
                    }
                }
            }
        }
        for (int i = 0; i < nums; ++i) {
            for (int j = 0; j < nums; ++j) {
                //更新裁剪边对应关系
                if (associatehoods[i][j] && (neighborhoods[i][j]) && (j != i)) {
                    const auto& hexa = hexas.polys[i];
                    const auto& hexb = hexas.polys[j];
                    int res = int(hex_neighbor(hexa.coord, hexb.coord));
                    if (res >= 0) {
                        const auto& h2pmap1 = hexa.h2pPointMap;
                        auto itr1 = h2pmap1.find(res); ///<当前线段起始点为六角网格网格顶点
                        if (itr1 != h2pmap1.end()) {
                            const auto& h2pmap2 = hexb.h2pPointMap;
                            const int inx = (res + 3 + 1) % 6; ///<关联线段终点为六角网格网格顶点
                            const auto& itr2 = h2pmap2.find(inx);
                            if (itr2 != h2pmap2.end()) {
                                const int sizeb = hexb.poly.size();
                                hexas.polys[i].edges[itr1->second].associate = j;
                                hexas.polys[j].edges[(itr2->second + sizeb - 1) % sizeb].associate = i;
                                associatehoods[i][j] = 0;
                                associatehoods[j][i] = 0;
                            }
                        }
                    }
                }
            }
        }
        //更新内外轮廓之间的网格相邻关系
        for (int i = 0; i < nums; ++i) {
            for (int j = 0; j < nums; ++j) {
                if (associatehoods[i][j] && (neighborhoods[i][j]) && (j != i)) {
                    const auto& hexa = hexas.polys[i];
                    const auto& hexb = hexas.polys[j];
                    int res = int(hex_neighbor(hexa.coord, hexb.coord));
                    if (res >= 0) {
                        const auto& h2pmap1 = hexa.h2pPointMap;
                        auto itr1 = h2pmap1.find(res); ///<当前线段起始点为六角网格网格顶点
                        if (itr1 != h2pmap1.end()) {
                            const auto& h2pmap2 = hexb.h2pPointMap;
                            const int inx = (res + 3) % 6; ///<关联线段始点为六角网格网格顶点
                            const auto & itr2 = h2pmap2.find(inx);
                            if (itr2 != h2pmap2.end()) {
                                const int sizeb = hexb.poly.size();
                                hexas.polys[i].edges[itr1->second].associate = j;
                                hexas.polys[j].edges[itr2->second].associate = i;
                                associatehoods[i][j] = 0;
                                associatehoods[j][i] = 0;
                            }
                        }
                    }
                    if (res >= 0) {
                        const auto& h2pmap1 = hexa.h2pPointMap;
                        auto itr1 = h2pmap1.find((res + 1) % 6); ///<当前线段终点点为六角网格网格顶点
                        if (itr1 != h2pmap1.end()) {
                            const auto& h2pmap2 = hexb.h2pPointMap;
                            const int inx = (res + 3 + 1) % 6; ///<关联线终点为六角网格网格顶点
                            const auto& itr2 = h2pmap2.find(inx);
                            if (itr2 != h2pmap2.end()) {
                                const int sizea = hexa.poly.size();
                                const int sizeb = hexb.poly.size();
                                hexas.polys[i].edges[(itr1->second + sizea - 1) % sizea].associate = j;
                                hexas.polys[j].edges[(itr2->second + sizeb - 1) % sizeb].associate = i;
                                associatehoods[i][j] = 0;
                                associatehoods[j][i] = 0;
                            }
                        }
                    }
                }
            }
        }
        for (int i = 0; i < nums; ++i) {
            const int size = hexas.polys[i].poly.size();
            for (int j = 0; j < size; ++j) {
                const int associate = hexas.polys[i].edges[j].associate;
                const int neighbor = hexas.polys[i].edges[j].neighbor;
                hexas.polys[i].edges[j].relate = associate < 0 ? neighbor : associate;
            }
        }
        const float cradius = hexas.side * param.ratio * 0.5f;
        const float cdelta = param.delta + 2 * cradius;
        const float cheight = param.height + cradius;
        float cmax = param.height + 2 * cradius;
        //更新六棱柱每个可加孔洞侧面最低点最高点坐标值
        for (auto& hexa : hexas.polys) {
            const int polysize = hexa.poly.size();
            for (int j = 0; j < polysize; ++j) {
                const int res = hexa.p2hEdgeMap[j];
                if (hexa.edges[j].neighbor >= 0 && (!hexa.edges[j].canAdd) && (res >= 0)) {
                    auto& edge = hexa.edges[j];
                    auto& next = hexa.edges[(j + 1) % polysize];
                    auto& oh = hexas.polys[hexa.edges[j].neighbor];
                    const auto& h2pmap = oh.h2pEdgeMap;
                    const int inx = (res + 3) % 6;
                    auto itr = h2pmap.find(inx);
                    if (itr != h2pmap.end()) {
                        const int ind = itr->second;
                        auto& oe = oh.edges[ind];
                        const int osize = oh.poly.size();
                        auto& onext = oh.edges[(ind + 1) % osize];
                        edge.lowHeight = std::max(edge.lowHeight, next.lowHeight);
                        edge.topHeight = std::min(edge.topHeight, next.topHeight);
                        oe.lowHeight = std::max(oe.lowHeight, onext.lowHeight);
                        oe.topHeight = std::min(oe.topHeight, onext.topHeight);
                        if (param.holeConnect) {
                            float lowerlimit = std::max(edge.lowHeight, oe.lowHeight);
                            float upperlimit = std::min(edge.topHeight, oe.topHeight);
                            float medium = (lowerlimit + upperlimit) * 0.5;
                            float scope = (upperlimit - lowerlimit) * 0.5;
                            if (medium < cheight && (cmax - medium) < scope) {
                                hexa.edges[j].canAdd = true;
                                oh.edges[ind].canAdd = true;
                            } else if (medium >= cheight) {
                                float appro = std::round((medium - cheight) / cdelta);
                                float dist = std::abs((medium - cheight) - cdelta * appro);
                                if (dist + cradius < scope) {
                                    hexa.edges[j].canAdd = true;
                                    oh.edges[ind].canAdd = true;
                                }
                            }
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
                if (edges[i].relate >= 0)
                    polys.emplace_back(hexas.polys[edges[i].relate].poly);
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
                int val = edges[dir].relate;
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
        int hexagonsize = 0, columnvertexs = 0;
        for (auto& hexa : hexas.polys) {
            hexa.startIndex = columnvertexs;
            const auto& poly = hexa.poly;
            for (int i = 0; i < poly.size(); ++i) {
                const auto& p = poly[i];
                points.emplace_back(trimesh::vec3(p.x, p.y, hexa.edges[i].lowHeight));
                ++columnvertexs, ++hexagonsize;
            }
        }

        for (const auto& hexa : hexas.polys) {
            const auto& poly = hexa.poly;
            for (int i = 0; i < poly.size(); ++i) {
                const auto& p = poly[i];
                points.emplace_back(trimesh::vec3(p.x, p.y, hexa.edges[i].topHeight));
                ++columnvertexs;
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
            ///底部网格中间矩形区域
            for (auto& hexa : hexas.polys) {
                const int size = hexa.poly.size();
                //底部完整网格中间矩形区域
                for (int i = 0; i < size; ++i) {
                    const int res = hexa.p2hEdgeMap[i]; ///<当前线段对应六角网格边长
                    if ((!hexa.edges[i].addRect) && (hexa.edges[i].neighbor >= 0) && (res >= 0)) {
                        auto& oh = hexas.polys[hexa.edges[i].neighbor];
                        const auto& h2pmap = oh.h2pEdgeMap;
                        const int inx = (res + 3) % 6; ///<临近线段对应六角网格的边长
                        auto itr = h2pmap.find(inx);
                        if (itr != h2pmap.end()) {
                            const int& start = hexa.startIndex;
                            const int& a = start + i;
                            const int& b = start + (i + 1) % size;

                            const int& ostart = oh.startIndex;
                            const int osize = oh.poly.size();
                            const int ind = itr->second;
                            const int& c = ostart + ind;
                            const int& d = ostart + (ind + 1) % osize;
                            faces.emplace_back(trimesh::ivec3(b, c, d));
                            faces.emplace_back(trimesh::ivec3(b, d, a));
                            hexa.edges[i].addRect = true;
                            oh.edges[ind].addRect = true;
                            bottomfacenums += 2;
                        }
                    }
                }
                //底部边界裁剪网格中间矩形区域
                for (int i = 0; i < size; ++i) {
                    const int res = hexa.p2hPointMap[i]; ///<当前线段起始点对应六角网格顶点
                    if ((!hexa.edges[i].addRect) && (hexa.edges[i].associate >= 0) && (res >= 0)) {
                        auto& oh = hexas.polys[hexa.edges[i].associate];
                        const auto& h2pmap = oh.h2pPointMap;
                        const int inx = (res + 3 + 1) % 6; ///<关联线段终点对应六角网格的边顶点
                        auto itr = h2pmap.find(inx);
                        if (itr != h2pmap.end()) {
                            const int& start = hexa.startIndex;
                            const int& a = start + i;
                            const int& b = start + (i + 1) % size;

                            const int& ostart = oh.startIndex;
                            const int osize = oh.poly.size();
                            const int ind = itr->second;
                            const int& d = ostart + ind;
                            const int& c = ostart + (ind + osize - 1) % osize;
                            faces.emplace_back(trimesh::ivec3(b, c, d));
                            faces.emplace_back(trimesh::ivec3(b, d, a));
                            hexa.edges[i].addRect = true;
                            oh.edges[(ind + osize - 1) % osize].addRect = true;
                            bottomfacenums += 2;
                        }
                    }
                }
            }
            //内外轮廓之间的两个裁剪网格
            for (auto& hexa : hexas.polys) {
                const int size = hexa.poly.size();
                for (int i = 0; i < size; ++i) {
                    const int res = hexa.p2hPointMap[i]; ///<当前线段始点对应六角网格顶点
                    if ((!hexa.edges[i].addRect) && (hexa.edges[i].associate >= 0) && (res >= 0)) {
                        auto& oh = hexas.polys[hexa.edges[i].associate];
                        const auto& h2pmap = oh.h2pPointMap;
                        const int inx = (res + 3) % 6; ///<始点
                        auto itr = h2pmap.find(inx);
                        if (itr != h2pmap.end()) {
                            const int& start = hexa.startIndex;
                            const int& a = start + i;
                            const int& b = start + (i + 1) % size;

                            const int& ostart = oh.startIndex;
                            const int osize = oh.poly.size();
                            const int ind = itr->second;
                            const int& c = ostart + ind;
                            const int& d = ostart + (ind + 1) % osize;
                            faces.emplace_back(trimesh::ivec3(b, c, d));
                            faces.emplace_back(trimesh::ivec3(b, d, a));
                            hexa.edges[i].addRect = true;
                            oh.edges[ind].addRect = true;
                            bottomfacenums += 2;
                        }
                    }
                }
                for (int i = 0; i < size; ++i) {
                    const int res = hexa.p2hPointMap[(i + 1) % size]; ///<当前线段终点对应六角网格顶点
                    if ((!hexa.edges[i].addRect) && (hexa.edges[i].associate >= 0) && (res >= 0)) {
                        auto& oh = hexas.polys[hexa.edges[i].associate];
                        const auto& h2pmap = oh.h2pPointMap;
                        const int inx = (res + 3) % 6; ///<终点
                        auto itr = h2pmap.find(inx);
                        if (itr != h2pmap.end()) {
                            const int& start = hexa.startIndex;
                            const int& a = start + i;
                            const int& b = start + (i + 1) % size;

                            const int& ostart = oh.startIndex;
                            const int osize = oh.poly.size();
                            const int ind = itr->second;
                            const int& c = ostart + (ind + osize - 1) % osize;
                            const int& d = ostart + ind;
                            faces.emplace_back(trimesh::ivec3(b, c, d));
                            faces.emplace_back(trimesh::ivec3(b, d, a));
                            hexa.edges[i].addRect = true;
                            oh.edges[(ind + osize - 1) % osize].addRect = true;
                            bottomfacenums += 2;
                        }
                    }
                }
            }
            ///底部三个网格中间部分
            for (auto& hexa : hexas.polys) {
                const int size = hexa.poly.size();
                for (int i = 0; i < size; ++i) {
                    const int nb = hexa.edges[i].relate;
                    const int nc = hexa.edges[(i + 1) % size].relate;
                    const int res = hexa.p2hPointMap[(i + 1) % size];
                    if ((!hexa.edges[(i + 1) % size].addTriangle) && (nb >= 0) && (nc >= 0) && (res >= 0)) {
                        auto& ohb = hexas.polys[nb];
                        auto& ohc = hexas.polys[nc];
                        const int& a = hexa.startIndex + (i + 1) % size; ///<线段终点

                        const auto& h2pmapb = ohb.h2pPointMap;
                        const int inxb = (res + 2) % 6; ///<线段始点
                        auto itrb = h2pmapb.find(inxb);
                        if (itrb != h2pmapb.end()) {
                            const int indb = itrb->second;
                            const int& b = ohb.startIndex + indb; 

                            const auto& h2pmapc = ohc.h2pPointMap;
                            const int inxc = (res + 3 + 1) % 6; ///<线段终点
                            auto itrc = h2pmapc.find(inxc);
                            if (itrc != h2pmapc.end()) {
                                const int indc = itrc->second;
                                const int sizec = ohc.poly.size();
                                const int& c = ohc.startIndex + indc;
                                faces.emplace_back(trimesh::ivec3(a, c, b));
                                hexa.edges[(i + 1) % size].addTriangle = true;
                                ohb.edges[indb].addTriangle = true;
                                ohc.edges[indc].addTriangle = true;
                                ++bottomfacenums;
                            }
                        }
                    }
                }
            }
            for (auto& hexa : hexas.polys) {
                const int size = hexa.poly.size();
                for (int i = 0; i < size; ++i) {
                    const int nb = hexa.edges[i].relate;
                    const int nc = hexa.edges[(i + 1) % size].relate;
                    const int res = hexa.p2hPointMap[(i + 1) % size];
                    if ((!hexa.edges[(i + 1) % size].addTriangle) && (nb >= 0) && (nc >= 0) && (res >= 0)) {
                        auto& ohb = hexas.polys[nb];
                        auto& ohc = hexas.polys[nc];
                        const int& a = hexa.startIndex + (i + 1) % size; ///<六角网格顶点

                        const auto& h2pmapb = ohb.h2pPointMap;
                        const int inxb = (res + 3) % 6; ///<下一个是六角网格顶点
                        auto itrb = h2pmapb.find(inxb);
                        if (itrb != h2pmapb.end()) {
                            const int sizeb = ohb.poly.size();
                            const int indb = (itrb->second + sizeb - 1) % sizeb;
                            const int& b = ohb.startIndex + indb;

                            const auto& h2pmapc = ohc.h2pPointMap;
                            const int inxc = (res + 3) % 6; ///<上一个是六角网格顶点
                            auto itrc = h2pmapc.find(inxc);
                            if (itrc != h2pmapc.end()) {
                                const int sizec = ohc.poly.size();
                                const int indc = (itrc->second + 1) % sizec;
                                const int& c = ohc.startIndex + indc;
                                faces.emplace_back(trimesh::ivec3(a, c, b));
                                hexa.edges[(i + 1) % size].addTriangle = true;
                                ohb.edges[indb].addTriangle = true;
                                ohc.edges[indc].addTriangle = true;
                                ++bottomfacenums;
                            }
                        }
                    }
                }
            }
        }
        ///添加连接孔洞的顶点坐标
        for (auto& hexa : hexas.polys) {
            const int size = hexa.poly.size();
            for (int i = 0; i < size; ++i) {
                const int res = hexa.p2hEdgeMap[i];
                if (hexa.edges[i].canAdd && (!hexa.edges[i].hasAdd) && (hexa.edges[i].neighbor >= 0) && (res >= 0)) {
                    auto& edge = hexa.edges[i];
                    auto& next = hexa.edges[(i + 1) % size];
                    auto& oh = hexas.polys[hexa.edges[i].neighbor];
                    const int inx = (res + 3) % 6;
                    const auto& h2pmap = oh.h2pEdgeMap;
                    auto itr = h2pmap.find(inx);
                    if (itr != h2pmap.end()) {
                        const int ind = itr->second;
                        const int osize = oh.poly.size();
                        auto& oe = oh.edges[ind];
                        auto& onext = oh.edges[(ind + 1) % osize];
                        float lowerlimit = std::max(edge.lowHeight, oe.lowHeight);
                        float upperlimit = std::min(edge.topHeight, oe.topHeight);
                        float medium = (lowerlimit + upperlimit) * 0.5;
                        float scope = (upperlimit - lowerlimit) * 0.5;
                        const auto& a = hexa.poly[i];
                        const auto& b = hexa.poly[(i + 1) % size];
                        const auto& c = oh.poly[ind];
                        const auto& d = oh.poly[(ind + 1) % osize];
                        if (medium < cheight && (cmax - medium) < scope) {
                            ///此时最高点下方有几个圆就可增加几个圆
                            int n = std::floor((upperlimit - cmax) / cdelta);
                            edge.holenums = n + 1, oe.holenums = n + 1;
                            edge.starts.reserve((size_t)n + 1);
                            oe.starts.reserve((size_t)n + 1);
                            const auto& ab = (a + b) * 0.5;
                            const auto& c1 = trimesh::vec3(ab.x, ab.y, cheight);
                            std::vector<int> corner1;
                            auto&& poly1 = traitPlanarCircle(c1, cradius, corner1, a - b, nslices);
                            edge.corners.swap(corner1);
                            edge.starts.emplace_back(points.size());
                            holeStarts.emplace_back(points.size());
                            points.insert(points.end(), poly1.begin(), poly1.end());

                            const auto& cd = (c + d) * 0.5;
                            const auto& c2 = trimesh::vec3(cd.x, cd.y, cheight);
                            std::vector<int> corner2;
                            auto&& poly2 = traitPlanarCircle(c2, cradius, corner2, c - d, nslices);
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
                                edge.pointIndexs.reserve(n);
                                edge.holeIndexs.reserve(n);
                                oe.pointIndexs.reserve(n);
                                oe.holeIndexs.reserve(n);
                            }
                            for (int j = 1; j <= n; ++j) {
                                const auto& z = cheight + float(2 * j - 1) * addz;
                                edge.pointIndexs.emplace_back(points.size());
                                edge.holeIndexs.emplace_back(j);
                                points.emplace_back(trimesh::vec3(a.x, a.y, z));
                                points.emplace_back(trimesh::vec3(b.x, b.y, z));
                                oe.pointIndexs.emplace_back(points.size());
                                oe.holeIndexs.emplace_back(j);
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
                            oh.edges[ind].hasAdd = true;
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
                                const auto& ab = (a + b) * 0.5;
                                const auto& c1 = trimesh::vec3(ab.x, ab.y, cheight + n1 * cdelta);
                                std::vector<int> corner1;
                                auto&& poly1 = traitPlanarCircle(c1, cradius, corner1, a - b, nslices);
                                edge.corners.swap(corner1);
                                edge.starts.emplace_back(points.size());
                                holeStarts.emplace_back(points.size());
                                points.insert(points.end(), poly1.begin(), poly1.end());

                                const auto& cd = (c + d) * 0.5;
                                const auto& c2 = trimesh::vec3(cd.x, cd.y, cheight + n1 * cdelta);
                                std::vector<int> corner2;
                                auto&& poly2 = traitPlanarCircle(c2, cradius, corner2, c - d, nslices);
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
                                    edge.pointIndexs.reserve(size_t(n) - 1);
                                    oe.pointIndexs.reserve(size_t(n) - 1);
                                    edge.holeIndexs.reserve(size_t(n) - 1);
                                    oe.holeIndexs.reserve(size_t(n) - 1);
                                }
                                edge.lowerHoles = n1;
                                oe.lowerHoles = n1;
                                for (int j = n1; j < n2 - 1; ++j) {
                                    const auto& z = cheight + float(2 * j + 1) * addz;
                                    edge.pointIndexs.emplace_back(points.size());
                                    edge.holeIndexs.emplace_back(j + 1);
                                    points.emplace_back(trimesh::vec3(a.x, a.y, z));
                                    points.emplace_back(trimesh::vec3(b.x, b.y, z));
                                    oe.pointIndexs.emplace_back(points.size());
                                    oe.holeIndexs.emplace_back(j + 1);
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
                                oh.edges[ind].hasAdd = true;
                            }
                        }
                    }
                }
            }
        }
        int noHoleEdges = hexagonsize - singleHoleEdges - mutiHoleEdges;
        int holeFaces = holesnum * singlenums;
        int rectfacenums = noHoleEdges * 2;
        int upperfacenums = hexagonsize;
        int allfacenums = holeFaces + rectfacenums + upperfacenums + bottomfacenums;
        faces.reserve(allfacenums);
        std::vector<int> topfaces;
        if (hexas.bSewTop) {
            for (int i = 0; i < hexas.polys.size(); ++i) {
                auto&& hexa = hexas.polys[i];
                const auto& poly = hexa.poly;
                int polynums = poly.size();
                const int& start = hexa.startIndex;
                int upstart = start + hexagonsize;
                ///六角网格顶部
                if (hexa.standard) {
                    for (int j = 1; j < polynums - 1; ++j) {
                        const int& cur = upstart + j;
                        const int& next = upstart + j + 1;
                        faces.emplace_back(trimesh::ivec3(upstart, next, cur));
                        topfaces.push_back(faces.size() - 1);
                    }
                } else {
                    int near = 0, far = 0;
                    const auto& center = hexa.center;
                    float mindist = trimesh::len2(poly[0] - center);
                    for (int j = 1; j < polynums; ++j) {
                        float d = trimesh::len2(poly[j] - hexa.center);
                        if (d < mindist) {
                            mindist = d;
                            near = j;
                        }
                    }
                    const trimesh::vec3& nearPt = poly[near];
                    const auto& dir = nearPt - center;
                    float maxdist = 0;
                    if (hexa.h2pPointMap.size() > 0) {
                        for (int j = 0; j < polynums; ++j) {
                            if (hexa.p2hPointMap[j] >= 0) {
                                float d = std::fabs((poly[j] - center) DOT dir);
                                if (d > maxdist) {
                                    maxdist = d;
                                    far = j;
                                }
                            }
                        }
                    } else {
                        for (int j = 0; j < polynums; ++j) {
                            float d = std::fabs((poly[j] - center) DOT dir);
                            if (d > maxdist) {
                                maxdist = d;
                                far = j;
                            }
                        }
                    }
                    int newupstart = start + hexagonsize + far;
                    for (int j = far + 1; j < polynums + far - 1; ++j) {
                        const int& cur = upstart + j % polynums;
                        const int& next = upstart + (j + 1) % polynums;
                        faces.emplace_back(trimesh::ivec3(newupstart, next, cur));
                        topfaces.push_back(faces.size() - 1);
                    }
                }
            }
        }
        hexas.topfaces.swap(topfaces);
        for (int i = 0; i < hexas.polys.size(); ++i) {
            const auto& hexa = hexas.polys[i];
            const auto& poly = hexa.poly;
            const int polysize = poly.size();
            const int& start = hexa.startIndex;
            //六角网格侧面
            for (int j = 0; j < polysize; ++j) {
                ///侧面下方两个顶点索引
                int a = start + j;
                int b = start + (j + 1) % polysize;
                ///侧面上方两个顶点索引
                int c = a + hexagonsize;
                int d = b + hexagonsize;
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
                    const auto& ledge = hexa.edges[(j + polysize - 1) % polysize];
                    const auto& nedge = hexa.edges[(j + 1) % polysize];
                    const auto& cIndexs = edge.holeIndexs;
                    const auto& lIndexs = ledge.holeIndexs;
                    const auto& nIndexs = nedge.holeIndexs;

                    if (edge.bmutihole) {
                        //第一个圆孔的第一个顶点的索引
                        std::vector<int> ldiffs, ndiffs;
                        std::set_difference(lIndexs.begin(), lIndexs.end(), cIndexs.begin(), cIndexs.end(), std::back_inserter(ldiffs));
                        std::set_difference(nIndexs.begin(), nIndexs.end(), cIndexs.begin(), cIndexs.end(), std::back_inserter(ndiffs));
                        int cstart = edge.starts.front();
                        c = edge.pointIndexs.front(); d = c + 1;
                        faces.emplace_back(trimesh::ivec3(cstart + r, c, d));
                        for (int k = 0; k < s; ++k) {
                            faces.emplace_back(trimesh::ivec3(cstart + k, d, cstart + k + 1));
                        }
                        //处理顶端b附近区域
                        if (!ndiffs.empty()) {
                            auto bf = std::find(ndiffs.begin(), ndiffs.end(), cIndexs.front() - 1);
                            if (bf != ndiffs.end()) {
                                int n0 = nedge.lowerHoles + 1;
                                int pos = std::distance(ndiffs.begin(), bf);
                                const auto& npoints = nedge.pointIndexs;
                                const int& b0 = npoints[(size_t)ndiffs[pos] - n0];
                                faces.emplace_back(trimesh::ivec3(b0, cstart + s, d));
                                for (size_t k = 0; k < pos; ++k) {
                                    faces.emplace_back(trimesh::ivec3(npoints[(size_t)ndiffs[k] - n0], cstart + s, npoints[(size_t)ndiffs[k + 1] - n0]));
                                }
                                faces.emplace_back(trimesh::ivec3(b, cstart + s, npoints.front()));
                            } else {
                                faces.emplace_back(trimesh::ivec3(b, cstart + s, d));
                            }
                        } else {
                            faces.emplace_back(trimesh::ivec3(b, cstart + s, d));
                        }
                        for (int k = s; k < t; ++k) {
                            faces.emplace_back(trimesh::ivec3(cstart + k, b, cstart + k + 1));
                        }
                        faces.emplace_back(trimesh::ivec3(b, a, cstart + t));
                        //处理顶端a附近区域
                        if (!ldiffs.empty()) {
                            auto af = std::find(ldiffs.begin(), ldiffs.end(), cIndexs.front() - 1);
                            if (af != ldiffs.end()) {
                                int n0 = ledge.lowerHoles + 1;
                                int pos = std::distance(ldiffs.begin(), af);
                                const auto& lpoints = ledge.pointIndexs;
                                faces.emplace_back(trimesh::ivec3(lpoints.front() + 1, cstart + t, a));
                                for (size_t k = 0; k < pos; ++k) {
                                    faces.emplace_back(trimesh::ivec3(lpoints[(size_t)ldiffs[k] - n0] + 1, cstart + t, lpoints[(size_t)ldiffs[k + 1] - n0]) + 1);
                                }
                                a = lpoints[(size_t)ldiffs[pos] - n0] + 1;
                            }
                        }
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
                            c = edge.pointIndexs[num]; d = c + 1;

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
                        c = start + j + hexagonsize;
                        d = start + (j + 1) % polysize + hexagonsize;
                        faces.emplace_back(trimesh::ivec3(cstart + r, c, d));
                        //处理顶端d附近区域
                        if (!ndiffs.empty()) {
                            auto df = std::find(ndiffs.begin(), ndiffs.end(), cIndexs.back() + 1);
                            if (df != ndiffs.end()) {
                                int n0 = nedge.lowerHoles + 1;
                                int pos = std::distance(ndiffs.begin(), df);
                                const auto& npoints = nedge.pointIndexs;
                                for (size_t k = pos; k < ndiffs.size() - 1; ++k) {
                                    faces.emplace_back(trimesh::ivec3(npoints[(size_t)ndiffs[k] - n0], cstart + r, npoints[(size_t)ndiffs[k + 1] - n0]));
                                }
                                faces.emplace_back(trimesh::ivec3(npoints[(size_t)ndiffs.back() - n0], cstart + r, d));
                                d = npoints[(size_t)ndiffs[pos] - n0];
                            }
                        }
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
                        //处理顶端c附近区域
                        if (!ldiffs.empty()) {
                            auto cf = std::find(ldiffs.begin(), ldiffs.end(), cIndexs.back() + 1);
                            if (cf != ldiffs.end()) {
                                int n0 = ledge.lowerHoles + 1;
                                int pos = std::distance(ldiffs.begin(), cf);
                                const auto& lpoints = ledge.pointIndexs;
                                faces.emplace_back(trimesh::ivec3(a, lpoints[(size_t)ldiffs[pos] - n0] + 1, cstart + u));
                                for (size_t k = pos; k < ldiffs.size() - 1; ++k) {
                                    faces.emplace_back(trimesh::ivec3(lpoints[(size_t)ldiffs[k] - n0] + 1, lpoints[(size_t)ldiffs[k + 1] - n0] + 1, cstart + r));
                                }
                                faces.emplace_back(trimesh::ivec3(lpoints.back() + 1, c, cstart + r));
                                c = lpoints[(size_t)ldiffs[pos] - n0] + 1;
                            } else {
                                faces.emplace_back(trimesh::ivec3(a, c, cstart + u));
                            }
                        } else {
                            faces.emplace_back(trimesh::ivec3(a, c, cstart + u));
                        }
                        for (int k = u; k < nslices; ++k) {
                            faces.emplace_back(trimesh::ivec3(c, cstart + (k + 1) % nslices, cstart + k));
                        }
                    } else {
                        ///只有单孔的侧面
                        int cstart = edge.starts.front();
                        faces.emplace_back(trimesh::ivec3(cstart + r, c, d));
                        //处理顶端d附近区域
                        if (!nIndexs.empty()) {
                            auto df = std::find(nIndexs.begin(), nIndexs.end(), 1);
                            if (df != nIndexs.end()) {
                                int pos = std::distance(nIndexs.begin(), df);
                                const auto& npoints = nedge.pointIndexs;
                                for (size_t k = pos; k < nIndexs.size() - 1; ++k) {
                                    faces.emplace_back(trimesh::ivec3(npoints[(size_t)nIndexs[k] - 1], cstart + r, npoints[(size_t)nIndexs[k + 1] - 1]));
                                }
                                faces.emplace_back(trimesh::ivec3(npoints[(size_t)nIndexs.back() - 1], cstart + r, d));
                                d = npoints[(size_t)nIndexs[pos] - 1];
                            }
                        }
                        for (int k = 0; k < s; ++k) {
                            faces.emplace_back(trimesh::ivec3(cstart + k, d, cstart + k + 1));
                        }
                        //b端不可能插入矩阵点，无需处理
                        faces.emplace_back(trimesh::ivec3(b, cstart + s, d));
                        for (int k = s; k < t; ++k) {
                            faces.emplace_back(trimesh::ivec3(cstart + k, b, cstart + k + 1));
                        }
                        faces.emplace_back(trimesh::ivec3(b, a, cstart + t));
                        //a端不可能插入矩阵点，无需处理
                        for (int k = t; k < u; ++k) {
                            faces.emplace_back(trimesh::ivec3(a, cstart + k + 1, cstart + k));
                        }
                        //处理顶端c附近区域
                        if (!lIndexs.empty()) {
                            auto cf = std::find(lIndexs.begin(), lIndexs.end(), 1);
                            if (cf != lIndexs.end()) {
                                int pos = std::distance(lIndexs.begin(), cf);
                                const auto& lpoints = ledge.pointIndexs;
                                faces.emplace_back(trimesh::ivec3(a, lpoints[(size_t)lIndexs[pos] - 1] + 1, cstart + u));
                                for (size_t k = pos; k < lIndexs.size() - 1; ++k) {
                                    faces.emplace_back(trimesh::ivec3(lpoints[(size_t)lIndexs[k] - 1] + 1, lpoints[(size_t)lIndexs[k + 1] - 1] + 1, cstart + r));
                                }
                                faces.emplace_back(trimesh::ivec3(lpoints.back() + 1, c, cstart + r));
                                c = lpoints[(size_t)lIndexs[pos] - 1] + 1;
                            } else {
                                faces.emplace_back(trimesh::ivec3(a, c, cstart + u));
                            }
                        } else {
                            faces.emplace_back(trimesh::ivec3(a, c, cstart + u));
                        }
                        for (int k = u; k < nslices; ++k) {
                            faces.emplace_back(trimesh::ivec3(c, cstart + (k + 1) % nslices, cstart + k));
                        }
                    }
                } else {
                    ///无孔的侧面
                    const auto& ledge = hexa.edges[(j + polysize - 1) % polysize];
                    const auto& nedge = hexa.edges[(j + 1) % polysize];
                    if (ledge.bmutihole && nedge.bmutihole) {
                        const auto& lIndexs = ledge.pointIndexs;
                        const auto& nIndexs = nedge.pointIndexs;
                        faces.emplace_back(trimesh::ivec3(b, a, lIndexs.front() + 1));
                        for (size_t k = 0; k < lIndexs.size() - 1; ++k) {
                            faces.emplace_back(trimesh::ivec3(b, lIndexs[k] + 1, lIndexs[k + 1] + 1));
                        }
                        faces.emplace_back(trimesh::ivec3(b, lIndexs.back() + 1, c));
                        faces.emplace_back(trimesh::ivec3(b, c, nIndexs.front()));
                        for (size_t k = 0; k < nIndexs.size() - 1; ++k) {
                            faces.emplace_back(trimesh::ivec3(nIndexs[k], c, nIndexs[k + 1]));
                        }
                        faces.emplace_back(trimesh::ivec3(nIndexs.back(), c, d));
                    } else if (ledge.bmutihole && (!nedge.bmutihole)) {
                        const auto& lIndexs = ledge.pointIndexs;
                        faces.emplace_back(trimesh::ivec3(b, a, lIndexs.front() + 1));
                        for (size_t k = 0; k < lIndexs.size() - 1; ++k) {
                            faces.emplace_back(trimesh::ivec3(b, lIndexs[k] + 1, lIndexs[k + 1] + 1));
                        }
                        faces.emplace_back(trimesh::ivec3(b, lIndexs.back() + 1, c));
                        faces.emplace_back(trimesh::ivec3(b, c, d));
                    } else if ((!ledge.bmutihole) && nedge.bmutihole) {
                        const auto& nIndexs = nedge.pointIndexs;
                        faces.emplace_back(trimesh::ivec3(b, a, c));
                        faces.emplace_back(trimesh::ivec3(b, c, nIndexs.front()));
                        for (size_t k = 0; k < nIndexs.size() - 1; ++k) {
                            faces.emplace_back(trimesh::ivec3(nIndexs[k], c, nIndexs[k + 1]));
                        }
                        faces.emplace_back(trimesh::ivec3(nIndexs.back(), c, d));
                    } else {
                        faces.emplace_back(trimesh::ivec3(b, a, c));
                        faces.emplace_back(trimesh::ivec3(b, c, d));
                    }
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
        dumplicateMesh(triMesh.get());
        return triMesh;
    }

    HexaPolygons generateEmbedHolesColumnar(trimesh::TriMesh* trimesh, const HoneyCombParam& honeyparams, ccglobal::Tracer* tracer, HoneyCombDebugger* debugger)
    {
        HexaPolygons hexas;
        //0.初始化Cmesh,并旋转
        std::shared_ptr<trimesh::TriMesh> result_trimesh(new trimesh::TriMesh());
        if (trimesh == nullptr) return hexas;
        CMesh cmesh(trimesh);
        //1.重新整理输入参数
        /*trimesh::vec3 dir=adjustHoneyCombParam(trimesh, honeyparams);
        cmesh.Rotate(dir, trimesh::vec3(0, 0, 1));*/
        //2.找到生成蜂窝指定的区域（自定义或者是用户自己指定）          
        if (honeyparams.faces.empty()) {
            //自定义最大底面朝向
            if (honeyparams.axisDir == trimesh::vec3(0, 0, 0)) {
                honeyLetterOpt letterOpts;
                //inputMesh.WriteSTLFile("inputmesh");
                //第1步，寻找底面（最大平面）朝向
                std::vector<int>bottomFaces;
                trimesh::vec3 dir = cmesh.FindBottomDirection(&bottomFaces);
                cmesh.Rotate(dir, trimesh::vec3(0, 0, -1));
                const auto & minPt = cmesh.mbox.min;
                cmesh.Translate(-minPt);
                letterOpts.bottom.resize(bottomFaces.size());
                letterOpts.bottom.assign(bottomFaces.begin(), bottomFaces.end());
                std::vector<int> honeyFaces;
                honeyFaces.reserve(cmesh.mfaces.size());
                for (int i = 0; i < cmesh.mfaces.size(); ++i) {
                    honeyFaces.emplace_back(i);
                }
                std::sort(bottomFaces.begin(), bottomFaces.end());
                std::vector<int> otherFaces(honeyFaces.size() - bottomFaces.size());
                std::set_difference(honeyFaces.begin(), honeyFaces.end(), bottomFaces.begin(), bottomFaces.end(), otherFaces.begin());
                letterOpts.others = std::move(otherFaces);
                //第2步，平移至xoy平面后底面整平
                cmesh.FlatBottomSurface(&bottomFaces);
                if (debugger) {
                    //显示底面区域
                    TriPolygons polygons;
                    const auto& inPoints = cmesh.mpoints;
                    const auto& inIndexs = cmesh.mfaces;
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
                //第3步，生成底面六角网格
                GenerateBottomHexagons(cmesh, honeyparams, letterOpts, debugger);
                trimesh::TriMesh&& mesh = cmesh.GetTriMesh();
                //mesh.write("result.ply");
                trimesh = &mesh;
                trimesh->need_bbox();
                int row = 200;
                int col = 200;
                std::vector<std::tuple<trimesh::point, trimesh::point, trimesh::point>> Upfaces;
                for (int i : letterOpts.others)
                    Upfaces.push_back(std::make_tuple(trimesh->vertices[trimesh->faces[i][0]], trimesh->vertices[trimesh->faces[i][1]],
                        trimesh->vertices[trimesh->faces[i][2]]));
                topomesh::SolidTriangle upST(&Upfaces, row, col, trimesh->bbox.max.x, trimesh->bbox.min.x, trimesh->bbox.max.y, trimesh->bbox.min.y);
                upST.work();

                hexas.side = letterOpts.side;

                /*std::random_device rd;
                std::mt19937 engine(rd());
                std::uniform_real_distribution<double> ldist(0.0, 1.0);
                std::uniform_real_distribution<double> udist(12, 15);  */
                trimesh::point max_xy = trimesh->bbox.max;
                trimesh::point min_xy = trimesh->bbox.min;
                float lengthx = (trimesh->bbox.max.x - trimesh->bbox.min.x) / (col * 1.f);
                float lengthy = (trimesh->bbox.max.y - trimesh->bbox.min.y) / (row * 1.f);
                trimesh::TriMesh * pointmesh = new trimesh::TriMesh();
                for (auto& hg : letterOpts.hexgons) {
                    std::vector<float> height;
                    std::vector<trimesh::ivec2> coord;
                    int max_xi = std::numeric_limits<int>::min();
                    int min_xi = std::numeric_limits<int>::max();
                    int max_yi = std::numeric_limits<int>::min();
                    int min_yi = std::numeric_limits<int>::max();
                    for (int i = 0; i < hg.poly.size(); i++) {
                        trimesh::point p = hg.poly[i] - min_xy;
                        int xi = p.x / lengthx;
                        int yi = p.y / lengthy;
                        xi = xi == col ? --xi : xi;
                        yi = yi == row ? --yi : yi;
                        if (xi > max_xi) max_xi = xi;
                        if (xi < min_xi) min_xi = xi;
                        if (yi > max_yi) max_yi = yi;
                        if (yi < min_yi) min_yi = yi;
                        float min_z = upST.getDataMinZCoord(xi, yi);
                        if (min_z != std::numeric_limits<float>::max()) {
                            min_z -= honeyparams.shellThickness;
                            height.push_back(min_z);
                            pointmesh->vertices.push_back(trimesh::point(p.x, p.y, min_z));
                        } else {
                            height.push_back(0.f);
                            pointmesh->vertices.push_back(trimesh::point(p.x, p.y, 0.f));
                        }
                        coord.push_back(trimesh::ivec2(xi, yi));
                    }
#if 0
                    for (int i = 0; i < height.size(); i++) {
                        if (height[i] == 0.f) continue;
                        int next = (i + 1) % height.size();
                        trimesh::point c = hg.poly[next] - hg.poly[i];
                        float ch = height[next] - height[i];
                        int cx = coord[next].x - coord[i].x;
                        int cy = coord[next].y - coord[i].y;
                        if (cx == 0 && cy == 0) continue;
                        cx = std::abs(cx);
                        cy = std::abs(cy);
                        if (cx > cy) {
                            cx += 1;
                            trimesh::point t = c / (cx * 1.f);
                            float th = ch / (cx * 1.f);
                            float max_c = std::numeric_limits<float>::max();
                            for (int j = 1; j < cx; j++) {
                                trimesh::point tt = hg.poly[i] + j * 1.f * t;
                                float min_z = upST.getDataMinZ(tt.x, tt.y) - honeyparams.shellThickness;
                                float temp_z = height[i] + j * 1.f * th;
                                if (temp_z > min_z) {
                                    if (max_c > (temp_z - min_z))
                                        max_c = (temp_z - min_z);
                                }
                            }
                            if (max_c != std::numeric_limits<float>::max()) {
                                if (height[next] > max_c)
                                    height[next] -= max_c;
                                if (height[i] > max_c)
                                    height[i] -= max_c;
                            }
                        } else {
                            cy += 1;
                            trimesh::point t = c / (cy * 1.f);
                            float th = ch / (cy * 1.f);
                            float max_c = std::numeric_limits<float>::max();
                            for (int j = 1; j < cy; j++) {
                                trimesh::point tt = hg.poly[i] + j * 1.f * t;
                                float min_z = upST.getDataMinZ(tt.x, tt.y) - honeyparams.shellThickness;
                                float temp_z = height[i] + j * 1.f * th;
                                if (temp_z > min_z) {
                                    if (max_c > (temp_z - min_z))
                                        max_c = (temp_z - min_z);
                                }
                            }
                            if (max_c != std::numeric_limits<float>::max()) {
                                if (height[next] > max_c)
                                    height[next] -= max_c;
                                if (height[i] > max_c)
                                    height[i] -= max_c;
                            }
                        }
                    }
#else
                    float last_z = std::numeric_limits<float>::max();
                    for (int i = min_yi; i <= max_yi; i++) {
                        for (int j = min_xi; j <= max_xi; j++) {
                            float min_z = upST.getDataMinZCoord(j, i);
                            if (min_z < last_z)
                                last_z = min_z;
                        }
                    }
                    if (last_z >= honeyparams.shellThickness)
                        last_z -= honeyparams.shellThickness;
                    last_z += 1.0f;
#endif
                    hg.edges.resize(hg.poly.size());
                    for (int i = 0; i < hg.edges.size(); i++) {
                        //hg.edges[i].lowHeight = ldist(engine);
                       // hg.edges[i].topHeight = height[i];
                        hg.edges[i].topHeight = last_z;
                    }
                    hexas.polys.push_back(hg);
                }
                std::vector<bool> deletefaces(trimesh->faces.size(), false);
                for (int f = 0; f < bottomFaces.size(); f++) {
                    deletefaces[bottomFaces[f]] = true;
                }

                trimesh::remove_faces(trimesh, deletefaces);
                trimesh::remove_unused_vertices(trimesh);
                //trimesh->write("trimesh.ply");
                //pointmesh->write("pointmesh.ply");
                
                return hexas;
            } else {
                trimesh::vec3 dir = honeyparams.axisDir;
                trimesh::apply_xform(trimesh, trimesh::xform::rot_into(dir, trimesh::vec3(0, 0, 1)));
                //Rotate(dir, trimesh::vec3(0, 0, 1));

                std::vector<int> botfaces;
                trimesh::TriMesh* newmesh = findOutlineOfDir(trimesh, botfaces);
                trimesh->need_bbox();
                std::vector<std::tuple<trimesh::point, trimesh::point, trimesh::point>> facecontianer;
                for (trimesh::TriMesh::Face& f : newmesh->faces)
                    facecontianer.push_back(std::make_tuple(newmesh->vertices[f[0]], newmesh->vertices[f[1]], newmesh->vertices[f[2]]));
                topomesh::SolidTriangle soildtri(&facecontianer, 100, 100, trimesh->bbox.max.x, trimesh->bbox.min.x, trimesh->bbox.max.y, trimesh->bbox.min.y);
                soildtri.work();
                return hexas;
            }
        } else {
            //用户指定方向，需要计算指定面片的轮廓
            trimesh::vec3 dir = honeyparams.axisDir;
            trimesh::apply_xform(trimesh, trimesh::xform::rot_into(dir, trimesh::vec3(0, 0, 1)));
            return hexas;
        }

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