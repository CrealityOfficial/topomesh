#pragma once
#include "Point.h"
#include "Polyline.h"
#include "Mesh.h"
#include "Matrix.h"
#include "Geometry.h"
#include "Hexagon.h"

namespace honeycomb {
    ////蜂窝结构参数
    //struct HoneyCombOpt {
    //    Polygon* polyline = nullptr; ///<蜂窝填充区域多边形边界，默认底面全部填充
    //    Vec3 axisDir = Vec3(0, 0, 1); ///<三维向量，蜂窝正六棱柱的旋转轴方向，默认垂直于底面
    //    Vec2 arrayDir = Vec2(1, 0); ///<二维向量，蜂窝底面正六边形的布局方向，默认顶点朝y轴方向
    //    double honeyCombRadius = 3; ///<蜂窝底面正六边形网格的边长
    //    double nestWidth = 1; ///<蜂窝底面正六边形与正六边形之间的壁厚
    //    double shellThickness = 2; ///<蜂窝顶部与模型顶部表面空间的抽壳厚度
    //};
    ///*
    //生成底面蜂窝正六边形结构
    //*/
    //void GenerateHoneyCombs(const TriMesh& mesh, TriMesh& resultMesh, const HoneyCombOpt& opt)
    //{

    //}
    void SaveEdgesToMesh(const std::vector<std::vector<Point>>& saveEdges, Mesh& mesh, double r = 0.1, size_t nslices = 5)
    {
        const size_t nums = saveEdges.size();
        auto& points = mesh.GetPoints();
        points.reserve(2 * nums * nslices);
        double delta = 2.0 * M_PI / nslices;
        Point z(0, 0, 1);
        for (size_t i = 0; i < nums; ++i) {
            const auto& a = saveEdges[i][0];
            const auto& b = saveEdges[i][1];
            const auto& n = (b - a).Normalized();
            auto & x = std::move(z.Cross(n));
            if (x.Norm() < EPS) {
                x = Point(1, 0, 0);
            }
            const auto& y = n.Cross(x);
            for (int j = 0; j < nslices; ++j) {
                const auto& theta = delta * j;
                const auto& p = b + x * r * std::cos(theta) + y * r * std::sin(theta);
                points.emplace_back(p);
            }
            for (int j = 0; j < nslices; ++j) {
                const auto& theta = delta * j;
                const auto& p = a + x * r * std::cos(theta) + y * r * std::sin(theta);
                points.emplace_back(p);
            }
        }
        mesh.IndexsReserve(2 * nums * nslices);
        for (size_t i = 0; i < nums; ++i) {
            for (size_t j = 0; j < nslices; ++j) {
                const auto& i0 = j + 2 * nslices * i;
                const auto& i1 = (j + 1) % nslices + 2 * nslices * i;
                const auto & j0 = i0 + nslices;
                const auto & j1 = i1 + nslices;
                mesh.AddFace(j0, i1, i0);
                mesh.AddFace(j0, j1, i1);
            }
        }
        mesh.GenerateFaceNormals();
    }

    void SavePointsToMesh(const std::vector<Point> & savePoints, Mesh & mesh, double radius = 0.1, size_t nrows = 4, size_t ncolumns = 5)
    {
        const size_t spherenums = savePoints.size();
        const size_t nums = nrows * ncolumns + 2;
        const size_t pointnums = nums * spherenums;
        const size_t facenums = 2 * (nums - 2) * spherenums;
        mesh.IndexsReserve(facenums);
        auto & points = mesh.GetPoints();
        points.reserve(pointnums);
        for (int k = 0; k < spherenums; ++k) {
            const auto& p = savePoints[k];
            points.emplace_back(p.x, p.y, p.z + radius);
            for (int i = 0; i < nrows; ++i) {
                const auto& phi = M_PI * (i + 1.0) / double(nrows + 1.0);
                const auto & z = radius * std::cos(phi);
                const auto & r = radius * std::sin(phi);
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& theta = 2.0 * M_PI * j / ncolumns;
                    const auto& x = r * std::cos(theta);
                    const auto& y = r * std::sin(theta);
                    points.emplace_back(p.x + x, p.y + y, p.z + z);
                }
            }
            points.emplace_back(p.x, p.y, p.z - radius);
            const auto & maxInx = points.size() - 1;
            const auto & v0 = k * nums;
            //上下底部两部分
            for (size_t i = 0; i < ncolumns; ++i) {
                const auto& i0 = i + 1 + v0;
                const auto& i1 = (i + 1) % ncolumns + 1 + v0;
                mesh.AddFace(v0, i0, i1);
                const auto & j0 = i0 + (nrows - 1) * ncolumns;
                const auto & j1 = i1 + (nrows - 1) * ncolumns;
                mesh.AddFace(j1, j0, maxInx);
            }
            //中间部分
            for (size_t i = 0; i < nrows - 1; ++i) {
                const auto& j0 = i * ncolumns + 1 + v0;
                const auto& j1 = (i + 1) * ncolumns + 1 + v0;
                for (size_t j = 0; j < ncolumns; ++j) {
                    const auto& i0 = j0 + j;
                    const auto& i1 = j0 + (j + 1) % ncolumns;
                    const auto & i2 = j1 + j;
                    const auto & i3 = j1 + (j + 1) % ncolumns;
                    mesh.AddFace(i2, i1, i0);
                    mesh.AddFace(i2, i3, i1);
                }
            }
        }
        mesh.GenerateFaceNormals();
    }

    void GenerateBottomHexagons()
    {
        using namespace honeycomb;
        Matrix3d mat = Matrix3d::GetMatrix(Point(73, -56, 89), Point(12, 103, 891));
        Mesh mesh;
        mesh.ReadFromSTL("种植牙模-实心.stl");
        mesh.GenerateFaceEdgeAdjacency();
        mesh.Rotate(mat);
        mesh.WriteSTLFile("旋转后模型");
        //std::vector<std::vector<int>>& neighbors = mesh.GetFaceNeighborFaces();
        Mesh bottomMesh;
        std::vector<int> resultFaces;
        /*resultFaces = mesh.SelectLargetPlanar(0.95);
        mesh.SaveFacesToMesh(resultFaces, bottomMesh);
        bottomMesh.WriteSTLFile("筛选底面");*/
        Point direc = mesh.FindBottomDirection(&resultFaces);
        mesh.Rotate(direc, Point(0, 0, -1));
        BoundBox bound = mesh.Bound();
        mesh.Translate(-bound.Min());
        mesh.FlatBottomSurface(&resultFaces);
        mesh.WriteSTLFile("平移至xoy平面模型");
        Mesh cutMesh;
        cutMesh.MiniCopy(mesh);
        cutMesh.DeleteFaces(resultFaces, true);
        cutMesh.WriteSTLFile("删除底面后剩余面片");
        std::vector<int> edges;
        cutMesh.SelectIndividualEdges(edges);
        //底面边界轮廓点索引
        std::vector<std::vector<int>>sequentials;
        cutMesh.GetSequentialPoints(edges, sequentials);
        Mesh pointMesh;
        cutMesh.SavePointsToMesh(sequentials[0], pointMesh, 0.1);
        pointMesh.WriteSTLFile("底面边界轮廓点");
        std::vector<Poly2d> polys;
        polys.reserve(sequentials.size());
        const auto & points = cutMesh.GetPoints();
        for (const auto& seq : sequentials) {
            std::vector<Point2d> border;
            border.reserve(seq.size());
            for (const auto& v : seq) {
                const auto& p = points[v];
                border.emplace_back(p.x, p.y);
            }
            Poly2d poly(border);
            polys.emplace_back(poly);
        }
        ExPoly2d boundarys(polys);
        Point2d centroid = boundarys.MassCentroid();
        Hexagon hex(Point2d(30, 12));
        const auto& border = hex.Border();
        bool inside = false;
        for (int i = 0; i < border.size(); ++i) {
            const auto& p = border[i];
            inside = boundarys.InExPoly2d(p);
            if (!inside)break;
        }
        hex.WriteTxtFile("中心六角网格");
        //Point2d centroid = poly.MassCentroid();
    }
}
