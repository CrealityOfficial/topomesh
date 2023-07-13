#pragma once
#include "Point.h"
#include "Polyline.h"
#include "Mesh.h"
#include "Matrix.h"
#include "Geometry.h"
#include "Hexagon.h"

namespace honeycomb {
    ////���ѽṹ����
    //struct HoneyCombOpt {
    //    Polygon* polyline = nullptr; ///<��������������α߽磬Ĭ�ϵ���ȫ�����
    //    Vec3 axisDir = Vec3(0, 0, 1); ///<��ά����������������������ת�᷽��Ĭ�ϴ�ֱ�ڵ���
    //    Vec2 arrayDir = Vec2(1, 0); ///<��ά���������ѵ����������εĲ��ַ���Ĭ�϶��㳯y�᷽��
    //    double honeyCombRadius = 3; ///<���ѵ���������������ı߳�
    //    double nestWidth = 1; ///<���ѵ���������������������֮��ıں�
    //    double shellThickness = 2; ///<���Ѷ�����ģ�Ͷ�������ռ�ĳ�Ǻ��
    //};
    ///*
    //���ɵ�������������νṹ
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
            auto x = std::move(z.Cross(n));
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
            //���µײ�������
            for (size_t i = 0; i < ncolumns; ++i) {
                const auto& i0 = i + 1 + v0;
                const auto& i1 = (i + 1) % ncolumns + 1 + v0;
                mesh.AddFace(v0, i0, i1);
                const auto & j0 = i0 + (nrows - 1) * ncolumns;
                const auto & j1 = i1 + (nrows - 1) * ncolumns;
                mesh.AddFace(j1, j0, maxInx);
            }
            //�м䲿��
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
        mesh.ReadFromSTL("GenerateBottomHexagons0.stl");
        mesh.GenerateFaceEdgeAdjacency();
        mesh.Rotate(mat);
        mesh.WriteSTLFile("GenerateBottomHexagons2");
        //std::vector<std::vector<int>>& neighbors = mesh.GetFaceNeighborFaces();
        Mesh bottomMesh;
        std::vector<int> resultFaces;
        /*resultFaces = mesh.SelectLargetPlanar(0.95);
        mesh.SaveFacesToMesh(resultFaces, bottomMesh);
        bottomMesh.WriteSTLFile("ɸѡ����");*/
        Point direc = mesh.FindBottomDirection(&resultFaces);
        mesh.Rotate(direc, Point(0, 0, -1));
        BoundBox bound = mesh.Bound();
        mesh.Translate(-bound.Min());
        mesh.FlatBottomSurface(&resultFaces);
        mesh.WriteSTLFile("GenerateBottomHexagons3");
        Mesh cutMesh;
        cutMesh.MiniCopy(mesh);
        cutMesh.DeleteFaces(resultFaces, true);
        cutMesh.WriteSTLFile("GenerateBottomHexagons4");
        std::vector<int> edges;
        cutMesh.SelectIndividualEdges(edges);
        //����߽�����������
        std::vector<std::vector<int>>sequentials;
        cutMesh.GetSequentialPoints(edges, sequentials);
        Mesh pointMesh;
        cutMesh.SavePointsToMesh(sequentials[0], pointMesh, 0.1);
        pointMesh.WriteSTLFile("GenerateBottomHexagons5");
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
        hex.WriteTxtFile("GenerateBottomHexagons6");
        //Point2d centroid = poly.MassCentroid();
    }
}
