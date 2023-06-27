#pragma once
#include "Point.h"
#include "Mesh.h"
namespace honeycomb {
    class Cylinder :public Mesh {
    public:
        Cylinder(Point& a, Point& b, double r = 1.0, size_t nslices = 20)
        {
            // 计算圆柱的轴方向，默认a->b
            Point n = (b - a).Normalized();
            Point z(0, 0, 1);
            Point x = z.Cross(n);
            if (x.Norm() < EPS) {
                x = Point(1, 0, 0);
            }
            Point y = n.Cross(x);
            auto& points = GetPoints();
            points.reserve(2 * nslices);
            double delta = 2.0 * M_PI / nslices;
            for (int i = 0; i < nslices; ++i) {
                const auto& theta = delta * i;
                const auto& p = b + x * r * std::cos(theta) + y * r * std::sin(theta);
                points.emplace_back(p);
            }
            for (int i = 0; i < nslices; ++i) {
                const auto& theta = delta * i;
                const auto& p = a + x * r * std::cos(theta) + y * r * std::sin(theta);
                points.emplace_back(p);
            }
            IndexsReserve(2 * nslices);
            for (size_t i = 0; i < nslices; ++i) {
                const auto& i0 = i;
                const auto& i1 = (i + 1) % nslices;
                const auto & j0 = i0 + nslices;
                const auto & j1 = i1 + nslices;
                AddFace(j0, i1, i0);
                AddFace(j0, j1, i1);
            }
            GenerateFaceNormals();
        }
        ~Cylinder() {}
    };

    class UVSphere : public Mesh {
    public:
        // nrows>=1, ncolumns>=2
        UVSphere(double radius = 1.0, int nrows = 20, int ncolumns = 40)
        {
            auto& points = GetPoints();
            int pointnums = nrows * ncolumns + 2;
            points.reserve(pointnums);
            points.emplace_back(0, 0, radius);
            for (int i = 0; i < nrows; ++i) {
                double phi = M_PI * (i + 1.0) / double(nrows + 1.0);
                const auto & z = radius * std::cos(phi);
                const auto & r = radius * std::sin(phi);
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& theta = 2.0 * M_PI * j / ncolumns;
                    const auto& x = r * std::cos(theta);
                    const auto& y = r * std::sin(theta);
                    points.emplace_back(x, y, z);
                }
            }
            points.emplace_back(0, 0, -radius);
            const int maxInx = points.size() - 1;
            const int facenums = 2 * nrows * ncolumns;
            IndexsReserve(facenums);
            for (int i = 0; i < ncolumns; ++i) {
                auto i0 = i + 1;
                auto i1 = (i + 1) % ncolumns + 1;
                AddFace(0, i0, i1);
                i0 = maxInx - ncolumns + i;
                i1 = maxInx - ncolumns + (i + 1) % ncolumns;
                AddFace(i1, i0, maxInx);
            }
            for (int i = 0; i < nrows - 1; ++i) {
                const auto& j0 = i * ncolumns + 1;
                const auto& j1 = (i + 1) * ncolumns + 1;
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& i0 = j0 + j;
                    const auto& i1 = j0 + (j + 1) % ncolumns;
                    const auto & i2 = j1 + j;
                    const auto & i3 = j1 + (j + 1) % ncolumns;
                    AddFace(i2, i1, i0);
                    AddFace(i2, i3, i1);
                }
            }
            GenerateFaceNormals();
        }

        ~UVSphere() {}
    };

    class UVOvum :public Mesh {
    public:
        //鸡蛋体，方程为(x^2/a^2+y^2/b^2)(1+kz)+z^2/c^2=1;
        UVOvum(double a = 1.0, double b = 2.0, double c = 3.0, double k = 0.1, int nrows = 20, int ncolumns = 40)
        {
            auto& points = GetPoints();
            int pointnums = nrows * ncolumns + 2;
            points.reserve(pointnums);
            points.emplace_back(0, 0, c);
            for (int i = 0; i < nrows; ++i) {
                double phi = M_PI * (i + 1.0) / double(nrows + 1.0);
                const auto & z = c * std::cos(phi);
                const auto & r = std::sqrt((c * c - z * z) / std::fabs(1 + k * z));
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& theta = 2.0 * M_PI * j / ncolumns;
                    const auto& x = a * r * std::cos(theta);
                    const auto& y = b * r * std::sin(theta);
                    points.emplace_back(x, y, z);
                }
            }
            points.emplace_back(0, 0, -c);
            const int maxInx = points.size() - 1;
            const int facenums = 2 * nrows * ncolumns;
            IndexsReserve(facenums);
            for (int i = 0; i < ncolumns; ++i) {
                auto i0 = i + 1;
                auto i1 = (i + 1) % ncolumns + 1;
                AddFace(0, i0, i1);
                i0 = maxInx - ncolumns + i;
                i1 = maxInx - ncolumns + (i + 1) % ncolumns;
                AddFace(i1, i0, maxInx);
            }
            for (int i = 0; i < nrows - 1; ++i) {
                const auto& j0 = i * ncolumns + 1;
                const auto& j1 = (i + 1) * ncolumns + 1;
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& i0 = j0 + j;
                    const auto& i1 = j0 + (j + 1) % ncolumns;
                    const auto & i2 = j1 + j;
                    const auto & i3 = j1 + (j + 1) % ncolumns;
                    AddFace(i2, i1, i0);
                    AddFace(i2, i3, i1);
                }
            }
            GenerateFaceNormals();
        }

        ~UVOvum() {}
    };

    class UVTorus : public Mesh {
    public:
        // 环面，没有上下天穹点，nrows>=2, ncolumns>=2
        UVTorus(double Radius = 3, double radius = 1.0, int nrows = 20, int ncolumns = 40)
        {
            auto& points = GetPoints();
            int pointnums = nrows * ncolumns;
            points.reserve(pointnums);
            for (int i = 0; i < nrows; ++i) {
                double phi = 2.0 * M_PI * double(i) / double(nrows);
                const auto& z = radius * std::cos(phi);
                const auto& r = Radius + radius * std::sin(phi);
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& theta = 2.0 * M_PI * double(j) / double(ncolumns);
                    const auto& x = r * std::cos(theta);
                    const auto& y = r * std::sin(theta);
                    points.emplace_back(x, y, z);
                }
            }
            const int facenums = 2 * nrows * ncolumns;
            IndexsReserve(facenums);
            for (int i = 0; i < nrows; ++i) {
                const auto& j0 = i * ncolumns;
                const auto& j1 = ((i + 1) % nrows) * ncolumns;
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& i0 = j0 + j;
                    const auto& i1 = j0 + (j + 1) % ncolumns;
                    const auto & i2 = j1 + j;
                    const auto & i3 = j1 + (j + 1) % ncolumns;
                    AddFace(i2, i1, i0);
                    AddFace(i2, i3, i1);
                }
            }
            GenerateFaceNormals();
        }

        ~UVTorus() {}
    };
}