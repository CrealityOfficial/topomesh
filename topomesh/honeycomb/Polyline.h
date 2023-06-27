#pragma once

#include <vector>
#include "Point.h"
#include "Mesh.h"
#include "Hexagon.h"
namespace honeycomb {
    class Poly2d {
    public:
        Poly2d() {}
        Poly2d(const std::vector<Point2d>& pts) :points(std::move(pts)) {}
        std::vector<Point2d> Points()const
        {
            return points;
        }
        void AddPoint(const Point2d& p)
        {
            points.emplace_back(std::move(p));
        }
        void AddPoints(const std::vector<Point2d>& ps)
        {
            points.insert(points.end(), ps.begin(), ps.end());
        }
        bool IsClosed()const
        {
            return points.front() == points.back();
        }
        void SetClosed()
        {
            if (IsClosed()) return;
            points.emplace_back(std::move(points.front()));
        }
        void Reverse()
        {
            std::reverse(points.begin(), points.end());
        }
        double Area()const
        {
            std::vector<Point2d> plist = std::move(points);
            if (!IsClosed()) {
                plist.emplace_back(std::move(plist.front()));
            }
            const auto& nums = plist.size();
            double area = 0;
            for (size_t i = 0; i < nums - 1; ++i) {
                const auto& p0 = plist[i];
                const auto& p1 = plist[i + 1];
                area += (p0.x * p1.y - p1.x * p0.y);
            }
            return area * 0.5;
        }
        size_t size()const
        {
            return points.size();
        }
        void Expand(const Point2d& factor)
        {
            for (auto& p : points) {
                p *= factor;
            }
        }
        void Shrink(const Point2d& factor)
        {
            for (auto& p : points) {
                p /= factor;
            }
        }
        //计算多边形的物质中心
        Point2d MassCentroid()const
        {
            double area = 0;
            Point2d c(0, 0);
            const size_t len = points.size();
            for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
                const auto& a = points[j];
                const auto& b = points[i];
                const auto& f = a.Cross(b);
                c += (a + b) * f;
                area += f * 3;
            }
            if (std::fabs(area) < EPS) return points.at(0);
            return c / area;
        }
        BoundBox2d Bound()const
        {
            constexpr double minValue = std::numeric_limits<double>::lowest();
            constexpr double maxValue = std::numeric_limits<double>::max();
            Point2d min(maxValue, maxValue);
            Point2d max(minValue, minValue);
            for (const auto& p : points) {
                if (p.x < min.x) min.x = p.x;
                if (p.y < min.y) min.y = p.y;
                if (p.x > max.x) max.x = p.x;
                if (p.y > max.y) max.y = p.y;
            }
            return BoundBox2d(min, max);
        }
        bool operator<(const Poly2d& poly) const
        {
            return std::fabs(Area()) < std::fabs(poly.Area());
        }
        ~Poly2d() {}
    private:
        std::vector<Point2d> points;
    };
    class Hole :public Poly2d {

    };
    class ExPoly2d {
    private:
        std::vector<Poly2d> polys;
    public:
        ExPoly2d(){}
        ExPoly2d(const ExPoly2d&plys):polys(plys.polys){
        }
        ExPoly2d(const std::vector<Poly2d>&plys):polys(plys){
            SortPoly2d();
        }
        std::vector<Poly2d> Polys() const
        {
            return polys;
        }
        void Add(const ExPoly2d& plys)
        {
            const auto& pos = plys.Polys();
            polys.insert(polys.end(), pos.begin(), pos.end());
            SortPoly2d();
        }
        void Copy(const ExPoly2d& plys)
        {
            const auto& otherPolys = plys.Polys();
            for (const auto& ps : otherPolys) {
                polys.emplace_back(std::move(ps));
            }
        }
        std::vector<Point2d> Points()const
        {
            std::vector<Point2d> points;
            points.reserve(size() * Outer().size());
            for (const auto& poly : polys) {
                const auto& plist = poly.Points();
                for (const auto& p : plist) {
                    points.emplace_back(std::move(p));
                }
            }
            return points;
        }
        std::vector<std::vector<Point2d>> Edges()const
        {
            std::vector<std::vector<Point2d>> points;
            points.reserve(size() * Outer().size());
            for (const auto& poly : polys) {
                const auto& plist = poly.Points();
                const auto& len = plist.size();
                for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
                    std::vector<Point2d> edge;
                    edge.reserve(2);
                    const auto& a = plist[j];
                    const auto& b = plist[i];
                    edge.emplace_back(a);
                    edge.emplace_back(b);
                    points.emplace_back(edge);
                }
                
            }
            return points;
        }
        bool IsValid()const
        {
            return polys.size() >= 1;
        }
        void Expand(const Point2d& factor)
        {
            for (auto& poly : polys) {
                poly.Expand(factor);
            }
        }
        void Shrink(const Point2d& factor)
        {
            for (auto& poly : polys) {
                poly.Shrink(factor);
            }
        }
        bool HasHoles()const
        {
            return polys.size() >= 2;
        }
        Poly2d Outer()const
        {
            return polys.front();
        }
        std::vector<Poly2d> Holes()const
        {
            const auto& nums = polys.size();
            std::vector<Poly2d> plys;
            plys.reserve(nums - 1);
            for (int i = 1; i < polys.size(); ++i) {
                plys.emplace_back(polys[i]);
            }
            return plys;
        }
        //按照轮廓大小作排序，默认降序
        void SortPoly2d(bool descend = true)
        {
            std::sort(polys.begin(), polys.end());
            if (descend) {
                std::reverse(polys.begin(), polys.end());
            }
            for (auto& poly : polys) {
                if (poly.Area() < 0) {
                    poly.Reverse();
                }
            }
        }
        BoundBox2d Bound()const
        {
            constexpr double minValue = std::numeric_limits<double>::lowest();
            constexpr double maxValue = std::numeric_limits<double>::max();
            Point2d min(maxValue, maxValue);
            Point2d max(minValue, minValue);
            for (const auto& poly : polys) {
                const auto& ps = poly.Points();
                for(const auto&p:ps){
                    if (p.x < min.x) min.x = p.x;
                    if (p.y < min.y) min.y = p.y;
                    if (p.x > max.x) max.x = p.x;
                    if (p.y > max.y) max.y = p.y;
                }
            }
            return BoundBox2d(min, max);
        }
        bool InExPoly2d(const Point2d& p)
        {
            bool inside = false;
            BoundBox2d& bound = Bound();
            Point2d min = bound.Min(), max = bound.Max();
            if (p.x <= min.x || p.x >= max.x || p.y <= min.y || p.y >= max.y) {
                return inside;
            }
            for (const auto& poly : polys) {
                const auto& points = poly.Points();
                const size_t len = points.size();
                for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
                    const auto& a = points[i];
                    const auto& b = points[j];
                    if (((a.y > p.y) != (b.y > p.y)) &&
                        ((p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x))) inside = !inside;
                }
            }
            return inside;
        }
        Point2d MassCentroid()const
        {
            double area = 0;
            Point2d c(0, 0);
            for (const auto& poly : polys) {
                const auto& plist = poly.Points();
                const size_t len = plist.size();
                for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
                    const auto& a = plist[j];
                    const auto& b = plist[i];
                    auto f = a.Cross(b);
                    /*c.x += (a.x + b.x) * f;
                    c.y += (a.y + b.y) * f;*/
                    c += (a + b) * f;
                    area += f * 3;
                }
            }
            if (std::fabs(area) < EPS) return polys.at(0).Points().at(0);
            return c / area;
        }
        void SetClosed()
        {
            for (auto& poly : polys) {
                poly.SetClosed();
            }
        }
        size_t size()const
        {
            return polys.size();
        }
        
        //正值向内偏移，负值向外偏移
        ExPoly2d Offset(const double offset = 0.2)const
        {
            std::vector<Poly2d> newPolys;
            newPolys.reserve(size());
            const auto& outer = Outer();
            const auto& holes = Holes();
            const auto& outPoints = outer.Points();
            const size_t outerSize = outer.size();
            std::vector<Point2d> newOuterPoints;
            newOuterPoints.reserve(outerSize);
            for (size_t i = 0; i < outerSize; ++i) {
                const auto& prev = outPoints[(i == 0) ? (outerSize - 1) : (i - 1)];
                const auto& next = outPoints[(i + 1 == outerSize) ? 0 : (i + 1)];
                const auto& curv = outPoints[i];
                const auto& dir1 = prev - curv;
                const auto& dir2 = next - curv;
                const auto& len1 = dir1.Norm();
                const auto& len2 = dir2.Norm();
                const auto& dir = (dir1 * len2 + dir2 * len1) / (len1 + len2);
                if (dir.Norm() < EPS) continue;
                const auto& nor = dir.Normalized();
                const auto& flag = dir2.Cross(nor) > 0 ? 1 : -1;
                const auto& vec = nor * flag * offset;
                const auto& newPt = curv + vec;
                newOuterPoints.emplace_back(std::move(newPt));
            }
            newPolys.emplace_back(newOuterPoints);
            for (auto& hole : holes) {
                const auto& holePoints = hole.Points();
                const auto& holeSize = hole.size();
                std::vector<Point2d> newHolePoints;
                newHolePoints.reserve(holeSize);
                for (size_t i = 0; i < holeSize; ++i) {
                    const auto& prev = holePoints[(i == 0) ? (holeSize - 1) : (i - 1)];
                    const auto & next = holePoints[(i + 1 == holeSize) ? 0 : (i + 1)];
                    const auto & curv = holePoints[i];
                    const auto& dir1 = prev - curv;
                    const auto& dir2 = next - curv;
                    const auto& len1 = dir1.Norm();
                    const auto& len2 = dir2.Norm();
                    const auto& dir = (dir1 * len2 + dir2 * len1) / (len1 + len2);
                    if (dir.Norm() < EPS) continue;
                    const auto & nor = dir.Normalized();
                    const auto & flag = dir2.Cross(nor) > 0 ? -1 : 1;
                    const auto & vec = nor * flag * offset;
                    const auto& newPt = curv + vec;
                    newHolePoints.emplace_back(std::move(newPt));
                }
                newPolys.emplace_back(newHolePoints);
            }
            ExPoly2d newExpoly2d(newPolys);
            return newExpoly2d;
        }
        std::vector<Point2d> GenerateCrossPoints(const Point2d& pt, double xdelta = 1.0)const
        {
            std::vector<Point2d> crossPoints;
            std::vector<double> leftXs, rightXs;
            bool inside = false;
            for (auto& poly : polys) {
                const auto& points = poly.Points();
                const size_t len = points.size();
                for (std::size_t i = 0, j = len - 1; i < len; j = i++) {
                    const auto& a = points[i];
                    const auto& b = points[j];
                    if (a.y == b.y) continue;
                    if (pt.y <= Min(a.y, b.y)) continue;
                    if (pt.y >= Max(a.y, b.y)) continue;
                    // 求交点的x坐标（由直线两点式方程转化而来）  
                    double x = (double)(pt.y - a.y) * (double)(b.x - a.x) / (double)(b.y - a.y) + a.x;
                    // 统计p1p2与p向右射线的交点及左射线的交点  
                    if (pt.x < x) {
                        rightXs.push_back(x);
                        inside = !inside;
                    } else {
                        leftXs.push_back(x);
                    }
                }
            }
            std::sort(leftXs.begin(), leftXs.end());
            std::sort(rightXs.begin(), rightXs.end());
            const size_t lns = leftXs.size();
            const size_t rns = rightXs.size();
            const double w = xdelta;
            const double y = pt.y;
            if (inside) {
                leftXs.insert(leftXs.end(), rightXs.begin(), rightXs.end());
                for (size_t i = 0; i < lns + rns; i += 2) {
                    const auto& xl = leftXs[i];
                    const auto& xr = leftXs[i + 1];
                    const auto& lnums = std::ceil(xl / w);
                    const auto& rnums = std::floor(xr / w);
                    for (int xnum = lnums; xnum <= rnums; ++xnum) {
                        crossPoints.emplace_back(std::move(Point2d(xnum * w, y)));
                    }
                }
            } else {
                // left part
                for (size_t i = 0; i < lns; i += 2) {
                    const auto& xl = leftXs[i];
                    const auto& xr = leftXs[i + 1];
                    const auto& lnums = std::ceil(xl / w);
                    const auto & rnums = std::floor(xr / w);
                    for (int xnum = lnums; xnum <= rnums; ++xnum) {
                        crossPoints.emplace_back(std::move(Point2d(xnum * w, y)));
                    }
                }
                // right part
                for (size_t i = 0; i < rns; i += 2) {
                    const auto& xl = rightXs[i];
                    const auto& xr = rightXs[i + 1];
                    const auto& lnums = std::ceil(xl / w);
                    const auto& rnums = std::floor(xr / w);
                    for (int xnum = lnums; xnum <= rnums; ++xnum) {
                        crossPoints.emplace_back(std::move(Point2d(xnum * w, y)));
                    }
                }
            }
            return crossPoints;
        }
        std::vector<std::vector<Point2d>> GenerateGridPoints(double xdelta = 1.0, double ydelta = 1.0)const
        {
            const BoundBox2d bound = Bound();
            const Point2d& max = bound.Max();
            const Point2d& min = bound.Min();
            const auto& dy = std::ceil(min.y / ydelta);
            const auto& uy = std::floor(max.y / ydelta);
            const auto& lx = std::ceil(min.x / xdelta);
            const auto& rx = std::floor(max.x / xdelta);
            const auto& xmedim = (max.x + min.x) / 2.0;
            const auto& xnums = std::ceil(xmedim / xdelta);
            const auto& x = xnums * xdelta;
            std::vector<std::vector<Point2d>> gridPoints;
            for (int ynum = uy; ynum >= dy; --ynum) {
                /*std::vector<Point2d> crossPts;
                for (int xnum = lx; xnum < rx; ++xnum) {
                    crossPts.emplace_back(xnum * xdelta, ynum * ydelta);
                }*/
                const std::vector<Point2d>& crossPts = GenerateCrossPoints(Point2d(x, ynum * ydelta), xdelta);
                //在多边形内部横纵轴分别均匀采样
                gridPoints.emplace_back(crossPts);
            }
            return gridPoints;
        }
        std::vector<Hexagon> GenerateHexagons(const HexagonOpt& opt, std::vector<std::vector<Point2d>>& borders) const
        {
            std::vector<Hexagon> hexagons;
            const double radius = opt.hexagonRadius;
            const double nestWidth = opt.nestWidth;
            const double thickness = opt.shellThickness;
            const double xdist = 3.0 / 2.0 * radius;
            const double ydist = SQRT3 / 2.0 * radius;
            const double side = radius - nestWidth / SQRT3;
            const double offdist = Min(xdist, thickness);
            ExPoly2d polys = Offset(xdist + thickness);
            if (true) {
                Mesh polyMesh;
                ExPoly2d tempPolys;
                tempPolys.Copy(polys);
                tempPolys.Add(*this);
                tempPolys.SavePolysToMesh(polyMesh);
                polyMesh.WriteSTLFile("向内偏移多边形");
            }
            const auto& gridPoints = polys.GenerateGridPoints(xdist, ydist);
            std::cout << xdist << " " << ydist <<" "<<radius << std::endl;
            const size_t nrows = gridPoints.size();
            hexagons.reserve(nrows);
            //计算六角网格对应边界包围盒坐标奇偶性质
            const BoundBox2d bound = Bound();
            const Point2d& min = bound.Min();
            for (size_t row = 0; row < nrows; row += 2) {
                const auto& rowPoints = gridPoints[row];
                const auto& ncols = rowPoints.size();
                for (size_t i = 0; i < ncols; ++i) {
                    const auto& cur = rowPoints[i];
                    const int& num = int((cur - min).x / xdist);
                    if (num % 2 == 0) {
                        const auto& center = cur;
                        Hexagon hexagon(center, side);
                        hexagons.emplace_back(std::move(hexagon));
                        const auto& edges = hexagon.BorderEdges();
                        borders.insert(borders.end(), edges.begin(), edges.end());
                    } else {
                        const auto& center = cur + Point2d(0, -ydist);
                        Hexagon hexagon(center, side);
                        hexagons.emplace_back(std::move(hexagon));
                        const auto& edges = hexagon.BorderEdges();
                        borders.insert(borders.end(), edges.begin(), edges.end());
                    }
                    
                }
            }
            return hexagons;
        }

        void SavePolysToMesh(Mesh& mesh, double r = 0.1, size_t nslices = 5)
        {
            const auto& edges = Edges();
            const auto& edgePoints = WrapToPoints(edges);
            const size_t nums = edges.size();
            auto& points = mesh.GetPoints();
            points.reserve(2 * nums * nslices);
            double delta = 2.0 * M_PI / nslices;
            Point z(0, 0, 1);
            for (size_t i = 0; i < nums; ++i) {
                const auto& a = edgePoints[i][0];
                const auto& b = edgePoints[i][1];
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
    };
    class Polyline {
    public:
        Polyline() {}
        Polyline(const std::vector<Point>& pts) :points(std::move(pts)) {}
        std::vector<Point>& Points()
        {
            return points;
        }
        void AddPoint(const Point& p)
        {
            points.emplace_back(std::move(p));
        }
        void AddPoints(const std::vector<Point>& ps)
        {
            points.insert(points.end(), ps.begin(), ps.end());
        }
        bool IsClosed()const
        {
            return points.front() == points.back();
        }
        void SetClosed()
        {
            if (IsClosed()) return;
            points.emplace_back(std::move(points.front()));
        }
        void Reverse()
        {
            std::reverse(points.begin(), points.end());
        }
        void PolyToMesh(Mesh& mesh, double r = 0.01, size_t nslices = 20)
        {
            const size_t nums = points.size();
            auto& plists = mesh.GetPoints();
            plists.reserve(2 * (nums - 1) * nslices);
            double delta = 2.0 * M_PI / nslices;
            Point z(0, 0, 1);
            for (size_t i = 0; i < nums - 1; ++i) {
                const auto& a = points[i];
                const auto& b = points[i + 1];
                const auto& n = (b - a).Normalized();
                auto & x = std::move(z.Cross(n));
                if (x.Norm() < EPS) {
                    x = Point(1, 0, 0);
                }
                const auto& y = n.Cross(x);
                for (int j = 0; j < nslices; ++j) {
                    const auto& theta = delta * j;
                    const auto& p = b + x * r * std::cos(theta) + y * r * std::sin(theta);
                    plists.emplace_back(p);
                }
                for (int j = 0; j < nslices; ++j) {
                    const auto& theta = delta * j;
                    const auto& p = a + x * r * std::cos(theta) + y * r * std::sin(theta);
                    plists.emplace_back(p);
                }
            }
            mesh.IndexsReserve(2 * (nums - 1) * nslices);
            for (size_t i = 0; i < nums - 1; ++i) {
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
        ~Polyline() {}
    private:
        std::vector<Point> points;
    };
}