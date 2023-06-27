#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#ifndef EPS
#define EPS 1E-8
#endif // !EPS

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif // !M_PI

#ifndef SQRT3
#define SQRT3 1.732050807568877
#endif // !SQRT3
namespace honeycomb {
    inline double Min(const double a, const double b)
    {
        return a < b ? a : b;
    }
    inline double Max(const double a, const double b)
    {
        return a > b ? a : b;
    }
    //二维空间点
    class Point2d {
    public:
        double x = 0, y = 0;
        Point2d() :x(0), y(0) {}
        Point2d(const Point2d& p)
        {
            x = p.x;
            y = p.y;
        }
        Point2d(std::vector<double>& datas)
        {
            const auto& size = datas.size();
            if (size <= 0) return;
            else if (size <= 1) {
                x = datas[0];
            } else {
                x = datas[0];
                y = datas[1];
            }
        }
        Point2d(double x0, double y0)
        {
            x = x0, y = y0;
        }
        inline double Dist(const Point2d& p)const
        {
            return (*this - p).Norm();
        }
        inline double Norm()const
        {
            return std::sqrt(x * x + y * y);
        }
        inline double Norm2()const
        {
            return x * x + y * y;
        }
        inline void Normalize()
        {
            const double norm = sqrt(x * x + y * y);
            if (norm < EPS) {
                return;
            } else {
                x /= norm;
                y /= norm;
            }
        }
        inline Point2d Normalized()const
        {
            const double norm = sqrt(x * x + y * y);
            if (norm < EPS) {
                return Point2d(0.0, 0.0);
            } else {
                return Point2d(x / norm, y / norm);
            }
        }
        inline double Dot(const Point2d& p)const
        {
            return x * p.x + y * p.y;
        }
        inline double Cross(const Point2d& p)const
        {
            return x * p.y - y * p.x;
        }
        inline Point2d operator-()const
        {
            return Point2d(-x, -y);
        }
        inline Point2d operator+(const Point2d& p)const
        {
            return Point2d(x + p.x, y + p.y);
        }
        inline Point2d operator-(const Point2d& p)const
        {
            return Point2d(x - p.x, y - p.y);
        }
        inline Point2d operator*(const double factor)const
        {
            return Point2d(x * factor, y * factor);
        }
        inline Point2d operator/(const double factor)const
        {
            if (std::fabs(factor) < EPS) {
                fputs("Mod Zero Error", stderr);
                exit(1);
            }
            return Point2d(x / factor, y / factor);
        }
        inline void operator*=(const Point2d& factor)
        {
            x *= factor.x;
            y *= factor.y;
        }
        inline void operator/=(const Point2d& factor)
        {
            x /= factor.x;
            y /= factor.y;
        }
        inline void operator+=(const Point2d& p)
        {
            x += p.x;
            y += p.y;
        }
        inline bool operator<(const Point2d& p) const
        {
            auto isEqual = [&](double a, double b, double eps = EPS) {
                return std::fabs(a - b) < eps;
            };
            return x < p.x || (isEqual(x, p.x) && y < p.y);
        }
        inline bool operator==(const Point2d& p) const
        {
            auto isEqual = [&](double a, double b, double eps = EPS) {
                return std::fabs(a - b) < eps;
            };
            return isEqual(x, p.x) && isEqual(y, p.y);
        }
        inline void Print()const
        {
            std::cout << "x=" << x << ",y=" << y << std::endl;
        }
        ~Point2d() {}
    };
    //三维空间点
    class Point {
    public:
        double x = 0, y = 0, z = 0;
        Point() :x(0), y(0), z(0) {}
        Point(const Point& p)
        {
            x = p.x;
            y = p.y;
            z = p.z;
        }
        
        Point(std::vector<double>& datas)
        {
            const auto& size = datas.size();
            if (size <= 0) return;
            else if (size <= 1) {
                x = datas[0];
            }
            else if (size <= 2) {
                x = datas[0];
                y = datas[1];
            } else {
                x = datas[0];
                y = datas[1];
                z = datas[2];
            }
        }
        Point(double x0, double y0, double z0)
        {
            x = x0, y = y0, z = z0;
        }

        inline double Dist(const Point& p)const
        {
            return (*this - p).Norm();
        }
        inline double Norm()const
        {
            return std::sqrt(x * x + y * y + z * z);
        }
        inline double Norm2()const
        {
            return x * x + y * y + z * z;
        }
        inline void Normalize()
        {
            const double norm = sqrt(x * x + y * y + z * z);
            if (norm < EPS) {
                return;
            } else {
                x /= norm;
                y /= norm;
                z /= norm;
            }
        }
        inline Point Normalized()const
        {
            const double norm = sqrt(x * x + y * y + z * z);
            if (norm < EPS) {
                return Point(0, 0, 0);
            } else {
                return Point(x / norm, y / norm, z / norm);
            }
        }
        inline double Dot(const Point & p)const
        {
            return x * p.x + y * p.y + z * p.z;
        }
        inline Point Cross(const Point & p)const
        {
            return Point(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
        }
        inline Point operator-()const
        {
            return Point(-x, -y, -z);
        }
        inline Point operator+(const Point & p)const
        {
            return Point(x + p.x, y + p.y, z + p.z);
        }
        inline Point operator-(const Point & p)const
        {
            return Point(x - p.x, y - p.y, z - p.z);
        }
        inline Point operator*(const double factor)const
        {
            return Point(x * factor, y * factor, z * factor);
        }
        inline Point operator/(const double factor)const
        {
            if (std::fabs(factor) < EPS) {
                fputs("Mod Zero Error", stderr);
                exit(1);
            }
            return Point(x / factor, y / factor, z / factor);
        }
        inline void operator+=(const Point& p)
        {
            x += p.x;
            y += p.y;
            z += p.z;
        }
        inline bool operator<(const Point & p) const
        {
            auto isEqual = [&](double a, double b, double eps = EPS) {
                return std::fabs(a - b) < eps;
            };
            return x < p.x || (isEqual(x, p.x) && y < p.y) || (isEqual(x, p.x) && isEqual(y, p.y) && z < p.z);
        }
        inline bool operator==(const Point & p) const
        {
            auto isEqual = [&](double a, double b, double eps = EPS) {
                return std::fabs(a - b) < eps;
            };
            return isEqual(x, p.x) && isEqual(y, p.y) && isEqual(z, p.z);
        }
        inline void Print()const
        {
            std::cout << "x=" << x << ",y=" << y << ",z=" << z << std::endl;
        }
        ~Point() {}
    };
    class BoundBox2d {
    private:
        Point2d min, max;
    public:
        BoundBox2d(const BoundBox2d& box)
        {
            min = box.min;
            max = box.max;
        }
        BoundBox2d(const Point2d& a, const Point2d& b) :min(a), max(b) {}
        Point2d Centroid() const
        {
            return (min + max) / 2.0;
        }
        Point2d Min()const
        {
            return min;
        }
        Point2d Max()const
        {
            return max;
        }
        double Diagonal()const
        {
            return (min - max).Norm();
        }
        void Translate(const Point2d & trans)
        {
            min += trans;
            max += trans;
        }
        ~BoundBox2d() {}
    };
    class BoundBox {
    public:
        BoundBox(const BoundBox& box)
        {
            min = box.min;
            max = box.max;
        }
        BoundBox(const Point& a,const Point&b):min(a),max(b){}
        Point Min() const
        {
            return min;
        }
        Point Max()const
        {
            return max;
        }
        Point Centroid() const
        {
            return (min + max) / 2.0;
        }
        double Diagonal()const
        {
            return (min - max).Norm();
        }
        void Translate(const Point& trans)
        {
            min += trans;
            max += trans;
        }
        ~BoundBox(){}
    private:
        Point min, max;
    };

    inline std::vector<Point> WrapToPoints(const std::vector<Point2d>& plist)
    {
        const size_t nums = plist.size();
        std::vector<Point> points;
        points.reserve(nums);
        for (int i = 0; i < nums; ++i) {
            const auto& p = plist[i];
            points.emplace_back(Point(p.x, p.y, 0));
        }
        return points;
    }
    inline std::vector<std::vector<Point>> WrapToPoints(const std::vector<std::vector<Point2d>>& plist)
    {
        const size_t nums = plist.size();
        std::vector<std::vector<Point>> points;
        points.reserve(nums);
        for (int i = 0; i < nums; ++i) {
            const auto& ps = WrapToPoints(plist[i]);
            points.emplace_back(ps);
        }
        return points;
    }
}