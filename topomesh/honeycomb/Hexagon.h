#pragma once
#include "Point.h"
namespace honeycomb {
    struct HexagonOpt {
        Point2d arrayDir = Point2d(1, 0); ///<二维向量，蜂窝底面正六边形的布局方向，默认边朝上
        double hexagonRadius = 1; ///<蜂窝底面正六边形网格的边长
        double nestWidth = 0.2; ///<蜂窝底面正六边形与正六边形之间的壁厚
        double shellThickness = 1; ///<蜂窝顶部与模型顶部表面空间的抽壳厚度
    };
    //默认六角网格边朝上
    class Hexagon {
    private:
        Point2d centroid;
        double radius = 1; ///<初始网格边长
        std::vector<Point2d> border;
    public:
        Hexagon() :centroid(Point2d(0, 0)) {
            const auto& theta = 2.0 * M_PI / 6;
            for (int i = 0; i < 6; ++i) {
                const auto& phi = double(i) * theta;
                const auto& p = centroid + Point2d(std::cos(phi), std::sin(phi)) * radius;
                border.emplace_back(std::move(p));
            }
        }
        Hexagon(const Point2d& c)
        {
            centroid = c;
            const auto& theta = 2.0 * M_PI / 6;
            for (int i = 0; i < 6; ++i) {
                const auto& phi = double(i) * theta;
                const auto& p = centroid + Point2d(std::cos(phi), std::sin(phi)) * radius;
                border.emplace_back(std::move(p));
            }
        }
        Hexagon(const double& d)
        {
            radius = d;
            const auto& theta = 2.0 * M_PI / 6;
            for (int i = 0; i < 6; ++i) {
                const auto& phi = double(i) * theta;
                const auto& p = centroid + Point2d(std::cos(phi), std::sin(phi)) * radius;
                border.emplace_back(std::move(p));
            }
        }
        Hexagon(const Point2d& c, const double& d)
        {
            centroid = c;
            radius = d;
            const auto& theta = 2.0 * M_PI / 6;
            for (int i = 0; i < 6; ++i) {
                const auto& phi = double(i) * theta;
                const auto& p = centroid + Point2d(std::cos(phi), std::sin(phi)) * radius;
                border.emplace_back(std::move(p));
            }
        }
        Hexagon(const Hexagon& hex)
        {
            centroid = hex.centroid;
            radius = hex.radius;
            const auto& theta = 2.0 * M_PI / 6;
            for (int i = 0; i < 6; ++i) {
                const auto& phi = double(i) * theta;
                const auto& p = centroid + Point2d(std::cos(phi), std::sin(phi)) * radius;
                border.emplace_back(std::move(p));
            }
        }
        double Radius() const
        {
            return radius;
        }
        Point2d Centroid() const
        {
            return centroid;
        }
        std::vector<Point2d> Border() const
        {
            return border;
        }
        std::vector<std::vector<Point2d>> BorderEdges()const
        {
            const size_t nums = border.size();
            std::vector<std::vector<Point2d>> edgelist;
            edgelist.reserve(nums);
            for (size_t i = 0; i < nums; ++i) {
                std::vector<Point2d> edge;
                edge.reserve(2);
                const auto& curv = border[i];
                const auto& next = border[i + 1 == nums ? 0 : i + 1];
                edge.emplace_back(std::move(curv));
                edge.emplace_back(std::move(next));
                edgelist.emplace_back(edge);
            }
            return edgelist;
        }
        bool WriteTxtFile(const char* filename)
        {
            std::ofstream fs(std::string(filename) + "_hexagon.txt");
            if (!fs) { fs.close(); return false; }
            fs << "六角网格中心点: \n" << (float)centroid.x << " " << (float)centroid.y << std::endl;
            fs << "六角网格边界点: " << std::endl;
            for (int i = 0; i < border.size(); ++i) {
                const auto& p = border[i];
                fs << (float)p.x << " " << (float)p.y << std::endl;
            }
            fs.close();
            return true;
        }
    };
}