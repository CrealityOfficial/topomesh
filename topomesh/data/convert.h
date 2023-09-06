#pragma once
#include <map>
#include "topomesh/interface/idata.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // !M_PI
#ifndef M_PI_2
#define M_PI_2 1.570796326794897
#endif // !M_PI_2

#ifndef EPS
#define EPS 1E-8F
#endif // !EPS

namespace topomesh
{
	struct CameraParam
	{	
		trimesh::ivec2 ScreenSize;  // hw
		trimesh::ivec2 p1;		
		trimesh::ivec2 p2;
		//
		float n;
		float f;		
		float t;
		float b;
		float l;
		float r;

		float fov,aspect;

		//
		trimesh::point pos;
		trimesh::point lookAt;
		trimesh::point right;
		trimesh::point up;
		trimesh::point dir;
	};

    struct HexaEdge {
        int neighbor = -1;
        bool canAdd = false;
        bool hasAdd = false;
        bool addRect = false;
        bool addTriangle = false;
        int holenums = -1;
        bool bmutihole = false;

        std::vector<int> starts;
        int lowerHoles = 0;
        std::vector<int> holeIndexs;
        std::vector<int> pointIndexs;
        std::vector<int> corners;
        float lowHeight = 0.f;
        float topHeight = 0.f;
    };
    struct Hexagon {
        trimesh::vec2 centroid;
        float radius = 1;
        std::vector<trimesh::vec2> border;
        Hexagon(const trimesh::vec2& c, const float& d)
        {
            centroid = c, radius = d;
            const auto& theta = 2.0 * M_PI / 6;
            for (int i = 0; i < 6; ++i) {
                const auto& phi = float(i) * theta;
                const auto& p = centroid + trimesh::vec2(std::cos(phi), std::sin(phi)) * radius;
                border.emplace_back(std::move(p));
            }
        }
    };
    struct HexaPolygon {
        bool standard = true;
        trimesh::vec3 center;
        TriPolygon poly;
        int startIndex = 0; ///< 六棱柱第一个点的索引
        trimesh::ivec3 coord; ///<三轴坐标系下的坐标
        std::vector<HexaEdge> edges;
        std::map<int, int> edgemap;
    };
    struct HexaPolygons {
        bool bSewTop = true; ///棱柱的顶部是否需要缝合
        bool bSewBottom = true; ///<棱柱的底部连接部分是否需要缝合
        float side = 0.0f; ///< 每个棱柱底面六角网格的边长
        std::vector<HexaPolygon> polys;
    };

    void translateTriPolygon(TriPolygon& poly, const trimesh::vec3& trans);
    TriPolygons convertFromHexaPolygons(const HexaPolygons& hexas);
    HexaPolygons convertFromTriPolygons(const TriPolygons& polys, const trimesh::ivec3& coords);
    trimesh::TriMesh SaveTriPolygonToMesh(const TriPolygon& poly, double r = 0.01, size_t nslices = 20);
    trimesh::TriMesh SaveTriPolygonsToMesh(const TriPolygons& polys, double r = 0.01, size_t nslices = 20);
}