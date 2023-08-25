#pragma once
#include "topomesh/interface/idata.h"

namespace topomesh
{
	struct CameraParam
	{	
		trimesh::ivec2 ScreenSize;  // hw
		trimesh::ivec2 p1;		
		trimesh::ivec2 p2;
		//�ڲ�
		float n;
		float f;		
		float t;
		float b;
		float l;
		float r;

		float fov,aspect;

		//���
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
        std::vector<int> addRects;
        std::vector<int> corners;
        bool blower = true;
        float lowHeight = 0.f;
        float topHeight = 0.f;
    };
    struct HexaPolygon {
        TriPolygon poly;
        int startIndex = 0; ///< 六棱柱第一个点的索引
        trimesh::ivec3 coord; ///<三轴坐标系下的坐标
        std::vector<HexaEdge> edges;
        HexaPolygon() : edges(6) {}
    };
    struct HexaPolygons {
        bool bSewTop = true; ///棱柱的顶部是否需要缝合
        bool bSewBottom = true; ///<棱柱的底部是否需要缝合
        float side = 0.0f; ///< 每个棱柱底面六角网格的边长
        std::vector<HexaPolygon> polys;
    };

    void translateTriPolygon(TriPolygon& poly, const trimesh::vec3& trans);
    TriPolygons convertFromHexaPolygons(const HexaPolygons& hexas);
    HexaPolygons convertFromTriPolygons(const TriPolygons& polys, const trimesh::ivec3& coords);
}