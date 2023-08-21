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

    struct HexaPolygon {
        TriPolygon poly;
        float radius = 0.0f;
        float ratio = 0.5f;
        float depth = 10.0f;
        float ctop = 0.0f;
        int startIndex = 0;
        trimesh::ivec3 coord;
        std::vector<int> neighbors;
        std::vector<bool> canAdds;
        std::vector<bool> hasAdds;
        std::vector<int> holeIndexs;
        std::vector<std::vector<int>> corners;
        HexaPolygon() : neighbors(6, -1), canAdds(6, false), hasAdds(6, false), holeIndexs(6, -1), corners(6) {}
    };
    typedef std::vector<HexaPolygon> HexaPolygons;
    TriPolygons convertFromHexaPolygons(const HexaPolygons& hexas);
    HexaPolygons convertFromTriPolygons(const TriPolygons& polys, const trimesh::ivec3& coords);
}