#pragma once
#include "topomesh/data/convert.h"
#include "ccglobal/tracer.h"
#include "topomesh/interface/idata.h"
#include <memory>

namespace topomesh {

	

	class HoneyCombDebugger
	{
	public:
		virtual void onGenerateBottomPolygons(const TriPolygons& polygons) = 0;
		virtual void onGenerateInfillPolygons(const TriPolygons& polygons) = 0;
	};

    trimesh::TriMesh* generateHoneyCombs(trimesh::TriMesh* trimesh = nullptr, const HoneyCombParam& honeyparams = HoneyCombParam(),
        ccglobal::Tracer* tracer = nullptr, HoneyCombDebugger* debugger = nullptr);


	trimesh::TriMesh* GenerateHoneyCombs(trimesh::TriMesh* trimesh = nullptr, const HoneyCombParam& honeyparams = HoneyCombParam(),
		ccglobal::Tracer* tracer = nullptr, HoneyCombDebugger* debugger = nullptr);

    struct HexagonArrayParam {
		trimesh::vec dir=trimesh::vec3(0,0,1);
        trimesh::vec3 pos; ///<阵列左上角起始点
        int nrows = 2;
        int ncols = 2;
        double radius = 2.0; ///<未收缩前蜂窝六边形边长
        double nestWidth = 0.3; ///<蜂窝六边形壁厚（向内收缩的距离的2倍）
    };
    HexaPolygons GenerateHexagons(const HexagonArrayParam& hexagonparams = HexagonArrayParam());
    void GenerateHexagonNeighbors(HexaPolygons& hexas, float cheight = 5.0);
	trimesh::vec3 adjustHoneyCombParam(trimesh::TriMesh* trimesh,const HoneyCombParam& honeyparams);
    TriPolygons traitCurrentPolygons(const HexaPolygons& hexas, int index);
    TriPolygons traitNeighborPolygons(const HexaPolygons& hexas, int index);
    TriPolygons traitDirctionPolygon(const HexaPolygons& hexas, int index, int dir);
    struct ColumnarHoleParam {
        int nslices = 17; ///<圆孔默认正17边形
        float cheight = 5.0f; ///<圆孔圆心高度
        float radius = 1.0f; ///<圆孔半径
    };
    TriPolygon traitPlanarCircle(const trimesh::vec3& c, float r, std::vector<int>& indexs, const trimesh::vec3& edgeDir = trimesh::vec3(0, 0, 1), int nums = 17);
    std::shared_ptr<trimesh::TriMesh> generateHolesColumnar(HexaPolygons& hexas, const ColumnarHoleParam& param);

	class MMeshT;
	void findNeighVertex(trimesh::TriMesh* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, int>>& vertex_distance);
	void findNeighVertex(MMeshT* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, float>>& vertex_distance);
	void innerHex(MMeshT* mesh, std::vector<std::vector<trimesh::vec2>>& poly, std::vector<int>& inFace, std::vector<int>& outFace,float len);

    struct honeyLetterOpt {
        std::vector<int>bottom; ///<底面大块平面的面片索引
        std::vector<int>others; ///<去掉底面后其余面片索引（已保留原模型对应的索引）
        /*
        hexagon: 每个六角网格结构体;
        radius: 包含收缩前的网格边长;
        borders: 收缩后轮廓顶点(不一定6个点，默认xoy坐标系逆时针排序);
        neighbors: 六个方向相邻的六角网格索引(6个方向不一定都有相邻，没有记-1);
        */
        struct hexagon {
            double radius = 1.0;
            std::vector<trimesh::vec3>borders;
            std::vector<int> neighbors;
            hexagon() : neighbors(6, -1) {}
        };
        //所有六角网格及邻居关系
        std::vector<hexagon>hexgons;
    };

    void findHoneyCombsCoord(trimesh::TriMesh* mesh,const  honeyLetterOpt& honeycombs,std::vector<std::vector<std::pair<float,float>>>& coord);
}