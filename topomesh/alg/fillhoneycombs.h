#pragma once
#include "topomesh/data/convert.h"
#include "ccglobal/tracer.h"

namespace topomesh {

	struct HoneyCombParam
	{
        double resolution = 1E-4; ///<多边形相交允许误差
        TriPolygon* polyline = nullptr; ///<自由编辑区域生成内蜂窝（默认底面全部填充）
        trimesh::vec3 axisDir = trimesh::vec3(0, 0, 1); ///<蜂窝多边形平面朝向（默认z轴正方向）
        trimesh::vec2 arrayDir = trimesh::vec2(1, 0); ///<蜂窝多边形平面布局（默认边朝上结构）
		double honeyCombRadius = 2.0; ///<未收缩前蜂窝六边形边长
		double nestWidth = 0.3; ///<蜂窝六边形壁厚（向内收缩的距离的2倍）
		double shellThickness = 1.0; ///<抽壳厚度

		//debug
		int step_return = 9999; // debug quick return
	};

	class HoneyCombDebugger
	{
	public:
		virtual void onGenerateBottomPolygons(const TriPolygons& polygons) = 0;
		virtual void onGenerateInfillPolygons(const TriPolygons& polygons) = 0;
	};

    trimesh::TriMesh* generateHoneyCombs(trimesh::TriMesh* trimesh = nullptr, const HoneyCombParam& honeyparams = HoneyCombParam(),
        ccglobal::Tracer* tracer = nullptr, HoneyCombDebugger* debugger = nullptr);

    struct HexagonArrayParam {
        trimesh::vec3 pos; ///<阵列左上角起始点
        int nrows = 2;
        int ncols = 2;
        double radius = 2.0; ///<未收缩前蜂窝六边形边长
        double nestWidth = 0.3; ///<蜂窝六边形壁厚（向内收缩的距离的2倍）
    };
    TriPolygons GenerateHexagons(const HexagonArrayParam& hexagonparams = HexagonArrayParam());


	class MMeshT;
	void GenerateHoneyCombs(const trimesh::TriMesh* mesh, trimesh::TriMesh& resultmesh, const TriPolygon& poly, trimesh::vec3 axisDir=trimesh::vec3(0,0,1),
		trimesh::vec2 arrayDir = trimesh::vec2(0,1), double honeyCombRadius=2.0f, double nestWidth=1.0f, double shellThickness=1.0f);
	void findNeighVertex(trimesh::TriMesh* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, int>>& vertex_distance);
	void findNeighVertex(MMeshT* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, float>>& vertex_distance);
	void innerHex(MMeshT* mesh, std::vector<std::vector<trimesh::vec2>>& poly, std::vector<int>& inFace, std::vector<int>& outFace,float len);
}