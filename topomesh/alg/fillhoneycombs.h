#pragma once
#include "topomesh/data/convert.h"
#include "ccglobal/tracer.h"
#include "topomesh/interface/idata.h"
#include "topomesh/data/CMesh.h"
#include <memory>

namespace topomesh {
	
    struct HoneyCombParam
    {
        int mode = 0; ///< 0 is shell, 1 is backfill.
        TriPolygon* polyline = nullptr; ///< input polygons boundary.
        trimesh::vec3 axisDir = trimesh::vec3(0, 0, 0); ///< default cylindrical orientation.
        //trimesh::vec3* axisDir = nullptr;
        bool isdelect = false; ///< 
        std::vector<int> faces; ///< faces ids of the input flat region.
        // for honeycomb
        double honeyCombRadius = 2.0; ///<to generate honeycomb radius.
        double nestWidth = 1.0; ///<to generate honeycomb nest width.
        double shellThickness = 0.6; ///<exshell thickness of honeycomb model.
        
        // for bottom boundary clip
        double resolution = 1E-4; ///< clipper polygon resolution.
        double keepHexagonRate = 0.1; ///<minimum area ratio of the hexagonal grid to retain.
        double keepHexagonArea = 2.0; ///<minimum area of the hexagonal grid to retain.
        trimesh::vec2 arrayDir = trimesh::vec2(1, 0); ///<default the hexagonal grids array orientation.

        // for circular holes
        bool holeConnect = true; ///< is holes connected.
        float cheight = 5.0f; ///< the lowest point height of the first layer holes. 
        float ratio = 0.5f; ///< the portion of holes diameter to the edge length of the hexagonal grids.
        float delta = 0.5f; ///< the distance between the two separated circular holes.
        int nslices = 17;   ///< number of edges of the positive polygon fitting the circular hole.
        
        //debug
        bool bKeepHexagon = false; ///< Hexagons is to keep or not.
        int step_return = 9999; ///< debug quick return
    };

	class HoneyCombDebugger
	{
	public:
		virtual void onGenerateBottomPolygons(const TriPolygons& polygons) = 0;
		virtual void onGenerateInfillPolygons(const TriPolygons& polygons) = 0;
	};

    trimesh::TriMesh* generateHoneyCombs(trimesh::TriMesh* trimesh = nullptr, const HoneyCombParam& honeyparams = HoneyCombParam(),
        ccglobal::Tracer* tracer = nullptr, HoneyCombDebugger* debugger = nullptr);


	std::shared_ptr<trimesh::TriMesh> GenerateHoneyCombs(trimesh::TriMesh* trimesh = nullptr, const HoneyCombParam& honeyparams = HoneyCombParam(),
		ccglobal::Tracer* tracer = nullptr, HoneyCombDebugger* debugger = nullptr);

    struct HexagonArrayParam {
		trimesh::vec dir=trimesh::vec3(0,0,1);
        trimesh::vec3 pos; ///<
        int nrows = 2;
        int ncols = 2;
        double radius = 2.0; ///<
        double nestWidth = 0.3; ///<
    };
    struct honeyLetterOpt {
        double side = 1.0;
        std::vector<int>bottom; ///<
        std::vector<int>others; ///<
        /*
        hexagon: 
        radius: 
        borders: 
        neighbors: 
        */
        std::vector<HexaPolygon>hexgons;
    };
    TriPolygons GetOpenMeshBoundarys(const trimesh::TriMesh& triMesh, HoneyCombDebugger* debugger = nullptr);
    void getTriMeshBoundarys(trimesh::TriMesh& trimesh, std::vector<std::vector<int>>& sequentials);
    void GenerateBottomHexagons(const CMesh& honeyMesh, const HoneyCombParam& honeyparams, honeyLetterOpt& letterOpts, HoneyCombDebugger* debugger = nullptr);
    void GenerateTriPolygonsHexagons(const TriPolygons& polys, const HoneyCombParam& honeyparams, honeyLetterOpt& letterOpts, HoneyCombDebugger* debugger = nullptr);
    HexaPolygons GenerateHexagonsGridArray(const HexagonArrayParam& hexagonparams = HexagonArrayParam());
    struct ColumnarHoleParam {
        int nslices = 17; ///<
        float height = 5.0f; ///<
        float ratio = 0.5f; ///<
        float delta = 1.0f; ///<
        bool holeConnect = true;
    };
    void GenerateHexagonNeighbors(HexaPolygons& hexas, const ColumnarHoleParam& param = ColumnarHoleParam());
	trimesh::vec3 adjustHoneyCombParam(trimesh::TriMesh* trimesh,const HoneyCombParam& honeyparams);
    TriPolygons traitCurrentPolygons(const HexaPolygons& hexas, int index);
    TriPolygons traitNeighborPolygons(const HexaPolygons& hexas, int index);
    TriPolygons traitDirctionPolygon(const HexaPolygons& hexas, int index, int dir);
    TriPolygon traitPlanarCircle(const trimesh::vec3& c, float r, std::vector<int>& indexs, const trimesh::vec3& edgeDir = trimesh::vec3(0, 0, 1), int nums = 17);
    std::shared_ptr<trimesh::TriMesh> generateHolesColumnar(HexaPolygons& hexas, const ColumnarHoleParam& param);
    HexaPolygons generateEmbedHolesColumnar(trimesh::TriMesh* trimesh = nullptr, const HoneyCombParam& honeyparams = HoneyCombParam(),
        ccglobal::Tracer* tracer = nullptr, HoneyCombDebugger* debugger = nullptr);

	class MMeshT;
	void findNeighVertex(trimesh::TriMesh* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, int>>& vertex_distance);
	void findNeighVertex(MMeshT* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, float>>& vertex_distance);
	void innerHex(MMeshT* mesh, std::vector<std::vector<trimesh::vec2>>& poly, std::vector<int>& inFace, std::vector<int>& outFace,float len);

    trimesh::TriMesh* findOutlineOfDir(trimesh::TriMesh* mesh,std::vector<int>& botfaces);   
    void JointBotMesh(trimesh::TriMesh* mesh, trimesh::TriMesh* newmesh, std::vector<int>& botfaces);
    void SelectInnerFaces(trimesh::TriMesh* mesh,int indicate,std::vector<int>& out);
    void SelectBorderFaces(trimesh::TriMesh* mesh, int indicate, std::vector<int>& out);
}