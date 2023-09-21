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
        double resolution = 1E-4; ///<������ཻ�������
        TriPolygon* polyline = nullptr; ///<���ɱ༭���������ڷ��ѣ�Ĭ�ϵ���ȫ����䣩
        trimesh::vec3 axisDir = trimesh::vec3(0, 0, 0); ///<���Ѷ����ƽ�泯��Ĭ��z��������
        //trimesh::vec3* axisDir = nullptr;
        trimesh::vec2 arrayDir = trimesh::vec2(1, 0); ///<���Ѷ����ƽ�沼�֣�Ĭ�ϱ߳��Ͻṹ��
        double honeyCombRadius = 2.0; ///< ���ɵķ��������α߳�
        double nestWidth = 1.0; ///<���������αں����������ľ����2����
        double shellThickness = 0.6; ///<��Ǻ��
        double keepHexagonRate = 0.1; ///���������������С�������
        double keepHexagonArea = 2.0; ///���������������С�������ֵ�ο���
        bool isdelect = false; //�Ƿ�ɾ���������
        std::vector<int> faces; //��ѡ��������棬���ΪĬ��-1��Ϊ�Զ��巽��
        bool holeConnect = true;
        float cheight = 5.0f; ///<Բ��Բ�ĸ߶�
        float ratio = 0.5f; ///<Բ��ֱ���������߳�����
        float delta = 0.5f; //Բ�׼��
        int nslices = 17;   //��϶���α���
        bool bKeepHexagon = false;
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
   
}