#pragma once
#include "topomesh/data/convert.h"
#include "ccglobal/tracer.h"
#include "topomesh/interface/idata.h"
#include "topomesh/data/CMesh.h"
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
        trimesh::vec3 pos; ///<�������Ͻ���ʼ��
        int nrows = 2;
        int ncols = 2;
        double radius = 2.0; ///<δ����ǰ���������α߳�
        double nestWidth = 0.3; ///<���������αں����������ľ����2����
    };
    HexaPolygons GenerateHexagonsGridArray(const HexagonArrayParam& hexagonparams = HexagonArrayParam());
    struct ColumnarHoleParam {
        int nslices = 17; ///<Բ��Ĭ����17����
        float height = 5.0f; ///<Բ��Բ�ĸ߶�
        float ratio = 0.5f; ///<Բ��ֱ���������߳��ı���
        float delta = 1.0f; ///<������Բ��Բ�ľ���
    };
    void GenerateHexagonNeighbors(HexaPolygons& hexas, const ColumnarHoleParam& param = ColumnarHoleParam());
	trimesh::vec3 adjustHoneyCombParam(trimesh::TriMesh* trimesh,const HoneyCombParam& honeyparams);
    TriPolygons traitCurrentPolygons(const HexaPolygons& hexas, int index);
    TriPolygons traitNeighborPolygons(const HexaPolygons& hexas, int index);
    TriPolygons traitDirctionPolygon(const HexaPolygons& hexas, int index, int dir);
    
    TriPolygon traitPlanarCircle(const trimesh::vec3& c, float r, std::vector<int>& indexs, const trimesh::vec3& edgeDir = trimesh::vec3(0, 0, 1), int nums = 17);
    std::shared_ptr<trimesh::TriMesh> generateHolesColumnar(HexaPolygons& hexas, const ColumnarHoleParam& param);

	class MMeshT;
	void findNeighVertex(trimesh::TriMesh* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, int>>& vertex_distance);
	void findNeighVertex(MMeshT* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, float>>& vertex_distance);
	void innerHex(MMeshT* mesh, std::vector<std::vector<trimesh::vec2>>& poly, std::vector<int>& inFace, std::vector<int>& outFace,float len);

    trimesh::TriMesh* findOutlineOfDir(trimesh::TriMesh* mesh);   

    struct honeyLetterOpt {
        std::vector<int>bottom; ///<������ƽ�����Ƭ����
        std::vector<int>others; ///<ȥ�������������Ƭ�������ѱ���ԭģ�Ͷ�Ӧ��������
        /*
        hexagon: ÿ����������ṹ��;
        radius: ��������ǰ������߳�;
        borders: ��������������(��һ��6���㣬Ĭ��xoy����ϵ��ʱ������);
        neighbors: �����������ڵ�������������(6������һ���������ڣ�û�м�-1);
        */
        struct hexagon {
            double radius = 1.0;
            std::vector<trimesh::vec3>borders;
            std::vector<int> neighbors;
            hexagon() : neighbors(6, -1) {}
        };
        //�������������ھӹ�ϵ
        std::vector<hexagon>hexgons;
    };
   
}