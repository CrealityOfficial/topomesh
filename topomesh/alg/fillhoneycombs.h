#pragma once
#include "topomesh/data/convert.h"
#include "ccglobal/tracer.h"

namespace topomesh {

	struct HoneyCombParam
	{
        TriPolygon* polyline = nullptr;
        trimesh::vec3 axisDir = trimesh::vec3(0, 0, 1);
        trimesh::vec2 arrayDir = trimesh::vec2(1, 0);
		double honeyCombRadius = 1.0;
		double nestWidth = 0.3;
		double shellThickness = 1.0;

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




	class MMeshT;
	void GenerateHoneyCombs(const trimesh::TriMesh* mesh, trimesh::TriMesh& resultmesh, const TriPolygon& poly, trimesh::vec3 axisDir=trimesh::vec3(0,0,1),
		trimesh::vec2 arrayDir = trimesh::vec2(0,1), double honeyCombRadius=2.0f, double nestWidth=1.0f, double shellThickness=1.0f);
	void findNeighVertex(trimesh::TriMesh* mesh, const std::vector<int>& upfaceid, const std::vector<int>& botfaceid, std::vector<std::pair<int, int>>& vertex_distance);
	void innerHex(MMeshT* mesh, std::vector<std::vector<trimesh::vec2>>& poly, std::vector<int>& inFace, std::vector<int>& outFace,float len);
}