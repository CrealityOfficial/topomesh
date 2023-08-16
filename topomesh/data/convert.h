#pragma once
#include "trimesh2/TriMesh.h"
#include <memory>

namespace topomesh
{
	struct LetterParam
	{		
		bool concave = true;
		float deep = 2.0f;
	};

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

	class LetterInput
	{
	public:
		LetterParam letterparam;

		void save(const std::string& fileName);
		void load(const std::string& fileName);
	};

	typedef std::vector<trimesh::vec3> TriPolygon;
	typedef std::vector<TriPolygon> TriPolygons;
    struct HexaPolygon {
        TriPolygon poly;
        trimesh::ivec3 coord;
        std::vector<int> neighbors;
        HexaPolygon() : neighbors(6, -1) {}
    };
    typedef std::vector<HexaPolygon> HexaPolygons;
    TriPolygons convertFromHexaPolygons(const HexaPolygons& hexas);
    HexaPolygons convertFromTriPolygons(const TriPolygons& polys, const trimesh::ivec3& coords);

	struct SimpleCamera
	{
		float f;
		float n;
		float fov;
		float aspect;

		trimesh::point pos;
		trimesh::point center;
		trimesh::point up;
	};

	class LetterDebugger
	{
	public:
		virtual ~LetterDebugger() {}

		virtual void onMeshProjected(const std::vector<trimesh::vec3>& triangles) = 0;
	};

	typedef std::vector<int> FacePatch;
	typedef std::vector<FacePatch> FacePatchs;

	typedef std::shared_ptr<trimesh::TriMesh> TopoTriMeshPtr;
}