#pragma once
#include "trimesh2/TriMesh.h"
#include "topomesh/data/mmesht.h"
#include "clipperxyz/clipper.hpp"

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
		MMeshT mesh;
		CameraParam cameraparam;
		LetterParam letterparam;
		ClipperLibXYZ::Paths paths;

		void save(const std::string& fileName);
		void save(trimesh::TriMesh* mesh, LetterParam& letterparam, CameraParam& cameraparam, ClipperLibXYZ::Paths& paths);
		void load(const std::string& fileName);
	};

	typedef std::vector<trimesh::vec3> TriPolygon;
	typedef std::vector<TriPolygon> TriPolygons;

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
}