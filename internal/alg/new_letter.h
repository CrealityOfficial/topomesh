#pragma once
#include "trimesh2/TriMesh.h"
#include "trimesh2/XForm.h"
#include "topomesh/interface.h"



namespace topomesh {





	//生成字体mesh接口 参数1.为字体轮廓、2.字体高度、3.返回每个字的坐标、4.返回每个字的点区间、5.返回每个字的BBX
	TOPOMESH_API trimesh::TriMesh* CreateFontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter, float height,
		std::vector<float>& word_location = std::vector<float>(), std::vector<int>& mesh_vertex_sizes = std::vector<int>(),
		std::vector<trimesh::vec3>& word_mesh_center=std::vector<trimesh::vec3>());

	//根据交互实现字体mesh的变换 参数1.目标mesh、2.字体mesh、3.鼠标点击faceid、4.鼠标点击的坐标、5.面的法线、
	//6.7.8 上一个接口的字体mesh的三个信息 9、暂定字体上方向 10.是否环绕
	TOPOMESH_API void MeshTransform(trimesh::TriMesh* traget_meshes, trimesh::TriMesh* font_mesh,int face_id,
		trimesh::vec3 location,trimesh::vec3 dir,std::vector<float>& word_location,std::vector<int>& mesh_vertex_sizes, std::vector<trimesh::vec3>& word_mesh_center,
		trimesh::vec3 up=trimesh::vec3(0,1,0), bool is_surround=false,float angle=0.f);
}