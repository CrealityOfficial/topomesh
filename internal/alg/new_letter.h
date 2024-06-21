#pragma once
#include "trimesh2/TriMesh.h"
#include "trimesh2/XForm.h"


namespace topomesh {
	trimesh::TriMesh* CreateFontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter, float height,
		std::vector<float>& word_location = std::vector<float>(), std::vector<int>& mesh_vertex_sizes = std::vector<int>(),
		std::vector<trimesh::vec3>& word_mesh_center=std::vector<trimesh::vec3>());

	void MeshTransform(trimesh::TriMesh* traget_meshes,trimesh::xform& font_xform, trimesh::TriMesh* font_mesh,int face_id,
		trimesh::vec3 location,trimesh::vec3 dir,std::vector<float>& word_location,std::vector<int>& mesh_vertex_sizes, std::vector<trimesh::vec3>& word_mesh_center,
		trimesh::vec3 up=trimesh::vec3(0,0,1), bool is_surround=false);
}