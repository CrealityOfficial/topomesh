#pragma once
#include "trimesh2/TriMesh.h"
#include "trimesh2/XForm.h"
#include "topomesh/interface.h"



namespace topomesh {





	//��������mesh�ӿ� ����1.Ϊ����������2.����߶ȡ�3.����ÿ���ֵ����ꡢ4.����ÿ���ֵĵ����䡢5.����ÿ���ֵ�BBX
	TOPOMESH_API trimesh::TriMesh* CreateFontMesh(const std::vector<std::vector<std::vector<trimesh::vec2>>>& letter, float height,
		std::vector<float>& word_location = std::vector<float>(), std::vector<int>& mesh_vertex_sizes = std::vector<int>(),
		std::vector<trimesh::vec3>& word_mesh_center=std::vector<trimesh::vec3>());

	//���ݽ���ʵ������mesh�ı任 ����1.Ŀ��mesh��2.����mesh��3.�����faceid��4.����������ꡢ5.��ķ��ߡ�
	//6.7.8 ��һ���ӿڵ�����mesh��������Ϣ 9���ݶ������Ϸ��� 10.�Ƿ���
	TOPOMESH_API void MeshTransform(trimesh::TriMesh* traget_meshes, trimesh::TriMesh* font_mesh,int face_id,
		trimesh::vec3 location,trimesh::vec3 dir,std::vector<float>& word_location,std::vector<int>& mesh_vertex_sizes, std::vector<trimesh::vec3>& word_mesh_center,
		trimesh::vec3 up=trimesh::vec3(0,1,0), bool is_surround=false,float angle=0.f);
}