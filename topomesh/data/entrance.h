#pragma once
#include "mmesht.h"

namespace topomesh {

	class BaseParameterClass {};

	class InternelData
	{
	public:
		InternelData() {};
		InternelData(trimesh::TriMesh* trimesh) :_mesh(trimesh){};
		~InternelData() { _mesh.clear(); };
	private:
		MMeshT _mesh;

	public:
		void chunkedMesh(int n); //�ֿ�
		void getChunkFace(int ni, std::vector<int>& faceindexs);
		void getVertex(const std::vector<int>& faceindexs, std::vector<std::tuple<trimesh::ivec3, trimesh::point, trimesh::point, trimesh::point>>& vertexindexs);
		trimesh::TriMesh* mmesht2trimesh(bool save = false); //����trimesh ������ٶ�mesh���� ��saveΪĬ��ֵfalse ��������ں������� save��Ϊtrue
		void loopSubdivsion(std::vector<int>& faceindexs, std::vector<std::tuple<int, trimesh::point>>& vertex,
			std::vector<std::tuple<int, trimesh::ivec3>>& face_vertex, bool is_move=false,int iteration=1);
		void SimpleSubdivsion(std::vector<int>& faceindexs);
		void SimpleRemeshing(const std::vector<int>& faceindexs, float thershold);
	};
}