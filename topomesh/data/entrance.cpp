#include "entrance.h"
#include "topomesh/alg/subdivision.h"
#include "topomesh/alg/segmentation.h"
#include "topomesh/alg/remesh.h"
#include "unordered_map"

namespace topomesh {
	trimesh::TriMesh* InternelData::mmesht2trimesh(bool save)
	{
		trimesh::TriMesh* newmesh = new trimesh::TriMesh();
		if (!save)
			this->_mesh.quickTransform(newmesh);
		else
			this->_mesh.mmesh2trimesh(newmesh);
		return newmesh;
	}

	void InternelData::loopSubdivsion(std::vector<int>& faceindexs, std::vector<std::tuple<int, trimesh::point>>& vertex,
		std::vector<std::tuple<int, trimesh::ivec3>>& face_vertex, bool is_move, int iteration)
	{
		topomesh::loopSubdivision(&this->_mesh, faceindexs, vertex, face_vertex,is_move, iteration);
		return;
	}

	void InternelData::SimpleSubdivsion(std::vector<int>& faceindexs)
	{
		topomesh::SimpleMidSubdivision( &this->_mesh,faceindexs);
		return;
	}

	void InternelData::chunkedMesh(int n)
	{
		if (!this->_mesh.is_BoundingBox()) this->_mesh.getBoundingBox();
		float x = (this->_mesh.boundingbox.max_x - this->_mesh.boundingbox.min_x)*1.f / n * 1.f;
		float y = (this->_mesh.boundingbox.max_y - this->_mesh.boundingbox.min_y)*1.f / n * 1.f;
		float z = (this->_mesh.boundingbox.max_z - this->_mesh.boundingbox.min_z)*1.f / n * 1.f;
		int Chunknum = std::pow(n, 3);		
		std::vector<std::vector<int>> chunkface(Chunknum);
		for (MMeshFace& f : this->_mesh.faces)
		{
			if (f.GetU() > 0) continue;
			trimesh::point c = (f.V0(0)->p + f.V1(0)->p + f.V2(0)->p) / 3.f;
			int xi = (c.x - this->_mesh.boundingbox.min_x) / x;
			int yi = (c.y - this->_mesh.boundingbox.min_y) / y;
			int zi = (c.z - this->_mesh.boundingbox.min_z) / z;
			int chunked = zi * std::pow(n, 2) + yi * n + xi + 1;
			f.SetU(chunked);
			chunkface[chunked-1].push_back(f.index);		
		}

		for (std::vector<int>& vec : chunkface)
		{
			MMeshT mt(&this->_mesh, vec);
			this->_ChunkMesh.emplace_back(std::move(mt));
		}

	}

	void InternelData::getChunkFaces(int ni ,std::vector<int>& faceindexs)
	{
		int deleteFnum = 0;
		for (MMeshFace& f : this->_mesh.faces)
		{
			if (f.IsD())
			{
				deleteFnum++; continue;
			}
			if (f.GetU() == ni)
				faceindexs.push_back(f.index-deleteFnum);
		}
	}

	void InternelData::getVertex(const std::vector<int>& faceindexs, std::vector<std::tuple<trimesh::ivec3, trimesh::point, trimesh::point, trimesh::point>>& vertexindexs)
	{
		for (int i : faceindexs)
		{
			MMeshFace& f = this->_mesh.faces[i];
			vertexindexs.push_back(std::make_tuple(trimesh::ivec3(f.V0(0)->index, f.V1(0)->index, f.V2(0)->index),
				f.V0(0)->p, f.V1(0)->p, f.V2(0)->p));
		}
	}

	void InternelData::SimpleRemeshing(const std::vector<int>& faceindexs, float thershold)
	{
		topomesh::SimpleRemeshing(&this->_mesh, faceindexs, thershold);
	}

	int InternelData::getFaceChunk(int faceindex)
	{
		return this->_mesh.faces[faceindex].GetU();
	}

	void InternelData::ChunkMeshSimpleRemshing(const std::vector<int>& Chunkid, const std::vector<std::vector<int>>& ChunkMeshfaceindexs, float thershold)
	{
		for (int i=0 ; i<Chunkid.size();i++)
		{
			topomesh::SimpleRemeshing(&this->_ChunkMesh[Chunkid[i]],ChunkMeshfaceindexs[i],thershold);
		}
	}

	void InternelData::QuickCombinationMesh()
	{
		this->_mesh.clear();
		int Accface = 0;
		//std::unordered_map<unsigned long long, int> map;
		for (MMeshT& m : this->_ChunkMesh)
		{
			int deleteVnum = 0;		
			for (MMeshVertex& v : m.vertices)
			{
				if (v.IsD())
				{
					deleteVnum++; continue;
				}
				v.index -= deleteVnum;
				this->_mesh.appendVertex(v.p);
				//unsigned long long hash_n = (int)(v.p.x * 99971) ^ (int)(v.p.y * 99989) ^ (int)(v.p.z * 99991);
				//map[hash_n] = this->_mesh.vertices.size() - 1;
			}
			for (MMeshFace& f : m.faces)
			{
				if (!f.IsD())
				{					
					this->_mesh.appendFace(f.connect_vertex[0]->index+ Accface,f.connect_vertex[1]->index+ Accface,f.connect_vertex[2]->index+ Accface);
				}
			}
			Accface = this->_mesh.vertices.size();
		}
	}
}