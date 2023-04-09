#ifndef TOPOMESH_MMESHT_1680853426715_H
#define TOPOMESH_MMESHT_1680853426715_H
#include "trimesh2/Vec.h"
#include "mmeshFace.h"
#include "mmeshVertex.h"

namespace topomesh
{
	class MMeshT
	{
	public:
		MMeshT() :VVadjacent(true), VFadjacent(true), FFadjacent(true) { faces.reserve(1024); vertices.reserve(3000); };//单个加入限制在1000以内
		MMeshT(const MMeshT& mt) {};
		MMeshT(MMeshT& mt) {};
		MMeshT(trimesh::TriMesh* currentMesh);
		MMeshT& operator=(MMeshT mt) {};
		virtual ~MMeshT() { vertices.clear(); faces.clear(); };

		std::vector<MMeshVertex> vertices;
		std::vector<MMeshFace> faces;

	public:
		static inline double det(trimesh::point& p0, trimesh::point& p1, trimesh::point& p2);
		inline double det(int faceIndex);
		inline double det(int VertexIndex1, int VertexIndex2, int VertexIndex3);		
		void mmesh2trimesh(trimesh::TriMesh* currentMesh);

		void appendVertex(MMeshVertex& v);
		void appendVertex(trimesh::point& v);

		//inline void appendFace(MMeshFace& f);
		void appendFace(MMeshVertex& v0, MMeshVertex& v1, MMeshVertex& v2);
		void appendFace(int i0, int i1, int i2);

		void deleteFace(MMeshFace& f);
		void deleteFace(int i);

		void deleteVertex(MMeshVertex& v);
		void deleteVertex(int i);

		void shrinkMesh();
		void getMeshBoundary();
		void getEdge(std::vector<trimesh::ivec2>& edge, bool is_select = false);
		void calculateCrossPoint(std::vector<trimesh::ivec2>& edge, std::vector<trimesh::point>& line, std::vector<trimesh::vec3>& tc);
	

		inline bool is_VVadjacent() { return VVadjacent; }
		inline bool is_VFadjacent() { return VFadjacent; }
		inline bool is_FFadjacent() { return FFadjacent; }

		inline int VN() const { return vn; }
		inline int FN() const { return fn; }

	private:
		int vn = 0;
		int fn = 0;
		bool VVadjacent = false;
		bool VFadjacent = false;
		bool FFadjacent = false;
	};

	double  getTotalArea(std::vector<trimesh::point>& inVertices);
}

#endif // TOPOMESH_MMESHT_1680853426715_H