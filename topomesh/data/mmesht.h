#ifndef TOPOMESH_MMESHT_1680853426715_H
#define TOPOMESH_MMESHT_1680853426715_H
#include "trimesh2/Vec.h"
#include "mmeshFace.h"
#include "mmeshVertex.h"
#include <map>

namespace topomesh
{
	class MMeshT
	{
	public:
		MMeshT() :VVadjacent(true), VFadjacent(true), FFadjacent(true) { faces.reserve(1024); vertices.reserve(3000); };//单个加入限制在1000以内
		MMeshT(const MMeshT& mt) {};
		MMeshT(MMeshT& mt) {};
		MMeshT(trimesh::TriMesh* currentMesh);
		MMeshT(trimesh::TriMesh* currentMesh,std::vector<int>& faces,std::map<int,int>& vmap, std::map<int, int>& fmap);
		MMeshT& operator=(MMeshT mt) {};
		virtual ~MMeshT() { vertices.clear(); faces.clear(); };

		std::vector<MMeshVertex> vertices;
		std::vector<MMeshFace> faces;

	public:
		static inline double det(trimesh::point& p0, trimesh::point& p1, trimesh::point& p2)
		{
			trimesh::vec3 a = p0 - p1;
			trimesh::vec3 b = p1 - p2;
			return sqrt(pow((a.y * b.z - a.z * b.y), 2) + pow((a.z * b.x - a.x * b.z), 2)
				+ pow((a.x * b.y - a.y * b.x), 2)) / 2.0f;
		}
		inline double det(int faceIndex)
		{
			return det(this->faces[faceIndex].V0(0)->index, this->faces[faceIndex].V0(1)->index, this->faces[faceIndex].V0(2)->index);
		}
		inline double det(int VertexIndex1, int VertexIndex2, int VertexIndex3)
		{
			return det(this->vertices[VertexIndex1].p, this->vertices[VertexIndex2].p, this->vertices[VertexIndex3].p);
		}
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
		void getMeshBoundaryFaces();
		void getEdge(std::vector<trimesh::ivec2>& edge, bool is_select = false);
		void calculateCrossPoint(std::vector<trimesh::ivec2>& edge, std::pair<trimesh::point,trimesh::point>& line, std::vector<std::pair<float,trimesh::ivec2>>& tc);
		void getFacesNormals();

		inline bool is_VVadjacent() { return VVadjacent; }
		inline bool is_VFadjacent() { return VFadjacent; }
		inline bool is_FFadjacent() { return FFadjacent; }
		inline bool is_VertexNormals() { return VertexNormals; }
		inline bool is_FaceNormals() { return FacesNormals; }

		inline void set_VVadjacent(bool b) { VVadjacent = b; }
		inline void set_VFadjacent(bool b) { VFadjacent = b; }
		inline void set_FFadjacent(bool b) { FFadjacent = b; }
		inline void set_VertexNormals(bool b) { VertexNormals = b; }
		inline void set_FacesNormals(bool b) { FacesNormals = b; }

		inline void clear() { vertices.clear(); faces.clear();
		}

		inline int VN() const { return vn; }
		inline int FN() const { return fn; }

	private:
		int vn = 0;
		int fn = 0;
		bool VVadjacent = false;
		bool VFadjacent = false;
		bool FFadjacent = false;
		bool VertexNormals = false;
		bool FacesNormals = false;
	};

	double  getTotalArea(std::vector<trimesh::point>& inVertices);
}

#endif // TOPOMESH_MMESHT_1680853426715_H