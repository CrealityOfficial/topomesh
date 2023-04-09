#pragma once
#include "trimesh2/Vec.h"
#include "trimesh2/TriMesh.h"
#include "mmeshVertex.h"


namespace topomesh
{
	//class MMeshVertex;

	class MMeshFace
	{
	public:
		MMeshFace() {};
		MMeshFace(trimesh::TriMesh::Face f) 
		{
			for (int i = 0; i < 3; i++)
			{
				vertex_index[i] = f.at(i);
			}
		};
		MMeshFace(MMeshVertex* v0, MMeshVertex* v1, MMeshVertex* v2)
		{		
			connect_vertex.push_back(v0); connect_vertex.push_back(v1); connect_vertex.push_back(v2);			
		}
		virtual ~MMeshFace() {};	
		int vertex_index[3] = { -1 };
		int connected_face_index[3];		
		int index = -1;
		bool operator==(MMeshFace* f) const { return this == f ? true : false; };

		trimesh::point normal;
		std::vector<MMeshVertex*> connect_vertex;
		std::vector<MMeshFace*> connect_face;

		std::vector<trimesh::vec2> uv_coord;
	private:
		enum faceflag
		{
			MF_DELETE = 0x00000001,
			MF_SELECT = 0x00000002,
			MF_BORDER = 0x00000004,
			MF_VISITED= 0x00000008,
			MF_ACCUMULATE   = 0x00000ff0
		};
		int flag=0;
		
	public:
		inline void SetD() { flag |= MF_DELETE; }
		inline bool IsD() { return (MF_DELETE & flag)!=0 ? 1 : 0; }
		inline void ClearD() { flag &= ~MF_DELETE; }

		inline void SetS() { flag |= MF_SELECT; }
		inline bool IsS() { return (MF_SELECT & flag) != 0 ? 1 : 0; }
		inline void ClearS() { flag &= ~MF_SELECT; }

		inline void SetB() { flag |= MF_BORDER; }
		inline bool IsB() { return (MF_BORDER & flag) != 0 ? 1 : 0; }
		inline void ClearB() { flag &= ~MF_BORDER; }

		inline void SetV() { flag |= MF_VISITED; }
		inline bool IsV() { return (MF_VISITED & flag) != 0 ? 1 : 0; }
		inline void ClearV() { flag &= ~MF_VISITED; }

		//0 - 255
		inline bool SetA() { int copy = flag & MF_ACCUMULATE; copy = copy >> 4; copy += 1; if (copy > 255) return false; copy = copy << 4;flag&=~MF_ACCUMULATE ; flag |= copy; return true; }
		inline bool IsA(int i) { int copy = flag & MF_ACCUMULATE; copy = copy >> 4; if (i == copy) return true; return false; }
		inline void ClearA() { flag &= ~MF_ACCUMULATE; }

		inline void open_uv() { uv_coord.reserve(8); }

		inline MMeshVertex* V0(int i) { return connect_vertex[(i + 0) % 3]; }
		inline MMeshVertex* V1(int i) { return connect_vertex[(i + 1) % 3]; }
		inline MMeshVertex* V2(int i) { return connect_vertex[(i + 2) % 3]; }
		inline MMeshVertex* V0(MMeshVertex* v){
			for (int i = 0; i < 3; i++)
			{
				if (v == this->connect_vertex[i])
					return this->connect_vertex[i];
			}
			return nullptr;
		}
		inline MMeshVertex* V1(MMeshVertex* v){ 
			for (int i = 0; i < 3; i++)
			{
				if (v == this->connect_vertex[i])
					return this->connect_vertex[(i+1)%3];
			}
			return nullptr;
		}
		inline MMeshVertex* V2(MMeshVertex* v) { 
			for (int i = 0; i < 3; i++)
			{
				if (v == this->connect_vertex[i])
					return this->connect_vertex[(i + 2) % 3];
			}
			return nullptr;
		}

		static inline bool antiClockWise(trimesh::point v0, trimesh::point v1, trimesh::point v2) {
			float f = v0.x * (v1.y * v2.z - v1.z * v2.y) + v0.y * (v1.z * v2.x - v1.x * v2.z) + v0.z * (v1.x * v2.y - v1.y * v2.x);
			if (f < 0)
				return true;
			return false;
		}

		static inline bool antiClockWise(MMeshVertex* v0, MMeshVertex* v1, MMeshVertex* v2) {
			return antiClockWise(v0->p, v1->p, v2->p);
		}
	};
}