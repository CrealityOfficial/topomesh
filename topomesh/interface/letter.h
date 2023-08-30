#ifndef TOPOMESH_LETTER_1692613164094_H
#define TOPOMESH_LETTER_1692613164094_H
#include "topomesh/interface/idata.h"
#include "cxnd/serial/serial.h"

namespace topomesh
{
	struct LetterParam
	{
		bool concave = true;
		float deep = 2.0f;

		//debug
		bool cacheInput = false;
		std::wstring fileName;
	};

	class LetterDebugger 
	{
	public:
		virtual ~LetterDebugger() {}

		virtual void onMeshProjected(const std::vector<trimesh::vec3>& triangles) = 0;
	};


	class TOPOMESH_API LetterInput : public cxnd::Serializeable
	{
	public:
		trimesh::TriMesh mesh;
		LetterParam param;
		SimpleCamera camera;
		std::vector<TriPolygons> polys;

		LetterInput();
		~LetterInput();

		trimesh::TriMesh* letter(LetterDebugger* debugger = nullptr, ccglobal::Tracer* tracer = nullptr);
	protected:
		int version() override;
		bool save(std::fstream& out, ccglobal::Tracer* tracer) override;
		bool load(std::fstream& in, int ver, ccglobal::Tracer* tracer) override;
	};

	TOPOMESH_API trimesh::TriMesh* letter(trimesh::TriMesh* mesh, const SimpleCamera& camera, 
		const LetterParam& param, const std::vector<TriPolygons>& polygons,
		LetterDebugger* debugger = nullptr, ccglobal::Tracer* tracer = nullptr);
}

#endif // TOPOMESH_LETTER_1692613164094_H