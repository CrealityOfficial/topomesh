#include "letter.h"
#include "topomesh/alg/letter.h"
#include "cxnd/serial/trimeshserial.h"
#include "ccglobal/profile.h"

namespace topomesh
{
	LetterInput::LetterInput()
	{

	}

	LetterInput::~LetterInput()
	{

	}

	int LetterInput::version()
	{
		return 0;
	}

	bool LetterInput::save(boost::nowide::fstream& out, ccglobal::Tracer* tracer)
	{
		cxnd::cxndSaveT(out, param.concave);
		cxnd::cxndSaveT(out, param.deep);
		cxnd::cxndSaveT(out, camera);
		cxnd::saveTrimesh(out, mesh);

		int size = (int)polys.size();
		cxnd::cxndSaveT(out, size);
		for (int i = 0; i < size; ++i)
			cxnd::savePolys(out, polys.at(i));

		return true;
	}

	bool LetterInput::load(boost::nowide::fstream& in, int ver, ccglobal::Tracer* tracer)
	{
		if (ver == 0)
		{
			cxnd::cxndLoadT(in, param.concave);
			cxnd::cxndLoadT(in, param.deep);
			cxnd::cxndLoadT(in, camera);
			cxnd::loadTrimesh(in, mesh);

			int size = 0;
			cxnd::cxndLoadT(in, size);
			if (size > 0)
			{
				polys.resize(size);
				for (int i = 0; i < size; ++i)
					cxnd::loadPolys(in, polys.at(i));
			}
			return true;
		}
		return false;
	}

	trimesh::TriMesh* LetterInput::letter(LetterDebugger* debugger, ccglobal::Tracer* tracer)
	{
		param.cacheInput = false;

		return topomesh::letter(&mesh, camera, param, polys, debugger, tracer);
	}

	trimesh::TriMesh* letter(trimesh::TriMesh* mesh, const SimpleCamera& camera,
		const LetterParam& param, const std::vector<TriPolygons>& polygons,
		LetterDebugger* debugger, ccglobal::Tracer* tracer) 
	{
		if (param.cacheInput && mesh)
		{
			LetterInput input;
			input.mesh = *mesh;
			input.camera = camera;
			input.param = param;
			input.polys = polygons;

			cxnd::cxndSave(input, param.fileName);
		}

		bool letterOpState = true;
		SYSTEM_TICK("letter");
		trimesh::TriMesh* result = letter(mesh, camera, param, polygons, letterOpState, debugger, tracer);
		SYSTEM_TICK("letter");
		if (letterOpState && result)
		{
			delete result;
			result = nullptr;
		}

		return result;
	}
}