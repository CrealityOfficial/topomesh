#ifndef TOPOMESH_HONEYCOMB_1692753685514_H
#define TOPOMESH_HONEYCOMB_1692753685514_H
#include "topomesh/interface/idata.h"

namespace topomesh
{
	struct CombParam
	{
		float width = 1.0f;
		float radius = 1.0f;
		float combShell = 1.0f;

		bool holeConnect = false;
		float holeHeight = 1.0f;
		float holeRadius = 1.0f;
		float holeGap = 1.0f;
	};

	TOPOMESH_API std::shared_ptr<trimesh::TriMesh> honeyCombGenerate(trimesh::TriMesh* trimesh, const CombParam& honeyparams = CombParam(),
		ccglobal::Tracer* tracer = nullptr);

	/// <summary>
	/// 0  
	/// 1
	/// </summary>
	/// <param name="param"></param>
	/// <returns></returns>
	TOPOMESH_API int checkParam(const CombParam& param);
}

#endif // TOPOMESH_HONEYCOMB_1692753685514_H