#ifndef TOPOMESH_HONEYCOMB_1692753685514_H
#define TOPOMESH_HONEYCOMB_1692753685514_H
#include "topomesh/interface/idata.h"

namespace topomesh
{
	struct CombParam
	{
		float width = 2.0f;
		float diameter = 6.0f;
		float combShell = 2.0f;

		bool holeConnect = true;
		float holeHeight = 1.0f;
		float holeDiameter = 2.5f;
		float holeGap = 1.0f;
	};

	TOPOMESH_API std::shared_ptr<trimesh::TriMesh> honeyCombGenerate(trimesh::TriMesh* trimesh, const CombParam& honeyparams = CombParam(),
		ccglobal::Tracer* tracer = nullptr);

	/// <summary>
	/// 0: parameters qualified.  
	/// 1: combShell is not qualified.
    /// 2: diameter is not qualified.
    /// 3: width is not qualified.
    /// 4: holeDiameter is not qualified.
    /// 5: holeGap is not qualified.
    /// 6: holeHeight is not qualified.
	/// </summary>
	/// <param name="param"></param>
	/// <returns></returns>
	TOPOMESH_API int checkParam(const CombParam& param);
}

#endif // TOPOMESH_HONEYCOMB_1692753685514_H