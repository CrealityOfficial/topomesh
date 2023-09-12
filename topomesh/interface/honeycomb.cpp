#include "honeycomb.h"
#include "topomesh/alg/fillhoneycombs.h"
namespace topomesh
{
	std::shared_ptr<trimesh::TriMesh> honeyCombGenerate(trimesh::TriMesh* trimesh, const CombParam& honeyparams,
		ccglobal::Tracer* tracer)
	{
        HoneyCombParam params;
        params.honeyCombRadius = honeyparams.diameter / 2.0;
        params.nestWidth = honeyparams.width;
        params.shellThickness = honeyparams.combShell;
        params.cheight = honeyparams.holeHeight;
        params.delta = honeyparams.holeGap;
        params.ratio = honeyparams.holeDiameter / honeyparams.diameter;
        params.holeConnect = honeyparams.holeConnect;
        std::shared_ptr<trimesh::TriMesh> mesh = GenerateHoneyCombs(trimesh, params, tracer);
		return mesh;
	}

    int checkParam(const CombParam& param)
    {
        if (param.combShell < 0.5f || param.combShell > 10.f) {
            return 1;
        }else if (param.diameter < 1.0f || param.diameter > 10.0f) {
            return 2;
        } else if (param.width < 1.0f || param.width > param.diameter) {
            return 3;
        } else if (param.holeDiameter < 1.0f || param.holeDiameter > param.diameter / 2.0f - 0.5f) {
            return 4;
        } else if (param.holeGap < 1.0f || param.holeGap > 10.0f) {
            return 5;
        } else if (param.holeHeight < 1.0f || param.holeHeight > 3.0f) {
            return 6;
        }
        return 0;
    }
}