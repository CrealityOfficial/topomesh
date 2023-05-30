#include "topomesh/data/mmesht.h"
#include "topomesh/alg/convert.h"
#include "topomesh/alg/letter.h"
#include "trimesh2/XForm.h"
#include "trimesh2/quaternion.h"
#include "Eigen/Dense"

#include "ccglobal/tracer.h"

namespace topomesh {

	bool ModleCutting(const std::vector<trimesh::TriMesh*>& inMesh, std::vector<trimesh::TriMesh*>& outMesh, const SimpleCamera& camera,
		const TriPolygon& paths, ccglobal::Tracer* tracer = nullptr);
	bool JudgeCloseOfPath(const TriPolygon& paths);
	bool JudgeMeshIsVaild(const trimesh::TriMesh* inMesh);
	void setMarkOfCorssPoint(std::vector<MMeshT*>& meshs, const TriPolygon& paths,const std::vector<std::pair<int,int>>& corssPoint );
	void splitMesh(MMeshT* mesh, std::vector<MMeshT>& outmesh);
}