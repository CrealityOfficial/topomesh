#include "topomesh/data/mmesht.h"
#include "topomesh/alg/convert.h"
#include "trimesh2/XForm.h"
#include "trimesh2/quaternion.h"
#include "Eigen/Dense"

#include "ccglobal/tracer.h"

namespace topomesh {

	bool ModleCutting(const std::vector<trimesh::TriMesh*>& inMesh, std::vector<trimesh::TriMesh*>& outMesh, const SimpleCamera& camera,
		const TriPolygon& paths, ccglobal::Tracer* tracer = nullptr);
}