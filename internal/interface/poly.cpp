#include "topomesh/interface/poly.h"
#include "internal/polygon/conv1.h"

namespace topomesh
{
	void simplifyPolygons(TriPolygons& polys)
	{
		ClipperLib::Paths paths;
		convertRaw(polys, paths);

		ClipperLib::SimplifyPolygons(paths);
		convertRaw(paths, polys);
	}
}