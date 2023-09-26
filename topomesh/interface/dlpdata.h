#ifndef TOPOMESH_CXSW_DLPDATA_1593762618888_H
#define TOPOMESH_CXSW_DLPDATA_1593762618888_H
#include "topomesh/interface.h"
#include "trimesh2/Vec.h"
#include <vector>

namespace topomesh
{
	typedef std::vector<trimesh::dvec2> DLPPoly;
	typedef std::vector<DLPPoly> DLPPolys;

	typedef std::vector<std::vector<trimesh::ivec2>> UmDLPPolys;

	struct LayerDetailInfo
	{
		double maxArea = 0.0;    //um * 2
		double totalArea = 0.0;  //um * 2
		double maxDistances = 0.0;  //um
	};

	struct SliceInfo
	{
		double volume = 0.0;   // mm * 3
		std::vector<LayerDetailInfo> details;
	};

	class DLPDataImpl;
	class TOPOMESH_API DLPData
	{
		friend class DLPSlicer;
	public:
		DLPData();
		virtual ~DLPData();

		int layers() const;
		bool isValid() const;
		void traitPolys(int layer, DLPPolys& polys) const;
		void traitUmPolys(int layer, UmDLPPolys& polys) const;
		std::vector<std::vector<trimesh::vec2>> traitPolys(int layer) const;

		void calculateVolumeAreas(SliceInfo& info) const;
		//ClipperLib::PolyTree* layerData(int layer) const;

	protected:
		DLPDataImpl* impl;
	};
}


#endif // TOPOMESH_CXSW_DLPDATA_1593762618888_H
