#ifndef TOPOMESH_IDATA_1692613164081_H
#define TOPOMESH_IDATA_1692613164081_H
#include "topomesh/interface.h"
#include "trimesh2/TriMesh.h"
#include "trimesh2/XForm.h"
#include "ccglobal/tracer.h"
#include <memory>

namespace topomesh
{
	typedef std::vector<trimesh::vec3> TriPolygon;
	typedef std::vector<TriPolygon> TriPolygons;

	struct SimpleCamera
	{
		float f;
		float n;
		float fov;
		float aspect;

		trimesh::point pos;
		trimesh::point center;
		trimesh::point up;
	};

	typedef std::vector<int> FacePatch;
	typedef std::vector<FacePatch> FacePatchs;

	typedef std::shared_ptr<trimesh::TriMesh> TopoTriMeshPtr;

	struct HoneyCombParam
	{
		double resolution = 1E-4; ///<多边形相交允许误差
		TriPolygon* polyline = nullptr; ///<自由编辑区域生成内蜂窝（默认底面全部填充）
		trimesh::vec3 axisDir = trimesh::vec3(0, 0, 0); ///<蜂窝多边形平面朝向（默认z轴正方向）
		//trimesh::vec3* axisDir = nullptr;
		trimesh::vec2 arrayDir = trimesh::vec2(1, 0); ///<蜂窝多边形平面布局（默认边朝上结构）
		double honeyCombRadius = 2.0; ///< 生成的蜂窝六边形边长
		double nestWidth = 1.0; ///<蜂窝六边形壁厚（向内收缩的距离的2倍）
		double shellThickness = 0.6; ///<抽壳厚度
		double keepHexagonRate = 0.1; ///允许保留网格的最小面积比例
        double keepHexagonArea = 2.0; ///允许保留网格的最小面积绝对值参考量
		bool isdelect = false; //是否删除标记区域
		std::vector<int> faces; //所选择的区域面，如果为默认-1则为自定义方向
        bool holeConnect = true;
        float cheight = 5.0f; ///<圆孔圆心高度
		float ratio = 0.5f; ///<圆孔直径相对网格边长比例
		float delta = 0.5f; //圆孔间隔
		int nslices = 17;   //拟合多边形边数
        bool bKeepHexagon = false;
		//debug
		int step_return = 9999; // debug quick return
	};
}

#endif // TOPOMESH_IDATA_1692613164081_H