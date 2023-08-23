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
		double resolution = 1E-4; ///<������ཻ�������
		TriPolygon* polyline = nullptr; ///<���ɱ༭���������ڷ��ѣ�Ĭ�ϵ���ȫ����䣩
		trimesh::vec3 axisDir = trimesh::vec3(0, 0, 1); ///<���Ѷ����ƽ�泯��Ĭ��z��������
		trimesh::vec2 arrayDir = trimesh::vec2(1, 0); ///<���Ѷ����ƽ�沼�֣�Ĭ�ϱ߳��Ͻṹ��
		double honeyCombRadius = 1.0; ///<δ����ǰ���������α߳�
		double nestWidth = 0.1; ///<���������αں����������ľ����2����
		double shellThickness = 0.1; ///<��Ǻ��
		double keepHexagonRate = 0.1; ///���������������С�������
		bool isdelect = false; //�Ƿ�ɾ���������
		std::vector<int> faces; //��ѡ��������棬���ΪĬ��-1��Ϊ�Զ��巽��
		float cheight = 5.0f; ///<Բ��Բ�ĸ߶�
		float ratio = 0.5f; ///<Բ�װ뾶
		float delta = 0.5f; //Բ�׼��
		int nslices = 17;   //��϶���α���
		//debug
		int step_return = 9999; // debug quick return
	};
}

#endif // TOPOMESH_IDATA_1692613164081_H