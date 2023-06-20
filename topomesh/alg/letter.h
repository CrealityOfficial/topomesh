#ifndef TOPOMESH_LETTER_1680853426716_H
#define TOPOMESH_LETTER_1680853426716_H
#include "mmesh/trimesh/polygon.h"
#include "mmesh/util/mnode.h"
#include "topomesh/data/mmesht.h"
#include "topomesh/alg/convert.h"
#include "convert.h"
#include "trimesh2/XForm.h"
#include "trimesh2/quaternion.h"
#include "Eigen/Dense"
#include "omp.h"
#include "cmath"
#include <vector>
#include <queue>
#include <numeric>

#include "ccglobal/tracer.h"

namespace topomesh
{	
	void lettering(MMeshT* mesh, const std::vector<ClipperLibXYZ::Paths>& paths,  CameraParam& camera, const LetterParam& Letter, std::vector<int>* faceindex = nullptr);
	void concaveOrConvexOfFaces(MMeshT* mt,std::vector<int>& faces, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix ,bool concave=false,float deep=2.0);
	void splitPoint(MMeshT* mt, MMeshVertex* v, trimesh::point ori);	
	void embedingAndCutting(MMeshT* mesh,const std::vector<std::vector<trimesh::vec2>>& lines,const std::vector<int>& facesIndex,bool is_close=true);
	void wordToWorldPoint(const CameraParam& camera,const LetterParam& letter, const std::vector<ClipperLibXYZ::Paths>& paths,std::vector<std::vector<trimesh::point>>& points);
	trimesh::point getWorldPoint(const CameraParam& camera, trimesh::ivec2 p);
	bool intersectionTriangle(MMeshT* mt,trimesh::point p,trimesh::point normal);
	void polygonInnerFaces(MMeshT* mt, std::vector<std::vector<std::vector<trimesh::vec2>>>& poly, std::vector<int>& infaceIndex, std::vector<int>& outfaceIndex);
	void getScreenWidthAndHeight(const CameraParam& camera,std::pair<float,float>&wh);
	void getViewMatrixAndProjectionMatrix(const CameraParam& camera,Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix);
	void loadCameraParam(CameraParam& camera);
	void getEmbedingPoint(std::vector<std::vector<trimesh::point>>& lines, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix, std::vector<std::vector<trimesh::vec2>>& poly);
	void unTransformationMesh(MMeshT* mesh, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix);
	void unTransformationMesh(trimesh::TriMesh* mesh, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix);
	void TransformationMesh(MMeshT* mesh, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix);
	void TransformationMesh(trimesh::TriMesh* mesh, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix);
	void fillTriangle(MMeshT* mesh, std::vector<int>& vindex);
	void fillTriangleForTraverse(MMeshT* mesh, std::vector<int>& vindex,bool is_rollback=false);
	void getMeshFaces(MMeshT* mesh, const std::vector<std::vector<trimesh::vec2>>& polygons, const CameraParam& camera, std::vector<int>& faces);
	void getMeshFaces(trimesh::TriMesh* mesh, const std::vector<std::vector<std::vector<trimesh::vec2>>>& polygons, const CameraParam& camera, std::vector<int>& faces,float threshold);
	void getDisCoverFaces(MMeshT* mesh, std::vector<int>& faces, std::map<int, int>& fmap);
	void mapping(MMeshT* mesh, trimesh::TriMesh* trimesh, std::map<int, int>& vmap, std::map<int, int>& fmap,bool is_thread=false);
	void fillholes(trimesh::TriMesh* mesh);
	void simpleCutting(MMeshT* mesh, const std::vector<std::vector<std::vector<trimesh::vec2>>>& polygons, std::vector<std::vector<int>>& faceindexs);
	trimesh::TriMesh* letter(trimesh::TriMesh* mesh, const SimpleCamera& camera, const LetterParam& Letter, const std::vector<TriPolygons>& polygons, bool& letterOpState,
		LetterDebugger* debugger = nullptr, ccglobal::Tracer* tracer = nullptr);
	trimesh::TriMesh* letter(trimesh::TriMesh* mesh, const SimpleCamera& camera, const LetterParam& Letter,const std::vector<TriPolygons>& polygons,
		LetterDebugger* debugger = nullptr, ccglobal::Tracer* tracer = nullptr);
}

#endif // TOPOMESH_LETTER_1680853426716_H