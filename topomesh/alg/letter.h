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
#include <vector>
#include <numeric>

namespace topomesh
{	
	void lettering(MMeshT* mesh, const std::vector<ClipperLibXYZ::Paths>& paths,  CameraParam& camera, const LetterParam& Letter, std::vector<int>* faceindex = nullptr);
	void concaveOrConvexOfFaces(MMeshT* mt,std::vector<int>& faces, bool concave=false);
	void splitPoint(MMeshT* mt, MMeshVertex* v, trimesh::point ori);	
	void embedingAndCutting(MMeshT* mesh, std::vector<trimesh::point>& lines, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix);
	void wordToWorldPoint(const CameraParam& camera,const LetterParam& letter, const std::vector<ClipperLibXYZ::Paths>& paths,std::vector<trimesh::point>& points);
	trimesh::point getWorldPoint(const CameraParam& camera, trimesh::ivec2 p);
	bool intersectionTriangle(MMeshT* mt,trimesh::point p,trimesh::point normal);
	void polygonInnerFaces(MMeshT* mt, std::vector<trimesh::vec3>& poly, std::vector<int>& faceIndex);
	void getScreenWidthAndHeight(const CameraParam& camera,std::pair<float,float>&wh);
	void getViewMatrixAndProjectionMatrix(const CameraParam& camera,Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix);
	void loadCameraParam(CameraParam& camera);
	void getEmbedingPoint(std::vector<trimesh::point>& lines, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix, std::vector<trimesh::vec3>& poly);
	void unTransformationMesh(MMeshT* mesh, Eigen::Matrix4f& ViewMatrix, Eigen::Matrix4f& ProjectionMatrix);
}

#endif // TOPOMESH_LETTER_1680853426716_H