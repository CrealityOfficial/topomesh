#ifndef TOPOMESH_CONVERT_1681265581987_H
#define TOPOMESH_CONVERT_1681265581987_H
#include "trimesh2/TriMesh.h"
#include "topomesh/data/mmesht.h"
#include "clipperxyz/clipper.hpp"

namespace topomesh
{
	struct LetterParam
	{
		bool concave = true;
		float deep = 2.0f;
	};

	struct Camera
	{

	};

	class LetterInput
	{
	public:
		MMeshT mesh;
		Camera camera;
		LetterParam param;
		ClipperLibXYZ::Paths paths;

		void save(const std::string& fileName);
		void load(const std::string& fileName);
	};
}

#endif // TOPOMESH_CONVERT_1681265581987_H