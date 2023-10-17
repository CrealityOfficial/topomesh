#include "volumeMesh.h"

namespace topomesh {
	float getMeshVolume(trimesh::TriMesh* mesh, std::vector<int>& faces)
	{
		float vol = 0.f;
		for (int fi = 0; fi < faces.size(); fi++)
		{
			trimesh::point v0 = mesh->vertices[mesh->faces[fi][0]];
			trimesh::point v1 = mesh->vertices[mesh->faces[fi][1]];
			trimesh::point v2 = mesh->vertices[mesh->faces[fi][2]];

			float v = (-v2.x * v1.y * v0.z + v1.x * v2.y * v0.z + v2.x * v0.y * v1.z - v0.x * v2.y * v1.z - v1.x * v0.y * v2.z + v0.x * v1.y * v2.z)/6.f;
			vol += v;
		}
		return vol;
	}


}