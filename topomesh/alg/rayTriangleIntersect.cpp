#include "rayTriangleIntersect.h"

namespace topomesh
{
	bool rayTriangleIntersect(const point& v0, const point& v1, const point& v2, const point& ori, const point& dir, point& result)
	{
		point S = ori - v0;
		point E1 = v1 - v0;
		point E2 = v2 - v1;
		point S1 = dir % E2;
		point S2 = S % E1;

		float S1E1 = S1 ^ E1;
		float t = (S2 ^ E2)*1.f / S1E1*1.f;
		float b1 = (S1 ^ S) * 1.f / S1E1 * 1.f;
		float b2 = (S2 ^ dir) * 1.f / S1E1 * 1.f;

		if (t > 0.f && b1 >= 0.f && b2 >= 0.f && (1 - b1 - b2) >= 0.f)
		{
			result = ori + t * dir;
			return true;
		}
		return false;
	}

	bool faceTriangleIntersect(const trimesh::point& va1, const trimesh::point& va2, const trimesh::point& va3,
		const trimesh::point& vb1, const trimesh::point& vb2, const trimesh::point& vb3,
		trimesh::point& cross1, trimesh::point& cross2)
	{
		return false;
	}
}