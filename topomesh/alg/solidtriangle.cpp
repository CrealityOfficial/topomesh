#include "solidtriangle.h"
#ifdef _WIN32
#include "omp.h"
#endif
#include <intrin.h>


#ifndef EPS
#define EPS 1e-6f
#endif 

namespace topomesh {
	SolidTriangle::SolidTriangle(const std::vector<std::tuple<trimesh::point, trimesh::point, trimesh::point>>* data,
		int row, int col,
		float bbox_max_x,float bbox_min_x,float bbox_max_y,float bbox_min_y)
		:_data(data),_row(row),_col(col),
		_bbox_max_x(bbox_max_x),
		_bbox_min_x(bbox_min_x),
		_bbox_max_y(bbox_max_y),
		_bbox_min_y(bbox_min_y)
	{
		_result.resize(_row, std::vector<float>(_col, std::numeric_limits<float>::max()));
		_length_x = (_bbox_max_x - _bbox_min_x) * 1.f / (_col * 1.f);
		_length_y = (_bbox_max_y - _bbox_min_y) * 1.f / (_row * 1.f);
	}

	void SolidTriangle::work()
	{
#ifdef _OPENMP
#pragma omp parallel for 
#endif 
		for (int i = 0; i < _data->size(); i++)
		{
			trimesh::point v0;
			trimesh::point v1;
			trimesh::point v2;
			std::tie(v0, v1, v2) = _data->at(i);			
			
			float xa[] = { v0.x,v1.x,v2.x };
			float ya[] = { v0.y,v1.y,v2.y };
			float za[] = { v0.z,v1.z,v2.z };
			std::sort(za, za + 3);
			std::sort(xa, xa + 3);
			std::sort(ya, ya + 3);
			float min_z = za[0];
			float min_x = xa[0];
			float max_x = xa[2];
			float min_y = ya[0];
			float max_y = ya[2];

			int min_xi = (min_x - _bbox_min_x) / _length_x;
			int min_yi = (min_y - _bbox_min_y) / _length_y;
			int max_xi = (max_x - _bbox_min_x) / _length_x;
			int max_yi = (max_y - _bbox_min_y) / _length_y;
			if (max_xi == _col) max_xi--;
			if (max_yi == _row) max_yi--;
			float I1 = v0.y - v1.y, I2 = v1.y - v2.y, I3 = v2.y - v0.y;
			float J1 = v1.x - v0.x, J2 = v2.x - v1.x, J3 = v0.x - v2.x;
			float F1 = (v0.x * v1.y - v0.y * v1.x);
			float F2 = (v1.x * v2.y - v1.y * v2.x);
			float F3 = (v2.x * v0.y - v2.y * v0.x);

			float start1 = I1 * (_bbox_min_x + min_xi * _length_x + _length_x / 4.f) + J1 * (_bbox_min_y + min_yi * _length_y + _length_y / 4.f) + F1;
			float start2 = I2 * (_bbox_min_x + min_xi * _length_x + _length_x / 4.f) + J2 * (_bbox_min_y + min_yi * _length_y + _length_y / 4.f) + F2;
			float start3 = I3 * (_bbox_min_x + min_xi * _length_x + _length_x / 4.f) + J3 * (_bbox_min_y + min_yi * _length_y + _length_y / 4.f) + F3;
			float CY1 = start1, CY2 = start2, CY3 = start3;
			I1 = I1 * _length_x; I2 = I2 * _length_x; I3 = I3 * _length_x;
			J1 = J1 * _length_y; J2 = J2 * _length_y; J3 = J3 * _length_y;
			for (int y = min_yi; y <= max_yi; y++)
			{
				float CX1 = CY1, CX2 = CY2, CX3 = CY3;
				for (int x = min_xi; x <= max_xi; x++)
				{
					if (CX1 >= -EPS && CX2 >= -EPS && CX3 >= -EPS)
#ifdef _OPENMP
#pragma omp critical
#endif 
						_result[x][y] = min_z;
					CX1 += I1; CX2 += I2; CX3 += I3;
				}
				CY1 += J1; CY2 += J2; CY3 += J3;
			}
			
		}
	}
}