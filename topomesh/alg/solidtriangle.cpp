#include "solidtriangle.h"

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
		_result.resize(_row, std::vector<float>(_col, std::numeric_limits<float>::min()));
		_length_x = (_bbox_max_x - _bbox_min_x) * 1.f / (_col * 1.f);
		_length_y = (_bbox_max_y - _bbox_min_y) * 1.f / (_row * 1.f);
	}

	void SolidTriangle::work()
	{
		for (int i = 0; i < _data->size(); i++)
		{
			trimesh::point v0;
			trimesh::point v1;
			trimesh::point v2;
			std::tie(v0, v1, v2) = _data->at(i);
			int xi0 = (v0.x - _bbox_min_x) / _length_x;
			int yi0 = (v0.y - _bbox_min_y) / _length_y;
			xi0 == _col ? xi0-- : xi0; yi0 == _col ? yi0-- : yi0;
			_result[xi0][yi0] = v0.z;
			int xi1 = (v1.x - _bbox_min_x) / _length_x;
			int yi1 = (v1.y - _bbox_min_y) / _length_y;
			xi1 == _col ? xi1-- : xi1; yi1 == _col ? yi1-- : yi1;
			_result[xi1][yi1] = v1.z;
			int xi2 = (v2.x - _bbox_min_x) / _length_x;
			int yi2 = (v2.y - _bbox_min_y) / _length_y;
			xi2 == _col ? xi2-- : xi2; yi2 == _col ? yi2-- : yi2;
			_result[xi2][yi2] = v2.z;
			//暂时用最低值，加快运行
			float mz[] = { v0.z,v1.z,v2.z };
			std::sort(mz, mz + 3);
			float min_z = mz[0];

			int xa[] = { xi0,xi1,xi2 };
			int ya[] = { yi0,yi1,yi2 };
			std::sort(xa, xa + 3);
			std::sort(ya, ya + 3);
			int min_x = xa[0];
			int max_x = xa[2];
			int min_y = ya[0];
			int max_y = ya[2];


			int I1 = v0.y - v1.y, I2 = v1.y - v2.y, I3 = v2.y - v0.y;
			int J1 = v1.x - v0.x, J2 = v2.x - v1.x, J3 = v0.x - v2.x;
			int F1 = (v0.x * v1.y - v0.y * v1.x);
			int F2 = (v1.x * v2.y - v1.y * v2.x);
			int F3 = (v2.x * v0.y - v2.y * v0.x);
			int CY1 = F1, CY2 = F2, CY3 = F3;
			for (int y = min_y; y <= max_y; y++)
			{
				int CX1 = CY1, CX2 = CY2, CX3 = CY3;
				for (int x = min_x; x <= max_x; x++)
				{
					if (CX1 > 0 && CX2 > 0 && CX3 > 0)
						_result[x][y] = min_z;
					CX1 += I1; CX2 += I2; CX3 += I3;
				}
				CY1 += J1; CY2 += J2; CY3 += J3;
			}
		}
	}
}