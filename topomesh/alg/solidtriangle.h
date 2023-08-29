#pragma once
#include "trimesh2/TriMesh.h"

namespace topomesh {
	class SolidTriangle {
	public :
		SolidTriangle() {};
		SolidTriangle(const std::vector<std::tuple<trimesh::point,trimesh::point,trimesh::point>>*,int,int,float,float,float,float);
		~SolidTriangle() {};


		void work();
		const std::vector<std::vector<float>>& getResult() { return _result; };
		float getData(int r, int c) { if (r >= _row || c >= _col) return std::numeric_limits<float>::max(); return _result[r][c]; };
		float getCoordData(float x, float y) { int xi = (x - _bbox_min_x) / (_col * 1.f); int yi = (y - _bbox_min_y) / (_row * 1.f); return _result[xi][yi]; };
	private:
		const std::vector<std::tuple<trimesh::point, trimesh::point, trimesh::point>>* _data=nullptr;
		int _row;
		int _col;
		float _bbox_min_x;
		float _bbox_min_y;
		float _bbox_max_x;
		float _bbox_max_y;
		float _length_x;
		float _length_y;
		std::vector<std::vector<float>> _result;
	};
}