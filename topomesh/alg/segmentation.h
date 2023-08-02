#pragma once
#include "topomesh/data/convert.h"
#include "topomesh/data/mmesht.h"

namespace topomesh {
	class SpectralClusteringCuts
	{
	public:
		SpectralClusteringCuts() {};
		SpectralClusteringCuts(float delte, float eta):_delta(delte), _eta(eta) {};
		SpectralClusteringCuts(MMeshT* mesh, float delte, float eta);
		~SpectralClusteringCuts() {};

		void BlockSpectralClusteringCuts(MMeshT* mesh);
		std::vector<std::vector<int>>* get_result() { return &result; }
		void test();
	private:
		std::vector<std::vector<int>> result;
		float _delta;
		float _eta;
	};
	
}