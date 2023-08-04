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
	
	struct SegmentationParam
	{

	};

	enum class SegmentatStrategy
	{
		ss_angle,
		ss_sdf,
		ss_num
	};

	class Segmentation
	{
	public:
		Segmentation(trimesh::TriMesh* mesh);
		~Segmentation();

		void autoSegment(int num);   // num  <= 0 auto

		//seed
		int createGroup();
		void removeGroup(int index);
		void addSeed2Group(int groupIndex, int index);
		void addSeeds2Group(int groupIndex, const std::vector<int>& indices);
		void removeGroupSeed(int groupIndex, int index);

		const FacePatchs& currentPatches();
	protected:
		FacePatchs m_patches;
	};
}