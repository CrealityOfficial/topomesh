#include "serial.h"
#include "ccglobal/serial/serial.h"

namespace topomesh
{
	void load(const std::string& fileName, TriPolygons& polys)
	{
		serialFunc f = [&polys](std::fstream& stream) {
			load(stream, polys);
		};
		serialLoad(fileName, f);
	}

	void save(const std::string& fileName, const TriPolygons& polys)
	{
		serialFunc f = [&polys](std::fstream& stream) {
			save(stream, polys);
		};
		serialSave(fileName, f);
	}

	void load(std::fstream& in, TriPolygons& polys)
	{
		int size = 0;
		loadT(in, size);
		if (size > 0)
		{
			polys.resize(size);
			for (int i = 0; i < size; ++i)
				loadVectorT(in, polys.at(i));
		}
	}

	void save(std::fstream& out, const TriPolygons& polys)
	{
		int size = (int)polys.size();
		saveT(out, size);
		for (int i = 0; i < size; ++i)
			saveVectorT(out, polys.at(i));
	}
}