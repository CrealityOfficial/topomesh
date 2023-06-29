#ifndef TOPOMESH_SERIAL_1688038279478_H
#define TOPOMESH_SERIAL_1688038279478_H
#include "topomesh/data/convert.h"
#include <fstream>

namespace topomesh
{
	void load(const std::string& fileName, TriPolygons& polys);
	void save(const std::string& fileName, const TriPolygons& polys);

	void load(std::fstream& in, TriPolygons& polys);
	void save(std::fstream& out, const TriPolygons& polys);
}

#endif // TOPOMESH_SERIAL_1688038279478_H