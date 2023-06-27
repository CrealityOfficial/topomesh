#ifndef CRCOMMON_HEADER_INTERFACE
#define CRCOMMON_HEADER_INTERFACE
#include "topomesh/math/interface.h"

#include "ccglobal/tracer.h"
#include "ccglobal/log.h"

#include "trimesh2/TriMesh.h"
#include "trimesh2/XForm.h"
#include "trimesh2/TriMesh_algo.h"

#include <memory>
#include <fstream>
#include <unordered_map>

typedef std::shared_ptr<trimesh::TriMesh> TriMeshPtr;

#include "topomesh/math/Settings.h"
typedef std::shared_ptr<crcommon::Settings> SettingsPtr;

typedef std::unordered_map<std::string, std::string> KValues;

template<class T>
void templateSave(const T& t, std::ofstream& out)
{
	out.write((const char*)&t, sizeof(T));
}

template<class T>
T templateLoad(std::ifstream& in)
{
	T t;
	in.read((char*)&t, sizeof(T));
	return t;
}

template<class T>
void templateSave(const std::vector<T>& ts, std::ofstream& out)
{
	int size = (int)ts.size();
	templateSave<int>(size, out);
	out.write((const char*)ts.data(), size * sizeof(T));
}

template<class T>
void templateLoad(std::vector<T>& ts, std::ifstream& in)
{
	int size = templateLoad<int>(in);
	if (size > 0)
	{
		ts.resize(size);
		in.read((char*)ts.data(), size * sizeof(T));
	}
}
#endif // CRCOMMON_HEADER_INTERFACE