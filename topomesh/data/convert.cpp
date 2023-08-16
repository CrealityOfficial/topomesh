#include "convert.h"

namespace topomesh
{
    TriPolygons convertFromHexaPolygons(const HexaPolygons& hexas)
    {
        TriPolygons polys;
        polys.reserve(hexas.size());
        for (const auto& h : hexas) {
            polys.emplace_back(std::move(h.poly));
        }
        return polys;
    }
    HexaPolygons convertFromTriPolygons(const TriPolygons& polys, const trimesh::ivec3& coords)
    {
        HexaPolygons hexas;
        hexas.reserve(polys.size());
        const int size = polys.size();
        if (polys.size() == coords.size()) {
            for (int i = 0; i < size; ++i) {
                hexas.emplace_back();
                auto& hexa = hexas.back();
                hexa.poly = polys[i];
                hexa.coord = coords[i];
            }
        } else {
            for (int i = 0; i < size; ++i) {
                hexas.emplace_back();
                auto& hexa = hexas.back();
                hexa.poly = polys[i];
            }
        }
        return hexas;
    }
}