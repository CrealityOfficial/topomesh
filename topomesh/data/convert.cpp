#include "convert.h"

namespace topomesh
{
    void translateTriPolygon(TriPolygon& poly, const trimesh::vec3& trans)
    {
        for (auto& p : poly) {
            p += trans;
        }
    }
    TriPolygons convertFromHexaPolygons(const HexaPolygons& hexas)
    {
        TriPolygons polys;
        polys.reserve(hexas.polys.size());
        for (const auto& h : hexas.polys) {
            polys.emplace_back(std::move(h.poly));
        }
        return polys;
    }
    HexaPolygons convertFromTriPolygons(const TriPolygons& polys, const trimesh::ivec3& coords)
    {
        HexaPolygons hexas;
        hexas.polys.reserve(polys.size());
        const int size = polys.size();
        if (polys.size() == coords.size()) {
            for (int i = 0; i < size; ++i) {
                hexas.polys.emplace_back();
                auto& hexa = hexas.polys.back();
                hexa.poly = polys[i];
                hexa.coord = coords[i];
            }
        } else {
            for (int i = 0; i < size; ++i) {
                hexas.polys.emplace_back();
                auto& hexa = hexas.polys.back();
                hexa.poly = polys[i];
            }
        }
        return hexas;
    }
}