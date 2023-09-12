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

    trimesh::TriMesh SaveTriPolygonToMesh(const TriPolygon& poly, double r, size_t nslices)
    {
        trimesh::TriMesh edgeMesh;
        const size_t polysize = poly.size();
        std::vector<trimesh::vec3> points;
        std::vector<trimesh::ivec3> faces;
        points.reserve(2 * polysize * nslices);
        faces.reserve(2 * polysize * nslices);
        double delta = 2.0 * M_PI / nslices;
        trimesh::vec3 z(0, 0, 1);
        for (size_t i = 0; i < polysize; ++i) {
            const auto& a = poly[i];
            const auto& b = poly[(i + 1) % polysize];
            const auto& n = trimesh::normalized(b - a);
            auto && x = std::move(z.cross(n));
            if (trimesh::len(x) < EPS) {
                x = trimesh::vec3(1, 0, 0);
            }
            const auto& y = n.cross(x);
            for (int j = 0; j < nslices; ++j) {
                const auto& theta = delta * j;
                const auto& p = b + x * r * std::cos(theta) + y * r * std::sin(theta);
                points.emplace_back(p);
            }
            for (int j = 0; j < nslices; ++j) {
                const auto& theta = delta * j;
                const auto& p = a + x * r * std::cos(theta) + y * r * std::sin(theta);
                points.emplace_back(p);
            }
        }
        for (size_t i = 0; i < polysize; ++i) {
            for (size_t j = 0; j < nslices; ++j) {
                const auto& i0 = j + 2 * nslices * i;
                const auto& i1 = (j + 1) % nslices + 2 * nslices * i;
                const auto & j0 = i0 + nslices;
                const auto & j1 = i1 + nslices;
                faces.emplace_back(j0, i1, i0);
                faces.emplace_back(j0, j1, i1);
            }
        }
        edgeMesh.faces.swap(faces);
        edgeMesh.vertices.swap(points);
        return edgeMesh;
    }
    trimesh::TriMesh SaveTriPolygonsToMesh(const TriPolygons & polys, double r, size_t nslices)
    {
        trimesh::TriMesh edgeMesh;
        size_t pointnums = 0;
        for (int i = 0; i < polys.size(); ++i) {
            pointnums += polys[i].size();
        }
        std::vector<trimesh::vec3> points;
        std::vector<trimesh::ivec3> faces;
        points.reserve(2 * pointnums * nslices);
        faces.reserve(2 * pointnums * nslices);
        double delta = 2.0 * M_PI / nslices;
        trimesh::vec3 z(0, 0, 1);
        int start = 0;
        for (size_t k = 0; k < polys.size(); ++k) {
            const auto& poly = polys[k];
            const int polysize = poly.size();
            for (size_t i = 0; i < polysize; ++i) {
                const auto& a = poly[i];
                const auto& b = poly[(i + 1) % polysize];
                const auto& n = trimesh::normalized(b - a);
                auto && x = std::move(z.cross(n));
                if (trimesh::len(x) < EPS) {
                    x = trimesh::vec3(1, 0, 0);
                }
                const auto& y = n.cross(x);
                for (int j = 0; j < nslices; ++j) {
                    const auto& theta = delta * j;
                    const auto& p = b + x * r * std::cos(theta) + y * r * std::sin(theta);
                    points.emplace_back(p);
                }
                for (int j = 0; j < nslices; ++j) {
                    const auto& theta = delta * j;
                    const auto& p = a + x * r * std::cos(theta) + y * r * std::sin(theta);
                    points.emplace_back(p);
                }
            }
            for (size_t i = 0; i < polysize; ++i) {
                for (size_t j = 0; j < nslices; ++j) {
                    const auto& i0 = start + j + 2 * nslices * i;
                    const auto& i1 = start + (j + 1) % nslices + 2 * nslices * i;
                    const auto & j0 = i0 + nslices;
                    const auto & j1 = i1 + nslices;
                    faces.emplace_back(j0, i1, i0);
                    faces.emplace_back(j0, j1, i1);
                }
            }
            start = points.size();
        }
        edgeMesh.faces.swap(faces);
        edgeMesh.vertices.swap(points);
        return edgeMesh;
    }
}