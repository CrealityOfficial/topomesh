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
            trimesh::normalize(x);
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
                auto&& x = std::move(z.cross(n));
                if (trimesh::len(x) < EPS) {
                    x = trimesh::vec3(1, 0, 0);
                }
                trimesh::normalize(x);
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
    trimesh::TriMesh SaveTriPolygonPointsToMesh(const TriPolygon& poly, double radius, size_t nrows, size_t ncolumns)
    {
        trimesh::TriMesh pointMesh;
        const size_t spherenums = poly.size();
        const size_t nums = nrows * ncolumns + 2;
        const size_t pointnums = nums * spherenums;
        const size_t facenums = 2 * (nums - 2) * spherenums;
        std::vector<trimesh::vec3> points;
        points.reserve(pointnums);
        std::vector<trimesh::ivec3> faces;
        faces.reserve(facenums);
        for (int k = 0; k < spherenums; ++k) {
            const auto& p = poly[k];
            points.emplace_back(p.x, p.y, p.z + radius);
            for (int i = 0; i < nrows; ++i) {
                const auto& phi = M_PI * (i + 1.0) / double(nrows + 1.0);
                const auto& z = radius * std::cos(phi);
                const auto& r = radius * std::sin(phi);
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& theta = 2.0 * M_PI * j / ncolumns;
                    const auto& x = r * std::cos(theta);
                    const auto& y = r * std::sin(theta);
                    points.emplace_back(p.x + x, p.y + y, p.z + z);
                }
            }
            points.emplace_back(p.x, p.y, p.z - radius);
            const auto& maxInx = points.size() - 1;
            const auto& v0 = k * nums;
            //上下底部两部分
            for (size_t i = 0; i < ncolumns; ++i) {
                const auto& i0 = i + 1 + v0;
                const auto& i1 = (i + 1) % ncolumns + 1 + v0;
                faces.emplace_back(trimesh::ivec3(v0, i0, i1));
                const auto& j0 = i0 + (nrows - 1) * ncolumns;
                const auto& j1 = i1 + (nrows - 1) * ncolumns;
                faces.emplace_back(trimesh::ivec3(j1, j0, maxInx));
            }
            //中间部分
            for (size_t i = 0; i < nrows - 1; ++i) {
                const auto& j0 = i * ncolumns + 1 + v0;
                const auto& j1 = (i + 1) * ncolumns + 1 + v0;
                for (size_t j = 0; j < ncolumns; ++j) {
                    const auto& i0 = j0 + j;
                    const auto& i1 = j0 + (j + 1) % ncolumns;
                    const auto& i2 = j1 + j;
                    const auto& i3 = j1 + (j + 1) % ncolumns;
                    faces.emplace_back(trimesh::ivec3(i2, i1, i0));
                    faces.emplace_back(trimesh::ivec3(i2, i3, i1));
                }
            }
        }
        pointMesh.vertices.swap(points);
        pointMesh.faces.swap(faces);
        return pointMesh;
    }
    trimesh::TriMesh SaveTriPolygonsPointsToMesh(const TriPolygons& polys, double radius, size_t nrows, size_t ncolumns)
    {
        trimesh::TriMesh pointMesh;
        TriPolygon poly;
        size_t spherenums = 0;
        for (const auto& pys : polys) {
            spherenums += pys.size();
        }
        poly.reserve(spherenums);
        for (const auto& pys : polys) {
            for (const auto& p : pys) {
                poly.emplace_back(p);
            }
        }
        return SaveTriPolygonPointsToMesh(poly, radius, nrows, ncolumns);
    }
}