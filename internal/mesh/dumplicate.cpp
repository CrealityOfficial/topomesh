#include "dumplicate.h"
#include "trimesh2/TriMesh_algo.h"
#include <map>

namespace topomesh
{
    bool testNeedfitMesh(trimesh::TriMesh* mesh, float& scale)
    {
        if (!mesh)
            return false;

        mesh->need_bbox();
        trimesh::vec3 size = mesh->bbox.size();

        bool needScale = false;
        scale = 1.0f;
        if (size.max() > 1000.0f)
        {
            needScale = true;
            scale = 100.0f / size.max();
        }
        else if (size.min() < 1.0f && size.min() > 0.00001f)
        {
            needScale = true;
            scale = 100.0f / size.min();
        }

        return needScale;
    }

    struct hash_vec3 {
        size_t operator()(const trimesh::vec3& v)const
        {
            return std::abs(v.x) * 10000.0f / 23.0f + std::abs(v.y) * 10000.0f / 19.0f + std::abs(v.z) * 10000.0f / 17.0f;
        }
    };

    struct hash_func1 {
        size_t operator()(const trimesh::vec3& v)const
        {
            return (int(v.x * 99971)) ^ (int(v.y * 99989)) ^ (int(v.z * 99991));
        }
    };

    struct hash_func2 {
        size_t operator()(const trimesh::vec3& v)const
        {          
            size_t r = (int(v.x * 99971)) ^ (int(v.y * 99989) << 1) ^ (int(v.z * 99991) << 2);
            return r;
        }
    };

    bool hashMesh(trimesh::TriMesh* mesh, std::function<size_t(const trimesh::vec3&)> hash_func, ccglobal::Tracer* tracer, float ratio = 0.3f)
    {
        if (!mesh)
            return false;

        float sValue = 1.0f;
        bool needScale = testNeedfitMesh(mesh, sValue);

        if (needScale)
            trimesh::apply_xform(mesh, trimesh::xform::scale(sValue));

        size_t vertexNum = mesh->vertices.size();

        struct equal_vec3 {
            bool operator()(const trimesh::vec3& v1, const trimesh::vec3& v2) const
            {
                return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
            }
        };

        typedef std::unordered_map<trimesh::vec3, size_t, decltype(hash_func), equal_vec3> unique_point;
        typedef typename unique_point::value_type unique_value;
        const int buckets = vertexNum * ratio + 1;
        unique_point points(buckets, hash_func);


        size_t faceNum = mesh->faces.size();

        if (vertexNum == 0 || faceNum == 0)
            return false;

        trimesh::TriMesh* optimizeMesh = new trimesh::TriMesh();
        bool interuptted = false;

        std::vector<int> vertexMapper;
        vertexMapper.resize(vertexNum, -1);

        if (tracer)
            tracer->formatMessage("dumplicateMesh %d", (int)vertexNum);

        size_t pVertex = vertexNum / 20;
        if (pVertex == 0)
            pVertex = vertexNum;

        for (size_t i = 0; i < vertexNum; ++i) {
            trimesh::vec3 p = mesh->vertices.at(i);
            auto it = points.find(p);
            if (it != points.end()) {
                int index = (*it).second;
                vertexMapper.at(i) = index;
            }
            else {
                int index = (int)points.size();
                points.insert(unique_value(p, index));

                vertexMapper.at(i) = index;
            }

            if (i % pVertex == 1) {
                if (tracer) {
                    tracer->formatMessage("dumplicateMesh %i", (int)i);
                    tracer->progress((float)i / (float)vertexNum);
                    if (tracer->interrupt()) {
                        interuptted = true;
                        break;
                    }
                }
            }
        }

        if (tracer)
            tracer->formatMessage("dumplicateMesh over %d", (int)points.size());

        if (interuptted) {
            delete optimizeMesh;
            return false;
        }
        trimesh::TriMesh* omesh = optimizeMesh;
        omesh->vertices.resize(points.size());
        for (auto it = points.begin(); it != points.end(); ++it) {
            omesh->vertices.at(it->second) = it->first;
        }

        omesh->faces = mesh->faces;
        for (size_t i = 0; i < faceNum; ++i) {
            trimesh::TriMesh::Face& of = omesh->faces.at(i);
            for (int j = 0; j < 3; ++j) {
                int index = of[j];
                of[j] = vertexMapper[index];
            }
        }

        mesh->vertices.swap(omesh->vertices);
        mesh->faces.swap(omesh->faces);
        mesh->flags.clear();

        if (needScale)
            trimesh::apply_xform(mesh, trimesh::xform::scale(1.0f / sValue));

        mesh->flags.clear();
        mesh->clear_bbox();
        mesh->need_bbox();

        delete omesh;
        return true;
    }

    bool dumplicateMesh(trimesh::TriMesh* mesh, ccglobal::Tracer* tracer, float ratio)
    {
        return hashMesh(mesh, hash_func1(), tracer, ratio);
    }
    bool mergeNearPoints(trimesh::TriMesh* mesh, ccglobal::Tracer* tracer, float eps)
    {
        if (!mesh)
            return false;

        float sValue = 1.0f;
        bool needScale = testNeedfitMesh(mesh, sValue);

        if (needScale)
            trimesh::apply_xform(mesh, trimesh::xform::scale(sValue));

        struct PointComparator {
            float epsilon = 1E-8F;
            PointComparator(float error) :epsilon(error) {}
            bool operator()(const trimesh::vec3& v1, const trimesh::vec3& v2) const
            {
                if (std::fabs(v1.x - v2.x) > epsilon)
                    return v1.x < v2.x;
                else if (std::fabs(v1.y - v2.y) > epsilon)
                    return v1.y < v2.y;
                else if (std::fabs(v1.z - v2.z) > epsilon)
                    return v1.z < v2.z;
                return false;
            }
        };
        PointComparator compare(eps);
        std::map<trimesh::vec3, int, PointComparator> points(compare);
        size_t vertexNum = mesh->vertices.size();
        size_t faceNum = mesh->faces.size();

        if (vertexNum == 0 || faceNum == 0)
            return false;

        trimesh::TriMesh * optimizeMesh = new trimesh::TriMesh();
        bool interuptted = false;

        std::vector<int> vertexMapper;
        vertexMapper.resize(vertexNum, -1);

        if (tracer)
            tracer->formatMessage("mergeNearPoints %d", (int)vertexNum);

        size_t pVertex = vertexNum / 20;
        if (pVertex == 0)
            pVertex = vertexNum;

        auto Equal = [&](const trimesh::dvec3 & v1, const trimesh::dvec3 & v2, double epsilon = 1E-8) {
            const auto& v = v1 - v2;
            return (std::fabs(v.x) < epsilon) && (std::fabs(v.y) < epsilon) && (std::fabs(v.z) < epsilon);
        };
        for (size_t i = 0; i < vertexNum; ++i) {
            trimesh::vec3 p = mesh->vertices.at(i);
            auto it = points.find(p);
            if (it != points.end()) {
                int index = (*it).second;
                vertexMapper.at(i) = index;
            } else {
                int index = (int)points.size();
                points.emplace(p, index);
                vertexMapper.at(i) = index;
            }

            if (i % pVertex == 1) {
                if (tracer) {
                    tracer->formatMessage("mergeNearPoints %i", (int)i);
                    tracer->progress((float)i / (float)vertexNum);
                    if (tracer->interrupt()) {
                        interuptted = true;
                        break;
                    }
                }
            }
        }

        if (tracer)
            tracer->formatMessage("mergeNearPoints over %d", (int)points.size());

        if (interuptted) {
            delete optimizeMesh;
            return false;
        }
        trimesh::TriMesh* omesh = optimizeMesh;
        omesh->vertices.resize(points.size());
        for (auto it = points.begin(); it != points.end(); ++it) {
            omesh->vertices.at(it->second) = it->first;
        }

        omesh->faces = mesh->faces;
        for (size_t i = 0; i < faceNum; ++i) {
            trimesh::TriMesh::Face& of = omesh->faces.at(i);
            for (int j = 0; j < 3; ++j) {
                int index = of[j];
                of[j] = vertexMapper[index];
            }
        }

        mesh->vertices.swap(omesh->vertices);
        mesh->faces.swap(omesh->faces);
        mesh->flags.clear();

        if (needScale)
            trimesh::apply_xform(mesh, trimesh::xform::scale(1.0f / sValue));

        mesh->flags.clear();
        mesh->clear_bbox();
        mesh->need_bbox();

        delete omesh;
        return true;
    }
}