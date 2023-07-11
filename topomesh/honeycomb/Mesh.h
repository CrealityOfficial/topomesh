#pragma once
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <condition_variable>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <numeric>

#include "Matrix.h"
//#include "Polyline.h"
namespace honeycomb {
    //暂时未用到
    /*class HalfEdge {
    public:
        int a, b; ///<逆时针排序
        HalfEdge():a(0),b(0){}
        HalfEdge(int a0,int b0):a(a0),b(b0){}
        inline bool operator<(const HalfEdge& e) const
        {
            return a < e.a;
        }
        ~HalfEdge() {}
    };*/
    class Edge {
    public:
        int a, b; ///< a < b
        Edge() :a(0), b(0) {}
        Edge(int a0, int b0) :a(a0), b(b0) {}
        inline bool operator<(const Edge& e) const
        {
            return a < e.a /*|| (a == e.a && b < e.b)*/;
        }
        ~Edge() {}
    };

    class Mesh {
    public:
        bool ReadFromSTL(const char* filename)
        {
            // only support ASCII or BINARY format
            // https://zhuanlan.zhihu.com/p/398208443?utm_id=0
            Clear();
            std::ifstream fid;
            std::string format, header, tail;
            char buffer[80];
            fid.open(filename, std::ios::in | std::ios::binary);
            fid.seekg(0, std::ios::end);
            size_t fidsize = fid.tellg(); // Check the size of the file
            if ((fidsize - 84) % 50 > 0)
                format = "ascii";
            else {
                fid.seekg(0, std::ios::beg);//go to the beginning of the file
                fid.read(buffer, 80);
                header = buffer;
                bool isSolid = header.find("solid") + 1;
                fid.seekg(-80, std::ios::end);          // go to the end of the file minus 80 characters
                fid.read(buffer, 80);
                tail = buffer;
                bool isEndSolid = tail.find("endsolid") + 1;
                if (isSolid && isEndSolid)
                    format = "ascii";
                else format = "binary";
            }
            fid.close();
            size_t facenums = 0;
            std::vector<Point> originPts;
            if (format == "ascii")//判断格式
            {
                // ReadASCII
                fid.open(filename, std::ios::in);
                // Read the STL header
                std::getline(fid, header);
                std::string linestr;
                size_t np1, np2;
                float x, y, z;
                while (std::getline(fid, linestr)) {
                    std::string alpha = linestr.substr(0, linestr.find_first_of(' '));
                    if (alpha == "facet") {
                        if (linestr.back() == '\r') {
                            linestr.pop_back();
                        }
                        // skip "facet normal"
                        linestr = linestr.substr(13);
                        np1 = linestr.find_first_of(' ');
                        np2 = linestr.find_last_of(' ');
                        x = std::stof(linestr.substr(0, np1));
                        y = std::stof(linestr.substr(np1, np2));
                        z = std::stof(linestr.substr(np2));
                        //sscanf_s(linestr.c_str(), "facet normal %f %f %f", &x, &y, &z);
                        normals_.emplace_back(x, y, z);
                        ++facenums;
                    } else if (alpha == "vertex") {
                        if (linestr.back() == '\r') {
                            linestr.pop_back();
                        }
                        // skip "vertex "
                        linestr = linestr.substr(7);
                        np1 = linestr.find_first_of(' ');
                        np2 = linestr.find_last_of(' ');
                        x = std::stof(linestr.substr(0, np1));
                        y = std::stof(linestr.substr(np1, np2));
                        z = std::stof(linestr.substr(np2));
                        //sscanf_s(linestr.c_str(), "vertex %f %f %f", &x, &y, &z);
                        originPts.emplace_back(x, y, z);
                    }
                }
                fid.close();
            } else if (format == "binary") {
                //ReadBinary;
                char objectname[80];
                fid.open(filename, std::ios::in | std::ios::binary);
                fid.read(reinterpret_cast<char*>(&objectname), sizeof(objectname));
                object_name_.assign(objectname);
                //fid.seekg(80, std::ios::beg);
                fid.read(reinterpret_cast<char*>(&facenums), sizeof(int));
                originPts.resize(3 * facenums);
                normals_.resize(facenums);
                for (size_t i = 0; i < facenums; ++i) {
                    float pdata[12];
                    fid.read(reinterpret_cast<char*>(&pdata), sizeof(float) * 12);
                    //读取法向量信息
                    normals_[i] = std::move(Point(pdata[0], pdata[1], pdata[2]));
                    //读取三个顶点坐标
                    for (size_t j = 0; j < 3; ++j) {
                        originPts[3 * i + j] = std::move(Point(pdata[3 * j + 3], pdata[3 * j + 4], pdata[3 * j + 5]));
                    }
                    char tailbuffer[2];
                    fid.read(reinterpret_cast<char*>(&tailbuffer), sizeof(tailbuffer));
                }
                fid.close();
            } else {
                std::cout << "format error!" << std::endl;
                return false;
            }
            std::cout << facenums << std::endl;
            std::cout << originPts.size() << std::endl;

            pointIndexMap_.rehash(facenums * 2);
            faces_.resize(facenums);
            indexs_.resize(facenums);
            for (size_t i = 0; i < facenums; ++i) {
                std::vector<int> faceIndex;
                faceIndex.reserve(3); // 预留空间
                for (size_t j = 0; j < 3; ++j) {
                    const auto& p = originPts[3 * i + j];
                    auto iter = pointIndexMap_.find(p);
                    if (iter != pointIndexMap_.end()) {
                        faceIndex.emplace_back(iter->second);
                    } else {
                        const int newIndex = pointIndexMap_.size();
                        pointIndexMap_.emplace(p, newIndex);
                        //uniquePoints.insert(p);
                        points_.emplace_back(std::move(p));
                        faceIndex.emplace_back(newIndex);
                    }
                }
                faces_[i] = std::move(i);
                indexs_[i] = std::move(faceIndex);
            }
            return true;
        }
        //bBinary=true,写入二进制文件，bBinary=false写入ASCII文件
        bool WriteSTLFile(const char* filename, bool bBinary = true)
        {
            int face_num = faces_.size();
            if (normals_.empty()) {
                GenerateFaceNormals();
            }
            int normal_num = normals_.size();
            if (face_num != normal_num) {
                std::cout << "data error!" << std::endl;
                return false;
            }
            //https://blog.csdn.net/kxh123456/article/details/105814498/
            if (bBinary) {
                //二进制写入文件
                std::ofstream fs(std::string(filename) + "_bin.stl", std::ios::binary | std::ios::trunc);
                if (!fs) { fs.close(); return false; }
                int intSize = sizeof(int);
                int floatSize = sizeof(float);
                // 文件头
                char fileHead[3];
                for (int i = 0; i < 3; ++i) {
                    fileHead[i] = ' ';
                }
                fs.write(fileHead, sizeof(char) * 3);
                // 附加信息
                char fileInfo[77];
                for (int i = 0; i < 77; ++i)
                    fileInfo[i] = ' ';
                fs.write(fileInfo, sizeof(char) * 77);
                // 面的个数
                fs.write((char*)(&face_num), intSize);
                // 点列表，面列表
                char a[2];
                int a_size = sizeof(char) * 2;
                for (int i = 0; i < face_num; ++i) {
                    int pIndex0 = indexs_[i][0];
                    int pIndex1 = indexs_[i][1];
                    int pIndex2 = indexs_[i][2];
                    float nx = normals_[i].x;
                    float ny = normals_[i].y;
                    float nz = normals_[i].z;
                    //保存法向量
                    fs.write((char*)(&nx), floatSize);
                    fs.write((char*)(&ny), floatSize);
                    fs.write((char*)(&nz), floatSize);
                    // 保存顶点
                    float p0x = points_[pIndex0].x;
                    float p0y = points_[pIndex0].y;
                    float p0z = points_[pIndex0].z;
                    fs.write((char*)(&p0x), floatSize);
                    fs.write((char*)(&p0y), floatSize);
                    fs.write((char*)(&p0z), floatSize);

                    float p1x = points_[pIndex1].x;
                    float p1y = points_[pIndex1].y;
                    float p1z = points_[pIndex1].z;
                    fs.write((char*)(&p1x), floatSize);
                    fs.write((char*)(&p1y), floatSize);
                    fs.write((char*)(&p1z), floatSize);

                    float p2x = points_[pIndex2].x;
                    float p2y = points_[pIndex2].y;
                    float p2z = points_[pIndex2].z;
                    fs.write((char*)(&p2x), floatSize);
                    fs.write((char*)(&p2y), floatSize);
                    fs.write((char*)(&p2z), floatSize);
                    fs.write(a, a_size);
                }
                fs.close();
                return true;
            } else {
                //points_:模型顶点
                //faces_:三角面片
                //normals_:法线
                if (face_num <= 1e6) {
                    std::ofstream fs(std::string(filename) + "_ast.stl");
                    if (!fs) { fs.close(); return false; }
                    fs << "solid WRAP" << std::endl;
                    for (const auto& f : faces_) {
                        fs << "facet normal ";
                        const auto& n = normals_[f];
                        fs << (float)n.x << " " << (float)n.y << " " << (float)n.z << std::endl;
                        fs << "outer loop" << std::endl;
                        const auto& vertexs = indexs_[f];
                        for (const auto& inx : vertexs) {
                            fs << "vertex ";
                            const auto& p = points_[inx];
                            fs << (float)p.x << " " << (float)p.y << " " << (float)p.z << std::endl;
                        }
                        fs << "end loop" << std::endl;
                        fs << "end facet" << std::endl;
                    }
                    fs << "endsolid WRAP";
                    fs.close();
                    return true;
                }
                //针对大模型分块读取，多线程加速读取
                const int block_size = 10000; // 每个块包含的三角形数量
                const int num_threads = 8;     // 使用的线程数
                const int num_blocks = (face_num + block_size - 1) / block_size;
                // 创建输出文件
                std::ofstream of(std::string(filename) + "_ast.stl", std::ios::trunc);
                of << "solid ASCII_STL\n";
                of.close();
                // 分块写入数据
                std::vector<std::thread> threads(num_threads);
                std::mutex mutex;
                std::condition_variable cv;
                std::atomic<int> block_counter(0);
                const auto & faces = faces_;
                const auto & normals = normals_;
                const auto & points = points_;
                const auto & indexs = indexs_;
                std::ofstream out(std::string(filename) + "_ast.stl", std::ios::app);
                for (int i = 0; i < num_threads; ++i) {
                    threads[i] = std::thread([&]() {
                        while (true) {
                            int block_index = block_counter.fetch_add(1);
                            if (block_index >= num_blocks) {
                                break;
                            }
                            int start = block_index * block_size;
                            int end = std::min(start + block_size, face_num);
                            std::string buffer;
                            for (int j = start; j < end; ++j) {
                                const auto& f = faces[j];
                                const auto& n = normals[f];
                                buffer += "facet normal " + std::to_string(n.x) + " " + std::to_string(n.y) + " " + std::to_string(n.z) + "\n";
                                buffer += "outer loop\n";
                                const auto& neighbors = indexs[f];
                                for (const auto& inx : neighbors) {
                                    const auto& p = points[inx];
                                    buffer += "vertex " + std::to_string(p.x) + " " + std::to_string(p.y) + " " + std::to_string(p.z) + "\n";
                                }
                                buffer += "endloop\n";
                                buffer += "endfacet\n";
                                if (j % 1000 == 0) {
                                    std::unique_lock<std::mutex> lock(mutex);
                                    out << buffer;
                                    buffer.clear();
                                    cv.notify_all();
                                }
                            }
                            std::unique_lock<std::mutex> lock(mutex);
                            out << buffer;
                            cv.notify_all();
                        }
                    });
                }
                // 等待所有线程完成
                for (auto& thread : threads) {
                    thread.join();
                }
                // 写入文件尾
                out << "endsolid ASCII_STL\n";
                out.close();
                return true;
            }
        }
        std::string GetName()const
        {
            return object_name_;
        }
        void SetName(std::string & name)
        {
            object_name_ = name;
        }
        void Clone(const Mesh& mesh)
        {
            *this = mesh;
        }
        void MiniCopy(const Mesh& mesh)
        {
            points_ = std::move(mesh.points_);
            faces_ = std::move(mesh.faces_);
            indexs_ = std::move(mesh.indexs_);
            normals_ = std::move(mesh.normals_);
        }

        void Merge(const Mesh& mesh)
        {
            return;
        }

        void PointsReserve(size_t num)
        {
            points_.reserve(num);
            pointIndexMap_.reserve(num);
        }

        void IndexsReserve(size_t num)
        {
            faces_.reserve(num);
            indexs_.reserve(num);
        }

        size_t AddPoint(const Point & p)
        {
            const auto& iter = pointIndexMap_.find(p);
            if (iter != pointIndexMap_.end()) {
                return iter->second;
            }
            size_t size = pointIndexMap_.size();
            points_.emplace_back(std::move(p));
            pointIndexMap_.emplace(std::move(p), size);
            return size;
        }

        size_t AddFace(int v0, int v1, int v2)
        {
            std::vector<int> vertexs;
            vertexs.reserve(3);
            vertexs.emplace_back(v0);
            vertexs.emplace_back(v1);
            vertexs.emplace_back(v2);
            const size_t faceIndex = indexs_.size();
            indexs_.emplace_back(std::move(vertexs));
            faces_.emplace_back(faceIndex);
            return faceIndex;
        }

        void GenerateFaceNormals(bool Normalized = true, bool calculateArea = false)
        {
            std::vector<double> areas;
            std::vector<Point> normals;
            const auto& faces = GetFaces();
            const auto& points = GetPoints();
            const auto& indexs = GetFaceVertexAdjacency();
            normals.resize(faces.size());
            std::transform(faces.begin(), faces.end(), normals.begin(), [&](const int& f) {
                const auto& neighbor = indexs[f];
                const auto& p0 = points[neighbor[0]];
                const auto& p1 = points[neighbor[1]];
                const auto& p2 = points[neighbor[2]];
                return (p2 - p1).Cross(p0 - p1);
            });
            //std::transform针对大量数据有加速效应
            /*normals.reserve(faces.size());
            for (const auto& f : faces) {
                const auto& neighbor = indexs[f];
                const auto& p0 = points[neighbor[0]];
                const auto& p1 = points[neighbor[1]];
                const auto& p2 = points[neighbor[2]];
                const auto& n = (p2 - p1).Cross(p0 - p1);
                normals.emplace_back(std::move(n));
            }*/
            if (calculateArea) {
                areas.resize(faces.size());
                std::transform(normals.begin(), normals.end(), areas.begin(), [&](const Point& n) {
                    return n.Norm() / 2.0;
                });
                areas_.swap(areas);
            }
            if (Normalized) {
                for (auto& nor : normals) {
                    nor.Normalize();
                }
            }
            normals_.swap(normals);
        }

        void GenerateFaceAreas()
        {
            const auto& faces = GetFaces();
            const auto& points = GetPoints();
            const auto& indexs = GetFaceVertexAdjacency();
            areas_.resize(faces.size());
            std::transform(faces.begin(), faces.end(), areas_.begin(), [&](const int& f) {
                const auto& neighbor = indexs[f];
                const auto& p0 = points[neighbor[0]];
                const auto& p1 = points[neighbor[1]];
                const auto& p2 = points[neighbor[2]];
                return (p2 - p1).Cross(p0 - p1).Norm() / 2.0;
            });
        }
        void GenerateFaceEdgeAdjacency(bool bGenerateEdgeFaceAdjacency = false, bool bGenerateEgdeLength = false)
        {
            const size_t facenums = faces_.size();
            //避免indexs被改变
            std::vector<std::vector<int>> indexs;
            indexs.reserve(facenums);
            indexs.assign(indexs_.begin(), indexs_.end());
            //定义边的哈希函数
            struct EdgeHash {
                size_t operator()(const Edge& e) const
                {
                    return int((e.a * 99989)) ^ (int(e.b * 99991) << 2);
                }
            };
            //判定两条边是否相同
            struct EdgeEqual {
                bool operator()(const Edge& e1, const Edge& e2) const
                {
                    return e1.a == e2.a && e1.b == e2.b;
                }
            };
            std::unordered_map<Edge, size_t, EdgeHash, EdgeEqual> edgeIndexMap_;
            edges_.reserve(3 * facenums);
            edgeIndexMap_.reserve(3 * facenums);
            faceEdgeAdjacency_.resize(facenums);
            for (size_t i = 0; i < facenums; ++i) {
                std::vector<int>elist;
                elist.reserve(3);
                auto& vertexs = indexs[i];
                std::sort(vertexs.begin(), vertexs.end());
                // 第1 2条边
                for (size_t j = 0; j < 2; ++j) {
                    Edge e(vertexs[j], vertexs[j + 1]);
                    const auto & itr = edgeIndexMap_.find(e);
                    if (itr != edgeIndexMap_.end()) {
                        elist.emplace_back(itr->second);
                    }
                    if (itr == edgeIndexMap_.end()) {
                        const size_t edgeIndex = edgeIndexMap_.size();
                        edgeIndexMap_.emplace(std::move(e), edgeIndex);
                        edges_.emplace_back(std::move(e));
                        elist.emplace_back(edgeIndex);
                    }
                }
                // 第3条边
                Edge e(vertexs[0], vertexs[2]);
                const auto& itr = edgeIndexMap_.find(e);
                if (itr != edgeIndexMap_.end()) {
                    elist.emplace_back(itr->second);
                }
                if (itr == edgeIndexMap_.end()) {
                    const size_t edgeIndex = edgeIndexMap_.size();
                    edgeIndexMap_.emplace(std::move(e), edgeIndex);
                    edges_.emplace_back(std::move(e));
                    elist.emplace_back(edgeIndex);
                }
                faceEdgeAdjacency_[i] = std::move(elist);
            }
            if (bGenerateEdgeFaceAdjacency) {
                const size_t edgenums = edges_.size();
                edgeFaceAdjacency_.resize(edgenums);
                for (size_t i = 0; i < edgenums; ++i) {
                    edgeFaceAdjacency_[i].reserve(2);
                }
                for (size_t i = 0; i < facenums; ++i) {
                    const auto& edges = std::move(faceEdgeAdjacency_[i]);
                    for (int j = 0; j < 3; ++j) {
                        edgeFaceAdjacency_[edges[j]].emplace_back(i);
                    }
                }
            }
            if (bGenerateEgdeLength) {
                const size_t edgenums = edges_.size();
                edgeLengths_.reserve(edgenums);
                for (const auto& e : edges_) {
                    const auto& d = (points_[e.a] - points_[e.b]).Norm();
                    edgeLengths_.emplace_back(d);
                }
            }
        }
        std::vector<Point> GetPoints() const
        {
            return points_;
        }
        std::vector<Point>& GetPoints()
        {
            return points_;
        }
        std::vector<Point>& GetFaceNormals()
        {
            return normals_;
        }
        std::vector<double>& GetFaceAreas()
        {
            return areas_;
        }
        std::vector<int>& GetFaces()
        {
            return faces_;
        }
        std::vector<Edge>& GetEdges()
        {
            return edges_;
        }
        std::vector<std::vector<int>> GetFaceVertexAdjacency() const
        {
            return indexs_;
        }
        std::vector<std::vector<int>>& GetFaceVertexAdjacency()
        {
            return indexs_;
        }
        std::vector<std::vector<int>>& GetFaceEdgeAdjacency()
        {
            return faceEdgeAdjacency_;
        }
        std::vector<std::vector<int>>& GetEdgeFaceAdjacency()
        {
            return edgeFaceAdjacency_;
        }

        void Clear()
        {
            std::vector<Point>().swap(points_);
            std::vector<int>().swap(faces_);
            std::vector<Point>().swap(normals_);
            std::vector<std::vector<int>>().swap(indexs_);
            std::vector<std::vector<int>>().swap(faceEdgeAdjacency_);
            std::vector<std::vector<int>>().swap(edgeFaceAdjacency_);
            std::vector<double>().swap(edgeLengths_);
            std::vector<Edge>().swap(edges_);
        }

        inline void Translate(const Point& trans)
        {
            for (auto& pt : points_) {
                pt += trans;
            }
        }

        inline void Rotate(const Point& axis, const double& angle)
        {
            const auto& mat = Matrix3d::GetMatrix(axis, angle);
            Rotate(mat);
        }

        inline void Rotate(const Matrix3d& mat)
        {
            for (auto& p : points_) {
                mat *= p;
            }
        }

        inline void Rotate(const Point& a, const Point& b)
        {
            const auto& mat = Matrix3d::GetMatrix(a, b);
            Rotate(mat);
        }
        BoundBox Bound()const
        {
            constexpr double minValue = std::numeric_limits<double>::lowest();
            constexpr double maxValue = std::numeric_limits<double>::max();
            Point min(maxValue, maxValue, maxValue);
            Point max(minValue, minValue, minValue);
            for (const auto& p : points_) {
                if (p.x < min.x) min.x = p.x;
                if (p.y < min.y) min.y = p.y;
                if (p.z < min.z) min.z = p.z;
                if (p.x > max.x) max.x = p.x;
                if (p.y > max.y) max.y = p.y;
                if (p.z > max.z) max.z = p.z;
            }
            return BoundBox(min, max);
        }
        inline int EdgeOppositePoint(int e, int f)const
        {
            const auto& neighbor = indexs_[f];
            const auto& edge = edges_[e];
            for (int i = 0; i < 3; ++i) {
                if ((neighbor[i] != edge.a) && (neighbor[i] != edge.b)) {
                    return neighbor[i];
                }
            }
            return -1;
        }

        //默认逆时针排序好的
        inline bool EdgeJoinEdge(int e1, int e2)const
        {
            return (edges_[e1].a == edges_[e2].b) || (edges_[e1].b == edges_[e2].a);
        }

        inline bool FaceJoinFace(int f1, int f2)const
        {
            const auto& indexs = GetFaceVertexAdjacency();
            const auto& fs1 = indexs[f1];
            const auto& fs2 = indexs[f2];
            /*for (auto& it1 = fs1.begin(); it1 != fs1.end(); ++it1) {
                for (auto& it2 = fs2.begin(); it2 != fs2.end(); ++it2) {
                    if (*it1 == *it2) {
                        return true;
                    }
                }
            }*/
            //下面这种循环遍历更快
            const int n1 = fs1.size() - 1;
            const int n2 = fs2.size() - 1;
            for (int i = n1; i >= 0; --i) {
                for (int j = n2; j >= 0; --j) {
                    if (fs1[i] == fs2[j]) {
                        return true;
                    }
                }
            }
            return false;
        }

        std::vector<int> GetFaceNeighbor(const int f)
        {
            const auto& faces = GetFaces();
            std::vector<int> neighbor;
            neighbor.reserve(15);
            const int maxInx = faces.size() - 1;
            for (int fa = maxInx; fa >= 0; --fa) {
                if (fa != f) {
                    if (FaceJoinFace(f, fa)) {
                        neighbor.emplace_back(fa);
                    }
                }
            }
            return neighbor;
        }

        std::vector<std::vector<int>> GetFaceNeighborFaces()
        {
            const auto& faces = GetFaces();
            const auto& points = GetPoints();
            const auto& indexs = GetFaceVertexAdjacency();//面->顶点索引
            const int facenums = faces.size();
            std::vector<std::vector<int>> vertexFaceMap;
            vertexFaceMap.resize(points.size());
            for (int f = 0; f < facenums; ++f) {
                const auto& fs = std::move(indexs[f]);
                vertexFaceMap[fs[0]].emplace_back(f);
                vertexFaceMap[fs[1]].emplace_back(f);
                vertexFaceMap[fs[2]].emplace_back(f);
            }
            std::vector<std::vector<int>>neighbors;
            neighbors.resize(facenums);
            //#pragma omp parallel for
            for (int f = 0; f < facenums; ++f) {
                std::vector<int> neighbor;
                neighbor.reserve(20);
                const auto& vertexs = std::move(indexs[f]);
                for (const auto& v : vertexs) {
                    const auto& fs = std::move(vertexFaceMap[v]);
                    //neighbor.reserve(neighbor.size() + std::distance(fs.begin(), fs.end()));
                    for (const auto& fa : fs) {
                        if (fa != f) {
                            neighbor.emplace_back(fa);
                        }
                    }
                }
                std::sort(neighbor.begin(), neighbor.end());
                auto last = std::unique(neighbor.begin(), neighbor.end());
                neighbor.erase(last, neighbor.end());
                neighbors[f] = std::move(neighbor);
            }
            return neighbors;
        }

        std::vector<int> SelectLargetPlanar(double threshold = 0.95)
        {
            GenerateFaceNormals(true, true);
            const auto& normals = GetFaceNormals();
            const auto& faces = GetFaces();
            const auto& areas = GetFaceAreas();
            const auto& neighbors = GetFaceNeighborFaces();
            const size_t facenums = faces.size();
            std::vector<bool> masks(facenums, true);
            std::vector<std::vector<int>> selectFaces;
            selectFaces.reserve(facenums);
            for (const auto& f : faces) {
                if (masks[f]) {
                    const auto& nf = normals[f];
                    std::vector<int>currentFaces;
                    std::queue<int>currentQueue;
                    currentQueue.emplace(f);
                    currentFaces.emplace_back(f);
                    masks[f] = false;
                    while (!currentQueue.empty()) {
                        auto& fr = currentQueue.front();
                        currentQueue.pop();
                        const auto& neighbor = neighbors[fr];
                        for (const auto& fa : neighbor) {
                            if (masks[fa]) {
                                const auto& na = normals[fa];
                                const auto& nr = normals[fr];
                                if (nr.Dot(na) > threshold && nf.Dot(na) > threshold) {
                                    currentQueue.emplace(fa);
                                    currentFaces.emplace_back(fa);
                                    masks[fa] = false;
                                }
                            }
                        }
                    }
                    selectFaces.emplace_back(std::move(currentFaces));
                }
            }

            // 找到最大的平面
            double maxArea = 0.0f;
            std::vector<int> resultFaces;
            for (auto& fs : selectFaces) {
                double currentArea = 0.0;
                for (const auto& f : fs) {
                    currentArea += areas[f];
                }
                if (currentArea > maxArea) {
                    maxArea = currentArea;
                    resultFaces.swap(fs);
                }
            }
            return resultFaces;
        }
        Point FindBottomDirection(std::vector<int>* bottomfaces = nullptr, double threshold = 0.95)
        {
            *bottomfaces = SelectLargetPlanar(threshold);
            const auto& areas = GetFaceAreas();
            const auto& normals = GetFaceNormals();
            Point normal(0, 0, 0);
            for (const auto& f : *bottomfaces) {
                normal += normals[f] * areas[f];
            }
            return normal.Normalized();
        }
        //模型放平xoy平面后，底面整平
        void FlatBottomSurface(std::vector<int>* bottomfaces = nullptr)
        {
            auto& points = GetPoints();
            const auto& indexs = GetFaceVertexAdjacency();
            const auto& faceIndexs = *bottomfaces;
            const auto& nums = bottomfaces->size();
            const auto& minz = Bound().Min().z;
            for (int i = 0; i < nums; ++i) {
                const auto& neighbor = indexs[faceIndexs[i]];
                for (int j = 0; j < 3; ++j) {
                    points[neighbor[j]].z = minz;
                }
            }
        }
        void SavePointsToMesh(std::vector<int>& pointIndexs, Mesh& mesh, double radius=0.01, size_t nrows = 20, size_t ncolumns = 20)
        {
            const size_t spherenums = pointIndexs.size();
            const size_t nums = nrows * ncolumns + 2;
            const size_t pointnums = nums * spherenums;
            const size_t facenums = 2 * (nums - 2) * spherenums;
            IndexsReserve(facenums);
            auto& points = mesh.GetPoints();
            points.reserve(pointnums);
            for (int k = 0; k < spherenums; ++k) {
                const auto& p = points_[pointIndexs[k]];
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
                    mesh.AddFace(v0, i0, i1);
                    const auto& j0 = i0 + (nrows - 1) * ncolumns;
                    const auto& j1 = i1 + (nrows - 1) * ncolumns;
                    mesh.AddFace(j1, j0, maxInx);
                }
                //中间部分
                for (size_t i = 0; i < nrows - 1; ++i) {
                    const auto& j0 = i * ncolumns + 1 + v0;
                    const auto& j1 = (i + 1) * ncolumns + 1 + v0;
                    for (size_t j = 0; j < ncolumns; ++j) {
                        const auto& i0 = j0 + j;
                        const auto& i1 = j0 + (j + 1) % ncolumns;
                        const auto & i2 = j1 + j;
                        const auto & i3 = j1 + (j + 1) % ncolumns;
                        mesh.AddFace(i2, i1, i0);
                        mesh.AddFace(i2, i3, i1);
                    }
                }
            }
            mesh.GenerateFaceNormals();
        }

        void SaveEdgesToMesh(std::vector<int>& edgeIndexs, Mesh& mesh, double r = 0.01, size_t nslices = 20)
        {
            const size_t nums = edgeIndexs.size();
            auto& points = mesh.GetPoints();
            points.reserve(2 * nums * nslices);
            double delta = 2.0 * M_PI / nslices;
            Point z(0, 0, 1);
            for (size_t i = 0; i < nums; ++i) {
                const auto& a = points_[edges_[edgeIndexs[i]].a];
                const auto& b = points_[edges_[edgeIndexs[i]].b];
                const auto& n = (b - a).Normalized();
                auto x = z.Cross(n);
                if (x.Norm() < EPS) {
                    x = Point(1, 0, 0);
                }
                const auto& y = n.Cross(x);
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
            mesh.IndexsReserve(2 * nums * nslices);
            for (size_t i = 0; i < nums; ++i) {
                for (size_t j = 0; j < nslices; ++j) {
                    const auto& i0 = j + 2 * nslices * i;
                    const auto& i1 = (j + 1) % nslices + 2 * nslices * i;
                    const auto & j0 = i0 + nslices;
                    const auto & j1 = i1 + nslices;
                    mesh.AddFace(j0, i1, i0);
                    mesh.AddFace(j0, j1, i1);
                }
            }
            mesh.GenerateFaceNormals();
        }

        void SaveFacesToMesh(std::vector<int>& faceIndexs, Mesh& faceMesh)
        {
            std::vector<int> selectFaces;
            const int facenums = faces_.size();
            for (const auto& f : faceIndexs) {
                if (f < facenums) {
                    selectFaces.emplace_back(f);
                }
            }
            faceIndexs.swap(selectFaces);
            const auto& points = GetPoints();
            const auto& indexs = GetFaceVertexAdjacency();
            for (const auto& f : faceIndexs) {
                const auto& vertexs = indexs[f];
                int v0 = faceMesh.AddPoint(points[vertexs[0]]);
                int v1 = faceMesh.AddPoint(points[vertexs[1]]);
                int v2 = faceMesh.AddPoint(points[vertexs[2]]);
                faceMesh.AddFace(v0, v1, v2);
            }
            faceMesh.GenerateFaceNormals();
        }

        void DeleteFaces(std::vector<int> & faces, bool bKeepPoints = false)
        {
            std::sort(faces.begin(), faces.end());
            std::vector<int> result(faces_.size() - faces.size());
            std::set_difference(faces_.begin(), faces_.end(), faces.begin(), faces.end(), result.begin());
            const size_t facenums = result.size();
            std::vector<int> newFaces(facenums);
            std::iota(newFaces.begin(), newFaces.end(), 0);
            if (bKeepPoints) {
                std::vector<std::vector<int>> otherIndexs;
                otherIndexs.resize(facenums);
                for (int i = 0; i < facenums; ++i) {
                    otherIndexs[i] = std::move(indexs_[result[i]]);
                }
                faces_.swap(newFaces);
                indexs_.swap(otherIndexs);
                GenerateFaceNormals();
                return;
            }
            std::vector<std::vector<int>> otherIndexs;
            otherIndexs.resize(facenums);
            for (int i = 0; i < facenums; ++i) {
                otherIndexs[i] = std::move(indexs_[result[i]]);
            }
            std::unordered_map<Point, size_t, PointHash, PointEqual> pointMap;
            pointMap.rehash(2 * facenums);
            std::vector<Point> newPoints;
            std::vector<std::vector<int>>newIndexs;
            newPoints.reserve(3 * facenums);
            newIndexs.reserve(facenums);
            for (const auto& vertexs : otherIndexs) {
                std::vector<int> faceIndex;
                faceIndex.reserve(3);
                for (const auto& v : vertexs) {
                    const auto& p = points_[v];
                    const auto& iter = pointMap.find(p);
                    if (iter != pointMap.end()) {
                        faceIndex.emplace_back(iter->second);
                    } else {
                        const int newIndex = pointMap.size();
                        newPoints.emplace_back(std::move(p));
                        pointMap.emplace(std::move(p), newIndex);
                        faceIndex.emplace_back(newIndex);
                    }
                }
                newIndexs.emplace_back(std::move(faceIndex));
            }
            points_.swap(newPoints);
            indexs_.swap(newIndexs);
            faces_.swap(newFaces);
            pointIndexMap_.swap(pointMap);
            GenerateFaceNormals();
            return;
        }

        void SelectIndividualEdges(std::vector<int>& edges, bool bCounterClockWise = false)
        {
            std::vector<int>().swap(edges);
            if (edgeFaceAdjacency_.empty()) {
                GenerateFaceEdgeAdjacency(true);
            }
            const size_t edgenums = edgeFaceAdjacency_.size();
            for (int i = 0; i < edgenums; ++i) {
                const auto& neighbor = edgeFaceAdjacency_[i];
                if (neighbor.size() == 1) {
                    edges.emplace_back(i);
                }
            }
            if (bCounterClockWise) {
                auto& elist = GetEdges();
                const auto& points = GetPoints();
                const auto& indexs = GetFaceVertexAdjacency();
                for (auto& e : edges) {
                    const auto& neighbor = edgeFaceAdjacency_[e];
                    const auto& f = neighbor.front();
                    const auto& v0 = indexs[f][0];
                    const auto& v1 = indexs[f][1];
                    const auto& v2 = indexs[f][2];
                    const auto& n = (points[v2] - points[v1]).Cross(points[v0] - points[v1]).Normalized();
                    const auto& a = elist[e].a;
                    const auto& b = elist[e].b;
                    const auto& c = EdgeOppositePoint(e, f);
                    const auto& d= (points[b] - points[a]).Cross(points[c] - points[a]).Normalized();
                    if (n.Dot(d) < 0) {
                        elist[e].a = elist[e].a + elist[e].b;
                        elist[e].b = elist[e].a - elist[e].b;
                        elist[e].a = elist[e].a - elist[e].b;
                    }
                }
            }
        }

        void GetSequentialPoints(std::vector<int>& edges, std::vector<std::vector<int>>& sequentials)
        {
            auto& elist = GetEdges();
            const size_t nums = edges.size();
            std::vector<bool> masks(elist.size(), false);
            // 所有边顶点按照逆时针排序
            if (true) {
                const auto& points = GetPoints();
                const auto& indexs = GetFaceVertexAdjacency();
                for (auto& e : edges) {
                    const auto& neighbor = edgeFaceAdjacency_[e];
                    const auto& f = neighbor.front();
                    const auto& v0 = indexs[f][0];
                    const auto& v1 = indexs[f][1];
                    const auto& v2 = indexs[f][2];
                    const auto& n = (points[v2] - points[v1]).Cross(points[v0] - points[v1]).Normalized();
                    const auto & a = elist[e].a;
                    const auto & b = elist[e].b;
                    const auto & c = EdgeOppositePoint(e, f);
                    const auto & d = (points[b] - points[a]).Cross(points[c] - points[a]).Normalized();
                    if (n.Dot(d) < 0) {
                        elist[e].a = elist[e].a + elist[e].b;
                        elist[e].b = elist[e].a - elist[e].b;
                        elist[e].a = elist[e].a - elist[e].b;
                        masks[e] = true;
                    }
                }
            }
            std::queue<int> Queues;
            for (int i = 0; i < nums; ++i) {
                Queues.push(edges[i]);
            }
            //边的若干有序序列
            std::vector<std::vector<int>> edgeRings;
            while (!Queues.empty()) {
                std::vector<int> current;
                current.reserve(nums);
                const auto& e = Queues.front();
                Queues.pop();
                current.emplace_back(e);
                int count = Queues.size();
                int times = 0;
                while (!Queues.empty()) {
                    const auto& front = current.front();
                    const auto& back = current.back();
                    const auto& ef = Queues.front();
                    if (elist[ef].b == elist[front].a) {
                        current.insert(current.begin(), ef);
                        Queues.pop();
                        times = 0;
                    } else if (elist[ef].a == elist[back].b) {
                        current.emplace_back(ef);
                        Queues.pop();
                        times = 0;
                    } else {
                        Queues.pop();
                        Queues.push(ef);
                        ++times;
                    }
                    if (elist[front].a == elist[back].b) {
                        break;
                    }
                    if (times > count) {
                        break;
                    }
                }
                edgeRings.emplace_back(current);
            }
            sequentials.reserve(edgeRings.size());
            //修改的边的顶点顺序需要还原
            for (auto& es : edgeRings) {
                std::vector<int> pointList;
                pointList.reserve(es.size());
                const auto& fr = elist[es.front()].a;
                const auto& ba = elist[es.back()].b;
                if (fr != ba) {
                    pointList.emplace_back(fr);
                }
                for (auto& e : es) {
                    pointList.emplace_back(elist[e].b);
                    if (masks[e]) {
                        elist[e].a = elist[e].a + elist[e].b;
                        elist[e].b = elist[e].a - elist[e].b;
                        elist[e].a = elist[e].a - elist[e].b;
                    }
                }
                sequentials.emplace_back(pointList);
            }
        }

        ~Mesh() {}

    private:
        std::string object_name_;
        std::vector<std::vector<int>> indexs_;
        std::vector<Point> points_;
        std::vector<int> faces_;
        std::vector<Edge> edges_;
        std::vector<Point> normals_;
        std::vector<double> areas_;
        std::vector<double> edgeLengths_;
        std::vector<std::vector<int>> faceEdgeAdjacency_;
        std::vector<std::vector<int>> edgeFaceAdjacency_;
        std::vector<std::vector<int>> vertexFaceAdjacency_;


        //定义点的哈希函数
        struct PointHash {
            size_t operator()(const Point& p) const
            {
                return (int(p.x * 99971)) ^ (int(p.y * 99989) << 2) ^ (int(p.z * 99991) << 3);
            }
        };
        //判定两个点是否相同
        struct PointEqual {
            bool operator()(const Point& p1, const Point& p2) const
            {
                auto isEqual = [&](double a, double b, double eps = EPS) {
                    return std::fabs(a - b) < eps;
                };
                return isEqual(p1.x, p2.x) && isEqual(p1.y, p2.y) && isEqual(p1.z, p2.z);
            }
        };
        //std::unordered_set<Point, PointHash, PointEqual> uniquePoints;
        std::unordered_map<Point, size_t, PointHash, PointEqual> pointIndexMap_;

    };

}