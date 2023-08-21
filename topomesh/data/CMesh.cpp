#include "CMesh.h"
#include <numeric>
#include <unordered_map>
#include <queue>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

#ifndef EPS
#define EPS 1E-8f
#endif // !EPS
namespace topomesh {
    CMesh::CMesh()
    {
        mbox.clear();
        ::std::vector<PPoint>().swap(mpoints);
        ::std::vector<FFace>().swap(mfaces);
        ::std::vector<EEdge>().swap(medges);
        ::std::vector<PPoint>().swap(mnorms);
        ::std::vector<float>().swap(mareas);
        ::std::vector<float>().swap(medgeLengths);
        ::std::vector<std::vector<int>>().swap(mfaceEdges);
        ::std::vector<std::vector<int>>().swap(medgeFaces);
    }
    CMesh::~CMesh()
    {
        Clear();
    }
    void CMesh::Clear()
    {
        mbox.clear();
        std::vector<PPoint>().swap(mpoints);
        std::vector<FFace>().swap(mfaces);
        std::vector<PPoint>().swap(mnorms);
        std::vector<std::vector<int>>().swap(mfaceEdges);
        std::vector<std::vector<int>>().swap(medgeFaces);
        std::vector<float>().swap(medgeLengths);
        std::vector<EEdge>().swap(medges);
        std::vector<float>().swap(mareas);
    }
    bool CMesh::ReadFromSTL(const char* filename)
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
        if (format == "ascii")//判断格式
        {
            // ReadASCII
            fid.open(filename, std::ios::in);
            // Read the STL header
            std::string header;
            std::getline(fid, header);
            std::string linestr;
            size_t np1, np2;
            float x, y, z;
            mfaces.reserve(facenums);
            mnorms.reserve(facenums);
            mpoints.reserve(3 * facenums);
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
                    mnorms.emplace_back(x, y, z);
                    mfaces.emplace_back();
                    const auto& index = 3 * facenums;
                    mfaces.back().x = index;
                    mfaces.back().y = index + 1;
                    mfaces.back().z = index + 2;
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
                    mpoints.emplace_back(std::move(PPoint(x, y, z)));
                }
            }
            fid.close();
        } else if (format == "binary") {
            //ReadBinary;
            char objectname[80];
            fid.open(filename, std::ios::in | std::ios::binary);
            fid.read(reinterpret_cast<char*>(&objectname), sizeof(objectname));
            //fid.seekg(80, std::ios::beg);
            fid.read(reinterpret_cast<char*>(&facenums), sizeof(int));
            mpoints.resize(3 * facenums);
            mnorms.resize(facenums);
            //mfaces.resize(facenums);
            for (size_t i = 0; i < facenums; ++i) {
                float pdata[12];
                fid.read(reinterpret_cast<char*>(&pdata), sizeof(float) * 12);
                //读取法向量信息
                float* pbuffer = pdata;
                for (size_t j = 0; j < 3; ++j) {
                    mnorms[i][j] = *pbuffer;
                    ++pbuffer;
                }
                //读取三个顶点坐标
                for (size_t j = 0; j < 3; ++j) {
                    for (size_t k = 0; k < 3; ++k) {
                        mpoints[3 * i + j][k] = *pbuffer;
                        ++pbuffer;
                    }
                }
                //const int& i0 = 3 * i;
                //mfaces[i] = std::move(std::vector<int>({ i0,i0 + 1,i0 + 2 }));
                char tailbuffer[2];
                fid.read(reinterpret_cast<char*>(&tailbuffer), sizeof(tailbuffer));
            }
            fid.close();
        } else {
            std::cout << "format error!" << std::endl;
            return false;
        }
        return true;
    }
    bool CMesh::WriteSTLFile(const char* filename, bool bBinary)
    {
        int face_num = mfaces.size();
        if (mnorms.empty()) {
            GenerateFaceNormals();
        }
        int normal_num = mnorms.size();
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
                int pIndex0 = mfaces[i][0];
                int pIndex1 = mfaces[i][1];
                int pIndex2 = mfaces[i][2];
                float nx = mnorms[i].x;
                float ny = mnorms[i].y;
                float nz = mnorms[i].z;
                //保存法向量
                fs.write((char*)(&nx), floatSize);
                fs.write((char*)(&ny), floatSize);
                fs.write((char*)(&nz), floatSize);
                // 保存顶点
                float p0x = mpoints[pIndex0].x;
                float p0y = mpoints[pIndex0].y;
                float p0z = mpoints[pIndex0].z;
                fs.write((char*)(&p0x), floatSize);
                fs.write((char*)(&p0y), floatSize);
                fs.write((char*)(&p0z), floatSize);

                float p1x = mpoints[pIndex1].x;
                float p1y = mpoints[pIndex1].y;
                float p1z = mpoints[pIndex1].z;
                fs.write((char*)(&p1x), floatSize);
                fs.write((char*)(&p1y), floatSize);
                fs.write((char*)(&p1z), floatSize);

                float p2x = mpoints[pIndex2].x;
                float p2y = mpoints[pIndex2].y;
                float p2z = mpoints[pIndex2].z;
                fs.write((char*)(&p2x), floatSize);
                fs.write((char*)(&p2y), floatSize);
                fs.write((char*)(&p2z), floatSize);
                fs.write(a, a_size);
            }
            fs.close();
            return true;
        } else {
            //mpoints:模型顶点
            //faces_:三角面片
            //mnorms:法线
            if (face_num <= 1e6) {
                std::ofstream fs(std::string(filename) + "_ast.stl");
                if (!fs) { fs.close(); return false; }
                fs << "solid WRAP" << std::endl;
                for (int f = 0; f < face_num; ++f) {
                    fs << "facet normal ";
                    const auto& n = mnorms[f];
                    fs << (float)n.x << " " << (float)n.y << " " << (float)n.z << std::endl;
                    fs << "outer loop" << std::endl;
                    const auto& vertexs = mfaces[f];
                    for (const auto& inx : vertexs) {
                        fs << "vertex ";
                        const auto& p = mpoints[inx];
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
            const auto & normals = mnorms;
            const auto & points = mpoints;
            const auto & indexs = mfaces;
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
                            const auto& n = normals[j];
                            buffer += "facet normal " + std::to_string(n.x) + " " + std::to_string(n.y) + " " + std::to_string(n.z) + "\n";
                            buffer += "outer loop\n";
                            const auto& neighbors = indexs[j];
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
    CMesh::CMesh(const trimesh::TriMesh* mesh)
    {
        trimesh::TriMesh trimesh = *mesh;
        mfaces.swap(trimesh.faces);
        mpoints.swap(trimesh.vertices);
        std::swap(mbox, trimesh.bbox);
    }
    trimesh::TriMesh CMesh::GetTriMesh() const
    {
        trimesh::TriMesh mesh;
        mesh.vertices = mpoints;
        mesh.faces = mfaces;
        mesh.bbox = mbox;
        return mesh;
    }
    void CMesh::Merge(const CMesh& mesh)
    {
    }
    void CMesh::Clone(const CMesh& mesh)
    {
        *this = mesh;
    }
    void CMesh::MiniCopy(const CMesh& mesh)
    {
        mbox = std::move(mesh.mbox);
        mfaces = std::move(mesh.mfaces);
        mpoints = std::move(mesh.mpoints);
    }
    void CMesh::Translate(const PPoint& trans)
    {
        for (auto& pt : mpoints) {
            pt += trans;
        }
    }
    void CMesh::Rotate(const trimesh::fxform& mat)
    {
        for (auto& pt : mpoints) {
            pt = mat * pt;
        }
    }
    void CMesh::Rotate(const PPoint& axis, const double& angle)
    {
        trimesh::fxform&& mat = trimesh::fxform::rot(angle, axis);
        Rotate(mat);
    }
    void CMesh::Rotate(const PPoint& dir1, const PPoint& dir2)
    {
        trimesh::fxform&& mat = trimesh::fxform::rot_into(dir1, dir2);
        Rotate(mat);
    }
    void CMesh::BuildFromBox(const BBox& box)
    {
        const trimesh::point& min = box.min;
        const trimesh::point& max = box.max;
        const auto& diag = max - min;
        AddPoint(min + trimesh::point(diag.x, 0, 0));
        AddPoint(min + trimesh::point(diag.x, diag.y, 0));
        AddPoint(min + trimesh::point(diag.x, diag.y, diag.z));
        AddPoint(min + trimesh::point(diag.x, 0, diag.z));
        AddPoint(min + trimesh::point(0, 0, 0));
        AddPoint(min + trimesh::point(0, diag.y, 0));
        AddPoint(min + trimesh::point(0, diag.y, diag.z));
        AddPoint(min + trimesh::point(0, 0, diag.z));
        AddFace(0, 1, 2);
        AddFace(0, 2, 3);
        AddFace(1, 5, 6);
        AddFace(1, 6, 2);
        AddFace(3, 2, 6);
        AddFace(3, 6, 7);
        AddFace(4, 0, 3);
        AddFace(4, 3, 7);
        AddFace(4, 5, 1);
        AddFace(4, 1, 0);
        AddFace(5, 4, 7);
        AddFace(5, 7, 6);
    }
    void CMesh::GenerateBoundBox()
    {
        if (!mbox.valid) {
            constexpr float max = std::numeric_limits<float>::max();
            constexpr float min = std::numeric_limits<float>::lowest();
            trimesh::vec3 a(max, max, max), b(min, min, min);
            for (const auto& p : mpoints) {
                if (a.x > p.x) a.x = p.x;
                if (a.y > p.y) a.y = p.y;
                if (a.z > p.z) a.z = p.z;
                if (b.x < p.x) b.x = p.x;
                if (b.y < p.y) b.y = p.y;
                if (b.z < p.z) b.z = p.z;
            }
            mbox = trimesh::box3(a, b);
            mbox.valid = true;
        }
    }
    int CMesh::AddPoint(const PPoint& p)
    {
        mpoints.emplace_back(std::move(p));
        return mpoints.size() - 1;
    }
    int CMesh::AddFace(int i0, int i1, int i2)
    {
        FFace f(i0, i1, i2);
        mfaces.emplace_back(std::move(f));
        return mfaces.size() - 1;
    }

    int CMesh::EdgeOppositePoint(int e, int f) const
    {
        const auto& neighbor = mfaces[f];
        const auto& edge = medges[e];
        for (int i = 0; i < 3; ++i) {
            if ((neighbor[i] != edge.a) && (neighbor[i] != edge.b)) {
                return neighbor[i];
            }
        }
        return -1;
    }

    void CMesh::GenerateFaceAreas(bool calculateAgain)
    {
        if (mareas.empty() || calculateAgain) {
            const int facenums = mfaces.size();
            mareas.resize(facenums);
            std::transform(mfaces.begin(), mfaces.end(), mareas.begin(), [&](const FFace & f) {
                const auto& p0 = mpoints[f[0]];
                const auto& p1 = mpoints[f[1]];
                const auto& p2 = mpoints[f[2]];
                return trimesh::len((p2 - p1).cross(p0 - p1)) / 2.0;
            });
        }
    }

    void CMesh::GenerateFaceNormals(bool Normalized, bool calculateArea)
    {
        const int facenums = mfaces.size();
        ::std::vector<PPoint>().swap(mnorms);
        mnorms.reserve(facenums);
        for (const auto& f : mfaces) {
            const auto& p0 = mpoints[f[0]];
            const auto& p1 = mpoints[f[1]];
            const auto& p2 = mpoints[f[2]];
            const auto& n = (p2 - p1).cross(p0 - p1);
            mnorms.emplace_back(std::move(n));
        }
        //std::transform针对大量数据有加速效应
        if (calculateArea) {
            mareas.reserve(facenums);
            for (const trimesh::vec3& n : mnorms) {
                mareas.emplace_back(trimesh::len(n) / 2.0);
            }
        }
        if (Normalized) {
            for (auto& n : mnorms) {
                trimesh::normalize(n);
            }
        }
    }

    std::vector<std::vector<int>> CMesh::GenerateFaceNeighborFaces()
    {
        const int facenums = mfaces.size();
        std::vector<std::vector<int>> vertexFaceMap;
        vertexFaceMap.resize(mpoints.size());
        for (int f = 0; f < facenums; ++f) {
            const auto& fs = std::move(mfaces[f]);
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
            const auto& vertexs = std::move(mfaces[f]);
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

    void CMesh::GenerateFaceEdgeAdjacency(bool bGenerateEdgeFaceAdjacency, bool bGenerateEgdeLength)
    {
        const size_t facenums = mfaces.size();
        medges.reserve(3 * facenums);
        std::vector<EdgeToFace> edgesmap;
        edgesmap.reserve(3 * facenums);
        for (size_t facet_idx = 0; facet_idx < facenums; ++facet_idx) {
            for (size_t j = 0; j < 3; ++j) {
                edgesmap.push_back({});
                EdgeToFace& e2f = edgesmap.back();
                e2f.vertex_low = mfaces[facet_idx][j];
                e2f.vertex_high = mfaces[facet_idx][(j + 1) % 3];
                e2f.at_face = facet_idx;
                e2f.which_edge = j + 1;
                if (e2f.vertex_low > e2f.vertex_high) {
                    std::swap(e2f.vertex_low, e2f.vertex_high);
                    e2f.which_edge = -e2f.which_edge;
                }
            }
        }
        std::sort(edgesmap.begin(), edgesmap.end());
        std::vector<std::vector<int>>fedges(facenums, std::vector<int>(3, -1));
        int num_edges = 0;
        for (size_t i = 0; i < edgesmap.size(); ++i) {
            EdgeToFace& efi = edgesmap[i];
            const auto& fi = efi.at_face;
            if (fi == -1)
                continue;
            size_t j;
            bool found = false;
            for (j = i + 1; j < edgesmap.size() && efi == edgesmap[j]; ++j)
                if (efi.which_edge * edgesmap[j].which_edge < 0 && edgesmap[j].at_face != -1) {
                    found = true;
                    break;
                }
            if (!found) {
                for (j = i + 1; j < edgesmap.size() && efi == edgesmap[j]; ++j)
                    if (edgesmap[j].at_face != -1) {
                        found = true;
                        break;
                    }
            }
            const auto& ei = efi.which_edge;
            fedges[fi][std::fabs(ei) - 1] = num_edges;
            medges.emplace_back(std::move(EEdge(efi.vertex_low, efi.vertex_high)));
            if (found) {
                EdgeToFace& efj = edgesmap[j];
                const auto& fj = efj.at_face;
                const auto& ej = efj.which_edge;
                fedges[fj][std::fabs(ej) - 1] = num_edges;
                medges.emplace_back(std::move(EEdge(efj.vertex_low, efj.vertex_high)));
                efj.at_face = -1;
            }
            ++num_edges;
        }
        auto itr = std::unique(medges.begin(), medges.end());
        medges.erase(itr, medges.end());
        if (bGenerateEdgeFaceAdjacency) {
            const size_t edgenums = medges.size();
            std::vector<std::vector<int>> efaces;
            efaces.resize(edgenums);
            for (int i = 0; i < edgenums; ++i) {
                efaces[i].reserve(2);
            }
            for (size_t i = 0; i < facenums; ++i) {
                const auto& es = fedges[i];
                for (const auto& e : es) {
                    efaces[e].emplace_back(i);
                }
            }
            medgeFaces.swap(efaces);
        }
        mfaceEdges.swap(fedges);
        if (bGenerateEgdeLength) {
            const size_t edgenums = medges.size();
            medgeLengths.reserve(edgenums);
            for (const auto& e : medges) {
                const auto& d = trimesh::len(mpoints[e.a] - mpoints[e.b]);
                medgeLengths.emplace_back(d);
            }
        }
    }

    void CMesh::GenerateFaceEdgeAdjacency2(bool bGenerateEdgeFaceAdjacency, bool bGenerateEgdeLength)
    {
        //避免indexs被改变
        std::vector<FFace> indexs(mfaces.begin(), mfaces.end());
        const size_t facenums = indexs.size();
        //定义边的哈希函数
        struct EEdgeHash {
            size_t operator()(const EEdge& e) const
            {
                return int((e.a * 99989)) ^ (int(e.b * 99991) << 2);
            }
        };
        //判定两条边是否相同
        struct EEdgeEqual {
            bool operator()(const EEdge& e1, const EEdge& e2) const
            {
                return e1.a == e2.a && e1.b == e2.b;
            }
        };
        std::unordered_map<EEdge, size_t, EEdgeHash, EEdgeEqual> edgeIndexMap;
        medges.reserve(3 * facenums);
        edgeIndexMap.reserve(3 * facenums);
        mfaceEdges.resize(facenums);
        for (size_t i = 0; i < facenums; ++i) {
            std::vector<int>elist;
            elist.reserve(3);
            auto& vertexs = indexs[i];
            std::sort(vertexs.begin(), vertexs.end());
            // 第1 2条边
            for (size_t j = 0; j < 2; ++j) {
                EEdge e(vertexs[j], vertexs[j + 1]);
                const auto & itr = edgeIndexMap.find(e);
                if (itr != edgeIndexMap.end()) {
                    elist.emplace_back(itr->second);
                }
                if (itr == edgeIndexMap.end()) {
                    const size_t edgeIndex = edgeIndexMap.size();
                    edgeIndexMap.emplace(std::move(e), edgeIndex);
                    medges.emplace_back(std::move(e));
                    elist.emplace_back(edgeIndex);
                }
            }
            // 第3条边
            EEdge e(vertexs[0], vertexs[2]);
            const auto& itr = edgeIndexMap.find(e);
            if (itr != edgeIndexMap.end()) {
                elist.emplace_back(itr->second);
            }
            if (itr == edgeIndexMap.end()) {
                const size_t edgeIndex = edgeIndexMap.size();
                edgeIndexMap.emplace(std::move(e), edgeIndex);
                medges.emplace_back(std::move(e));
                elist.emplace_back(edgeIndex);
            }
            mfaceEdges[i] = std::move(elist);
        }
        if (bGenerateEdgeFaceAdjacency) {
            const size_t edgenums = medges.size();
            medgeFaces.resize(edgenums);
            for (size_t i = 0; i < edgenums; ++i) {
                medgeFaces[i].reserve(2);
            }
            for (size_t i = 0; i < facenums; ++i) {
                const auto& es = std::move(mfaceEdges[i]);
                for (int j = 0; j < 3; ++j) {
                    medgeFaces[es[j]].emplace_back(i);
                }
            }
        }
        if (bGenerateEgdeLength) {
            const size_t edgenums = medges.size();
            medgeLengths.reserve(edgenums);
            for (const auto& e : medges) {
                const auto& d = trimesh::len(mpoints[e.a] - mpoints[e.b]);
                medgeLengths.emplace_back(d);
            }
        }
    }

    std::vector<int> CMesh::SelectLargetPlanar(float threshold)
    {
        GenerateFaceNormals(true, true);
        const size_t facenums = mfaces.size();
        std::vector<int> sequence(facenums);
        for (int i = 0; i < facenums; ++i) {
            sequence[i] = i;
        }
        trimesh::TriMesh trimesh = GetTriMesh();
        trimesh.need_across_edge();
        mffaces = trimesh.across_edge;
        std::vector<int> masks(facenums, 1);
        std::vector<std::vector<int>> selectFaces;
        selectFaces.reserve(facenums);
        for (const auto& f : sequence) {
            if (masks[f]) {
                const auto& nf = mnorms[f];
                std::vector<int>currentFaces;
                std::queue<int>currentQueue;
                currentQueue.emplace(f);
                currentFaces.emplace_back(f);
                masks[f] = false;
                while (!currentQueue.empty()) {
                    auto& fr = currentQueue.front();
                    currentQueue.pop();
                    const FFace& neighbor = mffaces[fr];
                    for (const auto& fa : neighbor) {
                        if (masks[fa]) {
                            const auto& na = mnorms[fa];
                            const auto& nr = mnorms[fr];
                            if ((nr DOT na) > threshold && (nf DOT na) > threshold) {
                                currentQueue.emplace(fa);
                                currentFaces.emplace_back(fa);
                                masks[fa] = 0;
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
                currentArea += mareas[f];
            }
            if (currentArea > maxArea) {
                maxArea = currentArea;
                resultFaces.swap(fs);
            }
        }
        return resultFaces;
    }

    void CMesh::FlatBottomSurface(std::vector<int>* bottomfaces)
    {
        GenerateBoundBox();
        const auto& minz = mbox.min.z;
        for (const auto& f : *bottomfaces) {
            const auto& neighbor = mfaces[f];
            for (int j = 0; j < 3; ++j) {
                mpoints[neighbor[j]].z = minz;
            }
        }
    }

    CMesh::PPoint CMesh::FindBottomDirection(std::vector<int>* bottomfaces, float threshold)
    {
        *bottomfaces = SelectLargetPlanar(threshold);
        PPoint normal(0, 0, 0);
        for (const auto& f : *bottomfaces) {
            normal += mnorms[f] * mareas[f];
        }
        return trimesh::normalized(normal);
    }

    void CMesh::DeleteFaces(std::vector<int>& faceIndexs, bool bKeepPoints)
    {
        std::sort(faceIndexs.begin(), faceIndexs.end());
        std::vector<int> sequence(mfaces.size());
        for (int i = 0; i < mfaces.size(); ++i) {
            sequence[i] = i;
        }
        std::vector<int> result(sequence.size() - faceIndexs.size());
        std::set_difference(sequence.begin(), sequence.end(), faceIndexs.begin(), faceIndexs.end(), result.begin());
        const size_t facenums = result.size();
        std::vector<FFace> otherFaces;
        otherFaces.resize(facenums);
        for (int i = 0; i < facenums; ++i) {
            otherFaces[i] = std::move(mfaces[result[i]]);
        }
        if (bKeepPoints) {
            mfaces = std::move(otherFaces);
            GenerateFaceNormals();
            return;
        }
        //定义点的哈希函数
        struct PPointHash {
            size_t operator()(const PPoint& p) const
            {
                return (int(p.x * 99971)) ^ (int(p.y * 99989) << 2) ^ (int(p.z * 99991) << 3);
            }
        };
        //判定两个点是否相同
        struct PPointEqual {
            bool operator()(const PPoint& p1, const PPoint& p2) const
            {
                auto isEqual = [&](float a, float b, float eps = EPS) {
                    return std::fabs(a - b) < eps;
                };
                return isEqual(p1.x, p2.x) && isEqual(p1.y, p2.y) && isEqual(p1.z, p2.z);
            }
        };
        std::unordered_map<PPoint, int, PPointHash, PPointEqual> pointMap;
        pointMap.rehash(2 * facenums);
        std::vector<PPoint> newPoints;
        std::vector<FFace>newFaces;
        newPoints.reserve(3 * facenums);
        newFaces.reserve(facenums);
        for (const auto& vs : otherFaces) {
            std::vector<int> f;
            f.reserve(3);
            for (const auto& v : vs) {
                const auto& p = mpoints[v];
                const auto& iter = pointMap.find(p);
                if (iter != pointMap.end()) {
                    f.emplace_back(iter->second);
                } else {
                    const int newIndex = pointMap.size();
                    newPoints.emplace_back(std::move(p));
                    pointMap.emplace(std::move(p), newIndex);
                    f.emplace_back(newIndex);
                }
            }
            newFaces.emplace_back(std::move(f));
        }
        mpoints.swap(newPoints);
        mfaces.swap(newFaces);
        GenerateFaceNormals();
        return;
    }

    void CMesh::SelectIndividualEdges(std::vector<int>& edgeIndexs, bool bCounterClockWise)
    {
        std::vector<int>().swap(edgeIndexs);
        if (medgeFaces.empty()) {
            GenerateFaceEdgeAdjacency(true);
        }
        const size_t edgenums = medgeFaces.size();
        for (int i = 0; i < edgenums; ++i) {
            const auto& neighbor = medgeFaces[i];
            if (neighbor.size() == 1) {
                edgeIndexs.emplace_back(i);
            }
        }
        if (bCounterClockWise) {
            for (auto& e : edgeIndexs) {
                const auto& neighbor = medgeFaces[e];
                const auto& f = neighbor.front();
                const auto& v0 = mfaces[f][0];
                const auto& v1 = mfaces[f][1];
                const auto& v2 = mfaces[f][2];
                const auto& n = trimesh::normalized((mpoints[v2] - mpoints[v1]).cross(mpoints[v0] - mpoints[v1]));
                const auto & a = medges[e].a;
                const auto & b = medges[e].b;
                const auto & c = EdgeOppositePoint(e, f);
                const auto & d = trimesh::normalized((mpoints[b] - mpoints[a]).cross(mpoints[c] - mpoints[a]));
                if ((n DOT d) < 0) {
                    medges[e].a = medges[e].a + medges[e].b;
                    medges[e].b = medges[e].a - medges[e].b;
                    medges[e].a = medges[e].a - medges[e].b;
                }
            }
        }
    }

    void CMesh::GetSequentialPoints(std::vector<int>& edgeIndexs, std::vector<std::vector<int>>& sequentials)
    {
        const size_t nums = edgeIndexs.size();
        std::vector<bool> masks(medges.size(), false);
        // 所有边顶点按照逆时针排序
        if (true) {
            for (auto& e : edgeIndexs) {
                const auto& neighbor = medgeFaces[e];
                const auto& f = neighbor.front();
                const auto& v0 = mfaces[f][0];
                const auto& v1 = mfaces[f][1];
                const auto& v2 = mfaces[f][2];
                const auto& n = trimesh::normalized((mpoints[v2] - mpoints[v1]).cross(mpoints[v0] - mpoints[v1]));
                const auto & a = medges[e].a;
                const auto & b = medges[e].b;
                const auto & c = EdgeOppositePoint(e, f);
                const auto & d = trimesh::normalized((mpoints[b] - mpoints[a]).cross(mpoints[c] - mpoints[a]));
                if ((n DOT d) < 0) {
                    medges[e].a = medges[e].a + medges[e].b;
                    medges[e].b = medges[e].a - medges[e].b;
                    medges[e].a = medges[e].a - medges[e].b;
                    masks[e] = true;
                }
            }
        }
        std::queue<int> Queues;
        for (int i = 0; i < nums; ++i) {
            Queues.push(edgeIndexs[i]);
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
                if (medges[ef].b == medges[front].a) {
                    current.insert(current.begin(), ef);
                    Queues.pop();
                    times = 0;
                } else if (medges[ef].a == medges[back].b) {
                    current.emplace_back(ef);
                    Queues.pop();
                    times = 0;
                } else {
                    Queues.pop();
                    Queues.push(ef);
                    ++times;
                }
                if (medges[front].a == medges[back].b) {
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
            const auto& fr = medges[es.front()].a;
            const auto& ba = medges[es.back()].b;
            if (fr != ba) {
                pointList.emplace_back(fr);
            }
            for (auto& e : es) {
                pointList.emplace_back(medges[e].b);
                if (masks[e]) {
                    medges[e].a = medges[e].a + medges[e].b;
                    medges[e].b = medges[e].a - medges[e].b;
                    medges[e].a = medges[e].a - medges[e].b;
                }
            }
            sequentials.emplace_back(pointList);
        }
    }

    void CMesh::SavePointsToMesh(std::vector<int>& pointIndexs, CMesh& mesh, double radius, size_t nrows, size_t ncolumns)
    {
        const size_t spherenums = pointIndexs.size();
        const size_t nums = nrows * ncolumns + 2;
        const size_t pointnums = nums * spherenums;
        const size_t facenums = 2 * (nums - 2) * spherenums;
        auto& points = mesh.mpoints;
        points.reserve(pointnums);
        for (int k = 0; k < spherenums; ++k) {
            const auto& p = mpoints[pointIndexs[k]];
            points.emplace_back(p.x, p.y, p.z + radius);
            for (int i = 0; i < nrows; ++i) {
                const auto& phi = M_PI * (i + 1.0) / double(nrows + 1.0);
                const auto & z = radius * std::cos(phi);
                const auto & r = radius * std::sin(phi);
                for (int j = 0; j < ncolumns; ++j) {
                    const auto& theta = 2.0 * M_PI * j / ncolumns;
                    const auto& x = r * std::cos(theta);
                    const auto& y = r * std::sin(theta);
                    points.emplace_back(p.x + x, p.y + y, p.z + z);
                }
            }
            points.emplace_back(p.x, p.y, p.z - radius);
            const auto & maxInx = points.size() - 1;
            const auto & v0 = k * nums;
            //上下底部两部分
            for (size_t i = 0; i < ncolumns; ++i) {
                const auto& i0 = i + 1 + v0;
                const auto& i1 = (i + 1) % ncolumns + 1 + v0;
                mesh.AddFace(v0, i0, i1);
                const auto & j0 = i0 + (nrows - 1) * ncolumns;
                const auto & j1 = i1 + (nrows - 1) * ncolumns;
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
    }
    void CMesh::SaveEdgesToMesh(std::vector<int>& edgeIndexs, CMesh& mesh, double r, size_t nslices)
    {
        const size_t nums = edgeIndexs.size();
        auto& points = mesh.mpoints;
        points.reserve(2 * nums * nslices);
        double delta = 2.0 * M_PI / nslices;
        trimesh::vec3 z(0, 0, 1);
        for (size_t i = 0; i < nums; ++i) {
            const auto& a = points[medges[edgeIndexs[i]].a];
            const auto& b = points[medges[edgeIndexs[i]].b];
            const auto& n = trimesh::normalized(b - a);
            auto&& x = std::move(z.cross(n));
            if (trimesh::len(x) < EPS) {
                x = PPoint(1, 0, 0);
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
    }
    void CMesh::SaveFacesToMesh(std::vector<int>& faceIndexs, CMesh& faceMesh)
    {
        std::vector<int> selectFaces;
        const int facenums = mfaces.size();
        for (const auto& f : faceIndexs) {
            if (f < facenums) {
                selectFaces.emplace_back(f);
            }
        }
        faceIndexs.swap(selectFaces);
        const auto& points = mpoints;
        const auto& indexs = mfaces;
        for (const auto& f : faceIndexs) {
            const auto& vertexs = indexs[f];
            int v0 = faceMesh.AddPoint(points[vertexs[0]]);
            int v1 = faceMesh.AddPoint(points[vertexs[1]]);
            int v2 = faceMesh.AddPoint(points[vertexs[2]]);
            faceMesh.AddFace(v0, v1, v2);
        }
    }
}

