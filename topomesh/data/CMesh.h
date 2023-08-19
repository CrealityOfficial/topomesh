#pragma once
#include "trimesh2/TriMesh.h"
#include "trimesh2/XForm.h"


namespace topomesh {
    struct EEdge {
        int a, b; ///< a < b
        EEdge(int a0, int b0) :a(a0), b(b0) {}
        inline bool operator<(const EEdge& e) const
        {
            return a < e.a || (a == e.a && b < e.b);
        };
        inline bool operator==(const EEdge & e) const
        {
            return a == e.a && b == e.b;
        };
    };
    struct EdgeToFace {
        int vertex_low;
        int vertex_high;
        int at_face;
        int which_edge;
        bool operator==(const EdgeToFace& e2f)const
        {
            return vertex_low == e2f.vertex_low && vertex_high == e2f.vertex_high;
        }
        bool operator<(const EdgeToFace& e2f)const
        {
            return vertex_low < e2f.vertex_low || (vertex_low == e2f.vertex_low && vertex_high < e2f.vertex_high);
        }
    };
    class CMesh{
    public:
        trimesh::TriMesh mtrimesh;
        typedef trimesh::Vec<3, float> PPoint;
        typedef trimesh::Vec<3, int> FFace;
        typedef trimesh::Box<3, float>BBox;
        BBox mbox;
        ::std::vector<PPoint> mpoints;
        ::std::vector<FFace> mfaces;
        ::std::vector<FFace> mffaces;
        ::std::vector<EEdge> medges;
        ::std::vector<PPoint> mnorms;
        ::std::vector<float> mareas;
        ::std::vector<float> medgeLengths;
        ::std::vector<std::vector<int>> mfaceEdges;
        ::std::vector<std::vector<int>> medgeFaces;
        CMesh();
        ~CMesh();
        void Clear();
        bool ReadFromSTL(const char* filename);
        void DuplicateSTL(double ratio = 0.3);
        //bBinary=true,д��������ļ���bBinary=falseд��ASCII�ļ�
        bool WriteSTLFile(const char* filename, bool bBinary = true);
        CMesh(trimesh::TriMesh* mesh);
        void Merge(const CMesh& mesh);
        void Clone(const CMesh& mesh);
        void MiniCopy(const CMesh& mesh);
        void Translate(const PPoint& trans);
        void Rotate(const trimesh::fxform& mat);
        void Rotate(const PPoint& axis, const double& angle);
        void Rotate(const PPoint& dir1, const PPoint& dir2);
        void BuildFromBox(const BBox& box);
        BBox Bound();
        int AddPoint(const PPoint& p);
        int AddFace(int i0, int i1, int i2);
        int EdgeOppositePoint(int e, int f)const;
        void GenerateFaceAreas(bool calculateAgain = false);
        void GenerateFaceNormals(bool Normalized = true, bool calculateArea = false);
        std::vector<std::vector<int>> GenerateFaceNeighborFaces();
        void GenerateFaceEdgeAdjacency(bool bGenerateEdgeFaceAdjacency = false, bool bGenerateEgdeLength = false);
        void GenerateFaceEdgeAdjacency2(bool bGenerateEdgeFaceAdjacency = false, bool bGenerateEgdeLength = false);

        std::vector<int> SelectLargetPlanar(float threshold = 0.95f);
        void FlatBottomSurface(std::vector<int>* bottomfaces = nullptr);
        PPoint FindBottomDirection(std::vector<int>* bottomfaces = nullptr, float threshold = 0.95f);
        void DeleteFaces(std::vector<int>& faceIndexs, bool bKeepPoints = false);
        void SelectIndividualEdges(std::vector<int>& edgeIndexs, bool bCounterClockWise = false);
        void GetSequentialPoints(std::vector<int>& edgeIndexs, std::vector<std::vector<int>>& sequentials);

        void SavePointsToMesh(std::vector<int>& pointIndexs, CMesh& mesh, double r = 0.01, size_t nrows = 20, size_t ncolumns = 20);
        void SaveEdgesToMesh(std::vector<int>& edgeIndexs, CMesh& mesh, double r = 0.01, size_t nslices = 20);
        void SaveFacesToMesh(std::vector<int>& faceIndexs, CMesh& faceMesh);

    };
}
