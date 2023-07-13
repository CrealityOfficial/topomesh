#pragma once

#include <vector>
#include <istream>
#include <iostream>
#include "Point.h"
namespace honeycomb {
    class Matrix2d {
    public:
        const int rows = 2;
        const int cols = 2;
        double ma[2][2] = { 0 };
        Matrix2d()
        {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    ma[i][j] = 0.0;
                }
            }
        }
        Matrix2d(const double& val)
        {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    ma[i][j] = val;
                }
            }
        }
        Matrix2d(std::vector<double>& arr)
        {
            size_t size = arr.size();
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    ma[i][j] = (rows * i + j < size) ? arr[rows * i + j] : 0;
                }
            }
        }
        Matrix2d(double arr[], int size = 4)
        {
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    ma[i][j] = arr[rows * i + j];
                }
            }
        }
        Matrix2d(double arr[][2], int row_nums = 2)
        {
            for (int i = 0; i < row_nums; ++i) {
                for (int j = 0; j < cols; ++j) {
                    ma[i][j] = arr[i][j];
                }
            }
        }
        Matrix2d(double a0, double a1, double a2, double a3)
        {
            ma[0][0] = a0;
            ma[0][1] = a1;
            ma[1][0] = a2;
            ma[1][1] = a3;
        }
        bool IsIdentify() const
        {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    if (i == j) {
                        if (std::fabs(ma[i][j] - 1) > EPS) {
                            return false;
                        }
                    } else {
                        if (std::fabs(ma[i][j]) > EPS) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
        static Matrix2d Identify()
        {
            return Matrix2d(1, 0, 0, 1);
        }
        inline Point2d Row(int row)const
        {
            return Point2d(ma[row][0], ma[row][1]);
        }
        inline Point2d Col(int col)const
        {
            return Point2d(ma[0][col], ma[1][col]);
        }
        //����ʽ,�˴�û��������ʽ
        inline double Minors(int row, int col)const
        {
            if (row < 0 || row >= rows) {
                std::abort();
            }
            if (col < 0 || col >= cols) {
                std::abort();
            }
            for (int i = 0; i < 2; ++i) {
                if (i == row) continue;
                for (int j = 0; j < 2; ++j) {
                    if (j == col) continue;
                    const auto & M = ma[i][j];
                    return M;
                }
            }
            return 0.0;
        }
        //��������ʽ
        inline double Cofactors(int row, int col) const{
            const auto& k = std::pow(-1, row + col);
            return Minors(row, col)* k;
        }
        inline double Determinant() const
        {
            return ma[0][0] * ma[1][1] - ma[0][1] * ma[1][0];
        }
        Matrix2d Inverse() const
        {
            const auto& d = Determinant();
            if (std::fabs(d) < EPS) {
                std::abort();
            } else {
                Matrix2d mat;
                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        mat.ma[i][j] = Cofactors(j, i) / d;
                    }
                }
                return mat;
            }
        }
        Matrix2d Transpose() const
        {
            Matrix2d mat;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    mat.ma[i][j] = ma[j][i];
                }
            }
            return mat;
        }
        Matrix2d operator*(const double k) const
        {
            Matrix2d mat;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    mat.ma[i][j] = ma[i][j] * k;
                }
            }
            return mat;
        }
        inline Point2d operator*(const Point2d& p) const
        {
            std::vector<double> datas;
            datas.reserve(rows);
            for (int i = 0; i < rows; ++i) {
                datas.emplace_back(Row(i).Dot(p));
            }
            return Point2d(datas);
        }
        Matrix2d operator*(const Matrix2d& m) const
        {
            Matrix2d mat;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    mat.ma[i][j] = ma[i][0] * m.ma[0][j] + ma[i][1] * m.ma[1][j];
                }
            }
            return mat;
        }
        static Matrix2d GetMatrix(double angle)
        {
            const auto& costh = std::cos(angle);
            const auto& sinth = std::sin(angle);
            return Matrix2d(costh, -sinth, sinth, costh);
        }
        static Matrix2d GetMatrix(Point2d& a, Point2d& b)
        {
            const auto& theta = std::cos(a.Dot(b) / (a.Norm() * b.Norm()));
            double flag = a.Cross(b) > 0 ? 1.0 : -1.0;
            return GetMatrix(flag * theta);
        }
    };
    class Matrix3d {
    public:
        int rows = 3;
        int cols = 3;
        double ma[3][3] = { 0.0 };
        Matrix3d() 
        {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    ma[i][j] = 0.0;
                }
            }
        }
        Matrix3d(const double& val)
        {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    ma[i][j] = val;
                }
            }
        }
        Matrix3d(std::vector<double>& arr)
        {
            size_t size = arr.size();
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    ma[i][j] = (rows * i + j < size) ? arr[rows * i + j] : 0;
                }
            }
        }
        Matrix3d(double arr[], int size = 9)
        {
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    ma[i][j] = arr[rows * i + j];
                }
            }
        }
        Matrix3d(double arr[][3], int row_nums = 3)
        {
            for (int i = 0; i < row_nums; ++i) {
                for (int j = 0; j < cols; ++j) {
                    ma[i][j] = arr[i][j];
                }
            }
        }
        Matrix3d(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8)
        {
            ma[0][0] = a0;
            ma[0][1] = a1;
            ma[0][2] = a2;
            ma[1][0] = a3;
            ma[1][1] = a4;
            ma[1][2] = a5;
            ma[2][0] = a6;
            ma[2][1] = a7;
            ma[2][2] = a8;
        }
        bool IsIdentify() const
        {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    if (i == j) {
                        if (std::fabs(ma[i][j] - 1) > EPS) {
                            return false;
                        }
                    } else {
                        if (std::fabs(ma[i][j]) > EPS) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
        static Matrix3d Identify()
        {
            return Matrix3d(1, 0, 0, 0, 1, 0, 0, 0, 1);
        }

        static Matrix3d GetMatrix(const Point& axis, const double& angle)
        {
            Point u = axis;
            u.Normalize();
            Matrix3d m;
            m.ma[0][0] = cos(angle) + u.x * u.x * (1 - cos(angle));
            m.ma[0][1] = u.x * u.y * (1 - cos(angle)) - u.z * sin(angle);
            m.ma[0][2] = u.y * sin(angle) + u.x * u.z * (1 - cos(angle));
            //m.ma[0][3] = 0.0;

            m.ma[1][0] = u.z * sin(angle) + u.x * u.y * (1 - cos(angle));
            m.ma[1][1] = cos(angle) + u.y * u.y * (1 - cos(angle));
            m.ma[1][2] = -u.x * sin(angle) + u.y * u.z * (1 - cos(angle));
            //ma[1][3] = 0.0;

            m.ma[2][0] = -u.y * sin(angle) + u.x * u.z * (1 - cos(angle));
            m.ma[2][1] = u.x * sin(angle) + u.y * u.z * (1 - cos(angle));
            m.ma[2][2] = cos(angle) + u.z * u.z * (1 - cos(angle));
            //m.ma[2][3] = 0.0;

            //m.ma[3][0] = 0.0;
            //m.ma[3][1] = 0.0;
            //m.ma[3][2] = 0.0;
            //m.ma[3][3] = 1.0;
            return m;
        }
        static Matrix3d GetMatrix(const Point& a, const Point& b)
        {
            double rotationAngle = std::acos(a.Dot(b) / (a.Norm() * b.Norm()));
            Point rotationAxis = a.Cross(b);
            return GetMatrix(rotationAxis, rotationAngle);
        }
        void Copy(const Matrix3d& m)
        {
            rows = m.rows, cols = m.cols;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    ma[i][j] = m.ma[i][j];
                }
            }
        }

        inline Point Row(int row)const
        {
            return Point(ma[row][0], ma[row][1], ma[row][2]);
        }
        inline Point Col(int col)const
        {
            return Point(ma[0][col], ma[1][col], ma[2][col]);
        }
        //����ʽ,�˴�û��������ʽ
        inline Matrix2d Minors(int row, int col)const
        {
            if (row < 0 || row >= rows) {
                std::abort();
            }
            if (col < 0 || col >= cols) {
                std::abort();
            }
            std::vector<double> datas;
            datas.reserve(4);
            for (int i = 0; i < rows; ++i) {
                if (i == row) continue;
                for (int j = 0; j < cols; ++j) {
                    if (j == col) continue;
                    const auto & M = ma[i][j];
                    datas.emplace_back(M);
                }
            }
            return Matrix2d(datas);
        }
        //��������ʽ
        inline double Cofactors(int row, int col) const
        {
            const auto& k = std::pow(-1, row + col);
            return Minors(row,col).Determinant()* k;
        }
        inline double Determinant() const{
            double ans = 0.0;
            for (int j = 0; j < cols; ++j) {
                ans += ma[0][j] * Cofactors(0, j);
            }
            return ans;
        }
        //��������ʽ����������
        //��˹��Ԫ����δʵ��
        Matrix3d Inverse()
        {
            const auto& d = Determinant();
            if (std::fabs(d) < EPS) {
                std::abort();
            } else {
                Matrix3d mat;
                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        mat.ma[i][j] = Cofactors(j, i) / d;
                    }
                }
                return mat;
            }
        }
        Matrix3d operator*(const double k) const
        {
            Matrix3d mat;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    mat.ma[i][j] = ma[i][j] * k;
                }
            }
            return mat;
        }
        inline Point operator*(const Point& p) const
        {
            std::vector<double> datas;
            datas.reserve(rows);
            for (int i = 0; i < rows; ++i) {
                datas.emplace_back(Row(i).Dot(p));
            }
            return Point(datas);
        }
        inline void operator*=(Point& p) const
        {
            const auto& x = Row(0).Dot(p);
            const auto& y = Row(1).Dot(p);
            const auto& z = Row(2).Dot(p);
            p.x = x;
            p.y = y;
            p.z = z;
        }
        Matrix3d operator*(const Matrix3d& m) const
        {
            Matrix3d mat;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    mat.ma[i][j] = ma[i][0] * m.ma[0][j] + ma[i][1] * m.ma[1][j] + ma[i][2] * m.ma[2][j];
                }
            }
            return mat;
        }
        friend std::istream& operator>>(std::istream& is, Matrix3d& m)
        {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    is >> m.ma[i][j];
                }
            }
            return is;
        }
        friend std::ostream& operator<<(std::ostream& os, const Matrix3d& m)
        {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    std::cout << m.ma[i][j] << " ";
                }
                std::cout << std::endl;
            }
            return os;
        }
        
        ~Matrix3d() {}
    };
    class Matrix4d {
    public:

    };
}