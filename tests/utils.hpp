#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <catch2/catch.hpp>

inline void make_first_positive(Eigen::VectorXi& v)
{
    for(uint32_t n = 0; n < v.size(); ++n)
    {
        if(v(n) == 0)
        {
            continue;
        }
        if(v(n) < 0)
        {
            v *= -1;
            return;
        }
        else
        {
            return;
        }
    }
}

struct hash_vector
{
    size_t operator()(const Eigen::VectorXi& v) const
    {
        std::size_t seed = v.size();
        for(uint32_t n = 0; n < v.size(); ++n)
        {
            seed ^= abs(v(n)) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
        }
        return seed;
    }
};

inline int powi(int base, unsigned int exp)
{
    int res = 1;
    while(exp != 0U)
    {
        if((exp & 1U) != 0U)
        {
            res *= base;
        }
        exp >>= 1U;
        base *= base;
    }
    return res;
}

template<typename UINT, template<typename> class Basis>
Eigen::MatrixXd basisMatrix(const Basis<UINT>& basis)
{
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(1U << basis.getN(), basis.getDim());
    for(unsigned int n = 0; n < basis.getDim(); ++n)
    {
        auto bvec = basis.basisVec(n);
        for(const auto p : bvec)
        {
            res(p.first, n) = p.second;
        }
    }
    return res;
}

template<typename Basis> Eigen::VectorXd flip(Basis&& basis, const Eigen::VectorXd& r)
{
    Eigen::VectorXd res(r.size());
    for(int i = 0; i < r.size(); i++)
    {
        res(basis.flip(i)) = r(i);
    }
    return res;
}

inline void TestBasisMatrix(const Eigen::MatrixXd& r)
{
    using Catch::WithinAbs;
    Eigen::MatrixXd id = Eigen::MatrixXd::Identity(r.cols(), r.cols());
    REQUIRE_THAT((r.transpose() * r - id).cwiseAbs().maxCoeff(), WithinAbs(0.0, 1e-8));
}

Eigen::SparseMatrix<double> getSX();
Eigen::SparseMatrix<std::complex<double>> getSY();
Eigen::SparseMatrix<double> getSZ();
Eigen::SparseMatrix<double> getSXXYY();
Eigen::SparseMatrix<double> getSXX();
Eigen::SparseMatrix<double> getSYY();
Eigen::SparseMatrix<double> getSZZ();
