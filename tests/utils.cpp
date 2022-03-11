#include "utils.hpp"

Eigen::SparseMatrix<double> getSXXYY()
{
    Eigen::SparseMatrix<double> res(4, 4);
    res.insert(1, 2) = 2.0;
    res.insert(2, 1) = 2.0;
    res.makeCompressed();
    return res;
}

Eigen::SparseMatrix<double> getSXX()
{
    Eigen::SparseMatrix<double> res(4, 4);
    res.insert(0, 3) = 1.0;
    res.insert(1, 2) = 1.0;
    res.insert(2, 1) = 1.0;
    res.insert(3, 0) = 1.0;
    res.makeCompressed();
    return res;
}

Eigen::SparseMatrix<double> getSYY()
{
    Eigen::SparseMatrix<double> res(4, 4);
    res.insert(0, 3) = -1.0;
    res.insert(1, 2) = 1.0;
    res.insert(2, 1) = 1.0;
    res.insert(3, 0) = -1.0;
    res.makeCompressed();
    return res;
}

Eigen::SparseMatrix<double> getSZZ()
{
    Eigen::SparseMatrix<double> res(4, 4);
    res.insert(0, 0) = 1.0;
    res.insert(1, 1) = -1.0;
    res.insert(2, 2) = -1.0;
    res.insert(3, 3) = 1.0;
    res.makeCompressed();
    return res;
}

Eigen::SparseMatrix<double> getSX()
{
    Eigen::SparseMatrix<double> res(2, 2);
    res.insert(0, 1) = 1.0;
    res.insert(1, 0) = 1.0;
    res.makeCompressed();
    return res;
}

Eigen::SparseMatrix<std::complex<double>> getSY()
{
    Eigen::SparseMatrix<std::complex<double>> res(2, 2);
    constexpr std::complex<double> I(0., 1.);
    res.insert(0, 1) = -I;
    res.insert(1, 0) = I;
    res.makeCompressed();
    return res;
}

Eigen::SparseMatrix<double> getSZ()
{
    Eigen::SparseMatrix<double> res(2, 2);
    res.insert(0, 0) = 1.0;
    res.insert(1, 1) = -1.0;
    res.makeCompressed();
    return res;
}
