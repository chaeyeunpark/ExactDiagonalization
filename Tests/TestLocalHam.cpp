#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 


#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include <EDP/LocalHamiltonian.hpp>
#include <EDP/ConstructSparseMat.hpp>


Eigen::SparseMatrix<double> getSX()
{
	Eigen::SparseMatrix<double> res(2,2);
	res.insert(0,1) = 1.0;
	res.insert(1,0) = 1.0;
	res.makeCompressed();
	return res;
}
Eigen::SparseMatirx<double> getSZ()
{
	Eigen::SparseMatrix<double> res(2,2);
	res.insert(0,0) = 1.0;
	res.insert(1,1) = -1.0;
	res.makeCompressed();
	return res;
}

Eigen::MatrixXd (int N, const Eigen::Vector3d& v)
{
}

Eigen::SparseMatrix<double> getSXXYY()
{
	Eigen::SparseMatrix<double> res(4,4);
	res.insert(1,2) = 2.0;
	res.insert(2,1) = 2.0;
	res.makeCompressed();
	return res;
}
Eigen::SparseMatrix<double> getSXX()
{
	Eigen::SparseMatrix<double> res(4,4);
	res.insert(0,3) = 1.0;
	res.insert(1,2) = 1.0;
	res.insert(2,1) = 1.0;
	res.insert(3,0) = 1.0;
	res.makeCompressed();
	return res;
}
Eigen::SparseMatrix<double> getSYY()
{
	Eigen::SparseMatrix<double> res(4,4);
	res.insert(0,3) = -1.0;
	res.insert(1,2) = 1.0;
	res.insert(2,1) = 1.0;
	res.insert(3,0) = -1.0;
	res.makeCompressed();
	return res;
}

Eigen::SparseMatrix<double> getSZZ()
{
	Eigen::SparseMatrix<double> res(4,4);
	res.insert(0,0) = 1.0;
	res.insert(1,1) = -1.0;
	res.insert(2,2) = -1.0;
	res.insert(3,3) = 1.0;
	res.makeCompressed();
	return res;
}

