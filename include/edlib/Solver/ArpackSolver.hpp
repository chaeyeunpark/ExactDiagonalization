#pragma once
#include "arpackdef.h"

#include <Eigen/Dense>
#include <random>
#include <vector>

extern "C"
{
    void dsaupd_(a_int*, char*, a_int*, char*, a_int*, double*, double*, a_int*, double*, a_int*,
                 a_int*, a_int*, double*, double*, a_int*, a_int*);

    void dseupd_(a_int* rvec, char* howmny, a_int* select, double* d, double* z, a_int* ldz,
                 double* sigma, char* bmat, a_int* n, char* which, a_int* nev, double* tol,
                 double* resid, a_int* ncv, double* v, a_int* ldv, a_int* iparam, a_int* inptr,
                 double* workd, double* workl, a_int* lworkl, a_int* info);
}

namespace edlib
{
enum class ErrorType
{
    NormalExit,
    NotConverged,
    IncorrectParams,
    Others
};

template<typename MatrixVectorOp> class ArpackSolver
{
private:
    MatrixVectorOp& op_;
    const a_int dim_;
    std::vector<double> d_;
    std::vector<double> z_;

public:
    ArpackSolver(MatrixVectorOp& op) : op_{op}, dim_{static_cast<a_int>(op.rows())}
    {
        assert(op.rows() == op.cols());
    }

    ErrorType solve(a_int nev, uint32_t max_iter = 1024, double tol = 1e-10)
    {
        a_int dim = dim_;
        if(dim <= 0 || nev <= 0)
        {
            return ErrorType::IncorrectParams;
        }

        a_int ido = 0;
        char bmat[] = "I";
        a_int n = dim;
        char which[] = "SA";
        std::vector<double> resid(dim);

        a_int ncv = 3 * nev;
        std::vector<double> v(ncv * dim);

        a_int ldv = dim;
        a_int iparam[11] = {
            1, // ishift
            0, // levec (not used)
            static_cast<a_int>(max_iter), // maxiter
            1, // nb
            0, // nconv
            0, // iupd (not used)
            1, // mode (1 is usual eigenvalue problem)
            0, // np
            0,
            0,
            0 // only for output
        };

        a_int ipntr[14];
        std::vector<double> workd(3 * dim, 0.);

        int lworkl = 3 * ncv * ncv + 6 * ncv;
        std::vector<double> workl(lworkl, 0);

        a_int info = 0;

        // first call
        dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid.data(), &ncv, v.data(), &ldv, iparam,
                ipntr, workd.data(), workl.data(), &lworkl, &info);

        while(ido == -1 || ido == 1)
        {
            op_.perform_op(workd.data() + ipntr[0] - 1, workd.data() + ipntr[1] - 1);

            dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid.data(), &ncv, v.data(), &ldv, iparam,
                    ipntr, workd.data(), workl.data(), &lworkl, &info);
        }

        if(info == 1 || iparam[4] != nev)
        {
            return ErrorType::NotConverged;
        }
        else if(info != 0)
        {
            return ErrorType::IncorrectParams;
        }

        a_int rvec = 1;

        d_.resize(nev + 1);
        z_.resize((dim + 1) * (nev + 1));
        a_int ldz = dim + 1;
        double sigma = 0.0;

        std::vector<a_int> select(ncv);
        std::fill(select.begin(), select.end(), 1);

        char howmny[] = "All";
        dseupd_(&rvec, howmny, select.data(), d_.data(), z_.data(), &ldz, &sigma, bmat, &n, which,
                &nev, &tol, resid.data(), &ncv, v.data(), &ldv, iparam, ipntr, workd.data(),
                workl.data(), &lworkl, &info);

        return ErrorType::NormalExit;
    }

    Eigen::Map<const Eigen::VectorXd> eigenvalues() const
    {
        return Eigen::Map<const Eigen::VectorXd>(d_.data(), d_.size() - 1);
    }

    Eigen::Map<const Eigen::MatrixXd, 0, Eigen::OuterStride<Eigen::Dynamic>> eigenvectors() const
    {
        return {z_.data(), dim_, static_cast<Eigen::Index>(d_.size() - 1), {dim_ + 1}};
    }
};
} // namespace edlib
