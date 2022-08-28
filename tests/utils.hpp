#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <catch2/catch_all.hpp>

template<typename T> struct remove_complex
{
    using type = T;
};
template<typename T> struct remove_complex<std::complex<T>>
{
    using type = T;
};
template<typename T> using remove_complex_t = typename remove_complex<T>::type;

template<typename T> struct is_complex : std::false_type
{ };

template<typename T> struct is_complex<std::complex<T>> : std::true_type
{ };

template<typename T> constexpr bool is_complex_v = is_complex<T>::value;

template<typename T> concept std_complex = is_complex_v<T>;

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

template<class T> struct is_std_vector
{
    constexpr static bool value = false;
};
template<class T, class Alloc> struct is_std_vector<std::vector<T, Alloc>>
{
    constexpr static bool value = true;
};
template<class T> concept std_vector = is_std_vector<T>::value;

template<typename T, typename AllocComp> struct VectorApprox : Catch::Matchers::MatcherGenericBase
{
private:
    const std::vector<T, AllocComp>& comp_;
    mutable Catch::Approx approx = Catch::Approx::custom();

public:
    explicit VectorApprox(const std::vector<T, AllocComp>& comp) : comp_{comp} { }

    std::string describe() const override
    {
        std::ostringstream ss;
        ss << "is approx to [";
        for(const auto& elt : comp_)
        {
            ss << elt << ", ";
        }
        ss << "]";
        return ss.str();
    }

    template<typename AllocMatch> bool match(const std::vector<T, AllocMatch>& v) const
    {
        if(comp_.size() != v.size())
        {
            return false;
        }

        if constexpr(is_complex_v<T>)
        {
            for(size_t i = 0; i < v.size(); i++)
            {
                if((std::real(comp_[i]) != approx(std::real(v[i])))
                   || (std::imag(comp_[i]) != approx(std::imag(v[i]))))
                {
                    return false;
                }
            }
        }
        else
        {
            for(size_t i = 0; i < v.size(); i++)
            {
                if(comp_[i] != approx(v[i]))
                {
                    return false;
                }
            }
        }
        return true;
    }

    VectorApprox& epsilon(remove_complex_t<T> new_eps)
    {
        approx.epsilon(new_eps);
        return *this;
    }

    VectorApprox& margin(remove_complex_t<T> new_margin)
    {
        approx.margin(new_margin);
        return *this;
    }
};

template<std_vector T> auto Approx(const T& comp) -> decltype(auto)
{
    return VectorApprox(comp);
}

template<typename UINT, template<typename> class Basis>
Eigen::MatrixXd basisMatrix(const Basis<UINT>& basis)
{
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(1U << basis.getN(), basis.getDim());
    for(unsigned int n = 0; n < basis.getDim(); ++n)
    {
        auto bvec = basis.basisVec(n);
        for(const auto& p : bvec)
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
    using Catch::Matchers::WithinAbs;
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
