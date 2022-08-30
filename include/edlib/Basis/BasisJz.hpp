#pragma once
#include <iterator>

namespace edlib
{
//! \ingroup Basis
//! Basis for U(1) symmetric subspace.
template<typename UINT> class BasisJz
{
private:
    uint32_t N_;
    uint32_t nup_;

public:
    struct BasisJzIterator : public std::iterator<std::forward_iterator_tag, UINT>
    {
        /*
        using iterator_category = std::forward_iterator_tag;
        using value_type = UINT;
        using different_type = void;
        using pointer = void;
        using reference = void;
        */

        UINT n_;

        BasisJzIterator() : n_{static_cast<UINT>(0)} { }

        explicit BasisJzIterator(UINT val) : n_(val) { }

        BasisJzIterator& operator++() // prefix
        {
            next();
            return *this;
        }
        BasisJzIterator operator++(int) // postfix
        {
            BasisJzIterator r(*this);
            next();
            return r;
        }
        void
        next() // https://stackoverflow.com/questions/8281951/bit-hack-to-generate-all-integers-with-a-given-number-of-1s
        {
            UINT t = n_ | (n_ - 1);
            // NOLINTNEXTLINE(hicpp-signed-bitwise)
            UINT w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(n_) + 1));
            n_ = w;
        }
        UINT operator*() const { return n_; }
        bool operator==(const BasisJzIterator& rhs) { return n_ == rhs.n_; }
        bool operator!=(const BasisJzIterator& rhs) { return n_ != rhs.n_; }
    };

    //! Construct a basis for the subspace. The dimension is \f$N \choose nup\f$.
    //! \param N number of total spins
    //! \param nup number of spin ups(\f$|\uparrow \rangle\f$)
    explicit BasisJz(uint32_t N, uint32_t nup) : N_(N), nup_(nup) { }

    BasisJzIterator begin() { return BasisJzIterator{(UINT(1) << UINT(nup_)) - 1}; }

    BasisJzIterator end()
    {
        BasisJzIterator r{((UINT(1) << UINT(nup_)) - 1) << (UINT(N_) - UINT(nup_))};
        r.next();
        return r;
    }

    [[nodiscard]] UINT size() const
    {
        auto res = static_cast<UINT>(1U);
        UINT k = nup_;
        // Since C(n, k) = C(n, n-k)
        if(k > N_ - k)
        {
            k = N_ - k;
        }

        // Calculate value of
        // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
        for(UINT i = 0U; i < k; ++i)
        {
            res *= (N_ - i);
            res /= (i + 1);
        }
        return res;
    }
};
} // namespace edlib
