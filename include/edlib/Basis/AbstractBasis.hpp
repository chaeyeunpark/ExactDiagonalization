#pragma once
#include <climits>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <tuple>
#include <vector>

namespace edlib
{
template<typename UINT> class AbstractBasis
{
private:
    const uint32_t N_;
    const UINT ups_;

public:
    explicit AbstractBasis(uint32_t N)
        : N_{N}, ups_{(~UINT(0)) >> (sizeof(UINT) * CHAR_BIT - N)} { }

    [[nodiscard]] inline uint32_t getN() const { return N_; }

    [[nodiscard]] inline UINT getUps() const { return ups_; }

    [[nodiscard]] UINT mask(std::initializer_list<uint32_t> pos) const
    {
        UINT s = 0;
        for(uint32_t p : pos)
        {
            s ^= (UINT(1) << p);
        }
        return s;
    }

    [[nodiscard]] virtual std::size_t getDim() const = 0;
    [[nodiscard]] virtual UINT getNthRep(uint32_t n) const = 0;
    [[nodiscard]] virtual std::pair<int, double> hamiltonianCoeff(UINT bsigma, int adix) const = 0;

    [[nodiscard]] virtual std::vector<std::pair<UINT, double>> basisVec(uint32_t n) const = 0;

    AbstractBasis(const AbstractBasis&) = delete;
    AbstractBasis(AbstractBasis&&) noexcept = default;

    AbstractBasis& operator=(const AbstractBasis&) = delete;
    AbstractBasis& operator=(AbstractBasis&&) noexcept = default;

    virtual ~AbstractBasis() = default;
};
} // namespace edlib
