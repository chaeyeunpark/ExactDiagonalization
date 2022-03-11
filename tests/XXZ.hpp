#pragma once
#include <Eigen/Eigen>

class XXZ
{
private:
    uint32_t n_;
    double J_;
    double Delta_;

public:
    XXZ(uint32_t n, double J, double Delta) : n_(n), J_(J), Delta_(Delta) { }

    template<class State> typename State::T operator()(const State& smp) const
    {
        typename State::T s = 0.0;
        // Nearest-neighbor
        for(uint32_t i = 0; i < n_; i++)
        {
            int zz = smp.sigmaAt(i) * smp.sigmaAt((i + 1) % n_);
            s += J_ * Delta_ * zz; // zz
            s += J_ * (1 - zz) * smp.ratio(i, (i + 1) % n_); // xx+yy
        }
        return s;
    }

    std::map<uint32_t, double> operator()(uint32_t col) const;
};
