#include <Basis/TIBasisZ2.hpp>
#include <iostream>

int main()
{
    using UINT = uint32_t;
    TIBasisZ2<UINT> basis(10, 0, true, -1);

    std::cout << "dim: " << basis.getDim() << std::endl;

    uint32_t idx = 8;
    std::cout << basis.getNthRep(idx) << std::endl;

    auto p_vec = basis.basisVec(idx);
    for(const auto& p : p_vec)
    {
        std::cout << p.first << "\t" << p.second << std::endl;
    }
}
