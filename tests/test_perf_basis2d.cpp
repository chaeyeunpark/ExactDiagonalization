#include <chrono>
#include <edlib/Basis/Basis2D.hpp>
#include <random>

int main()
{
    using namespace edlib;

    std::random_device rd;
    std::default_random_engine re{rd()};

    const uint32_t Lx = 5;
    const uint32_t Ly = 4;
    const uint32_t N = (Lx * Ly);

    auto start = std::chrono::high_resolution_clock::now();

    Basis2D<uint32_t> basis(Lx, Ly, 0, 0, false);
    auto end = std::chrono::high_resolution_clock::now();

    auto elapsed_miliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time to construct the basis: " << elapsed_miliseconds.count() << std::endl;

    std::uniform_int_distribution<uint32_t> basisDist(0, basis.getDim() - 1);
    std::uniform_int_distribution<uint32_t> sigmaDist(0, (1u << N) - 1);

    start = std::chrono::high_resolution_clock::now();
    const auto n_iter = 1'000'000u;
    for(uint32_t n = 0; n < n_iter; ++n)
    {
        volatile auto r = basis.hamiltonianCoeff(sigmaDist(rd), basisDist(rd));
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed_miliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time to call " << n_iter << " hamiltonianCoeff : " << elapsed_miliseconds.count()
              << std::endl;

    return 0;
}
