#pragma once

#if defined(INTEL_MKL_VERSION)
#include "NodeMV_MKL.hpp"
#else
#include "NodeMV_Defualt.hpp"
#endif
