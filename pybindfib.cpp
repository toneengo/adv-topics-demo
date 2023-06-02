#include <pybind11/pybind11.h>
#include "cfib.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pybindfib, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("pfib", &fib, "fibonacci");
}
