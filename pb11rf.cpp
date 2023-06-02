#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "multiDomainRootFinder.cpp"

namespace py = pybind11;
using ddfunc = std::function<double(double)>;

double pb_linear_root_finder(
    double start, double end,
    ddfunc const &f, ddfunc const &fprime,
    int verificationLevels=2, int subdivisionLevels=1,
    double dfsTerminationDistance= 1 / pow(2, DFS_TERMINATION_LIMIT),
    int maxNewtonIterations=250
) {
    int totalCalls = 0;
    return linearRootFinder(
        start, end,
        f, fprime,
        totalCalls,
        verificationLevels, subdivisionLevels,
        dfsTerminationDistance,
        maxNewtonIterations
    );
}

std::tuple<double, double> pb_newton(
    ddfunc const &f,
    double x0,
    ddfunc const &fprime,
    int useless = 0,
    double tol = 1e-8,
    int maxiter = 100,
    double rtol = 0.0,
    bool disp = false
) {
    int fcalls = 0;
    return newton(f, fprime, x0, fcalls, tol, maxiter);
}


std::pair<std::vector<int>, std::vector<double>> pb_dfs_multi_domain(
    double start, double end,
    std::vector<ddfunc> fs, std::vector<ddfunc> fprimes,
    int verificationLevels = 2, int subdivisionLevels = 1,
    double dfsTerminationDistance = 1 / pow(2, DFS_TERMINATION_LIMIT)
) {
    //std::vector<ddfunc> fs = bernsteinGenerator(4);
    //std::vector<ddfunc> fprimes = bernsteinPrimeGenerator(4);
    return dfsMultiDomain(
        start, end,
        fs, fprimes,
        verificationLevels, subdivisionLevels,
        dfsTerminationDistance,
        250
    );
}

PYBIND11_MODULE(pb11rf, m) {
    m.doc() = "pybind11 plugin"; // optional module docstring

    m.def("pb_linear_root_finder", &pb_linear_root_finder, "linear root finder",
            py::arg("start"), py::arg("end"),
            py::arg("f"), py::arg("fprime"),
            py::arg("verificationLevels") = 2, py::arg("subdivisionLevels") = 1,
            py::arg("dfsTerminationDistance") = 1 / pow(2, DFS_TERMINATION_LIMIT),
            py::arg("maxNewtonIterations") = 250
    );
    m.def("pb_newton", &pb_newton, "newton",
            py::arg("f"), py::arg("x0"), py::arg("fprime"),
            py::arg("useless") = 0, py::arg("tol") = 1.48e-8,
            py::arg("maxiter") = 100, py::arg("rtol") = 0.0,
            py::arg("disp") = false
    );
    m.def("pb_dfs_multi_domain", &pb_dfs_multi_domain, "dfs multi domain",
            py::arg("start"), py::arg("end"),
            py::arg("fs"), py::arg("fprimes"),
            py::arg("verificationLevels") = 2, py::arg("subdivisionLevels") = 1,
            py::arg("dfsTerminationDistance") = 1 / pow(2, DFS_TERMINATION_LIMIT)
    );
}
