// BioFVM external-solver cross-check driver for the Rust reaction_diffusion
// module (#408). Steady-state diffusion-decay D grad^2 c - k c = 0 (lambda =
// sqrt(D/k)) with Dirichlet point sources (c=1, re-clamped each step) and no-flux
// boundaries. Canonical BioFVM usage. Serial build (omp pragmas ignored); see the
// runner script for the clang/SDK/omp-stub incantation.
//   usage: biofvm_rd_check N h lambda dt nsteps  r0 c0 l0 [r1 c1 l1 ...]
#include "BioFVM.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace BioFVM;
int main(int argc, char** argv) {
    int N = atoi(argv[1]); double h = atof(argv[2]); double lambda = atof(argv[3]);
    double dt = atof(argv[4]); int nsteps = atoi(argv[5]);
    double D = lambda * lambda, k = 1.0;            // sqrt(D/k) = lambda
    Microenvironment M;
    M.set_density(0, "supply", "dimensionless");
    M.resize_space_uniform(0, N * h, 0, N * h, 0, N * h, h);
    M.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_3D;
    M.diffusion_coefficients[0] = D;
    M.decay_rates[0] = k;
    for (int n = 0; n < M.number_of_voxels(); n++) M.density_vector(n)[0] = 0.0;
    std::vector<int> src;
    for (int a = 6; a + 2 < argc; a += 3) {
        std::vector<double> p = {(atof(argv[a]) + 0.5) * h, (atof(argv[a+1]) + 0.5) * h, (atof(argv[a+2]) + 0.5) * h};
        src.push_back(M.nearest_voxel_index(p));
    }
    for (int s : src) M.density_vector(s)[0] = 1.0;
    for (int step = 0; step < nsteps; step++) {
        M.simulate_diffusion_decay(dt);
        for (int s : src) M.density_vector(s)[0] = 1.0;   // hard Dirichlet
    }
    FILE* f = fopen("biofvm_field.csv", "w"); fprintf(f, "x,y,z,c\n");
    for (int n = 0; n < M.number_of_voxels(); n++) {
        std::vector<double>& ctr = M.mesh.voxels[n].center;
        fprintf(f, "%.4f,%.4f,%.4f,%.8f\n", ctr[0], ctr[1], ctr[2], M.density_vector(n)[0]);
    }
    fclose(f);
    fprintf(stderr, "biofvm done N=%d voxels=%d sources=%zu dt=%g steps=%d\n", N, M.number_of_voxels(), src.size(), dt, nsteps);
    return 0;
}
