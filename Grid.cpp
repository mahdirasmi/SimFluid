#include "Grid.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdexcept>

using namespace std;

// ============================================================
//  CONSTRUCTOR 1 — UNIFORM
//  Equal cell size everywhere.
// ============================================================
Grid::Grid(int Nx_total, int Ny_total, double Ld, double Hd)
    : meshType(MeshType::UNIFORM),
      Nx{Nx_total}, Ny{Ny_total},
      L(Ld), H(Hd),
      First_cell_size(0.0), r_growth(1.0),
      dxx(Nx_total, 0.0), dyy(Ny_total, 0.0),
      x_c(Nx_total + 1, 0.0), y_c(Ny_total + 1, 0.0),
      dl{}, dH{}
{
    if (Nx_total <= 0 || Ny_total <= 0)
        throw invalid_argument("[Grid] Nx and Ny must be > 0.");
    if (L <= 0.0 || H <= 0.0)
        throw invalid_argument("[Grid] Domain dimensions L and H must be > 0.");

    compute_uniform_dxx();
    compute_uniform_dyy();
    build_coordinates();
}

// ============================================================
//  CONSTRUCTOR 2 — NON-UNIFORM / SQUARE CYLINDER
//  3-block stretching in x and y, refined around (Xc,Yc).
//  Middle block node count is auto-computed from FC.
// ============================================================
Grid::Grid(int Nx1, int Nx3,
           int Ny1, int Ny3,
           double Ld, double Hd, double FC,
           double Xc, double Yc, double D)
    : meshType(MeshType::NONUNIFORM_CYL),
      Nx{Nx1, 0, Nx3}, Ny{Ny1, 0, Ny3},
      L(Ld), H(Hd),
      First_cell_size(FC), r_growth(1.0),
      dl{ Xc - D/2.0, D, L - Xc - D/2.0 },
      dH{ Yc - D/2.0, D, H - Yc - D/2.0 },
      dxx{}, dyy{},
      x_c{}, y_c{}
{
    Middle_block_size();

    int total_nx = Nx[0] + Nx[1] + Nx[2];
    int total_ny = Ny[0] + Ny[1] + Ny[2];

    dxx.assign(total_nx, 0.0);
    dyy.assign(total_ny, 0.0);
    x_c.assign(total_nx + 1, 0.0);
    y_c.assign(total_ny + 1, 0.0);

    compute_dxx_cyl();
    compute_dyy_cyl();
    build_coordinates();
}

// ============================================================
//  CONSTRUCTOR 3 — NON-UNIFORM / SYMMETRIC
//  Cells grow from each wall toward the centre.
//  growth_rate > 1  (typical: 1.05 – 1.15)
// ============================================================
Grid::Grid(int Nx_total, int Ny_total,
           double Ld, double Hd,
           double growth_rate)
    : meshType(MeshType::NONUNIFORM_Symmetry),
      Nx{Nx_total}, Ny{Ny_total},
      L(Ld), H(Hd),
      First_cell_size(0.0), r_growth(growth_rate),
      dxx(Nx_total, 0.0), dyy(Ny_total, 0.0),
      x_c(Nx_total + 1, 0.0), y_c(Ny_total + 1, 0.0),
      dl{}, dH{}
{
    if (Nx_total <= 0 || Ny_total <= 0)
        throw invalid_argument("[Grid] Nx and Ny must be > 0.");
    if (L <= 0.0 || H <= 0.0)
        throw invalid_argument("[Grid] Domain dimensions L and H must be > 0.");
    if (growth_rate <= 1.0)
        throw invalid_argument("[Grid] growth_rate must be > 1.0 for wall refinement.");

    compute_dxx_symmetry();
    compute_dyy_symmetry();
    build_coordinates();
}

// ============================================================
//  Middle block size (cylinder constructor only)
// ============================================================
void Grid::Middle_block_size()
{
    Nx[1] = static_cast<int>(round(dl[1] / First_cell_size));
    Ny[1] = static_cast<int>(round(dH[1] / First_cell_size));
}

// ============================================================
//  Bisection solver for geometric-series growth rate r
// ============================================================
double Grid::solve_r_bisection(double Lseg, double dx1, int N)
{
    double lo = 1.0 + 1e-10;
    double hi = pow(Lseg / dx1, 1.0 / (N - 1));

    auto f = [&](double r) {
        return dx1 * (pow(r, N) - 1.0) / (r - 1.0) - Lseg;
    };

    for (int k = 0; k < 200; ++k) {
        double mid = 0.5 * (lo + hi);
        (f(mid) > 0.0) ? hi = mid : lo = mid;
    }
    return 0.5 * (lo + hi);
}

// ============================================================
//  Stretching function
// ============================================================
double Grid::stretch(double Lseg, int N, int i, double dx1)
{
    double r = solve_r_bisection(Lseg, dx1, N);
    return dx1 * pow(r, i);
}

// ============================================================
//  CYLINDER: DXX
// ============================================================
void Grid::compute_dxx_cyl()
{
    int N0 = Nx[0], N1 = Nx[1], N2 = Nx[2];
    int Ntot = N0 + N1 + N2;
    vector<double> w(Ntot, 0.0);
    double dx1 = First_cell_size;

    // Block 0 (left): finest cell at right end, grows leftward
    for (int i = 0; i < N0; ++i)
        w[N0 - 1 - i] = stretch(dl[0], N0, i, dx1);

    // Block 1 (middle): uniform at the finest size
    for (int i = 0; i < N1; ++i)
        w[N0 + i] = w[N0 - 1];

    // Block 2 (right): finest cell at left end, grows rightward
    for (int i = 0; i < N2; ++i)
        w[N0 + N1 + i] = stretch(dl[2], N2, i, dx1);

    // Normalise each outer block to its exact physical length
    double s0 = 0.0, s2 = 0.0;
    for (int i = 0;     i < N0;   ++i) s0 += w[i];
    for (int i = N0+N1; i < Ntot; ++i) s2 += w[i];

    for (int i = 0; i < N0; ++i)
        dxx[i] = dl[0] * w[i] / s0;
    for (int i = 0; i < N1; ++i)
        dxx[N0 + i] = dl[1] * w[N0 + i];   // middle already uniform
    for (int i = 0; i < N2; ++i)
        dxx[N0 + N1 + i] = dl[2] * w[N0 + N1 + i] / s2;
}

// ============================================================
//  CYLINDER: DYY
// ============================================================
void Grid::compute_dyy_cyl()
{
    int N0 = Ny[0], N1 = Ny[1], N2 = Ny[2];
    int Ntot = N0 + N1 + N2;
    vector<double> wy(Ntot, 0.0);
    double dy1 = First_cell_size;

    for (int j = 0; j < N0; ++j)
        wy[N0 - 1 - j] = stretch(dH[0], N0, j, dy1);

    for (int j = 0; j < N1; ++j)
        wy[N0 + j] = wy[N0 - 1];

    for (int j = 0; j < N2; ++j)
        wy[N0 + N1 + j] = stretch(dH[2], N2, j, dy1);

    double s0 = 0.0, s2 = 0.0;
    for (int j = 0;     j < N0;   ++j) s0 += wy[j];
    for (int j = N0+N1; j < Ntot; ++j) s2 += wy[j];

    for (int j = 0; j < N0; ++j)
        dyy[j] = dH[0] * wy[j] / s0;
    for (int j = 0; j < N1; ++j)
        dyy[N0 + j] = dH[1] * wy[N0 + j];
    for (int j = 0; j < N2; ++j)
        dyy[N0 + N1 + j] = dH[2] * wy[N0 + N1 + j] / s2;

    cout << "[Grid] Cylinder mesh — middle block: Nx=" << Nx[1]
         << "  Ny=" << Ny[1] << "\n";
}

// ============================================================
//  UNIFORM: DXX / DYY
// ============================================================
void Grid::compute_uniform_dxx()
{
    double cell = L / Nx[0];
    dxx.assign(Nx[0], cell);
}

void Grid::compute_uniform_dyy()
{
    double cell = H / Ny[0];
    dyy.assign(Ny[0], cell);
}

// ============================================================
//  SYMMETRIC: DXX
//  Cells shrink from the centre toward each wall.
// ============================================================
void Grid::compute_dxx_symmetry()
{
    int N    = Nx[0];
    int half = N / 2;
    vector<double> w(N, 0.0);
    double r = r_growth;

    // Left half: finest cell at the wall (i=0), grows toward centre
    for (int i = 0; i < half; ++i)
        w[i] = pow(r, i);

    // Right half: mirror
    for (int i = 0; i < half; ++i)
        w[half + half-1 - i] = pow(r, i);

    // Centre cell if odd total
    if (N % 2 != 0)
        w[half] = w[half - 1];

    // Normalise so sum == L
    double sum = 0.0;
    for (int i = 0; i < N; ++i) sum += w[i];
    for (int i = 0; i < N; ++i) dxx[i] = L * w[i] / sum;
}

// ============================================================
//  SYMMETRIC: DYY
// ============================================================
void Grid::compute_dyy_symmetry()
{
    int M    = Ny[0];
    int half = M / 2;
    vector<double> wy(M, 0.0);
    double r = r_growth;

    for (int j = 0; j < half; ++j)
        wy[j] = pow(r, j);

    for (int j = 0; j < half; ++j)
        wy[half + half-1- j] = pow(r, j);

    if (M % 2 != 0)
        wy[half] = wy[half - 1];

    double sum = 0.0;
    for (int j = 0; j < M; ++j) sum += wy[j];
    for (int j = 0; j < M; ++j) dyy[j] = H * wy[j] / sum;
}

// ============================================================
//  Build node coordinates from cell sizes
// ============================================================
void Grid::build_coordinates()
{
    for (size_t i = 1; i < x_c.size(); ++i)
        x_c[i] = x_c[i - 1] + dxx[i - 1];

    for (size_t j = 1; j < y_c.size(); ++j)
        y_c[j] = y_c[j - 1] + dyy[j - 1];
}

// ============================================================
//  Write Mesh.dat
// ============================================================
void Grid::mesh_generation()
{
    int Nx_total = 0, Ny_total = 0;
    switch (meshType) {
        case MeshType::NONUNIFORM_CYL:
            Nx_total = Nx[0] + Nx[1] + Nx[2];
            Ny_total = Ny[0] + Ny[1] + Ny[2];
            break;
        default:
            Nx_total = Nx[0];
            Ny_total = Ny[0];
            break;
    }

    ofstream f("Mesh.dat");
    for (int j = 0; j <= Ny_total; ++j) {
        for (int i = 0; i <= Nx_total; ++i)
            f << x_c[i] << '\t' << y_c[j] << '\n';
        f << '\n';
    }
    cout << "[Grid] Mesh.dat written — "
         << Nx_total << " x " << Ny_total << " cells.\n";
}

// ============================================================
//  Accessors
// ============================================================
double Grid::getX(int i) const { return x_c[i]; }
double Grid::getY(int j) const { return y_c[j]; }
double Grid::dx(int i)   const { return dxx[i]; }
double Grid::dy(int j)   const { return dyy[j]; }
int    Grid::getNx()     const { return static_cast<int>(x_c.size()) - 1; }
int    Grid::getNy()     const { return static_cast<int>(y_c.size()) - 1; }
double Grid::getD()      const { return 1.0; }
