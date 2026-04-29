#ifndef GRID_H
#define GRID_H

#include <vector>

class Grid {
public:

    // ── Mesh type ────────────────────────────────────────────────────────────
    enum class MeshType {
        UNIFORM,              // equal cell sizes everywhere
        NONUNIFORM_CYL,       // 3-block stretching around a square cylinder
        NONUNIFORM_Symmetry   // symmetric wall refinement 
    };

    // ── Constructor 1: Uniform ───────────────────────────────────────────────
    Grid(int Nx_total, int Ny_total,
         double Ld, double Hd);

    // ── Constructor 2: Non-uniform cylinder ─────────────────────────────────
    //   3-block stretching in x and y, refined around (Xc,Yc) with size D.
    //   Middle block node count is auto-computed from FC.
    Grid(int Nx1, int Nx3,
         int Ny1, int Ny3,
         double Ld, double Hd, double FC,
         double Xc, double Yc, double D);

    // ── Constructor 3: Non-uniform symmetric ────────────────────────────────
    //   Cells grow from each wall toward the centre.
    //   growth_rate > 1  
    Grid(int Nx_total, int Ny_total,
         double Ld, double Hd,
         double growth_rate);

    // ── Accessors (unchanged interface) ─────────────────────────────────────
    double getX(int i) const;
    double getY(int j) const;
    double dx(int i)   const;
    double dy(int j)   const;
    int    getNx()     const;
    int    getNy()     const;
    double getD()      const;

    void   mesh_generation();

public:
    MeshType meshType;

    std::vector<int> Nx, Ny;
    double L, H;
    double First_cell_size;
    double r_growth;

    std::vector<double> x_c, y_c;
    std::vector<double> dxx,  dyy;

    std::vector<double> dl;
    std::vector<double> dH;

private:
    void   Middle_block_size();
    double solve_r_bisection(double Lseg, double dx1, int N);
    double stretch(double Lseg, int N, int i, double dx1);

    // cylinder mesh
    void   compute_dxx_cyl();
    void   compute_dyy_cyl();

    // uniform mesh
    void   compute_uniform_dxx();
    void   compute_uniform_dyy();

    // symmetric mesh
    void   compute_dxx_symmetry();
    void   compute_dyy_symmetry();

    void   build_coordinates();
};

#endif
