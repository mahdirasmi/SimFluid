#include "Solver.h"
#include "Vectors.h"
#include "Cylinder.h"
#include "Property.h"
#include "Grid.h"
#include "BoundaryConditions.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <complex>
#include <functional>

using namespace std;

Solver::Solver(Grid &grid, Property &fluidprop, Vectors &vec, Cylinder &cyl,bool AT_flag, const BoundaryConditions& bc, bool cylinder_present, bool Obc)
    : Domain_grid(grid),
      fluid_property(fluidprop),
      field_vectors(vec),
      Square_cylinder(cyl),
      BC(bc),
      AT(AT_flag),
      has_cylinder(cylinder_present),
      open_BC(Obc)
{

    double rho    = fluid_property.getrho();
    double mu     = fluid_property.getviscosity();
    double min_dx = Domain_grid.dx(0);
    int    Nx     = Domain_grid.getNx();
    int    Ny     = Domain_grid.getNy();
    for (int i = 0; i < Nx; ++i) min_dx = std::min(min_dx, Domain_grid.dx(i));
    for (int j = 0; j < Ny; ++j) min_dx = std::min(min_dx, Domain_grid.dy(j));

    double nu = mu / rho;

    if (cylinder_present || Obc) {
    
        initialdt = 0.01;
    } else {
        
        initialdt = CFL_diff * min_dx * min_dx / nu;
        double Uref = fluid_property.getinlet_average();
        if (Uref > 1e-12)
            initialdt = std::min(initialdt, CFD_conv * min_dx / Uref);
        initialdt = std::min(initialdt, 0.01);
    }
}

double Solver::getinitialdt() const { return initialdt; }

/*****************  Boundary Condition  *****************/
void Solver::Boundary_condition(double dtt)
{
 
    apply_top_bc   (dtt);
    apply_bottom_bc(dtt);
    apply_left_bc (dtt);
    apply_right_bc(dtt);

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    for (int j = 1; j < Ny - 1; ++j)
    for (int i = 1; i < Nx; ++i)
    if (Square_cylinder.is_insideU(Domain_grid.getX(i), Domain_grid.getY(j), j, i))
        field_vectors.u[j][i] = 0;

    for (int j = 1; j < Ny; ++j)
    for (int i = 1; i < Nx - 1; ++i)
    if (Square_cylinder.is_insideV(Domain_grid.getX(i), Domain_grid.getY(j), j, i))
        field_vectors.v[j][i] = 0;
}

// ─────────────────────────────────────────────────────────────────────────────
//  LEFT wall 
// ─────────────────────────────────────────────────────────────────────────────
void Solver::apply_left_bc(double dtt)
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double dt = dtt;                        
    const WallBC& B = BC.left;

    switch (B.type) {

        case BCType::DIRICHLET:
       
            for (int j = 0; j < Ny; ++j) {
                field_vectors.u[j][0] = fluid_property.get_inlet_u(j);
                field_vectors.v[j][0] = B.value_v;
            }
            break;

        case BCType::NO_SLIP:
            for (int j = 0; j < Ny; ++j) {
                field_vectors.u[j][0] = 0.0;
                field_vectors.v[j][0] = 0.0;
            }
            break;

        case BCType::NEUMANN:
            for (int j = 0; j < Ny; ++j) {         
                field_vectors.u[j][0] = field_vectors.u[j][1];
                field_vectors.v[j][0] = field_vectors.v[j][1];
            }
            break;

        case BCType::SYMMETRY:
            for (int j = 0; j < Ny; ++j) {        
                field_vectors.u[j][0] = 0.0;
                field_vectors.v[j][0] = field_vectors.v[j][1];
            }
            break;

        case BCType::CONVECTIVE: {
            double U_C = fluid_property.getinlet_average();
            for (int j = 0; j < Ny; ++j) {         
                field_vectors.u_next[j][0] = field_vectors.u[j][0]
                    - U_C * dt / Domain_grid.dx(0)
                    * (field_vectors.u[j][1] - field_vectors.u[j][0]);
            }
            double dx_avg = 0.5 * (Domain_grid.dx(0) + Domain_grid.dx(1));
            for (int j = 0; j < Ny; ++j) {         
                field_vectors.v[j][0]      = 0.0;
                field_vectors.v_next[j][0] = field_vectors.v[j][0]
                    - U_C * dt / dx_avg
                    * (field_vectors.v[j][1] - field_vectors.v[j][0]);
            }
            break;                                 
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  RIGHT wall  
// ─────────────────────────────────────────────────────────────────────────────
void Solver::apply_right_bc(double dtt)
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double dt = dtt;                            
    const WallBC& B = BC.right;

    switch (B.type) {

        case BCType::DIRICHLET:
            for (int j = 0; j < Ny; ++j) {
                field_vectors.u[j][Nx]     = B.value_u;
                field_vectors.v[j][Nx - 1] = B.value_v;
            }
            break;

        case BCType::NO_SLIP:
            for (int j = 0; j < Ny; ++j) {
                field_vectors.u[j][Nx]     = 0.0;
                field_vectors.v[j][Nx - 1] = 0.0;
            }
            break;

        case BCType::NEUMANN:
            for (int j = 0; j < Ny; ++j) {
                field_vectors.u_next[j][Nx]     = field_vectors.u[j][Nx - 1];
                field_vectors.v_next[j][Nx - 1] = field_vectors.v[j][Nx - 2];
            }
            break;

        case BCType::SYMMETRY:
            for (int j = 0; j < Ny; ++j) {     
                field_vectors.u[j][Nx]     = 0.0;
                field_vectors.v[j][Nx - 1] = field_vectors.v[j][Nx - 2];
            }
            break;

        case BCType::CONVECTIVE: {
            double U_C = fluid_property.getinlet_average();
            for (int j = 0; j < Ny; ++j) {         
                field_vectors.u_next[j][Nx] = field_vectors.u[j][Nx]
                    - U_C * dt / Domain_grid.dx(Nx - 1)
                    * (field_vectors.u[j][Nx] - field_vectors.u[j][Nx - 1]);
            }
            double dx_avg = 0.5 * (Domain_grid.dx(Nx - 2) + Domain_grid.dx(Nx - 1)); 
            for (int j = 0; j < Ny; ++j) {        
                field_vectors.v_next[j][Nx - 1] = field_vectors.v[j][Nx - 1]
                    - U_C * dt / dx_avg
                    * (field_vectors.v[j][Nx - 1] - field_vectors.v[j][Nx - 2]);
            }
            break;                              
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  BOTTOM wall  (y = 0,  j = 0)
// ─────────────────────────────────────────────────────────────────────────────
void Solver::apply_bottom_bc(double dtt)
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double dt = dtt;                               
    const WallBC& B = BC.bottom;

    switch (B.type) {

        case BCType::DIRICHLET:
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u[0][i] = B.value_u;
            for (int i = 0; i < Nx; ++i)
                field_vectors.v[0][i] = B.value_v;
            break;

        case BCType::NO_SLIP:
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u[0][i] = 0.0;
            for (int i = 0; i < Nx; ++i)
                field_vectors.v[0][i] = 0.0;
            break;

        case BCType::NEUMANN:
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u[0][i] = field_vectors.u[1][i];   
            for (int i = 0; i < Nx; ++i)
                field_vectors.v[0][i] = field_vectors.v[1][i];   
            break;

        case BCType::SYMMETRY:
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u[0][i] = field_vectors.u[1][i];
            for (int i = 0; i < Nx; ++i)
                field_vectors.v[0][i] = 0.0;                  
            break;

        case BCType::CONVECTIVE: {
            double U_C = fluid_property.getinlet_average();
            double dy_avg = 0.5 * (Domain_grid.dy(0) + Domain_grid.dy(1));
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u_next[0][i] = field_vectors.u[0][i]
                    - U_C * dt / dy_avg
                    * (field_vectors.u[0][i] - field_vectors.u[1][i]); 
            for (int i = 0; i < Nx; ++i)
                field_vectors.v_next[0][i] = field_vectors.v[0][i]
                    - U_C * dt / Domain_grid.dy(0)
                    * (field_vectors.v[0][i] - field_vectors.v[1][i]); 
            break;                                
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  TOP wall  (y = H,  j = Ny-1)
// ─────────────────────────────────────────────────────────────────────────────
void Solver::apply_top_bc(double dtt)
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double dt = dtt;                         
    const WallBC& B = BC.top;

    switch (B.type) {

        case BCType::DIRICHLET:
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u[Ny - 1][i] = B.value_u;
            for (int i = 0; i < Nx; ++i)
                field_vectors.v[Ny][i] = B.value_v;
            break;

        case BCType::NO_SLIP:
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u[Ny - 1][i] = 0.0;
            for (int i = 0; i < Nx; ++i)
                field_vectors.v[Ny][i] = 0.0;
            break;

        case BCType::NEUMANN:
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u[Ny - 1][i] = field_vectors.u[Ny - 2][i];
            for (int i = 0; i < Nx; ++i)
                field_vectors.v[Ny][i] = field_vectors.v[Ny - 1][i];
            break;

        case BCType::SYMMETRY:
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u[Ny - 1][i] = field_vectors.u[Ny - 2][i];
            for (int i = 0; i < Nx; ++i)
                field_vectors.v[Ny][i] = 0.0;
            break;

        case BCType::CONVECTIVE: {
            double V_C = fluid_property.getinlet_average();
            double dy_avg = 0.5 * (Domain_grid.dy(Ny - 1) + Domain_grid.dy(Ny - 2));
            for (int i = 0; i <= Nx; ++i)
                field_vectors.u_next[Ny - 1][i] = field_vectors.u[Ny - 1][i]
                    - V_C * dt / dy_avg
                    * (field_vectors.u[Ny - 1][i] - field_vectors.u[Ny - 2][i]);
            for (int i = 0; i < Nx; ++i)
                field_vectors.v_next[Ny][i] = field_vectors.v[Ny][i]
                    - V_C * dt / Domain_grid.dy(Ny - 1)
                    * (field_vectors.v[Ny][i] - field_vectors.v[Ny - 1][i]);
            break;                                  // FIX: missing break
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────

/***************** Interior Velocity Update *****************/
void Solver::Interrior_Velocity(double timestep)
{

    double dt = timestep;
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double rho = fluid_property.getrho();
    double mu  = fluid_property.getviscosity();

    // Adams-Bashforth weights: Forward Euler on first step, AB2 afterwards
    double ab1 = 1.5, ab2 = -0.5;
    if (step_count <= 1) { ab1 = 1.0;  ab2 = 0.0; }  // Forward Euler for first two steps

    /********** u predictor **********/
    for (int j = 1; j < Ny - 1; ++j)
    {
        for (int i = 1; i < Nx; ++i)
        {

            double xu = Domain_grid.getX(i);
            double yu = Domain_grid.getY(j);

            if (Square_cylinder.is_insideU(xu, yu, j, i))
            {
                field_vectors.uupp[j][i] = 0;
                field_vectors.u[j][i]    = 0;
                continue;
            }

            double vn = (Square_cylinder.is_insideV(0, 0, j + 1, i - 1) || Square_cylinder.is_insideV(0, 0, j + 1, i)) ? 0.0
                      : 0.5 * (field_vectors.v[j + 1][i - 1] + field_vectors.v[j + 1][i]);

            double vs = (Square_cylinder.is_insideV(0, 0, j, i - 1) || Square_cylinder.is_insideV(0, 0, j, i)) ? 0.0
                      : 0.5 * (field_vectors.v[j][i - 1] + field_vectors.v[j][i]);

            double ue = Square_cylinder.is_insideU(0, 0, j, i + 1) ? 0.0 : 0.5 * (field_vectors.u[j][i + 1] + field_vectors.u[j][i]);
            double uw = Square_cylinder.is_insideU(0, 0, j, i - 1) ? 0.0 : 0.5 * (field_vectors.u[j][i - 1] + field_vectors.u[j][i]);
            double un = Square_cylinder.is_insideU(0, 0, j + 1, i) ? 0.0 : 0.5 * (field_vectors.u[j + 1][i] + field_vectors.u[j][i]);
            double us = Square_cylinder.is_insideU(0, 0, j - 1, i) ? 0.0 : 0.5 * (field_vectors.u[j - 1][i] + field_vectors.u[j][i]);

            double Ae = Domain_grid.dy(j);
            double Aw = Domain_grid.dy(j);
            double An = (Domain_grid.dx(i - 1) + Domain_grid.dx(i)) / 2;
            double As = An;

            double m_e = rho * ue * Ae;
            double m_w = rho * uw * Aw;
            double m_n = rho * vn * An;
            double m_s = rho * vs * As;

            double U_N = field_vectors.u[j + 1][i];
            double U_S = field_vectors.u[j - 1][i];
            double U_E = field_vectors.u[j][i + 1];
            double U_W = field_vectors.u[j][i - 1];
            double U_P = field_vectors.u[j][i];

            double d_pe = Domain_grid.dx(i);
            double d_pw = Domain_grid.dx(i - 1);
            double d_pn = (Domain_grid.dy(j) + Domain_grid.dy(j + 1)) / 2;
            double d_ps = (Domain_grid.dy(j) + Domain_grid.dy(j - 1)) / 2;

            double diff_e = mu * (U_E - U_P) * Ae / d_pe;
            double diff_w = mu * (U_P - U_W) * Aw / d_pw;
            double diff_n = mu * (U_N - U_P) * An / d_pn;
            double diff_s = mu * (U_P - U_S) * As / d_ps;

            double vn_n_1 = (Square_cylinder.is_insideV(0, 0, j + 1, i - 1) || Square_cylinder.is_insideV(0, 0, j + 1, i)) ? 0.0
                          : 0.5 * (field_vectors.v_old[j + 1][i - 1] + field_vectors.v_old[j + 1][i]);

            double vs_n_1 = (Square_cylinder.is_insideV(0, 0, j, i - 1) || Square_cylinder.is_insideV(0, 0, j, i)) ? 0.0
                          : 0.5 * (field_vectors.v_old[j][i - 1] + field_vectors.v_old[j][i]);

            double ue_n_1 = Square_cylinder.is_insideU(0, 0, j, i + 1) ? 0.0 : 0.5 * (field_vectors.u_old[j][i + 1] + field_vectors.u_old[j][i]);
            double uw_n_1 = Square_cylinder.is_insideU(0, 0, j, i - 1) ? 0.0 : 0.5 * (field_vectors.u_old[j][i - 1] + field_vectors.u_old[j][i]);
            double un_n_1 = Square_cylinder.is_insideU(0, 0, j + 1, i) ? 0.0 : 0.5 * (field_vectors.u_old[j + 1][i] + field_vectors.u_old[j][i]);
            double us_n_1 = Square_cylinder.is_insideU(0, 0, j - 1, i) ? 0.0 : 0.5 * (field_vectors.u_old[j - 1][i] + field_vectors.u_old[j][i]);

            double m_e_n_1 = rho * ue_n_1 * Ae;
            double m_w_n_1 = rho * uw_n_1 * Aw;
            double m_n_n_1 = rho * vn_n_1 * An;
            double m_s_n_1 = rho * vs_n_1 * As;

            double U_N_n_1 = field_vectors.u_old[j + 1][i];
            double U_S_n_1 = field_vectors.u_old[j - 1][i];
            double U_E_n_1 = field_vectors.u_old[j][i + 1];
            double U_W_n_1 = field_vectors.u_old[j][i - 1];
            double U_P_n_1 = field_vectors.u_old[j][i];

            double diff_e_n_1 = mu * (U_E_n_1 - U_P_n_1) * Ae / d_pe;
            double diff_w_n_1 = mu * (U_P_n_1 - U_W_n_1) * Aw / d_pw;
            double diff_n_n_1 = mu * (U_N_n_1 - U_P_n_1) * An / d_pn;
            double diff_s_n_1 = mu * (U_P_n_1 - U_S_n_1) * As / d_ps;

            double gamma = ((Domain_grid.dx(i - 1) + Domain_grid.dx(i)) / 2) * Domain_grid.dy(j);

            double Ru_x     = -(m_e * ue - m_w * uw + m_n * un - m_s * us)     + (diff_e - diff_w + diff_n  - diff_s);
            double Ru_x_n_1 = -(m_e_n_1 * ue_n_1 - m_w_n_1 * uw_n_1 + m_n_n_1 * un_n_1 - m_s_n_1 * us_n_1) + (diff_e_n_1 - diff_w_n_1 + diff_n_n_1 - diff_s_n_1);

            field_vectors.uupp[j][i] = U_P + (dt / (rho * gamma)) * (ab1 * Ru_x + ab2 * Ru_x_n_1);
        }
    }


    for (int j = 0; j < Ny; ++j)
        field_vectors.uupp[j][Nx] = field_vectors.u_next[j][Nx];

    for (int j = 1; j < Ny - 1; ++j)
    {
        for (int i = 1; i < Nx - 1; ++i)
        {

            double xv = Domain_grid.getX(i);
            double yv = Domain_grid.getY(j);

            if (Square_cylinder.is_insideV(xv, yv, j, i))
            {
                field_vectors.vvpp[j][i] = 0;
                field_vectors.v[j][i]    = 0;
                continue;
            }

            double ue = (Square_cylinder.is_insideU(0, 0, j, i + 1) || Square_cylinder.is_insideU(0, 0, j - 1, i + 1)) ? 0.0
                      : 0.5 * (field_vectors.u[j][i + 1] + field_vectors.u[j - 1][i + 1]);

            double uw = (Square_cylinder.is_insideU(0, 0, j, i) || Square_cylinder.is_insideU(0, 0, j - 1, i)) ? 0.0
                      : 0.5 * (field_vectors.u[j][i] + field_vectors.u[j - 1][i]);

            double vn = Square_cylinder.is_insideV(0, 0, j + 1, i) ? 0.0 : 0.5 * (field_vectors.v[j + 1][i] + field_vectors.v[j][i]);
            double vs = Square_cylinder.is_insideV(0, 0, j - 1, i) ? 0.0 : 0.5 * (field_vectors.v[j - 1][i] + field_vectors.v[j][i]);
            double ve = Square_cylinder.is_insideV(0, 0, j, i + 1) ? 0.0 : 0.5 * (field_vectors.v[j][i + 1] + field_vectors.v[j][i]);
            double vw = Square_cylinder.is_insideV(0, 0, j, i - 1) ? 0.0 : 0.5 * (field_vectors.v[j][i - 1] + field_vectors.v[j][i]);

            double Ae = (Domain_grid.dy(j) + Domain_grid.dy(j - 1)) / 2;
            double Aw = Ae;
            double An = Domain_grid.dx(i);
            double As = An;

            double m_e = rho * ue * Ae;
            double m_w = rho * uw * Aw;
            double m_n = rho * vn * An;
            double m_s = rho * vs * As;

            double V_N = field_vectors.v[j + 1][i];
            double V_S = field_vectors.v[j - 1][i];
            double V_E = field_vectors.v[j][i + 1];
            double V_W = field_vectors.v[j][i - 1];
            double V_P = field_vectors.v[j][i];

            double d_pe = (Domain_grid.dx(i) + Domain_grid.dx(i + 1)) / 2;
            double d_pw = (Domain_grid.dx(i) + Domain_grid.dx(i - 1)) / 2;
            double d_pn = Domain_grid.dy(j);
            double d_ps = Domain_grid.dy(j - 1);

            double diff_e = mu * (V_E - V_P) * Ae / d_pe;
            double diff_w = mu * (V_P - V_W) * Aw / d_pw;
            double diff_n = mu * (V_N - V_P) * An / d_pn;
            double diff_s = mu * (V_P - V_S) * As / d_ps;

            double ue_n_1 = (Square_cylinder.is_insideU(0, 0, j, i + 1) || Square_cylinder.is_insideU(0, 0, j - 1, i + 1)) ? 0.0
                          : 0.5 * (field_vectors.u_old[j][i + 1] + field_vectors.u_old[j - 1][i + 1]);

            double uw_n_1 = (Square_cylinder.is_insideU(0, 0, j, i) || Square_cylinder.is_insideU(0, 0, j - 1, i)) ? 0.0
                          : 0.5 * (field_vectors.u_old[j][i] + field_vectors.u_old[j - 1][i]);

            double vn_n_1 = Square_cylinder.is_insideV(0, 0, j + 1, i) ? 0.0 : 0.5 * (field_vectors.v_old[j + 1][i] + field_vectors.v_old[j][i]);
            double vs_n_1 = Square_cylinder.is_insideV(0, 0, j - 1, i) ? 0.0 : 0.5 * (field_vectors.v_old[j - 1][i] + field_vectors.v_old[j][i]);
            double ve_n_1 = Square_cylinder.is_insideV(0, 0, j, i + 1) ? 0.0 : 0.5 * (field_vectors.v_old[j][i + 1] + field_vectors.v_old[j][i]);
            double vw_n_1 = Square_cylinder.is_insideV(0, 0, j, i - 1) ? 0.0 : 0.5 * (field_vectors.v_old[j][i - 1] + field_vectors.v_old[j][i]);

            double m_e_n_1 = rho * ue_n_1 * Ae;
            double m_w_n_1 = rho * uw_n_1 * Ae;
            double m_n_n_1 = rho * vn_n_1 * An;
            double m_s_n_1 = rho * vs_n_1 * As;

            double V_N_n_1 = field_vectors.v_old[j + 1][i];
            double V_S_n_1 = field_vectors.v_old[j - 1][i];
            double V_E_n_1 = field_vectors.v_old[j][i + 1];
            double V_W_n_1 = field_vectors.v_old[j][i - 1];
            double V_P_n_1 = field_vectors.v_old[j][i];

            double diff_e_n_1 = mu * (V_E_n_1 - V_P_n_1) * Ae / d_pe;
            double diff_w_n_1 = mu * (V_P_n_1 - V_W_n_1) * Aw / d_pw;
            double diff_n_n_1 = mu * (V_N_n_1 - V_P_n_1) * An / d_pn;
            double diff_s_n_1 = mu * (V_P_n_1 - V_S_n_1) * As / d_ps;

            double gamma = Domain_grid.dx(i) * (Domain_grid.dy(j) + Domain_grid.dy(j - 1)) / 2;

            double Ru_y     = -(m_e * ve   - m_w * vw  + m_n * vn  - m_s * vs)     + (diff_e - diff_w + diff_n - diff_s);
            double Ru_y_n_1 = -(m_e_n_1 * ve_n_1 - m_w_n_1 * vw_n_1 + m_n_n_1 * vn_n_1 - m_s_n_1 * vs_n_1) + (diff_e_n_1 - diff_w_n_1 + diff_n_n_1 - diff_s_n_1);

            field_vectors.vvpp[j][i] = V_P + (dt / (rho * gamma)) * (ab1 * Ru_y + ab2 * Ru_y_n_1);
        }
    }
}

void Solver::pressure_boundary()
{
}

/********************** Gauss-Seidel ***************************/
void Solver::Guessiedel(double dtt)
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double rho = fluid_property.getrho();
    double dt  = dtt;

    for (int iter = 0; iter < MaxIter; ++iter)
    {
        double maxerror = 0;
    if (has_cylinder){
    for (int j = 1; j < Ny - 1; ++j)
    {
    for (int i = 1; i < Nx - 1; ++i)
    {
        bool left_face   = false;
        bool right_face  = false;
        bool bottom_face = false;
        bool top_face    = false;

        Square_cylinder.get_faces(0, 0, left_face, right_face, bottom_face, top_face, j, i);
          if (!(left_face || right_face || bottom_face || top_face))
                    continue;

                if (left_face)   field_vectors.P_next[j][i + 1] = field_vectors.P_next[j][i];
                if (right_face)  field_vectors.P_next[j][i - 1] = field_vectors.P_next[j][i];
                if (bottom_face) field_vectors.P_next[j + 1][i] = field_vectors.P_next[j][i];
                if (top_face)    field_vectors.P_next[j - 1][i] = field_vectors.P_next[j][i];
            }
        }
    }

    // ── Gauss-Seidel pressure sweep ──────────────────────────────────────
    
        if (has_cylinder|| open_BC) {
            // ── OPEN DOMAIN: interior sweep  ──────────────────
            for (int j = 1; j < Ny - 1; ++j)
            {
                for (int i = 1; i < Nx - 1; ++i)
                {
                    double x = Domain_grid.getX(i) + Domain_grid.dx(i) / 2;
                    double y = Domain_grid.getY(j) + Domain_grid.dy(j) / 2;

                    if (Square_cylinder.is_inside(x, y, j, i)) continue;

                    double Ae=Domain_grid.dy(j), Aw=Ae, An=Domain_grid.dx(i), As=An;
                    double d_pe=(Domain_grid.dx(i)+Domain_grid.dx(i+1))/2;
                    double d_pw=(Domain_grid.dx(i)+Domain_grid.dx(i-1))/2;
                    double d_pn=(Domain_grid.dy(j)+Domain_grid.dy(j+1))/2;
                    double d_ps=(Domain_grid.dy(j)+Domain_grid.dy(j-1))/2;
                    double aE=Ae/d_pe, aW=Aw/d_pw, aN=An/d_pn, aS=As/d_ps;
                    double ap=aE+aW+aN+aS;

                    double pN=field_vectors.P_next[j+1][i], pS=field_vectors.P_next[j-1][i];
                    double pE=field_vectors.P_next[j][i+1], pW=field_vectors.P_next[j][i-1];

                    double vpN=field_vectors.vvpp[j+1][i], vpS=field_vectors.vvpp[j][i];
                    double upE=field_vectors.uupp[j][i+1], upW=field_vectors.uupp[j][i];
                    double bp=-((rho*upE)*Ae-(rho*upW)*Aw+(rho*vpN)*An-(rho*vpS)*As)/dt;
                    double PP=(aN*pN+aS*pS+aE*pE+aW*pW+bp)/ap;
                  
                    double newP=field_vectors.P_next[j][i]*(1-omega)+omega*PP;

                    maxerror=max(maxerror,fabs(newP-field_vectors.P_next[j][i]));
                    field_vectors.P_next[j][i]=newP;
                }
            }
            // Open domain boundary pressures — Neumann top/bottom/left, P=0 outlet
            for (int i = 0; i < Nx; ++i) {
                field_vectors.P_next[0][i]      = field_vectors.P_next[1][i];
                field_vectors.P_next[Ny - 1][i] = field_vectors.P_next[Ny - 2][i];
            }
            for (int j = 0; j < Ny; ++j) {
                field_vectors.P_next[j][0]      = field_vectors.P_next[j][1];
                field_vectors.P_next[j][Nx - 1] = 0;   // P=0 at outlet
            }

        } else {

            // ── CLOSED DOMAIN: all-cells sweep with Neumann-flag coefficients ─
            //
            // Key differences from the open-domain path:
            //   1. No under-relaxation (omega=1)
            //   2. Convergence is checked on AVERAGE error over all cells, not  maximum error.
            
            
            double total_error_cav = 0.0;

            for (int j = 0; j < Ny; ++j)
            {
                for (int i = 0; i < Nx; ++i)
                {
                    int an=(j==Ny-1)?0:1, as=(j==0)?0:1;
                    int ae=(i==Nx-1)?0:1, aw=(i==0)?0:1;
                    int ap_count = an+as+ae+aw;

                    if (ap_count == 0) continue;

                    double Ae=Domain_grid.dy(j), Aw=Ae, An=Domain_grid.dx(i), As=An;
                    double d_pe=(ae)?(Domain_grid.dx(i)+Domain_grid.dx(i+1))/2:Domain_grid.dx(i);
                    double d_pw=(aw)?(Domain_grid.dx(i)+Domain_grid.dx(i-1))/2:Domain_grid.dx(i);
                    double d_pn=(an)?(Domain_grid.dy(j)+Domain_grid.dy(j+1))/2:Domain_grid.dy(j);
                    double d_ps=(as)?(Domain_grid.dy(j)+Domain_grid.dy(j-1))/2:Domain_grid.dy(j);

                    double aE=ae*Ae/d_pe, aW=aw*Aw/d_pw, aN=an*An/d_pn, aS=as*As/d_ps;
                    double aP=aE+aW+aN+aS;

                    // Neumann wall (dP/dn = 0)
                    double pN=(an)?field_vectors.P_next[j+1][i]:field_vectors.P_next[j][i];
                    double pS=(as)?field_vectors.P_next[j-1][i]:field_vectors.P_next[j][i];
                    double pE=(ae)?field_vectors.P_next[j][i+1]:field_vectors.P_next[j][i];
                    double pW=(aw)?field_vectors.P_next[j][i-1]:field_vectors.P_next[j][i];

                    double vpN=field_vectors.vvpp[j+1][i], vpS=field_vectors.vvpp[j][i];
                    double upE=field_vectors.uupp[j][i+1], upW=field_vectors.uupp[j][i];


                    double rhs = (rho * Ae * An / dt) * (vpN - vpS + upE - upW);
                    double bp=-((rho*upE)*Ae-(rho*upW)*Aw+(rho*vpN)*An-(rho*vpS)*As)/dt;

                   
                 
                    double PP = (aN*pN + aS*pS + aE*pE + aW*pW +bp) / aP;

                    double diff = fabs(PP - field_vectors.P_next[j][i]);
                    total_error_cav += diff;
                    maxerror = max(maxerror, diff);

                    field_vectors.P_next[j][i] = PP;
                }
            }

            field_vectors.P_next[0][0] = 0.0;
            double avg_error = total_error_cav / (Nx * Ny);
            field_vectors.P = field_vectors.P_next;
            report_error = avg_error;

            if (avg_error < Max_error)
                break;

            continue;
        }

        field_vectors.P  = field_vectors.P_next;
        report_error = maxerror;

        if (maxerror < Max_error)
            break;
    }
}

// Should be checked ..........?!!
void Solver::ConjugateGradient(double dtt)
{
    int    Nx  = Domain_grid.getNx();
    int    Ny  = Domain_grid.getNy();
    double rho = fluid_property.getrho();
    double dt  = dtt;


    const bool open_path = has_cylinder || open_BC;

    std::vector<std::pair<int,int>> fluid_cells;
    fluid_cells.reserve(Nx * Ny);

    if (open_path) {

        for (int j = 1; j < Ny - 1; ++j)
            for (int i = 1; i < Nx - 1; ++i) {
                double x = Domain_grid.getX(i) + 0.5*Domain_grid.dx(i);
                double y = Domain_grid.getY(j) + 0.5*Domain_grid.dy(j);
                if (!Square_cylinder.is_inside(x, y, j, i))
                    fluid_cells.push_back({j, i});
            }
    } else {
       
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                fluid_cells.push_back({j, i});
    }

    
    auto apply_boundary = [&](std::vector<std::vector<double>>& P)
    {
        if (open_path) {
        
            if (has_cylinder) {
                for (int j = 1; j < Ny - 1; ++j)
                    for (int i = 1; i < Nx - 1; ++i) {
                        bool lf=false, rf=false, bf=false, tf=false;
                        Square_cylinder.get_faces(0,0,lf,rf,bf,tf,j,i);
                        if (!(lf||rf||bf||tf)) continue;
                        if (lf) P[j][i+1] = P[j][i];
                        if (rf) P[j][i-1] = P[j][i];
                        if (bf) P[j+1][i] = P[j][i];
                        if (tf) P[j-1][i] = P[j][i];
                    }
            }
            // Domain boundary: Neumann top/bottom/left, P=0 right (outlet)
            for (int i = 0; i < Nx; ++i) {
                P[0][i]      = P[1][i];        // bottom Neumann
                P[Ny-1][i]   = P[Ny-2][i];     // top Neumann
            }
            for (int j = 0; j < Ny; ++j) {
                P[j][0]      = P[j][1];         // left Neumann
                P[j][Nx-1]   = 0.0;             // right: P=0 reference
            }
        } else {
            // Closed domain: Neumann on all four walls
            for (int i = 0; i < Nx; ++i) {
                P[0][i]      = P[1][i];
                P[Ny-1][i]   = P[Ny-2][i];
            }
            for (int j = 0; j < Ny; ++j) {
                P[j][0]      = P[j][1];
                P[j][Nx-1]   = P[j][Nx-2];
            }
            // Pressure anchor 
            P[0][0] = 0.0;
        }
    };


    auto stencil = [&](const std::vector<std::vector<double>>& p_field,
                             std::vector<std::vector<double>>& Ap_field)
    {
        for (auto& [j, i] : fluid_cells) {

            double Ae = Domain_grid.dy(j),  Aw = Ae;
            double An = Domain_grid.dx(i),  As = An;

            double aE, aW, aN, aS;

            if (open_path) {
                // All active cells are interior: neighbours always exist.
                // Exactly the same coefficient formula as the GS open sweep.
                double d_pe = (Domain_grid.dx(i) + Domain_grid.dx(i+1)) / 2;
                double d_pw = (Domain_grid.dx(i) + Domain_grid.dx(i-1)) / 2;
                double d_pn = (Domain_grid.dy(j) + Domain_grid.dy(j+1)) / 2;
                double d_ps = (Domain_grid.dy(j) + Domain_grid.dy(j-1)) / 2;
                aE = Ae/d_pe;  aW = Aw/d_pw;
                aN = An/d_pn;  aS = As/d_ps;
            } else {
                // Boundary cells use Neumann flags (0 at domain wall).
                int an=(j==Ny-1)?0:1, as=(j==0)?0:1;
                int ae=(i==Nx-1)?0:1, aw=(i==0)?0:1;
                double d_pe=(ae)?(Domain_grid.dx(i)+Domain_grid.dx(i+1))/2:Domain_grid.dx(i);
                double d_pw=(aw)?(Domain_grid.dx(i)+Domain_grid.dx(i-1))/2:Domain_grid.dx(i);
                double d_pn=(an)?(Domain_grid.dy(j)+Domain_grid.dy(j+1))/2:Domain_grid.dy(j);
                double d_ps=(as)?(Domain_grid.dy(j)+Domain_grid.dy(j-1))/2:Domain_grid.dy(j);
                aE=ae*Ae/d_pe; aW=aw*Aw/d_pw;
                aN=an*An/d_pn; aS=as*As/d_ps;
            }

            double ap = aE + aW + aN + aS;

            // Neighbour values — apply_boundary has already set ghost rows
            // so reading p_field[j±1][i] and p_field[j][i±1] is always safe.
            double pN = p_field[j+1][i], pS = p_field[j-1][i];
            double pE = p_field[j][i+1], pW = p_field[j][i-1];

            Ap_field[j][i] = ap*p_field[j][i] - aE*pE - aW*pW - aN*pN - aS*pS;
        }
    };

    // ── dot product over active cells only ───────────────────────────────────
    auto dot = [&](const std::vector<std::vector<double>>& A,
                   const std::vector<std::vector<double>>& B) -> double
    {
        double s = 0.0;
        for (auto& [j, i] : fluid_cells) s += A[j][i] * B[j][i];
        return s;
    };

    
    std::vector<std::vector<double>> b(Ny, std::vector<double>(Nx, 0.0));
    for (auto& [j, i] : fluid_cells) {
        double Ae = Domain_grid.dy(j), Aw = Ae;
        double An = Domain_grid.dx(i), As = An;
        double upE = field_vectors.uupp[j][i+1];   // valid: u has Nx+1 columns
        double upW = field_vectors.uupp[j][i];
        double vpN = field_vectors.vvpp[j+1][i];   // valid: v has Ny+1 rows
        double vpS = field_vectors.vvpp[j][i];
        b[j][i] = -((rho*upE)*Ae - (rho*upW)*Aw + (rho*vpN)*An - (rho*vpS)*As) / dt;
    }

    // ── CG iteration ─────────────────────────────────────────────────────────
    // Initialise from the current pressure guess (warm start).
    std::vector<std::vector<double>> p  = field_vectors.P_next;
    std::vector<std::vector<double>> Ap (Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> r  (Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> d  (Ny, std::vector<double>(Nx, 0.0));
    std::vector<std::vector<double>> Ad (Ny, std::vector<double>(Nx, 0.0));

    // Apply BCs to the initial guess, compute initial residual r = b - A·p
    apply_boundary(p);
    stencil(p, Ap);
    for (auto& [j, i] : fluid_cells) {
        r[j][i] = b[j][i] - Ap[j][i];
        d[j][i] = r[j][i];   // initial search direction = residual
    }

    double r_dot_r = dot(r, r);

    for (int iter = 0; iter < MaxIter; ++iter)
    {
            stencil(d, Ad);

        double d_dot_Ad = dot(d, Ad);
        if (fabs(d_dot_Ad) < 1e-30) break;

        double alpha = r_dot_r / d_dot_Ad;

        double maxerror = 0.0;
        for (auto& [j, i] : fluid_cells) {
            double p_old = p[j][i];
            p[j][i] += alpha * d[j][i];
            r[j][i] -= alpha * Ad[j][i];
            maxerror = max(maxerror, fabs(p[j][i] - p_old));
        }

        // Apply BCs to the updated solution
        apply_boundary(p);
        report_error = maxerror;

        // Convergence check — same criterion as GS open path (maxerror)
        if (maxerror < Max_error) break;

        double r_dot_r_new = dot(r, r);
        if (fabs(r_dot_r) < 1e-30) break;
        double beta = r_dot_r_new / r_dot_r;
        r_dot_r = r_dot_r_new;
        for (auto& [j, i] : fluid_cells) d[j][i] = r[j][i] + beta * d[j][i];
    }

    field_vectors.P_next = p;
    field_vectors.P      = p;
}

/********************* Velocity Correction *************************/
void Solver::correction_velocity(double dtt)
{

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double rho = fluid_property.getrho();
    double dt  = dtt;

    for (int j = 1; j < Ny - 1; ++j)
    {
        for (int i = 1; i < Nx; ++i)
        {
            double xu = Domain_grid.getX(i);
            double yu = Domain_grid.getY(j);

            if (Square_cylinder.is_insideU(xu, yu, j, i))
            {
                field_vectors.u_next[j][i] = 0;
                continue;
            }

            field_vectors.u_next[j][i] = field_vectors.uupp[j][i] - dt * (field_vectors.P_next[j][i] - field_vectors.P_next[j][i - 1]) / (rho * (Domain_grid.dx(i) + Domain_grid.dx(i - 1)) / 2);
        }
    }

    for (int j = 1; j < Ny; ++j)
    {
        for (int i = 1; i < Nx - 1; ++i)
        {
            double xv = Domain_grid.getX(i);
            double yv = Domain_grid.getY(j);

            if (Square_cylinder.is_insideV(xv, yv, j, i))
            {
                field_vectors.v_next[j][i] = 0;
                continue;
            }

            field_vectors.v_next[j][i] = field_vectors.vvpp[j][i] - dt * (field_vectors.P_next[j][i] - field_vectors.P_next[j - 1][i]) / (rho * (Domain_grid.dy(j) + Domain_grid.dy(j - 1)) / 2);
        }
    }
}

/********************* MAIN SOLVER LOOP ***************************/
void Solver::solve(double ultimate_time)
{

    double sampling_Cdp=0, Time_averaged_Cdp=0;
    double sampling_Cdf=0, Time_averaged_Cdf=0;
    double sampling_Cl=0,  Time_averaged_Cl=0;
    double Tavg = 0;

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    double rho = fluid_property.getrho();
    double mu  = fluid_property.getviscosity();
    double Un  = fluid_property.getinlet_average();
    double dt  = initialdt;

    // Asymmetric perturbation to trigger vortex shedding (cylinder only)

    if (has_cylinder) {
        for (int j = 1; j < Ny; ++j)
            for (int i = 1; i < Nx - 1; ++i) {
                if (Square_cylinder.is_insideV(0, 0, j, i)) continue;
                double yv = Domain_grid.getY(j);
                if (yv >= 9.5 && yv <= 10.5)
                    field_vectors.v[j][i] += 0.01 * Un * sin(M_PI * (yv - 9.5));
            }
    }

    while (flowtime < ultimate_time)
    {
        field_vectors.u_old = field_vectors.u;
        field_vectors.v_old = field_vectors.v;

        Boundary_condition(dt);

        field_vectors.uupp = field_vectors.u;
        field_vectors.vvpp = field_vectors.v;

        Interrior_Velocity(dt);
        pressure_boundary();

        if (pressure_solver == 1)
            ConjugateGradient(dt);
        else
            Guessiedel(dt);

        correction_velocity(dt);

         if(has_cylinder || open_BC){
        // Re-apply right convective BC after correction
        if (BC.right.type == BCType::CONVECTIVE)
            for (int j = 0; j < Ny; ++j)
            field_vectors.u_next[j][Nx] = field_vectors.u[j][Nx] - (field_vectors.u[j][Nx] - field_vectors.u[j][Nx - 1]) * (Un * dt / Domain_grid.dx(Nx - 1));

        
        if (BC.left.type == BCType::DIRICHLET) {
            for (int j = 0; j < Ny; ++j) {
                field_vectors.u_next[j][0] = fluid_property.get_inlet_u(j);
                field_vectors.v_next[j][0] = BC.left.value_v;
            }
        } else if (BC.left.type == BCType::NO_SLIP) {
            for (int j = 0; j < Ny; ++j) {
                field_vectors.u_next[j][0] = 0.0;
                field_vectors.v_next[j][0] = 0.0;
            }
        } else if (BC.left.type == BCType::NEUMANN) {
            for (int j = 0; j < Ny; ++j)
                field_vectors.u_next[j][0] = field_vectors.u_next[j][1];
        }
 
        if (BC.bottom.type == BCType::NO_SLIP) {
            for (int i = 0; i <= Nx; ++i) field_vectors.u_next[0][i] = 0.0;
            for (int i = 0; i <  Nx; ++i) field_vectors.v_next[0][i] = 0.0;
        } else if (BC.bottom.type == BCType::SYMMETRY) {
            for (int i = 0; i <= Nx; ++i) field_vectors.u_next[0][i] = field_vectors.u_next[1][i];
            for (int i = 0; i <  Nx; ++i) field_vectors.v_next[0][i] = 0.0;
        } else if (BC.bottom.type == BCType::DIRICHLET) {
            for (int i = 0; i <= Nx; ++i) field_vectors.u_next[0][i] = BC.bottom.value_u;
            for (int i = 0; i <  Nx; ++i) field_vectors.v_next[0][i] = BC.bottom.value_v;
        }

        if (BC.top.type == BCType::NO_SLIP) {
            for (int i = 0; i <= Nx; ++i) field_vectors.u_next[Ny-1][i] = 0.0;
            for (int i = 0; i <  Nx; ++i) field_vectors.v_next[Ny][i]   = 0.0;
        } else if (BC.top.type == BCType::SYMMETRY) {
            for (int i = 0; i <= Nx; ++i) field_vectors.u_next[Ny-1][i] = field_vectors.u_next[Ny-2][i];
            for (int i = 0; i <  Nx; ++i) field_vectors.v_next[Ny][i]   = 0.0;
        } else if (BC.top.type == BCType::DIRICHLET) {
            for (int i = 0; i <= Nx; ++i) field_vectors.u_next[Ny-1][i] = BC.top.value_u;
            for (int i = 0; i <  Nx; ++i) field_vectors.v_next[Ny][i]   = BC.top.value_v;
        }
    }


        if (!has_cylinder && !open_BC) {
      
            for (int i = 0; i <= Nx; ++i) {
                field_vectors.u_next[0][i]      = 0.0;          // bottom u = 0
                field_vectors.u_next[Ny - 1][i] = BC.top.value_u; // lid u = Uref
            }
            for (int j = 0; j < Ny; ++j) {
                field_vectors.u_next[j][0]  = 0.0;   // left wall
                field_vectors.u_next[j][Nx] = 0.0;   // right wall
            }
            for (int i = 0; i < Nx; ++i) {
                field_vectors.v_next[0][i]  = 0.0;   // bottom v = 0
                field_vectors.v_next[Ny][i] = 0.0;   // top v = 0
            }
            for (int j = 0; j <= Ny; ++j) {
                field_vectors.v_next[j][0]      = 0.0;  // left v = 0
                field_vectors.v_next[j][Nx - 1] = 0.0;  // right v = 0
            }
        }

        field_vectors.u = field_vectors.u_next;
        field_vectors.v = field_vectors.v_next;

        ++step_count;

        if (has_cylinder) {
            if (flowtime >= 50.0)
            {
                for (int j = 0; j < Ny; ++j)
                    for (int i = 0; i <= Nx; ++i)
                        field_vectors.sumu[j][i] += field_vectors.u[j][i] * dt;

                double CD  = Drag_pressure();
                double CDF = Drag_friction();
                double CL  = Lift_coefficient();

                sampling_Cdp += CD  * dt;
                sampling_Cdf += CDF * dt;
                sampling_Cl  += CL  * dt;
                Tavg         += dt;
                Tavg_stored   = Tavg;

                Time_averaged_Cdp = sampling_Cdp / Tavg;
                Time_averaged_Cdf = sampling_Cdf / Tavg;
                Time_averaged_Cl  = sampling_Cl  / Tavg;

                CL_history.push_back(CL);
                Time_history.push_back(flowtime);
                Time_period += dt;

                cout << " Averaged Cdp = " << Time_averaged_Cdp
                     << "\t Averaged Cdf = " << Time_averaged_Cdf
                     << "\t Averaged total Cd = " << Time_averaged_Cdp + Time_averaged_Cdf
                     << "\t Averaged Cl = " << Time_averaged_Cl << endl;
            }
            else
            {
                cout << " Instant Cdp = " << Drag_pressure()
                     << "\t Instant Cdf = " << Drag_friction()
                     << "\t Instant Total Cd = " << Drag_pressure() + Drag_friction()
                     << "\t Instant Cl = " << Lift_coefficient() << endl;
            }
        }

        // CFL adaptive time step
        double min_stable_dt = initialdt;
        double Max_velocity  = 0.0;

        for (int j = 1; j < Ny; ++j)
            for (int i = 1; i < Nx - 1; ++i) {
                double uvel = 0.5 * (field_vectors.u[j][i] + field_vectors.u[j][i + 1]);
                double vvel = 0.5 * (field_vectors.v[j][i] + field_vectors.v[j + 1][i]);
                double V_local = sqrt(uvel * uvel + vvel * vvel);
                Max_velocity = max(Max_velocity, V_local);
                if (V_local > 1e-4) {
                    double local_dt_conv = CFD_conv * min(Domain_grid.dx(i), Domain_grid.dy(j)) / V_local;
                    double local_dt_dif  = CFL_diff  * min(Domain_grid.dx(i), Domain_grid.dy(j)) / (mu / rho);
                    min_stable_dt = min(min_stable_dt, min(local_dt_conv, local_dt_dif));
                }
            }

        dt = (Max_velocity > 1e-4) ? min_stable_dt : initialdt;

        Monior(flowtime, dt);

        // ── Steady-state detection for cavity ────────────────────────────────
     
        if (!has_cylinder && !open_BC) {
            double max_du = 0.0;
            for (int jj = 1; jj < Ny - 1; ++jj)
                for (int ii = 1; ii < Nx; ++ii)
                    max_du = max(max_du, fabs(field_vectors.u[jj][ii] - field_vectors.u_old[jj][ii]));
            if (max_du < 3e-7 && flowtime > 0.1) {
                cout << "[Cavity] Steady state at t=" << flowtime
                     << "  max_du=" << max_du << endl;
                flowtime += dt;
                break;
            }
        }

        if (AT && flowtime > 160.0)
            export_paraview_Cylinder(flowtime);

        flowtime += dt;
    }

}

/********************* Drag — Pressure *************************/
double Solver::Drag_pressure()
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double rho = fluid_property.getrho();
    double Uin = fluid_property.getinlet_average();
    double D   = Domain_grid.getD();
    double Fx  = 0.0;

    for (int j = 1; j < Ny - 1; ++j)
        for (int i = 1; i < Nx - 1; ++i) {
            bool left_face=false, right_face=false, bottom_face=false, top_face=false;
            Square_cylinder.get_faces(Domain_grid.getX(i), Domain_grid.getY(j), left_face, right_face, bottom_face, top_face, j, i);
            if (!(left_face || right_face || top_face || bottom_face)) continue;
            double dS=0, nx=0;
            if (left_face)  { nx = -1.0; dS = Domain_grid.dy(j); }
            if (right_face) { nx =  1.0; dS = Domain_grid.dy(j); }
            if (left_face || right_face)
                Fx += (-field_vectors.P[j][i] * nx) * dS;
        }

    return Cdp = (2.0 * Fx) / (rho * Uin * Uin * D);
}

/********************* Drag — Friction *************************/
double Solver::Drag_friction()
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double rho = fluid_property.getrho();
    double mu  = fluid_property.getviscosity();
    double Uin = fluid_property.getinlet_average();
    double D   = Domain_grid.getD();
    double Fx  = 0.0;

    for (int j = 1; j < Ny - 1; ++j)
        for (int i = 1; i < Nx - 1; ++i) {
            bool left_face=false, right_face=false, bottom_face=false, top_face=false;
            Square_cylinder.get_faces(Domain_grid.getX(i), Domain_grid.getY(j), left_face, right_face, bottom_face, top_face, j, i);
            if (!(left_face || right_face || bottom_face || top_face)) continue;
            if (bottom_face || top_face) {
                double dy_half  = Domain_grid.dy(j) / 2;
                double dS       = Domain_grid.dx(i);
                double u_centre = 0.5 * (field_vectors.u[j][i] + field_vectors.u[j][i + 1]);
                Fx += mu * (u_centre / dy_half) * dS;
            }
        }

    return Cdf = (2.0 * Fx) / (rho * Uin * Uin * D);
}

/********************* Lift Coefficient *************************/
double Solver::Lift_coefficient()
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double rho = fluid_property.getrho();
    double Uin = fluid_property.getinlet_average();
    double D   = Domain_grid.getD();
    double Fy  = 0.0;

    for (int j = 1; j < Ny - 1; ++j)
        for (int i = 1; i < Nx - 1; ++i) {
            bool left_face=false, right_face=false, bottom_face=false, top_face=false;
            Square_cylinder.get_faces(Domain_grid.getX(i), Domain_grid.getY(j), left_face, right_face, bottom_face, top_face, j, i);
            if (!(left_face || right_face || bottom_face || top_face)) continue;
            double dS = Domain_grid.dx(i);
            if (bottom_face) Fy += (-field_vectors.P[j][i] * (-1.0)) * dS;
            if (top_face)    Fy += (-field_vectors.P[j][i] * ( 1.0)) * dS;
        }

    return Cl = (2.0 * Fy) / (rho * Uin * Uin * D);
}

/********************* Generate Report — Square Cylinder *************************/
void Solver::cylinder_reports()
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    // ── Full-field ParaView VTS ───────────────────────────────────────────────
    {
        std::ofstream vtk("cylinder_results.vts");
        vtk << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "<StructuredGrid WholeExtent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n"
            << "<Piece Extent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n"
            << "<Points>\n"
            << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j <= Ny; ++j)
            for (int i = 0; i <= Nx; ++i)
                vtk << Domain_grid.getX(i) << " " << Domain_grid.getY(j) << " 0.0\n";
        vtk << "</DataArray>\n</Points>\n<CellData>\n"
            << "<DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << field_vectors.P[j][i] << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"u_velocity\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << 0.5*(field_vectors.u[j][i] + field_vectors.u[j][i+1]) << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"v_velocity\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << 0.5*(field_vectors.v[j][i] + field_vectors.v[j+1][i]) << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"VelocityMagnitude\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                double uc = 0.5*(field_vectors.u[j][i] + field_vectors.u[j][i+1]);
                double vc = 0.5*(field_vectors.v[j][i] + field_vectors.v[j+1][i]);
                vtk << std::sqrt(uc*uc + vc*vc) << "\n";
            }
        vtk << "</DataArray>\n</CellData>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
        vtk.close();
    }

    // ── Wake centreline u (instantaneous + time-averaged) — CSV ──────────────
    {
        int jj = 0;
        double best = 1e9;
        for (int j = 0; j < Ny; ++j) {
            double yc = Domain_grid.getY(j) + 0.5*Domain_grid.dy(j);
            double d  = std::abs(yc - 10.0);
            if (d < best) { best = d; jj = j; }
        }
        double tavg = (Tavg_stored > 0.0) ? Tavg_stored : flowtime;
        std::ofstream f("cylinder_wake_centreline.csv");
        f << "x,u_instant,u_mean\n";
        for (int i = 0; i <= Nx; ++i) {
            double u_inst = field_vectors.u[jj][i];
            double u_mean = (tavg > 0.0) ? field_vectors.sumu[jj][i] / tavg : u_inst;
            f << Domain_grid.getX(i) << "," << u_inst << "," << u_mean << "\n";
        }
        f.close();
    }

    std::cout << "Cylinder output: cylinder_results.vts  cylinder_wake_centreline.csv\n";
}

/********************* Cavity Post-Processing *************************/
void Solver::cavity_report()
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double L    = Domain_grid.getX(Nx);
    double H    = Domain_grid.getY(Ny);
    double Uref = fluid_property.getinlet_average();

    // Column closest to x = L/2
    int i_half = 0;
    { double best = 1e30;
      for (int i = 0; i < Nx; ++i) {
          double xc = Domain_grid.getX(i) + 0.5*Domain_grid.dx(i);
          if (std::fabs(xc - 0.5*L) < best) { best = std::fabs(xc - 0.5*L); i_half = i; }
      }
    }
    // Row closest to y = H/2
    int j_half = 0;
    { double best = 1e30;
      for (int j = 0; j < Ny; ++j) {
          double yc = Domain_grid.getY(j) + 0.5*Domain_grid.dy(j);
          if (std::fabs(yc - 0.5*H) < best) { best = std::fabs(yc - 0.5*H); j_half = j; }
      }
    }

    // ── u(y) at x=L/2  (CSV — compare with Ghia et al. 1982) ────────────────
    {
        std::ofstream f("cavity_u_profile.csv");

        f << "y_over_L,u_over_Uref,y_m,u_ms\n";
        f << "0.0,0.0,0.0,0.0\n";   // bottom wall (y=0, u=0)
        for (int j = 0; j < Ny; ++j) {
            double yc  = Domain_grid.getY(j) + 0.5*Domain_grid.dy(j);
            double u_c = 0.5*(field_vectors.u[j][i_half] + field_vectors.u[j][i_half+1]);
            f << yc/H << "," << u_c/Uref << "," << yc << "," << u_c << "\n";
        }
        f << "1.0,1.0," << H << "," << Uref << "\n";  // top lid (y=H, u=Uref)
        f.close();
    }

    // ── v(x) at y=H/2  (CSV) ────────────────────────────────────────────────
    {
        std::ofstream f("cavity_v_profile.csv");
        f << "x_over_L,v_over_Uref,x_m,v_ms\n";
        f << "0.0,0.0,0.0,0.0\n";   // left wall (x=0, v=0)
        for (int i = 0; i < Nx; ++i) {
            double xc  = Domain_grid.getX(i) + 0.5*Domain_grid.dx(i);
            double v_c = 0.5*(field_vectors.v[j_half][i] + field_vectors.v[j_half+1][i]);
            f << xc/L << "," << v_c/Uref << "," << xc << "," << v_c << "\n";
        }
        f << "1.0,0.0," << L << ",0.0\n";  // right wall (x=L, v=0)
        f.close();
    }

    // ── ParaView VTS — full field (P, u, v, |V|) ─────────────────────────────
    {
        std::ofstream vtk("cavity_results.vts");
        vtk << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "<StructuredGrid WholeExtent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n"
            << "<Piece Extent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n"
            << "<Points>\n"
            << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j <= Ny; ++j)
            for (int i = 0; i <= Nx; ++i)
                vtk << Domain_grid.getX(i) << " " << Domain_grid.getY(j) << " 0.0\n";
        vtk << "</DataArray>\n</Points>\n"
            << "<CellData>\n"
            << "<DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << field_vectors.P[j][i] << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"u_velocity\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << 0.5*(field_vectors.u[j][i] + field_vectors.u[j][i+1]) << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"v_velocity\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << 0.5*(field_vectors.v[j][i] + field_vectors.v[j+1][i]) << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"VelocityMagnitude\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                double uc = 0.5*(field_vectors.u[j][i] + field_vectors.u[j][i+1]);
                double vc = 0.5*(field_vectors.v[j][i] + field_vectors.v[j+1][i]);
                vtk << std::sqrt(uc*uc + vc*vc) << "\n";
            }
        vtk << "</DataArray>\n</CellData>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
        vtk.close();
    }
    std::cout << "Cavity output: cavity_u_profile.csv  cavity_v_profile.csv  cavity_results.vts\n";
}
//============================================
void Solver::open_Bc_report()
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    double L    = Domain_grid.getX(Nx);
    double H    = Domain_grid.getY(Ny);
    double Uref = fluid_property.getinlet_average();


    // ── ParaView VTS — full field (P, u, v, |V|) ─────────────────────────────
    {
        std::ofstream vtk("Open_BC_results.vts");
        vtk << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "<StructuredGrid WholeExtent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n"
            << "<Piece Extent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n"
            << "<Points>\n"
            << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j <= Ny; ++j)
            for (int i = 0; i <= Nx; ++i)
                vtk << Domain_grid.getX(i) << " " << Domain_grid.getY(j) << " 0.0\n";
        vtk << "</DataArray>\n</Points>\n"
            << "<CellData>\n"
            << "<DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << field_vectors.P[j][i] << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"u_velocity\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << 0.5*(field_vectors.u[j][i] + field_vectors.u[j][i+1]) << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"v_velocity\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                vtk << 0.5*(field_vectors.v[j][i] + field_vectors.v[j+1][i]) << "\n";
        vtk << "</DataArray>\n"
            << "<DataArray type=\"Float64\" Name=\"VelocityMagnitude\" format=\"ascii\">\n";
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                double uc = 0.5*(field_vectors.u[j][i] + field_vectors.u[j][i+1]);
                double vc = 0.5*(field_vectors.v[j][i] + field_vectors.v[j+1][i]);
                vtk << std::sqrt(uc*uc + vc*vc) << "\n";
            }
        vtk << "</DataArray>\n</CellData>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
        vtk.close();
    }
    std::cout << "Open BC output: Open_BC_results.vts\n";
}


void Solver::Monior(double flowtime, double dt)
{
    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();
    cout << "time = " << setw(6) << flowtime << "\t dt= " << dt
         << "\t u = " << field_vectors.u[Ny/10][Nx/10]
         << "\t Error= " << report_error << endl;

}

/********************* ParaView Export *************************/
void Solver::export_paraview_Cylinder(double flowtime)
{
    static int last_written_second = -1;
    int current_second = static_cast<int>(flowtime);
    if (current_second <= last_written_second) return;
    last_written_second = current_second;

    int Nx = Domain_grid.getNx();
    int Ny = Domain_grid.getNy();

    std::ostringstream fname;
    fname << "flow_" << current_second << ".vts";
    std::ofstream vtk(fname.str());
    if (!vtk.is_open()) { std::cerr << "Error opening ParaView file\n"; return; }

    vtk << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<StructuredGrid WholeExtent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n"
        << "<Piece Extent=\"0 " << Nx << " 0 " << Ny << " 0 0\">\n"
        << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int j = 0; j <= Ny; ++j)
        for (int i = 0; i <= Nx; ++i)
            vtk << Domain_grid.getX(i) << " " << Domain_grid.getY(j) << " 0.0\n";
    vtk << "</DataArray>\n</Points>\n<CellData>\n"
        << "<DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">\n";
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
            vtk << field_vectors.P[j][i] << "\n";
    vtk << "</DataArray>\n"
        << "<DataArray type=\"Float64\" Name=\"u_velocity\" format=\"ascii\">\n";
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
            vtk << 0.5*(field_vectors.u[j][i] + field_vectors.u[j][i+1]) << "\n";
    vtk << "</DataArray>\n"
        << "<DataArray type=\"Float64\" Name=\"v_velocity\" format=\"ascii\">\n";
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
            vtk << 0.5*(field_vectors.v[j][i] + field_vectors.v[j+1][i]) << "\n";
    vtk << "</DataArray>\n"
        << "<DataArray type=\"Float64\" Name=\"VelocityMagnitude\" format=\"ascii\">\n";
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i) {
            double uc = 0.5*(field_vectors.u[j][i] + field_vectors.u[j][i+1]);
            double vc = 0.5*(field_vectors.v[j][i] + field_vectors.v[j+1][i]);
            vtk << std::sqrt(uc*uc + vc*vc) << "\n";
        }
    vtk << "</DataArray>\n</CellData>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
    std::cout << "ParaView export written: " << fname.str() << std::endl;
}

/********************* Strouhal Number *************************/
double Solver::compute_strouhal()
{
    std::ofstream CLfile("CL_history.txt");
    for (size_t i = 0; i < CL_history.size(); ++i)
        CLfile << Time_history[i] << " " << CL_history[i] << std::endl;
    CLfile.close();

    int Nsig = (int)CL_history.size();
    if (Nsig < 4 || (int)Time_history.size() < 4) return 0.0;

    double t0 = Time_history.front(), t1 = Time_history.back();
    double Twindow = t1 - t0;
    if (Twindow <= 0.0) return 0.0;

    int Nuni = 1;
    while (Nuni < Nsig) Nuni <<= 1;
    double dt_uni = Twindow / (Nuni - 1);

    std::vector<double> cl_uniform(Nuni);
    int src = 0;
    for (int k = 0; k < Nuni; ++k) {
        double t = t0 + k * dt_uni;
        while (src < Nsig - 2 && Time_history[src + 1] < t) ++src;
        double t_lo = Time_history[src];
        double t_hi = (src + 1 < Nsig) ? Time_history[src + 1] : t_lo;
        double frac = (t_hi > t_lo) ? (t - t_lo) / (t_hi - t_lo) : 0.0;
        cl_uniform[k] = CL_history[src] + frac * (CL_history[min(src+1, Nsig-1)] - CL_history[src]);
    }

    double mean = std::accumulate(cl_uniform.begin(), cl_uniform.end(), 0.0) / Nuni;
    for (auto &v : cl_uniform) v -= mean;

    std::vector<double> magnitude;
    FFT(cl_uniform, magnitude);

    int peak_index = 1;
    for (int i = 2; i < (int)magnitude.size(); ++i)
        if (magnitude[i] > magnitude[peak_index]) peak_index = i;

    double fs   = 1.0 / dt_uni;
    double freq = (double)peak_index * fs / Nuni;
    double D    = Domain_grid.getD();
    double U    = fluid_property.getinlet_average();
    double St   = freq * D / U;

    if (St < 0.05 || St > 0.5) {
        const double T_min_gap = 3.0;
        std::vector<double> cross_times;
        double t_prev = -1e9;
        for (int k = 1; k < Nsig; ++k)
            if (CL_history[k-1] < 0.0 && CL_history[k] >= 0.0) {
                double tc = Time_history[k];
                if (tc - t_prev >= T_min_gap) { cross_times.push_back(tc); t_prev = tc; }
            }
        int Nc = (int)cross_times.size();
        if (Nc < 2) return 0.0;
        int use = min(4, Nc - 1);
        double Tperiod = (cross_times[Nc-1] - cross_times[Nc-1-use]) / use;
        St = D / (U * Tperiod);
    }

    return St;
}

/********************* RMS *************************/
double Solver::RMS(const std::vector<double> &signal)
{
    double Mean = std::accumulate(signal.begin(), signal.end(), 0.0) / signal.size();
    double sum_squares = 0.0;
    for (const auto &val : signal) sum_squares += (val - Mean) * (val - Mean);
    return sqrt(sum_squares / signal.size());
}

/********************* FFT *************************/
void Solver::FFT(const std::vector<double> &signal, std::vector<double> &magnitude)
{
    int N = signal.size(), Nfft = 1;
    while (Nfft < N) Nfft <<= 1;

    std::vector<std::complex<double>> data(Nfft);
    for (int i = 0; i < N;    ++i) data[i] = signal[i];
    for (int i = N; i < Nfft; ++i) data[i] = 0.0;

    std::function<void(std::vector<std::complex<double>>&)> fft =
        [&](std::vector<std::complex<double>>& x) {
            int n = x.size();
            if (n <= 1) return;
            std::vector<std::complex<double>> even(n/2), odd(n/2);
            for (int i = 0; i < n/2; ++i) { even[i] = x[2*i]; odd[i] = x[2*i+1]; }
            fft(even); fft(odd);
            for (int k = 0; k < n/2; ++k) {
                std::complex<double> t = std::polar(1.0, -2*M_PI*k/n) * odd[k];
                x[k] = even[k] + t;  x[k + n/2] = even[k] - t;
            }
        };

    fft(data);
    magnitude.resize(Nfft / 2);
    for (int i = 0; i < Nfft/2; ++i) magnitude[i] = std::abs(data[i]);
}

/********************* pressure_face *************************/
double Solver::pressure_face(int i, int j)
{
    return field_vectors.P[j][i];
}
