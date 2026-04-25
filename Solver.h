#ifndef SOLVER_H
#define SOLVER_H

#include "Grid.h"
#include "Property.h"
#include "Cylinder.h"
#include "Vectors.h"
#include "BoundaryConditions.h"
#include <vector>

class Solver {

private:
    Grid&     Domain_grid;
    Property& fluid_property;
    Vectors&  field_vectors;
    Cylinder& Square_cylinder;

    BoundaryConditions BC;

public:
    double flowtime    = 0.0, initialdt;
    int    MaxIter     = 500;
    double Max_error   = 1e-7;
    double dtsum       = 0;
    double total_error = 0.0;
    double report_error = 0;
    double omega       = 0.7;
    double Cdf = 0.0, Cdp = 0.0, Cd = 0.0;
    double Cl  = 0.0;
    double CFD_conv = 0.15;
    double CFL_diff = 0.15;
    double Time_period = 0.0;
    double Tavg_stored = 0.0;
    int    step_count  = 0;
    int    TP = 0;
    bool   has_cylinder = false;
    bool   open_BC = false;
    std::vector<double> CL_history;
    std::vector<double> Time_history;

    int  pressure_solver = 0;   // 0 = Gauss-Seidel,  1 = Conjugate Gradient
    bool AT = false;

    // FIX: removed Ultime parameter — constructor in .cpp does not take it
    Solver(Grid& grid, Property& fluidprop, Vectors& vec_vec, Cylinder& cy,
           bool AT_flag, const BoundaryConditions& bc, bool cylinder_present,bool OBC);

    double getinitialdt() const;

    void   Boundary_condition       (double dtt);
    void   apply_left_bc           (double dtt);
    void   apply_right_bc          (double dtt);
    void   apply_bottom_bc          (double dtt);
    void   apply_top_bc             (double dtt);
    void   apply_pressure_domain_bc ();

    void   Interrior_Velocity  (double timestep);
    void   pressure_boundary   ();
    void   Guessiedel          (double dtt);
    void   ConjugateGradient   (double dtt);
    void   correction_velocity (double dtt);
    void   solve               (double ultimate_time);

    double Drag_pressure    ();
    double Drag_friction    ();
    double Lift_coefficient ();
    void   FFT(const std::vector<double>& signal, std::vector<double>& magnitude);
    double compute_strouhal ();
    double RMS(const std::vector<double>& signal);
    void   cylinder_reports  ();
    void   cavity_report    ();
    void open_Bc_report();
    double pressure_face    (int i, int j);
    void   Monior           (double t, double dt);
    void   export_paraview_Cylinder  (double flowtime);
};

#endif
