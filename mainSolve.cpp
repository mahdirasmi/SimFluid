#include "Grid.h"
#include "Property.h"
#include "Cylinder.h"
#include "Vectors.h"
#include "Solver.h"
#include "BoundaryConditions.h"
#include <iostream>
#include <stdexcept>


// A function to interactively ask user for BC type and values for each boundary
static WallBC ask_bc(const std::string& wall_name)
{
    std::cout << "\n  " << wall_name << ":\n"
              << "    1 - Dirichlet   (fixed velocity)\n"
              << "    2 - Neumann     (zero gradient)\n"
              << "    3 - Convective  (outlet)\n"
              << "    4 - No-slip     (solid wall, u=v=0)\n"
              << "    5 - Symmetry    (du/dn=0, v_normal=0)\n"
              << "  Choice: ";

    int c = 0;
    while (!(std::cin >> c) || c < 1 || c > 5) {
        std::cin.clear();
        std::cout << "  Enter 1-5: ";
    }

    WallBC w;
    switch (c) {
        case 1:
            w.type = BCType::DIRICHLET;
            std::cout << "    u value: "; std::cin >> w.value_u;
            std::cout << "    v value: "; std::cin >> w.value_v;
            break;
        case 2: w.type = BCType::NEUMANN;    break;
        case 3: w.type = BCType::CONVECTIVE; break;
        case 4: w.type = BCType::NO_SLIP;    break;
        case 5: w.type = BCType::SYMMETRY;   break;
    }

    std::cout << "  -> " << wall_name << " = " << bc_name(w.type) << "\n";
    return w;
}

// ────────────────────────────: Main Block ─Starts :─────────────────────────────────────────────────//
int main()
{

 std::cout << R"(

  ███████╗██╗███╗   ███╗███████╗██╗     ██╗   ██╗██╗██████╗ 
  ██╔════╝██║████╗ ████║██╔════╝██║     ██║   ██║██║██╔══██╗
  ███████╗██║██╔████╔██║█████╗  ██║     ██║   ██║██║██║  ██║
  ╚════██║██║██║╚██╔╝██║██╔══╝  ██║     ██║   ██║██║██║  ██║
  ███████║██║██║ ╚═╝ ██║██║     ███████╗╚██████╔╝██║██████╔╝
  ╚══════╝╚═╝╚═╝     ╚═╝╚═╝     ╚══════╝ ╚═════╝ ╚═╝╚═════╝ 

)" << std::endl;

    std::cout << "\nDeveloped by Mahdi.Rs\n"<<std::endl;

    // ================================== 1. Domain=======================================================
    double L = 0, H = 0;
    std::cout << "Domain length L: "; std::cin >> L;
    std::cout << "Domain height H: "; std::cin >> H;

    // ================================== 2. Reynolds number=======================================================
    double Re = 100.0;
    std::cout << "Reynolds number Re: "; std::cin >> Re;
    std::cout << "Density is set to 1.0 and viscosity is calculated accordingly.\n";

    // ================================== 3. Square cylinder=======================================================
    bool has_cylinder = false;
    std::cout << "Does the domain contain a square cylinder? (1=yes / 0=no): "; std::cin >> has_cylinder;

    double Xc = 0, Yc = 0, D = 1;
    if (has_cylinder) {
        std::cout << "  Characteristic length D: "; std::cin >> D;
        std::cout << "  Center Xc:    "; std::cin >> Xc;
        std::cout << "  Center Yc:    "; std::cin >> Yc;
    }

    // ================================== 4. Simulation time=======================================================
    double Ultime = 0;
    bool   AT     = false;
    std::cout << "Ultimate time: "; std::cin >> Ultime;
    std::cout << "Do you want Time averaging? (1=yes / 0=no): "; std::cin >> AT;

    // ================================== 5. Boundary conditions=======================================================
    std::cout << "\n--- Boundary Conditions ---\n";
    if (has_cylinder)
        std::cout << "  (Cylinder surface is automatically set to No-slip)\n";

    BoundaryConditions BC;
    // Ask user for BC type and values for each wall
    BC.left   = ask_bc("(left,   x=0)");
    BC.right  = ask_bc("(right,  x=L)");
    BC.bottom = ask_bc("(Bottom y=0)");
    BC.top    = ask_bc("(top    y=H)");


    if (BC.left.type == BCType::DIRICHLET) {
        std::cout << "\n  Inlet velocity profile:\n"
                  << "    1 - Uniform   (flat plug flow)\n"
                  << "    2 - Parabolic (Poiseuille — u_max = 1.5 * Uavg at centreline)\n"
                  << "  Choice: ";
        int pc = 1;
        while (!(std::cin >> pc) || pc < 1 || pc > 2) {
            std::cin.clear();
            std::cout << "  Enter 1 or 2: ";
        }
        BC.left.inlet_type = (pc == 2) ? InletType::PARABOLIC : InletType::UNIFORM;
        std::cout << "  -> Inlet profile = " << inlet_name(BC.left.inlet_type) << "\n";
    }

    bool open_BC =
        BC.left.type   == BCType::CONVECTIVE || BC.left.type   == BCType::NEUMANN ||
        BC.right.type  == BCType::CONVECTIVE || BC.right.type  == BCType::NEUMANN ||
        BC.top.type    == BCType::CONVECTIVE || BC.top.type    == BCType::NEUMANN ||
        BC.bottom.type == BCType::CONVECTIVE || BC.bottom.type == BCType::NEUMANN;
 
    // ================================== 6. Mesh=======================================================
    std::cout << "\n--- Mesh ---\n"
              << "  1 - Uniform\n"
              << "  2 - Non-uniform  (square cylinder, 3-block stretching)\n"
              << "  3 - Non-uniform  (symmetric wall refinement)\n"
              << "  Choice: ";
    int mesh_choice = 0;
    std::cin >> mesh_choice;

    // ================================== 7. Build grid=======================================================
    Grid* grid_ptr = nullptr;

    if (mesh_choice == 1) {
        int Nx_uni = 0, Ny_uni = 0;
        std::cout << "  Cells in x (Nx): "; std::cin >> Nx_uni;
        std::cout << "  Cells in y (Ny): "; std::cin >> Ny_uni;
        grid_ptr = new Grid(Nx_uni, Ny_uni, L, H);
    }
    else if (mesh_choice == 2) {
        double FCC = 0;
        int lx1 = 0, lx2 = 0, ly1 = 0, ly2 = 0;
        std::cout << "  First cell size (FCC): ";     std::cin >> FCC;
        std::cout << "  Nodes - left block  (lx1): "; std::cin >> lx1;
        std::cout << "  Nodes - right block (lx2): "; std::cin >> lx2;
        std::cout << "  Nodes - bottom block (ly1): "; std::cin >> ly1;
        std::cout << "  Nodes - top block   (ly2): "; std::cin >> ly2;
        grid_ptr = new Grid(lx1, lx2, ly1, ly2, L, H, FCC, Xc, Yc, D);
    }
    else {
        int Nx_sym = 0, Ny_sym = 0;
        double r = 1.05;
        std::cout << "  Cells in x (Nx): ";           std::cin >> Nx_sym;
        std::cout << "  Cells in y (Ny): ";           std::cin >> Ny_sym;
        std::cout << "  Growth rate (e.g. 1.05): ";   std::cin >> r;
        grid_ptr = new Grid(Nx_sym, Ny_sym, L, H, r);
    }

    Grid& grid = *grid_ptr;
    grid.mesh_generation();

    // ================================== 8. Objects======================================================= 
    Property property(1.0, Re, grid);
    Vectors  vectors(grid, property);

    // FIX: place cylinder far outside domain when not present so is_inside always returns false
    double cxc = has_cylinder ? Xc : -1e6;
    double cyc = has_cylinder ? Yc : -1e6;
    double cD  = has_cylinder ? D  : 1.0;
    Cylinder cylinder(cxc, cyc, grid, vectors, cD);

    vectors.initilization_u();

    // Build the inlet velocity profile now that the grid is constructed.
    // For open-flow cases (DIRICHLET left wall) the profile (uniform or
    // parabolic) is stored in Property and applied every time step by
    // apply_left_bc(). For all other cases (cavity, channel with symmetry
    // etc.) the default uniform profile is harmless.
    property.build_inlet_profile(BC.left.inlet_type);

    vectors.u_old = vectors.u;  vectors.uupp   = vectors.u;  vectors.u_next = vectors.u;
    vectors.v_old = vectors.v;  vectors.vvpp   = vectors.v;  vectors.v_next = vectors.v;

    // ================================== 9. Pressure solver=======================================================
    std::cout << "\nPressure solver  (0=Gauss-Seidel / 1=Conjugate Gradient): ";
    int ps = 0; std::cin >> ps;

    // ================================== 10. Summary=======================================================
    std::cout << "\n--- Summary ---\n"
              << "  Left   : " << bc_name(BC.left.type);
    if (BC.left.type == BCType::DIRICHLET)
        std::cout << "  [" << inlet_name(BC.left.inlet_type) << "]";
    std::cout << "\n"
              << "  Right  : " << bc_name(BC.right.type)  << "\n"
              << "  Bottom : " << bc_name(BC.bottom.type) << "\n"
              << "  Top    : " << bc_name(BC.top.type)    << "\n";
    if (has_cylinder)
        std::cout << "  Cylinder wall: No-slip (automatic)\n";

    // ================================== 11. Run=======================================================
    // FIX: removed Ultime from constructor call — Solver no longer takes it
    Solver solver(grid, property, vectors, cylinder, AT, BC, has_cylinder,open_BC);
    solver.pressure_solver = ps;

    solver.solve(Ultime);



    if (has_cylinder) {
    solver.cylinder_reports();
        double st  = solver.compute_strouhal();
        double rms = solver.RMS(solver.CL_history);
        std::cout << "Strouhal number:       " << st  << "\n";
        std::cout << "RMS lift coefficient:  " << rms << "\n";
    }

    if (!has_cylinder && !open_BC)
        solver.cavity_report();
    if (!has_cylinder )
        solver.open_Bc_report();

    std::cout << "\nSimulation complete.\n";
    delete grid_ptr;
    return 0;
}
