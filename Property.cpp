#include "Grid.h"
#include "Property.h"
#include "BoundaryConditions.h"
#include <iostream>
#include <cmath>

using namespace std;

Property::Property(double density, double Re_flow, Grid& g)
    : rho(density), Re(Re_flow), grid(g), inletVel(1.0),
      viscosity(density * 1.0 * 1.0 / Re_flow),
      inlet_profile(g.getNy(), 1.0),  // default: uniform u=1
      profile_type(InletType::UNIFORM)
{}

void Property::build_inlet_profile(InletType type)
{
    profile_type = type;
    int    Ny = grid.getNy();
    double H  = grid.getY(Ny);   // domain height

    if (type == InletType::UNIFORM) {
        // Every row gets the uniform average velocity
        for (int j = 0; j < Ny; ++j)
            inlet_profile[j] = inletVel;

    } else {
      
        double y_mid = 0.5 * H;
        double R     = 0.5 * H;   // half-channel height 
        for (int j = 0; j < Ny; ++j) {
            double yc = grid.getY(j) + 0.5 * grid.dy(j);  // cell centre
            double eta = (yc - y_mid) / R;                 // -1 at walls, 0 at centre
            inlet_profile[j] = 2 * inletVel * (1.0 - eta * eta);
        }
    }

    cout << "[Property] Inlet profile: " << inlet_name(type)
         << "  (Uavg = " << inletVel << ")\n";
}

