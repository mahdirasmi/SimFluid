#ifndef PROPERTY_H
#define PROPERTY_H

#include "Grid.h"
#include "BoundaryConditions.h"
#include <vector>

class Property {
private:
    double rho, viscosity, Re;
    double inletVel;
    Grid&  grid;


    std::vector<double> inlet_profile;
    InletType profile_type = InletType::UNIFORM;

public:
    Property(double density, double Re_flow, Grid& g);

 
    void build_inlet_profile(InletType type);

    double get_inlet_u(int j) const { return inlet_profile[j]; }

    double getinlet_average() const { return inletVel;  }
    double getrho()           const { return rho; }
    double getviscosity()     const { return viscosity; }
    InletType getInletType()  const { return profile_type;  }
};

#endif