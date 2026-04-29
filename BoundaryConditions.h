#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H


enum class BCType {
    DIRICHLET, NEUMANN, CONVECTIVE, NO_SLIP, SYMMETRY
};

inline const char* bc_name(BCType t)
{
    switch (t) {
        case BCType::DIRICHLET:  return "DIRICHLET";
        case BCType::NEUMANN:    return "NEUMANN";
        case BCType::CONVECTIVE: return "CONVECTIVE";
        case BCType::NO_SLIP:    return "NO_SLIP";
        case BCType::SYMMETRY:   return "SYMMETRY";
        default:                 return "UNKNOWN";
    }
}


enum class InletType {
    UNIFORM,    // u = uniform
    PARABOLIC   // u = 2 * Uavg * (1 - ((y - y_mid) / (H/2))^2)
               
};

inline const char* inlet_name(InletType t)
{
    switch (t) {
        case InletType::UNIFORM:   return "UNIFORM";
        case InletType::PARABOLIC: return "PARABOLIC (Poiseuille)";
        default:                   return "UNKNOWN";
    }
}

struct WallBC {
    BCType    type       = BCType::SYMMETRY;
    InletType inlet_type = InletType::UNIFORM;  // only used when type == DIRICHLET
    double    value_u    = 0.0;
    double    value_v    = 0.0;
    double    value_p    = 0.0;
};


struct BoundaryConditions{
    WallBC left;    // left  wall  (x = 0)
    WallBC right;   // right wall  (x = L)
    WallBC bottom;   // bottom wall (y = 0)
    WallBC top;      // top wall    (y = H)
};

#endif
