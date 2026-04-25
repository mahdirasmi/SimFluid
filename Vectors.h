#ifndef VECTOR_H
#define VECTOR_H

#include "Grid.h"
#include "Property.h"
#include "Cylinder.h"
#include <vector>
class Vectors {
public:
    Grid& grid;
    Property& fluid;


    int Nx, Ny; // Number of grid points in x and y directions
    std::vector<std::vector<double>> u, v, P, uupp, vvpp, u_next,v_next,P_next,u_old,v_old,sumu;

    Vectors(Grid& g,Property& f);
    void initilization_u();
};
#endif