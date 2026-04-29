#include "Grid.h"
#include "Property.h"
#include "Cylinder.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace std;

Cylinder::Cylinder(double xc, double yc, Grid& g, Vectors& v, double D_in)
    : xsc(xc), ysc(yc), grid(g), D(D_in), field_vectors(v)
{
    
    CYL_X0 = xsc - D / 2;
    CYL_X1 = xsc + D / 2;
    CYL_Y0 = ysc - D / 2;
    CYL_Y1 = ysc + D / 2;
}

// -----------------------------------------------------------------------
// is_inside — pressure cell masking
// -----------------------------------------------------------------------
bool Cylinder::is_inside(double /*xx*/, double /*yy*/, int j, int i) const
{
    double xp = grid.getX(i) + 0.5 * grid.dx(i);
    double yp = grid.getY(j) + 0.5 * grid.dy(j);
    return (xp > CYL_X0 + EPS && xp < CYL_X1 - EPS &&
            yp > CYL_Y0 + EPS && yp < CYL_Y1 - EPS);
}

// -----------------------------------------------------------------------
// is_insideU — u-face masking
// -----------------------------------------------------------------------
bool Cylinder::is_insideU(double /*xx*/, double /*yy*/, int j, int i) const
{
    double xu     = grid.getX(i);
    double yc_row = 0.5 * (grid.getY(j) + grid.getY(j + 1));
    return (xu     >= CYL_X0 - EPS && xu     <= CYL_X1 + EPS &&
            yc_row >= CYL_Y0 - EPS && yc_row <= CYL_Y1 + EPS);
}

// -----------------------------------------------------------------------
// is_insideV — v-face masking
// -----------------------------------------------------------------------
bool Cylinder::is_insideV(double /*xx*/, double /*yy*/, int j, int i) const
{
    double xc_col = 0.5 * (grid.getX(i) + grid.getX(i + 1));
    double yv     = grid.getY(j);
    return (xc_col >= CYL_X0 - EPS && xc_col <= CYL_X1 + EPS &&
            yv     >= CYL_Y0 - EPS && yv     <= CYL_Y1 + EPS);
}

// -----------------------------------------------------------------------
// get_faces for postprocessing
// -----------------------------------------------------------------------
void Cylinder::get_faces(double /*xx*/, double /*yy*/,bool& lf, bool& rf, bool& bf, bool& tf,int j, int i) const
{
    const double eps = 1e-4;
    lf = rf = bf = tf = false;

    if (fabs(grid.getX(i + 1) - CYL_X0) < eps &&
        grid.getY(j)   >= CYL_Y0 - eps &&
        grid.getY(j+1) <= CYL_Y1 + eps)
        lf = true;

    if (fabs(grid.getX(i) - CYL_X1) < eps &&
        grid.getY(j)   >= CYL_Y0 - eps &&
        grid.getY(j+1) <= CYL_Y1 + eps)
        rf = true;

    if (fabs(grid.getY(j+1) - CYL_Y0) < eps &&
        grid.getX(i)   >= CYL_X0 - eps &&
        grid.getX(i+1) <= CYL_X1 + eps)
        bf = true;

    if (fabs(grid.getY(j) - CYL_Y1) < eps &&
        grid.getX(i)   >= CYL_X0 - eps &&
        grid.getX(i+1) <= CYL_X1 + eps)
        tf = true;
}

bool Cylinder::is_inside2(double x, double y, int /*j*/, int /*i*/) const
{
    return (x >= xsc - D/2.0 - 0.001 && x <= xsc + D/2.0 &&
            y >= ysc - D/2.0 - 0.001 && y <= ysc + D/2.0);
}
