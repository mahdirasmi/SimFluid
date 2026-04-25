#ifndef CYLINDER_H
#define CYLINDER_H
#include "Grid.h"
#include "Property.h"
#include "Vectors.h"
#include <vector>

class Cylinder {
private:
    Grid&    grid;
    Vectors& field_vectors;
public:
    double xsc, ysc, D;
    double CYL_X0, CYL_X1, CYL_Y0, CYL_Y1;   // class members — assigned in constructor body
    static constexpr double EPS = 1e-5;

    Cylinder(double xc, double yc, Grid& g, Vectors& v, double D_in);
    bool is_inside (double xx, double yy, int j, int i) const;
    bool is_inside2(double xx, double yy, int j, int i) const;
    bool is_insideU(double xx, double yy, int j, int i) const;
    bool is_insideV(double xx, double yy, int j, int i) const;
    double getXs() const { return xsc; }
    double getYs() const { return ysc; }
    void get_faces(double xx, double yy,
                   bool& lf, bool& rf, bool& bf, bool& tf,
                   int j, int i) const;
};

#endif
