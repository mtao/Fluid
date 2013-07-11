#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <QtGui>

#include "pcgsolver/pcg_solver.h"

#include "mac/macgridfactory.h"

typedef Eigen::Vector2f Vec2f;

#define GRID_SIZE 50
#define GRID_WIDTH 1

#ifndef M_PI
const double M_PI = 3.1415926535897932384626433832795;
#endif

class Simulation
{
public:
    Simulation();
    void initialization(int grid_size, float grid_width);
    void step(float dt);
    void computeVelocity(float x_old, float y_old, float x, float y);

    std::vector<Vec2f> particles;
    float particle_radius;
    int n;
    float dx;
    Vec2f get_velocity(Vec2f x);

protected:
    void addExternalForces(float dt);
    void advect(float dt);
    void project(float dt);
    void advectParticles(float dt);
    Vec2f trace_rk2(Vec2f x, float dt);
    float linear_interpolation(Vec2f x, int flag);
    void linearSolver(int b, float a, float c);
    void setBoundary(int c);



private:
    // velocity field

    float u[GRID_SIZE+1][GRID_SIZE];
    float u_temp[GRID_SIZE+1][GRID_SIZE];
    float v[GRID_SIZE][GRID_SIZE+1];
    float v_temp[GRID_SIZE][GRID_SIZE+1];

/*
    Eigen::VectorXd pressure;
    Eigen::VectorXd rhs;
    Eigen::SparseMatrix<double> matrix;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
    std::vector<Eigen::Triplet<double> > triplets;
*/

    //Solver data
    PCGSolver<double> solver;
    SparseMatrixd matrix;
    std::vector<double> rhs;
    std::vector<double> pressure;

};

#endif // SIMULATION_H
