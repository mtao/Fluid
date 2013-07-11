#include "simulation.h"

Simulation::Simulation()
{

}

void Simulation::initialization(int grid_size, float grid_width){

    dx = (float) grid_width / grid_size;
    particle_radius = dx/sqrt(2.0f);
    n = grid_size;
/*
    Eigen::AlignedBox<double,2> m(Eigen::Vector2d(0,0),Eigen::Vector2d(grid_width,grid_width));
    MACGridFactory<double,2> factory(Eigen::Vector2i(n,n),m);

    auto&& temp = factory.create<mtao::UGrid>();
*/

    //The grid will be fined by an axis aligned bounding box...
    Eigen::AlignedBox<double,3> m(Eigen::Vector3d(0,0,0),Eigen::Vector3d(1,1,1));
    //In order to have all grids created be consistent wit heach other you need to make them
    //from a factory
    MACGridFactory<double,3> factory(Eigen::Vector3i(10,10,10),m);

    //After you have a factory you can create a grid rather easily!
    auto&& g = factory.create<mtao::UGrid>();

    /*
    u = factory.create<mtao::UGrid>();
    u_old = factory.create<mtao::UGrid>();
    v = factory.create<mtao::VGrid>();
    v_old = factory.create<mtao::VGrid>();

    std::fill(u.begin(), u.end(), 0);
    std::fill(u_old.begin(), u_old.end(), 0);
    std::fill(v.begin(), v.end(), 0);
    std::fill(v_old.begin(), v_old.end(), 0);
*/

    // u
    for (int i=0; i<=n; i++){
        for (int j=0; j<n; j++){
            u[i][j] = 0;
            u_temp[i][j] = 0;
        }
    }

    // v
    for (int i=0; i<n; i++){
        for (int j=0; j<=n; j++){
            v[i][j] = 0;
            v_temp[i][j] = 0;
        }
    }

    for (int i=0; i<n*n*5; i++){
        float x = float(rand() % 1000) / 1000 * n * dx;
        float y = float(rand() % 1000) / 1000 * n * dx;

        if(x < 0.5)
            particles.push_back(Vec2f(x,y));
    }

}

void Simulation::step(float dt){
    advectParticles(dt);

    addExternalForces(dt);
    advect(dt);
    project(dt);
}

void Simulation::addExternalForces(float dt){

    for (int i=0; i<n; i++){
        for (int j=0; j<=n; j++){
            //v[i][j] -= 0.2;
            /*if(i==20 && j==20)
                u[i][j] = 100;*/
        }
    }

}


void Simulation::advect(float dt){

    Vec2f x(0,0);

    // u
    for (int i=0; i<=n; i++){
        for (int j=0; j<n; j++){
            // position, origin is (0,0)
            Vec2f x(i*dx, (j + 0.5)*dx);

            // trace particle
            // second order Runge-Kutta method
            Vec2f p = trace_rk2(x, -dt);

            // linear interpolation
            u_temp[i][j] = get_velocity(p)[0];
        }
    }
    // v
    for (int i=0; i<n; i++){
        for (int j=0; j<=n; j++){
            // position, origin is (0,0)
            Vec2f x((i + 0.5)*dx, j*dx);

            // trace particle
            // second order Runge-Kutta method
            Vec2f p = trace_rk2(x, -dt);

            // linear interpolation
            v_temp[i][j] = get_velocity(p)[1];
        }
    }

    // u
    for (int i=0; i<=n; i++){
        for (int j=0; j<n; j++){
            u[i][j] = u_temp[i][j];
        }
    }

    // v
    for (int i=0; i<n; i++){
        for (int j=0; j<=n; j++){
            v[i][j] = v_temp[i][j];
        }
    }
}

Vec2f Simulation::trace_rk2(Vec2f x, float dt){
    Vec2f velocity = get_velocity(x);
    velocity = get_velocity(x + (0.5 * dt * velocity));
    return x + dt * velocity;
}

Vec2f Simulation::get_velocity(Vec2f x){

    float u_val = linear_interpolation(x / dx - Vec2f(0, 0.5), 0);
    float v_val = linear_interpolation(x / dx - Vec2f(0.5, 0), 1);

    return Vec2f(u_val, v_val);
}

float Simulation::linear_interpolation(Vec2f t, int flag){
    int i, j;
    // find the four points
    // round down
    i = floor(t[0]);
    j = floor(t[1]);

    // cutoff at boundary
    if(i < 0){
        i = 0;
        t[0] = i;
    }else if(i > n - 2){
        i = n - 2;
        t[0] = i+1;
    }
    if(j < 0){
        j = 0;
        t[1] = j;
    }else if(j > n - 2){
        j = n - 2;
        t[1] = j+1;
    }
    float r1, r2;
    // bilinear interpolation
    if(flag == 0){
        r1 = (i + 1 - t[0]) * u[i][j] + (t[0] - i) * u[i+1][j];
        r2 = (i + 1 - t[0]) * u[i][j+1] + (t[0] - i) * u[i+1][j+1];
    }else{
        r1 = (i + 1 - t[0]) * v[i][j] + (t[0] - i) * v[i+1][j];
        r2 = (i + 1 - t[0]) * v[i][j+1] + (t[0] - i) * v[i+1][j+1];
    }

    return (j + 1 - t[1])*r1 + (t[1] - j)*r2;
}

void Simulation::project(float dt){
    int size = n * n;

/*
    matrix.resize(n, n);
    matrix.setZero();
    rhs.resize(n);
    pressure.resize(n);

    float term = dx * dx;

    // setup matrix and rhs
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

            int index = i + n * j;
            int count = 0;
            rhs[index] = 0;
            pressure[index] = 0;
            //right
            if(i < n - 1){
                matrix.insert(index, index + 1) = -1;
                triplets.push_back(Eigen::Triplet<double>(index,index+1,-1));
                rhs[index] -= u[i+1][j] * dx;
                count++;
            }

            //left
            if(i > 0){
                matrix.insert(index, index - 1) = -1;
                triplets.push_back(Eigen::Triplet<double>(index,index-1, -1));
                rhs[index] += u[i][j] * dx;
                count++;
            }

            //top
            if(j < n - 1){
                matrix.insert(index, index + n) = -1;
                triplets.push_back(Eigen::Triplet<double>(index,index+n, -1));
                rhs[index] -= v[i][j+1] * dx;
                count++;
            }

            //bottom
            if(j > 0){
                matrix.insert(index, index - n) = -1;
                triplets.push_back(Eigen::Triplet<double>(index,index-n, -1));
                rhs[index] += v[i][j] * dx;
                count++;
            }

            matrix.insert(index, index) = count;
            triplets.push_back(Eigen::Triplet<double>(index,index, count));
        }
    }
*/
    /*
    // checking sum of matrix
    for (int k=0; k<matrix.outerSize(); ++k){
        float count = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix,k); it; ++it)
        {
            count += it.value();
        }
        if(count != 0)
            qDebug("problem");
    }

    // check symmetry
    for (int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(matrix.coeff(i, j) != matrix.coeff(j, i))
                qDebug("problem");
        }
    }

    */
    /*
    qDebug("done");

    solver.compute(matrix);
    pressure = solver.solve(rhs);


    Eigen::SparseMatrix<double> mat(n,n);
    mat.setFromTriplets(triplets.begin(),triplets.end());
    Eigen::Map<Eigen::VectorXd> e_pressure(pressure.data(),pressure.size());
    Eigen::Map<Eigen::VectorXd> e_rhs(rhs.data(),rhs.size());
    qDebug("%f", (mat*e_pressure - e_rhs).norm());

*/
    if(rhs.size() != size) {
       rhs.resize(size);
       pressure.resize(size);
       matrix.resize(size);
    }
    matrix.zero();
    for(int j=0; j<n; j++){
        for(int i=0; i<n; i++){

            int index = i + n * j;
            int count = 0;
            rhs[index] = 0;
            pressure[index] = 0;
            //right
            if(i < n - 1){
                matrix.add_to_element(index, index + 1, -1);
                rhs[index] -= u[i+1][j] * dx;
                count++;
            }

            //left
            if(i > 0){
                matrix.add_to_element(index, index - 1, -1);
                rhs[index] += u[i][j] * dx;
                count++;
            }

            //top
            if(j < n - 1){
                matrix.add_to_element(index, index + n, -1);
                rhs[index] -= v[i][j+1] * dx;
                count++;
            }

            //bottom
            if(j > 0){
                matrix.add_to_element(index, index - n, -1);
                rhs[index] += v[i][j] * dx;
                count++;
            }

            matrix.add_to_element(index, index, count);
       }
    }

    //Solve the system using Robert Bridson's incomplete Cholesky PCG solver

    double tolerance;
    int iterations;
    bool success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
    if(!success) {
        qDebug("WARNING: Pressure solve failed!************************************************\n");
    }
    //qDebug("%f", tolerance);

    // u
    for (int i=1; i<=n; i++){
        for (int j=0; j<n; j++){
            int index = i + n * j;
            u[i][j] -= (pressure[index] - pressure[index - 1]) / dx;
        }
    }
    for (int j=0; j<n; j++)
        u[0][j] = u[1][j];


    // v
    for (int i=0; i<n; i++){
        for (int j=1; j<n; j++){
            int index = i + n * j;
            v[i][j] -= (pressure[index] - pressure[index - n]) / dx;
        }
    }
    for (int i=0; i<n; i++)
        v[i][0] = v[i][1];
}

void Simulation::advectParticles(float dt){
    for(int i=0; i < particles.size(); i++){
        Vec2f x = particles[i];

        particles[i] = trace_rk2(x, dt);

        if( particles[i][0] < 0 || particles[i][0] > n * dx
                || particles[i][1] < 0 || particles[i][1] > n * dx)
            particles[i] = x;


    }
}


void Simulation::computeVelocity(float x_old, float y_old, float x, float y){

    qDebug("%f, %f, %f, %f", x_old, y_old, x, y);
    Vec2f start(x_old, y_old);
    Vec2f end(x, y);
    Vec2f dir = end - start;

    if(dir.norm() < 1e-14)
       return;

    int i = floor(x_old / dx);
    int j = floor(y_old / dx);

    u[i][j] += dir[0]*100;
    v[i][j] += dir[1]*100;
    qDebug("%f, %f", dir[0], dir[1]);
}
