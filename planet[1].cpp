//Planet.cpp
#include "planets.h"

using namespace std;

//Constructor with known values
Planet::Planet(const double& x_, const double& y_, const double& xdot_, const double& ydot_, const double& m_) { 
    x = x_;
    y = y_;
    xdot = xdot_;
    ydot = ydot_;
    m = m_; 
}

//Constructor with no values
Planet::Planet() {
    x = 0;
    y = 0;
    xdot = 0;
    ydot = 0;
    m = 0;
}

//Destructor
Planet::~Planet() {}

//Getter definitions
double Planet::GetX() {
    return x;
}
double Planet::GetY() {
    return y;
}
double Planet::GetXdot() {
    return xdot;
}
double Planet::GetYdot() {
    return ydot;
}
double Planet::GetM() {
    return m;
}

//Setter definitions
void Planet::SetX(double x_) {
    x = x_;
}
void Planet::SetY(double y_) {
    y = y_;
}
void Planet::SetXdot(double xdot_) {
    xdot = xdot_;
}
void Planet::SetYdot(double ydot_) {
    ydot = ydot_;
}
void Planet::SetM(double m_) {
    m = m_;
}

//Operator overload of = to make all variables of a planet equal to another plant's variables.
Planet& Planet::operator=(const Planet& planet2) {
    x = planet2.x;
    y = planet2.y;
    xdot = planet2.xdot;
    ydot = planet2.ydot;
    m = planet2.m;
    return *this;
}

//Operator overload of += to add a vector of positions and velocities to a planet's position and velocities.
Planet& Planet::operator+=(vector<double> a) {
    x += a[0];
    y += a[1];
    xdot += a[2];
    ydot += a[3];
    return *this;
}

//constant find function, finds the specific rk constant by evaluating the function derivative wrt to all planets, and then multiplying by step size.
vector<vector<double>> Planet::constantfind(vector<Planet>& planetstore,const double& no_lines ,const double& G, const double& h) {

    vector<vector<double>> derivatives(no_lines-1, vector<double>(4));

    vector<vector<double>>  k(no_lines - 1, vector<double>(4));

    for (int z = 0; z < no_lines - 1; ++z) {
        k[z] = { 0.0,0.0,0.0,0.0 };
        for (int j = 0; j < no_lines - 1; ++j) {
            if (j != z) {
                double dx = planetstore[j].x - planetstore[z].x;
                double dy = planetstore[j].y - planetstore[z].y;
                double d12 = sqrt(dx * dx + dy * dy);

                //lambda is the multiplier that represents gravitational pull  
                double lambda = G * planetstore[j].m / (pow(d12, 3));
                derivatives[z] = { planetstore[z].xdot, planetstore[z].ydot, lambda * dx, lambda * dy };
                vector<double> tempk = {0.0, 0.0, k[z][2], k[z][3]}; // to add on, only to velocities, from rk logic
                k[z] = tempk + derivatives[z] * h;
            }
        }
    }
    return k;
}
//update pos function, evaluates the specific yn at which the function needs to be calculated, takes in a condition, so that it works with all the constants
vector<Planet> Planet::update_pos(vector<Planet>& planetstore, const vector<vector<double>>& k, const double& no_lines, const double& condition) {
    vector<Planet> pos_update(planetstore.size());
    pos_update = planetstore;

    for (int z = 0; z < no_lines - 1; ++z) {
        pos_update[z] += k[z] * (condition/ 2.0);
    }
    return pos_update;
}

//rk4 member function, computes the effects on a planet, by another planet. Must be provided with G and h, 
//which are constants of the problem. Also must be provided with planets vector. The function returns the same planets vector, with the change to its dynamics applied for that timestep
 vector<Planet> Planet::rk4(vector<Planet>& planets, const double& no_lines, const double& G, const double& h) {
    
     //declare a vector of planets, it will serve to be a copy of the planets at the start of each timestep
     vector<Planet> planetstore(planets.size());

     //generate a copy of the planet, so that if there are multiple planets they can all refer to the original starting position
     planetstore = planets;

     vector<vector<double>>  k1(no_lines - 1, vector<double>(4)); //this vector is a 4xno planets matrix, which holds all rk k1 values
     vector<vector<double>>  k2(no_lines - 1, vector<double>(4));
     vector<vector<double>>  k3(no_lines - 1, vector<double>(4));
     vector<vector<double>>  k4(no_lines - 1, vector<double>(4));

     vector<Planet> pos_update(planetstore.size());

     k1 = planetstore[0].constantfind(planetstore, no_lines, G, h);
     pos_update = planetstore[0].update_pos(planetstore, k1, no_lines, 1.0);
     k2 = planetstore[0].constantfind(pos_update, no_lines, G, h);
     pos_update = planetstore[0].update_pos(planetstore, k2, no_lines, 1.0);
     k3 = planetstore[0].constantfind(pos_update, no_lines, G, h);
     pos_update = planetstore[0].update_pos(planetstore, k3, no_lines, 2.0);
     k4 = planetstore[0].constantfind(pos_update, no_lines, G, h);
     for (int z = 0; z < no_lines - 1; ++z) {
         planets[z] += (k1[z] + k2[z] * 2.0 + k3[z] * 2.0 + k4[z]) / 6.0;
     }
    return planets;
};