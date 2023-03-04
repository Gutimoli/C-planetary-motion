//Planet.h
#ifndef PLANETS_H
#define PLANETS_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>

using namespace std;

/**
 * @class Planet
 * @brief This class serves to represent a planet
 *
 * The class stores a planet's positions (x,y) and velocities (xdot,ydot) in a 2D cartesian system,
 * it also stores the planet's mass (m). All of these variables remain private. The class has defined
 * setter and getter functions, so that the variables can still be accessed and modified. Operators
 * +=, and = are overloaded. += adds a vector of positions and velocities to the planet, hence
 * updating its position and velocities. = sets a planet's variables equal to another planets
 * variables. There are also a member functions, which help to compute the change in dynamics of a 
 * set of planets as determined by a 4th order Runge-Kutta numerical scheme.
 *
 * @author Daniel Gutierrez
 * @date January 21, 2023
 */

//declare namespace for robustness
namespace Universe {
    class Planet {

    public:

        //Constructor with known values
        Planet(const double& x_, const double& y_, const double& xdot_, const double& ydot_, const double& m_);

        //Constructor with no values
        Planet();

        //Destructor
        ~Planet();

        //Getter functions
        double GetX();
        double GetY();
        double GetXdot();
        double GetYdot();
        double GetM();

        //Setter functions
        void SetX(double x_);
        void SetY(double y_);
        void SetXdot(double xdot_);
        void SetYdot(double ydot_);
        void SetM(double m_);

        //Operator overload of = to make all variables of a planet equal to another plant's variables.
        Planet& operator=(const Planet& planet);

        //Operator overload of += to add a vector of positions and velocities to a planet's position and velocities.
        Planet& operator+=(std::vector<double> a);

        //rk4 member function, computes the effects on a planet, by another planet. Must be provided with G and h, 
        //which are constants of the problem.Also provide it with planets The function returns the same planet, with the change to its constants applied
        vector<Planet> rk4(vector<Planet>& planets, const double& no_lines, const double& G, const double& h);

        //Finds the rk specific constant
        vector<vector<double>> constantfind(vector<Planet>&, const double& no_lines, const double& G, const double& h);

        //position update function
        vector<Planet> update_pos(vector<Planet>& planetstore, const vector<vector<double>>& k, const double& no_lines, const  double& conditi);

    private:

        double x; //position x of planet
        double y; //position y of planet
        double xdot; //velocity in x direction of planet
        double ydot; //velocity in y direction of planet
        double m; //mass of planet
    };
}

using namespace Universe;

//Simple, yet useful vector overloading
//Adding two vectors element wise
template <typename U>
std::vector<U> operator+(std::vector<U> v1, const std::vector<U>& v2) {
    std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), [] (const U& a, const U& b) {return a + b; });
    return v1;
}

//Multiplying element wise a vector by a scalar
template <typename U>
std::vector<U> operator*(std::vector<U> v1, const double& k) {
    std::transform(v1.begin(), v1.end(), v1.begin(), [k] (const U& a) {return a * k; });
    return v1;
}
//divide element wise by a scalar
template <typename U>
std::vector<U> operator/(std::vector<U> v1, const double& k) {
    std::transform(v1.begin(), v1.end(), v1.begin(), [k](const U& a) {return a / k; });
    return v1;
}
#endif
