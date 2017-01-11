#ifndef UTILITY_H
#define UTILITY_H
#include <cmath>
#include <vector>
#include <iostream>

double product4(std::vector<double> const& A, std::vector<double> const& B);
void print4vec(std::vector<double> const& A);
 
//=======================Rotation==============================================================
// Rotation can take care of 
//---------------------General Euler angles----------------------------------------------------
// This rotation operation takes 3 euler angles and rotate (new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
// the new frame is achieved by the old frame from Z(1_=alpha), X(2_=beta) and Z(3_=gamma)
// R = Z3T*X2T*Z1T
std::vector<double> rotate_ByEuler(std::vector<double> const& A, double & alpha, double & beta, double & gamma);

std::vector<double> rotate_back_from_D(std::vector<double> const& A, double & Dx, double & Dy, double & Dz);
// Rotation around ith axis, return vector components in the new frame, passive
// dir = 1(x), 2(y), 3(z)
std::vector<double> rotate_ByAxis(std::vector<double> const& A, double & alpha, unsigned int & dir);

//=======================Boost==============================================================

//---------------------General Boost (beta_x, beta_y, beta_z)-----------------------------------------
// This boost operation takes 3 velocity (vx, vy, vz, of the new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
std::vector<double> boost4_By3(std::vector<double> const& A, std::vector<double> const& v);

// This boost operation takes 4 velocity (u0, u1, u2, u3 of the new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
std::vector<double> boost4_By4(std::vector<double> const& A, std::vector<double> const& u);
// Boost along ith axis, return vector components in the new frame, passive
// dir = 1(x), 2(y), 3(z)
std::vector<double> boost4_ByAxis(std::vector<double> const& A, double const & vd, unsigned int & dir);

// boost two 4-vectors to their center of mass frame:
void go_to_CoM(std::vector<double> const& Pcom,
			   std::vector<double> const& A, std::vector<double> const& B,
			   std::vector<double> & Ap, std::vector<double> & Bp);

struct particle{
	std::vector<double> p;
	std::vector<double> x;
	double t_last;
	double weight;
	bool freeze;
	double T_dec;
	std::vector<double> v_dec;
	particle() : p(4), x(4), t_last(0.0), weight(1.0), freeze(false), T_dec(-1.0), v_dec(3){
	}
};

#endif
