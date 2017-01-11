#include "utility.h"

double product4(std::vector<double> const& A, std::vector<double> const& B){
	return A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
}

void print4vec(std::vector<double> const& A){
	std::cout << A[0] << "\t" << A[1] << "\t" << A[2] << "\t" << A[3] << std::endl;
}
 
//=======================Rotation==============================================================
// Rotation can take care of 
//---------------------General Euler angles----------------------------------------------------
// This rotation operation takes 3 euler angles and rotate (new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
// the new frame is achieved by the old frame from Z(1_=alpha), X(2_=beta) and Z(3_=gamma)
// R = Z3T*X2T*Z1T
std::vector<double> rotate_ByEuler(std::vector<double> const& A, double & alpha, double & beta, double & gamma){
	double c1 = std::cos(alpha), s1 = std::sin(alpha);
	double c2 = std::cos(beta), s2 = std::sin(beta);
	double c3 = std::cos(gamma), s3 = std::sin(gamma);
	size_t offset = 0;
	std::vector<double> Ap;
	if (A.size() == 4) {offset = 1; Ap.resize(4); Ap[0] = A[0]; }
	else {Ap.resize(3);}
	size_t ix = offset, iy = offset+1, iz = offset+2;

	Ap[ix] = (c1*c3-c2*s1*s3)*A[ix] 	+ (c3*s1+c1*c2*s3)*A[iy] 	+ s2*s3*A[iz];
	Ap[iy] = (-c1*s3-c2*c3*s1)*A[ix] 	+ (c1*c2*c3-s1*s3)*A[iy] 	+ c3*s2*A[iz];
	Ap[iz] = s1*s2*A[ix] 				+ (-c1*s2)*A[iy] 			+ c2*A[iz];
	return Ap;
}

std::vector<double> rotate_back_from_D(std::vector<double> const& A, double & Dx, double & Dy, double & Dz){
	double Dperp = std::sqrt(Dx*Dx + Dy*Dy);
	double D = std::sqrt(Dperp*Dperp + Dz*Dz);
	double c2 = Dz/D, s2 = -Dperp/D;
	double c3 = Dx/Dperp, s3 = -Dy/Dperp;
	size_t offset = 0;
	std::vector<double> Ap;
	if (A.size() == 4) {offset = 1; Ap.resize(4); Ap[0] = A[0];}
	else {Ap.resize(3);}
	size_t ix = offset, iy = offset+1, iz = offset+2;

	Ap[ix] = c3*A[ix] 	+ c2*s3*A[iy] 	+ s2*s3*A[iz];
	Ap[iy] = -s3*A[ix] 	+ c2*c3*A[iy] 	+ c3*s2*A[iz];
	Ap[iz] =            - s2*A[iy] 	    + c2*A[iz];
	return Ap;
}

// Rotation around ith axis, return vector components in the new frame, passive
// dir = 1(x), 2(y), 3(z)
std::vector<double> rotate_ByAxis(std::vector<double> const& A, double & alpha, unsigned int & dir){
	dir -= 1;
	double c1 = std::cos(alpha), s1 = std::sin(alpha);
	size_t offset = 0;
	std::vector<double> Ap;
	if (A.size() == 4) {offset = 1; Ap.resize(4); Ap[0] = A[0];}
	else {Ap.resize(3);}
	size_t i1 = offset + dir, i2 = offset + ((dir+1)%3), i3 = offset + ((dir+2)%3);
	
	Ap[i1] = A[i1];
	Ap[i2] = c1*A[i2] + s1*A[i3];
	Ap[i3] = -s1*A[i2] + c1*A[i3];
	return Ap;
}

//=======================Boost==============================================================

//---------------------General Boost (beta_x, beta_y, beta_z)-----------------------------------------
// This boost operation takes 3 velocity (vx, vy, vz, of the new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
std::vector<double> boost4_By3(std::vector<double> const& A, std::vector<double> const& v){
	double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	double absv = std::sqrt(v2)+1e-32;
	double nx = v[0]/absv, ny = v[1]/absv, nz = v[2]/absv;
	double gamma = 1./sqrt(1. - v2 + 1e-32);
	double gb = gamma - 1.;
	double gb_vecn_dot_vecA = gb*(nx*A[1] + ny*A[2] + nz*A[3]);
	double gammaA0 = gamma*A[0];
	std::vector<double> Ap(4);
	Ap[0] = gamma*(A[0] - v[0]*A[1] - v[1]*A[2] - v[2]*A[3]);
	Ap[1] = -gammaA0*v[0] + A[1] + nx*gb_vecn_dot_vecA;
	Ap[2] = -gammaA0*v[1] + A[2] + ny*gb_vecn_dot_vecA;
	Ap[3] = -gammaA0*v[2] + A[3] + nz*gb_vecn_dot_vecA;
	return Ap;
}

// This boost operation takes 4 velocity (u0, u1, u2, u3 of the new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
std::vector<double> boost4_By4(std::vector<double> const& A, std::vector<double> const& u){
	double gamma = u[0];
	double gb = gamma - 1., gv = std::sqrt(gamma*gamma - 1.) + 1e-32;
	double nx = u[1]/gv, ny = u[2]/gv, nz = u[3]/gv;
	std::vector<double> Ap(4);
	Ap[0] = gamma*A[0]	 -  u[1]*A[1] 		 - u[2]*A[2]		   -u[3]*A[3];
	Ap[1] = (-u[1])*A[0] + (1. + gb*nx*nx)*A[1] + gb*nx*ny*A[2] 	   + gb*nx*nz*A[3];
	Ap[2] = (-u[2])*A[0] + gb*ny*nx*A[1] 		+ (1. + gb*ny*ny)*A[2] + gb*ny*nz*A[3];
	Ap[3] = (-u[3])*A[0] + gb*nz*nx*A[1] 		+ gb*nz*ny*A[2] 	   + (1. + gb*nz*nz)*A[3];
	return Ap;
}

// Boost along ith axis, return vector components in the new frame, passive
// dir = 1(x), 2(y), 3(z)
std::vector<double> boost4_ByAxis(std::vector<double> const& A, double const & vd, unsigned int & dir){
	dir -= 1;
	double gamma = 1./std::sqrt(1. - vd*vd);
	size_t i1 = 1+dir, i2 = 1+(dir+1)%3, i3 = 1+(dir+2)%3;
	std::vector<double> Ap(4);
	Ap[0] = gamma*A[0] - gamma*vd*A[i1];
	Ap[i1] = -gamma*vd*A[0] + gamma*A[i1];
	Ap[i2] = A[i2];
	Ap[i3] = A[i3];
	return Ap;
}

// boost two 4-vectors to their center of mass frame:
void go_to_CoM(std::vector<double> const& Pcom,
			   std::vector<double> const& A, std::vector<double> const& B,
			   std::vector<double> & Ap, std::vector<double> & Bp){
	// center of mass velocity relative to the present frame
	std::vector<double> vcom(3);
	vcom[0] = Pcom[1]/Pcom[0]; vcom[1] =  Pcom[2]/Pcom[0]; vcom[2] =  Pcom[3]/Pcom[0];
	boost4_By3(A, vcom);
	boost4_By3(B, vcom);
}

