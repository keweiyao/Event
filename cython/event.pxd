from libcpp cimport bool
from libcpp.vector cimport vector

#------------Import c++ utility function for rotation and transformation-----------
cdef extern from "../src/utility.h":
	cdef double product4(const vector[double] & A, const vector[double] & B)
	cdef void print4vec(const vector[double] & A)
	cdef vector[double] rotate_ByEuler(const vector[double] & A, double & alpha, double & beta, double & gamma)
	cdef vector[double] rotate_ByAxis(const vector[double] & A, double & alpha, unsigned int & dir)
	cdef vector[double] rotate_back_from_D(const vector[double] & A, double & Dx, double & Dy, double & Dz)
	cdef vector[double] boost4_By3(const vector[double] & A, const vector[double] & v)
	cdef vector[double] boost4_By4(const vector[double] & A, const vector[double] & u)
	cdef vector[double] boost4_ByAxis(const vector[double] & A, const double & vd, unsigned int & dir)
	cdef void go_to_CoM(const vector[double] & Pcom,
			   const vector[double] & A, const vector[double] & B,
			   vector[double] & Ap, vector[double] & Bp);

