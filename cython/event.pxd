from libcpp cimport bool
from libcpp.vector cimport vector

#------------Import c++ utility function for rotation and transformation-----------
cdef extern from "../src/utility.h":
	cdef double product4(const vector[double] & A, const vector[double] & B)
	cdef void print4vec(vector[double] & Ap, const vector[double] & A)
	cdef void rotate_ByEuler(vector[double] & Ap, const vector[double] & A, double & alpha, double & beta, double & gamma)
	cdef void  rotate_ByAxis(vector[double] & Ap, const vector[double] & A, double & alpha, unsigned int & dir)
	cdef void  rotate_back_from_D(vector[double] & Ap, const vector[double] & A, double & Dx, double & Dy, double & Dz)
	cdef void  boost4_By3(vector[double] & Ap, const vector[double] & A, const vector[double] & v)
	cdef void  boost4_By4(vector[double] & Ap, const vector[double] & A, const vector[double] & u)
	cdef void  boost4_ByAxis(vector[double] & Ap, const vector[double] & A, const double & vd, unsigned int & dir)
	cdef void go_to_CoM(const vector[double] & Pcom,
			   const vector[double] & A, const vector[double] & B,
			   vector[double] & Ap, vector[double] & Bp);

