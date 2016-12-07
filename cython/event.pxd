from libcpp cimport bool
from libcpp.vector cimport vector

#------------Import c++ utility function for rotation and transformation-----------
cdef extern from "../src/utility.h":
	cdef double product4(const vector[double] & A, const vector[double] & B)
	cdef void print4vec(const vector[double] & A)
	cdef void rotate_Euler(const vector[double] & A, vector[double] & Ap, double alpha, double beta, double gamma)
	cdef void rotate_axis(const vector[double] & A, vector[double] & Ap, double alpha, unsigned int dir)
	cdef void boost_by3(const vector[double] & A, vector[double] & Ap, const vector[double] v)
	cdef void boost_by4(const vector[double] & A, vector[double] & Ap, const vector[double] u)
	cdef void boost_axis(const vector[double] & A, vector[double] & Ap, const double vd, unsigned int dir)
	cdef void go_to_CoM(const vector[double] & Pcom,
			   const vector[double] & A, const vector[double] & B,
			   vector[double] & Ap, vector[double] & Bp);

