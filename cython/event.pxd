from libcpp cimport bool
from libcpp.vector cimport vector

#------------Import c++ utility function for rotation and transformation-----------
cdef extern from "../src/utility.h":
	cdef double product4(const vector[double] & A, const vector[double] & B) 
	cdef void  rotate_back_from_D(vector[double] & Ap, const vector[double] & A, double & Dx, double & Dy, double & Dz) 
	cdef void  boost4_By3(vector[double] & Ap, const vector[double] & A, const vector[double] & v) 
	cdef void  boost4_By3_back(vector[double] & Ap, const vector[double] & A, const vector[double] & v) 


