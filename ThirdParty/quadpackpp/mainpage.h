/** \mainpage quadpack++
 
 <b>quadpack++</b> implements adaptive quadrature routines from the 
 <a href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> for an 
 arbitrary floating-point type using C++ templates. The only requirement is 
 that standard arithmetic operators, including mixed arithmetic with standard 
 types such as "int", be overloaded for the floating-point type. The core 
 routines are implemented in header files, requiring nothing to be built or 
 linked against.
 
 The library is intended for use with multiple-precision calculations. One 
 package supporting this is <a href="http://www.mpfr.org/">MPFR</a>, whose 
 website has links to wrappers that define a multiple-precision floating-point 
 type as a C++ class.

 <h2><a href="http://sourceforge.net/projects/quadpackpp/develop">
 Sourceforge repository</a></h2>

 <h2>Design principles</h2>

 Source code for the quadrature routines is taken from the GSL Integration 
 package. All adaptive routines (QAG, QAGS, QAGI, etc.) in this package require 
 a "workspace" structure to be allocated which is then passed as an argument 
 to the routine and subsequent helper routines. quadpack++ uses an object 
 oriented approach: \ref Workspace is a class, and the adaptive 
 quadrature routine(s) such as QAG are member functions.
 
 Memory management is handled by by class constructors and destructors, and 
 numerical difficulties encountered in the adaptive quadrature routines throw 
 an exception. 

 <h2>Gauss-Kronrod rules</h2>

 Gauss-Kronrod quadrature makes robust error estimation and control in the 
 QUADPACK/GSL routines possible. In the original packages, routines of order 
 2m+1 for only m=7, 10, 15, 20, 25, and 30 are provided, and the abscissae 
 and weights for each are tabulated as constants in the source codes. To 
 accommodate calculations to arbitrary precision, the \ref GaussKronrod class 
 constructor computes abscissae and weights for a single rule to a precision 
 commensurate with "machine epsilon" for the floating-point type with which
 the template parameter is instantiated. 

 The weights and abscissae are available through member functions of 
 \ref GaussKronrod, albeit in the compact QUADPACK storage format. In addition 
 to numerical quadrature, therefore, this class can be used to compute 
 abscissae and weights for rules of order 2m+1 for practically any m and 
 to a degree of numerical precision restricted only by that of the template 
 parameter.

 <h2>Machine parameters</h2>
 
 The \ref Machar class constructor computes "machine" parameters needed by methods 
 in the \ref GaussKronrod class and the Workspace class, the essential ones being 
 machine epsilon, overflow limit, and underflow limit. It is based on the 
 MACHAR program described in
 - William Cody, <em>Algorithm 665: MACHAR, a subroutine to dynamically 
 determine machine parameters</em>, ACM Trans. Math. Software, Vol 14 (1988) 
 pages 303-311.
 */
 