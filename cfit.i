%module cfit
%{
#include "cfit.h"
%}

%include "stl.i"
%template(_string_list) std::vector< std::string >;
%include "std_vector.i"
%include "std_string.i"
%include "cfit.h"



