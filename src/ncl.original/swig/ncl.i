// minimal wrapper around NCL
/* Setting modulename to "ncl" with directors enabled */
%module(directors="0") nclwrapper
%{
 #include "ncl.h"
%}


/* Redirecting to C++ exception found running wrapper code */
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}

/* Allowing synonymous employment of certain STL classes */
typedef std::string string;
typedef std::ostream ostream;
typedef std::istream istream;
typedef long file_pos;

/* Enabling default typemaps for std::strings */
%include std_string.i

/* Ignored members of NCL */
%ignore NxsBlock::CloneBlock;
%ignore	NxsCharactersBlock::GetStateSymbolIndex;


/* Conditional processing based on language of wrappers being produced */


%include "std_vector.i"

namespace std {
   %template(vectori) vector< int >;
   %template(vectorvi) vector< vector<int> >;
   %template(vectord) vector< double >;
   %template(vectorst) vector< string >;

class exception {
	public:
		exception();
		exception(const exception& rhs);
		virtual ~exception();
		virtual const char *what(void);
};

};


%include nxsexception.h
%include nxstoken.h
%include nxsblock.h
%include nxsreader.h
%include nxstaxablock.h
%include nxstreesblock.h
%include nxscharactersblock.h
%include nxsassumptionsblock.h
%include nxsdatablock.h
%include nxsdistancesblock.h
%include nxssetreader.h
%include nxspublicblocks.h
%include nxsmultiformat.h

