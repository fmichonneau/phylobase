#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(NCL, 2.0.04, paul.lewis@uconn.edu)
AC_CONFIG_HEADERS([config.h])

AC_PREREQ(2.57)

#	Directory that contains install-sh and other auxiliary files
AC_CONFIG_AUX_DIR([config])

################################################################################
#	According to (http://www.mail-archive.com/autoconf@gnu.org/msg14232.html)
#		this macro should be after AC_INIT but before AM_INIT_AUTOMAKE
################################################################################
AC_CONFIG_MACRO_DIR(config/m4)

AM_INIT_AUTOMAKE([1.6 foreign dist-zip tar-ustar filename-length-max=299])

AC_CONFIG_SRCDIR([ncl/ncl.h])

# Checks for programs.
AC_PROG_LN_S
AC_PROG_INSTALL

################################################################################
# calling AC_PROG_CXX resets CXXFLAGS, we use our own flags set in the 
#	the AX_CXXFLAGS macro below.   
# So we cache the argument to configure
#	here in ARG_CXX_FLAGS (so that our macro does not override them)
ARG_CXX_FLAGS="$CXXFLAGS"
################################################################################
ARG_CXX_FLAGS="$CXXFLAGS"
AC_PROG_CXX
AC_PROG_RANLIB
AC_PROG_LIBTOOL
CXXFLAGS="$ARG_CXX_FLAGS"
AX_CXXFLAGS
AC_SUBST(CXXFLAGS)


# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MKTIME
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memmove strchr strtol])

################################################################################
# the install-check requires python with subprocess to actually perform a check
# these tests should not cause failure of configure if python is not found.
################################################################################
AC_PATH_PROG([PYTHON], [python], [], [$PATH])
AC_PYTHON_MODULE([subprocess], [], [])  
AC_SUBST(PYTHON)

if test -z "$PYTHON"
then
	MAYBE_TEST_DIR=
else
	MAYBE_TEST_DIR="test"
fi

AC_SUBST(MAYBE_TEST_DIR)

AC_PROG_MAKE_SET

AC_CONFIG_FILES([					\
	Makefile 						\
	example/Makefile				\
	example/ncltest/Makefile		\
	example/nclsimplest/Makefile	\
	example/basiccmdline/Makefile	\
	example/normalizer/Makefile		\
	ncl/Makefile					\
	test/Makefile
	])
AC_OUTPUT
