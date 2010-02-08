#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(NCL, 2.1.09, paul.lewis@uconn.edu)
AC_CONFIG_HEADERS([config.h])

AC_PREREQ(2.57)

#	Directory that contains install-sh and other auxiliary files
AC_CONFIG_AUX_DIR([config])

################################################################################
#	According to (http://www.mail-archive.com/autoconf@gnu.org/msg14232.html)
#		this macro should be after AC_INIT but before AM_INIT_AUTOMAKE
################################################################################
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.6 foreign dist-zip tar-ustar filename-length-max=299])

AC_CONFIG_SRCDIR([ChangeLog])

# Checks for programs.
AC_PROG_LN_S
AC_PROG_INSTALL
AC_PROG_LIBTOOL

################################################################################
# calling AC_PROG_CXX resets CXXFLAGS, we use our own flags set in the 
#	the AX_CXXFLAGS macro below.   
# So we cache the argument to configure
#	here in ARG_CXX_FLAGS (so that our macro does not override them)
ARG_CXX_FLAGS="$CXXFLAGS"
################################################################################
ARG_CXX_FLAGS="$CXXFLAGS"
AC_PROG_CXX
#Ranlib handled by check for libtool
#AC_PROG_RANLIB
CXXFLAGS="$ARG_CXX_FLAGS"
AX_CXXFLAGS


################################################################################
# the install-check requires python with subprocess to actually perform a check
# these tests should not cause failure of configure if python is not found.
################################################################################
#AC_PATH_PROG([SWIG], [swig], [], [$PATH])

#if test -z "$SWIG"
#then
#	if test -z "$USE_SWIG"
#	then
#		MAYBE_SWIG_DIR=
#	else
#		AC_MSG_ERROR([--with-swig was specified, but swig was not found on the path])
#	fi
#else
#	if test -z "$USE_SWIG"
#	then
#		MAYBE_SWIG_DIR=
#	else
#		MAYBE_SWIG_DIR="swig"
#	fi
#fi
#
#PYTHON_SDK_DIR="$PYTHON_SDK_DIR"
#if test -z "$PYTHON_SDK_DIR"
#then
#	if ! test -z "$MAYBE_SWIG_DIR"
#	then
#		AC_MSG_WARN([Currently the swig compilation requires that you specify the PYTHON_SDK_DIR environmental variable. It is very likely that the swig compilation will fail without this variable. The variable should be set to the directory that contains the python installation with subdirectories bin, include, and lib ])
#	fi
#fi
#
#PYTHON_NAME_WITH_VERSION="$PYTHON_NAME_WITH_VERSION"
#if test -z "$PYTHON_NAME_WITH_VERSION"
#then
#	if ! test -z "$MAYBE_SWIG_DIR"
#	then
#		AC_MSG_WARN([Currently the swig compilation requires that you specify the PYTHON_NAME_WITH_VERSION environmental variable. It is very likely that the swig compilation will fail without this variable. It should be set to your python version, such as python2.6 ])
#	fi
#fi

#MAYBE_SWIG_DIR=
#USE_SWIG=
#AC_ARG_WITH(
#	[swig],
#	AC_HELP_STRING(
#		[--with-swig-python=yes],
#		[Compiles NCL for use with SWIG bindings for use with languages other than C++ - currently only python examples are included]
#		),
#	[
#	if test "$withval" = "yes"; then
#		AM_PATH_PYTHON
#		AC_PROG_SWIG(1.3.17)
#		SWIG_ENABLE_CXX
#		SWIG_PYTHON
#		MAYBE_SWIG_DIR="swig"
#	fi
#	])


AC_SUBST(CXXFLAGS)
#AC_CXX_RTTI

CPPFLAGS="-I\$(top_srcdir) $CPPFLAGS"
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

AC_ARG_WITH(
	[constfuncs],
	AC_HELP_STRING(
		[--with-constfuncs=yes],
		[Defines the NCL_CONST_FUNCS macro so functions that should be const are declared as such]
		),
	[
	if test "$withval" = "yes"; then
		CPPFLAGS="-DNCL_CONST_FUNCS $CPPFLAGS"
	fi
	])




AC_SUBST(MAYBE_TEST_DIR)
AC_SUBST(MAYBE_SWIG_DIR)
#AC_SUBST(PYTHON_SDK_DIR)
#AC_SUBST(PYTHON_NAME_WITH_VERSION)

AC_PROG_MAKE_SET

AC_CONFIG_FILES([					\
	Makefile 						\
	example/Makefile				\
	example/basiccmdline/Makefile	\
	example/basicfactory/Makefile	\
	example/gapcode/Makefile		\
	example/nclsimplest/Makefile	\
	example/ncltest/Makefile		\
	example/normalizer/Makefile		\
	example/patristic/Makefile		\
	example/splitsinfile/Makefile	\
	example/translate/Makefile   	\
	ncl/Makefile					\
	test/Makefile                    \
	nclv2.1.pc
	])
AC_OUTPUT
