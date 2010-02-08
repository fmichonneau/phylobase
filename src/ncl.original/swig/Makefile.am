BUILT_SOURCES = $(srcdir)/ncl_wrap.cpp
SWIG_FILES = ncl.i

pkgpython_PYTHON = nclwrapper.py
pkgpyexec_LTLIBRARIES = _nclwrapper.la
_nclwrap_la_SOURCES = $(srcdir)/ncl_wrap.cpp $(SWIG_SOURCES)
_nclwrap_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/ncl
_nclwrap_la_LDFLAGS = -module
_nclwrap_la_LIBADD = ../ncl/libncl.la

ncl_wrap.cpp : $(SWIG_FILES)
	 $(SWIG) $(SWIG_PYTHON_OPT) -DNCL_COULD_BE_CONST=const -DNCL_CONST_FUNCS -I$(top_srcdir) -I$(top_srcdir)/ncl -o $@ $<
