LDADD       = @top_builddir@/ncl/libncl.la
AM_CPPFLAGS = -I@top_srcdir@/ncl
bin_PROGRAMS = normalizer validator
normalizer_SOURCES = normalizer.cpp normalizer.h

validator_SOURCES = normalizer.cpp normalizer.h

validator_CPPFLAGS = $(AM_CPPFLAGS) -DJUST_VALIDATE_NEXUS
