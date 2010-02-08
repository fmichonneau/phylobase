LDADD       = @top_builddir@/ncl/libncl.la
AM_CPPFLAGS = -I@top_srcdir@/ncl
bin_PROGRAMS = NEXUSnormalizer NEXUSvalidator NEXUSinspector NEX_us2ml NEXUSunion NCLconverter


NEXUSnormalizer_SOURCES = normalizer.cpp normalizer.h

NEXUSvalidator_SOURCES = normalizer.cpp normalizer.h
NEXUSvalidator_CPPFLAGS = $(AM_CPPFLAGS) -DJUST_VALIDATE_NEXUS

NEXUSinspector_SOURCES = normalizer.cpp normalizer.h
NEXUSinspector_CPPFLAGS = $(AM_CPPFLAGS) -DJUST_REPORT_NEXUS

NEX_us2ml_SOURCES = normalizer.cpp normalizer.h us2ml.cpp
NEX_us2ml_CPPFLAGS = $(AM_CPPFLAGS) -DTO_NEXML_CONVERTER

NEXUSunion_SOURCES = normalizer.cpp normalizer.h
NEXUSunion_CPPFLAGS = $(AM_CPPFLAGS) -DMULTIFILE_NEXUS_READER

NCLconverter_SOURCES = normalizer.cpp normalizer.h us2ml.cpp converter.cpp
NCLconverter_CPPFLAGS = $(AM_CPPFLAGS) -DNCL_CONVERTER_APP
