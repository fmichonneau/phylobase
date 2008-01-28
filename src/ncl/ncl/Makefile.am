pkglib_LTLIBRARIES = libncl.la
libncl_ladir = $(includedir)/ncl
libncl_la_LDFLAGS = -release @PACKAGE_VERSION@


NCL_HEADER_FILES = \
	ncl.h \
	nxsassumptionsblock.h \
	nxsblock.cpp \
	nxsblock.h \
	nxscharactersblock.h \
	nxsdatablock.h \
	nxsdefs.h \
	nxsdiscretedatum.h \
	nxsdiscretematrix.h \
	nxsdistancedatum.h \
	nxsdistancesblock.h \
	nxsexception.h \
	nxsindent.h \
	nxsreader.h \
	nxssetreader.h \
	nxsstring.h \
	nxstaxablock.h \
	nxstoken.h \
	nxstreesblock.h \
	nxsunalignedblock.h \
	nxsunalignedblock.inl

NCL_SOURCE_FILES = \
	nxsassumptionsblock.cpp \
	nxsblock.cpp \
	nxscharactersblock.cpp \
	nxsdatablock.cpp \
	nxsdiscretedatum.cpp \
	nxsdiscretematrix.cpp \
	nxsdistancedatum.cpp \
	nxsdistancesblock.cpp \
	nxsexception.cpp \
	nxsreader.cpp \
	nxssetreader.cpp \
	nxsstring.cpp \
	nxstaxablock.cpp \
	nxstoken.cpp \
	nxstreesblock.cpp \
	nxsunalignedblock.cpp

libncl_la_HEADERS = $(NCL_HEADER_FILES)
libncl_la_SOURCES = $(NCL_SOURCE_FILES)

#lib_LIBRARIES = libncl.a
#libncl_adir = $(includedir)/ncl
#libncl_a_HEADERS = $(NCL_HEADER_FILES)
#libncl_a_SOURCES = $(NCL_SOURCE_FILES)
	
