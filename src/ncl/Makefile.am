SUBDIRS = ncl example $(MAYBE_TEST_DIR)

EXTRA_DIST = data html test ncl/nxsemptyblock.h ncl/nxsemptyblock.cpp

dist-hook:
	find "$(distdir)/data" -depth -and -name .svn -and -type d -and -exec rm -rf {} \; 
	find "$(distdir)/html" -depth -and -name .svn -and -type d -and -exec rm -rf {} \; 
	find "$(distdir)/test" -depth -and -name .svn -and -type d -and -exec rm -rf {} \; 

