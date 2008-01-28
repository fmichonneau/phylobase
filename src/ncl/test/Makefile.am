installcheck-local:
	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/normalizer $(srcdir)/NTSValidIn $(srcdir)/NTSValidOut
	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/normalizer $(top_srcdir)/data/characters.nex $(srcdir)/data
	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/normalizer $(top_srcdir)/data/sample.tre $(srcdir)/data
