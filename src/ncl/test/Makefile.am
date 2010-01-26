installcheck-local:
	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/NEXUSnormalizer $(top_srcdir)/data/characters.nex $(srcdir)/data
	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/NEXUSnormalizer $(top_srcdir)/data/sample.tre $(srcdir)/data
	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/NEXUSnormalizer $(srcdir)/NTSValidIn $(srcdir)/NTSValidOut
	$(PYTHON) $(srcdir)/roundTripNCLTest.py -e $(bindir)/NEXUSnormalizer $(srcdir)/ExternalValidIn $(srcdir)/ExternalValidOut
	$(PYTHON) $(srcdir)/roundTripNCLTest.py -i -e $(bindir)/NEXUSnormalizer $(srcdir)/ExternalInvalidIn
	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/ncltest $(srcdir)/OldValidIn $(srcdir)/OldValidOut
	$(PYTHON) $(srcdir)/roundTripNCLTest.py -o $(bindir)/NEXUStranslate  $(srcdir)/TranslateIn $(srcdir)/TranslateOut
	$(PYTHON) $(srcdir)/roundTripNCLTest.py -a $(bindir)/NEXUSnormalizer 
