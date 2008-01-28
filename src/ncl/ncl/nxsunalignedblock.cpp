//	Copyright (C) 2007 Paul O. Lewis
//
//	This file is part of NCL (Nexus Class Library) version 2.0.
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#include "ncl.h"
#include "nxsunalignedblock.h"

//@POL Note: This file is not yet ready for use (Paul Lewis, 19-May-2007)

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `id' to "UNALIGNED", `taxa' to `tb', `assumptionsBlock' to `ab', `ntax' and `ntaxTotal' to 0, `newtaxa' 
|	and `respectingCase' to false, `labels' to true, `datatype' to `NxsUnalignedBlock::standard', `missing' to '?', and
|	`taxonPos' and `activeTaxon' to NULL. The `equates' map and `umatrix' vector are both cleared. The ResetSymbols 
|	function is called to reset the `symbols' data member. Assumes that `tb' and `ab' point to valid NxsTaxaBlock and 
|	NxsAssumptionsBlock objects, respectively.
*/
NxsUnalignedBlock::NxsUnalignedBlock(
  NxsTaxaBlock * tb,			/* is the taxa block object to consult for taxon labels */
  NxsAssumptionsBlock * ab)		/* is the assumptions block object to consult for exclusion sets */
  : NxsBlock()
	{
	assert(tb != NULL);
	assert(ab != NULL);

	taxa				= tb;
	assumptionsBlock	= ab;
	id					= "UNALIGNED";

	// These need to be initialized to NULL so Reset member function will not try to delete them
	taxonPos			= NULL;
	activeTaxon			= NULL;
	symbols				= NULL;

	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes any memory allocated to the arrays `symbols', `taxonPos' and `activeTaxon'. Flushes the containers 
|	`equates', and `umatrix'.
*/
NxsUnalignedBlock::~NxsUnalignedBlock()
	{
	Reset();

	equates.clear();

	if (symbols != NULL)
		delete [] symbols;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns NxsUnalignedBlock object to the state it was in when first created. See NxsUnalignedBlock constructor for
|	details.
*/
void NxsUnalignedBlock::Reset()
	{ 
	NxsBlock::Reset();
	ntax				= 0;
	ntaxTotal			= 0;
	newtaxa				= false;
	respectingCase		= false;
	labels				= true;
	datatype			= NxsUnalignedBlock::standard;
	missing				= '?';

	ResetSymbols();	// also resets equates

	umatrix.clear();

	if (taxonPos != NULL)
		{
		delete [] taxonPos;
		taxonPos = NULL;
		}

	if (activeTaxon != NULL)
		{
		delete [] activeTaxon;
		activeTaxon = NULL;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets standard symbol set after a change in `datatype' is made. Also flushes equates list and installs standard 
|	equate macros for the current `datatype'.
*/
void NxsUnalignedBlock::ResetSymbols()
	{
	// Symbols might be NULL (if a different NxsUnalignedBlock has consumed this one
	//
	if (symbols == NULL)
		{
		symbols = new char[NCL_MAX_STATES+1];
		symbols[0] = '0';
		symbols[1] = '1';
		symbols[2] = '\0';
		}

	switch(datatype)
		{
		case NxsUnalignedBlock::dna:
			strcpy(symbols, "ACGT");
			break;

		case NxsUnalignedBlock::rna:
			strcpy(symbols, "ACGU");
			break;

		case NxsUnalignedBlock::nucleotide:
			strcpy(symbols, "ACGT");
			break;

		case NxsUnalignedBlock::protein:
			strcpy(symbols, "ACDEFGHIKLMNPQRSTVWY*");
			break;

		default:
			strcpy(symbols, "01");
		}

	// Setup standard equates
	//
	equates.clear();
	this->equates = GetDefaultEquates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes (i.e., excludes from further analyses) taxa whose indices are contained in the set `delset'. The taxon 
|	indices refer to those in the TAXA block, not current indices (originals will equal current ones if number of taxa 
|	in TAXA block equals number of taxa in MATRIX command). Returns the number of taxa actually deleted (some may have 
|	already been deleted)
*/
unsigned NxsUnalignedBlock::ApplyDelset(
  NxsUnsignedSet & delset)	/* set of taxon indices to delete in range [0..`ntaxTotal') */
	{
	assert(activeTaxon != NULL);
	assert(taxonPos != NULL);

	unsigned num_deleted = 0;
	unsigned k;

	NxsUnsignedSet::const_iterator i;
	for (i = delset.begin(); i != delset.end(); ++i)
		{
		k = taxonPos[*i];
		if (k == UINT_MAX)
			continue;

		// k less than UINT_MAX means data was supplied for this taxon and therefore it can be deleted
		if (activeTaxon[k] == true)
			num_deleted++;
		activeTaxon[k] = false;
		}

	return num_deleted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores (i.e., includes in further analyses) taxa whose indices are contained in the set `restoreset'. The taxon 
|	indices refer to original taxon indices, not current indices (originals will equal current ones if number of taxa 
|	in TAXA block equals number of taxa in MATRIX command).
*/
unsigned NxsUnalignedBlock::ApplyRestoreset(
  NxsUnsignedSet & restoreset)	/* set of taxon indices to restore in range [0..`ntaxTotal') */
	{
	assert(activeTaxon != NULL);
	assert(taxonPos != NULL);

	unsigned num_restored = 0;
	unsigned k;

	NxsUnsignedSet::const_iterator i;
	for (i = restoreset.begin(); i != restoreset.end(); i++)
		{
		k = taxonPos[*i];
		if (k == UINT_MAX)
			continue;

		// k less than UINT_MAX means data was supplied for this taxon and therefore it can be restored
		if (activeTaxon[k] == false)
			num_restored++;
		activeTaxon[k] = true;
		}

	return num_restored;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides a dump of the contents of the `umatrix' variable. Useful for testing whether data is being read as 
|	expected. If `marginText' is NULL, output is flush left. If each line of output should be prefaced with 
|	a tab character, specify "\t" for `marginText'.
*/
void NxsUnalignedBlock::DebugShowMatrix(
  ostream & out,		/* is the output stream on which to print */
  char * marginText)	/* is text printed first on each line */
	{
	assert(taxonPos != NULL);

	unsigned i, k;
	unsigned width = taxa->GetMaxTaxonLabelLength();
	unsigned first_taxon = UINT_MAX;

	for (i = 0; i < ntaxTotal; ++i)
		{
		if (taxonPos[i] == UINT_MAX)
			{
			// No data was provided for the taxon at index i in the TAXA block
			continue;
			}
		else
			{
			// Data was provided for the taxon at index i in the TAXA block
			// This data is in row taxonPos[i] of umatrix
			if (first_taxon == UINT_MAX)
				first_taxon = i;

			if (marginText != NULL)
				out << marginText;

			NxsString currTaxonLabel = taxa->GetTaxonLabel(i);
				out << currTaxonLabel;

			// Print out enough spaces to even up the left edge of the matrix output
			unsigned currTaxonLabelLen = (unsigned)currTaxonLabel.size();
			unsigned diff = width - currTaxonLabelLen;
			for (k = 0; k < diff + 5; k++)
				out << ' ';
			}

		unsigned row = taxonPos[i];
		UnalignedVect::const_iterator j = umatrix[row].begin();
		for (; j != umatrix[row].end(); ++j)
			{
			out << FormatState(*j);
			}

		out << endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string containing a formatted representation of the state `x'. 
*/
std::string NxsUnalignedBlock::FormatState(
  NxsDiscreteDatum x)		/* is the element of `umatrix' to format */
  const
	{
	// x.states could be one of these:
	//	     ?              NULL
	//	     G              [1][2]
	//	(AG) polymorphic    [2][0][2][1]
	//	{AG} ambiguous      [2][0][2][0]
	//
	// Note that gaps are not allowed, so x.states cannot be this:
	//	     -              [0]

	std::string s;
	if (x.states == NULL)
		s += missing;
	else
		{
		// First element of x is the number of states
		unsigned nsymbols = (unsigned)strlen(symbols);
		unsigned n = x.states[0];
		if (n == 1)
			{
			// unambiguous state
			unsigned v = x.states[1];
			assert(v < nsymbols);
			s += symbols[v];
			}
		else
			{
			// ambiguity or polymorphism
			bool is_polymorphic = (bool)(x.states[n+1] == 1);
			if (is_polymorphic)
				{
				// polymorphic
				s += '(';
				}
			else
				{
				// ambiguous
				s += '{';
				}
			for (unsigned k = 0; k < n; k++)
				{
				unsigned v = x.states[1+k];
				assert(v < nsymbols);
				s += symbols[v];
				}
			if (is_polymorphic)
				s += ')';
			else
				s += '}';
			}
		}

	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a count of the number of taxa for which `activeTaxon' array reports true.
*/
unsigned NxsUnalignedBlock::GetNumActiveTaxa()
	{
	unsigned num_active_taxa = 0;
	for (unsigned i = 0; i < ntax; i++)
		{
		if (activeTaxon[i] == false)
			continue;
		num_active_taxa++;
		}

	return num_active_taxa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	For row `i' in `umatrix', returns the index of the corresponding taxon in the TAXA block. The return value will be 
|	equal to `i' unless data was not provided for some taxa listed in a preceding TAXA block or the taxa were presented
|	in a different order in the TAXA block vs. the UNALIGNED MATRIX command.
*/
unsigned NxsUnalignedBlock::GetOrigTaxonIndex(
  unsigned i)	/* is the row in `umatrix' in range [0..`ntax') */
	{
	assert(taxonPos != NULL);

	unsigned k = 0;
	for (; k < ntaxTotal; ++k)
		{
		if (taxonPos[k] == i)
			break;
		}
	assert(k < ntaxTotal);	// every row in umatrix should correspond to a taxon in the TAXA block
	return k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `ch' can be found in the `symbols' array. The value of `respectingCase' is used to determine 
|	whether or not the search should be case sensitive. Assumes `symbols' is non-NULL.
*/
bool NxsUnalignedBlock::IsInSymbols(
  char ch)	/* the symbol character to search for */
	{
	assert(symbols != NULL);
	unsigned symbolsLength = (unsigned)strlen(symbols);
	bool found = false;
	for (unsigned i = 0; i < symbolsLength; i++)
		{
		char char_in_symbols = (respectingCase ? symbols[i] : (char)toupper(symbols[i]));
		char char_in_question = (respectingCase ? ch : (char)toupper(ch));
		if (char_in_symbols != char_in_question)
			continue;
		found = true;
		break;
		}

	return found;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when DIMENSIONS command needs to be parsed from within the UNALIGNED block. Deals with everything after the 
|	token DIMENSIONS up to and including the semicolon that terminates the DIMENSIONS command.
*/
void NxsUnalignedBlock::HandleDimensions(
  NxsToken & token)			/* the token used to read from `in' */
	{
	for (;;)
		{
		token.GetNextToken();

		if (token.Equals("NEWTAXA"))
			{
			newtaxa = true;
			}
		else if (token.Equals("NTAX")) 
			{
			newtaxa = true;
			DemandEquals(token, "after NTAX in DIMENSIONS command");
			ntax = DemandPositiveInt(token, "NTAX");
			ntaxTotal = (newtaxa ? ntax : taxa->GetNumTaxonLabels());
			}
		else if (token.Equals(";"))
			{
			break;
			}
		}

	if (ntaxTotal < ntax)
		{
		errormsg = "NTAX";
		errormsg += " in ";
		errormsg += id;
		errormsg += " block must be less than or equal to NTAX in TAXA block";
		errormsg += "\nNote: one circumstance that can cause this error is ";
		errormsg += "\nforgetting to specify ";
		errormsg += "NTAX";
		errormsg += " in DIMENSIONS command when ";
		errormsg += "\na TAXA block has not been provided";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	if (newtaxa)
		taxa->Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when the END or ENDBLOCK command needs to be parsed from within the UNALIGNED block. Does two things: 
|~
|	o checks to make sure the next token in the data file is a semicolon
|	o eliminates character labels and character state labels for characters that have been eliminated
|~
*/
void NxsUnalignedBlock::HandleEndblock(
  NxsToken & token)		/* the token used to read from `in' */
	{
	DemandEndSemicolon(token, "END or ENDBLOCK");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when FORMAT command needs to be parsed from within the DIMENSIONS block. Deals with everything after the 
|	token FORMAT up to and including the semicolon that terminates the FORMAT command.
*/
void NxsUnalignedBlock::HandleFormat(
  NxsToken & token)	/* is the token used to read from `in' */
	{
	bool standardDataTypeAssumed = false;
	bool ignoreCaseAssumed = false;

	for (;;)
		{
		token.GetNextToken();

		if (token.Equals("DATATYPE"))
			{
			DemandEquals(token, "after keyword DATATYPE");
			// This should be one of the following: STANDARD, DNA, RNA, NUCLEOTIDE or PROTEIN
			token.GetNextToken();

			if (token.Equals("STANDARD"))
				datatype = standard;
			else if (token.Equals("DNA"))
				datatype = dna;
			else if (token.Equals("RNA"))
				datatype = rna;
			else if (token.Equals("NUCLEOTIDE"))
				datatype = nucleotide;
			else if (token.Equals("PROTEIN"))
				datatype = protein;
			else
				{
				errormsg = token.GetToken();
				errormsg += " is not a valid DATATYPE within a ";
				errormsg += id;
				errormsg += " block";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			if (standardDataTypeAssumed && datatype != standard)
				{
				errormsg = "DATATYPE must be specified first in FORMAT command";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			ResetSymbols();
			}
		else if (token.Equals("RESPECTCASE"))
			{
			if (ignoreCaseAssumed)
				{
				errormsg = "RESPECTCASE must be specified before MISSING and SYMBOLS in FORMAT command";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			standardDataTypeAssumed = true;
			respectingCase = true;
			}
		else if (token.Equals("MISSING"))
			{
			DemandEquals(token, "after keyword MISSING");
			// This should be the missing data symbol (single character)
			token.GetNextToken();

			if (token.GetTokenLength() != 1)
				{
				errormsg = "MISSING symbol should be a single character, but ";
				errormsg += token.GetToken();
				errormsg += " was specified";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			else if (token.IsPunctuationToken() && !token.IsPlusMinusToken())
				{
				errormsg = "MISSING symbol specified cannot be a punctuation token (";
				errormsg += token.GetToken();
				errormsg += " was specified)";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			else if (token.IsWhitespaceToken())
				{
				errormsg = "MISSING symbol specified cannot be a whitespace character (";
				errormsg += token.GetToken();
				errormsg += " was specified)";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			missing = token.GetToken()[0];

			ignoreCaseAssumed = true;
			standardDataTypeAssumed = true;
			}
		else if (token.Equals("SYMBOLS"))
			{
			int numDefStates;
			unsigned maxNewStates;
			switch(datatype)
				{
				case NxsUnalignedBlock::dna:
				case NxsUnalignedBlock::rna:
				case NxsUnalignedBlock::nucleotide:
					numDefStates = 4;
					maxNewStates = NCL_MAX_STATES-4;
					break;
				case NxsUnalignedBlock::protein:
					numDefStates = 21;
					maxNewStates = NCL_MAX_STATES-21;
					break;
				default:
					numDefStates = 0; // replace symbols list for standard datatype
					symbols[0] = '\0';
					maxNewStates = NCL_MAX_STATES;
				}
			DemandEquals(token, "after keyword SYMBOLS");

			// This should be the symbols list
			token.SetLabileFlagBit(NxsToken::doubleQuotedToken);
			token.GetNextToken();

			token.StripWhitespace();
			unsigned numNewSymbols = token.GetTokenLength();

			if (numNewSymbols > maxNewStates)
				{
				errormsg = "SYMBOLS defines ";
				errormsg += numNewSymbols;
				errormsg += " new states but only ";
				errormsg += maxNewStates;
				errormsg += " new states allowed for this DATATYPE";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			NxsString t = token.GetToken();
			unsigned tlen = (unsigned)t.size();

			// Check to make sure user has not used any symbols already in the
			// default symbols list for this data type
			for (unsigned i = 0; i < tlen; i++)
				{
				if (IsInSymbols(t[i]))
					{
					errormsg = "The character ";
					errormsg += t[i];
					errormsg += " defined in SYMBOLS has already been predefined for this DATATYPE";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
				}

			// If we've made it this far, go ahead and add the user-defined
			// symbols to the end of the list of predefined symbols
			strcpy(symbols + numDefStates, t.c_str());

			ignoreCaseAssumed = true;
			standardDataTypeAssumed = true;
			}

		else if (token.Equals("EQUATE"))
			{
			DemandEquals(token, "after keyword EQUATE");

			// This should be a double-quote character
			token.GetNextToken();

			if (!token.Equals("\""))
				{
				errormsg = "Expecting '\"' after keyword EQUATE but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			// Loop until second double-quote character is encountered
			for (;;)
				{
				token.GetNextToken();
				if (token.Equals("\""))
					break;

				// If token is not a double-quote character, then it must be the equate symbol (i.e., the 
				// character to be replaced in the data matrix)
				if (token.GetTokenLength() != 1)
					{
					errormsg = "Expecting single-character EQUATE symbol but found ";
					errormsg += token.GetToken();
					errormsg += " instead";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}

				// Check for bad choice of equate symbol
				NxsString t = token.GetToken();
				char ch = t[0];
				bool badEquateSymbol = false;

				// The character '^' cannot be an equate symbol
				if (ch == '^')
					badEquateSymbol = true;

				// Equate symbols cannot be punctuation (except for + and -)
				if (token.IsPunctuationToken() && !token.IsPlusMinusToken())
					badEquateSymbol = true;

				// Equate symbols cannot be same as matchchar, missing, or gap
				if (ch == missing)
					badEquateSymbol = true;

				// Equate symbols cannot be one of the state symbols currently defined
				if (IsInSymbols(ch))
					badEquateSymbol = true;

				if (badEquateSymbol)
					{
					errormsg = "EQUATE symbol specified (";
					errormsg += token.GetToken();
					errormsg += ") is not valid; must not be same as missing, \nmatchchar, gap, state symbols, or any of the following: ()[]{}/\\,;:=*'\"`<>^";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}

				NxsString k = token.GetToken();

				DemandEquals(token, "in EQUATE definition");
				
				// This should be the token to be substituted in for the equate symbol
				token.SetLabileFlagBit(NxsToken::parentheticalToken);
				token.SetLabileFlagBit(NxsToken::curlyBracketedToken);
				token.GetNextToken();
				NxsString v = token.GetToken();

				// Add the new equate association to the equates list
				equates[k] = v;
				}

			standardDataTypeAssumed = true;
			}
		else if (token.Equals("LABELS"))
			{
			labels = true;
			standardDataTypeAssumed = true;
			}
		else if (token.Equals("NOLABELS"))
			{
			labels = false;
			standardDataTypeAssumed = true;
			}
		else if (token.Equals(";"))
			{
			break;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called from HandleMatrix function to read in the next state. Returns true if next token encountered is a comma, 
|	false otherwise. A comma signals the end of data for the current taxon in an UNALIGNED block.
*/
bool NxsUnalignedBlock::HandleNextState(
  NxsToken & token,			/* is the token used to read from `in' */
  unsigned i,				/* is the row in range [0..ntax) (used for error reporting only) */
  UnalignedVect & new_row)	/* is the container for storing new state */
	{
	// This should be the next state for taxon i
	token.SetLabileFlagBit(NxsToken::parentheticalToken);
	token.SetLabileFlagBit(NxsToken::curlyBracketedToken);
	token.SetLabileFlagBit(NxsToken::singleCharacterToken);
	token.GetNextToken();

	// Make sure we didn't run out of file
	if (token.AtEOF())
		{
		errormsg = "Unexpected end of file encountered";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	// If we didn't run out of file, there is no reason why we should have a zero-length token on our hands
	assert(token.GetTokenLength() > 0);

	// See if any equate macros apply
	NxsString skey = NxsString(token.GetToken(true)); // equates should always respect case
	NxsStringMap::iterator p = equates.find(skey);
	if (p != equates.end())
		{
		NxsString sval = (*p).second;
		token.ReplaceToken(sval.c_str());
		}

	// Handle case of single-character state symbol or comma/semicolon terminator
	if (token.GetTokenLength() == 1)
		{
		char ch = token.GetToken()[0];

		// Check for end of data signal (comma)
		if (ch == ',')
			{
			return true;
			}

		// Check for end of data signal (semicolon)
		if (ch == ';')
			{
			if (i + 1 != ntax)
				{
				errormsg = "Semicolon terminating the MATRIX statement of the UNALIGNED block was found while reading data for taxon ";
				errormsg += (i+1);
				errormsg += " of ";
				errormsg += ntax;
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			return true;
			}

		// Check for missing data symbol
		if (ch == missing)
			{
			NxsDiscreteDatum d;
			// d.states could be one of these:
			//	     ?              NULL          <-- this one applies here
			//	     G              [1][2]
			//	(AG) polymorphic    [2][0][2][1]
			//	{AG} ambiguous      [2][0][2][0]
			d.states = NULL;
			new_row.push_back(d);
			}

		// Look up the position of this state in the symbols array
		else
			{
			int p = PositionInSymbols(ch);
			if (p < 0)
				{
				errormsg = "State specified (";
				errormsg += token.GetToken();
				errormsg += ") for taxon ";
				errormsg += (i+1);
				errormsg += ", character ";
				errormsg += (unsigned)new_row.size() + 1;
				errormsg += ", not found in list of valid symbols";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			NxsDiscreteDatum d;
			// d.states could be one of these:
			//	     ?              NULL
			//	     G              [1][2]       <-- this one applies here
			//	(AG) polymorphic    [2][0][2][1]
			//	{AG} ambiguous      [2][0][2][0]
			d.states = new unsigned[2];
			d.states[0] = 1;
			d.states[1] = (unsigned)p;
			new_row.push_back(d);
			}
		}	// if (token.GetTokenLength() == 1)

	// Handle case of state sets
	else
		{
		// Token should be in one of the following forms: LEFT_SQUIGGLYacgRIGHT_SQUIGGLY LEFT_SQUIGGLYa~gRIGHT_SQUIGGLY LEFT_SQUIGGLYa c gRIGHT_SQUIGGLY (acg) (a~g) (a c g) 
		//
		NxsString t = token.GetToken();
		unsigned tlen = (unsigned)t.size();
		unsigned poly = (t[0] == '(');

		// Perform some sanity checks in debug builds
		assert(poly || t[0] == '{');
		assert((poly && t[tlen-1] == ')') || (!poly && t[tlen-1] == '}'));

		unsigned first_nonblank = 1;
		while (t[first_nonblank] == ' ' || t[first_nonblank] == '\t')
			first_nonblank++;

		unsigned last_nonblank = tlen - 2;
		while (t[last_nonblank] == ' ' || t[last_nonblank] == '\t')
			last_nonblank--;

		if (t[first_nonblank] == '~' || t[last_nonblank] == '~')
			{
			errormsg = token.GetToken();
			errormsg += " does not represent a valid range of states";
			throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

		// Starting right after the initial '(' or '{', record all states found in the vector stateset
		// After reaching the terminating ')' or '}', create the NxsDiscreteDatum and store it in new_row
		// The vector stateset is needed because we need to know how many states are present before allocating
		// memory for the NxsDiscreteDatum states member.
		unsigned k = 1;
		bool tildeFound = false;
		int p_before_tilde = -1;
		int p_after_tilde = -1;
		std::vector<unsigned> stateset; // used to store states until we're sure we've read them all
		for (;;)
			{
			if (t[k] == ')' || t[k] == '}')
				break;

			if (t[k] == ' ' || t[k] == '\t')
				{
				k++;
				continue;
				}

			// t[k] should be either '~' or one of the state symbols
			if (t[k] == '~')
				{
				tildeFound = true;
				}
			else
				{
				// Add state symbol and record if it is the first or last one in case we encounter a tilde
				if (tildeFound)
					{
					// Check to make sure t[k] is a valid symbol and that its position in symbols list
					// is after the symbol that appeared before the tilde
					p_after_tilde = PositionInSymbols(t[k]);
					if (p_after_tilde < 0)
						{
						errormsg = "State specified (";
						errormsg += t[k];
						errormsg += ") for taxon ";
						errormsg += (i+1);
						errormsg += ", character ";
						errormsg += (unsigned)new_row.size() + 1;
						errormsg += ", not found in list of valid symbols";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}
					if (p_after_tilde < p_before_tilde)
						{
						errormsg = "Final state (";
						errormsg += t[k];
						errormsg += ") comes before starting state (";
						errormsg += symbols[p_before_tilde];
						errormsg += ") in stateset range defined for taxon ";
						errormsg += (i+1);
						errormsg += ", character ";
						errormsg += (unsigned)new_row.size() + 1;
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}

					// Add all states from p_before_tilde + 1 to p_after_tilde then set tildeFound to false again
					for (unsigned pos = p_before_tilde + 1; pos <= (unsigned)p_after_tilde; ++pos)
						{
						stateset.push_back(pos);
						}

					tildeFound = false;
					p_before_tilde = -1;
					p_after_tilde = -1;
					}
				else
					{
					p_before_tilde = PositionInSymbols(t[k]);
					if (p_before_tilde < 0)
						{
						errormsg = "State specified (";
						errormsg += t[k];
						errormsg += ") for taxon ";
						errormsg += (i+1);
						errormsg += ", character ";
						errormsg += (unsigned)new_row.size() + 1;
						errormsg += ", not found in list of valid symbols";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}
					stateset.push_back(p_before_tilde);
					}
				} // if (t[k] == '~') ... else ... loop
			k++;
			} // for (;;) loop

		NxsDiscreteDatum d;
		// d.states could be one of these:
		//	     ?              NULL
		//	     G              [1][2]
		//	(AG) polymorphic    [2][0][2][1] <-- this one applies if poly true
		//	{AG} ambiguous      [2][0][2][0] <-- this one applies if poly false
		unsigned num_states = (unsigned)stateset.size();
		assert(num_states > 1);
		std::sort(stateset.begin(), stateset.end());
		d.states = new unsigned[num_states + 2];
		d.states[0] = num_states;
		k = 1;
		for (std::vector<unsigned>::const_iterator v = stateset.begin(); v != stateset.end(); ++v)
			{
			d.states[k++] = (*v);
			}
		d.states[num_states + 1] = (poly ? 1 : 0);
		new_row.push_back(d);
		}	// if (token.GetTokenLength() == 1) ... else ...

	// Normal situation is to return false, indicating that we have not yet reached the end of the 
	// data for the current taxon.
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when MATRIX command needs to be parsed from within the UNALIGNED block. Deals with everything after the 
|	token MATRIX up to and including the semicolon that terminates the MATRIX command.
*/
void NxsUnalignedBlock::HandleMatrix(
  NxsToken & token)	/* is the token used to read from `in' */
	{
	unsigned i;

	if (ntax == 0)
		{
		errormsg = "Must precede ";
		errormsg += id;
		errormsg += " block with a TAXA block or specify NEWTAXA and NTAX in the DIMENSIONS command";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	if (ntaxTotal == 0)
		ntaxTotal = taxa->GetNumTaxonLabels();

	// Start over with an empty umatrix, which is a vector of vectors of NxsDiscreteDatum objects
	umatrix.clear();

	// Allocate memory for (and initialize) the array activeTaxon. All taxa are initially active.
	//
	activeTaxon = new bool[ntax];
	for (i = 0; i < ntax; i++)
		activeTaxon[i] = true;

	// The value of ntaxTotal equals the number of taxa specified in the TAXA block, whereas ntax equals the number of
	// taxa specified in the DIMENSIONS command of the UNALIGNED block. These two numbers will be identical unless 
	// some taxa were left out of the MATRIX command of the UNALIGNED block, in which case ntax < ntaxTotal.
	if (taxonPos != NULL)
		delete [] taxonPos;
	taxonPos = new unsigned[ntaxTotal];
	for (i = 0; i < ntaxTotal; i++)
		taxonPos[i] = UINT_MAX;

	// Call an implementation function to handle most of the work
	HandleMatrixImpl(token);

	// If we've gotten this far, presumably it is safe to tell the ASSUMPTIONS block that were ready to take on the 
	// responsibility of being the current character-containing block (to be consulted if taxa are deleted or restored)
	//@POL what if someone tries to include/exclude characters? Need to find a way to prevent that if UNALIGNED block
	// is the current character-containing block
	//@POL oops, SetCallback requires a reference to a NxsCharactersBlock
	//assumptionsBlock->SetCallback(this);

	// The terminating semicolon at the end of the matrix command has already been read
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the HandleMatrix function, this function does most of the work of reading in the contents of the MATRIX 
|	command.
*/
void NxsUnalignedBlock::HandleMatrixImpl(
  NxsToken & token)	/* is the token used to read from `in' */
	{
	assert(taxonPos != NULL);

	unsigned i;
	//************************************************
	//******** Beginning of loop through taxa ********
	//************************************************
	for (i = 0; i < ntax; i++)
		{
		if (labels)
			{
			// This should be the taxon label
			token.GetNextToken();

			if (newtaxa)
				{
				// This section:
				// - labels provided
				// - previous TAXA block, if any, will be obliterated and replaced with information herein

				// Check for duplicate taxon names
				if (taxa->IsAlreadyDefined(token.GetToken()))
					{
					errormsg = "Data for this taxon (";
					errormsg += token.GetToken();
					errormsg += ") has already been saved";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}

				// Taxon labels have not yet been stored in the taxa block because taxa->Reset() was called
				// when the NEWTAXA subcommand was processed, thus obliterating any previous taxon information.
				taxa->AddTaxonLabel(token.GetToken());

				// Order of occurrence of taxa in TAXA block will be the same as umatrix because a new TAXA block
				// is being built up from the taxon names in this UNALIGNED block. The taxonPos vector
				// lets us look up which row in umatrix corresponds to any particular taxon in the taxa block.
				// That is, taxonPos[taxa_block_index] = umatrix_row.
				taxonPos[i] = i;
				}	// if (newtaxa)

			else
				{
				// This section:
				// - labels provided
				// - TAXA block provided or has been created already

				// Cannot assume taxon in same position in taxa block. Set up taxonPos array so that we can look up
				// the correct row in umatrix for any given taxon. That is, taxonPos[taxa_block_index] = umatrix_row.
				unsigned positionInTaxaBlock;
				try
					{
					positionInTaxaBlock = taxa->FindTaxon(token.GetToken());
					}
				catch(NxsTaxaBlock::NxsX_NoSuchTaxon)
					{
					if (token.Equals(";") && i == 0)
						{
						errormsg = "Unexpected semicolon (';') found while reading the UNALIGNED block MATRIX command";
						}
					else
						{
						errormsg = "Could not find taxon named ";
						errormsg += token.GetToken();
						errormsg += " among stored taxon labels";
						}
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}

				// Make sure user has not duplicated data for a single taxon
				if (taxonPos[positionInTaxaBlock] != UINT_MAX)
					{
					errormsg = "Data for this taxon (";
					errormsg += token.GetToken();
					errormsg += ") has already been saved";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}

				taxonPos[positionInTaxaBlock] = i;
				}	// if (newtaxa) ... else
			}	// if (labels)

		else
			{
			// No taxon labels provided in the MATRIX command, so we must assume that the positions of taxa are
			// the same here and in the taxa block
			taxonPos[i] = i;
			}	// if (labels) ... else

		//******************************************************
		//******** Beginning of loop through characters ********
		//******************************************************
		UnalignedVect new_row;
		bool done = false;
		while (!done)
			{
			done = HandleNextState(token, i, new_row);
			}
		umatrix.push_back(new_row);
		} // for (i = 0; i < ntax; i++)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when TAXLABELS command needs to be parsed from within the UNALIGNED block. Deals with everything after the 
|	token TAXLABELS up to and including the semicolon that terminates the TAXLABELS command.
*/
void NxsUnalignedBlock::HandleTaxlabels(
  NxsToken & token)	/* the token used to read from `in' */
	{
	if (!newtaxa)
		{
		errormsg = "NEWTAXA must have been specified in DIMENSIONS command to use the TAXLABELS command in a ";
		errormsg += id;
		errormsg += " block";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	for (;;)
		{
		token.GetNextToken();

		// Token should either be ';' or the name of a taxon
		if (token.Equals(";"))
			{
			break;
			}
		else
			{
			// Check to make sure user is not trying to read in more taxon labels than there are taxa
			if (taxa->GetNumTaxonLabels() > ntaxTotal)
				{
				errormsg = "Number of taxon labels exceeds NTAX specified in DIMENSIONS command";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			taxa->AddTaxonLabel(token.GetToken());
			}
		}

	// OPEN ISSUE: Some may object to setting newtaxa to false here, because then the fact that new taxa were 
	// specified in this UNALIGNED block rather than in a preceding TAXA block is lost. This will only be 
	// important if we wish to recreate the original data file, which I don't anticipate anyone doing with
	// this code (too difficult to remember all comments, the order of blocks in the file, etc.)

	newtaxa = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns position of `ch' in `symbols' array. The value of `respectingCase' is used to determine whether the search 
|	should be case sensitive or not. Assumes `symbols' is non-NULL. Returns UINT_MAX if `ch' is not found in `symbols'.
*/
unsigned NxsUnalignedBlock::PositionInSymbols(
  char ch)	/* the symbol character to search for */
	{
	assert(symbols != NULL);
	unsigned symbolsLength = (unsigned)strlen(symbols);
	bool found = false;
	unsigned i;
	for (i = 0; i < symbolsLength; i++)
		{
		char char_in_symbols	= (respectingCase	? symbols[i]	: (char)toupper(symbols[i]));
		char char_in_question	= (respectingCase	? ch			: (char)toupper(ch));
		if (char_in_symbols != char_in_question)
			continue;
		found = true;
		break;
		}
	return (found ? i : UINT_MAX);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the ability to read everything following the block name (which is read by the NxsReader 
|	object) to the END or ENDBLOCK statement. Characters are read from the input stream `in'. Overrides the abstract 
|	virtual function in the base class.
*/
void NxsUnalignedBlock::Read(
  NxsToken & token)	/* is the token used to read from `in' */
	{
	isEmpty = false;
	isUserSupplied = true;

	// This should be the semicolon after the block name
	token.GetNextToken(); 
	if (!token.Equals(";"))
		{
		errormsg = "Expecting ';' after ";
		errormsg += id;
		errormsg += " block name, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	ntax = taxa->GetNumTaxonLabels();

	for (;;)
		{
		token.GetNextToken();

		if (token.Equals("DIMENSIONS"))
			{
			HandleDimensions(token);
			}
		else if (token.Equals("FORMAT"))
			{
			HandleFormat(token);
			}
		else if (token.Equals("TAXLABELS"))
			{
			HandleTaxlabels(token);
			}
		else if (token.Equals("MATRIX"))
			{
			HandleMatrix(token);
			}
		else if (token.Equals("END"))
			{
			HandleEndblock(token);
			break;
			}
		else if (token.Equals("ENDBLOCK"))
			{
			HandleEndblock(token);
			break;
			}
		else
			{
			SkippingCommand(token.GetToken());

			do
				{
				token.GetNextToken();
				}
			while (!token.AtEOF() && !token.Equals(";"));

			if (token.AtEOF())
				{
				errormsg = "Unexpected end of file encountered";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			}	// else
		}	// for (;;)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function outputs a brief report of the contents of this UNALIGNED block. Overrides the abstract virtual 
|	function in the base class.
*/
void NxsUnalignedBlock::Report(
  ostream & out)	/* is the output stream to which to write the report */
	{
	out << endl;
	out << id << " block contains ";
	if (ntax == 0)
		out << "no taxa";
	else if (ntax == 1)
		out << "one taxon";
	else
		out << ntax << " taxa";
	out << endl;

	out << "  Data type is \"" << this->GetDatatypeName() << "\"" << endl;
	
	if (respectingCase)
		out << "  Respecting case" << endl;
	else
		out << "  Ignoring case" << endl;

	if (labels)
		out << "  Taxon labels were provided on left side of matrix" << endl;
	else
		out << "  No taxon labels were provided on left side of matrix" << endl;

	out << "  Missing data symbol is '" << missing << '\'' << endl;
	out << "  Valid symbols are: " << symbols << endl;

	int numEquateMacros = (int)equates.size();
	if (numEquateMacros > 0)
		{
		out << "  Equate macros in effect:" << endl;
		typedef NxsStringMap::const_iterator CI;
		for (CI i = equates.begin(); i != equates.end(); ++i)
			{
			out << "    " << (*i).first << " = " << (*i).second << endl;
			}
		}
	else
		out << "  No equate macros have been defined" << endl;

	out << "  The following taxa have been deleted:" << endl;
	unsigned nx = 0;
	for (unsigned k = 0; k < ntax; k++)
		{
		if (activeTaxon[k])
			continue;
		out << "    " << (k+1) << endl;
		nx++;
		}

	if (nx == 0)
		out << "    (no taxa deleted)" << endl;

	out << "  Data matrix:" << endl;
	DebugShowMatrix(out, "    ");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Writes out the information in this block in NEXUS format to the specified ostream.
*/
void NxsUnalignedBlock::WriteAsNexus(
  std::ostream & out)	/* is the output stream on which to write */
  const 
	{
	out << "BEGIN UNALIGNED;\n\t";
	this->WriteTitleCommand(out);
	out << '\t';
	if (this->newtaxa)
		{
		out << "DIMENSIONS NewTaxa NTax=" << this->ntax << ";\n\t";
		}
	this->WriteFormatCommand(out);
	
	if (this->newtaxa)
		{
		out << "\n\t";
		this->taxa->WriteTaxLabelsCommand(out);
		}

	out << "\n\t";
	this->WriteMatrixCommand(out);

	out << "\nEND;\n";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Writes out the information in the MATRIX command in NEXUS format to the specified ostream.
*/
void NxsUnalignedBlock::WriteMatrixCommand(
  std::ostream & out)	/* is the output stream on which to print the matrix */
  const
	{
	assert(taxonPos != NULL);

	unsigned width = taxa->GetMaxTaxonLabelLength();
	out << "Matrix";
	
	//std::vector<unsigned> origIndexVec = this->GetOrigMatrixIndicesToWrite();
		
	for (unsigned i = 0; i < this->ntaxTotal; ++i)
		{
		// Grab taxon name from taxa block. Taxa may not have been presented in the matrix in the same order as they
		// were stored in the taxa block, so use taxonPos array to determine which name goes with which matrix row.
		// That is, taxonPos[taxa_block_index] = umatrix_row. If taxonPos[i] = UINT_MAX, it means that no data exists
		// for the taxon with index i in the TAXA block.
		if (this->taxonPos[i] != UINT_MAX)
			{
			out << "\n\t\t";
			NxsString nm = taxa->GetTaxonLabel(i);
			std::string s = nm.c_str();
			const std::string currTaxonLabel = NxsToken::EscapeString(taxa->GetTaxonLabel(i));
			out << currTaxonLabel;

			// Print out enough spaces to even up the left edge of the matrix output
			unsigned currTaxonLabelLen = (unsigned)currTaxonLabel.size();
			unsigned diff = width - currTaxonLabelLen;
			for (unsigned k = 0; k < diff + 5; k++)
				out << ' ';

			// Print out states for all characters store for taxon i
			// Remember taxon i in the TAXA block corresponds to row taxonPos[i] in umatrix
			unsigned row = taxonPos[i];
			UnalignedVect::const_iterator j = umatrix[row].begin();
			for (; j != umatrix[row].end(); ++j)
				{
				out << FormatState(*j);
				}
			}
		}
	out << ';';
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Writes out the information in the FORMAT command in NEXUS format to the specified ostream.
*/
void NxsUnalignedBlock::WriteFormatCommand(std::ostream &out) const
	{
	// Output name of command
	out << "FORMAT";
	
	// Output DATATYPE statement
	out << " Datatype=" << this->GetDatatypeName();

	// Output RESPECTCASE statement
	if (this->respectingCase)
		out << " RespectCase"; 

	// Output MISSING statement
	if (this->missing != '?')
		{
		out << " Missing=";
		out << this->missing;
		}

	// Output SYMBOLS statement
	// Determine length of symbols array based on current datatype
	unsigned numDefStates = 4;
	if (this->datatype == NxsUnalignedBlock::protein)
		numDefStates = 21;
	else if (this->datatype == NxsUnalignedBlock::standard)	
		numDefStates = 0;

	// Print out a symbols statement if non-default symbols have been added
	// Only include symbols in the symbols statement that are not default symbols
	// for the selected datatype
	unsigned nSym = (unsigned)strlen(this->symbols);
	if (nSym > numDefStates)
		out << " Symbols=\"" << (this->symbols + numDefStates) <<"\"";

	// Output EQUATES statement
	// Print out an equates statement if any equates beyond the default ones
	// for the selected datatype have been defined. The following code iterates
	// through all defined equates, copying those not found in the default equates
	// map into a new map. This new map is printed out in the form of an equates
	// statement if it is non-empty.
	const NxsStringMap defEquates = this->GetDefaultEquates();
	NxsStringMap toWrite;
	const NxsStringMap::const_iterator notFound = defEquates.end();
	NxsStringMap::const_iterator inDefEquates;
	for (NxsStringMap::const_iterator i = equates.begin(); i != equates.end(); ++i)
		{
		const NxsString key =  (*i).first;
		const NxsString val =  (*i).second;
		inDefEquates = defEquates.find(key);
		if (inDefEquates == notFound || inDefEquates->second != val)
			toWrite[key] = val;
		}
	if (toWrite.size() > 0)
		{
		out << " Equate=\"";		
		for (NxsStringMap::const_iterator j = toWrite.begin(); j != toWrite.end(); ++j)
			out << ' ' << j->first << '=' << j->second;
		out <<"\"";
		}

	// Output LABELS statement
	if (labels)
		out << " LABELS";
	else
		out << " NOLABELS";

	// Output terminating semicolon
	out << ";";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates and returns a NxsStringMap containing all default equates for the current datatype. For dna, rna and 
|	nucleotide datatypes, equate macros are created for the standard ambiguity codes R, Y, M, K, S, W, H, B, V, D, N and
|	X. For proteins, only two equate macros are created, for B and Z.
*/
NxsStringMap NxsUnalignedBlock::GetDefaultEquates() const
	{
	NxsStringMap defEquates;
	const DataTypesEnum dt = this->datatype;
	if (dt == NxsUnalignedBlock::dna || dt == NxsUnalignedBlock::rna || dt == NxsUnalignedBlock::nucleotide)
		{
		//@POL shouldn't the rna datatype involve Us rather than Ts?
		defEquates[ NxsString("R") ] = NxsString("{AG}");
		defEquates[ NxsString("Y") ] = NxsString("{CT}");
		defEquates[ NxsString("M") ] = NxsString("{AC}");
		defEquates[ NxsString("K") ] = NxsString("{GT}");
		defEquates[ NxsString("S") ] = NxsString("{CG}");
		defEquates[ NxsString("W") ] = NxsString("{AT}");
		defEquates[ NxsString("H") ] = NxsString("{ACT}");
		defEquates[ NxsString("B") ] = NxsString("{CGT}");
		defEquates[ NxsString("V") ] = NxsString("{ACG}");
		defEquates[ NxsString("D") ] = NxsString("{AGT}");
		defEquates[ NxsString("N") ] = NxsString("{ACGT}");
		defEquates[ NxsString("X") ] = NxsString("{ACGT}");
		}
	else if (dt == NxsUnalignedBlock::protein)
		{
		defEquates[ NxsString("B") ] = NxsString("{DN}");
		defEquates[ NxsString("Z") ] = NxsString("{EQ}");
		}
	return defEquates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a C-string containing the name of the current datatype.
*/
const char * NxsUnalignedBlock::GetDatatypeName() const
	{
	switch(datatype)
		{
		case NxsUnalignedBlock::dna:
			return "DNA";
		case NxsUnalignedBlock::rna:
			return "RNA";
		case NxsUnalignedBlock::nucleotide:
			return "Nucleotide";
		case NxsUnalignedBlock::protein:
			return "Protein";
		default:
			return "Standard";
		}
	}
