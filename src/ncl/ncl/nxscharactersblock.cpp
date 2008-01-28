//	Copyright (C) 1999-2003 Paul O. Lewis
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
/**
 * This file includes contributions by Brian O'Meara. August 2005.
 * These changes include the ability to parse continuous data types.
 */

#include <iomanip>
#include "ncl.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `id' to "CHARACTERS", `taxa' to `tb', `assumptionsBlock' to `ab', `ntax', `ntaxTotal', `nchar' and 
|	`ncharTotal' to 0, `newchar' to true, `newtaxa', `interleaving', `transposing', `respectingCase', `tokens' and 
|	`formerly_datablock' to false, `datatype' to `NxsCharactersBlock::standard', `missing' to '?', `gap' and `matchchar'
|	to '\0', and `matrix', `charPos', `taxonPos', `activeTaxon', and `activeChar' to NULL. The ResetSymbols member 
|	function is called to reset the `symbols' data member. Assumes that `tb' and `ab' point to valid NxsTaxaBlock and 
|	NxsAssumptionsBlock objects, respectively.
*/
NxsCharactersBlock::NxsCharactersBlock(
  NxsTaxaBlock *tb,			/* the taxa block object to consult for taxon labels */
  NxsAssumptionsBlock *ab)	/* the assumptions block object to consult for exclusion sets */
  : NxsBlock()
	{
	assert(tb != NULL);
	assert(ab != NULL);

	taxa				= tb;
	assumptionsBlock	= ab;
	id					= "CHARACTERS";

	// These need to be initialized to NULL so Reset member function will not try to delete them
	//
	matrix				= NULL;
	charPos				= NULL;
	taxonPos			= NULL;
	activeTaxon			= NULL;
	activeChar			= NULL;
	symbols				= NULL;

	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes any memory allocated to the arrays `symbols', `charPos', `taxonPos', `activeChar', and `activeTaxon'. 
|	Flushes the containers `charLabels', `eliminated', and `deleted'. Also deletes memory allocated to `matrix'.
*/
NxsCharactersBlock::~NxsCharactersBlock()
	{
	Reset();

	if (symbols != NULL)
		delete [] symbols;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes (i.e., excludes from further analyses) taxa whose indices are contained in the set `delset'. The taxon 
|	indices refer to original taxon indices, not current indices (originals will equal current ones if number of taxa 
|	in TAXA block equals number of taxa in MATRIX command). Returns the number of taxa actually deleted (some may have 
|	already been deleted)
*/
unsigned NxsCharactersBlock::ApplyDelset(
  NxsUnsignedSet &delset)	/* set of taxon indices to delete in range [0..`ntaxTotal') */
	{
	assert(activeTaxon != NULL);
	assert(taxonPos != NULL);

	unsigned num_deleted = 0;
	unsigned k;

	NxsUnsignedSet::const_iterator i;
	for (i = delset.begin(); i != delset.end(); i++)
		{
		k = taxonPos[*i];
		if (k == UINT_MAX)
			continue;

		// k equal to UINT_MAX means data was supplied for
		// this taxon and therefore it can be deleted
		//
		if (activeTaxon[k] == true)
			num_deleted++;
		activeTaxon[k] = false;
		}

	return num_deleted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Excludes characters whose indices are contained in the set `exset'. The indices supplied should refer to the 
|	original character indices, not current character indices. Returns number of characters actually excluded (some 
|	may have already been excluded).
*/
unsigned NxsCharactersBlock::ApplyExset(
  NxsUnsignedSet &exset)	/* set of character indices to exclude in range [0..`ncharTotal') */
	{
	assert(activeChar != NULL);
	assert(charPos != NULL);

	int num_excluded = 0;
	unsigned k;

	NxsUnsignedSet::const_iterator i;
	for (i = exset.begin(); i != exset.end(); i++)
		{
		k = charPos[*i];
		if (k == UINT_MAX)
			continue;

		// k equal to UINT_MAX means character was not eliminated
		// and therefore can be excluded
		//
		if (activeChar[k] == true)
			num_excluded++;
		activeChar[k] = false;
		}

	return num_excluded;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Includes characters whose indices are contained in the set `inset'. The indices supplied should refer to the 
|	original character indices, not current character indices.
*/
unsigned NxsCharactersBlock::ApplyIncludeset(
  NxsUnsignedSet &inset)	/* set of character indices to include in range [0..`ncharTotal') */
	{
	assert(activeChar != NULL);
	assert(charPos != NULL);

	unsigned num_included = 0;
	unsigned k;

	NxsUnsignedSet::const_iterator i;
	for (i = inset.begin(); i != inset.end(); i++)
		{
		k = charPos[*i];
		if (k == UINT_MAX)
			continue;

		// k equal to UINT_MAX means character was not eliminated
		// and therefore can be excluded
		//
		if (activeChar[k] == false)
			num_included++;
		activeChar[k] = true;
		}

	return num_included;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores (i.e., includes in further analyses) taxa whose indices are contained in the set `restoreset'. The taxon 
|	indices refer to original taxon indices, not current indices (originals will equal current ones if number of taxa 
|	in TAXA block equals number of taxa in MATRIX command).
*/
unsigned NxsCharactersBlock::ApplyRestoreset(
  NxsUnsignedSet &restoreset)	/* set of taxon indices to restore in range [0..`ntaxTotal') */
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

		// k equal to UINT_MAX means data was supplied for
		// this taxon and therefore it can be restored
		//
		if (activeTaxon[k] == false)
			num_restored++;
		activeTaxon[k] = true;
		}

	return num_restored;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use to allocate memory for (and initialize) `charPos' array, which keeps track of the original character index in 
|	cases where characters have been eliminated. This function is called by HandleEliminate in response to encountering 
|	an ELIMINATE command in the data file, and this is probably the only place where BuildCharPosArray should be called 
|	with `check_eliminated' true. BuildCharPosArray is also called in HandleMatrix, HandleCharstatelabels, 
|	HandleStatelabels, and HandleCharlabels.
*/
void NxsCharactersBlock::BuildCharPosArray(
  bool check_eliminated)	/* if true, eliminated set has something in it and should be consulted (default is false) */
	{
	assert(charPos == NULL);

	charPos = new unsigned[ncharTotal];

	unsigned k = 0;
	for (unsigned j = 0; j < ncharTotal; j++)
		{
		if (check_eliminated && IsEliminated(j))
			charPos[j] = UINT_MAX;
		else
			charPos[j] = k++;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Converts a character label to a 1-offset number corresponding to the character's position within `charLabels'. This
|	method overrides the virtual function of the same name in the NxsBlock base class. If `s' is not a valid character 
|	label, returns the value 0.
*/
unsigned NxsCharactersBlock::CharLabelToNumber(
  NxsString s)	/* the character label to convert */
	{
	NxsStringVector::const_iterator iter = find(charLabels.begin(), charLabels.end(), s);

	unsigned k = 1;
	if (iter != charLabels.end())
		k += (unsigned)(iter - charLabels.begin());
	else
		k = 0;

	return k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Transfers all data from `other' to this object, leaving `other' completely empty. Used to convert a NxsDataBlock 
|	object to a NxsCharactersBlock object in programs where it is desirable to just have a NxsCharactersBlock for 
|	storage but also allow users to enter the information in the form of the deprecated NxsDataBlock. This function 
|	does not make a copy of such things as the data matrix, instead just transferring the pointer to that object from 
|	other to this. This is whay it was named Consume rather than CopyFrom.
*/
void NxsCharactersBlock::Consume(
  NxsCharactersBlock &other)	/* NxsCharactersBlock object from which to copy */
	{
	ntax				= other.ntax;
	ntaxTotal			= other.ntaxTotal;
	nchar				= other.nchar;
	ncharTotal			= other.ncharTotal;

	newtaxa				= other.newtaxa;
	newchar				= other.newchar;

	formerly_datablock	= true;
	respectingCase		= other.respectingCase;
	transposing			= other.transposing;
	interleaving		= other.interleaving;
	tokens				= other.tokens;
	labels				= other.labels;

	missing				= other.missing;
	gap					= other.gap;
	matchchar			= other.matchchar;

	datatype			= other.datatype;

	if (symbols != NULL)
		delete [] symbols;
	symbols				= other.symbols;
	other.symbols		= NULL;

	if (charPos != NULL)
		delete [] charPos;
	charPos				= other.charPos;
	other.charPos		= NULL;

	if (taxonPos != NULL)
	delete [] taxonPos;
	taxonPos			= other.taxonPos;
	other.taxonPos		= NULL;

	if (activeChar != NULL)
	delete [] activeChar;
	activeChar			= other.activeChar;
	other.activeChar	= NULL;

	if (activeTaxon != NULL)
	delete [] activeTaxon;
	activeTaxon			= other.activeTaxon;
	other.activeTaxon	= NULL;

	if (matrix != NULL)
	delete matrix;
	matrix				= other.matrix;
	other.matrix		= NULL;

	equates.clear();
	int size = (int)other.equates.size();
	if (size > 0)
		{
		NxsStringMap::const_iterator i;
		for (i = other.equates.begin(); i != other.equates.end(); ++i)
			equates[(*i).first] = (*i).second;
		other.equates.clear();
		}

	eliminated.clear();
	size = (int)eliminated.size();
	if (size > 0)
		{
		NxsUnsignedSet::const_iterator i;
		for (i = other.eliminated.begin(); i != other.eliminated.end(); i++)
			eliminated.insert(*i);
		other.eliminated.clear();
		}

	charLabels.clear();
	size = (int)charLabels.size();
	if (size > 0)
		{
		NxsStringVector::const_iterator i;
		for (i = other.charLabels.begin(); i != other.charLabels.end(); i++)
			charLabels.push_back((*i));
		other.charLabels.clear();
		}

	charStates.clear();
	size = (int)charStates.size();
	if (size > 0)
		{
		NxsStringVectorMap::const_iterator i;
		for (i = other.charStates.begin(); i != other.charStates.end(); i++)
			charStates[ (*i).first ] = (*i).second;
		other.charStates.clear();
		}

	isEmpty = false;
	isUserSupplied = other.isUserSupplied;
	items = other.items;
	continuousCharMatrix = other.continuousCharMatrix;
	statesFormat = other.statesFormat;
	other.Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides a dump of the contents of the `matrix' variable. Useful for testing whether data is being read as 
|	expected. If marginText is NULL, matrix output is placed flush left. If each line of output should be prefaced with 
|	a tab character, specify "\t" for `marginText'.
*/
void NxsCharactersBlock::DebugShowMatrix(
  ostream &out,			/* output stream on which to print matrix */
  bool use_matchchar,	/* if true, matchchar symbol used; otherwise, states shown for all taxa */
  char *marginText)		/* for printing first on each line */
	{
	assert(charPos != NULL);
	assert(taxonPos != NULL);

	unsigned i, k;
	unsigned width = taxa->GetMaxTaxonLabelLength();
	unsigned first_taxon = UINT_MAX;

	for (i = 0; i < ntaxTotal; i++)
		{
		// Grab taxon name from taxa block. Taxa may not have been presented in the matrix in the same order
		// as they were stored in the taxa block, so use taxonPos array to spit them out in the order they 
		// appeared in the TAXA command. If the taxonPos cell is UINT_MAX, then that means there is no row of
		// the data matrix corresponding to that taxon.
		//
		if (taxonPos[i] == UINT_MAX)
			continue;
		else
			{
			if (first_taxon == UINT_MAX)
				first_taxon = i;

			if (marginText != NULL)
				out << marginText;

			NxsString currTaxonLabel = taxa->GetTaxonLabel(taxonPos[i]);
				out << currTaxonLabel;

			// Print out enough spaces to even up the left edge of the matrix output
			//
			unsigned currTaxonLabelLen = (unsigned)currTaxonLabel.size();
			unsigned diff = width - currTaxonLabelLen;
			for (k = 0; k < diff+5; k++)
				out << ' ';
			}

		for (unsigned currChar = 0; currChar < ncharTotal; currChar++)
			{
			unsigned j = charPos[currChar];
			if (j == UINT_MAX)
				continue;
			ShowStateLabels(out, i, j, (use_matchchar ? first_taxon : UINT_MAX));
			}

		out << endl;
		}	// for (i = 0; i < ntaxTotal; i++)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the maximum observed number of states for any character. Note: this function is rather slow, as it must 
|	walk through each row of each column, adding the states encountered to a set,  then finally returning the size of 
|	the set. Thus, if this function is called often, it would be advisable to initialize an array using this function, 
|	then refer to the array subsequently. 
*/
unsigned NxsCharactersBlock::GetMaxObsNumStates()
	{
	unsigned max = 2;
	for (unsigned j = 0; j < nchar; j++)
		{
		unsigned ns = GetObsNumStates(j);
		if (ns <= max)
			continue;
		max = ns;
		}

	return max;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a count of the number of characters for which `activeChar' array reports true.
*/
unsigned NxsCharactersBlock::GetNumActiveChar()
	{
	unsigned num_active_char = 0;
	for (unsigned i = 0; i < nchar; i++)
		{
		if (activeChar[i])
			num_active_char++;
		}

	return num_active_char;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a count of the number of taxa for which `activeTaxon' array reports true.
*/
unsigned NxsCharactersBlock::GetNumActiveTaxa()
	{
	unsigned num_active_taxa = 0;
	for (unsigned i = 0; i < ntax; i++)
		{
		if (activeTaxon[i])
			num_active_taxa++;
		}

	return num_active_taxa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the original character index in the range [0..`ncharTotal'). Will be equal to `j' unless some characters 
|	were eliminated.
*/
unsigned NxsCharactersBlock::GetOrigCharIndex(
  unsigned j) const /* the character in range [0..`nchar') */
	{
	unsigned k = j;
	while (k < ncharTotal && charPos[k] < j)
		k++;

	assert(k < ncharTotal);
	return k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the original taxon index in the range [0..`ntaxTotal'). Will be equal to `i' unless data was not provided 
|	for some taxa listed in a preceding TAXA block.
*/
unsigned NxsCharactersBlock::GetOrigTaxonIndex(
  unsigned i)	/* the taxon in range [0..`ntax') */
	{
	assert(taxonPos != NULL);

	unsigned k = i;
	while (k < ntaxTotal && taxonPos[k] < i)
		k++;

	assert(k < ntaxTotal);
	return k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns label for character state `j' at character `i', if a label has been specified. If no label was specified, 
|	returns string containing a single blank (i.e., " ").
*/
NxsString NxsCharactersBlock::GetStateLabel(
  unsigned i,	/* the locus in range [0..`nchar') */
  unsigned j)	/* the 0-offset index of the state of interest */
	{
	NxsString s = " ";
	NxsStringVectorMap::const_iterator cib = charStates.find(i);
	if (cib != charStates.end() && static_cast<unsigned>(j) < (*cib).second.size())
		{
		s = (*cib).second[j];
		}

	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if character number `origCharIndex' was eliminated, false otherwise. Returns false immediately if 
|	`eliminated' set is empty.
*/
bool NxsCharactersBlock::IsEliminated(
  unsigned origCharIndex)	/* the character in question */
	{
	// Note: it is tempting to try to streamline this method by just looking up character j in charPos to see if it
	// has been eliminated, but this temptation should be resisted because this function is used in setting up
	// charPos in the first place!

	if (eliminated.empty())
		return false;

	NxsUnsignedSet::const_iterator found = eliminated.find(origCharIndex);
	if (found == eliminated.end())
		return false;

	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `ch' can be found in the `symbols' array. The value of `respectingCase' is used to determine 
|	whether or not the search should be case sensitive. Assumes `symbols' is non-NULL.
*/
bool NxsCharactersBlock::IsInSymbols(
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
|	Called when CHARLABELS command needs to be parsed from within the DIMENSIONS block. Deals with everything after 
|	the token CHARLABELS up to and including the semicolon that terminates the CHARLABELS command. If an ELIMINATE 
|	command has been processed, labels for eliminated characters will not be stored.
*/
void NxsCharactersBlock::HandleCharlabels(
  NxsToken &token)	/* the token used to read from `in' */
	{
	unsigned num_labels_read = 0;
	charLabels.clear();

	if (charPos == NULL)
		BuildCharPosArray();

	for (;;)
		{
		token.GetNextToken();

		// Token should either be ';' or the name of a character (an isolated '_' character is 
		// converted automatically by token.GetNextToken() into a space, which is then stored
		// as the character label)
		//
		if (token.Equals(";"))
			{
			break;
			}
		else
			{
			num_labels_read++;

			// Check to make sure user is not trying to read in more character labels than 
			// there are characters
			//
			if (num_labels_read > ncharTotal)
				GenerateNxsException(token, "Number of character labels exceeds NCHAR specified in DIMENSIONS command");

			if (!IsEliminated(num_labels_read - 1))
				charLabels.push_back(token.GetToken());
			}
		}

	newchar = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when CHARSTATELABELS command needs to be parsed from within the CHARACTERS block. Deals with everything 
|	after the token CHARSTATELABELS up to and including the semicolon that terminates the CHARSTATELABELS command. 
|	Resulting `charLabels' vector will store labels only for characters that have not been eliminated, and likewise for 
|	`charStates'. Specifically, `charStates[0]' refers to the vector of character state labels for the first 
|	non-eliminated character.
*/
void NxsCharactersBlock::HandleCharstatelabels(
  NxsToken &token)	/* the token used to read from `in' */
	{
	unsigned currChar = 0;
	bool semicolonFoundInInnerLoop = false;
	bool tokenAlreadyRead = false;
	bool save = true;

	charStates.clear();
	charLabels.clear();

	if (charPos == NULL)
		BuildCharPosArray();

	for (;;)
		{
		save = true;

		if (semicolonFoundInInnerLoop)
			break;

		if (tokenAlreadyRead)
			tokenAlreadyRead = false;
		else
			token.GetNextToken();

		if (token.Equals(";"))
			break;

		// Token should be the character number; create a new association
		//
		int sn = atoi(token.GetToken().c_str());
		unsigned n = (unsigned)sn;
		if (sn < 1 || n > ncharTotal || n <= currChar)
			{
			errormsg = "Invalid character number (";
			errormsg += token.GetToken();
			errormsg += ") found in CHARSTATELABELS command (either out of range or not interpretable as an integer)";
			throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

		// If n is not the next character after currChar, need to add some dummy
		// labels to charLabels list
		//
		while (n - currChar > 1) 
			{
			currChar++;
			if (!IsEliminated(currChar - 1))
				charLabels.push_back(" ");
			}

		// If n refers to a character that has been eliminated, go through the motions of
		// reading in the information but don't actually save any of it
		//
		currChar++;
		assert(n == currChar);
		if (IsEliminated(currChar-1))
			save = false;

		token.GetNextToken();

		// Token should be the character label
		//
		if (save) 
			charLabels.push_back(token.GetToken());

		token.GetNextToken();

		// Token should be a slash character if state labels were provided for this character; otherwise, 
		// token should be one of the following:
		// 1) the comma separating information for different characters, in which case we read in the 
		//	  next token (which should be the next character number)
		// 2) the semicolon indicating the end of the command
		//
		if (!token.Equals("/"))
			{
			if (!token.Equals(",") && !token.Equals(";"))
				{
				errormsg = "Expecting a comma or semicolon here, but found (";
				errormsg += token.GetToken();
				errormsg += ") instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			if (token.Equals(","))
				token.GetNextToken();
			tokenAlreadyRead = true;
			continue;
			}

		// Now create a new association for the character states list

		for (;;)
			{
			token.GetNextToken();

			if (token.Equals(";"))
				{
				semicolonFoundInInnerLoop = true;
				break;
				}

			if (token.Equals(","))
				{
				break;
				}

			if (save)
				{
				if (datatype == continuous)
					GenerateNxsException(token, "State Labels cannot be specified when the datatype is continuous");

				// Token should be a character state label; add it to the list
				//
				NxsString cslabel = token.GetToken();
				unsigned k = GetCharPos(n - 1);
				charStates[k].push_back(cslabel);
				}

			} // inner for (;;) loop (grabbing state labels for character n)
		} // outer for (;;) loop

	newchar = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when DIMENSIONS command needs to be parsed from within the CHARACTERS block. Deals with everything after 
|	the token DIMENSIONS up to and including the semicolon that terminates the DIMENSIONs command. `newtaxaLabel', 
|	`ntaxLabel' and `ncharLabel' are simply "NEWTAXA", "NTAX" and "NCHAR" for this class, but may be different for 
|	derived classes that use `newtaxa', `ntax' and `nchar' for other things (e.g., ntax is number of populations in 
|	an ALLELES block)
*/
void NxsCharactersBlock::HandleDimensions(
  NxsToken &token,			/* the token used to read from `in' */
  NxsString newtaxaLabel,	/* the label used in data file for `newtaxa' */
  NxsString ntaxLabel,		/* the label used in data file for `ntax' */
  NxsString ncharLabel)		/* the label used in data file for `nchar' */
	{
	for (;;)
		{
		token.GetNextToken();

		if (token.Equals(newtaxaLabel))
			{
			newtaxa = true;
			}
		else if (token.Equals(ntaxLabel)) 
			{
			DemandEquals(token, "DIMENSIONS command");
			ntax = DemandPositiveInt(token, ntaxLabel.c_str());
			if (newtaxa)
				ntaxTotal = ntax;
			else
				{
				ntaxTotal = taxa->GetNumTaxonLabels();
				if (ntaxTotal < ntax)
					{
					errormsg = ntaxLabel;
					errormsg += " in ";
					errormsg += id;
					errormsg += " block must be less than or equal to NTAX in TAXA block\nNote: one circumstance that can cause this error is \nforgetting to specify ";
					errormsg += ntaxLabel;
					errormsg += " in DIMENSIONS command when \na TAXA block has not been provided";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}
				}
			}
		else if (token.Equals(ncharLabel)) 
			{
			DemandEquals(token, "in DIMENSIONS command");
			nchar = DemandPositiveInt(token, ncharLabel.c_str());
			ncharTotal = nchar;
			}
		else if (token.Equals(";"))
			{
			break;
			}
		}

	if (newtaxa)
		taxa->Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when ELIMINATE command needs to be parsed from within the CHARACTERS block. Deals with everything after the 
|	token ELIMINATE up to and including the semicolon that terminates the ELIMINATE command. Any character numbers 
|	or ranges of character numbers specified are stored in the NxsUnsignedSet `eliminated', which remains empty until 
|	an ELIMINATE command is processed. Note that like all sets the character ranges are adjusted so that their offset 
|	is 0. For example, given "eliminate 4-7;" in the data file, the eliminate array would contain the values 3, 4, 5 
|	and 6 (not 4, 5, 6 and 7). It is assumed that the ELIMINATE command comes before character labels and/or character 
|	state labels have been specified; an error message is generated if the user attempts to use ELIMINATE after a 
|	CHARLABELS, CHARSTATELABELS, or STATELABELS command.
*/
void NxsCharactersBlock::HandleEliminate(
  NxsToken &token)	/* the token used to read from `in' */
	{
	// Construct an object of type NxsSetReader, then call its run function
	// to store the set in the eliminated set
	//
	NxsSetReader(token, ncharTotal, eliminated, *this, NxsSetReader::charset).Run();
	assert(eliminated.size() <= ncharTotal);
	nchar = ncharTotal - eliminated.size();

	if (nchar != ncharTotal && (charLabels.size() > 0 || charStates.size() > 0)) 
		GenerateNxsException(token, "The ELIMINATE command must appear before character\n(or character state) labels are specified");

	if (charPos != NULL) 
		GenerateNxsException(token, "Only one ELIMINATE command is allowed, and it must appear before the MATRIX command");

	BuildCharPosArray(true);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when the END or ENDBLOCK command needs to be parsed from within the CHARACTERS block. Does two things: 
|~
|	o checks to make sure the next token in the data file is a semicolon
|	o eliminates character labels and character state labels for characters that have been eliminated
|~
*/
void NxsCharactersBlock::HandleEndblock(
  NxsToken &token,		/* the token used to read from `in' */
  NxsString charToken)	/* */
	{
	DemandEndSemicolon(token, "END or ENDBLOCK");

	if (charLabels.empty() && !charStates.empty())
		{
		// Make up labels for characters since user has provided labels
		// for character states; that way, we know that charLabels
		// and charStates are either both empty or both full
		//
		for (unsigned k = 0; k < ncharTotal; k++)
			{
			NxsString nm = charToken;
			nm += " ";
			nm += (k+1);
			charLabels.push_back(nm.c_str());
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when FORMAT command needs to be parsed from within the DIMENSIONS block. Deals with everything after the 
|	token FORMAT up to and including the semicolon that terminates the FORMAT command.
*/
void NxsCharactersBlock::HandleFormat(
  NxsToken &token)	/* the token used to read from `in' */
	{
	bool standardDataTypeAssumed = false;
	bool ignoreCaseAssumed = false;

	for (;;)
		{
		token.GetNextToken();

		if (token.Equals("DATATYPE"))
			{
			DemandEquals(token, "after keyword DATATYPE");
			
			// This should be one of the following: STANDARD, DNA, RNA, NUCLEOTIDE, PROTEIN, or CONTINUOUS
			//
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
			else if (token.Equals("CONTINUOUS"))
				{
				datatype = continuous;
				statesFormat = INDIVIDUALS;
				items = std::vector<std::string>(1, std::string("AVERAGE"));
				}
			else
				{
				errormsg = token.GetToken();
				errormsg += " is not a valid DATATYPE within a ";
				errormsg += id;
				errormsg += " block";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			if (standardDataTypeAssumed && datatype != standard)
				GenerateNxsException(token, "DATATYPE must be specified first in FORMAT command");

			ResetSymbols();

			if (datatype == continuous)
				tokens = true;
			}

		else if (token.Equals("RESPECTCASE"))
			{
			if (ignoreCaseAssumed)
				GenerateNxsException(token, "RESPECTCASE must be specified before MISSING, GAP, SYMBOLS, and MATCHCHAR in FORMAT command");
			standardDataTypeAssumed = true;
			respectingCase = true;
			}

		else if (token.Equals("MISSING"))
			{
			DemandEquals(token, "after keyword MISSING");

			// This should be the missing data symbol (single character)
			//
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

		else if (token.Equals("GAP"))
			{
			DemandEquals(token, "after keyword GAP");

			// This should be the gap symbol (single character)
			//
			token.GetNextToken();

			if (token.GetTokenLength() != 1)
				{
				errormsg = "GAP symbol should be a single character, but ";
				errormsg += token.GetToken();
				errormsg += " was specified";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			else if (token.IsPunctuationToken() && !token.IsPlusMinusToken())
				{
				errormsg = "GAP symbol specified cannot be a punctuation token (";
				errormsg += token.GetToken();
				errormsg += " was specified)";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			else if (token.IsWhitespaceToken())
				{
				errormsg = "GAP symbol specified cannot be a whitespace character (";
				errormsg += token.GetToken();
				errormsg += " was specified)";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			gap = token.GetToken()[0];

			ignoreCaseAssumed = true;
			standardDataTypeAssumed = true;
			}

		else if (token.Equals("SYMBOLS"))
			{
			if (datatype == NxsCharactersBlock::continuous)
				GenerateNxsException(token, "SYMBOLS subcommand not allowed for DATATYPE=CONTINUOUS");

			int numDefStates;
			unsigned maxNewStates;
			switch(datatype)
				{
				case NxsCharactersBlock::dna:
				case NxsCharactersBlock::rna:
				case NxsCharactersBlock::nucleotide:
					numDefStates = 4;
					maxNewStates = NCL_MAX_STATES-4;
					break;

				case NxsCharactersBlock::protein:
					numDefStates = 21;
					maxNewStates = NCL_MAX_STATES-21;
					break;

				default:
					numDefStates = 0; // replace symbols list for standard datatype
					symbols[0] = '\0';
					maxNewStates = NCL_MAX_STATES;
				}

			DemandEquals(token, " after keyword SYMBOLS");

			// This should be the symbols list
			//
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
			//
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
			//
			strcpy(symbols + numDefStates, t.c_str());

			ignoreCaseAssumed = true;
			standardDataTypeAssumed = true;
			}

		else if (token.Equals("EQUATE"))
			{
			if (datatype == NxsCharactersBlock::continuous)
				GenerateNxsException(token, "EQUATE subcommand not allowed for DATATYPE=CONTINUOUS");

			DemandEquals(token, "after keyword EQUATE");
			
			// This should be a double-quote character
			//
			token.GetNextToken();

			if (!token.Equals("\""))
				{
				errormsg = "Expecting '\"' after keyword EQUATE but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			// Loop until second double-quote character is encountered
			//
			for (;;)
				{
				token.GetNextToken();
				if (token.Equals("\""))
					break;

				// If token is not a double-quote character, then it must be the equate symbol (i.e., the 
				// character to be replaced in the data matrix)
				//
				if (token.GetTokenLength() != 1)
					{
					errormsg = "Expecting single-character EQUATE symbol but found ";
					errormsg += token.GetToken();
					errormsg += " instead";
					throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
					}

				// Check for bad choice of equate symbol
				//
				NxsString t = token.GetToken();
				char ch = t[0];
				bool badEquateSymbol = false;

				// The character '^' cannot be an equate symbol
				//
				if (ch == '^')
					badEquateSymbol = true;

				// Equate symbols cannot be punctuation (except for + and -)
				//
				if (token.IsPunctuationToken() && !token.IsPlusMinusToken())
					badEquateSymbol = true;

				// Equate symbols cannot be same as matchchar, missing, or gap
				//
				if (ch == missing || ch == matchchar || ch == gap)
					badEquateSymbol = true;

				// Equate symbols cannot be one of the state symbols currently defined
				//
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
				DemandEquals(token, " in EQUATE definition");

				// This should be the token to be substituted in for the equate symbol
				//
				token.SetLabileFlagBit(NxsToken::parentheticalToken);
				token.SetLabileFlagBit(NxsToken::curlyBracketedToken);
				token.GetNextToken();
				NxsString v = token.GetToken();

				// Add the new equate association to the equates list
				//
				equates[k] = v;
				}

			standardDataTypeAssumed = true;
			}

		else if (token.Equals("MATCHCHAR"))
			{
			DemandEquals(token, "after keyword MATCHCHAR");

			// This should be the matchchar symbol (single character)
			//
			token.GetNextToken();

			if (token.GetTokenLength() != 1)
				{
				errormsg = "MATCHCHAR symbol should be a single character, but ";
				errormsg += token.GetToken();
				errormsg += " was specified";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			else if (token.IsPunctuationToken() && !token.IsPlusMinusToken())
				{
				errormsg = "MATCHCHAR symbol specified cannot be a punctuation token (";
				errormsg += token.GetToken();
				errormsg += " was specified) ";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			else if (token.IsWhitespaceToken())
				{
				errormsg = "MATCHCHAR symbol specified cannot be a whitespace character (";
				errormsg += token.GetToken();
				errormsg += " was specified)";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			matchchar = token.GetToken()[0];

			ignoreCaseAssumed = true;
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

		else if (token.Equals("TRANSPOSE"))
			{
			transposing = true;
			standardDataTypeAssumed = true;
			}

		else if (token.Equals("INTERLEAVE"))
			{
			interleaving = true;
			standardDataTypeAssumed = true;
			}

		else if (token.Equals("ITEMS"))
			{
			DemandEquals(token, "after keyword ITEMS");
			items.clear();
			// This should be STATES (no other item is supported at this time)
			//
			token.GetNextToken();
			if (datatype == NxsCharactersBlock::continuous)
				{
				NxsString s;
				if (token.Equals("("))
					{
					token.GetNextToken();
					while (!token.Equals(")"))
						{
						s = token.GetToken();
						s.ToUpper();
						items.push_back(std::string(s.c_str()));
						token.GetNextToken();
						}
					}
				else
					{
					s = token.GetToken();
					s.ToUpper();
					items.push_back(std::string(s.c_str()));
					}
				}
			else 
				{
				if (!token.Equals("STATES"))
					GenerateNxsException(token, "Sorry, only ITEMS=STATES is supported for discrete datatypes at this time");
				items = std::vector<std::string>(1, std::string("STATES"));
				}
			standardDataTypeAssumed = true;
			}

		else if (token.Equals("STATESFORMAT"))
			{
			DemandEquals(token, "after keyword STATESFORMAT");

			// This should be STATESPRESENT (no other statesformat is supported at this time)
			//
			token.GetNextToken();

			if (token.Equals("STATESPRESENT"))
				statesFormat = STATES_PRESENT;
			else
				{
				if (datatype == NxsCharactersBlock::continuous)
					{
					if (token.Equals("INDIVIDUALS"))
						statesFormat = INDIVIDUALS;
					else
						GenerateNxsException(token, "Sorry, only STATESFORMAT=STATESPRESENT or STATESFORMAT=INDIVIDUALS are supported for continuous datatypes at this time");
					}
				else
					GenerateNxsException(token, "Sorry, only STATESFORMAT=STATESPRESENT supported for discrete datatypes at this time");
				}
			standardDataTypeAssumed = true;
			}

		else if (token.Equals("TOKENS"))
			{
			tokens = true;
			standardDataTypeAssumed = true;
			}

		else if (token.Equals("NOTOKENS"))
			{
			if (datatype == NxsCharactersBlock::continuous)
				GenerateNxsException(token, "NOTOKENS is not allowed for the CONTINUOUS datatype");
			tokens = false;
			standardDataTypeAssumed = true;
			}

		else if (token.Equals(";"))
			break;
		}

	// Perform some last checks before leaving the FORMAT command
	//
	if (!tokens && datatype == continuous)
		GenerateNxsException(token, "TOKENS must be defined for DATATYPE=CONTINUOUS");
	if (tokens && (datatype == dna || datatype == rna || datatype == nucleotide))
		GenerateNxsException(token, "TOKENS not allowed for the DATATYPEs DNA, RNA, or NUCLEOTIDE");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called from HandleStdMatrix or HandleTransposedMatrix function to read in the next state. Always returns true 
|	except in the special case of an interleaved matrix, in which case it returns false if a newline character is 
|	encountered before the next token.
*/
bool NxsCharactersBlock::HandleNextState(
  NxsToken &token,	/* the token used to read from `in' */
  unsigned i,		/* the taxon index, in range [0..`ntax') */
  unsigned j)		/* the character index, in range [0..`nchar') */
	{
	if (interleaving)
		token.SetLabileFlagBit(NxsToken::newlineIsToken);
	if (datatype == NxsCharactersBlock::continuous)
		{
		token.SetLabileFlagBit(NxsToken::hyphenNotPunctuation);
		std::vector<double> v;
		std::vector<int> scored;
		token.GetNextToken();
		NxsString t;
		if (interleaving && token.AtEOL())
			return false;
		if (token.Equals("("))
			{
			token.SetLabileFlagBit(NxsToken::hyphenNotPunctuation);		   
			token.GetNextToken();
			if (token.AtEOF())
				GenerateNxsException(token, "Unexpected end of file encountered");
			while (!token.Equals(")"))
				{
				t = token.GetToken();
				if (t.length() == 1 && (t[0] == missing || t[0] == gap))
					{
					v.push_back(DBL_MAX);
					scored.push_back(0);
					}
				else if (t.length() == 1 && t[0] == matchchar)
					{
					v.push_back(DBL_MAX);
					scored.push_back(2);
					}
				else if (!t.IsADouble())
					GenerateUnexpectedTokenNxsException(token, "a number");
				else
					{
					v.push_back(t.ConvertToDouble());
					scored.push_back(1);
					}
				token.SetLabileFlagBit(NxsToken::hyphenNotPunctuation);		   
				token.GetNextToken();
				if (token.AtEOF())
					GenerateNxsException(token, "Unexpected end of file encountered");
				}
			}
		else
			{
			t = token.GetToken();
			if (t.length() == 1 && (t[0] == missing || t[0] == gap))
				{
				v.push_back(DBL_MAX);
				scored.push_back(0);
				}
			else if (t.length() == 1 && t[0] == matchchar)
				{
				v.push_back(DBL_MAX);
				scored.push_back(2);
				}
			else if (!t.IsADouble())
				GenerateUnexpectedTokenNxsException(token, "a number");
			else
				{
				v.push_back(t.ConvertToDouble());
				scored.push_back(1);
				}
			}
		unsigned n_read = v.size();
		if (n_read < items.size())
			{
			errormsg = "For each cell of the MATRIX a value for each of the ";
			errormsg += items.size();
			errormsg += " ITEMS listed in the FORMAT command is expected.\nOnly ";
			errormsg += n_read;
			errormsg += " values read.";
			GenerateNxsException(token);
			}
		// We've read in the state now, so if this character has been eliminated, we don't want to go any further with it
		//
		if (j == UINT_MAX)
			return true;

		if (i > continuousCharMatrix.size())
			GenerateNxsException(token, "Internal Error: taxon index out of range in continuousCharMatrix.");
		ContinuousCharRow & row = continuousCharMatrix[i];
		if (j > row.size())
			GenerateNxsException(token, "Internal Error: character index out of range in continuousCharMatrix.");

		ContinuousCharCell & cell = row[j];
		cell.clear();

		std::vector<std::string >::const_iterator itemIt = items.begin();
		std::string key;
		unsigned curr_ind_in_v = 0;
		for (; itemIt != items.end(); ++itemIt, ++curr_ind_in_v)
			{
			key = *itemIt;
			if (scored[curr_ind_in_v] == 1)
				cell[key] = vector<double>(1, v[curr_ind_in_v]);
			else if (scored[curr_ind_in_v] == 0)
				cell[key] = vector<double>();
			else
				{
				if (i == 0)
					GenerateNxsException(token, "MATCHCHAR cannot be used in the first taxon");
				const vector<double> & first_taxon_vector = continuousCharMatrix[0][j][key];
				if (first_taxon_vector.empty())
					GenerateNxsException(token, "First taxon does not have a value to copy, but a MATCHCHAR was found.");
				else
					cell[key] = vector<double>(first_taxon_vector[0]);
				}
			}
		unsigned curr_ind_mapped = 1;
		if (!key.empty() && curr_ind_in_v < n_read)
			{
			vector<double> & curr_cell_vector = cell[key];
			for (; curr_ind_in_v < n_read; ++curr_ind_in_v, ++curr_ind_mapped)
				{
				if (scored[curr_ind_in_v] == 1)
					curr_cell_vector.push_back(v[curr_ind_in_v]);
				else if (scored[curr_ind_in_v] != 0)
					curr_cell_vector.push_back(DBL_MAX);
				else
					{
					if (i == 0)
						GenerateNxsException(token, "MATCHCHAR cannot be used in the first taxon");
					const vector<double> & first_taxon_vector = continuousCharMatrix[0][j][key];
					if (first_taxon_vector.size() < curr_ind_mapped+1)
						GenerateNxsException(token, "First taxon does not have a value to copy, but a MATCHCHAR was found.");
					else
						curr_cell_vector.push_back(first_taxon_vector[curr_ind_mapped]);
					}
				}
			}
		return true;
		}
	
	// This should be the state for taxon i and character j
	//
	if (!tokens)
		{
		token.SetLabileFlagBit(NxsToken::parentheticalToken);
		token.SetLabileFlagBit(NxsToken::curlyBracketedToken);
		token.SetLabileFlagBit(NxsToken::singleCharacterToken);
		}


	token.GetNextToken();

	if (interleaving && token.AtEOL())
		return false;

	// Make sure we didn't run out of file
	//
	if (token.AtEOF())
		GenerateNxsException(token, "Unexpected end of file encountered");

	// If we didn't run out of file, there is no reason why we should have a zero-length token on our hands
	//
	assert(token.GetTokenLength() > 0);

	// We've read in the state now, so if this character has been eliminated, we don't want to go any further with it
	//
	if (j == UINT_MAX)
		return true;

	// See if any equate macros apply
	//
	NxsString skey = NxsString(token.GetToken(true)); // equates should always respect case

	NxsStringMap::iterator p = equates.find(skey);
	if (p != equates.end())
		{
		NxsString sval = (*p).second;
		token.ReplaceToken(sval.c_str());
		}

	// Handle case of single-character state symbol
	//
	if (!tokens && token.GetTokenLength() == 1)
		{
		char ch = token.GetToken()[0];

		// Check for missing data symbol
		//
		if (ch == missing)
			{
			matrix->SetMissing(i, j);
			}

		// Check for matchchar symbol
		//
		else if (matchchar != '\0' && ch == matchchar)
			{
			matrix->CopyStatesFromFirstTaxon(i, j);
			}

		// Check for gap symbol
		//
		else if (gap != '\0' && ch == gap)
			{
			matrix->SetGap(i, j);
			}

		// Look up the position of this state in the symbols array
		//
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
				errormsg += (j+1);
				errormsg += ", not found in list of valid symbols";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}
			matrix->AddState(i, j, p);
			matrix->SetPolymorphic(i, j, 0);
			}
		}	// if (!tokens && token.GetTokenLength() == 1)

	// Handle case of state sets when tokens is not in effect
	//
	else if (!tokens && token.GetTokenLength() > 1)
		{
		// Token should be in one of the following forms: LEFT_SQUIGGLYacgRIGHT_SQUIGGLY LEFT_SQUIGGLYa~gRIGHT_SQUIGGLY LEFT_SQUIGGLYa c gRIGHT_SQUIGGLY (acg) (a~g) (a c g) 
		//
		NxsString t = token.GetToken();
		unsigned tlen = (unsigned)t.size();
		unsigned poly = (t[0] == '(');
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

		unsigned k = 1;
		char *pFirst = symbols;
		bool tildeFound = false;
		for (;;)
			{
			if (t[k] == ')' || t[k] == '}')
				break;

			if (t[k] == ' ' || t[k] == '\t' || t[k] == ',')
				{
				k++;
				continue;
				}

			// t[k] should be either '~' or one of the state symbols
			//
			if (t[k] == '~')
				{
				tildeFound = true;
				}
			else
				{
				// Add state symbol and record if it is the first or last one in case we encounter a tilde
				//
				if (tildeFound)
					{
					// Add all states from firstState to t[k] then set tildeFound to false again
					//
					pFirst++;
					while (*pFirst != '\0' && *pFirst != t[k])
						{
						int p = PositionInSymbols(*pFirst);
						if (p < 0)
							{
							errormsg = "State specified (";
							errormsg += *pFirst;
							errormsg += ") for taxon ";
							errormsg += (i+1);
							errormsg += ", character ";
							errormsg += (j+1);
							errormsg += ", not found in list of valid symbols";
							throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
							}
						matrix->AddState(i, j, p);
						pFirst++;
						}

					tildeFound = false;
					}
				else
					{
					int p = PositionInSymbols(t[k]);
					if (p < 0)
						{
						errormsg = "State specified (";
						errormsg += t[k];
						errormsg += ") for taxon ";
						errormsg += (i+1);
						errormsg += ", character ";
						errormsg += (j+1);
						errormsg += ", not found in list of valid symbols";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}
					pFirst = (symbols + p);
					matrix->AddState(i, j, p);
					}

				} // if (t[k] == '~') ... else ... loop

			k++;
			} // for (;;) loop

		matrix->SetPolymorphic(i, j, poly);
		}	// if (!tokens && token.GetTokenLength() == 1) ... else if (!tokens && token.GetTokenLength() > 1)

	// Handle case in which TOKENS was specified in the FORMAT command
	//
	else
		{
		// Token should be in one of the following forms: "LEFT_SQUIGGLY"  "a"	"bb"
		//
		int polymorphism = token.Equals("(");
		int uncertainty	 = token.Equals("LEFT_SQUIGGLY");

		if (!uncertainty && !polymorphism)
			{
			int k = HandleTokenState(token, j);
			matrix->AddState(i, j, k);
			}

		else	// either uncertainty or polymorphism
			{
			bool tildeFound = false;
			unsigned first = UINT_MAX;
			unsigned last;
			for (;;)
				{
				// OPEN ISSUE: What about newlines if interleaving? I'm assuming
				// that the newline must come between characters to count.

				token.SetLabileFlagBit(NxsToken::tildeIsPunctuation);
				token.GetNextToken();

				if (polymorphism && token.Equals(")"))
					{
					if (tildeFound)
						GenerateNxsException(token, "Range of states still being specified when ')' encountered");
					break;
					}

				else if (uncertainty && token.Equals("RIGHT_SQUIGGLY"))
					{
					if (tildeFound)
						GenerateNxsException(token, "Range of states still being specified when '}' encountered");
					break;
					}

				else if (token.Equals("~"))
					{
					if (first == UINT_MAX)
						GenerateNxsException(token, "Tilde character ('~') cannot precede token indicating beginning of range");
					tildeFound = true;
					}

				else if (tildeFound)
					{
					// Add all states from first+1 to last, then reset tildeFound to false
					//
					last = HandleTokenState(token, j);

					if (last <= first)
						{
						errormsg = "Last state in specified range (";
						errormsg += token.GetToken();
						errormsg += ") must be greater than the first";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}

					for (unsigned k = first+1; k <= last; k++)
						matrix->AddState(i, j, k);

					tildeFound = false;
					first = UINT_MAX;
					}

				else
					{
					// Add current state, then set first to that state's value
					// State's value is its position within the list of states
					// for that character
					//
					first = HandleTokenState(token, j);
					matrix->AddState(i, j, first);
					}
				}	// if (!uncertainty && !polymorphism) ... else 

			if (polymorphism)
				matrix->SetPolymorphic(i, j, 1);
			}	// if (!uncertainty && !polymorphism) ... else

		}	// if (!tokens && token.GetTokenLength() == 1) ... else if (!tokens && token.GetTokenLength() > 1) ... else 

	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called from HandleNextState to read in the next state when TOKENS was specified. Looks up state in character 
|	states listed for the character to make sure it is a valid state, and returns state's value (0, 1, 2, ...). Note: 
|	does NOT handle adding the state's value to matrix. Save the return value (call it k) and use the following command
|	to add it to matrix: matrix->AddState(i, j, k);
*/
unsigned NxsCharactersBlock::HandleTokenState(
  NxsToken &token,	/* the token used to read from `in' */
  unsigned j)		/* the character index, in range [0..`nchar') */
	{
	// Token should be one of the character states listed for character j in charStates
	//
	if (charStates.find(j) == charStates.end())
		{
		errormsg = "No states were defined for character ";
		errormsg += (1 + GetOrigCharIndex(j));
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	// TO DO: this section is very UGLY - need to find some cleaner way of comparing
	// the token NxsString to the strings representing valid characters states
	// in the NxsStringVector associated with character j
	//
	NxsStringVectorMap::const_iterator bagIter	= charStates.find(j);
	NxsStringVector::const_iterator ci_begin	= (*bagIter).second.begin();
	NxsStringVector::const_iterator ci_end		= (*bagIter).second.end();
	NxsString t									= token.GetToken(respectingCase);
	NxsStringVector::const_iterator cit;
	if (respectingCase)
		cit = find(ci_begin, ci_end, t);
	else
		cit = find_if (ci_begin, ci_end, bind2nd(NxsStringEqual(), t));

	if (cit == ci_end)
		{
		errormsg = "Character state ";
		errormsg += t;
		errormsg += " not defined for character ";
		errormsg += (1 + GetOrigCharIndex(j));
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	// Ok, the state has been identified, so return the state's internal representation. That is, 
	// if the list of state labels was "small medium large" and "small" was specified in the data file,
	// state saved in matrix would be 0 (it would be 1 if "medium" were specified in the data file, 
	// and 2 if "large" were specified in the data file).
	//
	unsigned k = (unsigned)(cit - ci_begin);
	return k;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called from HandleMatrix function to read in a standard (i.e., non-transposed) matrix. Interleaving, if 
|	applicable, is dealt with herein.
*/
void NxsCharactersBlock::HandleStdMatrix(
  NxsToken &token)	/* the token used to read from `in' */
	{
	assert(charPos != NULL);
	assert(taxonPos != NULL);

	unsigned i, j, currChar;
	unsigned firstChar = 0;
	unsigned lastChar = ncharTotal;
	unsigned nextFirst = 0;
	unsigned page = 0;

	for (;;)
		{
		//************************************************
		//******** Beginning of loop through taxa ********
		//************************************************

		for (i = 0; i < ntax; i++)
			{
			if (labels)
				{
				// This should be the taxon label
				//
				token.GetNextToken();

				if (page == 0 && newtaxa)
					{
					// This section:
					// - labels provided
					// - on first (or only) interleave page
					// - no previous TAXA block

					// Check for duplicate taxon names
					//
					if (taxa->IsAlreadyDefined(token.GetToken()))
						{
						errormsg = "Data for this taxon (";
						errormsg += token.GetToken();
						errormsg += ") has already been saved";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}

					// Labels provided and not already stored in the taxa block with
					// the TAXLABELS command; taxa->Reset() and taxa->SetNTax() have
					// were already called, however, when the NTAX subcommand was
					// processed.
					//
					taxa->AddTaxonLabel(token.GetToken());

					// Order of occurrence in TAXA block same as row in matrix
					//
					taxonPos[i] = i;

					}	// if (page == 0 && newtaxa)

				else
					{
					// This section:
					// - labels provided
					// - TAXA block provided or has been created already
					// - may be on any (interleave) page					

					// Cannot assume taxon in same position in
					// taxa block. Set up taxonPos array so that we can look up
					// the correct row in matrix for any given taxon
					//
					unsigned positionInTaxaBlock;
					try
						{
						positionInTaxaBlock = taxa->FindTaxon(token.GetToken());
						}
					catch(NxsTaxaBlock::NxsX_NoSuchTaxon)
						{
						if (token.Equals(";") && i == 0)
							{
							errormsg = "Unexpected ; (after only ";
							errormsg += currChar;
							errormsg += " characters were read)";
							}
						else
							{
							errormsg = "Could not find taxon named ";
							errormsg += token.GetToken();
							errormsg += " among stored taxon labels";
							}
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}

					if (page == 0)
						{
						// Make sure user has not duplicated data for a single taxon
						//
						if (taxonPos[positionInTaxaBlock] != UINT_MAX)
							{
							errormsg = "Data for this taxon (";
							errormsg += token.GetToken();
							errormsg += ") has already been saved";
							throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
							}

						// Make sure user has kept same relative ordering of taxa in both the TAXA
						// block and the CHARACTERS block
						//
						if (positionInTaxaBlock != i)
							GenerateNxsException(token, "Relative order of taxa must be the same in both the TAXA and CHARACTERS blocks");
						taxonPos[i] = positionInTaxaBlock;
						}	// if (page == 0)

					else
						{
						// Make sure user has kept the ordering of taxa the same from one interleave page to the next
						//
						if (taxonPos[positionInTaxaBlock] != i)
							GenerateNxsException(token, "Ordering of taxa must be identical to that in first interleave page");
						}	// // if (page == 0) ... else
					}	// if (page == 0 && newtaxa) ... else
				}	// if (labels)

			else
				{
				// No labels provided, assume taxon position same as in taxa block
				//
				if (page == 0)
					taxonPos[i] = i;
				}	// if (labels) ... else

			//******************************************************
			//******** Beginning of loop through characters ********
			//******************************************************

			for (currChar = firstChar; currChar < lastChar; currChar++)
				{
				// It is possible that character currChar has been eliminated, in which case we need to 
				// go through the motions of reading in the data but we don't store it. The variable j 
				// will be our guide when it comes time to store data since j will be UINT_MAX for
				// characters that were eliminated and will be set to the correct row for characters 
				// that have not been eliminated.
				//
				j = charPos[currChar];

				// ok will be false only if a newline character is encountered before character j processed
				//
				bool ok = HandleNextState(token, i, j);
				if (interleaving && !ok)
					{
					if (lastChar < ncharTotal && j != lastChar)
						GenerateNxsException(token, "Each line within an interleave page must comprise the same number of characters");

					// currChar should be firstChar in next go around
					//
					nextFirst = currChar;

					// Set lastChar to currChar so that we can check to make sure the remaining lines 
					// in this interleave page end at the same place
					//
					lastChar = currChar;

					// Since j is now equal to lastChar, we are done with this innermost loop
					}
				} // for (currChar = firstChar; currChar < lastChar; currChar++)

			} // for (i = 0; i < ntax; i++)

		firstChar = nextFirst;
		lastChar = ncharTotal;

		// If currChar equals ncharTotal, then we've just finished reading the last interleave page 
		// and thus should break from the outer loop. Note that if we are not interleaving, this will 
		// still work since lastChar is initialized to ncharTotal and never changed
		//
		if (currChar == ncharTotal)
			break;

		page++;
		} // for (;;)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called from HandleMatrix function to read in a transposed matrix. Interleaving, if applicable, is dealt with herein.
*/
void NxsCharactersBlock::HandleTransposedMatrix(
  NxsToken &token)	/* the token used to read from in */
	{
	assert(charPos != NULL);
	assert(taxonPos != NULL);

	unsigned i, j, currChar;
	unsigned firstTaxon = 0;
	unsigned lastTaxon = ntaxTotal;
	unsigned nextFirst = 0;
	unsigned page = 0;

	for (;;)
		{
		//******************************************************
		//******** Beginning of loop through characters ********
		//******************************************************

		for (currChar = 0; currChar < ncharTotal; currChar++)
			{
			// It is possible that character currChar has been eliminated, in
			// which case we need to go through the motions of reading in the
			// data but we don't store it.	The variable j will be our guide
			// when it comes time to store data since j will be UINT_MAX for
			// characters that were eliminated and will be set to the correct
			// row for characters that have not been eliminated.
			//
			j = charPos[currChar];

			if (labels)
				{
				// This should be the character label
				//
				token.GetNextToken();

				if (page == 0 && newchar)
					{
					// Check for duplicate character names
					//
					NxsString s = token.GetToken();
					NxsStringVector::const_iterator iter = find(charLabels.begin(), charLabels.end(), s);

					bool charLabelFound = (iter != charLabels.end());
					if (charLabelFound)
						{
						errormsg = "Data for this character (";
						errormsg += token.GetToken();
						errormsg += ") has already been saved";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}

					// Labels provided, need to add them to charLabels list. We're not supposed to save 
					// anything for this character since it has been eliminated, but the labels must be 
					// saved for purposes of numbering. Otherwise a more complicated system would be needed
					// wherein an association is set up between character number and character label. Since
					// this is not done in the case of taxa that are effectively eliminated when they are 
					// included in the TAXA block but not in the CHARACTERS MATRIX command, I see no reason
					// to not save the full character labels here even for those that have been eliminated.
					// Also, for interleaved matrices, it is necessary to have the full labels saved somewhere
					// so that it is possible to detect characters out of order or duplicated.
					//
					charLabels.push_back(token.GetToken());
					}	// if (page == 0 && newchar)

				else // either not first interleaved page or character labels not previously defined
					{
					NxsString s = token.GetToken();
					NxsStringVector::const_iterator iter = find(charLabels.begin(), charLabels.end(), s);

					if (iter == charLabels.end())
						{
						errormsg = "Could not find character named ";
						errormsg += token.GetToken();
						errormsg += " among stored character labels";
						throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
						}

					unsigned positionInCharLabelsList = (unsigned)(iter - charLabels.begin());

					// Make sure user has not duplicated data for a single character or changed the order 
					// in which characters appear in different interleave pages
					//
					if (positionInCharLabelsList != currChar)
						{
						if (page == 0)
							{
							errormsg = "Data for this character (";
							errormsg += token.GetToken();
							errormsg += ") has already been saved";
							throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
							}
						else
							{
							errormsg = "Ordering of characters must be identical to that in first interleave page";
							throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
							}
						}
					}	// if (page == 0 && newchar) ... else
				} // if (labels)

			//************************************************
			//******** Beginning of loop through taxa ********
			//************************************************

			for (i = firstTaxon; i < lastTaxon; i++)
				{
				if (page == 0)
					{
					// We are forced to assume that the user did not leave out any
					// taxa, since without taxon labels in the matrix we would
					// have no way of detecting which were left out; thus,
					// ntax == ntaxTotal in this case.	Order of occurrence in
					// TAXA block is the same as the row in matrix.
					//
					taxonPos[i] = i;
					}

				// ok will be 0 only if a newline character is encountered before
				// taxon i processed
				//
				bool ok = HandleNextState(token, i, j);
				if (interleaving && !ok)
					{
					if (lastTaxon < ntaxTotal && i != lastTaxon)
						GenerateNxsException(token, "Each line within an interleave page must comprise the same number of taxa");

					// i should be firstChar in next go around
					//
					nextFirst = i;

					// Set lastTaxon to i so that we can check to make sure the
					// remaining lines in this interleave page end at the same place
					lastTaxon = i;

					// Since i is now equal to lastTaxon, we are done with this innermost loop
					}	// if (interleaving && !ok)
				} // for (i = firstTaxon; i < lastTaxon; i++)

			} // for (currChar = 0; currChar < ncharTotal; currChar++)

		firstTaxon = nextFirst;
		lastTaxon = ntaxTotal;

		// If i equals ncharTotal, then we've just finished reading the last
		// interleave page and thus should break from the outer loop
		// Note that if we are not interleaving, this will still work since
		// lastTaxon is initialized to ntaxTotal and never changed
		//
		if (i == ntaxTotal)
			break;

		page++;
		} // for (;;)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when MATRIX command needs to be parsed from within the CHARACTERS block. Deals with everything after the 
|	token MATRIX up to and including the semicolon that terminates the MATRIX command.
*/
void NxsCharactersBlock::HandleMatrix(
  NxsToken &token)	/* the token used to read from `in' */
	{
	unsigned i, j;

	if (ntax == 0)
		{
		errormsg = "Must precede ";
		errormsg += id;
		errormsg += " block with a TAXA block or specify NEWTAXA and NTAX in the DIMENSIONS command";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}

	if (ntaxTotal == 0)
		ntaxTotal = taxa->GetNumTaxonLabels();

	if (matrix != NULL)
		{
		delete matrix;
		matrix = NULL;
		}
	continuousCharMatrix.clear();

	// Allocate memory for (and initialize) the arrays activeTaxon and activeChar.
	// All characters and all taxa are initially active.
	//
	activeTaxon = new bool[ntax];
	for (i = 0; i < ntax; i++)
		activeTaxon[i] = true;

	activeChar = new bool[nchar];
	for (j = 0; j < nchar; j++)
		activeChar[j] = true;

	// The value of ncharTotal is normally identical to the value of nchar specified
	// in the CHARACTERS block DIMENSIONS command.	If an ELIMINATE command is
	// processed, however, nchar < ncharTotal.	Note that the ELIMINATE command
	// will have already been read by now, and the eliminated character numbers
	// will be stored in the NxsUnsignedSet eliminated.
	//
	// Note that if an ELIMINATE command has been read, charPos will have already
	// been created; thus, we only need to allocate and initialize charPos if user
	// did not specify an ELIMINATE command
	//
	if (charPos == NULL)
		BuildCharPosArray();

	// The value of ntaxTotal equals the number of taxa specified in the
	// TAXA block, whereas ntax equals the number of taxa specified in
	// the DIMENSIONS command of the CHARACTERS block.	These two numbers
	// will be identical unless some taxa were left out of the MATRIX
	// command of the CHARACTERS block, in which case ntax < ntaxTotal.
	//
	if (taxonPos != NULL)
		delete [] taxonPos;
	taxonPos = new unsigned[ntaxTotal];

	for (i = 0; i < ntaxTotal; i++)
		taxonPos[i] = UINT_MAX;

	if (datatype == NxsCharactersBlock::continuous)
		{
		ContinuousCharRow row(nchar);
		continuousCharMatrix.assign(ntax, row);
		}
	else 
		{
		matrix = new NxsDiscreteMatrix(ntax, nchar);
		}
	if (transposing)
		HandleTransposedMatrix(token);
	else
		HandleStdMatrix(token);
	DemandEndSemicolon(token, "MATRIX");
	// If we've gotten this far, presumably it is safe to
	// tell the ASSUMPTIONS block that were ready to take on
	// the responsibility of being the current character-containing
	// block (to be consulted if characters are excluded or included
	// or if taxa are deleted or restored)
	//
	assumptionsBlock->SetCallback(this);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when STATELABELS command needs to be parsed from within the DIMENSIONS block. Deals with everything after 
|	the token STATELABELS up to and including the semicolon that terminates the STATELABELS command. Note that the 
|	numbers of states are shifted back one before being stored so that the character numbers in the NxsStringVectorMap 
|	objects are 0-offset rather than being 1-offset as in the NxsReader data file.
*/
void NxsCharactersBlock::HandleStatelabels(
  NxsToken &token)	/* the token used to read from `in' */
	{
	bool semicolonFoundInInnerLoop = false;
	if (datatype == continuous)
		GenerateNxsException(token, "STATELABELS cannot be specified when the datatype is continuous");
	charStates.clear();

	if (charPos == NULL)
		BuildCharPosArray();

	for (;;)
		{
		if (semicolonFoundInInnerLoop)
			break;

		token.GetNextToken();

		if (token.Equals(";"))
			break;

		// Token should be the character number; create a new association
		//
		unsigned n = atoi(token.GetToken().c_str());

		if (n < 1 || n > ncharTotal)
			{
			errormsg = "Invalid character number (";
			errormsg += token.GetToken();
			errormsg += ") found in STATELABELS command (either out of range or not interpretable as an integer)";
			throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
			}

		for (;;)
			{
			token.GetNextToken();

			if (token.Equals(";"))
				{
				semicolonFoundInInnerLoop = true;
				break;
				}

			if (token.Equals(","))
				break;

			// Token should be a character state label; add it to the list
			//
			if (!IsEliminated(n - 1))
				{
				unsigned k = GetCharPos(n - 1);
				charStates[k].push_back(token.GetToken());
				}

			} // for (;;)
		} // for (;;)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when TAXLABELS command needs to be parsed from within the CHARACTERS block. Deals with everything after the 
|	token TAXLABELS up to and including the semicolon that terminates the TAXLABELS command.
*/
void NxsCharactersBlock::HandleTaxlabels(
  NxsToken &token)	/* the token used to read from `in' */
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
		//
		if (token.Equals(";"))
			{
			break;
			}
		else
			{
			// Check to make sure user is not trying to read in more
			// taxon labels than there are taxa
			//
			if (taxa->GetNumTaxonLabels() > ntaxTotal)
				GenerateNxsException(token, "Number of taxon labels exceeds NTAX specified in DIMENSIONS command");

			taxa->AddTaxonLabel(token.GetToken());
			}
		}

	// OPEN ISSUE: Some may object to setting newtaxa to false here, because then the fact that new taxa were 
	// specified in this CHARACTERS block rather than in a preceding TAXA block is lost. This will only be 
	// important if we wish to recreate the original data file, which I don't anticipate anyone doing with
	// this code (too difficult to remember all comments, the order of blocks in the file, etc.)

	newtaxa = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns position of `ch' in `symbols' array. The value of `respectingCase' is used to determine whether the search 
|	should be case sensitive or not. Assumes `symbols' is non-NULL. Returns UINT_MAX if `ch' is not found in `symbols'.
*/
unsigned NxsCharactersBlock::PositionInSymbols(
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
void NxsCharactersBlock::Read(
  NxsToken &token)	/* the token used to read from `in' */
	{
	isEmpty = false;
	isUserSupplied = true;

	NxsString s;
	s = "BEGIN ";
	s += id;
	DemandEndSemicolon(token, s.c_str());
	
	ntax = taxa->GetNumTaxonLabels();

	for (;;)
		{
		token.GetNextToken();

		if (token.Equals("DIMENSIONS"))
			{
			HandleDimensions(token, "NEWTAXA", "NTAX", "NCHAR");
			}
		else if (token.Equals("FORMAT"))
			{
			HandleFormat(token);
			}
		else if (token.Equals("ELIMINATE"))
			{
			HandleEliminate(token);
			}
		else if (token.Equals("TAXLABELS"))
			{
			HandleTaxlabels(token);
			}
		else if (token.Equals("CHARSTATELABELS"))
			{
			HandleCharstatelabels(token);
			}
		else if (token.Equals("CHARLABELS"))
			{
			HandleCharlabels(token);
			}
		else if (token.Equals("STATELABELS"))
			{
			HandleStatelabels(token);
			}
		else if (token.Equals("MATRIX"))
			{
			HandleMatrix(token);
			}
		else if (token.Equals("END"))
			{
			HandleEndblock(token, "Character");
			break;
			}
		else if (token.Equals("ENDBLOCK"))
			{
			HandleEndblock(token, "Character");
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
				GenerateNxsException(token, "Unexpected end of file encountered");
			}	// else
		}	// for (;;)
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function outputs a brief report of the contents of this CHARACTERS block. Overrides the abstract virtual 
|	function in the base class.
*/
void NxsCharactersBlock::Report(
  ostream &out) /* the output stream to which to write the report */
	{
	out << endl;
	out << id << " block contains ";
	if (ntax == 0)
		out << "no taxa";
	else if (ntax == 1)
		out << "one taxon";
	else
		out << ntax << " taxa";
	out << " and ";
	if (nchar == 0)
		out << "no characters";
	else if (nchar == 1)
		out << "one character";
	else
		out << nchar << " characters";
	out << endl;

	if (formerly_datablock)
		{
		out << "Originally read in as a DATA block" << endl;
		out << endl;
		}

	out << "  Data type is \"" << this->GetDatatypeName() << "\"" << endl;
	
	if (respectingCase)
		out << "  Respecting case" << endl;
	else
		out << "  Ignoring case" << endl;

	if (tokens)
		out << "  Multicharacter tokens allowed in data matrix" << endl;
	else
		out << "  Data matrix entries are expected to be single symbols" << endl;

	if (labels && transposing)
		out << "  Character labels are expected on left side of matrix" << endl;
	else if (labels && !transposing)
		out << "  Taxon labels are expected on left side of matrix" << endl;
	else
		out << "  No labels are expected on left side of matrix" << endl;

	if (charLabels.size() > 0)
		{
		out << "  Character and character state labels:" << endl;
		for (unsigned k = 0; k < nchar; k++) 
			{
			if (charLabels[k].length() == 0)
				out << '\t' << (1 + GetOrigCharIndex(k)) << '\t' << "(no label provided for this character)" << endl;
			else
				out << '\t' << (1 + GetOrigCharIndex(k)) << '\t' << charLabels[k] << endl;

			// Output state labels if any are defined for this character
			//
			NxsStringVectorMap::const_iterator cib = charStates.find(k);
			if (cib != charStates.end())
				{
				int ns = (int)(*cib).second.size();
				for (int m = 0; m < ns; m++)
					{
					out << "\t\t" << (*cib).second[m] << endl;
					}
				}
			}
		}

	if (transposing && interleaving)
		out << "  Matrix transposed and interleaved" << endl;
	else if (transposing && !interleaving)
		out << "  Matrix transposed but not interleaved" << endl;
	else if (!transposing && interleaving)
		out << "  Matrix interleaved but not transposed" << endl;
	else
		out << "  Matrix neither transposed nor interleaved" << endl;

	out << "  Missing data symbol is '" << missing << '\'' << endl;

	if (matchchar != '\0')
		out << "  Match character is '" << matchchar << '\'' << endl;
	else
		out << "  No match character specified" << endl;

	if (gap != '\0')
		out << "  Gap character specified is '" << gap << '\'' << endl;
	else
		out << "  No gap character specified" << endl;

	out << "  Valid symbols are: " << symbols << endl;

	int numEquateMacros = (int)equates.size();
	if (numEquateMacros > 0)
		{
		out << "  Equate macros in effect:" << endl;
		typedef NxsStringMap::const_iterator CI;
		for (CI i = equates.begin(); i != equates.end(); ++i)
			{
			out << "	" << (*i).first << " = " << (*i).second << endl;
			}
		}
	else
		out << "  No equate macros have been defined" << endl;

	if (ncharTotal == nchar)
		out << "  No characters were eliminated" << endl;
	else
		{
		out << "  The following characters were eliminated:" << endl;
		NxsUnsignedSet::const_iterator k;
		for (k = eliminated.begin(); k != eliminated.end(); k++)
			{
			out << "	" << ((*k)+1) << endl;
			}
		}

	out << "  The following characters have been excluded:" << endl;
	unsigned k;
	unsigned nx = 0;
	for (k = 0; k < nchar; k++)
		{
		if (activeChar[k])
			continue;
		out << "	" << (k+1) << endl;
		nx++;
		}

	if (nx == 0)
		out << "	(no characters excluded)" << endl;

	out << "  The following taxa have been deleted:" << endl;
	nx = 0;
	for (k = 0; k < ntax; k++)
		{
		if (activeTaxon[k])
			continue;
		out << "	" << (k+1) << endl;
		nx++;
		}

	if (nx == 0)
		out << "	(no taxa deleted)" << endl;

	out << "  Data matrix:" << endl;
	DebugShowMatrix(out, false, "	 ");
	}

void NxsCharactersBlock::WriteAsNexus(std::ostream &out) const
	{
	writeAllChars = true; // include eliminated chars
	
	out << "BEGIN CHARACTERS;\n";
	this->WriteTitleCommand(out);

	out << "\tDIMENSIONS";
	if (this->newtaxa)
		out << " NewTaxa NTax=" << this->ntax; 
	out << " NChar=" << (this->writeAllChars ? this->nchar : this->ncharTotal)<< ";";

	out << "\n\t";
	this->WriteEliminateCommand(out);

	if (this->newtaxa)
		{
		out << "\n\t";
		this->taxa->WriteTaxLabelsCommand(out);
		}

	out << "\n\t";
	this->WriteCharStateLabelsCommand(out);

	out << "\n\t";
	this->WriteFormatCommand(out);
	
	out << "\n";
	this->WriteMatrixCommand(out);

	out << "\nEND;\n";
	}

std::vector<unsigned> NxsCharactersBlock::GetOrigIndexVector() const
	{
	std::vector<unsigned> v(this->nchar);
	for (unsigned i = 0; i < this->nchar; ++i)
		v[i] = this->GetOrigCharIndex(i);
	return v;
	}

std::vector<unsigned> NxsCharactersBlock::GetOrigMatrixIndicesToWrite() const
	{
	if (writeAllChars)
		{
		std::vector<unsigned> origIndexVec;
		for (unsigned i = 0; i < this->ncharTotal; ++i)
			origIndexVec.push_back(i);			//@TEMP = GetOrigIndexVector(this->nCharTotal);
		return origIndexVec;
		}
	return GetOrigIndexVector();
	}

void NxsCharactersBlock::WriteEliminateCommand(
  ostream &out) const /* output stream on which to print matrix */
	{
	if (this->eliminated.empty() || !this->writeAllChars)
		return;
	out << "ELIMINATE";
	for (NxsUnsignedSet::const_iterator u = this->eliminated.begin(); u != this->eliminated.end(); ++u)
		out << (1 + *u);
	out << ";";
	}

void NxsCharactersBlock::WriteMatrixCommand(
  ostream &out) const /* output stream on which to print matrix */
	{
	assert(charPos != NULL);
	assert(taxonPos != NULL);

	unsigned width = taxa->GetMaxTaxonLabelLength();
	unsigned first_taxon = UINT_MAX;
	const bool use_matchchar = (this->matchchar != '\0');
	out << "Matrix";
	
	std::vector<unsigned> origIndexVec = this->GetOrigMatrixIndicesToWrite();
	if (datatype == continuous)
		out << setprecision(10);
	for (unsigned i = 0; i < this->ntaxTotal; i++)
		{
		out << '\n';
		// Grab taxon name from taxa block. Taxa may not have been presented in the matrix in the same order
		// as they were stored in the taxa block, so use taxonPos array to spit them out in the order they 
		// appeared in the TAXA command. If the taxonPos cell is UINT_MAX, then that means there is no row of
		// the data matrix corresponding to that taxon.
		//
		if (this->taxonPos[i] != UINT_MAX)
			{
			if (first_taxon == UINT_MAX)
				first_taxon = i;

			const std::string currTaxonLabel = NxsToken::EscapeString(taxa->GetTaxonLabel(taxonPos[i]));
			out << currTaxonLabel;

			// Print out enough spaces to even up the left edge of the matrix output
			//
			unsigned currTaxonLabelLen = (unsigned)currTaxonLabel.size();
			unsigned diff = width - currTaxonLabelLen;
			for (unsigned k = 0; k < diff+5; k++)
				out << ' ';
			WriteStatesForMatrixRow(out, i, (use_matchchar ? first_taxon : UINT_MAX), origIndexVec);
			}
		}
	out << ';';
	if (datatype == continuous)
		out << setprecision(6);
	}

void NxsCharactersBlock::WriteStatesForMatrixRow(
  ostream &out,				/* the output stream on which to write */
  unsigned currTaxonIndex,	/* the taxon, in range [0..`ntax') */
  unsigned first_taxon,		/* the index of the first taxon (if UINT_MAX, don't use matchchar) */
  const std::vector<unsigned> &oiv) const /* the original indices of all non eliminated chars */
	{
	NxsCharactersBlock * t = const_cast<NxsCharactersBlock *>(this);
	if (this->tokens)
		{
		for (std::vector<unsigned>::const_iterator oIt = oiv.begin(); oIt != oiv.end(); ++oIt)
			{
			const unsigned j = *oIt; // should invert this so we can walk down charPos
			out << ' ';
			if (datatype == continuous)
				{
				t->ShowStateLabels(out, currTaxonIndex, j, first_taxon);
				}
			else
				{
				unsigned n = this->matrix->GetNumStates(currTaxonIndex, j);
				if (n == 0 && this->matrix->IsGap(currTaxonIndex, j))
					out << this->gap;
				else if (n == 0 && this->matrix->IsMissing(currTaxonIndex, j))
					out << this->missing;
				else if (n == 1) 
					{
					unsigned s = matrix->GetState(currTaxonIndex, j);
					const bool use_matchchar = ((currTaxonIndex > first_taxon) && (matrix->GetState(first_taxon, j) == s));
					if (use_matchchar)
						out << this->matchchar;
					else
						{
						NxsStringVectorMap::const_iterator ci = charStates.find(j);
						if (ci == charStates.end())
							out << s ;
						else
							out << NxsToken::EscapeString((*ci).second[s]);
						}
					}
				else 
					{
					out <<	(matrix->IsPolymorphic(currTaxonIndex, j) ? '(' : '{');
					for (unsigned k = 0; k < n; k++)
						{
						unsigned s = matrix->GetState(currTaxonIndex, j, k);
						NxsStringVectorMap::const_iterator ci = charStates.find(j);
						if (ci == charStates.end())
							out << s;
						else
							out << NxsToken::EscapeString((*ci).second[s]);
						}
					out <<	(matrix->IsPolymorphic(currTaxonIndex, j) ? ')' : '}');
					}
				}
			}
		}	
	else
		{
		for (std::vector<unsigned>::const_iterator oIt = oiv.begin(); oIt != oiv.end(); ++oIt)
			{
			const unsigned j = *oIt; // should invert this so we can walk down charPos
			if (currTaxonIndex > first_taxon) 
				{
				char s[NCL_MAX_STATES + 3];
				t->WriteStates(matrix->GetDiscreteDatum(currTaxonIndex, j), s, NCL_MAX_STATES + 3);
	
				char ss[NCL_MAX_STATES + 3];
				t->WriteStates(matrix->GetDiscreteDatum(first_taxon, j), ss, NCL_MAX_STATES + 3);
	
				if (strcmp(s, ss) == 0)
					out << this->matchchar;
				else
					t->ShowStates(out, currTaxonIndex, j);
				}
			else 
				t->ShowStates(out, currTaxonIndex, j);
			}
		}
	}
	

void NxsCharactersBlock::WriteCharLabelsCommand(std::ostream &out) const
	{
	const std::vector<unsigned> origIndexVec = this->GetOrigMatrixIndicesToWrite();
	bool isFirst = true;
	const unsigned nLabels = this->charLabels.size();
	for (std::vector<unsigned>::const_iterator oIt = origIndexVec.begin(); oIt != origIndexVec.end(); ++oIt)
		{
		const unsigned k = *oIt; // should invert this so we can walk down charPos
		if (k >= nLabels)
			break;
		const NxsString cLabel = this->charLabels[k];
		if (cLabel.length() > 0)
			{
			out <<	(isFirst ? "CHARLABELS " : "") ;
			out << ' ' << NxsToken::EscapeString(cLabel);
			isFirst = false;
			}
		}
	if (!isFirst)
		out << ';';
	}

void NxsCharactersBlock::WriteCharStateLabelsCommand(std::ostream &out) const
	{
	const std::vector<unsigned> origIndexVec = this->GetOrigMatrixIndicesToWrite();

	if (charStates.empty())
		{
		this->WriteCharLabelsCommand(out);
		return;
		}
	const NxsString mtString;
	const unsigned nLabels = this->charLabels.size();
	bool isFirst = true;
	for (std::vector<unsigned>::const_iterator oIt = origIndexVec.begin(); oIt != origIndexVec.end(); ++oIt)
		{
		const unsigned k = *oIt; // should invert this so we can walk down charPos
		const NxsString cLabel =  (k < nLabels ? this->charLabels[k] : mtString);
		const NxsStringVectorMap::const_iterator cib = this->charStates.find(k);
		if (cib != this->charStates.end())
			{
			const unsigned n = 1 + this->GetOrigCharIndex(k);
			const NxsStringVector & stateLabelsVec = cib->second;
			unsigned ns = (int)stateLabelsVec.size();
			out <<	(isFirst ? "CHARSTATELABELS " : ", ") ;
			out << n << ' ';
			if (cLabel.length() > 0)
				out << NxsToken::EscapeString(cLabel);
				
			out << " / ";
			for (unsigned m = 0; m < ns; m++)
				out << " " << NxsToken::EscapeString(stateLabelsVec[m]);
			isFirst = false;
			}
		else if (cLabel.length() > 0)
			{
			out <<	(isFirst ? "CHARSTATELABELS " : ", ") ;
			unsigned p = 1 + this->GetOrigCharIndex(k);
			out << p << NxsToken::EscapeString(cLabel);
			isFirst = false;
			}
		}
	if (!isFirst)
		out << ';';
	}

void NxsCharactersBlock::WriteFormatCommand(std::ostream &out) const
	{
	out << "FORMAT Datatype=" << this->GetDatatypeName();
	if (this->respectingCase)
		out << " RespectCase"; 
	if (this->missing != '?')
		{
		out << " Missing=";
		out << this->missing;
		}
	if (this->gap != '\0')
		{
		out << "  Gap=";
		out << this->gap;
		}
	if (this->datatype != NxsCharactersBlock::continuous)
		{
		unsigned numDefStates = 4;
		if (this->datatype == NxsCharactersBlock::protein)
			numDefStates = 21;
		else if (this->datatype == NxsCharactersBlock::standard)	
			numDefStates = 0;
		unsigned nSym = strlen(this->symbols);
		if (nSym > numDefStates)
			out << " Symbols=\"" << (this->symbols + numDefStates) <<"\"";
		}
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
	if (this->matchchar != '\0')
		{
		out << "  MatchChar=";
		out << this->matchchar;
		}
	if (this->datatype == continuous)
		{
		out << " Items = (";
		for (vector<std::string>::const_iterator iIt = items.begin(); iIt != items.end(); ++iIt)
			out << *iIt << ' '; 
		out << ")";
		if (this->statesFormat == STATES_PRESENT)
			out << " StatesFormat=StatesPresent";
		}
	else 
		{
		if (this->statesFormat == INDIVIDUALS)
			out << " StatesFormat=Individuals";
		}
	out << ";";
	}

NxsStringMap NxsCharactersBlock::GetDefaultEquates() const
	{
	NxsStringMap defEquates;
	const DataTypesEnum dt = this->datatype;
	if (dt == NxsCharactersBlock::dna || dt == NxsCharactersBlock::rna || dt == NxsCharactersBlock::nucleotide)
		{
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
	else if (dt == NxsCharactersBlock::protein)
		{
		defEquates[ NxsString("B") ] = NxsString("{DN}");
		defEquates[ NxsString("Z") ] = NxsString("{EQ}");
		}
	return defEquates;
	}

const char *  NxsCharactersBlock::GetDatatypeName() const
	{
	switch(datatype)
		{
		case NxsCharactersBlock::dna:
			return "DNA";
		case NxsCharactersBlock::rna:
			return "RNA";
		case NxsCharactersBlock::nucleotide:
			return "Nucleotide";
		case NxsCharactersBlock::protein:
			return "Protein";
		case NxsCharactersBlock::continuous:
			return "Continuous";
		default:
			return "Standard";
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns NxsCharactersBlock object to the state it was in when first created.
*/
void NxsCharactersBlock::Reset()
	{ 
	NxsBlock::Reset();
	ntax				= 0;
	ntaxTotal			= 0;
	nchar				= 0;
	ncharTotal			= 0;
	newchar				= true;
	newtaxa				= false;
	interleaving		= false;
	transposing			= false;
	respectingCase		= false;
	formerly_datablock	= false;
	labels				= true;
	tokens				= false;
	datatype			= NxsCharactersBlock::standard;
	missing				= '?';
	gap					= '\0';
	matchchar			= '\0';
	matrix				= NULL;
	charPos				= NULL;
	taxonPos			= NULL;
	activeTaxon			= NULL;
	activeChar			= NULL;
	symbols				= NULL;

	ResetSymbols();

	charLabels.clear();
	charStates.clear();
	equates.clear();
	eliminated.clear();

	if (matrix != NULL)
		{
		delete matrix;
		matrix = NULL;
		}

	if (charPos != NULL)
		{
		delete [] charPos;
		charPos = NULL;
		}

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

	if (activeChar != NULL)
		{
		delete [] activeChar;
		activeChar = NULL;
		}
	continuousCharMatrix.clear();
	items = std::vector<std::string>(1, std::string("STATES"));
	statesFormat = STATES_PRESENT;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets standard symbol set after a change in `datatype' is made. Also flushes equates list and installs standard 
|	equate macros for the current `datatype'.
*/
void NxsCharactersBlock::ResetSymbols()
	{
	// Symbols might be NULL (if a different NxsCharactersBlock has consumed this one
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
		case NxsCharactersBlock::dna:
			strcpy(symbols, "ACGT");
			break;

		case NxsCharactersBlock::rna:
			strcpy(symbols, "ACGU");
			break;

		case NxsCharactersBlock::nucleotide:
			strcpy(symbols, "ACGT");
			break;

		case NxsCharactersBlock::protein:
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
|	Looks up the state(s) at row `i', column `j' of matrix and writes it (or them) to out. If there is uncertainty or 
|	polymorphism, the list of states is surrounded by the appropriate set of symbols (i.e., parentheses for 
|	polymorphism, curly brackets for uncertainty). If TOKENS was specified, the output takes the form of the defined 
|	state labels; otherwise, the correct symbol is looked up in `symbols' and output.
*/
void NxsCharactersBlock::ShowStateLabels(
  ostream &out,				/* the output stream on which to write */
  unsigned i,				/* the taxon, in range [0..`ntax') */
  unsigned j,				/* the character, in range [0..`nchar') */
  unsigned first_taxon)		/* the index of the first taxon (if UINT_MAX, don't use matchchar) */
	{
	if (datatype == continuous)
		{
		const ContinuousCharCell & cell = continuousCharMatrix.at(i).at(j);
		std::vector<std::string>::const_iterator itemIt = items.begin();
		out	 << '(';

		for (; itemIt != items.end(); ++itemIt)
			{
			ContinuousCharCell::const_iterator cit = cell.find(*itemIt);
			if (cit == cell.end())
				{
				out << missing << ' ';
				}
			else
				{
				if (cit->second.empty())
					{
					out << missing << ' ';
					}
				else
					{
					vector<double>::const_iterator vIt = cit->second.begin();
					for(; vIt != cit->second.end(); ++vIt)
						{
						if (*vIt == DBL_MAX)
							out << missing << ' ';
						else
							out << *vIt << ' ';
						}
					}
					
				}
			}
		out	 << ") ";
		return;
		}
	if (tokens)
		{
		unsigned n = matrix->GetNumStates(i, j);
		if (n == 0 && matrix->IsGap(i, j))
			out << gap;
		else if (n == 0 && matrix->IsMissing(i, j))
			out << missing;
		else if (n == 1) 
			{
			unsigned s = matrix->GetState(i, j);
			const bool use_matchchar = ((i > first_taxon) && (matrix->GetState(first_taxon, j) == s));
			if (use_matchchar)
				out << "  .";
			else
				{
				NxsStringVectorMap::const_iterator ci = charStates.find(j);

				// OPEN ISSUE: need to eliminate state labels for characters that have
				// been eliminated
				//
				if (ci == charStates.end())
					out << "  " << s << "[<-no label found]";
				else
					{
					// Show label at index number s in NxsStringVector at ci
					//
					out << "  " << (*ci).second[s];
					}
				}	// if (use_matchchar) ... else
			}	// if (n == 0 && matrix->IsGap(i, j)) ... else if (n == 0 && matrix->IsMissing(i, j)) ... else
		else 
			{
			// TODO: handle matchchar possibility here too
			//
			if (matrix->IsPolymorphic(i, j))
				out << "  (";
			else
				out << "  {";
			for (unsigned k = 0; k < n; k++)
				{
				unsigned s = matrix->GetState(i, j, k);
				NxsStringVectorMap::const_iterator ci = charStates.find(j);
				if (ci == charStates.end())
					out << "  " << s << "[<-no label found]";
				else
					{
					// Show label at index number s in NxsStringVector at ci
					//
					out << "  " << (*ci).second[s];
					}
				}
			if (matrix->IsPolymorphic(i, j))
				out << ')';
			else
				out << '}';
			}

		}	// if (tokens)

	else
		{
		if (i > first_taxon) 
			{
			char s[NCL_MAX_STATES + 3];
			WriteStates(matrix->GetDiscreteDatum(i, j), s, NCL_MAX_STATES + 3);

			char ss[NCL_MAX_STATES + 3];
			WriteStates(matrix->GetDiscreteDatum(first_taxon, j), ss, NCL_MAX_STATES + 3);

			if (strcmp(s, ss) == 0)
				out << '.';
			else
				ShowStates(out, i, j);
			}
		else 
			ShowStates(out, i, j);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Writes out the state (or states) stored in this NxsDiscreteDatum object to the buffer `s' using the symbols array 
|	to do the necessary translation of the numeric state values to state symbols. In the case of polymorphism or 
|	uncertainty, the list of states will be surrounded by brackets or parentheses (respectively). Assumes `s' is 
|	non-NULL and long enough to hold everything printed.
*/
void NxsCharactersBlock::WriteStates(
  NxsDiscreteDatum &d,	/* the datum to be queried */
  char *s,				/* the buffer to which to print */
  unsigned slen)		/* the length of the buffer `s' */
	{
	assert(s != NULL);
	assert(slen > 1);

	if (matrix->IsMissing(d))
		{
		s[0] = missing;
		s[1] = '\0';
		}
	else if (matrix->IsGap(d))
		{
		s[0] = gap;
		s[1] = '\0';
		}
	else
		{
		assert(symbols != NULL);
		unsigned symbolListLen = (unsigned)strlen(symbols);

		unsigned numStates = matrix->GetNumStates(d);
		unsigned numCharsNeeded = numStates;
		if (numStates > 1)
			numCharsNeeded += 2;
		assert(slen > numCharsNeeded);

		if (numStates == 1)
			{
			unsigned v = matrix->GetState(d);
			assert(v < symbolListLen);
			s[0] = symbols[v];
			s[1] = '\0';
			}

		else
			{
			// numStates must be greater than 1
			//
			unsigned i = 0;
			if (matrix->IsPolymorphic(d))
				s[i++] = '(';
			else
				s[i++] = '{';
			for (unsigned k = 0; k < numStates; k++)
				{
				unsigned v = matrix->GetState(d, k);
				assert(v < symbolListLen);
				s[i++] = symbols[v];
				s[i] = '\0';
				}
			if (matrix->IsPolymorphic(d))
				s[i++] = ')';
			else
				s[i++] = '}';
			s[i] = '\0';
			}
		}
	}
