//  Copyright (C) 2007-2008 Brian O'Meara & Derrick Zwickl
//  A modification of the BasicCommandLine file of the NCL (see below)
//  to use for loading trees and data from Nexus into R. Licensing as below.

//	Copyright (C) 1999-2002 Paul O. Lewis
//
//	This file is part of NCL (Nexus Class Library).
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

#include <ncl/ncl.h>
#include "NCLInterface.h"

#include <iostream>
#include <sstream>
#include <cassert>



std::string BASICCMDLINE::TestRunning() {
	string output="The NCL has been loaded";
	return output;
}

#if defined(NEW_NCL_INTERFACE)

MultiFormatReader * instantiateReader(int gStrictLevel);



void BASICCMDLINE::AppendRContent(
  NxsString & nexuscharacters,
  NxsCharactersBlock &characters,
  bool allchar,
  bool polymorphictomissing,
  bool levelsall,
  bool returnlabels) const
{

	int nchartoreturn = 0;
	NxsTaxaBlockAPI * taxa = characters.GetTaxaBlockPtr();
	if (!taxa)
		return;
	int ntax = taxa->GetNumTaxonLabels();

	if (!characters.IsEmpty())
		{
		nexuscharacters+="|";
		nexuscharacters+=characters.GetNameOfDatatype(characters.GetDataType());
		nexuscharacters+="|";

		if (std::string("Standard")==characters.GetNameOfDatatype(characters.GetDataType())) {
			nexuscharacters+="data.frame(";

			if (allchar) {
				nchartoreturn=characters.GetNCharTotal();
			}
			else {
				nchartoreturn=characters.GetNumIncludedChars();
			}
			for (int character=0; character<nchartoreturn; character++) { //We only pass the non-eliminated chars
				NxsString charlabel=characters.GetCharLabel(character);
				if (character>0) {
					nexuscharacters+=", ";
				}
				if (charlabel.length()>1) {
					nexuscharacters+="'";
					nexuscharacters+=charlabel;
					nexuscharacters+="'";
				}
				else {
					nexuscharacters+="standard_char_";
					nexuscharacters+=character+1;
				}
				nexuscharacters+=" = ";

				nexuscharacters+="factor(c(";
				for (int taxon=0;taxon<ntax;taxon++) {
					int statenumber=characters.GetInternalRepresentation(taxon,character,0);

					if(characters.IsMissingState(taxon,character)) {
						nexuscharacters+="NA";
					}
					else if (characters.GetNumStates(taxon,character)>1) {
						if(polymorphictomissing) {
							nexuscharacters+="NA";
						}
						else {
							nexuscharacters+='"';
							nexuscharacters+='{';
							for (unsigned int k=0;k<characters.GetNumStates(taxon,character);k++) {
								nexuscharacters+=characters.GetInternalRepresentation(taxon,character,k);	
								if (k+1<characters.GetNumStates(taxon,character)) {
									nexuscharacters+=',';
								}
							}
							nexuscharacters+='}';
							nexuscharacters+='"';
						}
					}
					else {
						nexuscharacters+='"';
						nexuscharacters+=statenumber;
						nexuscharacters+='"';
					}
					if (taxon+1<ntax) {
						nexuscharacters+=',';
					}
				}
				nexuscharacters+=')';
				if (returnlabels) {
					NxsString labels=", labels=c(";
					unsigned int totallabellength=0;
					for (unsigned int l=0;l<characters.GetObsNumStates(character); l++) {
						labels+='"';
						labels+= characters.GetStateLabel(character,l);
						totallabellength+=(characters.GetStateLabel(character,l)).length();
						labels+='"';
						if (l+1<characters.GetObsNumStates(character)) {
								labels+=',';
						}
					}
					labels+=')';
					if (totallabellength>characters.GetObsNumStates(character)) {
						nexuscharacters+=labels;
					}
				}
				nexuscharacters+=")\n";
			}
			nexuscharacters+=", row.names=c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				nexuscharacters+=RemoveUnderscoresAndSpaces(characters.GetTaxonLabel(taxon));
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}
			nexuscharacters+="))";
		}
		else if (std::string("DNA")==characters.GetNameOfDatatype(characters.GetDataType()) ||
		         std::string("RNA")==characters.GetNameOfDatatype(characters.GetDataType()) ||
			 std::string("Nucleotide")==characters.GetNameOfDatatype(characters.GetDataType())) {

			nexuscharacters+="data.frame(";
			if (2==characters.GetDataType()) {
				nexuscharacters+="dna_alignment_1=c(";
			}
			else if (3==characters.GetDataType()) {
				nexuscharacters+="rna_alignment_1=c(";
			}
			else {
				nexuscharacters+="nucleotide_alignment_1=c(";
			}

			if (allchar) {
				nchartoreturn=characters.GetNCharTotal();
			}
			else {
				nchartoreturn=characters.GetNumIncludedChars();
			}
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				for (int character=0; character<nchartoreturn; character++) {
					int numstates=characters.GetNumStates(taxon,character);
					if (characters.IsGapState(taxon,character)) {
						nexuscharacters+="-";
					}
					else if (characters.IsMissingState(taxon,character)) {
						nexuscharacters+="?";
					}
					else if ( numstates == 1) {
						nexuscharacters+=characters.GetState(taxon,character,0);
					}
					else {
						bool hasA=false;
						bool hasT=false;
						bool hasG=false;
						bool hasC=false;
						for (int i=0;i<numstates;i++) {
							unsigned currentstate=characters.GetState(taxon,character,i);
							if (currentstate=='A') {
								hasA=true;
							}
							else if (currentstate=='T') {
								hasT=true;
							}
							else if (currentstate=='G') {
								hasG=true;
							}
							else if (currentstate=='C') {
								hasC=true;
							}
							else {
								cout<<"Error: currentstate "<<currentstate<<" does not match ATGC"<<endl;
							}
						}
						if (hasA && hasG && !hasT && !hasC) {
							nexuscharacters+="R";
						}
						else if (!hasA && !hasG && hasC && hasT) {
							nexuscharacters+="Y";
						}
						else if (hasA && !hasG && hasC && !hasT) {
							nexuscharacters+="M";
						}
						else if (!hasA && hasG && !hasC && hasT) {
							nexuscharacters+="K";
						}
						else if (!hasA && hasG && hasC && !hasT) {
							nexuscharacters+="S";
						}
						else if (hasA && !hasG && !hasC && hasT) {
							nexuscharacters+="W";
						}
						else if (hasA && !hasG && hasC && hasT) {
							nexuscharacters+="H";
						}
						else if (!hasA && hasG && hasC && hasT) {
							nexuscharacters+="B";
						}
						else if (hasA && hasG && hasC && !hasT) {
							nexuscharacters+="V";
						}
						else if (hasA && hasG && !hasC && hasT) {
							nexuscharacters+="D";
						}
						else if (hasA && hasG && hasC && hasT) {
							nexuscharacters+="N";
						}
						else {
							nexuscharacters+="N";
						}
					}
				}
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}

			nexuscharacters+="), row.names=c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				nexuscharacters+=RemoveUnderscoresAndSpaces(characters.GetTaxonLabel(taxon));
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}
			nexuscharacters+="), stringsAsFactors=FALSE)";

		}
		else if (std::string("Protein")==characters.GetNameOfDatatype(characters.GetDataType())) {
			nexuscharacters+="data.frame(";
			nexuscharacters+="aa_alignment_1=c(";


			if (allchar) {
				nchartoreturn=characters.GetNCharTotal();
			}
			else {
				nchartoreturn=characters.GetNumIncludedChars();
			}
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				for (int character=0; character<nchartoreturn; character++) {
					int numstates=characters.GetNumStates(taxon,character);
					if (characters.IsGapState(taxon,character)) {
						nexuscharacters+="-";
					}
					else if (characters.IsMissingState(taxon,character)) {
						nexuscharacters+="?";
					}
					else if ( numstates == 1) {
						nexuscharacters+=characters.GetState(taxon,character,0);
					}
					else {
						nexuscharacters+="X";
					}
				}
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}

			nexuscharacters+="), row.names=c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				nexuscharacters+=RemoveUnderscoresAndSpaces(characters.GetTaxonLabel(taxon));
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}
			nexuscharacters+="), stringsAsFactors=FALSE)";
		}
		else if (std::string("Continuous")==characters.GetNameOfDatatype(characters.GetDataType())) { 
			nexuscharacters+="data.frame(";

			if (allchar) {
				nchartoreturn=characters.GetNCharTotal();
			}
			else {
				nchartoreturn=characters.GetNumIncludedChars();
			}
			for (int character=0; character<nchartoreturn; character++) { //We only pass the non-eliminated chars
				NxsString charlabel=characters.GetCharLabel(character);
				if (character>0) {
					nexuscharacters+=", ";
				}
				if (charlabel.length()>1) {
					nexuscharacters+="'";
					nexuscharacters+=charlabel;
					nexuscharacters+="'";
				}
				else {
					nexuscharacters+="continuous_char_";
					nexuscharacters+=character+1;
				}
				nexuscharacters+=" = ";

				nexuscharacters+="c(";
				for (int taxon=0;taxon<ntax;taxon++) {
					double state=characters.GetSimpleContinuousValue(taxon,character);
					//cout<<"State at "<<taxon+1<<" char "<<character+1<<" = "<<state<<endl;
					if (state==DBL_MAX) {
							nexuscharacters+="NA";
					}
					else {					    
					    char buffer[100];
					    sprintf(buffer, "%.10f", state);
					    nexuscharacters+=buffer; 
					}

					if (taxon+1<ntax) {
						nexuscharacters+=',';
					}
				}
				nexuscharacters+=")";
			}
			nexuscharacters+=", row.names=c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				nexuscharacters+=RemoveUnderscoresAndSpaces(characters.GetTaxonLabel(taxon));
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}
			nexuscharacters+="))";	      		
		}
		else {
			std::string message="Matrix loaded but datatype: ";
			message+=characters.GetNameOfDatatype(characters.GetDataType());
			message+=" is not supported (yet)"; //"Error: character matrix loaded, but does not match any category (dna, standard, etc.)";
			errorMessage(message);
		}
		nexuscharacters=RemoveUnderscoresAndSpaces(nexuscharacters);
	}

}

void BASICCMDLINE::AppendRContent(
  NxsString & nexustrees,
  NxsTreesBlock &treesB) const
{
	NxsTaxaBlockAPI * taxa = treesB.GetTaxaBlockPtr();
	if (!taxa)
		return;
	if (!treesB.IsEmpty())
	{

		nexustrees+= "\nBEGIN TREES;\n";
		for (unsigned k = 0; k < treesB.GetNumTrees(); k++)
		{
			NxsString s = treesB.GetTreeName(k);
			s.BlanksToUnderscores();
			nexustrees+="\tTREE ";
			nexustrees+=s;
			nexustrees+=" = ";
			if (treesB.IsRootedTree(k))
				nexustrees+="[&R]";
			else
				nexustrees+="[&U]";
			nexustrees+=treesB.GetTranslatedTreeDescription(k);
			nexustrees+=";\n";
		}
		nexustrees+="END;\n";
	}
}

void BASICCMDLINE::AppendRContent(
  NxsString & nexusdistances,
  NxsDistancesBlock &distancesB) const
{
	NxsTaxaBlockAPI * taxa = distancesB.GetTaxaBlockPtr();
	if (!taxa)
		return;
	int ntax = taxa->GetNumTaxonLabels();

	if (!distancesB.IsEmpty())
	{ //fill cols first, first col has taxon 1, first row has taxon 2 (no diags)
		nexusdistances+="\ndistances <- structure(c(";
		vector<double> distancevector;
		for (int col=0;col<ntax-1;col++) {
			for (int row=col+1;row<ntax;row++) {
				distancevector.push_back(distancesB.GetDistance(row,col));
			}
		}
		for (unsigned int i=0;i<distancevector.size();i++) {
			nexusdistances+=distancevector.at(i);
			if (i+1<distancevector.size()) {
				nexusdistances+=',';
			}
		}
		nexusdistances+="), Size = ";
		nexusdistances+=ntax;
		nexusdistances+="L, Labels = c(";
		for (int taxon=0; taxon<ntax;taxon++) {
			nexusdistances+='"';
			nexusdistances+=RemoveUnderscoresAndSpaces(taxa->GetTaxonLabel(taxon));
			nexusdistances+='"';
			if (taxon+1<ntax) {
				nexusdistances+=", ";
			}
		}
		nexusdistances+="), Upper = FALSE, Diag = FALSE, class = \"dist\")\n";
	}
}


void BASICCMDLINE::RReturnCharacters(NxsString & outputstring, bool allchar, bool polymorphictomissing, bool levelsall, bool returnlabels)
{
	if (!nexusReader)
		return;
	const unsigned nTaxaBlocks = nexusReader->GetNumTaxaBlocks();
	for (unsigned t = 0; t < nTaxaBlocks; ++t) {
		const NxsTaxaBlock * tb = nexusReader->GetTaxaBlock(t);
		const unsigned nCharBlocks = nexusReader->GetNumCharactersBlocks(tb);
		if (nCharBlocks == 0)
			continue;
		for (unsigned i = 0; i < nCharBlocks; ++i) {
			NxsCharactersBlock * cb = nexusReader->GetCharactersBlock(tb, i);
			if (cb) {
				AppendRContent(outputstring, *cb, allchar, polymorphictomissing, levelsall, returnlabels);
			}
		}
	}
}


void BASICCMDLINE::RReturnTrees(NxsString & outputstring)
{
	if (!nexusReader)
		return;
	const unsigned nTaxaBlocks = nexusReader->GetNumTaxaBlocks();
	for (unsigned t = 0; t < nTaxaBlocks; ++t) {
		const NxsTaxaBlock * tb = nexusReader->GetTaxaBlock(t);
		const unsigned nTreesBlocks = nexusReader->GetNumTreesBlocks(tb);
		if (nTreesBlocks == 0)
			continue;
		for (unsigned i = 0; i < nTreesBlocks; ++i) {
			NxsTreesBlock * treesb = nexusReader->GetTreesBlock(tb, i);
			if (treesb) {
				AppendRContent(outputstring, *treesb);
			}
		}
	}
}

void BASICCMDLINE::RReturnDistances(NxsString & outputstring)
{
	if (!nexusReader)
		return;
	const unsigned nTaxaBlocks = nexusReader->GetNumTaxaBlocks();
	for (unsigned t = 0; t < nTaxaBlocks; ++t) {
		const NxsTaxaBlock * tb = nexusReader->GetTaxaBlock(t);
		const unsigned nDistancesBlocks = nexusReader->GetNumDistancesBlocks(tb);
		if (nDistancesBlocks == 0)
			continue;
		for (unsigned i = 0; i < nDistancesBlocks; ++i) {
			NxsDistancesBlock * distancesb = nexusReader->GetDistancesBlock(tb, i);
			if (distancesb) {
				AppendRContent(outputstring, *distancesb);
			}
		}
	}
}



void BASICCMDLINE::writeMessage(MessageLevel level, const std::string & m)  const {
	if (level >= int(this->verboseLevel)) {
		std::cerr << m << std::endl;
	}
}


MultiFormatReader * instantiateReader(int gStrictLevel)
{
	MultiFormatReader * nexusReader = new MultiFormatReader(-1, NxsReader::WARNINGS_TO_STDERR);


	if (gStrictLevel != 2)
		nexusReader->SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
	NxsCharactersBlock * charsB = nexusReader->GetCharactersBlockTemplate();
	NxsDataBlock * dataB = nexusReader->GetDataBlockTemplate();
	charsB->SetAllowAugmentingOfSequenceSymbols(true);
	dataB->SetAllowAugmentingOfSequenceSymbols(true);

	NxsTreesBlock * treesB = nexusReader->GetTreesBlockTemplate();
	assert(treesB);
	if (gStrictLevel < 2)
		treesB->SetAllowImplicitNames(true);
	if (gStrictLevel < 2)
		{
		NxsStoreTokensBlockReader *storerB =  nexusReader->GetUnknownBlockTemplate();
		assert(storerB);
		storerB->SetTolerateEOFInBlock(true);
		}
	nexusReader->conversionOutputRecord.addNumbersToDisambiguateNames = true;
	return nexusReader;
}


void BASICCMDLINE::Clear()
{
	if (nexusReader) {
		MultiFormatReader * nr = nexusReader;
		nexusReader = 0L;
		try {
			nr->DeleteBlocksFromFactories();
		}
		catch (...) {
			delete nr;
			throw;
		}
		delete nr;
	}

}



////////////////////////////////////////////////////////////////////////////////
// Creates a NxsReader, and tries to read the file `filename`.  If the
//	read succeeds, then processContent will be called.
////////////////////////////////////////////////////////////////////////////////
void BASICCMDLINE::processFilepath(const char * filename) // enum indicating the file format to expect.
{
	if (! (this->nexusReader)) {
		nexusReader = instantiateReader(this->gStrictLevel);
	}

	nexusReader->SetWarningOutputLevel(this->warningLevel);

	std::string m("Reading ");
	if (filename)
		m.append(filename);
	this->statusMessage(m);

	nexusReader->ReadFilepath(filename, fmt);

	m.assign("Finished reading ");
	if (filename)
		m.append(filename);
	statusMessage(m);


}







#else // #if defined(NEW_NCL_INTERFACE)



/*----------------------------------------------------------------------------------------------------------------------
 |	The constructor simply passes along `i' to the base class constructor. Nothing else needs to be done.
 */
MyNexusToken::MyNexusToken(
						   istream &i)	/* the input file stream attached to the NEXUS file to be read */
: NxsToken(i)
{
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Overrides the NxsToken::OutputComment virtual function (which does nothing) to display output comments [!comments
 |	like this one beginning with an exclamation point]. The output comment contained in `msg' is simply sent to the
 |	standard output stream cout.
 */
void MyNexusToken::OutputComment(
								 const NxsString &msg)	/* the output comment to be displayed */
{
	cout << msg << endl;
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Initializes the `id' data member to "BASICCMDLINE" and calls the FactoryDefaults member function to perform the
 |	remaining initializations. The data member `'' is set to NULL so that memory will be allocated for it in
 |	FactoryDefaults.
 */
BASICCMDLINE::BASICCMDLINE()
{
	id				= "BASICCMDLINE";

	// Make sure all data members that are pointers are initialized to NULL!
	// Failure to do this will result in problems because functions such as
	// FactoryDefaults() will try to delete an object if it is non-NULL.
	//
	taxa			= NULL;
	trees			= NULL;
	assumptions		= NULL;
	distances		= NULL;
	characters		= NULL;
	data			= NULL;
	next_command	= NULL;

	FactoryDefaults();

	taxa			= new NxsTaxaBlock();
	trees			= new NxsTreesBlock(taxa);
	assumptions		= new NxsAssumptionsBlock(taxa);
	characters		= new NxsCharactersBlock(taxa, assumptions);
	distances		= new NxsDistancesBlock(taxa);
	data			= new NxsDataBlock(taxa, assumptions);

	Add(taxa);
	Add(trees);
	Add(assumptions);
	Add(characters);
	Add(distances);
	Add(data);
	Add(this);

}

/*----------------------------------------------------------------------------------------------------------------------
 |	Closes `logf' if it is open and deletes memory allocated to `next_command'.
 */
BASICCMDLINE::~BASICCMDLINE()
{
	assert(next_command != NULL);
	delete [] next_command;

	if (logf_open)
		logf.close();
}

/*----------------------------------------------------------------------------------------------------------------------
 |	The code here is identical to the base class version (simply returns 0), so the code here should either be
 |	modified or this derived version eliminated altogether. Under what circumstances would you need to modify the
 |	default code, you ask? This function should be modified to something meaningful if this derived class needs to
 |	construct and run a NxsSetReader object to read a set involving characters. The NxsSetReader object may need to
 |	use this function to look up a character label encountered in the set. A class that overrides this method should
 |	return the character index in the range [1..`nchar']; i.e., add one to the 0-offset index.
 */
unsigned BASICCMDLINE::CharLabelToNumber(
										 NxsString)	/* the character label to be translated to character number */
{
	return 0;
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called by the NxsReader object when a block named `blockName' is entered. Allows program to notify user of
 |	progress in parsing the NEXUS file. Also gives program the opportunity to ask user if it is ok to purge data
 |	currently contained in this block. If user is asked whether existing data should be deleted, and the answer comes
 |	back no, then then return false, otherwise return true. Overrides pure virtual function in class NxsReader.
 */
bool BASICCMDLINE::EnteringBlock(
								 NxsString blockName)	/* the name of the block just entered */
{
	message = "Reading ";
	message += blockName;
	message += " block...";
	PrintMessage();

	return true;
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called by the NxsReader object when exiting a block named `blockName'. Allows program to notify user of progress
 |	in parsing the NEXUS file. Virtual function that overrides the pure virtual function in the base class NxsReader.
 */
void BASICCMDLINE::ExitingBlock(
								NxsString )	/* the name of the block just exited */
{
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Sets all data members to their factory default settings: `inf_open', `logf_open' and `quit_now' are set to false;
 |	`message' to the null string, and the pointers `data', `characters', `assumptions', `taxa' and `trees'
 |	are all set to NULL. The C-string `next_command' is allocated COMMAND_MAXLEN + 1 bytes if it is currently NULL,
 |	and its first byte is set to the null character to create an empty `next_command' string.
 */
void BASICCMDLINE::FactoryDefaults()
{
	isEmpty = true;
	inf_open = false;
	logf_open = false;
	quit_now = false;
	message.clear();

	if (trees != NULL)
	{
		Detach(trees);
		delete trees;
		trees = NULL;
	}

	if (taxa != NULL)
	{
		Detach(taxa);
		delete taxa;
		taxa = NULL;
	}

	if (assumptions != NULL)
	{
		Detach(assumptions);
		delete assumptions;
		assumptions = NULL;
	}

	if (distances != NULL)
	{
		Detach(distances);
		delete distances;
		distances = NULL;
	}

	if (characters != NULL)
	{
		Detach(characters);
		delete characters;
		characters = NULL;
	}

	if (data != NULL)
	{
		Detach(data);
		delete data;
		data = NULL;
	}

	if (next_command == NULL)
		next_command = new char[COMMAND_MAXLEN + 1];
	next_command[0] = '\0';
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Returns true if file named `fn' already exists, false otherwise.
 */
bool BASICCMDLINE::FileExists(
							  const char *fn)	/* the name of the file to check */
{
	bool exists = false;

	FILE *fp = fopen(fn, "r");
	if (fp != NULL)
	{
		fclose(fp);
		exists = true;
	}

	return exists;
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called whenever a file name needs to be read from either the command line or a file. Expects next token to be "="
 |	followed by the token representing the file name. Call this function after, say, the keyword "file" has been read
 |	in the following LOG command:
 |>
 |	log file=doofus.txt start replace;
 |>
 |	Note that this function will read only the "=doofus.txt " leaving "start replace;" in the stream for reading at
 |	a later time.
 */
NxsString BASICCMDLINE::GetFileName(
									NxsToken &token)	/* the token used to read from `in' */
{
	// Eat the equals sign
	//
	token.GetNextToken();

	if (!token.Equals("="))
	{
		errormsg = "Expecting an equals sign, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

	// Now get the filename itself
	//
	token.GetNextToken();

	return token.GetToken();
}

void BASICCMDLINE::HandleReturnData(
									NxsToken &token)	/*  the token used to read from `in' */
{
	// Get the semicolon following END or ENDBLOCK token
	//
	NxsString null=ReturnDataForR(false, true, false);
	token.GetNextToken();

	if (!token.Equals(";"))
	{
		errormsg = "Expecting ';' to terminate the RETURNDATA command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}
}

void BASICCMDLINE::RReturnCharacters(NxsString & nexuscharacters, bool allchar, bool polymorphictomissing, bool levelsall) {

	int nchartoreturn=0;
	int ntax = taxa->GetNumTaxonLabels();
	if (!characters->IsEmpty())
	{
		//characters->Report(cerr);
		if (1==characters->GetDataType()) { //standard datatype
		//if((characters->GetDatatypeName())=="standard") {
			//nexuscharacters=characters->GetDatatypeName();
			nexuscharacters+="data.frame(";

			if (allchar) {
				nchartoreturn=characters->GetNCharTotal();
			}
			else {
				nchartoreturn=characters->GetNumIncludedChars();
			}
			for (int character=0; character<nchartoreturn; character++) { //We only pass the non-eliminated chars
				NxsString charlabel=characters->GetCharLabel(character);
				if (character>0) {
					nexuscharacters+=", ";
				}
				if (charlabel.length()>1) {
					nexuscharacters+="'";
					nexuscharacters+=charlabel;
					nexuscharacters+="'";
				}
				else {
					nexuscharacters+="standard_char_";
					nexuscharacters+=character+1;
				}
				nexuscharacters+=" = ";

				nexuscharacters+="factor(c(";
				for (int taxon=0;taxon<ntax;taxon++) {
					int statenumber=characters->GetInternalRepresentation(taxon,character,0);

					if(characters->IsMissingState(taxon,character)) {
						nexuscharacters+="NA";
					}
					else if (characters->GetNumStates(taxon,character)>1) {
						if(polymorphictomissing) {
							nexuscharacters+="NA";
						}
						else {
							nexuscharacters+='{';
							for (unsigned int k=0;k<characters->GetNumStates(taxon,character);k++) {
								nexuscharacters+=characters->GetInternalRepresentation(taxon,character,0);
								if (k+1<characters->GetNumStates(taxon,character)) {
									nexuscharacters+=',';
								}
							}
							nexuscharacters+='}';
						}
					}
					else {
						nexuscharacters+=statenumber;
					}
					if (taxon+1<ntax) {
						nexuscharacters+=',';
					}
				}
				nexuscharacters+=')';
				if (levelsall) {
					nexuscharacters+=", levels=c(";
					for (unsigned int l=0;l<characters->GetMaxObsNumStates(); l++) {
						nexuscharacters+=l;
						if (l+1<characters->GetMaxObsNumStates()) {
							nexuscharacters+=',';
						}
					}
					nexuscharacters+=')';
				}
				else {

					NxsString levels=", levels=c(";
					NxsString labels=", labels=c(";
					unsigned int totallabellength=0;
					for (unsigned int l=0;l<characters->GetObsNumStates(character); l++) {
						labels+='"';
						labels+= characters->GetStateLabel(character,l);
						totallabellength+=(characters->GetStateLabel(character,l)).length();
						labels+='"';
						levels+= l;
						if (l+1<characters->GetObsNumStates(character)) {
							labels+=',';
							levels+=',';
						}
					}
					levels+=')';
					labels+=')';
					//cout<<"labels.length="<<labels.length()<<endl<<"levels.length="<<levels.length()<<endl<<"total label length="<<totallabellength<<endl;
					if (totallabellength>characters->GetObsNumStates(character)) {
						nexuscharacters+=levels;
						nexuscharacters+=labels;
					}
				}
				nexuscharacters+=")\n";
			}
			nexuscharacters+=", row.names=c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				nexuscharacters+=RemoveUnderscoresAndSpaces(characters->GetTaxonLabel(taxon));
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}
			nexuscharacters+="))";
		}
		else if (2==characters->GetDataType() || 3==characters->GetDataType() || 4==characters->GetDataType()) { //dna, rna, nucleotide
	//	if((characters->GetDatatypeName())=="dna") {

			nexuscharacters+="data.frame(";
			if (2==characters->GetDataType()) {
				nexuscharacters+="dna_alignment_1=c(";
			}
			else if (3==characters->GetDataType()) {
				nexuscharacters+="rna_alignment_1=c(";
			}
			else {
				nexuscharacters+="nucleotide_alignment_1=c(";
			}

			if (allchar) {
				nchartoreturn=characters->GetNCharTotal();
			}
			else {
				nchartoreturn=characters->GetNumIncludedChars();
			}
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				for (int character=0; character<nchartoreturn; character++) {
					int numstates=characters->GetNumStates(taxon,character);
					if (characters->IsGapState(taxon,character)) {
						nexuscharacters+="-";
					}
					else if (characters->IsMissingState(taxon,character)) {
						nexuscharacters+="?";
					}
					else if ( numstates == 1) {
						nexuscharacters+=characters->GetState(taxon,character,0);
					}
					else {
						bool hasA=false;
						bool hasT=false;
						bool hasG=false;
						bool hasC=false;
						for (int i=0;i<numstates;i++) {
							unsigned currentstate=characters->GetState(taxon,character,i);
							if (currentstate=='A') {
								hasA=true;
							}
							else if (currentstate=='T') {
								hasT=true;
							}
							else if (currentstate=='G') {
								hasG=true;
							}
							else if (currentstate=='C') {
								hasC=true;
							}
							else {
								cout<<"Error: currentstate "<<currentstate<<" does not match ATGC"<<endl;
							}
						}
						if (hasA && hasG && !hasT && !hasC) {
							nexuscharacters+="R";
						}
						else if (!hasA && !hasG && hasC && hasT) {
							nexuscharacters+="Y";
						}
						else if (hasA && !hasG && hasC && !hasT) {
							nexuscharacters+="M";
						}
						else if (!hasA && hasG && !hasC && hasT) {
							nexuscharacters+="K";
						}
						else if (!hasA && hasG && hasC && !hasT) {
							nexuscharacters+="S";
						}
						else if (hasA && !hasG && !hasC && hasT) {
							nexuscharacters+="W";
						}
						else if (hasA && !hasG && hasC && hasT) {
							nexuscharacters+="H";
						}
						else if (!hasA && hasG && hasC && hasT) {
							nexuscharacters+="B";
						}
						else if (hasA && hasG && hasC && !hasT) {
							nexuscharacters+="V";
						}
						else if (hasA && hasG && !hasC && hasT) {
							nexuscharacters+="D";
						}
						else if (hasA && hasG && hasC && hasT) {
							nexuscharacters+="N";
						}
						else {
							nexuscharacters+="N";
						}
					}
				}
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}

			nexuscharacters+="), row.names=c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				nexuscharacters+=RemoveUnderscoresAndSpaces(characters->GetTaxonLabel(taxon));
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}
			nexuscharacters+="), stringsAsFactors=FALSE)";

		}
		else if (5==characters->GetDataType()) { //protein
			nexuscharacters+="data.frame(";
			nexuscharacters+="aa_alignment_1=c(";


			if (allchar) {
				nchartoreturn=characters->GetNCharTotal();
			}
			else {
				nchartoreturn=characters->GetNumIncludedChars();
			}
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				for (int character=0; character<nchartoreturn; character++) {
					int numstates=characters->GetNumStates(taxon,character);
					if (characters->IsGapState(taxon,character)) {
						nexuscharacters+="-";
					}
					else if (characters->IsMissingState(taxon,character)) {
						nexuscharacters+="?";
					}
					else if ( numstates == 1) {
						nexuscharacters+=characters->GetState(taxon,character,0);
					}
					else {
						nexuscharacters+="X";
					}
				}
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}

			nexuscharacters+="), row.names=c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				nexuscharacters+=RemoveUnderscoresAndSpaces(characters->GetTaxonLabel(taxon));
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}
			nexuscharacters+="), stringsAsFactors=FALSE)";
		}
		else if (6==characters->GetDataType()) { //continuousnexuscharacters+="data.frame(";
			nexuscharacters+="data.frame(";

			if (allchar) {
				nchartoreturn=characters->GetNCharTotal();
			}
			else {
				nchartoreturn=characters->GetNumIncludedChars();
			}
			for (int character=0; character<nchartoreturn; character++) { //We only pass the non-eliminated chars
				NxsString charlabel=characters->GetCharLabel(character);
				if (character>0) {
					nexuscharacters+=", ";
				}
				if (charlabel.length()>1) {
					nexuscharacters+="'";
					nexuscharacters+=charlabel;
					nexuscharacters+="'";
				}
				else {
					nexuscharacters+="continuous_char_";
					nexuscharacters+=character+1;
				}
				nexuscharacters+=" = ";

				nexuscharacters+="c(";
				for (int taxon=0;taxon<ntax;taxon++) {
					double state=characters->GetSimpleContinuousValue(taxon,character);
					//cout<<"State at "<<taxon+1<<" char "<<character+1<<" = "<<state<<endl;
					if (state==DBL_MAX) {
							nexuscharacters+="NA";
					}
					else {
						nexuscharacters+=state;
					}

					if (taxon+1<ntax) {
						nexuscharacters+=',';
					}
				}
				nexuscharacters+=")";
			}
			nexuscharacters+=", row.names=c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				nexuscharacters+='"';
				nexuscharacters+=RemoveUnderscoresAndSpaces(characters->GetTaxonLabel(taxon));
				nexuscharacters+='"';
				if (taxon+1<ntax) {
					nexuscharacters+=',';
				}
			}
			nexuscharacters+="))";
			//message="Warning: Continuous characters do not work";
			//PrintMessage();
		}
		else {
				message="Error: character matrix loaded, but does not match any category (dna, standard, etc.)";
			PrintMessage();
		}
		nexuscharacters=RemoveUnderscoresAndSpaces(nexuscharacters);
	}
}

void BASICCMDLINE::RReturnTrees(NxsString & nexustrees) {
	if (!trees->IsEmpty())
	{

		nexustrees+= "\nBEGIN TREES;\n";
		for (unsigned k = 0; k < trees->GetNumTrees(); k++)
		{
			NxsString s = trees->GetTreeName(k);
			s.BlanksToUnderscores();
			nexustrees+="\tTREE ";
			nexustrees+=s;
			nexustrees+=" = ";
			if (trees->IsRootedTree(k))
				nexustrees+="[&R]";
			else
				nexustrees+="[&U]";
			nexustrees+=trees->GetTranslatedTreeDescription(k);
			nexustrees+=";\n";
		}
		nexustrees+="END;\n";
	}

}

void BASICCMDLINE::RReturnDistances(NxsString  & nexusdistances) {
	int ntax = taxa->GetNumTaxonLabels();

	if (!distances->IsEmpty())
	{ //fill cols first, first col has taxon 1, first row has taxon 2 (no diags)
		nexusdistances+="\ndistances <- structure(c(";
		vector<double> distancevector;
		for (int col=0;col<ntax-1;col++) {
			for (int row=col+1;row<ntax;row++) {
				distancevector.push_back(distances->GetDistance(row,col));
			}
		}
		for (unsigned int i=0;i<distancevector.size();i++) {
			nexusdistances+=distancevector.at(i);
			if (i+1<distancevector.size()) {
				nexusdistances+=',';
			}
		}
		nexusdistances+="), Size = ";
		nexusdistances+=ntax;
		nexusdistances+="L, Labels = c(";
		for (int taxon=0; taxon<ntax;taxon++) {
			nexusdistances+='"';
			nexusdistances+=RemoveUnderscoresAndSpaces(taxa->GetTaxonLabel(taxon));
			nexusdistances+='"';
			if (taxon+1<ntax) {
				nexusdistances+=", ";
			}
		}
		nexusdistances+="), Upper = FALSE, Diag = FALSE, class = \"dist\")\n";
	}
}

//Break into separate functions,Input is string reference rather than return string
NxsString BASICCMDLINE::ReturnDataForR(bool allchar, bool polymorphictomissing, bool levelsall) {
	//allchar: return even eliminated characters if true
	//polymorphictomissing: convert polymorphic observations to missing if true
	//levelsall: if true, takes as the number of states the max across all chars
	//Create one factor per char.
	//Factor can have info about whether or not a char is ordered, state labels, missing states
	//(i.e., if a particular char has just states 0 and 2, but the true set of possible states is 0,1,2,3,
	//you can specify this.
	//Create a vector of character names if present
	//Create a vector of taxon names
	//Make a data frame of all this
	//Give the data frame to R
	int nchartoreturn=0;
	int ntax = taxa->GetNumTaxonLabels();
	NxsString outputforR = "";
	if (!characters->IsEmpty())
	{
		outputforR+=characters->GetDatatypeName();
		outputforR+=" <- data.frame(taxa=c(";

		for (int taxon=0;taxon<ntax;taxon++) {
			outputforR+='"';
			outputforR+=RemoveUnderscoresAndSpaces(characters->GetTaxonLabel(taxon));
			outputforR+='"';
			if (taxon+1<ntax) {
				outputforR+=',';
			}
		}
		outputforR+=')';
		if (allchar) {
			nchartoreturn=characters->GetNCharTotal();
		}
		else {
			nchartoreturn=characters->GetNumIncludedChars();
		}
		for (int character=0; character<nchartoreturn; character++) { //We only pass the non-eliminated chars
			NxsString charlabel=characters->GetCharLabel(character);
			outputforR+=", ";
			if (charlabel.length()>1) {
				outputforR+="'";
				outputforR+=charlabel;
				outputforR+="'";
			}
			else {
				outputforR+="char";
				outputforR+=character+1;
			}
			outputforR+=" = ";

			outputforR+="factor(c(";
			for (int taxon=0;taxon<ntax;taxon++) {
				int statenumber=characters->GetInternalRepresentation(taxon,character,0);

				if(characters->IsMissingState(taxon,character)) {
					outputforR+="<NA>";
				}
				else if (characters->GetNumStates(taxon,character)>1) {
					if(polymorphictomissing) {
						outputforR+="<NA>";
					}
					else {
						outputforR+='{';
						for (unsigned int k=0;k<characters->GetNumStates(taxon,character);k++) {
							outputforR+=characters->GetInternalRepresentation(taxon,character,0);
							if (k+1<characters->GetNumStates(taxon,character)) {
								outputforR+=',';
							}
						}
						outputforR+='}';
					}
				}
				else {
					outputforR+=statenumber;
				}
				if (taxon+1<ntax) {
					outputforR+=',';
				}
			}
			outputforR+=')';
			if (levelsall) {
				outputforR+=", levels=c(";
				for (unsigned int l=0;l<characters->GetMaxObsNumStates(); l++) {
					outputforR+=l;
					if (l+1<characters->GetMaxObsNumStates()) {
						outputforR+=',';
					}
				}
				outputforR+=')';
			}
			else {

				NxsString levels=", levels=c(";
				NxsString labels=", labels=c(";

				for (unsigned int l=0;l<characters->GetObsNumStates(character); l++) {
					labels+= characters->GetStateLabel(character,l);
					levels+= l;
					if (l+1<characters->GetObsNumStates(character)) {
						labels+=',';
						levels+=',';
					}
				}
				levels+=')';
				labels+=')';
				if (labels.length()>levels.length()) {
					outputforR+=levels;
					outputforR+=labels;
				}
			}
			outputforR+=")\n";
		}
		outputforR+=")\n";
	}
	if (!taxa->IsEmpty())
	{
	}

	if (!trees->IsEmpty())
	{

		outputforR+= "\nBEGIN TREES;\n";
		for (unsigned k = 0; k < trees->GetNumTrees(); k++)
		{
			NxsString s = trees->GetTreeName(k);
			s.BlanksToUnderscores();
			outputforR+="\tTREE ";
			outputforR+=s;
			outputforR+=" = ";
			if (trees->IsRootedTree(k))
				outputforR+="[&R]";
			else
				outputforR+="[&U]";
			outputforR+=trees->GetTranslatedTreeDescription(k);
			outputforR+=";\n";
		}
		outputforR+="END;\n";

	}

	if (!assumptions->IsEmpty())
	{
	}

	if (!distances->IsEmpty())
	{ //fill cols first, first col has taxon 1, first row has taxon 2 (no diags)
		outputforR+="\ndistances <- structure(c(";
		vector<double> distancevector;
		for (int col=0;col<ntax-1;col++) {
			for (int row=col+1;row<ntax;row++) {
				distancevector.push_back(distances->GetDistance(row,col));
			}
		}
		for (unsigned int i=0;i<distancevector.size();i++) {
			outputforR+=distancevector.at(i);
			if (i+1<distancevector.size()) {
				outputforR+=',';
			}
		}
		outputforR+="), Size = ";
		outputforR+=ntax;
		outputforR+="L, Labels = c(";
		for (int taxon=0; taxon<ntax;taxon++) {
			outputforR+='"';
			outputforR+=RemoveUnderscoresAndSpaces(taxa->GetTaxonLabel(taxon));
			outputforR+='"';
			if (taxon+1<ntax) {
				outputforR+=", ";
			}
		}
		outputforR+="), Upper = FALSE, Diag = FALSE, class = \"dist\")\n";
	}

	if (!data->IsEmpty())
	{
	}

	//message=outputforR;
	//PrintMessage();
	return outputforR;
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called when the END or ENDBLOCK command needs to be parsed from within the BASICCMDLINE block. Basically just
 |	checks to make sure the next token in the data file is a semicolon.
 */
void BASICCMDLINE::HandleEndblock(
								  NxsToken &token)	/*  the token used to read from `in' */
{
	// Get the semicolon following END or ENDBLOCK token
	//
	token.GetNextToken();

	if (!token.Equals(";"))
	{
		errormsg = "Expecting ';' to terminate the END or ENDBLOCK command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Handles everything after the EXECUTE keyword and the terminating semicolon. Purges all blocks before executing
 |	file specified, and no warning is given of this.
 */
void BASICCMDLINE::HandleExecute(
								 NxsToken &token)	/* the token used to read from `in' */
{
	// Issuing the EXECUTE command from within a file is a no-no (at least in this program)
	//
	if (inf_open)
	{
		errormsg = "Cannot issue execute command from within a BASICCMDLINE block";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

	// Get the file name to execute (note: if filename contains underscores, these will be
	// automatically converted to spaces; user should surround such filenames with single quotes)
	//
	token.GetNextToken();

	NxsString fn = token.GetToken();

	// Get the semicolon terminating the EXECUTE command
	//
	token.GetNextToken();

	if (!token.Equals(";"))
	{
		errormsg = "Expecting ';' to terminate the EXECUTE command, but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}
	ReadFilePath(fn.c_str());

}

void BASICCMDLINE::ReadFilePath(const char * fn)
{
	if (FileExists(fn))
	{
		cerr << endl;
		cerr << "Opening " << fn << "..." << endl;

		PurgeBlocks();

		ifstream inf(fn, ios::binary | ios::in);

		inf_open = true;

		MyNexusToken ftoken(inf);

		try
		{
			Execute(ftoken);
		}
		catch(NxsException x)
		{
			NexusError(errormsg, x.pos, x.line, x.col);
		}

		if (inf_open)
			inf.close();
		inf_open = false;

		// Users are allowed to put DATA blocks in their NEXUS files, but internally the data is always
		// stored in a NxsCharacterBlock object.
		//
		if (characters->IsEmpty() && !data->IsEmpty())
		{
			data->TransferTo(*characters);
		}

	}
	else
	{
		cerr << endl;
		cerr << "Oops! Could not find specified file: " << fn << endl;
	}
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called when the HELP command needs to be parsed from within the BASICCMDLINE block.
 */
void BASICCMDLINE::HandleHelp(
							  NxsToken &token)	/* the token used to read from `in' */
{
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
	for (;;)
	{
		token.GetNextToken();

		if (token.Equals(";"))
		{
			break;
		}
		else
		{
			errormsg = "Unexpected keyword (";
			errormsg += token.GetToken();
			errormsg += ") encountered reading HELP command";
			throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}
	}

	message = "\nExamples of use of available commands:";
	message += "\n  help                     -> shows this message";
	message += "\n  log file=mylog.txt start -> opens log file named mylog.txt";
	message += "\n  log stop                 -> closes current log file";
	message += "\n  exe mydata.nex           -> executes nexus file mydata.nex";
	message += "\n  show                     -> reports on blocks currently stored";
	message += "\n  quit                     -> terminates application";
	message += "\n";
	PrintMessage();
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called when the HELP command needs to be parsed from within the BASICCMDLINE block.
 */
void BASICCMDLINE::HandleShow(
							  NxsToken &token)	/* the token used to read from `in' */
{
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
	for (;;)
	{
		token.GetNextToken();

		if (token.Equals(";"))
			break;
		else
		{
			errormsg = "Unexpected keyword (";
			errormsg += token.GetToken();
			errormsg += ") encountered reading HELP command";
			throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}
	}

	message = "\nNexus blocks currently stored:";
	PrintMessage();

	if (!taxa->IsEmpty())
	{
		cerr << "\n  TAXA block found" << endl;
		taxa->Report(cerr);
		if (logf_open)
			taxa->Report(logf);
	}

	if (!trees->IsEmpty())
	{
		cerr << "\n  TREES block found" << endl;
		trees->Report(cerr);
		if (logf_open)
			trees->Report(logf);
	}

	if (!assumptions->IsEmpty())
	{
		cerr << "\n  ASSUMPTIONS block found" << endl;
		assumptions->Report(cerr);
		if (logf_open)
			assumptions->Report(logf);
	}

	if (!distances->IsEmpty())
	{
		cerr << "\n  DISTANCES block found" << endl;
		distances->Report(cerr);
		if (logf_open)
			distances->Report(logf);
	}

	if (!characters->IsEmpty())
	{
		cerr << "\n  CHARACTERS block found" << endl;
		characters->Report(cerr);
		if (logf_open)
			characters->Report(logf);
	}

	if (!data->IsEmpty())
	{
		cerr << "\n  DATA block found" << endl;
		data->Report(cerr);
		if (logf_open)
			data->Report(logf);
	}
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called when the LOG command needs to be parsed from within the BASICCMDLINE block.
 */
void BASICCMDLINE::HandleLog(
							 NxsToken &token)	/* the token used to read from `in' */
{
	bool starting = false;
	bool stopping = false;
	bool appending = false;
	bool replacing = false;
	bool name_provided = false;
	NxsString logfname;

	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
	for (;;)
	{
		token.GetNextToken();

		if (token.Equals(";"))
		{
			break;
		}
		else if (token.Abbreviation("STOp"))
		{
			stopping = true;
		}
		else if (token.Abbreviation("STArt"))
		{
			starting = true;
		}
		else if (token.Abbreviation("Replace"))
		{
			replacing = true;
		}
		else if (token.Abbreviation("Append"))
		{
			appending = true;
		}
		else if (token.Abbreviation("File"))
		{
			logfname = GetFileName(token);
			name_provided = true;
		}
		else
		{
			errormsg = "Unexpected keyword (";
			errormsg += token.GetToken();
			errormsg += ") encountered reading LOG command";
			throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}
	}

	// Check for incompatible combinations of keywords
	//
	if (stopping && (starting || appending || replacing || name_provided))
	{
		errormsg = "Cannot specify STOP with any of the following START, APPEND, REPLACE, FILE";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

	if (appending && replacing)
	{
		errormsg = "Cannot specify APPEND and REPLACE at the same time";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

	if (logf_open && (starting || name_provided || appending || replacing))
	{
		errormsg = "Cannot start log file since log file is already open";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

	// Is user closing an open log file?
	//
	if (stopping)
	{
		logf.close();
		logf_open = false;

		message = "\nLog file closed";
		PrintMessage();

		return;
	}

	// If this far, must be attempting to open a log file
	//
	if (!name_provided)
	{
		errormsg = "Must provide a file name when opening a log file\n";
		errormsg += "e.g., log file=doofus.txt start replace;";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

	if (appending)
	{
		logf_open = true;
		logf.open(logfname.c_str(), ios::out | ios::app);

		message = "\nAppending to log file ";
		message += logfname;
		PrintMessage();
	}

	else if (replacing)
	{
		logf_open = true;
		logf.open(logfname.c_str());

		message = "\nReplacing log file ";
		message += logfname;
		PrintMessage();
	}

	else
	{
		bool exists = FileExists(logfname.c_str());
		bool userok = true;
		if (exists && !UserQuery("Ok to replace?", "Log file specified already exists", BASICCMDLINE::UserQueryEnum(BASICCMDLINE::uq_yes | BASICCMDLINE::uq_no)))
			userok = false;

		if (userok)
		{
			logf_open = true;
			logf.open(logfname.c_str());
		}

		if (exists && userok)
		{
			message = "\nReplacing log file ";
			message += logfname;
		}

		else if (userok)
		{
			message = "\nLog file ";
			message += logfname;
			message += " opened";
		}

		else
		{
			message = "\nLog command aborted";
		}

		PrintMessage();
	}
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Accepts a string in the form of a BASICCMDLINE block containing one command and processes it just like a real
 |	BASICCMDLINE block in a NEXUS data file.
 */
void BASICCMDLINE::HandleNextCommand()
{
	std::istringstream cmdin(next_command);

	MyNexusToken token(cmdin);
	try
	{
		Read(token);
	}
	catch(NxsException x)
	{
		NexusError(errormsg, x.pos, x.line, x.col);
	}
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called when an error is encountered in a NEXUS file. Allows program to give user details of the error as well as
 |	the precise location of the error.
 */
void BASICCMDLINE::NexusError(
							  NxsString msg,	/* the error message */
							  file_pos ,		/* the point in the NEXUS file where the error occurred */
							  long line,		/* the line in the NEXUS file where the error occurred */
							  long col)			/* the column in the NEXUS file where the error occurred */
{
	message = "\n";
	message += msg;
	PrintMessage();

	if (inf_open)
	{
		message = "Line:   ";
		message += line;
		PrintMessage();

		message = "Column: ";
		message += col;
		PrintMessage();
	}
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Begins with the command just entered by the user, which is stored in the data member `next_command', adds a
 |	semicolon (if the user failed to supply one), and then adds the string "end;" so the whole bundle looks like a
 |	very short BASICCMDLINE block. This is then passed to HandleNextCommand, which processes it just like a real
 |	BASICCMDLINE block in a NEXUS data file.
 */
void BASICCMDLINE::PreprocessNextCommand()
{
	// If user failed to add the terminating semicolon, we'll do it now. We will also remove the line feed
	// at the end and add the command "end;" to the end of the line (see explanation below).
	//
	unsigned len = (unsigned)strlen(next_command);
	assert(len > 0);

	// Remove any whitespace characters from end of string entered by user
	//
	unsigned i = len;
	while (i > 0 && (next_command[i-1] == ' ' || next_command[i-1] == '\t' || next_command[i-1] == '\n'))
		i--;

	// If character at position i - 1 is a semicolon, put '\0' terminator at position i;
	// otherwise, put a semicolon at position i and terminator at i + 1
	//
	if (next_command[i-1] != ';')
	{
		next_command[i] = ';';
		i++;
	}
	assert(i <= COMMAND_MAXLEN);
	next_command[i] = '\0';

	// Now add a semicolon at the beginning and terminate with an "END;" command
	// so that we can pretend this is simply a very short private NEXUS block
	// containing only one command.  This allows us to simply use the Read
	// function we inherited from the base class BstBase to process the command.
	//
	len = (unsigned)strlen(next_command);
	assert(len < COMMAND_MAXLEN-2);
	NxsString tmp = ";";
	tmp += next_command;
	tmp += "end;";
	strcpy(next_command, tmp.c_str());
}

/*----------------------------------------------------------------------------------------------------------------------
 |	All output is funneled through here. Writes string currently stored in `message' (a NxsString data member) to the
 |	output file stream, if open, and also to the console via cerr. Places a newline after the string if `linefeed' is
 |	true.
 */
void BASICCMDLINE::PrintMessage(
								bool linefeed)	/* if true, places newline character after message */
{
	cerr << message;
	if (linefeed)
		cerr << endl;

	if (logf_open)
	{
		logf << message;
		if (linefeed)
			logf << endl;
	}
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Detaches all blocks, deletes them, creates new blocks, and finally adds the new blocks. Call this function if
 |	you want to be sure that there is no data currently stored in any blocks.
 */
void BASICCMDLINE::PurgeBlocks()
{
	if (blockList != NULL)
	{
		Detach(taxa);
		Detach(trees);
		Detach(assumptions);
		Detach(distances);
		Detach(characters);
		Detach(data);
	}

	delete taxa;
	delete trees;
	delete assumptions;
	delete distances;
	delete characters;
	delete data;

	taxa		= new NxsTaxaBlock();
	trees		= new NxsTreesBlock(taxa);
	assumptions	= new NxsAssumptionsBlock(taxa);
	distances	= new NxsDistancesBlock(taxa);
	characters	= new NxsCharactersBlock(taxa, assumptions);
	data		= new NxsDataBlock(taxa, assumptions);

	Add(taxa);
	Add(trees);
	Add(assumptions);
	Add(distances);
	Add(characters);
	Add(data);
}

/*----------------------------------------------------------------------------------------------------------------------
 |	This function provides the ability to read everything following the block name (which is read by the NxsReader
 |	object) to the END or ENDBLOCK statement. Characters are read from the input stream `in'. Overrides the virtual
 |	function in the base class.
 */
void BASICCMDLINE::Read(
						NxsToken &token)	/* the token used to read from `in' */
{
	isEmpty = false;

	// This should be the semicolon after the block name
	//
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

	for (;;)
	{
		token.GetNextToken();

		if (token.Abbreviation("ENdblock"))
		{
			HandleEndblock(token);
			break;
		}
		else if (token.Abbreviation("Help"))
		{
			HandleHelp(token);
		}
		else if (token.Abbreviation("Log"))
		{
			HandleLog(token);
		}
		else if (token.Abbreviation("EXecute"))
		{
			HandleExecute(token);
		}
		else if (token.Abbreviation("Show"))
		{
			HandleShow(token);
		}
		else if (token.Abbreviation("REturn"))
		{
			HandleReturnData(token);
			//NxsString null=ReturnDataForR(false, true, false);
		}
		else if (token.Abbreviation("Quit"))
		{
			quit_now = true;

			message = "\nBASICCMDLINE says goodbye\n";
			PrintMessage();

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
		}
	}
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Overrides the NxsBlock virtual function. This function does nothing because the BASICCMDLINE block is simply a
 |	private command block and does not store any data.
 */
void BASICCMDLINE::Reset()
{
	NxsBlock::Reset();
}

/*----------------------------------------------------------------------------------------------------------------------
 |	This function outputs a brief report of the contents of this BASICCMDLINE block. Overrides the virtual function
 |	in the NxsBlock base class.
 */
void BASICCMDLINE::Report(
						  ostream &out)	/* the output stream to which to write the report */
{
	message.clear();
	PrintMessage();
	out << message << '\n';
	message = id;
	message += " block contains...";
	PrintMessage();
	out << message << '\n';
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Runs the command line interpreter, allowing BASICCMDLINE to interact with user. Typically, this is the only
 |	function called in main after a BASICCMDLINE object is created. If `infile_name' is non-NULL, the first command
 |	executed by the command interpreter will be "EXECUTE `infile_name'".
 */
void BASICCMDLINE::Run(
					   char *infile_name)	/* the name of the NEXUS data file to execute (can be NULL) */
{

	if (infile_name != NULL)
	{
		strcpy(next_command, "exe ");
		strncat(next_command, infile_name, 252);
		PreprocessNextCommand();
		HandleNextCommand();
	}

	quit_now = false;
	while (!quit_now)
	{
		cerr << endl;
		cerr << "BASICCMDLINE> ";
		//cin.getline(next_command, COMMAND_MAXLEN);
		unsigned i = 0;
		for(;;)
		{
			int ch = cin.get();
			if (i > COMMAND_MAXLEN)
			{
				cerr << endl;
				cerr << "Error: the length of any one command cannot exceed ";
				cerr << COMMAND_MAXLEN << " characters" << endl;
				break;
			}
			else if (ch == 10 || ch == 13)
				break;
			next_command[i++] = (char)ch;
			next_command[i] = '\0';
		}
		PreprocessNextCommand();
		HandleNextCommand();
	}
}

void BASICCMDLINE::Initialize(
					   char *infile_name,	/* the name of the NEXUS data file to execute (can be NULL) */
					   int verboseLevel
{
	ReadFilePath(infile_name);
}


/*----------------------------------------------------------------------------------------------------------------------
 |	Called when program does not recognize a block name encountered in a NEXUS file. Virtual function that overrides
 |	the virtual function in the base class NxsReader.
 */
void BASICCMDLINE::SkippingBlock(
								 NxsString blockName)	/* the unrecognized block name */
{
	message = "Skipping unknown block (";
	message += blockName;
	message += ")";
	PrintMessage();
}

/*----------------------------------------------------------------------------------------------------------------------
 |	This function is called when an unknown command named `commandName' is about to be skipped. This version of the
 |	function (which is identical to the base class version) does nothing (i.e., no warning is issued that a command
 |	was unrecognized). Modify this virtual function to provide such warnings to the user (or eliminate it altogether
 |	since the base class version already does what this does).
 */
void BASICCMDLINE::SkippingCommand(
								   NxsString commandName)	/* the name of the command being skipped */
{
	message = "Skipping unknown command (";
	message += commandName;
	message += ")";
	PrintMessage();
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Called by the NxsReader object when skipping a block named blockName that has been disabled. Allows program to
 |	notify user of progress in parsing the NEXUS file. Virtual function that overrides the virtual function in the
 |	base class NxsReader.
 */
void BASICCMDLINE::SkippingDisabledBlock(
										 NxsString )	/* the name of the block just exited */
{
}

/*----------------------------------------------------------------------------------------------------------------------
 |	The code here is identical to the base class version (simply returns 0), so the code here should either be modified
 |	or this derived version eliminated altogether. Under what circumstances would you need to modify the default code,
 |	you ask? This function should be modified to something meaningful if this derived class needs to construct and run
 |	a NxsSetReader object to read a set involving taxa. The NxsSetReader object may need to use this function to look
 |	up a taxon label encountered in the set. A class that overrides this method should return the taxon index in the
 |	range [1..ntax]; i.e., add one to the 0-offset index.
 */
unsigned BASICCMDLINE::TaxonLabelToNumber(
										  NxsString )	/* the taxon label to be translated to a taxon number */
{
	return 0;
}

/*----------------------------------------------------------------------------------------------------------------------
 |	Returns true if response is either "ok" or "yes", and returns false if response is either "no" or "cancel".
 |	This is a general query function that can handle many situations. The possible responses are enumerated in
 |	BASICCMDLINE::UserQueryEnum: uq_cancel, uq_ok, uq_yes, and uq_no. Not yet fully implemented: only handles uq_ok
 |	alone or the (uq_yes | uq_no) combination.
 */
bool BASICCMDLINE::UserQuery(
							 NxsString mb_message,						/* the question posed to the user */
							 NxsString mb_title,						/* the title of the message box */
							 BASICCMDLINE::UserQueryEnum mb_choices)	/* bit combination of uq_xx values indicating which buttons to show */
{
	const bool yes_no			= (mb_choices == (BASICCMDLINE::uq_yes | BASICCMDLINE::uq_no));
	const bool ok_only		= (mb_choices == BASICCMDLINE::uq_ok);
	assert(ok_only || yes_no); // Still working on other choices

	if (ok_only)
	{
		cerr << endl;
		cerr << mb_title << endl;
		cerr << "  " << mb_message;
		cerr << " (press return to acknowledge) ";
		cin.getline(next_command, COMMAND_MAXLEN);
		return true;
	}
	cerr << endl;
	cerr << mb_title << endl;
	cerr << "  " << mb_message;
	cerr << " (y/n) ";

	cin.getline(next_command, COMMAND_MAXLEN);

	// This could be made much simpler by just checking first letter: if 'y' then
	// assume yes, if 'n' assume no.
	//
	bool yep  = (next_command[0] == 'y' && next_command[1] == '\0');
	bool nope = (next_command[0] == 'n' && next_command[1] == '\0');

	while (!yep && !nope)
	{
		cerr << endl;
		cerr << "Must answer by typing either y or n and then pressing the Enter key" << endl;
		cerr << endl;
		cerr << mb_title << endl;
		cerr << "  " << mb_message;
		cerr << " (y/n) ";

		cin.getline(next_command, COMMAND_MAXLEN);
		yep  = (next_command[0] == 'y' && next_command[1] == '\0');
		nope = (next_command[0] == 'n' && next_command[1] == '\0');
	}

	return yep;
}


#endif // #if defined(NEW_NCL_INTERFACE)

NxsString BASICCMDLINE::RemoveUnderscoresAndSpaces(NxsString input) const
{
	while (input.find( "_", 0 ) != string::npos ) {
		input.erase((input.find( "_", 0 )),1);
	}
	while (input.find( " ", 0 ) != string::npos ) {
		input.erase((input.find( " ", 0 )),1);
	}
	input+="";
	return input;
}


int main(int argc, char *argv[])
{
	BASICCMDLINE basiccmdline;
	for (int i = 1; i < argc; ++i) {
		basiccmdline.Initialize(argv[i]);
	}
	NxsString t;
	basiccmdline.RReturnCharacters(t, true, false, false, false);
	if (!t.empty())
		cout << t << '\n';
	t.clear();
	basiccmdline.RReturnTrees(t);
	if (!t.empty())
		cout << t << '\n';
	t.clear();
	basiccmdline.RReturnDistances(t);
	if (!t.empty())
		cout << t << '\n';
	return 0;
}

