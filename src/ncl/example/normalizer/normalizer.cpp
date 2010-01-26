//	Copyright (C) 2007-2008 Mark T. Holder
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

/*******************************************************************************
 * This file contains code for 4 executables:
 *		NEXUSnormalizer, NEXUSvalidator, NEXUSinspector, and NEX_us2ml
 *	with conditional compilation used to determine the behavior.
 *
 *		* NEXUSnormalizer - writes a NEXUS version of the file with consistent
 *			ordering of blocks and commands. Ideally 2 equivalent files will 
 *			produce the same normalized output. This version of tthe program is
 *			less ambitious. The goal is to be able to run (for any valid NEXUS 
 *			in.nex file):
 *				$ NEXUSnormalizer in.nex > outOrig.nex
 *				$ NEXUSnormalizer outOrig.nex > outSecond.nex
 *				$ diff outOrig.nex outSecond.nex
 *			and find no differences.
 *		* NEXUSvalidator - reports errors and warnings to stderr. Invalid files
 *			cause exit with a non-zero return code
 *		* NEXUSinspector - writes a brief report of every block parsed
 *		* NEXUS_us2ml - writes a nexml version of the input (partially 
 *			implemented, note that the code to write nexml is in us2ml.cpp).
 * See the processFilepath() function for an example of how to deal with NCL
 *	to read a file using the new MultiFormatReader class. When the file
 *	is correctly read, the processContent() function is called.
 * 
 * All other code has to do with reading command line arguments and other
 * 	user-interface concerns.
 */
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"
#include <cassert>
//#include "ncl/nxscxxdiscretematrix.h"

#if defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
	void	writeAsNexml(PublicNexusReader & nexusReader, ostream & os);
#endif
#if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
	enum ExportFormatEnum
		{
		NEXUS_EXPORT_FORMAT,
		PHYLIP_EXPORT_FORMAT,
		RELAXED_PHYLIP_EXPORT_FORMAT,
		FASTA_EXPORT_FORMAT,
		NEXML_EXPORT_FORMAT,
		UNSUPPORTED_EXPORT_FORMAT
		};
	ExportFormatEnum gExportFormat = NEXML_EXPORT_FORMAT;
	std::string gExportPrefix("out");
	ExportFormatEnum readExportFormatName(const std::string &);
	void exportData(PublicNexusReader & nexusReader, MultiFormatReader::DataFormatType f, long interleavLen, std::string prefix);
	
	

ExportFormatEnum readExportFormatName(const std::string & s)
{
	const char * gExportFormatNames[] = {   "nexus",
											"phylip",
											"relaxedphylip",
											"fasta",
											"nexml"
											};
	const unsigned gNumExportFormats = 5;

	NxsString l(s.c_str());
	NxsString::to_lower(l);
	int ind = NxsString::index_in_array(l, gExportFormatNames, gNumExportFormats);
	if (ind < 0)
		return UNSUPPORTED_EXPORT_FORMAT;
	return ExportFormatEnum(ind);
}

#endif

bool gAltNexus = false;

void writeAsNexus(PublicNexusReader & nexusReader, ostream & os);

long gStrictLevel = 2;
bool gTreesViaInMemoryStruct = true;
long gInterleaveLen = -1;
bool blocksReadInValidation = false;

enum ProcessActionsEnum
	{
	REPORT_BLOCKS,
	OUTPUT_NORMALIZED_NEXUS,
	OUTPUT_ANY_FORMAT,
	OUTPUT_NEXML,
	VALIDATE_ONLY
	};
	
	
void processContent(PublicNexusReader & nexusReader, ostream *os, ProcessActionsEnum currentAction);
MultiFormatReader * instantiateReader();

#	if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
	MultiFormatReader * gNexusReader = NULL;
#	endif


void reportNexusStats(const PublicNexusReader & nexusReader, ostream *os)
{
	if (!os)
		return;

	const unsigned nTaxaBlocks = nexusReader.GetNumTaxaBlocks();
	*os <<  nTaxaBlocks << " taxa block(s) read.\n";
	for (unsigned t = 0; t < nTaxaBlocks; ++t)
		{
		NxsTaxaBlock * tb = nexusReader.GetTaxaBlock(t);
		*os << "Taxa block #" << t + 1 << ".\n";
		tb->Report(*os);
		const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(tb);
		*os <<  nCharBlocks << " Characters/Data block(s) read that link to this Taxa block.\n";
		for (unsigned i = 0; i < nCharBlocks; ++i)
			{
			NxsCharactersBlock * cb = nexusReader.GetCharactersBlock(tb, i);

			//NxsCXXDiscreteMatrix mat(*cb, true);

			*os << "Character block #" << i + 1 << " for this Taxa block.\n";
			cb->Report(*os);
			const unsigned nAssumpBlocks = nexusReader.GetNumAssumptionsBlocks(cb);
			*os <<  nAssumpBlocks << " Assumptions block(s) read that link to this Characters block.\n";
			for (unsigned j= 0; j < nAssumpBlocks; ++j)
				{
				NxsAssumptionsBlock * ab = nexusReader.GetAssumptionsBlock(cb, j);
				*os << "Assumptions block #" << j + 1 << " for this Characters block.\n";
				ab->Report(*os);
				}
			}
		const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(tb);
		*os <<  nTreesBlocks << " Trees/Data block(s) read that link to this Taxa block.\n";
		for (unsigned i = 0; i < nTreesBlocks; ++i)
			{
			NxsTreesBlock * cb = nexusReader.GetTreesBlock(tb, i);
			*os << "Trees block #" << i + 1 << " for this Taxa block.\n";
			cb->Report(*os);
			const unsigned nAssumpBlocks = nexusReader.GetNumAssumptionsBlocks(cb);
			*os <<  nAssumpBlocks << " Assumptions block(s) read that link to this Trees block.\n";
			for (unsigned j= 0; j < nAssumpBlocks; ++j)
				{
				NxsAssumptionsBlock * ab = nexusReader.GetAssumptionsBlock(cb, j);
				*os << "Assumptions block #" << j + 1 << " for this Trees block.\n";
				ab->Report(*os);
				}
			}
		const unsigned nAssumpBlocks = nexusReader.GetNumAssumptionsBlocks(tb);
		*os <<  nAssumpBlocks << " Assumptions block(s) read that link to this Taxa block.\n";
		for (unsigned j= 0; j < nAssumpBlocks; ++j)
			{
			NxsAssumptionsBlock * ab = nexusReader.GetAssumptionsBlock(tb, j);
			*os << "Assumptions block #" << j + 1 << " for this Taxa block.\n";
			ab->Report(*os);
			}
		const unsigned nDistancesBlocks = nexusReader.GetNumDistancesBlocks(tb);
		*os <<  nDistancesBlocks << " Distances block(s) read that link to this Taxa block.\n";
		for (unsigned j= 0; j < nDistancesBlocks; ++j)
			{
			NxsDistancesBlock * ab = nexusReader.GetDistancesBlock(tb, j);
			*os << "Distances block #" << j + 1 << " for this Taxa block.\n";
			ab->Report(*os);
			}
		const unsigned nUnalignedBlocks = nexusReader.GetNumUnalignedBlocks(tb);
		*os <<  nUnalignedBlocks << " Unaligned block(s) read that link to this Taxa block.\n";
		for (unsigned j= 0; j < nUnalignedBlocks; ++j)
			{
			NxsUnalignedBlock * ab = nexusReader.GetUnalignedBlock(tb, j);
			*os << "Unaligned block #" << j + 1 << " for this Taxa block.\n";
			ab->Report(*os);
			}
		*os << "\n\n";
		}
	const unsigned nUnknown = nexusReader.GetNumUnknownBlocks();
	*os <<  nUnknown << " private block(s) read.\n";
	for (unsigned t = 0; t < nUnknown; ++t)
		{
		NxsStoreTokensBlockReader * ub = nexusReader.GetUnknownBlock(t);
		*os << "Private block #" << t + 1 << " is a " << ub->GetID() << " block.\n";
		}
}


void writeAsNexus(PublicNexusReader & nexusReader, ostream & os)
{
	BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();
	os << "#NEXUS\n";
	for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt)
		{
		NxsBlock * b = *bIt;
		if (b)
			b->WriteAsNexus(os);
		}
}
////////////////////////////////////////////////////////////////////////////////
// Takes NxsReader that has successfully read a file, and processes the
//	information stored in the reader. 
//
// The caller is responsibel for calling DeleteBlocksFromFactories() to clean
//	up (if the reader uses the factory API).
////////////////////////////////////////////////////////////////////////////////
void processContent(PublicNexusReader & nexusReader, ostream *os, ProcessActionsEnum currentAction)
	{
	BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();

	if (currentAction == REPORT_BLOCKS)
		reportNexusStats(nexusReader, os);
	else if (OUTPUT_NORMALIZED_NEXUS == currentAction && os)
		{
		writeAsNexus(nexusReader, *os);
		}
	else if (OUTPUT_NEXML == currentAction && os)
		{
#		if defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
			writeAsNexml(nexusReader, *os);
#		else
			cerr << "Error nexml conversion not implemented\n";
			exit(1);
#		endif			
		}
	else if (OUTPUT_ANY_FORMAT == currentAction && os)
		{
#		if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
			std::string fullExportPrefix;
			MultiFormatReader::DataFormatType f = MultiFormatReader::NEXUS_FORMAT;
			if (gExportFormat == NEXUS_EXPORT_FORMAT)
				exportData(nexusReader, MultiFormatReader::NEXUS_FORMAT, gInterleaveLen, gExportPrefix);
			else if (gExportFormat == NEXML_EXPORT_FORMAT)
				exportData(nexusReader, MultiFormatReader::NEXML_FORMAT, gInterleaveLen, gExportPrefix);
			else if (gExportFormat == PHYLIP_EXPORT_FORMAT) {
				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".dna");
				f = (gInterleaveLen < 0 ? MultiFormatReader::PHYLIP_DNA_FORMAT : MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".rna");
				f = (gInterleaveLen < 0 ? MultiFormatReader::PHYLIP_RNA_FORMAT : MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".aa");
				f = (gInterleaveLen < 0 ? MultiFormatReader::PHYLIP_AA_FORMAT : MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix);

				//fullExportPrefix = gExportPrefix;
				//fullExportPrefix.append(".discrete");
				//f = (gInterleaveLen < 0 ? MultiFormatReader::PHYLIP_DISC_FORMAT : MultiFormatReader::INTERLEAVED_PHYLIP_DISC_FORMAT);
				//exportData(nexusReader, f, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				exportData(nexusReader, MultiFormatReader::PHYLIP_TREE_FORMAT, gInterleaveLen, fullExportPrefix);
			}
			else if (gExportFormat == RELAXED_PHYLIP_EXPORT_FORMAT) {
				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".dna");
				f = (gInterleaveLen < 0 ? MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT : MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".rna");
				f = (gInterleaveLen < 0 ? MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT : MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".aa");
				f = (gInterleaveLen < 0 ? MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT : MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT);
				exportData(nexusReader, f, gInterleaveLen, fullExportPrefix);

				//fullExportPrefix = gExportPrefix;
				//fullExportPrefix.append(".discrete");
				//f = (gInterleaveLen < 0 ? MultiFormatReader::RELAXED_PHYLIP_DISC_FORMAT : MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT);
				//exportData(nexusReader, f, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				exportData(nexusReader, MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT, gInterleaveLen, fullExportPrefix);
			}
			else if (gExportFormat == FASTA_EXPORT_FORMAT) {
				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".dna");
				exportData(nexusReader,  MultiFormatReader::FASTA_DNA_FORMAT, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".rna");
				exportData(nexusReader,  MultiFormatReader::FASTA_RNA_FORMAT, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				fullExportPrefix.append(".aa");
				exportData(nexusReader,  MultiFormatReader::FASTA_AA_FORMAT, gInterleaveLen, fullExportPrefix);

				fullExportPrefix = gExportPrefix;
				exportData(nexusReader, MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT, gInterleaveLen, fullExportPrefix);
			}
			else {
				cerr << "Unsupported export format requested.\n";
				exit(1);
			}
#		else
			cerr << "Exporting to any format is implemented by this executable.\n";
			exit(1);
#		endif			
		}
	else if (VALIDATE_ONLY == currentAction)
		{
		if (!blocks.empty())
			blocksReadInValidation = true;
		}
	}


MultiFormatReader * instantiateReader()
{
	MultiFormatReader * nexusReader = new MultiFormatReader(-1, NxsReader::WARNINGS_TO_STDERR);
	if (gStrictLevel != 2)
		nexusReader->SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
	NxsCharactersBlock * charsB = nexusReader->GetCharactersBlockTemplate();
	NxsDataBlock * dataB = nexusReader->GetDataBlockTemplate();
	charsB->SetAllowAugmentingOfSequenceSymbols(true);
	dataB->SetAllowAugmentingOfSequenceSymbols(true);
	if (gInterleaveLen > 0)
		{
		assert(charsB);
		charsB->SetWriteInterleaveLen(gInterleaveLen);
		dataB->SetWriteInterleaveLen(gInterleaveLen);
		}
	
	NxsTreesBlock * treesB = nexusReader->GetTreesBlockTemplate();
	assert(treesB);
	if (gStrictLevel < 2)
		treesB->SetAllowImplicitNames(true);
	treesB->SetWriteFromNodeEdgeDataStructure(gTreesViaInMemoryStruct);
	if (gAltNexus)
		treesB->setWriteTranslateTable(false);
	if (gStrictLevel < 2)
		{
		NxsStoreTokensBlockReader *storerB =  nexusReader->GetUnknownBlockTemplate();
		assert(storerB);
		storerB->SetTolerateEOFInBlock(true);  
		}
	nexusReader->conversionOutputRecord.addNumbersToDisambiguateNames = true;
	return nexusReader;
}

////////////////////////////////////////////////////////////////////////////////
// Creates a NxsReader, and tries to read the file `filename`.  If the
//	read succeeds, then processContent will be called.
////////////////////////////////////////////////////////////////////////////////
void processFilepath(
	const char * filename, // name of the file to be read
	ostream *os, // output stream to use (NULL for no output). Not that cerr is used to report errors.
	MultiFormatReader::DataFormatType fmt, // enum indicating the file format to expect.
	ProcessActionsEnum currentAction) // enum that is passed on to processContent to indicate what should be done with the content of the file.
	{
	assert(filename);
	try
		{
		MultiFormatReader * nexusReader;
#		if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
			nexusReader = gNexusReader;
#		else
			nexusReader = instantiateReader();
#		endif

		cerr << "Executing" <<endl;
		try {
			nexusReader->ReadFilepath(filename, fmt);
#			if !defined(MULTIFILE_NEXUS_READER) ||  !MULTIFILE_NEXUS_READER
				processContent(*nexusReader, os, currentAction);
#			endif
			}
		catch(...) 
			{
			nexusReader->DeleteBlocksFromFactories();
			throw;
			}
#		if !defined(MULTIFILE_NEXUS_READER) ||  !MULTIFILE_NEXUS_READER
			nexusReader->DeleteBlocksFromFactories();
			delete nexusReader;
#		endif
		}
	catch (const NxsException &x)
		{
		cerr << "Error:\n " << x.msg << endl;
		if (x.line > 0 || x.pos > 0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
		}
	}

void readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt, ProcessActionsEnum currentAction)
	{
	cerr << "[Reading " << filename << "	 ]" << endl;
	try {
		ostream * outStream = 0L;
		if (currentAction != VALIDATE_ONLY)
			outStream = &cout;
		processFilepath(filename, outStream, fmt, currentAction);
		}
	catch (...) 
		{
		cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
		exit(1);
		}
	}	

void readFilesListedIsFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt, ProcessActionsEnum currentAction)
	{
	ifstream masterStream(masterFilepath, ios::binary);
	if (masterStream.bad())
		{
		cerr << "Could not open " << masterFilepath << "." << endl;
		exit(3);
		}
	char filename[1024];
	while ((!masterStream.eof())  && masterStream.good())
		{
		masterStream.getline(filename, 1024);
		if (strlen(filename) > 0 && filename[0] != '#')
			readFilepathAsNEXUS(filename, fmt, currentAction);
		}
	}

void printHelp(ostream & out)
	{
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		out << "NEXUSvalidator takes reads a file and exits with a success (return code 0) if the file is valid.\n";
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		out << "NEXUSinspector takes reads a file and writes a report of the content to standard out.\n";
#	elif defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
		out << "NCLconverter takes reads a file and writes a report of the content to a file prefix (specified with the -o flag) in the chosen output format (specified with the -e flag).\n";
#	else
#		if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
			out << "NEXUSunion reads a series of NEXUS file and writes the union of all of their content to standard out (using the NEXUSnormalizer conventions of indentation and syntax).\n";
#		else
			out << "NEXUSnormalizer takes reads a file and rewrites the file to standard out with consistent indentation and syntax.\n";
#		endif
# 	endif
	out << "\nThe most common usage is simply:\n    NEXUSnormalizer <path to NEXUS file>\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
	out << "    -l<path> reads a file and treats each line of the file as a path to NEXUS file\n\n";
	out << "    -a output AltNexus (no translation table in trees)\n\n";
	out << "    -s<non-negative integer> controls the NEXUS strictness level.\n";
	out << "        the default level is equivalent to -s2 invoking the program with \n";
	out << "        -s3 or a higher number will convert some warnings into fatal errors.\n";
	out << "        Running with -s1 will cause the parser to accept dangerous constructs,\n";
	out << "        and running with -s0 will cause the parser make every attempt to finish\n";
	out << "        parsing the file (warning about very serious errors).\n\n";
	out << "        Note that when -s0 strictness level is used, and the parser fails to\n";
	out << "        finish, it will often be the result of an earlier error than the \n";
	out << "        error that is reported in the last message.\n";
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		//pass
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		//pass
#	elif defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
		//pass
#	else
		out << "    -i<number> specifies the length of the interleaved pages to create\n";
#	endif
	out << "    -f<format> specifies the input file format expected:\n";
	out << "            -fnexus     NEXUS (this is also the default)\n";
	out << "            -faafasta   Amino acid data in fasta\n";
	out << "            -fdnafasta  DNA data in fasta\n";
	out << "            -frnafasta  RNA data in fasta\n";
	out << "        The complete list of format names that can follow the -f flag is:\n";
	std::vector<std::string> fmtNames =  MultiFormatReader::getFormatNames();
	for (std::vector<std::string>::const_iterator n = fmtNames.begin(); n != fmtNames.end(); ++n)
		{
		out << "            "<< *n << "\n";
		}
#	if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
		out << "    -e<format> specifies the output file format expected:\n";
		out << "            -enexus  \"normalized\" NEXUS output\n";
		out << "            -efasta  Character data in fasta (could result in multiple output files)\n";
		out << "            -ephylip  Trees and character data in phylip (could result in multiple output files)\n";
		out << "            -erelaxedphylip  Trees and character data in relaxed phylip (could result in multiple output files)\n";
		out << "            -enexml  nexml output (this is also the default)\n";
		out << "    -o<fn> specifies the output prefix.  An appropriate suffix and extension are added\n";
#	endif
	}

int main(int argc, char *argv[])
	{
	NxsReader::setNCLCatchesSignals(true);
	MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		ProcessActionsEnum currentAction = VALIDATE_ONLY;
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		ProcessActionsEnum currentAction = REPORT_BLOCKS;
#	elif defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
		ProcessActionsEnum currentAction = OUTPUT_ANY_FORMAT;
#	elif defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
		ProcessActionsEnum currentAction = OUTPUT_NEXML;
#	else
		ProcessActionsEnum currentAction = OUTPUT_NORMALIZED_NEXUS;
# 	endif

	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 2 || filepath[0] != '-')
			continue;
		if (filepath[1] == 'h')
			printHelp(cout);
		else if (filepath[1] == 's')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gStrictLevel)))
				{
				cerr << "Expecting an integer after -s\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		//pass
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		//pass
#	elif defined(TO_NEXML_CONVERTER) && TO_NEXML_CONVERTER
		//pass
#	else
		else if (filepath[1] == 'i')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gInterleaveLen)) || gInterleaveLen < 1)
				{
				cerr << "Expecting a positive integer after -i\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}		
		else if (filepath[1] == 'a')
			{
			if ((slen != 2))
				{
				cerr << "Not expecting a value after -a\n" << endl;
				printHelp(cerr);
				return 2;
				}
			gAltNexus = true;
			}		
#	endif
#	if defined(NCL_CONVERTER_APP) && NCL_CONVERTER_APP
		else if (filepath[1] == 'e') {
			if (slen > 2)
				{
				std::string efmtName(filepath + 2, slen - 2);
				gExportFormat = readExportFormatName(efmtName);
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting a format after -e\n" << endl;
				printHelp(cerr);
				return 2;
				}
		}
		else if (filepath[1] == 'o') {
			if (slen > 2)
				{
				std::string oname(filepath + 2, slen - 2);
				gExportPrefix = oname;
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting an output file prefix after -o\n" << endl;
				printHelp(cerr);
				return 2;
				}
		}
#	endif
		else if (filepath[1] == 'f')
			{
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2)
				{
				std::string fmtName(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(fmtName);
				if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
					{
					cerr << "Unknowon format \"" << fmtName << "\" after -f\n" << endl;
					printHelp(cerr);
					return 3;
					}
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting a format after -f\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		}

#	if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
		gNexusReader = instantiateReader();
		gNexusReader->cullIdenticalTaxaBlocks(true);
#	endif
	bool readfile = false;
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 1)
			continue;
		if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l')
			{
			readfile = true;
			readFilesListedIsFile(filepath+2, f, currentAction);
			}
		else if (filepath[0] != '-')
			{
			readfile = true;
			readFilepathAsNEXUS(filepath, f, currentAction);
			}
		}
#	if defined(MULTIFILE_NEXUS_READER) && MULTIFILE_NEXUS_READER
		if (gNexusReader)
			{
			processContent(*gNexusReader, &std::cout, OUTPUT_NORMALIZED_NEXUS);
			gNexusReader->DeleteBlocksFromFactories();
			delete gNexusReader;
			}
#	endif

	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
#	if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
		if (blocksReadInValidation)
			return  0;
		std::cerr << "No blocks read\n";
		return 1;
#	else
		return 0;
#	endif
	}

