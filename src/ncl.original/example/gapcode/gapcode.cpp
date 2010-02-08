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

long gStrictLevel = 2;
long gInterleaveLen = -1;
	
void processContent(PublicNexusReader & nexusReader, ostream *out);


void vecToListOfIntWeights(const std::vector<int> &, NxsTransformationManager::ListOfIntWeights & iw);

void vecToListOfDblWeights(const std::vector<double> &dblWts, NxsTransformationManager::ListOfDblWeights & dw);


void vecToListOfIntWeights(const std::vector<int> & v, NxsTransformationManager::ListOfIntWeights & iw)
	{
	std::map<int, NxsUnsignedSet> m;
	unsigned index = 0;
	for (std::vector<int>::const_iterator vIt = v.begin(); vIt != v.end(); ++vIt, ++index)
		{
		std::map<int, NxsUnsignedSet>::iterator el = m.find(*vIt);
		if (el == m.end())
			{
			std::set<unsigned> s;
			s.insert(index);
			m[*vIt] = s;
			}
		else
			el->second.insert(index);
		}
	for (std::map<int, NxsUnsignedSet>::const_iterator mIt = m.begin(); mIt != m.end(); ++mIt)
		iw.push_back(NxsTransformationManager::IntWeightToIndexSet(mIt->first, mIt->second));
	}

void vecToListOfDblWeights(const std::vector<double> & v, NxsTransformationManager::ListOfDblWeights & iw)
	{
	std::map<double, NxsUnsignedSet> m;
	unsigned index = 0;
	for (std::vector<double>::const_iterator vIt = v.begin(); vIt != v.end(); ++vIt, ++index)
		{
		std::map<double, NxsUnsignedSet>::iterator el = m.find(*vIt);
		if (el == m.end())
			{
			std::set<unsigned> s;
			s.insert(index);
			m[*vIt] = s;
			}
		else
			el->second.insert(index);
		}
	for (std::map<double, NxsUnsignedSet>::const_iterator mIt = m.begin(); mIt != m.end(); ++mIt)
		iw.push_back(NxsTransformationManager::DblWeightToIndexSet(mIt->first, mIt->second));
	}


void writeCharactersAsGapped(ostream & out, NxsCharactersBlock * cb);

////////////////////////////////////////////////////////////////////////////////
// Takes NxsReader that has successfully read a file, and processes the
//	information stored in the reader. 
//
// The caller is responsibel for calling DeleteBlocksFromFactories() to clean
//	up (if the reader uses the factory API).
////////////////////////////////////////////////////////////////////////////////
void processContent(PublicNexusReader & nexusReader, ostream *out)
	{
	if (!out)
		return;
	BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();

	*out << "#NEXUS\n";
	for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt)
		{
		NxsBlock * b = *bIt;
		if (b->GetID() == "TAXA")
			b->WriteAsNexus(*out);
		else if (b->GetID() == "CHARACTERS" || b->GetID() == "DATA" ) 
			{
			NxsCharactersBlock * charactersBPtr = (NxsCharactersBlock *) b;
			writeCharactersAsGapped(*out, charactersBPtr);
			}
		}
	}

////////////////////////////////////////////////////////////////////////////////
// Creates a NxsReader, and tries to read the file `filename`.  If the
//	read succeeds, then processContent will be called.
////////////////////////////////////////////////////////////////////////////////
void processFilepath(
	const char * filename, // name of the file to be read
	ostream *out, // output stream to use (NULL for no output). Not that cerr is used to report errors.
	MultiFormatReader::DataFormatType fmt) // enum indicating the file format to expect.
	{
	assert(filename);
	try
		{
		MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
		if (gStrictLevel != 2)
			nexusReader.SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
			NxsCharactersBlock * charsB = nexusReader.GetCharactersBlockTemplate();
			NxsDataBlock * dataB = nexusReader.GetDataBlockTemplate();
			charsB->SetAllowAugmentingOfSequenceSymbols(true);
			dataB->SetAllowAugmentingOfSequenceSymbols(true);
		if (gInterleaveLen > 0)
			{
			assert(charsB);
			charsB->SetWriteInterleaveLen(gInterleaveLen);
			dataB->SetWriteInterleaveLen(gInterleaveLen);
			}
		NxsTreesBlock * treesB = nexusReader.GetTreesBlockTemplate();
		assert(treesB);
		if (gStrictLevel < 2)
			treesB->SetAllowImplicitNames(true);
		if (gStrictLevel < 2)
			{
			NxsStoreTokensBlockReader *storerB =  nexusReader.GetUnknownBlockTemplate();
			assert(storerB);
			storerB->SetTolerateEOFInBlock(true);  
			}
		cerr << "Executing" <<endl;
		try {
			nexusReader.ReadFilepath(filename, fmt);
			processContent(nexusReader, out);
			}
		catch(...) 
			{
			nexusReader.DeleteBlocksFromFactories();
			throw;
			}
		nexusReader.DeleteBlocksFromFactories();
		}
	catch (const NxsException &x)
		{
		cerr << "Error:\n " << x.msg << endl;
		if (x.line >=0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
		}
	}

void readFilepathAsNEXUS(const char *filename, MultiFormatReader::DataFormatType fmt)
	{
	cerr << "[Reading " << filename << "	 ]" << endl;
	try {
		ostream * outStream = 0L;
		outStream = &cout;
		processFilepath(filename, outStream, fmt);
		}
	catch (...) 
		{
		cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
		exit(1);
		}
	}	

void readFilesListedIsFile(const char *masterFilepath, MultiFormatReader::DataFormatType fmt)
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
			readFilepathAsNEXUS(filename, fmt);
		}
	}

void printHelp(ostream & out)
	{
	out << "NEXUSgapcode takes reads a file and rewrites the characters blocks with gaps coded as missing and an additional character block that displayse the presence/absence of a base for all gapped columns in the alignment.\n";
	out << "\nThe most common usage is simply:\n    NEXUSgapcode <path to NEXUS file>\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
	out << "    -l<path> reads a file and treats each line of the file as a path to NEXUS file\n\n";
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
		//passs
#	elif defined(JUST_REPORT_NEXUS) && JUST_REPORT_NEXUS
		//passs
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
	}

int main(int argc, char *argv[])
	{
	MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);
	
	bool readfile = false;
	for (int i = 1; i < argc; ++i)
		{
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
			printHelp(cout);
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 's')
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
#	else
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'i')
			{
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gInterleaveLen)) || gInterleaveLen < 1)
				{
				cerr << "Expecting a positive integer after -i\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}		
#	endif
		else if (slen > 1 && filepath[0] == '-' && filepath[1] == 'f')
			{
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2)
				{
				std::string fmtName(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(fmtName);
				}
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT)
				{
				cerr << "Expecting a format after -f\n" << endl;
				printHelp(cerr);
				return 2;
				}
			}
		else if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l')
			{
			readfile = true;
			readFilesListedIsFile(filepath+2, f);
			}
		else
			{
			readfile = true;
			readFilepathAsNEXUS(filepath, f);
			}
		}
	if (!readfile)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
		}
	return 0;
	}


void writeCoreAssumptions(ostream &out, NxsCharactersBlock * cb, const char * newTitle)
{
	std::vector<int>			intWts;
	std::vector<double>			dblWts;
	NxsUnsignedSet activeExSet = cb->GetExcludedIndexSet();
	
	const NxsTransformationManager &tm = cb->GetNxsTransformationManagerRef();
	intWts = tm.GetDefaultIntWeights();
	if (intWts.empty())
		dblWts = tm.GetDefaultDoubleWeights();

	if (!(activeExSet.empty() && intWts.empty() && dblWts.empty()))
		{
		out << "BEGIN ASSUMPTIONS; \n    LINK CHARACTERS = ";
		out << NxsString::GetEscaped(newTitle) << " ;\n";
		if (!activeExSet.empty())
			{
			NxsString exsetName;
			exsetName << newTitle;
			exsetName.append("ExSet");
			NxsUnsignedSetMap m;
			m[exsetName] = activeExSet;
			NxsWriteSetCommand("EXSET", m, out, exsetName.c_str());;
			}
		if (!(intWts.empty() && dblWts.empty()))
			{
			NxsTransformationManager &cbntm = cb->GetNxsTransformationManagerRef();
			const std::string &wtSetName =  cbntm.GetDefaultWeightSetName(); 
			NxsTransformationManager ntm;
			if (!intWts.empty())
				{
				NxsTransformationManager::ListOfIntWeights iw;
				vecToListOfIntWeights(intWts, iw);
				ntm.AddIntWeightSet(wtSetName, iw, true);
				}
			else
				{
				NxsTransformationManager::ListOfDblWeights dw;
				vecToListOfDblWeights(dblWts, dw);
				ntm.AddRealWeightSet(wtSetName, dw, true);
				}
			ntm.WriteWtSet(out);
			}
		out << "END;\n";
		}
}

void writeCharactersAsGapped(ostream &out, NxsCharactersBlock * cb)
{
	if (!cb )
		return;
	NxsTaxaBlockAPI *tb = cb->GetTaxaBlockPtr(NULL);
	if (tb == NULL)
		return;
	NxsCharactersBlock::DataTypesEnum dt = cb->GetDataType();
	const char g = cb->GetGapSymbol();
	const std::string baseTitle = cb->GetTitle();
	if (tb == NULL
		|| dt == NxsCharactersBlock::standard 
		|| dt == NxsCharactersBlock::continuous 
		|| dt == NxsCharactersBlock::mixed
		|| g == '\0')
		{
		cb->WriteAsNexus(out);
		writeCoreAssumptions(out, cb, baseTitle.c_str());

		return;
		}
	bool isAutogen = cb->IsAutoGeneratedTitle();
	std::string newTitle = baseTitle;
	newTitle.append("GapsAsMissing");
	cb->SetTitle(newTitle, isAutogen);
	
	std::set<unsigned> gappedColumns;

	out << "BEGIN CHARACTERS;\n";
	cb->WriteBasicBlockCommands(out);

	const unsigned ntaxTotal = tb->GetNTax();
	out << "    DIMENSIONS";
	if (tb)
		{
		const unsigned wod = cb->GetNTaxWithData();
		if (wod > 0 && wod != ntaxTotal)
			out << " NTax=" << wod;
		}
	const unsigned nc = cb->GetNCharTotal();
	out << " NChar=" << nc << ";\n";
	cb->WriteEliminateCommand(out);
	cb->SetGapSymbol('\0');
	cb->WriteFormatCommand(out);
	cb->SetGapSymbol(g);

	cb->WriteCharStateLabelsCommand(out);

	std::vector<int>			intWts;
	std::vector<double>			dblWts;
	NxsUnsignedSet activeExSet = cb->GetExcludedIndexSet();
	
	const NxsTransformationManager &tm = cb->GetNxsTransformationManagerRef();
	intWts = tm.GetDefaultIntWeights();
	if (intWts.empty())
		dblWts = tm.GetDefaultDoubleWeights();

	unsigned width = tb->GetMaxTaxonLabelLength();
	out << "Matrix\n";
	unsigned begCharInd = 0;
	unsigned endCharInd = nc;
	for (unsigned i = 0; i < ntaxTotal; i++)
		{
		if (cb->TaxonIndHasData(i))
			{
			const std::string currTaxonLabel = NxsString::GetEscaped(tb->GetTaxonLabel(i));
			out << currTaxonLabel;
			unsigned currTaxonLabelLen = (unsigned)currTaxonLabel.size();
			unsigned diff = width - currTaxonLabelLen;
			for (unsigned k = 0; k < diff + 5; k++)
				out << ' ';

			const NxsDiscreteStateRow & row = cb->GetDiscreteMatrixRow(i);
			if (!row.empty())
				{
				const NxsDiscreteDatatypeMapper * dm = cb->GetDatatypeMapperForChar(0);
				if (dm == NULL)
					throw NxsNCLAPIException("No DatatypeMapper in WriteStatesForTaxonAsNexus");
				if (cb->IsMixedType())
					throw NxsNCLAPIException("Mixed datatypes are not supported by the NEXUSgapcode program.");
				if (cb->IsTokens())
					{
					for (unsigned charInd = begCharInd; charInd < endCharInd; ++charInd)
						{
						int sc = row[charInd];
						out << ' ';
						if (sc == NXS_GAP_STATE_CODE)
							out << dm->GetMissingSymbol();
						else 
							{
							NxsString sl = cb->GetStateLabel(charInd, sc);
							if (sl == " ")
								{
								NxsString errormsg = "Writing character state ";
								errormsg << 1 + sc << " for character " << 1+charInd << ", but no appropriate chararcter label or symbol was found.";
								throw NxsNCLAPIException(errormsg);
								}
							else
								out  << NxsString::GetEscaped(sl);
							}
						}
					}
				else
					{
					std::vector<int>::const_iterator endIt = row.begin() + begCharInd;
					std::vector<int>::const_iterator begIt = endIt;
					if (endCharInd == row.size())
						endIt = row.end();
					else
						endIt += endCharInd - begCharInd;
					unsigned j = begCharInd;
					for (; begIt != endIt; ++begIt, ++j)
						{
						int c;
						if (*begIt == NXS_GAP_STATE_CODE)
							{
							c = NXS_MISSING_CODE;
							gappedColumns.insert(j);
							}
						else
							c = *begIt;
						dm->WriteStateCodeAsNexusString(out, c, true);
						}
					}
				}
			out << '\n';
			}
		}
	out << ";\n";
	cb->WriteSkippedCommands(out);
	out << "END;\n";
	
	
	writeCoreAssumptions(out, cb, newTitle.c_str());
	const unsigned nGappedCols = gappedColumns.size();
	if (nGappedCols > 0)
		{
		newTitle = baseTitle;
		newTitle.append("GapsAsBinary");
		cb->SetTitle(newTitle, isAutogen);
			out << "BEGIN CHARACTERS;\n";
		cb->WriteBasicBlockCommands(out);
	
		out << "    DIMENSIONS";
		if (tb)
			{
			const unsigned wod = cb->GetNTaxWithData();
			if (wod > 0 && wod != ntaxTotal)
				out << " NTax=" << wod;
			}
		out << " NChar=" << nGappedCols << ";\n";
		out << " CharStateLabels " ;
		unsigned currChNumber = 1;
		std::set<unsigned>::iterator gcIt = gappedColumns.begin(); 
		out << currChNumber++ << " col_" << (1 + *gcIt);
		for (++gcIt ; gcIt != gappedColumns.end(); ++gcIt)
			out << ",\n    " << currChNumber++ << " col_" << (1 + *gcIt);
		out << " ;\n" ;
		out << "Format Datatype = Standard Symbols=\"01\" missing = '?' ;\n";
	
		out << "Matrix\n";
		for (unsigned i = 0; i < ntaxTotal; i++)
			{
			if (cb->TaxonIndHasData(i))
				{
				const std::string currTaxonLabel = NxsString::GetEscaped(tb->GetTaxonLabel(i));
				out << currTaxonLabel;
				unsigned currTaxonLabelLen = (unsigned)currTaxonLabel.size();
				unsigned diff = width - currTaxonLabelLen;
				for (unsigned k = 0; k < diff + 5; k++)
					out << ' ';
	
				const NxsDiscreteStateRow & row = cb->GetDiscreteMatrixRow(i);
				for (std::set<unsigned>::iterator cIt = gappedColumns.begin(); cIt != gappedColumns.end(); ++cIt)
					{
					int sc = row[*cIt];
					if (sc == NXS_GAP_STATE_CODE)
						out << '0';
					else 
						out << '1';
					}
				out << '\n';
				}
			}
		out << ";\n";
		cb->WriteSkippedCommands(out);
		out << "END;\n";
		}

	std::vector<int>			gapintWts;
	std::vector<double>			gapdblWts;
	std::set<unsigned>			gapactiveExSet;
	if (!activeExSet.empty())
		{
		unsigned gapind = 0;
		for (std::set<unsigned>::iterator cIt = gappedColumns.begin(); cIt != gappedColumns.end(); ++cIt, ++gapind)
			{
			if (activeExSet.find(*cIt) != activeExSet.end())
				gapactiveExSet.insert(gapind);
			}
		}
	if (!intWts.empty())
		{
		for (std::set<unsigned>::iterator cIt = gappedColumns.begin(); cIt != gappedColumns.end(); ++cIt)
			gapintWts.push_back(intWts[*cIt]);
		}
	if (!dblWts.empty())
		{
		for (std::set<unsigned>::iterator cIt = gappedColumns.begin(); cIt != gappedColumns.end(); ++cIt)
			gapdblWts.push_back(dblWts[*cIt]);
		}

	if (!(gapactiveExSet.empty() && gapintWts.empty() && gapdblWts.empty()))
		{
		out << "BEGIN ASSUMPTIONS; \n    LINK CHARACTERS = ";
		out << NxsString::GetEscaped(newTitle) << " ;\n";
		if (!gapactiveExSet.empty())
			{
			NxsString exsetName;
			exsetName << newTitle;
			exsetName.append("GapExSet");
			NxsUnsignedSetMap m;
			m[exsetName] = gapactiveExSet;
			NxsWriteSetCommand("EXSET", m, out, exsetName.c_str());;
			}
		if (!(gapintWts.empty() && gapdblWts.empty()))
			{
			NxsTransformationManager &cbntm = cb->GetNxsTransformationManagerRef();
			std::string wtSetName =  cbntm.GetDefaultWeightSetName(); 
			wtSetName.append("GapWtSet");
			NxsTransformationManager ntm;
			if (!gapintWts.empty())
				{
				NxsTransformationManager::ListOfIntWeights iw;
				vecToListOfIntWeights(gapintWts, iw);
				ntm.AddIntWeightSet(wtSetName, iw, true);
				}
			else
				{
				NxsTransformationManager::ListOfDblWeights dw;
				vecToListOfDblWeights(gapdblWts, dw);
				ntm.AddRealWeightSet(wtSetName, dw, true);
				}
			ntm.WriteWtSet(out);
			}
		out << "END;\n";
		}


	cb->SetTitle(baseTitle, isAutogen);
}
