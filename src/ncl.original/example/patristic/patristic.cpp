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

/*******************************************************************************
 *	This file contains code for an executable that takes path to a NEXUS file
 *		as a command line argument and writes a "normalized" version of the 
 *		blocks.	 This is useful for testing. 
 *
 *	Ideally 2 equivalent files will produce the same normalized output. This 
 *		version of tthe program is less ambitious. The goal is to be able to run 
 *		(for any valid NEXUS in.nex file):
 *			$ normalizer in.nex > outOrig.nex
 *			$ normalizer outOrig.nex > outSecond.nex
 *			$ diff outOrig.nex outSecond.nex
 *		and find no differences.
 */
#include "ncl/ncl.h"
#include <iostream>
#include <climits>
#include <cassert>
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"


static long gStrictLevel = 1;
static long gTreeNum = -1;
enum PATRISTIC_RUN_MODE
	{
		NO_OUTPUT             = 0,
		ONLY_PATRISTIC        = 1,
		ONLY_TO_MRCA          = 2,
		PATRISTIC_AND_TO_MRCA = 3
		};
static long gMode = ONLY_PATRISTIC;

void filepathToPatristic(const char * filename, std::ostream * os);
bool printPatristic(const NxsFullTreeDescription &treeDesc, std::ostream * os, std::ostream * errStr, NxsTaxaBlockAPI &taxaB, const unsigned treeN);

bool printPatristic(const NxsFullTreeDescription &treeDesc, std::ostream * os, std::ostream * errStr, NxsTaxaBlockAPI &taxaB, const unsigned treeN)
{
	const unsigned ntax = taxaB.GetNTax();
	if (treeDesc.AllEdgesHaveLengths()) {
		for (int flag = 1; flag < 3; ++flag) {	
			if (flag&gMode) {
				const bool toMRCA = (flag == ONLY_TO_MRCA);
				if (errStr)
					{
					if (toMRCA)
						*errStr << "The to-MRCA matrix:" << std::endl;
					else
						*errStr << "The patristic distance matrix:" << std::endl;
					}
				NxsSimpleTree tree(treeDesc, 1, 1.0);
				if (!treeDesc.EdgeLengthsAreAllIntegers()) {
					NxsRealStepMatrix::DblMatrix dm = tree.GetDblPathDistances(toMRCA);
					const unsigned n = (ntax < dm.size() ? ntax : dm.size());
					if (n > 0 && os) {
						for (unsigned i = 0; i < n; ++i)
							*os << '\t' << NxsString::GetEscaped(taxaB.GetTaxonLabel(i));
						for (unsigned i = 0; i < n; ++i) {
							*os << '\n' << NxsString::GetEscaped(taxaB.GetTaxonLabel(i));
							const std::vector<double> & row = dm[i];
							for (unsigned j= 0; j < n; ++j) {
								const double & d = row[j];
								*os << '\t';
								if (i == j)
									*os << '-';
								else if (d == DBL_MAX)
									*os << 'i';
								else 
									*os << d;
							}
						}
						*os << '\n';
					}
				}
				else {
					NxsIntStepMatrix::IntMatrix im = tree.GetIntPathDistances(toMRCA);
					const unsigned n = (ntax < im.size() ? ntax : im.size());
					if (n > 0 && os) {
						for (unsigned i = 0; i < n; ++i)
							*os << '\t' << NxsString::GetEscaped(taxaB.GetTaxonLabel(i));
						for (unsigned i = 0; i < n; ++i) {
							*os << '\n' << NxsString::GetEscaped(taxaB.GetTaxonLabel(i));
							const std::vector<int> & row = im[i];
							for (unsigned j= 0; j < n; ++j) {
								const int & d = row[j];
								*os << '\t';
								if (i == j)
									*os << '-';
								else if (d == INT_MAX)
									*os << 'i';
								else 
									*os << d;
							}
						}
						*os << '\n';
					}
				}
			}
		}
	}
	else if (errStr) {
		if (treeDesc.SomeEdgesHaveLengths())
			*errStr << "Not all of the branches in tree " << treeN + 1 << " have edge lengths.\n";
		else
			*errStr << "Tree " << treeN + 1 << " does not have edge lengths.\n";
	}
	return false;
}

void filepathToPatristic(const char * filename, ostream *os)
{
	assert(filename);
	BlockReaderList blocks;
	try{
		ExceptionRaisingNxsReader nexusReader(NxsReader::WARNINGS_TO_STDERR);
		if (gStrictLevel != 2)
			nexusReader.SetWarningToErrorThreshold((int)NxsReader::FATAL_WARNING + 1 - (int) gStrictLevel);
		ifstream inf(filename, ios::binary);
		if (!inf.good()) {
			NxsString err;
			err << "Could not parse the file \"" << filename <<"\"";
			nexusReader.NexusError(err, 0, -1, -1);
		}
		cerr << "Creating token" <<endl;
		NxsToken token(inf);	
		NxsCloneBlockFactory factory;

		nexusReader.AddFactory(&factory);
		NxsCharactersBlock charsB(NULL, NULL);
		charsB.SetCreateImpliedBlock(true);
		charsB.SetImplementsLinkAPI(true);
		charsB.SetSupportMixedDatatype(true);
		
		NxsDataBlock dataB(NULL, NULL);
		dataB.SetCreateImpliedBlock(true);
		dataB.SetImplementsLinkAPI(true);
		dataB.SetSupportMixedDatatype(true);
		
		NxsDistancesBlock distancesB(NULL);
		distancesB.SetCreateImpliedBlock(true);
		distancesB.SetImplementsLinkAPI(true);
		
		NxsTaxaBlock taxaB;
		taxaB.SetImplementsLinkAPI(false);
		
		NxsTreesBlock treesB(NULL);
		treesB.SetCreateImpliedBlock(true);
		treesB.SetImplementsLinkAPI(true);
		treesB.SetProcessAllTreesDuringParse(gTreeNum < 0);
		if (gStrictLevel < 2)
			treesB.SetAllowImplicitNames(true);
		treesB.SetWriteFromNodeEdgeDataStructure(true);
		
		NxsUnalignedBlock unalignedB(NULL);
		unalignedB.SetCreateImpliedBlock(true);
		unalignedB.SetImplementsLinkAPI(true);
		
		factory.AddPrototype(&charsB, "CHARACTERS");
		factory.AddPrototype(&dataB, "DATA");
		factory.AddPrototype(&distancesB);
		factory.AddPrototype(&taxaB);
		factory.AddPrototype(&treesB);
		factory.AddPrototype(&unalignedB);

		try {
			cerr << "Executing" <<endl;
			nexusReader.Execute(token);
		}
		catch(...) {
			nexusReader.RemoveFactory(&factory);
			throw;
		}
		nexusReader.RemoveFactory(&factory);
		blocks = nexusReader.GetUsedBlocksInOrder();
		if (os) {
			int prevTrees = 0;
			for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt) {
				NxsBlock * b = *bIt;
				if (b && b->GetID() == "TREES") {
					//*os << "TREES block found" << std::endl;
					
					NxsTreesBlock * treesBPtr = (NxsTreesBlock *) b;
					NxsTaxaBlockAPI * taxaBPtr =  treesBPtr->GetTaxaBlockPtr(NULL);
					if (!taxaBPtr)
						throw NxsException("Trees block is not connected to a taxa block -- I don\'t know how that happened");
					const int nTreesThisBlock = (int) treesBPtr->GetNumTrees();
					if (gTreeNum < 0) {
						for (unsigned i = 0; i < (unsigned) nTreesThisBlock; ++i) {
							const NxsFullTreeDescription & treeDesc = treesBPtr->GetFullTreeDescription(i);
							printPatristic(treeDesc, os, &(std::cerr), *taxaBPtr, i + prevTrees);
							*os << std::endl;
						}
					}
					else if (prevTrees + nTreesThisBlock > gTreeNum) {
						const NxsFullTreeDescription & treeDesc = treesBPtr->GetFullTreeDescription(gTreeNum - prevTrees);
						printPatristic(treeDesc, os, &(std::cerr), *taxaBPtr, gTreeNum);
						break;
					}
				}
			}
		}
		for (BlockReaderList::const_iterator bIt = blocks.begin(); bIt != blocks.end(); ++bIt) {
			NxsBlock * b = *bIt;
			if (b) 
				delete b;
		}
	}
	catch (const NxsException &x) {
		cerr << "Error:\n " << x.msg << endl;
		if (x.line >=0)
			cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << endl;
		exit(2);
	}
}

void readFilepathAsNEXUS(const char *filename) {
	cerr << "[Reading " << filename << "	 ]" << endl;
	try {
		ostream * outStream = &std::cout;
		filepathToPatristic(filename, outStream);
	}
	catch (...) {
		cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
		exit(1);
	}
}	

void readFilesListedIsFile(const char *masterFilepath)
{
	ifstream masterStream(masterFilepath);
	if (masterStream.bad()) {
		cerr << "Could not open " << masterFilepath << "." << endl;
		exit(3);
	}
	char filename[1024];
	while ((!masterStream.eof())  && masterStream.good()) {
		masterStream.getline(filename, 1024);
		if (strlen(filename) > 0 && filename[0] != '#')
			readFilepathAsNEXUS(filename);
	}
}

void printHelp(ostream & out) 
{
	out << "patristicmat reads a NEXUS file.\n\nFor every tree that it encounters with branch lengths, it prints a tab-separated table of patristic distances";
	out << "\nThe most common usage is simply:\n    patristicmat <path to NEXUS file>\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -h on the command line shows this help message\n\n";
	out << "    -t<non-negative integer> specifies the tree to use (overrides the behavior of producing a matrix for every tree).\n";
	out << "    -m enables display of the distance from leaf i to the mrca of i and j in row i and col j.\n";
	out << "    -n disables the  display of the patristic distance matrix.\n";
}

int main(int argc, char *argv[])
{
	//sleep(10);
	if (argc < 2) {
		cerr << "Expecting the path to NEXUS file as the only command line argument!\n" << endl;
		printHelp(cerr);
		return 1;
	}
	for (int i = 1; i < argc; ++i) {
		const char * filepath = argv[i];
		
		if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
			printHelp(cout);
		else if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 't') {
			if (!NxsString::to_long(filepath + 2, &gTreeNum) || gTreeNum < 1) {
				cerr << "Expecting a positive integer after -t\n" << endl;
				printHelp(cerr);
				return 2;
			}
			--gTreeNum;
		}
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'm')
			gMode |= ONLY_TO_MRCA;
		else if (strlen(filepath) == 2 && filepath[0] == '-' && filepath[1] == 'n')
			gMode ^= ONLY_PATRISTIC;
		else if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l')
			readFilesListedIsFile(filepath+2);
		else
			readFilepathAsNEXUS(filepath);
	}
	return 0;
}

