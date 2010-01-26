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

#include "ncl/ncl.h"
#include "ncltest.h"
#include <cassert>
/*******************************************************************************
 * 	This file uses the version 2.0 API for querying characters blocks to write
 * 		out a NEXUS version of the input file.
 *
 *	It is intended to be a mechanism for regression testing to protect that API
 *
 *	The normalizer program is intended to provide more robust round-trip 
 */
long gStrictLevel = 2;


class NormalizingReader : public NxsReader
	{
	public:
		ostream * outf;
		NormalizingReader(ostream * outStream)
			:NxsReader(),
			outf(outStream),
			storingBlocks(false),
			supportingSharedImpliedBlocks(false)
			{
			}

		~NormalizingReader();
	void Clear();
	void ExecuteStarting() {}
	void ExecuteStopping() {}

	bool EnteringBlock(NxsString )
		{
		// Returning true means it is ok to delete any data associated with previous blocks of this type
		return true;	
		}
	
	void SkippingBlock(NxsString blockName)
		{
		if (outf != 0L)
			{
			*outf << "[!Skipping unknown block (" << blockName << ")...]\n";
			outf->flush();
			}
		}

	void SkippingDisabledBlock(NxsString blockName) 
		{
		if (outf != 0L)
			{
			*outf << "[!Skipping disabled block (" << blockName << ")...]\n";
			outf->flush();
			}
		}

	void PostBlockReadingHook(NxsBlock & block);
	void NexusError(NxsString msg, file_pos pos, long line, long col)
		{
		cerr << "\nError found at line " << line << ", column " << col ;
		cerr << " (file position " << pos << "):\n" << msg << endl;

		if (outf != 0L)
			{
			*outf << "\nError found at line " << line << ", column " << col ;
			*outf << " (file position " << pos << "):\n" << msg;
			outf->flush();
			}
		exit(2);
		}
	private:
		
		NxsBlockMapper 		blockMapper;
		VecBlockPtr			blocksToDelete;
		std::set<BlockUniqueID> writtenBlocks;
		/*  need to make this true for a full normalizer that can rearrange block order.
			The storingBlocks=true branch of the code is only half-written.
			Do NOT change the default in the ctor and expect it to work!!!
		*/
		bool				storingBlocks; 
		/*	The code has only been tested with supportingSharedImpliedBlocks=false */
		bool				supportingSharedImpliedBlocks;
	};

class NormalizingToken : public NxsToken
	{
	public:
		NormalizingToken(istream &is, ostream * os)
			:NxsToken(is),
			out(os)
			{
			}
		void OutputComment(const NxsString & msg)
			{
			if (out == NULL)
				return;
			*out << "[!" << msg << "]\n";
			out->flush();
			}
	private:
		ostream * out;
	};

NormalizingReader::~NormalizingReader()
	{
	Clear();
	}
	
void NormalizingReader::Clear()
	{
	writtenBlocks.clear();
	blockMapper.clear();
	for (VecBlockPtr::iterator it = blocksToDelete.begin(); it != blocksToDelete.end(); ++it)
		delete *it;
	}

void WriteCharactersBlock(ostream & out, NxsCharactersBlock *ncb)
	{
	const unsigned nt = ncb->GetNTax();
	const unsigned nc = ncb->GetNChar();
	int status;
	NxsTaxaBlockAPI* tb = ncb->GetTaxaBlockPtr(&status);
	assert(tb);
	out << "BEGIN DATA;\n";
	out << "    DIMENSIONS NTax="<< ncb->GetNTaxWithData() << " NChar=" << ncb->GetNumActiveChar() << ";\n";
	out << "    FORMAT Datatype = "<< ncb->GetDatatypeName() << " missing=" << ncb->GetMissingSymbol();
	const char gap = ncb->GetGapSymbol();
	if (gap != '\0')
		out << " gap=" << gap;
	out << ";\nMatrix\n";
	for (unsigned i = 0; i < nt; ++i)
		{
		if (ncb->TaxonIndHasData(i))
			{
			out << NxsToken::EscapeString(tb->GetTaxonLabel(i)) << ' ';
			for (unsigned j = 0; j < nc; ++j)
				{
				const unsigned ns = ncb->GetNumStates(i, j);
				if (ns > 1)
					{
					out << '{';
					for (unsigned k = 0; k < ns; ++k)
						out << ncb->GetState(i, j, k);
					out << '}';
					}
				else
					out << ncb->GetState(i, j, 0);
				}
			out << '\n';
			}
		}
	out << ";\nEND;\n";
	}


void processNexusBlock(ostream & out, NxsBlock *nb)
	{
	NxsString id = nb->GetID();
	if (NxsString::case_insensitive_equals(id.c_str(), "DATA")
		|| NxsString::case_insensitive_equals(id.c_str(), "CHARACTERS"))
		{
		NxsCharactersBlock * ncb = static_cast<NxsCharactersBlock *>(nb);
		WriteCharactersBlock(out, ncb);
		}
	else if (NxsString::case_insensitive_equals(id.c_str(), "SETS")
		|| NxsString::case_insensitive_equals(id.c_str(), "ASSUMPTIONS"))
		{}
	else
		nb->WriteAsNexus(out);
	}
void NormalizingReader::PostBlockReadingHook(NxsBlock & block)
	{
	//cerr << "PostBlockReadingHook for " << GetBlockIDTitleString(block) << endl;

	NxsBlock * toStore;
	if (this->storingBlocks)
		{
		try {
			toStore = block.CloneBlock(blockMapper);
			}
		catch (...)
			{
			cerr << "CloneBlock of " << GetBlockIDTitleString(block) << " failed with an exception.  Only Clonable blocks can be normalized." << endl;
			throw;
			}
		blocksToDelete.push_back(toStore);
		blockMapper[block.GetInstanceIdentifier()] = toStore;
		}
	else
		toStore = &block;
	/*
	VecBlockPtr impliedBlocks = toStore->GetImpliedBlocks();
	for (VecBlockPtr::const_iterator ibIt = impliedBlocks.begin(); ibIt != impliedBlocks.end(); ++ibIt)
		{
		const BlockUniqueID currBlockID = (*ibIt)->GetInstanceIdentifier();
		if ((!this->supportingSharedImpliedBlocks) || writtenBlocks.find(currBlockID) == writtenBlocks.end())
			{
			if (this->storingBlocks)
				blocksToDelete.push_back(*ibIt);
			try {
				if (outf != 0L)
					processNexusBlock(*outf, *ibIt);
				}
			catch (...)
				{
				cerr << block.GetInstanceName() << "raised an exception when writing as NEXUS." << endl;
				throw;
				}
			if (this->supportingSharedImpliedBlocks) 
				writtenBlocks.insert(currBlockID);
			}
		}
	*/
	const BlockUniqueID toStoreBlockID = toStore->GetInstanceIdentifier();
	try {
		if (outf != 0L)
			processNexusBlock(*outf, toStore);
		}
	catch (...)
		{
		cerr << GetBlockIDTitleString(block) << " raised an exception when writing as NEXUS." << endl;
		throw;
		}
	writtenBlocks.insert(toStoreBlockID);
	}

void filepathToNormalizedNEXUS(const char * filename, ostream *os)
	{
	ifstream inf(filename, ios::binary);
	toNormalizedNEXUS(inf, os);
	}


void readFilesListedIsFile(const char *masterFilepath, ostream *out)
{
	ifstream masterStream(masterFilepath);
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
			filepathToNormalizedNEXUS(filename, out);
	}
}

void toNormalizedNEXUS(ifstream & inf, ostream *os)
	{
	NxsTaxaBlock taxa;
	NxsTreesBlock trees(&taxa);
	NxsAssumptionsBlock assumptions(&taxa);
	NxsCharactersBlock character(&taxa, &assumptions);
	NxsDataBlock data(&taxa, &assumptions);
	NxsDistancesBlock distances(&taxa);
	NxsUnalignedBlock unaligned(&taxa);

	NormalizingReader nexus(os);
	nexus.Add(&taxa);
	nexus.Add(&trees);
	nexus.Add(&assumptions);
	nexus.Add(&character);
	nexus.Add(&data);
	nexus.Add(&distances);
	nexus.Add(&unaligned);

	if (os)
		{
		*os << "#NEXUS\n";
		}
	NormalizingToken token(inf, os);
	
	nexus.Execute(token);
	}

void printHelp(ostream & out)
	{
	out << "NCLtest NEXUS normalizer takes reads a NEXUS file and rewrites the file to standard out with consistent indentation and syntax.\n";
	out << "\nThe most common usage is simply:\n    NEXUSnormalizer <path to NEXUS file>\n";
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
	}

int main(int argc, char *argv[])
	{
	if (argc < 2)
		{
		cerr << "Expecting the path to NEXUS file as the only command line argument" << endl;
		return 1;
		}
	for (int i = 1; i < argc; ++i)
		{
		char * filename = argv[1];
		cerr << "[Reading " << filename << "     ]" << endl;
		try {
#			if defined(JUST_VALIDATE_NEXUS) && JUST_VALIDATE_NEXUS
				ostream * outStream = 0L;
#			else	
				ostream * outStream = &cout;
#			endif	
			
			const char * filepath = argv[i];
			if (strlen(filepath) > 1 && filepath[0] == '-' && filepath[1] == 'h')
				printHelp(cout);
			else if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 's')
				{
				if (!NxsString::to_long(filepath + 2, &gStrictLevel))
					{
					cerr << "Expecting an integer after -s\n" << endl;
					printHelp(cerr);
					return 2;
					}
				}
			else if (strlen(filepath) > 2 && filepath[0] == '-' && filepath[1] == 'l')
					readFilesListedIsFile(filepath + 2, outStream);
			else
				filepathToNormalizedNEXUS(argv[i], outStream);
			}
		catch (...) 
			{
			cerr << "Normalizing of " << filename << " failed (with an exception)" << endl;
			exit(1);
			}
		}
	return 0;
	}

