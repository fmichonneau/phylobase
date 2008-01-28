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
 * 	This file contains code for an executable that takes path to a NEXUS file
 *  	as a command line argument and writes a "normalized" version of the 
 *		blocks.  This is useful for testing. 
 *
 *	Ideally 2 equivalent files will produce the same normalized output. This 
 *		version of tthe program is less ambitious. The goal is to be able to run 
 *		(for any valid NEXUS in.nex file):
 *			$ normalizer in.nex > outOrig.nex
 *			$ normalizer outOrig.nex > outSecond.nex
 *  		$ diff outOrig.nex outSecond.nex
 *		and find no differences.
 */
#include "ncl.h"
#include "normalizer.h"
#include "nxsblock.h"


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
					(*ibIt)->WriteAsNexus(*outf);
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
	const BlockUniqueID toStoreBlockID = toStore->GetInstanceIdentifier();
	try {
		if (outf != 0L)
			toStore->WriteAsNexus(*outf);
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
	
void toNormalizedNEXUS(ifstream & inf, ostream *os)
	{
	NxsTaxaBlock taxa;
	NxsTreesBlock trees(&taxa);
	NxsAssumptionsBlock assumptions(&taxa);
	NxsCharactersBlock character(&taxa, &assumptions);
	NxsDataBlock data(&taxa, &assumptions);
	NxsDistancesBlock distances(&taxa);
	NxsUnalignedBlock unaligned(&taxa, &assumptions);

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

