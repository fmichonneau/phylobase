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

#include <cassert>
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"
#include "ncl/nxsmultiformat.h"


void writeAsNexml(PublicNexusReader & nexusReader, ostream & os);
void writeAsNexus(PublicNexusReader & nexusReader, ostream & os);


void exportData(PublicNexusReader & nexusReader, MultiFormatReader::DataFormatType f, long interleaveLen, std::string prefix);
void exportCharacters(PublicNexusReader & nexusReader, MultiFormatReader::DataFormatType f, long interleaveLen, std::string prefix);
void exportTrees(PublicNexusReader & nexusReader, MultiFormatReader::DataFormatType f, std::string prefix);
const char * getFileExtension( MultiFormatReader::DataFormatType f);
std::vector<NxsNameToNameTrans> nameTranslationDict(const std::vector<std::string> & origTaxa, MultiFormatReader::DataFormatType f);
std::string getLegalTaxonName(const std::string & origName, const std::set<std::string> & used,  MultiFormatReader::DataFormatType f);
std::string getLegalPhylipTaxonName(const std::string & origName, const std::set<std::string> & used);
std::string getLegalRelaxedPhylipTaxonName(const std::string & origName, const std::set<std::string> & used);
bool formLegalPhylipName(const std::string & origName, const NxsString & numericExtension, const std::set<std::string> & used, std::string & toReturn);
bool formLegalRelaxedPhylipName(const std::string & origName, const NxsString & numericExtension, const std::set<std::string> & used, std::string & toReturn);

std::string purgeIllegalCharactersFromPhylipName(const std::string &origName);
std::string purgeIllegalCharactersFromRelaxedPhylipName(const std::string &origName);


const char * getFileExtension( MultiFormatReader::DataFormatType f) {
	if (f == MultiFormatReader::NEXUS_FORMAT)
		return ".nex";
	else if (f == MultiFormatReader::FASTA_DNA_FORMAT
			|| f == MultiFormatReader::FASTA_AA_FORMAT
			|| f == MultiFormatReader::FASTA_RNA_FORMAT)
		return ".fasta";
	else if (f == MultiFormatReader::PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::PHYLIP_DISC_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DISC_FORMAT
			|| f == MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::RELAXED_PHYLIP_DISC_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT)
		return ".phy";
	else if (f == MultiFormatReader::PHYLIP_TREE_FORMAT
			|| f == MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT)
		return ".tre";
	else if (f == MultiFormatReader::NEXML_FORMAT)
		return ".xml";
	else {
		throw NxsException("export requested for unsupported format");
	}
}

std::string purgeIllegalCharactersFromPhylipName(const std::string &origName)
{
	std::string cleanedName;
	cleanedName.reserve(origName.length());
	for (std::string::const_iterator it = origName.begin(); it != origName.end(); ++it)
		{
		if (strchr("()[];,:\n", *it) == NULL)
			cleanedName.append(1, *it);
		}
	return cleanedName;
}

std::string purgeIllegalCharactersFromRelaxedPhylipName(const std::string &origName)
{
	std::string cleanedName;
	cleanedName.reserve(origName.length());
	for (std::string::const_iterator it = origName.begin(); it != origName.end(); ++it)
		{
		if (strchr(" \n\t", *it) != NULL)
			cleanedName.append(1, '_');
		else if (strchr("()[];,:.", *it) == NULL) // RAxML (at least in the past) has not been happy with .'s in names -- since it is one of the targetted "relaxed phylip" programs, we will suppress .'s
			cleanedName.append(1, *it);
		}
	return cleanedName;
}

bool formLegalPhylipName(const std::string & origName, const NxsString & numericExtension, const std::set<std::string> & used, std::string & toReturn)
{
	const unsigned MAX_PHYLIP_NAME_LENGTH = 10;
	unsigned postLen = numericExtension.length();
	if (postLen > MAX_PHYLIP_NAME_LENGTH)
		throw NxsException("Number of duplicate names exceed the capacity of our poorly thought out mechanism for avoiding name clashes");
	const unsigned unPaddedLen = postLen + origName.length();
	if (unPaddedLen <= MAX_PHYLIP_NAME_LENGTH)
		toReturn = origName;
	else if (postLen == MAX_PHYLIP_NAME_LENGTH)
		toReturn.clear();
	else
		toReturn.assign(origName.c_str(), MAX_PHYLIP_NAME_LENGTH - postLen);
	toReturn.append(numericExtension.c_str());
	if (unPaddedLen < MAX_PHYLIP_NAME_LENGTH)
		toReturn.append(MAX_PHYLIP_NAME_LENGTH - unPaddedLen, ' ');
	const std::string cap = NxsString::get_upper(toReturn);
	return (used.find(cap) == used.end());
}

bool formLegalRelaxedPhylipName(const std::string & origName, const NxsString & numericExtension, const std::set<std::string> & used, std::string & toReturn)
{
	toReturn = origName;
	toReturn.append(numericExtension.c_str());
	const std::string cap = NxsString::get_upper(toReturn);
	return (used.find(cap) == used.end());
}

std::string getLegalPhylipTaxonName(const std::string & origName, const std::set<std::string> & used)
{
	NxsString numericExtension;
	std::string toReturn;
	const std::string cleanedName(purgeIllegalCharactersFromPhylipName(origName));
	
	for (unsigned i = 1;; ++i)
		{
		if (formLegalPhylipName(cleanedName, numericExtension, used, toReturn))
			return toReturn;
		numericExtension.clear();
		numericExtension << i;
		}
}

std::string getLegalRelaxedPhylipTaxonName(const std::string & origName, const std::set<std::string> & used)
{
	NxsString numericExtension;
	std::string toReturn;
	const std::string cleanedName(purgeIllegalCharactersFromRelaxedPhylipName(origName));
	
	for (unsigned i = 1;; ++i)
		{
		if (formLegalRelaxedPhylipName(cleanedName, numericExtension, used, toReturn))
			return toReturn;
		numericExtension.clear();
		numericExtension << i;
		}
}

std::string getLegalTaxonName(const std::string & origName, const std::set<std::string> & used,  MultiFormatReader::DataFormatType f)
{
	if (f == MultiFormatReader::NEXUS_FORMAT 
		|| f == MultiFormatReader::FASTA_DNA_FORMAT
		|| f == MultiFormatReader::FASTA_AA_FORMAT
		|| f == MultiFormatReader::FASTA_RNA_FORMAT
		|| f == MultiFormatReader::NEXML_FORMAT)
		return origName;
	if (f == MultiFormatReader::PHYLIP_DNA_FORMAT
		|| f == MultiFormatReader::PHYLIP_RNA_FORMAT
		|| f == MultiFormatReader::PHYLIP_AA_FORMAT
		|| f == MultiFormatReader::PHYLIP_DISC_FORMAT
		|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT
		|| f == MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT
		|| f == MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT
		|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DISC_FORMAT
		|| f == MultiFormatReader::PHYLIP_TREE_FORMAT) {
		return getLegalPhylipTaxonName(origName, used);
	}
	if (f == MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT
		|| f == MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT
		|| f == MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT
		|| f == MultiFormatReader::RELAXED_PHYLIP_DISC_FORMAT
		|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT
		|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT
		|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT
		|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT
		|| f == MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT) {
		return getLegalRelaxedPhylipTaxonName(origName, used);
	}
	throw NxsException("export requested for unsupported format");
}

/// returns an empty vector if no translation of names is needed or a vector
///	 of the internal name to the getLegaldName
std::vector<NxsNameToNameTrans> nameTranslationDict(const std::vector<std::string> & origTaxa, MultiFormatReader::DataFormatType f)
{
	bool transNeeded = false;
	std::vector<NxsNameToNameTrans> translationTable;
	std::set<std::string> usedNames;
	for (std::vector<std::string>::const_iterator origIt = origTaxa.begin(); origIt != origTaxa.end(); ++origIt) {
		std::string e = getLegalTaxonName(*origIt, usedNames, f);
		if (e != *origIt)
			transNeeded = true;
		translationTable.push_back(NxsNameToNameTrans(*origIt, e));
		NxsString::to_upper(e);
		usedNames.insert(e);
	}
	if (transNeeded)
		return translationTable;
	return std::vector<NxsNameToNameTrans>();
}

bool IsPhylipType(MultiFormatReader::DataFormatType f)
{
	return (f == MultiFormatReader::PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::PHYLIP_DISC_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DISC_FORMAT);
}

bool IsRelaxedPhylipType(MultiFormatReader::DataFormatType f)
{
	return (f == MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::RELAXED_PHYLIP_DISC_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT);
}

bool IsInterleaveType(MultiFormatReader::DataFormatType f)
{
	return (f == MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DISC_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT
			|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT);
}	

bool IsFastaType(MultiFormatReader::DataFormatType f)
{
	return (f == MultiFormatReader::FASTA_DNA_FORMAT
			|| f == MultiFormatReader::FASTA_RNA_FORMAT
			|| f == MultiFormatReader::FASTA_AA_FORMAT);
}

void writeCharactersBlockToStream(
  const NxsCharactersBlock & cb, 
  ostream & outf, 
  const std::vector<std::string> & taxaNames, 
  MultiFormatReader::DataFormatType f, 
  long interleaveLen)
{
	const unsigned nt = taxaNames.size();
	const unsigned nc = cb.GetNChar();
	unsigned nCharsToWrite;
	unsigned seqStartColumn = 0;
	
	if (IsRelaxedPhylipType(f))
		{
		for (unsigned i = 0; i < nt; ++i)
			{
			const std::string & name = taxaNames[i];
			if (name.length() > seqStartColumn)
				seqStartColumn = name.length();
			}
		seqStartColumn += 1;
		}
	
				
	if (IsPhylipType(f) || IsRelaxedPhylipType(f))
		{
		std::string sep;

		outf << nt << ' ' << nc << '\n';

		if (IsInterleaveType(f) && interleaveLen > 0) 
			{
			std::vector<std::string> storedSeqs;
			storedSeqs.reserve(nt);
			std::string *sp;
			nCharsToWrite = (nc > (unsigned) interleaveLen ? (unsigned)interleaveLen : nc);
			for (unsigned i = 0; i < nt; ++i)
				{
				const std::string & name = taxaNames[i];
				storedSeqs.push_back(cb.GetMatrixRowAsStr(i));
				sp = &(storedSeqs[i]);
				if (IsRelaxedPhylipType(f))
					{
					sep.clear();
					sep.append(seqStartColumn - name.length(), ' ');
					}
				outf << name << sep << sp->substr(0, nCharsToWrite) << '\n';
				}
			
			for (unsigned currIndex = (unsigned)interleaveLen; currIndex < nc; currIndex += (unsigned)interleaveLen) 
				{
				outf << '\n';
				nCharsToWrite = ((nc - currIndex) > (unsigned)interleaveLen ? (unsigned)interleaveLen : (nc - currIndex));
				for (unsigned i = 0; i < nt; ++i)
					{
					sp = &(storedSeqs[i]);
					outf << sp->substr(currIndex, nCharsToWrite) << '\n';
					}
				}
			}
		else
			{
			if (interleaveLen > 0) 
				{
				// not interleaved, but wrapping at interleaveLen
				for (unsigned i = 0; i < nt; ++i)
					{
					nCharsToWrite = (nc > (unsigned)interleaveLen ? (unsigned)interleaveLen : nc);
					const std::string & name = taxaNames[i];
					std::string seq = cb.GetMatrixRowAsStr(i);
					if (IsRelaxedPhylipType(f))
						{
						sep.clear();
						sep.append(seqStartColumn - name.length(), ' ');
						}
					outf << name << sep << seq.substr(0, nCharsToWrite) << '\n';
					for (unsigned currIndex = (unsigned)interleaveLen; currIndex < nc; currIndex += (unsigned)interleaveLen) 
						{
						nCharsToWrite = ((nc - currIndex) > (unsigned)interleaveLen ? (unsigned)interleaveLen : (nc - currIndex));
						outf << seq.substr(currIndex, nCharsToWrite) << '\n';
						}
					}
				}
			else 
				{
				// not interleaved, and not wrapping
				for (unsigned i = 0; i < nt; ++i)
					{
					const std::string & name = taxaNames[i];
					std::string seq = cb.GetMatrixRowAsStr(i);
					if (IsRelaxedPhylipType(f))
						{
						sep.clear();
						sep.append(seqStartColumn - name.length(), ' ');
						}
					outf << name << sep << seq << '\n';
					}
				}
			}
		}
	else if (IsFastaType(f))
		{
		if (interleaveLen < 1)
			interleaveLen = 60; // default FASTA line length
		for (unsigned i = 0; i < nt; ++i)
			{
			nCharsToWrite = (nc > (unsigned)interleaveLen ? (unsigned)interleaveLen : nc);
			const std::string & name = taxaNames[i];
			std::string seq = cb.GetMatrixRowAsStr(i);
			outf << '>' << name << '\n' << seq.substr(0, nCharsToWrite) << '\n';
			for (unsigned currIndex = (unsigned)interleaveLen; currIndex < nc; currIndex += (unsigned)interleaveLen) 
				{
				nCharsToWrite = ((nc - currIndex) > (unsigned)interleaveLen ? (unsigned)interleaveLen : (nc - currIndex));
				outf << seq.substr(currIndex, nCharsToWrite) << '\n';
				}
			}
		}
	else 
		{
		throw NxsException("writeCharactersBlockToStream requested for unsupported format");
		}
}

void exportCharacters(
  PublicNexusReader & nexusReader, 
  MultiFormatReader::DataFormatType f, 
  long interleaveLen,
  std::string prefix)
{	
	const unsigned nTaxaBlocks = nexusReader.GetNumTaxaBlocks();
	for (unsigned t = 0; t < nTaxaBlocks; ++t)
		{
		const NxsTaxaBlock * tb = nexusReader.GetTaxaBlock(t);
		const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(tb);
		if (nCharBlocks == 0)
			continue;
		
		NxsString tbSpecificPrefix;
		if (t > 0)
			tbSpecificPrefix << (1 + t);
		tbSpecificPrefix << prefix;
		
		typedef std::pair<const NxsCharactersBlock *, std::string> PairCBAndString;
		typedef std::vector< PairCBAndString > VecPairCBAndString;
		VecPairCBAndString cbToWrite;
		for (unsigned i = 0; i < nCharBlocks; ++i)
			{
			NxsString fn = tbSpecificPrefix;
			if (i > 0)
				fn << (1 + i);
			fn << getFileExtension(f);

			const NxsCharactersBlock * cb = nexusReader.GetCharactersBlock(tb, i);
			
			NxsCharactersBlock::DataTypesEnum dt = cb->GetDataType();
			bool writeBlock = false;
			if ((dt == NxsCharactersBlock::standard)
				&& (f == MultiFormatReader::PHYLIP_DISC_FORMAT
					|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DISC_FORMAT
					|| f == MultiFormatReader::RELAXED_PHYLIP_DISC_FORMAT
					|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT))
				writeBlock = true;
			else if ((dt == NxsCharactersBlock::dna || dt == NxsCharactersBlock::nucleotide)
				&& (f == MultiFormatReader::PHYLIP_DNA_FORMAT
					|| f == MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT
					|| f == MultiFormatReader::FASTA_DNA_FORMAT
					|| f == MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT
					|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT))
				writeBlock = true;
			else if ((dt == NxsCharactersBlock::rna)
				&& (f == MultiFormatReader::PHYLIP_RNA_FORMAT
					|| f == MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT
					|| f == MultiFormatReader::FASTA_RNA_FORMAT
					|| f == MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT
					|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT))
				writeBlock = true;
			else if ((dt == NxsCharactersBlock::protein)
				&& (f == MultiFormatReader::PHYLIP_AA_FORMAT
					|| f == MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT
					|| f == MultiFormatReader::FASTA_AA_FORMAT
					|| f == MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT
					|| f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT))
				writeBlock = true;
				
			if (writeBlock)
				cbToWrite.push_back(std::pair<const NxsCharactersBlock *, std::string>(cb, fn));
			}
		if (!cbToWrite.empty())
			{
			std::vector<std::string> namesToPrint(tb->GetAllLabels());
			std::vector<NxsNameToNameTrans> nameTrans = nameTranslationDict(namesToPrint, f);
			if (!nameTrans.empty())
				{
				namesToPrint.clear();
				for (std::vector<NxsNameToNameTrans>::const_iterator nIt = nameTrans.begin(); nIt != nameTrans.end(); ++nIt)
					namesToPrint.push_back(nIt->second);
				nexusReader.conversionOutputRecord.writeNameTranslation(nameTrans, tb);
				}

			for (VecPairCBAndString::const_iterator vIt = cbToWrite.begin(); vIt != cbToWrite.end(); ++vIt)
				{
				const NxsCharactersBlock * cbP = vIt->first;
				std::string fn = vIt->second;
				ofstream outf;
				outf.open(fn.c_str());
				if (!outf.good())
					{
					NxsString msg; 
					msg << "Could not open the file " << fn;
					throw NxsException(msg);
					}
				std::cerr << "Writing " << fn << '\n';
				writeCharactersBlockToStream(*cbP, outf, namesToPrint, f, interleaveLen);
				}
			}

		}
}

void exportTrees(
  PublicNexusReader & nexusReader, 
  MultiFormatReader::DataFormatType f, 
  std::string prefix)
{	
	const unsigned nTaxaBlocks = nexusReader.GetNumTaxaBlocks();
	for (unsigned t = 0; t < nTaxaBlocks; ++t)
		{
		const NxsTaxaBlock * tb = nexusReader.GetTaxaBlock(t);
		const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(tb);
		if (nTreesBlocks == 0)
			continue;

		NxsString tbSpecificPrefix;
		if (t > 0)
			tbSpecificPrefix << (1 + t);
		tbSpecificPrefix << prefix;
		std::vector<std::string> namesToPrint(tb->GetAllLabels());
		std::vector<NxsNameToNameTrans> nameTrans = nameTranslationDict(namesToPrint, f);
		if (!nameTrans.empty())
			{
			namesToPrint.clear();
			for (std::vector<NxsNameToNameTrans>::const_iterator nIt = nameTrans.begin(); nIt != nameTrans.end(); ++nIt)
				namesToPrint.push_back(nIt->second);
			nexusReader.conversionOutputRecord.writeNameTranslation(nameTrans, tb);
			}
			
		for (unsigned i = 0; i < nTreesBlocks; ++i)
			{
			NxsString fn = tbSpecificPrefix;
			if (i > 0)
				fn << (1 + i);
			fn << getFileExtension(f);

			const NxsTreesBlock * trb = nexusReader.GetTreesBlock(tb, i);
			
			ofstream outf;
			outf.open(fn.c_str());
			if (!outf.good())
				{
				NxsString msg; 
				msg << "Could not open the file " << fn << " to write character data";
				throw NxsException(msg);
				}
			std::cerr << "Writing " << fn << '\n';
			trb->ProcessAllTrees();
			for (unsigned j = 0; j < trb->GetNumTrees(); ++j) 
				{
				const NxsFullTreeDescription & ftd = trb->GetFullTreeDescription(j);
				NxsSimpleTree tree(ftd, -1, -1.0);
				std::vector<NxsSimpleNode *> & leaves = tree.GetLeavesRef();
				for (std::vector<NxsSimpleNode *>::const_iterator leafIt = leaves.begin(); leafIt != leaves.end(); ++leafIt) 
					{
					NxsSimpleNode * leaf = *leafIt;
					assert(leaf);
					const std::string name = namesToPrint[leaf->GetTaxonIndex()];
					std::cerr << "Setting name=" << name << '\n';
					leaf->SetName(name);
					}
				tree.WriteAsNewick(outf, true, true, false, tb);
				outf << ";\n";
				}
			}
		}
}

void exportData(
  PublicNexusReader & nexusReader, 
  MultiFormatReader::DataFormatType f,
  long interleaveLen, 
  std::string prefix)
{
	std::string fullName = prefix;
	if (f == MultiFormatReader::NEXUS_FORMAT) {
		// hack-alert: we don't need to pass in the interleaveLen because it is set as a global in the calling code 
		fullName.append(".nex");
		std::ofstream nexOut(fullName.c_str());
		std::cerr << "Writing " << fullName << '\n';
		writeAsNexus(nexusReader, nexOut);
		nexOut.close();
	}
	else if (f == MultiFormatReader::FASTA_DNA_FORMAT
			 || f == MultiFormatReader::FASTA_AA_FORMAT	
			 || f == MultiFormatReader::FASTA_RNA_FORMAT	
			 || f == MultiFormatReader::PHYLIP_DNA_FORMAT	
			 || f == MultiFormatReader::PHYLIP_RNA_FORMAT	
			 || f == MultiFormatReader::PHYLIP_AA_FORMAT	
//			 || f == MultiFormatReader::PHYLIP_DISC_FORMAT	
			 || f == MultiFormatReader::INTERLEAVED_PHYLIP_DNA_FORMAT	
			 || f == MultiFormatReader::INTERLEAVED_PHYLIP_RNA_FORMAT	
			 || f == MultiFormatReader::INTERLEAVED_PHYLIP_AA_FORMAT	
//			 || f == MultiFormatReader::INTERLEAVED_PHYLIP_DISC_FORMAT	
			 || f == MultiFormatReader::RELAXED_PHYLIP_DNA_FORMAT	
			 || f == MultiFormatReader::RELAXED_PHYLIP_RNA_FORMAT	
			 || f == MultiFormatReader::RELAXED_PHYLIP_AA_FORMAT	
//			 || f == MultiFormatReader::RELAXED_PHYLIP_DISC_FORMAT	
			 || f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DNA_FORMAT	
			 || f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_RNA_FORMAT	
//			 || f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_DISC_FORMAT
			 || f == MultiFormatReader::INTERLEAVED_RELAXED_PHYLIP_AA_FORMAT)
		exportCharacters(nexusReader, f, interleaveLen, prefix);
	else if (f == MultiFormatReader::PHYLIP_TREE_FORMAT
			 || f == MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT)
		exportTrees(nexusReader, f, prefix);
	else if (f == MultiFormatReader::NEXML_FORMAT) {
		fullName.append(".xml");
		std::ofstream nexOut(fullName.c_str());
		std::cerr << "Writing " << fullName << '\n';
		writeAsNexml(nexusReader, nexOut);
		nexOut.close();
	}
	else {
		throw NxsException("export requested for unsupported format");
	}
}
