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

#include <cstdio>
#include <cassert>
#include "ncl/ncl.h"
#include "ncl/nxsblock.h"
#include "ncl/nxspublicblocks.h"

#if defined HAVE_CONFIG
#	include "config.h"
#endif
#if ! defined(NCL_NAME_AND_VERSION)
#	define NCL_NAME_AND_VERSION "NCL unknown version"
#endif

void writeAsNexml(PublicNexusReader & nexusReader, ostream & os);


class NexmlIDStrorer;

typedef std::pair<const NxsDiscreteDatatypeMapper *, std::vector<std::string> > MapperStateLabelVec;
template <typename T>
std::string getOrGenId(std::pair<T, unsigned> & p, std::map<T, std::string> & toId, std::map<std::string, T> & toObj, const std::string & pref);

typedef std::pair<const NxsTaxaBlock *, unsigned> TaxaBlockPtrIndPair;
typedef std::pair<const NxsCharactersBlock *, unsigned> CharBlockPtrIndPair;
typedef std::pair<const NxsTreesBlock *, unsigned> TreesBlockPtrIndPair;
typedef std::pair<MapperStateLabelVec, unsigned> MapperStateLabelVecIndPair;

void writeOTUS(ostream & os, const NxsTaxaBlock *, const std::vector<const NxsAssumptionsBlock *> & assumps, NexmlIDStrorer &, unsigned );
void writeCharacters(ostream & os, const NxsCharactersBlock *, const std::vector<const NxsAssumptionsBlock *> & assumps, NexmlIDStrorer &, unsigned);
void writeTrees(ostream & os, const NxsTreesBlock *, const std::vector<const NxsAssumptionsBlock *> & assumps, NexmlIDStrorer &, unsigned);

typedef std::pair<std::string, std::string> AttributeData;
typedef std::vector<AttributeData> AttributeDataVec;

void writeAttribute(std::ostream & out, const AttributeData & aIt);




std::string generateID(std::string prefix, unsigned n)
	{
	char tmp[81];
	sprintf(tmp, "%u", n);
	prefix.append(tmp);
	return prefix;
	}


template <typename T>
std::string getOrGenId(std::pair<T, unsigned> & p, std::map<T, std::string> & toId, std::map<std::string, T> & toObj, const std::string & pref)
{
	typedef typename std::map<T, std::string>::const_iterator ToIDIterator;
	T & t = p.first;
	ToIDIterator f = toId.find(t);
	if (f == toId.end())
		{
		unsigned index = p.second;
		if (index == UINT_MAX)
			{
			cerr << "could not find ID for NEXUS object";
			exit(3);
			}
		std::string id;
		do 
			{
			id = generateID(pref, index++);
			}
		while (toObj.find(id) != toObj.end());
		toId[t] = id;
		toObj[id] = t;
		return id;
		}
	return f->second;
}


class NexmlIDStrorer
	{
	public:
		
		std::string getID(TaxaBlockPtrIndPair taxa)
			{
			const std::string pref("t");
			return getOrGenId<const NxsTaxaBlock *>(taxa, taxaBToID, idToTaxaB, pref);
			}
		std::string getID(CharBlockPtrIndPair chars)
			{
			const std::string pref("c");
			return getOrGenId<const NxsCharactersBlock *>(chars, charsBToID, idToCharsB, pref);
			}
		std::string getID(TreesBlockPtrIndPair trees)
			{
			const std::string pref("g");
			return getOrGenId<const NxsTreesBlock *>(trees, treesBToID, idToTreesB, pref);
			}
		std::string getID(TaxaBlockPtrIndPair taxa, unsigned taxonInd)
			{
			const std::string pref("t");
			std::string p =  getOrGenId<const NxsTaxaBlock *>(taxa, taxaBToID, idToTaxaB, pref);
			p.append(1, 'n');
			return generateID(p, taxonInd);
			}
		std::string getCharID(CharBlockPtrIndPair chars, unsigned charInd)
			{
			const std::string pref("c");
			std::string p =  getOrGenId<const NxsCharactersBlock *>(chars, charsBToID, idToCharsB, pref);
			p.append(1, 'n');
			return generateID(p, charInd);
			}
		std::string getID(CharBlockPtrIndPair chars, unsigned charInd)
			{
			const std::string pref("r");
			std::string p =  getOrGenId<const NxsCharactersBlock *>(chars, charsBToID, idToCharsB, pref);
			p.append(1, 'n');
			return generateID(p, charInd);
			}
		std::string getID(TreesBlockPtrIndPair trees, unsigned treeInd)
			{
			const std::string pref("g");
			std::string p =  getOrGenId<const NxsTreesBlock *>(trees, treesBToID, idToTreesB, pref);
			p.append(1, 'n');
			return generateID(p, treeInd);
			}
		std::string getID(MapperStateLabelVecIndPair m, unsigned sIndex)
			{
			const std::string pref("s");
			std::string p =  getOrGenId<MapperStateLabelVec>(m, mapperToID, idToMapper, pref);
			p.append(1, 'n');
			return generateID(p, sIndex);
			}
			
		
	protected:
		std::map<const NxsTaxaBlock *, std::string> taxaBToID;
		std::map<std::string, const NxsTaxaBlock *> idToTaxaB;
		std::map<const NxsCharactersBlock *, std::string> charsBToID;
		std::map<std::string, const NxsCharactersBlock *> idToCharsB;
		std::map<const NxsTreesBlock *, std::string> treesBToID;
		std::map<std::string, const NxsTreesBlock *> idToTreesB;
		std::map<MapperStateLabelVec, std::string> mapperToID;
		std::map<std::string, MapperStateLabelVec> idToMapper;
	};





inline void writeAttribute(ostream & out, const AttributeData & aIt)
	{
	out << ' ' << aIt.first << '=';
	writeAttributeValue(out, aIt.second);
	}


class XMLElement
{
	public:
	
	XMLElement(const char *name, ostream &outstream, bool hasSubElements, const char *indentStr, const AttributeDataVec *ovec=NULL)
		:out(outstream),
		elementName(name)
		,containsElements(hasSubElements)
		,indentation(indentStr)
		{
		if (ovec)
			this->open(*ovec);
		}
		
	void open()
		{
		const std::vector<AttributeData> atts;
		this->open(atts);
		}

	void open(const std::vector<AttributeData> &atts)
		{
		out << indentation << "<" << this->elementName;
		std::vector<AttributeData>::const_iterator aIt = atts.begin();
		for (; aIt != atts.end(); ++aIt)
			{
			writeAttribute(out, *aIt);
			}
		if (containsElements)
			out << ">\n";
		else
			out << "/>\n";
		
		}
	virtual ~XMLElement()
		{
		if (containsElements)
			out << indentation << "</" << this->elementName << ">\n";
		}
	protected:
		ostream & out;
		const std::string elementName;
		bool containsElements;
		const char *indentation;
};

const char * getNexmlCharPref(NxsCharactersBlock::DataTypesEnum dt)
{
	if (dt == NxsCharactersBlock::standard)
		return "nex:Standard";
	if (dt == NxsCharactersBlock::dna)
		return "nex:Dna";
	if (dt == NxsCharactersBlock::rna)
		return "nex:Rna";
	if (dt == NxsCharactersBlock::protein)
		return "nex:Protein";
	if (dt == NxsCharactersBlock::continuous)
		return "nex:Continuous";
	cerr << "Mixed and Nucleotide data (int " << int(dt) <<") type not supported for nexml output\n";
	exit(2);
}

std::string getNexmlCharSeqType(NxsCharactersBlock::DataTypesEnum dt)
{
	std::string p(getNexmlCharPref(dt));
	p.append("Seqs");
	return p;
}

std::string getNexmlCharCellsType(NxsCharactersBlock::DataTypesEnum dt)
{
	std::string p(getNexmlCharPref(dt));
	p.append("Cells");
	return p;
}

class IDLabelledElement: public XMLElement
{
	public:
		IDLabelledElement(const char *elN, std::string id, std::string titleStr, ostream &outstream, bool contains, const char *indent, const AttributeDataVec *ovec=NULL)
			:XMLElement(elN, outstream, contains, indent)
			{
			AttributeDataVec v;
			v.push_back(AttributeData("id", id));
			if (!titleStr.empty())
				v.push_back(AttributeData("label", titleStr));
			if (ovec)
				v.insert(v.end(), ovec->begin(), ovec->end());
			XMLElement::open(v);
			}
};

class OtusElement: public IDLabelledElement
{
	public:
		OtusElement(std::string id, std::string titleStr, ostream &outstream, const AttributeDataVec *ovec=NULL)
			:IDLabelledElement("otus", id, titleStr, outstream, true, "  ", ovec)
			{}
};

class OtuElement: public IDLabelledElement
{
	public:
		OtuElement(std::string id, std::string titleStr, ostream &outstream, const AttributeDataVec *ovec=NULL)
			:IDLabelledElement("otu", id, titleStr, outstream, false, "    ", ovec)
			{}			
};

class OtherObjLinkedElement : public XMLElement
{
	public:
		OtherObjLinkedElement(const char  * elN, std::string id, std::string titleStr, const char *otherAttN, std::string otherID, ostream &outstream, bool contains, const char * indent, const AttributeDataVec * att)
			:XMLElement(elN, outstream, contains, indent)
			{
			AttributeDataVec v;
			v.push_back(AttributeData("id", id));
			v.push_back(AttributeData(otherAttN, otherID));
			if (!titleStr.empty())
				v.push_back(AttributeData("label", titleStr));
			if (att)
				v.insert(v.end(), att->begin(), att->end());
			XMLElement::open(v);
			}
};

class OTULinkedElement: public OtherObjLinkedElement
{
	public:
		OTULinkedElement(const char  * elN, std::string id, std::string titleStr, std::string taxaBlockID, ostream &outstream, bool contains, const char * indent, const AttributeDataVec * att)
			:OtherObjLinkedElement(elN, id, titleStr, "otu", taxaBlockID, outstream, contains, indent, att)
			{}
};

class OTUSLinkedElement: public OtherObjLinkedElement
{
	public:
		OTUSLinkedElement(const char  * elN, std::string id, std::string titleStr, std::string taxaBlockID, ostream &outstream, bool contains, const char * indent, const AttributeDataVec * att)
			:OtherObjLinkedElement(elN, id, titleStr, "otus", taxaBlockID, outstream, contains, indent, att)
			{}
};
class CharactersElement: public OTUSLinkedElement
{
	public:
		CharactersElement(std::string id, std::string titleStr, std::string taxaBlockID, ostream &outstream, const AttributeDataVec * att)
			:OTUSLinkedElement("characters", id, titleStr, taxaBlockID, outstream, true, "  ", att)
			{}			
};


void writeAsNexml(PublicNexusReader & nexusReader, ostream & os)
{
	std::string version("0.8");
	os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	os << "<nex:nexml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.nexml.org/1.0 ../xsd/nexml.xsd\" xmlns:nex=\"http://www.nexml.org/1.0\" xmlns=\"http://www.nexml.org/1.0\" ";
	os << " xmlns:xml=\"http://www.w3.org/XML/1998/namespace\"  version=\"" << version << "\" generator=\"NEX_us2ml " << NCL_NAME_AND_VERSION << "\">\n";
	
	const unsigned nTaxaBlocks = nexusReader.GetNumTaxaBlocks();
	NexmlIDStrorer memo;
	unsigned nCharBlocksRead = 0;
	unsigned nTreeBlocksRead = 0;
	
	for (unsigned t = 0; t < nTaxaBlocks; ++t)
		{
		const NxsTaxaBlock * tb = nexusReader.GetTaxaBlock(t);

		std::vector<const NxsAssumptionsBlock *> assumps;
		for (unsigned j= 0; j < nexusReader.GetNumAssumptionsBlocks(tb); ++j)
			assumps.push_back(nexusReader.GetAssumptionsBlock(tb, j));

		writeOTUS(os, tb, assumps, memo, t);

		const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(tb);
		for (unsigned i = 0; i < nCharBlocks; ++i)
			{
			NxsCharactersBlock * cb = nexusReader.GetCharactersBlock(tb, i);

			assumps.clear();
			for (unsigned j= 0; j < nexusReader.GetNumAssumptionsBlocks(cb); ++j)
				assumps.push_back(nexusReader.GetAssumptionsBlock(cb, j));
			
			writeCharacters(os, cb, assumps, memo, nCharBlocksRead++);
			}

		const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(tb);
		for (unsigned i = 0; i < nTreesBlocks; ++i)
			{
			NxsTreesBlock * cb = nexusReader.GetTreesBlock(tb, i);

			assumps.clear();
			for (unsigned j= 0; j < nexusReader.GetNumAssumptionsBlocks(cb); ++j)
				assumps.push_back(nexusReader.GetAssumptionsBlock(cb, j));

			writeTrees(os, cb, assumps, memo, nTreeBlocksRead++);
			}
		}
	
	os << "</nex:nexml>\n";
}


void writeOTUS(ostream & os, const NxsTaxaBlock *taxa, const std::vector<const NxsAssumptionsBlock *> & , NexmlIDStrorer &memo, unsigned index)
{
	if (taxa == NULL)
		return;
	TaxaBlockPtrIndPair tbp(taxa, index);
	std::string taxaBlockID = memo.getID(tbp);
	std::string title = taxa->GetTitle();
	OtusElement otus(taxaBlockID, title, os);

	const std::vector<std::string> labels = taxa->GetAllLabels();
	std::vector<std::string>::const_iterator labelIt = labels.begin();
	unsigned n = 0;
	for (; labelIt != labels.end(); ++labelIt, ++n)
		{
		std::string taxonId = memo.getID(tbp, n);
		OtuElement o(taxonId, *labelIt, os);
		}
}


void writeCharLabels(ostream & os, const NxsCharactersBlock *cb, NexmlIDStrorer &memo, unsigned index, std::vector<std::string> * indToStatesID)
{
	const std::string emptyString;
	CharBlockPtrIndPair cbp(cb, index);
	const unsigned nchars = cb->GetNumChar();
	for (unsigned i = 0; i < nchars; ++i)
		{
		std::string label = cb->GetCharLabel(i);
		if (true ) // <char  elements required now!
			{
			std::string charId = memo.getCharID(cbp, i);
			if (indToStatesID)
				{
				std::string labelToShow;
				if (label != " ")
					labelToShow = label;
				std::string statesID = (*indToStatesID)[i];
				AttributeDataVec v(1, AttributeData("states", statesID));
				IDLabelledElement c("char", charId, labelToShow, os, false, "      ", &v);
				}
			else if (label != " ") 
				IDLabelledElement c2("char", charId, label, os, false, "      ");
			else  
				IDLabelledElement c3("char", charId, emptyString, os, false, "      ");
			}
		}
}

std::vector<std::string> getStateLabesVec(const NxsCharactersBlock *cb, unsigned charIndex, const NxsDiscreteDatatypeMapper *mapper)
{
	const unsigned nst = mapper->GetNumStates();
	unsigned i = 0;
	std::vector<std::string> v;
	const std::string emptyString;
	for (unsigned n = 0; n < nst; ++n)
		{
		std::string l = cb->GetStateLabel(charIndex, n);
		if (l != " ")
			{
			if (i > 0)
				{
				for (unsigned j = 0; j < i; ++j)
					v.push_back(emptyString);
				i = 0;
				}
			v.push_back(l);
			}
		else
			i++;
		}
	return v;
}


const char * gDigits = "0123456789";
const char * gLetters = "abcdefghijklmnopqrstuvwxyz";

////////////////////////////////////////////////////////////////////////////////
// Creates a new states element and returns its ID.  
// 	the id's of the each state within the states element will be the states
//	element's ID + 's' + the state index.
//	except that that:
//		the missing state's id will have the last state code's index 
//			(mapper->GetNumStateCodes() - 1).
//		the gap state's id will have the next-to-last state code's index 
//			(mapper->GetNumStateCodes() - 2).
////////////////////////////////////////////////////////////////////////////////
std::string writeStatesElement(const MapperStateLabelVec & mslp, ostream &os, NexmlIDStrorer &memo, unsigned charIndex, unsigned statesIndex, bool useNexusSymbols)
{
	MapperStateLabelVecIndPair m(mslp, charIndex);
	const std::string emptyString;
	const std::string id = memo.getID(m, statesIndex);
	IDLabelledElement statesElement("states", id, emptyString, os, true, "        ", NULL);
	const NxsDiscreteDatatypeMapper * mapper = mslp.first;
	std::vector<std::string> symbols;
	
	// figure out what symbols to use (in almost all cases this will be the NEXUS symbols
	//	but it is possible that (if the matrix was entered in TOKENS mode) there 
	//	will not be NEXUS symbols for all states -- in this case we'll just numbre states.
	
	const unsigned nsc = mapper->GetNumStateCodes();
	const bool hasGap = (mapper->GetNumStatesIncludingGap() > mapper->GetNumStates());
	const int endNum = (hasGap ? nsc - 2: nsc - 1); // minus one for the missing state and one more for the gap state (if used)

	
	int unnumberedCutoff = 10;
	if (useNexusSymbols)
		{
		for (int i = 0; useNexusSymbols && i < endNum; ++i)
			{
			//cerr << i << '\n';
			std::string s = mapper->StateCodeToNexusString(i, false);
			if (s.length() != 1)
				useNexusSymbols = false;
			else
				symbols.push_back(string(1, s[0]));
			}
		if (useNexusSymbols)
			symbols.push_back(string(1, mapper->GetGapSymbol()));
		if (useNexusSymbols)
			symbols.push_back(string(1, mapper->GetMissingSymbol()));
		unnumberedCutoff = 36;
		}
	if (!useNexusSymbols)
		{
		unnumberedCutoff = 10;
		symbols.clear();
		if (nsc <= (unsigned) unnumberedCutoff)
			{
			for (int i = 0; i < 10 && (int) symbols.size() < endNum; ++i)
				symbols.push_back(std::string(1, gDigits[i]));
			for (int i = 0; i < 26 && (int) symbols.size() < endNum; ++i)
				symbols.push_back(std::string(1, gLetters[i]));
			}
		else
			{
			for (int i = 0; i < endNum; ++i)
				symbols.push_back(generateID(emptyString, i));
			}
		/*
		if (hasGap)
			symbols.push_back("-");
		symbols.push_back("?");
		*/
		if (hasGap)
			symbols.push_back(generateID(emptyString, symbols.size()));
		symbols.push_back(generateID(emptyString, symbols.size()));
		}
	/// now we write the "fundamental" states
	const unsigned nst = mapper->GetNumStates();
	const std::vector<std::string> & sl = mslp.second;
	std::string stateIDpref = id;
	stateIDpref.append(1, 's');
	//std::string symStr(symbols.begin(), symbols.end());
	//std::cerr << "symbols = " << symStr << '\n';
	for (unsigned i = 0; i < nst; ++i)
		{
		const std::string stateID = generateID(stateIDpref, i);
		const std::string label = (i < sl.size() ? sl[i] : emptyString);
		AttributeDataVec v(1, AttributeData("symbol", symbols.at(i)));
		IDLabelledElement c("state", stateID, label, os, false, "          ", &v);
		}
	/// now we write the state sets:
	std::string gapStateID;
	if (hasGap)
		gapStateID = generateID(stateIDpref, nsc - 2);
	std::string missingStateID = generateID(stateIDpref, nsc - 1);
	/// now we deal with the gap "state"
	if (hasGap)
		{
		const std::string label("Gap");
		AttributeDataVec v(1, AttributeData("symbol", symbols[nsc - 2]));
		IDLabelledElement c("state", gapStateID, label, os, false, "          ", &v);
		}
	/// now we deal with the gap "state"
	if (hasGap)
		{
		const std::string label("Missing");
		AttributeDataVec v(1, AttributeData("symbol", symbols[nsc - 1]));
		IDLabelledElement c("state", missingStateID, label, os, false, "          ", &v);
		}
	for (int polyuncertain = 0; polyuncertain < 2; ++polyuncertain)
		{
		for (int i = nst; i < endNum; ++i)
			{
			const bool isPoly = mapper->IsPolymorphic(i);
			if ((isPoly && (polyuncertain == 0)) || ((!isPoly) && (polyuncertain == 1)))
				{
				const char * elName = (isPoly ? "polymorphic_state_set" : "uncertain_state_set");
				const std::string stateID = generateID(stateIDpref, i);
				AttributeDataVec v(1, AttributeData("symbol", symbols[i]));
				IDLabelledElement stateSetElement(elName, stateID, emptyString, os, true, "          ", &v);
				const std::set<int>	 & ss = mapper->GetStateSetForCode(i);
				for (std::set<int>::const_iterator subStateIt = ss.begin(); subStateIt != ss.end(); ++subStateIt)
					{
					const int subStateCode = *subStateIt;
					string subStateID;
					if (subStateCode < 0)
						{
						if (subStateCode == NXS_GAP_STATE_CODE)
							subStateID = gapStateID;
						else
							{
							cerr << "unexpected negative state code\n";
							exit(4);
							}
						}
					else
						subStateID = generateID(stateIDpref, (unsigned) subStateCode);
					AttributeDataVec v2(1, AttributeData("state", subStateID));
					XMLElement ss2("member", os, false, "            ", &v2);
					}
				}
			}
		}
	return id;
}


void writeAllStatesElements(
	ostream & os, 
	const NxsCharactersBlock *cb, 
	NexmlIDStrorer &memo, 
	unsigned index, 
	std::vector<std::string> & statesIDVec,
	const bool useNexusSymbols)
{
	const unsigned nchars = cb->GetNumChar();
	typedef std::map<MapperStateLabelVec, std::string> MSLToId;
	MSLToId mappersUsed;
	for (unsigned i = 0; i < nchars; ++i)
		{
		const NxsDiscreteDatatypeMapper * mapper = cb->GetDatatypeMapperForChar(i);
		std::vector<std::string> lv = getStateLabesVec(cb, i, mapper);
		MapperStateLabelVec mslPair(mapper, lv);
		MSLToId::const_iterator prev = mappersUsed.find(mslPair);
		const bool wse = (prev == mappersUsed.end());
		std::string sid;
		if (wse)
			{
			sid = writeStatesElement(mslPair, os, memo, index, mappersUsed.size(), useNexusSymbols);
			mappersUsed[mslPair] = sid;
			}
		else
			sid = prev->second;
		statesIDVec.push_back(sid);
		}
}
void writeCharacters(ostream & os, const NxsCharactersBlock *cb , const std::vector<const NxsAssumptionsBlock *> & , NexmlIDStrorer &memo, unsigned index)
{
	if (!cb)
		return;
	NxsTaxaBlock * taxa  = dynamic_cast<NxsTaxaBlock *>(cb->GetTaxaBlockPtr(NULL));
	if (!taxa)
		return;
	TaxaBlockPtrIndPair tbp(taxa, UINT_MAX);
	std::string taxaBlockID = memo.getID(tbp);
	CharBlockPtrIndPair cbp(cb, index);
	std::string charBlockID = memo.getID(cbp);
	std::string title = cb->GetTitle();
	
	NxsCharactersBlock::DataTypesEnum dt = cb->GetDataType();
	const unsigned nchars = cb->GetNumChar();
		
	if (dt == NxsCharactersBlock::standard || cb->GetNumUserEquates() > 0)
		{
		if (dt != NxsCharactersBlock::standard)
			{
			cerr << "Warning: user defined equates causing the coercion of " << getNexmlCharCellsType(dt) << " type to nex:StandardCells.\n";
			dt = NxsCharactersBlock::standard;
			}

		std::string dtStr = getNexmlCharCellsType(dt);
		AttributeDataVec atts(1, AttributeData("xsi:type", dtStr));
		CharactersElement charEl(charBlockID, title, taxaBlockID, os, &atts);
		std::vector<std::string> statesIDVec;
		if (true)
			{
			XMLElement format("format", os, true, "    ");
			format.open();
			writeAllStatesElements(os, cb, memo, index, statesIDVec, false);
			writeCharLabels(os, cb, memo, index, &statesIDVec);
			}
		if (true)
			{
			XMLElement mat("matrix", os, true, "    ");
			mat.open();
			
			const std::vector<std::string> labels = taxa->GetAllLabels();
			std::vector<std::string>::const_iterator labelIt = labels.begin();
			unsigned n = 0;
			const std::string emptyString;
			for (; labelIt != labels.end(); ++labelIt, ++n)
				{
				if (cb->TaxonIndHasData(n))
					{
					std::string rowId = memo.getID(cbp, n);
					std::string otuId = memo.getID(tbp, n);
					OTULinkedElement row("row", rowId,  emptyString, otuId, os, true, "      ", NULL);
					AttributeDataVec csAtts;
					csAtts.push_back(AttributeData("char",emptyString));
					csAtts.push_back(AttributeData("state",emptyString));
					AttributeData & charAttribute = csAtts[0];
					AttributeData & stateAttribute = csAtts[1];
					const NxsDiscreteStateRow & dataRow =  cb->GetDiscreteMatrixRow(n);
					for (unsigned k = 0; k < nchars; ++k)
						{
						charAttribute.second = memo.getCharID(cbp, k);
						const int sc = dataRow[k];
						unsigned nexmlStatesIndex = 0;
						if (sc >= 0)
							nexmlStatesIndex = (unsigned) sc; 
						else
							{
							const NxsDiscreteDatatypeMapper * mapper = cb->GetDatatypeMapperForChar(k);
							const unsigned nsc = mapper->GetNumStateCodes();
							if (sc == NXS_GAP_STATE_CODE)
								nexmlStatesIndex = nsc - 2; 
							else if (sc == NXS_MISSING_CODE)
								nexmlStatesIndex = nsc - 1; 
							else
								{
								cerr << "Unknown state code " << sc << '\n';
								exit(5);
								}
							}
						std::string stateIDpref = statesIDVec[k];
						stateIDpref.append(1, 's');
						stateAttribute.second = generateID(stateIDpref, nexmlStatesIndex);
						XMLElement("cell", os, false, "        ", &csAtts);
						}
					}
				}
			}
		}
	else
		{
		std::string dtStr = getNexmlCharSeqType(dt);
		AttributeDataVec atts(1, AttributeData("xsi:type", dtStr));
		CharactersElement charEl(charBlockID, title, taxaBlockID, os, &atts);
		if (true) // cb->HasCharLabels()) // new xsd requires format
			{
			XMLElement format("format", os, true, "    ");
			format.open();
			std::vector<std::string> statesIDVec;
			writeAllStatesElements(os, cb, memo, index, statesIDVec, true);
			writeCharLabels(os, cb, memo, index, &statesIDVec);
			}
		if (true)
			{
			XMLElement mat("matrix", os, true, "    ");
			mat.open();
			
			const std::vector<std::string> labels = taxa->GetAllLabels();
			std::vector<std::string>::const_iterator labelIt = labels.begin();
			unsigned n = 0;
			const std::string emptyString;
			for (; labelIt != labels.end(); ++labelIt, ++n)
				{
				if (cb->TaxonIndHasData(n))
					{
					std::string rowId = memo.getID(cbp, n);
					std::string otuId = memo.getID(tbp, n);
					OTULinkedElement row("row", rowId,  emptyString, otuId, os, true, "      ", NULL);
					if (true)
						{
						os << "        <seq>";
						cb->WriteStatesForTaxonAsNexus(os, n, 0, nchars);
						os << "</seq>\n";
						}
					}
				}
			}
		}
}

std::string writeSimpleNode(ostream & os, const NxsSimpleNode &nd, const TaxaBlockPtrIndPair & taxa, unsigned nodeIndex, NexmlIDStrorer &memo, AttributeDataVec*oatts)
{
	AttributeDataVec v;
	std::string prefix("n");
	unsigned otuInd = nd.GetTaxonIndex();
	std::string otuID;
	std::string label;
	std::string id = generateID(prefix, nodeIndex);
	if (otuInd != UINT_MAX)
		v.push_back(AttributeData("otu", memo.getID(taxa, otuInd)));
	else
		label = nd.GetName();
	if (oatts)
		v.insert(v.end(), oatts->begin(), oatts->end());
	IDLabelledElement nodeEl ("node", id, label, os, false, "      ", &v);
	return id;
}

std::string writeSimpleEdge(ostream & os, const NxsSimpleNode *nd, std::map<const NxsSimpleNode *, std::string>  & ndToIdMap, bool edgesAsIntegers)
{
	const NxsSimpleEdge edge = nd->GetEdgeToParent();
	bool defEdgeLen = edge.EdgeLenIsDefaultValue();
	assert(edge.GetChild() == nd);
	std::string eid(1, 'e');
	const std::string & nid = ndToIdMap[nd];
	eid.append(nid);
	NxsString lenstring;
	AttributeDataVec v;
	if (edgesAsIntegers)
		lenstring << edge.GetIntEdgeLen();
	else
		lenstring << edge.GetDblEdgeLen();
	if (!defEdgeLen)
		v.push_back(AttributeData("length", lenstring));
	v.push_back(AttributeData("target", nid));
	const NxsSimpleNode * par = edge.GetParent();
	assert(par);
	assert(ndToIdMap.find(par) != ndToIdMap.end());
	v.push_back(AttributeData("source", ndToIdMap[par]));
	IDLabelledElement edgeEl("edge", eid, std::string(), os, false, "      ", &v);
	return eid;
}
void writeTrees(ostream & os, const NxsTreesBlock *tb, const std::vector<const NxsAssumptionsBlock *> & , NexmlIDStrorer &memo, unsigned index)
{
	if (!tb)
		return;
	const NxsTaxaBlock * taxa  = dynamic_cast<const NxsTaxaBlock *>(tb->GetTaxaBlockPtr(NULL));
	if (!taxa)
		return;
	TaxaBlockPtrIndPair tbp(taxa, UINT_MAX);
	std::string taxaBlockID = memo.getID(tbp);
	TreesBlockPtrIndPair treesbp(tb, index);
	std::string treesBlockID = memo.getID(treesbp);
	std::string title = tb->GetTitle();
	const unsigned ntrees = tb->GetNumTrees();
	if (ntrees == 0)
		return;
	OTUSLinkedElement treesEl("trees", treesBlockID, title, taxaBlockID, os, true, "  ", NULL);
	tb->ProcessAllTrees();
	for (unsigned treen = 0; treen < ntrees; ++treen)
		{
		const NxsFullTreeDescription &ftd = tb->GetFullTreeDescription(treen);
		const bool edgesAsIntegers = ftd.EdgeLengthsAreAllIntegers();
		const char * treeType = (edgesAsIntegers ?  "nex:IntTree": "nex:FloatTree" );
		std::string id = memo.getID(treesbp, treen);
		AttributeDataVec treeAtts(1, AttributeData("xsi:type", std::string(treeType)));
		IDLabelledElement treeEl("tree", id, ftd.GetName(), os, true, "    ", &treeAtts);
		NxsSimpleTree tree(ftd, INT_MAX, DBL_MAX);
		std::vector<const NxsSimpleNode *> preorder = tree.GetPreorderTraversal();
		std::vector<const NxsSimpleNode *>::const_iterator ndIt = preorder.begin();
		std::map<const NxsSimpleNode *, std::string> nodeToIDMap;
		unsigned nodeIndex = 0;
		if (ndIt != preorder.end())
			{
			AttributeDataVec rootAtts;
			string rv(ftd.IsRooted() ? "true" : "false");
			rootAtts.push_back(AttributeData("root", rv));
			const NxsSimpleNode * nd = *ndIt;
			nodeToIDMap[nd] = writeSimpleNode(os, *nd, tbp, nodeIndex++, memo, &rootAtts);
			++ndIt;
			for (; ndIt != preorder.end(); ++ndIt)
				{
				nd = *ndIt;
				nodeToIDMap[nd] = writeSimpleNode(os, *nd, tbp, nodeIndex++, memo, NULL);
				}
			}
		ndIt = preorder.begin();
		nodeIndex = 0;
		if (ndIt != preorder.end())
			{
			const NxsSimpleNode * nd = *ndIt;
			const NxsSimpleEdge edge = nd->GetEdgeToParent();
			bool defEdgeLen = edge.EdgeLenIsDefaultValue();
			if (!defEdgeLen)
				{
				std::string eid(1, 'e');
				const std::string & nid = nodeToIDMap[nd];
				eid.append(nid);
				NxsString lenstring;
				AttributeDataVec v;
				if (edgesAsIntegers)
					lenstring << edge.GetIntEdgeLen();
				else
					lenstring << edge.GetDblEdgeLen();
				v.push_back(AttributeData("length", lenstring));
				v.push_back(AttributeData("target", nid));
				IDLabelledElement edgeEl("rootedge", eid, std::string(), os, false, "      ", &v);
				}
			++ndIt;
			for (; ndIt != preorder.end(); ++ndIt)
				{
				nd = *ndIt;
				writeSimpleEdge(os, nd, nodeToIDMap, edgesAsIntegers);
				}
			}
		}

}
