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

#include "ncl.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes id to "TAXA" and ntax to 0.
*/
NxsTaxaBlock::NxsTaxaBlock()
  : NxsBlock()
	{
	ntax	= 0;
	id		= "TAXA";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Erases taxonLabels vector.
*/
NxsTaxaBlock::~NxsTaxaBlock()
	{
	taxonLabels.erase(taxonLabels.begin(), taxonLabels.end());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the ability to read everything following the block name (which is read by the NxsReader 
|	object) to the end or endblock statement. Characters are read from the input stream in. Overrides the abstract 
|	virtual function in the base class. 
*/
void NxsTaxaBlock::Read(
  NxsToken &token)	/* the token used to read from in */
	{
	ntax				= 0;
	unsigned nominal_ntax	= 0;
	isEmpty				= false;
	isUserSupplied		= true;

	DemandEndSemicolon(token, "BEGIN TAXA");

	for (;;)
		{
		token.GetNextToken();

		if (token.Equals("DIMENSIONS"))
			{
			// This should be the NTAX keyword
			//
			token.GetNextToken(); 

			if (!token.Equals("NTAX"))
				{
				errormsg = "Expecting NTAX keyword, but found ";
				errormsg += token.GetToken();
				errormsg += " instead";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			DemandEquals(token, "after NTAX");
			nominal_ntax = DemandPositiveInt(token, "NTAX");
			DemandEndSemicolon(token, "DIMENSIONS");
			}	// if (token.Equals("DIMENSIONS"))

		else if (token.Equals("TAXLABELS")) 
			{
			if (nominal_ntax <= 0) 
				{
				errormsg = "NTAX must be specified before TAXLABELS command";
				throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
				}

			for (unsigned i = 0; i < nominal_ntax; i++)
				{
				token.GetNextToken();
				//@pol should check to make sure this is not punctuation
				AddTaxonLabel(token.GetToken());
				}
			DemandEndSemicolon(token, "TAXLABELS");
			}	// if (token.Equals("TAXLABELS")) 

		else if (token.Equals("END") || token.Equals("ENDBLOCK"))
			{
			DemandEndSemicolon(token, "ENDBLOCK");
			break;
			}	// if (token.Equals("END") || token.Equals("ENDBLOCK"))

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
			}	// token not END, ENDBLOCK, TAXLABELS, or DIMENSIONS
		}	// GetNextToken loop
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function outputs a brief report of the contents of this taxa block. Overrides the abstract virtual function in
|	the base class.
*/
void NxsTaxaBlock::Report(
  ostream &out)	/* the output stream to which to write the report */
	{
	out << endl;
	out << id << " block contains ";

	if (ntax == 0)
		{
		out << "no taxa" << endl;
		}
	else if (ntax == 1)
		out << "one taxon" << endl;
	else
		out << ntax << " taxa" << endl;

	if (ntax == 0)
		return;

	for (unsigned k = 0; k < ntax; k++)
		{
		out << '\t' << (k+1) << '\t' << taxonLabels[k] << endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Writes contents of this block in NEXUS format to `out'.
*/
void NxsTaxaBlock::WriteAsNexus(std::ostream &out) const
	{
	out << "BEGIN TAXA;\n";
	this->WriteTitleCommand(out);
	out << "\tDIMENSIONS NTax = " << this->ntax << ";\n\t";
	this->WriteTaxLabelsCommand(out);
	out << "\nEND;\n";
	}

void NxsTaxaBlock::WriteTaxLabelsCommand(std::ostream &out) const
	{
	const unsigned nLabels = this->GetNumTaxonLabels();
	if (nLabels > 0)
		{
		out << "TAXLABELS";
		for (unsigned k = 0; k < this->ntax; ++k)
			out << ' ' << NxsToken::EscapeString(this->taxonLabels[k]) ;
		out << ";";
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Flushes taxonLabels and sets ntax to 0 in preparation for reading a new TAXA block.
*/
void NxsTaxaBlock::Reset()
	{
	NxsBlock::Reset();
	ntax			= 0;
	taxonLabels.clear();
	needsQuotes.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds taxon label 's' to end of list of taxon labels and increments ntax by 1. Returns index of taxon label just 
|	added.
*/
unsigned NxsTaxaBlock::AddTaxonLabel(
  NxsString s)	/* the taxon label to add */
	{
	isEmpty = false;
	if (s.QuotesNeeded())
		needsQuotes.push_back(true);
	else
		needsQuotes.push_back(false);
	
	taxonLabels.push_back(s);
	ntax++;
	return (ntax-1);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Changes the label for taxon 'i' to 's'.
*/
void NxsTaxaBlock::ChangeTaxonLabel(
  unsigned i,	/* the taxon label number to change */
  NxsString s)	/* the string used to replace label i */
	{
	assert(i < (unsigned)taxonLabels.size());

	if (s.QuotesNeeded())
		needsQuotes[i] = true;
	else
		needsQuotes[i] = false;

	taxonLabels[i] = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the length of the longest taxon label stored. Useful for formatting purposes in outputting the data matrix 
|	(i.e., you want the left edge of the matrix to line up).
*/
unsigned NxsTaxaBlock::GetMaxTaxonLabelLength()
	{
	assert(ntax == (unsigned)taxonLabels.size());

	unsigned maxlen = 0;
	for (unsigned i = 0; i < ntax; i++)
		{
		unsigned thislen = (unsigned)taxonLabels[i].size();
		if (thislen > maxlen)
			maxlen = thislen;
		}
	return maxlen;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the label for taxon 'i'.
*/
NxsString NxsTaxaBlock::GetTaxonLabel(
  unsigned i)	/* the taxon label number to return */
	{
	assert(i < (unsigned)taxonLabels.size());

	return taxonLabels[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if taxonLabels[i] contains embedded spaces and thus should be surrounded by single quotes if output is
|	NEXUS format.
*/
bool NxsTaxaBlock::NeedsQuotes(
  unsigned i)	/* the taxon label number in question */
	{
	assert(i < (unsigned)taxonLabels.size());

	return needsQuotes[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if taxon label equal to 's' can be found in the taxonLabels list, and returns false otherwise.
*/
bool NxsTaxaBlock::IsAlreadyDefined(
  NxsString s)	/* the s to attempt to find in the taxonLabels list */
	{
	NxsStringVector::const_iterator iter = find(taxonLabels.begin(), taxonLabels.end(), s);
	bool taxonLabelFound = (iter != taxonLabels.end());
	return taxonLabelFound;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns index of taxon named 's' in taxonLabels list. If taxon named 's' cannot be found, or if there are no 
|	labels currently stored in the taxonLabels list, throws NxsX_NoSuchTaxon exception.
*/
unsigned NxsTaxaBlock::FindTaxon(
  const NxsString &s) const /* the string to attempt to find in the taxonLabels list */
	{
	unsigned k = 0;
	NxsStringVector::const_iterator i;
	for (i = taxonLabels.begin(); i != taxonLabels.end(); ++i)
		{
		if (*i == s)
			return k;
		k++;
		}
	throw NxsTaxaBlock::NxsX_NoSuchTaxon();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of taxon labels currently stored.
*/
unsigned NxsTaxaBlock::GetNumTaxonLabels() const
	{
	return (unsigned)taxonLabels.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets ntax to n.
*/
void NxsTaxaBlock::SetNtax(
  unsigned n)	/* the number of taxa */
	{
	ntax = n;
	}
