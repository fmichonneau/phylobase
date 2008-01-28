//	Copyright (C) 2007 Paul O. Lewis
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

#ifndef NCL_NXSUNALIGNEDBLOCK_INL
#define NCL_NXSUNALIGNEDBLOCK_INL

//@POL Note: This file is not yet ready for use (Paul Lewis, 19-May-2007)

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes taxon whose index in the TAXA block is `i'. If taxon has already been deleted, or if no data is stored for
|	this taxon in this UNALIGNED block, this function has no effect.
*/
inline void NxsUnalignedBlock::DeleteTaxon(
  unsigned i)	/* index of taxon to delete in range [0..`ntaxTotal') */
	{
	assert(i < ntaxTotal);
	unsigned row = taxonPos[i];
	if (row < ntax)
		activeTaxon[row] = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Restores taxon whose index in the TAXA block is `i'. If taxon is already active, or if no data is stored for this
|	taxon in this UNALIGNED block, this function has no effect.
*/
inline void NxsUnalignedBlock::RestoreTaxon(
  unsigned i)	/* index of taxon to restore in range [0..`ntaxTotal') */
	{
	assert(i < ntaxTotal);
	unsigned row = taxonPos[i];
	if (row < ntax)
		activeTaxon[row] = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if taxon `i' is active. If taxon `i' has been deleted, returns false. Also returns false if there is
|	no data for this taxon stored in `umatrix'. Note that `i' should be the index of the taxon in the TAXA block. This 
|	may be equal to this taxon's row in `umatrix', but not necessarily, because data may have been provided (in this 
|	UNALIGNED block) for fewer taxa than appear in the TAXA block.
*/
inline bool NxsUnalignedBlock::IsActiveTaxon(
  unsigned i)	/* the taxon in question, in the range [0..`ntaxTotal') */
	{
	assert(i < ntaxTotal);
	unsigned row = taxonPos[i];
	if (row > ntax)
		return false;
	else
		return activeTaxon[row];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if taxon number `i' has been deleted, false otherwise. This function simply calls IsActiveTaxon, so
|	please also see the documentation for the NxsUnalignedBlock::IsActiveTaxon function.
*/
inline bool NxsUnalignedBlock::IsDeleted(
  unsigned i)	/* the taxon in question, in the range [0..`ntaxTotal') */
	{
	return !IsActiveTaxon(i);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `activeTaxon' data member (pointer to first element of the `activeTaxon' array). Access to this protected 
|	data member is necessary in certain circumstances, such as when a NxsUnalignedBlock object is stored in another 
|	class, and that other class needs direct access to the `activeTaxon' array even though it is not derived from 
|	NxsUnalignedBlock.
*/
inline bool * NxsUnalignedBlock::GetActiveTaxonArray()
	{
	return activeTaxon;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `umatrix' row corresponding to taxon having `taxon_index' in the TAXA block. These may be different if some 
|	taxa were listed in the TAXA block but not in the UNALIGNED block. The parameter `taxon_index' is assumed to range
|	from 0 to `ntaxTotal' - 1.
*/
inline unsigned NxsUnalignedBlock::GetTaxPos(
  unsigned taxon_index)	/* is the index of the taxon in the TAXA block, should be in the range [0..`ntaxTotal') */
	{
	assert(taxonPos != NULL);
	assert(taxon_index != UINT_MAX);
	assert(taxon_index < ntaxTotal);

	return taxonPos[taxon_index];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `datatype' as an unsigned integer. If you want the name of the datatype, you should call 
|	NxsUnalignedBlock::GetDatatypeName instead.
*/
inline unsigned NxsUnalignedBlock::GetDataType()
	{
	return (unsigned)datatype;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns internal representation of the state for taxon `i', character `j', as a vector of integer values. In the 
|	normal situation, there is only one state with no uncertainty or polymorphism and the vector returned will contain
|	only a single, positive integer value. If there are multiple states, the vector will contain k values, where k
|	equals the number of states plus 1. The first value in the vector will be either 0 (indicating ambiguity) or 1
|	(meaning polymorphism), and the remaining values will be positive integers each of which is an index into the 
|	symbols array. In the case of missing data, an empty vector will be returned. In an UNALIGNED block, there are a
|	different number of characters for each taxon. Use NumCharsForTaxon before calling this function to make sure you
|	do not ask for a character beyond the end. If that happens, a NxsUnalignedBlock::NxsX_NoSuchCharacter exception 
|	will be thrown. If no data is stored for taxon `i' in this UNALIGNED block, a NxsUnalignedBlock::NxsX_NoDataForTaxon
|	exception will be thrown, with the exception object storing the offending taxon index in its public data member
|	`taxon_index'.
*/
inline NxsIntVector NxsUnalignedBlock::GetInternalRepresentation(
  unsigned i,	/* is the index of the taxon in the TAXA block in range [0..`ntaxTotal') */
  unsigned j)	/* is the character index (greater than or equal to 0) */
	{
	assert(ntax > 0);
	assert(ntaxTotal > 0);
	assert(i != UINT_MAX);
	assert(i < ntaxTotal);
	unsigned row = taxonPos[i];

	if (row >= ntax)
		{
		throw NxsUnalignedBlock::NxsX_NoDataForTaxon(i);
		}

	// Reminder of what the states data member of each element of umatrix looks like:
	//		     ?              NULL
	//		     G              [1][2]
	//		(AG) polymorphic    [2][0][2][1]
	//		{AG} ambiguous      [2][0][2][0]

	NxsIntVector v;

	UnalignedVect & rowseq = umatrix[row];
	unsigned nchars = (unsigned)rowseq.size();

	// Check to see if user is asking for a character that is beyond the last one 
	if (j > nchars - 1)
		{
		throw NxsUnalignedBlock::NxsX_NoSuchCharacter();
		}

	// Check to see if missing data
	NxsDiscreteDatum & d = umatrix[row][j];
	bool missing = (d.states == NULL);
	if (missing)
		return v;

	// Fill v with the states for this character
	unsigned nstates = d.states[0];
	if (nstates == 1)
		{
		v.push_back(d.states[1]);
		}
	else
		{
		bool poly = (d.states[nstates + 1] == 1);
		v.push_back(poly ? 1 : 0);
		for (unsigned k = 1; k <= nstates; ++k)
			v.push_back(d.states[k]);
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of characters stored for taxon whose index in the TAXA block is `i'. In an UNALIGNED block, 
|	each taxon can have a different number of characters, and this function can be used to find out how many characters
|	are stored for any particular taxon. Note that `i' should be the index of the taxon of interest as it appears in
|	the TAXA block. Because there may be fewer taxa in this UNALIGNED block (`ntax') than there are in the TAXA block
|	(`ntaxTotal'), it is possible that no data were stored for the taxon having index `i', in which case a
|	NxsUnalignedBlock::NxsX_NoDataForTaxon exception is thrown.
*/
inline unsigned NxsUnalignedBlock::NumCharsForTaxon(
  unsigned i)	/* is the index of the taxon in range [0..`ntaxTotal') */
	{
	assert(ntax > 0);
	assert(ntaxTotal > 0);
	assert(i != UINT_MAX);
	assert(i < ntaxTotal);
	unsigned row = taxonPos[i];

	if (row >= ntax)
		{
		throw NxsUnalignedBlock::NxsX_NoDataForTaxon(i);
		}

	// Reminder of what the states data member of each element of umatrix looks like:
	//		     ?              NULL
	//		     G              [1][2]
	//		(AG) polymorphic    [2][0][2][1]
	//		{AG} ambiguous      [2][0][2][0]

	return (unsigned)umatrix[row].size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the missing data symbol currently in effect. If no missing data symbol specified, returns '\0'.
*/
inline char NxsUnalignedBlock::GetMissingSymbol()
	{
	return missing;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ntax', which is the number of rows in `umatrix' (which in turn is the number of taxa for which
|	data is stored). This is not necessarily the same as `ntaxTotal', which is the number of taxa defined in the TAXA
|	block (use the member function GetNTaxTotal to query for that).
*/
inline unsigned NxsUnalignedBlock::GetNTax()
	{
	return ntax;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ntaxTotal', which is the number of taxa defined in the TAXA block. This is not necessarily the
|	same as `ntax', which is the number of taxa for which data were stored in this UNALIGNED block (use the member 
|	function GetNTax to query for that).
*/
inline unsigned NxsUnalignedBlock::GetNTaxTotal()
	{
	return ntaxTotal;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of stored equate associations.
*/
inline unsigned NxsUnalignedBlock::GetNumEquates()
	{
	return (unsigned)equates.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of actual rows in `matrix'. This number is equal to `ntax', and hence this function is identical
|	to GetNTax. Note that `ntax' can be smaller than `ntaxTotal' since the user did not have to provide data for all 
|	taxa specified in the TAXA block.
*/
inline unsigned NxsUnalignedBlock::GetNumMatrixRows()
	{
	return ntax;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of states for taxon `i', character `j'. If `j' is equal to or greater than the number of 
|	characters for taxon `i', returns UINT_MAX. If there is missing data, the return value is 0, otherwise a positive
|	integer will be returned. An alternative is to use the function GetInternalRepresentation to obtain a vector of all
|	states, and the size of that vector could be used to determine both the number and the identity of the states. If
|	no data was stored for the taxon having index i in the UNALIGNED block, a NxsUnalignedBlock::NxsX_NoDataForTaxon 
|	exception is thrown.
*/
inline unsigned NxsUnalignedBlock::GetNumStates(
  unsigned i,	/* the taxon in range [0..`ntaxTotal') */
  unsigned j)	/* the character in range [0..`nchar') */
	{
	assert(ntax > 0);
	assert(ntaxTotal > 0);
	assert(i != UINT_MAX);
	assert(i < ntaxTotal);
	unsigned row = taxonPos[i];

	if (row >= ntax)
		{
		throw NxsUnalignedBlock::NxsX_NoDataForTaxon(i);
		}

	// Reminder of what the states data member of each element of umatrix looks like:
	//		     ?              NULL
	//		     G              [1][2]
	//		(AG) polymorphic    [2][0][2][1]
	//		{AG} ambiguous      [2][0][2][0]

	UnalignedVect & rowseq = umatrix[row];
	unsigned nchars = (unsigned)rowseq.size();

	// Check to see if user is asking for a character that is beyond the last one 
	if (j > nchars - 1)
		{
		return UINT_MAX;
		}

	// Check to see if missing data
	NxsDiscreteDatum & d = umatrix[row][j];
	bool missing = (d.states == NULL);
	if (missing)
		return 0;

	return d.states[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the original taxon number (used in the NEXUS data file) in the range [1..`ntaxTotal']. Will be equal to 
|	`i' + 1 unless data was not provided for some taxa listed in a preceding TAXA block.
*/
inline unsigned NxsUnalignedBlock::GetOrigTaxonNumber(
  unsigned i)	/* the character in range [0..`ntax') */
	{
	if (ntax == ntaxTotal)
		return 1 + i;	// use the fast method if possible
	else
		return (1 + GetOrigTaxonIndex(i));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns data member `symbols'. Warning: returned value may be NULL.
*/
inline char * NxsUnalignedBlock::GetSymbols()
	{
	return symbols;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns label for taxon `i'. Note that `i' should be the index of the taxon in the TAXA block. This may be equal to
|	this taxon's row in `umatrix', but not necessarily, because data may have been provided (in this UNALIGNED block) 
|	for fewer taxa than appear in the TAXA block.
*/
inline NxsString NxsUnalignedBlock::GetTaxonLabel(
  unsigned i)	/* the taxon index in range [0..`ntaxTotal') */
	{
	assert(i < ntaxTotal);
	NxsString s = taxa->GetTaxonLabel(i);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if LABELS was specified in the FORMAT command, false otherwise.
*/
inline bool NxsUnalignedBlock::IsLabels()
	{
	return labels;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the state at taxon `i', character `j' is the missing state, false otherwise. Throws NxsException if
|	`j' is too large (i.e. specifies a character beyond the last character for `umatrix' row `i'). Calls 
|	NxsUnalignedBlock::GetInternalRepresentation, so unless all you need is information about missing data, it is more
|	efficient to simply call GetInternalRepresentation and see if the returned vector is empty. Note that `i' should be
|	the index of the taxon in the TAXA block. If data for that taxon has not been stored in this UNALIGNED block, then
|	a NxsUnalignedBlock::NxsX_NoDataForTaxon exception will be thrown by GetInternalRepresentation.
*/
inline bool NxsUnalignedBlock::IsMissingState(
  unsigned i,	/* the taxon, in range [0..`ntaxTotal') */
  unsigned j)	/* the character, in range [0..infinity) */
	{
	NxsIntVector v = GetInternalRepresentation(i, j);
	return v.empty();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if taxon `i' is polymorphic for character `j', false otherwise. Throws NxsException if `j' is too large
|	(i.e. specifies a character beyond the last character for `umatrix' row `i'). Calls 
|	NxsUnalignedBlock::GetInternalRepresentation, so unless all you need is information about polymorphism, it is more
|	efficient to simply call GetInternalRepresentation and extract the information you need from the returned vector.
|	Note that `i' should be the index of the taxon in the TAXA block. If data for that taxon has not been stored in this
|	UNALIGNED block, then a NxsUnalignedBlock::NxsX_NoDataForTaxon exception will be thrown by 
|	GetInternalRepresentation.
*/
inline bool NxsUnalignedBlock::IsPolymorphic(
  unsigned i,	/* the taxon in range [0..`ntaxTotal') */
  unsigned j)	/* the character in range [0..infinity) */
	{
	NxsIntVector v = GetInternalRepresentation(i, j);
	unsigned sz = (unsigned)v.size();
	bool poly = false;
	if (sz > 1)
		poly = (v[0] == 1);
	return poly;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if RESPECTCASE was specified in the FORMAT command, false otherwise.
*/
inline bool NxsUnalignedBlock::IsRespectCase()
	{
	return respectingCase;
	}

#endif
