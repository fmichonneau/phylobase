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

#ifndef NCL_NXSUNALIGNEDBLOCK_H
#define NCL_NXSUNALIGNEDBLOCK_H

//@POL Note: This file is not yet ready for use (Paul Lewis, 19-May-2007)

class NxsTaxaBlock;
class NxsAssumptionsBlock;

/*----------------------------------------------------------------------------------------------------------------------
|	This class handles reading and storage for the NEXUS block UNALIGNED. It overrides the member functions Read and 
|	Reset, which are abstract virtual functions in the base class NxsBlock.
|>
|	Below is a table showing the correspondence between the elements of an UNALIGNED block in a NEXUS file and the
|	variables and member functions of the NxsUnalignedBlock class that can be used to access each piece of information
|	stored. Items in parenthesis should be viewed as "see also" items.
|>
|	NEXUS         Command        Data           Member
|   Command       Atribute       Member         Functions
|	---------------------------------------------------------------------
|	DIMENSIONS    NEWTAXA        newtaxa
|	
|	              NTAX           ntax           GetNTax
|	
|	FORMAT        DATATYPE       datatype       GetDataType
|	
|	              RESPECTCASE    respectingCase IsRespectCase
|	
|	              MISSING        missing        GetMissingSymbol
|	
|	              SYMBOLS        symbols        GetSymbols
|	
|	              EQUATE         equates        GetEquateKey
|	                                            GetEquateValue
|	                                            GetNumEquates
|	
|	              (NO)LABELS     labels         IsLabels
|	
|	TAXLABELS                    taxonLabels    GetTaxonLabels

|	MATRIX                       matrix         GetState
|	                                            GetInternalRepresentation
|	                                            GetNumStates
|	                                            GetNumMatrixRows
|	                                            IsPolymorphic
|>
*/
class NxsUnalignedBlock
  : public NxsBlock
	{
	friend class NxsAssumptionsBlock;

	typedef std::vector<NxsDiscreteDatum>	UnalignedVect;
	typedef std::vector<UnalignedVect>		UnalignedMatrix;

	public:

		class NxsX_NoSuchCharacter {};	/* thrown if a function is called with an index to a character beyond the last character for a particular taxon */

		class NxsX_NoDataForTaxon
			{
			public:
				NxsX_NoDataForTaxon(unsigned i) : taxon_index(i) {}
				unsigned taxon_index;
			};	/* thrown if a function is called with an index to a taxon for which no data is stored */

		enum DataTypesEnum		/* values used to represent different basic types of data stored in an UNALIGNED block, and used with the data member `datatype' */
			{
			standard = 1,		/* indicates `matrix' holds characters with arbitrarily-assigned, discrete states, such as discrete morphological data */
			dna,				/* indicates `matrix' holds DNA sequences (states A, C, G, T) */
			rna,				/* indicates `matrix' holds RNA sequences (states A, C, G, U) */
			nucleotide,			/* indicates `matrix' holds nucleotide sequences */
			protein				/* indicates `matrix' holds amino acid sequences */
			};

								NxsUnalignedBlock(NxsTaxaBlock * tb, NxsAssumptionsBlock * ab);
		virtual					~NxsUnalignedBlock();

		//These are not needed and should be deleted as soon as the code for this block is finished
		////unsigned				ApplyExset(NxsUnsignedSet &exset);
		////unsigned				ApplyIncludeset(NxsUnsignedSet &inset);
		////unsigned				GetCharPos(unsigned origCharIndex);
		////unsigned				GetNChar();
		////unsigned				GetNCharTotal();
		////unsigned				GetNumActiveChar();
		////unsigned				GetNumEliminated();
		////unsigned				GetNumMatrixCols();
		////unsigned				GetOrigCharIndex(unsigned j) const;
		////unsigned				GetOrigCharNumber(unsigned j);
		////char					GetState(unsigned i, unsigned j, unsigned k = 0);
		////char					GetGapSymbol();
		////char					GetMatchcharSymbol();
		////bool					IsGapState(unsigned i, unsigned j);
		////bool					IsInterleave();
		////bool					IsTranspose();
		////bool					IsEliminated(unsigned origCharIndex);
		////void					ExcludeCharacter(unsigned i);
		////void					IncludeCharacter(unsigned i);
		////bool					IsActiveChar(unsigned j);
		////bool					IsExcluded(unsigned j);
		////bool					*GetActiveCharArray();
		////NxsString				GetCharLabel(unsigned i);
		////virtual unsigned		CharLabelToNumber(NxsString s);
		////void					WriteCharLabelsCommand(std::ostream &out) const;
		////void					WriteCharStateLabelsCommand(std::ostream &out) const;
		////void					WriteEliminateCommand(std::ostream &out) const;
		////virtual unsigned		GetMaxObsNumStates();
		////virtual unsigned		GetObsNumStates(unsigned j);
		////bool					IsTokens();
		////void					Consume(NxsUnalignedBlock &other);
		////NxsString				GetStateLabel(unsigned i, unsigned j);
		////virtual unsigned		TaxonLabelToNumber(NxsString s);

		void					ShowStateLabels(ostream & out, NxsDiscreteDatum s);
		unsigned				ApplyDelset(NxsUnsignedSet & delset);
		unsigned				ApplyRestoreset(NxsUnsignedSet & restoreset);
		unsigned				GetTaxPos(unsigned origTaxonIndex);
		unsigned				GetDataType();
		NxsIntVector			GetInternalRepresentation(unsigned i, unsigned j);
		unsigned				GetNTax();
		unsigned				GetNTaxTotal();
		unsigned				GetNumActiveTaxa();
		unsigned				GetNumEquates();
		unsigned				GetNumMatrixRows();
		unsigned				GetNumStates(unsigned i, unsigned j);
		unsigned				GetOrigTaxonIndex(unsigned j);
		unsigned				GetOrigTaxonNumber(unsigned j);
		unsigned				NumCharsForTaxon(unsigned i);
		char					GetMissingSymbol();
		bool					IsLabels();
		bool					IsMissingState(unsigned i, unsigned j);
		bool					IsPolymorphic(unsigned i, unsigned j);
		bool					IsRespectCase();
		void					DeleteTaxon(unsigned i);
		void					RestoreTaxon(unsigned i);
		bool					IsActiveTaxon(unsigned i);
		bool					IsDeleted(unsigned i);
		unsigned				GetStateSymbolIndex(unsigned i, unsigned j, unsigned k = 0);	// added by mth for standard data types
		char *					GetSymbols();
		bool *					GetActiveTaxonArray();
		NxsString				GetTaxonLabel(unsigned i);
		virtual void			DebugShowMatrix(ostream & out, char * marginText = NULL);
		virtual void			Report(ostream & out);
		virtual void			Reset();

		void					WriteAsNexus(std::ostream & out) const;
		void					WriteFormatCommand(std::ostream & out) const;
		void					WriteMatrixCommand(std::ostream & out) const;
		const char *			GetDatatypeName() const;

	protected:
	
		//These are not needed and should be deleted as soon as the code for this block is finished
		////std::vector<unsigned>	GetOrigIndexVector() const;
		////std::vector<unsigned>	GetOrigMatrixIndicesToWrite() const;
		////void 					WriteStatesForMatrixRow(std::ostream &out, unsigned taxon, unsigned first_taxon, const std::vector<unsigned> &origIndices) const;
		////void					BuildCharPosArray(bool check_eliminated = false);
		////void					HandleCharlabels(NxsToken &token);
		////void					HandleCharstatelabels(NxsToken &token);
		////void					HandleEliminate(NxsToken &token);
		////virtual void			HandleTransposedMatrix(NxsToken &token);
		////void					HandleStatelabels(NxsToken &token);
		////virtual void			HandleStdMatrix(NxsToken &token);
		////virtual unsigned		HandleTokenState(NxsToken & token, unsigned c);
		////void					WriteStates(NxsDiscreteDatum &d, char *s, unsigned slen);
		////void					ShowStates(ostream & out, unsigned i, unsigned j);
		//unsigned				nchar;				/* number of columns in matrix (same as `ncharTotal' unless some characters were eliminated, in which case `ncharTotal' > `nchar') */
		//unsigned				ncharTotal;			/* total number of characters (same as `nchar' unless some characters were eliminated, in which case `ncharTotal' > `nchar') */
		//bool					newchar;			/* true unless CHARLABELS or CHARSTATELABELS command read */
		//bool					transposing;		/* indicates matrix will be in transposed format */
		//bool					interleaving;		/* indicates matrix will be in interleaved format */
		//char					gap;				/* gap symbol for use with molecular data */
		//char					matchchar;			/* match symbol to use in matrix */
		//NxsUnsignedSet			eliminated;			/* array of (0-offset) character numbers that have been eliminated (will remain empty if no ELIMINATE command encountered) */
		//bool					*activeChar;		/* `activeChar[i]' true if character `i' not excluded; `i' is in range [0..`nchar') */
		//NxsStringVector			charLabels;			/* storage for character labels (if provided) */
		//NxsStringVectorMap		charStates;			/* storage for character state labels (if provided) */
		//mutable bool 			writeAllChars;
		//unsigned				*charPos;			/* maps character numbers in the data file to column numbers in matrix (necessary if some characters have been eliminated) */
		//NxsDiscreteMatrix *		matrix;			/* storage for discrete data */

		NxsStringMap			GetDefaultEquates() const;
		bool					IsInSymbols(char ch);
		void					HandleDimensions(NxsToken & token);
		void					HandleEndblock(NxsToken & token);
		virtual void			HandleFormat(NxsToken & token);
		virtual void			HandleMatrix(NxsToken & token);
		void					HandleMatrixImpl(NxsToken & token);
		virtual bool			HandleNextState(NxsToken & token, unsigned i, UnalignedVect & new_row);
		virtual void			Read(NxsToken & token);
		unsigned				PositionInSymbols(char ch);
		void					HandleTaxlabels(NxsToken & token);
		void					ResetSymbols();
		std::string				FormatState(NxsDiscreteDatum x) const;

		NxsTaxaBlock *			taxa;				/* pointer to the TAXA block in which taxon labels are stored */
		NxsAssumptionsBlock *	assumptionsBlock;	/* pointer to the ASSUMPTIONS block in which exsets, taxsets and charsets are stored */
		unsigned				ntax;				/* number of rows in umatrix (same as `ntaxTotal' unless fewer taxa appeared in UNALIGNED MATRIX command than were specified in the TAXA block, in which case `ntaxTotal' > `ntax') */
		unsigned				ntaxTotal;			/* number of taxa (same as `ntax' unless fewer taxa appeared in UNALIGNED MATRIX command than were specified in the TAXA block, in which case `ntaxTotal' > `ntax') */
		bool					newtaxa;			/* true if NEWTAXA keyword encountered in DIMENSIONS command */
		bool					respectingCase;		/* if true, RESPECTCASE keyword specified in FORMAT command */
		bool					tokens;				/* if false, data matrix entries must be single symbols; if true, multicharacter entries are allows */
		bool					labels;				/* indicates whether or not labels will appear on left side of matrix */
		char					missing;			/* missing data symbol */
		char *					symbols;			/* list of valid character state symbols */
		NxsStringMap			equates;			/* list of associations defined by EQUATE attribute of FORMAT command */
		unsigned *				taxonPos;			/* `taxonPos[i]' is the row number in `umatrix' corresponding to the taxon at index i in the TAXA block (necessary if fewer taxa appear in UNALIGNED block MATRIX command than are specified in the TAXA block). This array is `ntaxTotal' elements long, where `ntaxTotal' is the number of taxa in the TAXA block */
		bool *					activeTaxon;		/* `activeTaxon[i]' true if taxon `i' not deleted; `i' is in range [0..`ntax') */
		UnalignedMatrix			umatrix;			/* storage for unaligned discrete data */

	private:
		
		DataTypesEnum			datatype;			/* flag variable (see datatypes enum) */
	};

#include "nxsunalignedblock.inl"

#endif
