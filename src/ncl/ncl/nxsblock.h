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
#ifndef NCL_NXSBLOCK_H
#define NCL_NXSBLOCK_H

#include <vector>

class NxsReader;
class NxsBlock;
typedef std::vector<NxsBlock *> VecBlockPtr;
typedef std::vector<const NxsBlock *> VecConstBlockPtr;
typedef std::pair<const NxsBlock *, std::string> BlockUniqueID;
typedef std::map<BlockUniqueID, NxsBlock *> NxsBlockMapper;

std::string GetBlockIDTitleString(NxsBlock &);
/*----------------------------------------------------------------------------------------------------------------------
|	This is the base class from which all block classes are derived. A NxsBlock-derived class encapsulates a Nexus block
|	(e.g. DATA block, TREES block, etc.). The abstract virtual function Read must be overridden for each derived class 
|	to provide the ability to read everything following the block name (which is read by the NxsReader object) to the 
|	end or endblock statement. Derived classes must provide their own data storage and access functions. The abstract
|	virtual function Report must be overridden to provide some feedback to user on contents of block. The abstract
|	virtual function Reset must be overridden to empty the block of all its contents, restoring it to its 
|	just-constructed state.
*/
class NxsBlock
	{
	friend class NxsReader;

	public:
							NxsBlock();
		virtual				~NxsBlock();

		void				SetNexus(NxsReader *nxsptr);

		NxsString			GetID();
		bool				IsEmpty();

		void				Enable();
		void				Disable();
		bool				IsEnabled();
		bool				IsUserSupplied();

		virtual unsigned	CharLabelToNumber(NxsString s);
		virtual unsigned	TaxonLabelToNumber(NxsString s);

		virtual void		SkippingCommand(NxsString commandName);

		virtual void		Report(std::ostream &out);
		virtual void		Reset();

		NxsString			errormsg;			/* workspace for creating error messages */


		BlockUniqueID		GetInstanceIdentifier() const
			{
			return BlockUniqueID(this, GetInstanceName());
			}
			
		const std::string  &GetInstanceName() const
			{
			return title;
			}
			
		virtual NxsBlock			*CloneBlock(NxsBlockMapper &memo) const;
		virtual void				WriteAsNexus(std::ostream &out) const;
		virtual VecBlockPtr			GetImpliedBlocks();
		virtual VecConstBlockPtr	GetImpliedBlocksConst() const;
		virtual void 				WriteTitleCommand(std::ostream &out) const;


	protected:
		void				DemandEndSemicolon(NxsToken &token, const char *contextString);
		void				DemandEquals(NxsToken &token, const char *contextString);
		unsigned 			DemandPositiveInt(NxsToken &token, const char *contextString);
		void				GenerateNxsException(NxsToken &token, const char *message = NULL);
		void				GenerateUnexpectedTokenNxsException(NxsToken &token, const char *expected = NULL);
		bool				isEmpty;			/* true if this object is currently storing data */
		bool				isEnabled;			/* true if this block is currently ebabled */
		bool				isUserSupplied;		/* true if this object has been read from a file; false otherwise */
		NxsReader			*nexus;				/* pointer to the Nexus file reader object */
		NxsBlock			*next;				/* pointer to next block in list */
		NxsString			id;					/* holds name of block (e.g., "DATA", "TREES", etc.) */
		std::string			title;				/* holds the title of the block empty by default */
				
		virtual void		Read(NxsToken &token);
	};
	
	
inline std::string GetBlockIDTitleString(NxsBlock &b)
	{
	const std::string &t = b.GetInstanceName();
	std::string r = b.GetID();
	r.append(" block");
	if (t.length() > 0)
		{
		r.append(" (");
		r.append(t);
		r.append(")");
		}
	return r;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Advances the token, and returns the unsigned int that the token represents
|
| 	Sets errormsg and raises a NxsException on failure.
|	`contextString` is used in error messages:
|		"${contextString} must be a number greater than 0"
*/
inline unsigned NxsBlock::DemandPositiveInt(NxsToken &token, const char *contextString)
	{
	return NxsToken::DemandPositiveInt(token, this->errormsg, contextString);
	}

inline void NxsBlock::DemandEndSemicolon(NxsToken &token, const char *contextString)
	{
	NxsToken::DemandEndSemicolon(token, this->errormsg, contextString);
	}
#endif


