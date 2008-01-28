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
|	Initializes all pointer data members to NULL, and all bool data members to true except isUserSupplied, which is
|	initialized to false.
*/
NxsBlock::NxsBlock()
	{
	next			= NULL;
	nexus			= NULL;
	isEmpty			= true;
	isEnabled		= true;
	isUserSupplied	= false;

	id.clear();
	errormsg.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Advances the token, and raise an exception if it is not an equals sign.
|
| 	Sets errormsg and raises a NxsException on failure.
|	`contextString` is used in error messages:
|		"Expecting '=' ${contextString} but found..."
*/
void NxsBlock::DemandEquals(NxsToken &token, const char *contextString)
	{
	token.GetNextToken();
	if (!token.Equals("="))
		{
		errormsg = "Expecting '=' ";
		if (contextString)
			errormsg.append(contextString);
		errormsg += " but found ";
		errormsg += token.GetToken();
		errormsg += " instead";
		throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
		}
	}

void NxsBlock::GenerateNxsException(NxsToken &token, const char *message)
	{
	if (message)
		errormsg = message;
	throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}

void NxsBlock::GenerateUnexpectedTokenNxsException(NxsToken &token, const char *expected)
	{
	errormsg = "Unexpected token";
	if (expected)
		{
		errormsg += ". Expecting ";
		errormsg += expected;
		errormsg += ", but found: ";
		}
	else
		{
		errormsg += ": ";
		}
	errormsg += token.GetToken();
	throw NxsException(errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn());
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to be done.
*/
NxsBlock::~NxsBlock()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version simply returns 0 but a derived class should override this function if it needs to construct
|	and run a NxsSetReader object to read a set involving characters. The NxsSetReader object may need to use this 
|	function to look up a character label encountered in the set. A class that overrides this method should return the
|	character index in the range [1..nchar].
*/
unsigned NxsBlock::CharLabelToNumber(
  NxsString)	/* the character label to be translated to the character's number */
	{
	return 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of isEnabled to false. A NxsBlock can be disabled (by calling this method) if blocks of that type
|	are to be skipped during execution of the NEXUS file. If a disabled block is encountered, the virtual
|	NxsReader::SkippingDisabledBlock function is called, giving your application the opportunity to inform the user
|	that a block was skipped.
*/
void NxsBlock::Disable()
	{
	isEnabled = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of isEnabled to true. A NxsBlock can be disabled (by calling Disable) if blocks of that type are to
|	be skipped during execution of the NEXUS file. If a disabled block is encountered, the virtual 
|	NxsReader::SkippingDisabledBlock function is called, giving your application the opportunity to inform the user
|	that a block was skipped.
*/
void NxsBlock::Enable()
	{
	isEnabled = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of isEnabled, which can be controlled through use of the Enable and Disable member functions. A 
|	NxsBlock should be disabled if blocks of that type are to be skipped during execution of the NEXUS file. If a 
|	disabled block is encountered, the virtual NxsReader::SkippingDisabledBlock function is called, giving your 
|	application the opportunity to inform the user that a block was skipped.
*/
bool NxsBlock::IsEnabled()
	{
	return isEnabled;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of isUserSupplied, which is true if and only if this block's Read function is called to process a 
|	block of this type appearing in a data file. This is useful because in some cases, a block object may be created 
|	internally (e.g. a NxsTaxaBlock may be populated using taxon names provided in a DATA block), and such blocks do 
|	not require permission from the user to delete data stored therein.
*/
bool NxsBlock::IsUserSupplied()
	{
	return isUserSupplied;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if Read function has not been called since the last Reset. This base class version simply returns the 
|	value of the data member isEmpty. If you derive a new block class from NxsBlock, be sure to set isEmpty to true in 
|	your Reset function and isEmpty to false in your Read function.
*/
bool NxsBlock::IsEmpty()
	{
	return isEmpty;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the id NxsString.
*/
NxsString NxsBlock::GetID()
	{
	return id;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This virtual function must be overridden for each derived class to provide the ability to read everything following
|	the block name (which is read by the NxsReader object) to the end or endblock statement. Characters are read from 
|	the input stream 'in'. Note that to get output comments displayed, you must derive a class from NxsToken, override 
|	the member function OutputComment to display a supplied comment, and then pass a reference to an object of the 
|	derived class to this function.
*/
void NxsBlock::Read(
  NxsToken &)	/* the NxsToken to use for reading block */
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This virtual function should be overridden for each derived class to completely reset the block object in 
|	preparation for reading in another block of this type. This function is called by the NxsReader object just prior to
|	calling the block object's Read function.
*/
void NxsBlock::Reset()
	{
	title = std::string();
	// Reset base class data members that could have changed
	//
	errormsg.clear();
	isEnabled      = true;
	isEmpty        = true;
	isUserSupplied = false;


	}

/*----------------------------------------------------------------------------------------------------------------------
|	This virtual function provides a brief report of the contents of the block.
*/
void NxsBlock::Report(
  ostream &)	/* the output stream to which the report is sent */
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the nexus data member of the NxsBlock object to 'nxsptr'.
*/
void NxsBlock::SetNexus(
  NxsReader *nxsptr)	/* pointer to a NxsReader object */
	{
	nexus = nxsptr;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	This function is called when an unknown command named commandName is about to be skipped. This version of the 
|	function does nothing (i.e., no warning is issued that a command was unrecognized). Override this virtual function 
|	in a derived class to provide such warnings to the user.
*/
void NxsBlock::SkippingCommand(
  NxsString )	/* the name of the command being skipped */
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version simply returns 0, but a derived class should override this function if it needs to construct
|	and run a NxsSetReader object to read a set involving taxa. The NxsSetReader object may need to use this function to
|	look up a taxon label encountered in the set. A class that overrides this method should return the taxon index in
|	the range [1..ntax].
*/
unsigned NxsBlock::TaxonLabelToNumber(
  NxsString )	/* the taxon label to be translated to a taxon number */
	{
	return 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of Blocks that were created by the reading in of this block (the prototypical case is the taxa block
|	that is implied by a data block).
*/
VecBlockPtr	NxsBlock::GetImpliedBlocks()
	{
	return VecBlockPtr();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of Blocks that were created by the reading in of this block (the prototypical case is the taxa block
|	that is implied by a data block).
*/
VecConstBlockPtr NxsBlock::GetImpliedBlocksConst() const
	{
	return VecConstBlockPtr();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Should returns a  new instance (deep copy) of the same type of block with the same state.
|	Note the NxsReader field should not be cloned (it should be aliased).
|
|	NxsBlock version throws NxsUnimplementedException (in future versions of NCL this will be a pure virtual.
|
|	NxsBlocks are expected to clone their linked blocks, but memo is passed in to avoid double cloning of shared references.
|	memo is an mapper of an old block to a new instance (used when groups of blocks are being cloned).
*/
NxsBlock * NxsBlock::CloneBlock(
  NxsBlockMapper & /// memo is an mapper of an old block to a new instance (used when groups of blocks are being cloned)
  ) const
	{
	throw NxsUnimplementedException();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Should serialize the content of the block as NEXUS.
|	NxsBlock version throws NxsUnimplementedException (in future versions of NCL this will be a pure virtual.
*/
void NxsBlock::WriteAsNexus(std::ostream &) const
	{
	throw NxsUnimplementedException();
	}

void NxsBlock::WriteTitleCommand(std::ostream &out) const
	{
	const std::string &t = this->GetInstanceName();
	if (t.length() > 0)
		out << "TITLE " << NxsToken::EscapeString(t) << ';';
	}


