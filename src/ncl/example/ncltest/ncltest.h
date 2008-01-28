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

/*----------------------------------------------------------------------------------------------------------------------
|	MyToken is derived from NxsToken in order to provide a way to report printable NEXUS comments (the square-bracket
|	comments in NEXUS files that begin with an exclamation point and are intended to be printed in the outout).
*/
class MyToken
  : public NxsToken
	{
	public:
		MyToken(istream &is, ostream &os);
		void OutputComment(const NxsString &msg);

	private:
		ostream& out;	/* output file stream */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	This class represents a NEXUS file reader that supervises the interpretation of a NEXUS-formatted data file. It 
|	overloads several important virtual base class functions, including EnteringBlock and ExitingBlock (which allow it
|	to inform the user when particular blocks are being entered or exited), ExecuteStarting and ExecuteStopping (which
|	allow any initialization before a NEXUS file is read or cleanup afterwards, respectively), SkippingBlock and 
|	SkippingDisabledBlock (which allow user to be informed when unknown NEXUS blocks are encountered in the data file
|	or when NEXUS blocks are encountered that are known to this application but are being skipped because these blocks
|	were temporarily disabled), OutputComment (which determines how user-specified printable comments in NEXUS files are
|	displayed), and NexusError (which determines how errors encountered in the NEXUS data file are to be reported).
*/
class MyNexusFileReader
  : public NxsReader
	{
	public:
					MyNexusFileReader(char *infname, char *outfname);
					~MyNexusFileReader();

		bool		EnteringBlock(NxsString blockName);
		void		ExitingBlock(NxsString blockName);
		
		void		SkippingBlock(NxsString blockName);
		void		SkippingDisabledBlock(NxsString blockName);

		bool		FileExists(char *fn);
		void		OutputComment(const NxsString &s);
		void		DebugReportBlock(NxsBlock &nexusBlock);
		void		NexusError(NxsString msg, file_pos pos, unsigned line, unsigned col);

		bool		infile_exists;	/* true if input file is determined to exist, false otherwise */
		ifstream	inf;			/* input file stream used for reading NEXUS file */
		ofstream	outf;			/* output file stream used for reporting progress */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	This derived version of NxsCharactersBlock is necessary in order to provide for the use of an ostream to give 
|	feedback to the user and report on the information contained in any CHARACTERS blocks found. It overloads the base 
|	class virtual function SkippingCommand.
*/
class MyCharactersBlock
  : public NxsCharactersBlock
	{
	public:
				MyCharactersBlock(NxsTaxaBlock *tb, NxsAssumptionsBlock *ab, ostream &o);
		void	SkippingCommand(NxsString s);

	private:
		ostream &outf;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	This derived version of NxsDataBlock is necessary in order to provide for the use of an ostream to give 
|	feedback to the user and report on the information contained in any DATA blocks found. It overloads the base class
|	virtual function SkippingCommand.
*/
class MyDataBlock
  : public NxsDataBlock
	{
	public:
				MyDataBlock(NxsTaxaBlock *tb, NxsAssumptionsBlock *ab, ostream &o);
		void	SkippingCommand(NxsString s);

	private:
		ostream &outf;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	This derived version of NxsAssumptionsBlock was created in order to provide feedback to the user on commands that
|	are not yet impemented (much of the NxsAssumptionsBlock is not yet implemented, and the NxsAssumptionsBlock version
|	of SkippingCommand is simply the one it inherited from NxsBlock, which does nothing).
*/
class MyAssumptionsBlock
  : public NxsAssumptionsBlock
	{
	public:
				MyAssumptionsBlock(NxsTaxaBlock *tb, ostream &o);
		void	SkippingCommand(NxsString s);

	private:
		ostream &outf;
	};

