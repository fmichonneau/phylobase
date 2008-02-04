//  Copyright (C) 2007-2008 Brian O'Meara & Derrick Zwickl
//  A modification of the BasicCommandLine file of the NCL (see below)
//  to use for loading trees and data from Nexus into R. Licensing as below.

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

#ifndef NCL_BASICCMDLINE_H
#define NCL_BASICCMDLINE_H

#define COMMAND_MAXLEN  255

/*----------------------------------------------------------------------------------------------------------------------
|	BASICCMDLINE provides a template for creating a program that reads NEXUS data files and provides a basic command 
|	line. After compiling BASICCMDLINE, you will already have a program that understands the following commands, either 
|	typed in at the console or provided in a BASICCMDLINE block in a NEXUS data file (exception is the execute command,
|	which can only be entered at the console). Keywords in the descriptions below are given in uppercase, however the
|	commands themselves are case-insensitive. Lower-case indicates a parameter supplied by the user (e.g., "filename" 
|	would be replaced by the actual name of the file). Square brackets indicate optional keywords or subcommands.
|>
|	EXECUTE filename;
|	
|	LOG [options];
|	
|	  Option         Action
|	  ------------------------------------------------------
|	  FILE=filename  specifies name of log file to start
|	  START          indicates logging is to be started
|	  STOP           indicates logging is to be stopped
|	  APPEND         append to log file if it already exists
|	  REPLACE        replace log file without asking
|	
|	QUIT;
|>
|	See the Read function for details and to add other commands.
|	
|	To change the name of the program (which is also the prompt name and the name of the program's private NEXUS 
|	block), replace all occurrences of BASICCMDLINE with the name of your program (also search for the string 
|	"basiccmdline" and replace with an appropriate string at each occurrence).
|	
|	This class handles reading and storage for the NxsReader block BASICCMDLINE. It also serves as the main class for 
|	the program BASICCMDLINE, acting as both a NxsReader object (in order to be capable of parsing data files) as well 
|	as a NxsBlock object (in order to be able to process commands in a BASICCMDLINE block). 
|
|	Acting as a NxsBlock, it overrides the member functions Read and Reset, which are virtual functions in the base 
|	class NxsBlock. Acting as a NxsReader object, it overrides the member functions EnteringBlock, SkippingBlock, and 
|	NexusError.
|	
|	Adding a new data member? Don't forget to:
|~
|	o Describe it in the class header comment at the top of "basiccmdline.h"
|	o Initialize it (unless it is self-initializing) in the constructor and reinitialize it in the Reset function
|	o Describe the initial state in the constructor documentation
|	o Delete memory allocated to it in both the destructor and Reset function
|	o Report it in some way in the Report function
|~
*/
class BASICCMDLINE
  : public NxsBlock,
  public NxsReader
	{
	public:

		enum UserQueryEnum		/* enumeration used with UserQuery member function to specify which choices to provide the user */
			{
			uq_cancel = 0x01,	/* provide opportunity to cancel */
			uq_ok	  = 0x02,	/* provide opportunity to answer ok */
			uq_yes	  = 0x04,	/* provide opportunity to answer yes */
			uq_no 	  = 0x08	/* provide opportunity to answer no */
			};

							BASICCMDLINE();
		virtual				~BASICCMDLINE();

		bool				EnteringBlock(NxsString blockName);
		void				ExitingBlock(NxsString blockName);
		void				ExecuteStarting();
		void				ExecuteStopping();
		void				OutputComment(const NxsString &msg);
		void				HandleNextCommand();
		void				NexusError(NxsString msg, file_pos pos, long line, long col);
		void				PreprocessNextCommand();
		void				PrintMessage(bool linefeed = true);
		virtual void		Report(ostream &out);
		void				Run(char *infile_name);
		void				Initialize(char *infile_name);
		void				SkippingBlock(NxsString blockName);
		void				SkippingCommand(NxsString commandName);
		void				SkippingDisabledBlock(NxsString blockName);
		virtual bool		UserQuery(NxsString mb_message, NxsString mb_title, BASICCMDLINE::UserQueryEnum mb_choices = BASICCMDLINE::uq_ok);
		void				HandleReturnData(NxsToken& token);
		NxsString			ReturnDataForR(bool allchar, bool polymorphictomissing, bool levelsall);
		NxsString			RemoveUnderscoresAndSpaces(NxsString input);
		void				RReturnCharacters(NxsString & nexuscharacters, bool allchar, bool polymorphictomissing, bool levelsall);
		void				RReturnTrees(NxsString & nexustrees);
		void				RReturnDistances(NxsString & nexusdistances);
		string				TestRunning();
		
	protected:

		bool				inf_open;			/* true if `inf' is currently open */
		bool				logf_open;			/* true if `logf' is currently open */
		bool				quit_now;			/* set to false at beginning of Run and turns true only when QUIT command processed */
		ofstream			logf;				/* the log file output stream */
		NxsString			message;			/* workspace for composing output strings */
		NxsTreesBlock		*trees;				/* pointer to NxsTreesBlock object */
		NxsTaxaBlock		*taxa;				/* pointer to NxsTaxaBlock object */
		NxsAssumptionsBlock	*assumptions;		/* pointer to NxsAssumptionsBlock object */
		NxsDistancesBlock	*distances;			/* pointer to NxsDistancesBlock object */
		NxsCharactersBlock	*characters;		/* pointer to NxsCharactersBlock object */
		NxsDataBlock		*data;				/* pointer to NxsDataBlock object */
		char				*next_command;		/* workspace for processing next command entered interactively by user */
		unsigned			CharLabelToNumber(NxsString s);
		bool				FileExists(const char* fn);
		NxsString			GetFileName(NxsToken& token);
		void				FactoryDefaults();
		void				HandleEndblock(NxsToken& token);
		void				HandleShow(NxsToken& token);
		void				HandleHelp(NxsToken& token);
		void				HandleLog(NxsToken& token);
		void				HandleExecute(NxsToken& token);
		void				PurgeBlocks();
		virtual void		Read(NxsToken& token);
		virtual void		Reset();
		unsigned			TaxonLabelToNumber(NxsString s);
	};

/*----------------------------------------------------------------------------------------------------------------------
|	The MyNexusToken class provides a NxsToken-derived object that can display output comments as it encounters them.
|	The virtual function NxsToken::OutputComment is overridden in this class for this purpose.
*/
class MyNexusToken
  : public NxsToken
	{
	public:
				MyNexusToken(istream &i);

		void	OutputComment(const NxsString &msg);
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Will be called by NxsReader::Execute after the initial "#NEXUS" keyword is found in a NEXUS file but before other
|	tokens are read. Add code here if you need to do any initializations prior to encountering any NEXUS blocks in a
|	NEXUS data file.
*/
inline void BASICCMDLINE::ExecuteStarting()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Will be called by NxsReader::Execute just before it exits after reading to the end of a NEXUS data file (or until
|	encountering a LEAVE command between NEXUS blocks. Add code here if you need to clean up any memory allocated in
|	ExecuteStarting.
*/
inline void BASICCMDLINE::ExecuteStopping()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called if an "output comment" is encountered in a NEXUS data file. An output comment is a comment [text enclosed in
|	square brackets] that begins with an exclamation point. [!This is an example of a NEXUS output comment]. Output
|	comments are supposed to be displayed when encountered. Modify this function's body to display output comments, 
|	which are made available as they are encountered via the `msg' argument.
*/
inline void	BASICCMDLINE::OutputComment(const NxsString &)
	{
	}



#endif

