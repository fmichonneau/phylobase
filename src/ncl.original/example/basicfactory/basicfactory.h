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

#ifndef NCL_BASICFACTORY_H
#define NCL_BASICFACTORY_H

#define COMMAND_MAXLEN  255
#include "ncl/nxsassumptionsblock.h"
#include "ncl/nxscharactersblock.h"
#include "ncl/nxsdatablock.h"
#include "ncl/nxsdistancesblock.h"
#include "ncl/nxstaxablock.h"
#include "ncl/nxstreesblock.h"
#include "ncl/nxsunalignedblock.h"

/*----------------------------------------------------------------------------------------------------------------------
|	BASICFACTORY provides a template for creating a program that reads NEXUS data files using the factory API.
|	See BASICCMDLINE for a full description of the basic features, and compare this source to that BASICCMDLINE's
|	source to see how the factory API differs from NCL's classic, recycling API
*/
class BASICFACTORY
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

								BASICFACTORY();
		virtual					~BASICFACTORY();

		bool					EnteringBlock(NxsString blockName);
		void					ExitingBlock(NxsString blockName);
		void					ExecuteStarting();
		void					ExecuteStopping();
		void					OutputComment(const NxsString &msg);
		void					HandleNextCommand();
		void					NexusError(NxsString msg, file_pos pos, long line, long col);
		void					PreprocessNextCommand();
		void					PrintMessage(bool linefeed = true) NCL_COULD_BE_CONST ;
		virtual void			Report(ostream &out) NCL_COULD_BE_CONST ;
		void					Run(char *infile_name);
		void					SkippingBlock(NxsString blockName);
		void					SkippingCommand(NxsString commandName);
		void					SkippingDisabledBlock(NxsString blockName);
		virtual bool			UserQuery(NxsString mb_message, NxsString mb_title, BASICFACTORY::UserQueryEnum mb_choices = BASICFACTORY::uq_ok);

	protected:

		bool					inf_open;			/* true if `inf' is currently open */
		bool					logf_open;			/* true if `logf' is currently open */
		bool					quit_now;			/* set to false at beginning of Run and turns true only when QUIT command processed */
		mutable ofstream		logf;				/* the log file output stream */
		mutable NxsString		message;			/* workspace for composing output strings */
		char *					next_command;		/* workspace for processing next command entered interactively by user */

		unsigned				CharLabelToNumber(NxsString s) NCL_COULD_BE_CONST;
		bool					FileExists(const char* fn);
		NxsString				GetFileName(NxsToken& token);
		void					FactoryDefaults();
		void					HandleEndblock(NxsToken& token);
		void					HandleShow(NxsToken& token);
		void					HandleHelp(NxsToken& token);
		void					HandleLog(NxsToken& token);
		void					HandleExecute(NxsToken& token);
		void					PurgeBlocks();
		virtual void			Read(NxsToken& token);
		virtual void			Reset();
		unsigned				TaxonLabelToNumber(NxsString s) const;
		
		
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
inline void BASICFACTORY::ExecuteStarting()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Will be called by NxsReader::Execute just before it exits after reading to the end of a NEXUS data file (or until
|	encountering a LEAVE command between NEXUS blocks. Add code here if you need to clean up any memory allocated in
|	ExecuteStarting.
*/
inline void BASICFACTORY::ExecuteStopping()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called if an "output comment" is encountered in a NEXUS data file. An output comment is a comment [text enclosed in
|	square brackets] that begins with an exclamation point. [!This is an example of a NEXUS output comment]. Output
|	comments are supposed to be displayed when encountered. Modify this function's body to display output comments, 
|	which are made available as they are encountered via the `msg' argument.
*/
inline void	BASICFACTORY::OutputComment(const NxsString &)
	{
	}

#endif

