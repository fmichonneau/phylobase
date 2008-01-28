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

#include "ncl.h"
#include "ncltest.h"

/*----------------------------------------------------------------------------------------------------------------------
|	Takes an ostream reference in addition to the istream reference required by the base class NxsToken. The istream 
|	'is' is passed along to the base class constructor, and 'os' is used to initialize the data member out.
*/
MyToken::MyToken(
  istream &is,	/* reference to input file stream attached to the NEXUS data file to be read */
  ostream &os)	/* reference to the output stream used for outputting results */
  :NxsToken(is),
  out(os)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides the base class version. When a NEXUS printable comment is encountered in a NEXUS file, the comment is 
|	passed to the virtual function OutputComment in its argument 'msg'. This version outputs 'msg' directly to both the
|	standard output (cout) and the output file (out).
*/
void MyToken::OutputComment(
  const NxsString &msg)	/* reference to string containing the printable comment embedded in the NEXUS data file */
	{
	cout << msg << endl;
	out << msg << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes infile_exists to true if a file by the name 'infname' exists and to false otherwise. If a file named
|	'infname' exists, proceeds to open the ifstream data member inf to that file and the ofstream data member outf to
|	the file specified by the 'outfname' argument. Does not check for the existance of the output file beforehand, so
|	will overwrite a file if one by than name already exists. In a real application, use the FileExists member function
|	to check for the existance of the output file as well and warn the user before opening it.
*/
MyNexusFileReader::MyNexusFileReader(
  char *infname,	/* the name of the input file (should be a NEXUS-formatted data file) */
  char *outfname)	/* the name of a file to send output to */
 : NxsReader() 
	{
	infile_exists = true;
	if (!FileExists(infname))
		infile_exists = false;

	if (infile_exists)
		{
		inf.open(infname, ios::in | ios::binary);
		outf.open(outfname);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor simply closes the input file stream inf and the output file stream outf.
*/
MyNexusFileReader::~MyNexusFileReader() 
	{
	inf.close();
	outf.close();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. Simply reports to both standard output and the output file stream outf that
|	a block named 'blockName' was just entered by the NEXUS file reader. Always returns true.
*/
bool MyNexusFileReader::EnteringBlock(
  NxsString blockName)	/* the name of the block just entered */
	{
	cout << "Reading \"" << blockName << "\" block..." << endl;
	outf << "Reading \"" << blockName << "\" block..." << endl;
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. Simply reports to both standard output and the output file stream outf that
|	the NEXUS file reader just finished with a block named 'blockName'.
*/
void MyNexusFileReader::ExitingBlock(
  NxsString blockName)	/* the name of the block being exited */
	{
	cout << "Finished with \"" << blockName << "\" block." << endl;
	outf << "Finished with \"" << blockName << "\" block." << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. Simply outputs the string 's' to both the standard output and the output
|	stream outf. This function will be called when a printable NEXUS comment is found in a NEXUS file. The contents of
|	the comment are provided in the string argument 's'.
*/
void MyNexusFileReader::OutputComment(
  const NxsString &s)	/* the contents of the printable NEXUS comment */
	{
	cout << endl;
	cout << s << endl;
	outf << endl;
	outf << s << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	A simple function to test for the existance of a file named 'fn'. Returns true if an attempt to open the file 
|	succeeds, and false otherwise. No attempt is made to diagnose failure to open the file, so this function may return
|	false even when the file exists if the file could not be opened for some other reason (e.g. user does not have
|	read permissions for that file).
*/
bool MyNexusFileReader::FileExists(
  char *fn)	/* the name of the file to test (path should be included if you do not want to assume the current working directory) */
	{
	if (fn == NULL) 
		return false;

	bool exists = false;

	FILE *f = fopen(fn, "r");
	if (f != NULL)
		{
		exists = true;
		fclose(f);
		f = NULL;
		}

	return exists;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. Simply reports to both the standard output and to the output stream outf
|	that a block named 'blockName' is being skipped because this NEXUS file reader does not know how to interpret blocks
|	by this name. This situation is common if the NEXUS file has private blocks, e.g. a PAUP block containing commands
|	understood only by the program PAUP.
*/
void MyNexusFileReader::SkippingBlock(
  NxsString blockName)	/* the name of the block being skipped */
	{
	cout << "Skipping unknown block (" << blockName << ")..." << endl;
	outf << "Skipping unknown block (" << blockName << ")..." << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. Simply reports to both the standard output and to the output stream outf
|	that a block named 'blockName' is being skipped because this block, though known to the NEXUS reader, is currently
|	disabled. This situation arises when reading a treefile for instance (the application may disable all blocks except
|	NxsTreesBlock because only TREES blocks expected in a treefile. Overriding this function gives the application a
|	chance to report to the user that NEXUS blocks were found in the file that were unexpected under the circumstances.
*/
void MyNexusFileReader::SkippingDisabledBlock(
  NxsString blockName)	/* the name of the disabled block being skipped */
	{
	cout << "Skipping disabled block (" << blockName << ")..." << endl;
	outf << "Skipping disabled block (" << blockName << ")..." << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. If a block of the type referenced by 'nexusBlock' has been read (i.e. the
|	IsEmpty member function returns false), this function calls the Report member function of 'nexusBlock' to dump the
|	contents of the block object to the output file stream outf.
*/
void MyNexusFileReader::DebugReportBlock(
  NxsBlock &nexusBlock)	/* reference to the NxsBlock to be dumped if not empty */
	{
	if (!nexusBlock.IsEmpty())
		{
		outf << endl;
		outf << "********** Contents of the " << nexusBlock.GetID() << " block **********" << endl;
		nexusBlock.Report(outf);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. This function simply reports the error message contained in the string
|	'msg' and the file position, line and column (contained in 'pos', 'line' and 'col', respectively) to both the 
|	standard error stream and the output file stream outf. It then aborts the application by invoking exit(0).
*/
void MyNexusFileReader::NexusError(
  NxsString msg,	/* reference to the error message string */
  file_pos pos,		/* the file position where the error was discovered */
  unsigned line,		/* the line where the error was discovered */
  unsigned col)			/* the column where the error was discovered */
	{
	cerr << endl;
	cerr << "Error found at line " << line;
	cerr << ", column " << col;
	cerr << " (file position " << pos << "):" << endl;
	cerr << msg << endl;

	outf << endl;
	outf << "Error found at line " << line;
	outf << ", column " << col;
	outf << " (file position " << pos << "):" << endl;
	outf << msg << endl;

	cerr << endl;
	cerr << "Press return to quit..." << endl;
	cin.get();

	exit(0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Class derived from NxsCharactersBlock so that an output file stream can be used for reporting unknown commands that
|	are being skipped. Overrides SkippingCommand function for this purpose.
*/
MyCharactersBlock::MyCharactersBlock(
  NxsTaxaBlock *tb,			/* pointer to the NxsTaxaBlock object to use to store taxon labels */
  NxsAssumptionsBlock *ab,	/* pointer to the NxsAssumptionsBlock object to use for storing charsets, taxsets, etc. */
  ostream &o)				/* reference to the ostream object for use in reporting commands that are being skipped */
  :NxsCharactersBlock(tb, ab),
  outf(o)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. Simply reports the command being skipped (name of command is in 's') both 
|	to the standard output stream as well as to the output stream outf.
*/
void MyCharactersBlock::SkippingCommand(
  NxsString s)	/* the name of the command being skipped */
	{
	cout << "Skipping unknown command (" << s << ")..." << endl;
	outf << "Skipping unknown command (" << s << ")..." << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Class derived from NxsDataBlock so that an output file stream can be used for reporting unknown commands that are 
|	being skipped. Overrides SkippingCommand function for this purpose.
*/
MyDataBlock::MyDataBlock(
  NxsTaxaBlock *tb,			/* pointer to the NxsTaxaBlock object to use to store taxon labels */
  NxsAssumptionsBlock *ab,	/* pointer to the NxsAssumptionsBlock object to use for storing charsets, taxsets, etc. */
  ostream &o)				/* reference to the ostream object for use in reporting commands that are being skipped */
  :NxsDataBlock(tb, ab),
  outf(o)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. Simply reports the command being skipped (name of command is in 's') both 
|	to the standard output stream as well as to the output stream outf.
*/
void MyDataBlock::SkippingCommand(
  NxsString s)	/* the name of the command being skipped */
	{
	cout << "Skipping unknown command (" << s << ")..." << endl;
	outf << "Skipping unknown command (" << s << ")..." << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Class derived from NxsAssumptionsBlock so that an output file stream can be used for reporting unknown commands 
|	that are being skipped. Overrides SkippingCommand function for this purpose.
*/
MyAssumptionsBlock::MyAssumptionsBlock(
  NxsTaxaBlock *tb,			/* pointer to the NxsTaxaBlock object to use to store taxon labels */
  ostream &o)				/* reference to the ostream object for use in reporting commands that are being skipped */
  :NxsAssumptionsBlock(tb),
   outf(o)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class virtual function. Simply reports the command being skipped (name of command is in 's') both 
|	to the standard output stream as well as to the output stream outf.
*/
void MyAssumptionsBlock::SkippingCommand(
  NxsString s)	/* the name of the command being skipped */
	{
	cout << "Skipping unknown command (" << s << ")..." << endl;
	outf << "Skipping unknown command (" << s << ")..." << endl;
	}

int main(int argc, char *argv[])
	{
	char infname[256];
	char outfname[256];
	strcpy(infname, "characters.nex");
	strcpy(outfname, "output.txt");
	if (argc > 1)
		strcpy(infname, argv[1]);
	if (argc > 2)
		strcpy(outfname, argv[2]);

	MyNexusFileReader nexus(infname, outfname);
	if (!nexus.infile_exists)
		{
		ofstream outf(outfname);
		outf << "Error: specified input file (";
		outf << infname;
		outf << ") does not exist." << endl;
		outf.close();

		cerr << "Error: specified input file (";
		cerr << infname;
		cerr << ") does not exist." << endl;

		cerr << endl;
		cerr << "Press return to quit..." << endl;
		cin.get();

		return 0;
		}

	NxsTaxaBlock		*taxa			= new NxsTaxaBlock();
	MyAssumptionsBlock	*assumptions	= new MyAssumptionsBlock(taxa, nexus.outf);
	NxsTreesBlock		*trees			= new NxsTreesBlock(taxa);
	MyCharactersBlock	*characters		= new MyCharactersBlock(taxa, assumptions, nexus.outf);
	MyDataBlock			*data			= new MyDataBlock(taxa, assumptions, nexus.outf);
	NxsDistancesBlock	*distances		= new NxsDistancesBlock(taxa);

	nexus.Add(taxa);
	nexus.Add(assumptions);
	nexus.Add(trees);
	nexus.Add(characters);
	nexus.Add(data);
	nexus.Add(distances);

	MyToken token(nexus.inf, nexus.outf);
	nexus.Execute(token);

	nexus.outf << endl;
	nexus.outf << "Reports follow on all blocks currently in memory." << endl;

	if (!taxa->IsEmpty())
		taxa->Report(nexus.outf);
	if (!trees->IsEmpty())
		trees->Report(nexus.outf);
	if (!characters->IsEmpty())
		characters->Report(nexus.outf);
	if (!data->IsEmpty())
		data->Report(nexus.outf);
	if (!distances->IsEmpty())
		distances->Report(nexus.outf);
	if (!assumptions->IsEmpty())
		assumptions->Report(nexus.outf);

	NxsStringVector z;
	nexus.outf << endl;
	if (assumptions->GetNumTaxSets() > 0)
		{
		assumptions->GetTaxSetNames(z);
		NxsStringVector::const_iterator zi;
		for (zi = z.begin(); zi != z.end(); zi++)
			{
			nexus.outf << "Taxset " << (*zi) << ": ";
			NxsUnsignedSet &x = assumptions->GetTaxSet(*zi);
			NxsUnsignedSet::const_iterator xi;
			for (xi = x.begin(); xi != x.end(); xi++)
				nexus.outf << (*xi + 1) << " ";
			nexus.outf << endl;
			}
		}

	if (assumptions->GetNumCharSets() > 0)
		{
		assumptions->GetCharSetNames(z);
		NxsStringVector::const_iterator zi;
		for (zi = z.begin(); zi != z.end(); zi++)
			{
			nexus.outf << "Charset " << (*zi) << ": ";
			NxsUnsignedSet &x = assumptions->GetCharSet(*zi);
			NxsUnsignedSet::const_iterator xi;
			for (xi = x.begin(); xi != x.end(); xi++)
				nexus.outf << (*xi + 1) << " ";
			nexus.outf << endl;
			}
		}

	if (assumptions->GetNumExSets() > 0)
		{
		assumptions->GetExSetNames(z);
		NxsStringVector::const_iterator zi;
		for (zi = z.begin(); zi != z.end(); zi++)
			{
			nexus.outf << "Exset " << (*zi) << ": ";
			NxsUnsignedSet& x = assumptions->GetExSet(*zi);
			NxsUnsignedSet::const_iterator xi;
			for (xi = x.begin(); xi != x.end(); xi++)
				nexus.outf << (*xi + 1) << " ";
			nexus.outf << endl;
			}
		}

	//cerr << endl;
	//cerr << "Press return to quit..." << endl;
	//cin.get();

	return 0;
	}
