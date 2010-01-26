#include "ncl/ncl.h"

NxsTaxaBlock	*taxa	= NULL;
NxsTreesBlock	*trees	= NULL;

class MyReader : public NxsReader
	{
	public:
		ifstream inf;
		ofstream outf;

		MyReader(char *infname, char *outfname) : NxsReader()
			{
			inf.open(infname, ios::binary);
			outf.open(outfname);
			}

		~MyReader()
			{
			inf.close();
			outf.close();
			}

	void ExecuteStarting() {}
	void ExecuteStopping() {}

	bool EnteringBlock(NxsString blockName)
		{
		cout << "Reading \"" << blockName << "\" block..." << endl;
		outf << "Reading \"" << blockName << "\" block..." << endl;

		// Returning true means it is ok to delete any data associated with 
		// blocks of this type read in previously
		//
		return true;	
		}

	void SkippingBlock(NxsString blockName)
		{
		cout << "Skipping unknown block (" << blockName << ")..." << endl;
		outf << "Skipping unknown block (" << blockName << ")..." << endl;
		}

	void SkippingDisabledBlock(NxsString ) 
		{
		}

	void OutputComment(const NxsString &msg)
		{
		outf << msg;
		}

	void NexusError(NxsString msg, file_pos pos, long line, long col)
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

		exit(0);
		}
	};

class MyToken : public NxsToken
	{
	public:

		MyToken(istream &is, ostream &os)
			:NxsToken(is),
			out(os)
			{
			}

		void OutputComment(const NxsString &msg)
			{
			cout << msg << endl;
			out << msg << endl;
			}

	private:
		ostream &out;
	};

int main(int , char *argv[])
	{
	taxa = new NxsTaxaBlock();
	trees = new NxsTreesBlock(taxa);

	MyReader nexus(argv[1], argv[2]);
	nexus.Add(taxa);
	nexus.Add(trees);

	MyToken token(nexus.inf, nexus.outf);
	nexus.Execute(token);

	taxa->Report(nexus.outf);
	trees->Report(nexus.outf);

	return 0;
	}

