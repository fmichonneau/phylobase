// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-

//  Copyright (C) 2007-2008 Brian O'Meara & Derrick Zwickl

#include <Rcpp.h>

#include <ncl/ncl.h>
#include "NCLInterface.h"

using namespace std;

// //This function receives a list of parameters from R, which can then be extracted below
// RcppExport SEXP ReadWithNCL(SEXP params) {

//     SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
//     char* exceptionMesg=NULL;
//     try {

// 		// Get parameters in params - only 1 is gotten now
// 		RcppParams rparam(params);
// #		if defined(FILENAME_AS_NEXUS)
// 			string filename = "'";
// 			filename+=rparam.getStringValue("filename");
// 			filename+="'";
// #		else
// 			string filename = rparam.getStringValue("filename");
// #		endif

// 		BASICCMDLINE reader;

// 		//this is where the reader would be passed the filename to read
// 		//reader.Run(filename.c_str()); //Will not compile
// 		//reader.Run(NULL);
// 		reader.Initialize(const_cast < char* > (filename.c_str()));

// 		string filenameString = "I was told to read a file named ";
// 		filenameString += filename;

// 		//Various calls to the reader can be made here to fill the various strings
// 		//string treeString = "This was passed in the tree string";
// 		NxsString treeStringNxs;
// 		reader.RReturnTrees(treeStringNxs);
// 		string treeString=treeStringNxs.c_str();

// 		string otherString = "";

// 		//string discreteString = "This was passed in the discrete string";

// 		NxsString characterStringNxs;
// 		reader.RReturnCharacters(characterStringNxs,false, true, false, false);
// 		//string characterString=characterStringNxs.c_str();

// 		string continuousString = "This was passed in the continuous string";
// 		//reader.GetContinuousString();

// 		//Just to test to see if reader is running
// 		string testString = reader.TestRunning();

// 		// Build result set to be returned as a list to R.
// 		RcppResultSet rs;

// 		//if the various strings are nonempty, add them to the ResultSet
// 		if(filenameString.length() > 0)
// 			rs.add("filenamestring", filenameString);
// 		if(treeString.length() > 0)
// 			rs.add("treestring", treeString);
// 		if(otherString.length() > 0)
// 			rs.add("otherstring", otherString);
// 		if(characterStringNxs.length() > 0)
// 			rs.add("characterstring", characterStringNxs.c_str());
// 		if(continuousString.length() > 0)
// 			rs.add("continuousstring", continuousString);
// 		if(testString.length() > 0)
// 			rs.add("teststring", testString);


// 		// Get the list to be returned to R.
// 		rl = rs.getReturnList();

//     	} catch(std::exception& ex) {
// 			exceptionMesg = copyMessageToR(ex.what());
//     	} catch(...) {
// 			exceptionMesg = copyMessageToR("unknown reason");
//     	}

//     if(exceptionMesg != NULL)
// 	Rf_error(exceptionMesg);

//     return rl;

// }

RcppExport SEXP ReadTreesWithNCL(SEXP params) {

    SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
    char* exceptionMesg=NULL;
    try {


		// Get parameters in params - only 1 is gotten now
		RcppParams rparam(params);
#		if defined(FILENAME_AS_NEXUS)
			string filename = "'";
			filename+=rparam.getStringValue("filename");
			filename+="'";
#		else
			string filename = rparam.getStringValue("filename");
#		endif

		BASICCMDLINE reader;

		//this is where the reader would be passed the filename to read
		//reader.Run(filename.c_str()); //Will not compile
		//reader.Run(NULL);
		reader.Initialize(const_cast < char* > (filename.c_str()));

		NxsString treeStringNxs;
		reader.RReturnTrees(treeStringNxs);
		string treeString=treeStringNxs.c_str();

		// Build result set to be returned as a list to R.
		RcppResultSet rs;

		//if the various strings are nonempty, add them to the ResultSet
		if(treeString.length() > 0)
			rs.add("treestring", treeString);


		// Get the list to be returned to R.
		rl = rs.getReturnList();

	} catch(std::exception& ex) {
		exceptionMesg = copyMessageToR(ex.what());
	} catch(...) {
		exceptionMesg = copyMessageToR("unknown reason");
	}

    if(exceptionMesg != NULL)
		Rf_error(exceptionMesg);

    return rl;

}

// old API version -- version for the new API below
// RcppExport SEXP ReadCharsWithNCL(SEXP params) {

//     SEXP  rl = R_NilValue;
//     char* exceptionMesg = NULL;
//     try {

// 	// Get parameters in params - only 1 is gotten now
// 	RcppParams rparam(params);
// #	if defined(FILENAME_AS_NEXUS)
// 	    string filename = "'";
// 	    filename+=rparam.getStringValue("filename");
// 	    filename+="'";
// #	else
// 	    string filename = rparam.getStringValue("filename");
// #	endif

// 	bool allchar = rparam.getBoolValue("allchar");
// 	bool levelsall=rparam.getBoolValue("levelsall");
// 	bool polymorphictomissing=rparam.getBoolValue("polymorphictomissing");
// 	bool returnlabels=rparam.getBoolValue("returnlabels");

// 	BASICCMDLINE reader;

// 	//this is where the reader would be passed the filename to read
// 	//reader.Run(filename.c_str()); //Will not compile
// 	//reader.Run(NULL);
// 	reader.Initialize(const_cast < char* > (filename.c_str()));

// 	NxsString charStringNxs;
// 	reader.RReturnCharacters(charStringNxs,allchar, polymorphictomissing, levelsall, returnlabels);
// 	string charString = charStringNxs.c_str();
// 	std::cout << "charString is: " << charString << endl;

// 	// Build result set to be returned as a list to R.
// 	RcppResultSet rs;

// 	//if the various strings are nonempty, add them to the ResultSet
// 	if(charString.length() > 0)
// 	    rs.add("charstring", charString);

// 	// Get the list to be returned to R.
// 	rl = rs.getReturnList();

//     } catch(std::exception& ex) {
// 	exceptionMesg = copyMessageToR(ex.what());
//     } catch(...) {
// 	exceptionMesg = copyMessageToR("unknown reason");
//     }

//     if(exceptionMesg != NULL)
// 	Rf_error(exceptionMesg);

//     return rl;

// }

RcppExport SEXP ReadCharsWithNCL(SEXP params) {

    Rcpp::List list(params);
    #if defined(FILENAME_AS_NEXUS)
        string filename = "'" + list["filename"] + "'";
    #else
	string filename = list["filename"];
    #endif
    
    BASICCMDLINE reader;
    reader.Initialize(const_cast < char* > (filename.c_str()));
    
    NxsString charStringNxs;
    reader.RReturnCharacters(charStringNxs,
			     list["allchar"], 			// boolean flags 
			     list["polymorphictomissing"],
			     list["levelsall"], 
			     list["returnlabels"]);
    string charString = charStringNxs.c_str();

    return Rcpp::List::create(Rcpp::Named("charstring") = charString);
}


