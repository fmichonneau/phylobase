//  Copyright (C) 2007-2008 Brian O'Meara & Derrick Zwickl
//  A modification of the RcppExample.cp file of the R/C++ interface class
//  library (see below) to use for loading trees and data from Nexus into R. 
//  Licensing as below.


// RcppExample.cpp: Part of the R/C++ interface class library, Version 5.0
//
// Copyright (C) 2005-2006 Dominick Samperi
//
// This library is free software; you can redistribute it and/or modify it 
// under the terms of the GNU Lesser General Public License as published by 
// the Free Software Foundation; either version 2.1 of the License, or (at 
// your option) any later version.
//
// This library is distributed in the hope that it will be useful, but 
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public License 
// along with this library; if not, write to the Free Software Foundation, 
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 

#include <fstream>

#include <Rcpp.h>

//TODO - figure out where/why length is getting defined as Rf_length so that
//this isn't necessary to compile
// DE: This comes from Rinternals.h via the R_NO_REMAP used in Rcpp.h -- I have
//     found that some libraries fail to build when R defines length, error, ... 
//     so I prefer to use the safer and explicit Rf_ prefixes as eg Rf_error
//     used here three times
#define  Rf_length length

#include <ncl/ncl.h>
#include "NCLInterface.h"

using namespace std;

//This function receives a list of parameters from
//R, which can then be extracted below
RcppExport SEXP ReadWithNCL(SEXP params) {

    SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
    char* exceptionMesg=NULL;
    try {


		// Get parameters in params - only 1 is gotten now
		RcppParams rparam(params);
		string filename = "'";
		filename+=rparam.getStringValue("filename");
		filename+="'";
		
		BASICCMDLINE reader;
		
		//this is where the reader would be passed the filename to read
		//reader.Run(filename.c_str()); //Will not compile
		//reader.Run(NULL);
		reader.Initialize(const_cast < char* > (filename.c_str()));
	
		string filenameString = "I was told to read a file named ";
		filenameString += filename;
	
		//Various calls to the reader can be made here to fill the various strings
		//string treeString = "This was passed in the tree string";
		NxsString treeStringNxs;
		reader.RReturnTrees(treeStringNxs);
		string treeString=treeStringNxs.c_str();
	
		string otherString = "";
	
		//string discreteString = "This was passed in the discrete string";
		
		NxsString characterStringNxs;
		reader.RReturnCharacters(characterStringNxs,false, true, false);
		//string characterString=characterStringNxs.c_str();
		
		string continuousString = "This was passed in the continuous string";
		//reader.GetContinuousString();
	
		//Just to test to see if reader is running
		string testString = reader.TestRunning();
	
		// Build result set to be returned as a list to R.
		RcppResultSet rs;
		
		//if the various strings are nonempty, add them to the ResultSet
		if(filenameString.length() > 0)
			rs.add("filenamestring", filenameString);
		if(treeString.length() > 0)
			rs.add("treestring", treeString);
		if(otherString.length() > 0)
			rs.add("otherstring", otherString);
		if(characterStringNxs.length() > 0)
			rs.add("characterstring", characterStringNxs.c_str());
		if(continuousString.length() > 0)
			rs.add("continuousstring", continuousString);
		if(testString.length() > 0)
			rs.add("teststring", testString);

		
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

RcppExport SEXP ReadTreesWithNCL(SEXP params) {
	
    SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
    char* exceptionMesg=NULL;
    try {
		
		
		// Get parameters in params - only 1 is gotten now
		RcppParams rparam(params);
		string filename = "'";
		filename+=rparam.getStringValue("filename");
		filename+="'";
		
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

RcppExport SEXP ReadCharsWithNCL(SEXP params) {
	
    SEXP  rl=R_NilValue; // Use this when there is nothing to be returned.
    char* exceptionMesg=NULL;
    try {
		
		
		// Get parameters in params - only 1 is gotten now
		RcppParams rparam(params);
		string filename = "'";
		filename+=rparam.getStringValue("filename");
		filename+="'";
		
		bool allchar = rparam.getBoolValue("allchar");
		bool levelsall=rparam.getBoolValue("levelsall");
		bool polymorphictomissing=rparam.getBoolValue("polymorphictomissing");
		
		BASICCMDLINE reader;
		
		//this is where the reader would be passed the filename to read
		//reader.Run(filename.c_str()); //Will not compile
		//reader.Run(NULL);
		reader.Initialize(const_cast < char* > (filename.c_str()));
		
		NxsString charStringNxs;
		reader.RReturnCharacters(charStringNxs,allchar, polymorphictomissing, levelsall);
		string charString=charStringNxs.c_str();
		
		// Build result set to be returned as a list to R.
		RcppResultSet rs;
		
		//if the various strings are nonempty, add them to the ResultSet
		if(charString.length() > 0)
			rs.add("charstring", charString);
		
		
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



