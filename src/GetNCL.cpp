#include <Rcpp.h>
#include "ncl/nxsmultiformat.h"

//#define NEW_TREE_RETURN_TYPE

NxsString contData(NxsCharactersBlock& charBlock, NxsString& charString, 
		   const int& eachChar, const int& nTax) {
    for (int taxon=0; taxon < nTax; ++taxon) {
	double state=charBlock.GetSimpleContinuousValue(taxon,eachChar);			
	if (state==DBL_MAX) {
	    charString+="NA";
	}
	else {					    
	    char buffer[100];
	    sprintf(buffer, "%.10f", state);
	    charString+=buffer; 
	}
	
	if (taxon+1 < nTax) {
	    charString+=',';
	}
    }
    return charString;
}


NxsString stdData(NxsCharactersBlock& charBlock, NxsString& charString, const int& eachChar,
		  const int& nTax, bool polyconvert) {
    for (int taxon=0; taxon<nTax; ++taxon) {

	int stateNumber=charBlock.GetInternalRepresentation(taxon, eachChar, 0);

	if(charBlock.IsMissingState(taxon, eachChar)) {
	    charString+="NA";
	}
	else if (charBlock.GetNumStates(taxon, eachChar)>1) {
	    if(polyconvert) {
		charString+="NA";
	    }
	    else {
		charString+='"';
		charString+='{';
		for (unsigned int k=0; k < charBlock.GetNumStates(taxon, eachChar); ++k) {
		    charString += charBlock.GetInternalRepresentation(taxon, eachChar, k);	
		    if (k+1 < charBlock.GetNumStates(taxon, eachChar)) {
			charString+=',';
		    }
		}
		charString+='}';
		charString+='"';
	    }
	}
	else {
	    charString+='"';
	    charString+=stateNumber;
	    charString+='"';
	}
	if (taxon+1 < nTax) {
	    charString+=',';
	}
    }
    return charString;
}


extern "C" SEXP GetNCL(SEXP params, SEXP paramsVecR) {

    Rcpp::List list(params);
    Rcpp::LogicalVector paramsVec(paramsVecR);

    bool charall = paramsVec[0];
    bool polyconvert = paramsVec[1];
    bool levelsUnif = paramsVec[2];
    bool returnTrees = paramsVec[3];
    bool returnData = paramsVec[4];

    int nCharToReturn = 0;

    std::vector<std::string> dataTypes;      //vector of datatypes for each character block
    std::vector<int> nbCharacters;           //number of characters for each character block
    std::vector<std::string> dataChr;        //characters
    std::vector<std::string> charLabels;     //labels for the characters
    std::vector<std::string> stateLabels;    //labels for the states
    std::vector<int> nbStates;               //number of states for each character (for Standard datatype)
#   if defined (NEW_TREE_RETURN_TYPE)
        std::vector<std::string> taxonLabelVector; //Index of the parent. 0 means no parent.
        std::vector<unsigned> parentVector; //Index of the parent. 0 means no parent.
        std::vector<double> branchLengthVector; 
#   else
        std::vector<std::string> trees;          //vector of Newick strings holding the names
#   endif
    std::vector<std::string> treeNames;      //vector of tree names
    std::vector<std::string> taxaNames;      //vector of taxa names
    std::string errorMsg;                    //error message

    std::vector<bool> test(3);
    test[0] = charall;
    test[1] = polyconvert;
    test[2] = levelsUnif;

#   if defined(FILENAME_AS_NEXUS)
    std::string filename = "'" + list["fileName"] + "'";
#   else
    std::string filename = list["fileName"];
#   endif

    MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);

    /* make NCL less strict */
    NxsTreesBlock * treesB = nexusReader.GetTreesBlockTemplate();
    treesB->SetAllowImplicitNames(true);
    nexusReader.cullIdenticalTaxaBlocks(true);
    /* End of making NCL less strict */

    MultiFormatReader::DataFormatType fileFormat =  MultiFormatReader::NEXUS_FORMAT;
    std::string fileFormatString = list["fileFormat"];
    if (!fileFormatString.empty())
        {
        fileFormat = MultiFormatReader::formatNameToCode(fileFormatString);
        if (fileFormat == MultiFormatReader::UNSUPPORTED_FORMAT)
            {
            std::string m = "Unsupported format \"";
            m.append(fileFormatString);
            m.append("\"");
            Rcpp::List res = Rcpp::List::create(Rcpp::Named("ErrorMsg") = m);
	        return res;
            }
        }
/* 

fileFormatString should be one of these: 	"nexus",
								"dnafasta",
								"aafasta",
								"rnafasta",
								"dnaphylip",
								"rnaphylip",
								"aaphylip",
								"discretephylip",
								"dnaphylipinterleaved",
								"rnaphylipinterleaved",
								"aaphylipinterleaved",
								"discretephylipinterleaved",
								"dnarelaxedphylip",
								"rnarelaxedphylip",
								"aarelaxedphylip",
								"discreterelaxedphylip",
								"dnarelaxedphylipinterleaved",
								"rnarelaxedphylipinterleaved",
								"aarelaxedphylipinterleaved",
								"discreterelaxedphylipinterleaved",
								"dnaaln",
								"rnaaln",
								"aaaln",
								"phyliptree",
								"relaxedphyliptree",
								"nexml",
								"dnafin",
								"aafin",
								"rnafin"
							};
							*/
    try {
	nexusReader.ReadFilepath(const_cast < char* > (filename.c_str()), fileFormat);  
    }
    catch (NxsException &x) {
	errorMsg = x.msg;
	Rcpp::List res = Rcpp::List::create(Rcpp::Named("ErrorMsg") = errorMsg);
	return res;
    }

    const unsigned nTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    for (unsigned t = 0; t < nTaxaBlocks; ++t) {
	/* Get blocks */
	const NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(t);
	const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
	const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);

	int nTax = taxaBlock->GetNumTaxonLabels();
	
	/* Get taxa names */
	for (int j=0; j < nTax; ++j) {	
	    taxaNames.push_back (taxaBlock->GetTaxonLabel(j));
	}

#   if defined (NEW_TREE_RETURN_TYPE)
        taxonLabelVector.reserve(nTax);
        parentVector.reserve(2*nTax);
        branchLengthVector.reserve(2*nTax);
#   endif
	/* Get trees */
	if (returnTrees) {
	    if (nTreesBlocks == 0) {
		continue;
	    }
	    for (unsigned i = 0; i < nTreesBlocks; ++i) {
		NxsTreesBlock* treeBlock = nexusReader.GetTreesBlock(taxaBlock, i);
		const unsigned nTrees = treeBlock->GetNumTrees();
		if (nTrees > 0) {
		    for (unsigned k = 0; k < nTrees; k++) {
#           if defined(NEW_TREE_RETURN_TYPE)
                taxonLabelVector.clear();
                parentVector.clear();
                branchLengthVector.clear();
                
                const NxsFullTreeDescription & ftd = treeBlock->GetFullTreeDescription(k); 
                treeNames.push_back(ftd.GetName());
                NxsSimpleTree simpleTree(ftd, -1, -1.0);
                std::vector<const NxsSimpleNode *> ndVector =  simpleTree.GetPreorderTraversal();
                unsigned internalNdIndex = nTax;
                for (std::vector<const NxsSimpleNode *>::const_iterator ndIt = ndVector.begin(); ndIt != ndVector.end(); ++ndIt)
                    {
                    NxsSimpleNode * nd = (NxsSimpleNode *) *ndIt;
                    unsigned nodeIndex;
                    if (nd->IsTip())
                        {
                        nodeIndex = nd->GetTaxonIndex();
                        taxonLabelVector.push_back(taxaNames[nodeIndex]);
                        std::cout << " leaf node # = " <<  nodeIndex << '\n';
                        }
                    else
                        {
                        nodeIndex = internalNdIndex++;
                        nd->SetTaxonIndex(nodeIndex);
                        std::cout << " internal node # = " << nd->GetTaxonIndex()  << '\n';
                        }
                    if (parentVector.size() < nodeIndex + 1)
                        {
                        parentVector.resize(nodeIndex + 1);
                        }
                    if (branchLengthVector.size() < nodeIndex + 1)
                        {
                        branchLengthVector.resize(nodeIndex + 1);
                        }
                    NxsSimpleEdge edge = nd->GetEdgeToParent();
    
                    NxsSimpleNode * par = 0L;
                    par = (NxsSimpleNode *) edge.GetParent();
                    if (par != 0L)
                        {
                        parentVector[nodeIndex] = 1 + par->GetTaxonIndex();
                        branchLengthVector[nodeIndex] = edge.GetDblEdgeLen();
                        }
                    else
                        {
                        parentVector[nodeIndex] = 0;
                        branchLengthVector[nodeIndex] = -1.0;
                        }
                    }
                std::cout << "Parents = [";
                for (std::vector<unsigned>::const_iterator nIt = parentVector.begin(); nIt != parentVector.end(); ++nIt)
                    {
                    std::cout << *nIt << ", ";				
                    }
                std::cout << "]\nbranch lengths = [";
                for (std::vector<double>::const_iterator nIt = branchLengthVector.begin(); nIt != branchLengthVector.end(); ++nIt)
                    {
                    std::cout << *nIt << ", ";				
                    }
                std::cout << "]\n";




#           else
    			NxsString trNm = treeBlock->GetTreeName(k);
	    		treeNames.push_back(trNm);
    			NxsString ts = treeBlock->GetTreeDescription(k);
    			trees.push_back (ts);
#           endif
		    }
		}
		else {
		    continue;
		}
	    }
	}
	
	/* Get data */
	if (returnData) {
	    for (unsigned k = 0; k < nCharBlocks; ++k) {
		NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, k);
		
		if (nCharBlocks == 0) {
		    continue;
		}
		else {
		    NxsString dtType = charBlock->GetNameOfDatatype(charBlock->GetDataType());
		    dataTypes.push_back(dtType);
		
		    if (charall) {
			nCharToReturn=charBlock->GetNCharTotal();
		    }
		    else {
			nCharToReturn=charBlock->GetNumIncludedChars();
		    }
		    nbCharacters.push_back (nCharToReturn);
		    for (int eachChar=0; eachChar < nCharToReturn; ++eachChar) { //We only pass the non-eliminated chars
			NxsString charLabel=charBlock->GetCharLabel(eachChar);
			if (charLabel.length()>1) {
			    charLabels.push_back (charLabel);
			}
			else {
			    charLabels.push_back ("standard_char"); //FIXME: needs to fixed for sequence data
			}
			
			NxsString tmpCharString;
			if (std::string("Continuous") == dtType) {
			    tmpCharString = contData(*charBlock, tmpCharString, eachChar, nTax);
			    nbStates.push_back (0);			    
			}
			else {
			    if (std::string("Standard") == dtType) {			    
				tmpCharString = stdData(*charBlock, tmpCharString, eachChar, nTax,
							polyconvert);
				unsigned int nCharStates = charBlock->GetNumObsStates(eachChar, false);
				nbStates.push_back (nCharStates);
				for (unsigned int l=0; l < nCharStates; ++l) {
				    NxsString label = charBlock->GetStateLabel(eachChar, l);
				    stateLabels.push_back (label);
				}
			    }
			    else {
				if (std::string("DNA") == dtType) {
				    for (int taxon=0; taxon < nTax; ++taxon) {
					for (int eachChar=0; eachChar < nCharToReturn; ++eachChar) {
					    unsigned int nCharStates = charBlock->GetNumStates(taxon, eachChar);
					    if (charBlock->IsGapState(taxon, eachChar)) {
						tmpCharString += "-";
					    }
					    else {
						if (charBlock->IsMissingState(taxon, eachChar)) {
						    tmpCharString += "?";
						}
						else {
						    if (nCharStates == 1) {
							tmpCharString += charBlock->GetState(taxon, eachChar, 0);
						    }
						    else {
							tmpCharString += "?"; //FIXME
						    }			    
						}
					    }
					}
				    }
				}
				else { // other type of data not yet supported
				    tmpCharString = "";
				    nbStates.push_back (0);
				    stateLabels.push_back (std::string(""));
				}
			    }
			}
			std::string charString = "c(" + tmpCharString + ");";
			dataChr.push_back (charString);
		    }				
		}  
	    }
	}
    }

    /* Prepare list to return */
    Rcpp::List res = Rcpp::List::create(Rcpp::Named("taxaNames") = taxaNames,
					Rcpp::Named("treeNames") = treeNames,
#               if defined (NEW_TREE_RETURN_TYPE)
                    Rcpp::Named("parentVector") = parentVector,
                    Rcpp::Named("branchLengthVector") = branchLengthVector,
#               else
					Rcpp::Named("trees") = trees,
#               endif
                    Rcpp::Named("dataTypes") = dataTypes,
					Rcpp::Named("nbCharacters") = nbCharacters,
					Rcpp::Named("charLabels") = charLabels,
					Rcpp::Named("nbStates") = nbStates,
					Rcpp::Named("stateLabels") = stateLabels,
					Rcpp::Named("dataChr") = dataChr,
					Rcpp::Named("Test") = test);
    return res;				
}
