#!/usr/bin/env python
import sys
import subprocess
import filecmp
import itertools

if len(sys.argv) < 4:
    sys.exit('''Expecting the path to a normalizer that takes:
     1. the path to NEXUS file as a command-line argument, interprets the content, and writes  what it understands (as NEXUS) to standard output.
     2. the path to the parent of the inputs files to test
     3. the path to the parent of the expected output.
     ''')
import os
normalizer = os.path.abspath(sys.argv[1])
if not os.path.exists(normalizer):
    sys.exit(normalizer + " does not exist")
    
fileNamePatterns = ['*.nex', '*.tre' ]

inArgPath = os.path.abspath(sys.argv[2])
if not os.path.exists(inArgPath):
    sys.exit(inArgPath + " does not exist")
if os.path.isdir(inArgPath):
    import glob
    inParent = inArgPath
    inFiles = []
    for pat in fileNamePatterns:
         inFiles.extend([os.path.basename(i) for i in glob.glob(os.path.join(inArgPath, pat))])
else:
    inParent = os.path.dirname(inArgPath)
    inFiles = [os.path.basename(inArgPath)]

if not inFiles:
    print >>sys.stderr, "No files to test (looking for " + ' '.join(fileNamePatterns) + " files)"
    sys.exit(0)

outParent = os.path.abspath(sys.argv[3])
if not os.path.exists(outParent):
    sys.exit(outParent + " does not exist")
if not os.path.isdir(outParent):
    sys.exit(outParent + " must be a directory")
if os.path.samefile(outParent, inParent):
    sys.exit("the input and output directories cannot be the same")

for f in inFiles:
    expectedOut = os.path.join(outParent, f)
    if not os.path.exists(expectedOut):
        sys.exit("Expected output file " + expectedOut + " does not exist.")

    tf = os.path.join(".roundTripNCLOut.nex")

    inFile = os.path.join(inParent, f)
    tFileObj = open(tf, 'w')
    retCode = subprocess.call([normalizer, inFile], stdout=tFileObj)
    tFileObj.close()
    if retCode != 0:
        sys.exit("Call to " + normalizer + " failed for " +inFile)
    else:
        e = file(expectedOut, 'rU')
        o = file(tf, 'rU')
        for n, lineTuple in enumerate(itertools.izip(e, o)):
            eLine, oLine = lineTuple
            if eLine != oLine:
                basicMsg = "Call to " +  os.path.basename(normalizer) + " failed to produce expected output for " + inFile + ".\n"
                fileInfo = "\nThe output files:\n " + tf + "\nand the expected output:\n " + expectedOut + "\ndiffer at line " + str(n) + '\n'
                diff = "Got:\n" + oLine + "\nExpected:\n" + eLine
                sys.exit(basicMsg + fileInfo + diff)
        e.close()
        o.close()
        print >>sys.stderr, sys.argv[0], ": Call to", normalizer, "succeeded for", inFile
    if os.path.exists(tf):
        os.remove(tf)
        
