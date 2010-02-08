#!/usr/bin/env python
import sys
import os
import subprocess
import filecmp
import itertools
import cStringIO
import shutil
                
def runTest(inArgPath, 
            outParent, 
            compareOut=True, 
            copyOutput=False, 
            strictLevel=2, absentOnlyMode=False,
            copyFile=False,
            external=False,
            invalid=False,
            parseOutput=True):
    if invalid:
        compareOut = False
    if external:
        fileNamePatterns = ["*.dat", "*.phy", "*.fasta", "*.txt", "*.fas", "*.nex", "*.tre", "*.aln"]
    else:
        fileNamePatterns = ["*.nex", "*.tre" ]
    if not os.path.exists(inArgPath):
        sys.exit("input file " + inArgPath + " does not exist")
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
        print >>sys.stderr, "No files to test (looking for " + " ".join(fileNamePatterns) + " files)"
        sys.exit(0)
    
    if compareOut:
        if not os.path.exists(outParent):
            sys.exit(outParent + " does not exist")
        if not os.path.isdir(outParent):
            sys.exit(outParent + " must be a directory")
    #if os.path.samefile(outParent, inParent):
    #    sys.exit("the input and output directories cannot be the same")
    
    for f in inFiles:
        pref = os.path.join(".roundTripNCLOut.nex")
        
        invokedInFile = os.path.join(inParent, f)

        toCheck = [invokedInFile]
        if compareOut:
            expectedOut = os.path.join(outParent, f)
            if not copyOutput and not os.path.exists(expectedOut):
                sys.exit("Expected output file " + expectedOut + " does not exist.")
            if parseOutput:
                toCheck.append(expectedOut)
            

        if absentOnlyMode and os.path.exists(expectedOut):
            continue
        
        for pRep, inFile in enumerate(toCheck):
            if pRep == 0:
                n = open(".roundTripNCLInputName.txt", "w")
                n.write('%s\n' % os.path.abspath(inFile))
                n.close()
                if not invalid:
                    n = open(".roundTripNCLExpectedOutName.txt", "w")
                    n.write('%s\n' % os.path.abspath(expectedOut))
                    n.close()
                if copyFile:
                    shutil.copy2(inFile, ".roundTripNCLInput.nex")
            inFile = os.path.abspath(inFile)
            tf = pref
            if os.path.exists(tf):
                os.remove(tf)
            tFileObj = open(tf, "w")
            if external and pRep == 0:
                fileFormat = f.split('_')[0]
                invocation = [normalizer, "-s%d" % strictLevel, "-f%s" % fileFormat, inFile]
            else:
                invocation = [normalizer, "-s%d" % strictLevel, inFile]
            sys.stderr.write('"%s"\n' % '" "'.join(invocation))
            retCode = subprocess.call(invocation, stdout=tFileObj)
            tFileObj.close()
            if invalid:
                if retCode == 0:
                    sys.exit("Call to " + normalizer + " accepted the invalid file " + inFile)
            elif retCode != 0:
                sys.exit("Call to " + normalizer + " failed for "+ inFile)
            if copyOutput and pRep == 0:
                shutil.copy2(tf, expectedOut)
            elif compareOut:
                e = file(expectedOut, "rU") 
                o = file(tf, "rU")
                eit = iter(e)
                oit = iter(o)
                n = 1
                while True:
                    try:
                        eLine = eit.next()
                        try:
                            oLine = oit.next()
                        except StopIteration:
                            oLine = ""
                    except StopIteration:
                        eLine = ""
                        try:
                            oLine = oit.next()
                        except StopIteration:
                            break
                    if eLine != oLine:
                        basicMsg = "Call to " +  os.path.basename(normalizer) + " failed to produce expected output for " + inFile + ".\n"
                        fileInfo = "\nThe output files:\n " + tf + "\nand the expected output:\n " + expectedOut + "\ndiffer at line " + str(n) + "\n"
                        diff = "Got:\n" + oLine + "\nExpected:\n" + eLine
                        sys.exit(basicMsg + fileInfo + diff)
                e.close()
                o.close()
                print >>sys.stderr, sys.argv[0], ": Call to", normalizer, "succeeded for", inFile
            if os.path.exists(tf):
                os.remove(tf)


from optparse import OptionParser
usage = "%prog [options] <path to normalizer> <input dir> <expected output dir>"    
parser = OptionParser(usage=usage, add_help_option=True, version = 0.2, description="Tool for using NCL's normalizer to test NCL functionality by comparing round trip output to expected output")
parser.add_option("-a", "--auto",
              dest="auto",
              default=False,
              action="store_true",
              help="Get <input dir> and <expected output dir> batch file (default is ~/.ncl_round_triprc)")   
parser.add_option("-b", "--batch",
              dest="batch",
              default="~/.ncl_round_triprc",
              help="The path to a file with and pairs of lines indicating <input dir> and <expected output dir>")   
parser.add_option("-p", "--parse",
              dest="parseOnly",
              default=False,
              action="store_true",
              help="Parse only (do not compare output from \"normalizer\" to expected output)")   
parser.add_option("-c", "--copy",
              dest="copy",
              default=False,
              action="store_true",
              help="Copy file to ./.roundTripNCLInput.nex before parsing")   
parser.add_option("-f", "--force-equal",
              dest="force",
              default=False,
              action="store_true",
              help="Copies the generated output to the expected output loctation (does not test validity -- this is used for creating an expected output when you are confident that the normalizer works for a set of files)")   
parser.add_option("-m", "--missing",
              dest="absent",
              default=False,
              action="store_true",
              help="Only run if the expected output is missing normalizer.")
parser.add_option("-e", "--external",
              dest="external",
              default=False,
              action="store_true",
              help="Treat the file as an external (non-NEXUS) format.  The format name should be the prefix to the file name, separated from the rest of the name by an underscore.")
parser.add_option("-i", "--invalid",
              dest="invalid",
              default=False,
              action="store_true",
              help="Demand that the parser rejects the file(s) because they are invalid.")
parser.add_option("-o", "--output-invalid",
              dest="parseOutput",
              default=True,
              action="store_false",
              help="Do not check the output file by a second round-trip test")
parser.add_option("-s", "--strictness",
              dest="strict",
              default=2,
              type="int",
              help="Strictness argument passed in as a -s argument to the normalizer.")   

(options, args) = parser.parse_args()
if options.batch != "~/.ncl_round_triprc":
    options.auto = True
if options.absent:
    options.force = True
if options.auto:
    if len(args) != 1:
        sys.exit("Expecting the path to the normalizer as the only argument when running in auto or batch mode")
    resourceFilePath = os.path.expanduser(options.batch)
    if not os.path.exists(resourceFilePath):
        m = "%s does not exist." % resourceFilePath
        if options.batch != "~/.ncl_round_triprc":
            sys.exit(m)
        print(m + " Test stepped skipped.")
        sys.exit(0)
    resourceFileStream = open(resourceFilePath, "rU")
else:
    if len(args) != 3:
        if not options.invalid:
            sys.exit("""Expecting 3 arguments:
             1. the path to NEXUS file as a command-line argument, interprets the content, and writes  what it understands (as NEXUS) to standard output.
             2. the path to the parent of the inputs files to test
             3. the path to the parent of the expected output.
             """)
        elif len(args) != 2:
            sys.exit("""Expecting 2 arguments:
             1. the path to NEXUS file as a command-line argument, interprets the content, and writes  what it understands (as NEXUS) to standard output.
             2. the path to the parent of the inputs files to test
             """)
    if not options.invalid:
        fs = "%s\n%s\n" % (args[1], args[2])
    else:
        fs = "%s\n" % (args[1])
    resourceFileStream = cStringIO.StringIO(fs)
normalizer = os.path.abspath(args[0])
if not os.path.exists(normalizer):
    sys.exit(normalizer + " does not exist")

lineIter = iter(resourceFileStream)
inputParentPath = "#"
try:
    while True:
        inputParentPath = "#"
        while len(inputParentPath) < 1 or inputParentPath.strip().startswith("#"):
            inputParentPath = lineIter.next().strip()
        if options.invalid:
            outputParentPath = None
        else:
            outputParentPath = "#"
            while  len(outputParentPath) < 1 or outputParentPath.strip().startswith("#"):
                outputParentPath = lineIter.next().strip()
        #print "next round: ", inputParentPath,"\n ", outputParentPath
        runTest(inputParentPath, outputParentPath, not options.parseOnly, options.force, options.strict, options.absent, options.copy, options.external, options.invalid, options.parseOutput)
        
except StopIteration:
    if len(inputParentPath) > 0 and inputParentPath != "#":
        sys.exit("No matching output path found for %s" % inputParentPath)
