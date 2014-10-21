#!/usr/bin/env python
#combines reduce.py, aligning, replace.py and concat.py into one script
#created 31.07.2014
#changes 06.09.2014
#reslp

# known issues:
# can't change aligner settings - in align
# NNs at the beginning of alignment - in concat - DONE
# windows compatibility not checked
# problem if file list is given! doesn't work properly. - in reduce - DONE
# if concat only is used it works only for single lined sequences! - DONE
# input dir has to be specified when concat is used!

import argparse
import glob
import os
import platform
from subprocess import Popen, PIPE
import sys

parser = argparse.ArgumentParser(description="%(prog)s (build: Aug 14, 2014) will create concatened alignments from seperate single locus files")
parser.add_argument("-t", dest="taxonfile", help="specify place of file containing desired taxon set")
parser.add_argument("-d", dest="directory", help="specify directory with single locus files that should be used")
parser.add_argument("-v", action="version", help="outputs version of concat script", version="%(prog)s Version 0.1 (06.09.14)")
parser.add_argument("-c", dest="concat",action="store_true", help="concatenate sequences")
parser.add_argument("-r", dest="reduce",action="store_true", help="reduces sequences to desired set of taxa")
parser.add_argument("-a", dest="align",action="store_true", help="align set of sequences. used aligner: mafft")
parser.add_argument("-p", dest="replace",action="store_true", help="replace - with ? at beginning and end of aligned sequences")
parser.add_argument("-N", dest="NNN",action="store_true", help="adds Ns between loci to seperate them in concatenated alignment")
#parser.add_argument("-aligner", dest="align_param",action="store_true", help="Arguments that should be sent to the aligner (not yet implemented)")
parser.add_argument("-clean", dest="clean",action="store_true", help="performs a clean run: removes all intermediate files")
#parser.add_argument("-trim", dest="trim",action="store_true", help="reduce alignment to include only sites in at least number of sequences (not yet implemented)")
parser.add_argument("-i", dest="input_file", nargs="*", type=str, help="names & path to input files that should be used. Seperate filenames with space. e.g.: \"its.fas lsu.fas\"")

#-------------------------------------------------------------- set things up
#check for number of commandline arguments, print help if none given
if len(sys.argv)<2:
   parser.print_help()
   sys.exit()
Args = parser.parse_args()
print Args.input_file
print Args.taxonfile
#print Args
#parser.print_help()
#sys.exit()


if Args.input_file == None and Args.directory == None:
	print "Need input files. Use either -i oder -d. Try also -h for more information."
	sys.exit()
if Args.taxonfile == None:
	print "Need file containing taxon IDs. Use -t. Try also -h for more information."
	sys.exit()


#reading Taxa from taxonfile
TaxonFile = open(Args.taxonfile, "r")
TaxonList = []
for Line in TaxonFile:
	TaxonList.append(Line.strip("\n"))
TaxonListOutput = TaxonList [:] #create Taxon List for Output
TaxonFile.close()

#get working directory:
WD = os.getcwd()
WD = str(WD)
print "Working directory: ", WD

#sys.exit()
	
if Args.input_file != None and Args.directory != None:
	print "Can't use both input files (-i flag) and directory (-d flag). Please use only one option."
	sys.exit()
	
#create directory for temp files
if not os.path.exists("tmp"):
    os.makedirs("tmp")

TMPD = os.path.join(WD,"tmp")
print TMPD
#sys.exit()

#getting input files ready

if Args.input_file != None:
	#print Args.input_file
	SeqFileList = Args.input_file
	OnlyFileName = []
	OnlyDirName = []
	for Name in SeqFileList:
		OnlyDirName.append(os.path.dirname(Name))
		OnlyFileName.append(os.path.basename(Name))
	print OnlyFileName	
	print OnlyDirName
	#sys.exit()

if Args.directory != None: #open directory with raw files
	OnlyFileName = []
	OnlyDirName = []
	Directory = Args.directory
	os.chdir(Directory)
	OnlyFileName = glob.glob("*.fas")
	for File in OnlyFileName:
		OnlyDirName.append(Directory)
	os.chdir(WD)
	
OS = platform.system()
print OS
	
#sys.exit()
#print TaxonList
#-------------------------------------------------------------- Functions 
def reduce(): #reduces alignments to desired set of taxa
	file=0
	TotTax = len(TaxonList)
	for RawFile, RawDir in zip(OnlyFileName, OnlyDirName):
		file=1
		Filename = os.path.join(TMPD, (RawFile + "_reduced"))
		Outfile = open(Filename,"w")
		Sequenzfile = open(os.path.join(RawDir,RawFile), "U")
		#create a blank Sequencelist for given number of taxa
		SequenceList = [""]*len(TaxonList)
		#Open and read Files
		Found = 0
		for Taxon in TaxonList:
			for Line in Sequenzfile:
				if Line.startswith(">"+Taxon):
					#print Line.strip("\n")
					Outfile.write(Line)
					Found = 1
					continue
				if Found == 1:
					if Line.find(">") == -1:
						#print Line.strip("\n")
						Outfile.write(Line.strip("\n"))
					else:
						Outfile.write("\n")
						Found = 0						
			Sequenzfile.seek(0)
		Outfile.close()
	if file == 0:
		print "Reduce: No files found in specified directory"

def align(): # aligns sequence files
	file = 0
	if OS == "Darwin" or OS == "Linux":
		for RedFile in glob.glob(os.path.join(TMPD,"*_reduced")):
			print RedFile
			file = 1
			try:
				process = Popen(["mafft","--auto",RedFile], stdout=PIPE)
				(output, err) = process.communicate()
				exit_code = process.wait()
			except:
				print "mafft executable not found. Is it in your path?"
				sys.exit()
			Outfile = open(RedFile+"_aligned","w")
			Outfile.write(output)
			Outfile.close()
		if file == 0:
			print "Align: No \"*_reduced\" files found in temp directory", TMPD
		return
	
	if OS == "Windows":
		for RedFile in glob.glob(TMPD+"*_reduced"):
			print RedFile
			file = 1
			process = Popen(["mafft","--auto",RedFile], stdout=PIPE)
			(output, err) = process.communicate()
			exit_code = process.wait()
			Outfile = open(RedFile+"_aligned","w")
			Outfile.write(output)
			Outfile.close()
		if file == 0:
			print "Align: No \"*_reduced\" files found in temp directory", TMPD

def replace():	# replaces - with ? at the beginning and end of fasta alignments
	def replace_begin (Sequenz_0):
		zahl = 0
		for Base in Sequenz_0:
			if Base == "-":
				Sequenz_0[zahl] = "?"
			else:
				break
			zahl += 1
		return "".join(Sequenz_0)
	# end replace_begin
		
	#replace all - from the end of the sequence
	def replace_end (Sequenz_0):
		zahl2 = -1
		for Base in Sequenz_0:
			if Sequenz_0[len(Sequenz_0) + zahl2] == "-":
				Sequenz_0[len(Sequenz_0) + zahl2] = "?"	
			else:
				break
			zahl2 -= 1
		return "".join(Sequenz_0)
	#end replace_end
	file = 0
	for AlFile in glob.glob(os.path.join(TMPD,"*_aligned")):
		file = 1
		TaxonList = []	
		File = open(AlFile, "U")
		for Line in File:
			if Line[0] == ">":
				TaxonList.append(Line)
		
		MaxSeq = len(TaxonList)
		SequenceList = [""] * MaxSeq
		File.seek(0)
		SeqNumber = -1	
		for Line in File:
			if Line[0] != ">":
				SequenceList[SeqNumber] += Line
			else:
				SeqNumber += 1
		SeqNumber = 0
		for Single in SequenceList:
			SequenceList[SeqNumber] = Single.replace("\n", "")
			SequenceList[SeqNumber] = replace_begin(list(SequenceList[SeqNumber]))
			SequenceList[SeqNumber] = replace_end(list(SequenceList[SeqNumber]))
			TaxonList[SeqNumber] = TaxonList[SeqNumber].replace("\n","")
			SeqNumber += 1
		
		Outfile = open(AlFile+"_replaced", "w")
		for Taxon in range(0,MaxSeq):
			#print TaxonList[Taxon]
			Outfile.write(TaxonList[Taxon]+"\n")
			#print "seqlist:"
			#print SequenceList
			Outfile.write(SequenceList[Taxon])
			Outfile.write("\n")
		Outfile.close()
		File.close()
	if file == 0:
		print "Replace: No \"*_aligned\" files found in specified directory", TMPD
		
def concat():
	TaxonListOutput = TaxonList [:] #Taxon List for Output

	def add_missing(SList):#add ? marks when sequence is missing
		Index = 0
		LongestItem = max(SList, key=len)	
		for Item in SList:
			if len(Item) <= len(LongestItem):
				SList[Index] += "?" * (len(LongestItem) - len(Item))
				if Args.NNN == True:
					if SList[Index] != "": #avoid NN at the beginning of alignment
						SList[Index] += "N" * 30
			Index += 1
		return SList
	# end of add_missing()
		
	def add_to_taxon(WhichTaxon): #add 0 for missing locus to name
		TaxonListOutput[WhichTaxon] += "O"
		return
	# end of add_to_taxon
	#get total number of taxa
	TotTax = len(TaxonList)	
	#create a blank Sequencelist for given number of taxa
	SequenceList = [""]*len(TaxonList)	
	file = 0
	
	#get whole sequence / to deal with interleaved file format
	def read_sequence(SeqFile):
		Seq = ""
		#SeqFile.next()
		for Line in SeqFile:
			if ">" not in Line:
				Seq += Line.strip("\n")
			if ">" in Line:
				break

		return Seq

	for ReplFile in glob.glob(os.path.join(TMPD,"*_replaced")):
		file = 1
		#print ReplFile 
		TaxonNum = 0
		LineIndex = 0
		LineNum = 0	
		Sequenzfile = open(ReplFile, "U")
		TaxonNum = 0
		for Element in TaxonList:
			Found = 0
			for Line in Sequenzfile:
				if Element in Line:
					SequenceList[TaxonNum] += read_sequence(Sequenzfile)
					#SequenceList[TaxonNum] += Sequenzfile.next().strip("\n") #works only for sequences in single lines!!!!
					TaxonListOutput[TaxonNum] += "X"
					Found = 1
					break					
				LineIndex += 1
			if Found == 0:
				add_to_taxon(TaxonNum)
			LineIndex = 0
			Sequenzfile.seek(0)	
			TaxonNum +=1
		SequenceList = add_missing(SequenceList)
	
	
		
	if file == 0:
		print "Concat: No \"*_replaced\" files found in specified directory", TMPD
		return

	for i in range(0, TotTax): #remove last 30 NNs (this needs work!)
		SequenceList[i] = SequenceList[i][:-30]
		
		
	# write output
	Outfile = open("concat.fas", "w")
	for i in range(0, TotTax):
		#print ">" + TaxonListOutput[i]
		Outfile.write(">" + TaxonListOutput[i] + "\n")
		#print SequenceList[i]
		Outfile.write(SequenceList[i]+"\n")
		

#-------------------------------------------------------------- MAIN	
if Args.reduce == True:
	print "-r specified: Reducing taxa"
	reduce()
if Args.align == True:
	print "-a specified: Aligning taxa"
	align()
if Args.replace == True:
	print "-p specified: Replacing gaps"
	replace()
if Args.concat == True:
	print "-c specified: Concatenate"
	concat()
if Args.reduce == False and Args.align == False and Args.replace == False and Args.concat == False:
	print "nothing specified: Will reduce, replace, align and concat..."
	reduce()
	align()
	replace()
	concat()
if Args.clean == True:
	print "cleaning up..."
	for DelFile in glob.glob(TMPD+"*_reduced*"):
		print DelFile
		os.remove(DelFile)
print "done"
			