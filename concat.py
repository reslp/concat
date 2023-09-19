#!/usr/bin/env python
# script to produce concatenated alignments and different alignment statistics
# originally created 31.07.2014, extended since then
# written by Philipp Resl

# known issues:
# can't change aligner settings - in align
# names in ID file and sequence files have to be identical

import argparse
import glob
import os
import platform
import time
from subprocess import Popen, PIPE
import sys
from collections import defaultdict
import copy

ver = "0.32"
parser = argparse.ArgumentParser(description="%(prog)s (build: Aug 2, 2023) will create concatened alignments from seperate single locus files and produce different kinds of alignment statistics.")
parser.add_argument("-t", dest="taxonfile", help="specify place of file containing desired taxon set")
parser.add_argument("-d", dest="directory", help="specify directory with single locus files that should be used")
parser.add_argument("-i", dest="input_file", nargs="*", type=str, help="names & path to input files that should be used. Seperate filenames with space. e.g.: \"its.fas lsu.fas\"")
parser.add_argument("-o", dest="o", help="specify output directory")
parser.add_argument("-l", dest="minlen", default = 1, type=int, help="Minimum length of sequence to keep (used in reduce). Default: 1")
parser.add_argument("-v", action="version", help="outputs version of concat script", version="%(prog)s Version "+str(ver)+" (08.09.20)")
#parser.add_argument("--add-to-name", dest="add", action="store_true", help="Append X to species name if sequences is present. O if it is missing")
parser.add_argument("-m", "--runmode", dest="runmode",action="store", help="Specify runmode. Possible options: replace, reduce, align, concat, all. Runmode all will ignore -o")
parser.add_argument("-N", dest="NNN",action="store_true", help="adds Ns between loci to seperate them in concatenated alignment. Default: false")
parser.add_argument("-M", dest="MMM",default="?", help="Character to add for missing sequences (concat) and positions at beginning and end of sequences (replace). Default: ?")
parser.add_argument("--biopython", dest="biopython",action="store_true", default=False, help="Optional: Use approach based on biopython. biopython v1.77 needs to be installed. Currently works only for --runmode concat.")
parser.add_argument("--statistics", dest="partition",action="store_true", default=False, help="Optional: Output start and end position of partitions in concatenated alignment. Needs --biopython")
parser.add_argument("--seqtype", dest="seqtype", default="nu", help="Specify sequence type. Possible values: nu and aa. default = nu. This is needed to correctly identify parsimony informative sites.")
parser.add_argument("--noseq", dest="noseq", action="store_true", default=False, help="Can be used in combination with --statistics. When this is specified only the statistics file but NO sequence file will be produced.")

#parser.add_argument("-aligner", dest="align_param",action="store_true", help="Arguments that should be sent to the aligner (not yet implemented)")
#parser.add_argument("-clean", dest="clean",action="store_true", help="performs a clean run: removes all intermediate files")
#parser.add_argument("-trim", dest="trim",action="store_true", help="reduce alignment to include only sites in at least number of sequences (not yet implemented)")
Args = parser.parse_args()

def now():
	return time.strftime("%Y-%m-%d %H:%M") + " -"

# very simple check if file is in fasta format
def check_fasta(filename):
	if not os.path.isabs(filename):
		filename = os.path.abspath(filename)
	if os.stat(filename).st_size == 0:
		print(now(), "File", filename, "is empty. This file will be skipped.")
		return None
	f = open(filename, "r")
	if not f.readline().startswith(">"):
		print(now(), "File", filename, "is probably not in FASTA format. This file will be skipped.")
		return None
	f.close()
	return filename
		
	
def get_input_files(which_type, file_info):
	if which_type == "files":
		SeqFileList = file_info
		OnlyFileName = [check_fasta(f) for f in SeqFileList if check_fasta(f) != None] 
		return sorted(OnlyFileName, key=lambda i: os.path.splitext(os.path.basename(i))[0])
	if which_type == "dir":
		OnlyFileName = [check_fasta(f) for f in glob.glob(file_info+"/*") if check_fasta(f) != None]
		return sorted(OnlyFileName, key=lambda i: os.path.splitext(os.path.basename(i))[0])


def reduce(taxon_list, file_list, outdir, WD): #reduces alignments to desired set of taxa
	if outdir == None or outdir == "":
		print(now(), "(Reduce) No output directory specified. Will create output files in folder reduce.", file=sys.stderr)
		outdir = "reduce"
	path = os.path.join(WD, outdir)
	if not os.path.exists(path):
		os.makedirs(path)
	TMPD = path
	TaxonList = taxon_list
	file=0
	TotTax = len(TaxonList)
	for filename in file_list:
		file=1
		Filename = os.path.join(TMPD, (os.path.basename(filename)))
		Outfile = open(Filename,"w")
		Sequenzfile = open(os.path.join(WD,filename), "r")
		#create a blank Sequencelist for given number of taxa
		SequenceList = [""]*len(TaxonList)
		#Open and read Files
		Found = 0
		for Taxon in TaxonList:
			Sequence = ""
			SeqName = ""
			for Line in Sequenzfile:
				if Line.startswith(">"+Taxon):
					SeqName = Line.strip("\n")
					Found = 1
					continue
				if Found == 1:
					if Line.find(">") == -1:
						Sequence += Line.strip("\n")
					else:
						Found = 0
			if len(Sequence) > Args.minlen:
				Outfile.write(SeqName+"\n")
				Outfile.write(Sequence+"\n")
			else:
				if SeqName != "":
					print (now(), "(reduce) %s: Sequence %s is too short (<200 positions) and was not included" % (filename,SeqName), file=sys.stderr)
			Sequenzfile.seek(0)
		Outfile.close()
	if file == 0:
		print (now(), "(reduce) No files provided.")

def align(files_list, outdir, WD, OS): # aligns sequence files
	if outdir == None or outdir == "":
		print(now(), "(align) No output directory specified. Will place output files in folder align.", file=sys.stderr)
		outdir = "align"
	path = os.path.join(WD, outdir)
	if not os.path.exists(path):
		os.makedirs(path)
	file = 0
	if OS == "Darwin" or OS == "Linux":
		for RedFile in files_list:
			if os.stat(RedFile).st_size == 0:
				print(now(), "(align) File seems to be empty:", RedFile, ". Will skip.", file=sys.stderr)
			file = 1
			try:
				print(now(), "(align) Aligning file with mafft:", RedFile, file=sys.stderr)
				process = Popen(["mafft","--auto",RedFile], stdout=PIPE, stderr=PIPE)
				(output, err) = process.communicate()
				exit_code = process.wait()
			except:
				print (now(), "(align) mafft executable not found. Is it in your path?", file=sys.stderr)
				sys.exit(1)
			Outfile = open(os.path.join(path, os.path.basename(RedFile)),"wb")
			Outfile.write(output)
			Outfile.close()
		if file == 0:
			print (now(), "(align) no files given.", file=sys.stderr)
		return


def replace(file_list, outdir, WD):	# replaces - with ? at the beginning and end of fasta alignments
	if outdir == None or outdir == "":
		print(now(), "(replace) No output directory specified. Will place output files in folder replace.", file=sys.stderr)
		outdir = "replace"
	path = os.path.join(WD, outdir)
	if not os.path.exists(path):
		os.makedirs(path)
	def replace_begin (Sequenz_0):
		zahl = 0
		for Base in Sequenz_0:
			if Base == "-":
				Sequenz_0[zahl] = Args.MMM
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
				Sequenz_0[len(Sequenz_0) + zahl2] = Args.MMM
			else:
				break
			zahl2 -= 1
		return "".join(Sequenz_0)
	#end replace_end
	file = 0
	for AlFile in file_list:
		file = 1
		TaxonList = []
		File = open(AlFile, "r")
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

		Outfile = open(os.path.join(WD, outdir, os.path.basename(AlFile)), "w")
		for Taxon in range(0,MaxSeq):
			#print TaxonList[Taxon]
			if len(SequenceList[Taxon])-SequenceList[Taxon].count("?")-SequenceList[Taxon].count("-") >= Args.minlen:
				Outfile.write(TaxonList[Taxon]+"\n")
				#print "seqlist:"
				#print SequenceList
				Outfile.write(SequenceList[Taxon])
				Outfile.write("\n")
			else:
				print (now(), "(replace) %s was not added because it is too short." % TaxonList[Taxon])
		Outfile.close()
		File.close()
	if file == 0:
		print (now(), "(replace) No files provided.")


def add_missing(seqdict):
	longest = max(seqdict, key=lambda k: len(seqdict[k]))
	length = len(seqdict[longest]) 
	#print(length)
	for taxon in seqdict.keys():
		if len(seqdict[taxon]) < length:
			seqdict[taxon] += Args.MMM * (length-len(seqdict[taxon]))
		if Args.NNN == True:
			if Args.seqtype == "nu":
				seqdict[taxon] += "N" * 5
			if Args.seqtype == "aa":
				seqdict[taxon] += "X" * 5
	return seqdict


def get_all_taxa(file_list):
	taxon_list = []
	print(now(), "(get_all_taxa): No IDs specified, will use sequence names.", file= sys.stderr)
	for seqfile in file_list:
		taxlist = []
		for line in open(seqfile, "r"):
			if line.startswith(">"):
				current_taxon = line.strip().split(">")[-1]
			else:
				continue
			taxlist.append(current_taxon)
		for taxon in set(taxlist):
			if taxlist.count(taxon) > 1:
				print(now(), "(get_all_taxa):WARNING: Duplicated taxon: ",taxon," in file ",seqfile,  file= sys.stderr)
			taxon_list += list(set(taxlist))
	return taxon_list

def concat(taxon_list, file_list, outdir, WD):
	if not file_list:
		print(now(), "(concat): No files specified", file= sys.stderr)
		return
	if taxon_list == None:
		taxon_list=get_all_taxa(file_list)
	if outdir == None or outdir == "":
		print(now(), "(concat) No output directory specified. Will place output files in wd.", file=sys.stderr)
		outdir = ""
	path = os.path.join(WD, outdir)
	if not os.path.exists(path):
		os.makedirs(path)
	print(now(), "(concat) Output is set to", path, file=sys.stderr)
	TaxonList = taxon_list
	TaxonListOutput = TaxonList [:] #Taxon List for Output
	print(now(), "(concat): Concatenating", len(file_list), "files", file= sys.stderr)
	concat_dict = dict.fromkeys(taxon_list, "")
	for seqfile in file_list:
		#print(seqfile)
		current_taxon = ""
		for line in open(seqfile, "r"):
			if line.startswith(">"):
				current_taxon = line.strip().split(">")[-1]
			else:
				concat_dict[current_taxon] += line.strip("\n")
		concat_dict = add_missing(concat_dict)
	Outfile = open(os.path.join(path,"concat.fas"), "w")
	for taxon in taxon_list:
		Outfile.write(">" + taxon + "\n")
		Outfile.write(concat_dict[taxon]+"\n")
		
def bio_concat(taxon_list, file_list, outdir, WD):
	# this currently only works with fasta files!
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import UnknownSeq, Seq
	from Bio.Align import MultipleSeqAlignment 
	if outdir == None or outdir == "":
		print(now(), "(concat) No output directory specified. Will place output files in wd.", file=sys.stderr)
		outdir = ""
	path = os.path.join(WD, outdir)
	if not os.path.exists(path):
		os.makedirs(path)
	print(now(), "(concat) Output is set to", path, file=sys.stderr)
	if taxon_list == None:
		taxon_list=get_all_taxa(file_list)
	concat_dict = dict.fromkeys(taxon_list, "")
	alignments = []
	for seqfile in file_list:
		alignment = []
		for record in SeqIO.parse(seqfile, "fasta"):
			if record.id in taxon_list:
				alignment.append(record)
		#print(len(alignment))
		alignments.append(MultipleSeqAlignment(alignment))
	new_alignments =[] 
	position = 0
	
	#create temp dictionary
	new_alignment = defaultdict(list)
	alignment_infos = []
	for alignment,seqfile in zip(alignments,file_list):
		
		length = alignment.get_alignment_length()
		#gather alignment information if --statistics flag is specified:
		if Args.partition:
			nspecies = len(alignment)
			nparsimony = get_parsimony_sites(alignment)
			nvariable = get_variable_sites(alignment)
			nfixed = get_fixed_sites(alignment)
			rcv = get_rcv_value(alignment)
			alignment_infos.append((seqfile, position + 1, position + length, length, nspecies, nparsimony, nvariable, nfixed, rcv))
			position += length
		
		labels_in_alignment = set(seq.id for seq in alignment)
		# get sequence names missing in alignment:
		missing = set(taxon_list) - labels_in_alignment

		# now create sequences for missing labels:
		for seqid in missing:
			new_seq = UnknownSeq(length, alphabet=alignment._alphabet)
			new_alignment[seqid].append(str(new_seq))
		# also add all other sequences to temporary dict
		for seq in alignment:
			new_alignment[seq.id].append(str(seq.seq))
	
	if not Args.noseq:
		Outfile = open(os.path.join(path,"concat.fas"), "w")
		msa = MultipleSeqAlignment(SeqRecord(Seq(''.join(sequence), alphabet=alignments[0]._alphabet), id=seqid) for (seqid, sequence) in new_alignment.items())
		for seq in msa:
			Outfile.write(">" + seq.id + "\n")
			Outfile.write(str(seq.seq) + "\n")
		Outfile.close()
	
	if Args.partition:
		Outfile = open(os.path.join(path,"statistics.txt"), "w")
		print(now(), "(concat) --statistics specified. Will create alignment statistics file.", file=sys.stderr)
		Outfile.write("alignment\tstart\tend\tlength\tnseqs\tnparsimony\tnvariable\tnfixed\trcv\n")
		for info in alignment_infos:
			out = str(info[0].split("/")[-1])+"\t"+str(info[1])+"\t"+str(info[2])+"\t"+str(info[3])+"\t"+str(info[4])+"\t"+str(info[5])+"\t"+str(info[6])+"\t"+str(info[7])+"\t"+str(info[8])+"\n"
			Outfile.write(out)
		Outfile.close()
		#print(alignment_infos)

def get_rcv_value(alignment):
	from Bio.Align import MultipleSeqAlignment
	from Bio.Seq import Seq
	# calculate rcv (relative composition variability) according to
	# Phillips & Penny (2003) - Molecular phylogenetics and Evolution - https://doi.org/10.1016/S1055-7903(03)00057-5
	# remove fixed positions from alignment:
	pos_to_keep = []
	for i in range(0, alignment.get_alignment_length()):
		if len(set(alignment[:, i])) == 1: # check if alignment column contains a fixed letter
			next
		else:
			pos_to_keep.append(i)	
	shortened_records_list = []
	for record in alignment:
		tmpseq = "".join([record.seq[i] for i in pos_to_keep])
		tmprecord = copy.copy(record) #create a shallow copy of alignment for sanity.
		tmprecord.seq = Seq(tmpseq, alignment._alphabet)
		shortened_records_list.append(tmprecord)
	shortalignment = MultipleSeqAlignment(shortened_records_list)

	ntaxa = len(shortalignment)

	# calculate averages for each occuring letter:
	combinedseq = ""
	for seq in shortalignment:
		combinedseq += str(seq.seq)
	if Args.seqtype == "aa": # remove gaps and ambiguous positions
		combinedseq = combinedseq.replace("-", "").replace("?","").replace("X", "")
	else:
		combinedseq = combinedseq.replace("-", "").replace("?","").replace("N", "")
	av = {}	
	for letter in set(combinedseq):
		av[letter] = combinedseq.count(letter) / ntaxa	
	
	# calculate rcv values and combine them
	rcv_values = []
	for seq in shortalignment:
		rcv = 0
		for letter in av.keys():
			rcv += abs(seq.seq.count(letter) - av[letter])
		rcv_values.append(rcv)

	if (shortalignment.get_alignment_length() > 0): #check contains variable sites to calculate rcv
		final_rcv = round(sum(rcv_values) / (ntaxa * shortalignment.get_alignment_length()), 4)
	else: # if not set rcv to 1
		final_rcv = 1

	return final_rcv

def get_parsimony_sites(alignment):
	# calculate parsimony informative sites:
	# considered sites have at least two states with at least two occurences.
	npars = 0
	for i in range(0, alignment.get_alignment_length()):
		unique_chars = set(alignment[:, i])
		#skip positions which contain ambiguous positions.
		#handle DNA and AA sequences differently here:
		if Args.seqtype == "aa":
			if "X" in unique_chars or "-" in unique_chars or "?" in unique_chars:
				continue
		else:
			if "N" in unique_chars or "-" in unique_chars or "?" in unique_chars:
				continue
		counts = [pos for pos in unique_chars if alignment[:, i].count(pos) >=2]
		ncounts = len(counts)
		if ncounts >= 2:
			npars += 1
	return npars	
	
def get_variable_sites(alignment):	
	# number of variable sites
	nvar = 0
	for i in range(0, alignment.get_alignment_length()):
		if len(set(alignment[:, i])) >= 2:
			nvar += 1
	return nvar
	
def get_fixed_sites(alignment):
	# number of fixed sites
	nfixed = 0
	for i in range(0, alignment.get_alignment_length()):
		unique_chars = set(alignment[:, i])
		#skip positions which contain ambiguous positions.
		#handle DNA and AA sequences differently here:
		if Args.seqtype == "aa":
			if "X" in unique_chars or "-" in unique_chars or "?" in unique_chars:
				continue
		else:
			if "N" in unique_chars or "-" in unique_chars or "?" in unique_chars:
				continue
		if len(unique_chars) == 1:
			nfixed += 1
	return nfixed
		

def check_sequence_input(input_i, input_d):
	if input_i == None and input_d == None:
		print(now(), "Need input files. Use either -i oder -d. Try also -h for more information.", file=sys.stderr)
		return
	if input_i != None and input_d != None:
                print (now(), "Can't use both input files (-i flag) and directory (-d flag). Please use only one option.", file=sys.stderr)
                return
	print(now(), "Checking sequence input files...", file=sys.stderr)
	if input_i != None:
                which_type = "files"
                file_info = input_i
	elif input_d != None:
                which_type = "dir"
                file_info = input_d
        #print(which_type)
	input_file_list = get_input_files(which_type, file_info)
	print(now(), "Found", len(input_file_list), "alignment files in FASTA format")
	return input_file_list

def check_taxid_file(taxonfile):
	if taxonfile == None:
		return
	print(now(), "Reading taxon file...", file=sys.stderr)
	TaxonFile = open(taxonfile, "r")
	taxon_list = []
	for Line in TaxonFile:
		if Line.strip("\n") != "": #only append to list of line is not empty
			taxon_list.append(Line.strip("\n"))
	TaxonFile.close()
	if len(taxon_list) == 0:
		print(now(), "Taxon file", TaxonFile, "is empty.", file=sys.stderr)
		sys.exit(1)
	print(now(), "Taxon file containes", len(taxon_list), "taxa", file=sys.stderr)

	return taxon_list 

if __name__ == "__main__":
	print(now(), "---- Welcome to concat v%s ----" % ver, file=sys.stderr)
	if len(sys.argv)<2:
		parser.print_help()
		sys.exit(1)
	WD = os.getcwd()
	WD = str(WD)
	OS = platform.system()
	print (now(), "Running on platform:",OS, file=sys.stderr)
	print (now(), "Working directory: ", WD, file=sys.stderr)
	print (now(), "Assuming", Args.seqtype, "sequences.", file=sys.stderr)
	
	outdir = Args.o
	
	if Args.runmode == "reduce":
		print (now(), "Runmode: reduce. Reducing sequence files to specified taxa...", file=sys.stderr)
		input_file_list = check_sequence_input(Args.input_file, Args.directory)
		taxon_list = check_taxid_file(Args.taxonfile)
		if input_file_list and taxon_list:
			reduce(taxon_list, input_file_list, outdir, WD)
		else:
			print(now(), "Taxon id file (-t) or sequence files (-i or -d) are missing. Please specify")
			sys.exit(1)
	if Args.runmode == "align":
		print (now(), "Runmode: aligning. Align sequences...", file=sys.stderr)
		input_file_list = check_sequence_input(Args.input_file, Args.directory)
		if input_file_list:
			align(input_file_list, outdir, WD, OS)
		else:
			print(now(), "Taxon id file (-t) or sequence files (-i or -d) are missing. Please specify")
			sys.exit(1)
	if Args.runmode == "replace":
		print (now(), "Runmode: replace. Replace gaps at ends of sequences...", file=sys.stderr)
		input_file_list = check_sequence_input(Args.input_file, Args.directory)
		if input_file_list:
			replace(taxon_list, input_file_list, outdir, WD)
		else:
			print(now(), "Taxon id file (-t) or sequence files (-i or -d) are missing. Please specify")
			sys.exit(1)
	if Args.runmode == "concat":
		print (now(), "Runmode: concat. Concatenating sequences...", file=sys.stderr)
		input_file_list = check_sequence_input(Args.input_file, Args.directory)
		if not input_file_list:
			print(now(), "No input files found. Check your specified -i or -d")
			sys.exit(1)
		taxon_list = check_taxid_file(Args.taxonfile)
		if Args.biopython:
			print (now(), "Using --biopython method", file=sys.stderr)
			if Args.noseq:
				print (now(), "--noseq specified. Will not write sequence file.", file=sys.stderr)
			bio_concat(taxon_list, input_file_list, outdir, WD)
		else:
			print (now(), "Using standard method", file=sys.stderr)
			concat(taxon_list, input_file_list, outdir, WD)
	if Args.runmode == "all":
		print (now(), "Runmode: all. Will run reduce, align, replace and concat.", file=sys.stderr)
		print (now(), "Running reduce...", file=sys.stderr)
		input_file_list = check_sequence_input(Args.input_file, Args.directory)
		taxon_list = check_taxid_file(Args.taxonfile)
		if input_file_list == None or taxon_list == None:
			sys.exit(1)
		reduce(taxon_list, input_file_list, "", WD)
		print (now(), "Running align...", file=sys.stderr)
		input_file_list = get_input_files("dir", "reduce")
		align(input_file_list, "align", WD, OS)
		print (now(), "Running replace...", file=sys.stderr)
		input_file_list = get_input_files("dir", "align")
		replace(input_file_list, "replace", WD)
		print (now(), "Running concat...", file=sys.stderr)
		input_file_list = get_input_files("dir", "replace")
		concat(taxon_list, input_file_list,None, WD)
	print (now(), "done", file=sys.stderr)
