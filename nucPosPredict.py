import numpy as np;
import time
import scipy.io
import sys
import math

def getNucleosomeRegions(data, out):
	arr = []
	Enriched_Regions = []
	Depleted_Regions = []

	# Each row is a chromosome, elements are tuples of start and end positions
	for i in range(16):
		Enriched_Regions.append([])
		Depleted_Regions.append([])
	string = ""
	j = 0;

	# nucleosome forming / depleting threshold
	threshold_enriched = 0.75;
	threshold_depleted = -0.75;

	e_count = 0;
	curr_chr = 0;
	curr_pos = 0;
	# e_count_loc = 0;
	d_count = 0;
	consec_length = 50;

	for lines in data:
		line=lines.split("\n")[0].split('\t')
		# arr.append([])
		if line[0] == "Chromosome":
			print("Skipped Header")
			continue;
		# position = line[0]+":"+(line[1)]
		# arr.append({'pos': position, 'score':float(line[2])})
		if int(line[0]) != curr_chr:
			print("Processing Chr ");
			print(line[0]);
			e_count = 0;
			e_start = 0;
			e_end = 0;
			d_count = 0;
			d_start = 0;
			d_end = 0;
			curr_chr = int(line[0]);

		# Enriched
		if float(line[2]) >= threshold_enriched:
			if e_count == 0:
				e_start = int(line[1]);
			e_count += 1;
			
			if d_count >= consec_length:
				d_end = int(line[1]);
				Depleted_Regions[curr_chr-1].append((d_start, d_end));
			d_count = 0;

		# Neither
		elif float(line[2]) < threshold_enriched and float(line[2]) > threshold_depleted:

			# Found enriched region
			if e_count >= consec_length:
				e_end = int(line[1]);
				Enriched_Regions[curr_chr-1].append((e_start, e_end));
			e_count = 0;

			# Found depleted region
			if d_count >= consec_length:
				d_end = int(line[1]);
				Depleted_Regions[curr_chr-1].append((d_start, d_end));
			d_count = 0;

		# Depleted
		elif float(line[2]) <= threshold_depleted:
			if e_count >= consec_length:
				e_end = int(line[1]);
				Enriched_Regions[curr_chr-1].append((e_start, e_end));
			e_count = 0;

			if d_count == 0:
				d_start = int(line[1]);
			d_count += 1;

	print("Fetching Enriched/Depleted Region Seqeuences \n");
	string += "Enriched_Regions \n";
	pt("Enriched_Regions", 0);
	for chromosome in range(16):

		string += "chromosome: \n";
		pt("chromosome: ", 0);
		string += str(chromosome+1) + "\n";
		pt(chromosome+1, 0);

		for reg in Enriched_Regions[chromosome]:
			reg_length = reg[1] - reg[0] + 1;
			seq = fetchSeqence(chromosome, reg[0], reg[1]);
			pt(f'{reg} {reg_length} {seq}', 0);
			string += str(reg) + " " + str(reg_length) + " " + str(seq) + "\n";

			

	string += "Depleted_Regions \n";
	pt("Depleted_Regions", 0)
	for chromosome in range(16):

		string += "chromosome: \n";
		pt("chromosome: ", 0);
		string += str(chromosome+1) + "\n";
		pt(chromosome+1, 0);

		for reg in Depleted_Regions[chromosome]:
			reg_length = reg[1] - reg[0] + 1;
			seq = fetchSeqence(chromosome, reg[0], reg[1]);
			pt(f'{reg} {reg_length} {seq}', 0);
			string += str(reg) + " " + str(reg_length) + " " + str(seq) + "\n";
	pt(string, 0)
	# write out the sequence in fasta format
	f = open(out, 'w+');
	f.write(string)
	f.write("\n")
	f.close()
	return Enriched_Regions, Depleted_Regions;

def pt(content, debug):
	if debug == 1:
		print(content)

def fetchSeqence(chromosome, start, end):
	from Bio import Entrez, SeqIO
	GI = ["330443391",
	"330443482",
	"330443489",
	"330443520",
	"330443531",
	"330443543",
	"330443578",
	"330443590",
	"330443595",
	"330443638",
	"330443667",
	"330443681",
	"330443688",
	"330443715",
	"330443743",
	"330443753"]
	Entrez.email = "yibeijia@usc.edu"
	handle = Entrez.efetch(db="nucleotide", 
		id=GI[chromosome], 
		rettype="fasta", 
		strand=1, 
		seq_start=start, 
		seq_stop=end)
	record = SeqIO.read(handle, "fasta")
	handle.close()
	return record.seq;

def readInputAsArray(fileName):
	with open(fileName, 'r') as myfile:
		data = myfile.readlines()
	# print(data)
	return data

def readInputAsString(fileName):
	with open(fileName, 'r') as myfile:
		data=myfile.read().replace('\n', ' ')
		# print 
	return data

# 147 rows X 4 columns
def oneHotEncode(seq):
    import numpy as np
    seq2=list()
    mapping = {"A":[1., 0., 0., 0.], "C": [0., 1., 0., 0.], "G": [0., 0., 1., 0.], "T":[0., 0., 0., 1.]};
    for i in seq:
    	seq2.append(mapping[i]  if i in mapping.keys() else [0., 0., 0., 0.]);
    return seq2;
    # return np.array(seq2);


def encodeNucSeq(data, total_sections, section):

	enrichSeqCount = 0;
	depleteSeqCount = 0;
	depleteLineStart = 0;
	enrich = True;

	# Count # of enriched sequences
	total_lines = 0;
	total_sections = int(total_sections);
	section = int(section);
	for lines in data:
		line=lines.split("\n")[0].split(' ')
		total_lines += 1;
		if line[0] == "Enriched_Regions":
			# print("in enriched region\n")
			enrich = True;
			
		if line[0] == "Depleted_Regions":
			# print("in depleted region\n")
			enrich = False;
			depleteLineStart = total_lines;
		if line[0] == "":
			continue;
		if line[0][0] == "(":
			if enrich:
				enrichSeqCount += 1;
			else:
				depleteSeqCount += 1;
	print(f"Total # of regions is: {total_sections}\n");			
	print(f"Current regions is: {section}\n");			
	print(f"Total # of enriched sequences is: {enrichSeqCount}\n");
	print(f"Total # of Depleted sequences is: {depleteSeqCount}\n");
	print(f"Total # of lines is: {total_lines}\n");

	# one hot encode the DNA
	allEnrichSeqArr=[];
	allDepleteSeqArr=[];
	enrich = True;
	
	tic = time.perf_counter();
	section_length = math.floor(total_lines / total_sections);
	start_line = (section-1) * section_length;
	end_line = section * section_length-1;
	
	if (end_line > total_lines):
		end_line = total_lines-1;
	print(f"section_length is: {section_length}\n")
	print(f"start_line is: {start_line}\n")
	print(f"end_line is: {end_line}\n")
	print(f"depleteLineStart is: {depleteLineStart}\n")
	
	line_count = 0;
	for lines in data:
		# debug
		# print(f"line_count is: {line_count}\n")
		if(line_count < start_line):
			line_count += 1;
			continue;
		if(line_count > end_line):
			# print("larger")
			break;
		line_count += 1;
		# print(f"line_count is: {line_count}\n")
		# debug
		line=lines.split("\n")[0].split(' ')

		if line_count < depleteLineStart:
			# print("in enriched region\n")
			enrich = True;
			
		if line_count >= depleteLineStart:
			# print("in depleted region\n")
			enrich = False;
		if line[0] == "":
			continue;

		if line[0][0] == "(":
			curr_seq_length = int(line[2]);
			curr_seq =line[3];
			encodedDNAArr = [];

			if enrich:
				## for seqeunces that are shorter
				if(len(curr_seq) < 147):	
					while len(curr_seq) < 147:
						curr_seq = curr_seq + 'N';
					encodedDNAArr = oneHotEncode(curr_seq);
					allEnrichSeqArr+=[encodedDNAArr];
				else:
					## for seqeunces that are long enough	
					for start in range(curr_seq_length-146):
						curr_start = start;
						curr_end = 146 + start;
						curr_seq = line[3][curr_start: curr_end+1];
						encodedDNAArr = oneHotEncode(curr_seq);
						allEnrichSeqArr+=[encodedDNAArr];

			else:
				## for seqeunces that are shorter
				if(len(curr_seq) < 147):	
					while len(curr_seq) < 147:
						curr_seq = curr_seq + 'N';
					encodedDNAArr = oneHotEncode(curr_seq);
					allDepleteSeqArr+=[encodedDNAArr];
				else:
					## for seqeunces that are long enough	
					for start in range(curr_seq_length-146):
						curr_start = start;
						curr_end = 146 + start;
						curr_seq = line[3][curr_start: curr_end+1];
						encodedDNAArr = oneHotEncode(curr_seq);
						allDepleteSeqArr += [encodedDNAArr];
		else:
			print("header")
	# print(allEnrichSeqArr);
	toc = time.perf_counter();
	print(f"Getting the DNA encoded took {toc - tic:0.4f} seconds");
	nrow_enrich = len(allEnrichSeqArr);
	# ncol_enrich = -1;
	# ncol_enrich2 = -1;
	# if (len(allEnrichSeqArr) > 0):
	# 	ncol_enrich = len(allEnrichSeqArr[0]);
	# 	if (len(allEnrichSeqArr[0]) > 0):
	# 		ncol_enrich2 = len(allEnrichSeqArr[0][0]);
	# print(f" allEnrichSeqArr.shape: {nrow_enrich} x {ncol_enrich} x {ncol_enrich2}");

	nrow_deplete = len(allDepleteSeqArr)
	# ncol_deplete = -1;
	# ncol_deplete2 = -1;
	# if (len(allDepleteSeqArr) > 0):
	# 	ncol_deplete = len(allDepleteSeqArr[0])
	# 	if (len(allDepleteSeqArr[0]) > 0):
	# 		ncol_deplete2 = len(allDepleteSeqArr[0][0]);
	# print(f" allDepleteSeqArr.shape: {nrow_deplete} x {ncol_deplete} x {ncol_deplete2}");

	# split files if they get too large
	if (nrow_enrich > 100000000):
		num_sub_array = 10;
		sub_size = math.floor(nrow_enrich / num_sub_array);
		for index in range(0, num_sub_array):
			index_str = index;
			sub_arr_start = sub_size*index;
			sub_arr_end = sub_size*(index+1);
			if (sub_arr_end > nrow_enrich):
				sub_arr_end = nrow_enrich;
			np_allEnrichSeqArr = np.array(allEnrichSeqArr[sub_arr_start:sub_arr_end])
			np_allDepleteSeqArr = np.array(allDepleteSeqArr)
			# print(np_allDepleteSeqArr)
			print(f"1 np_allEnrichSeqArr.shape {np_allEnrichSeqArr.shape}, np_allDepleteSeqArr.shape {np_allDepleteSeqArr.shape}");
			Data = {"EnrichedData": np_allEnrichSeqArr, "DepletedData": np_allDepleteSeqArr};
			# mat_filename = '/scratch2/yibeijia/data/nucleosome_occupancy_' + str(section) + '_' + str(index_str) + '.mat';
			mat_filename = '/Users/yibeijia/Downloads/nucleosome_occupancy/train_test_data/nucleosome_occupancy_' + str(section) + '_' + str(index_str) + '.mat';
			scipy.io.savemat(mat_filename, Data,  do_compression=True);

	# split files if they get too large
	elif (nrow_deplete > 100000000):
		num_sub_array = 20;
		sub_size = math.floor(nrow_deplete / num_sub_array);
		for index in range(0, num_sub_array):
			index_str = index;
			sub_arr_start = sub_size*index;
			sub_arr_end = sub_size*(index+1);
			if (sub_arr_end > nrow_deplete):
				sub_arr_end = nrow_enrich;
			np_allEnrichSeqArr = np.array(allEnrichSeqArr)
			np_allDepleteSeqArr = np.array(allDepleteSeqArr[sub_arr_start:sub_arr_end])
			# print(np_allDepleteSeqArr)
			print(f"2 np_allEnrichSeqArr.shape {np_allEnrichSeqArr.shape}, np_allDepleteSeqArr.shape {np_allDepleteSeqArr.shape}");
			Data = {"EnrichedData": np_allEnrichSeqArr, "DepletedData": np_allDepleteSeqArr};

			# mat_filename = '/scratch2/yibeijia/data/nucleosome_occupancy_' + str(section) + '_' + str(index_str) + '.mat';
			mat_filename = '/Users/yibeijia/Downloads/nucleosome_occupancy/train_test_data/nucleosome_occupancy_' + str(section) + '_' + str(index_str) + '.mat';
			scipy.io.savemat(mat_filename, Data,  do_compression=True);

	else:
		np_allEnrichSeqArr = np.array(allEnrichSeqArr)
		np_allDepleteSeqArr = np.array(allDepleteSeqArr)
		print(f"3 np_allEnrichSeqArr.shape {np_allEnrichSeqArr.shape}, np_allDepleteSeqArr.shape {np_allDepleteSeqArr.shape}");

		Data = {"EnrichedData": np_allEnrichSeqArr, "DepletedData": np_allDepleteSeqArr};
		# mat_filename = '/scratch2/yibeijia/data/nucleosome_occupancy_' + str(section) + '.mat';
		mat_filename = '/Users/yibeijia/Downloads/nucleosome_occupancy/train_test_data/nucleosome_occupancy_' + str(section) + '.mat';
		scipy.io.savemat(mat_filename, Data,  do_compression=True);


def main():
	##### ########## ########## ########## ########## #####
	#Get the sequencnes in depleted ir enriched regions
	# data = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/GSE13622_RAW/normalizedMap/GSM351491_InVitro_normalized.tab')
	# print(getNucleosomeRegions(data,'regions_out.txt'))
	# tic = time.perf_counter();
	# getNucleosomeRegions(data,'InVitro_regions_out.txt');
	# toc = time.perf_counter();
	# print(f"Getting the nucleosome regions took {toc - tic:0.4f} seconds");

	# getNucleosomeSeq(Enriched_Regions)
	# getNucleosomeSeq(Enriched_Regions)

	##### ########## ########## ########## ########## #####
	##### Get the sequencnes in depleted ir enriched regions
	# data = readInputAsArray('/project/rohs_108/yibeijia/nucleosome_occupancy/InVitro_regions_out.txt');
	data = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/InVitro_regions_out.txt')
	encodeNucSeq(data, sys.argv[1], sys.argv[2])
	# dna="ACGTAC";
	# encodedDna=oneHotEncode(dna)
	# print("dna\n",list(dna))
	# print('encoded dna\n')
main();
