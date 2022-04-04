import numpy as np;
import time
import scipy.io
import sys
import math
import re

# ####################### ####################### ####################### ######################
# This function outputs all the nucleosome positions for enriched/depleted regions
# Cutoff is determined by threshold_enriched and threshold_depleted
# Output sequence has minimum 50 bp in length
def processYeastData(data):
	Enriched_Regions = []
	Depleted_Regions = []
	num_chromsome=16;

	# nucleosome forming / depleting threshold
	threshold_enriched = 0.75;
	threshold_depleted = -0.75;

	e_count = 0;
	curr_chr = 0;
	curr_pos = 0;
	d_count = 0;

	# minimum consecutive length
	consec_length = 50;
	for i in range(num_chromsome):
		Enriched_Regions.append([])
		Depleted_Regions.append([])
		string = ""
		j = 0;
	for lines in data:
		line=lines.split("\n")[0].split('\t')
		if line[0] == "Chromosome":
			print("Skipped Header")
			continue;
		# Find Enriched/Depleted Region
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
	for chromosome in range(num_chromsome):

		string += "chromosome: \n";
		pt("chromosome: ", 0);
		string += str(chromosome+1) + "\n";
		pt(chromosome+1, 0);

		for reg in Enriched_Regions[chromosome]:
			reg_length = reg[1] - reg[0] + 1;
			seq = fetchSeqence(chromosome, reg[0], reg[1], "yeast");
			pt(f'{reg} {reg_length} {seq}', 0);
			string += str(reg) + " " + str(reg_length) + " " + str(seq) + "\n";
	string += "Depleted_Regions \n";
	pt("Depleted_Regions", 0)
	for chromosome in range(num_chromsome):

		string += "chromosome: \n";
		pt("chromosome: ", 0);
		string += str(chromosome+1) + "\n";
		pt(chromosome+1, 0);

		for reg in Depleted_Regions[chromosome]:
			reg_length = reg[1] - reg[0] + 1;
			seq = fetchSeqence(chromosome, reg[0], reg[1], "yeast");
			string += str(reg) + " " + str(reg_length) + " " + str(seq) + "\n";
	pt(string, 0)	
	return string

def processOtherData(data, species):
	print(f"species: {species}")
	# line format is [chr, s_pos, e_pos, seq]
	Depleted_Regions = []
	Enriched_Regions = []
	All_chromosomes=[]
	d_start=0
	d_end=0
	e_start=0
	e_end=0
	e_start_prev=0
	e_end_prev=0
	# minimum consecutive length
	consec_length = 50;

	start = True
	string = ""
	print("Fetching Enriched/Depleted Region Seqeuences \n");
	string += "Enriched_Regions \n";
	pt("Enriched_Regions", 0);

	for i in range(len(data)):
		lines = data[i]
		line=lines.split("\n")[0].split()
		chrom_start_end = line[0].split('_')
		chrom_start_end[0] = chrom_start_end[0].split('chr')[1]
		line=chrom_start_end + [line[1]]
		
		if line[0] not in All_chromosomes: # get all chrom from data
			string += "chromosome: \n";
			string += line[0] + "\n";

			All_chromosomes.append(line[0])
			if len(Enriched_Regions) >0:
				Enriched_Regions[len(Enriched_Regions)-1].sort(key=lambda y: y[0])
			Enriched_Regions.append([])

			Depleted_Regions.append([])
			print(All_chromosomes)
			pt(line[0], 0);

		e_start=int(float(line[1])) # avoid scientific notation
		e_end=int(float(line[2]))
		Enriched_Regions[len(All_chromosomes)-1].append((e_start, e_end));
		
		reg_length = e_end - e_start
		string += "(" + str(line[1]) + ", "+ str(line[2])+ ") " + str(reg_length) + " " + line[3] + "\n";

		if start:
			e_start_prev=e_start
			e_end_prev=e_end
			start = False
		else: # get depleted regions that are between enriched regions
			if e_start - e_end_prev > consec_length:
				d_start = e_end_prev + 1
				d_end = e_start - 1
				Depleted_Regions[len(All_chromosomes)-1].append((d_start, d_end));
			e_start_prev=e_start
			e_end_prev=e_end

# 	print(">>>>>>>>>>> Getting Depleted Seqs")
# 	string += "Depleted_Regions \n";
# 	for chromosome in range(len(All_chromosomes)):

# 		string += "chromosome: \n";
# 		pt("chromosome: ", 0);
# 		string += All_chromosomes[chromosome] + "\n";

# 		for reg in Depleted_Regions[chromosome]:
# 			reg_length = reg[1] - reg[0] + 1;
# 			seq = fetchSeqence(All_chromosomes[chromosome], reg[0], reg[1],species);
# 			string += str(reg) + " " + str(reg_length) + " " + str(seq) + "\n";

	return string;



def getNucleosomeRegions(data, out, species):
	
	num_chromsome=0
	# process yeast dataset
	if species == "yeast": 
		string= processYeastData(data,species)

	# process other dataset
	else: 
		string= processOtherData(data,species)
	# return;
	# write out the sequence in fasta format
	f = open(out, 'w+');
	f.write(string)
	f.write("\n")
	f.close()
	return ;

# ####################### ####################### ####################### ######################
# This function outputs all the nucleosome positions for neutral regions
# Cutoff is determined by threshold_enriched and threshold_depleted
# Output sequence has minimum 50 bp in length
def getNucleosomeNeutralRegions(data, out):
	arr = []
	Neutral_Regions = []

	# Each row is a chromosome, elements are tuples of start and end positions
	for i in range(16):
		Neutral_Regions.append([])
	string = ""
	j = 0;

	# cut off for neutral sequences
	threshold_enriched = 0.5;
	threshold_depleted = -0.5;

	curr_chr = 0;
	curr_pos = 0;

	# minimum consecutive length
	consec_length = 100;

	for lines in data:
		line=lines.split("\n")[0].split('\t')
		if line[0] == "Chromosome":
			print("Skipped Header")
			continue;
		# Find Neutral positioning sequences
		if int(line[0]) != curr_chr:
			print("Processing Chr ");
			print(line[0]);
			count = 0;
			start = 0;
			end = 0;
			curr_chr = int(line[0]);
		if float(line[2]) <= threshold_enriched and float(line[2]) >= threshold_depleted:

			# Found neutral region
			if count == 0:
				start = int(line[1]);
			count += 1;
		else:
			if count >= consec_length:
				end = int(line[1]);
				Neutral_Regions[curr_chr-1].append((start, end));
			count = 0;
	num_seq = 0;
	for chromosome in range(len(Neutral_Regions)):
		print(f"chromosome {chromosome+1}  has {len(Neutral_Regions[chromosome])} seqs")
		num_seq += len(Neutral_Regions[chromosome])

	print(f"# of Sequence in Neutral Region: {num_seq}")
	print("Fetching Neutral Region Seqeuences \n");
	string += "Neutral_Regions \n";
	for chromosome in range(16):
		print(f"Fetching for chromosome {chromosome+1} ")
		string += "chromosome: \n";
		string += str(chromosome+1) + "\n";
		for reg in Neutral_Regions[chromosome]:
			reg_length = reg[1] - reg[0] + 1;
			seq = fetchSeqence(chromosome, reg[0], reg[1]);
			string += str(reg) + " " + str(reg_length) + " " + str(seq) + "\n";

	# write out the sequence in fasta format
	f = open(out, 'w+');
	f.write(string)
	f.write("\n")
	f.close()
	return Neutral_Regions;

def pt(content, debug):
	if debug == 1:
		print(content)

#mdenv
def fetchSeqence(chromosome, start, end, species):
	from Bio import Entrez, SeqIO
	if species == "yeast":
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
	elif species == "worm":
		import pysam

		# print(f"fetching chr{chromosome}, {start}, {end}")

		fa_file="/Users/yibeijia/Downloads/data/chr" + chromosome+".fa"
		ref = pysam.FastaFile(fa_file)
		seq = ref.fetch('chr'+chromosome, start, end)
	return seq;

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


# ####################### ####################### ####################### ######################
# This function outputs all the one hot encoded nucleosome sequences
# And it outputs all sequences into a numpy array 
def encodeNucSeq(data, total_sections, section, species):

	enrichSeqCount = depleteSeqCount = neutralSeqCount = 0;
	depleteLineStart = -1
	neutralLineStart = -1;
	enrich = False;
	neutral = False;
	deplete = False;

	# Count # of enriched sequences
	total_lines = 0;
	total_sections = int(total_sections);
	section = int(section);
	for lines in data:
		line=lines.split("\n")[0].split(' ')
		total_lines += 1;
		if line[0] == "Enriched_Regions":
			print("in enriched region\n")
			enrich = True;

		elif line[0] == "Depleted_Regions":
			print("in depleted region\n")
			enrich = False;
			deplete = True
			depleteLineStart = total_lines;

		elif line[0] == "Neutral_Regions":
			print("in neutral region\n")
			neutral = True;
			enrich = False;
			neutralLineStart = total_lines;

		if line[0] == "":
			continue;
		if line[0][0] == "(":
			if enrich:
				enrichSeqCount += 1;
			elif deplete and enrich==False and neutral==False:
				depleteSeqCount += 1;
			else:
				neutralSeqCount += 1;
	print(f"Total # of regions is: {total_sections}\n");			
	print(f"Current regions is: {section}\n");			
	print(f"Total # of enriched sequences is: {enrichSeqCount}\n");
	print(f"Total # of Depleted sequences is: {depleteSeqCount}\n");
	print(f"Total # of Depleted sequences is: {neutralSeqCount}\n");
	print(f"Total # of lines is: {total_lines}\n");

	# one hot encode the DNA
	allEnrichSeqArr=[];
	allDepleteSeqArr=[];
	allNeutralSeqArr=[];
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
	print(f"neutralLineStart is: {neutralLineStart}\n")
	
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

		if line_count < depleteLineStart or depleteLineStart<0:
			enrich = True;
			
		if line_count >= depleteLineStart and depleteLineStart >0:
			enrich = False;
			deplete = True;

		if line_count >= neutralLineStart and neutralLineStart >0:
			deplete = False;	
			neutral = True;	

		if line[0] == "":
			continue;

		if line[0][0] == "(":
			curr_seq_length = int(line[2]);
			curr_seq =line[3];
			encodedDNAArr = [];

			## for seqeunces that are shorter
			if(len(curr_seq) < 147):	
				while len(curr_seq) < 147:
					curr_seq = curr_seq + 'N';
			## for seqeunces that are long enough	
			else:
				for start in range(curr_seq_length-146):
					curr_start = start;
					curr_end = 146 + start;
					curr_seq = line[3][curr_start: curr_end+1];

			if enrich:
				# print("in enrich")
				encodedDNAArr = oneHotEncode(curr_seq);
				allEnrichSeqArr+=[encodedDNAArr];

				
			elif deplete:
				# print("in deplete")
				encodedDNAArr = oneHotEncode(curr_seq);
				allDepleteSeqArr+=[encodedDNAArr];

			elif neutral:
				# print("in neutral")
				encodedDNAArr = oneHotEncode(curr_seq);
				allNeutralSeqArr+=[encodedDNAArr];
		else:
			print("header")
	# print(allEnrichSeqArr);
	toc = time.perf_counter();
	print(f"Getting the DNA encoded took {toc - tic:0.4f} seconds");
	print(len(allEnrichSeqArr))
	print(len(allEnrichSeqArr[0]))
	print(len(allEnrichSeqArr[0][0]))
	print(len(allEnrichSeqArr[1]))
	print(len(allEnrichSeqArr[1][0]))

	nrow_enrich = len(allEnrichSeqArr);
	nrow_deplete = len(allDepleteSeqArr)
	nrow_neutral = len(allNeutralSeqArr)
	np_allEnrichSeqArr= np.array([])
	np_allDepleteSeqArr= np.array([])
	np_allNeutralSeqArr= np.array([])
	print(f"nrow_enrich: {nrow_enrich}")
	print(f"nrow_deplete: {nrow_deplete}")
	print(f"nrow_neutral: {nrow_neutral}")

	# split files if they get too large
	if nrow_enrich!=0:
		print("creating enriched arr")

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
				# np_allDepleteSeqArr = np.array(allDepleteSeqArr)
				# print(np_allDepleteSeqArr)
		else:
			np_allEnrichSeqArr = np.array(allEnrichSeqArr)

	# split files if they get too large
	if nrow_deplete!=0:
		print("creating depleted arr")

		if (nrow_deplete > 100000000):
			num_sub_array = 20;
			sub_size = math.floor(nrow_deplete / num_sub_array);
			for index in range(0, num_sub_array):
				index_str = index;
				sub_arr_start = sub_size*index;
				sub_arr_end = sub_size*(index+1);
				if (sub_arr_end > nrow_deplete):
					sub_arr_end = nrow_enrich;
				# np_allEnrichSeqArr = np.array(allEnrichSeqArr)
				np_allDepleteSeqArr = np.array(allDepleteSeqArr[sub_arr_start:sub_arr_end])
				# print(np_allDepleteSeqArr)
		else:
			print("deplete")
			print(len(allDepleteSeqArr))
			# print(len(allDepleteSeqArr[0]))
			# print(len(allDepleteSeqArr[0][0]))
			np_allDepleteSeqArr = np.array(allDepleteSeqArr)

	if nrow_neutral!=0:
		print("creating neutral arr")
		if (nrow_neutral > 100000000):
			num_sub_array = 20;
			sub_size = math.floor(nrow_neutral / num_sub_array);
			for index in range(0, num_sub_array):
				index_str = index;
				sub_arr_start = sub_size*index;
				sub_arr_end = sub_size*(index+1);
				if (sub_arr_end > nrow_deplete):
					sub_arr_end = nrow_enrich;
				# np_allEnrichSeqArr = np.array(allEnrichSeqArr)
				np_allNeutralSeqArr = np.array(allNeutralSeqArr[sub_arr_start:sub_arr_end])
				# print(np_allDepleteSeqArr)	
		else:
			np_allNeutralSeqArr = np.array(allNeutralSeqArr)

	print(f"np_allEnrichSeqArr.shape {np_allEnrichSeqArr.shape}, np_allDepleteSeqArr.shape {np_allDepleteSeqArr.shape}, np_allNeutralArr.shape {np_allNeutralSeqArr.shape}");
	# Data = {"EnrichedData": np_allEnrichSeqArr, "DepletedData": np_allDepleteSeqArr, "NeutralData": np_allNeutralSeqArr};
	Data = {"EnrichedData": np_allEnrichSeqArr}
	# mat_filename = '/scratch2/yibeijia/data/nucleosome_occupancy_' + str(section) + '.mat';
	#mat_filename = '/Users/yibeijia/Downloads/nucleosome_occupancy/train_test_data/nucleosome_occupancy_' + str(section) + '.mat';
	# print(Data)
	# mat_filename = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/sampleSeqs.mat';
	# mat_filename = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/'+species+'Seqs_'+str(section) + '.mat';
	mat_filename = '/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/'+species+'Seqs_'+str(section) + '.mat';
	# print(Data).mat'
	scipy.io.savemat(mat_filename, Data,  do_compression=True);

def main():
	##### ########## ########## ########## ########## #####
	#conda activate mdenv

	#Get the sequencnes in depleted ir enriched regions
	#data = readInputAsArray('/Users/yibeijia/Downloads/data/GSE13622_RAW/GSM351491_InVitro_normalized.tab')
	# data = readInputAsArray('/project/rohs_108/yibeijia/nucleosome_occupancy/data/GSM351491_InVitro_normalized.tab')
	# data = readInputAsArray('/Users/yibeijia/Downloads/data/NBS/size.147.plusOrMinus1bp.worm.GSM514735.clean.seq.bed.uniq')
	# data = readInputAsArray('/Users/yibeijia/Downloads/data/NBS/size.147.plusOrMinus1bp.human.GSM907783.clean.seq.uniq')
	# getNucleosomeRegions(data,'/Users/yibeijia/Downloads/data/human_enriched_regions_out.txt', "human")

	# data = readInputAsArray('/Users/yibeijia/Downloads/data/NBS/size.147.plusOrMinus1bp.fly.Kingston.EE.clean.seq.bed.uniq')
	# getNucleosomeRegions(data,'/Users/yibeijia/Downloads/data/fly_enriched_regions_out.txt', "fly")

	# tic = time.perf_counter();
	# getNucleosomeNeutralRegions(data,'/Users/yibeijia/Downloads/data/neutral_regions_out.txt')
	# getNucleosomeNeutralRegions(data,'/project/rohs_102/share/nucleosome_occupancy_data/neutral_regions_out.txt')
	# getNucleosomeNeutralRegions(data,'/project/rohs_102/share/nucleosome_occupancy_data/neutral_regions_out.txt')
	# getNucleosomeRegions(data,'InVitro_regions_out.txt');
	# toc = time.perf_counter();
	# print(f"Getting the nucleosome regions took {toc - tic:0.4f} seconds");

	##### ########## ########## ########## ########## #####
	##### Get the sequencnes in depleted ir enriched regions
	# data = readInputAsArray('/Users/yibeijia/Downloads/data/human_enriched_regions_out.txt')
	data = readInputAsArray('/project/rohs_108/yibeijia/data/human_enriched_regions_out.txt')
	encodeNucSeq(data, sys.argv[1], sys.argv[2],'human')
	# data = readInputAsArray('/Users/yibeijia/Downloads/data/fly_enriched_regions_out.txt')
	# encodeNucSeq(data, sys.argv[1], sys.argv[2],'fly')
	# data = readInputAsArray('/project/rohs_108/yibeijia/nucleosome_occupancy/InVitro_regions_out.txt');
	# data = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/InVitro_regions_out.txt')
	# data = readInputAsArray('/Users/yibeijia/Downloads/data/worm_enriched_regions_out.txt')
	# data = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/sampleSeqs.txt')
	# encodeNucSeq(data, sys.argv[1], sys.argv[2],'worm')
	# encodeNucSeq(data, 1, 1)
main();
