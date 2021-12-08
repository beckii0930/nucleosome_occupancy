import numpy as np;
import time
import scipy.io
import sys
import math

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
	allEnrichSeqArr='';
	allDepleteSeqArr='';
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
					# encodedDNAArr = oneHotEncode(curr_seq);
					allEnrichSeqArr+= (curr_seq + '\n');
				else:
					## for seqeunces that are long enough	
					for start in range(curr_seq_length-146):
						curr_start = start;
						curr_end = 146 + start;
						curr_seq = line[3][curr_start: curr_end+1];
						# encodedDNAArr = oneHotEncode(curr_seq);
						allEnrichSeqArr+=(curr_seq + '\n');

			else:
				## for seqeunces that are shorter
				if(len(curr_seq) < 147):	
					while len(curr_seq) < 147:
						curr_seq = curr_seq + 'N';
					# encodedDNAArr = oneHotEncode(curr_seq);
					allDepleteSeqArr+=(curr_seq + '\n');
				else:
					## for seqeunces that are long enough	
					for start in range(curr_seq_length-146):
						curr_start = start;
						curr_end = 146 + start;
						curr_seq = line[3][curr_start: curr_end+1];
						# encodedDNAArr = oneHotEncode(curr_seq);
						allDepleteSeqArr += (curr_seq + '\n');
		else:
			print("header")
	# print(allEnrichSeqArr);
	toc = time.perf_counter();
	print(f"Getting the DNA encoded took {toc - tic:0.4f} seconds");
	nrow_enrich = len(allEnrichSeqArr);
	nrow_deplete = len(allDepleteSeqArr);

	num_sub_array = 10;
	sub_size = math.floor(nrow_enrich / num_sub_array);
	for index in range(0, num_sub_array):
		index_str = index;
		sub_arr_start = sub_size*index;
		sub_arr_end = sub_size*(index+1);
		if (sub_arr_end > nrow_enrich):
			sub_arr_end = nrow_enrich;
		# np_allEnrichSeqArr = np.array(allEnrichSeqArr[sub_arr_start:sub_arr_end])
		# np_allDepleteSeqArr = np.array(allDepleteSeqArr)
		# print(np_allDepleteSeqArr)
		print(f"1 np_allEnrichSeqArr.shape {np_allEnrichSeqArr.shape}, np_allDepleteSeqArr.shape {np_allDepleteSeqArr.shape}");
		# Data = {"EnrichedData": np_allEnrichSeqArr, "DepletedData": np_allDepleteSeqArr};
		out_filename_d = '/project/rohs_102/share/nucleosome_occupancy_depleted_'+ str(index_str) + '.txt';
		out_filename_e = '/project/rohs_102/share/nucleosome_occupancy_enriched_'+ str(index_str) + '.txt';
		#out_filename = '/Users/yibeijia/Downloads/nucleosome_occupancy/train_test_data/nucleosome_occupancy_' + str(section) + '_' + str(index_str) + '.mat';
		file1 = open(out_filename_d, "w+")
		file1.write(allDepleteSeqArr)
		file1.close()
		file2 = open(out_filename_e, "w+")
		file2.write(allEnrichSeqArr)
		file2.close()




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
	data = readInputAsArray('/project/rohs_108/yibeijia/nucleosome_occupancy/InVitro_regions_out.txt');
	#data = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/InVitro_regions_out.txt')
	encodeNucSeq(data, sys.argv[1], sys.argv[2])
	# dna="ACGTAC";
	# encodedDna=oneHotEncode(dna)
	# print("dna\n",list(dna))
	# print('encoded dna\n')
main();
