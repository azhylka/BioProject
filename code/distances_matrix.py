import csv
import numpy
import math

def parse_line( line ):
	strain_end = line.find('\t')
	strain = line[:strain_end]
	parsed_line = [parse_int(x) for x in line[strain_end + 1:].split("\t")]
	return {strain : parsed_line}

def parse_int( entry ): 	# in case of missing entry ('-')
	if '-' in entry:
		return 0
	else:
		return int(entry)

def read_genome(source):
	genome_map = {}
	genome_snps = []
	source_file = open(source)
	genome_snps = source_file.readlines()[2:] #skip header and annotations
	source_file.close()
	for line in genome_snps:		
		genome_map.update(parse_line(line))	     
	return genome_map

def hamming_distance(first_code, second_code):	    
	if len(first_code) != len(second_code):
		raise ValueError("Undefined for sequences of unequal length")
	return sum(bit1 != bit2 for bit1, bit2 in zip(first_code, second_code))

def normalize_distance(distance):
	if distance == 0:
		return distance
	else:
		return math.log(distance)

def build_distance_map(genome_map):
	distances_map = {}
	for strain, snps in genome_map.iteritems():
		for strain2, snps2 in genome_map.iteritems():				
#			distances_map[strain, strain2] = normalize_distance(hamming_distance(snps, snps2))
			distances_map[strain, strain2] = hamming_distance(snps, snps2)
	return distances_map

def get_sub_map(original_map, subset):
	return {(key1, key2) : original_map[key1, key2] for key1 in subset for key2 in subset}


def export_to_csv(doublekey_map, header, filename):
	csv_file = open(filename, 'wb')
	row = ['\t']
	row[:0] = header
	for name in row:
		csv_file.write(','+name)	
	for row_name in header:
		csv_file.write('\n'+row_name)
		for column_name in header:
			csv_file.write(','+str(doublekey_map[row_name, column_name]))			
	csv_file.close()

def main():
	genome_map = read_genome("apr15.snps.matrix")
	strains = genome_map.keys()
	distance_map = build_distance_map(genome_map)
	control_snps_subset = ['G55417', 'G55419', 'G55424', 'G55426', 'G55432',
			'G55438', 'G55448', 'G55455', 'G55461', 'G55464',
			'G55469', 'G55474', 'G55478', 'G55481', 'G55500',
			'G55505', 'G55507', 'G55510', 'G55512', 'G55514',
			'G55522', 'G55524', 'G55526', 'G55543', 'G55546',
			'G55553', 'G55560']
	control_snps_distances = get_sub_map(distance_map, control_snps_subset)
	export_to_csv(distance_map, sorted(strains), 'resources/distance_table.csv')
	export_to_csv(control_snps_distances, sorted(control_snps_subset), 'resources/control_distance_table.csv')


if __name__ == "__main__":
	main()
