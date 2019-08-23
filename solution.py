
amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
amino_acids = list(amino_acids) + ['-']

def input_sequences(file):
	sequences = []
	with open(file) as f:
		for line in f:
			line = line.strip('\n')
			line = line[30: ]
			#print(line + '\n', len(line))
			sequences.append(line)
	return sequences

def consensus_1a(sequences, amino_acids, pseudocount = 0):
	sequence_length = len(sequences[0])
	profile_matrix = {}
	for acid in amino_acids:
		profile_matrix[acid] = [pseudocount for i in range(sequence_length)]
	
	for i in range(len(sequences)):
		seq = sequences[i]
		for j in range(len(seq)):
			#print(i, j)
			profile_matrix[seq[j]][j] += 1

	for aa in profile_matrix:
		l = profile_matrix[aa]
		for i in range(len(l)):
			if pseudocount > 0:
				l[i] /= (len(sequences)*2)
			else:
				l[i] /= len(sequences)

	consensus_string = ''
	for i in range(sequence_length):
		l = []
		for aa in profile_matrix:
			l.append(profile_matrix[aa][i])
		index = l.index(max(l))
		consensus_string += amino_acids[index]

	return consensus_string, len(consensus_string)


file = 'PF00167_full.txt'
sequences = input_sequences(file)
#print(sequences[87][258: 261])
print(consensus_1a(sequences, amino_acids))
