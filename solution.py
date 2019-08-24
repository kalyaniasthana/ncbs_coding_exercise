
amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
amino_acids = list(amino_acids) + ['-']

def input_sequences(file):
	sequences = []
	sequences_with_name = []
	with open(file) as f:
		for line in f:
			line = line.strip('\n')
			sequences_with_name.append(line)
			line = line[30: ]
			#print(line + '\n', len(line))
			sequences.append(line)
	return sequences, sequences_with_name

def second_largest(numbers):
    count = 0
    m1 = m2 = float('-inf')
    for x in numbers:
        count += 1
        if x > m2:
            if x >= m1:
                m1, m2 = x, m1            
            else:
                m2 = x
    return m2 if count >= 2 else None

def consensus(x, sequences, sequences_with_name, amino_acids, pseudocount = 0):
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


	if x == 7:
		#identifying the position in the alignment which has the most dashes
		max_value = max(profile_matrix['-'])
		if max_value == 1:
			max_value = second_largest(profile_matrix['-'])

		positions = []
		for i in range(len(profile_matrix['-'])):
			if profile_matrix['-'][i] == max_value:
				positions.append(i)
		#return positions, max_value

		bad_sequence_numbers = []
		for i in range(len(sequences)):
			for position in positions:
				if sequences[i][position] != '-':
					if i not in bad_sequence_numbers:
						bad_sequence_numbers.append(i)

		bad_sequences = []
		for number in bad_sequence_numbers:
			bad_sequences.append(sequences_with_name[number][0: 30])

		return bad_sequences

	elif x == 1 or x == 2:
		for i in range(sequence_length):
			l = []
			for aa in profile_matrix:
				l.append(profile_matrix[aa][i])
			index = l.index(max(l))
			consensus_string += amino_acids[index]

		if x == 1:
			return consensus_string
		else:
			return len(consensus_string)

	elif x == 3 or x == 4:
		for i in range(sequence_length):
			l = []
			for aa in profile_matrix:
				l.append(profile_matrix[aa][i])
			index = l.index(max(l))
			if amino_acids[index] == '-':
				continue
			consensus_string += amino_acids[index]

		if x == 3:
			return consensus_string
		else:
			return len(consensus_string)

	elif x == 5 or 6:
		for i in range(sequence_length):
			l = []
			for aa in profile_matrix:
				l.append(profile_matrix[aa][i])
			index = l.index(max(l))
			if amino_acids[index] == '-':
				if l[index] < 0.5:
					index = l.index(second_largest(l))
				elif l[index] >= 0.5:
					continue
			consensus_string += amino_acids[index]

		if x == 5:
			return consensus_string
		else:
			return len(consensus_string)


file = 'PF00167_full.txt'
sequences, sequences_with_name = input_sequences(file)
#print(sequences[87][258: 261])
print('1a. - Enter 1, 1b. - Enter 2, 2a. - Enter 3, 2b. - Enter 4, 3a. - Enter 5, 3b. - Enter 6, 4. - Enter 7')
x = int(input())
#print(consensus_1a(x, sequences, amino_acids))
#print(consensus_2a(sequences, amino_acids))
#print(consensus_3a(sequences, amino_acids))
#print(bad_sequences(sequences, sequences_with_name, amino_acids))
print(consensus(x, sequences, sequences_with_name, amino_acids))
