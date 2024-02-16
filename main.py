
def read_input_file(filename):
    with open(filename, 'r') as file:
        sequence = file.read()
        sequence_list = [char for char in sequence if char != ' ' and char != '\n']
    return sequence_list




def get_basecounts(list):
    '''


    :param list:
    :return: the count of each nucleotid
    '''
    count = {}
    for base in list:
        if base in count:
            count[base] += 1
        else:
            count[base] = 1
    return count



def get_complement(seq, reverse=False):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_sequence = ''.join(complement_dict[base] for base in seq)
    if reverse:
        complement_sequence = complement_sequence[::-1]
    return complement_sequence

def test_rev_comp():
    assert get_complement('ATCG')=='TAGC', "complement test"
    assert get_complement('ATCG', True)=='CGAT', "reverse complement test"
    print("Tests passed")

def calculate_gc_content(list):
    '''
    calculate the gc content
    :param list:
    :return:
    '''
    gc = [char for char in list if char == "G" or char =="C"]
    content = (len(gc) / len(list))*100
    return content



def test_gc_content():
    '''
    test the gc content
    :return:
    '''
    assert calculate_gc_content('ggggaaaaaaaatttatatatcgcc')==0.32, "gc_content test"
    print("Tests passed")

def detect_gc_islands(sequence, window_size, gc_threshold):
    '''

    :param sequence:
    :param window_size:
    :param gc_threshold:
    :return: the island of gc
    '''
    gc_islands = []
    for i in range(0, len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        gc_content = calculate_gc_content(window)
        if gc_content >= gc_threshold:
            start = i
            end = i + window_size - 1
            gc_islands.append((start, end, gc_content))
    return gc_islands


def main(input_files, base_count=False, reverse_complement=False, GC_content=False, number_of_islands=False, output_file=None):
    if len(input_files) < 3:
        print("please input 3 or more files")
        return

    sequences = []
    for file in input_files:
        sequence = read_input_file(file)
        sequences.append(sequence)

    if base_count:
        for i, sequence in enumerate(sequences, start=1):
            base_counts = get_basecounts(sequence)
            print(f"Base counts for Sequence {i}: {base_counts}")

    if reverse_complement:
        for i, sequence in enumerate(sequences, start=1):
            rc_sequence = get_complement(sequence, reverse=True)
            print(f"Reverse complement of Sequence {i}: {rc_sequence[:50]}")

    if GC_content:
        for i, sequence in enumerate(sequences, start=1):
            gc_content = calculate_gc_content(sequence)
            print(f"GC content of Sequence {i}: {gc_content:.2f}")

    if number_of_islands:
        window_size = 200
        gc_threshold = 70
        for i, sequence in enumerate(sequences, start=1):
            gc_islands = detect_gc_islands(sequence, window_size, gc_threshold)
            print(f"Number of GC islands in Sequence {i}: {len(gc_islands)}")
            for start, end, gc_content in gc_islands:
                print(f"GC island position: Start: {start}, End: {end}, GC Content: {gc_content:.2f}")

    if output_file:
        with open(output_file, 'w') as f:
            for operation in [base_count, reverse_complement, GC_content, number_of_islands]:
                if operation:
                    f.write(f"Operation: {operation}\n")
                    for i, sequence in enumerate(sequences, start=1):
                        if base_count:
                            base_counts = get_basecounts(sequence)
                            f.write(f"Base counts for Sequence {i}: {base_counts}\n")
                        if reverse_complement:
                            rc_sequence = get_complement(sequence, reverse=True)
                            f.write(f"Reverse complement of Sequence {i}: {rc_sequence[:50]}\n")
                        if GC_content:
                            gc_content = calculate_gc_content(sequence)
                            f.write(f"GC content of Sequence {i}: {gc_content:.2f}\n")
                        if number_of_islands:
                            window_size = 200
                            gc_threshold = 70
                            gc_islands = detect_gc_islands(sequence, window_size, gc_threshold)
                            f.write(f"Number of GC islands in Sequence {i}: {len(gc_islands)}\n")
                            for start, end, gc_content in gc_islands:
                                f.write(f"GC island position: Start: {start}, End: {end}, GC Content: {gc_content:.2f}\n")


if __name__ == "__main__":
    main(["DNA_lambdavirus.fasta", "redmt_DNA_algae.fasta", "redmt_DNA_greentea.fasta"],True,True,True,True,"output_file")
    #print(detect_gc_islands(read_input_file("DNA_lambdavirus.fasta"), 50, 0.4))





