'''
This is a solved problem from Rosalind website on bioinformatics problems

Title: Identifying Maximal Repeats
ID: MREP

Given: A DNA string s of length at most 1 kbp.
Return: A list containing all maximal repeats of s having length at least 20.
'''


def extend_right(seq, rep_seq, rep_num):

    """    Input: DNA sequence (seq) as a string, a repeating substring of the DNA sequence (rep_seq) and the number of
    times the substring repeats itself on the DNA sequence (rep_num)

    Output: the repeating substring (rep_seq) but one base extended to the right if possible to do so without eliminating
    a repeat (otherwise will return None) and a dictionary where the keys are the possible one base extensions of the
    repeating substring and their respective values are the number of times they repeat in the DNA sequence"""

    right_ext_dict = {rep_seq + 'A': 0,
                      rep_seq + 'T': 0,
                      rep_seq + 'G': 0,
                      rep_seq + 'C': 0}  # extended repeating sequences to the right of the given repeating sequence

    for ext_rep_seq in right_ext_dict:
        right_ext_dict[ext_rep_seq] = seq.count(ext_rep_seq)
        # if there were an agreement between all after the extension, there should be one ext_rep_seq with the same
        # repetition as the rep_seq.
        if right_ext_dict[ext_rep_seq] == rep_num:
            return ext_rep_seq, right_ext_dict
    return None, right_ext_dict


def extend_left(seq, rep_seq, rep_num):

    """    Input: DNA sequence (seq) as a string, a repeating substring of the DNA sequence (rep_seq) and the number of
    times the substring repeats itself on the DNA sequence (rep_num)

    Output: the repeating substring (rep_seq) but one base extended to the left if possible to do so without eliminating
    a repeat (otherwise will return None) and a dictionary where the keys are the possible one base extensions of the
    repeating substring and their respective values are the number of times they repeat in the DNA sequence"""

    left_ext_dict = {'A' + rep_seq: 0,
                     'T' + rep_seq: 0,
                     'G' + rep_seq: 0,
                     'C' + rep_seq: 0}  # extended repeating sequences to the left of the given repeating sequence
    
    for ext_rep_seq in left_ext_dict:
        left_ext_dict[ext_rep_seq] = seq.count(ext_rep_seq)
        # if there were an agreement between all after the extension, there should be one ext_rep_seq with the same
        # repetition as the rep_seq.
        if left_ext_dict[ext_rep_seq] == rep_num:
            ext_rep_seq += ext_rep_seq[0]
            return ext_rep_seq, left_ext_dict
    return None, left_ext_dict


def maximal_repeats(seq, rep_seq='', rep_num=0):
    """    Input: a DNA sequence as a string (seq)
    Output: all of the repeating sequences in the given DNA that cannot be extended to the right with agreement between
    all of them (the repeating sequences)"""

    right_ext_seq, right_ext_dict = extend_right(seq, rep_seq, rep_num)
    max_rep_list = []  # all of the repeating sequences found
    extended_rep_seq = rep_seq

    # if cannot be extended to both left and right with agreement, than it is a maximal repeat
    if right_ext_seq is None and rep_seq != '':
        left_ext_seq = extend_left(seq, rep_seq, rep_num)[0]
        if left_ext_seq is None:
            max_rep_list.append(rep_seq)

    for new_rep_seq in right_ext_dict:
        # if the new sequence after the extension repeats itself, call the function to extend it further to the right
        if right_ext_dict[new_rep_seq] > 1:
            max_rep_list += maximal_repeats(seq, new_rep_seq, right_ext_dict[new_rep_seq])

    return max_rep_list


# the following will print all of the maximal repeats in the defined 'seq' that are bigger or equal to length of 20.
seq = 'TAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTATTATATAGAGATAGAATGGGTCCAGAGTTTTGTAATTTCCATGGGTCCAGAGTTTTGTAATTTAT'
max_reps = maximal_repeats(seq)

print([max_rep for max_rep in max_reps if len(max_rep) >= 20])
