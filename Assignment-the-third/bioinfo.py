#!/usr/bin/env python

# Author: <Daraudom Nhem <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.1"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = 'ACTG'
RNA_bases = 'ACUG'

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33
    

def qual_score(phred_score: str) -> float:
    '''Given a phred string, convert each character into a quality score (Phred 33)'''
    n = len(phred_score)
    total_score = 0
    for x in phred_score:
        total_score += convert_phred(x)

    return total_score/n


def validate_base_seq(seq, RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    if len(seq) == 0: return False
    valid_bases = RNA_bases if RNAflag else DNA_bases
    return all([base in valid_bases for base in seq.upper()])

def gc_content(seq):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    gc_count = 0
    if set(seq).issubset(set(DNA_bases)) or set(seq).issubset(set(RNA_bases)):
        for base in seq.upper():
            if base == 'G' or base == 'C': gc_count += 1
    else:
        raise AssertionError("Not a DNA String")
    return gc_count/len(seq)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    n = len(lst)
    if n%2 == 0:
        return (lst[n//2] + lst[(n//2) - 1])/2
    else:
        return lst[n//2]

def oneline_fasta(input_file, output_file):
    '''Take in an input fasta file and convert it to a one line fasta.
    The format is >header then an one line sequence.
    '''
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        header_record = True
        for line in fin:
            if header_record:
                fout.write(line)
                header_record = False
            elif line.startswith('>'):
                fout.write(f'\n{line}')
            else:
                fout.write(line.strip())
    

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    # test_functions for qual_Score
    assert qual_score("I") == 40.0, "Expected score of 40.0 for 'I'"
    assert qual_score("IIII") == 40.0, "Expected score of 40.0 for 'IIII'"
    assert qual_score("ABCDE") == 34.0, "Expected score of 34.0 for 'ABCDE'"
    assert qual_score("!!!!") == 0.0, "Expected score of 0.0 for '!!!!'"
    print("Successfully passed Quality Score test!")

    # test_functions for validate_base_seq
    assert validate_base_seq("") == False, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

    # test_functions for gc_content
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    print("correctly calculated GC content")

    # test_functions for calc_median
    assert calc_median([9,15]) == 12
    assert calc_median([25,36,47,58]) == 41.5
    assert calc_median([2,21,101,13,41,1,31,1,100]) == 41
    assert calc_median([72]) == 72
    assert calc_median([50,100]) == 75
    print("correctly calculated medians")


