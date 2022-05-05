# Ortholog Finder, final project for BIOENG 140L
# written by Sara Smith

# notes section ########################################################
#sys.argv[i]
#print("\nArguments passed:", end = " ")
#for i in range(1, n):
#    print(sys.argv[i], end = " ")
#https://www.geeksforgeeks.org/how-to-use-sys-argv-in-python/
########################################################################

## Imports
import sys #to read commandline arguments
#from Bio.Blast import NCBIWWW
## read commandline arguments

def print_help():
    help_message = 
    """Ortholog Finder
    version 1
    Arguments:
    -h :print help menu
    -i :input sequence in fasta format (no default)
    -org :input organism as string (no default)
    -n :number of orthologs to return (default 1)
    -oligo :if included, generate suggested oligos
    -out :output file name (default 'Ortholog_Finder_Output' in current directory)
    """
    print (help_message)
    sys.exit("exiting the program")

def parse_arguments() :
    #set defaults
    phelp = False
    seq_file = None
    orgainism = "testing parse... organism"
    number = 1
    oligo = False
    out = "./Ortholog_Finder_Output"

    i = 1
    while i < len(sys.argv):
        token = sys.argv[i]
        if token == "-h":
            phelp = True
            i += 1
        elif token == "-i":
            seq_file = sys.argv[i+1]
            print("printing seq_file", seq_file)
            i += 2
        elif token == "-org":
            organism = sys.argv[i+1]
            i += 2
        elif token == "-n":
            number = sys.argv[i+1]
            i += 2
        elif token == "-oligo":
            oligo = True
            i += 1
        elif token == "-out":
            out = sys.argv[i+1]
            i += 2
        else:
            print(sys.argv[i], "is an unrecognized token. use -h to see the list of accepted arguments")
            break
    return phelp, seq_file, orgainism, number, oligo, out
    


def blast(input_seq_file):
    sequence_data = open(input_seq_file).read() 
    result_handle = NCBIWWW.qblast("tblastn", "nt", sequence_data)
    with open('results.xml', 'w') as save_file:
        blast_results = result_handle.read()
        save_file.write(blast_results)
    pass

def generate_oligos():
    pass


# main part of the program ##########################################################
help_message, input_seq_file, input_organism, num, generate_oligos, out_directory = parse_arguments()
#testing parse
print(help_message, input_seq_file, input_organism, num, generate_oligos, out_directory)
if help_message:
    print_help()









