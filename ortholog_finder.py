# Ortholog Finder, final project for BIOENG 140L
# written by Sara Smith

## notes section #######################################################

########################################################################

## Imports #############################################################
import sys #to read commandline arguments
from Bio.Blast import NCBIWWW # to access NCBI web server
from Bio.Blast import NCBIXML # to parse results
from Bio import SeqIO
from Bio import Entrez
########################################################################

## function define #####################################################
def print_help():
    help_message = """Ortholog Finder
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

    #update with user input
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
    if seq_file == None:
        print("missing required argument")
        print_help()
    return phelp, seq_file, orgainism, number, oligo, out
    


def blast(input_seq_file):
    sequence_data = open(input_seq_file).read()
    result_handle = NCBIWWW.qblast("tblastn", "nt", sequence_data, hitlist_size = 1, entrez_query = '(txid562[ORGN])')
    with open('results.xml', 'w') as save_file:
        blast_results = result_handle.read()
        save_file.write(blast_results)

def generate_oligos(output_seq_file):
    pass

def get_sequence(accession_number, start, end):
    Entrez.email = "sara.smith@berkeley.edu" #change this to your email
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    ortholog = record.seq[(start-1):end]
    pre = record.seq[(start-201):start]
    post = record.seq[(end-1):(end+200)]
    return ortholog, pre, post
    #//TODO:more or less working, could have off by one errors 
    #takes about 30 seconds to run



#####################################################################################

## workaround ######### https://stackoverflow.com/questions/27835619/urllib-and-ssl-certificate-verify-failed-error#comment72358900_42334357, https://github.com/biopython/biopython/issues/3671
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
########

# main part of the program ##########################################################
help_message, input_seq_file, input_organism, num, generate_oligos, out_directory = parse_arguments()
if help_message:
    print_help()
#blast(input_seq_file)

E_VALUE_THRESH = 1e-20 
for record in NCBIXML.parse(open("results.xml")): 
     if record.alignments: 
        print("\n") 
        print("query: %s" % record.query[:100]) 
        for align in record.alignments: 
           for hsp in align.hsps: 
              if hsp.expect < E_VALUE_THRESH: 
                 print("match: %s " % align.title[:100])

print("here")
print(get_sequence("LR133996.1", 2987872, 2989731))










