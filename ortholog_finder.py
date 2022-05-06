## notes section #######################################################

# Ortholog Finder, final project for BIOENG 140L
# written by Sara Smith

########################################################################

## Imports #############################################################
import sys #to read commandline arguments
import os # helps with creating output directories
from Bio.Blast import NCBIWWW # to access NCBI web server
from Bio.Blast import NCBIXML # to parse results
from Bio import SeqIO
from Bio import Entrez # to access sequences
########################################################################

## function define #####################################################
def print_help():
    help_message = """Ortholog Finder
    version 1
    Arguments:
    -h :print help menu
    -i :input sequence in fasta format (no default)
    -orgs :file of taxonomy number of organisms to search (default e.coli, 562)
    -oligo :if included, generate suggested oligos
    -out :output directory name (default 'Ortholog_Finder_Output' in current directory)
    """
    print (help_message)
    sys.exit("exiting the program")

def parse_arguments() :
    #set defaults
    phelp = False
    seq_file = None
    organisms = "default_organisms.txt"
    oligo = False
    out = "./Ortholog_Finder_Output/"

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
        elif token == "-orgs":
            organisms = sys.argv[i+1]
            i += 2
        elif token == "-oligo":
            oligo = True
            i += 1
        elif token == "-out":
            out = sys.argv[i+1]
            i += 2
        else:
            print(sys.argv[i], "is an unrecognized token. use -h to see the list of accepted arguments")
            i += 1
            continue
    if seq_file == None:
        print("missing required argument")
        print_help()
    return phelp, seq_file, organisms, oligo, out
    


def blast(input_seq_file, organism, out):
    sequence_data = open(input_seq_file).read()
    result_handle = NCBIWWW.qblast("tblastn", "nt", sequence_data, alignments = 1, hitlist_size = 1, entrez_query = '(txid{}[ORGN])'.format(organism))
    filename = "{}{}/results.xml".format(out, organism)
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w') as save_file:
        blast_results = result_handle.read()
        save_file.write(blast_results)
    #this can be quite slow. it starts by polling NCBI after 20 seconds, then polls once per minute
    
    #gets just the highest scoring alignment. Future improvements to this program could consider more than one
    for record in  NCBIXML.parse(open("{}{}/results.xml".format(out, organism))):
        acc = None
        start = None
        end = None
        if record.alignments: 
            align = record.alignments[0]
            acc = align.accession
            hsp = align.hsps[0]
            start = hsp.sbjct_start
            end = hsp.sbjct_end
    return acc, start, end


def get_sequence(accession_number, start, end):
    Entrez.email = "sara.smith@berkeley.edu" #change this to your email
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    segment = record.seq[(start-201):end+200]
    start_codon = segment[200:203]
    if not (start_codon == "ATG" or start_codon == "GTG" or start_codon == "TTG"):
        segment = segment.reverse_complement()
        if not (start_codon == "ATG" or start_codon == "GTG" or start_codon == "TTG"):
            segment = segment.reverse_complement()
    pre = segment[0:200]
    ortholog = segment[200:len(segment)-200]
    post = segment[(len(segment)-201):]

    return ortholog, pre, post
    #//TODO:more or less working, could have off by one errors 
    #takes about 30 seconds to run



#####################################################################################

## workaround ######### https://stackoverflow.com/questions/27835619/urllib-and-ssl-certificate-verify-failed-error#comment72358900_42334357, https://github.com/biopython/biopython/issues/3671
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
#####################################################################################

# main part of the program ##########################################################
print("Starting Program")
print("Parsing arguments.....")
help_message, input_seq_file, input_organism, generate_oligos, out_directory = parse_arguments()
if help_message:
    print_help()

#loop through organisms
organisms = open(input_organism, 'r')
for organism in organisms.readlines():
    print("Starting blast on {}.....".format(organism.strip()))
    acc, start, end = blast(input_seq_file, organism.strip(), out_directory)  
    print("Getting sequences for {}.....".format(organism.strip()))
    orth, pre, post = get_sequence(acc, start, end)
    print(orth, pre, post)


    #generate oligos
    if generate_oligos:
        print("Generating oligos for {}.....".format(organism.strip()))
        fwd = "ccaaaGGTCTCaCCCT" + pre[:20]
        rev =  "gctagGGTCTCgGAAC" + post[:20].reverse_complement()

    #write sequences to a fasta file 
    print("Writing file for {}.....".format(organism.strip()))       
    filename = "{}{}/result_sequences.fasta".format(out_directory, organism.strip())
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w') as save_file:
        save_file.write('>ortholog\n')
        save_file.write('{}\n'.format(orth))
        save_file.write('>pre_200\n')
        save_file.write('>{}\n'.format(pre))
        save_file.write('>post_200\n')
        save_file.write('>{}\n'.format(post))
        if generate_oligos:
            save_file.write('>forward_oligo\n')
            save_file.write('>{}\n'.format(fwd))
            save_file.write('>reverse_oligo\n')
            save_file.write('>{}\n'.format(rev))

    #check for BSA1 sites and start codons //TODO: make the enzyme a parameter
    print("Checking {} for restriction sites.....".format(organism.strip()))
    filename = "{}{}/warnings.txt".format(out_directory, organism.strip())
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w') as warn_file:
        if "GGTCTC" in orth:
            warn_file.write("WARNING: Ortholog contains a BSA1 site\n")
        if "CCAGAG" in orth:
            warn_file.write("WARNING: Ortholog contains a BSA1 site\n")
        if "GGTCTC" in pre:
            warn_file.write("WARNING: Pre-200 sequence contains a BSA1 site\n")
        if "CCAGAG" in pre:
            warn_file.write("WARNING: Pre-200 sequence contains a BSA1 site\n")
        if "GGTCTC" in post:
            warn_file.write("WARNING: Post-200 sequence contains a BSA1 site\n")
        if "CCAGAG" in post:
            warn_file.write("WARNING: Post-200 sequence contains a BSA1 site\n")
        start_codon = orth[:3]
        if not (start_codon == "ATG" or start_codon == "GTG" or start_codon == "TTG"):
            warn_file.write("WARNING: Alignment does not start with a start codon\n")
    


sys.exit("complete")








