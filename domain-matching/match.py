#Identifies proteins in a given genome containing the protein domain containing a certain profile
import sys
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

def rpstblastn(input_file_name, output_file_name, database_path):
    subprocess.run(["rpstblastn", "-query", input_file_name, "-out", output_file_name, "-db", database_path])

def read_output(filename):
    with open(filename) as f:
        for record in NCBIXML.parse(f):
            if record.alignments:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        print(hsp)
                    
def convert_to_protein_fasta(filename, profile_name):
    record = SeqIO.read(filename + ".gbk", "genbank")
    proteins = []
    index = 0
    for feature in filter(lambda f: f.type == "gene", record.features):
        gene = feature.extract(record).seq
        protein = SeqRecord(gene.transcribe().translate())
        protein.id = str(index)
        index += 1
        proteins.append(protein)
        print(protein)

    SeqIO.write(proteins, filename + ".fsa", "fasta")

def make_protein_blast(filename):
    record = SeqIO.read(filename + ".gbk", "genbank")
    proteins = []
    index = 0
    for feature in filter(lambda f: f.type == "gene", record.features):
        gene = feature.extract(record).seq
        protein = SeqRecord(gene.transcribe().translate())
        protein.id = str(index)
        index += 1
        proteins.append(protein)
        print(protein)

    SeqIO.write(proteins, filename + ".fsa", "fasta")

def extract_upstream(filenames, number):
    record = SeqIO.read(filename + ".gbk", "genbank")
    records = []
    index = 0
    for filename in filenames:
        for feature in filter(lambda f: f.type == "gene", record.features):
            end = feature.location.start
            if start > 0:
                start = max(end - number, 0)
                upstream = record.seq[start..end]
                upstream.id = str(index)
                index += 1
                records.append(Record)
                print(upstream)

    SeqIO.write(records, "output.fsa", "fasta")

read_output("hits")
