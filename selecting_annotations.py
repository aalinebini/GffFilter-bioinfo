import sys
import getopt

from Bio import SeqIO
import pandas as pd
import re
import csv

class Selecting():
    """A class to select the gene annotations in the .gff file by the fasta file

    Returns:
        [gff] -- A GFF with gene annotations selected
    """

    def __init__(self):
        self.list_fasta = None
        self.result = None

    def reading_fasta(self, fasta_path):
            """this method loads a fasta sequence and create a list with
            the name and proteinId

            Args:
                fasta_path ([str]): the fasta file path
            """
            
            arq_fasta = SeqIO.parse(fasta_path, "fasta")
            
            self.list_fasta = list()
            
            regex_proteinID = re.compile(r"\|[0-9]+\|")
            regex_name = re.compile(r"\|[A-Za-z0-9_.#]+$")

            for seq_record in arq_fasta:
                
                match_protein = regex_proteinID.search(seq_record.id)
                proteinID = seq_record.id[(match_protein.start()+1):(match_protein.end()-1)]
                
                match_name = regex_name.search(seq_record.id)
                name = seq_record.id[(match_name.start()+1):(match_name.end())]

                self.list_fasta.append(str(name + ',' + proteinID))

    def comparisons(self, gff_path):
        """this method select the gff sequences when there's a match
            between name and proteinID of fasta list (reading_fasta) and gff file

            Args:
                gff_path ([str]): gff file path
        """

        self.result = pd.read_csv(gff_path, sep='\t', header=None)

        self.result = self.result.rename(columns={0:'SeqId', 1:'Source', 2:'type', 3:'start', 4:'end', 5:'score', 6:'strand', 7:'phase', 8:'attributes'})

        self.result = self.result.join(self.result['attributes'].str.split(';', expand=True))
        self.result.rename(columns={0:'name', 1:'proteinId', 2:'exonNumber', 3:'product_name'}, inplace=True)
        self.result.drop('exonNumber', inplace=True, axis=1)
        self.result.drop('product_name', inplace=True, axis=1)

        self.result['name'] = self.result['name'].str.extract('"(.+)"')
        self.result['proteinId'] = self.result['proteinId'].str.extract(r'(proteinId [0-9]+)')
        self.result['proteinId'] = self.result['proteinId'].str.extract(r'([0-9]+)')

        self.result.dropna(subset=['proteinId'], inplace=True)

        self.result['nameProt'] = self.result['name'] + ',' + self.result['proteinId']

        array_bool = self.result.nameProt.isin(self.list_fasta)

        self.result = self.result[array_bool]

        self.result.drop('name', inplace=True, axis=1)
        self.result.drop('proteinId', inplace=True, axis=1)
        self.result.drop('nameProt', inplace=True, axis=1)

                
    def save_to_GFF(self, output):
        """A method to save the ordenated sequence into a .gff file

        Args:
            output ([gff]): name to save the gff file
        """

        self.result.to_csv(output+'.gff', sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)


if __name__ == "__main__":

    try:
        OPTS, ARGS = getopt.getopt(sys.argv[1:], 'f:g:o:h', ['fasta_path', 'gff_path', 'output', 'help'])

    except getopt.GetoptError as err:
        print(err)
        sys.exit(1)

    SELECTING = Selecting()
    FASTA_PATH = None
    GFF_PATH = None
    OUTPUT = None

    for opt, arg in OPTS:

        if opt in ('-h', '--help'):
            print('''
            This program select the gene annotations in the .gff file by the fasta file.
            selecting_annotations.py -f --fasta_path -g --gff_path -o --output -h --help
            ''')
            sys.exit(2)

        elif opt in ('-f', '--fasta_path'):
            if re.match('.+fasta$|.+txt$', arg):
                FASTA_PATH = arg
            else:
                print("-f isn't a fasta or txt file")
                sys.exit(3)

        elif opt in ('-g', '--gff_path'):
            if re.match('.+gff$|.+txt$', arg):
                GFF_PATH = arg
            else:
                print("-g isn't a gff or txt file")
                sys.exit(4)

        elif opt in ('-o', '--output'):
            if re.match('.+', arg):
                OUTPUT = arg
            else:
                print("-o missing the output name")
                sys.exit(5)

    if FASTA_PATH and GFF_PATH and OUTPUT:
        SELECTING.reading_fasta(FASTA_PATH)
        SELECTING.comparisons(GFF_PATH)
        SELECTING.save_to_GFF(OUTPUT)

    else:
        print('missing arguments -f, -g or -o')
        sys.exit(6)
