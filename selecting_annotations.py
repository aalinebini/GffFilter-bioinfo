import sys
import getopt
import re
from Bio import SeqIO

class Selecting():
    """A class to select the gene annotations in the .gff file by the fasta file

    Returns:
        [gff] -- A GFF with gene annotations selected
    """

    def __init__(self):
        self.dict_fasta = None

    def reading_fasta(self, fasta_path):
        """this method loads a fasta sequence and create a dict with
        the name and proteinId ('name':proteinId)

        Args:
            fasta_path ([str]): the fasta file path
        """
        
        arq_fasta = SeqIO.parse(fasta_path, "fasta")
        
        self.dict_fasta = dict()
        
        regex_proteinID = re.compile(r"\|[0-9]+\|")
        regex_name = re.compile(r"\|[A-Za-z0-9_.#]+$")

        for seq_record in arq_fasta:
            
            match_protein = regex_proteinID.search(seq_record.id)
            proteinID = seq_record.id[(match_protein.start()+1):(match_protein.end()-1)]
            
            match_name = regex_name.search(seq_record.id)
            name = seq_record.id[(match_name.start()+1):(match_name.end())]

            self.dict_fasta[name] = proteinID
    
    def comparisons(self, gff_path, output):
        """this method select the gff sequences when there's a match
        between name and proteinID of dict (fasta file) and gff file

        Args:
            gff_path ([str]): gff file path
            output ([str]): output name
        """
    
        arq_gff = open(gff_path,'r')
        new_arq_gff = open(output, 'w')
        
        regex_name = re.compile('"(.+)"')
        regex_proteinID = re.compile('[0-9]+')
        
        for line in arq_gff:
                    
            line_analysis = line
            
            if 'proteinId' in line_analysis:

                line_analysis = line_analysis.split('\t')[8].split(';')

                try:
                    name = line_analysis[0]
                    name = str(regex_name.findall(name)).replace("'", "").strip('[]')
                    
                    proteinID = line_analysis[1]
                    proteinID = str(regex_proteinID.findall(proteinID)).replace("'", "").strip('[]')
                
                    if name in self.dict_fasta.keys():
                        if self.dict_fasta[name] == proteinID:
                            new_arq_gff.write(line)

                except:
                    continue
                
            else:
                continue
        
        arq_gff.close()
        new_arq_gff.close()


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
            if re.match('.+gff$', arg):
                OUTPUT = arg
            else:
                print("-o isn't a gff file output")
                sys.exit(5)

    if FASTA_PATH and GFF_PATH and OUTPUT:
        SELECTING.reading_fasta(FASTA_PATH)
        SELECTING.comparisons(GFF_PATH, OUTPUT)

    else:
        print('missing arguments -f, -g or -o')
        sys.exit(6)
