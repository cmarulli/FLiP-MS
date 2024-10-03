import os
import csv

class Peptide:
    """Class for handling peptide data"""
    
    def __init__(self, infos):
        self.protein = infos[0]
        self.peptide = infos[1][1:-1]
        self.F = infos[2]
        self.p = infos[3]
        self.q = float(infos[4])
        
        self.MSError = float(infos[5])
        self.MSGroup = float(infos[6])
        self.mean100K = float(infos[7])
        self.mean50K = float(infos[8])
        self.mean30K = float(infos[9])
        self.mean10K = float(infos[10])
        self.tryptics = list()
        self.tryptics.append(infos[11][1:-1])
        self.L2FC = float(infos[12])
        
        self.effect_size = float(infos[13])
        self.difference_mean = infos[14]
        self.expected_at_interface = infos[15]
        if self.expected_at_interface=='TRUE':
            self.sign = True
        else:
            self.sign = False
        self.positions = list()
        self.positions.append(int(infos[16]))
        self.sequence = infos[17]
        self.fully_tryptic = infos[18]
    
    ###__repr__(self)####add this
    def show(self):
        print('protein: ', self.protein, 'peptide: ', self.peptide,
              'q: ', self.q, 'is expected at interface: ', self.expected_at_interface)
    

class PeptideFamily:
    """Class for handling family peptide data"""
    
    def __init__(self, peptide, q, protein, idxs):
        self.peptide = peptide
        self.q = q
        self.protein = protein
        self.idxs = idxs
    
    ###__repr__(self)####add this
    def show(self):
        print('peptide: ', self.peptide, 'q: ', self.q, 'protein: ', self.protein, 'peptides mapped: ', self.idxs)

########################################################################################################
#########################################  FUNCTIONS  ##################################################
########################################################################################################


def load_peptides(data_path):
    '''returns a list of objects peptide'''

    FLiPR_path = os.path.join(data_path, 'FLiPRFraction_anova_results.csv')
    with open(FLiPR_path) as file:
        csv_reader = csv.reader(file)
        peptides = list()
        skip = 0
        for row in csv_reader:
            if skip==0:
                skip=1
                continue
            peptides.append(Peptide(row[0].split(';')))

    print('loaded {} peptides'.format(len(peptides)))
    

    return peptides


def group_peptides(peptides):
    '''group the peptides based on the protein they belong to. Returns a dictionary uniprot_id:list_of_peptides'''

    clusters = dict()
    for idx,p in enumerate(peptides):
        if p.protein not in clusters.keys():
            clusters[p.protein]=list()
        clusters[p.protein].append(idx)

    print('number of proteins {}'.format(len(clusters.keys())))
    
    return clusters




def store_sequences(clusters, yeast_sequences_file):
    '''stores the sequences contained in clusters into a dictionary'''
    sequences = dict()
    sequence_id = ''
    with open(yeast_sequences_file, 'r') as file:
        fasta_file = file.readlines()
        for row in fasta_file:
            if row[0]=='>':
                sequence_id = row.split('|')[1]
                if sequence_id in clusters.keys():
                    sequences[sequence_id] = ''
                    
            else:
                if sequence_id in sequences.keys():
                    sequences[sequence_id]+=row.strip()
    
    print('stored {} sequences'.format(len(sequences)))
                
    return sequences


###FUNCTIONS FOR COMPUTING PROTEOME COVERAGE####


def find_all_peptide_positions(sequence, pep):
    positions = list()
    x1 = 0
    x2 = 0
    while x1!=-1:
        x1 = sequence[x2:].find(pep)
        if x1!=-1:
            x2 = x2+x1+1
            positions.append(x2)
                
    return positions


def add_positions(sequence, pep, positions_list):
    positions = find_all_peptide_positions(sequence, pep)
    for pos in positions:
        for x in range(pos, pos+len(pep)):
            positions_list.append(x)


def compute_coverage(clusters, sequences, peptides):
    '''create a dictionary with the fraction of each sequences observed within the peptides'''
    coverage = dict()

    for prot in clusters:
        positions_list = list()
        for pep_idx in clusters[prot]:
            add_positions(sequences[prot], peptides[pep_idx].peptide, positions_list)
    
        positions_list = set(positions_list)
        sequence_length = len(sequences[prot])
        coverage[prot] = len(positions_list)/sequence_length
    
    return coverage