# Real-time PCR primer design
# define the criteria for primers in the class
import math
from Bio import pairwise2
from colorama import Fore, Style

class primerRequisites:
    
    def __init__(self, primer, length):
        self.primer = primer
        self.length = length
   # define GC content
    def primer_GC(self):
        GC_number = self.primer.count('G') + self.primer.count('C')
        GC_content = (GC_number / self.length) * 100
        if 40 < GC_content < 60:
            return True
  # define Tm      
    def primer_Tm(self):
        GC_number = self.primer.count('G') + self.primer.count('C')
        Tm = 64.9 + 41 * (GC_number -16.4) / self.length
        if  58 < Tm < 65:
            return True
  # define to avoid 3'GC clamp
    def primer_GCclamp(self):
        last5bases = self.primer[(self.length - 6): self.length]
        last5GC = last5bases.count('G') + last5bases.count('C')
        if last5GC < 3:
            return True
  # define repeatie nucleotides  
    def primer_nucleotiderepeat(self):
        if all (['AAAA' not in self.primer, 'TTTT' not in self.primer, 'CCCC' not in self.primer, 'GGGG' not in self.primer]):
            return True
    # define to avoid repeative di-nucleotides
    def primer_dinucleotidesrepeat(self):
        if any(['ATATATAT' in self.primer, 'ACACACAC' in self.primer, 'AGAGAGAG' in self.primer,
                'TATATATA' in self.primer, 'TCTCTCTC' in self.primer, 'TGTGTGTG' in self.primer,
                'CACACACA' in self.primer, 'CTCTCTCT' in self.primer, 'CGCGCGCG' in self.primer,
                'GAGAGAGA' in self.primer, 'GTGTGTGT' in self.primer, 'GCGCGCGC' in self.primer]):
            return False
        else:
            return True
    # define to remove primers has matches inside
    def homo_hairpin(self):
        revPrimer = self.primer[::-1]
        # get the reverse complementary sequence of the primer
        reverse_nucleotide = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        revComple = ''
        for i in range(self.length):
            revComple += reverse_nucleotide[revPrimer[i]]
        # do the pariwise alignment in primer itself
        alignments = pairwise2.align.globalxx(self.primer, revComple)
        # restrain the matched number less than 5
        if alignments[0][2] <= (self.length // 2):
            return True

# define second class to get the reverse complementary primer, and the second class inherit all attributes from first class
class primerParing(primerRequisites):
    
    def __init__(self, primer, length, primerpair):
        super().__init__(primer, length) # Alternatively: primerRequisites.__init__(self, primer, length)
        self.primerpair = primerpair
    # define to remove paired primers that have more that 2C Tm difference
    def primerTm_differ(self):
        GC_number1 = self.primer.count('G') + self.primer.count('C')
        GC_number2 = self.primerpair.count('G') + self.primerpair.count('C')
        Tm1 = 64.9 + 41 * (GC_number1 - 16.4) / self.length
        Tm2 = 64.9 + 41 * (GC_number2 - 16.4) / len(self.primerpair)
        if math.fabs(Tm1 - Tm2) < 2:
            return True
    # define to remove paired primers that have matches
    def heter_hairpin(self):
        # get the reverse complementary sequence
        reverse_nucleotide = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        revPrimer1 = self.primer[::-1]
        revPrimer2 = self.primerpair[::-1]
        revComple1 = ''
        revComple2 = ''
        for i in range(self.length):
            revComple1 += reverse_nucleotide[revPrimer1[i]]
        for i in range(len(self.primerpair)):
            revComple2 += reverse_nucleotide[revPrimer2[i]]
        # do the pariwise alignment between primer pairs
        alignment1 = pairwise2.align.globalxx(self.primer, revComple2)
        alignment2 = pairwise2.align.globalxx(revComple1, self.primerpair)
        # restrain the matched number less than half lenght of primers
        if alignment1[0][2] <= (self.length // 2) and alignment2[0][2] <= (len(self.primerpair) //2):
            return True

    
# Find all potential primers meeting all requisites for real-time PCR primers
def primerFinding (codingsequence):
    allprimers = []
    codinglength = len(codingsequence)
    for i in range(20,27):
        for j in range(0, codinglength):
            if j < (codinglength - i):
                primer = codingsequence[j: (j+i+1)]
                primer = ''.join(primer)
                a = primerRequisites(primer, i)
                if a.primer_GC() and a.primer_Tm() and a.primer_GCclamp() and a.primer_nucleotiderepeat() and a.primer_dinucleotidesrepeat() and a.homo_hairpin():
                    allprimers.append(primer)

   # paring primers
    # get the index of all primers on coding sequence
    allprimers_index = [codingsequence.index(i) for i in allprimers]
    # make a ditionary to link the index with primers
    allprimers_dict = dict(zip(allprimers_index, allprimers)) 
    # alternative method to make a dictionary between two lists: allprimers_dict = {k:v for k, v in zip(allprimers_index, allprimers)}
    
    # pair primers to let the amplicon size within 90-150
    allpairedprimers = []
    x = 0
    for i in range(len(allprimers_index)):
        for j in range(i+1, len(allprimers_index)):
            if 90 < (allprimers_index[j] - allprimers_index[i]) < 150:
                k = allprimers_dict[allprimers_index[i]]
                v = allprimers_dict[allprimers_index[j]]
                w = primerParing(k, len(k), v)
                # rule out all paired primers has more than 2C Tm difference
                if w.primerTm_differ() and w.heter_hairpin():
                    x += 1
                    oneprimer = Fore.RED + k + Style.RESET_ALL
                    theotherprimer = Fore.RED + v + Style.RESET_ALL
                    primerReve = v[::-1]
                    # get the reverse complementary sequence of the primer
                    reverse_nucleotide = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                    primerR = ''
                    for y in range(len(primerReve)):
                        primerR += reverse_nucleotide[primerReve[y]]
                    p1 = codingsequence.index(k) # label the first primer with red
                    p2 = codingsequence.index(v) # label the other primer with red
                    # print the whole sequence with primers in red color
                    print('primer pair:', x, ', from', p1, 'to', p2)
                    print(codingsequence[0:p1] + oneprimer + codingsequence[(p1+len(k)):p2] + theotherprimer + codingsequence[(p2+len(v)):codinglength])
                    # write the positions and forward and reverse primers to a text file
                    with open('all primers.txt', 'a') as primerlist:
                        primerlist.write('primer pair ' + str(x) + ':' + '\t' + str(p1) + '-->' + str(p2))
                        primerlist.write('\n' + 'Forward primer: ' + k + '\t'*2 + 'Reverse primer: ' + primerR + '\n')
                        primerlist.write('\n')

    primerlist.close()
# input a sequence
codingsequence = input("Input your DNA sequence here: ")
codingsequence = codingsequence.upper()
primerFinding(codingsequence)

# ****To do the blast in NCBI****
print('*'*100)
print()
print('....All forward and reverse primers have been outputed to \'all primer.txt\'....')
x = input('Would you want to do the blast for your primer: <Y or N, enter to quit>')
while x == 'Y' or x =='y':
    from Bio.Blast import NCBIWWW

    seq_for_blast = input('Enter your primer for Blast: ')
    result_handle = NCBIWWW.qblast('blastn', 'nt', seq_for_blast)

    with open('my_blast.xml', 'w') as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()

    result_handle = open('my_blast.xml')

    from Bio.Blast import NCBIXML
    blast_record = NCBIXML.read(result_handle)

    E_VALUE_THRESH = 0.01
    t = 0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                t += 1
                print('****Hits', t, '****')
                print('sequnce: ', alignment.title)
                print('length', alignment.length)
                print('e_evalue: ', hsp.expect)
                print(hsp.query[0:100] + '...')
                print(hsp.match[0:100] + '...')
                print(hsp.sbjct[0:100] + '...')
    print('*'*100)
    print('Blast is done!')
    print()
    x = input('Would you want to do another blast for your primer: <Y or N, enter to quit>')
else:
    print()
    print()
    print('Thank you for using the program. Have a good day!')
