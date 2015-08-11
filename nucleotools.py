"""A module for handling simple string methods involving nucleotides in bioinformatics problems"""
import re
class DNACodonDictionary(dict):
    """A class which acts as a dictionary with DNA codons acting as a key for the amino acids produced upon translation"""
    def __init__(self):
        self.__dict__ = {}
        ##T__
        self.__dict__["TTT"] = "F"
        self.__dict__["TTC"] = "F"
        self.__dict__["TTA"] = "L"
        self.__dict__["TTG"] = "L"
        
        self.__dict__["TCT"] = "S"
        self.__dict__["TCC"] = "S"
        self.__dict__["TCA"] = "S"
        self.__dict__["TCG"] = "S"

        self.__dict__["TAT"] = "Y"
        self.__dict__["TAC"] = "Y"
        self.__dict__["TAA"] = False
        self.__dict__["TAG"] = False

        self.__dict__["TGT"] = "C"
        self.__dict__["TGC"] = "C"
        self.__dict__["TGA"] = False
        self.__dict__["TGG"] = "W"
        ##C__
        self.__dict__["CTT"] = "L"
        self.__dict__["CTC"] = "L"
        self.__dict__["CTA"] = "L"
        self.__dict__["CTG"] = "L"

        self.__dict__["CCT"] = "P"
        self.__dict__["CCC"] = "P"
        self.__dict__["CCA"] = "P"
        self.__dict__["CCG"] = "P"

        self.__dict__["CAT"] = "H"
        self.__dict__["CAC"] = "H"
        self.__dict__["CAA"] = "Q"
        self.__dict__["CAG"] = "Q"

        self.__dict__["CGT"] = "R"
        self.__dict__["CGC"] = "R"
        self.__dict__["CGA"] = "R"
        self.__dict__["CGG"] = "R"
        ##A__
        self.__dict__["ATT"] = "I"
        self.__dict__["ATC"] = "I"
        self.__dict__["ATA"] = "I"
        self.__dict__["ATG"] = "M"

        self.__dict__["ACT"] = "T"
        self.__dict__["ACC"] = "T"
        self.__dict__["ACA"] = "T"
        self.__dict__["ACG"] = "T"

        self.__dict__["AAT"] = "N"
        self.__dict__["AAC"] = "N"
        self.__dict__["AAA"] = "K"
        self.__dict__["AAG"] = "K"

        self.__dict__["AGT"] = "S"
        self.__dict__["AGC"] = "S"
        self.__dict__["AGA"] = "R"
        self.__dict__["AGG"] = "R"
        ##G__
        self.__dict__["GTT"] = "V"
        self.__dict__["GTC"] = "V"
        self.__dict__["GTA"] = "V"
        self.__dict__["GTG"] = "V"

        self.__dict__["GCT"] = "A"
        self.__dict__["GCC"] = "A"
        self.__dict__["GCA"] = "A"
        self.__dict__["GCG"] = "A"

        self.__dict__["GAT"] = "D"
        self.__dict__["GAC"] = "D"
        self.__dict__["GAA"] = "E"
        self.__dict__["GAG"] = "E"

        self.__dict__["GGT"] = "G"
        self.__dict__["GGC"] = "G"
        self.__dict__["GGA"] = "G"
        self.__dict__["GGG"] = "G"

    def __getitem__(self, key):
        return self.__dict__[key]

    def values(self):
        return self.__dict__.values()

    def __len__(self): 
        return len(self.__dict__)

    def __repr__(self): 
        return repr(self.__dict__)

    def __contains__(self, item):
        return item in self.__dict__

    def __iter__(self):
        return iter(self.__dict__)

class RNACodonDictionary(DNACodonDictionary):
    """A class which acts as a dictionary with RNA codons acting as a key for the amino acids produced upon translation"""
    def __init__(self):
        for key in DNACodonDictionary():
            self.__dict__[re.sub("T","U",key)] = DNACodonDictionary()[key]


class DNATools:
    "tools for translating, transcribing, complementing or reverse-complementing DNA strands"
    def __init__(self):
        self.complement = {}
        self.complement["A"]="T"
        self.complement["T"]="A"
        self.complement["C"]="G"
        self.complement["G"]="C"
        self.startCodon = re.compile("(?=(ATG))")
        self.tranDict = DNACodonDictionary()
        

    def __complement__(self, sequence):
        returny = ""
        for character in sequence:
            returny+=self.complement[character]
        return returny

    def reverseComplement(self, sequence):
        "returns the reverse complement of parameter sequence"
        return self.__complement__(reversed(sequence))

    def translateStrandNoFrame(self, sequence):
        "translates parameter sequence to polypeptide sequence from the beginning of parameter sequence until the end of parameter string or until a stop codon is reached"
        returnString = ""
        for i in range(3, len(sequence), 3):
            if not self.tranDict[sequence[i-3:i]]:
                break
            returnString+=self.tranDict[sequence[i-3:i]]
        return returnString

    def translateRCompNoFrame(self, sequence):
        "translates the reverse complement of parameter sequence to polypeptide sequence from the beginning of parameter sequence until the end of parameter string or until a stop codon is reached"
        return self.translateStrandNoFrame(self.reverseComplement(sequence))

    def translateStrandORF(self, sequence):
        "finds all polypeptide sequences produced from translating the parameter nucleotide strand beginning at a start codon and ending at a stop codon"
        starts = re.finditer(self.startCodon,sequence)
        aminoSequences = []
        for start in starts:
            aminoSequence = ""
            for i in range(start.start()+3,len(sequence),3):
                amino = self.tranDict[sequence[i-3:i]]
                if not amino:
                    aminoSequences.append(aminoSequence)
                    break
                aminoSequence+=amino
        return list(set(aminoSequences))

    def translateRCompORF(self, sequence):
        "finds all polypeptide sequences produced from translating the reverse complement of parameter nucleotide strand beginning at a start codon and ending at a stop codon"
        return self.translateStrandORF(self.reverseComplement(sequence))

    def translateAllORF(self, sequence):
        "finds all polypeptide sequences produced from translating the parameter nucleotide strand and its reverse complement beginning at a start codon and ending at a stop codon"
        returnList = []
        for aminoSeq in self.translateStrandORF(sequence):
            returnList.append(aminoSeq)
        for aminoSeq in self.translateRCompORF(sequence):
            returnList.append(aminoSeq)
        return list(set(returnList))

    def transcribe(self, sequence):
        "transcribes DNA sequence to equivalent RNA sequence"
        return re.sub("T","U",sequence)


class RNATools(DNATools):
    "tools for translating, transcribing, complementing or reverse-complementing DNA strands"
    def __init__(self):
        self.complement = {}
        self.complement["A"]="U"
        self.complement["U"]="A"
        self.complement["C"]="G"
        self.complement["G"]="C"
        self.startCodon = re.compile("(?=(AUG))")
        self.tranDict = RNACodonDictionary()

    def transcribe(self, sequence):
        "disabled in RNA tools"
        raise Exception('No such method for RNA tools')



