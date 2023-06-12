# Gaurav Sablok
# Thulani Makhalanyane
# Senior Postdoctoral Fellow
# Faculty of Natural and Agricultural Sciences
# Room 7-35, Agricultural Sciences Building
# University of Pretoria, Private Bag X20
# Hatfield 0028, South Africa
"""
A python function for the normalization of the functional 
gene overlap and gene numbers
across a large number of genomes
Mathematical formula from there and i coded it
Reference : Jurdzinski et al., Sci. Adv. 9, eadg2059 (2023)
"""
import os
import pandas as pd
class KONormalize(self):
    def __init__(self, filename, genome1, genome2):
        self.file = filename
        self.path = os.path.join(os.getcwd(), self.filename)
        self.datafile = pd.read_csv(self.path, sep = ",")
        self.genome1 = genome1
        self.genome2 = genome2
        print(f'the path to the supplied file is {self.path}')
        print(f'the name of the supplied first genome is {self.genome1}')
        print(f'the name of the supplied second genome is {self.genome2}')
    def getcommonKO(self):
        """Find the common shared KO between the genomes"""
        self.datafile = pd.read_csv(self.path, sep = ",")
        self.KO1 = self.datafile["self.genome1"].tolist()
        self.KO2 = self.datafile["self.genome2"].tolist()
        self.commKO = [self.__KO1[i] for i in range(len(set(self.KO1))) if self.KO1[i] in set(self.KO2)]
        return self.commKO

    def getgeneOverlap(self):
        """Obtain the gene overlap between the two genomes"""
        self.datafile = pd.read_csv(self.path, sep = ",")
        self.KO1 = self.datafile["self.genome1"].tolist()
        self.KO2 = self.datafile["self.genome2"].tolist()
        self.shared_KO_length = len([self.KO1[i] for i in range(len(set(self.KO1))) \
                                              if self.KO1[i] in set(self.KO2)])
        self.complete_KO_length = (len(set(self.KO1))+len(set(self.KO2)))//2
        self.gene_overlap = self.shared_KO_length//self.complete_KO_length
        return self.gene_overlap

    def uniqueKO(self):
        """Obtain the unique KO between the two genomes
        self.commKO1 are those that are present in the first genome
        self.commKO2 are those that are present in the second genome
        """
        self.datafile = pd.read_csv(self.path, sep = ",")
        self.KO3 = self.datafile["self.genome1"].tolist()
        self.KO4 = self.datafile["self.genome2"].tolist()
        self.commKO1 = [self.KO3[i] for i in range(len(set(self.KO3))) \
                       if self.KO3[i] not in set(self.KO4)]
        self.commKO2 = [self.KO4[i] for i in range(len(set(self.KO4))) \
                       if self.KO4[i] not in set(self.KO3)]
        return self.commKO1, self.commKO2
