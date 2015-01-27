#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Standard library imports
import ConfigParser
import optparse

# Third party imports
import HTSeq

# Local imports
from pyDNA.Utilities import is_readable_file

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Quade(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#
    
    def __init__ (self, conf_file = None):
        
        """
        Initialization function, parse options in command line and configuration file and verify
        their values. All self.variables are initialized explicitly in init/
        conf_file arg is mandatory in interactive interpreter and import, else optparse will parse
        command line arguments
        """

        # Define fundamental variables
        self.version = "quade 0.1"
        self.usage = "Usage: %prog -c Conf.txt"

        # Use Conf file if interactive interpreter else parse arguments with optparse
        if conf_file:
            self.conf = conf_file
        else:
            self.conf = self._optparser()
        
        # Parse the configuration file and verify the values of variables
        try:
            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=True)
            cp.read(self.conf)

            # Get Values from the conf file
            self.outdir = cp.get("general", "outdir") if cp.get("general", "outdir") else "./"
            self.mean_qual = cp.getint("quality", "mean_qual")
            self.minimal_qual = cp.getint("quality", "minimal_qual")
            self.qual_scale = cp.get("quality", "qual_scale")
            self.seq_R1 = cp.get("fastq", "seq_R1").split()
            self.seq_R2 = cp.get("fastq", "seq_R2").split()
            self.index_R1 = cp.get("fastq", "index_R1").split()
            self.index_R2 = cp.get("fastq", "index_R2").split()

            # For each index verify its existence and store start and end positions if needed
            self.index1_start = cp.getint("index", "index1_start")
            self.index1_end = cp.getint("index", "index1_end")
            
            self.index2 = cp.getboolean("index", "index2")
            self.molecular1 = cp.getboolean("index", "molecular1")
            self.molecular2 = cp.getboolean("index", "molecular2")
            if self.index2:
                self.index2_start = cp.getint("index", "index2_start")
                self.index2_end = cp.getint("index", "index2_end")
            if self.molecular1:
                self.molecular1_start = cp.getint("index", "molecular1_start")
                self.molecular1_end = cp.getint("index", "molecular1_end")
            if self.molecular2:
                self.molecular2_start = cp.getint("index", "molecular2_start")
                self.molecular2_end = cp.getint("index", "molecular2_end")

            # Samples are a special case since the number of sections is variable
            # Iterate only on sections starting by "sample" and create a autoreferenced Sample object
            for sample in [i for i in cp.sections() if i.startswith("sample")]:
                if self.index2:
                    Sample(name = cp.get(sample, "name"), index1 = cp.get(sample, "index1_seq"), index2 = cp.get(sample, "index2_seq"))
                else:
                    Sample(name = cp.get(sample, "name"), index1 = cp.get(sample, "index1_seq"))

            # Values are tested in a private function
            self._test_values()

        # Handle the many possible errors occurring during conf file parsing or variable test
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as E:
            print ("Option or section missing. Report to the template configuration file\n" + E.message)
            exit(1)
        except (ValueError, AssertionError) as E:
            print ("One of the value in the configuration file is not correct\n" + E.message)
            exit(1)

    def __repr__(self):
        msg = "QUADE CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        Main function of the script
        """
        print(repr(self))
        
        if self.index2:
            in_list = self.reader_double()
            #out_list = self.worker_double(in_list)
        else:
            in_list = self.reader_simple()
            #out_list = self.worker_simple(in_list)
        
        #self.writer(out_list)
        
        for read_set in in_list:
            for read in read_set:
                print read.seq
            print("")
    
    def reader_double (self):
        """
        Read fatsq files in case of double indexing = 4 files
        """
        # Init empty list to store read sets
        in_list = []
        
        # Iterate over fastq chunks for sequence and index reads
        for S1, S2, I1, I2 in zip(self.seq_R1, self.seq_R2, self.index_R1, self.index_R2):
            
            # Init HTSeq fastq reader generator
            S1_gen = HTSeq.FastqReader(S1, qual_scale=self.qual_scale)
            S2_gen = HTSeq.FastqReader(S2, qual_scale=self.qual_scale)
            I1_gen = HTSeq.FastqReader(I1, qual_scale=self.qual_scale)
            I2_gen = HTSeq.FastqReader(I2, qual_scale=self.qual_scale)
            
            # Iterate over read in fastq files and put each set in a list
            for s1, s2 ,i1 ,i2 in zip (S1_gen, S2_gen, I1_gen, I2_gen):
                in_list.append([s1, s2 ,i1 ,i2])
                
        return in_list
    
    def reader_simple (self):
        """
        Read fatsq files in case of simple indexing = 3 files
        """
        # Init empty list to store read sets
        in_list = []
        
        # Iterate over fastq chunks for sequence and index reads
        for S1, S2, I1 in zip(self.seq_R1, self.seq_R2, self.index_R1):
            
            # Init HTSeq fastq reader generator
            S1_gen = HTSeq.FastqReader(S1, qual_scale=self.qual_scale)
            S2_gen = HTSeq.FastqReader(S2, qual_scale=self.qual_scale)
            I1_gen = HTSeq.FastqReader(I1, qual_scale=self.qual_scale)
            
            # Iterate over read in fastq files and put each set in a list
            for s1, s2 ,i1 in zip (S1_gen, S2_gen, I1_gen):
                in_list.append([s1, s2 ,i1])
                
        return in_list
        

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _optparser(self):
        """
        Parse CLI and return a valid configuration file path
        """
        # Define parser usage, options
        optparser = optparse.OptionParser(usage = self.usage, version = self.version)
        optparser.add_option('-c', dest="conf_file", help= "Path to the configuration file")

        # Parse arguments
        options, args = optparser.parse_args()

        # Verify conf file
        if not options.conf_file:
            optparser.print_help()
            optparser.error("incorrect number of arguments")

        return options.conf_file

    def _test_values(self):
        """
        Test the validity of options in the configuration file
        """
        # Verify values from the quality section
        assert self.qual_scale in ["solexa-old", "solexa-old", 'phred'], "Authorized values for quality_scale : solexa-old, solexa-old, phred"
        assert 0 <= self.mean_qual <= 40, "Authorized values for mean_qual : 0 to 40"
        assert 0 <= self.minimal_qual <= 40, "Authorized values for minimal_qual : 0 to 40"

        # Verify values of the fastq section and readability of fastq files
        if self.index2:
            assert len(self.seq_R1) == len(self.seq_R2) == len(self.index_R1) == len(self.index_R2) > 0,\
            "seq_R1, seq_R2, index_R1 and index_R2 are mandatory and have to contain the same number of files"
            self._verify_file (self.seq_R1 + self.seq_R2 + self.index_R1 + self.index_R2)
        else:
            assert len(self.seq_R1) == len(self.seq_R2) == len(self.index_R1) > 0,\
            "seq_R1, seq_R2 and index_R1 are mandatory and have to contain the same number of files"
            self._verify_file (self.seq_R1 + self.seq_R2 + self.index_R1)

        # Verify values of the index section
        assert self.index1_start >= 0, "Autorized values for index1_start : >= 0"
        assert self.index1_end >= self.index1_start, "Autorized values for index1_end : > index1_start"
        if self.index2:
            assert self.index2_start >= 0, "Autorized values for index2_start : >= 0"
            assert self.index2_end >= self.index2_start, "Autorized values for index2_end : > index2_start"
        if self.molecular1:
            assert self.molecular1_start >= 0, "Autorized values for molecular1_start : >= 0"
            assert self.molecular1_end >= self.molecular1_start, "Autorized values for molecular_end : > molecular1_start"
        if self.molecular2:
            assert self.molecular2_start >= 0, "Autorized values for molecular2_start : >= 0"
            assert self.molecular2_end >= self.molecular2_start, "Autorized values for molecular2_end : > molecular2_start"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sample(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    Instances = [] # Class field used for instance tracking
    ID_count = 1
    DNA = ["A","T","C","G"]

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def next_ID (self):
        current = self.ID_count
        self.ID_count +=1
        return current
    
    @ classmethod
    def index_coresp (self, index1, index2=""):
        pass
        # Return the sample name if index corespond to a sample 
        # Else return Undetermined
    
    @ classmethod
    def all_get (self, key):
        # Return a list of a self variable from all sample objects
        return [sample[key] for sample in self.Instances]

    def __init__ (self, name, index1, index2=""):
        
        # Store object variables
        print ("Creating Sample object {}".format(name))
               
        # Name uniqness
        assert name not in self.all_get("name"), "{} : Name is not unique".format(name)
        self.name = name
        
        # Index combination uniqueness is verified and DNA sequence
        index = (index1+index2).upper()
        assert index not in self.all_get("index"), "{} : Index is not unique".format(name)
        assert not [base for base in index if base not in self.DNA],\
            "{} : Non canonical DNA base in index".format(name)
        self.index = index
        
        # Attibute an id and append sample object to the Instance list
        self.ID = self.next_ID()
        self.Instances.append(self)
        
    # Fundamental class functions str and repr
    def __repr__(self):
        msg = "SAMPLE CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
    
    def __getitem__(self, key):
        return self.__dict__[key]
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    quade = Quade()
    quade()
