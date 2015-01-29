#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Standard library imports
import ConfigParser
import optparse
import sys
import gzip

# Third party imports
import HTSeq

# Local imports
from pyDNA.Utilities import is_readable_file, rm_blank, mkdir

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
        print("INITIALIZE QUADE")

        # Define fundamental variables
        self.version = "quade 0.1"
        self.usage = "Usage: %prog -c Conf.txt"

        # Use Conf file if interactive interpreter else parse arguments with optparse
        self.conf = conf_file if conf_file else self._optparser()

        # Parse the configuration file and verify the values of variables
        try:
            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=True)
            cp.read(self.conf)

            # Quality section
            self.mean_qual = cp.getint("quality", "mean_qual")
            self.minimal_qual = cp.getint("quality", "minimal_qual")
            self.qual_scale = cp.get("quality", "qual_scale")

            # Boolean flag of subindex presence
            self.index1 = True
            self.index2 = cp.getboolean("index", "index2")
            self.molecular1 = self.index1 and cp.getboolean("index", "molecular1")
            self.molecular2 = self.index2 and cp.getboolean("index", "molecular2")

            # Positions of subindex
            self.index1_pos = [cp.getint("index", "index1_start"), cp.getint("index", "index1_end")] if self.index1 else [0,0]
            self.index2_pos = [cp.getint("index", "index2_start"), cp.getint("index", "index2_end")] if self.index2 else [0,0]
            self.molecular1_pos = [cp.getint("index", "molecular1_start"), cp.getint("index", "molecular1_end")] if self.molecular1 else [0,0]
            self.molecular2_pos = [cp.getint("index", "molecular2_start"), cp.getint("index", "molecular2_end")] if self.molecular2 else [0,0]

            # List of fastqfiles
            self.seq_R1 = cp.get("fastq", "seq_R1").split()
            self.seq_R2 = cp.get("fastq", "seq_R2").split()
            self.index_R1 = cp.get("fastq", "index_R1").split() if self.index1 else []
            self.index_R2 = cp.get("fastq", "index_R2").split() if self.index2 else []

            # Samples are a special case since the number of sections is variable
            # Iterate only on sections starting by "sample"
            for sample in [i for i in cp.sections() if i.startswith("sample")]:

                # Fuse index if needed
                index_seq = cp.get(sample, "index1_seq")+(cp.get(sample, "index2_seq") if self.index2 else "")
                # Create a autoreferenced Sample_identifier object
                Sample_identifier(name = cp.get(sample, "name"), index = index_seq)

            # Values are tested in a private function
            self._test_values()

        # Handle the many possible errors occurring during conf file parsing or variable test
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as E:
            print ("Option or section missing. Report to the template configuration file\n" + E.message)
            sys.exit(1)
        except (ValueError, AssertionError) as E:
            print ("One of the value in the configuration file is not correct\n" + E.message)
            sys.exit(1)

    def __repr__(self):
        msg = "QUADE CLASS\n\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        for sample in Sample_identifier.Instances:
            msg+="\t{}\n".format(sample.__repr__())
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        Main function of the script
        """
        print("READ FASTQ FILES")
        in_list = self.reader()
        print("PROCESS INDEX AND IDENTIFY SAMPLES")
        out_list = self.worker(in_list)
        print("WRITE DEMULTIPLEXED FASTQ")
        self.writer(out_list)

        #for s1, s2, sample_name, quality_passed in out_list:
            #if not quality_passed and sample_name != "Undetermined":
                #print ("{}\t {}...".format(s1.name, s1.seq))
                #print ("{}\t {}...".format(s2.name, s2.seq))
                #print ("{}\t {}\n".format(sample_name, quality_passed))

    def reader (self):
        """
        Read fatsq files
        """
        # Init empty list to store read sets
        in_list = []

        # For double indexing
        if self.index2:
            # Iterate over fastq chunks for sequence and index reads
            for S1, S2, I1, I2 in zip(self.seq_R1, self.seq_R2, self.index_R1, self.index_R2):

                # Init HTSeq fastq reader generator
                S1_gen = HTSeq.FastqReader(S1, qual_scale=self.qual_scale)
                S2_gen = HTSeq.FastqReader(S2, qual_scale=self.qual_scale)
                I1_gen = HTSeq.FastqReader(I1, qual_scale=self.qual_scale)
                I2_gen = HTSeq.FastqReader(I2, qual_scale=self.qual_scale)

                # Iterate over read in fastq files and put each set in a list
                for s1, s2, i1, i2 in zip (S1_gen, S2_gen, I1_gen, I2_gen):
                    in_list.append([s1, s2 ,i1 ,i2])

        # For simple indexing
        else:
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

    def worker (self, in_list):
        """
        filter and asign sample to a set of fastq reads
        """
        # Init empty list to store read sets
        out_list = []

        # For double indexing
        if self.index2:
            for s1, s2, i1, i2 in in_list:

                # Extract index and molecular sequences from index reads and merge them
                index = i1[self.index1_pos[0]:self.index1_pos[1]]+i2[self.index2_pos[0]:self.index2_pos[1]]
                molecular = i1[self.molecular1_pos[0]:self.molecular1_pos[1]]+i2[self.molecular2_pos[0]:self.molecular2_pos[1]]

                # Identify sample correspondance and verify index quality
                sample_name = Sample_identifier.index_coresp(index)
                if not self._quality_filter(i1) or not self._quality_filter(i2):
                    sample_name+="_fail_qual"

                # Suffix read name with index and molecular sequence
                suffix = ":"+index+":"+molecular+":" if molecular else ":"+index+":"
                s1.name += suffix
                s2.name += suffix

                out_list.append([s1, s2, sample_name])

        # For simple indexing
        else:
            for s1, s2, i1 in in_list:

                # Extract index and molecular sequences from index read
                index = i1[self.index1_pos[0]:self.index1_pos[1]]
                molecular = i1[self.molecular1_pos[0]:self.molecular1_pos[1]]

                # Identify sample correspondance and verify index quality
                sample_name = Sample_identifier.index_coresp(index)
                if not self._quality_filter(i1):
                    sample_name+="_fail_qual"

                # Suffix read name with index and molecular sequence
                suffix = ":"+index+":"+molecular+":" if molecular else ":"+index+":"
                s1.name += suffix
                s2.name += suffix

                out_list.append([s1, s2, sample_name])

        return out_list

    def writer (self, out_list):
        """
        """

        for s1, s2, sample_name in out_list:
            with gzip.open ("./{}_R1.fastq.gz".format(sample_name), 'ab') as fastq_file:
                fastq_file.write("@{}\n{}\n+\n{}\n".format(s1.name, s1.seq, s1.qualstr))
            with gzip.open ("./{}_R2.fastq.gz".format(sample_name), 'ab') as fastq_file:
                fastq_file.write("@{}\n{}\n+\n{}\n".format(s2.name, s2.seq, s2.qualstr))


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
            is_readable_file (self.seq_R1 + self.seq_R2 + self.index_R1 + self.index_R2)
        else:
            assert len(self.seq_R1) == len(self.seq_R2) == len(self.index_R1) > 0,\
            "seq_R1, seq_R2 and index_R1 are mandatory and have to contain the same number of files"
            is_readable_file (self.seq_R1 + self.seq_R2 + self.index_R1)

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


    def _quality_filter(self, fastq):

        # test if all base are > minimal quality
        for val in fastq.qual:
            if val < self.minimal_qual:
                return False

        # test if index mean qual > mean qual
        if sum(fastq.qual)/len(fastq) < self.mean_qual:
            return False

        return True


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sample_identifier(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    Instances = [] # Class field used for instance tracking
    Index_seq_to_name = {}
    DNA = ["A","T","C","G"]

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def index_coresp (self, index):
        return Index_seq_access.get(index, "Undetermined")

    @ classmethod
    def all_get (self, key):
        # Return a list of a self variable from all sample objects
        return [sample[key] for sample in self.Instances]

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, index):

        # Create self variables
        self.name = name
        self.index = index.upper()

        # test uniqueness of name and index and DNA composition of index
        assert self.name not in self.all_get("name"), "{} : Name is not unique".format(self.name)
        assert self.index not in self.all_get("index"), "{} : Index is not unique".format(self.name)
        assert self._is_dna(index), "{} : Non canonical DNA base in index".format(self.name)

        # Append sample object to the class list of Instance
        self.Instances.append(self)

        # Add sample object to a class dict with index seq as a key
        self.Index_seq_to_name[self.index] = self.name

    # Fundamental class functions str and repr
    def __repr__(self):
        return "Sample ID:{}\tName : {}\tINDEX : {}".format(self.ID, self.name, self.index)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    # Function allowing to acces objet self values with key (ex sample[id])
    def __getitem__(self, key):
        return self.__dict__[key]

    #~~~~~~~PRIVATE METHODS~~~~~~~#
    def _is_dna (self, sequence):
        for base in sequence:
            if base not in self.DNA:
                return False
        return True

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#class Sample_aggregator(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Will received sample_name id = responssible for bufferized writting of fastq


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    quade = Quade()
    quade()
