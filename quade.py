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
            self.index1_pos = [cp.getint("index", "index1_start")-1, cp.getint("index", "index1_end")] if self.index1 else [0,0]
            self.index2_pos = [cp.getint("index", "index2_start")-1, cp.getint("index", "index2_end")] if self.index2 else [0,0]
            self.molecular1_pos = [cp.getint("index", "molecular1_start")-1, cp.getint("index", "molecular1_end")] if self.molecular1 else [0,0]
            self.molecular2_pos = [cp.getint("index", "molecular2_start")-1, cp.getint("index", "molecular2_end")] if self.molecular2 else [0,0]

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
                index = i1.seq[self.index1_pos[0]:self.index1_pos[1]]+i2.seq[self.index2_pos[0]:self.index2_pos[1]]
                molecular = i1.seq[self.molecular1_pos[0]:self.molecular1_pos[1]]+i2.seq[self.molecular2_pos[0]:self.molecular2_pos[1]]

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
                index = i1.seq[self.index1_pos[0]:self.index1_pos[1]]
                molecular = i1.seq[self.molecular1_pos[0]:self.molecular1_pos[1]]

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
        Simple and naive version of fastq writer
        """
        sample_dict = {}
        for s1, s2, sample_name in out_list:

            # Use of a exception handling to limit the number of test
            try:
                # Add new reads to the sample
                sample_dict[sample_name](s1, s2)

            except KeyError:
                # Create a new sample_aggretor  and add new reads to the sample
                sample_dict[sample_name] = Sample_writer(sample_name)
                sample_dict[sample_name](s1, s2)

        for sample in sample_dict.values():
            sample.flush()

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
        for start, end in [self.index1_pos, self.index2_pos, self.molecular1_pos, self.molecular2_pos]:
            assert start >= 0
            assert end >= start

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
        return self.Index_seq_to_name.get(index, "Undetermined")

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
        return "SAMPLE\tName : {}\tINDEX : {}".format(self.name, self.index)

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
class Sample_writer(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__ (self, name):

        # Define self variables
        self.name = name

        # Counters
        self.total_pair = 0
        self.current_pair = 0

        # Str aggregators
        self.R1_buffer = ""
        self.R2_buffer = ""

        # File names
        self.R1_file = "./{}_R1.fastq.gz".format(self.name)
        self.R2_file = "./{}_R2.fastq.gz".format(self.name)

        # Init empty files
        with gzip.open (self.R1_file, 'wb'):
            pass
        with gzip.open (self.R2_file, 'wb'):
            pass


    # Fundamental class functions str and repr
    #def __repr__(self):
        #return "SAMPLE\tName : {}\tINDEX : {}".format(self.name, self.index)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    # Function allowing to acces objet self values with key (ex sample[id])
    def __getitem__(self, key):
        return self.__dict__[key]

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self, read1, read2):

        self.R1_buffer += "@{}\n{}\n+\n{}\n".format(read1.name, read1.seq, read1.qualstr)
        self.R2_buffer += "@{}\n{}\n+\n{}\n".format(read2.name, read2.seq, read2.qualstr)

        self.total_pair += 1
        self.current_pair += 1

        # If buffers contains more that 20 sequences, write in file and reset the counter
        if self.current_pair >= 20:
            self.flush()
            self.current_pair = 0

    def flush(self):

        # Append the content of the buffers in file and reset the buffers
        with gzip.open (self.R1_file, 'ab') as fastq_file:
            fastq_file.write(self.R1_buffer)
            self.R1_buffer = ""

        with gzip.open (self.R2_file, 'ab') as fastq_file:
            fastq_file.write(self.R2_buffer)
            self.R2_buffer = ""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    quade = Quade()
    quade()
