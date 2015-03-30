#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    Quade
@brief      Main file of the program
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

try:
    # Standard library imports
    import ConfigParser
    import optparse
    import sys
    import os
    import gzip
    from time import time
    from datetime import datetime

    # Third party imports
    import numpy
    import HTSeq

    # Local imports
    from Sample import Sample
    from Conf_file import write_example_conf

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Quade(object):
    """
    Fastq file demultiplexer, handling double indexing, molecular indexing and filtering based
    on index quality
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "Quade 0.2"
    USAGE = "Usage: %prog -c Conf.txt [-i -h]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        # Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)
        optparser.add_option('-c', dest="conf_file",
            help= "Path to the configuration file [Mandatory]")
        optparser.add_option('-i', dest="init_conf", action='store_true',
            help= "Generate an example configuration file and exit [Facultative]")

        # Parse arguments
        options, args = optparser.parse_args()

        return Quade(options.conf_file, options.init_conf)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, conf_file=None, init_conf=None):
        """
        Initialization function, parse options from configuration file and verify their values.
        All self.variables are initialized explicitly in init.
        """

        # Create a example conf file if needed
        if init_conf:
            print("Create an example configuration file in the current folder")
            write_example_conf()
            sys.exit(0)

        print("Initialize Quade")
        # Parse the configuration file and verify the values of variables
        try:

            #verify if conf file was
            assert conf_file, "A path to the configuration file is mandatory"
            self.conf = conf_file

            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=True)
            cp.read(self.conf)

            # Quality section
            self.minimal_qual = cp.getint("quality", "minimal_qual")
            self.qual_scale = cp.get("quality", "qual_scale")

            # Boolean flag of subindex presence
            self.idx1 = True
            self.idx2 = cp.getboolean("index", "index2")
            self.mol1 = self.idx1 and cp.getboolean("index", "molecular1")
            self.mol2 = self.idx2 and cp.getboolean("index", "molecular2")

            # Positions of subindex
            self.idx1_pos = {
                "start":cp.getint("index", "index1_start")-1,
                "end":cp.getint("index", "index1_end")}
            self.idx2_pos = {
                "start":0 if not self.idx2 else cp.getint("index", "index2_start")-1,
                "end":0 if not self.idx2 else cp.getint("index", "index2_end")}
            self.mol1_pos = {
                "start":0 if not self.mol1 else cp.getint("index", "molecular1_start")-1,
                "end":0 if not self.mol1 else cp.getint("index", "molecular1_end")}
            self.mol2_pos = {
                "start":0 if not self.mol2 else cp.getint("index", "molecular2_start")-1,
                "end":0 if not self.mol2 else cp.getint("index", "molecular2_end")}

            # List of fastq files
            self.seq_R1 = cp.get("fastq", "seq_R1").split()
            self.seq_R2 = cp.get("fastq", "seq_R2").split()
            self.index_R1 = cp.get("fastq", "index_R1").split()
            self.index_R2 = [] if not self.idx2 else cp.get("fastq", "index_R2").split()

            # Init the class Sample with generic values (unconventional initialization)
            Sample.CLASS_INIT(
                write_undetermined = cp.getboolean("output", "write_undetermined"),
                write_pass = cp.getboolean("output", "write_pass"),
                write_fail = cp.getboolean("output", "write_fail"),
                min_qual = self.minimal_qual)

            # Samples are a special case, since the number of sections is variable
            # Iterate only on sections starting by "sample"
            for sample in [i for i in cp.sections() if i.startswith("sample")]:

                # Create a autoreferenced Sample object and fuse index if needed
                if self.idx2:
                    Sample(name=cp.get(sample, "name"), index=(cp.get(sample, "index1_seq")+cp.get(sample, "index2_seq")))
                else:
                    Sample(name=cp.get(sample, "name"), index_seq=cp.get(sample, "index1_seq"))

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
        return (msg)

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """ Main function of the script """

        start_time = time()

        print ("Start parsing files")
        # For double indexing
        if self.idx2:
            self.double_index_parser()
        # For simple indexing
        else:
            self.simple_index_parser()

        # Flush remaining content in sample buffers
        Sample.FLUSH_ALL()

        # Write a report
        print ("Generate_a csv report")
        with open ("Quade_report.csv", "wb") as report:
            report.write ("Program {}\tDate {}\n\n".format(self.VERSION,str(datetime.today())))
            for descr, value in Sample.REPORT():
                report.write ("{}\t{}\n".format(descr, value))

        print ("Done in {}s".format(round(time()-start_time, 3)))
        return(0)

    def double_index_parser (self):

        # Iterate over fastq chunks for sequence and index reads
        for n, (R1, R2, I1, I2) in enumerate (zip (self.seq_R1, self.seq_R2, self.index_R1, self.index_R2)):

            print("Parsing chunk {}/{}".format(n+1, len(self.seq_R1)))

            # Init HTSeq fastq reader generator
            R1_gen = HTSeq.FastqReader(R1, qual_scale=self.qual_scale)
            R2_gen = HTSeq.FastqReader(R2, qual_scale=self.qual_scale)
            I1_gen = HTSeq.FastqReader(I1, qual_scale=self.qual_scale)
            I2_gen = HTSeq.FastqReader(I2, qual_scale=self.qual_scale)

            # Iterate over read in fastq files and put each set in a list
            for read1, read2, index1, index2 in zip (R1_gen, R2_gen, I1_gen, I2_gen):

                # Extract index and molecular sequences from index reads and merge
                index = HTSeq.SequenceWithQualities(
                    seq = index1.seq[self.idx1_pos["start"]:self.idx1_pos["end"]] + index2.seq[self.idx2_pos["start"]:self.idx2_pos["end"]],
                    name = "fused_index",
                    qualstr = index1.qualstr[self.idx1_pos["start"]:self.idx1_pos["end"]] + index2.qualstr[self.idx2_pos["start"]:self.idx2_pos["end"]],
                    qualscale = self.qual_scale)

                molecular = HTSeq.SequenceWithQualities(
                    seq = index1.seq[self.mol1_pos["start"]:self.mol1_pos["end"]] + index2.seq[self.mol2_pos["start"]:self.mol2_pos["end"]],
                    name = "fused_molecular",
                    qualstr = index1.qualstr[self.mol1_pos["start"]:self.mol1_pos["end"]] + index2.qualstr[self.mol2_pos["start"]:self.mol2_pos["end"]],
                    qualscale = self.qual_scale)

                # Identify sample and verify index quality
                Sample.FINDER (read1,read2, index, molecular)

    def simple_index_parser (self):

        # Iterate over fastq chunks for sequence and index reads
        for n, (R1, R2, I1) in enumerate (zip(self.seq_R1, self.seq_R2, self.index_R1)):

            print("Parsing chunk {}/{}".format(n+1, len(self.seq_R1)))

            # Init HTSeq fastq reader generator
            R1_gen = HTSeq.FastqReader(R1, qual_scale=self.qual_scale)
            R2_gen = HTSeq.FastqReader(R2, qual_scale=self.qual_scale)
            I1_gen = HTSeq.FastqReader(I1, qual_scale=self.qual_scale)

            # Iterate over read in fastq files and put each set in a list
            for read1, read2 ,index1 in zip (R1_gen, R2_gen, I1_gen):

                # Extract index and molecular sequences from index reads
                index = HTSeq.SequenceWithQualities(
                    seq = index1.seq[self.idx1_pos["start"]:self.idx1_pos["end"]],
                    name = "index",
                    qualstr = index1.qualstr[self.idx1_pos["start"]:self.idx1_pos["end"]],
                    qualscale = self.qual_scale)

                molecular = HTSeq.SequenceWithQualities(
                    seq = index1.seq[self.mol1_pos["start"]:self.mol1_pos["end"]],
                    name = "molecular",
                    qualstr = index1.qualstr[self.mol1_pos["start"]:self.mol1_pos["end"]],
                    qualscale = self.qual_scale)

                # Identify sample and verify index quality
                Sample.FINDER (read1,read2, index, molecular)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _test_values(self):
        """ Test the validity of options in the configuration file """

        # Verify values from the quality section
        assert self.qual_scale in ["solexa", "solexa-old", 'phred'], "Authorized values for quality_scale : solexa, solexa-old, phred"
        assert 0 <= self.minimal_qual <= 40, "Authorized values for minimal_qual : 0 to 40"

        # Verify values of the fastq section and readability of fastq files
        if self.idx2:
            assert len(self.seq_R1) == len(self.seq_R2) == len(self.index_R1) == len(self.index_R2) > 0,\
            "seq_R1, seq_R2, index_R1 and index_R2 are mandatory and have to contain the same number of files"
            self._is_readable_file (self.seq_R1 + self.seq_R2 + self.index_R1 + self.index_R2)
        else:
            assert len(self.seq_R1) == len(self.seq_R2) == len(self.index_R1) > 0,\
            "seq_R1, seq_R2 and index_R1 are mandatory and have to contain the same number of files"
            self._is_readable_file (self.seq_R1 + self.seq_R2 + self.index_R1)

        # Verify values of the index section
        for pos in [self.idx1_pos, self.idx2_pos, self.mol1_pos, self.mol2_pos]:
            assert pos["start"] >= 0
            assert pos["end"] >= pos["start"]


    def _is_readable_file (self, file_list):
        """ Verify the readability of a file or list of file """

        for fp in file_list:
            if not os.access(fp, os.R_OK):
                raise ValueError ("{} is not a valid file".format(fp))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    quade = Quade.class_init()
    quade()
