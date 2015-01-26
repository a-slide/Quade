#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Imports
import ConfigParser
import optparse
import os

# Definition of the main class
class Quade(object):

    # Class init method
    def __init__ (self): # conf_file arg is for interactive interpreter and import
        """
        Initialization function, parse options in command line and configuration file and verify
        their values. All self.variables are initialized explicitly in init
        """

        #### DEFINE FUNDAMENTAL VARIABLES ####

        self.version = "quade 0.1"
        self.usage = "Usage: %prog -c Conf.txt"

        #### PARSE COMMAND LINE ARGUMENTS ####

        self.conf = self._optparser()

        try:
            #### PARSE CONF FILE ####

            # Define a configuration file parser object and load the configuration file
            cp = ConfigParser.RawConfigParser(allow_no_value=True)
            cp.read(self.conf)

            # Get Values from the conf file
            self.outdir = cp.get("general", "outdir") if cp.get("general", "outdir") else "./"
            self.mean_qual = cp.getint("quality", "mean_qual")
            self.minimal_qual = cp.getint("quality", "minimal_qual")
            self.quality_scale = cp.get("quality", "quality_scale")
            self.seq_R1 = cp.get("fastq", "seq_R1").split()
            self.seq_R2 = cp.get("fastq", "seq_R2").split()
            self.index_R1 = cp.get("fastq", "index_R1").split()
            self.index_R2 = cp.get("fastq", "index_R2").split()

            # For each index verify its existence and store start and end positions if needed
            self.index1_start = cp.getint("index", "index1_start")
            self.index1_end = cp.getint("index", "index1_end")

            self.index2 = cp.getboolean("index", "index2")
            if self.index2:
                self.index2_start = cp.getint("index", "index2_start")
                self.index2_end = cp.getint("index", "index2_end")

            self.molecular1 = cp.getboolean("index", "molecular1")
            if self.molecular1:
                self.molecular1_start = cp.getint("index", "molecular1_start")
                self.molecular1_end = cp.getint("index", "molecular1_end")

            self.molecular2 = cp.getboolean("index", "molecular2")
            if self.molecular2:
                self.molecular2_start = cp.getint("index", "molecular2_start")
                self.molecular2_end = cp.getint("index", "molecular2_end")

            # Sample section(s) is(are) a special case since the number of sample sections is variable
            self.sample_dict = {}
            # Iterate only on sections starting by "sample"
            for sample in [i for i in cp.sections() if i.startswith("sample")]:

                # Extract sample name and fuse index together
                sample_name = cp.get(sample, "name")
                fused_index = (cp.get(sample, "index1_seq")+cp.get(sample, "index2_seq")).upper()

                # Verify if the sample name is uniq
                if sample_name in self.sample_dict:
                    raise ValueError("All sample name have to be uniq")
                else:
                    self.sample_dict[sample_name] = fused_index

            ### TEST VALUES

            self._test_values()

        #### HANDLE ERRORS ####

        # Handle the possible exceptions occurring during conf file parsing or testing variable values
        except (ConfigParser.NoOptionError, ConfigParser.NoSectionError) as E:
            print ("Option or section missing. Report to the template configuration file\n" + E.message)
            exit(1)
        except (ValueError, AssertionError) as E:
            print ("One of the value in the configuration file is not correct\n" + E.message)
            exit(1)

    # Add fundamental class functions str and repr
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


    ## Main function for future development (__call__ allow object to be callable)
    def __call__(self):
        print(repr(self))


        for R1,R2,R3 in zip(R1_list, R2_list, R3_list):
            R1_gen = HTSeq.FastqReader(R1, qual_scale="phred")
            R2_gen = HTSeq.FastqReader(R2, qual_scale="phred")
            R3_gen = HTSeq.FastqReader(R3, qual_scale="phred")
            R4_gen = HTSeq.FastqReader(R3, qual_scale="phred")
            for r1,r2,r3,r4 in zip (R1_gen, R2_gen, R3_gen, R4_gen):
                print (r1.name+"\n"+r2.name+"\n"+r3.name+"\n"+r4.name+"\n\n")





    ## Private methods

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
        assert self.quality_scale in ["solexa-old", "solexa-old", 'phred'], "Authorized values for quality_scale : solexa-old, solexa-old, phred"
        assert 0 <= self.mean_qual <= 40, "Authorized values for mean_qual : 0 to 40"
        assert 0 <= self.minimal_qual <= 40, "Authorized values for minimal_qual : 0 to 40"

        # Verify values of the fastq section and readability of fastq files
        if self.index2:
            assert len(self.seq_R1) == len(self.seq_R2) == len(self.index_R1) == len(self.index_R2) > 0,\
            "seq_R1, seq_R2, index_R1 and index_R2 are mandatory and have to contain the same number of files"
            self._verify_files (self.seq_R1 + self.seq_R2 + self.index_R1 + self.index_R2)
        else:
            assert len(self.seq_R1) == len(self.seq_R2) == len(self.index_R1) > 0,\
            "seq_R1, seq_R2 and index_R1 are mandatory and have to contain the same number of files"
            self._verify_files (self.seq_R1 + self.seq_R2 + self.index_R1)

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

        # Verify values of the sample sections
        for name, seq in self.sample_dict.items():
            for base in seq:
                if not base in ["A","T","C","G"]:
                    raise ValueError( "{} contains {} in index sequence".format(name, base))


    def _verify_files (self, file_list):
        """
        Verify the readability of a list of file
        """
        for fp in file_list :
            if not os.access(fp, os.R_OK):
                raise ValueError ("{} is not a valid file".format(fp))


#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

    quade = Quade()
    quade()
