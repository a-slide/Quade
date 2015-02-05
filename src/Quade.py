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

    # Third party imports
    import numpy
    import HTSeq

except ImportError as E:
    print (E)
    print ("Please verify your dependencies. See Readme for more informations\n")
    exit()


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
        print("Initialize quade")

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

            # Output option sections
            self.write_undetermined = cp.getboolean("output", "write_undetermined")
            self.write_fail = cp.getboolean("output", "write_fail")
            self.write_report = cp.getboolean("output", "write_report")

            # Samples are a special case, since the number of sections is variable
            # Iterate only on sections starting by "sample"
            for sample in [i for i in cp.sections() if i.startswith("sample")]:

                # Fuse index if needed
                if self.idx2:
                   index_seq = (cp.get(sample, "index1_seq")+(cp.get(sample, "index2_seq"))).upper()
                else:
                    index_seq = (cp.get(sample, "index1_seq")).upper()

                # Create a autoreferenced Sample object
                Sample(
                    name=cp.get(sample, "name"),
                    index=index_seq,
                    write_pass=True ,
                    write_fail=self.write_fail,
                    minimal_qual=self.minimal_qual)
                    
            # Create a Sample object for undetermined reads
            Sample (
                name="Undetermined",
                index="N",
                write_pass=True if self.write_undetermined else False,
                write_fail=self.write_fail if self.write_undetermined else False,
                minimal_qual=0)

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
        """
        Main function of the script
        """
        # For double indexing
        if self.idx2:
            self.double_index_parser()
        # For simple indexing
        else:
            self.simple_index_parser()

        # Flush remaining content in sample buffers
        Sample.FLUSH_ALL()

        if self.write_report:
            print ("Generate_report")
            Sample.WRITE_REPORT()

        print ("Done")
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
                
            # Update report
            if self.write_report:
                Sample.WRITE_REPORT()

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
                
            # Update report
            if self.write_report:
                Sample.WRITE_REPORT()


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
        """
        Verify the readability of a file or list of file
        """
        for fp in file_list:
            if not os.access(fp, os.R_OK):
                raise ValueError ("{} is not a valid file".format(fp))
                

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sample(object):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    TOTAL = TOTAL_FAIL = TOTAL_PASS = 0
    BUFFER_SIZE = 20
    NAME_TO_SAMPLE = {}
    INDEX_TO_SAMPLE = {}
    DNA = ["A","T","C","G", "N"]

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def increment_TOTAL(self):
        self.TOTAL += 1

    @ classmethod
    def increment_TOTAL_FAIL(self):
        self.TOTAL_FAIL += 1

    @ classmethod
    def increment_TOTAL_PASS(self):
        self.TOTAL_PASS += 1

    @ classmethod
    def FINDER (self, read1, read2, index, molecular=""):
        """
        Try to find the sample corresponding to the index sequence, Undetermined by default
        """
        # Look up in the sample dict access by INDEX sequence
        sample = self.INDEX_TO_SAMPLE.get(index.seq, self.NAME_TO_SAMPLE["Undetermined"])
        
        # Call the Sample object by the __call__ method
        sample(read1, read2, index, molecular)

    @ classmethod
    def FLUSH_ALL (self):
        for sample in self.NAME_TO_SAMPLE.values():
            if sample.write_pass:
                sample._flush_to_file(file = sample.R1_pass_file, buffer = sample.R1_pass_buffer)
                sample._flush_to_file(file = sample.R2_pass_file, buffer = sample.R2_pass_buffer)
            if sample.write_fail:
                sample._flush_to_file(file = sample.R1_fail_file, buffer = sample.R1_fail_buffer)
                sample._flush_to_file(file = sample.R2_fail_file, buffer = sample.R2_fail_buffer)

    @ classmethod
    def WRITE_REPORT (self):
        with open ("Quade_report.txt", "wb") as report:
            report.write ("Total pair analysed\t{}\n".format(self.TOTAL))
            report.write ("Total pair pass quality\t{} ({}%)\n".format(self.TOTAL_PASS, self.TOTAL_PASS*100/self.TOTAL))
            report.write ("Total pair fail quality\t{} ({}%)\n".format(self.TOTAL_FAIL, self.TOTAL_FAIL*100/self.TOTAL))
            report.write ("Number of samples\t{}\n".format(len(self.NAME_TO_SAMPLE.values())-1))
            sample_list = [sample for sample in self.NAME_TO_SAMPLE.keys() if sample != "Undetermined"]
            sample_list.sort()
            for sample in sample_list:
                report.write ("\nSample Name\t{}\n".format(self.NAME_TO_SAMPLE[sample].name))
                report.write ("\tIndex sequence\t{}\n".format(self.NAME_TO_SAMPLE[sample].index))
                if self.NAME_TO_SAMPLE[sample].total == 0:
                    report.write ("\tNo read found\n") 
                else:
                    report.write ("\tTotal pair analysed\t{} ({}% all)\n".format(
                        self.NAME_TO_SAMPLE[sample].total,
                        self.NAME_TO_SAMPLE[sample].total * 100 / self.TOTAL))
                    report.write ("\tTotal pair pass quality\t{} ({}% sample)\n".format(
                        self.NAME_TO_SAMPLE[sample].total_pass,
                        self.NAME_TO_SAMPLE[sample].total_pass * 100 / self.NAME_TO_SAMPLE[sample].total))
                    report.write ("\tTotal pair fail quality\t{} ({}% sample)\n".format(
                        self.NAME_TO_SAMPLE[sample].total_fail,
                        self.NAME_TO_SAMPLE[sample].total_fail * 100 / self.NAME_TO_SAMPLE[sample].total))
            report.write ("\nUndetermined reads\n")
            report.write ("\tTotal pair analysed\t{} ({}% all)\n".format(
                self.NAME_TO_SAMPLE["Undetermined"].total,
                self.NAME_TO_SAMPLE["Undetermined"].total * 100 / self.TOTAL))

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, index, write_pass=True, write_fail=False, minimal_qual=0):

        # Create self variables
        self.name = name
        self.index = index 
        self.write_pass = write_pass
        self.write_fail = write_fail if minimal_qual > 0 else False # write fail if quality filtering only
        self.minimal_qual = minimal_qual

        # test uniqueness of name and index and DNA composition of index
        assert self.name not in self.NAME_TO_SAMPLE, "{} : Name is not unique".format(self.name)
        assert self.index not in self.INDEX_TO_SAMPLE, "{} : Index is not unique".format(self.name)
        assert self._is_dna(index), "{} : Non canonical DNA base in index".format(self.name)

        # Counters
        self.total = self.total_pass = self.total_fail = 0

        # Str Buffers and file to write seq
        if self.write_pass:
            self.current_pass = 0
            self.R1_pass_buffer = ""
            self.R2_pass_buffer = ""
            self.R1_pass_file = "./{}_R1_pass.fastq.gz".format(self.name)
            self.R2_pass_file = "./{}_R2_pass.fastq.gz".format(self.name)
            self._init_file (self.R1_pass_file)
            self._init_file (self.R2_pass_file)

        if self.write_fail:
            self.current_fail = 0
            self.R1_fail_buffer = ""
            self.R2_fail_buffer = ""
            self.R1_fail_file = "./{}_R1_fail.fastq.gz".format(self.name)
            self.R2_fail_file = "./{}_R2_fail.fastq.gz".format(self.name)
            self._init_file (self.R1_fail_file)
            self._init_file (self.R2_fail_file)

        # Add sample object to class dictionaries with name or index as a keys
        self.INDEX_TO_SAMPLE[self.index] = self
        self.NAME_TO_SAMPLE[self.name] = self

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

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self, read1, read2, index, molecular=""):

        # Increment class and object total read counter
        self.increment_TOTAL()
        self.total += 1

        # If the index quality is good
        if self._pass_quality(index):

            # Increment counters
            self.increment_TOTAL_PASS()
            self.total_pass += 1

            if self.write_pass:
                self.current_pass += 1
                # Update buffers
                self.R1_pass_buffer += self._fastq_str (read1, index, molecular)
                self.R2_pass_buffer += self._fastq_str (read2, index, molecular)

                if self.current_pass >= self.BUFFER_SIZE:
                    # Empty buffer in file and reset counter and buffers
                    self._flush_to_file (file = self.R1_pass_file, buffer = self.R1_pass_buffer)
                    self._flush_to_file (file = self.R2_pass_file, buffer = self.R2_pass_buffer)
                    self.current_pass = 0
                    self.R1_pass_buffer = ""
                    self.R2_pass_buffer = ""

        # Else if the quality is not high enough
        else:

            # Increment counters
            self.increment_TOTAL_FAIL()
            self.total_fail += 1

            if self.write_fail:
                self.current_fail += 1
                # Update buffers
                self.R1_fail_buffer += self._fastq_str (read1, index, molecular)
                self.R2_fail_buffer += self._fastq_str (read2, index, molecular)

                if self.current_fail >= self.BUFFER_SIZE:

                    # Empty buffer in file and reset counter and buffers
                    self._flush_to_file (file = self.R1_fail_file, buffer = self.R1_fail_buffer)
                    self._flush_to_file (file = self.R2_fail_file, buffer = self.R2_fail_buffer)
                    self.current_fail = 0
                    self.R1_fail_buffer = ""
                    self.R2_fail_buffer = ""


    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _is_dna (self, sequence):
        """
        Verify if the sequence is DNA
        """
        for base in sequence:
            if base not in self.DNA:
                return False
        return True

    def _pass_quality(self, seq):
        """
        Test minimal quality of a phred quality Seq
        """
        return min(seq.qual) >= self.minimal_qual

    def _fastq_str (self, read, index, molecular=""):
        """
        Generate a fastq str of read from read, index and molecular Seq Record
        """
        if molecular:
            return "@{}:{}:{}\n{}\n+\n{}\n".format(read.name.split()[0], index.seq, molecular.seq, read.seq, read.qualstr)
        else:
            return "@{}:{}\n{}\n+\n{}\n".format(read.name.split()[0], index.seq, read.seq, read.qualstr)

    def _init_file (self, file):
        """
        Initialize an empty file = erase existing file
        """
        with gzip.open (file, "wb") as fastq_file:
            pass

    def _flush_to_file (self, file, buffer):
        """
        Flush str buffer in file 
        """
        with gzip.open (file, "ab") as fastq_file:
            fastq_file.write(buffer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    quade = Quade()
    quade()
