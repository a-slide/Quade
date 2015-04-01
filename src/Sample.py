# -*- coding: utf-8 -*-

"""
@package    Quade
@brief      Helper class for Quade to represent Samples
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Local imports
from FastqWriter import FastqWriter

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Sample(object):
    """
    This class store and verify the sample names and indexes at the instantiation of Quade. Sample
    object are auto referenced by the class collections NAME_TO_SAMPLE,INDEX_TO_SAMPLE and
    SAMPLE_LIST. The class method CLASS_INIT can be used to define which fatsq files will be
    written thereafter. In a second time the class method FINDER is able to attribute a pair of
    reads with an index to the corresponding sample where it will be written in a file if required
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    # Counters for the overall number of read at class level
    TOTAL = FAIL_QUAL = PASS_QUAL = UNDETERMINED = 0
    # Flag to be modified with CLASS_INIT to write fastq files
    WRITE_UNDETERMINED = WRITE_PASS = WRITE_FAIL = False
    # Class dictionaries to index Sample object either by index sequence or by name
    NAME_TO_SAMPLE = {}
    INDEX_TO_SAMPLE = {}
    SAMPLE_LIST = []
    # Definition of canonical DNA sequences
    DNA = ["A","T","C","G","N"]
    # Minimal numeric phred
    MIN_QUAL = 0
    # Init a class owned FastqWriter for the undetermined reads
    UNDETERMINED_WRITER = FastqWriter(name="Undetermined")

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def CLASS_INIT (self, write_undetermined=False, write_pass=False, write_fail=False, min_qual=0):
        """Generic class method to define all class fields at once"""
        self.WRITE_UNDETERMINED = write_undetermined
        self.WRITE_PASS = write_pass
        self.WRITE_FAIL = write_fail
        self.MIN_QUAL = min_qual

    @ classmethod
    def FINDER (self, read1, read2, index, molecular=""):
        """
        Try to find the sample corresponding to the index sequence, Undetermined by default
        """
        # Increment pair counter
        self.TOTAL += 1

        # Check if the index sequence correspond to a sample index
        if index.seq.upper() in self.INDEX_TO_SAMPLE:
            # Extract the sample object
            sample = self.INDEX_TO_SAMPLE[index.seq.upper()]

            # If the quality of the index is high enough + increment counters
            if min(index.qual) >= self.MIN_QUAL:
                self.PASS_QUAL +=1
                sample.pass_qual += 1
                # Write only if required
                if self.WRITE_PASS:
                    sample.pass_writer(read1, read2, index, molecular)

            # In case the quality filter is not passed + increment counters
            else:
                self.FAIL_QUAL += 1
                sample.fail_qual += 1
                # Write only if required
                if self.WRITE_FAIL:
                    sample.fail_writer(read1, read2, index, molecular)

        # If the index sequence do not correspond to a sample
        else:
            # Increment counters
            self.UNDETERMINED += 1
            # Write only if required
            if self.WRITE_UNDETERMINED:
                self.UNDETERMINED_WRITER(read1, read2, index, molecular)

    @ classmethod
    def FLUSH_ALL (self):
        """Flush all fastq buffers at once"""
        for sample in self.NAME_TO_SAMPLE.values():
            if sample.pass_writer.counter > 0:
                sample.pass_writer.flush_buffers()
            if sample.fail_writer.counter > 0:
                sample.fail_writer.flush_buffers()
        if self.UNDETERMINED_WRITER.counter > 0:
            self.UNDETERMINED_WRITER.flush_buffers()

    @ classmethod
    def REPORT (self):
        """Output a report under the form of a list"""
        report = []
        report.append(["Total pair", self.TOTAL])
        report.append(["Pair pass quality", self.PASS_QUAL])
        report.append(["Pair fail quality", self.FAIL_QUAL])
        report.append(["Pair Undetermined", self.UNDETERMINED])
        if self.TOTAL-self.UNDETERMINED > 0:
            report.append(["Percent Pair pass quality (wo Undetermined)", self.PASS_QUAL*100/(self.TOTAL-self.UNDETERMINED)])
            report.append(["Percent Pair fail quality (wo Undetermined)", self.FAIL_QUAL*100/(self.TOTAL-self.UNDETERMINED)])
            report.append(["Percent Pair Undetermined", self.UNDETERMINED*100/self.TOTAL])

        for sample in self.SAMPLE_LIST:
            report.append([" ", " "])
            report.append(["Sample Name", sample.name])
            report.append(["Total pair", sample.total])
            report.append(["Pair pass quality", sample.pass_qual])
            report.append(["Pair fail quality", sample.fail_qual])
            if sample.total > 0:
                report.append(["Percent of total pair", sample.total*100/(self.TOTAL-self.UNDETERMINED)])
                report.append(["Percent Pair pass quality", sample.pass_qual*100/sample.total])
                report.append(["Percent Pair fail quality", sample.fail_qual*100/sample.total])

        return report

    #~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, name, index):

        # Create self variables
        self.name = name
        self.index = index.upper()

        # test uniqueness of name and index and DNA composition of index
        assert self.name not in self.NAME_TO_SAMPLE, "{} : Name is not unique".format(self.name)
        assert self.index not in self.INDEX_TO_SAMPLE, "{} : Index is not unique".format(self.name)
        assert self._is_dna(index), "{} : Non canonical DNA base in index".format(self.name)

        # Counters
        self.pass_qual = self.fail_qual = 0

        # fastq_writer objects to manage the writing in fastq files
        self.pass_writer = FastqWriter(name = "{}_pass".format(self.name))
        self.fail_writer = FastqWriter(name = "{}_fail".format(self.name))

        # Add sample object to class dictionaries with name or index as a keys
        self.INDEX_TO_SAMPLE[self.index] = self
        self.NAME_TO_SAMPLE[self.name] = self
        self.SAMPLE_LIST.append(self)

    @property
    def total(self):
        return (self.pass_qual+self.fail_qual)

    # Fundamental class functions str and repr
    def __str__(self):
        msg = "SAMPLE CLASS\n"
        for key, value in self.__dict__.items():
            msg+="\t{}\t{}\n".format(key, value)
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _is_dna (self, sequence):
        for base in sequence:
            if base not in self.DNA:
                return False
        return True
