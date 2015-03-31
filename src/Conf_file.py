#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    Conf_file
@brief      Contain the template of the empty configuration file
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

def write_example_conf():

    with open ("Conf_file.txt", 'wb') as fp:
        fp.write ("""
###################################################################################################
#                                   QUADE CONFIGURATION FILE                                      #
###################################################################################################
# Values can by customized with users values, but the file template must remain unchanged,
# otherwise the program will not be able to load default values.
#

###################################################################################################
[quality]
# The quality encoding of your sequence have to be Illumina 1.8+ Phred+33. The program do not
# manage the other encoding scales

# Minimal quality for one base of the index to consider a read pair valid. 0 if no filtering
# required. (INTEGER)
minimal_qual : 25

###################################################################################################
[fastq]

# Path to the fastq files containing non demultiplexed sequences. Since fastq splited in several
# chunks of data a list of fastq can be given for each categories, in the same order in all
# categories. seq_R1 and seq_R2 are for the fastq files containing the insert sequencing reads,
# index_R1 is for the first index read and index_R1 for the second index read (required if double
# indexing only). Usually for simple indexing : seq_R1 = R1, seq_R2 = R3 and index_R1 = R2. For
# double indexing: seq_R1 = R1, seq_R2 = R4, index_R1 = R2 and index_R2 = R3

seq_R1 :   ../dataset/C1_R1.fastq.gz  ../dataset/C2_R1.fastq.gz  ../dataset/C3_R1.fastq.gz
seq_R2 :   ../dataset/C1_R4.fastq.gz  ../dataset/C2_R4.fastq.gz  ../dataset/C3_R4.fastq.gz
index_R1 : ../dataset/C1_R2.fastq.gz  ../dataset/C2_R2.fastq.gz  ../dataset/C3_R2.fastq.gz
index_R2 : ../dataset/C1_R3.fastq.gz  ../dataset/C2_R3.fastq.gz  ../dataset/C3_R3.fastq.gz

###################################################################################################
[index]

# At the exception of index 1 which is mandatory, indicate by a boolean if the index is to be used.
# index 2 is used only if sample are double indexed, molecular1 if a molecular barcoding was
# performed and molecular2 in case of double molecular indexing (BOOLEAN)
index2 : True
molecular1 : True
molecular2 : True

# Indicate the start and end positions of each index and molecular barcode within in index read
# using 1 base coordinates. If an index is not in usage (from the previous section), the positions
# of the index are not required(INTEGERS)
index1_start : 1
index1_end : 4
index2_start : 1
index2_end : 4
molecular1_start : 4
molecular1_end : 6
molecular2_start : 4
molecular2_end : 6

###################################################################################################
[output]

# Write fastq files for each sample containing reads passing the quality filter **
write_pass : True

# Write fastq files for each sample containing reads failing the quality filter
write_fail : True

# Write fastq files containing reads whose sample is undetermined
write_undetermined : True

###################################################################################################
# SAMPLE DEFINITIONS

# It is possible to include as many independant sample as required by duplicating a entire sample
# section. All sample name and index have to be unique and or to be organized as follow :
#   [sampleX] = Sample identifier section, where X is the sample id number starting from 1 for the
#   first sample and incrementing by 1 for each additional sample
#   name = Unique identifier that will be used to prefix the read files (STRING)
#   index1_seq = Index 1 DNA sequences associated with the sample (STRING)
#   index2_seq = Similar to index 1, only in case of double indexing (STRING)

[sample1]
name : S1
index1_seq : ACAG
index2_seq : ACAG

[sample2]
name : S2
index1_seq : CTTG
index2_seq : CTTG""")
