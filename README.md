# Quade

**Fastq file demultiplexer, handling double indexing, molecular indexing and filtering based on index quality**

[see GitHub Page](http://a-slide.github.io/Quade) 

**Creation : 2015/01/07**

**Last update : 2015/02/06** 

## Motivation

Quade is a **python2.7** object oriented script developed to demultiplex sample from mixed fastq files.

Specific features:

* The program parse chunks of non demultiplexed raw paired end fastq files (converted in fastq by CASAVA but not demultiplexed, see [convert-bcl-to-fastq](https://gist.github.com/brantfaircloth/3125885))
* Reads are attributed to a specific sample defined by index sequences. Simple and double indexing are supported.
* Reads can be filtered to exclude index read with bad quality base (to avoid sample cross contamination due to errors in index read)
* The program handle stochastic molecular index. Sample index and molecular index are appended at the end of sequence name for further use.

![Double index](https://raw.githubusercontent.com/a-slide/Quade/master/doc/img/modified_illumina_adapters.png)

## Principle

1. A configuration file containing all program parameters (including sample/index association) is parsed and thoroughly verified for validity
2. Non demultiplexed fastq files are read chunks by chunks, with reads from sample read1, sample read2, index read1 [and index read2] in different files.
3. sample barcode and molecular index sequence are extracted from index read. In the case of double indexing sample barcode is the result of the fusion between sample barcode from index1 and index2. The same applies for molecular index 
4. Depending of the index (or fused index) sequence, sample reads1 and read2 are attributed to a defined sample, or to the undetermined category if appropriate.
5. The phred quality of the index read is checked so as every positions is above a defined minimal value. If so reads are written in a sample **pass** file, else reads are written in a sample **fail** file.
6. reads names are written as follow : *@Name_in_original_fastq:SAMPLE_BARCODE[:MOLECULAR_BARCODE]*
7. a distribution report is generated

![Quade](https://raw.githubusercontent.com/a-slide/Quade/master/doc/img/quade-io.png)

## Dependencies

The program was developed under Linux Mint 17 and require a python 2.7 environment.
The following dependencies are required for proper program execution:

* python package [numpy](http://www.numpy.org/) 1.7.1+
* python package [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/index.html#) 0.6.1p1+

If you have pip already installed, enter the following line to install packages: ```sudo pip install numpy HTSeq```

## Get and install Quade

* Clone the repository or download the archive ```git clone https://github.com/a-slide/Quade.git```

* Enter the src folder of the program folder and make the main script executable ```sudo chmod u+x Quade.py```

* Alternatively, you can compile Quade as a C executable with cython by using the makefile ```make```. Be sure that you have cython installed and added to your PATH

* Finally, Add Quade.py (or Quade if you choose the cython alternative)  in your PATH

## Usage

In the folder where fastq files will be output
    
    Usage: Quade -c Conf.txt
    
    Options:
      --version     show program's version number and exit
      -h, --help    show this help message and exit
      -c CONF_FILE  Path to the configuration file
  
An example configuration file is provided with the program, with detailed information for each sections

## Additional information

See Blog post [Starting a new NGS project with python : Example with a fastq demultiplexer](http://a-slide.github.io/blog/fastq_demultiplexer)

## Authors and Contact

Adrien Leger - 2014

* <adrien.leger@gmail.com> - <adrien.leger@inserm.fr> - <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089](http://www.atlantic-gene-therapies.fr/)