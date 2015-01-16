# Imports
import ConfigParser
from sys import argv

# Definition of the main class
class Main(object):

    # Class init method
    def __init__ (self): # conf_file arg is only for interactive interpreter. Will be remove in final program

        # Define a configuration file parser object
        self.conf = ConfigParser.RawConfigParser(allow_no_value=True)

        # Indicate the path to the configuration file
        self.conf.read(argv[1])
        print(self.conf.sections())

        # In a try block due to possible errors while reading the file
        try:
            #### PARSE CONF FILE ####

            # Get Values from the required section in the appropriate type
            self.outdir = self.conf.get("general", "outdir") # for string values
            self.mean_qual = self.conf.getint("quality", "mean_qual") # for integer values
            self.minimal_qual = self.conf.getint("quality", "minimal_qual")
            self.quality_scale = self.conf.get("quality", "quality_scale")

            # for each index verify its existence and store start and end positions if needed
            self.index1 = self.conf.getboolean("index", "index1")  # for boolean values
            if self.index1:
                self.index1_start = self.conf.getint("index", "index1_start")
                self.index1_end = self.conf.getint("index", "index1_end")
            self.index2 = self.conf.getboolean("index", "index2")
            if self.index2:
                self.index2_start = self.conf.getint("index", "index2_start")
                self.index2_end = self.conf.getint("index", "index2_end")
            self.molecular1 = self.conf.getboolean("index", "molecular1")
            if self.molecular1:
                self.molecular1_start = self.conf.getint("index", "molecular1_start")
                self.molecular1_end = self.conf.getint("index", "molecular1_end")
            self.molecular2 = self.conf.getboolean("index", "molecular1")
            if self.molecular2:
                self.molecular2_start = self.conf.getint("index", "molecular2_start")
                self.molecular2_end = self.conf.getint("index", "molecular2_end")

            # Sample section are a special case since the number of sample section is variable
            # Only section with a name starting by "sample" will be parsed
            self.sample_list = []
            for sample in [i for i in parser.sections() if i.startswith("sample")]:
                self.sample_list.append({
                    'name'   : self.conf.get(sample, "name"),
                    'index1_seq' :  self.conf.get(sample, "index1_seq"),
                    'index2_seq' :  self.conf.get(sample, "index2_seq")})

            #### TEST VALUES ####

            # Now all variable can be tested and an assertion error will be raised in case of invalid value
            self.outdir = self.outdir if self.outdir else "./"
            self.quality_scale = self.conf.get("quality", "quality_scale")

            assert self.mean_qual >= 0 and self.mean_qual <= 40, "Autorized values for mean_qual : 0 to 40"
            assert self.mean_qual >= 0 and self.mean_qual <= 40, "Autorized values for minimal_qual : 0 to 40"
            if self.index1:
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


        # Handle the possible exception by printing simple message for users.
        except ConfigParser.NoOptionError as E:
            print (E)
            print ("An option is missing in the configuration file")
            print ("Please report to the descriptions in the configuration file\n")
            exit()
        except ConfigParser.NoSectionError as E:
            print (E)
            print ("An section is missing in the configuration file")
            print ("Please report to the descriptions in the configuration file\n")
            exit()
        except (ValueError, AssertionError) as E:
            print (E)
            print ("One of the value in the configuration file is not in the correct format")
            print ("Please report to the descriptions in the configuration file\n")
            exit()

    def __call__()
        print (self.__dict__)


#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#
if __name__ == '__main__':

    main = Main()
    main()
