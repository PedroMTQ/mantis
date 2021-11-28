RequirementsNotMet='Installation check not passed! Make sure you\'ve setup the databases and your system meets all the requirements!'
NoValidFiles='No valid files to annotate'
TargetFileNotFound='Target file does not exist!\n'
InvalidTargetFile='You did not insert a valid target file!\n'
InvalidFastaFormat='Fasta format is not valid!\n'
MissingMetadataFile='Metadata file not found. Did you run setup_databases before?\n'
InstallationCheckNotPassed='Installation check not passed! Make sure you\'ve setup the databases and your system meets all the requirements!'
CythonNotCompiled= 'Cython has not been correctly compiled! Please go to mantis/source/ and run python utils.py'
BadNumberWorkers='You should not be seeing this, please contact the developer. Invalid number of workers in '
ConnectionError='Could not connect to url:\n'
InvalidTranslation='Invalid residues for translation. Please make sure you provided a CDS of DNA or RNA in the target file:\n'
InvalidGFFVersion='No valid GFF version found!'
InvalidNOGType='Your config file does not contain a valid database type for NOG (i.e. <dmnd> or <hmm>)!'

class RequirementsNotMet(Exception):
    def __init__(self, *args):
        if args:
            self.message = args
        else:
            self.message = None

    def __str__(self):
        if self.message:
            res = 'Requirements not met:\n'
            res += '\n'.join(i for i in self.message)
            return res
        else:
            return 'Requirements not met'


