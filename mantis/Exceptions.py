NoValidFiles='No valid files to annotate'
TargetFileNotFound='Target file does not exist!\n'
InvalidTargetFile='You did not insert a valid target file!\n'
InvalidFastaFormat='Fasta format is not valid!\n'
InstallationCheckNotPassed='Installation check not passed! Make sure you\'ve setup the databases and your system meets all the requirements!'
SQLCheckNotPassed='SQL check not passed! Make sure you\'ve setup the databases and your system meets all the requirements!'
CythonNotCompiled= 'Cython has not been correctly compiled! Please run:\n'
BadNumberWorkers='You should not be seeing this, please contact the developer. Invalid number of workers in '
ConnectionError='Could not connect to url:\n'
InvalidTranslation='Invalid residues for translation. Please make sure you provided a CDS of DNA or RNA in the target file:\n'
InvalidGFFVersion='No valid GFF version found!'
InvalidNOGType='Your config file does not contain a valid database type for NOG (i.e. <dmnd> or <hmm>)!'




class SQLCheckNotPassed(Exception):
    def __init__(self):
        self.message = 'SQL check not passed!'
    def __str__(self):
        return self.message



class InstallationCheckNotPassed(Exception):
    def __init__(self):
        self.message = 'Installation check not passed! Make sure you\'ve setup the databases'
    def __str__(self):
        return self.message


