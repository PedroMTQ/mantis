class RequirementsNotMet(Exception):
    def __init__(self,*args):
        if args: self.message=args
        else: self.message=None

    def __str__(self):
        if self.message:
            res='Requirements not met:\n'
            res+='\n'.join(i for i in self.message)
            return res
        else:
            return 'Requirements not met'

class MissingCondaEnvironment(Exception):
    def __str__(self):
        return 'Conda environment is missing from configuration file!'

class InvalidFunction(Exception):
    def __str__(self):
        return 'This function is not available in this module!'

class InvalidTargetFile(Exception):
    def __str__(self):
        return 'You did not insert a valid target file! Make sure it follows the same format as the file provided in tests/test_file.tsv'

class InstallationCheckNotPassed(Exception):
    def __str__(self):
        return 'Installation check not passed! Make sure you\'ve setup the databases and your system meets all the requirements!'

class CythonNotCompiled(Exception):
    def __str__(self):
        return 'Cython has not been correctly compiled! Please go to mantis/source/ and run python utils.py'

