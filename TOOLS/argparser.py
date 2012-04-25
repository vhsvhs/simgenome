import sys

class ArgParser:
    def __init__(self, cliArgs):
        self.args = cliArgs
        self.params = {}
    
    # use this method to grab REQUIRED command-line parameters
    # flag = the command line flag.  For example: "--model".
    def getArg(self, flag):
        if self.args.__contains__(flag):
            i = self.args.index(flag)
            return self.args[i+1]
        else:
            message = "You did not specify the required command-line argument: " + flag 
            print "(ERROR) " + message
            raise AssertionError
            exit(1)
            
    # this is somewhat a hack:
    def setArg(self, flag, content):
        if self.args.__contains__(flag):
            i = self.args.index(flag)
            self.args[i] = flag
            self.args[i+1] = content
        else:
            self.args.append(flag)
            self.args.append(content)
            
    # use this method to grab OPTIONAL command-line parameters.
    def getOptionalArg(self, flag):
        if self.args.__contains__(flag):
            i = self.args.index(flag)
            return self.args[i+1]
        else:
            return False
    
    def doesContainArg(self, flag):
        return self.args.__contains__(flag)
    
    # use this method to grab OPTIONAL command-line toggles (boolean on/off switches)
    def getOptionalToggle(self, flag):
        if self.args.__contains__(flag):
            return True
        else:
            return False
    