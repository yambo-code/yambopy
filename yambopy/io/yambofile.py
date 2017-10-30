import yamboparser

class YamboFile(yamboparser.YamboFile):
    """
    Wrapper around the Yambofile class of yamboparser
    """
    def from_json(filename):
        """
        Read a local JSON file already parsed with yambopy
        """
        pass
                
    def from_file(filename):
        """
        Read a file and find what type it is
        """
        #detect the type
        #store the data
        pass

    def write_json(self,filename):
        """
        Write a json file with the data for this file
        """
        pass

class YamboGW(YamboFile):
    """
    Provide functions specific of GW calculations
    """
    def plot_gw(self):
        pass

class YamboEPS(YamboFile):
    """
    Provide functions specific of BSE calculations
    """
    def plot_eps(self):
        pass

