from yambopy.tools.duck import isstring
import yamboparser

class YamboFile(yamboparser.YamboFile):
    """
    Wrapper around the Yambofile class of yamboparser
    """
    def from_dict(filename):
        """
        intialize from a dicitonary
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

    @staticmethod
    def has_tag(filename,tags):
        """check if the filename has a tag in its name"""
        if isstring(tags):
            tags = (tags,)
        return any([tag in filename for tag in tags])

    @staticmethod
    def is_output(filename):
        """check if the file is output file"""
        return filename.startswith('o.') or filename.startswith('o-')

    @staticmethod
    def is_log(filename):
        """check if the file is log file"""
        return filename.startswith('l.') or filename.startswith('l-')

    @staticmethod
    def is_report(filename):
        """check if the file is report file"""
        return filename.startswith('r-')

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

