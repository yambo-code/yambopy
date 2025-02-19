import os
import re
from subprocess import Popen, PIPE
from yambopy.env import yambopyenv
from yambopy.tools.duck import isstring

class YamboIO(object):
    """
    This class provides input/output capabilities to yambopy.

    Currently: called by other classes via the msg function, 
               it opens/closes a LOG file and fills it with the messages 

    Example of use:

        .. code-block:: python

            yf = YamboIO(print_to_shell=False)
            yf.IO_start()
            yf.msg('message to print')
            yf.IO_close()
    """
    def __init__(self,out_name='yambopy.log',out_path='.',print_to_shell=True):
        self.fpath = '%s/%s'%(out_path,out_name)
        self.print_to_shell = print_to_shell
        self.s = ""

    def IO_start(self):
        self.f = open(self.fpath,'w')

    def IO_close(self):
        self.f.close()

    def msg(self,msg_string):
        self.f.write('%s\n'%msg_string)
        self.s += '%s\n'%msg_string
        if self.print_to_shell: print(msg_string)   

    def __str__(self):
        s = self.s
        return s
