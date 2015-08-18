'''
Several scripts useful for other modules in ISA.

Author: Devin N. Taylor, UW-Madison
'''

import os
import sys
import errno

def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise
