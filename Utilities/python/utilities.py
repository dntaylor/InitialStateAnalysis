'''
Several scripts useful for other modules in ISA.

Author: Devin N. Taylor, UW-Madison
'''

import os
import sys
import errno
import hashlib
import re
import logging
import argparse

logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise

def get_parser(desc):
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('analysis', type=str, choices=['Z','WZ','WZ_W','WZ_Dijet','Hpp2l','Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('channel', type=str, choices=['Z','WZ','FakeRate','TT','Hpp2l','Hpp3l','Hpp4l','LowMass','Charge'], help='Channel to run for given analysis')
    parser.add_argument('period', type=int, choices=[8,13], help='Energy (TeV)')
    parser.add_argument('-l','--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')

    return parser


def hashfile(path, blocksize = 65536):
    '''Hash a file'''
    afile = open(path, 'rb')
    hasher = hashlib.md5()
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    afile.close()
    return hasher.hexdigest()

def hashstring(string):
    return hashlib.md5(string).hexdigest()

def hashlist(alist):
    return hashlib.md5(str(alist)).hexdigest()

def recurseList(index,items):
    result = []
    nextindex = index
    for i,item in enumerate(items):
        if i<nextindex: continue # not considered in this sub list
        if item=='(':            # new sub list
            subresult, nextindex = recurseList(i+1,items)
            result += [subresult]
        elif item==')':          # end of sub list
            return (result, i+1)
        else:                    # just a new item for the list
            result += [item]
    return result

def sortList(items):
    logicalPositions = [i for i,x in enumerate(items) if x=='&&' or x=='||']
    sortedElements = sorted([sortList(x) if isinstance(x,(list)) else x for i,x in enumerate(items) if i not in logicalPositions])
    for i,item in enumerate(items):
        if i not in logicalPositions:
            items[i] = sortedElements.pop(0)
    return items

def split(data):
    '''Take a string, split it up into distinct logical segments'''
    if data.count('(') != data.count(')'):
        logger.error('Unmatched parentheses in %s' % data)
        return []
    # split at ( and )
    items = re.split('(\(|\)|&&|\|\|)',data)
    # remove empty elements
    items = [x for x in items if x]
    # recombine functional forms (e.g. a(b)<c)
    prevItem = ''
    newItems = []
    functionDepth = 0
    for i,currItem in enumerate(items):
       if currItem in ['&&', '||']:         # logical break
           newItems += [prevItem]
           prevItem = currItem
       elif currItem in ['(']:
           if i==0:                         # start of cut
               prevItem = currItem
           elif prevItem in ['&&', '||', '('] and not functionDepth:   # this is a real start
               newItems += [prevItem]
               prevItem = currItem
           else:                            # this is the opening of a function
               prevItem += currItem
               functionDepth += 1
       elif currItem in [')']:
           if functionDepth:                # end of function
               prevItem += currItem
               functionDepth -= 1
           else:                            # real end
               newItems += [prevItem]
               prevItem = currItem
       elif functionDepth:                  # add to prev
           prevItem += currItem
       else:                                # something new
           if prevItem in ['&&', '||', '(', ')']:
               newItems += [prevItem]
               prevItem = currItem
           else:                            # continue item
               prevItem += currItem
    newItems += [prevItem]
    items = newItems
    # nest parantheses
    nestedItems = recurseList(0,items)
    # sort
    sortedItems = sortList(nestedItems)
    return sortedItems

def combine(data):
    '''Take the split list and combine it'''
    return reduce(lambda x,y: x+'('+combine(y)+')' if isinstance(y,(list)) else x+y, data, '')

def order(data):
    '''Order my cut to always be the same'''
    data = split(data)
    return combine(data)

def hashcut(cut):
    '''Hash a cut, doing basic error checking'''
    cut = cut.replace(' ','')
    # put all to bitwise
    cut = re.sub('&&+','&',cut)
    cut = re.sub('\|\|+','|',cut)
    # put all to logical
    cut = re.sub('&','&&',cut)
    cut = re.sub('\|','||',cut)
    # order the cuts
    cut = order(cut)
    # hash the cut
    return (cut, hashlib.md5(cut).hexdigest())

def hashscalefactor(scalefactor):
    scalefactor = '*'.join(sorted(scalefactor.split('*')))
    return (scalefactor, hashlib.md5(scalefactor).hexdigest())
