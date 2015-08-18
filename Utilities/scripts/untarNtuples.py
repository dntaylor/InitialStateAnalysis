#!/usr/bin/env python

import os
import sys
import pwd
import glob
import argparse


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Untar ntuples from submit directory")

    parser.add_argument('jobName', type=str, nargs='+', help='Job name for submission')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    userDir = '/nfs_scratch/%s' % pwd.getpwuid(os.getuid())[0]
    baseDirs = [os.path.basename(fname)
                for string in args.jobName
                for fname in glob.glob("%s/%s" % (userDir, os.path.basename(string)))]

    submitDirs = ['%s/%s' % (userDir,jobName) for jobName in baseDirs]
    for submitDir in submitDirs:
        if not os.path.exists(submitDir):
            print "%s: Submit directory does not exist" %submitDir
            continue

        if os.path.exists('%s/results.tar.gz' % submitDir):
            print "Extracting %s"
            os.system('tar -zxf %s/results.tar.gz' % submitDir)
        else:
            sample_names = [os.path.basename(dir)
                            for dir in glob.glob("%s/*" % submitDir)]
            for sample in sample_names:
                if os.path.exists('%s/%s/results.tar.gz' % (submitDir,sample)):
                    print "Extracting %s" % sample
                    os.system('tar -zxf %s/%s/results.tar.gz' % (submitDir,sample))
                else:
                    print '%s: no results.tar.gz found' % sample

    return 0


if __name__ == "__main__":
    main()

