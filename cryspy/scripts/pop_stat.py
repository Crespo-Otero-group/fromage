#!/usr/bin/env python
"""
Give basic information on a Gaussian population analysis

Usage:

pop_stat.py mol.log

"""
import argparse
import sys
import numpy as np

from fromage.io import read_file as rf


def main(in_log, kind):
    charges, energy = rf.read_g_char(in_log, pop=kind)
    maxi = max(charges)
    mini = min(charges)
    avg = sum([abs(i) for i in charges]) / len(charges)

    print("{:}{:20.6f}".format("Maximum charge", maxi))
    print("{:}{:20.6f}".format("Minimum charge", mini))
    print("{:}{:10.6f}".format("Averagre absolute charge", avg))
    print("{:}{:20.6f}".format("Std. deviation", np.std(charges)))
    print("{:}{:21.6f}".format("Energy (a.u.)", energy))


if __name__ == '__main__':
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument("in_log", help="Input .log file with RESP analysis",
                        default="gaussian.log")
    parser.add_argument("-k", "--kind", help="Kind of population, mulliken or esp",
                        default="esp", type=str)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    main(args.in_log, args.kind)
