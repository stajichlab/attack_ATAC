#!/usr/bin/env python3

import gzip, csv, os, sys, argparse


parser = argparse.ArgumentParser(description='gather mosdepth coverage.')
parser.add_argument('-i','--indir', default="coverage",
                    help='input directory to read mosdepth coverage files from')
parser.add_argument('-o','--outdir', default="summary_results",
                    help='outdir directory to write mosdepth coverage summary')
args = parser.parse_args()

os.listdir()
