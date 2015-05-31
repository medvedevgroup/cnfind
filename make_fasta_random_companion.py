import sys

if len(sys.argv) <= 1:
    print "will generate random fasta files corresponding to each chromosome, to make the --nogc option work" 
    print "usage: [path to companion package]"
    exit(1)


import glob
cpath = sys.argv[1]

import os
allchr = cpath + "/allchr.txt"
if not os.path.isfile(allchr):
    print cpath + " does not look like a companion package"
    exit(1)

#get chromosomes list
chrs = []
for line in open(allchr):
    chr = line.strip()
    chrs += [chr]

random_folder =  cpath + "/fasta_files_folder/random/" 
if not os.path.exists(random_folder):
        os.makedirs(random_folder)

cnfind_dir = os.path.dirname(os.path.realpath(__file__))
import subprocess
for chr in chrs:
    cmd = cnfind_dir + "/make_fasta_random < " + cpath + "/fasta_files_folder/" + chr + ".fa > " + random_folder + chr + ".fa"
    print "calling: " + cmd
    subprocess.call([cmd], shell=True)

