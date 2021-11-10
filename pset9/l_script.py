#! /usr/bin/env python3

# w09-naive.py
#   The analysis done by Lestrade et al.
#   Also serves as an example of parsing the input data.
# 
# Usage:
#   ./w09-naive.py <infile>
#
# Estimate the number of reads assigned to each isoform by a naive
# method.  Assign each read to its cognate segment, trusting the
# mapping and ignoring accuracy; then assign to isoform
# uniformly. (i.e. if 3 isoforms share the same segment, assign 1/3
# count to each isoform.)
#

import numpy as np
import string
import sys
import re

# Parse the input file.
#
with open(sys.argv[1]) as f:
    #   The first line is "The <n> transcripts of the sand mouse Arc locus"
    line  = f.readline()
    match = re.search(r'^The (\d+) transcripts', line)
    T     = int(match.group(1))

    # The next T lines are 
    #   <Arcn>  <true_tau> <L> <structure>
    # tau's may be present, or obscured ("xxxxx")
    tau       = np.zeros(T)
    L         = np.zeros(T).astype(int)
    tau_known = True   # until we see otherwise
    for i in range(T):
        fields    = f.readline().split()
        if fields[1] == "xxxxx":
            tau_known = False
        else:
            tau[i] = float(fields[1])
        L[i]      = int(fields[2])

    # after a blank line,
    # 'The <n> read sequences':
    line  = f.readline()
    line  = f.readline()
    match = re.search(r'The (\d+) read sequences', line)
    N     = int(match.group(1))

    # the next T lines are 
    #  <read a-j> <count>
    r = np.zeros(T).astype(int)
    for k in range(T):
        fields = f.readline().split()
        r[k]   = fields[1]


S = T    # S = R = T : there are T transcripts (Arc1..Arc10), S segments (A..J), R reads (a..j)
R = T
Slabel   = list(string.ascii_uppercase)[:S]               # ['A'..'J']        : the upper case labels for Arc locus segments 
Tlabel   = [ "Arc{}".format(d) for d in range(1,T+1) ]    # ['Arc1'..'Arc10'] : the labels for Arc transcript isoforms
Rlabel   = list(string.ascii_lowercase)[:T]               # ['a'..'j']        : lower case labels for reads


# Count how often each segment A..J is used in the isoforms i
# We'll use that to split observed read counts across the isoforms
# that they might have come from.
#
segusage = np.zeros(S).astype(int)
for i in range(T):
    for j in range(i,i+L[i]): 
        segusage[j%S] += 1


# Naive analysis:
#
c  = np.zeros(T)
for i in range(T):
    for k in range(i,i+L[i]):
        c[i] += (1.0 / float(segusage[k%S])) * float(r[k%S])  # For each read k, assume read k-> segment j,
                                                              # and assign 1/usage count to each transcript
                                                              # that contains segment j.
Z       = np.sum(c)
est_nu  = np.divide(c, Z)       # nucleotide abundance

print(c)

est_tau = np.divide(est_nu, L)  # convert to TPM, transcript abundance
est_tau = np.divide(est_tau, np.sum(est_tau))

# Print a table of the resulting estimates for tau
for i in range(T):
    print ("{0:10s}  {1:5.3f}".format(Tlabel[i], est_tau[i]))

