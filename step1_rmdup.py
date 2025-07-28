#!/usr/bin/python
#
# Split raw fastqs by i7 index
# Identify/remove PCR duplicates using i5 molecular barcode
#

import argparse
import sys

# Parse command line options
parser = argparse.ArgumentParser(description='Use molecular barcodes on i5 index to remove PCR duplicates. Fastq files should be gzipped and of the form:\n  prefix.[R1,R2,i5,i7].suffix')
parser.add_argument('-p', '--prefix', required=True, help='Prefix for fastq, up until, e.g., OPF2_5_A1D12')
parser.add_argument('-s', '--suffix', required=True, help='Fastq suffix, e.g., fastq')
args = parser.parse_args()

# Set prefix and suffix from command-line arguments
prefix = args.prefix
suffix = args.suffix

# Initialize variables
total = 0  # Total reads
dupes = 0  # PCR duplicates
uniqseqs = {}  # Store potential PCR duplicates (i5seq -> set of readtest)

# Debugging information
print(f"Opening files with prefix: {prefix}")
print(f"Suffix used: {suffix}")

# Open files using 'with' to ensure they close properly
with open("/home/ps2267/ycga_work/hybrid_necrosis/raw_data/Sample_PS_P2_L1_2024/" + prefix + "_R1." + suffix, 'r') as fwdin, \
     open("/home/ps2267/ycga_work/hybrid_necrosis/raw_data/Sample_PS_P2_L1_2024/" + prefix + "_R2." + suffix, 'r') as revin, \
     open("/home/ps2267/ycga_work/hybrid_necrosis/raw_data/Sample_PS_P2_L1_2024/" + prefix + "_i5." + suffix, 'r') as i5in, \
     open("/home/ps2267/ycga_work/hybrid_necrosis/raw_data/Sample_PS_P2_L1_2024/" + prefix + "_i7." + suffix, 'r') as i7in, \
     open(prefix + '.rmdup_R1.' + suffix, 'w') as fwdout, \
     open(prefix + '.rmdup_R2.' + suffix, 'w') as revout, \
     open(prefix + '.rmdup_i5.' + suffix, 'w') as i5out, \
     open(prefix + '.rmdup_i7.' + suffix, 'w') as i7out, \
     open(prefix + ".rmdup.log", 'w') as logout:

    sys.stderr.write("\n")

    while True:
        if total % 10000 == 0:  # Output log every 10,000 reads
            sys.stderr.write(f"{prefix} -- clusters examined: {total}; PCR dupes ID'd: {dupes}\r")
            logout.write(f"{prefix} -- clusters examined: {total}; PCR dupes ID'd: {dupes}\n")

        # Read i7 index
        i7header = i7in.readline()
        if not i7header:  # End of file
            sys.stderr.write(f"{prefix} -- clusters examined: {total}; PCR dupes ID'd: {dupes}\r")
            logout.write(f"{prefix} -- clusters examined: {total}; PCR dupes ID'd: {dupes}\n")
            break

        i7seq   = i7in.readline().strip()
        i7plus  = i7in.readline()
        i7phred = i7in.readline()
        total   += 1

        # Read i5 index, forward, and reverse reads
        i5header = i5in.readline()
        i5seq    = i5in.readline().strip()
        i5plus   = i5in.readline()
        i5phred  = i5in.readline()

        r1header = fwdin.readline()
        r1seq    = fwdin.readline().strip()
        r1plus   = fwdin.readline()
        r1phred  = fwdin.readline()

        r2header = revin.readline()
        r2seq    = revin.readline().strip()
        r2plus   = revin.readline()
        r2phred  = revin.readline()

        # Check for duplicates using first 10 bp of forward and reverse reads
        r1test = r1seq[:10]
        r2test = r2seq[:10]
        readtest = r1test + r2test

        if i5seq not in uniqseqs:
            uniqseqs[i5seq] = set([readtest])
            fwdout.write(r1header)
            fwdout.write(r1seq + '\n')
            fwdout.write(r1plus)
            fwdout.write(r1phred)

            revout.write(r2header)
            revout.write(r2seq + '\n')
            revout.write(r2plus)
            revout.write(r2phred)

            i5out.write(i5header)
            i5out.write(i5seq + '\n')
            i5out.write(i5plus)
            i5out.write(i5phred)

            i7out.write(i7header)
            i7out.write(i7seq + '\n')
            i7out.write(i7plus)
            i7out.write(i7phred)

        else:
            if readtest not in uniqseqs[i5seq]:
                uniqseqs[i5seq].add(readtest)
                fwdout.write(r1header)
                fwdout.write(r1seq + '\n')
                fwdout.write(r1plus)
                fwdout.write(r1phred)

                revout.write(r2header)
                revout.write(r2seq + '\n')
                revout.write(r2plus)
                revout.write(r2phred)

                i5out.write(i5header)
                i5out.write(i5seq + '\n')
                i5out.write(i5plus)
                i5out.write(i5phred)

                i7out.write(i7header)
                i7out.write(i7seq + '\n')
                i7out.write(i7plus)
                i7out.write(i7phred)

            else:
                dupes += 1

    # Write final log
    logout.write(f"\n\nClusters examined for sample {prefix}\n")
    logout.write(f"Total clusters examined: {total}\n")
    logout.write(f"PCR duplicates ID'd:     {dupes}\n")

    sys.stderr.write(f"\n\nClusters examined for sample {prefix}\n")
    sys.stderr.write(f"Total clusters examined: {total}\n")
    sys.stderr.write(f"PCR duplicates ID'd:     {dupes}\n")



