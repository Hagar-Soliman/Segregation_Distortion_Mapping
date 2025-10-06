#-------------------------------------------------------------------------------
# JKK 1st program applied to vcf on recombinants
# 
#
#-------------------------------------------------------------------------------

#from __future__ import division
import gzip
import sys
import os

library = "1A" #change the library ID!!!!!!!!!!!!

def main():
    #edit file paths!!!!!!!!!!!!!!
    input_path = f"/home/hks25/palmer_scratch/SNP_calling/{library}_VCF/{library}.vcf.gz"
    output_prefix = f"/home/hks25/palmer_scratch/GOOGA/{library}_output/"
    os.makedirs(output_prefix, exist_ok=True)
    
    src = gzip.open(input_path, "rt")  #read in compressed VCF file 
    out2 = open(output_prefix + "read_depth.all.txt", "w")
    out2a = open(output_prefix + "read_depth.included.txt", "w")
    out1 = open(output_prefix + "Slim.vcf", "w")
    out5 = open(output_prefix + "SNPsPerScaff.txt", "w")

    # user defined
    Min_MQ_score = 20
    Min_lines_called = 50
    minQ = 0.2
    maxQ = 0.8

    LineID = []
    sc_count = {}

    g_snps = [0, 0, 0]
    mbsnps = 0

    # consider every line in input file
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t')

        if len(cols) < 2:  # headline
            out1.write(line)

        elif cols[0] == "#CHROM":
            out1.write(line)
            for i in range(len(cols)):
                print(i, cols[i])
                if i > 8:
                    LineID.append(cols[i])

            Number_RILs = len(cols) - 9
            print("Samples ", Number_RILs)
        else:
            scx = cols[0]  # scaffold ID
            position = int(cols[1])
            ref_base = cols[3]
            alt_base = cols[4]
            if len(alt_base) > 1:  #Multi-nucleotide polymorphism (MNP)
                mbsnps += 1

            else:
                # extract MQ score
                infoSNP = cols[7].split(";")
                u78 = 2
                while u78 > 0:
                    mq = infoSNP[u78].split("=")
                    if mq[0] == "MQ":
                        mq_score = float(mq[1])
                        u78 = -1
                    else:
                        u78 += 1

                if mq_score >= Min_MQ_score:
                    g_snps[0] += 1
                    gc = [0, 0, 0.0]

                    # Get the index of the FORMAT field
                    format_fields = cols[8].split(":")
                    ad_index = format_fields.index("AD")

                    # Process sample data
                    for j in range(9, 9 + Number_RILs):
                        bases = cols[j].split(":")
                        if bases[0] != "./.":
                            ad_field = bases[ad_index]
                            alt_ref = ad_field.split(",")
                            cx1 = int(alt_ref[1])
                            cx2 = int(alt_ref[0])
                            if (cx1 + cx2) > 0:
                                gc[0] += 1
                                gc[1] += (cx1 + cx2)
                                gc[2] += float(cx2) / float(cx1 + cx2)

                    if gc[0] >= Min_lines_called:
                        g_snps[1] += 1
                        qR = gc[2] / float(gc[0])
                        out2.write(cols[0] + '\t' + cols[1] + '\t' + str(gc[0]) + '\t' + str(float(gc[1]) / float(gc[0])) + '\t' + str(qR) + '\n')
                        if qR <= maxQ and qR >= minQ:
                            g_snps[2] += 1
                            out1.write(line)
                            out2a.write(cols[0] + '\t' + cols[1] + '\t' + str(gc[0]) + '\t' + str(float(gc[1]) / float(gc[0])) + '\t' + str(qR) + '\n')
                            try:
                                sc_count[scx][0] += 1
                                sc_count[scx][1] = position
                            except KeyError:
                                sc_count[scx] = [1, position]

        if line_idx % 100000 == 0 and line_idx > 0:
            print(line_idx, cols[0], cols[1])

    print("Snps in Slim: passed mq, passed min n, passed p min ", g_snps)
    print("Multi-basers thrown ", mbsnps)
    for zz in sc_count:
        out5.write(zz + "\t" + str(sc_count[zz][0]) + "\t" + str(sc_count[zz][1]) + "\n")                        

if __name__ == "__main__":
    main()
