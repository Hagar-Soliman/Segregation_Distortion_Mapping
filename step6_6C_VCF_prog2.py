#VERY IMPORTANTLY- line 62 denotes the line number in your vcf file that are your parents (use the output from prog1 to change the script)

# you will have to run this script two times since this is a diagnostic script. first run produced info about the depth to decide on the cut offs, change the cut off to produce the SNP.limited.txt file when you run the script again
#note: it is okay not to be super strict here since we will do depth filtering in step9B at a window level

RADTAGsize=50
plants=195

library="6C"

src = open(f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/Slim.vcf", "r", newline=None)
out2 = open(f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/SNPs.all.txt", "w")
out2a = open(f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/SNPs.limited.txt", "w")
outD = open(f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/DepthInfo.txt", "w")
outAFD = open(f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/spectrum.txt", "w")

final_cat = 50
Dstat = [[0, 0.0] for j in range(final_cat + 1)]
Estat = [[] for j in range(final_cat + 1)]

for j in range(final_cat + 1):
    for k in range(j + 1):
        Estat[j].append(0)

snpsOut = 0
LineID = []

maxMedianDepth = 15.0  # set first to high value, inspect spectrum, reset to create real SNP.limited

includeX = [0, 0]

last_scaff = ''
bestsnp = ''
lastpos = 0
cp = 0

for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

    if len(cols) < 2:
        pass

    elif len(cols[0].split("_")) == 2:
        scID = int(cols[0].split("_")[1])
        if scID <= 17:
            scaff = cols[0]
            position = int(cols[1])

            ref_base = cols[3]
            alt_base = cols[4]
            if line_idx % 10000 == 0:
                print(scaff, position)

            if len(alt_base) > 1:
                print("why are there multiple bases here?")

            cx = 0
            cy = 0.0
            cz = 0
            datums = {}
            dlist = []
            for j in range(9, 9 + plants):
                geno = cols[j].split(":")
                if geno[0] != "./." and j != 9 and j != 10:  # j=9 = IM62, j=10 = CCC9 (parents excluded)
                    cx += 1

                    ref, alt = geno[3].split(",")

                    if (int(ref) + int(alt)) > 0:
                        datums[j] = [int(ref), int(alt)]
                        cy += float(datums[j][0]) / float(datums[j][0] + datums[j][1])
                        cz += (datums[j][0] + datums[j][1])
                        dlist.append(datums[j][0] + datums[j][1])

            dlist.sort()
            scored_snps = len(dlist)
            if scored_snps > 0:
                outD.write(scaff + '\t' + cols[1] + '\t' + str(cz) + '\t'
                    + str(cx) + '\t' + str(cz / float(cx)) + '\t'
                    + str(cy / float(cx)) + '\t' + str(dlist[int(scored_snps / 2)]) + '\n')

            if dlist[int(scored_snps / 2)] <= maxMedianDepth:
                snpsOut += 1
                includeX[0] += 1
                out2.write(scaff + '\t' + cols[1] + '\t')

                if scaff != last_scaff or (position - lastpos) > RADTAGsize:
                    out2a.write(bestsnp)
                    cp = scored_snps
                    bestsnp = scaff + '\t' + cols[1] + '\n'
                elif scored_snps > cp:
                    cp = scored_snps
                    bestsnp = scaff + '\t' + cols[1] + '\n'
                last_scaff = scaff
                lastpos = position

                linesp = datums.keys()
                cc = [0, 0, 0]
                for z in linesp:

                    if datums[z][0] + datums[z][1] >= 5:
                        if datums[z][0] >= 10 * datums[z][1]:
                            cc[0] += 1
                        elif datums[z][1] >= 10 * datums[z][0]:
                            cc[2] += 1
                        else:
                            cc[1] += 1
                    n = datums[z][0] + datums[z][1]
                    q = float(datums[z][0]) / float(n)

                    if n < final_cat:
                        Dstat[n][0] += 1
                        Dstat[n][1] += (q - 0.5) ** 2.0
                        Estat[n][datums[z][0]] += 1
                    else:
                        Dstat[final_cat][0] += 1
                        Dstat[final_cat][1] += (q - 0.5) ** 2.0
                        rval = int(float(final_cat) * q)
                        Estat[final_cat][rval] += 1

                out2.write(str(cc[0]) + '\t' + str(cc[1]) + '\t' + str(cc[2]) + '\n')

            else:
                includeX[1] += 1

out2a.write(bestsnp)  # last snp

for j in range(2, final_cat + 1):
    for k in range(j + 1):
        outAFD.write(str(j) + '\t' + str(k) + '\t' + str(Estat[j][k]) + '\n')

print("Snps output", snpsOut, "in/ex clude", includeX)

out2a.close()
out2.close()
