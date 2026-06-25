import os

plants = 195
library = "6C"

base_dir = f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/"
output_dir = f"/home/hks25/palmer_scratch/googaV3/refAlt/{library}_refAlt/"
os.makedirs(output_dir, exist_ok=True)
src = open(f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/Slim.vcf", "r")
in1 = open(f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/SNPs.limited.txt", "r")
outT = open(f"/home/hks25/palmer_scratch/googaV3/VCF/{library}_VCF/{library}_output/snps_perscaff.txt", "w")

Relevant_snps = {}
scafftots = {}

for line_idx, line in enumerate(in1):
    cols = line.replace('\n', '').split('\t')
    Relevant_snps[cols[0] + "_" + cols[1]] = 1
    sc = cols[0]
    try:
        scafftots[sc] += 1
    except KeyError:
        scafftots[sc] = 1

for kk in scafftots:
    outT.write(kk + '\t' + str(scafftots[kk]) + '\n')

LineID = []

from collections import defaultdict

file_cache = defaultdict(lambda: None)

for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

    if len(cols) < 2:
        pass

    elif cols[0] == "#CHROM":
        for i in range(len(cols)):
            if i > 8:
                LineID.append(cols[i])
        plants = len(cols) - 9
        print("samps ", plants)

    else:
        scaff = cols[0]
        position = int(cols[1])
        try:
            Relevant_snps[scaff + '_' + cols[1]]
            ref_base = cols[3]
            alt_base = cols[4]
            if line_idx % 10000 == 0:
                print(scaff, position)

            for j in range(9, 9 + plants):
                all_scaffold_filename = output_dir + "all." + scaff + ".txt"
                line_id_filename = output_dir + os.path.basename(LineID[j - 9]) + ".txt"

                if file_cache[all_scaffold_filename] is None:
                    file_cache[all_scaffold_filename] = open(all_scaffold_filename, 'a')

                if file_cache[line_id_filename] is None:
                    file_cache[line_id_filename] = open(line_id_filename, 'a')

                current_file0 = file_cache[all_scaffold_filename]
                current_file1 = file_cache[line_id_filename]

                geno = cols[j].split(":")
                if geno[0] != "./.":
                    ref, alt = geno[3].split(",")
                else:
                    ref = "0"
                    alt = "0"

                current_file1.write(scaff + '\t' + cols[1] + '\t' + ref + '\t' + alt + '\n')
                current_file0.write(scaff + '\t' + cols[1] + '\t' + ref + '\t' + alt + '\n')

        except KeyError:
            pass

src.close()
in1.close()
outT.close()

for file in file_cache.values():
    if file is not None:
        file.close()
