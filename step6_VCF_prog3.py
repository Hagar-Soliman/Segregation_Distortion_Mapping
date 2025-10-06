plants = 387
library = "1A"

base_dir = f"/home/hks25/palmer_scratch/GOOGA/{library}_output/"
src = open(f"/home/hks25/palmer_scratch/GOOGA/{library}_output/Slim.vcf", "r")
in1 = open(f"/home/hks25/palmer_scratch/GOOGA/{library}_output/SNPs.limited.txt", "r")
outT = open(f"/home/hks25/palmer_scratch/GOOGA/{library}_output/snps_perscaff.txt", "w")

#Relevant_snps is used to store relevant SNPs, and scafftots is used to store the total number of SNPs in each genomic scaffold
Relevant_snps = {}
scafftots = {}


for line_idx, line in enumerate(in1):
    #Remove newline characters and split each line using tab delimiters
    cols = line.replace('\n', '').split('\t')
    #Combine scaffold and position to generate a unique key, and set its value to 1
    Relevant_snps[cols[0] + "_" + cols[1]] = 1
    sc = cols[0]
    try:
        #If the scaffold already exists in the scafftots dictionary, increment its count by 1
        scafftots[sc] += 1
    except KeyError:
        #If the scaffold does not exist in the scafftots dictionary, initialize its count to 1
        scafftots[sc] = 1

#Write the data from the scafftots dictionary to the outT file
for kk in scafftots:
    outT.write(kk + '\t' + str(scafftots[kk]) + '\n')

#Initialize a list to store sample IDs
LineID = []

#Create a default dictionary to cache file objects
from collections import defaultdict

file_cache = defaultdict(lambda: None)

#Read the src file line by line
for line_idx, line in enumerate(src):
    #Remove newline characters and split each line using tab delimiters
    cols = line.replace('\n', '').split('\t')

    #Skip the line if the number of columns is less than 2 (this is usually the file header)
    if len(cols) < 2:
        pass

    #If the first column is "#CHROM", this line is the header; process sample IDs
    elif cols[0] == "#CHROM":
        for i in range(len(cols)):
            if i > 8:
                LineID.append(cols[i])
        
        #The number of samples equals the number of columns minus 9
        plants = len(cols) - 9
        print("samps ", plants)

    #Process other lines (i.e., actual data rows)
    else:
        scaff = cols[0]
        position = int(cols[1])
        try:
            #Check whether the current line corresponds to a relevant SNP
            Relevant_snps[scaff + '_' + cols[1]]
            ref_base = cols[3]
            alt_base = cols[4]
            if line_idx % 10000 == 0:
                print(scaff, position)

            for j in range(9, 9 + plants):
                #Optimization: use a caching mechanism to avoid frequent opening and closing of files
                all_scaffold_filename = "all." + scaff + ".txt"
                line_id_filename = LineID[j - 9] + ".txt"
                
                if file_cache[all_scaffold_filename] is None:
                    file_cache[all_scaffold_filename] = open(all_scaffold_filename, 'a')
                
                if file_cache[line_id_filename] is None:
                    file_cache[line_id_filename] = open(line_id_filename, 'a')
                
                current_file0 = file_cache[all_scaffold_filename]
                current_file1 = file_cache[line_id_filename]
                
                geno = cols[j].split(":")
                #Check whether genotype information is missing
                if geno[0] != "./.":
                    ref, alt = geno[3].split(",")
                else:
                    ref = "0"
                    alt = "0"

                #Write the results to a file
                current_file1.write(scaff + '\t' + cols[1] + '\t' + ref + '\t' + alt + '\n')
                current_file0.write(scaff + '\t' + cols[1] + '\t' + ref + '\t' + alt + '\n')

        except KeyError:
            pass

#close opned files
src.close()
in1.close()
outT.close()

# Close all cached file handles
for file in file_cache.values():
    if file is not None:
        file.close()
