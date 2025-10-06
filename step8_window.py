from math import exp,log
import sys



### Main Program

LineID = sys.argv[1]  # Aim3.8.7.1.bam.m.txt
name = sys.argv[2]  # 8.7.1
output_directory = sys.argv[3]
out3 = open(f"{output_directory}/Genostats.{name}.txt", "w")
out2 = open(f"{output_directory}/Genotypes.{name}.txt", "w")
out1 = open(f"{output_directory}/window.{name}.txt", "w")
src = open(LineID, "r")

# Add headers
out1.write("Scaffold\tWindowStart\tNumSNPs\tREF_Reads\tALT_Reads\n")
out2.write("Sample\tScaffold\tWindowStart\tGenotype\n")
out3.write("Sample\tAA_Count\tAB_Count\tBB_Count\tNN_Count\n")


blacklist={}
if len(sys.argv)>4:
	inb =open(sys.argv[4], "r")
	for line_idx, line in enumerate(inb):
		cols = line.replace('\n', '').split('\t') 
		blacklist[cols[0] + "_" + cols[1]] = 1


MinCount=5
hetdev=0.2
WindowSize=100000
MinSNPs_perwindow=1

snplist=[]
obs=[0,0] 

plant=0
sID=0

ghu=[0,0,0,0]
for line_idx, line in enumerate(src):
	cols = line.replace('\n', '').split('\t') 
	#1	588	0	1
	scaff=cols[0]
	pos = int(cols[1])
	window=int(pos/WindowSize)
	if line_idx==0:
		Cscaff=cols[0]
		Cwin=int(pos/WindowSize)
		cx=[0,0,0]
	else:
		if scaff != Cscaff or window != Cwin: # output last window
			ntot=cx[0]+cx[1]
			out1.write(Cscaff+'\t'+str(Cwin*WindowSize)+'\t'+str(cx[2])+'\t'+str(cx[0])+'\t'+str(cx[1])+'\n')
			if ntot>=MinCount and cx[2]>=MinSNPs_perwindow: #rule set for hard calls
				pR=float(cx[0])/float(ntot)
				if pR>0.95:
					out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tAA\n') 
					ghu[0]+=1
				elif pR<0.05:
					out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tBB\n')
					ghu[2]+=1 
				elif pR>=(0.5-hetdev) and pR<=(0.5+hetdev):
					out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tAB\n') 
					ghu[1]+=1
				else:
					out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tNN\n') 
					ghu[3]+=1

			else:
				out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tNN\n') 
				ghu[3]+=1

			Cscaff=cols[0]
			Cwin=window
			cx=[0,0,0]
			


	try:
		amIbad=blacklist[cols[0]+"_"+cols[1]]
	except KeyError: # OK
		cx[0]+=int(cols[2]) #IM alleles in window
		cx[1]+=int(cols[3]) #Alt alleles in window
		if int(cols[2])+int(cols[3])>0:
			cx[2]+=1


ntot = cx[0] + cx[1]
out1.write(Cscaff+'\t'+str(Cwin*WindowSize)+'\t'+str(cx[2])+'\t'+str(cx[0])+'\t'+str(cx[1])+'\n') 
if ntot>=MinCount and cx[2]>=MinSNPs_perwindow: #rule set for hard calls
	pR=float(cx[0])/float(ntot)
	if pR>0.95:
		out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tAA\n') 
		ghu[0]+=1
	elif pR<0.05:
		out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tBB\n')
		ghu[2]+=1 
	elif pR>=(0.5-hetdev) and pR<=(0.5+hetdev):
		out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tAB\n') 
		ghu[1]+=1
	else:
		out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tNN\n') 
		ghu[3]+=1

else:
	out2.write(name+'\t'+Cscaff+'\t'+str(Cwin*WindowSize)+'\tNN\n') 
	ghu[3]+=1

out3.write(name+'\t'+str(ghu[0])+'\t'+str(ghu[1])+'\t'+str(ghu[2])+'\t'+str(ghu[3])+'\n') 

out1.close()
out2.close()
out3.close()
src.close()
if len(sys.argv) > 4:
    inb.close()
