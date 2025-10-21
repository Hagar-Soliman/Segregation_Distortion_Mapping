#HS: this script was edited in 10/21/2025 for Map 3A and to try to solve the math domain error I keep getting 

# likelihoods based on window calls as input
# all plants for a scaff together in all*txt
# reports likelihood of observed states in forward_backward

#v2:: bounded optimization of error rates
#v3:: mask entire markers based on Geno.summary


from scipy import optimize
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from math import exp,log
import sys

plantID=sys.argv[1]

#genotyping error probs
zy=0.00001 # edge for bounds
rbp = 0.1/1000000.0  # recombination rate per bp (morgans / megabase)


def calc_v0(e_rates):
    def scipy_ln_like0(x):
        # Compute log-likelihood; if LL returns NaN or non-finite, return a large positive value
        val = LL(x)
        if val is None:
            return 1e300
        try:
            # catch NaN or infinite values
            if val != val or val == float('inf') or val == float('-inf'):
                return 1e300
        except Exception:
            return 1e300
        return -val

    bounds = [(zy, 0.5), (zy, 0.5), (zy, 1.0 - zy)]
    best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like0, e_rates, approx_grad=True, bounds=bounds)
    solution = list(best)
    ln_l = -scipy_ln_like0(solution)
    solution.append(ln_l)
    return solution



def foward_backward(obs, states, start_p, transition_probability, er):
    # Forward-backward with numeric guards
    epsilon = 1e-200

    L = len(obs)
    alpha = [dict() for _ in range(L)]
    lnFactor = 0.0

    # initialize alpha[0] and ensure non-zero values, then normalize
    for y in states:
        alpha[0][y] = max(start_p.get(y, 0.0) * emission_probability(y, obs[0], er), epsilon)

    # normalize alpha[0]
    s0 = alpha[0]['AA'] + alpha[0]['AB'] + alpha[0]['BB']
    s0 = max(s0, epsilon)
    for y in states:
        alpha[0][y] /= s0
    lnFactor += log(s0)

    # forward pass
    for t in range(1, L):
        for y in states:
            acc = 0.0
            for y0 in states:
                # multiply previous scaled alpha, transition and emission
                acc += alpha[t - 1][y0] * transition_probability[t - 1][y0][y] * emission_probability(y, obs[t], er)
            # guard acc
            alpha[t][y] = max(acc, 0.0)

        # compute normalizer; use epsilon if necessary
        normalizer = max(alpha[t]['AA'], alpha[t]['AB'], alpha[t]['BB'], epsilon)
        # if normalizer not finite, return a safe failure
        if not (normalizer == normalizer and normalizer != float('inf')):
            # return early with an indicator that LL is invalid
            return alpha, None, float('nan')

        lnFactor += log(normalizer)
        for y in states:
            alpha[t][y] = alpha[t][y] / normalizer

    # final sum and LLobs
    final_sum = alpha[L - 1]['AA'] + alpha[L - 1]['AB'] + alpha[L - 1]['BB']
    final_sum = max(final_sum, epsilon)
    LLobs = lnFactor + log(final_sum)

    # backward pass
    beta = [dict() for _ in range(L)]
    for y in states:
        beta[L - 1][y] = 1.0

    for t in range(L - 2, -1, -1):
        for y in states:
            acc = 0.0
            for y0 in states:
                acc += beta[t + 1][y0] * transition_probability[t][y][y0] * emission_probability(y0, obs[t + 1], er)
            beta[t][y] = max(acc, 0.0)

        normalizer = max(beta[t]['AA'], beta[t]['AB'], beta[t]['BB'], epsilon)
        if not (normalizer == normalizer and normalizer != float('inf')):
            return alpha, None, float('nan')
        for y in states:
            beta[t][y] = beta[t][y] / normalizer

    return alpha, beta, LLobs


def emission_probability(genotype, calledG, x):  # cc [ AA, AB, BB, NN ]
    # clamp probabilities to avoid 0 or >1 values
    epsilon = 1e-200

    e1 = x[0]  # probability of sequencing error to het
    e2 = x[1]  # probability of sequencing error to homozygous opposite
    beta = x[2]  # probability of true het miscall

    # default fallback
    prob = epsilon

    if calledG == 'NN':
        return 1.0

    elif calledG == 'AA':
        if genotype == 'AA':
            prob = 1 - e1 - e2
        elif genotype == 'AB':
            prob = beta / 2.0
        elif genotype == 'BB':
            prob = e2 + 0.01

    elif calledG == 'AB':
        if genotype == 'AA' or genotype == 'BB':
            prob = e1
        elif genotype == 'AB':
            prob = 1 - beta

    elif calledG == 'BB':
        if genotype == 'AA':
            prob = e2
        elif genotype == 'AB':
            prob = beta / 2.0
        elif genotype == 'BB':
            prob = 1 - e1 - e2

    # clamp into safe numeric range
    prob = max(min(prob, 1.0), epsilon)
    return prob

 
def LL(x):

	Total_LL=0.0
	for v1s in v1scaffs:
		total_snps=v1scaffs[v1s] # updated for each scaff
		# transition probs a global
		transition_probability=[{} for j in xrange(total_snps-1)] # global that is updated within LL(x)
		for x1 in xrange(total_snps-1): # v1 scaff
			dist=abs(Position[plantID][v1s][x1+1]-Position[plantID][v1s][x1])
			r = rbp*float(dist)
			transition_probability[x1] ={'BB': {'BB': (1 - r), 'AB': r, 'AA': 0.0}, 'AB': {'BB': r, 'AB': (1 - r), 'AA': r}, 'AA': {'BB': 0.0, 'AB': r, 'AA': (1 - r)} }
	
		if Gcalls[v1s]>0:
			fprbs,rprbs,llx=foward_backward(obsA[plantID][v1s],states,start_probability,transition_probability,x)
			#print v1s,Gcalls[v1s],"LL= ",llx
		#print "forward ",fprbs
		#print "backward ",rprbs
		#postProb=[{} for j in range(len(obsA[plantID][v1s]))] # forward:: alpha[j][X] is probability that true genotye is X at marker j (starts at 0)

		#for j in range(len(fprbs)):
		#	denom=0.0
		#	for y in states: 
		#		denom+=(fprbs[j][y]*rprbs[j][y])
		#	for y in states: 
		#		postProb[j][y]=(fprbs[j][y]*rprbs[j][y])/denom

		#print postProb
		else:
			llx=0.0

		Total_LL+=llx
	#print x,Total_LL
	return Total_LL


####################################################################################################
### Main Program

states = ('AA','AB','BB')
start_probability = {'AA':0.1,'AB':0.4,'BB':0.5}


inZ = open("bad.marks.txt","rU")
badmark={}
for line_idx, line in enumerate(inZ):
	cols = line.replace('\n', '').split('\t') 

# 103a	100000	
	key=cols[0]+"_"+cols[1]
	badmark[key]=1


Position={}
obsA={}
v1scaffs={}
Gcalls={}
cscaff=''
calls_total=0
src  =open("g.Genotypes."+plantID+".bam.txt", "rU")
for line_idx, line in enumerate(src):
	cols = line.replace('\n', '').split('\t') 
# isg480	1	400000	AB
	key=cols[1]+"_"+cols[2]
	try:
		uck=badmark[key]
		#print "suppressing ", key
	except KeyError:

		if plantID!=cols[0].replace('.bam', ''):
			print "Whoa",plantID,cols[0]
		if line_idx==0:
			Position[plantID]={}
			obsA[plantID]={}


		if cols[1] !=cscaff: # new scaff
			Position[plantID][cols[1]]=[]
			obsA[plantID][cols[1]]=[]
			cscaff=cols[1]
			v1scaffs[cols[1]]=0
			Gcalls[cols[1]]=0

		Position[plantID][cols[1]].append(int(cols[2]))
		obsA[plantID][cols[1]].append(cols[3])
		v1scaffs[cols[1]]+=1 # will need to be updated if you do more than one plant in a run
		if cols[3] != 'NN':
			Gcalls[cols[1]]+=1
			calls_total+=1

#initial values for e1,e2,beta
e_rates=[0.01, 0.01,0.01]

zsol= calc_v0(e_rates) 

print plantID,calls_total,zsol[0],zsol[1],zsol[2],zsol[3]



