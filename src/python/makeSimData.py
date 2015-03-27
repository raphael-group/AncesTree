#!/usr/bin/python

#Load required modules
import sys, os, argparse, math, random
import numpy as np

from convertData import collapse_mat, get_sample_labels, get_mut_labels, print_new_file

# Load and parse arguments
def parse_arguments():
	"""
	Parse command line arguments

	Returns:

	"""

	description = "Create simulated data for AncesTree Algorithm."
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("-m","--NUM_MUT",required=False, metavar="NUM_MUT",
		type=int, default=20, help="The number of mutations to simulate.")
	parser.add_argument("-n","--NUM_CLONE",required=False, metavar="NUM_CLONE",
		type=int, default=20, help="The number of clones to simulate.")
	parser.add_argument("-c", "--COV", required=False, metavar="COV",
		type=int, default=100, help="The coverage to simulate.")
	parser.add_argument("-s","--SAMPLES",required=False, metavar="SAMPLES",
		type=int, default=4, help="The number of samples to simulate.")
	parser.add_argument("--MAX", required=False, metavar="MAX",
		type=int, default=4, help="The maximum number of clones to mix.")
	parser.add_argument("--MIN", required=False, metavar="MIN",
		type=int, default=1, help="The minimum number of clones to mix.")
	parser.add_argument("-o","--OUT_DIR", required=False, metavar="OUT_DIR",
		default="./", help="The directory to write results to.")
	parser.add_argument("-d","--DRAWS", required=False, metavar="DRAWS",
		type=int, default=1, help="The number of datasets to make.")
	parser.add_argument("-r","--RAND", required=False, metavar="RAND",
		type=int, default=47, help="Random seed.")
	parser.add_argument("-p","--PERFECT", required=False,
		default=False, action = "store_true")
	args = parser.parse_args()

	num=args.NUM_MUT
	clones=args.NUM_CLONE
	if (clones > num): clones = num
	cov=args.COV
	samples=args.SAMPLES
	minNum=args.MIN
	maxNum=args.MAX
	directory=args.OUT_DIR
	draws=args.DRAWS
	rand=args.RAND
	perfect=args.PERFECT

	print "NUM MUT: " + str(num)
	print "NUM CLONES: " + str(clones)
	print "COV: " + str(cov)
	print "SAMPLES: " + str(samples)
	print "MIN: " + str(minNum)
	print "MAX: " + str(maxNum)
	print "OUT_DIR: " + directory
	print "DRAWS: " + str(draws)
	print "RAND SEED: " + str(rand)
	print "PERFECT DATA: " + str(perfect)


	return num, clones, cov, samples, minNum, maxNum, directory, draws, rand, perfect


def make_sim(num, samples, minNum, maxNum):
	"""
	Makes a simulated dataset.
	"""
	B, parents=generate_B(num)
	U=generate_U(num,samples, minNum, maxNum)

	return U, B

def make_sim_all_leaves(num, clones, samples, minNum, maxNum):
	"""
		Makes a simulation where all leaves are sampled.
	"""
	maxMixings = maxNum*samples

	valid = False
	while not valid:
		B, parents=generate_B_with_cluster(num,clones)
		leaves = get_leaves(parents, clones)

		#if maxMixings >= len(leaves): valid = True
		if maxMixings >= 2*len(leaves): valid = True

	U = generate_U_with_leaves(clones, samples, minNum, maxNum, leaves)

	return U, B


def get_leaves(parents, num):
	leaves = range(num)

	for val in parents.keys():
		par = parents[val]
		if par in leaves: leaves.remove(par)

	return leaves


def sample_reads(num, samples, cov, U, B):
	"""
	Creates R (ref) and V (variant) read counts
	"""
	T = [[np.random.poisson(cov) for i in xrange(num)] for i in xrange(samples)]
	F = 0.5 * np.dot(U,B)


	R = [[0 for i in xrange(num)] for i in xrange(samples)]
	V = [[0 for i in xrange(num)] for i in xrange(samples)]

	for i in xrange(samples):
		for j in xrange(num):
			v = np.random.binomial(T[i][j],F[i][j])
			r = T[i][j] - v
			R[i][j] = r
			V[i][j] = v

	return R, V

def sample_perfect_reads(num, samples, cov, U, B):
	"""
	Creates R(ref) and V (variant) read counts for perfect data.
	"""
	T = [[cov for i in xrange(num)] for i in xrange(samples)]
	F = 0.5 * np.dot(U,B)

	R = [[0 for i in xrange(num)] for i in xrange(samples)]
	V = [[0 for i in xrange(num)] for i in xrange(samples)]

	for i in xrange(samples):
		for j in xrange(num):
			v = int(round(T[i][j]*F[i][j]))
			r = T[i][j] - v
			R[i][j] = r
			V[i][j] = v

	return R, V



def generate_B(num):
	"""
	Generate a true matrix B.
	"""
	B = [[0 for i in xrange(num)] for i in xrange(num)]
	order=range(num)
	random.shuffle(order)

	parents = dict() #keep track of the parents of each "mutation"
	
	#Initialize root
	B[order[0]][order[0]] = 1
	parents[order[0]] = -1

	#Build B by randomly picking a parent node from the tree
	for row in [i+1 for i in xrange(num-1)]:
		
		p = random.choice(xrange(row))
		B[order[row]][:] = B[order[p]][:]
		B[order[row]][order[row]] = 1 #add a new one based on the order
		parents[order[row]] = order[p]

	return B, parents


def generate_B_with_cluster(num,clones):
	"""
	Generate a true matrix B with the specified number of mutations and clones (induces clustering).
	"""
	B = [[0 for i in xrange(num)] for i in xrange(clones)]

	order=range(clones)
	random.shuffle(order)
	clone_map = get_clone_mut_map(num, clones)
	

	parents = dict() #keep track of the parents of each "clone"
	
	#Initialize root
	for mut in clone_map[order[0]]:
		B[order[0]][mut] = 1

	parents[order[0]] = -1

	#Build B by randomly picking a parent node from the tree
	for row in [i+1 for i in xrange(clones-1)]:
		
		p = random.choice(xrange(row))
		B[order[row]][:] = B[order[p]][:]
		for mut in clone_map[order[row]]:
			B[order[row]][mut] = 1 #add a new one based on the order

		parents[order[row]] = order[p]

	return B, parents


def get_clone_mut_map(num, clones):
	"""
	Returns a random map from clones to mutations (effectively clustering)
	"""
	map = dict()

	breaks = [0]
	poss_starts = [x+1 for x in xrange(num-1)]

	for c in xrange(clones-1):
		val = random.choice(poss_starts)
		breaks.append(val)
		poss_starts.remove(val)
	breaks.append(num)

	breaks.sort()
	for c in xrange(clones):
		cur = range(breaks[c],breaks[c+1])
		map[c] = cur

	return map


def generate_U(num, samples, minNum, maxNum):
	"""
		Generate a U matrix
	"""
	U = [[0 for i in xrange(num)] for i in xrange(samples)]

	#update row by row
	for i in xrange(samples):
		numMix = random.randint(minNum, maxNum)
		draw = sample_simplex(numMix+1) # random draw from the simplex

		indices = [j for j in xrange(num)] #takes random indices to be mixed
		random.shuffle(indices)
		for j in range(numMix):

			U[i][indices[j]] = draw[j]
	
	return U

def generate_U_with_leaves(num, samples, minNum, maxNum, leaves):
	"""
		Generates a U matrix that uses the specified leaves.
		Assumes that the number of leaves is such that they can all
		be used given the current constraints.
	"""
	valid = False
	while not valid:
		U = generate_U(num, samples, minNum, maxNum)
		#valid = check_valid_U(U,leaves)
		valid = check_valid_U_twice(U,leaves)

	return U


def check_valid_U(U, leaves):
	valid = True
	samples = len(U)
	for leaf in leaves:
		s = 0
		for i in xrange(samples): s += U[i][leaf]

		if s == 0 : return False

	return valid

def check_valid_U_twice(U, leaves):
	valid = True
	samples = len(U)
	for leaf in leaves:
		s = 0
		for i in xrange(samples): 

			if U[i][leaf] > 0: s += 1

		if s < 2 : return False

	return valid

	
def sample_simplex(num):
	"""
	Make a uniform draw from the simplex. (Uses fact that draws are equivalent to partitions of the unit interval)
	"""
	X = [0 for i in xrange(num+1)]
	for i in [x+1 for x in xrange(num-1)]:
		X[i] = random.random()
	X[num] = 1
	X.sort()

	Z = [0 for i in xrange(num)]
	for i in xrange(num):
		Z[i] = X[i+1]-X[i]

	return Z


def save_input(R, V, directory, i):
	filename = "sim_" + str(i) + ".input"
	path = os.path.join(directory,filename)

	samples = len(R)
	num = len(R[0])

	f = open(path, 'w')
	header = ['gene_id']
	for i in xrange(samples):
		#header.append('ref_'+ str(i))
		#header.append('var_'+ str(i))
        # MEK: got rid of ref_/var_ prefix
		header.append(str(i))
		header.append(str(i))

	f.write("\t".join(header) + "\n")

	for i in xrange(num):
		row = ['Mut_'+ str(i)]
		for j in xrange(samples):
			row.append(str(R[j][i]))
			row.append(str(V[j][i]))
		f.write("\t".join(row) + "\n")

	f.close()


def save_true(U,B,F,directory, i):
	filename = "sim_" + str(i) + ".true"
	path = os.path.join(directory,filename)

	samples = len(U)
	clones = len(U[0])
	muts = len(B[0])

	f = open(path,'w')

	f.write(str(samples) +"\n")
	f.write(str(clones)+ "\n")

	#Write U
	for j in xrange(samples):
		f.write(" ".join([str(x) for x in U[j]]) + "\n")

	f.write("\n")

	f.write(str(clones) +"\n")
	f.write(str(muts)+ "\n")

	for j in xrange(clones):
		f.write(" ".join([str(x) for x in B[j]]) + "\n")

	f.write("\n")

	f.write(str(samples) +"\n")
	f.write(str(clones)+ "\n")

	for j in xrange(samples):
		f.write(" ".join([str(x) for x in F[j]]) + "\n")

	f.close()



# Main method
def main():

	#Load arguments
	num, clones, cov, samples, minNum, maxNum, directory, draws, seed, perfect = parse_arguments()

	random.seed(seed)

	#Create the datasets
	for i in xrange(draws):
		U, B = make_sim_all_leaves(num, clones, samples, minNum, maxNum)
		F = 0.5 * np.dot(U,B)	
		
		if perfect:
			R,V = sample_perfect_reads(num, samples, cov, U, B)
		else:
			R,V = sample_reads(num, samples,cov, U, B)

		col_B, B_labels = collapse_mat(B)
		#col_F, F_labels = collapse_mat(F)
		col_F = 0.5*np.dot(U,col_B)

		sample_labels = get_sample_labels(F)
		mut_labels = get_mut_labels(F)

		filename = "sim_" + str(i) + ".true"
		path = os.path.join(directory,filename)
		print_new_file(path, U, B, F, col_B, B_labels, col_F, sample_labels, mut_labels)

		save_input(R,V,directory,i)

		#save_true(U,B,F,directory, i)



if __name__ == '__main__':
	main()
