#USAGE :  python python_mdanalysis_skeleton.py [config file name]

#DEPENDCIES : numpy, MDAnalysis, math

#CONFIG FILE FORMAT:
#   TopFile = [topology file name (prmtop file)]
#   TrajFile = [trajectory file name (mdcrd file)]
#   OutFile = [output data file name]

import sys
import os
import numpy as np
import math
import MDAnalysis
from scipy import interpolate

kB = 0.001987 # kcal/mol/K
thresh = 1E-4

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global top_file, file_update, sim_in_file, atom_data_root, ib_lambda, n_iter, param_out_file, kT, n_start
	f = open(cfg_file)
	file_update = ""
        n_start = 0
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			# check value
			if option.lower()=='topfile':
				top_file = value
			elif option.lower()=='paramfile':
				param_out_file = value
			elif option.lower()=='siminfile':
				sim_in_file = value
			elif option.lower()=='atomisticdata':
				atom_data_root = value
			elif option.lower()=='lambda':
				ib_lambda = float(value)
			elif option.lower()=='iterations':
				n_iter = int(value)
			elif option.lower()=='restartnum':
				n_start = int(value)
			elif option.lower()=='fileupdate':
				file_update = value
			elif option.lower()=='temperature':
				T = float(value)
			else :
				print "Option:", option, " is not recognized"
	f.close()
	if file_update != "true" and file_update != "false":
		file_update = "true"
		print "file_update being assigned default value of true"
	kT = kB*T
        n_iter = n_start + n_iter

def ParsePsfFile(psf_file):
	global n_uniq_atom_types, atom_types, uniq_bond_atom_types, n_uniq_bonds, bonds, n_bonds
	f = open(psf_file)
	atom_flag = bond_flag = "false"
	atom = bond = angle = dihedral = 0
	bond_count = 0
	atom_types = []
	bonds = []
	uniq_bond_atom_types = []
	angles = []
	uniq_angle_atom_types = []
	dihedrals = []
	uniq_dihedral_atom_types = []
	n_uniq_atom_types = 0
	n_uniq_bonds = 0
	n_uniq_angles = 0
	n_uniq_dihedrals = 0
	for line in f:
		if atom_flag == "true" and atom < n_atoms:
			atom_types.append(line[29:33].strip())
			same = "false"
			if atom > 0:
				for atom2 in range(atom):
					if atom_types[atom] == atom_types[atom2]:
						same = "true"
						break
			if same == "false":
				n_uniq_atom_types += 1
			atom += 1
		elif bond_flag == "true" and bond < n_bonds:
			for line_bond in range(4):
				bonds.append([])
				bonds[bond].append(line[line_bond*16:line_bond*16+8])
				bonds[bond].append(line[line_bond*16+8:line_bond*16+16])
				same = "false"
				if bond > 0:
					for bond2 in range(bond):
						if (atom_types[int(bonds[bond][0])-1] == atom_types[int(bonds[bond2][0])-1] and atom_types[int(bonds[bond][1])-1] == atom_types[int(bonds[bond2][1])-1]) or (atom_types[int(bonds[bond][0])-1] == atom_types[int(bonds[bond2][1])-1] and atom_types[int(bonds[bond][1])-1] == atom_types[int(bonds[bond2][0])-1]):
							same = "true"
							uniq_bond_num = bonds[bond2][2]
							break
				if same == "false":
					uniq_bond_atom_types.append([])
					uniq_bond_atom_types[n_uniq_bonds].append(atom_types[int(bonds[bond][0])-1])
					uniq_bond_atom_types[n_uniq_bonds].append(atom_types[int(bonds[bond][1])-1])
					uniq_bond_num = str(n_uniq_bonds)
					n_uniq_bonds += 1
				bonds[bond].append(uniq_bond_num)
				bond += 1
				bond_count += 1
				if bond == n_bonds:
					break
		if line[9:15] == "!NATOM":
			n_atoms = int(line[0:9])
			atom_flag = "true"
		elif line[9:15] == "!NBOND":
			n_bonds = int(line[0:9])
			bond_flag = "true"
		if line[9:15] == "!NTHET":
                        break

	f.close()
	bonds = np.asmatrix(bonds,dtype=int)

# add to bond distance histograms 
def ComputeBondHists(atom_positions, bonds, bond_hists, bond_min, bond_delta):
	# get sizes etc
	n_bonds = bonds.shape[0]
	n_bins = bond_hists.shape[1]
	for bond in range(n_bonds):

		dist = math.sqrt(ComputeDist2(atom_positions[bonds[bond,0]-1,:],atom_positions[bonds[bond,1]-1,:]))
		dist_bin = int((dist-bond_min)/bond_delta)
		if dist_bin >= 0 and dist_bin < n_bins:
			bond_hists[bonds[bond,2],dist_bin] += 1

# average bond distance histograms (convert them to probability densities)
def ReadBondHists(bond_hists, bond_min, bond_delta, uniq_bond_atom_types):
	global file_update

	n_bonds = bond_hists.shape[0]
	n_bins = bond_hists.shape[1]

	for bond in range(n_bonds):
		filename = "../../rna_params_v02/bonds/"+uniq_bond_atom_types[bond][0].strip()+"_"+uniq_bond_atom_types[bond][1].strip()+".hist"
		filename2 = "../../rna_params_v02/bonds/"+uniq_bond_atom_types[bond][1].strip()+"_"+uniq_bond_atom_types[bond][0].strip()+".hist"
		if os.path.isfile(filename):
			if file_update == "true":
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					bond_hists[bond,count] += float(val)
					count += 1
				inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			if file_update == "true":
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					bond_hists[bond,count] += float(val)
					count += 1
				inp.close()
		else :
			if uniq_bond_atom_types[bond][0][1] == "1":
				type1 = uniq_bond_atom_types[bond][0][0] + str(2)
			else:
				type1 = uniq_bond_atom_types[bond][0][0] + str(1)
			if uniq_bond_atom_types[bond][1][1] == "1":
				type2 = uniq_bond_atom_types[bond][1][0] + str(2)
			else:
				type2 = uniq_bond_atom_types[bond][1][0] + str(1)
			filename = "../../rna_params_v02/bonds/"+type1+"_"+type2+".hist"
			filename2 = "../../rna_params_v02/bonds/"+type2+"_"+type1+".hist"
			if os.path.isfile(filename):
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					bond_hists[bond,count] = float(val)
					count += 1
				inp.close()
			elif os.path.isfile(filename2):
				filename = filename2
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					bond_hists[bond,count] = float(val)
					count += 1
				inp.close()
			else:
				print "Did not find prameters for:", uniq_bond_atom_types[bond]
				sys.exit()


# compute the dihedral between four points
def ComputeDih(r1,r2,r3,r4):

	r12 = np.empty(3,dtype=float)
	r23 = np.empty(3,dtype=float)
	r34 = np.empty(3,dtype=float)

	# compute two vectors (r1-r2 and r3-r2) making sure they are in same box
	for j in range(0,3):
		temp1 = r1[j]-r2[j]
		temp2 = r2[j]-r3[j]
		temp3 = r3[j]-r4[j]
		r12[j] = temp1
		r23[j] = temp2
		r34[j] = temp3
	
	A = np.cross(r12,r23)
	A /= math.sqrt(np.dot(A,A))
	B = np.cross(r23,r34)
	B /= math.sqrt(np.dot(B,B))
	C = np.cross(r23,A)
	C /= math.sqrt(np.dot(C,C))
	cos_phi = np.dot(A,B)
	sin_phi = np.dot(C,B)
	
	phi = -math.atan2(sin_phi,cos_phi)*180.0/3.1415926535
	return phi

# compute the angle between three points
def ComputeAng(r1,r2,r3):

	r21 = np.empty(3,dtype=float)
	r23 = np.empty(3,dtype=float)

	# compute two vectors (r1-r2 and r3-r2) making sure they are in same box
	for j in range(0,3):
		temp1 = r1[j]-r2[j]
		temp2 = r3[j]-r2[j]
		r21[j] = temp1
		r23[j] = temp2
	
	theta = math.acos(np.dot(r21,r23)/(math.sqrt(np.dot(r21,r21))*math.sqrt(np.dot(r23,r23))))*180.0/3.1415926535
	return theta

# compute the distance between two points taking into account periodic boundary conditions
def ComputeDist2(r1,r2):
	dist2 = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		dist2 += temp*temp

	return dist2;

# compute new CG probability distributions from CG trajectory
def AnalyzeCGTraj(top_file, traj_file, cg_bond_potentials):
	global bonds, bond_min, bond_delta, ib_iter
        # initialize histogram arrays
	n_bonds = cg_bond_potentials.shape[0]
	n_bond_bins = cg_bond_potentials.shape[1]
        cg_bond_hists = np.zeros((n_bonds,n_bond_bins),dtype=float)
        cg_bond_forces = np.zeros((n_bonds,n_bond_bins),dtype=float)
	# initiate MDAnalysis coordinate universe
	coord = MDAnalysis.Universe(top_file, traj_file)
	# make an atom selection
	sel = coord.select_atoms("all")
	n_frames = coord.trajectory.n_frames
	print "Number of frames in trajectory file: ", n_frames
	# declare array that will contain the selected atom positions at each time frame
	positions = np.empty((sel.n_atoms,3),dtype=float)

	# Loop through trajectory
	for ts in coord.trajectory:

		# track progress
		if ts.frame % 10 == 0:
			print "Frame: ", ts.frame, " of ", n_frames
		
		# save selected positions into array
		positions = sel.positions
		
		# add to bond histograms
		ComputeBondHists(positions, bonds, cg_bond_hists, bond_min, bond_delta)
        filename = "bond" + str(ib_iter) + "_out.hist"
        WriteBondHists(cg_bond_hists, bond_min, bond_delta, filename)
        CreateBondPotentials(cg_bond_hists, cg_bond_potentials, cg_bond_forces, bond_min, bond_delta)

def ReadAtomHists(atom_bond_potentials):
	global n_uniq_bonds, n_bond_bins, bond_min, bond_delta, param_out_file
        # allocate hist and force arrays
        atom_bond_hists = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float) 
        atom_bond_forces = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float) 
        # read hists
	ReadBondHists(atom_bond_hists, bond_min, bond_delta, uniq_bond_atom_types)
        # open param file
	param_out = open(param_out_file, 'w')
        # create potentials from hists
        CreateBondPotentials(atom_bond_hists, atom_bond_potentials, atom_bond_forces, bond_min, bond_delta)
        WriteBondPotentials(atom_bond_potentials, atom_bond_forces, bond_min, bond_delta,"bond0.ib","bond0.ener")
        WriteBondHists(atom_bond_hists, bond_min, bond_delta, "bond0_in.hist")
	param_out.close()


# average bond distance histograms (convert them to probability densities)
def CreateBondPotentials(bond_hists, bond_potentials, bond_forces, bond_min, bond_delta):
	global kT

	n_bonds = bond_hists.shape[0]
	n_bins = bond_hists.shape[1]
	coeff_mat = np.empty((n_bins,3),dtype=float)
        # zero out array
        bond_potentials.fill(0.0)

        # make coefficient matrix
	for i in range(n_bins):
		x = bond_min+(i+0.5)*bond_delta
		coeff_mat[i,0] = 1.0
		coeff_mat[i,1] = x
		coeff_mat[i,2] = x*x

	for bond in range(n_bonds):
		# create probability density
		bond_hists[bond,:] /= (np.sum(bond_hists[bond,:])*bond_delta)
		# convert to energies
		min_val = 0
		count = 0
		for i in range(n_bins):
			x = bond_min+(i+0.5)*bond_delta
                        x2 = x*x
			if bond_hists[bond,i]/x2 > thresh :
				bond_potentials[bond,i] = -kT*math.log(bond_hists[bond,i]/x2)

		# find the minimum energy
		min_val = np.amin(bond_potentials[bond,:])
                # find start
                for i in range (n_bins):
                        if bond_potentials[bond,i] != 0:
				start = i
                                break
                # find stop
                for i in range(n_bins-1,-1,-1):
			if bond_potentials[bond,i] != 0:
				stop = i + 1 # because python is non-inclusive of upperbound
                                break

		# now smooth the potential using a cubic spline fit	
                tck = interpolate.splrep(coeff_mat[start:stop,1],bond_potentials[bond,start:stop]-min_val,s=3)
                bond_potentials[bond,start:stop] = interpolate.splev(coeff_mat[start:stop,1], tck, der=0)

		# fit previous section if we have enough data
		if (stop-start) > 4:
		    # fit function to parabola to extrapolate to short and long distances
		    k, rss, rank, sv = np.linalg.lstsq(coeff_mat[start:stop],bond_potentials[bond,start:stop])
                    # place fit data into potential
	    	    for j in range(start):
		        bond_potentials[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
		    for j in range(stop,n_bins):
			bond_potentials[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]


		# now smooth the total  potential using a cubic spline fit	
		tck = interpolate.splrep(coeff_mat[:,1],bond_potentials[bond,:]-min_val,s=5)
                bond_potentials[bond,:] = interpolate.splev(coeff_mat[:,1], tck, der=0)
                bond_forces[bond,:] = -interpolate.splev(coeff_mat[:,1], tck, der=1)
               
# write bond potentials and forces to files
def WriteBondPotentials(bond_potentials, bond_forces, bond_min, bond_delta, ib_out,ener_out):

	n_bonds = bond_potentials.shape[0]
	n_bins = bond_potentials.shape[1]
	
	out = open(ib_out,'w')
	for bond in range(n_bonds):
		out.write("\nBOND_%s\n" % (str(bond+1)))
		out.write("N %5d\n\n" % (n_bins))
		for i in range(n_bins):
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,bond_min+(i+0.5)*bond_delta, bond_potentials[bond,i],bond_forces[bond,i]))
        # close ib file
        out.close()

        # write file containing all bond energies in easily graphable format
	temp = open(ener_out, 'w')
	for i in range(n_bins):
		temp.write("%20.5f" % (bond_min+(i+0.5)*bond_delta))
		for bond in range(n_bonds):
			temp.write("%20.7f" % (bond_potentials[bond,i]))
		temp.write("\n")
	temp.close()

# write bond hists
def WriteBondHists(bond_hists, bond_min, bond_delta, outfile):

	n_bonds = bond_hists.shape[0]
	n_bins = bond_hists.shape[1]
        norm = np.empty(n_bonds,dtype=float)

        for bond in range(n_bonds):
            norm[bond] = np.sum(bond_hists[bond,:])*bond_delta

        # write file containing all bond energies in easily graphable format
	temp = open(outfile, 'w')
	for i in range(n_bins):
		temp.write("%20.5f" % (bond_min+(i+0.5)*bond_delta))
		for bond in range(n_bonds):
			temp.write("%20.7f" % (bond_hists[bond,i]/norm[bond]))
		temp.write("\n")
	temp.close()


def WritePotentials(potentials, x_min, x_delta, chars, ib_out,ener_out):

	n_potentials = potentials.shape[0]
	n_bins = potentials.shape[1]

	out = open(ib_out,'w')
	for potential in range(n_potentials):
		out.write("\n%s_%s\n" % (chars,str(potential+1)))
		out.write("N %5d\n\n" % (n_bins))
		for i in range(n_bins):
			# compute force
                        if i > 0 and i < n_bins-1:
				force = -(potentials[potential,i+1]-potentials[potential,i-1])/(2*x_delta)
                        elif i == 0:
				force = -(potentials[potential,i+1]-potentials[potential,i])/(x_delta)
			elif i == n_bins-1:
				force = -(potentials[potential,i]-potentials[potential,i-1])/(x_delta)
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,x_min+(i+0.5)*x_delta, potentials[potential,i],force))
	out.close()

	temp = open(ener_out, 'w')
	for i in range(n_bins):
		temp.write("%20.5f" % (x_min+(i+0.5)*x_delta))
		for potential in range(n_potentials):
			temp.write("%20.5f" % (potentials[potential,i]))
		temp.write("\n")
	temp.close()


# MM_ANGLES
#
def UpdateBonds(ib_iter, atom_bond_potentials, prev_cg_bond_potentials, cg_bond_potentials, param_out):
	global file_update, kT, ib_lambda, bond_delta, bond_min
        
	n_bonds = cg_bond_potentials.shape[0]
	n_bins = cg_bond_potentials.shape[1]

	cg_bond_forces = np.zeros((n_bonds,n_bins),dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	ib_out = "bond" + str(ib_iter+1) + ".ib"
	ener_out = "bond" + str(ib_iter+1) + ".ener"

        for bond in range(n_bonds):
	    for i in range(n_bins):
	        cg_bond_potentials[bond,i] = prev_cg_bond_potentials[bond,i] + ib_lambda * (atom_bond_potentials[bond,i] - cg_bond_potentials[bond,i])

        # make coefficient matrix
	for i in range(n_bins):
		x_mat[i] = bond_min+(i+0.5)*bond_delta
        # let's smooth the potentials!
        for bond in range(n_bonds):
            # find minimum
            min_val = np.amin(cg_bond_potentials[bond,i])
            # now smooth the total  potential using a cubic spline fit	
            tck = interpolate.splrep(x_mat,cg_bond_potentials[bond,:]-min_val,s=5)
            cg_bond_potentials[bond,:] = interpolate.splev(x_mat, tck, der=0)
            cg_bond_forces[bond,:] = -interpolate.splev(x_mat, tck, der=1)
        # now to write them to file
        WriteBondPotentials(cg_bond_potentials,cg_bond_forces,bond_min, bond_delta,ib_out,ener_out)

# MM MM MM

def UpdatePotentials(ib_iter, atom_bond_potentials, prev_cg_bond_potentials, cg_bond_potentials):
	global file_update, kT, ib_lambda, param_out_file, bond_delta

	param_out = open(param_out_file, 'w')
	UpdateBonds(ib_iter, atom_bond_potentials, prev_cg_bond_potentials, cg_bond_potentials, param_out)
        param_out.close()

# read restart potentials
def ReadRestartPotentials(prev_cg_bond_potentials,ib_iter):
	n_bonds = prev_cg_bond_potentials.shape[0]
	n_bins = prev_cg_bond_potentials.shape[1]

#	ib_out = "bond" + str(ib_iter+1) + ".ib"
	ener_in = "bond" + str(ib_iter) + ".ener"
	temp = open(ener_in, 'r')
        bin_count = 0
        for line in temp:
		for bond in range(n_bonds):
                    prev_cg_bond_potentials[bond,bin_count] = float(line[20*(bond+1):20*(bond+2)])
                bin_count += 1
	temp.close()


# run CG simulation using lammps
def RunCGSim(ib_iter, equil_in_file, run_in_file):
	
	new_equil_file = "equil.iter" + str(ib_iter) + ".in"
	new_run_file = "run.iter" + str(ib_iter) + ".in"
	equil_log_file = "equil.iter" + str(ib_iter) + ".log"
	run_log_file = "run.iter" + str(ib_iter) + ".log"
        command0 = "cp bond" + str(ib_iter) + ".ib bonds.ib"
	command1 = "sed -e s/NUM/" + str(ib_iter) + "/ < " + equil_in_file + " > " +  new_equil_file
	command2 = "mpirun -np 3 lmp_mac_mpi -i " + new_equil_file + " > " + equil_log_file
	command3 = "sed -e s/NUM/" + str(ib_iter) + "/ < " + run_in_file + " > " +  new_run_file
	command4 = "mpirun -np 3 lmp_mac_mpi -i " + new_run_file + " > " + run_log_file
	
	print command0
	os.system(command0)
	print command1
	os.system(command1)
	print command2
	os.system(command2)
	print command3
	os.system(command3)
	print command4
	os.system(command4)
	print "Done with CG simulation for iteration", ib_iter

#############################################################################################################################################################
####################################################              MAIN PROGRAM               ################################################################
#############################################################################################################################################################

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Topology file:", top_file

# parse psf file
ParsePsfFile(top_file)

print "Number of unique atom types:", n_uniq_atom_types
print "Number of bonds:", n_bonds
print "Number of unique bonds:", n_uniq_bonds

# declare bond, angle and dihedral histograms
bond_min = 0.0
bond_max = 24.0
bond_delta = 0.25
n_bond_bins  = int((bond_max-bond_min)/bond_delta)
atom_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
cg_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
prev_cg_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)

# read atomistic probability distributions and generate initial (inverse Boltzmann) potentials
ReadAtomHists(atom_bond_potentials)
# set previous potentials to atom potentials
prev_cg_bond_potentials[:,:] = atom_bond_potentials[:,:]

# read restrart if needed
if n_start > 0:
    ReadRestartPotentials(prev_cg_bond_potentials,n_start)

equil_in_file = "equil.TEMPLATE.in"
run_in_file = "run.TEMPLATE.in"

# loop through IB iterations
for ib_iter in range(n_start, n_iter):
	traj_file = "run.iter" + str(ib_iter) + ".out.dcd"
	
	# run CG simulation
	RunCGSim(ib_iter, equil_in_file, run_in_file)

	# compute CG probability distributions
	AnalyzeCGTraj(top_file, traj_file,cg_bond_potentials)

	# update potentials
	UpdatePotentials(ib_iter, atom_bond_potentials, prev_cg_bond_potentials, cg_bond_potentials)
    
        # update previous potentials
        prev_cg_bond_potentials[:,:] = cg_bond_potentials[:,:]