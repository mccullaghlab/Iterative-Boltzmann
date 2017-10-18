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
from scipy.interpolate import interp1d

kB = 0.001987 # kcal/mol/K
thresh = 1E-3 # this will make it pretty smooth
ang_thresh = 1E-3 # this will make it pretty smooth
dih_thresh = 1E-3 # this will make it pretty smooth
degrees_to_radians = 3.1415926535/180.0
# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global top_file, sim_in_file, atom_data_root, ib_lambda, n_iter, param_out_file, kT, n_start
	f = open(cfg_file)
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
			elif option.lower()=='temperature':
				T = float(value)
			else :
				print "Option:", option, " is not recognized"
	f.close()
	kT = kB*T
        n_iter = n_start + n_iter
# MMMMMMMM

def ParsePsfFile(psf_file):
	global n_uniq_atom_types, atom_types, uniq_bond_atom_types, n_uniq_bonds, bonds, n_bonds, uniq_angle_atom_types, n_uniq_angles, angles, n_angles,  uniq_dihedral_atom_types, n_uniq_dihedrals, dihedrals, n_dihedrals
	f = open(psf_file)
	atom_flag = bond_flag = angle_flag = dihedral_flag = "false"
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
		elif angle_flag == "true" and angle < n_angles:
			for line_angle in range(3):
				angles.append([])
				angles[angle].append(line[line_angle*24:line_angle*24+8])
				angles[angle].append(line[line_angle*24+8:line_angle*24+16])
				angles[angle].append(line[line_angle*24+16:line_angle*24+24])
				same = "false"
				if angle > 0:
					for angle2 in range(angle):
						if (atom_types[int(angles[angle][0])-1] == atom_types[int(angles[angle2][0])-1] and atom_types[int(angles[angle][1])-1] == atom_types[int(angles[angle2][1])-1] and atom_types[int(angles[angle][2])-1] == atom_types[int(angles[angle2][2])-1]) or (atom_types[int(angles[angle][0])-1] == atom_types[int(angles[angle2][2])-1] and atom_types[int(angles[angle][1])-1] == atom_types[int(angles[angle2][1])-1] and atom_types[int(angles[angle][2])-1] == atom_types[int(angles[angle2][0])-1]):
							same = "true"
							uniq_angle_num = angles[angle2][3]
							break
				if same == "false":
					uniq_angle_atom_types.append([])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][0])-1])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][1])-1])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][2])-1])
					uniq_angle_num = str(n_uniq_angles)
					n_uniq_angles += 1
				angles[angle].append(uniq_angle_num)
				angle += 1
				if angle == n_angles:
					break
		elif dihedral_flag == "true" and dihedral < n_dihedrals:
			for line_dihedral in range(2):
				dihedrals.append([])
				dihedrals[dihedral].append(line[line_dihedral*32:line_dihedral*32+8])
				dihedrals[dihedral].append(line[line_dihedral*32+8:line_dihedral*32+16])
				dihedrals[dihedral].append(line[line_dihedral*32+16:line_dihedral*32+24])
				dihedrals[dihedral].append(line[line_dihedral*32+24:line_dihedral*32+32])
				same = "false"
				if dihedral > 0:
					for dihedral2 in range(dihedral):
						if (atom_types[int(dihedrals[dihedral][0])-1] == atom_types[int(dihedrals[dihedral2][0])-1] and atom_types[int(dihedrals[dihedral][1])-1] == atom_types[int(dihedrals[dihedral2][1])-1] and atom_types[int(dihedrals[dihedral][2])-1] == atom_types[int(dihedrals[dihedral2][2])-1] and atom_types[int(dihedrals[dihedral][3])-1] == atom_types[int(dihedrals[dihedral2][3])-1]) or (atom_types[int(dihedrals[dihedral][0])-1] == atom_types[int(dihedrals[dihedral2][3])-1] and atom_types[int(dihedrals[dihedral][1])-1] == atom_types[int(dihedrals[dihedral2][2])-1] and atom_types[int(dihedrals[dihedral][2])-1] == atom_types[int(dihedrals[dihedral2][1])-1] and atom_types[int(dihedrals[dihedral][3])-1] == atom_types[int(dihedrals[dihedral2][0])-1]):
							same = "true"
							uniq_dihedral_num = dihedrals[dihedral2][4]
							break
				if same == "false":
					uniq_dihedral_atom_types.append([])
					uniq_dihedral_atom_types[n_uniq_dihedrals].append(atom_types[int(dihedrals[dihedral][0])-1])
					uniq_dihedral_atom_types[n_uniq_dihedrals].append(atom_types[int(dihedrals[dihedral][1])-1])
					uniq_dihedral_atom_types[n_uniq_dihedrals].append(atom_types[int(dihedrals[dihedral][2])-1])
					uniq_dihedral_atom_types[n_uniq_dihedrals].append(atom_types[int(dihedrals[dihedral][3])-1])
					uniq_dihedral_num = str(n_uniq_dihedrals)
					n_uniq_dihedrals += 1
				dihedrals[dihedral].append(uniq_dihedral_num)
				dihedral += 1
				if dihedral == n_dihedrals:
					break

		if line[9:15] == "!NATOM":
			n_atoms = int(line[0:9])
			atom_flag = "true"
		elif line[9:15] == "!NBOND":
			n_bonds = int(line[0:9])
			bond_flag = "true"
		elif line[9:16] == "!NTHETA":
			n_angles = int(line[0:9])
			angle_flag = "true"
		elif line[9:14] == "!NPHI":
			n_dihedrals = int(line[0:9])
			dihedral_flag = "true"

	f.close()
	bonds = np.asmatrix(bonds,dtype=int)
	angles = np.asmatrix(angles,dtype=int)
	dihedrals = np.asmatrix(dihedrals,dtype=int)

def ComputeDihedralHists(atom_positions, dihedrals, dihedral_hists, dihedral_min, dihedral_delta):
	# get size info
	n_dihedrals = dihedrals.shape[0]
	n_bins = dihedral_hists.shape[1]
	for dihedral in range(n_dihedrals):

		dih = ComputeDih(atom_positions[dihedrals[dihedral,0]-1,:],atom_positions[dihedrals[dihedral,1]-1,:],atom_positions[dihedrals[dihedral,2]-1,:],atom_positions[dihedrals[dihedral,3]-1,:])
		dih_bin = int((dih-dihedral_min)/dihedral_delta)
		if dih_bin >= 0 and dih_bin < n_bins:
			dihedral_hists[dihedrals[dihedral,4],dih_bin] += 1

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
        global atom_data_root
	n_bonds = bond_hists.shape[0]
	n_bins = bond_hists.shape[1]
        x_mat = np.empty(n_bins,dtype=float)
        # these are the histogram details for the input files
        input_bond_min = 0.0
        input_bond_max = 24.0
        input_bond_delta = 0.25
        input_bond_bins  = int((input_bond_max-input_bond_min)/input_bond_delta)
        input_bond_hists = np.empty((n_bonds, input_bond_bins),dtype=float)
        input_x_mat = np.empty(input_bond_bins,dtype=float)
        #
        for i in range(n_bins):
            x_mat[i] = bond_min+(i+0.5)*bond_delta
        for i in range(input_bond_bins):
            input_x_mat[i] = input_bond_min+(i+0.5)*input_bond_delta
    
	for bond in range(n_bonds):
		filename = atom_data_root+"/bonds/"+uniq_bond_atom_types[bond][0].strip()+"_"+uniq_bond_atom_types[bond][1].strip()+".hist"
		filename2 = atom_data_root+"/bonds/"+uniq_bond_atom_types[bond][1].strip()+"_"+uniq_bond_atom_types[bond][0].strip()+".hist"
		if os.path.isfile(filename):
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				input_bond_hists[bond,count] += float(val)
				count += 1
			inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				input_bond_hists[bond,count] += float(val)
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
			filename = atom_data_root+"/bonds/"+type1+"_"+type2+".hist"
			filename2 = atom_data_root+"/bonds/"+type2+"_"+type1+".hist"
			if os.path.isfile(filename):
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					input_bond_hists[bond,count] = float(val)
					count += 1
				inp.close()
			elif os.path.isfile(filename2):
				filename = filename2
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					input_bond_hists[bond,count] = float(val)
					count += 1
				inp.close()
			else:
				print "Did not find prameters for:", uniq_bond_atom_types[bond]
				sys.exit()
                # interpolate so that we can make grid spacing arbitrary
                tck = interpolate.splrep(input_x_mat,input_bond_hists[bond,:])
                bond_hists[bond,:] = interpolate.splev(x_mat,tck,der=0)

# average angle distance histograms (convert them to probability densities)
def ReadAngleHists(angle_hists, angle_min, angle_delta, uniq_angle_atom_types):
        global atom_data_root
	n_angles = angle_hists.shape[0]
	n_bins = angle_hists.shape[1]
        x_mat = np.empty(n_bins,dtype=float)
        # these are the histogram details for the input files
        input_angle_min = 0.0
        input_angle_max = 180.0
        input_angle_delta = 2.0
        input_angle_bins  = int((input_angle_max-input_angle_min)/input_angle_delta)+1
        input_angle_hists = np.empty((n_angles, input_angle_bins),dtype=float)
        input_x_mat = np.empty(input_angle_bins,dtype=float)
        #
        for i in range(n_bins):
            x_mat[i] = angle_min+(i+0.5)*angle_delta
        for i in range(input_angle_bins):
            input_x_mat[i] = input_angle_min+(i+0.5)*input_angle_delta
    

	for angle in range(n_angles):
		filename = atom_data_root+"/angs/"+uniq_angle_atom_types[angle][0].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][2].strip()+".hist"
		filename2 = atom_data_root+"/angs/"+uniq_angle_atom_types[angle][2].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][0].strip()+".hist"
		if os.path.isfile(filename):
		    inp = open(filename,'r')
		    count = 0
		    for line in inp:
		        junk, val = line.split()
			input_angle_hists[angle,count] = float(val)
			count += 1
		    inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				input_angle_hists[angle,count] = float(val)
				count += 1
			inp.close()
		else :
			if uniq_angle_atom_types[angle][0][1] == "1":
				type1 = uniq_angle_atom_types[angle][0][0] + str(2)
			else:
				type1 = uniq_angle_atom_types[angle][0][0] + str(1)
			if uniq_angle_atom_types[angle][1][1] == "1":
				type2 = uniq_angle_atom_types[angle][1][0] + str(2)
			else:
				type2 = uniq_angle_atom_types[angle][1][0] + str(1)
			if uniq_angle_atom_types[angle][2][1] == "1":
				type3 = uniq_angle_atom_types[angle][2][0] + str(2)
			else:
				type3 = uniq_angle_atom_types[angle][2][0] + str(1)
			filename = atom_data_root+"/angs/"+type1+"_"+type2+"_"+type3+".hist"
			filename2 = atom_data_root+"/angs/"+type3+"_"+type2+"_"+type1+".hist"
			if os.path.isfile(filename):
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					input_angle_hists[angle,count] = float(val)
					count += 1
				inp.close()
			elif os.path.isfile(filename2):
				filename = filename2
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					input_angle_hists[angle,count] = float(val)
					count += 1
				inp.close()
			else:
				print "Did not find data for angle", uniq_angle_atom_types[angle]
				sys.exit()
                # interpolate so that we can make grid spacing arbitrary
                tck = interpolate.splrep(input_x_mat,input_angle_hists[angle,:])
                angle_hists[angle,:] = interpolate.splev(x_mat,tck,der=0)

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
def AnalyzeCGTraj(top_file, traj_file, cg_bond_potentials, cg_bond_start_stop, cg_angle_potential, cg_angle_start_stop, cg_dihedral_potential, cg_dihedral_start_stop):
	global bonds, bond_min, bond_delta, angles, angle_min, angle_delta, dihedrals, dihedral_min, dihedral_delta, ib_iter
        # initialize bond arrays
	n_bonds = cg_bond_potentials.shape[0]
	n_bond_bins = cg_bond_potentials.shape[1]
        cg_bond_hists = np.zeros((n_bonds,n_bond_bins),dtype=float)
        cg_bond_forces = np.zeros((n_bonds,n_bond_bins),dtype=float)
        # initialize angle arrays
	n_angles = cg_angle_potentials.shape[0]
	n_angle_bins = cg_angle_potentials.shape[1]
        cg_angle_hists = np.zeros((n_angles,n_angle_bins),dtype=float)
        cg_angle_forces = np.zeros((n_angles,n_angle_bins),dtype=float)
        # initialize dihedral arrays
	n_dihedrals = cg_dihedral_potentials.shape[0]
	n_dihedral_bins = cg_dihedral_potentials.shape[1]
        cg_dihedral_hists = np.zeros((n_dihedrals,n_dihedral_bins),dtype=float)
        cg_dihedral_forces = np.zeros((n_dihedrals,n_dihedral_bins),dtype=float)
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
		# add to angle histograms
		ComputeAngleHists(positions, angles, cg_angle_hists, angle_min, angle_delta)
		# add to dihedral histograms
		ComputeDihedralHists(positions, dihedrals, cg_dihedral_hists, dihedral_min, dihedral_delta)
        filename = "bond" + str(ib_iter) + "_out.hist"
        WriteBondHists(cg_bond_hists, bond_min, bond_delta, filename)
        filename = "angle" + str(ib_iter) + "_out.hist"
        WriteHists(cg_angle_hists,angle_min, angle_delta, filename)
        filename = "dihedral" + str(ib_iter) + "_out.hist"
        WriteHists(cg_dihedral_hists,dihedral_min, dihedral_delta, filename)
        CreateBondPotentials(cg_bond_hists, cg_bond_start_stop, cg_bond_potentials, cg_bond_forces, bond_min, bond_delta, ib_iter)
        CreateAnglePotentials(cg_angle_hists, cg_angle_start_stop, cg_angle_potentials, cg_angle_forces, angle_min, angle_delta, ib_iter)
        CreateDihedralPotentials(cg_dihedral_hists, cg_dihedral_start_stop, cg_dihedral_potentials, cg_dihedral_forces, dihedral_min, dihedral_delta, ib_iter)

def ReadAtomHists(atom_bond_potentials, atom_bond_start_stop, atom_angle_potentials, atom_angle_start_stop, atom_dihedral_potentials, atom_dihedral_start_stop):
	global n_uniq_bonds, n_bond_bins, bond_min, bond_delta, uniq_bond_atom_types, n_uniq_angles, n_angle_bins, angle_min, angle_delta, uniq_angle_atom_types, n_uniq_dihedrals, n_dihedral_bins, dihedral_min, dihedral_delta, uniq_dihedral_atom_types, param_out_file
        # allocate hist and force arrays
        atom_bond_hists = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float) 
        atom_bond_forces = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float) 
        atom_angle_hists = np.zeros((n_uniq_angles,n_angle_bins),dtype=float) 
        atom_angle_forces = np.zeros((n_uniq_angles,n_angle_bins),dtype=float) 
        atom_dihedral_hists = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float) 
        atom_dihedral_forces = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float) 
        # read hists
	ReadBondHists(atom_bond_hists, bond_min, bond_delta, uniq_bond_atom_types)
	ReadAngleHists(atom_angle_hists, angle_min, angle_delta, uniq_angle_atom_types)
	ReadDihedralHists(atom_dihedral_hists, dihedral_min, dihedral_delta, uniq_dihedral_atom_types)
        # write param file
	param_out = open(param_out_file, 'w')
        for bond in range(n_uniq_bonds):
		param_out.write("bond_coeff %2d %s BOND_%s\n" %(bond+1, "bonds.ib", str(bond+1).strip()))
        for angle in range(n_uniq_angles):
		param_out.write("angle_coeff %2d %s ANG_%s\n" %(angle+1, "angles.ib", str(angle+1).strip()))
        for dihedral in range(n_uniq_dihedrals):
		param_out.write("dihedral_coeff %2d %s DIH_%s\n" %(dihedral+1, "dihedrals.ib", str(dihedral+1).strip()))
        param_out.write("pair_coeff  *  * 0.1  10.4  10.4\n")
        param_out.close()
        # create potentials from hists
        CreateBondPotentials(atom_bond_hists, atom_bond_start_stop, atom_bond_potentials, atom_bond_forces, bond_min, bond_delta,0)
        WriteBondPotentials(atom_bond_potentials, atom_bond_forces, atom_bond_start_stop, bond_min, bond_delta,"bond0.ib","bond0.ener","bond0.start_stop")
        WriteBondHists(atom_bond_hists, bond_min, bond_delta, "bond0_in.hist")
	CreateAnglePotentials(atom_angle_hists, atom_angle_start_stop, atom_angle_potentials, atom_angle_forces, angle_min, angle_delta, 0)
        WriteHists(atom_angle_hists,angle_min, angle_delta,"angle0_in.hist")
        WriteAnglePotentials(atom_angle_potentials, atom_angle_forces, atom_angle_start_stop, angle_min, angle_delta,"angle0.ib","angle0.ener","angle0.start_stop", "ANG")
	CreateDihedralPotentials(atom_dihedral_hists, atom_dihedral_start_stop, atom_dihedral_potentials, atom_dihedral_forces, dihedral_min, dihedral_delta, 0)
        WriteAnglePotentials(atom_dihedral_potentials, atom_dihedral_forces, atom_dihedral_start_stop, dihedral_min, dihedral_delta,"dihedral0.ib","dihedral0.ener","dihedral0.start_stop", "DIH")
        WriteHists(atom_dihedral_hists,dihedral_min, dihedral_delta,"dihedral0_in.hist")


# average bond distance histograms (convert them to probability densities)
def CreateBondPotentials(bond_hists, bond_start_stop, bond_potentials, bond_forces, bond_min, bond_delta, ib_iter):
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
			if bond_hists[bond,i]/x2 > thresh :  # 1E-3 seems to work pretty well here
				bond_potentials[bond,i] = -kT*math.log(bond_hists[bond,i]/x2)

                # find start
                for i in range (n_bins):
                        if bond_potentials[bond,i] > thresh:
				start = i
                                break
                # find stop
                for i in range(n_bins-1,-1,-1):
			if bond_potentials[bond,i] > thresh:
				stop = i + 1 # because python is non-inclusive of upperbound
                                break
                bond_start_stop[bond,0] = start
                bond_start_stop[bond,1] = stop
		# find the minimum energy
		min_val = np.amin(bond_potentials[bond,start:stop])

		# now smooth the potential using a cubic spline fit	
                tck = interpolate.splrep(coeff_mat[start:stop,1],bond_potentials[bond,start:stop]-min_val)
                bond_potentials[bond,start:stop] = interpolate.splev(coeff_mat[start:stop,1], tck, der=0)

		# fit previous section if we have enough data
		if (stop-start) > 4:
		    # fit function to parabola to extrapolate to short and long distances
		    k, rss, rank, sv = np.linalg.lstsq(coeff_mat[start:start+4],bond_potentials[bond,start:start+4])
                    # check to make sure the potential is increasing
                    if k[2] > 0:
                        # place fit data into potential
	    	        for j in range(start):
		            bond_potentials[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                    # keep increasing size of fit region until k[2] > 0
                    else:
                        print "k[2] < 0 for decreasing r of bond ", bond+1, " at iteration", ib_iter
                        fit_size = 5
                        while k[2] < 0 and fit_size <= stop-start+1:
		            k, rss, rank, sv = np.linalg.lstsq(coeff_mat[start:start+fit_size],bond_potentials[bond,start:start+fit_size])
                            fit_size += 1
                        if fit_size > stop-start+1:
                            print "ERROR: Cannot fit data to increasing function... BOMBING OUT"
                            sys.exit()
                        print "k[2] > 0 for decreasing r of bond ", bond+1, " for fit_size=", fit_size
                        # place fit data into potential
	    	        for j in range(start):
		            bond_potentials[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
		    # fit function to parabola to extrapolate to short and long distances
		    k, rss, rank, sv = np.linalg.lstsq(coeff_mat[stop-4:stop],bond_potentials[bond,stop-4:stop])
                    # check to make sure the potential is increasing
                    if k[2] > 0:
                        # place fit data into potential
		        for j in range(stop,n_bins):
			    bond_potentials[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                    # keep increasing fit region until ...
                    else:
                        print "k[2] < 0 for increasing r of bond ", bond+1, " at iteration", ib_iter
                        fit_size = 5
                        while k[2] < 0 and fit_size <= stop-start+1:
		            k, rss, rank, sv = np.linalg.lstsq(coeff_mat[stop-fit_size:stop],bond_potentials[bond,stop-fit_size:stop])
                            fit_size += 1
                        if fit_size > stop-start+1:
                            print "ERROR: Cannot fit data to increasing function... BOMBING OUT"
                        print "k[2] > 0 for increasing r of bond ", bond+1, " for fit_size=", fit_size
                        # place fit data into potential
		        for j in range(stop,n_bins):
			    bond_potentials[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                else:
                    print "ERROR: Not enough data to fit for bond=", bond, " at iteration ", ib_iter
                    print "Bombing out..."
                    sys.exit()

		# now smooth the total  potential using a cubic spline fit	
		tck = interpolate.splrep(coeff_mat[:,1],bond_potentials[bond,:])
                bond_potentials[bond,:] = interpolate.splev(coeff_mat[:,1], tck, der=0)
                bond_forces[bond,:] = -interpolate.splev(coeff_mat[:,1], tck, der=1)
               
# write bond potentials and forces to files
def WriteBondPotentials(bond_potentials, bond_forces, bond_start_stop, bond_min, bond_delta, ib_out,ener_out, start_stop_out):

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
        # write potential start stop
	start_stop_file = open(start_stop_out, 'w')
        for bond in range(n_bonds):
            start_stop_file.write("%20d%20d\n" % (bond_start_stop[bond,0], bond_start_stop[bond,1]))
        start_stop_file.close()

def WriteHists(hists, x_min, x_delta, hist_out):

	n_hists = hists.shape[0]
	n_bins = hists.shape[1]

        # normalize histograms 
        norm = np.empty(n_bonds,dtype=float)
        for hist in range(n_hists):
            norm[hist] = np.sum(hists[hist,:])*x_delta
        # write historgram file
	hist_file = open(hist_out, 'w')
	for i in range(n_bins):
		hist_file.write("%20.5f" % (x_min+(i+0.5)*x_delta))
		for hist in range(n_hists):
			hist_file.write("%20.7f" % (hists[hist,i]/norm[hist]))
		hist_file.write("\n")
	hist_file.close()

def WriteAnglePotentials(potentials, forces, start_stop, x_min, x_delta, ib_out, ener_out, start_stop_out, char):
	n_hists = potentials.shape[0]
	n_bins = potentials.shape[1]
        # write LAMMPS potential file
	ib_file = open(ib_out,'w')
	for hist in range(n_hists):
		ib_file.write("\n%s_%s\n" % (char,str(hist+1)))
		ib_file.write("N %5d\n\n" % (n_bins))
		for i in range(n_bins):
			ib_file.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,x_min+(i)*x_delta, potentials[hist,i],forces[hist,i])) # angles must be from 0 to 180 and dihedrals -180 to 180
        ib_file.close()
        # write file containing all energies in easily graphable format
	ener_file = open(ener_out, 'w')
	for i in range(n_bins):
		ener_file.write("%20.5f" % (x_min+(i+0.5)*x_delta))
		for hist in range(n_hists):
			ener_file.write("%20.7f" % (potentials[hist,i]))
		ener_file.write("\n")
	ener_file.close()
        # write potential start stop
	start_stop_file = open(start_stop_out, 'w')
        for hist in range(n_hists):
            start_stop_file.write("%20d%20d\n" % (start_stop[hist,0], start_stop[hist,1]))
        start_stop_file.close()

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


# MM_ANGLES

def CreateAnglePotentials(angle_hists, angle_start_stop, angle_potentials, angle_forces, angle_min, angle_delta, ib_iter):
	global kT, degrees_to_radians

	n_angles = angle_hists.shape[0]
	n_bins = angle_hists.shape[1]
        # declare and fill coefficient matrices
	coeff_mat = np.empty((n_bins,3),dtype=float)
        init_x = np.empty((3,3),dtype=float)
	x_mat = np.empty((n_bins),dtype=float)
	sinx = np.empty((n_bins),dtype=float)
	for i in range(n_bins):
		x = angle_min+(i+0.5)*angle_delta
		coeff_mat[i,0] = 1.0
		coeff_mat[i,1] = x
		x_mat[i] = x
		coeff_mat[i,2] = x*x
                sinx[i] = math.sin(x*degrees_to_radians)
	# first check to see if we have two copies of same angle
	for angle in range(n_angles):

		# convert to probability density
		angle_hists[angle,:] /= (np.sum(angle_hists[angle,:])*angle_delta)
		# convert to energies
		min_val = 0
		count = 0
		for i in range(n_bins):
			if angle_hists[angle,i]/sinx[i] > ang_thresh:
				angle_potentials[angle,i] = -kT*math.log(angle_hists[angle,i]/sinx[i])
			else:
				angle_potentials[angle,i] = 0.0

                # find start
                for i in range (n_bins):
                        if angle_potentials[angle,i] > ang_thresh:
				start = i
                                break
                        else:
                            angle_potentials[angle,i] = 0.0
                # find stop
                for i in range(n_bins-1,-1,-1):
			if angle_potentials[angle,i] > ang_thresh:
				stop = i + 1 # because python is non-inclusive of upperbound
                                break
                        else:
                            angle_potentials[angle,i] = 0.0
                # smooth the potential and compute forces
                indice, = angle_potentials[angle,:].nonzero()
                tck = interpolate.splrep(x_mat[indice],angle_potentials[angle,indice])
                angle_potentials[angle,start:stop] = interpolate.splev(x_mat[start:stop], tck, der=0)
                angle_forces[angle,start:stop] = -interpolate.splev(x_mat[start:stop], tck, der=1)
                # trim the fat
                if start > 0 or stop < n_bins:
                    # make sure start is sloped appropriately
                    while angle_forces[angle,start] < 0 or angle_forces[angle,start+1] < 0 or angle_forces[angle,start+2] < 0:
                        start += 1
                    # make sure end is sloped appropriately
                    while angle_forces[angle,stop-1] > 0 or angle_forces[angle,stop-2] > 0 or angle_forces[angle,stop-3] > 0:
                        stop -= 1
                # save start stop
                angle_start_stop[angle,0] = start
                angle_start_stop[angle,1] = stop

		# fit unsampled regions
		hole_flag = "false"
                fit_hole = "false"
                for i in range(start,stop):

			if angle_potentials[angle,i] < ang_thresh and hole_flag == "false":
				hole_start = i
				hole_flag = "true"
                        elif angle_potentials[angle,i] > ang_thresh and hole_flag == "true":
				hole_stop = i
                                fit_hole = "true"
                                hole_flag = "false"
			
			# fit hole if identified
                        if fit_hole == "true":
                                if (hole_stop - hole_start) > 3:
                                    # make sure legs are sloped appropriately
                                    while hole_start >= 0 and (angle_forces[angle,hole_start] > 0 or angle_forces[angle,hole_start-1] > 0 or angle_forces[angle,hole_start-2] > 0):
                                        hole_start -= 1
                                    while hole_stop < n_bins and (angle_forces[angle,hole_stop] < 0 or angle_forces[angle,hole_stop+1] < 0 or angle_forces[angle,hole_stop+2] < 0):
                                        hole_stop += 1
                                    if hole_stop == n_bins or hole_start < 0: # have to bridge
                                        break
                                    data_y = np.append(angle_potentials[angle,hole_start-3:hole_start],angle_potentials[angle,hole_stop:hole_stop+3])
                                    data_x = np.append(coeff_mat[hole_start-3:hole_start], coeff_mat[hole_stop:hole_stop+3], axis=0)
				    k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
				    # fill in hole
                                    for j in range(hole_start,hole_stop):
					angle_potentials[angle,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                                else:
				    # smooth region with spline interpolation
                                    tck = interpolate.splrep(x_mat[hole_start-3:hole_stop+3],angle_potentials[angle,hole_start-3:hole_stop+3])
                                    angle_potentials[angle,hole_start-3:hole_stop+3] = interpolate.splev(x_mat[hole_start-3:hole_stop+3], tck, der=0)
                                # hole is filled
				fit_hole = "false"

		# connect stop and start using periodicity of angle
		# fit remaining missing data
		for j in range(init_x.shape[0]):
                        init_x[j,0] = 1.0
			init_x[j,1] = coeff_mat[start+j,1] + 180.0
			init_x[j,2] = init_x[j,1]*init_x[j,1]
                data_y = np.append(angle_potentials[angle,stop-3:stop], angle_potentials[angle,start:start+3])
                data_x = np.append(coeff_mat[stop-3:stop], init_x[:], axis=0)
		k, rss, rank, s = np.linalg.lstsq(data_x, data_y)

		for j in range(stop,n_bins):
			angle_potentials[angle,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
		for j in range(start):
			x = angle_min+(j+0.5)*angle_delta+180.0
			x2 = x*x
			angle_potentials[angle,j] = k[0] + k[1]*x + k[2]*x2

		# now smooth the potential using spline
                min_val = np.amin(angle_potentials[angle,:])
                tck = interpolate.splrep(x_mat,angle_potentials[angle,:]-min_val)
                angle_potentials[angle,:] = interpolate.splev(x_mat, tck, der=0)
                angle_forces[angle,:] = -interpolate.splev(x_mat, tck, der=1)

def UpdateAngles(ib_iter, atom_angle_potentials, atom_angle_start_stop, prev_cg_angle_potentials, prev_cg_angle_start_stop, cg_angle_potentials, cg_angle_start_stop):
	global kT, ib_lambda, angle_delta, angle_min
        
	n_angles = cg_angle_potentials.shape[0]
	n_bins = cg_angle_potentials.shape[1]

	cg_angle_forces = np.zeros((n_angles,n_bins),dtype=float)
        cg_angle_potential_temp = np.zeros(n_bins,dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	ib_out = "angle" + str(ib_iter+1) + ".ib"
	ener_out = "angle" + str(ib_iter+1) + ".ener"
	start_stop_out = "angle" + str(ib_iter+1) + ".start_stop"
        # make coefficient matrix
	for i in range(n_bins):
		x_mat[i] = angle_min+(i+0.5)*angle_delta

        for angle in range(n_angles):
            # get lowest min and largest max
            current_min = min(atom_angle_start_stop[angle,0],cg_angle_start_stop[angle,0],prev_cg_angle_start_stop[angle,0])
            current_max = max(atom_angle_start_stop[angle,1],cg_angle_start_stop[angle,1],prev_cg_angle_start_stop[angle,1])
            # perform the iterative Boltzmann procedure
            for i in range(current_min,current_max):
	        cg_angle_potential_temp[i] = prev_cg_angle_potentials[angle,i] + ib_lambda * (atom_angle_potentials[angle,i] - cg_angle_potentials[angle,i])
            # now add asymptotic behavior and smooth 
            AngleAsymptoticBehavior(cg_angle_potential_temp,cg_angle_forces[angle,:],x_mat,current_min,current_max)
            cg_angle_potentials[angle,:] = cg_angle_potential_temp[:] 
            cg_angle_start_stop[angle,0] = current_min
            cg_angle_start_stop[angle,1] = current_max
        # now to write them to file
        WriteAnglePotentials(cg_angle_potentials,cg_angle_forces,cg_angle_start_stop, angle_min, angle_delta,ib_out,ener_out,start_stop_out,"ANG")

def AngleAsymptoticBehavior(angle_potential,angle_force,x_mat,start,stop):

	n_bins = angle_potential.shape[0]
        # declare and fill coefficient matrices
	coeff_mat = np.empty((n_bins,3),dtype=float)
        init_x = np.empty((3,3),dtype=float)
	for i in range(n_bins):
		coeff_mat[i,0] = 1.0
		coeff_mat[i,1] = x_mat[i]
		coeff_mat[i,2] = x_mat[i]*x_mat[i]

        # trim the fat
        if start > 0 or stop < n_bins:
            # make sure start is sloped appropriately
            while angle_force[start] < 0 or angle_force[start+1] < 0 or angle_force[start+2] < 0:
                angle_potential[start] = 0.0
                start += 1
            # make sure end is sloped appropriately
            while angle_force[stop-1] > 0 or angle_force[stop-2] > 0 or angle_force[stop-3] > 0:
                angle_potential[stop] = 0.0
                stop -= 1
        # smooth the potential and compute forces
        indice, = angle_potential.nonzero()
        tck = interpolate.splrep(x_mat[indice],angle_potential[indice])
        angle_potential[start:stop] = interpolate.splev(x_mat[start:stop], tck, der=0)
        angle_force[start:stop] = -interpolate.splev(x_mat[start:stop], tck, der=1)
	# fit remaining missing data
	for j in range(init_x.shape[0]):
                init_x[j,0] = 1.0
		init_x[j,1] = coeff_mat[start+j,1] + 180.0
		init_x[j,2] = init_x[j,1]*init_x[j,1]
        data_y = np.append(angle_potential[stop-3:stop], angle_potential[start:start+3])
        data_x = np.append(coeff_mat[stop-3:stop], init_x[:], axis=0)
	k, rss, rank, s = np.linalg.lstsq(data_x, data_y)

	for j in range(stop,n_bins):
		angle_potential[j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
	for j in range(start):
		x = angle_min+(j+0.5)*angle_delta+180.0
		x2 = x*x
		angle_potential[j] = k[0] + k[1]*x + k[2]*x2

	# now smooth the potential using spline
        min_val = np.amin(angle_potential)
        tck = interpolate.splrep(x_mat,angle_potential-min_val)
        angle_potential[:] = interpolate.splev(x_mat, tck, der=0)
        angle_force[:] = -interpolate.splev(x_mat, tck, der=1)

# add to angle distance histograms 
def ComputeAngleHists(atom_positions, angles, angle_hists, angle_min, delta_angle):
	# get size info
	n_angles = angles.shape[0]
	n_bins = angle_hists.shape[1]
	for angle in range(n_angles):

		ang = ComputeAng(atom_positions[angles[angle,0]-1,:],atom_positions[angles[angle,1]-1,:],atom_positions[angles[angle,2]-1,:])
		ang_bin = int((ang-angle_min)/delta_angle)
		if ang_bin >= 0 and ang_bin < n_bins:
			angle_hists[angles[angle,3],ang_bin] += 1
# MMMMMMMM

# Dihedrals

# average dihedral distance histograms (convert them to probability densities)
def ReadDihedralHists(dihedral_hists, dihedral_min, delta_dihedral, uniq_dihedral_atom_types):
        global kT, atom_data_root

	n_dihedrals = dihedral_hists.shape[0]
	n_bins = dihedral_hists.shape[1]
        x_mat = np.empty(n_bins,dtype=float)

        input_dihedral_min = -179.0 
        input_dihedral_max = 180.0 # inclusive
        input_dihedral_delta = 1.0
        input_dihedral_bins  = int((input_dihedral_max-input_dihedral_min)/input_dihedral_delta) + 1
        input_dihedral_hists = np.zeros((n_dihedrals,input_dihedral_bins),dtype=float)
        input_x_mat = np.empty(input_dihedral_bins,dtype=float)
        for i in range(n_bins):
            x_mat[i] = dihedral_min+(i+0.5)*dihedral_delta
        for i in range(input_dihedral_bins):
            input_x_mat[i] = input_dihedral_min+(i+0.5)*input_dihedral_delta


	for dihedral in range(n_dihedrals):
		filename = atom_data_root+"/dihs/"+uniq_dihedral_atom_types[dihedral][0].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][3].strip()+".hist"
		filename2 = atom_data_root+"/dihs/"+uniq_dihedral_atom_types[dihedral][3].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][0].strip()+".hist"
		if os.path.isfile(filename):
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				input_dihedral_hists[dihedral,count] += float(val)
				count += 1
			inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				input_dihedral_hists[dihedral,count] += float(val)
				count += 1
			inp.close()
		else:
			if uniq_dihedral_atom_types[dihedral][0][1] == "1":
				type1 = uniq_dihedral_atom_types[dihedral][0][0] + str(2)
			else:
				type1 = uniq_dihedral_atom_types[dihedral][0][0] + str(1)
			if uniq_dihedral_atom_types[dihedral][1][1] == "1":
				type2 = uniq_dihedral_atom_types[dihedral][1][0] + str(2)
			else:
				type2 = uniq_dihedral_atom_types[dihedral][1][0] + str(1)
			if uniq_dihedral_atom_types[dihedral][2][1] == "1":
				type3 = uniq_dihedral_atom_types[dihedral][2][0] + str(2)
			else:
				type3 = uniq_dihedral_atom_types[dihedral][2][0] + str(1)
			if uniq_dihedral_atom_types[dihedral][3][1] == "1":
				type4 = uniq_dihedral_atom_types[dihedral][3][0] + str(2)
			else:
				type4 = uniq_dihedral_atom_types[dihedral][3][0] + str(1)
			filename = atom_data_root+"/dihs/"+type1+"_"+type2+"_"+type3+"_"+type4+".hist"
			filename2 = atom_data_root+"/dihs/"+type4+"_"+type3+"_"+type2+"_"+type1+".hist"
			if os.path.isfile(filename):
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					input_dihedral_hists[dihedral,count] = float(val)
					count += 1
				inp.close()
			elif os.path.isfile(filename2):
				filename = filename2
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					input_dihedral_hists[dihedral,count] = float(val)
					count += 1
				inp.close()
			else:
				print "Did not find data for dihedral", uniq_dihedral_atom_types[dihedral]
				sys.exit()
                # normalize
                input_dihedral_hists[dihedral,:] /= (np.sum(input_dihedral_hists[dihedral,:])*input_dihedral_delta)
                # interpolate so that we can make grid spacing arbitrary
                tck = interpolate.splrep(input_x_mat,input_dihedral_hists[dihedral,:])
                dihedral_hists[dihedral,:] = interpolate.splev(x_mat,tck,der=0)


def CreateDihedralPotentials(dihedral_hists, dihedral_start_stop, dihedral_potentials, dihedral_forces, dihedral_min, dihedral_delta, ib_iter):
	global kT, degrees_to_radians

	n_dihedrals = dihedral_hists.shape[0]
	n_bins = dihedral_hists.shape[1]
        # declare and fill coefficient matrices
	coeff_mat = np.empty((n_bins,3),dtype=float)
        init_x = np.empty((3,3),dtype=float)
	x_mat = np.empty((n_bins),dtype=float)
	for i in range(n_bins):
		x = dihedral_min+(i+0.5)*dihedral_delta
		coeff_mat[i,0] = 1.0
		coeff_mat[i,1] = x
		x_mat[i] = x
		coeff_mat[i,2] = x*x
	# first check to see if we have two copies of same dihedral
	for dihedral in range(n_dihedrals):

		# convert to probability density
		dihedral_hists[dihedral,:] /= (np.sum(dihedral_hists[dihedral,:])*dihedral_delta)
		# convert to energies
		min_val = 0
		count = 0
		for i in range(n_bins):
			if dihedral_hists[dihedral,i] > dih_thresh:
				dihedral_potentials[dihedral,i] = -kT*math.log(dihedral_hists[dihedral,i])
			else:
				dihedral_potentials[dihedral,i] = 0.0

                # find start
                for i in range (n_bins):
                        if dihedral_potentials[dihedral,i] > dih_thresh:
				start = i
                                break
                        else:
                            dihedral_potentials[dihedral,i] = 0.0
                # find stop
                for i in range(n_bins-1,-1,-1):
			if dihedral_potentials[dihedral,i] > dih_thresh:
				stop = i + 1 # because python is non-inclusive of upperbound
                                break
                        else:
                            dihedral_potentials[dihedral,i] = 0.0
                # smooth the potential and compute forces
                indice, = dihedral_potentials[dihedral,:].nonzero()
                tck = interpolate.splrep(x_mat[indice],dihedral_potentials[dihedral,indice],s=1)
                dihedral_potentials[dihedral,start:stop] = interpolate.splev(x_mat[start:stop], tck, der=0)
                dihedral_forces[dihedral,start:stop] = -interpolate.splev(x_mat[start:stop], tck, der=1)
                # trim the fat
                if start > 0 or stop < n_bins:
                    # make sure start is sloped appropriately
                    while dihedral_forces[dihedral,start] < 0 or dihedral_forces[dihedral,start+1] < 0 or dihedral_forces[dihedral,start+2] < 0:
                        start += 1
                    # make sure end is sloped appropriately
                    while dihedral_forces[dihedral,stop-1] > 0 or dihedral_forces[dihedral,stop-2] > 0 or dihedral_forces[dihedral,stop-3] > 0:
                        stop -= 1
                # save start stop
                dihedral_start_stop[dihedral,0] = start
                dihedral_start_stop[dihedral,1] = stop

		# fit unsampled regions
		hole_flag = "false"
                fit_hole = "false"
                for i in range(start,stop):

			if dihedral_potentials[dihedral,i] < dih_thresh and hole_flag == "false":
				hole_start = i
				hole_flag = "true"
                        elif dihedral_potentials[dihedral,i] > dih_thresh and hole_flag == "true":
				hole_stop = i
                                fit_hole = "true"
                                hole_flag = "false"
			
			# fit hole if identified
                        if fit_hole == "true":
                                if (hole_stop - hole_start) > 3:
                                    # make sure legs are sloped appropriately
                                    while hole_start > 2 and (dihedral_forces[dihedral,hole_start] > 0 or dihedral_forces[dihedral,hole_start-1] > 0 or dihedral_forces[dihedral,hole_start-2] > 0):
                                        hole_start -= 1
                                    while hole_stop < n_bins-3 and (dihedral_forces[dihedral,hole_stop] < 0 or dihedral_forces[dihedral,hole_stop+1] < 0 or dihedral_forces[dihedral,hole_stop+2] < 0):
                                        hole_stop += 1
                                    if hole_stop == n_bins-3 or hole_start ==2:
                                        break
                                    data_y = np.append(dihedral_potentials[dihedral,hole_start-3:hole_start],dihedral_potentials[dihedral,hole_stop:hole_stop+3])
                                    data_x = np.append(coeff_mat[hole_start-3:hole_start], coeff_mat[hole_stop:hole_stop+3], axis=0)
				    k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
				    # fill in hole
                                    for j in range(hole_start,hole_stop):
					dihedral_potentials[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                                else:
				    # smooth region with linear interpolation
                                    k, rss, rank, s = np.linalg.lstsq(coeff_mat[hole_start-3:hole_stop+3,0:2],dihedral_potentials[dihedral,hole_start-3:hole_stop+3])
				    # fill in hole
                                    for j in range(hole_start,hole_stop):
					dihedral_potentials[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] 

#                                    tck = interpolate.splrep(x_mat[hole_start-3:hole_stop+3],dihedral_potentials[dihedral,hole_start-3:hole_stop+3])
#                                    dihedral_potentials[dihedral,hole_start-3:hole_stop+3] = interpolate.splev(x_mat[hole_start-3:hole_stop+3], tck, der=0)
                                # hole is filled
				fit_hole = "false"

		# connect stop and start using periodicity of dihedral
		# fit remaining missing data
		for j in range(init_x.shape[0]):
                        init_x[j,0] = 1.0
			init_x[j,1] = coeff_mat[start+j,1] + 360.0
			init_x[j,2] = init_x[j,1]*init_x[j,1]
                data_y = np.append(dihedral_potentials[dihedral,stop-3:stop], dihedral_potentials[dihedral,start:start+3])
                data_x = np.append(coeff_mat[stop-3:stop], init_x[:], axis=0)
		k, rss, rank, s = np.linalg.lstsq(data_x, data_y)

		for j in range(stop,n_bins):
			dihedral_potentials[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
		for j in range(start):
			x = dihedral_min+(j+0.5)*dihedral_delta+360.0
			x2 = x*x
			dihedral_potentials[dihedral,j] = k[0] + k[1]*x + k[2]*x2

		# now smooth the potential using spline
                min_val = np.amin(dihedral_potentials[dihedral,:])
                tck = interpolate.splrep(x_mat,dihedral_potentials[dihedral,:]-min_val)
                dihedral_potentials[dihedral,:] = interpolate.splev(x_mat, tck, der=0)
                dihedral_forces[dihedral,:] = -interpolate.splev(x_mat, tck, der=1)

# add to dihedral distance histograms 
def ComputeDihedralHists(atom_positions, dihedrals, dihedral_hists, dihedral_min, delta_dihedral):
	# get size info
	n_dihedrals = dihedrals.shape[0]
	n_bins = dihedral_hists.shape[1]
	for dihedral in range(n_dihedrals):

		dih = ComputeDih(atom_positions[dihedrals[dihedral,0]-1,:],atom_positions[dihedrals[dihedral,1]-1,:],atom_positions[dihedrals[dihedral,2]-1,:],atom_positions[dihedrals[dihedral,3]-1,:])
		dih_bin = int((dih-dihedral_min)/delta_dihedral)
		if dih_bin >= 0 and dih_bin < n_bins:
			dihedral_hists[dihedrals[dihedral,4],dih_bin] += 1


def UpdateDihedrals(ib_iter, atom_dihedral_potentials, atom_dihedral_start_stop, prev_cg_dihedral_potentials, prev_cg_dihedral_start_stop, cg_dihedral_potentials, cg_dihedral_start_stop):
	global kT, ib_lambda, dihedral_delta, dihedral_min
        
	n_dihedrals = cg_dihedral_potentials.shape[0]
	n_bins = cg_dihedral_potentials.shape[1]

	cg_dihedral_forces = np.zeros((n_dihedrals,n_bins),dtype=float)
        cg_dihedral_potential_temp = np.zeros(n_bins,dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	ib_out = "dihedral" + str(ib_iter+1) + ".ib"
	ener_out = "dihedral" + str(ib_iter+1) + ".ener"
	start_stop_out = "dihedral" + str(ib_iter+1) + ".start_stop"
        # make coefficient matrix
	for i in range(n_bins):
		x_mat[i] = dihedral_min+(i+0.5)*dihedral_delta

        for dihedral in range(n_dihedrals):
            # get lowest min and largest max
            current_min = min(atom_dihedral_start_stop[dihedral,0],cg_dihedral_start_stop[dihedral,0],prev_cg_dihedral_start_stop[dihedral,0])
            current_max = max(atom_dihedral_start_stop[dihedral,1],cg_dihedral_start_stop[dihedral,1],prev_cg_dihedral_start_stop[dihedral,1])
            # perform the iterative Boltzmann procedure
            for i in range(current_min,current_max):
	        cg_dihedral_potential_temp[i] = prev_cg_dihedral_potentials[dihedral,i] + ib_lambda * (atom_dihedral_potentials[dihedral,i] - cg_dihedral_potentials[dihedral,i])
            # now add asymptotic behavior and smooth 
            DihedralAsymptoticBehavior(cg_dihedral_potential_temp,cg_dihedral_forces[dihedral,:],x_mat,current_min,current_max)
            cg_dihedral_potentials[dihedral,:] = cg_dihedral_potential_temp[:] 
            cg_dihedral_start_stop[dihedral,0] = current_min
            cg_dihedral_start_stop[dihedral,1] = current_max
        # now to write them to file
        WriteAnglePotentials(cg_dihedral_potentials,cg_dihedral_forces,cg_dihedral_start_stop, dihedral_min, dihedral_delta,ib_out,ener_out,start_stop_out,"DIH")

def DihedralAsymptoticBehavior(dihedral_potential,dihedral_force,x_mat,start,stop):

	n_bins = dihedral_potential.shape[0]
        # declare and fill coefficient matrices
	coeff_mat = np.empty((n_bins,3),dtype=float)
        init_x = np.empty((3,3),dtype=float)
	for i in range(n_bins):
		coeff_mat[i,0] = 1.0
		coeff_mat[i,1] = x_mat[i]
		coeff_mat[i,2] = x_mat[i]*x_mat[i]

        # trim the fat
        if start > 0 or stop < n_bins:
            # make sure start is sloped appropriately
            while dihedral_force[start] < 0 or dihedral_force[start+1] < 0 or dihedral_force[start+2] < 0:
                dihedral_potential[start] = 0.0
                start += 1
            # make sure end is sloped appropriately
            while dihedral_force[stop-1] > 0 or dihedral_force[stop-2] > 0 or dihedral_force[stop-3] > 0:
                dihedral_potential[stop] = 0.0
                stop -= 1
        # smooth the potential and compute forces
        indice, = dihedral_potential.nonzero()
        tck = interpolate.splrep(x_mat[indice],dihedral_potential[indice])
        dihedral_potential[start:stop] = interpolate.splev(x_mat[start:stop], tck, der=0)
        dihedral_force[start:stop] = -interpolate.splev(x_mat[start:stop], tck, der=1)
	# fit remaining missing data
	for j in range(init_x.shape[0]):
                init_x[j,0] = 1.0
		init_x[j,1] = coeff_mat[start+j,1] + 360.0
		init_x[j,2] = init_x[j,1]*init_x[j,1]
        data_y = np.append(dihedral_potential[stop-3:stop], dihedral_potential[start:start+3])
        data_x = np.append(coeff_mat[stop-3:stop], init_x[:], axis=0)
	k, rss, rank, s = np.linalg.lstsq(data_x, data_y)

	for j in range(stop,n_bins):
		dihedral_potential[j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
	for j in range(start):
		x = dihedral_min+(j+0.5)*dihedral_delta+360.0
		x2 = x*x
		dihedral_potential[j] = k[0] + k[1]*x + k[2]*x2

	# now smooth the potential using spline
        min_val = np.amin(dihedral_potential)
        tck = interpolate.splrep(x_mat,dihedral_potential-min_val)
        dihedral_potential[:] = interpolate.splev(x_mat, tck, der=0)
        dihedral_force[:] = -interpolate.splev(x_mat, tck, der=1)

# End Dihedrals

# BONDS
def UpdateBonds(ib_iter, atom_bond_potentials, atom_bond_start_stop, prev_cg_bond_potentials, prev_cg_bond_start_stop, cg_bond_potentials, cg_bond_start_stop):
	global kT, ib_lambda, bond_delta, bond_min
        
	n_bonds = cg_bond_potentials.shape[0]
	n_bins = cg_bond_potentials.shape[1]

	cg_bond_forces = np.zeros((n_bonds,n_bins),dtype=float)
        cg_bond_potential_temp = np.zeros(n_bins,dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	ib_out = "bond" + str(ib_iter+1) + ".ib"
	ener_out = "bond" + str(ib_iter+1) + ".ener"
	start_stop_out = "bond" + str(ib_iter+1) + ".start_stop"
        # make coefficient matrix
	for i in range(n_bins):
		x_mat[i] = bond_min+(i+0.5)*bond_delta

        for bond in range(n_bonds):
            # get lowest min
            current_min = min(atom_bond_start_stop[bond,0],cg_bond_start_stop[bond,0],prev_cg_bond_start_stop[bond,0])
            current_max = max(atom_bond_start_stop[bond,1],cg_bond_start_stop[bond,1],prev_cg_bond_start_stop[bond,1])
            # perform the iterative Boltzmann procedure
            for i in range(current_min,current_max):
	        cg_bond_potential_temp[i] = prev_cg_bond_potentials[bond,i] + ib_lambda * (atom_bond_potentials[bond,i] - cg_bond_potentials[bond,i])
            # now add asymptotic behavior and smooth 
            CheckAsymptoticBehavior(cg_bond_potential_temp,cg_bond_forces[bond,:],x_mat,current_min,current_max)
            cg_bond_potentials[bond,:] = cg_bond_potential_temp[:] 
            cg_bond_start_stop[bond,0] = current_min
            cg_bond_start_stop[bond,1] = current_max
        # now to write them to file
        WriteBondPotentials(cg_bond_potentials,cg_bond_forces,cg_bond_start_stop, bond_min, bond_delta,ib_out,ener_out,start_stop_out)

# routine to check potential is behaving appropriately in asymptotic limits and correct if not
def CheckAsymptoticBehavior(potential,force,x,start,stop):

	n_bins = potential.shape[0]
	coeff_mat = np.empty((n_bins,3),dtype=float)

        # make coefficient matrix
	for i in range(n_bins):
		coeff_mat[i,0] = 1.0
		coeff_mat[i,1] = x[i]
		coeff_mat[i,2] = x[i]*x[i]

	# fit beginning section if necessary
        if start > 0:
	    # fit function to parabola to extrapolate to short and long distances
            k, rss, rank, sv = np.linalg.lstsq(coeff_mat[start:start+4],potential[start:start+4])
            # check to make sure the potential is increasing
            if k[2] > 0:
                # place fit data into potential
	        for j in range(start):
	            potential[j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
            # keep increasing size of fit region until k[2] > 0
            else:
                fit_size = 5
                while k[2] < 0:
	            k, rss, rank, sv = np.linalg.lstsq(coeff_mat[start:start+fit_size],potential[start:start+fit_size])
                    fit_size += 1
                # place fit data into potential
	        for j in range(start):
		    potential[j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
	# fit end section if necessary
        if stop < n_bins-1:
	    # fit function to parabola to extrapolate to short and long distances
	    k, rss, rank, sv = np.linalg.lstsq(coeff_mat[stop-4:stop],potential[stop-4:stop])
            # check to make sure the potential is increasing
            if k[2] > 0:
                # place fit data into potential
	        for j in range(stop,n_bins):
	            potential[j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
            # keep increasing size of fit region until k[2] > 0
            else:
                fit_size = 5
                while k[2] < 0:
	            k, rss, rank, sv = np.linalg.lstsq(coeff_mat[stop-fit_size:stop],potential[stop-fit_size:stop])
                    fit_size += 1
                # place fit data into potential
	        for j in range(stop, n_bins):
		    potential[j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
        min_val = np.amin(potential) 
        # now smooth the total  potential using a cubic spline fit	
        tck = interpolate.splrep(x,potential-min_val)
        potential[:] = interpolate.splev(x, tck, der=0)
        force[:] = -interpolate.splev(x, tck, der=1)


# read restart potentials
def ReadRestartPotentials(prev_cg_potentials, prev_cg_start_stop, x_min, x_delta, ib_iter, char):
	n_potentials = prev_cg_potentials.shape[0]
	n_bins = prev_cg_potentials.shape[1]
        x_mat = np.empty(n_bins,dtype=float)
        # populate x values into array
        for i in range(n_bins):
            x_mat[i] = x_min + (i+0.5) * x_delta
        # read potential
	ener_in = char + str(ib_iter) + ".ener"
#	pot_in = open(ener_in, 'r')
#        bin_count = 0
#        input_cg_potentials = []
#        input_x_mat = []
#        for line in pot_in:
#                input_x_mat.append(line[0:20])
#                input_cg_potentials.append([])
#		for potential in range(n_potentials):
#                    input_cg_potentials[potential].append(line[20*(potential+1):20*(potential+2)])
#	pot_in.close()
#        input_cg_potentials = np.asmatrix(input_cg_potentials,dtype=float)
#        input_x_mat = np.asarray(input_x_mat,dtype=float)
#
        data = np.loadtxt(ener_in)
        print data.shape

        # interpolation is used so that we can change grid spacing if necessary
        for potential in range(n_potentials):
                # interpolate so that we can make grid spacing arbitrary
                tck = interpolate.splrep(data[:,0],data[:,potential+1])
                prev_cg_potentials[potential,:] = interpolate.splev(x_mat,tck,der=0)

        # read potential start stop
	start_stop_in = char + str(ib_iter) + ".start_stop"
        prev_cg_start_stop[:,:] = np.loadtxt(start_stop_in)
#	start_stop_file = open(start_stop_in, 'r')
#        potential_count = 0
#        start_stop_list = []
#        for line in start_stop_file:
#            start_stop_list.append([])
#            start_stop_list[potential_count].append(line[0:20])
#            start_stop_list[potential_count].append(line[21:40])
#            potential_count += 1
#        start_stop_file.close()
#        prev_cg_start_stop = np.asmatrix(start_stop_list,dtype=int)



# run CG simulation using lammps
def RunCGSim(ib_iter, equil_in_file, run_in_file):
	
	new_equil_file = "equil.iter" + str(ib_iter) + ".in"
	new_run_file = "run.iter" + str(ib_iter) + ".in"
	equil_log_file = "equil.iter" + str(ib_iter) + ".log"
	run_log_file = "run.iter" + str(ib_iter) + ".log"
        # copy current ib files to active ib files (read in by TEMPLATE cfg files)
        command0 = "cp bond" + str(ib_iter) + ".ib bonds.ib"
        command1 = "cp angle" + str(ib_iter) + ".ib angles.ib"
        command2 = "cp dihedral" + str(ib_iter) + ".ib dihedrals.ib"
        # update equil file and run
	command3 = "sed -e s/NUM/" + str(ib_iter) + "/ < " + equil_in_file + " > " +  new_equil_file
	command4 = "mpirun -np 3 lmp_mac_mpi -i " + new_equil_file + " > " + equil_log_file
        # update run file and run
	command5 = "sed -e s/NUM/" + str(ib_iter) + "/ < " + run_in_file + " > " +  new_run_file
	command6 = "mpirun -np 3 lmp_mac_mpi -i " + new_run_file + " > " + run_log_file
	
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
	print command5
	os.system(command5)
	print command6
	os.system(command6)
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
print "Number of angles:", n_angles
print "Number of unique angles:", n_uniq_angles
print "Number of dihedrals:", n_dihedrals
print "Number of unique dihedrals:", n_uniq_dihedrals

# declare bond, angle and dihedral potentials
bond_min = 0.0
bond_max = 24.0
bond_delta = 0.1 
n_bond_bins  = int((bond_max-bond_min)/bond_delta)
atom_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
atom_bond_start_stop = np.empty((n_uniq_bonds,2),dtype=int)
cg_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
cg_bond_start_stop = np.empty((n_uniq_bonds,2),dtype=int)
prev_cg_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
prev_cg_bond_start_stop = np.empty((n_uniq_bonds,2),dtype=int)
# angles
angle_min = 0.0
angle_max = 180 # inclusive
angle_delta = 1.0 
n_angle_bins  = int((angle_max-angle_min)/angle_delta) + 1
atom_angle_potentials = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
atom_angle_start_stop = np.empty((n_uniq_angles,2),dtype=int)
cg_angle_potentials = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
cg_angle_start_stop = np.empty((n_uniq_angles,2),dtype=int)
prev_cg_angle_potentials = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
prev_cg_angle_start_stop = np.empty((n_uniq_angles,2),dtype=int)
# dihedrals
dihedral_min = -179.0 
dihedral_max = 180.0 # inclusive
dihedral_delta = 1.0
n_dihedral_bins  = int((dihedral_max-dihedral_min)/dihedral_delta) + 1
atom_dihedral_potentials = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
atom_dihedral_start_stop = np.empty((n_uniq_dihedrals,2),dtype=int)
cg_dihedral_potentials = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
cg_dihedral_start_stop = np.zeros((n_uniq_dihedrals,2),dtype=int)
prev_cg_dihedral_potentials = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
prev_cg_dihedral_start_stop = np.zeros((n_uniq_dihedrals,2),dtype=int)

#MMMMMMMM

# read atomistic probability distributions and generate initial (inverse Boltzmann) potentials
ReadAtomHists(atom_bond_potentials,atom_bond_start_stop,atom_angle_potentials, atom_angle_start_stop, atom_dihedral_potentials, atom_dihedral_start_stop)
# set previous potentials to atom potentials
prev_cg_bond_potentials[:,:] = atom_bond_potentials[:,:]
prev_cg_bond_start_stop[:,:] = atom_bond_start_stop[:,:]
prev_cg_angle_potentials[:,:] = atom_angle_potentials[:,:]
prev_cg_angle_start_stop[:,:] = atom_angle_start_stop[:,:]
prev_cg_dihedral_potentials[:,:] = atom_dihedral_potentials[:,:]
prev_cg_dihedral_start_stop[:,:] = atom_dihedral_start_stop[:,:]

# read restrart if needed
if n_start > 0:
    ReadRestartPotentials(prev_cg_bond_potentials, prev_cg_bond_start_stop, bond_min, bond_delta, n_start, "bond")
    ReadRestartPotentials(prev_cg_angle_potentials, prev_cg_angle_start_stop, angle_min, angle_delta, n_start, "angle")
    ReadRestartPotentials(prev_cg_dihedral_potentials, prev_cg_dihedral_start_stop, dihedral_min, dihedral_delta, n_start, "dihedral")

equil_in_file = "equil.TEMPLATE.in"
run_in_file = "run.TEMPLATE.in"

delta_iter = n_iter - n_start
# run bond IB
for ib_iter in range(n_start, n_iter):
	traj_file = "run.iter" + str(ib_iter) + ".out.dcd"
    	# run CG simulation
	RunCGSim(ib_iter, equil_in_file, run_in_file)

	# compute CG probability distributions
	AnalyzeCGTraj(top_file, traj_file,cg_bond_potentials, cg_bond_start_stop, cg_angle_potentials, cg_angle_start_stop, cg_dihedral_potentials, cg_dihedral_start_stop)

	# update potentials
	UpdateBonds(ib_iter, atom_bond_potentials, atom_bond_start_stop, prev_cg_bond_potentials, prev_cg_bond_start_stop, cg_bond_potentials, cg_bond_start_stop)
    
        # update previous potentials
        prev_cg_bond_potentials[:,:] = cg_bond_potentials[:,:]
        prev_cg_bond_start_stop[:,:] = cg_bond_start_stop[:,:]
        prev_cg_angle_potentials[:,:] = cg_angle_potentials[:,:]
        prev_cg_angle_start_stop[:,:] = cg_angle_start_stop[:,:]
        prev_cg_dihedral_potentials[:,:] = cg_dihedral_potentials[:,:]
        prev_cg_dihedral_start_stop[:,:] = cg_dihedral_start_stop[:,:]
        # mv ib files for angles and dihedrals
        command0 = "cp angle" + str(ib_iter) + ".ib angle" + str(ib_iter+1) + ".ib"
	os.system(command0)
        command0 = "cp angle" + str(ib_iter) + ".ener angle" + str(ib_iter+1) + ".ener"
	os.system(command0)
        command0 = "cp angle" + str(ib_iter) + ".start_stop angle" + str(ib_iter+1) + ".start_stop"
	os.system(command0)
        command0 = "cp dihedral" + str(ib_iter) + ".ib dihedral" + str(ib_iter+1) + ".ib"
	os.system(command0)
        command0 = "cp dihedral" + str(ib_iter) + ".ener dihedral" + str(ib_iter+1) + ".ener"
	os.system(command0)
        command0 = "cp dihedral" + str(ib_iter) + ".start_stop dihedral" + str(ib_iter+1) + ".start_stop"
	os.system(command0)

# run angle IB
for ib_iter in range(n_iter, n_iter+delta_iter):
	traj_file = "run.iter" + str(ib_iter) + ".out.dcd"
    	# run CG simulation
	RunCGSim(ib_iter, equil_in_file, run_in_file)

	# compute CG probability distributions
	AnalyzeCGTraj(top_file, traj_file,cg_bond_potentials, cg_bond_start_stop, cg_angle_potentials, cg_angle_start_stop, cg_dihedral_potentials, cg_dihedral_start_stop)

	# update potentials
	UpdateAngles(ib_iter, atom_angle_potentials, atom_angle_start_stop, prev_cg_angle_potentials, prev_cg_angle_start_stop, cg_angle_potentials, cg_angle_start_stop)
    
        # update previous potentials
        prev_cg_bond_potentials[:,:] = cg_bond_potentials[:,:]
        prev_cg_bond_start_stop[:,:] = cg_bond_start_stop[:,:]
        prev_cg_angle_potentials[:,:] = cg_angle_potentials[:,:]
        prev_cg_angle_start_stop[:,:] = cg_angle_start_stop[:,:]
        prev_cg_dihedral_potentials[:,:] = cg_dihedral_potentials[:,:]
        prev_cg_dihedral_start_stop[:,:] = cg_dihedral_start_stop[:,:]
        # mv ib files for angles and dihedrals
        command0 = "cp bond" + str(ib_iter) + ".ib bond" + str(ib_iter+1) + ".ib"
	os.system(command0)
        command0 = "cp bond" + str(ib_iter) + ".ener bond" + str(ib_iter+1) + ".ener"
	os.system(command0)
        command0 = "cp bond" + str(ib_iter) + ".start_stop bond" + str(ib_iter+1) + ".start_stop"
	os.system(command0)
        command0 = "cp dihedral" + str(ib_iter) + ".ib dihedral" + str(ib_iter+1) + ".ib"
	os.system(command0)
        command0 = "cp dihedral" + str(ib_iter) + ".ener dihedral" + str(ib_iter+1) + ".ener"
	os.system(command0)
        command0 = "cp dihedral" + str(ib_iter) + ".start_stop dihedral" + str(ib_iter+1) + ".start_stop"
	os.system(command0)


# run dihedral IB
for ib_iter in range(n_iter+delta_iter,n_iter+2*delta_iter):
	traj_file = "run.iter" + str(ib_iter) + ".out.dcd"
    	# run CG simulation
	RunCGSim(ib_iter, equil_in_file, run_in_file)

	# compute CG probability distributions
	AnalyzeCGTraj(top_file, traj_file,cg_bond_potentials, cg_bond_start_stop, cg_angle_potentials, cg_angle_start_stop, cg_dihedral_potentials, cg_dihedral_start_stop)

	# update potentials
	UpdateDihedrals(ib_iter, atom_dihedral_potentials, atom_dihedral_start_stop, prev_cg_dihedral_potentials, prev_cg_dihedral_start_stop, cg_dihedral_potentials, cg_dihedral_start_stop)
    
        # update previous potentials
        prev_cg_bond_potentials[:,:] = cg_bond_potentials[:,:]
        prev_cg_bond_start_stop[:,:] = cg_bond_start_stop[:,:]
        prev_cg_angle_potentials[:,:] = cg_angle_potentials[:,:]
        prev_cg_angle_start_stop[:,:] = cg_angle_start_stop[:,:]
        prev_cg_dihedral_potentials[:,:] = cg_dihedral_potentials[:,:]
        prev_cg_dihedral_start_stop[:,:] = cg_dihedral_start_stop[:,:]
        # mv ib files for angles and dihedrals
        command0 = "cp angle" + str(ib_iter) + ".ib angle" + str(ib_iter+1) + ".ib"
	os.system(command0)
        command0 = "cp angle" + str(ib_iter) + ".ener angle" + str(ib_iter+1) + ".ener"
	os.system(command0)
        command0 = "cp angle" + str(ib_iter) + ".start_stop angle" + str(ib_iter+1) + ".start_stop"
	os.system(command0)
        command0 = "cp bond" + str(ib_iter) + ".ib bond" + str(ib_iter+1) + ".ib"
	os.system(command0)
        command0 = "cp bond" + str(ib_iter) + ".ener bond" + str(ib_iter+1) + ".ener"
	os.system(command0)
        command0 = "cp bond" + str(ib_iter) + ".start_stop bond" + str(ib_iter+1) + ".start_stop"
	os.system(command0)

