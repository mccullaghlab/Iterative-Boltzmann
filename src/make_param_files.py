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
bond_factor = 0.2
angle_factor = 0.2
dihedral_factor = 0.30

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global top_file, out_file, ff_data_root
	f = open(cfg_file)
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
			elif option.lower()=='outfile':
				out_file = value
			elif option.lower()=='ff_data_root':
				ff_data_root = value
			else :
				print "Option:", option, " is not recognized"
	f.close()

def ParsePsfFile(psf_file):
	global n_uniq_atom_types, atom_types, uniq_bond_atom_types, n_uniq_bonds, bonds, n_bonds, uniq_angle_atom_types, n_uniq_angles, angles, n_angles,  uniq_dihedral_atom_types, n_uniq_dihedrals, dihedrals, n_dihedrals
	f = open(psf_file)
	atom_flag = bond_flag = angle_flag = dihedral_flag = "false"
	atom = bond = angle = dihedral = 0
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
					for bond2 in range(n_uniq_bonds):
						if (atom_types[int(bonds[bond][0])-1] == uniq_bond_atom_types[bond2][0] and atom_types[int(bonds[bond][1])-1] == uniq_bond_atom_types[bond2][1]) or (atom_types[int(bonds[bond][0])-1] == uniq_bond_atom_types[bond2][1] and atom_types[int(bonds[bond][1])-1] == uniq_bond_atom_types[bond2][0]):
							same = "true"
							uniq_bond_num = bonds[bond2][2]
							break
				if same == "false":
					uniq_bond_atom_types.append([])
					uniq_bond_atom_types[n_uniq_bonds].append(atom_types[int(bonds[bond][0])-1])
					uniq_bond_atom_types[n_uniq_bonds].append(atom_types[int(bonds[bond][1])-1])
					uniq_bond_num = str(n_uniq_bonds)
					print n_uniq_bonds+1, uniq_bond_atom_types[n_uniq_bonds]
					n_uniq_bonds += 1
				bonds[bond].append(uniq_bond_num)
				bond += 1
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
					for angle2 in range(n_uniq_angles):
						if (atom_types[int(angles[angle][0])-1] == uniq_angle_atom_types[angle2][0] and atom_types[int(angles[angle][1])-1] == uniq_angle_atom_types[angle2][1] and atom_types[int(angles[angle][2])-1] == uniq_angle_atom_types[angle2][2]) or (atom_types[int(angles[angle][0])-1] == uniq_angle_atom_types[angle2][2] and atom_types[int(angles[angle][1])-1] == uniq_angle_atom_types[angle2][1] and atom_types[int(angles[angle][2])-1] == uniq_angle_atom_types[angle2][0]):
							same = "true"
							uniq_angle_num = angles[angle2][3]
							break
				if same == "false":
					uniq_angle_atom_types.append([])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][0])-1])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][1])-1])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][2])-1])
#					print n_uniq_angles+1, uniq_angle_atom_types[n_uniq_angles]
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
					for dihedral2 in range(n_uniq_dihedrals):
						if (atom_types[int(dihedrals[dihedral][0])-1] == uniq_dihedral_atom_types[dihedral2][0] and atom_types[int(dihedrals[dihedral][1])-1] == uniq_dihedral_atom_types[dihedral2][1] and atom_types[int(dihedrals[dihedral][2])-1] == uniq_dihedral_atom_types[dihedral2][2] and atom_types[int(dihedrals[dihedral][3])-1] == uniq_dihedral_atom_types[dihedral2][3]) or (atom_types[int(dihedrals[dihedral][0])-1] == uniq_dihedral_atom_types[dihedral2][3] and atom_types[int(dihedrals[dihedral][1])-1] == uniq_dihedral_atom_types[dihedral2][2] and atom_types[int(dihedrals[dihedral][2])-1] == uniq_dihedral_atom_types[dihedral2][1] and atom_types[int(dihedrals[dihedral][3])-1] == uniq_dihedral_atom_types[dihedral2][0]):
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
				if dihedral%1000 == 0:
					print "Dihedral ", dihedral, "of ", n_dihedrals
				sys.stdout.flush()
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

# average dihedral distance histograms (convert them to probability densities)
def WriteDihedralPotentials(dihedral_potentials, dihedral_min, dihedral_delta, ib_out, param_out):

	n_dihedrals = dihedral_potentials.shape[0]
	n_bins = dihedral_potentials.shape[1]

	force = np.empty(n_bins,dtype=float)
	x_mat = np.empty(n_bins,dtype=float)

	out = open(ib_out,'w')
	# first check to see if we have two copies of same dihedral
	for dihedral in range(n_dihedrals):
		for i in range(n_bins):
			x = dihedral_min+(i+0.5)*dihedral_delta
			x_mat[i] = x
		# write title to output 
#		out.write("\n%5d\n" % (dihedral+1))
		out.write("\nDIH_%s\n" % (str(dihedral+1)))
		out.write("N %5d DEGREES\n\n" % (n_bins))

		# compute forces
		tck = interpolate.splrep(x_mat,dihedral_potentials[dihedral,:])
		force[:] = -interpolate.splev(x_mat, tck, der=1)
		for i in range(n_bins):
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,dihedral_min+(i)*dihedral_delta, dihedra_potentials[dihedral,i],force[i]))
		param_out.write("dihedral_coeff %2d %s DIH_%s\n" %(dihedral+1, ib_out, str(dihedral+1).strip()))

	out.close()
# MM_DIHEDRALS

# average angle distance histograms (convert them to probability densities)
def WriteAnglePotentials(angle_potentials, angle_min, angle_delta, ib_out, param_out):

	n_angles = angle_potentials.shape[0]
	n_bins = angle_potentials.shape[1]

	force = np.empty(n_bins,dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	out = open(ib_out,'w')
	# first check to see if we have two copies of same angle
	for angle in range(n_angles):
		for i in range(n_bins):
			x = angle_min+(i+0.5)*angle_delta
			x_mat[i] = x
		out.write("\nANG_%s\n" % (str(angle+1)))
		out.write("N %5d\n\n" % (n_bins))

		# first connect point
                tck = interpolate.splrep(x_mat,angle_potentials[angle,:])
		force[:] = -interpolate.splev(x_mat, tck, der=1)
		for i in range(n_bins):
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,angle_min+(i)*angle_delta, angle_potentials[angle,i],force[i]))
		param_out.write("angle_coeff %2d %s ANG_%s\n" %(angle+1, ib_out, str(angle+1).strip()))
	out.close()

# MM_ANGLES

# average bond distance histograms (convert them to probability densities)
def WriteBondPotentials(bond_potentials, bond_min, bond_delta, ib_out,ener_out,param_out):

	n_bonds = bond_potentials.shape[0]
	n_bins = bond_potentials.shape[1]
	force = np.empty(n_bins,dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	
	out = open(ib_out,'w')
	for bond in range(n_bonds):
		for i in range(n_bins):
			x = bond_min+(i+0.5)*bond_delta
			x_mat[i] = x 
		# write title to output 
		out.write("\nBOND_%s\n" % (str(bond+1)))
		out.write("N %5d\n\n" % (n_bins))
                tck = interpolate.splrep(x_mat,bond_potentials[bond,:])
                force[:] = -interpolate.splev(x_mat, tck, der=1)
		for i in range(n_bins):
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,bond_min+(i+0.5)*bond_delta, bond_potentials[bond,i],force[i]))
		param_out.write("bond_coeff %2d %s BOND_%s\n" %(bond+1, ib_out, str(bond+1).strip()))
	out.close()


# average dihedral distance histograms (convert them to probability densities)
def ReadDihedralPotentials(dihedral_potentials, dihedral_min, dihedral_delta, uniq_dihedral_atom_types):

	n_dihedrals = dihedral_potentials.shape[0]
	n_bins = dihedral_potentials.shape[1]


	for dihedral in range(n_dihedrals):
		filename = ff_data_root+"/dihs/"+uniq_dihedral_atom_types[dihedral][0].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][3].strip()+".ib"
		filename2 = ff_data_root+"/dihs/"+uniq_dihedral_atom_types[dihedral][3].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][0].strip()+".ib"
		if os.path.isfile(filename):
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					dihedral_potentials[dihedral,count] += float(val)
					count += 1
				inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					dihedral_potentials[dihedral,count] += float(val)
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
			filename = ff_data_root+"/dihs/"+type1+"_"+type2+"_"+type3+"_"+type4+".ib"
			filename2 = ff_data_root+"/dihs/"+type4+"_"+type3+"_"+type2+"_"+type1+".ib"
			if os.path.isfile(filename):
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					dihedral_potentials[dihedral,count] = float(val)
					count += 1
				inp.close()
			elif os.path.isfile(filename2):
				filename = filename2
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					dihedral_potentials[dihedral,count] = float(val)
					count += 1
				inp.close()
			else:
				print "Did not find data for dihedral", uniq_dihedral_atom_types[dihedral]
				sys.exit()


# add to angle distance histograms 
def ComputeAngleHists(atom_positions, angles, angle_potentials, angle_min, angle_delta):
	# get size info
	n_angles = angles.shape[0]
	n_bins = angle_potentials.shape[1]
	for angle in range(n_angles):

		ang = ComputeAng(atom_positions[angles[angle,0]-1,:],atom_positions[angles[angle,1]-1,:],atom_positions[angles[angle,2]-1,:])
		ang_bin = int((ang-angle_min)/angle_delta)
		if ang_bin >= 0 and ang_bin < n_bins:
			angle_potentials[angles[angle,3],ang_bin] += 1

# average angle distance histograms (convert them to probability densities)
def ReadAnglePotentials(angle_potentials, angle_min, angle_delta, uniq_angle_atom_types):

	n_angles = angle_potentials.shape[0]
	n_bins = angle_potentials.shape[1]


	for angle in range(n_angles):
		filename = ff_data_root+"/angs/"+uniq_angle_atom_types[angle][0].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][2].strip()+".ib"
		filename2 = ff_data_root+"/angs/"+uniq_angle_atom_types[angle][2].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][0].strip()+".ib"
		if os.path.isfile(filename):
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				angle_potentials[angle,count] = float(val)
				count += 1
			inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				angle_potentials[angle,count] = float(val)
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
			filename = ff_data_root+"/angs/"+type1+"_"+type2+"_"+type3+".ib"
			filename2 = ff_data_root+"/angs/"+type3+"_"+type2+"_"+type1+".ib"
			if os.path.isfile(filename):
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					angle_potentials[angle,count] = float(val)
					count += 1
				inp.close()
			elif os.path.isfile(filename2):
				filename = filename2
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					angle_potentials[angle,count] = float(val)
					count += 1
				inp.close()
			else:
				print "Did not find data for angle", uniq_angle_atom_types[angle]
				sys.exit()


# average bond distance histograms (convert them to probability densities)
def ReadBondPotentials(bond_potentials, bond_min, bond_delta, uniq_bond_atom_types):

	n_bonds = bond_potentials.shape[0]
	n_bins = bond_potentials.shape[1]

	for bond in range(n_bonds):
		filename = ff_data_root+"/bonds/"+uniq_bond_atom_types[bond][0].strip()+"_"+uniq_bond_atom_types[bond][1].strip()+".ib"
		filename2 = ff_data_root+"/bonds/"+uniq_bond_atom_types[bond][1].strip()+"_"+uniq_bond_atom_types[bond][0].strip()+".ib"
		if os.path.isfile(filename):
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				bond_potentials[bond,count] += float(val)
				count += 1
			inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				bond_potentials[bond,count] += float(val)
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
			filename = ff_data_root+"/bonds/"+type1+"_"+type2+".ib"
			filename2 = ff_data_root+"/bonds/"+type2+"_"+type1+".ib"
			if os.path.isfile(filename):
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					bond_potentials[bond,count] = float(val)
					count += 1
				inp.close()
			elif os.path.isfile(filename2):
				filename = filename2
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					bond_potentials[bond,count] = float(val)
					count += 1
				inp.close()
			else:
				print "Did not find parameters for:", uniq_bond_atom_types[bond]
				sys.exit()

#############################################################################################################################################################
####################################################              MAIN PROGRAM               ################################################################
#############################################################################################################################################################

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Topology file:", top_file
print "Output data file:", out_file

# parse psf file
ParsePsfFile(top_file)

print "Number of unique atom types:", n_uniq_atom_types
print "Number of bonds:", n_bonds
print "Number of unique bonds:", n_uniq_bonds
print "Number of angles:", n_angles
print "Number of unique angles:", n_uniq_angles
print "Number of dihedrals:", n_dihedrals
print "Number of unique dihedrals:", n_uniq_dihedrals

# declare bond, angle and dihedral histograms
bond_min = 0.0
bond_max = 24.0
bond_delta = 0.25
n_bond_bins  = int((bond_max-bond_min)/bond_delta)
bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
angle_min = 0.0
angle_max = 180 # inclusive
angle_delta = 1.0 
n_angle_bins  = int((angle_max-angle_min)/angle_delta) + 1
angle_potentials = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
dihedral_min = -179.0 
dihedral_max = 180.0 # inclusive
dihedral_delta = 1.0
n_dihedral_bins  = int((dihedral_max-dihedral_min)/dihedral_delta) + 1
dihedral_potentials = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)

param_out = open(out_file,'w')

ReadBondPotentials(bond_potentials, bond_min, bond_delta, uniq_bond_atom_types)
ReadAnglePotentials(angle_potentials, angle_min, angle_delta, uniq_angle_atom_types)
ReadDihedralPotentials(dihedral_potentials, dihedral_min, dihedral_delta, uniq_dihedral_atom_types)
# scale 
#bond_potentials *= bond_factor
#angle_potentials *= angle_factor
#dihedral_potentials *= dihedral_factor

WriteBondPotentials(bond_potentials, bond_min, bond_delta, "bonds.ib",param_out)
WriteAnglePotentials(angle_potentials, angle_min, angle_delta, "angles.ib",param_out)
WriteDihedralPotentials(dihedral_potentials, dihedral_min, dihedral_delta, "dihedrals.ib",param_out)
param_out.close()
