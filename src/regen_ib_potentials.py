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
	global ff_out_root, n_iter, top_file
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
			elif option.lower()=='ffroot':
				ff_out_root = value
			elif option.lower()=='iteration':
				n_iter = int(value)
			else :
				print "Option:", option, " is not recognized"
	f.close()
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
					for bond2 in range(n_uniq_bonds):
		                                # check to see if the opposite strand parameters also exist and add to histogram
		                                if uniq_bond_atom_types[bond2][0][1] == "1":
			                            type1 = uniq_bond_atom_types[bond2][0][0] + str(2)
		                                else:
			                            type1 = uniq_bond_atom_types[bond2][0][0] + str(1)
		                                if uniq_bond_atom_types[bond2][1][1] == "1":
			                            type2 = uniq_bond_atom_types[bond2][1][0] + str(2)
		                                else:
			                            type2 = uniq_bond_atom_types[bond2][1][0] + str(1)
                                                # check if we already have this bond
						if (atom_types[int(bonds[bond][0])-1] == uniq_bond_atom_types[bond2][0] and atom_types[int(bonds[bond][1])-1] == uniq_bond_atom_types[bond2][1]) or (atom_types[int(bonds[bond][0])-1] == uniq_bond_atom_types[bond2][1] and atom_types[int(bonds[bond][1])-1] == uniq_bond_atom_types[bond2][0]):
							same = "true"
							uniq_bond_num = bond2
							break
						elif (atom_types[int(bonds[bond][0])-1] == type1 and atom_types[int(bonds[bond][1])-1] == type2) or (atom_types[int(bonds[bond][0])-1] == type2 and atom_types[int(bonds[bond][1])-1] == type1):
							same = "true"
							uniq_bond_num = bond2
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
					for angle2 in range(n_uniq_angles):
		                                # check to see if the opposite strand parameters also exist and add to histogram
		                                if uniq_angle_atom_types[angle2][0][1] == "1":
			                            type1 = uniq_angle_atom_types[angle2][0][0] + str(2)
		                                else:
			                            type1 = uniq_angle_atom_types[angle2][0][0] + str(1)
		                                if uniq_angle_atom_types[angle2][1][1] == "1":
			                            type2 = uniq_angle_atom_types[angle2][1][0] + str(2)
		                                else:
			                            type2 = uniq_angle_atom_types[angle2][1][0] + str(1)
		                                if uniq_angle_atom_types[angle2][2][1] == "1":
			                            type3 = uniq_angle_atom_types[angle2][2][0] + str(2)
		                                else:
			                            type3 = uniq_angle_atom_types[angle2][2][0] + str(1)
						if (atom_types[int(angles[angle][0])-1] == uniq_angle_atom_types[angle2][0] and atom_types[int(angles[angle][1])-1] == uniq_angle_atom_types[angle2][1] and atom_types[int(angles[angle][2])-1] == uniq_angle_atom_types[angle2][2]) or (atom_types[int(angles[angle][0])-1] == uniq_angle_atom_types[angle2][2] and atom_types[int(angles[angle][1])-1] == uniq_angle_atom_types[angle2][1] and atom_types[int(angles[angle][2])-1] == uniq_angle_atom_types[angle2][0]):
							same = "true"
							uniq_angle_num = angle2
							break
						elif (atom_types[int(angles[angle][0])-1] == type1 and atom_types[int(angles[angle][1])-1] == type2 and atom_types[int(angles[angle][2])-1] == type3) or (atom_types[int(angles[angle][0])-1] == type3 and atom_types[int(angles[angle][1])-1] == type2 and atom_types[int(angles[angle][2])-1] == type1):
							same = "true"
							uniq_angle_num = angle2
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
					for dihedral2 in range(n_uniq_dihedrals):
		                                # check to see if opposite signed strands exist and add to current histogram
		                                if uniq_dihedral_atom_types[dihedral2][0][1] == "1":
			                            type1 = uniq_dihedral_atom_types[dihedral2][0][0] + str(2)
		                                else:
			                            type1 = uniq_dihedral_atom_types[dihedral2][0][0] + str(1)
		                                if uniq_dihedral_atom_types[dihedral2][1][1] == "1":
			                            type2 = uniq_dihedral_atom_types[dihedral2][1][0] + str(2)
		                                else:
			                            type2 = uniq_dihedral_atom_types[dihedral2][1][0] + str(1)
		                                if uniq_dihedral_atom_types[dihedral2][2][1] == "1":
			                            type3 = uniq_dihedral_atom_types[dihedral2][2][0] + str(2)
		                                else:
			                            type3 = uniq_dihedral_atom_types[dihedral2][2][0] + str(1)
		                                if uniq_dihedral_atom_types[dihedral2][3][1] == "1":
			                            type4 = uniq_dihedral_atom_types[dihedral2][3][0] + str(2)
		                                else:
			                            type4 = uniq_dihedral_atom_types[dihedral2][3][0] + str(1)
						if (atom_types[int(dihedrals[dihedral][0])-1] == uniq_dihedral_atom_types[dihedral2][0] and atom_types[int(dihedrals[dihedral][1])-1] == uniq_dihedral_atom_types[dihedral2][1] and atom_types[int(dihedrals[dihedral][2])-1] == uniq_dihedral_atom_types[dihedral2][2] and atom_types[int(dihedrals[dihedral][3])-1] == uniq_dihedral_atom_types[dihedral2][3]) or (atom_types[int(dihedrals[dihedral][0])-1] == uniq_dihedral_atom_types[dihedral2][3] and atom_types[int(dihedrals[dihedral][1])-1] == uniq_dihedral_atom_types[dihedral2][2] and atom_types[int(dihedrals[dihedral][2])-1] == uniq_dihedral_atom_types[dihedral2][1] and atom_types[int(dihedrals[dihedral][3])-1] == uniq_dihedral_atom_types[dihedral2][0]):
							same = "true"
							uniq_dihedral_num = dihedral2
							break
						elif (atom_types[int(dihedrals[dihedral][0])-1] == type1 and atom_types[int(dihedrals[dihedral][1])-1] == type2 and atom_types[int(dihedrals[dihedral][2])-1] == type3 and atom_types[int(dihedrals[dihedral][3])-1] == type4) or (atom_types[int(dihedrals[dihedral][0])-1] == type4 and atom_types[int(dihedrals[dihedral][1])-1] == type3 and atom_types[int(dihedrals[dihedral][2])-1] == type2 and atom_types[int(dihedrals[dihedral][3])-1] == type1):
							same = "true"
							uniq_dihedral_num = dihedral2
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




## MM write potential energy files
def WriteBondPotentials(cg_bond_potentials, bond_min, bond_delta, uniq_bond_atom_types):
        global ff_out_root
	n_bonds = cg_bond_potentials.shape[0]
	n_bins = cg_bond_potentials.shape[1]

	for bond in range(n_bonds):
		filename = ff_out_root+"/bonds/"+uniq_bond_atom_types[bond][0].strip()+"_"+uniq_bond_atom_types[bond][1].strip()+".ib"
		out = open(filename,'w')
		for i in range(n_bins):
			out.write("%10.5f %20.2f\n" % (bond_min+(i+0.5)*bond_delta, cg_bond_potentials[bond,i]))
		out.close()

def WriteAnglePotentials(cg_angle_potentials, angle_min, angle_delta, uniq_angle_atom_types):
        global ff_out_root
	n_angles = cg_angle_potentials.shape[0]
	n_bins = cg_angle_potentials.shape[1]

	for angle in range(n_angles):
		filename = ff_out_root+"/angles/"+uniq_angle_atom_types[angle][0].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][2].strip()+".ib"
		out = open(filename,'w')
		for i in range(n_bins):
			out.write("%10.5f %20.2f\n" % (angle_min+(i)*angle_delta, cg_angle_potentials[angle,i]))
		out.close()

def WriteDihedralPotentials(cg_dihedral_potentials, dihedral_min, dihedral_delta, uniq_dihedral_atom_types):
        global ff_out_root
	n_dihedrals = cg_dihedral_potentials.shape[0]
	n_bins = cg_dihedral_potentials.shape[1]

	for dihedral in range(n_dihedrals):
		filename = ff_out_root+"/dihs/"+uniq_dihedral_atom_types[dihedral][0].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][3].strip()+".ib"
		out = open(filename,'w')
		for i in range(n_bins):
			out.write("%10.5f %20.2f\n" % (dihedral_min+(i)*dihedral_delta, cg_dihedral_potentials[dihedral,i]))
		out.close()

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
bond_delta = 0.01 
n_bond_bins  = int((bond_max-bond_min)/bond_delta)
cg_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
cg_bond_start_stop = np.empty((n_uniq_bonds,2),dtype=int)
# angles
angle_min = 0.0
angle_max = 180 # inclusive
angle_delta = 1.0 
n_angle_bins  = int((angle_max-angle_min)/angle_delta) + 1
cg_angle_potentials = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
cg_angle_start_stop = np.empty((n_uniq_angles,2),dtype=int)
# dihedrals
dihedral_min = -179.0 
dihedral_max = 180.0 # inclusive
dihedral_delta = 1.0
n_dihedral_bins  = int((dihedral_max-dihedral_min)/dihedral_delta) + 1
cg_dihedral_potentials = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
cg_dihedral_start_stop = np.zeros((n_uniq_dihedrals,2),dtype=int)

#MMMMMMMM

# read the potentials
ReadRestartPotentials(cg_bond_potentials, cg_bond_start_stop, bond_min, bond_delta, n_iter, "bond")
ReadRestartPotentials(cg_angle_potentials, cg_angle_start_stop, angle_min, angle_delta, n_iter, "angle")
ReadRestartPotentials(cg_dihedral_potentials, cg_dihedral_start_stop, dihedral_min, dihedral_delta, n_iter, "dihedral")

# print the potential files
WriteBondPotentials(cg_bond_potentials, bond_min, bond_delta, uniq_bond_atom_types)
WriteAnglePotentials(cg_angle_potentials, angle_min, angle_delta, uniq_angle_atom_types)
WriteDihedralPotentials(cg_dihedral_potentials, dihedral_min, dihedral_delta, uniq_dihedral_atom_types)

