import MDAnalysis as mda
from MDAnalysis.analysis import distances

import numpy as np
import matplotlib.pyplot as plt

import warnings

import MDAnalysis as mda
import prolif as plf
import sys
import pandas as pd
import matplotlib.pyplot as plt
from prolif.plotting.network import LigNetwork

#Import structure and trajectory
u = mda.Universe(sys.argv[1], sys.argv[2],continuous = False, guess_bonds= True, vdwradii={"H": 1.05, "O": 1.48, 'Cl': 1.75, 'Na': 2.27})

# create selections for the ligands and protein
prot = u.atoms.select_atoms("protein")
ion = u.atoms.select_atoms("resname CO3")
#bindingsite_ca = u.select_atoms('name CA and (resid 356 )') 
SOD_1 = u.atoms.select_atoms("resid 572")
SOD_2 = u.atoms.select_atoms("resid 573")
#after aligning two structures, D831 is D356
#use the resi within 6 of CO3 from the gmxMMPBSA then

#####distance between CO3 and D831
dists =[]
for ts in u.trajectory:
    ion = u.atoms.select_atoms("resname CO3")
    bindingsite_ca = u.select_atoms('resid 356') #only D831 here
    ###maybe center of mass?
    ion_com = ion.center_of_mass(compound='residues')
    bindingsite_ca_com = bindingsite_ca.center_of_mass(compound='residues')
    dist_arr = distances.distance_array(ion_com, # reference
                                    bindingsite_ca_com, # configuration
                                    box=u.dimensions)
    dists.append(dist_arr)

narr_CO3 = np.array(dists)
#maybe not the best way to do it
narr_simplified_CO3 = np.squeeze(narr_CO3)

num_frames = len(narr_simplified_CO3)
total_time = 250  # Total time in nanoseconds
time_per_frame = total_time / num_frames
#times
times = np.arange(0, num_frames) * time_per_frame
# Create a plot
plt.figure(figsize=(12, 5))
plt.plot(times,narr_simplified_CO3)  # Line plot with markers
plt.title('Distance between CO3 (COM) and D831 (COM)')
plt.xlabel('Time (ns)')
plt.ylabel('Distance(Å)')
plt.grid(True)  # Add grid for better readability
plt.savefig('distance_between_CO3_D831.png')


#####distance between SOD1 and D831
dists_SOD =[]
for ts in u.trajectory:
    SOD_1 = u.atoms.select_atoms("resid 572")
    bindingsite_ca = u.select_atoms('resid 356') #only D831 here
    ###maybe center of mass?
    SOD_1_com = SOD_1.center_of_mass(compound='residues')
    bindingsite_ca_com = bindingsite_ca.center_of_mass(compound='residues')
    dist_arr = distances.distance_array(SOD_1_com, # reference
                                    bindingsite_ca_com, # configuration
                                    box=u.dimensions)
    dists_SOD.append(dist_arr)

narr_SOD1 = np.array(dists_SOD)
#maybe not the best way to do it
narr_simplified_SOD1 = np.squeeze(narr_SOD1)

num_frames = len(narr_simplified_SOD1)
total_time = 300  # Total time in nanoseconds
time_per_frame = total_time / num_frames
#times
times = np.arange(0, num_frames) * time_per_frame
# Create a plot
plt.figure(figsize=(12, 5))
plt.plot(times,narr_simplified_SOD1)  # Line plot with markers
plt.title('Distance between SOD 1 (COM) and D831 (COM)')
plt.xlabel('Time (ns)')
plt.ylabel('Distance(Å)')
plt.grid(True)  # Add grid for better readability
plt.savefig('distance_between_SOD1_D831.png')


#####distance between SOD2 and D831
dists_SOD2 =[]
for ts in u.trajectory:
    SOD_2 = u.atoms.select_atoms("resid 573")
    bindingsite_ca = u.select_atoms('resid 356') #only D831 here
    ###maybe center of mass?
    SOD_2_com = SOD_2.center_of_mass(compound='residues')
    bindingsite_ca_com = bindingsite_ca.center_of_mass(compound='residues')
    dist_arr = distances.distance_array(SOD_2_com, # reference
                                    bindingsite_ca_com, # configuration
                                    box=u.dimensions)
    dists_SOD2.append(dist_arr)

narr_SOD2 = np.array(dists_SOD2)
#maybe not the best way to do it
narr_simplified_SOD2 = np.squeeze(narr_SOD2)

num_frames = len(narr_simplified_SOD2)
total_time = 300  # Total time in nanoseconds
time_per_frame = total_time / num_frames
#times
times = np.arange(0, num_frames) * time_per_frame
# Create a plot
plt.figure(figsize=(12, 5))
plt.plot(times,narr_simplified_SOD2)  # Line plot with markers
plt.title('Distance between SOD 2 (COM) and D831 (COM)')
plt.xlabel('Time (ns)')
plt.ylabel('Distance(Å)')
plt.grid(True)  # Add grid for better readability
plt.savefig('distance_between_SOD2_D831.png')

####combine the three plots and use default y-axis
# Create a single figure for all plots
plt.figure(figsize=(12, 5))
# Plotting each dataset
plt.plot(times, narr_simplified_CO3, color='blue', label='CO3-D831')  
plt.plot(times, narr_simplified_SOD1, color='red', label='SOD1-D831')  # Red for SOD1
plt.plot(times, narr_simplified_SOD2, color='green', label='SOD2-D831')  # Green for SOD2
# Adding titles and labels
plt.title('Distance between CO3, SOD1, SOD2 (COM) and D831 (COM) over Time')
plt.xlabel('Time (ns)')
plt.ylabel('Distance (Å)')
plt.grid(True)
# Setting the y-axis limit
plt.ylim(0, 20)  # This sets the upper limit of y-axis to 20 Å
# Addinga legend
plt.legend()
#plt.show()
# Optionally save the plot
plt.savefig('combined_distances_20A.png')