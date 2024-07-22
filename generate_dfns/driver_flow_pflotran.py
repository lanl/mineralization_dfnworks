import os, sys
from pydfnworks import *
from write_pflotran_files_vol_flow import *
from scipy.stats.qmc import LatinHypercube
import numpy as np

import shutil
import random
import pandas as pd 

def get_backbone(G,H):
        
    bb_nodes = list(np.sort(list(H.nodes())))
    bb_nodes.remove('s')
    bb_nodes.remove('t')

        
    new_bb = []
    for n in bb_nodes:
        try:
            new_bb.append(int(n)-1)
        except:
            new_bb.append(int(float(n))-1)
            
    bb_nodes = new_bb
    return bb_nodes 


def load_mas_file(filename):
    try:
        with open(filename,"r") as fp:
            line = fp.readline()
        
        keynames = line.split(",")
        keynames = [key.strip('\n') for key in keynames]
        keynames = [key.strip('"') for key in keynames]
        keynames = [key.strip(' "') for key in keynames]
        dictionary = dict.fromkeys(keynames, None)
        data = np.genfromtxt(filename,skip_header = 1)

        for i,key in enumerate(keynames):
            dictionary[key] = data[:,i]
        df = pd.DataFrame.from_dict(dictionary)
        return df
    except:
        exit()

home = os.getcwd()
ijob = int(sys.argv[1])

########### LHS SAMPLING ###########

# Num of samples should fill the space appropriately.
num_of_experiments = 500

# Get the input info:
iexp = ijob #% num_of_experiments

# Set the number of samples, and expose the parameter samples (which must be
# the same across simulations.
global latin_hypercube_sample

# Set the min and max values for the parameters we wish to investigate.
global min_inflow_rate
global max_inflow_rate

global min_diffusion_coef
global max_diffusion_coef

global min_porosity
global max_porosity

global min_gypsum_rate_const
global max_gypsum_rate_const

global min_calcite_rate_constant 
global max_calcite_rate_constant 

# Values chosen based on SKB report 10-52 Table 6.75
# TODO: Confirm with Jeffrey that this table is for Crystalline rock.
min_inflow_rate = np.log10(1e-9) ## m^3/s 
max_inflow_rate = np.log10(1e1)

min_diffusion_coef = np.log10(1e-12)
max_diffusion_coef = np.log10(1e-6)

min_porosity = 0.05
max_porosity = 0.50

# Gypsum dissolution rate, mol/m2/sec:
min_gypsum_rate_const = np.log10(1.2e-8)
max_gypsum_rate_const = np.log10(1.62e-3)

# Calcite dissolution rate, mol/m2/sec:
min_calcite_rate_constant = np.log10(2.0e-9)
max_calcite_rate_constant = np.log10(5.0e-5)


min_gypsum_surface_area = 486
max_gypsum_surface_area = 4670

min_calcite_surface_area = np.log10( 83.8)
max_calcite_surface_area = np.log10(2.3e4)

# Get the simulation-wide LHS
latin_hypercube_sampler = LatinHypercube(d=7, seed=13)
latin_hypercube_sample = latin_hypercube_sampler.random(n=num_of_experiments)
# THEN set the seed (in case scipy sets it for the whole file)
random.seed(ijob* 100)

# Grab parameter values for this experiment:
curr_inflow_rate = latin_hypercube_sample[iexp,0]*(max_inflow_rate - min_inflow_rate) + min_inflow_rate
curr_inflow_rate = 10*(10**curr_inflow_rate)

curr_diffusion_coef = latin_hypercube_sample[iexp,1]*(max_diffusion_coef - min_diffusion_coef) + min_diffusion_coef
curr_diffusion_coef = 10**curr_diffusion_coef

curr_porosity = latin_hypercube_sample[iexp,2]*(max_porosity - min_porosity) + min_porosity

curr_gypsum_rate_constant= latin_hypercube_sample[iexp,3]*(max_gypsum_rate_const - min_gypsum_rate_const) + min_gypsum_rate_const
curr_gypsum_rate_constant = 10**curr_gypsum_rate_constant

curr_calcite_rate_constant = latin_hypercube_sample[iexp,4]*(max_calcite_rate_constant - min_calcite_rate_constant) + min_calcite_rate_constant
curr_calcite_rate_constant = 10**curr_calcite_rate_constant  

curr_gypsum_surface_area = latin_hypercube_sample[iexp,5]*(max_gypsum_surface_area - min_gypsum_surface_area) + min_gypsum_surface_area 

curr_calcite_surface_area = latin_hypercube_sample[iexp,6]*(max_calcite_surface_area - min_calcite_surface_area) + min_calcite_surface_area 
curr_calcite_surface_area = 10**curr_calcite_surface_area


print("Sampling Complete")

#print(f"curr_inflow_rate: {curr_inflow_rate:0.5e}")
#print(f"curr_diffusion_coef: {curr_diffusion_coef:0.2e}")
#print(f"curr_porosity: {curr_porosity:0.2f}")
#print(f"curr_gypsum_rate_constant: {curr_gypsum_rate_constant:0.2e}")
#print(f"curr_calcite_rate_constant: {curr_calcite_rate_constant:0.2e}")


########### LHS SAMPLING ###########

# inflow_pressure = 2e6
# diffusion_coef = 1e-9
# gypsum_rate_const =1.2e-8 
# calcite_rate_constant = 2.43e-5 

src_path = os.getcwd() 
#### make network 
path = "/lclscratch/jhyman/rtm_LHS/"
path = os.getcwd() 
jobname =  f"{path}/sample_x{ijob:03d}"
DFN = DFNWORKS(jobname,
               ncpu=10)

DFN.params['domainSize']['value'] = [10,10,10]
DFN.params['h']['value'] = 0.2
DFN.params['orientationOption']['value'] = 1
DFN.params['stopCondition']['value'] = 1
DFN.params['seed']['value'] = ijob # *np.random.randint(1,num_of_experiments )

DFN.params['boundaryFaces']['value'] = [1, 1, 0, 0, 0, 0]
DFN.params['rFram']['value'] = True

DFN.add_fracture_family(shape="rect",
                        distribution="constant",
                        kappa=0.1,
                        aspect=1.0,
                        trend=0.0,
                        plunge=0.0,
                        constant=2.0,
                        p32 = 1.5,
                        hy_variable='aperture',
                        hy_function='constant',
                        hy_params={
                            "mu": 1e-4,
                        })

DFN.make_working_directory(delete=True)
DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()

G = DFN.create_graph('fracture', 'left', 'right')
H = DFN.current_flow_threshold(G,thrs = 1e-16)

bb_nodes = get_backbone(G,H)

backbone_volume = np.sum(DFN.surface_area[bb_nodes] * DFN.aperture[bb_nodes]) 
backbone_p32 = np.sum(DFN.surface_area[bb_nodes]) / (DFN.domain['x'] * DFN.domain['y'] * DFN.domain['z'])

initial_dfn_volume = np.sum(DFN.surface_area * DFN.aperture) 
p32 = np.sum(DFN.surface_area) / (DFN.domain['x'] * DFN.domain['y'] * DFN.domain['z'])

print(f"backbone volume: {backbone_volume}")
print(f"total volume: {initial_dfn_volume}")

DFN.mesh_network(uniform_mesh = True)
DFN.to_pickle()

## Control Parameters
flow_filename = f"dfn_flow_{ijob:02d}.in"
flow_file_path = src_path + os.sep + flow_filename
rxn_filename = f"dfn_rxn_{ijob:02d}.in"
rxn_file_path = src_path+ os.sep +   rxn_filename 

write_pflotran_steady_flow(flow_file_path, curr_inflow_rate, curr_diffusion_coef, curr_porosity)
write_pflotran_rxn(rxn_file_path, ijob, curr_inflow_rate, curr_diffusion_coef, curr_porosity, curr_gypsum_rate_constant, curr_calcite_rate_constant, curr_gypsum_surface_area, curr_calcite_surface_area)

DFN.lagrit2pflotran(boundary_cell_area = 1)
fin = open("boundary_bottom.ex","r")
fin.readline()
line = fin.readline()
fin.close()

fout = open("pinned.ex","w")
fout.write("CONNECTIONS 1\n")
fout.write(line)
fout.close()

DFN.ncpu = 10
DFN.dfnFlow_file = flow_file_path
DFN.local_dfnFlow_file = flow_filename
DFN.pflotran()
DFN.pflotran_cleanup()
DFN.parse_pflotran_vtk_python()

DFN.ncpu = 10
DFN.dfnFlow_file = rxn_file_path
DFN.local_dfnFlow_file = rxn_filename
DFN.pflotran()
DFN.pflotran_cleanup()
DFN.parse_pflotran_vtk_python()

print("* processing run *") 


filename = f"dfn_rxn_{ijob:02d}-mas.dat"
data = load_mas_file(filename)
if data["Time [y]"].to_numpy()[-1] > 0.9:

    shutil.copyfile(filename, home + os.sep + "mas_files_50y" + os. sep + filename)
    data = load_mas_file(filename)
    gypsum_total_mass = data['Region All Gypsum Total Mass [mol]'].array[-1]
    calcite_total_mass = data['Region All Calcite Total Mass [mol]'].array[-1]
    
    sample_info_filename = f'SA_params_{ijob:02d}.txt'
    with open(sample_info_filename, 'w') as fp:
        fp.write(f"dfn_seed: {DFN.params['seed']['value']}\n")
        fp.write(f"curr_inflow_rate: {np.log10(curr_inflow_rate):0.12e}\n")
        fp.write(f"curr_diffusion_coef: {np.log10(curr_diffusion_coef):0.12e}\n")
        fp.write(f"curr_porosity: {curr_porosity:0.2f}\n")
        fp.write(f"curr_gypsum_rate_constant: {np.log10(curr_gypsum_rate_constant):0.12e}\n")
        fp.write(f"curr_calcite_rate_constant: {np.log10(curr_calcite_rate_constant):0.12e}\n")
        fp.write(f"curr_gypsum_surface_area: {curr_gypsum_surface_area:0.12e}\n")
        fp.write(f"curr_calcite_surface_area: {np.log10(curr_calcite_surface_area):0.12e}\n")
        fp.write(f"dfn_p32: {p32}\n")
        fp.write(f"dfn_volume: {initial_dfn_volume}\n")
        fp.write(f"backbone_volume: {backbone_volume}\n")
        fp.write(f"backbone_p32: {backbone_p32}\n")
        fp.write(f"final_gypsum: {gypsum_total_mass:0.12e}\n")
        fp.write(f"final_calcite: {calcite_total_mass:0.12e}\n")
        
    print("* run complete *") 
    
    shutil.copyfile(sample_info_filename, home + os.sep + "output_files_50y" + os. sep + sample_info_filename)
else:
    print("* run failed to reach final time *")
    print(f'final time {data["Time [y]"].to_numpy()[-1]}') 

print("* cleaning up run")
os.chdir(home)
# shutil.rmtree(DFN.jobname)
print(f"* done with {DFN.jobname}") 





