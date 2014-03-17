"""
This file writes all of the materials data (multi-group nuclear 
cross-sections) for the LRA diffusion
benchmark problem to an HDF5 file. The script uses the h5py Python package
to interact with the HDF5 file format. This may be a good example for those
wishing ot write their nuclear data to an HDF5 file to import using the
OpenMOC 'materialize' Python module.
"""


# Create a Python dictionary to store LRA multi-group cross-sections
dataset = {}

# Set global parameters
dataset['Energy Groups'] = 2
dataset['Materials'] = {}

lra_materials = dataset['Materials']

###############################################################################
################################   region 1    ################################
###############################################################################

# Create a subdictionary for region 1 material data
lra_materials['region_1'] = {}

lra_materials['region_1']['FunctionalVariables'] = ['Temperature']

lra_materials['region_1']['Temperature'] = 300.0

lra_materials['region_1']['FunctionalTemperature'] = ['Absorption XS']

lra_materials['region_1']['Absorption XS'] = [0.008252, 0.1003]

lra_materials['region_1']['Total XS'] = [0.2656, 1.5798]

lra_materials['region_1']['Scattering XS'] = [0.232022, 0.02533, 0.00, 1.479479]

lra_materials['region_1']['Fission XS'] = [0.0018938272, 0.044897119]

lra_materials['region_1']['Nu Fission XS'] = [0.004602, 0.1091]

lra_materials['region_1']['Chi'] = [1.0, 0.0]

lra_materials['region_1']['Diffusion Coefficient'] = [1.255, 0.211]

lra_materials['region_1']['Buckling'] = [1e-4, 1e-4]

lra_materials['region_1']['Gamma'] = [3.034e-3, 0.0]


###############################################################################
################################   region 2    ################################
###############################################################################

# Create a subdictionary for region 2 material data
lra_materials['region_2'] = {}

lra_materials['region_2']['FunctionalVariables'] = ['Temperature']

lra_materials['region_2']['Temperature'] = 300.0

lra_materials['region_2']['FunctionalTemperature'] = ['Absorption XS']

lra_materials['region_2']['Absorption XS'] = [0.007181, 0.07047]

lra_materials['region_2']['Total XS'] = [0.2629, 1.7525]

lra_materials['region_2']['Scattering XS'] = [0.22803, 0.02767, 0.00, 1.682071]

lra_materials['region_2']['Fission XS'] = [0.0018967078, 0.035699588]

lra_materials['region_2']['Nu Fission XS'] = [0.004609, 0.08675]

lra_materials['region_2']['Chi'] = [1.0, 0.0]

lra_materials['region_2']['Diffusion Coefficient'] = [1.268, 0.1902]

lra_materials['region_2']['Buckling'] = [1e-4, 1e-4]

lra_materials['region_2']['Gamma'] = [3.034e-3, 0.0]


###############################################################################
################################   region 3    ################################
###############################################################################

# Create a subdictionary for region 3 material data
lra_materials['region_3'] = {}

lra_materials['region_3']['FunctionalVariables'] = ['Temperature']

lra_materials['region_3']['Temperature'] = 300.0

lra_materials['region_3']['FunctionalTemperature'] = ['Absorption XS']

lra_materials['region_3']['Absorption XS'] = [0.008002, 0.08344]

lra_materials['region_3']['Total XS'] = [0.2648, 1.5941]

lra_materials['region_3']['Scattering XS'] = [0.230588, 0.02617, 0.00, 1.510694]

lra_materials['region_3']['Fission XS'] = [0.0019189300, 0.042016461]

lra_materials['region_3']['Nu Fission XS'] = [0.004663, 0.1021]

lra_materials['region_3']['Chi'] = [1.0, 0.0]

lra_materials['region_3']['Diffusion Coefficient'] = [1.259, 0.2091]

lra_materials['region_3']['Buckling'] = [1e-4, 1e-4]

lra_materials['region_3']['Gamma'] = [3.034e-3, 0.0]


###############################################################################
################################   region 4    ################################
###############################################################################

# Create a subdictionary for region 4 material data
lra_materials['region_4'] = {}

lra_materials['region_4']['FunctionalVariables'] = ['Temperature']

lra_materials['region_4']['Temperature'] = 300.0

lra_materials['region_4']['FunctionalTemperature'] = ['Absorption XS']

lra_materials['region_4']['Absorption XS'] = [0.008002, 0.073324]

lra_materials['region_4']['Total XS'] = [0.2648, 1.5941]

lra_materials['region_4']['Scattering XS'] = [0.230588, 0.02617, 0.00, 1.520810]

lra_materials['region_4']['Fission XS'] = [0.0019189300, 0.042016461]

lra_materials['region_4']['Nu Fission XS'] = [0.004663, 0.1021]

lra_materials['region_4']['Chi'] = [1.0, 0.0]

lra_materials['region_4']['Diffusion Coefficient'] = [1.259, 0.2091]

lra_materials['region_4']['Buckling'] = [1e-4, 1e-4]

lra_materials['region_4']['Gamma'] = [3.034e-3, 0.0]


###############################################################################
################################   region 5    ################################
###############################################################################

# Create a subdictionary for region 5 material data
# bladed region
lra_materials['region_5'] = {}

lra_materials['region_5']['FunctionalVariables'] = ['Temperature', 'Time']

lra_materials['region_5']['Temperature'] = 300.0

lra_materials['region_5']['FunctionalTemperature'] = ['Absorption XS']

lra_materials['region_5']['FunctionalTime'] = ['Absorption XS']

lra_materials['region_5']['Time'] = [0.0, 2.0, 3.0]

lra_materials['region_5']['Absorption XS'] = [[0.008002, 0.08344], [0.008002, 0.073324], [0.008002, 0.073324]]

lra_materials['region_5']['Total XS'] = [0.2648, 1.5941]

lra_materials['region_5']['Scattering XS'] = [0.230588, 0.02617, 0.00, 1.510694]

lra_materials['region_5']['Fission XS'] = [0.0019189300, 0.042016461]

lra_materials['region_5']['Nu Fission XS'] = [0.004663, 0.1021]

lra_materials['region_5']['Chi'] = [1.0, 0.0]

lra_materials['region_5']['Diffusion Coefficient'] = [1.259, 0.2091]

lra_materials['region_5']['Buckling'] = [1e-4, 1e-4]

lra_materials['region_5']['Gamma'] = [3.034e-3, 0.0]

###############################################################################
################################   region 6    ################################
###############################################################################

# Create a subdictionary for region 6 material data
lra_materials['region_6'] = {}

lra_materials['region_6']['Absorption XS'] = [0.0006034, 0.01911]

lra_materials['region_6']['Total XS'] = [0.2652, 2.0938]

lra_materials['region_6']['Scattering XS'] = [0.217039, 0.04754, 0.00, 2.074692]

lra_materials['region_6']['Fission XS'] = [0.0, 0.0]

lra_materials['region_6']['Nu Fission XS'] = [0.0, 0.0]

lra_materials['region_6']['Chi'] = [1.0, 0.0]

lra_materials['region_6']['Diffusion Coefficient'] = [1.257, 0.1592]

lra_materials['region_6']['Buckling'] = [1e-4, 1e-4]
