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

dataset['Energy Groups'] = 2
dataset['Materials'] = {}

twigl_materials = dataset['Materials']

###############################################################################
################################   region 1    ################################
###############################################################################

# Create a subdictionary for region 1 material data
twigl_materials['region_1'] = {}

twigl_materials['region_1']['Time'] = [0.0, 0.2, 0.5]

twigl_materials['region_1']['FunctionalVariables'] = ['Time']

twigl_materials['region_1']['FunctionalTime'] = ['Absorption XS']

twigl_materials['region_1']['Absorption XS'] = [[0.01, 0.15], [0.01, 0.1465], [0.01, 0.1465]]

twigl_materials['region_1']['Scattering XS'] = [0.2181, 0.01, 0.00, 0.68333]

twigl_materials['region_1']['Fission XS'] = [0.002, 0.05]

twigl_materials['region_1']['Nu Fission XS'] = [0.007, 0.2]

twigl_materials['region_1']['Chi'] = [1.0, 0.0]

twigl_materials['region_1']['Diffusion Coefficient'] = [1.4, 0.4]

twigl_materials['region_1']['Gamma'] = [0.0, 0.0]

###############################################################################
################################   region 2    ################################
###############################################################################

# Create a subdictionary for region 2 material data
twigl_materials['region_2'] = {}

twigl_materials['region_2']['FunctionalVariables'] = ['Temperature']

twigl_materials['region_2']['Temperature'] = 300.0

twigl_materials['region_2']['FunctionalTemperature'] = ['Absorption XS']

twigl_materials['region_2']['Absorption XS'] = [0.01, 0.15]

twigl_materials['region_2']['Scattering XS'] = [0.2181, 0.01, 0.00, 0.68333]

twigl_materials['region_2']['Fission XS'] = [0.002, 0.05]

twigl_materials['region_2']['Nu Fission XS'] = [0.007, 0.2]

twigl_materials['region_2']['Chi'] = [1.0, 0.0]

twigl_materials['region_2']['Diffusion Coefficient'] = [1.4, 0.4]

twigl_materials['region_2']['Gamma'] = [0.0, 0.0]

###############################################################################
################################   region 3    ################################
###############################################################################

# Create a subdictionary for region 3 material data
twigl_materials['region_3'] = {}

twigl_materials['region_3']['FunctionalVariables'] = ['Temperature']

twigl_materials['region_3']['Temperature'] = 300.0

twigl_materials['region_3']['FunctionalTemperature'] = ['Absorption XS']

twigl_materials['region_3']['Absorption XS'] = [0.008, 0.05]

twigl_materials['region_3']['Scattering XS'] = [0.23841, 0.01, 0.00, 0.61667]

twigl_materials['region_3']['Fission XS'] = [0.002, 0.05]

twigl_materials['region_3']['Nu Fission XS'] = [0.003, 0.06]

twigl_materials['region_3']['Chi'] = [1.0, 0.0]

twigl_materials['region_3']['Diffusion Coefficient'] = [1.3, 0.5]

twigl_materials['region_3']['Gamma'] = [0.0, 0.0]
