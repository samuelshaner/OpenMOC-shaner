##
# @file materialize.py
# @package openmoc.materialize
# @brief The materialize module provides utility functions to read and write
#        multi-group materials cross-section data from HDF5 binary file format.
# @author William Boyd (wboyd@mit.edu)
# @date April 23, 2013

from openmoc import *
from log import *


##
# @brief
# @param filename
# @return a list of materials for OpenMOC
def materialize(filename):

    xs_types = ['Total XS', 'Absorption XS', 'Scattering XS', \
                    'Fission XS', 'Nu Fission XS', 'Chi',\
                    'Diffusion Coefficient']
    materials = {}

    # Check that the filename is a string
    if not isinstance(filename, str):
        py_printf('ERROR', 'Unable to materialize using filename %s ' + \
                      'since it is not a string', str(filename))


    ############################################################################
    #                               HDF5 DATA FILES
    ############################################################################
    if filename.endswith('.hdf5'):

        import h5py
        import numpy as np

        # Create a h5py file handle for the file
        f = h5py.File(filename)

        # Check that the file has an 'energy groups' attribute
        if not 'Energy Groups' in f.attrs:
            py_printf('ERROR', 'Unable to materialize file %s since it ' + \
                          'does not contain an \'Energy Groups\' attribute', \
                          filename)
    
        num_groups = f.attrs['Energy Groups']
        # Check that the number of energy groups is an integer
        if not isinstance(num_groups, int):
            py_printf('ERROR', 'Unable to materialize file %s since the ' + \
                          'number of energy groups %s is not an integer', \
                          filename, str(num_groups))

        material_names = list(f)

        # Loop over each material and 
        for name in material_names:

            py_printf('INFO', 'Importing material %s', str(name))

            new_material = Material(material_id())
            new_material.setNumEnergyGroups(int(num_groups))

            # Retrieve and load the cross-section data into the material object


            if 'Total XS' in f[name]:
                new_material.setSigmaT(f[name]['Total XS'][...])

            if 'Diffusion Coefficient' in f[name]:
                new_material.setDifCoef(f[name]['Diffusion Coefficient'][...])

            if 'Buckling' in f[name]:
                new_material.setBuckling(f[name]['Buckling'][...])
            
            if 'Absorption XS' in f[name]:
                new_material.setSigmaA(f[name]['Absorption XS'][...])
            
            if 'Scattering XS' in f[name]:
                new_material.setSigmaS(f[name]['Scattering XS'][...])
            
            if 'Fission XS' in f[name]:
                new_material.setSigmaF(f[name]['Fission XS'][...])
            
            if 'Nu Fission XS' in f[name]:
                new_material.setNuSigmaF(f[name]['Nu Fission XS'][...])
            
            if 'Chi' in f[name]:
                new_material.setChi(f[name]['Chi'][...])

            # Add this material to the list
            materials[name] = new_material


    ############################################################################
    #                      PYTHON DICTIONARY DATA FILES
    ############################################################################
    elif filename.endswith('.py'):

        import imp
        data = imp.load_source(filename, filename).dataset

        # Check that the file has an 'energy groups' attribute
        if not 'Energy Groups' in data.keys():
            py_printf('ERROR', 'Unable to materialize file %s since it ' + \
                      'does not contain an \'Energy Groups\' attribute', \
                      filename)
    
        num_groups = data['Energy Groups']

        # Check that the number of energy groups is an integer
        if not isinstance(num_groups, int):
            py_printf('ERROR', 'Unable to materialize file %s since the ' + \
                          'number of energy groups %s is not an integer', \
                          filename, str(num_groups))

        data = data['Materials']
        material_names = data.keys()

        # Loop over each material and 
        for name in material_names:

            py_printf('INFO', 'Importing material %s', str(name))

            # create material object
            if 'FunctionalVariables' in data[name].keys():
                new_material = FunctionalMaterial(material_id())
                
                if 'Time' in data[name]['FunctionalVariables']:
                    new_material.setNumEnergyGroups(int(num_groups), int(len(data[name]['Time'])))
                else:
                    new_material.setNumEnergyGroups(int(num_groups), 1)
            else:
                new_material = Material(material_id())
                new_material.setNumEnergyGroups(int(num_groups))
            
            if 'Diffusion Coefficient' in data[name].keys():
                new_material.setDifCoef(data[name]['Diffusion Coefficient'])

            if 'Buckling' in data[name].keys():
                new_material.setBuckling(data[name]['Buckling'])
            
            if 'FunctionalVariables' in data[name].keys():
                if 'Temperature' in data[name]['FunctionalVariables']:
                    if 'Absorption XS' in data[name]['FunctionalTemperature']:
                        new_material.sigmaAFuncTemp(True)

            # set time
            if 'FunctionalVariables' in data[name].keys():
                if 'Time' in data[name]['FunctionalVariables']:
                    if 'Time' in data[name].keys():
                        new_material.setTime(data[name]['Time'])
                    else:
                        py_printf('ERROR', 'If time is a functional variable,' + \
                                      'you must pass in an array of time values!')

            # set sigmaA
            if 'FunctionalVariables' in data[name].keys():
                if 'Time' in data[name]['FunctionalVariables']:
                    if 'Absorption XS' in data[name]['FunctionalTime']:
                        new_material.sigmaAFuncTime(True)
                        new_material.setSigmaATime(data[name]['Absorption XS'])
                    else:
                        new_material.setSigmaA(data[name]['Absorption XS'])
                else:
                    new_material.setSigmaA(data[name]['Absorption XS'])
            else:
                new_material.setSigmaA(data[name]['Absorption XS'])

            # set sigmaS
            if 'FunctionalVariables' in data[name].keys():
                if 'Time' in data[name]['FunctionalVariables']:
                    if 'Scattering XS' in data[name]['FunctionalTime']:
                        new_material.sigmaSFuncTime(True)
                        new_material.setSigmaSTime(data[name]['Scattering XS'])
                    else:
                        new_material.setSigmaS(data[name]['Scattering XS'])
                else:
                    new_material.setSigmaS(data[name]['Scattering XS'])
            else:
                new_material.setSigmaS(data[name]['Scattering XS'])

            if 'Fission XS' in data[name].keys():
                new_material.setSigmaF(data[name]['Fission XS'])

            if 'Nu Fission XS' in data[name].keys():
                new_material.setNuSigmaF(data[name]['Nu Fission XS'])

            if 'Chi' in data[name].keys():
                new_material.setChi(data[name]['Chi'])
            
            if 'Gamma' in data[name].keys():
                new_material.setGamma(data[name]['Gamma'])

            # compute sigma_t
            new_material.computeSigmaT()
              
            # Add this material to the list
            materials[name] = new_material


    ############################################################################
    #                      UNSUPPORTED DATA FILE TYPES
    ############################################################################
    else:
        py_printf('ERROR', 'Unable to materialize using filename %s ' + \
                      'since it has an unkown extension. Supported ' + \
                      'extension types are .hdf5 and .py', filename)



    # Return the list of materials
    return materials
