from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

tolerance = 1E-10
log.setLogLevel('INFO')

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data...')

materials = materialize.materialize('TWIGL-materials.py')

region1 = materials['region_1'].getId()
region2 = materials['region_2'].getId()
region3 = materials['region_3'].getId()

###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

planes = []
planes.append(XPlane(x=-40.0))
planes.append(XPlane(x=40.0))
planes.append(YPlane(y=-40.0))
planes.append(YPlane(y=40.0))
planes[0].setBoundaryType(REFLECTIVE)
planes[1].setBoundaryType(ZERO_FLUX)
planes[2].setBoundaryType(REFLECTIVE)
planes[3].setBoundaryType(ZERO_FLUX)

###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(universe=1, material=region1))
cells.append(CellBasic(universe=2, material=region2))
cells.append(CellBasic(universe=3, material=region3))
cells.append(CellFill(universe=21, universe_fill=31))
cells.append(CellFill(universe=22, universe_fill=32))
cells.append(CellFill(universe=23, universe_fill=33))
cells.append(CellFill(universe=0, universe_fill=4))

cells[3].addSurface(halfspace=+1, surface=planes[0])
cells[3].addSurface(halfspace=-1, surface=planes[1])
cells[3].addSurface(halfspace=+1, surface=planes[2])
cells[3].addSurface(halfspace=-1, surface=planes[3])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating LRA lattice...')

assembly1 = Lattice(id=31, width_x=1.0, width_y=1.0)
assembly1.setLatticeCells([[1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1, 1, 1, 1]])

assembly2 = Lattice(id=32, width_x=1.0, width_y=1.0)
assembly2.setLatticeCells([[2, 2, 2, 2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2, 2, 2, 2]])

assembly3 = Lattice(id=33, width_x=1.0, width_y=1.0)
assembly3.setLatticeCells([[3, 3, 3, 3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3, 3, 3, 3]])


lattice = Lattice(id=4, width_x=8.0, width_y=8.0)
lattice.setLatticeCells([[23, 23, 23, 23, 23, 23, 23, 23, 23, 23],
                         [23, 23, 23, 23, 23, 23, 23, 23, 23, 23],
                         [23, 23, 23, 23, 23, 23, 23, 23, 23, 23],
                         [22, 22, 22, 21, 21, 21, 21, 23, 23, 23],
                         [22, 22, 22, 21, 21, 21, 21, 23, 23, 23],
                         [22, 22, 22, 21, 21, 21, 21, 23, 23, 23],
                         [22, 22, 22, 21, 21, 21, 21, 23, 23, 23],
                         [23, 23, 23, 22, 22, 22, 22, 23, 23, 23],
                         [23, 23, 23, 22, 22, 22, 22, 23, 23, 23],
                         [23, 23, 23, 22, 22, 22, 22, 23, 23, 23]])


###############################################################################
###########################   Creating Cmfd Mesh   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd mesh...')

mesh = Mesh(DIFFUSION)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry(mesh)
for material in materials.values(): geometry.addMaterial(material)
for cell in cells: geometry.addCell(cell)
geometry.addLattice(assembly1)
geometry.addLattice(assembly2)
geometry.addLattice(assembly3)
geometry.addLattice(lattice)
geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the Cmfd module   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd...')

cmfd = Cmfd(geometry)
cmfd.setOmega(1.9)
cmfd.computeKeff()

log.py_printf('NORMAL', 'k_eff = %f', cmfd.getKeff())

###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plotMaterials(geometry, gridsize=500)
#plotter.plotCells(geometry, gridsize=500)
#plotter.plotFlatSourceRegions(geometry, gridsize=500)
plotter.plotMeshFluxes(mesh, energy_groups=[1,2])

log.py_printf('TITLE', 'Finished')

