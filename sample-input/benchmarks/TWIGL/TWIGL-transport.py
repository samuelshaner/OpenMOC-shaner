from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

num_threads = options.num_omp_threads
track_spacing = options.track_spacing
num_azim = options.num_azim
tolerance = options.tolerance
max_iters = options.max_iters
relax_factor = options.relax_factor
acceleration = options.acceleration
mesh_level = options.mesh_level
log.setLogLevel('NORMAL')

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
planes[1].setBoundaryType(VACUUM)
planes[2].setBoundaryType(REFLECTIVE)
planes[3].setBoundaryType(VACUUM)

###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(universe=1, material=region1))
cells.append(CellBasic(universe=2, material=region2))
cells.append(CellBasic(universe=3, material=region3))
cells.append(CellFill(universe=0, universe_fill=4))

cells[3].addSurface(halfspace=+1, surface=planes[0])
cells[3].addSurface(halfspace=-1, surface=planes[1])
cells[3].addSurface(halfspace=+1, surface=planes[2])
cells[3].addSurface(halfspace=-1, surface=planes[3])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating LRA lattice...')

lattice = Lattice(id=4, width_x=8.0, width_y=8.0)
lattice.setLatticeCells([[3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                         [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                         [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                         [2, 2, 2, 1, 1, 1, 1, 3, 3, 3],
                         [2, 2, 2, 1, 1, 1, 1, 3, 3, 3],
                         [2, 2, 2, 1, 1, 1, 1, 3, 3, 3],
                         [2, 2, 2, 1, 1, 1, 1, 3, 3, 3],
                         [3, 3, 3, 2, 2, 2, 2, 3, 3, 3],
                         [3, 3, 3, 2, 2, 2, 2, 3, 3, 3],
                         [3, 3, 3, 2, 2, 2, 2, 3, 3, 3]])


###############################################################################
###########################   Creating Cmfd Mesh   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd mesh...')

mesh = Mesh(MOC, acceleration, relax_factor, mesh_level)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry(mesh)
for material in materials.values(): geometry.addMaterial(material)
for cell in cells: geometry.addCell(cell)
geometry.addLattice(lattice)
geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.generateTracks()

###############################################################################
########################   Creating the Cmfd module   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd...')

cmfd = Cmfd(geometry)

###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = CPUSolver(geometry, track_generator, cmfd)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)
solver.printTimerReport()

###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plotMaterials(geometry, gridsize=500)
#plotter.plotCells(geometry, gridsize=500)
#plotter.plotFlatSourceRegions(geometry, gridsize=500)
plotter.plotMeshFluxes(mesh, energy_groups=[1,2])

log.py_printf('TITLE', 'Finished')

