from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

log.setLogLevel('NORMAL')
dt_moc = 1e-2
dt_cmfd = 1e-3
num_threads = options.num_omp_threads

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

mesh = Mesh(DIFFUSION)

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
########################   Creating the Cmfd module   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd...')

cmfd = Cmfd(geometry, 1.e-10)
cmfd.setNumThreads(num_threads)
cmfd.setOmega(1.5)

tcmfd = Tcmfd(geometry, 1.e-10)
tcmfd.setOmega(1.5)
tcmfd.setLambda([0.08])
tcmfd.setBeta([0.0075])
tcmfd.setVelocity([1.e7, 2.e5])

log.py_printf('NORMAL', 'Creating transient solver...')

transientSolver = TransientSolver(geometry, tcmfd, cmfd)
transientSolver.setKappa(1.0)
transientSolver.setAlpha(0.0)
transientSolver.setNu(2.43)
transientSolver.setDtMOC(dt_cmfd)
transientSolver.setDtCMFD(dt_cmfd)
transientSolver.setStartTime(0.0)
transientSolver.setEndTime(0.5)
transientSolver.setNumDelayGroups(1)
transientSolver.setTransientMethod('MAF')
transientSolver.setPowerInit(1.0e-6)

transientSolver.solveInitialState()

for t in range(int(0.5/dt_cmfd)):
   transientSolver.solveOuterStep()

###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plotMaterials(geometry, gridsize=500)
#plotter.plotCells(geometry, gridsize=500)
#plotter.plotFlatSourceRegions(geometry, gridsize=500)
plotter.plotMeshFluxes(mesh, energy_groups=[1,2])

log.py_printf('TITLE', 'Finished')

