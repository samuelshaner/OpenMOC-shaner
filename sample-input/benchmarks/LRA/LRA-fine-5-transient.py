from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

log.setLogLevel('NORMAL')
dt_cmfd = 1e-2;
num_threads = options.num_omp_threads

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from py...')

materials = materialize.materialize('LRA-materials-transient-ss.py')

region1 = materials['region_1'].getId()
region2 = materials['region_2'].getId()
region3 = materials['region_3'].getId()
region4 = materials['region_4'].getId()
region5 = materials['region_5'].getId()
region6 = materials['region_6'].getId()

###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

planes = []
planes.append(XPlane(x=-82.5))
planes.append(XPlane(x=82.5))
planes.append(YPlane(y=-82.5))
planes.append(YPlane(y=82.5))
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
cells.append(CellBasic(universe=4, material=region4))
cells.append(CellBasic(universe=5, material=region5))
cells.append(CellBasic(universe=6, material=region6))
cells.append(CellFill(universe=21, universe_fill=31))
cells.append(CellFill(universe=22, universe_fill=32))
cells.append(CellFill(universe=23, universe_fill=33))
cells.append(CellFill(universe=24, universe_fill=34))
cells.append(CellFill(universe=25, universe_fill=35))
cells.append(CellFill(universe=26, universe_fill=36))
cells.append(CellFill(universe=0, universe_fill=7))

cells[12].addSurface(halfspace=+1, surface=planes[0])
cells[12].addSurface(halfspace=-1, surface=planes[1])
cells[12].addSurface(halfspace=+1, surface=planes[2])
cells[12].addSurface(halfspace=-1, surface=planes[3])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating LRA lattice...')

assembly1 = Lattice(id=31, width_x=3.0, width_y=3.0)
assembly1.setLatticeCells([[1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 1]])

assembly2 = Lattice(id=32, width_x=3.0, width_y=3.0)
assembly2.setLatticeCells([[2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2],
                           [2, 2, 2, 2, 2]])


assembly3 = Lattice(id=33, width_x=3.0, width_y=3.0)
assembly3.setLatticeCells([[3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3],
                           [3, 3, 3, 3, 3]])


assembly4 = Lattice(id=34, width_x=3.0, width_y=3.0)
assembly4.setLatticeCells([[4, 4, 4, 4, 4],
                           [4, 4, 4, 4, 4],
                           [4, 4, 4, 4, 4],
                           [4, 4, 4, 4, 4],
                           [4, 4, 4, 4, 4]])


assembly5 = Lattice(id=35, width_x=3.0, width_y=3.0)
assembly5.setLatticeCells([[5, 5, 5, 5, 5],
                           [5, 5, 5, 5, 5],
                           [5, 5, 5, 5, 5],
                           [5, 5, 5, 5, 5],
                           [5, 5, 5, 5, 5]])


assembly6 = Lattice(id=36, width_x=3.0, width_y=3.0)
assembly6.setLatticeCells([[6, 6, 6, 6, 6],
                           [6, 6, 6, 6, 6],
                           [6, 6, 6, 6, 6],
                           [6, 6, 6, 6, 6],
                           [6, 6, 6, 6, 6]])

core = Lattice(id=7, width_x=15.0, width_y=15.0)
core.setLatticeCells([[26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26],
                         [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26],
                         [23, 23, 23, 23, 23, 23, 23, 26, 26, 26, 26],
                         [23, 23, 23, 23, 23, 23, 23, 24, 26, 26, 26],
                         [22, 21, 21, 21, 21, 22, 22, 25, 25, 26, 26],
                         [22, 21, 21, 21, 21, 22, 22, 25, 25, 26, 26],
                         [21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26],
                         [21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26],
                         [21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26],
                         [21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26],
                         [22, 21, 21, 21, 21, 22, 22, 23, 23, 26, 26]])


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
geometry.addLattice(assembly4)
geometry.addLattice(assembly5)
geometry.addLattice(assembly6)
geometry.addLattice(core)

geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the Cmfd module   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating cmfd...')

cmfd  = Cmfd(geometry)
cmfd.setNumThreads(num_threads)
cmfd.setOmega(1.75)


tcmfd = Tcmfd(geometry)
tcmfd.setOmega(1.75)
tcmfd.setLambda([0.0654, 1.35])
tcmfd.setBeta([0.0054, 0.001087])
tcmfd.setVelocity([3e7, 3e5])

log.py_printf('NORMAL', 'Creating transient solver...')

transientSolver = TransientSolver(geometry, tcmfd, cmfd)
transientSolver.setKappa(3.204e-11)
transientSolver.setAlpha(3.83e-11)
transientSolver.setNu(2.43)
transientSolver.setDtMOC(dt_cmfd)
transientSolver.setDtCMFD(dt_cmfd)
transientSolver.setStartTime(0.0)
transientSolver.setEndTime(3.0)
transientSolver.setNumDelayGroups(2)
transientSolver.setTransientMethod('ADIABATIC')
transientSolver.setPowerInit(1.e-6)

transientSolver.solveInitialState()

powers = [1e-6]
times  = [0.0]
temps  = [300.0]


for t in range(int(3.0/dt_cmfd)):
   transientSolver.solveOuterStep()
   powers.append(transientSolver.getPower())
   temps.append(transientSolver.getTemp())
   times.append(transientSolver.getTime())

   if (abs(t*dt_cmfd-0.1) < 1.0e-6):
      plotter.plotPrecursors(geometry, 0, gridsize=500)
      plotter.plotPrecursors(geometry, 1, gridsize=500)
      plotter.plotMeshFluxes(mesh, energy_groups=[1,2])
   elif (abs(t*dt_cmfd-1.45) < 1.0e-6):
      plotter.plotPrecursors(geometry, 0, gridsize=500)
      plotter.plotPrecursors(geometry, 1, gridsize=500)
      plotter.plotMeshFluxes(mesh, energy_groups=[1,2])
      

plt.figure()
plt.plot(times, powers)
plt.xlabel('time (s)')
plt.ylabel('Power (W/cc)')
plt.yscale('log')
plt.savefig('plots/powers.png')
   
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax2.plot(times, temps, 'r')
#plt.plot(times, temps)
plt.xlabel('time (s)')
ax2.set_ylabel('Average Temperature (C)')
plt.savefig('plots/temps.png')

###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plotPrecursors(geometry, 0, gridsize=500)
plotter.plotPrecursors(geometry, 1, gridsize=500)
#plotter.plotMaterials(geometry, gridsize=500)
#plotter.plotCells(geometry, gridsize=500)
#plotter.plotFlatSourceRegions(geometry, gridsize=500)
plotter.plotMeshFluxes(mesh, energy_groups=[1,2])

log.py_printf('TITLE', 'Finished')

