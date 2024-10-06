MSM_MD is a Molecular Dynamics code for CMU 12-623 (Molecular Simulation of Materials).

Developed by Alex Echols
Prof. Jerry Wang
Updated 10/7/24

Versions
--------
V1 (9/16/24): Initial version of MSM_MD. Satisfies the following
	Works in non-dimensional LJ units.
	Reads in the initial positions from an input file.
	Initializes the particle velocities to zero.
	Integrates the equations of motion in three dimensions using the velocity Verlet
	   scheme and a non-dimensional time step of 0.002.
	Calculates force and potential energy using the LJ potential.
	Calculates kinetic energy.
	Generates data that can be visualized using either OVITO or VMD.

	USAGE
	-----
	MSM_MD is used through the command line. Currently the following flags allow certain paramaters to be set
		-in : Path of the input file to read from
		-log : Path of the log file (this is different from the dumps). Defaults to ./msm_md.log
		-dump : Path of the dump file. Defaults tp ./msm_md.dump. Currently dump frequency is 1/timestep
		-n : Number of timesteps to run for. Default is 0

	Ex. ./MSM_MD.exe -in ../data/10.txt -dump ../data/10_dump.dump -n 100

	KNOWN ISSUES
	------------
	Logging is currently not done at the correct time for some items, causing them to appear in std out but not the log file

V2 (10/7/24): Second version of MSM_MD. Satisfies the following
	Works in non-dimensional LJ units.
	Reads in the initial positions from an input file.
	Initializes the particle velocities to zero.
	Integrates the equations of motion in three dimensions using the velocity Verlet
	   scheme and a non-dimensional time step of 0.002.
	Calculates force and potential energy using the LJ potential.
	Calculates kinetic energy.
	Generates data that can be visualized using either OVITO or VMD.
	Randomly initializes the particle velocities in a manner that gives the system
		zero total momentum.
	Implements the continuous-force / continuous-energy cutoff scheme with a 
		(dimensionless) cutoff of 2.5.
	Calculates the instantaneous temperature and pressure.
	Applies periodic boundary conditions with the nearest-image convention, so
		that you can model bulk systems

	USAGE
	-----
	MSM_MD is used through the command line. Currently the following flags allow certain paramaters to be set
		-in : Path of the input file to read from
		-log : Path of the log file (this is different from the dumps). Defaults to ./msm_md.log
		-dump : Path of the dump file. Defaults tp ./msm_md.dump. Currently dump frequency is 1/timestep
		-n : Number of timesteps to run for. Default is 0

	Ex. ./MSM_MD.exe -in ../data/10.txt -dump ../data/10_dump.dump -n 100

	KNOWN ISSUES
	------------
	Logging is currently not done at the correct time for some items, causing them to appear in std out but not the log file
	Unit Conversion is not handled properly
	All simulation parameters need to be hardcoded (so tool must be rebuilt on each trial)

TODO
----
- Timer
	- Average time of a timestep/timesteps per unit time
	- Total simulation time

- Fix Logger
- INPUT FILES OMG
- Better unit conversion