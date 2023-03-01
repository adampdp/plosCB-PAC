# plosCB-PAC


Model files in the NetPyNE simulation environment for the paper:

Theta-gamma phase amplitude coupling in a hippocampal CA1 microcircuit

by Adam Ponzi, Salvador Dura-Bernal, and Michele Migliore

Usage:

Download the files. 

Compile the mod files on the command line with 'nrnivmodl mechanisms'.

Install the NetPyNE simulation environment (http://netpyne.org/install.html)

The directory 't42_data_p2' contains configuration files for the large network simulations. Note results in the paper for simulations with 120 or 480 PYR will require a supercomputer.

The directory 'scanzres' contains configuration files for the small simulations to recreate the results in Fig1.

To run the code using mpi-exec: (1) Remove or rename the directory t42_data (which is subsequently recreated when the code is run). (2) Set the numcores parameter in t42_batch.py to the required number of cores. (3) On the command line : 'nrniv -python t42_batch.py'  (4) Result files are in the new directory t42_data.

To run the code directly: (1) On the command line: 'nrniv -python t42_init.py'. In this case new directories are not created. Results files are in the working directory.

On completion of the small simulation the file maketrace.py in the 'scanzres' directory can be run to recreate the results in Fig1C. Set the flag asBatch = True if the simulation was run using mpi-exec.

On completion of large simulations the file analysis2.py in the 'Fig2data' can be run to recreate single network simulation statistical analysis as those in Fig2 of the paper. This directory also contains the simulation data for the simulation shown in Fig2 of the paper. 

Questions on how to use this model should be directed to adam.ponzi@ibf.cnr.it

