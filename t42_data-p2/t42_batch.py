#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:34:42 2021

@author: adam
"""
from netpyne import specs
from netpyne.batch import Batch

#froam neuron import h
#h.nrnmpi_init()

def batchTauWeight():

    params = specs.ODict()
    
    seedbase = 576667
    
    #params['artifperpyr'] = [130 + x*4 for x in range(20)] 
    params['artifperpyr'] = [140]
    params['seedval'] = list(range(0 + seedbase, 100 + seedbase, 100)) 
    
    b = Batch(params=params, cfgFile='t42_cfg.py', netParamsFile='t42_netParams.py',)

 
    b.batchLabel = 't42'
    b.saveFolder = 't42_data'
    b.method = 'grid'
    
    doslurm = False
    if doslurm:
    
        b.runCfg = {    'type': 'hpc_slurm',
                        #'type': 'mpi_bulletin',
                            'mpiCommand': 'srun',
                            'custom': '#SBATCH --constraint=mc\n#SBATCH --partition=normal',
                            'allocation': 'ich011',
                            'nodes': 2,
                            'coresPerNode': 35,
                            'script': 't42_init.py',
                            'walltime': '24:00:00',
                            'skip': True}
    else:
   

    
        b.runCfg = {  'type': 'mpi_direct',
                'cores': 10,
                'mpiCommand': 'mpiexec --use-hwthread-cpus',
                        'script': 't42_init.py',
                        'skip': True}
    
    

    # Run batch simulations
    b.run()
    #h.quit()
# Main code
if __name__ == '__main__':
        batchTauWeight()
        import sys
        sys.exit()
