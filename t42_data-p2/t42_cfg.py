#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:41:05 2021

@author: adam
"""

from netpyne import specs

## Population parameters
cfg = specs.SimConfig()					# object of class SimConfig to store simulation configuration

cfg.duration = 1000
cfg.starttime = 0
cfg.seedval = 42

cfg.pyrpopsize = 20
cfg.pcscalenum = 20 
cfg.pvbcpopsize = 2
cfg.pvscalenum = 1
cfg.olmpopsize = 2
cfg.olmscalenum = 1 

cfg.connectPC2PC = False # False
cfg.connectPVBC2PVBC = False # False
cfg.connectPVBCPC = True # False 
cfg.connectOLMPC = True # False
cfg.connectPCPVBC = True # False
cfg.connectPCOLM = True # False

# =============================================================================
# cfg.connectPC2PC =  False
# cfg.connectPVBC2PVBC =  False
# cfg.connectPVBCPC = False 
# cfg.connectOLMPC =  False
# cfg.connectPCPVBC =  False
# cfg.connectPCOLM =  False
# =============================================================================

cfg.PVBCsomaDist = 50
cfg.OLMsomaDist = 250


cfg.pc_olm_use  = 0.07 
cfg.olmdepfact = 38
cfg.olmfacfact = 470
cfg.pvbcdep  = 110
cfg.pvbcfac = 0 
cfg.olm_pc_gaba_tau = 11.8
cfg.olm2pcDep = 1770
cfg.olm2pcFac = 6
cfg.pv_pc_gaba_tau_fact = 1
cfg.pvbc2pcDep = 965
cfg.pvbc2pcFac = 8.6
cfg.pc_olm_sec = 'basal'

cfg.pc_olm_conprob = 0.35
cfg.pc_olm_synfact = 1
cfg.pc_pc_conprob = 0.35
cfg.pc_pc_synfact = 1
cfg.pc_pv_conprob = 0.35
cfg.pc_pv_synfact  = 1
cfg.olm_pc_conprob = 0.4
cfg.olm_pc_synfact  = 1
cfg.pv_pc_conprob = 0.25
cfg.pv_pc_synfact = 1
cfg.pvbc_pvbc_conprob = 0.35
cfg.pvbc_pvbc_synfact = 1

cfg.pc_olm_hibound = 0.7 
cfg.pc_olm_lowbound = 0.5 
cfg.pc_olm_wei = 0.5


cfg.pc_pv_wei = 1

cfg.olm_pc_wei = 0.25
cfg.olm_pc_lobound = 4.1
cfg.olm_pc_hibound = 5.5



cfg.pv_pc_wei = 1
    


############################

cfg.artifperpyr = 100

cfg.artifpyrpars = {}
cfg.artifpyrpars['doartif'] = True
cfg.artifpyrpars['namestr'] = '_pyrart'
cfg.artifpyrpars['artif_starttime'] = 0
cfg.artifpyrpars['artif_duration'] = cfg.duration
cfg.artifpyrpars['pyrpopsize'] = cfg.pyrpopsize

cfg.artifpyrpars['artifpyrat'] =  int(cfg.artifperpyr/2)        
cfg.artifpyrpars['inp_uptime_start'] = 60000
cfg.artifpyrpars['inp_downtime_start'] = 60000
cfg.artifpyrpars['artifperpyr'] = cfg.artifperpyr

cfg.artifpyrpars['gshape'] = 1

meanCA3isi = 688
gscale = meanCA3isi/(cfg.artifpyrpars['gshape'])

cfg.artifpyrpars['lo_gscale_u'] = gscale
cfg.artifpyrpars['hi_gscale_u'] = gscale
cfg.artifpyrpars['lo_gscale_d'] = gscale
cfg.artifpyrpars['hi_gscale_d'] = gscale

cfg.artifpyrpars['scalesep'] = 0
cfg.artifpyrpars['grefract'] = 0
cfg.artifpyrpars['field_delay'] = 0
cfg.artifpyrpars['call_shift_fac'] = 0

cfg.artifpyrpars['npartifwei_low'] = 0.55
cfg.artifpyrpars['npartifwei_hi'] = 0.65
cfg.artifpyrpars['artifsynmech'] = 'PC-PCzero'
cfg.artifpyrpars['artifsynfact'] = 6




#############################

cfg.cvode_active = True
cfg.dt = 0.025
cfg.hParams = {'v_init': -65, 'celsius': 34} #, 'clamp_resist': 0.001}
cfg.verbose = False

cfg.distributeSynsUniformly = False
cfg.connRandomSecFromList = True

cfg.recordStep = 1 			# Step size in ms to save data (eg. V traces, LFP, etc)
cfg.savePickle = False 		# Save params, network and sim output to pickle file
cfg.saveFileStep = 1000 # step size in ms to save data to disk



cfg.analysis['plotRaster'] = { #'saveData': 'raster_data_' + str(runnum) + '.json', 
                              'include': ['artif_pyr', 'PYR_pop', 'PVBC_pop', 'OLM_pop'],
                                    'marker': 'o',
                                    'saveFig': True, 
                                    'showFig': False, 
                                    'markerSize': 6}


cfg.recordTraces['V_soma'] = {'sec':'soma_0','loc':0.5,'var':'v'}

 
    
cfg.analysis['plotTraces'] = { 
                                        'include': list(range(4)), 
                                      
                                        'saveFig': True, 
                                        'showFig': False,
                                 
                                        } # Plot cell traces 

    
cfg.saveJson = True
cfg.seeds = {'conn': cfg.seedval + 7515, 'stim': cfg.seedval + 84331, 'loc': cfg.seedval + 943}