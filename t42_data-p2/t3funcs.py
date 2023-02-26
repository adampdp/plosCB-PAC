#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 14:16:58 2021

@author: adam
"""


import numpy as np
import math
import random
from neuron import h
from scipy.signal import chirp



class spkGenGam:
    
    def __init__(self, gshape, gscale, grefract):
        self.gshape = gshape
        self.gscale = gscale
        self.grefract = grefract
    
        
    def call(self):
        spkTimes = (np.cumsum(np.random.gamma(self.gshape, self.gscale, 20000) 
                              + self.grefract)) - 50000 
        spkTimes = spkTimes[spkTimes>0]
        return spkTimes





class spkGenTel3:
    
    def __init__(self, sg1, sg2, onset, uptime, downtime, fulltime, shift, 
                 call_shift, soffsetstart = -10000):
        if isinstance(uptime,(list,np.ndarray)):
            uptime1 = uptime
            downtime1 = downtime
        else:
            uptime1 = [uptime]*10000
            downtime1 = [downtime]*10000
        self.spgt1 = spkGenTel2(sg1, onset, uptime1, downtime1, fulltime, shift, call_shift)
        self.spgt2 = spkGenTel2(sg2, onset + uptime1[0], downtime1, uptime1[1:], fulltime, shift, call_shift)
        # self.spks1 = np.append(self.spks1, spks2)
        
        
        
    def call(self):
        
        return np.append(self.spgt1.call(),self.spgt2.call())
        
class spkGenTel2:
    
    def __init__(self, sg1, onset, uptime, downtime, fulltime, shift, call_shift):
        
        self.sg1 = sg1
        self.onset = onset
        
        if isinstance(uptime,(list,np.ndarray)):
           self.uptime1 = uptime
           self.downtime1 = downtime
        else:
            self.uptime1 = [uptime]*1000
            self.downtime1 = [downtime]*1000
            
        self.fulltime = fulltime - onset
        self.shift = shift
        self.call_shift = call_shift
        self.call_shift_time = 0
        self.soffsetstart = -10*(sum(self.uptime1)/len(self.uptime1) 
                                     + sum(self.downtime1)/len(self.downtime1))
        
        
    def call(self):
        spktimes = self.sg1.call()
        spks = np.array([])
        tcount = 0
    
     
        soffset = self.shift + self.call_shift_time
        self.call_shift_time += self.call_shift
        stime = self.soffsetstart + soffset
        etime = stime + self.uptime1[tcount]
        while etime<self.fulltime:
            if etime > 0:
                if stime > 0:
                   mask = np.logical_and(spktimes>=stime, spktimes<etime)
                else:
                   mask = np.logical_and(spktimes>=0, spktimes<etime)
                spks = np.append(spks,spktimes[mask])
            stime = etime + self.downtime1[tcount]
            tcount += 1
            etime = stime + self.uptime1[tcount]
        else:
            if stime<self.fulltime:
                if stime > 0:
                   mask = np.logical_and(spktimes>=stime, spktimes<self.fulltime)
                else:
                   mask = np.logical_and(spktimes>=0, spktimes<self.fulltime) 
                spks = np.append(spks,spktimes[mask])
        
    
        spks = spks + self.onset
        return spks 



class ArtifPop:
    
    def __init__(self):
        self.artifconnlist = []
        self.cellsList = []
        self.cellcount = 0
        
    
    def create(self, inps_per_cell, pyrrat, cell_label, 
               cellconlist, vrate, spkgenfcn, starttime, fcnduration):
 
   
        for i in range(len(cellconlist)):
            x1 = int(inps_per_cell[i]/pyrrat)
      
            for j1 in range(x1):
                artifspkTimes2 = np.array([])
                for l1 in range(int(pyrrat)):
              
                    
                    v1 = spkgenfcn.call()
                    v1 = v1[v1>0] + starttime
                  
                    artifspkTimes2 =  np.append(artifspkTimes2,v1)
                
                
                artifspkTimes2 = np.array(sorted(artifspkTimes2)) 
                artifspkTimes = artifspkTimes2[artifspkTimes2<fcnduration]
                artifspkTimes3 = []
                for i3 in artifspkTimes:
                    if i3 not in artifspkTimes3:
                        artifspkTimes3.append(i3)
                
                self.cellsList.append({'cellLabel': cell_label + str(self.cellcount), 
                                  'spkTimes': artifspkTimes3})
                for i10 in cellconlist[i]:
                    self.artifconnlist.append([self.cellcount,i10])    
                self.cellcount+=1
              
    def makeListDrive(self,pars):
        
        inp_uptime_start = pars['inp_uptime_start']
        inp_downtime_start = pars['inp_downtime_start']
        call_shift_fac = pars['call_shift_fac'] 
        artifperpyr = pars['artifperpyr']
        pyrpopsize= pars['pyrpopsize']
        gshape = pars['gshape']
        lo_gscale_u = pars['lo_gscale_u']
        hi_gscale_u = pars['hi_gscale_u']
        lo_gscale_d = pars['lo_gscale_d']
        hi_gscale_d = pars['hi_gscale_d']
        grefract = pars['grefract']
        
        uptimelist = [inp_uptime_start]*10000
        downtimelist = [inp_downtime_start]*10000
    
        call_shift = call_shift_fac * (
         (inp_uptime_start + inp_downtime_start)/artifperpyr)
        
        for pcind in range(int(pyrpopsize)):
           
            inpspercell = [int(artifperpyr)]
            
            meanisi = np.random.uniform(lo_gscale_u, hi_gscale_u)
            gscale = meanisi/gshape
            sg1pc = spkGenGam(gshape, gscale, grefract)
            
            meanisi = np.random.uniform(lo_gscale_d, hi_gscale_d)
            gscale = meanisi/gshape
            sg2pc = spkGenGam(gshape, gscale, grefract)
            
            
         
            inp_onset = pcind*pars['field_delay'] 
            
            sgstep = spkGenTel3(sg1pc, sg2pc, pars['artif_starttime'], 
                                    uptimelist, downtimelist, 
                                    pars['artif_duration'], 
                                    inp_onset, call_shift)
            
            sgstep.soffsetstart = -10*(sum(uptimelist)/len(uptimelist) 
                                                 + sum(downtimelist)/len(downtimelist))
            
        
        
            cellconlist = [[pcind]]
            self.create(inpspercell, pars['artifpyrat'], pars['namestr'], 
                                  cellconlist, 0, sgstep, 0, pars['artif_duration'])     


