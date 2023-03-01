#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:49:12 2023

@author: adam
"""

import sys

#import myfuncs as mf
#from scipy import signal
import json
import matplotlib.pyplot as plt
#import matplotlib.ticker as mticker
#from matplotlib.ticker import StrMethodFormatter

import numpy as np
from sklearn.linear_model import LinearRegression
#from scipy.sparse.linalg import eigs
from json.decoder import JSONDecodeError

from pactools import Comodulogram #, REFERENCES
from scipy import signal
#import pickle
from scipy.signal import butter, lfilter
#import itertools



plt.rcParams.update({'font.sans-serif':'Helvetica'})

def MovAv(xmat, length, axis = 0):
    
    x = [slice(None)]*xmat.ndim
  
    s1 = list(xmat.shape)
    s1[axis] = s1[axis] - length + 1
    s1 = tuple(s1)
    
    xmat2 = np.zeros(s1)
  
    x2 = [slice(None)]*xmat2.ndim
    
    for i in list(range(xmat.shape[axis] - length + 1)):
        x[axis] = slice(i,i+length,1)
    
        x2[axis] = i
        xmat2[tuple(x2)] = np.mean(xmat[tuple(x)], axis = axis)
    return xmat2


class welchSig:
    
    def __init__(self, vals, fs, psdnperseg, psdnoverlap, avlen = 0):

          self.Fxx, self.Pxx = signal.welch(vals, fs=fs, 
                           window='hanning', 
                           nperseg=psdnperseg, noverlap=psdnoverlap)
         
          if avlen>0:
              self.Fxx  = MovAv(np.array(self.Fxx),avlen)
              self.Pxx  = MovAv(np.array(self.Pxx),avlen)
             
           
    def getPow(self, f1):
        
        sfmask = np.logical_and(self.Fxx>=f1[0], self.Fxx<f1[1])    
        return self.Pxx[sfmask].mean(axis = 0)



class mRates:

    def __init__(self, starttime, stoptime, binsize, belowStop = True):
        
        self.binsize = binsize
        self.starttime = starttime
        self.stoptime = stoptime
        self.binlist = np.arange(self.starttime, self.stoptime, self.binsize) + self.binsize/2
        self.numbins = self.binlist.shape[0] #int((stoptime - starttime)/binsize)
        self.rates = np.zeros(self.numbins)
        self.belowStop = belowStop
        
    def addVal(self, val):
        
        if hasattr(val, '__iter__'):
            for x in val:
               self.addVal(x)
        else:       
            if (val<self.stoptime):
               ind = int((val - self.starttime)/self.binsize)
               self.rates[ind] = self.rates[ind] + 1
            elif self.belowStop:   
               sys.exit('mRates')

    
    def binList(self):
        return self.binlist
    
# =============================================================================
#     def addArrayVal(self, val):
#         for x in val:
#             self.addVal(x)
# =============================================================================
     
    def ratesNormed(self):        
        return self.rates/self.binsize
    
    def reset(self):
        self.rates = self.rates * 0


class getRaster:
    
    def __init__(self,data,popstr):
      
        if isinstance(popstr, str):
           cellinds = data['net']['pops'][popstr]['cellGids']
        elif isinstance(popstr, list):
           cellinds = popstr 
        else:
           cellinds = [popstr] 
           
        rasttimes = data['simData']['spkt']
        rastinds = data['simData']['spkid']

     
        self.times = np.array([rasttimes[ind] for ind, val in enumerate(rastinds) if val in cellinds])
        self.vals =  np.array([val for val in rastinds if val in cellinds])
        
        
        
class getRates:
    
    def __init__(self,data,data_data,popstr,ratebinsize, cellind = None):
   
        if popstr is None and cellind is not None:
           cellinds = [cellind]
        else: 
           cellinds = data_data['net']['pops'][popstr]['cellGids']                   
        rasttimes = data_data['simData']['spkt']
        rastinds = data_data['simData']['spkid']
        duration = data['simConfig']['duration']
        if 'starttime' in data['simConfig']:
            starttime_orig = data['simConfig']['starttime']
        else:
            starttime_orig = 0
        starttime = starttime_orig
        #tottimelen = (duration - starttime)/1000
        
        mr1 = mRates(starttime,duration,ratebinsize)
        
        nprastimes = np.array(rasttimes)
        nprastinds = np.array(rastinds)
        nprastinds = nprastinds[nprastimes>starttime]
        nprastimes = nprastimes[nprastimes>starttime]
        
        for tcellind1 in cellinds:
       
            nprastimes1 = nprastimes[nprastinds == tcellind1]
            mr1.addVal(nprastimes1)     
     
        self.vals = mr1.ratesNormed()*1000.0/len(cellinds)
        self.times = mr1.binList()



def removedoubzero(x,y,smallval = 0.00000001):
    
    c1 = (x > smallval) & (y > smallval)
    return x[c1], y[c1]


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

class pacMethod:
    
   def __init__(self, signal, low_fr_lo, low_fr_hi, fs, outdirbase, pngname, method,  
                 low_fq_width = 1):
   
        self.cmap = 'inferno' 
        low_fq_range = np.linspace(low_fr_lo, low_fr_hi, 50)
       
        fig, axs = plt.subplots()
      
        self.estimator = Comodulogram(fs=fs, low_fq_range=low_fq_range,
                                     low_fq_width=low_fq_width, method=method,
                                     progress_bar=False)
        self.estimator.fit(signal)
        self.estimator.plot(titles=None, axs=[axs], cmap = self.cmap)
        axs.set_title(pngname)
        plt.tight_layout()
        plt.savefig(outdirbase + pngname + '.png')
       
        
        #plt.close()

class avspctlist:
    
     def __init__(self):
          self.sxxlist = []
        
    
     def add(self, tlist, st, sxx, lolen = 0, hilen = 1000000):
         
         for tind, tval in enumerate(tlist):
             if tind + 1 < len(tlist):
                 tlen = tlist[tind+1] - tval
                 if tlen>=lolen and tlen<hilen:
                     mask = np.logical_and(st>=tval, st<tlist[tind+1])         
                     self.sxxlist.append(sxx[:,mask])
     
     def setMaxMin(self):       
            
        
         self.maxlen = 0
         self.minlen = 100000
         for sval in self.sxxlist:
             if sval.shape[1]>self.maxlen:
                self.maxlen = sval.shape[1]
             if sval.shape[1]<self.minlen:
                self.minlen = sval.shape[1]   
         
     
         if self.maxlen == 0:
             sys.exit('maxlen avspctlist')
     
     
     def plotMinScale(self, outdir, nstr, label, inc = 1, start = 0):
         
         
         self.setMaxMin()
         self.totarr = np.zeros((self.sxxlist[0].shape[0], self.minlen))
         totadd = np.zeros_like(self.totarr)
         sumarr = 0
         for sval in self.sxxlist:
            ind = 0 
            for i in sval:
                totadd[ind] = moving_average(i, sval.shape[1] - self.minlen + 1)
                ind += 1
            self.totarr = self.totarr + totadd
            sumarr = sumarr + 1
         self.totarr = self.totarr/sumarr
         
         plotFig(np.arange(start,self.totarr.shape[1]*inc,inc), 
                    self.totarr[0,:], outdir, 'theta_averaged' + nstr, label = label) 
         
def maxlist(x, t):
    difflist = np.sign( x[1:] - x[:-1] )
    diffprods = difflist[:-1]*difflist[1:]
    return [t[ind + 1] for ind, val in enumerate(diffprods) if val < -0.9]         
  
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a      
            
            
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y            
            

class butterPass:

    def __init__(self,thetapeakfreq,fs,g3vals,g3times):
        
        oddbutter = 1
        butter_order = 5
        thetaoff = 2
        fperiod = 1/thetapeakfreq
        lowcut = max(1,thetapeakfreq - thetaoff)
        highcut = thetapeakfreq + thetaoff
                
                
        butter_y = butter_bandpass_filter(g3vals, lowcut, highcut, fs, order = butter_order)
                            
        butter_max = maxlist(butter_y, g3times)
        self.btl = np.array(butter_max)[oddbutter::2]/1000.0
        
        avspctoff = 0.1   
        self.avspctlolen = (1-avspctoff)*fperiod #- 30/1000
        self.avspcthilen = (1+avspctoff)*fperiod #+ 30/1000       
        self.btlarrmat = np.zeros((1,butter_y.shape[0]))
 
   
    def call(self,g3vals,g3times,popstr,outdir,nstr):
    
        self.btlarrmat[0,:] = g3vals 
                            
                            
        btlspct =  avspctlist()
       
       
        btlspct.add(self.btl, g3times/1000.0, self.btlarrmat,self.avspctlolen,self.avspcthilen)                 
        btlspct.plotMinScale(outdir, nstr, label = popstr, inc = ratebinsize/1000.0)

            


class powTSCalc:
    
    def __init__(self, ws1, powranges, spctthetapeakspread, spctgammapeakspread):
        
      
        wmmask = np.logical_and(ws1.Fxx>=powranges[1][0], ws1.Fxx<powranges[1][1])
        logammapeakfreq = ws1.Fxx[wmmask][np.argmax(ws1.Pxx[wmmask])]          
                    
        wmmask = np.logical_and(ws1.Fxx>=powranges[0][0], ws1.Fxx<powranges[0][1])
        thetapeakfreq = ws1.Fxx[wmmask][np.argmax(ws1.Pxx[wmmask])] 
        
        
        self.spctpowrange = []
        self.spctpowrange.append([thetapeakfreq - spctthetapeakspread, 
                        thetapeakfreq + spctthetapeakspread])
        self.spctpowrange.append([logammapeakfreq - spctgammapeakspread, 
                        logammapeakfreq + spctgammapeakspread])
        
        self.spctpowlabel = ['theta','low gamma']
        
        
    def call(self,s122,outdirbase,nstr):     
       
        powts = []
        for powrange in self.spctpowrange:                    
            powts.append(s122.getPow(powrange[0],powrange[1]))
           
       
        powts[0], powts[1] = removedoubzero(powts[0],powts[1])
        
        powts[0] = (powts[0] - np.mean(powts[0]))/np.std(powts[0])
        powts[1] = (powts[1] - np.mean(powts[1]))/np.std(powts[1])
        
        
        plt.figure('powrange')    
        for powt1,label in zip(powts,self.spctpowlabel):
            plt.plot(s122.st, powt1, label=label)
            
      
        plt.title('power variation' + nstr)
        plt.xlabel('time')
        plt.ylabel('standardized power')
        plt.legend()
        plt.savefig(outdirbase +'powrange' + nstr + '.png')    
        
        plt.figure('powrange xy')
        plt.plot(powts[0], powts[1])
        plt.title('power variation x-y' + nstr)
        plt.xlabel('standardized power ' + self.spctpowlabel[0])
        plt.ylabel('standardized power ' + self.spctpowlabel[1])
        plt.legend()
        plt.savefig(outdirbase +'powrange_xy' + nstr + '.png')    
        
        spctcc = np.corrcoef(np.stack((powts[0],powts[1])))  
        print('spctcc: ', spctcc)            
        regr = LinearRegression().fit(np.array(powts[0]).reshape(-1,1),
                                      np.array(powts[1]))            
        regscr = regr.score(powts[0].reshape(-1,1),powts[1])
        print('regress: ', regscr, regr.coef_)    
                
                 
def plotFig(xvals, yvals, outdir, pngname, label):
    

    plt.figure(pngname)
    plt.plot(xvals, yvals, label = label)
    plt.xlabel('time (secs)')
    plt.ylabel('firing rate (Hz)')
    plt.title(pngname)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir + pngname + '.png', transparent=True)
   # plt.close()


def plotMesh(t1, f1, sxx, outdir, pngname, cmap = 'inferno' , 
             xlabel = None, ylabel = None):
    
   
    plt.figure(pngname)
    plt.pcolormesh(t1, f1, sxx, shading='gouraud', cmap = cmap )
    plt.colorbar()
    plt.xlabel('time (secs)')
    plt.ylabel('frequency')
    plt.title(pngname)
    plt.tight_layout()
    plt.savefig(outdir + pngname + '.png', transparent=True)
   # plt.close()


class plotSpct:
    
        
    def __init__(self,vals,fs,spctnperseg,spctnoverlap, lof, hif, outdir, pngname, doPlot = True):
  
    
      self.sf, self.st, self.sxx = signal.spectrogram(vals, fs,                
      nperseg = spctnperseg, noverlap = spctnoverlap)
      
      if doPlot:
          sfmask = np.logical_and(self.sf>=lof, self.sf<hif)   
          plotMesh(self.st, self.sf[sfmask], self.sxx[sfmask,:], outdir, pngname) 
   
    def getPow(self,lof,hif):
      sfmask = np.logical_and(self.sf>=lof, self.sf<hif)    
      return self.sxx[sfmask,:].mean(axis = 0)


#dirname_in = 't42_data-m57'
dirname_in = './'
outdirbase = dirname_in + '/analysis_results/'
#outdirbase = dirname_in + '/rasters/test/'
par1 = 11
par2 = 0
nstr =  '_' + str(par1) +  '_' + str(par2)
filename11 = 't42'
popstrs = ['OLM_pop','PYR_pop','PVBC_pop']

g3tstart = 2000
ratebinsize = 4
fs = int(1000/ratebinsize)
psdnperseg= int(256*8/ratebinsize)
psdnoverlap= int(psdnperseg - 1)

spctnperseg = int(128*2/(ratebinsize))
spctnoverlap = int(spctnperseg - 1)

powranges = []
powranges.append([1,8])
powranges.append([15,50])
spctthetapeakspread = 1
spctgammapeakspread = 8

comod_powranges = []
comod_powranges.append([1,8])
comod_powranges.append([15,50])
comod_low, comod_hi = 1, 10


filename1 = dirname_in + '/' + filename11 + nstr + '_cfg.json'
filename1_data = dirname_in + '/' + filename11 + nstr + '_data.json'




do_cellavspect = True

with open(filename1) as f, open(filename1_data) as f_data:    
  try:
                
    data = json.load(f)
    data_data = json.load(f_data)
    
    reduceData = False
    if reduceData:
        if 'avgRate' in data_data['simData']:
            del data_data['simData']['avgRate']
        if 'V_soma' in data_data['simData']:
            del data_data['simData']['V_soma']
        if 'cells' in data_data['net']:
            del data_data['net']['cells']
        filename1_data_out = dirname_in + '/' + filename11 + nstr + '_reduced_data.json'    
        with open(filename1_data_out, 'w') as f:
             json.dump(data_data, f)
    
    g3 = getRates(data,data_data,'OLM_pop',ratebinsize)
    g3vals =  g3.vals[g3.times>g3tstart]
    g3times = g3.times[g3.times>g3tstart]
    ws1 = welchSig(g3vals, fs, psdnperseg, psdnoverlap)
    wmmask = np.logical_and(ws1.Fxx>=powranges[0][0], ws1.Fxx<powranges[0][1])
    thetapeakfreq = ws1.Fxx[wmmask][np.argmax(ws1.Pxx[wmmask])] 
    bp1 = butterPass(thetapeakfreq,fs,g3vals,g3times)
    
    
    for popstr in popstrs:
        targcells = data_data['net']['pops'][popstr]['cellGids']
        
        
        if do_cellavspect:
            ws1cell = []
            for tracecellnum in targcells:
                g3cell = getRates(data,data_data,None,ratebinsize,tracecellnum)
                g3cellvals =  g3cell.vals[g3cell.times>g3tstart]
                                
                                    
                ws1 = welchSig(g3cellvals, fs, psdnperseg, psdnoverlap)
                ws1cell.append(ws1.Pxx)
            psdcellav = np.array(ws1cell).mean(axis = 0)
            plt.figure('psdcell')
            plt.plot(ws1.Fxx, psdcellav, label=popstr)
        
        grast = getRaster(data_data,popstr)
        plt.figure('raster')
        rast_lotime, rast_hitime = 4000, 5000
        rast_topcell = 100
        tslice1 =  np.where(np.logical_and(grast.times>rast_lotime, grast.times<rast_hitime))
        grt = grast.times[tslice1]
        grv = grast.vals[tslice1]
        cslice1 = np.where(grv < rast_topcell)
        grt =  grt[cslice1]
        grv =  grv[cslice1]
        plt.plot(grt, grv, '.')
     

            
        g3 = getRates(data,data_data,popstr,ratebinsize)
        g3vals =  g3.vals[g3.times>g3tstart]
        g3times = g3.times[g3.times>g3tstart]
        
        bp1.call(g3vals,g3times,popstr,outdirbase,nstr)
        
        ws1 = welchSig(g3vals, fs, psdnperseg, psdnoverlap)
        plt.figure('psd')
        plt.plot(ws1.Fxx, ws1.Pxx, label=popstr)
       
        spct = plotSpct(g3vals,fs,spctnperseg,spctnoverlap, 0, 125,
                             outdirbase, 'spct' + nstr + '_' + popstr)
        
        
        if popstr == 'OLM_pop':
           powTSC1 = powTSCalc(ws1, powranges, spctthetapeakspread, spctgammapeakspread)
           powTSC1.call(spct,outdirbase,nstr)
        
      
        

        pm1 = pacMethod(g3vals, comod_low, comod_hi, fs, outdirbase, 
                        'pac' + nstr + '_' + popstr, 'penny')
 
        
        hifr = pm1.estimator.high_fq_range
        lofr = pm1.estimator.low_fq_range
        lofmask = np.logical_and(lofr>=comod_powranges[0][0], lofr<comod_powranges[0][1])
        hifmask = np.logical_and(hifr>=comod_powranges[1][0], hifr<comod_powranges[1][1])
        comod_orig = pm1.estimator.comod_
        comod_array = comod_orig[lofmask,:] 
        comod_array = comod_array[:,hifmask]
        comod_mean = comod_array.mean(axis=1).mean(axis=0)
        comod_max_value = np.max(comod_orig)
        print(popstr + ' pac values: ', comod_mean, comod_max_value) 
        
        
    plt.figure('raster') 
    plt.title('raster' + nstr)
    plt.xlabel('time (ms)')
    plt.ylabel('cell index #')
    plt.savefig(outdirbase +'raster' + nstr +'.png')
            
    plt.figure('psd')
    plt.title('population psd' + nstr)
    plt.yscale("log")
    plt.xlabel('frequency (Hz)')
    plt.ylabel('power')
    plt.legend()
    plt.savefig(outdirbase +'psd' + nstr + '.png')
    
    if do_cellavspect:
        plt.figure('psdcell')
        plt.title('mean cell psd' + nstr)
        plt.yscale("log")
        plt.xlabel('frequency (Hz)')
        plt.ylabel('power')
        plt.legend()
        plt.savefig(outdirbase +'psdcell' + nstr + '.png')
    
   
  except JSONDecodeError:
    print("JSONDecodeError")             