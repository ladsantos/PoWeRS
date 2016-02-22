#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from scipy.ndimage import convolve1d
import plotter

"""
This code is used to create a MOOG's driver file and run MOOGSILENT.
It is optimized for single-line analysis. Usage for a large region
of the spectrum, with many lines, is not advised.
"""


class driver(object):

    """
    Used to create the driver file for MOOG.
    """


    def __init__(self, params, **kwargs):

        # Output files
        self.standard_out = 'vm_long.out'
        self.summary_out = 'vm_li.out'
        self.smoothed_out = 'vm_smooth.out'

        # Input files
        self.model_in = 'star.mod'
        self.lines_in = 'lines.dat'
        self.observed_in = 'spectrum.dat'

        # Synthesis parameters
        self.syn_start = params[0,0]
        self.syn_end = params[1,0]
        self.wl_start = params[0,1]
        self.wl_end = params[1,1]
        self.step = params[2,0]
        self.opac = params[2,1]
        self.wlshift = params[9,1]
        self.vshift = params[9,0]
        self.yshifta = params[10,0]
        self.yshiftm = params[10,1]
        self.gauss = params[11,0]
        self.lorentz = params[11,1]
        self.dark = params[12,0]
        self.macrov = params[13,0]
        self.rotv = params[13,1]
        self.N = int(params[14,0])
        self.Z = np.array([params[k+15,0] for k in range(self.N)])
        self.abunds = np.array([params[k+15,1] for k in range(self.N)])

        # MOOG synth options
        self.atm = params[3,0]
        self.mol = params[4,0]
        self.tru = params[5,0]
        self.lin = params[6,0]
        self.flu = params[7,0]
        self.dam = params[8,0]


    def create_batch(self):

        """
        Writes the MOOG driver file batch.par
        """

        # Creating batch.par file
        with open('batch.par','w') as f:
            f.truncate()
            f.write('synth\n')
            f.write('standard_out  %s\n' % self.standard_out)
            f.write('summary_out  %s\n' % self.summary_out)
            f.write('smoothed_out  %s\n' % self.smoothed_out)
            f.write('model_in  %s\n' % self.model_in)
            f.write('lines_in  %s\n' % self.lines_in)
            f.write('observed_in  %s\n' % self.observed_in)
            f.write('atmosphere  %i\n' % self.atm)
            f.write('molecules  %i\n' % self.mol)
            f.write('trudamp  %i\n' % self.tru)
            f.write('lines    %i\n' % self.lin)
            f.write('flux/int  %i\n' % self.flu)
            f.write('damping    %i\n' % self.dam)
            f.write('freeform  1\n')
            f.write('plot    3\n')
            f.write('abundances  %i 1\n' % self.N)
            for k in range(self.N):
                f.write('   %i %f\n' % (self.Z[k], self.abunds[k]))
            f.write('isotopes   0  1\n')
            f.write('synlimits\n')
            f.write(' %.1f %.1f %.2f %.1f\n' % (self.syn_start, self.syn_end,
                                                self.step, self.opac))
            f.write('obspectrum  5\n')
            f.write('plotpars  1\n')
            f.write(' %.2f %.2f 0.5 1.05\n' % (self.wl_start, self.wl_end))
            f.write(' %.4f  %.4f  %.3f  %.3f\n' % (self.vshift, self.wlshift,
                                                   self.yshifta, self.yshiftm))
            f.write(' r  %.3f  0.0  %.1f  %.2f  %.1f' % (self.gauss,
                                                          #self.rotv,
                                                          self.dark,
                                                          self.macrov,
                                                          self.lorentz))


class run(object):

    """
    Used to run MOOG silent.
    """


    def __init__(self, **kwargs):

        # Do you want the silent version?
        if ('silent' in kwargs):
            self.silent = kwargs['silent']
        else:
            self.silent = False

        # Save file with a name instead of plotting on window
        if ('save' in kwargs):
            self.save = kwargs['save']
        else:
            self.save = 'window'

        # Write star name at the plot?
        if ('star_name' in kwargs):
            self.name = kwargs['star_name']
        else:
            self.name = 'Unnamed star'


        self.params = np.genfromtxt('params.txt',usecols=(1,2), skip_header=1,
                                    missing_values='None', filling_values=0)
        self.d = driver(self.params).create_batch()
        os.system('MOOGSILENT > moog.log 2>&1')

        # Performing rotational broadening
        self.synth_wl = np.loadtxt('vm_smooth.out',skiprows=2,usecols=(0,))
        self.synth_flux = np.loadtxt('vm_smooth.out',skiprows=2,usecols=(1,))
        self.obs_flux = np.loadtxt('spectrum.dat', usecols=(1,))*self.params[10,1]
        self.obs_wl = np.loadtxt('spectrum.dat', usecols=(0,))+self.params[9,1]
        self.synth_wl, self.synth_flux, self.obs_wl, self.obs_flux = \
            self.smart_cut(self.synth_wl, self.synth_flux, self.obs_wl, self.obs_flux)

        self.lines = np.loadtxt('lines.dat', usecols=(0,), skiprows=1)
        self.wl_0 = self.lines[np.where(
            np.abs(self.lines-self.params[0,0])==np.min(np.abs(self.lines-self.params[0,0]))
            )[0][0]]
        self.c = 2.998E5 # km/s
        self.synth_v = self.c*(self.synth_wl-self.wl_0)/self.wl_0
        self.obs_v = self.c*(self.obs_wl-self.wl_0)/self.wl_0
        self.vsini_profile = self.rot_prof(self.params[13,1], self.synth_v, self.params[12,0])
        self.conv_flux = convolve1d(self.synth_flux, self.vsini_profile)
        #self.ind = np.where(self.conv_flux == np.min(self.conv_flux[30:70]))[0][0]
        #self.ind0 = np.where(self.synth_flux == np.min(self.synth_flux))[0][0]
        #self.shift = self.synth_wl[self.ind0]-self.synth_wl[self.ind]
        #self.M = len(self.conv_flux_0)
        #self.what = -4
        #self.conv_flux = self.conv_flux_0[self.M/4-self.what:3*self.M/4-self.what-1]
        with open('vm_smooth.out', 'w') as self.f:
            self.f.write('Smoothed spectrum\n')
            self.f.write('Do not mind this useless line\n')
            for i in range(len(self.conv_flux)):
                self.f.write('  %.3f     %.5f\n' % (self.synth_wl[i], self.conv_flux[i]/max(self.conv_flux)))

        # Plotting
        if self.silent == False:
            self.p = plotter.line(self.params).plot(mode=self.save,
                                                    star_name=self.name)

    def rot_prof(self, vsini, vz, epsilon):
        self.N = len(vz)
        self.prof = np.zeros(self.N,float)
        for i in range(self.N):
            if np.abs(vz[i]) < vsini:
                self.prof[i] = (2.*(1.-epsilon)*(1-(vz[i]/vsini)**2)**0.5+\
                    0.5*np.pi*epsilon*(1-(vz[i]/vsini)**2))/(np.pi*vsini*(1.-epsilon/3))
            else:
                self.prof[i] = 0.0
        return self.prof

    def smart_cut(self, wl, flux, obs_wl, obs_flux):
        self.ind0 = np.where(flux == min(flux))[0][0]
        #self.cut = 40
        #wl = wl[self.ind0-self.cut:self.ind0+self.cut+1]
        #flux = flux[self.ind0-self.cut:self.ind0+self.cut+1]
        #obs_flux = obs_flux[self.ind0-self.cut:self.ind0+self.cut+1]
        #obs_wl = obs_wl[self.ind0-self.cut:self.ind0+self.cut+1]
        self.N = len(wl)
        if self.ind0 < (self.N-1)/2:
            if (self.ind0+1) % 2 == 0:
                wl = wl[1:2*self.ind0]
                flux = flux[1:2*self.ind0]
                obs_flux = obs_flux[1:2*self.ind0]
                obs_wl = obs_wl[1:2*self.ind0]
            else:
                wl = wl[0:2*self.ind0+1]
                flux = flux[0:2*self.ind0+1]
                obs_flux = obs_flux[0:2*self.ind0+1]
                obs_wl = obs_wl[0:2*self.ind0+1]
        elif self.ind0 > (self.N-1)/2:
            if (self.ind0+1) % 2 == 0:
                wl = wl[2*(self.ind0-(self.N-1)/2)+1:-1]
                flux = flux[2*(self.ind0-(self.N-1)/2)+1:-1]
                obs_flux = obs_flux[2*(self.ind0-(self.N-1)/2)+1:-1]
                obs_wl = obs_wl[2*(self.ind0-(self.N-1)/2)+1:-1]
            else:
                wl = wl[2*(self.ind0-(self.N-1)/2):]
                flux = flux[2*(self.ind0-(self.N-1)/2):]
                obs_flux = obs_flux[2*(self.ind0-(self.N-1)/2):]
                obs_wl = obs_wl[2*(self.ind0-(self.N-1)/2):]
        return wl, flux, obs_wl, obs_flux
