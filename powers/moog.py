#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from scipy.ndimage import convolve1d
import matplotlib.pyplot as plt

"""
This code is used to create a MOOG's driver file and run MOOGSILENT.
It is optimized for single-line analysis. Usage for a large region
of the spectrum, with many lines, is not advised.
"""


class Driver(object):

    """
    The MOOG driver object.
    
    Parameters
    ----------
    synth_interval : sequence
        The synthesis wavelength interval lower and upper limits in angstrons.
        Example: (6000, 6100).

    abunds : ``numpy.array``
        The atomic number (first column) and the abundance (second column) of
        the elements to be synthetisized.
        Example: numpy.array([[26, -0.05], [32, 0.01]])

    step: float, optional
        The wavelength step size of the synthesis. Default is 0.1.

    opac: float, optional
        Wavelength point to consider opacity contributions from neighboring
        transitions. Default is 2.0.

    wl_shift: float, optional
        Wavelength shift to be applied to the observed spectrum. Default is 0.

    v_shift: float, optional
        Doppler velocity shift to be applied to the observed spectrum. Default
        is 0.

    y_shift_add: float, optional
        Additive shift to be applied to the observed spectrum. Default is 0.

    y_shift_mult: float, optional
        Multiplicative factor to be applied to the observed spectrum. Default
        is 1.0 (no modification).

    gauss: float, optional
        Value of the 1 sigma dispersion of the Gaussian smoothing to be applied
        to the synthetic spectrum. Default is 0.

    lorentz: float, optional
        Default is 0.

    eps: float, optional
        Limb darkening coefficient. Default is 0.6.

    macro_v: float, optional
        Macroturbulence velocity of the star. Default is 0.

    vsini: float, optional
        The projected rotational velocity of the star. Default is 0.

    model_in: str, optional
        Name of the atmosphere model input file. Default is 'star.mod'.

    linelist_in: str, optional
        Name of the line list input file. Default is 'lines.dat'.

    observed_in: str, optional
        Name of the input file containing the observed spectrum. Default is
        'spectrum.dat'.

    std_out: str, optional
        Name of the standard output file. Default is 'vm_long.out'.

    summary_out: str, optional
        Default is 'vm_li.out'.

    smoothed_out: str, optional
        Name of the output file containing the smoothed synthetic spectrum.
        Default is 'vm_smooth.out'.

    atmosphere: int, optional
        Default is 1.

    molecules: int, optional
        Default is 1.

    trudamp: int, optional
        Default is 1.

    lines: int, optional
        Default is 1.

    flux: int, optional
        Default is 0.

    damping: int, optional
        Default is 0.

    star_name: str, optional
        Self-explanatory. Default is 'Unnamed star'.
    """

    def __init__(self, synth_interval, abunds, step=0.01, opac=2.0, 
                 wl_shift=0.0, v_shift=0.0, y_shift_add=0.0, y_shift_mult=1.0, 
                 gauss=0.0, lorentz=0.0, eps=0.6, macro_v=0.0, vsini=0.0,
                 model_in='star.mod', linelist_in='lines.dat',
                 observed_in='spectrum.dat', std_out='vm_long.out',
                 summary_out='vm_li.out', smoothed_out='vm_smooth.out', 
                 atmosphere=1, molecules=1, trudamp=1, lines=1, flux=0,
                 damping=0, star_name='Unnamed star', plot=True, savefig=False):

        self.name = star_name
        self.plot_switch = plot
        self.savefig = savefig
        # Output files
        self.standard_out = std_out
        self.summary_out = summary_out
        self.smoothed_out = smoothed_out

        # Input files
        self.model_in = model_in
        self.lines_in = linelist_in
        self.observed_in = observed_in

        # Synthesis parameters
        self.syn_start = synth_interval[0]
        self.syn_end = synth_interval[1]
        self.wl_start = synth_interval[0]
        self.wl_end = synth_interval[1]
        self.step = step
        self.opac = opac
        self.wl_shift = wl_shift
        if int(v_shift) != 0:
            raise NotImplementedError('Doppler shift in the observed spectrum'
                                      'is not implemented yet.')
        self.v_shift = v_shift
        self.y_shift_add = y_shift_add
        self.y_shift_mult = y_shift_mult
        self.gauss = gauss
        self.lorentz = lorentz
        self.dark = eps
        self.macro_v = macro_v
        self.vsini = vsini
        self.N, self.n_cols = np.shape(abunds)
        assert(self.n_cols == 2), 'Number of columns in `abunds` must be 2.'
        if self.N == 1:
            self.Z = abunds[0]
            self.abunds = abunds[1]
        elif self.N > 1:
            self.Z = abunds[:, 0]
            self.abunds = abunds[:, 1]

        # MOOG synth options
        self.atm = atmosphere
        self.mol = molecules
        self.tru = trudamp
        self.lin = lines
        self.flu = flux
        self.dam = damping

        # Reading the observed spectrum
        if isinstance(observed_in, str):
            self.obs_wl = np.loadtxt(observed_in, usecols=(0,)) + wl_shift
            self.obs_flux = np.loadtxt(observed_in, usecols=(1,)) * \
                            y_shift_mult + y_shift_add
        elif observed_in is None:
            self.observed_in = observed_in
        else:
            raise TypeError('observed_in must be ``str`` or ``None``.')

    def create_batch(self):
        """
        Writes the MOOG driver file batch.par
        """
        # Creating batch.par file
        with open('batch.par', 'w') as f:
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
            f.write(' %.4f  %.4f  %.3f  %.3f\n' % (self.v_shift, self.wl_shift,
                                                   self.y_shift_add,
                                                   self.y_shift_mult))
            f.write(' r  %.3f  0.0  %.1f  %.2f  %.1f' % (self.gauss,
                                                         self.dark,
                                                         self.macro_v,
                                                         self.lorentz))

    def run(self):

        """
        Used to run MOOG silent.
        """
        self.create_batch()

        # Running MOOGSILENT
        os.system('MOOGSILENT > moog.log 2>&1')
        os.remove('batch.par')

        # Grabing the synthetic and observed spectrum
        synth_wl = np.loadtxt('vm_smooth.out', skiprows=2, usecols=(0,))
        synth_flux = np.loadtxt('vm_smooth.out', skiprows=2, usecols=(1,))

        if self.vsini > 0.0:
            # Applying smart_cut() before rotational convolution
            synth_wl, synth_flux, self.obs_wl, self.obs_flux = \
                self.smart_cut(synth_wl, synth_flux, self.obs_wl, self.obs_flux)

            # Finding the line info corresponding to the wavelength in question
            lines = np.loadtxt('lines.dat', usecols=(0,), skiprows=1)
            wl_0 = lines[np.where(np.abs(lines - self.syn_start) ==
                                  np.min(np.abs(lines - self.syn_start)))[0][0]]

            # Wavelength to velocity space
            c = 2.998E5  # km/s
            synth_v = c * (synth_wl - wl_0) / wl_0

            # The rotational profile
            vsini_profile = self.rot_prof(synth_v)

            # Convolving the non-rotating spectrum with the rotational profile
            conv_flux = convolve1d(synth_flux, vsini_profile)
            with open('vm_smooth.out', 'w') as f:
                f.write('Smoothed spectrum\n')
                f.write('Do not mind this useless line\n')
                for i in range(len(conv_flux)):
                    f.write('  %.3f     %.5f\n' % (synth_wl[i], conv_flux[i] /
                                                   max(conv_flux)))

        # Plotting/saving a plot
        self.plot(synth_wl, synth_flux)

    def rot_prof(self, vz):
        
        """
        This function creates a rotational profile based on Gray (2005).

        Parameters
        ----------
        vz : ``numpy.array``
            The Doppler velocities from the spectral line center.

        Returns
        -------
        profile : ``numpy.array``
            The rotational profile.
        """
        
        n = len(vz)
        profile = np.zeros(n, float)
        for i in range(n):
            if np.abs(vz[i]) < self.vsini:
                profile[i] = (2. * (1. - self.dark) * (1. - (vz[i] /
                                                             self.vsini) ** 2)
                              ** 0.5 + 0.5 * np.pi * self.dark *
                              (1. - (vz[i] / self.vsini) ** 2)) / \
                             (np.pi * self.vsini * (1. - self.dark / 3))
            else:
                profile[i] = 0.0
        return profile

    @staticmethod
    def smart_cut(wl, flux, obs_wl, obs_flux):
        
        """
        smart_cut() is used to prepare the synthetic spectrum for a convolution
        with the rotational profile.
        """
        
        ind0 = np.where(flux == min(flux))[0][0]
        n = len(wl)
        if ind0 < (n - 1) / 2:
            if (ind0 + 1) % 2 == 0:
                wl = wl[1:2 * ind0]
                flux = flux[1:2 * ind0]
                obs_flux = obs_flux[1:2 * ind0]
                obs_wl = obs_wl[1:2 * ind0]
            else:
                wl = wl[0:2 * ind0 + 1]
                flux = flux[0:2 * ind0 + 1]
                obs_flux = obs_flux[0:2 * ind0 + 1]
                obs_wl = obs_wl[0:2 * ind0 + 1]
        elif ind0 > (n - 1) / 2:
            if (ind0 + 1) % 2 == 0:
                wl = wl[2*(ind0 - (n - 1) / 2) + 1:-1]
                flux = flux[2 * (ind0 - (n - 1) / 2) + 1:-1]
                obs_flux = obs_flux[2 * (ind0 - (n - 1) / 2) + 1:-1]
                obs_wl = obs_wl[2 * (ind0 - (n - 1) / 2) + 1:-1]
            else:
                wl = wl[2 * (ind0 - (n - 1) / 2):]
                flux = flux[2 * (ind0 - (n - 1) / 2):]
                obs_flux = obs_flux[2 * (ind0 - (n - 1) / 2):]
                obs_wl = obs_wl[2 * (ind0 - (n - 1) / 2):]
        return wl, flux, obs_wl, obs_flux

    def plot(self, synth_wl, synth_flux):
        """

        Parameters
        ----------
        synth_wl :

        synth_flux :
        """
        # TODO: also print abundances on plot
        fig, ax = plt.subplots()
        ax.plot(self.obs_wl, self.obs_flux, '.')
        ax.plot(synth_wl, synth_flux)
        ax.ticklabel_format(useOffset=False)
        plt.title(r'%s, $v \sin{i} = %.2f$, $v_{macro} = %.2f$, '
                  r'$FWHM_{gauss} = %.3f$' % (self.name, self.vsini,
                                              self.macro_v, self.gauss))
        plt.xlabel(r'$\lambda$ ($\AA$)')
        plt.ylabel(r'$F_\lambda$ $d\lambda$')
        plt.xlim(self.obs_wl[0], self.obs_wl[-1])
        plt.tight_layout()
        if self.savefig is not None:
            print('Saving plot as %s' % self.savefig)
            plt.savefig(self.savefig)
        if self.plot_switch is True:
            plt.show()
