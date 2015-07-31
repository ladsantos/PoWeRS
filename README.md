# pwoogs

pwoogs is a Python wrapper for MOOG synth, an application widely used in astronomy for synthesis of spectral lines. I recently started working on this after being fed up with the usual interface and plotting of MOOG. The code is very bare bones at the moment, but I'll be working on improving it, adding more features and a proper documentation as time goes on. Now it has a rotation estimator!

Dependencies
------------

Unfortunately, pwoogs depends on a full installation of MOOGSILENT, which does require superMONGO (SM), a very old and clunky closed-source plotting routine that we don't even use. I know, this is kind of bullshit, but in the future I'll try to figure out a way to avoid an installation of SM. For now, this is what you need installed on your machine:

* MOOG silent (therefore, also superMONGO, Fortran and C compilers)
* numpy
* matplotlib
* pip (optional: only for [un-]installation)
* crepe (optional: only if you plan to use the parameter estimation routines)

If you have problems installing SM or MOOG, which you most certainly will, take a look at the documentation of pwoogs. There is a bunch of tips there.

(Un-)Installation
------------

You might need superuser privileges to perform (un-)installation. If you have python installed through Anaconda, you won't need them.

Once you have the dependencies installed, you should declare the path of the MOOGSILENT installation to your system. One way to do this, for example, on Debian-based systems, is to open your .bashrc profile and put this line inside it:

    PATH=/whatever/path/you/installed/moogjul2014:$PATH

In order to install pwoogs, download it and run the setup:
 
    python setup.py install

After that, running pwoogs should be as simple as:

    from pwoogs import moog
    m = moog.run()
    
You can run pwoogs in silent mode (no plot) by using (contrary to the original MOOG, you don't need to compile a different version of the program! Horray!)

    m = moog.run(silent=True)

You can uninstall pwoogs simply by entering the following line on your terminal:

    pip uninstall pwoogs
    
Latest additions and corrections
------------

* estimate does not optimize for the wavelength shift anymore
* estimate.rotation now accepts a user-defined Gaussian smoothing and has more optional parameters
* estimate.rotation now requires a line-list file (MOOG standard)
* estimate.rotation now has a new, tunable performance function
* estimate.rotation.find() now has a silent option

Known issues
------------

* The code produce no or wrong results when there are too many points on the spectrum (a few tens of thousands, maybe?). This is a problem with MOOG, specifically, and there is not much we can do about it, except spoon-feeding it a nice small array of points.
* The cross-entropy method can converge too fast sometimes. This is a particular trait of this method. When this happens, you'll either get a visually bad fit (local minimum?), or a RuntimeWarning or IndexError by crepe, which essentially means there happened a division by zero when evaluating the new values for the optimization. The best way to deal with this is to use an appropriate smoothing factor (alpha and beta on crepe - see its [still incomplete] manual or try different values).

To-do list
------------

* Documentation for estimate
* A silent version of estimate
* Better error handling
* Undo the need of having some files already created (for instance, batch.par)

Version
------------

pwoogs 0.1 (Alpha) build 1507
