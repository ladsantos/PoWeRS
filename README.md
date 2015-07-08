# pwoogs
pwoogs is a Python wrapper for MOOG synth, an application widely used in astronomy for synthesis of spectral lines. I recently started working on this after being fed up with the usual interface and plotting of MOOG. The code is very bare at the moment, but I'll be working on improving it, adding more features and a proper documentation as time goes on. For now, you can check the source and figure out how to use it.

Running pwoogs is as simple as:

    from pwoogs import moog
    m = moog.run()
