These scripts were developed to analyze the debris distribution in LEO and model the effects of a high-powered laser satellite: https://drive.google.com/open?id=1R32u1WvZRNiiA9C3yZFgQ1g8wczhjzQl

Stephen Lisius wrote the ablate.m and main.m Matlab scripts, which modeled the delta-V due to ablation of Aluminum. These gave gross estimates that were used in the main simulation.

David Vallado's Matlab scripts (fundarg.m, iau80in.m, nut80.dat, nutation.m, precess.m, and teme2eci.m) were used to rotate from TEME to J2000. They are provided: https://celestrak.com/software/vallado-sw.php

In order to run the simulation:
`python debSim.py -p`

These scripts require a 64-bit Python interpreter, the Matlab API for Python, and STK. Other required external modules are listed in _requirements.txt_