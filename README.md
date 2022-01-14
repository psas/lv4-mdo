# Liquid Engine Analysis

Current location for LV4 Optimization (MDO) algorithms and work.

The archived LV4_Optimization white paper and notebook is essential reading for background info and conceptual explanation. My code follows directly from it and some other archived work.

Note, Openrocket only works on Java 8 and below. Check out the readme in rocket_farm.

# Active Development
* [LV4 MDO (3rd rev)](http://nbviewer.jupyter.org/github/psas/liquid-engine-analysis/blob/master/Simulation_and_Optimization/Multidisciplinary_Design_Optimization.ipynb)
* [Canonical LV4 description](http://nbviewer.ipython.org/github/psas/liquid-engine-analysis/blob/master/LV4_canonical/lv4_optimizer_output.txt)

# Archive
* [LV4 fuel tank FEA](http://nbviewer.ipython.org/github/psas/liquid-engine-analysis/blob/master/archive/AirframeFEA/LV4FuelTankParameters.ipynb)
* [Centifugal Pump Estimations](http://nbviewer.ipython.org/github/psas/liquid-engine-analysis/blob/master/archive/electric_pump_calcs/pump_sizing.ipynb)
* [LV4 MDO (1st rev)](http://nbviewer.ipython.org/github/psas/liquid-engine-analysis/blob/master/archive/LV4_Optimization.ipynb)
* [LV4 Development Notes](http://nbviewer.ipython.org/github/psas/liquid-engine-analysis/blob/master/archive/rocket_notes)
* [Aerobee analysis](http://nbviewer.ipython.org/github/psas/liquid-engine-analysis/blob/master/archive/aerobee-150-reconstruction/AJ11-26.ipynb)
* Misc. abandonware

# Python Environment
1) Install gfortran
   Required by [RocketCEA](https://rocketcea.readthedocs.io/en/latest/quickstart.html) 
    Ubuntu: apt-get install gfortran
    Windows: Dual-boot with Ubuntu :)
2) Install [anaconda](https://www.anaconda.com/products/individual) if you don't have it. Save yourself some time and use anaconda for python
3) Use conda to install the environment
   conda env create -f lv4_mdo_p3p7.yml
   If you need to update this later use:
   conda env update --file lv4_mdo_p3p7.yml --prune
4) Load the python environment:
   conda activate lv4_mdo_p3p7

# Detailed Manual Installation
This MDO has been successfully installed and compiled on Ubuntu and MacOS. 

### Python Environment
* numpy
* scipy
* matplotlib

The following required modules are problematic to install on Windows.

Some of them require a fortran compiler.

Each one will have more details below.

* pyhwm2014
* nrlmsise00 
* rocketCEA
* rbfopt

### pyhwm2014 installation
This module may be installed with `pip install pyhwm2014` but is not recommended.

Instead, download the code from [their GitHub repository](https://github.com/rilma/pyHWM14) and install it manually from the source code installation directions.

However, change the last command from `pip install -e . --process-dependency-links` to `pip install -e ./` as shown below:
```
git clone https://github.com/rilma/pyHWM14.git
cd pyHWM14
pip -q install coveralls
pip install numpy
pip install -e ./
```

### nrlmsise00 installation
Install with `pip install nrlmsise00`.

Requires fortran compiler.


### rocketCEA installation
Install with `pip install rocketcea`.

Requires fortran compiler.


### rbfopt installation
Install with `pip install rbfopt`.

This module requires Bonmin, which is distrubuted under the EPL (Eclipse Public License). 

Bonmin is also dependent on other third party code which are distributed under different licenses than Bonmin.

Check out Bonmin's pages on [Getting Started with Bonmin](https://projects.coin-or.org/Bonmin/wiki/GettingStarted) and [Third Party required code](https://projects.coin-or.org/Bonmin/wiki/ThirdParty)

