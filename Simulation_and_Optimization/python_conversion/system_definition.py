


# sets global variables




#from hanging_threads import start_monitoring
#start_monitoring(seconds_frozen=10, test_interval=100)
#import faulthandler
#faulthandler.enable()
#openrocket_interface imports
from math import pi, log, sqrt, exp, cos, sin, radians
import os
import contextlib
from sys import platform as _platform
import xml.etree.ElementTree as ET # xml library
from zipfile import ZipFile

#trajectory imports
from types import SimpleNamespace
import numpy as np
import math
import csv
import copy
from datetime import datetime
import nrlmsise00, pyhwm2014

#optimizer imports
from scipy.optimize import minimize, differential_evolution, shgo, Bounds
import pylab
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc

import rbfopt




########################
### propellant properties
########################



global PROPELLANT_SET

PROPELLANT_SET = True

if PROPELLANT_SET:
    from propellant_optimization import propellant_optimizer

    ipa_wt, of_ratio, p_ch, Tc, MW, gamma, propellant_string = propellant_optimizer(2413166)



########################
### constants
########################

# Physics, chemistry, and materials
G_N      = 9.80665  # kg.m/s^2     Standard gravity
R_UNIV   = 8314.46261815324 # universal gas constant, J/ (K * kmol)
M_PER_IN = 0.0254 # meters per inch




################################################
### materials
################################################

def Material(name, rho, mm=None, mu=None, Sy=None, Su=None, p_v=None):
    return { 'name': name,
             'rho': rho, # kg/m^3       Density
             'mu': mu, # Ns/m^2 Dynamic Viscosity
             'mm': mm, # g/mol Molar Mass
             'Sy': Sy, # Pa           Yield strength
             'Su': Su, # Pa           ultimate tensile strength
             'p_v': p_v} # Pa operating P_v

# Materials
# https://www.aircraftspruce.com/catalog/cmpages/anh4120honeycomb01-01574.php?clickkey=5444217
# note: i've calculated ~395 kg/m^3 for uniform airframe density based on our density/thickness estimates.
#       however, measurement of LV3.1 module says module is around 130 kg/m^3. Not sure why these don't agree.
NOMEX      = Material('Nomex', 48.06)
CRYOGEL    = Material('Cryogel', 160)
#FIBERGLASS = Material('Fiberglass', 1850, Sy=0.2068e9) # don't know where it came from
FIBERGLASS = Material('Fiberglass', 2460, Sy=0.2068e9) # from airframe team CAD
ALUM       = Material('Aluminum 6061-T6', 2800.0, Sy=0.270e9, Su=0.31e9)
#CFIBER     = Material('Carbon Fiber', 1550.0, Sy=0.450e9) # density from internet
CFIBER     = Material('Carbon Fiber', 1990.0, Sy=0.450e9) # what airframe team uses for CAD
LOX        = Material('LOX', 1141.0, mu=0.000009, p_v=8000) # kg/m^3  Density of LOX
#IPA        = Material('IPA/H20', 849.28) # kg/m^3  Density of 64.8% IPA / 35.2% H20
IPA        = Material('IPA', 786) # kg/m^3  Density of 64.8% IPA / 35.2% H20
H20        = Material('H20', 999.8) # kg/m^3  Density of 64.8% IPA / 35.2% H20
FUEL       = Material('IPA/H20', 0.648* IPA['rho'] + 0.352 * H20['rho'],
                     mu=0.00192, p_v=8840) # kg/m^3  Density of 64.8% IPA / 35.2% H20
# Nitrogen characteristics
# https://github.com/psas/reaction-control/blob/master/pubs/AIAA%20RCS%20Manuscript_FINAL2.pdf
N2_TEMP = 298.15 # K, holding this constant is sketchy but easy
N2_MM   = 28.01 # nitrogen molecular mass [g/mol]
N2_KE   = 1.4 # nitrogen specific heat ratio

# Helium characteristics
HE4_TEMP = 298.15 # K
HE4_MM   = 4.003 # g/mol
HE4_KE   = 1.66

ENG_P_CH = 783649.1830 # chamber pressure, PSI

MANUAL_PROP_CALCS = False
if not PROPELLANT_SET:
    if MANUAL_PROP_CALCS:
        # combustion gas properties ke, Re, T_ch, determined from CEArun
        # with chamber pressure=350 psi, fuel temp=419.15 K,
        #      lox temp=90 K, OF=1.3 for fuel = 64.8% IPA (2propanol) / 35.2% H20
        OF       = 1.3        # O/F ratio, this is somewhat arbitrary but CEA says its good.
        ENG_T_CH = 3097.82 # chamber temperature, K
        ENG_KE   = 1.1251 # specific heat ratio, propellant (aka gammas)
        ENG_MM   = 23.196 # molar mass
    else:
        IPA_WT, OF, ENG_P_CH, ENG_T_CH, ENG_MM, ENG_KE, _ = propellant_optimizer(ENG_P_CH)
        FUEL       = Material('IPA/H20', IPA_WT/100* IPA['rho'] + (100 - IPA_WT)/100 * H20['rho'],
                             mu=0.00192, p_v=8840) # kg/m^3  Density of 64.8% IPA / 35.2% H20






################################################
### launch parameters
################################################

# Launch constants
# Vertical Launch at Alkali Lake
LAUNCH_SITE_ALT = 1299 # m, altitude of launch site above sea level, from freemaptools.com/elevation-finder.htm
LAUNCH_TOWER    = 9.8 # launch rail height in m
LAUNCH_SITE_LOC = [42.977691975376736, -120.02818928644571] # dec deg N, E from google maps
AZ_PERTURB = 1
EL_PERTURB = 0.3421796303215663



################################################
### isogrid parameters
################################################

TANK_IN_DIA     = 11.5 * M_PER_IN # propellant tank inner diameter, m
TANK_IN_RAD     = TANK_IN_DIA * 0.5 # propellant tank inner radius, m
TANK_THICK      = 0.25 * M_PER_IN # propellant tank thickness, m
WELD_TENS_STR   = 1.655e8 / 4 # Pa, alum 6061-T6 welded tensile strength
TANK_OD         = TANK_IN_DIA + 2 * TANK_THICK
SKIN_T          = 0.05 * M_PER_IN #0.00127 # m
INSULTN_THKNS   = TANK_THICK - SKIN_T
RIB_DEPTH       = INSULTN_THKNS - M_PER_IN * 0.05
MIN_RIB_T       = 0.00127 # m
DELTA           = RIB_DEPTH / SKIN_T # rib depth to skin thickness ratio, pg. 2.0.008

RIB_T           = 0.00127 # m
NUM_RADL_DVSNS  = 24 # number of triangles



################################################
### rocket parameters
################################################

CPLNG_RING_THK     = 0.00635 # m
INNER_CF_THK       = 0.015 * M_PER_IN # m
NOMEX_THK          = 0.125 * M_PER_IN # m
OUTER_CF_THK       = 0.04 * M_PER_IN # m
AIRFRAME_THICKNESS = INNER_CF_THK + NOMEX_THK + OUTER_CF_THK #0.00508 # m, (0.2 in)
AIRFRM_IN_RAD   = TANK_OD * 0.5 - AIRFRAME_THICKNESS # rocket inner radius, m
THROTTLE_WINDOW = (100., 500.) # lower and upper bounds of drag force (N) for throttling
MIN_THROTTLE    = 1. # the internet says between 60 - 70% is doable without ruining our lives.
NOSETIP         = 0.3533 #kg, estimated for aluminum nose tip. Add more as needed for stability.
BALLAST         = 2



################################################
### fin geometry
################################################

FIN_ROOT        = 30 * M_PER_IN #0.7 # Root length, m
FIN_TIP         = 13 * M_PER_IN #0.45 # tip length, m
FIN_SEMISPAN    = 16 * M_PER_IN #0.4 # fin span/height, m
FIN_SWEEP_ANGLE = np.radians(40) #np.radians(70) # sweep angle, degrees
FIN_THICKNESS   = 0.125 * M_PER_IN #0.003175 # m, fin thickness from LV3.1
FIN_ROOT_HEIGHT = 3.605 * M_PER_IN #0.082525 # height of bottom of fin root above base of rocket
FIN_BRACKET     = 0.56382 # kg




################################################
### rcs parameters
################################################

# low priority, find constraints for these and optimize
# https://steelheadcomposites.com/composite-pressure-vessels/
RCS_CONTROL = False # whether sims have RCS enabled
RCS_MDOT    = 0.03784/3.15 # kg/s RCS nozzle mass flow rate
RCS_P_E     = 101300 # Pa, RCS nozzle exit pressure
RCS_P_CH    = 689476 # Pa, RCS chamber pressure
MAX_N2_TANK_P   = 6.895e7 # Pa, N2 tank pressure (10k PSI)
N2_TANK_OR  = 0.0825 # m, tank outer radius

################################################
### upper subsystem module dimensions
################################################
# if these change you will have to make changes in structure.ipynb
CPLNG_RING_THK     = 0.00635 # m
INNER_CF_THK       = 0.015 * M_PER_IN # m
NOMEX_THK          = 0.125 * M_PER_IN # m
OUTER_CF_THK       = 0.04 * M_PER_IN # m

COUPLING_RING_HT   = 1.4 * M_PER_IN # m, space between composite modules (not actually used)
HALF_CPL_RING      = 0.5 * COUPLING_RING_HT # m, half of assembly for ring, for metal-metal
THREE_QTR_CPL_RING = 0.75 * COUPLING_RING_HT # m, half of assembly for ring, for metal-composite
CON_NOSE_L         = 43.005 * M_PER_IN #1.1
CYL_NOSE_L         = 1 * M_PER_IN
NOSE_L             = CON_NOSE_L + CYL_NOSE_L       # m, with 1 foot cylinder at end
ERS_L              = 4.5 * M_PER_IN  # m 
RCS_L              = 3.305 * M_PER_IN  # m 
CAM_L              = 3.305 * M_PER_IN  # m
AV_L               = 15 * M_PER_IN # m
N2_L               = 23.5 * M_PER_IN # m 
PAS_L              = 3.305 * M_PER_IN  # m , height metal passthru
MICRO_L            = 9.5 * M_PER_IN # m, height micromodule between tanks
FIN_CAN_L          = 35.5 * M_PER_IN # m, height of fin can module
THRST_PLT          = 0.383 * M_PER_IN # m, height of thrust plate
ENG_CLEARANCE      = 4.59 * M_PER_IN # m, height of thrust plate



################################################
### upper subsystem module dimensions
################################################
# m_' = mass, 'l_' = length
# when you have time, these masses should transfer out of structures.ipynb and go into openrocket_interface.ipynb
ULLAGE             = 1.1          # percentage of length added to a tank to account for not filling
L_FEED             = 0.4572       # m, this is 18"
L_EMS              = 0.1016       # m, this is 4" 
L_ENGINE           = 0.300        # m

################################################
### piping
################################################

PIPE_THK      = 0.065 * M_PER_IN
N2_PIPE_OD    = 0.5 * M_PER_IN
FUEL_PIPE_OD  = 0.75 * M_PER_IN *1.5
N2_PIPE_IR    = N2_PIPE_OD * 0.5 - PIPE_THK
FUEL_PIPE_IR  = FUEL_PIPE_OD * 0.5 - PIPE_THK
N2_TO_ENG_L   = 173.725 * M_PER_IN
N2_TO_LOX_L   = 107.52 * M_PER_IN
N2_TO_IPA_L   = 34.2 * M_PER_IN
IPA_TO_ENG_L  = 80.5 * M_PER_IN



################################################
### EFS parameters
################################################

# EFS System parameters
USE_EFS = False
BAT_SPEC_POW       = 0.4 * 125 * 60 # W/kg, assuming 60 s operation time  (conservative est.)
#https://hobbyking.com/en_us/turnigy-heavy-duty-5000mah-7s-60c-lipo-pack-w-xt90.html?queryID=&objectID=74820&indexName=hbk_live_magento_en_us_products
MOT_SPEC_POW       = 0.4 * 9800 / (2.53 + 0.406) # W/kg, motor + ESC

#https://hobbyking.com/en_us/turnigy-rotomax-150cc-size-brushless-outrunner-motor.html?queryID=&objectID=47150&indexName=hbk_live_magento_en_us_products_hbk_price_stock_2_group_0_desc
#https://hobbyking.com/en_us/turnigy-fatboy-v2-300a-esc-4-15s-opto.html?queryID=&objectID=46320&indexName=hbk_live_magento_en_us_products_hbk_price_stock_2_group_0_desc
LOX_TANK_P         = 1859818.5060 # Pa, lox tank pressure
IPA_TANK_P         = 1530286.4229 # Pa, ipa tank pressure
D_PIPE             = FUEL_PIPE_IR * 2 #0.0157  # m Plumbing Pipe Diameter
A_PIPE             = np.pi/4 * D_PIPE**2 # m^2 Cross Sectional Area of Plumbing Pipe
EPSILON_PIPE       = 1.5 *10**(-6) # m Drawn Tubing Relative Roughness
PLUMBING_L_F       = 2.692 # m Length of Straight Pipe Section IPA
PLUMBING_L_O       = 1.046 # m Length of Straight Pipe Section LOX
UNDER_N2_M         = 6 # kg, FLIPS stuff here
ABOVE_FUEL_M       = 6 # kg, FLIPS stuff here
BETWEEN_TANKS_M    = 8 # kg, FLIPS stuff here
PLUMBING_M         = 11 # kg, assumed mass of plumbing in fin can
K_L                = 0.3 # Loss Coefficient for Regular Flanged 90
PUMP_EFF           = 0.6 # pump efficiency
U_SS               = 7000 # suction sp. speed
FRIC_F             = 0.025 # Friction Factor for Isopropyl Alcohol
FRIC_O             = 0.0149 # Friction Factor for Liquid Oxygen
#h_LPintle          = -689.5*1000 # Pa Pressure Loss from Pintle (estimated)
DELP_REGEN         = 0 #1.379 * 10**6 # Pa Guessed pressure loss from regenerative cooling channels (200 PSI)
DELP_INJ_F         = 344.8588707*1000 # Pa Experimental Pressure Loss across Pintle
DELP_INJ_O         = 689.5*1000 # Pa Experimental Pressure Loss across Pintle

LFETS_PIPE_AREA    = np.pi * FUEL_PIPE_IR**2

REGEN_F_COEFF = 0.07 # regen channel fric coeff (per Emilio)
REGEN_D       = 3e-3  # m, regen channel hydraulic diameter
REGEN_L       = 333e-3 #m, regen channel total length
REGEN_N       = 75 # number of regenerative channels
REGEN_MULT    = REGEN_F_COEFF * REGEN_L / (REGEN_D * 2)



################################################
### recovery system
################################################

# drag coefficients of parachutes
C_D_DROGUE = 0.97 # from rocketman
C_D_MAIN = 2.2 # from rocketman






################################################
### optimization parameters
################################################



DELTA      = 10**(-3)  # a guess at a good margin for design "convergence"
MU_0       = 0.0025    # barrier parameter, this value lets altitudes approach lower bound pretty quickly
RHO_0      = 5         # penalty parameter, i'm still playing with this value.
DT         = 0.025      # change starting time-step for trajectory simulation
ITERATIONS = 1         # number of escalating iteration sequences








################################################
### initial design guess
################################################



# be sure that you start with a feasible design, otherwise the problem will be ill-conditioned
M_PROP = 237.0335   # propellant mass (kg)
MDOT   = 6.8475    # Propellant mass flow rate (kg/s)
P_E    = 106525.5088  # Exit Pressure (Pa)

# initial design vector
X0 = np.array([M_PROP, MDOT, P_E, 0, FIN_ROOT, FIN_TIP, FIN_SWEEP_ANGLE, FIN_SEMISPAN, FIN_THICKNESS])





################################################
### optimization constraints
################################################



CONS_IMPLS   = 889600                    # maximum impulse, N s
CONS_AOA     = 8.                       # maximum angle of attack
CONS_MASS    = 450.                      # GLOW constraint, kg, somewhat arbitrary
CONS_LS      = 22.                       # min launch speed from 60' tower constraint, m/s
CONS_TWR     = 2.                        # TWR constraint
CONS_S_CRIT  = 0.35                      # Critical pressure ratio constraint
CONS_ACCEL   = 15.                       # Max acceleration constraint, g's
CONS_LD      = 25.                       # L/D ratio constraint, slightly arbitrary
CONS_ALT     = 105000. #* 0.5                   # Min altitude constraint, m
CONS_THRUST  = 10000                      # max ground-level thrust, N
CONS_CEILING = 150000.                   # base-11 maximum apogee requirement, m
CONS_STBLTY  = 2.0                       # minimum in flight stability margin caliber
#CONS_EFS     = 11000                      # maximum EFS pump power, W
CONS_EFS     = 0  # EFS is being removed
CONS_V_LFETS = 9.144                    # maximum fluid velocity in test stand
CONS_TANK_MIN = 689476 # Pa, minimum tank pressure (so FLIPS can use regulators)




################################################
### openrocket helper functions and definitions
################################################


RKT_PREFIX = "../rocket_farm/" # the rockets live on a cute little farm upstate where they frolic in fields

## Utility Functions
# unpack rocket template temporarily
def unzip():
    with ZipFile('../LV4_canonical/template.ork') as myzip:
        myzip.extract('rocket.ork')

# package our new rocket and remove temporary template
def zipit(index):
    with ZipFile(RKT_PREFIX+'psas_rocket_'+index+'.ork', 'w') as myzip:
        myzip.write('rocket.ork')
    if 'linux' in _platform:
        os.system('rm rocket.ork')
    elif "darwin" in _platform:
        os.system('rm rocket.ork')
    elif "win" in _platform:
        os.system('del rocket.ork')

# pulls ALL file references from given directory
def all_files(directory):
    for path, dirs, files in os.walk(directory):
        for f in sorted(files):
            yield os.path.join(path, f)

# counts how many rockets are in our directory and then increments by 1
def get_index():
    ork_files = [f for f in all_files(RKT_PREFIX)
                   if f.endswith('.ork')]
    return len(ork_files) + 1




################################################
### quaternion helper functions
################################################


# Consider that there are two tanks, and we will want to divide total mass flow rate and propellant mass
# oxygen first, fuel second
def proportion(amount, OF):
    stuff_o = amount * OF/(1 + OF)
    stuff_f = amount * 1/(1 + OF)
    return stuff_o, stuff_f

# this is hamilton's quaternion product
def product(a, b):
    v = b[0] * a[1:] + a[0] * b[1:] + np.cross(a[1:], b[1:])
    return np.array([a[0] * b[0] - np.dot(a[1:], b[1:]), v[0], v[1], v[2]])

# this is the inverse of a unit quat
def conjugate(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

# this rotates a vector with a fixed frame
def sandwich(q, v):
    return product(q, product(np.array([0, v[0], v[1], v[2]]), conjugate(q)))[1:]

# this rotates a frame with a fixed vector
def frame_rotation(q, v):
    return sandwich(conjugate(q), v)

# this constrains a quat to S3 or vector to S2
def normalize(q):
    norm = np.linalg.norm(q)
    norm = norm if norm !=0 else 1
    return q / norm

def eulerangle_to_quat(RA, dec, orientation):
    '''Encodes star tracker's attitude representation as a quaternion in
    3-2-1 order (yaw, pitch, roll). This quaternion transforms the inertial frame to the body frame.

    :params: right ascenscion, declination, roll.
    :returns: Quaternion representation of attitude.'''
    RA = RA / 2
    dec = dec / 2
    ortn = orientation / 2
    c_phi = np.cos(ortn)
    s_phi = np.sin(ortn)
    c_theta = np.cos(dec)
    s_theta = np.sin(dec)
    c_psi = np.cos(RA)
    s_psi = np.sin(RA)
    return np.array([c_phi * c_theta * c_psi  +  s_phi * s_theta * s_psi,
                     s_phi * c_theta * c_psi  -  c_phi * s_theta * s_psi,
                     c_phi * s_theta * c_psi  +  s_phi * c_theta * s_psi,
                     c_phi * c_theta * s_psi  -  s_phi * s_theta * c_psi
                     ])
