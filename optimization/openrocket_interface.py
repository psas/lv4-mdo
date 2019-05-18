#!/usr/bin/env python
# could improve:
# subsystem estimates
# fudge factors
# thrust estimates and curve
# idiot-proofing
# extendable file io
# make engine file with an xml library

from __future__ import print_function #not sure if i even need this, for python2 to use python3 style printing
from math import pi, log, sqrt
import os
from sys import platform as _platform
import xml.etree.ElementTree as ET # xml library
from zipfile import ZipFile

# Liquid motor variables (Change these!)
OF = 1.3      # O/F ratio, this is somewhat arbitrary but optimal
ullage = 1.1 # percentage of length added to a tank to account for not filling
loss_factor = 1.    # if < 1, then assume thrust is less than ideal (percentage)
factor_of_safety = 2   # factor of safety
mass_realism_coefficient = 2 #fudge factor for design mass, includes contribution of tank structural lugs, feed system, stress concentrations, welds, slosh baffles etc.

rkt_prefix = "./rocket_farm/" # this is where rockets live

# Physics
g_0 = 9.80665  # kg.m/s^2     Standard gravity

# Chemistry, note we are using LOX and IPA, not Ethanol
rho_lox = 1141.0   # kg/m^3  Density of LOX
rho_eth = 852.3   # kg/m^3  Density of Ethanol (with 70% H2O)
rho_ipa = 849.28   # kg/m^3 Density of 64.8% IPA / 35.2% H20

# Tank Materials
Al    = { 'rho': 2800.0,   # kg/m^3       Density
          'Sy':    0.270e9} # Pa          Yield strength
Steel = { 'rho': 7830.0,   # kg/m^3       Density
          'Sy':    0.250e9} # Pa          Yield strength
CF    = { 'rho': 1550.0,   # kg/m^3       Density
          'Sy':    0.450e9} # Pa          Yield strength

# Engine system dimensions
Tank     =  CF     # Choose from above table, ignore steel
gaps     =  0.050  # m ###CHECK ME
m_plumb  =  5.0    # kg
m_feed   =  10.0   # kg
l_feed   =  0.4572 # m, this is 18"
m_ems    =  1      # kg
l_ems    =  0.1016 # m, this is 4" 
m_engine =  3.0    # kg
l_engine =  0.300  # m ###CHECK ME
airframe_offset = 0.007 # m, empirical from lv3
dist_after_f = 36 * 0.0254 # m (converted from in), places fuel tank before avionics and N2
dist_after_o = 0.0 # m, oxygen tank still in front of feedsys, ems, engine (duh)


## Utility Functions
# NAR letter code
def nar_code(thrust, burn_time):
    impulse = thrust*burn_time
    print("")
    print("Total impulse:    %0.0f N.s" % impulse)

    nar_i = int(log(impulse/2.5)/log(2))
    nar_percent = impulse/(2.5*2**(nar_i+1))
    print('NAR:              "%s" (%0.0f%%)' % (chr(66+nar_i), nar_percent*100))
    return chr(66+nar_i), nar_percent*100

# unpack rocket template
def unzip():
    with ZipFile('psas_rocket.ork') as myzip:
        myzip.extract('rocket.ork')

# package our new rocket
def zipit(index):
    with ZipFile(rkt_prefix+'psas_rocket_'+index+'.ork', 'w') as myzip:
        myzip.write('rocket.ork')
    if 'linux' in _platform:
        os.system('rm rocket.ork')
    elif "darwin" in _platform:
        os.system('rm rocket.ork')
    elif "win" in _platform:
        os.system('del rocket.ork')

# pulls all file references from given directory
def all_files(directory):
    for path, dirs, files in os.walk(directory):
        for f in sorted(files):
            yield os.path.join(path, f)

# counts how many rockets are in directory and then increments by 1
def get_index():
    ork_files = [f for f in all_files(rkt_prefix)
               if f.endswith('.ork')]
    previous_iteration_index = len(ork_files)
    return previous_iteration_index + 1

# Consider that there are two tanks, and we will want to divide total mass flow rate and propellant mass
# generic division, oxygen first, fuel second
def proportion(amount):
    stuff_o = amount / (1 + (1/OF))
    stuff_f = amount / (1 + OF)
    return stuff_o, stuff_f


## Individual tank Geometry and Physics
# Radius of the tanks, converts in to m, 0.014 is empirical offset from airframe
def tank_r(total_dia):
    return (total_dia/2 * 0.0254) - airframe_offset

# Radius of airframe, m
def body_r(eng_r):
    return eng_r + airframe_offset

# length of tank based on its radius, and mass and density of contents
def tank_length(m, rho, r):
    return m / (rho * pi * r**2)

# turn total propellant mass and airframe diameter into two tanks of propellants
def split_tanks(prop_mass, total_dia):
    m_o, m_f = proportion(prop_mass)
    r = tank_r(total_dia)
    l_o = tank_length(m_o, rho_lox, r)
    l_f = tank_length(m_f, rho_ipa, r)
    l_o *= ullage # add 10% for ullage
    l_f *= ullage # add 10% for ullage
    return r, l_o, l_f

# tank thickness ###CHECK ME
def tank_thickness(tank, r):
    #P_i = 3.042e6  # Tank pressure in Pa, assuming pressure fed with regulator (roughly 441 psig) ###OLD ASSUMPTION
    P_i = 689475.7 # Tank pressure in Pa (~100 PSI), assuming pressurized by N2 and ramped up later by EFS
    radius_o = r   # outer radius, meters
    design_stress = tank['Sy']/factor_of_safety
    radius_i = sqrt(design_stress * (radius_o**2) / ((2*P_i) + design_stress)) # inner radius
    thickness = radius_o - radius_i
    return thickness

# Tank Mass
def tank_mass(l, tank, r):
    s_area = 2*pi*r*(l + r) #surface area of tank
    mass = s_area * tank_thickness(tank, r) * tank['rho'] * mass_realism_coefficient
    return mass

# bulkhead mass
def bulkhead(r):
    thick = 0.00635 # m (0.25")
    perforation = 0.4 # percentage of holiness
    return perforation * Al['rho'] * thick * pi * r**2 
    
    
## System level functions
# Returns mass of each tank
def tank_builder(r, l_o, l_f):
    return tank_mass(l_o, Al, r), tank_mass(l_f, Tank, r)

# Total length of engine system
def system_length(l_o, l_f):
    return sum([l_o, l_f, l_feed, l_ems, l_engine, 4*gaps])

# For openrocket thrust Curve, since the mass and cm change over time.
# total dry mass of engine system
def system_mass(r, l_o, l_f):
    bulkheads = bulkhead(body_r(r))
    t1, t2 = tank_builder(r, l_o, l_f)
    return sum([m_engine, m_plumb, m_ems, m_feed, t1, t2, 4*bulkheads])

# dry center of mass of engine system, note these are positions not lengths
def dry_c_of_m(r, l_o, l_f):
    m_tank_o, m_tank_f = tank_builder(r, l_o, l_f)
    bulkheads = bulkhead(body_r(r))
    
    # fuel tank cm is in the center of the tank
    cm_tank_f  = l_f / 2.0
    # including gaps since they have bulkheads now
    cm_gap1    = l_f + gaps/2.
    # next tank down (lox) has a gap.
    cm_tank_o  = l_f + gaps + dist_after_f + (l_o/2.0)
    # next gap
    cm_gap2    = l_f + gaps + dist_after_f + l_o + gaps/2.
    # now feedsystem
    cm_feed    = l_f + gaps + dist_after_f + l_o + gaps + dist_after_o + (l_feed/2.0)
    #next gap
    cm_gap3    = l_f + gaps + dist_after_f + l_o + gaps + dist_after_o + l_feed + gaps/2.
    #ems
    cm_ems     = l_f + gaps + dist_after_f + l_o + gaps + dist_after_o + l_feed + gaps + (l_ems/2.)
    #last gap
    cm_gap4    = l_f + gaps + dist_after_f + l_o + gaps + dist_after_o + l_feed + gaps + l_ems + gaps/2.
    # finally the engine
    cm_engine  = l_f + gaps + dist_after_f + l_o + gaps + dist_after_o + l_feed + gaps + l_ems + gaps + (l_engine/2.0)
    
    # sum cm
    dry_cm = sum([cm_tank_f*m_tank_f, cm_tank_o*m_tank_o, cm_feed*m_feed, cm_engine*m_engine, \
        cm_gap1*bulkheads, cm_gap2*bulkheads, cm_gap3*bulkheads, cm_gap4*bulkheads])
    
    dry_cm /= system_mass(r, l_o, l_f)
    return dry_cm

# Total mass of propellants at a time
def dprop_mass(total_propellant, mdot, t):
    return total_propellant - (mdot*t)

# Total center of mass (including propellants) at a time
def c_of_m(prop_mass, total_dia, mdot, t):
    M_o_0, M_f_0 = proportion(prop_mass) # initial propellant masses
    mdot_o, mdot_f = proportion(mdot) # mass flow rates
    r, l_o, l_f = split_tanks(prop_mass, total_dia) # geometry
    
    # dry stuff
    dry_cm = dry_c_of_m(r, l_o, l_f) 
    dry_mass = system_mass(r, l_o, l_f)
    
    # wet stuff
    m_o = dprop_mass(M_o_0, mdot_o, t)
    m_f = dprop_mass(M_f_0, mdot_f, t)
    #accounts for gravity as propellant is spent correctly
    cm_f_l = l_f - tank_length(m_f, rho_ipa, r)/2.0
    cm_o_l = l_f + gaps + dist_after_f + l_o - tank_length(m_o, rho_lox, r)/2.0
    cm_prop = ((m_f*cm_f_l) + (m_o*cm_o_l)) / (m_o + m_f +0.000001) # decimal to avoid divide by 0
    cm = ((dry_cm*dry_mass) + (cm_prop*(m_f + m_o))) / (dry_mass + m_f + m_o)
    return cm


### File IO functions
# Takes rocket properties from optimized trajectory and creates a list of all relevant properties
def print_characteristics(mdot, prop_mass, r, l_o, l_f, index, res_text):
    # Mass flow for each propllent
    mdot_o, mdot_f = proportion(mdot)
    res_text.append("\nOx flow: . . . . . . . %7.3f kg/s" % mdot_o)
    res_text.append("\nFuel flow:             %7.3f kg/s" % mdot_f)
    # Propellent Mass for each propllent
    mprop_o, mprop_f = proportion(prop_mass)
    res_text.append("\nOx mass:               %5.1f kg" % mprop_o)
    res_text.append("\nFuel mass: . . . . . . %5.1f kg" % mprop_f)
    # dimensions of each tank
    res_text.append("\nTank diameters:        %7.3f m" % (r*2))
    res_text.append("\nOx tank length + ullage: . . . .%7.3f m" % l_o)
    res_text.append("\nFuel tank length + ullage:      %7.3f m" % l_f)
    # Tank thickness for each tank (mm)
    thickness_o = tank_thickness(Al, r)
    thickness_f = tank_thickness(Tank, r)
    res_text.append("\nOx tank thickness:        %5.1f mm" % (thickness_o*1000))
    res_text.append("\nFuel tank thickness:        %5.1f mm" % (thickness_f*1000))
    # Mass of each tank
    m_tank_o = tank_mass(l_o, Al, r)
    m_tank_f = tank_mass(l_f, Tank, r)
    res_text.append("\nOx tank mass: . . . . .%5.1f kg" % m_tank_o)
    res_text.append("\nFuel tank mass:        %5.1f kg" % m_tank_f)
    
    # create a file with all this info in it
    with open(rkt_prefix+'psas_rocket_'+index+'_traj.txt', 'w') as traj:
        for line in res_text:
            traj.write(line)
        traj.close()

# create a rocket file for our engine's dimensions and characteristics
def update_body(index, eng_r, l_o, l_f):
    unzip() # unpack template
    with open('rocket.ork', 'rb') as xml_file:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for child in root.iter():
            #set radius, this node propagates itself in openrocket
            for kid in child.iterfind('aftradius'):
                kid.text = str(body_r(eng_r))
            for kid in child.iterfind("*[name='Fuel Tank']"):
                kid.find('length').text = str(l_f) # set fuel tank length
            for kid in child.iterfind("*[name='LOX Tank']"):
                kid.find('length').text = str(l_o) # set lox tank length
        tree.write('rocket.ork')
    zipit(index) # repack template
    print("Template rocket updated!")

# when a father engine and a mother engine love each other very much...
def make_engine(mdot, prop_mass, total_dia, Thrust, Burn_time, Isp, res_text):
    index = str(get_index())
    r, l_o, l_f = split_tanks(prop_mass, total_dia)
    
    print_characteristics(mdot, prop_mass, r, l_o, l_f, index, res_text)
    
    for F in Thrust:
        F *= loss_factor # allows us to introduce inefficiency later if we want
        
    length = system_length(l_o, l_f) + dist_after_f + dist_after_o
    dry_mass = system_mass(r, l_o, l_f)
    
    n = len(Thrust)
    peak = max(Thrust)
    average = float(sum(Thrust) / n)
    
    file_head = """<engine-database>
  <engine-list>
    <engine  mfg="PSAS" code="{code}" Type="Liquid" dia="{diameter}" len="{length}"
    initWt="{total_mass}" propWt="{M_prop}" delays="0" auto-calc-mass="0" auto-calc-cg="0"
    avgThrust="{a_thrust}" peakThrust="{p_thrust}" throatDia="0." exitDia="0." Itot="{impulse}"
    burn-time="{burn_time}" massFrac="{m_frac}" Isp="{Isp}" tDiv="10" tStep="-1." tFix="1"
    FDiv="10" FStep="-1." FFix="1" mDiv="10" mStep="-1." mFix="1" cgDiv="10"
    cgStep="-1." cgFix="1">
    <comments>Optimized engine</comments>
    <data>
""".format(**{'code': 'PSAS'+index,
                  'diameter': r*2*1000,
                  'length': length*1000,
                  'total_mass': (prop_mass + dry_mass)*1000,
                  'M_prop': prop_mass*1000,
                  'a_thrust': average,
                  'p_thrust': peak,
                  'burn_time': Burn_time,
                  'm_frac': dry_mass/(dry_mass + prop_mass),
                  'impulse': average*Burn_time,
                  'Isp': Isp,
        })
    
    data = [] # this is going to be our thrust curve! may be slightly inaccurate if/when altitudes vary between trajectory.py and openrocket
    resolution = float(Burn_time/(n-1)) # sec per step
    for i in range(n):
        t = i * resolution # sec
        data.append('     <eng-data  t="{t}" f="{thrust}" m="{mass}" cg="{cg}"/>\n'.format(**{
            't': t,
            'thrust': Thrust[i],
            'mass': (dry_mass + dprop_mass(prop_mass, mdot, t)) * 1000,
            'cg': c_of_m(prop_mass, total_dia, mdot, t) * 1000,
        }))
    
    file_tail = """
    </data>
  </engine>
</engine-list>
</engine-database>"""
    
    # we're gonna put the engine file in the default location
    prefix = "./"
    if 'linux' in _platform:
        home = os.path.expanduser("~")
        prefix =  os.path.join(home, '.openrocket/ThrustCurves/')
    elif "darwin" in _platform:
        home = os.path.expanduser("~")
        prefix =  os.path.join(home, 'Library/Application Support/OpenRocket/ThrustCurves/')
    elif "win" in _platform:
        home = os.getenv("APPDATA")
        prefix = os.path.join(home, "OpenRocket/ThrustCurves/")
    
    # now write the file, great job!
    with open(os.path.join(prefix, 'psas_motor_'+index+'.rse'), 'w') as eng:
        eng.write(file_head)
        for d in data:
            eng.write(d)
        eng.write(file_tail)
        eng.close()
        
    update_body(index, r, l_o, l_f) # make a rocket to correspond with our new engine
    print("Reopen OpenRocket to run simulation.")
    
