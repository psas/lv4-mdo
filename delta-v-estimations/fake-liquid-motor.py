#!/usr/bin/env python
from __future__ import print_function
from math import pi, log
import os
from sys import platform as _platform

# Assumptions
Isp      =   230.0      # s       Specific Impulse
OF       =     1.6      #         O/F ratio
r        =     0.0762   # m       Radius of the tanks (ID of rocket)

# Physics
g_0      =     9.80665  # kg.m/s^2     Standard gravity

# Tank Material
Al    = { 'rho': 2800.0,   # kg/m^3       Density
          'Ftu':    0.214} # GPa          Ultimate strength
Steel = { 'rho': 7830.0,   # kg/m^3       Density
          'Ftu':    0.862} # GPa          Ultimate strength
CF    = { 'rho': 1550.0,   # kg/m^3       Density
          'Ftu':    0.895} # GPa          Ultimate strength

# Chemestry
rho_lox  =   1141.0   # kg/m^3  Desity of LOX
rho_eth  =    852.3   # kg/m^3  Desity of Ethanol (with 70% H2O)

# Extra
m_engine =      2.0   # kg
l_engine =      0.300 # m
m_plumb  =      5.0   # kg
l_plumb  =      0.350 # m
gaps     =      0.100 # m

# Variables (Change these!)
Thrust    =  2500.0    # N       Thrust of engine
Burn_time =    45.0    # s       Duration of the burn
Tank      = Al         # Choose from above table


# ### Mass and Flow
# 
# Given the assumptions above we can solve the mass flow rate through the system:

# Total mass flow rate
mdot = Thrust / (g_0*Isp)

# Mass flow for each propllent
mdot_o = mdot / (1 + (1/OF))
mdot_f = mdot / (1 + OF)

# Propellent Mass
M_o = Burn_time * mdot_o
M_f = Burn_time * mdot_f
M_prop = M_o + M_f


print("Total Propellent Mass: %5.1f kg" % M_prop)
print("Ox mass:               %5.1f kg" % M_o)
print("Fuel mass: . . . . . . %5.1f kg" % M_f)
print("Mass flow:             %7.3f kg/s" % mdot)
print("Ox flow: . . . . . . . %7.3f kg/s" % mdot_o)
print("Fuel flow:             %7.3f kg/s" % mdot_f)

# ### Tank Geometry
def tank_length(m, rho):
    l = m / (rho*pi*r*r)
    return l

l_o = tank_length(M_o, rho_lox)
l_o += l_o*0.1 # add 10% for ullage
l_f = tank_length(M_f, rho_eth)
l_f += l_f*0.1 # add 10% for ullage
length = sum([l_o, gaps, l_f, gaps, l_plumb, gaps, l_engine])

print("Ox tank length: . . . .%7.3f m" % l_o)
print("Fuel tank length:      %7.3f m" % l_f)
print("System length: . . . . %7.3f m" % length)


# ### Tank Mass

def tank_mass(l):
    area = 2*pi*r*l + 2*pi*r*r
    thickness = 9.448e-4 / Tank['Ftu']
    print("Tank thickness:        %5.1f mm" % (thickness*1000))
    material = area * thickness
    mass = material * Tank['rho']
    return mass

m_tank_o = tank_mass(l_o)
m_tank_f = tank_mass(l_f)
print("Ox tank mass: . . . . .%5.1f kg" % m_tank_o)
print("Fuel tank mass:        %5.1f kg" % m_tank_f)


dry_mass = sum([m_engine, m_plumb, m_tank_o, m_tank_f])
print("Dry Mass: . . . . . . .%5.1f kg" % dry_mass)

# lox tank cm is in the center of the tank
cm_tank_o = l_o / 2.0

# next tank down has a gap.
cm_tank_f = l_o + gaps + (l_f/2.0)

# next down
cm_plumb = l_o + gaps + l_f + gaps +  (l_plumb/2.0)

# finally
cm_engine = l_o + gaps + l_f + gaps + l_plumb + gaps + (l_engine/2.0)

# sum cm
dry_cm = sum([cm_tank_o*m_tank_o, cm_tank_f*m_tank_f, cm_plumb*m_plumb, cm_engine*m_engine])
dry_cm /= dry_mass

#print "Dry CM:     %0.3f m" % dry_cm


# ## Thrust Curve
# 
# The mass and cm change over time.

# In[7]:

def mass(t):
    m = M_prop - (mdot*t)
    return m

def cm(t):
    m_o = M_o - (mdot_o*t)
    m_f = M_f - (mdot_f*t)
    lcm_o = l_o - tank_length(m_o, rho_lox)/2.0
    lcm_f = l_o + gaps + l_f - tank_length(m_f, rho_eth)/2.0
    cm_prop = ((m_o*lcm_o) + (m_f*lcm_f))/(m_o + m_f+0.000001)
    cm = ((dry_cm*dry_mass) + (cm_prop*(m_f+m_o)))/(dry_mass+m_f+m_o)
    return cm


# In[8]:

# NAR letter code
impulse = Thrust*Burn_time
print("")
print("Total impulse:    %0.0f N.s" % impulse)

nar_i = int(log(impulse/2.5)/log(2))
nar_percent = impulse/(2.5*2**(nar_i+1))
print('NAR:              "%s" (%0.0f%%)' % (chr(66+nar_i), nar_percent*100))

file_head = """<engine-database>
  <engine-list>
    <engine  mfg="PSAS" code="P10000-BS" Type="Liquid" dia="{diameter}" len="{length}"
    initWt="{total_mass}" propWt="{M_prop}" delays="0" auto-calc-mass="0" auto-calc-cg="0"
    avgThrust="{thrust}" peakThrust="{thrust}" throatDia="0." exitDia="0." Itot="{impulse}"
    burn-time="{burn_time}" massFrac="{m_frac}" Isp="{Isp}" tDiv="10" tStep="-1." tFix="1"
    FDiv="10" FStep="-1." FFix="1" mDiv="10" mStep="-1." mFix="1" cgDiv="10"
    cgStep="-1." cgFix="1">
    <comments>Made up</comments>
    <data>
""".format(**{'diameter': r*2*1000,
              'length': length*1000,
              'total_mass': (M_prop+dry_mass)*1000,
              'M_prop': M_prop*1000,
              'thrust': Thrust,
              'burn_time': Burn_time,
              'm_frac': dry_mass/(dry_mass+M_prop)*1000,
              'impulse': Thrust*Burn_time,
              'Isp': Isp,
    })

data = []
n = 100
res = Burn_time/float(n-1)
for i in range(n):
    t = i * res
    data.append('     <eng-data  t="{t}" f="{thrust}" m="{mass}" cg="{cg}"/>\n'.format(**{
        't': t,
        'thrust': Thrust,
        'mass': dry_mass*1000 + mass(t)*1000,
        'cg': cm(t)*1000,
    }))

file_tail = """
    </data>
  </engine>
</engine-list>
</engine-database>"""

prefix = "./"
if 'linux' in _platform:
    home = os.path.expanduser("~")
    prefix =  os.path.join(home, '.openrocket/ThrustCurves/')
elif _platform == "darwin":
    home = os.path.expanduser("~")
    prefix =  os.path.join(home, 'Library/Application Support/OpenRocket/')
elif "win" in _platform:
    prefix = os.path.join(os.getenv("APPDATA"), "OpenRocket/ThrustCurves/")

with open(os.path.join(prefix, 'psas_motor.rse'), 'w') as eng:
    eng.write(file_head)
    for d in data:
        eng.write(d)
    eng.write(file_tail)

print("Reopen OpenRocket to run simulation")
