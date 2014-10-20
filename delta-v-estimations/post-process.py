#!/usr/bin/env python
from __future__ import print_function
import csv
from math import log, fabs

times = []
t_bo = 0
altitudes = []
velocitys = []
masses = []
with open('sim.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        time = row[0]
        if time[0] == '#':
            if 'BURNOUT' in time:
                t_bo = time.split('t=')[1].split(' ')[0]
                t_bo = float(t_bo)
            continue

        times.append(float(time))
        altitudes.append(float(row[ 1]))
        velocitys.append(float(row[ 2]))
        masses.append(   float(row[19]))

i_t_bo = times.index(t_bo)

Isp      =   230.0      # s       Specific Impulse
g_0      =     9.80665  # kg.m/s/s

m_r = masses[0]/masses[-1]
dv_pure =  g_0 * Isp * log(m_r)
dv_grav = (g_0 * Isp * log(m_r)) - (g_0 * t_bo)

print("Altitude:        %10.1f [km]"  % (max(altitudes)/1000))
print("Burnout V:       %10.1f [m/s]" % velocitys[i_t_bo])
print("GLOW:            %10.1f [kg]"  % masses[0])
print("Empty Mass:      %10.1f [kg]"  % masses[-1])
print("Mass Ratio:      %11.2f"       % m_r)
print("Ideal dV:        %10.1f [m/s]" % dv_pure)
print("Gravity Loss:    %10.1f [m/s]" % fabs(dv_grav - dv_pure))
print("Drag Loss:       %10.1f [m/s]" % fabs(dv_grav - velocitys[i_t_bo]))
