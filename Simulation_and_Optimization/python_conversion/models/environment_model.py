
from datetime import datetime
import numpy as np
import math
import nrlmsise00
import pyhwm2014

from .quaternions import *



class Environment():
    def __init__(self, aero_model, ref_ground_wind, latitude, longitude, height):
        self.ka        = 1.4                   # Ratio of specific heats, air  
        self.Ra        = 287.058                 # Avg. specific gas constant (dry air)
        self.longitude = np.radians(longitude)
        self.latitude  = np.radians(latitude) # degrees north, launch site
        self.launch_pt = np.array([latitude, longitude])
        self.up = np.array([np.cos(self.longitude) * np.cos(self.latitude),
                            np.sin(self.longitude) * np.cos(self.latitude),
                            np.sin(self.latitude)])
        
        self.date_time = datetime(2021, 12, 10, 9, 0, 0)
        
        self.a = 6378137.0 # m, semimajor axis, defining parameter of wgs 84
        self.flattening_inv = 298.257223563 # 1/f, defining parameter along with semimajor axis
        self.mu_earth  = 3.986004418 * 10**14 # m^3/s^2, earth gravitational parameter
        
        # derived parameters
        self.f = 1 / self.flattening_inv
        self.b = 6356752.3142 # m, semiminor axis
        self.e = 0.081819190842622 # eccentricity approximation
        self.earth_rad = 6.378137e6 # m, earth mean equitorial radius
        self.earth_mass = 5.9733328e24 # kg, mass of earth including atmosphere
        
        # more derived parameters
        self.e2 = self.f * (2 - self.f)
        self.a1 = self.a * self.e2
        self.a2 = self.a1 * self.a1
        self.a3 = self.a1 * self.e2 * 0.5
        self.a4 = 2.5 * self.a2
        self.a5 = self.a1 + self.a3
        self.a6 = 1 - self.e2
        
        # ref markely & crassidis
        self.J2 = 1.08262668355e-3
        self.J3 = -2.53265648533e-6
        self.J4 = -1.61962159137e-6
        
        self.initial_position    = self.geodetic_to_ECEF(self.latitude, self.longitude, height)
        self.initial_height      = np.linalg.norm(self.initial_position)
        
        r = np.cross(np.array([0,0,1]), self.up)
        r = r / np.linalg.norm(r)
        theta = np.arccos(np.clip(np.dot(np.array([0,0,1]), self.up), -1, 1))
        v = np.sin(theta/2) * r
        self.initial_orientation = np.array([np.cos(theta/2), v[0], v[1], v[2]])
        
        self.erot   = 7.2921150e-5 # Sidearial Earth rotation rate, for coriolis accel
        
        self.aero_model = aero_model
        
    # coriolis acceleration
    # note, pretty sure openrocket does this wrong by adding a negative sign to v[0]
    # https://github.com/openrocket/openrocket/blob/b5cde10824440a1e125067b7d149873deec424b4/core/src/net/sf/openrocket/util/GeodeticComputationStrategy.java
    # double check the sign, low priority
    def coriolis(self, v):
        acc = -2*self.erot*np.cross(np.array([0, 0, 1]), v)
        return acc # acceleration
    
    # gravitational acceleration with up to J4 zonal terms
    # refer to Markely and Crassidis
    def g_accel(self, position):
        r = np.linalg.norm(position)
        xoverr = position[0] / r
        yoverr = position[1] / r
        zoverr = position[2] / r
        zoverrsquare = zoverr**2
        zoverrcube = zoverr*zoverrsquare
        zoverr4th = zoverrsquare**2
        muoverrsquared = self.mu_earth / r**2
        rearthoverr = self.earth_rad / r
        a   = - position / r
        aj2 = -3/2 * self.J2 * rearthoverr**2 * np.array([(1 - 5 * zoverrsquare) * xoverr,
                                                          (1 - 5 * zoverrsquare) * yoverr,
                                                          (3 - 5 * zoverrsquare) * zoverr])
        aj3 = -1/2 * self.J3 * rearthoverr**3 * np.array([5 * (7 * zoverrcube - 3 * zoverr) * xoverr,
                                                          5 * (7 * zoverrcube - 3 * zoverr) * yoverr,
                                                          3 * (10 * zoverrsquare - 35/3 * zoverr4th - 1)])
        aj4 = -5/8 * self.J4 * rearthoverr**4 * np.array([(3 - 42 * zoverrsquare + 63 * zoverr4th) * xoverr,
                                                          (3 - 42 * zoverrsquare + 63 * zoverr4th) * yoverr,
                                                          -(15 - 70 * zoverrsquare + 63 * zoverr4th) * zoverr])
        return muoverrsquared * (a + aj2 + aj3 + aj4)
        
    # gravitational torque. inputs in inertial frame, outputs in body frame
    def T_gg(self, x, attitude_BI, J):
        r = np.linalg.norm(x)
        n = frame_rotation(attitude_BI, - x / r)
        return (3 * self.mu_earth / r**3) * np.cross(n, J.dot(n)) # kg * m^2 / s^2
 
    # pushes a 2D vector in spherical coords into a 3D vector in cartesian coords, for wind model
    def pushforward(self, lat, lon):
        ca = np.cos(lat)
        sa = np.sin(lat)
        co = np.cos(lon)
        so = np.sin(lon)
        
        d_phi = np.array([[-sa*co, -ca*so],
                          [-sa*so, ca*co],
                          [ca, 0]])
        return d_phi
    
    # atmospheric model and wind model
    # input: ECEF position
    # output: pressure, density, temperature, viscosity, and a wind speed vector
    def atmo(self, x, gust_mag=0, gust_dir=0):
        f107a = 150
        f107 = 150
        ap = 4
        
        lat, long, h = self.ECEF_to_geodetic(x)
        h = h/1000
        
        time = self.date_time
        year = time.year
        doy = int(time.strftime("%j"))
        iyd = int(str(year)[-2:] + str(doy))
        sec = (time.hour * 3600.
                    + time.minute * 60.
                    + time.second
                    + time.microsecond * 1e-6)
        lst = sec / 3600. + long / 15.0

        densities, temps = nrlmsise00.msise_model(self.date_time, h, lat, long, f107a, f107, ap, lst)
        winds = pyhwm2014.hwm14.hwm14(iyd=iyd, sec=sec, alt=h, glat=lat, glon=long, stl=lst, f107a=f107a, f107=f107, ap=[-1, ap])
        rho = densities[5] * 1000 #  mass density, kg/m^3
        T_a = temps[1] # K, should be neutral temp (not exospheric)
        
        p_a = rho * T_a * self.Ra
        dyn_visc = 1.458e-6 * T_a**(3/2) / (T_a + 110.4)
        mu = dyn_visc / rho
        
        
        d_phi = self.pushforward(lat, long)
        v_wind = d_phi.dot(np.array(winds))
        if gust_mag != 0:
            v_wind   += normalize(v_wind) * gust_mag
        
        return np.array([p_a, rho, T_a, mu]), v_wind
    
    # transformation from WGS-84 geodetic coordinates to ECEF geocentric coordinates
    def geodetic_to_ECEF(self, lat, long, h):
        N = self.a / np.sqrt(1 - (self.e*np.sin(lat))**2)
        temp = (N + h) * np.cos(lat)
        x = temp * np.cos(long)
        y = temp * np.sin(long)
        z = (N*(self.b / self.a)**2 + h) * np.sin(lat) # from wgs def
        return np.array([x, y, z])
    
    # source: https://possiblywrong.wordpress.com/2014/02/14/when-approximate-is-better-than-exact/
    def ECEF_to_geodetic(self, ecef):
        """Convert ECEF (meters) to LLA (radians and meters).
        """
        # Olson, D. K., Converting Earth-Centered, Earth-Fixed Coordinates to
        # Geodetic Coordinates, IEEE Transactions on Aerospace and Electronic
        # Systems, 32 (1996) 473-476.
        w = math.sqrt(ecef[0] * ecef[0] + ecef[1] * ecef[1])
        z = ecef[2]
        zp = abs(z)
        w2 = w * w
        r2 = z * z + w2
        r  = math.sqrt(r2)
        s2 = z * z / r2
        c2 = w2 / r2
        u = self.a2 / r
        v = self.a3 - self.a4 / r
        if c2 > 0.3:
            s = (zp / r) * (1 + c2 * (self.a1 + u + s2 * v) / r)
            lat = math.asin(s)
            ss = s * s
            c = math.sqrt(1 - ss)
        else:
            c = (w / r) * (1 - s2 * (self.a5 - u - c2 * v) / r)
            lat = math.acos(c)
            ss = 1 - c * c
            s = math.sqrt(ss)
        g = 1 - self.e2 * ss
        rg = self.a / math.sqrt(g)
        rf = self.a6 * rg
        u = w - rg * c
        v = zp - rf * s
        f = c * u + s * v
        m = c * v - s * u
        p = m / (rf / g + f)
        lat = lat + p
        if z < 0:
            lat = -lat
        return (np.degrees(lat), np.degrees(math.atan2(ecef[1], ecef[0])), f + m * p / 2)
