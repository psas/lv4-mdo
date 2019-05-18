# could improve:
# drag coefficients
# trajectory efficiency

from math import sqrt, pi, exp, log, cos, sin, radians
import numpy as np
import csv
import copy
import openrocket_interface as openrkt

# A simple forward Euler integration for rocket trajectories
# could range-kutta or some other integration, but this is sufficient for now

# Physical constants
g_0 = openrkt.g_0  # Gravitational acceleration [m/s^2]
launch_tower = 60. # launch rail height in ft
air_cutoff = 25.   # Pa, ambient pressure at which we no longer care about air

# upper subsystem dimensions
nose_l = 1.25      # m
ers_l = 6 * 0.0254 # m (converted from in)
rcs_l = 6 * 0.0254 # m (converted from in)
av_l = 18 * 0.0254 # m (converted from in)
n2_l = 18 * 0.0254 # m (converted from in)

shirt_l = sum([nose_l, ers_l, rcs_l, av_l, n2_l]) # m, total of above. (what about gaps??)

# helper function for unit conversion
def unit_conv(p_ch, dia, p_e):
    p_ch = p_ch*6894.76    # Chamber pressure, convert psi to Pa
    # LV4 design variables
    dia = dia*0.0254       # Convert in. to m
    p_e = p_e*1000         # Exit pressure, convert kPa to Pa
    return p_ch, dia, p_e

# this class is for conceptual clarity, to keep rocket internals organized
class EngineSys:
    def __init__(self, mdot, p_e, p_ch, T_ch, ke, Re):
        self.f_init = False
        self.mdot = mdot
        self.p_e = p_e
        self.p_ch = p_ch
        self.T_ch = T_ch
        self.ke = ke
        self.Re = Re
        self.F = [self.thrust(self.alt[0])] # be aware that this line will crash if Environment class isn't used
    
    # calculates engine thrust
    def thrust(self, x=0):
        p_ch = self.p_ch
        T_ch = self.T_ch
        p_e = self.p_e
        ke = self.ke
        Re = self.Re
        mdot = self.mdot
        p_a = self.std_at(x)[0]                            # Ambient air pressure [Pa]
        if self.f_init == False:
            p_t = p_ch * (1 + (ke - 1)/2)**(-ke /(ke - 1)) # Throat pressure      [Pa]
            T_t = T_ch*(1 / (1 + (ke - 1)/2))              # Throat temperature   [K]
            self.A_t = (mdot / p_t) * sqrt(Re * T_t / ke)  # Throat area          [m^2]
            A_t = self.A_t
            self.A_e = A_t * (2 / (ke + 1))**(1 / (ke - 1)) * (p_ch / p_e)**(1/ke) \
                * 1/sqrt((ke + 1) / (ke - 1) * (1 - (p_e / p_ch)**((ke - 1) / ke))) 
            A_e = self.A_e                                 # Exit area [m^2]
            ex = A_e / A_t                                 # Expansion ratio
        
            alpha_t = [14, 11, 10, 9]                      # Lookup table of divergence angles, assuming 80% bell length
            ex_t = [5, 10, 15, 20]                         # Lookup table of expansion ratios from alpha_t
            alpha= np.interp(ex, ex_t, alpha_t)
            lam = 0.5 * (1 + cos(alpha * pi/180))          # Thrust cosine loss correction,
            # even in extreme cases this is definitely not an O(1) effect 
            #lam=1                                        # could use this instead of previous 4 lines, idk how much it matters
            self.Ve = lam * sqrt(2*ke / (ke - 1) * Re * T_ch * \
                (1 - (p_e/p_ch)**((ke - 1)/ke)))           # Exhaust velocity       [m/s]
            self.f_init = True
        return mdot*self.Ve + (p_e - p_a)*self.A_e         # Thrust force, ignoring that isp increases w/ p_ch [N]

# this class is for conceptual clarity to keep rocket externals and internals organized
class Rocket(EngineSys):
    def __init__(self, L, dia, mdot, p_e, p_ch, T_ch, ke, Re, x_init):
        # rocket trajectory variables
        self.t = [0.]
        self.alt = [x_init]
        self.v = [0.]
        self.a = [0.]
        
        EngineSys.__init__(self, mdot, p_e, p_ch, T_ch, ke, Re)            # rocket's got an engine don't it
        self.L = L                                                         # propellant tubes length (m)
        self.dia = dia                                                     # rocket diameter (m)
        self.A = pi * (dia/2)**2                                           # Airframe frontal area projected onto a circle of diameter variable dia
        # relevant measurements
        self.m_prop = [sum(self.propellant_mass())]                        # propellant mass
        radius, l_o, l_f = openrkt.split_tanks(self.m_prop[0], dia/0.0254) # tank lengths, split_tanks expects in, not m
        self.total_length = shirt_l + openrkt.system_length(l_o, l_f)      # length of entire rocket
        m_eng_sys = openrkt.system_mass(radius, l_o, l_f)                  # dry mass of whole engine system
        self.m_dry = self.dry_mass(m_eng_sys)                              # Total dry mass of rocket
        self.m = [self.m_dry + self.m_prop[0]]                             # GLOW
        
    # Propellant masses
    def propellant_mass(self):
        A=self.A
        L=self.L
        OF=openrkt.OF
        rho_eth = openrkt.rho_eth           # Density, ethanol fuel [kg/m^3]
        rho_ipa = openrkt.rho_ipa           # Density of 64.8% IPA / 35.2% H20 [kg/m^3]
        rho_lox = openrkt.rho_lox           # Density, lox          [kg/m^3]
        
        L_lox = L / (1 + rho_lox / (rho_ipa*OF))
        m_lox = rho_lox * L_lox * A         # Oxidizer mass         [kg]
        m_alc = rho_ipa * (L - L_lox) * A   # Fuel mass             [kg]
        return m_lox, m_alc                 # Propellant Mass       [kg]    
    
    # Subsystem dry masses, needs sanity check from dirty ME's
    def dry_mass(self, engine_sys_mass):
        m_ringsclamps = (.466 + 1) * 7
        m_nosecone = 17                                                                     # nosecone weight    [kg]
        m_recovery = 4                                                                      # Recovery system mass [kg]
        m_payload = 4                                                                       # Payload mass         [kg]
        m_avionics = 3.3                                                                    # Avionics mass        [kg]
        m_n2      = 4                                                                       # mass of n2 tank [kg]
        m_airframe = 0.00125 * ((self.total_length - .23) *pi*self.dia) * openrkt.CF['rho'] # Airframe mass        [kg]
        m_fins = 6.2                                                                        # total fin mass [kg], estimate from openrocket
        return (m_ringsclamps + m_nosecone + m_recovery + m_payload + m_avionics + m_n2 \
            + engine_sys_mass + m_airframe + m_fins)

# this class is for conceptual clarity to keep our rocket externals organized
class Environment(Rocket):
    def __init__(self, L, dia, mdot, p_e, p_ch, T_ch, ke, Re, x_init):
        Rocket.__init__(self, L, dia, mdot, p_e, p_ch, T_ch, ke, Re, x_init)
        # Drag coefficient look up
        self.C_d_t = []
        self.Ma_t = []
        #csv_file = open('lv4.csv')     # use estimated drag data from openrocket for one of our lv4 designs
        csv_file = open('aerobee.csv') # Use aerobee 150 drag data, perhaps a pessimistic estimate, but has the right shape and approximate range
        drag_data = csv.reader(csv_file, delimiter=',')
        for row in drag_data:
            self.Ma_t.append(float(row[0]))
            self.C_d_t.append(float(row[1]))
        # air
        self.g = [self.gravity(self.alt[0])]
        self.ka = 1.4                   # Ratio of specific heats, air  
        self.Ra = 287.1                 # Avg. specific gas constant (dry air)
        p_a, rho, T_a = self.std_at(self.alt[0])
        self.p_a = [p_a]
        self.rho = [rho]
        self.T_a = [T_a]
        
        # drag
        self.D = [0.]
        self.q = [0.]
        self.Ma = [0.]
    # all your atmospheric needs are here, probably should sanity check, matches openrocket at least
    def std_at(self, h):                # U.S. 1976 Standard Atmosphere
        if h < 11000:
            T = 15.04 - 0.00649*h
            p = 101.29*((T + 273.1)/288.08)**5.256
    
        elif 11000 <= h and h <25000:
            T = -56.46
            p = 22.65*exp(1.73 - 0.000157*h)
    
        else:
            T = -131.21 + 0.00299*h
            p = 2.488 * ((T + 273.1)/216.6)**(-11.388)
    
        p_a = p*1000                 # Ambient air pressure [Pa]
        rho = p/(0.2869*(T + 273.1)) # Ambient air density [kg/m^3]
        T_a = T + 273.1              # Ambient air temperature [K]
        return p_a, rho, T_a      
    
    # calculates drag force
    ###CHECK ME
    def drag(self, x, v):
        C_d_t=self.C_d_t
        Ma_t=self.Ma_t
        A=self.A
        # Check Knudsen number and switch drag models (e.g. rarified gas dyn vs. quadratic drag)
        # I'm not sure what the above comment means for us, so I ignored it
        
        (p_a, rho, T_a) = self.std_at(x)
        
        #C_d_t = [0.15, 0.15, 0.3, 0.45, 0.25, 0.2, 0.175, .15, .15] # V2 rocket drag coefficient lookup table
        #Ma_t = [0, 0.6, 1.0, 1.1, 2, 3, 4, 5, 5.6]                  # V2 rocket Mach number lookup table
        
        Ma = v / sqrt(self.ka * self.Ra * T_a)

        C_d = np.interp(Ma, Ma_t, C_d_t)                            # Drag coefficient function
        q = 0.5 * rho * v**2                                        # Dynamic pressure [Pa]
        D = q * C_d * A                                             # Drag force       [N]
        return D, q, Ma
    
    # https://www.sensorsone.com/local-gravity-calculator/
    # International Gravity Formula IGF) 1980
    # Geodetic Reference System 1980 (GRS80)
    # Free Air Correction (FAC)
    # i don't necessarily trust this, but i will research another time. for now it's approximately ok
    def gravity(self, h):
        lat = radians(32) # degrees north, launch site
        igf = 9.780327 * (1 + 0.0053024 * sin(lat)**2 - 0.0000058 * sin(2*lat)**2 )
        fac = -3.086 * 10**(-6) * h
        return igf + fac

# initialize trajectory sim
# Note combustion gas properties ke, Re, T_ch, etc, determined from CEA
# CEA run with chamber pressure=350 psi, fuel temp=419.15 K, lox temp=90 K, OF=1.3 for fuel = 64.8% IPA / 35.2% H20
def init(L, mdot, dia, p_e, x_init=0, p_ch=350, T_ch=3153.08, ke=1.1241, Re=361.6088):
    p_ch, dia, p_e = unit_conv(p_ch, dia, p_e) # gently change units
    simulation = Environment(L, dia, mdot, p_e, p_ch, T_ch, ke, Re, x_init)
    
    # launch site
    simulation.off_rail = False
    simulation.rail_height = x_init + (launch_tower * 0.3048) # 60 ft, convert to m, scale for site altitude
    simulation.launch_speed = 0.                              # in case rocket fails to launch
    
    # engine cut-off
    simulation.motor_burnout = False
    
    # this might be useful one day if we want to (or can) think about throttling...
    simulation.air_is_irrelevant = False
    return simulation

# convenience function, hide the eyebleed
def step_forward(sim, next_step):
    sim.t.append(next_step.t[-1])
    sim.alt.append(next_step.alt[-1])
    sim.v.append(next_step.v[-1])
    sim.a.append(next_step.a[-1])
    sim.g.append(next_step.g[-1])
    sim.F.append(next_step.F[-1])
    sim.D.append(next_step.D[-1])
    sim.m.append(next_step.m[-1])
    sim.m_prop.append(next_step.m_prop[-1])
    sim.q.append(next_step.q[-1])
    sim.Ma.append(next_step.Ma[-1])
    sim.p_a.append(next_step.p_a[-1])
    sim.rho.append(next_step.rho[-1])
    sim.T_a.append(next_step.T_a[-1])
    return sim

# abstracting some of the logic of which regime the rocket is in
def phase_driver(sim, next_step):
    # Check if the propellant tanks are currently non-empty
    if sim.m_prop[-1] > 0:                                              # fuelled phase
        
        if (sim.off_rail == False) and (sim.alt[-1] >= sim.rail_height):  # check if we just left launch rail to capture current speed
            sim.launch_speed = sim.v[-1]
            sim.off_rail = True
        next_step.m_prop.append(sim.m_prop[-1] - sim.mdot*sim.dt)       # burn fuel
        next_step.F.append(sim.thrust(sim.alt[-1]))                       # calculate thrust at current height
    elif sim.m_prop[-1] <= 0:                                           # coasting phase
        if (sim.motor_burnout == False) and (len(sim.m_prop) >= 2) and (sim.m_prop[-2] > 0):       # no propellant now, did we have it last moment?
            sim.motor_burnout = True
            # here we may want to eventually calculate CoM and CoP to get stability
            # it would be annoying, but possibly useful
            #this is mostly just a hook for potential usage
        next_step.m_prop.append(0)
        next_step.F.append(0)
    next_step.m.append(sim.m_dry + sim.m_prop[-1])
    return sim, next_step

def air_exists(sim, next_step):
    if (sim.air_is_irrelevant == False) and (sim.p_a[-1] <= air_cutoff): # air used to matter, but does it not anymore?
        sim.air_is_irrelevant == True # this is just a hook for later usage ;)
    # update air properties
    p_a_next, rho_next, T_a_next = sim.std_at(sim.alt[-1])
    next_step.p_a.append(p_a_next)
    next_step.rho.append(rho_next)
    next_step.T_a.append(T_a_next)
    
    # calculate current drag force    
    D_next, q_next, Ma_next = sim.drag(sim.alt[-1], sim.v[-1])
    next_step.D.append(D_next)
    next_step.q.append(q_next)
    next_step.Ma.append(Ma_next)
    return sim, next_step

# calculate trajectory!! we would like this to be more sophisticated (4th order runge-kutta) but that is low priority
def do_physics(sim, next_step):
    next_step.t.append(sim.t[-1] + sim.dt)
    next_step.g.append(sim.gravity(sim.alt[-1]))
    next_step.a.append((sim.F[-1] - sim.D[-1]) / sim.m[-1] - sim.g[-1])
    next_step.v.append(sim.v[-1] + sim.dt*sim.a[-1])
    next_step.alt.append(sim.alt[-1] + sim.dt*sim.v[-1])
    return sim, next_step

# calculates trajectory of a design, this is important to have right
# one day, might be cool to allow for varying OF ratios but that would be VERY complicated
# dt is Time step [s]
def trajectory(L, mdot, dia, p_e, x_init=0, dt=.1):
    sim = init(L, mdot, dia, p_e, x_init)
    sim.dt = dt
    next_step = copy.deepcopy(sim) # this is to prevent interference between steps
    # now step our rocket through time until we reach apogee
    while (sim.v[-1] > 0) or (sim.off_rail == False):
        sim, next_step = phase_driver(sim, next_step)
        
        sim, next_step = air_exists(sim, next_step)

        sim, next_step = do_physics(sim, next_step)
        
        sim = step_forward(sim, next_step)
    
    # this will run after we reach apogee to do last-minute calculations and conversions
    sim.TWR = sim.a[1] / sim.g[1]          # Thrust-to-weight ratio constraint
    sim.ex = sim.A_e / sim.A_t        # expansion ratio
    sim.S_crit = sim.p_e / sim.p_a[0] # Sommerfield criterion constraint
    
    sim.alt = np.array(sim.alt)
    sim.v = np.array(sim.v)
    sim.a = np.array(sim.a)
    sim.F = np.array(sim.F)
    sim.D = np.array(sim.D)
    sim.q = np.array(sim.q)
    sim.p_a = np.array(sim.p_a)
    sim.massratio = (sim.m_prop[0] + sim.m_dry) / sim.m_dry # Mass ratio
    sim.dV1 = sim.Ve * log(sim.massratio) / 1000                # Tsiolkovsky's bane (delta-V)
    # ends our suffering and returns all the stats of the launch
    return sim
