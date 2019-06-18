# could improve:
# drag coefficients (barrowman??)
# gravity

from math import sqrt, pi, exp, log, cos, sin, radians
import numpy as np
import csv
import copy
import openrocket_interface as openrkt

# A numerical integration for rocket trajectories

# Physical constants
g_n = openrkt.g_n  # Gravitational acceleration [m/s^2]
launch_tower = 60. # launch rail height in ft
air_cutoff   = 25. # Pa, ambient pressure at which we no longer care about air

# upper subsystem dimensions
nose_l = 1.50        # m
ers_l  = 6 * 0.0254  # m (converted from in)
rcs_l  = 6 * 0.0254  # m (converted from in)
av_l   = 18 * 0.0254 # m (converted from in)
n2_l   = 18 * 0.0254 # m (converted from in)

shirt_l = sum([nose_l, ers_l, rcs_l, av_l, n2_l]) # m, total of above. (what about gaps??)

# helper function for unit conversions
def unit_conv(p_ch, dia, p_e):
    p_ch = p_ch*6894.76    # Chamber pressure, convert psi to Pa
    # LV4 design variables
    dia = dia*0.0254       # Convert in. to m
    p_e = p_e*1000         # Exit pressure, convert kPa to Pa
    return p_ch, dia, p_e

# this class is for conceptual clarity, to keep rocket internals organized
class EngineSys:
    # let the model keep track of its own properties
    def __init__(self, mdot, p_e, p_ch, T_ch, ke, Re):
        self.f_init = False
        self.mdot = mdot
        self.p_e = p_e
        self.p_ch = p_ch
        self.T_ch = T_ch
        self.ke = ke
        self.Re = Re
    
    # calculates engine thrust at a height
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
            
            self.ex = A_e / A_t                                 # Expansion ratio
            ex = self.ex
            
            alpha_t = [14, 11, 10, 9]                      # Lookup table of divergence angles, assuming 80% bell length
            ex_t = [5, 10, 15, 20]                         # Lookup table of expansion ratios from alpha_t
            alpha= np.interp(ex, ex_t, alpha_t)
            lam = 0.5 * (1 + cos(alpha * pi/180))          # Thrust cosine loss correction,
            # even in extreme cases this is definitely not an O(1) effect 
            #lam=1    # could use this instead of previous 4 lines, idk how much it matters
            
            self.Ve = lam * sqrt(2*ke / (ke - 1) * Re * T_ch * \
                (1 - (p_e/p_ch)**((ke - 1)/ke)))           # Exhaust velocity       [m/s]
            
            self.f_init = True # save time by only running this block once
        return mdot*self.Ve + (p_e - p_a)*self.A_e         # Thrust force, ignoring that isp increases w/ p_ch [N]

# this class is for conceptual clarity to keep rocket externals and internals organized
class Rocket(EngineSys):
    # let the model keep track of its own properties
    def __init__(self, L, dia, mdot, p_e, p_ch, T_ch, ke, Re, x_init):
        # rocket trajectory variables
        self.t = [0.]
        self.alt = [x_init]
        self.v = [0.]
        self.a = [0.]
        
        # rocket's got an engine don't it
        EngineSys.__init__(self, mdot, p_e, p_ch, T_ch, ke, Re) 
                   
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
        m_ringsclamps = (.466 + 1) * 7                    # weights of rings and clamps [kg]
        m_nosecone = 9.25                                   # nosecone weight [kg]
        m_recovery = 4                                    # Recovery system mass [kg]
        m_payload = 4                                     # Payload mass  [kg]
        m_avionics = 3.3                                  # Avionics mass  [kg]
        m_n2      = 4                                     # mass of n2 tank [kg]
        m_airframe = 0.00125 * openrkt.CF['rho'] * \
            ((self.total_length - 0.0) *pi*self.dia)      # Airframe mass [kg]
        m_fins = 8.35                                      # total fin mass [kg], estimate from openrocket
        return (m_ringsclamps + m_nosecone + m_recovery \
            + m_payload + m_avionics + m_n2 + \
            engine_sys_mass + m_airframe + m_fins)

# this class is for conceptual clarity to keep our rocket externals organized
class Environment(Rocket):
    # let the model keep track of its own properties
    def __init__(self, L, dia, mdot, p_e, p_ch, T_ch, ke, Re, x_init):
        # an environment has a rocket in it, ideally
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
        self.g = [g_n] #[self.gravity(self.alt[0])]
        self.ka = 1.4                   # Ratio of specific heats, air  
        self.Ra = 287.1                 # Avg. specific gas constant (dry air)
        p_a, rho, T_a = self.std_at(self.alt[0])
        
        self.p_a = [p_a]
        self.rho = [rho]
        self.T_a = [T_a]
        
        # drag and thrust
        self.D = [0.]
        self.q = [0.]
        self.Ma = [0.]
        self.F = [self.thrust(self.alt[0])]
        
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
        
        Ma = v / sqrt(self.ka * self.Ra * T_a)

        C_d = np.interp(Ma, Ma_t, C_d_t)             # Drag coefficient function
        q = 0.5 * rho * v**2                         # Dynamic pressure [Pa]
        D = q * C_d * A                              # Drag force       [N]
        return D, q, Ma
    
    # https://www.sensorsone.com/local-gravity-calculator/
    # International Gravity Formula IGF) 1980
    # Geodetic Reference System 1980 (GRS80)
    # Free Air Correction (FAC)
    def gravity(self, h):
        lat = radians(32) # degrees north, launch site
        igf = 9.780327 * (1 + 0.0053024 * sin(lat)**2 - 0.0000058 * sin(2*lat)**2 )
        fac = -3.086 * 10**(-6) * h
        return igf + fac

# initialize trajectory sim
# Note combustion gas properties ke, Re, T_ch, etc, determined from CEA
# CEA run with chamber pressure=350 psi, fuel temp=419.15 K, 
#              lox temp=90 K, OF=1.3 for fuel = 64.8% IPA / 35.2% H20
def init(L, mdot, dia, p_e, x_init=0, p_ch=350, T_ch=3153.08, ke=1.1241, Re=361.6088):
    p_ch, dia, p_e = unit_conv(p_ch, dia, p_e) # gently change units
    #initialize sim
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

# convenience function, hide the eye-bleed!
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
    if sim.m_prop[-1] > 0: # fuelled phase
        
        # check if we just left launch rail to capture current speed
        if (sim.off_rail == False) and (sim.alt[-1] >= sim.rail_height): 
            #print((sim.alt[-1]+sim.alt[-2])/2, (sim.v[-1]+sim.v[-2])/2) # debug
            # take average because otherwise it over/undershoots
            sim.launch_speed = (sim.v[-1]+sim.v[-2])/2
            sim.off_rail = True

        next_step.m_prop.append(sim.m_prop[-1] - sim.mdot*sim.dt) # burn fuel
        next_step.F.append(sim.thrust(sim.alt[-1]))               # calculate thrust at current height
    
    elif sim.m_prop[-1] <= 0: # coasting phase

        # no propellant now, did we have it last moment (if there was a last moment)?
        if (sim.motor_burnout == False) and (len(sim.m_prop) >= 2) and (sim.m_prop[-2] > 0):
            sim.F_index = len(sim.m_prop) # need to calculate burn time
            sim.motor_burnout = True
            # here we may want to eventually calculate CoM and CoP to get stability
            # it would be annoying, but possibly useful
            #this is mostly just a hook for potential usage

        next_step.m_prop.append(0) # no fuel
        next_step.F.append(0)      # no thrust

    next_step.m.append(sim.m_dry + sim.m_prop[-1]) # update mass
    return sim, next_step

# this is a helper function
def air_exists(sim, next_step):
    # air used to matter, but does it not anymore?
    if (sim.air_is_irrelevant == False) and (sim.p_a[-1] <= air_cutoff):
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

# helper func for runge kutta, this is our ODE :)
# v = dx/dt, a = dv/dt
def f(t, x, v, F, D, m, g):
    v = v
    a = (F - D)/m - g
    return np.array([v, a])

# this is kinda ugly-looking but should give much better efficiency
# truly a proper integrator for us sophisticated folk
# it would probably be most sane to update F, D, m, and g in the same way as x, v, a
def rungekutta(sim, next_step):
    k1 = sim.dt * f(sim.t[-1],
                    sim.alt[-1],
                    sim.v[-1],
                    sim.F[-1],
                    sim.D[-1],
                    sim.m[-1],
                    sim.g[-1])
    k2 = sim.dt * f(sim.t[-1] + sim.dt/2,
                    sim.alt[-1] + k1[0]*sim.dt/2,
                    sim.v[-1] + k1[1]*sim.dt/2,
                    sim.F[-1],
                    sim.D[-1],
                    sim.m[-1],
                    sim.g[-1])
    k3 = sim.dt * f(sim.t[-1] + sim.dt/2,
                    sim.alt[-1] + k2[0]*sim.dt/2,
                    sim.v[-1] + k2[1]*sim.dt/2,
                    sim.F[-1],
                    sim.D[-1],
                    sim.m[-1],
                    sim.g[-1])
    k4 = sim.dt * f(sim.t[-1] + sim.dt,
                    sim.alt[-1] + k3[0]*sim.dt,
                    sim.v[-1] + k3[1]*sim.dt,
                    sim.F[-1],
                    sim.D[-1],
                    sim.m[-1],
                    sim.g[-1])
    
    next_step.alt.append(sim.alt[-1] + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6 )
    next_step.v.append(sim.v[-1] + (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6 )
    return next_step

# deprecated function for basic monkey integration
def forwardeuler(sim, next_step):
    next_step.v.append(sim.v[-1] + sim.dt*sim.a[-1])
    next_step.alt.append(sim.alt[-1] + sim.dt*sim.v[-1])
    return next_step

# calculate trajectory!!
def do_physics(sim, next_step):
    next_step.t.append(sim.t[-1] + sim.dt)
    #next_step.g.append(g_n) # constant gravity is for noobs?
    next_step.g.append(sim.gravity(sim.alt[-1]))
    next_step.a.append((sim.F[-1] - sim.D[-1]) / sim.m[-1] - sim.g[-1])
    
    next_step = rungekutta(sim, next_step)
    return sim, next_step

# calculates trajectory of a design, this is important to have right
# one day, might be cool to allow for varying OF ratios but that would be VERY complicated
# and would require lots of dumb chemistry
# dt is Time step [s]
# I am suspicious of our state transition from step n to n+1,
# and need to take some time to think carefully about the succession of moments
def trajectory(L, mdot, dia, p_e, x_init=0, dt=.1):
    sim = init(L, mdot, dia, p_e, x_init) # get the puppy rolling
    sim.dt = dt
    next_step = copy.deepcopy(sim) # this is to prevent interference between steps
    
    # now step our rocket through time until we reach apogee
    while (sim.v[-1] > 0) or (sim.off_rail == False):
        sim, next_step = phase_driver(sim, next_step)
        
        sim, next_step = air_exists(sim, next_step)

        sim, next_step = do_physics(sim, next_step)
        
        sim = step_forward(sim, next_step)
    
    # DONE!
    # this will run after we reach apogee to do last-minute calculations and conversions
    sim.TWR = sim.a[1] / g_n                  # Thrust-to-weight ratio constraint
    sim.S_crit = sim.p_e / sim.p_a[0]         # Sommerfield criterion constraint
    sim.massratio = (sim.m_prop[0] + sim.m_dry) / sim.m_dry # Mass ratio
    sim.dV1 = sim.Ve * log(sim.massratio) / 1000            # Tsiolkovsky's bane (delta-V)
    
    # numpy arrays are nicer for some reasons
    sim.alt = np.array(sim.alt)
    sim.v = np.array(sim.v)
    sim.a = np.array(sim.a)
    sim.F = np.array(sim.F)
    sim.D = np.array(sim.D)
    sim.q = np.array(sim.q)
    sim.p_a = np.array(sim.p_a)
    
    # ends our suffering and returns all the stats of the launch
    return sim
