# could improve:
# dry mass estimates
# thrust, and drag coefficients
# trajectory realism

from math import sqrt, pi, exp, log, cos
import numpy as np
import csv

# A simple forward Euler integration for rocket trajectories
# could range-kutta or some other integration, but this is sufficient for now

# Subsystem masses, needs sanity check from dirty ME's
def dry_mass(tankmass):
    m_nosecone = 30 + 0.1                        # nosecone weight    [kg]
    m_recovery = 4 + 0.1 + 1.2 + 0.1                         # Recovery system mass [kg]
    m_payload = 0.05 + 0.1                          # Payload mass         [kg]
    m_avionics = 3 + 0.1                       # Avionics mass        [kg]
    m_tankage = tankmass # Tank mass Estimation [kg] # avg 30 kg
    m_engine = 3 + 0.1 + 0.1                         # Engine mass          [kg] #CHECK THIS
    m_feedsys = 5 + 10                         # Feed system mass     [kg] #CHECK THIS, i'm including plumbing
    m_airframe  = 4.3 + 0.99 + 0.99 + 0.45 + 0.99 + 3.6 + 2.4 + 1.5                # Airframe mass        [kg]
    m_fins = 4.2                              # total fin mass [kg] #estimate from openrocket
    return (m_nosecone + m_recovery + m_payload + m_avionics + m_tankage 
        + m_engine + m_feedsys + m_airframe + m_fins)   # total Dry mass

# Propellant masses
def propellant_mass(A, L, OF=1.3):
    rho_alc = 852.3             # Density, ethanol fuel [kg/m^3]
    rho_ipa = 786     # kg/m^3 Density of Isopropyl Alcohol
    rho_lox = 1141.0            # Density, lox          [kg/m^3]
    
    L_lox = L/(rho_lox/(rho_ipa*OF) + 1)
    m_lox = rho_lox*L_lox*A     # Oxidizer mass         [kg]
    m_alc = rho_ipa*(L-L_lox)*A # Fuel mass             [kg]
    return m_lox, m_alc         # Propellant Mass       [kg]

# all your atmospheric needs are here, probably should sanity check
def std_at(h):                  # U.S. 1976 Standard Atmosphere
    if h < 11000:
        T = 15.04 - 0.00649*h
        p = 101.29*((T + 273.1)/288.08)**5.256

    elif 11000 <= h and h <25000:
        T = -56.46
        p = 22.65*exp(1.73 - 0.000157*h)

    else:
        T = -131.21 + 0.00299*h
        p = 2.488 * ((T + 273.1)/216.6)**(-11.388)

    rho = p/(0.2869*(T + 273.1)) # Ambient air density [kg/m^3]
    p_a = p*1000                 # Ambient air pressure [Pa]
    T_a = T + 273.1              # Ambient air temperature [K]
    return p_a, rho, T_a

# calculates engine thrust
def thrust(x, p_ch, T_ch, p_e, ke, Re, mdot):
    if p_e < 0.0001: #this is a kludge because optimizer divided by 0 once
        p_e = 0.01
    
    p_a = std_at(x)[0]                          # Ambient air pressure [Pa]
    p_t = p_ch*(1 + (ke - 1)/2)**(-ke/(ke - 1)) # Throat pressure      [Pa]
    T_t = T_ch*(1/(1 + (ke - 1)/2))             # Throat temperature   [K]
    A_t = (mdot / p_t)*sqrt(Re*T_t/ke)          # Throat area          [m^2]
    
    A_e = A_t*(2/(ke + 1))**(1/(ke - 1))*(p_ch/p_e)**(1/ke) * 1/sqrt((ke + 1)/(ke - 1)*(1 - (p_e/p_ch)**((ke - 1)/ke))) # Exit area [m^2]
    ex = A_e/A_t              # Expansion ratio
    
    alpha_t = [14, 11, 10, 9] # Lookup table of divergence angles, assuming 80% bell length
    ex_t = [5, 10, 15, 20]    # Lookup table of expansion ratios from alpha_t
    alpha= np.interp(ex, ex_t, alpha_t)
    lam = 0.5*(1 + cos(alpha *pi/180)) # Thrust cosine loss correction, even in extreme cases this is definitely not an O(1) effect 
    #lam=1 # could use this instead of previous 4 lines, idk how much it matters
    
    Ve = lam*sqrt(2*ke/(ke - 1)*Re*T_ch*(1 - (p_e/p_ch)**((ke - 1)/ke))) # Exhaust velocity                                  [m/s]
    F = mdot*Ve + (p_e - p_a)*A_e                                        # Thrust force, ignoring that isp increases w/ p_ch [N]
    return F, A_t, A_e, Ve

# calculates drag force
def drag(x, v, A, Ma, C_d_t, Ma_t):
    # Check Knudsen number and switch drag models (e.g. rarified gas dyn vs. quadratic drag)
    # I'm not sure what the above comment means for us, so I ignored it
    
    (p_a, rho, T_a) = std_at(x)
    
    #C_d_t = [0.15, 0.15, 0.3, 0.45, 0.25, 0.2, 0.175, .15, .15] # V2 rocket drag coefficient lookup table
    #Ma_t = [0, 0.6, 1.0, 1.1, 2, 3, 4, 5, 5.6]                  # V2 rocket Mach number lookup table
    
    C_d = np.interp(Ma, Ma_t, C_d_t)                            # Drag coefficient function
    q = 0.5 * rho * v**2                                        # Dynamic pressure [Pa]
    D = q * C_d * A                                             # Drag force       [N]
    return D, q

# calculates trajectory of a design, this is important to have right
# at some point, sanity check the constants we're assuming
# the default values of ke and Re in the function definition are for ethanol, theoretically
def trajectory(L, mdot, dia, p_e, p_ch=350, T_ch=3500, ke=1.3, Re=349, x_init=0, dt=.1, tankmass=30., propmass=0):
    # Note combustion gas properties ke, Re, T_ch, etc, determined from CEA
    # dt is Time step                  [s]
    ke = 1.1477 # from CEA for IPA, gammas. OF = 1.3, T = 3500
    Re = 472.7355 # from CEA for IPA, R/ molar weight (M)
    # Physical constants
    g_0 = 9.80665 # Gravitational acceleration [m/s^2]
    ka = 1.4   # Ratio of specific heats, air  
    Ra = 287.1 # Avg. specific gas constant (dry air)
    
    # Design constant
    p_ch = p_ch*6894.76    # Chamber pressure, convert psi to Pa
    
    # launch site constants
    off_rail = False
    rail_height = 60. * 0.3048 # 60 ft, convert to m
    launch_speed = 0. # kludge because optimizer crashed once
    
    # engine cut-off variables
    motor_burnout = False
    
    # LV4 design variables
    L = L                  # length of both tanks of propellant, in m
    mdot = mdot            # Mass flow rate [kg/s]
    dia = dia*0.0254       # Convert in. to m, sorry
    p_e = p_e*1000         # Exit pressure, convert kPa to Pa
    
    # relevant measurements
    A = pi*(dia/2)**2      # Airframe frontal area projected onto a circle of diameter variable dia
    m_dry = dry_mass(tankmass) # Dry mass, call from function dry_mass()

    # Initial conditions
    # rocket
    x = [x_init]
    v = [0.]
    a = [0.]
    t = [0.]
    
    # air
    rho = [std_at(x[-1])[1]]
    p_a = [std_at(x[-1])[0]]
    T_a = [std_at(x[-1])[2]]
    
    # propellant mass, else option is to dictate an amount, don't know if that's helpful
    if propmass==0:
        m_prop = [sum(propellant_mass(A, L))]
    else:
        m_prop = [propmass]
    m = [m_dry + m_prop[-1]] # GLOW
    
    (F, A_t, A_e, Ve) = thrust(x[-1], p_ch, T_ch, p_e, ke, Re, mdot)
    F = [F]
    D = [0.]
    Ma = [0.]
    q = [0.]
    
    r = (m_prop[0] + m_dry)/m_dry # Mass ratio
    dV1 = Ve*log(r)/1000          # Tsiolkovsky's bane (delta-V)
    mdot_old = mdot
    
    # Drag coefficient look up
    C_d_t = []
    Ma_t = []
    f = open('CD_sustainer_poweron.csv') # Use aerobee 150 drag data, perhaps a pessimistic estimate, but has the right shape and approximate range
    aerobee_cd_data = csv.reader(f, delimiter=',')
    for row in aerobee_cd_data:
        C_d_t.append(float(row[1]))
        Ma_t.append(float(row[0]))

    # now step our rocket through time until apogee
    while True:
        # air
        p_a.append(std_at(x[-1])[0])
        rho.append(std_at(x[-1])[1])
        T_a.append(std_at(x[-1])[2])
        
        # Check if the propellant tanks are empty
        if m_prop[-1] > 0: # fuelled phase
            (Fr, A_t, A_e, Ve) = thrust(x[-1], p_ch, T_ch, p_e, ke, Re, mdot)
            F.append(Fr)
            m_prop.append(m_prop[-1] - mdot*dt)
            mdot_old = mdot
            
            # check if we left launch rail
            if (off_rail == False) and (x[-1] >= rail_height):
                launch_speed = v[-1]
                off_rail = True
                
        else: # coasting phase
            if (motor_burnout == False) and (m_prop[-2] > 0): # no propellant now, did we have it last moment?
                motor_burnout = True
                # here we want to eventually calculate CoM and CoP to get stability
                # we will spit stability out of traj func
                #this is mostly just a hook for later usage
            Ve = thrust(x[-1], p_ch, T_ch, p_e, ke, Re, mdot_old)[3]
            F.append(0)
            mdot = 0
            m_prop[-1] = 0

        # kind of an simple method but it works
        # could improve significantly, but not essential
        q.append(drag(x[-1], v[-1], A, Ma[-1], C_d_t, Ma_t)[1])
        D.append(drag(x[-1], v[-1], A, Ma[-1], C_d_t, Ma_t)[0])
        a.append((F[-1] - D[-1])/m[-1] - g_0)
        v.append(a[-1]*dt + v[-1])
        x.append(v[-1]*dt + x[-1]) 
        Ma.append(v[-1]/sqrt(ka*Ra*T_a[-1]))
        t.append(t[-1] + dt)
        m.append(m_dry + m_prop[-1])
        TWR = a[1]/g_0      # Thrust-to-weight ratio constraint
        ex = A_e/A_t
        S_crit = p_e/p_a[0] # Sommerfield criterion constraint
        
        # check to see if we reached apogee
        if v[-1] <= 0:
            x = np.array(x)
            a = np.array(a)
            F = np.array(F)
            D = np.array(D)
            q = np.array(q)
            # ends our suffering and returns all the stats of the launch
            return x, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop, launch_speed
