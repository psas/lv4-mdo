from math import sqrt, pi, exp, log, cos
import numpy as np

class sim: #A simple forward Euler integration for rocket trajectories
    def dry_mass(L, dia):
        m_avionics = 3.3
        m_recovery = 4
        m_payload = 2
        m_tankage = 20.88683068354522*L*dia*pi
        m_engine = 2
        m_feedsys = 20
        m_airframe  = 6
        return (m_avionics + m_recovery + m_payload + m_tankage 
        + m_engine + m_feedsys + m_airframe)
        
    def propellant_mass(A, L, OF=1.3):
        rho_alc = 852.3 #density of ethanol fuel [kg/m^3]
        rho_lox = 1141.0 #density of lox [kg/m^3]
        L_lox = L/(rho_lox/(rho_alc*OF) + 1)
        m_lox = rho_lox*L_lox*A #oxidizer mass [kg]
        m_alc = rho_alc*(L-L_lox)*A #fuel mass [kg]
        return m_alc + m_lox
    
    def std_at(h): #U.S. 1976 Standard Atmosphere
        if h < 11000:
            T = 15.04 - 0.00649*h
            p = 101.29*((T + 273.1)/288.08)**5.256
        
        elif 11000 <= h and h <25000:
            T = -56.46
            p = 22.65*exp(1.73 - 0.000157*h)
        
        else:
            T = -131.21 + 0.00299*h
            p = 2.488 * ((T + 273.1)/216.6)**(-11.388)
            
        rho = p/(0.2869*(T + 273.1)) #ambient air density
        p_a = p*1000 #ambient air pressure
        T_a = T + 273.1 #ambient air temperature
        return p_a, rho, T_a
        
    def thrust(x, p_ch, T_ch, p_e, ke, Re, mdot):
        p_a = sim.std_at(x)[0] 
        p_t = p_ch*(1 + (ke - 1)/2)**(-ke/(ke - 1)) #Throat pressure
        T_t = T_ch*(1/(1 + (ke - 1)/2)) #Throat temperature
        A_t = (mdot / p_t)*sqrt(Re*T_t/ke) #Throat area
        A_e = A_t*(2/(ke + 1))**(1/(ke - 1))*(p_ch/p_e)**(1/ke) * 1/sqrt((ke + 1)/(ke - 1)*(1 - (p_e/p_ch)**((ke - 1)/ke))) #Exit area
        ex = A_e/A_t
        alpha_t = [14, 11, 10, 9] #lookup table of divergence angles, assuming 80% bell
        ex_t = [5, 10, 15, 20]
        alpha= np.interp(ex, ex_t, alpha_t)
        lam = 0.5*(1 + cos(alpha *pi/180)) #thrust cosine loss correction, even in extreme cases this is definitely not an O(1) effect 
        Ve = lam*sqrt(2*ke/(ke - 1)*Re*T_ch*(1 - (p_e/p_ch)**((ke - 1)/ke))) #exhaust velocity
        F = mdot*Ve + (p_e - p_a)*A_e  #Thrust force, ignoring that isp increases w/ p_ch
        return F, A_t, A_e, Ve
    
    def drag(x, v, A, Ma):
        #check Knudsen number and switch drag models (e.g. rarified gas dyn vs. quadratic drag)
        (p_a, rho, T_a) = sim.std_at(x)
        C_d_t = [0.15, 0.15, 0.3, 0.45, 0.25, 0.2, 0.175, .15, .15] #super cheesy lookup table for drag coefficients (historical data for V-2)
        Ma_t = [0, 0.6, 1.0, 1.1, 2, 3, 4, 5, 5.6]
        C_d = np.interp(Ma, Ma_t, C_d_t)
        
        """
        drag coefficient, placeholder scalar constant should be a function of Mach number 
        e.g. run some simulations in Star and curve fit via Cd = a*Ma**b + C_d0
        """
        
        q = 0.5 * rho * v**2 #dyanmic pressure [Pa]
        D = q * C_d * A #drag force [N]
        return D, q
        
    def trajectory(L, mdot, dia, p_e, p_ch=350, T_ch=3500, ke=1.3, Re=349, x_init=0):
        #Note combustion gas properties ke, Re, T_ch, etc, determined from CEA
        
        #physical constants
        g_0 = 9.81 #gravitational acceleration
        dt = 1 #time step
        ka = 1.4 #Ratio of specific heats, air
        Ra = 287.1 #avg. specific gas constant (dry air)
    
        #LV design variables
        dia = dia*0.0254 #convert in to m
        A = pi*(dia/2)**2 #airframe frontal area projected onto a circle of raduis r
        m_dry = sim.dry_mass(L, A) #dry mass, call from function
        mdot = mdot #mass flow rate
        p_ch = p_ch*6894.76 #chamber pressure, convert psi to Pa
        p_e = p_e*1000 #exit pressure, convert kPa to Pa
    
    
        #initial conditions
        x = [x_init]
        v = [0]
        a = [0]
        t = [0]
        rho = [sim.std_at(x[-1])[1]]
        p_a = [sim.std_at(x[-1])[0]]
        T_a = [sim.std_at(x[-1])[2]]
        m_prop = [sim.propellant_mass(A, L)]
        m = [m_dry + m_prop[-1]]
        (F, A_t, A_e, Ve) = sim.thrust(x[-1], p_ch, T_ch, p_e, ke, Re, mdot)
        F = [F]
        D = [0]
        Ma = [0]
        q = [0]
        r = (m_prop[0] + m_dry)/m_dry #mass ratio
        dV1 = Ve*log(r)/1000 #Tsiolkovsky's bane
        
        while True:
            p_a.append(sim.std_at(x[-1])[0])
            rho.append(sim.std_at(x[-1])[1])
            T_a.append(sim.std_at(x[-1])[2])
            if m_prop[-1] > 0:
                (Fr, A_t, A_e, Ve) = sim.thrust(x[-1], p_ch, T_ch, p_e, ke, Re, mdot)
                F.append(Fr)
                m_prop.append(m_prop[-1] - mdot*dt)
                mdot_old = mdot
            else:
                Ve = sim.thrust(x[-1], p_ch, T_ch, p_e, ke, Re, mdot_old)[3]
                F.append(0)
                mdot = 0
                m_prop[-1] = 0
            q.append(sim.drag(x[-1], v[-1], A, Ma[-1])[1])
            D.append(sim.drag(x[-1], v[-1], A, Ma[-1])[0])
            a.append((F[-1] - D[-1])/m[-1] - g_0)
            v.append(a[-1]*dt + v[-1])
            x.append(v[-1]*dt + x[-1]) 
            Ma.append(v[-1]/sqrt(ka*Ra*T_a[-1]))
            t.append(t[-1] + dt)
            m.append(m_dry + m_prop[-1])
            TWR = a[1]/g_0 #constraint
            ex = A_e/A_t
            S_crit = p_e/p_a[0] #constraint
            if v[-1] <= 0:
                x = np.array(x)
                a = np.array(a)
                F = np.array(F)
                D = np.array(D)
                q = np.array(q)
                return L, dia, x, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop
 