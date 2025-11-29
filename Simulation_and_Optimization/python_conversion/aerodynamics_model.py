
from system_definition import *
import scipy.interpolate
import numpy as np

project_normal = np.identity(3) - np.outer(np.array([0,0,1]), np.array([0,0,1]))

class AeroModel():
    def __init__(self, diameter, body_len, nose_len, fin):
        # kinematic viscosity of air
        self.mu                     = 1.5 * 10**-5 # m^2/s
        #self.surface_roughness      = 2 * 10**-6 # m, perhaps too ideal, from openrocket
        #self.surface_roughness      = 1.79e-4 # m, untreated from lv3
        self.surface_roughness      = 3.3e-6 # m, treated from lv3
        self.reference_len          = diameter
        self.body_len               = body_len
        self.nose_len               = nose_len
        self.A_ref                  = np.pi * self.reference_len**2 / 4
        self.height                 = fin.semispan # fin semi-span
        self.root                   = fin.root
        self.fin_lead               = body_len - self.root
        self.tip                    = fin.tip
        self.thickness              = fin.thickness
        self.sweep_angle            = fin.sweep_angle
        self.fin_area               = fin.volume / self.thickness
        self.MAC_length             = 2/3 * (self.root**2 + self.tip**2 + self.root * self.tip) / (self.root + self.tip)
        self.MACspan                = self.MAC_span()
        self.MAC_lead               = self.effective_MAC_LE()
        self.finAR                  = self.fin_aspect()
        # Commercial Airplane Design Principles, eq 5.11 for midchord sweep
        self.midchord_sweep         = np.arctan(self.sweep_angle - (2 / self.finAR) * (1 - self.tip / self.root) / (1 + self.tip / self.root))
        self.fin_poly               = self.coefficients(self.finAR)
        self.fin_pressure_mult      = np.cos(self.sweep_angle)**2 * 4 * self.thickness * self.height / self.A_ref
        self.bodytail               = self.bodytail_interference()
        self.finshape               = self.roll_damp_fin_shape()
        self.fin_interps            = self.init_interpolators()
        self.CDfric_coef            = self.frictcoef()
        self.averagechord           = (self.root + self.tip) / 2
        self.aspectratio            = self.height / self.averagechord
        self.partial_flutter_factor = 39.3 * (self.aspectratio * self.averagechord /self.thickness)**3 / (self.aspectratio + 2)
        self.taper_ratio            = self.tip / self.root
        self.flutter_factor         = self.partial_flutter_factor * (self.taper_ratio + 1) / (2 * 101325)
        self.tumbling_drag          = 0.5 * (1.2 + (1.42 * 1.41 * self.fin_area +
                                       0.56 * self.reference_len * self.body_len) / self.A_ref) # eq 3.99 of OR tech doc
        self.ka                     = 1.4                   # Ratio of specific heats, air  
        self.Ra                     = 287.058                 # Avg. specific gas constant (dry air)
    
    # this was an attempt to put constraints on fin designs, can't remember where I found it...
    def flutter_velocity(self, sound_speed, p_a):
        G_E = 2.6e10 # effective shear modulus, alum, pascals
        return sound_speed * np.sqrt(G_E / (self.flutter_factor*p_a))

    def fin_aspect(self):
        return 2*self.height**2/self.fin_area # multiplying by 2 because openrocket does, technically that may be wrong
    
    # fin-body interference coefficient, from section 5-3.4.2 of MIL-HDBK-762
    def bodytail_interference(self):
        b = (2 * self.height) + self.reference_len
        d = self.reference_len
        r = self.reference_len / 2
        return ((2/np.pi)/((1 - d / b)**2)) * ((1+(d**4 / b**4))*(0.5 * np.arctan(0.5*(b/r - d/b)) + np.pi/4) - (d**2 / b**2)*((b/d - d/b) + 2*np.arctan(d/b)))

    # subsonic fin shape term, eq 3.70
    def roll_damp_fin_shape(self):
        radius = self.reference_len / 2
        term1 = (self.root + self.tip) * self.height * radius**2 / 2
        term2 = (self.root + 2 * self.tip) * self.height**2 * radius/3
        term3 = (self.root + 3 * self.tip) * self.height**3 / 12
        return term1 + term2 + term3

    # skin friction drag coefficient, assuming nose cone is cylinder for now
    # func of body fineness, fin thickness, mean aero chord, wet area of body, wet area of fins (both sides incl)
    def frictcoef(self):
        return ((1 + 1/(2*(self.body_len/self.reference_len)))* (np.pi * self.reference_len * self.body_len) +
                      (1 + 2*self.thickness/self.MAC_length)* ((2*self.fin_area+self.tip*self.thickness)*4)) / self.A_ref
    
    # reciprocal of Prandtl factor P "which corrects subsonic force coefficients for compressible flow"
    def Beta(self, mach):
        return np.sqrt(abs(mach**2 - 1))
    
    # v0 is the free-stream velocity of the rocket, mu is the kinematic viscosity of air
    def Reynold(self, v0):
        return self.body_len*v0/self.mu

    #fin interpolator, yikes copypasta
    def coefficients(self, ar):
        denom = (1- 3.4641*ar)**2
        poly=[]
        poly.append((9.16049 * (-0.588838 + ar) * (-0.20624 + ar)) / denom)
        poly.append((-31.6049 * (-0.705375 + ar) * (-0.198476 + ar)) / denom)
        poly.append((55.3086 * (-0.711482 + ar) * (-0.196772 + ar)) / denom)
        poly.append((-39.5062 * (-0.72074 + ar) * (-0.194245 + ar)) / denom)
        poly.append((12.8395 * (-0.725688 + ar) * (-0.19292 + ar)) / denom)
        poly.append((-1.58025 * (-0.728769 + ar) * (-0.192105 + ar)) / denom)
        return poly

    # coefficients for single fin supersonic normal force coefficient derivative
    # we will follow their approach and precompute some coefficients
    # and we will return linear interpolators 
    # as you can see, I shoved a bunch of interpolators in here. pleaaase organize.
    def init_interpolators(self):
        gamma = 1.4 # iirc this was an arbitrary parameter that was picked to fit a curve or something
        Ma0 = 1.5
        def K1(beta):
            return 2.0/beta
        def K2(beta, mach):
            return ((gamma + 1) * mach**4 - 4 * beta**2) / (4 * beta**4)
        def K3(beta, mach):
            return ((gamma + 1) * mach**8 + (
                    2 * gamma**2 - 7 * gamma - 5) * mach**6
                    + 10 * (gamma + 1) * mach**4 + 8) / (6 * beta**7)
        Ma, k1, k2, k3 = [Ma0], [], [], []
        while Ma[-1] < 5.0:
            beta = self.Beta(Ma[-1])
            k1.append(K1(beta))
            k2.append(K2(beta, Ma[-1]))
            k3.append(K3(beta, Ma[-1]))
            Ma.append(Ma[-1] + 0.1)
        trash = Ma.pop()
        k1_intrp = lambda M: np.interp(M, Ma, k1)
        k2_intrp = lambda M: np.interp(M, Ma, k2)
        k3_intrp = lambda M: np.interp(M, Ma, k3)
        
        # sneaking a second interpolator creation in here
        self.CD_interp_angle = 17 * np.pi / 180
        start = (0,1)
        end = (self.CD_interp_angle, 1.3)
        self.CD_interp_below = scipy.interpolate.CubicHermiteSpline([start[0], end[0]],
                                                 [start[1], end[1]],
                                                 [0, 0],
                                                 extrapolate=False)
        start = (self.CD_interp_angle, 1.3)
        end = (np.pi/2, 0)
        self.CD_interp_above = scipy.interpolate.CubicHermiteSpline([start[0], end[0]],
                                                 [start[1], end[1]],
                                                 [0, 0],
                                                 extrapolate=False)
        
        # sneaking some values for a third interpolator here
        self.subsonic = 0.9
        self.supersonic = 1.5
        self.supersonic_b = (self.supersonic**2 - 1)**1.5
        self.sq = np.sqrt(1 +
                          (1 - self.subsonic**2) * (self.height**2 /
                                                    (self.fin_area * np.cos(self.midchord_sweep)))**2)
        self.subV = 2 * np.pi * self.height**2 / self.A_ref / (1 + self.sq)
        self.subD_coeff = 2 * np.pi * self.height**6 / ((self.fin_area * np.cos(self.midchord_sweep))**2 
                                                     * self.A_ref * self.sq * (1 + self.sq)**2)
        self.superD = - self.fin_area/self.A_ref * 2 * self.supersonic / self.supersonic_b
        self.superV = None
        
        # sneaking in a 4th interpolator
        self.CN_a1_sub = lambda b: 2*np.pi*self.height**2/self.A_ref/(1+np.sqrt(
                    1 + (self.height**2*b/(self.fin_area*np.cos(self.midchord_sweep)))**2))
        
        # oh, a 5th and 6th interpolator. for von Karman/LD-Haack nosecone pressure drag as function of mach
        # https://github.com/openrocket/openrocket/blob/unstable/core/src/net/sf/openrocket/aerodynamics/barrowman/SymmetricComponentCalc.java
        self.nose_extrap = np.log(self.nose_len / self.A_ref + 1) / np.log(4)
        self.vonKarmanInterpolator = scipy.interpolate.interp1d(
                                    [0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 2.0, 3.0],
                                    [0, 0.010, 0.027, 0.055, 0.070, 0.081, 0.095, 0.097, 0.091, 0.083],
                                    kind='linear', copy=False, fill_value="extrapolate")
        # this next block is maybe unnessary, just look at the first value
        # only leaving it for now in case we wanna different nosecone
        minValue = self.vonKarmanInterpolator.__call__(0.9)
        if minValue > 0.001:
            minDeriv = (self.vonKarmanInterpolator.__call__(0.9 + 0.01) - minValue) / 0.01
            a = minValue - 0.8 * np.sin(0.00)**2
            b = minDeriv / a
            self.nose_cone_pressure_subsonic = lambda M: a * M**b + 0.8 * np.sin(0.00)**2
        else:
            self.nose_cone_pressure_subsonic = lambda M: 0
        return (k1_intrp, k1_intrp, k1_intrp)

    # assuming CP constant, fin for mach < 0.5
    # for trapezoid fins with tip parallel to root
    def MAC_span(self):
        return (self.height/3) * (self.root + 2*self.tip)/(self.root + self.tip)
    def effective_MAC_LE(self):
        def Z_t():
            return np.tan(self.sweep_angle)*self.height
        return (Z_t()/3)*(self.root + 2*self.tip)/(self.root + self.tip)

    # correction for fin CP changing with speed
    # for mach above 2, and below is interpolated
    # by a 5th order overfit polynomial lol. RHS = right hand side of equation
    def fin_CoP(self, mach):
        if mach <=0.7: #subsonic, 1/4 of MAC
            RHS = 0.25
        elif mach >= 2: #supersonic, 1/2 of mac
            beta = self.Beta(mach)
            RHS = (self.finAR*beta-0.67)/(2*self.finAR*beta-1)
        else: #transonic, oh god
            x=1.0
            RHS=0
            for coeff in self.fin_poly:
                RHS += coeff*x
                x *= mach
        return self.fin_lead + self.MAC_lead + self.MAC_length*RHS

    #Tactical Missile Design by Fleeman, ignore variance with mach
    # this is from slender body and crossflow theories
    def bodynose_CoP(self, alpha):
        sin2a = np.sin(alpha)**2
        return 0.63*self.nose_len*(1-sin2a) + 0.5*self.body_len*sin2a
    
    # pitch damping coefficient, from openrocket    
    def pitch_damping_coefficient(self, pitch, yaw, v0, c_o_g):
        top_length  = self.body_len - c_o_g
        bot_length  = c_o_g
        l           = top_length**4 + bot_length**4
        fin_length  = c_o_g - self.root + self.MAC_lead + self.MAC_length*0.25
        body_moment = 0.275 * l * self.reference_len
        fin_moments = 16 * 0.6 * self.fin_area * fin_length**3
        pitch_term  = 0 if v0 == 0 else (pitch / v0)**2
        yaw_term    = 0 if v0 == 0 else (yaw / v0)**2
        mult        = (body_moment + fin_moments) / (self.A_ref * self.reference_len)
        return mult * pitch_term, mult * yaw_term
    
    # robert galejs' correction term for cylindrical body section lift
    # reference area is cross-section, planform is diam*length
    # technically this comes from pg 3-11 of Fluid Dynamic Drag, Hoerner
    def C_N_alpha_lift(self, alpha):
        K = 1.1 # fudge factor
        return K*(self.body_len-self.nose_len)*self.reference_len/self.A_ref * np.sin(alpha) * np.sinc(alpha/np.pi)

    # nosecone normal force coefficients
    def C_N_alpha_nose(self, alpha):
        return 2*np.sinc(alpha/np.pi)

    def normal_force_coefficient(self, alpha, mach):
        def C_N_alpha_fins(alpha, mach):
            # single fin subsonic normal force coefficient derivative
            def C_N_alpha_1_sub(beta):
                return self.CN_a1_sub(beta)

            # single fin supersonic normal force coefficient derivative
            # Barrowman used a third-order (!) expansion of Busemann theory
            # openrocket simplified his method by assuming fins are flat plates
            # and ignoring correction of fin-tip mach cones
            def C_N_alpha_1_super(alpha, mach):
                return self.fin_area/self.A_ref * (self.fin_interps[0](mach)
                                      + self.fin_interps[1](mach) * alpha
                                      + self.fin_interps[2](mach) * alpha**2)

            # black magic, don't ask. single fin transonic
            # https://github.com/openrocket/openrocket/blob/519379cf56ff7dcbcce237bf3ccf324208e8ae42/core/src/net/sf/openrocket/aerodynamics/barrowman/FinSetCalc.java#L423
            def C_N_alpha_1_trans(alpha, mach):
                if not self.superV:
                    self.superV = lambda a: C_N_alpha_1_super(a, self.supersonic)
                # cubic hermite spline sometimes led to infinite values
                #trans_interp = scipy.interpolate.CubicHermiteSpline([self.subsonic, self.supersonic],
                #                                     [self.subV, self.superV(alpha)],
                #                                     [self.subD_coeff*mach, self.superD])
                # this seems to be an appropriate fix so far, within 5% usually. don't know which is correct
                # so might be worth comparing to the exact polynomial construction from link above
                trans_interp = scipy.interpolate.BPoly.from_derivatives([self.subsonic, self.supersonic],
                                                     [[self.subV, self.subD_coeff*mach,0],
                                                     [self.superV(alpha), self.superD,0]])
                return trans_interp.__call__(mach)
            
            if mach <= self.subsonic:
                CNa1 = C_N_alpha_1_sub(self.Beta(mach))
            elif mach > self.supersonic:
                CNa1 = C_N_alpha_1_super(alpha, mach)
            else: #cry
                CNa1 = C_N_alpha_1_trans(alpha, mach)
            return CNa1 * 2

        #fin-body interference
        def C_N_alpha_tail(CNa_fins):
            return CNa_fins * self.bodytail
        
        return C_N_alpha_tail(C_N_alpha_fins(alpha, mach))        

    def roll_damping_coefficient(self, omega, v0, Mach):
        # subsonic roll damping from 4 fins, eq 3.69
        # v0 is free-stream velocity, omega is roll rate
        def C_ld_subsonic(Mach, v0, omega):
            CNa0 = 2*np.pi / self.Beta(Mach)
            v0 = v0 if v0 != 0 else 1
            return 4*CNa0*omega*self.finshape/(self.A_ref*self.reference_len*v0)

        def C_ld_supersonic(mach, omega, v0):
            # local pressure coefficient of a strip of fin
            def C_P_i(eta, K):
                return (K[0] * eta
                        + K[1] * eta**2
                        + K[2] * eta**3)
            def eta_i(distance, coeff):# eq 3.62, small angle approx
                return distance*coeff # reasonable up to 20 deg
            def width(y): #trapezoid, y is height
                return y*self.tip/self.height + (self.height-y)*self.root/self.height

            vel_coeff = omega/v0 if v0 != 0 else omega
            eval_K = [k(mach) for k in self.fin_interps]
            mult = 4/(self.A_ref*self.reference_len)

            num_terms = 5
            del_y = self.height/num_terms

            terms = sum([C_P_i(eta_i(
                                    i*del_y, vel_coeff),
                               eval_K) *i*del_y * width(i*del_y) * del_y
                for i in range(num_terms+1)])
            return mult*terms
        if omega == 0: return 0
        if Mach <= 1:
            return C_ld_subsonic(Mach, v0, omega)
        else:
            return C_ld_supersonic(Mach, omega, v0)

    # drag coefficient at zero angle of attack
    def C_D_0(self, R, M):
        # skin friction coefficient function of reynolds number and mach
        # turbulent flow, in BarrowmanCalculator.java there is a continuous interpolation of correction around Mach 1
        # but that probably isn't needed here. could precalc a little here?
        def C_f(R, M):
            correction1, correction2 = 1, 1
            if M < 1.1:# compressibility correction
                correction1 = 1 - 0.1*M**2
            if M > 0.9:
                correction2 = (1 + 0.15 * M**2)**-0.58

            if M < 0.9:
                correction = correction1
            elif M < 1.1:
                correction = correction2*(M-.9)/.2 + correction1*(1.1-M)/.2
            else:
                correction = correction2

            if M < 0.9:# compressibility correction
                roughcorrection = 1 - 0.1*M**2
            elif M > 1.1:
                roughcorrection = 1/(1 + 0.18 * M**2)
            else:
                c1 = 1 - 0.1*M**2
                c2 = 1/(1 + 0.18 * M**2)
                roughcorrection = c2 * (M - 0.9) / 0.2 + c1 * (1.1 - M) / 0.2

            if R < 10**4: # velocity < 1 m/s
                Cf = 1.48 * 10**-2
            else: #assuming smooth surface
                Cf = correction * (1.5* np.log(R) - 5.6)**-2

            roughlimited = roughcorrection*0.032*(self.surface_roughness / self.body_len)**0.2
            Cf = max(Cf, roughlimited)
            return Cf

        # skin friction drag coefficient
        # func of body fineness, fin thickness, mean aero chord, wet area of body, wet area of fins (both sides incl)
        def C_D_frict(Cfc):
            return Cfc * self.CDfric_coef

        def stagnation_pressure(M): # eq B.1 and B.2, blunt cylinder stagnation pressure
                if M <= 1:
                    return 0.85 * (1 + M**2 / 4 + M**4 / 40)
                else:
                    return 0.85 * (1.84 - 0.76 / M**2 + 0.166 / M**4 + 0.035 / M**6)

        # from openrocket tech doc, eqs 3.95-6
        # assuming 3 tubes and a button
        def parasitic_drag(stag_p):
            pipe_ref_A       = 0.000506707 # assuming half inch radius tube
            launch_lug_ref_A = 0.004 # m^2
            CD_pipes         = 3 * stag_p * pipe_ref_A / self.A_ref
            CD_lug           = stag_p * launch_lug_ref_A / self.A_ref
            return CD_pipes, CD_lug

        # nosecone pressure drag coefficient
        # assumes nosecone smoothly transitions
        def C_D_nose_press(M, stag_p):
            if M < 0.9:
                return self.nose_cone_pressure_subsonic(M) # this should be right
            else:
                # ref eq B.9
                # the two if blocks are a hack because this occasionally returns an invalid value, mysteriously
                if stag_p <= 0: return 0
                term = self.vonKarmanInterpolator.__call__(M) / stag_p
                if term <= 0: return 0
                return stag_p * (term)**self.nose_extrap

        # fin pressure drag coefficient, eq 3.89, 3.91
        # assumes leading edge rounded and trailing edge tapered
        def C_D_fin_press(M):
            if M < 0.9:
                LE = (1 - M**2)**-0.417 - 1
            elif M < 1:
                LE = 1 - 1.785*(M - 0.9)
            else:
                LE = 1.214 - 0.502 / M**2 + 0.1095 / M**4
            return LE * self.fin_pressure_mult # later check ref areas again

        # base drag coefficient
        # eq 3.94, at some point need to think about whether this only applies to coasting phase
        def C_D_base(M):
            if M < 1:
                return 0.12 + 0.13*M**2
            else:
                return 0.25/M

        Cfc = C_f(R, M)
        stag_p = stagnation_pressure(M)
        CD_f = C_D_frict(Cfc)
        CD_p1 = C_D_nose_press(M, stag_p)
        CD_p2 = C_D_fin_press(M)
        CD_b = C_D_base(M)
        CD_pipes, CD_lug = parasitic_drag(stag_p)
        return [CD_f, CD_p1, CD_p2, CD_b, CD_pipes, CD_lug]
    
    # drag coefficient adjusted for angle of attack
    def C_D_axial(self, alpha, CD0):
        if alpha < 0.001: return CD0
        alpha_ = np.pi - alpha if alpha > np.pi/2 else alpha
        
        if alpha_ <= self.CD_interp_angle:
            mult = self.CD_interp_below.__call__(alpha_)
        else:
            mult = self.CD_interp_above.__call__(alpha_)
            
        if alpha < np.pi/2:
            return mult * CD0
        else:
            return -mult * CD0
        
    def physics(self, alpha, Mach, v0, omega, mu, c_o_g):
        roll = omega[2]
        pitch= omega[0]
        yaw = omega[1]
        self.mu = mu
        CoPs = [self.bodynose_CoP(alpha), self.fin_CoP(Mach)]
        coeffs = [self.C_N_alpha_lift(alpha)+ self.C_N_alpha_nose(alpha), self.normal_force_coefficient(alpha, Mach)]
        CoP = sum([CoPs[i]* CN_alpha for i, CN_alpha in enumerate(coeffs)]) / sum(coeffs)
        
        # forces
        CN   = alpha *sum(coeffs)
        CD0  = sum(self.C_D_0(self.Reynold(v0), Mach))
        CDax = self.C_D_axial(alpha, CD0)
        
        # moments
        C_damp_p, C_damp_y = self.pitch_damping_coefficient(pitch, yaw, v0, c_o_g)
        Cld_times_d = self.roll_damping_coefficient(roll, v0, Mach) * self.reference_len
        return CoP, CN, CDax, Cld_times_d, C_damp_p* self.reference_len, C_damp_y* self.reference_len # remember that CoP is referenced from nosecone, not from engine
        
    # calculates all the relevant aerodynamic things at a given instant for the simulator
    def aero(self, state, rkt, air, wind):
        x, q, v, w, _        = state[1:]
        roll              = w[2]
        p_a, rho, T_a, mu = air
        
        # Check Knudsen number and consider other drag models (e.g. rarified gas dyn vs. quadratic drag), low priority
        v0                = np.linalg.norm(v - wind)
        v_hat             = normalize(v - wind)
        v_body            = frame_rotation(q, v_hat)
        
        alpha             = np.arccos(np.clip(np.dot(v_body, np.array([0,0,1])), -1, 1))
        sound_speed       = np.sqrt(self.ka * self.Ra * T_a)
        fin_flutter       = self.flutter_velocity(sound_speed, p_a) * 0.85 # 15% factor of safety
        Ma                = v0 / sound_speed # Mach
        dyn_press         = 0.5 * rho * v0**2         # Dynamic pressure [Pa]
        
        # at some point include an option for the old look up table here, med priority
        # also med priority, but a simple planform calculation for CoP could be useful for low-fidelity model
        # low priority, consider "Active Control Stabilization of High Power Rocket" for drag equations
        if rkt.descending:
            CoP, CN, CDax, Cld_times_d, C_damp_p, C_damp_y = 0, 0, 0, 0, 0, 0
        else:
            CoP, CN, CDax, Cld_times_d, C_damp_p, C_damp_y = self.physics(alpha, Ma, v0, w, mu, rkt.CoM[2])
        
        CoP               = np.array([0,0, rkt.length - CoP]) # CoP calculated nose ref, but rest of sim is base ref
        
        direction         = -normalize(project_normal.dot(v_body)) # normalize might be unnecessary
        norm_force        = CN * direction
        
        if not rkt.off_tower:
            norm_force        *= 0
            Cld_times_d        = 0
            C_damp_p, C_damp_y = 0, 0
            
        pitch_moment      = np.cross(CoP - rkt.CoM, norm_force)
        pitch_damp        = np.array([C_damp_p * np.sign(w[0]), C_damp_y * np.sign(w[1]), 0])
        mult              = dyn_press * rkt.frontal_area
        
        torque_body       = mult * (np.array([0, 0, -Cld_times_d]) + pitch_moment + pitch_damp) 
        
        if rkt.descending:
            chute = C_D_MAIN * 12 if rkt.main_chute else C_D_DROGUE * 3
            drag  = self.tumbling_drag * rkt.frontal_area if rkt.free_fall else chute
            force_body    = -dyn_press * drag * v_hat
        else:
            force_body    = mult * (np.array([0, 0, -CDax]) + norm_force)
            
        return np.array([force_body, torque_body, v0, dyn_press, Ma, alpha, CoP[2], fin_flutter, CN, CDax], dtype=object)
