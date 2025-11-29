


# hey, welcome to one of the silliest parts of the script

from system_definition import *
from efs_system_requirements import power_req


###
# Geometry
###

import copy

class Geometry:
    """ This class is for determining generic geometric properties of 
    uniformly dense and symmetrical objects, namely, how much space the material occupies and how mass is distributed
    """
    def __init__(self, shape):
        self.shape = shape
        
    def volume(self, axial, radial, thickness): 
        """Calculates the volume occupied, not the volume contained"""
        if self.shape == 'Point':
            return 0
        elif self.shape == 'Shell':
            return axial * math.pi * ((radial + thickness)**2 - radial**2)
        elif self.shape == 'Blob':
            return (4 * math.pi / 3) * (axial/2) * radial**2
        elif self.shape == 'Fin':
            return thickness[0] * thickness[1] * (axial + radial) / 2 # these names are not descriptive here
        elif self.shape == 'Cone':
            return (1/3) * math.pi * (axial * radial**2 - (axial - thickness) * (radial - thickness)**2)
        
    def moment(self, mass, params):
        """Calculates the moment of inertia
        
        Each moment of inertia is centered on the component's center of 
        mass in the component's frame of reference.
        The parallel axis theorem should be used to 
        find moments in superstructures
        """
        if self.shape == 'Point': # input position relative to some frame
            return mass * np.diag([np.linalg.norm(np.cross(np.array([1,0,0]), params))**2,
                                   np.linalg.norm(np.cross(np.array([0,1,0]), params))**2,
                                   np.linalg.norm(np.cross(np.array([0,0,1]), params))**2])
        elif self.shape == 'Shell':
            return (mass / 12) * np.diag([3 * (params[1]**2 + params[2]**2) + params[0]**2,
                                          3 * (params[1]**2 + params[2]**2) + params[0]**2,
                                          6 * (params[1]**2 + params[2]**2)])
        elif self.shape == 'Blob':
            return (mass / 5) * np.diag([(params[(i + 1) % 3]**2 + params[(i + 2) % 3]**2)
                                         for i in range(3)])
        elif self.shape == 'Fin':
            return (mass/12) * np.diag([(params[(i + 1) % 3]**2 + params[(i + 2) % 3]**2)
                                         for i in range(3)])
        elif self.shape == 'Cone':
            return (3*mass/5) * np.diag([params[0]**2 + params[1]**2 / 4,
                                         params[0]**2 + params[1]**2 / 4,
                                         params[1]**2 / 2])


class Component:
    """This class is for organizing the properties of individual components of the rocket"""
    def __init__(self, name, material, length, radius, thickness, shape, mass=0):
        self.name      = name
        self.material  = material
        self.length    = length
        self.radius    = radius
        self.thickness = thickness
        self.geo       = Geometry(shape)
        self.volume    = self.geo.volume(length, radius, thickness)
        
        # Calculate mass if it is not prescribed
        if mass == 0:
            self.mass = self.volume * material['rho'] if shape != 'Point' else material # for point-masses
        else:
            self.mass = mass
            
        # If mass is prescribed, then moment should also be prescribed
        # Prescribed moment has not been implemented yet
        if shape == 'Point':
            self.moment = self.geo.moment(self.mass, np.array([0,0, length/2]))
        elif shape == 'Shell':
            self.moment = self.geo.moment(self.mass, [radius, radius+thickness, length])
        elif shape == 'Blob':
            self.moment = self.geo.moment(self.mass, [radius, radius, length])
        elif shape == 'Cone':
            self.moment = self.geo.moment(self.mass, [length, radius])
            
    
    def center_of_mass(self, coords_rear):
        """Sets center of mass based on where they are created
        
        coords_rear = coordinates of component from base of rocket
        """
        if self.geo.shape == 'Cone':
            self.CoM = coords_rear + np.array([0,0, self.length/4])
        else:
            self.CoM = coords_rear + np.array([0, 0, self.length/2])

# fins require special calculations to determine their center of mass relative to the rocket
# because they are not located down the axis of symmetry
# at some point, consolidate most of the fin-related calculations from the aero model into this class, low priority
class Fin(Component):
    def __init__(self, name, material, root, tip, sweep_angle, semispan, thickness):
        self.name        = name
        self.material    = material
        self.root        = root
        self.tip         = tip
        self.sweep_angle = sweep_angle
        self.semispan    = semispan
        self.thickness   = thickness
        self.geo         = Geometry('Fin')
        self.volume      = self.geo.volume(root, tip, [semispan, thickness])
        self.mass        = self.volume * material['rho'] + FIN_BRACKET
        if name == 'Front' or name == 'Back':
            self.moment  = self.geo.moment(self.mass, [thickness, semispan, root])
        else:
            self.moment  = self.geo.moment(self.mass, [semispan, thickness, root])
        
    def center_of_mass(self, coords_rear):
        leg_b             = self.semispan / np.cos(self.sweep_angle)
        c                 = np.sqrt(leg_b**2 - self.semispan**2)
        d                 = self.tip + c - self.root
        leg_a             = np.sqrt(d**2 + self.semispan**2)
        self.sweep_length = c
        x                 = 0
        y                 = (self.semispan/3) * (self.root + 2*self.tip)/(self.root + self.tip)
        z                 = self.root/2 + (2*self.tip + self.root)*(leg_a**2 - leg_b**2)/(3*(self.root + self.tip))
        if self.name == 'Front':
            self.CoM      = coords_rear + np.array([x, y, z])
        elif self.name == 'Back':
            self.CoM      = coords_rear + np.array([x, -y, z])
        elif self.name == 'Left':
            self.CoM      = coords_rear + np.array([y, x, z])
        elif self.name == 'Right':
            self.CoM      = coords_rear + np.array([-y, x, z])

# engines and rcs nozzles don't have special geometry but they have many dynamic properties in the simulation
# this class exists because there isn't much difference between hot gas and cold gas thrusters
# and it will make simulation code more readable when we call thrust function
class Engine(Component):
    def __init__(self, mdot, p_e, p_ch, T_ch, ke, mm, throttle_window, min_throttle, is_RCS):
        if is_RCS:
            Component.__init__(self, 'Gas Jet', ALUM, 0.14, 0.07, 0, 'Blob')
        else:
            Component.__init__(self, 'Engine', ALUM, 0.3, .1, .011, 'Shell')
            
        self.thrust_vector        = np.array([0, 0, 1]) # as in, points straight up to the nose
        self.mdot                 = mdot
        self.p_e                  = p_e
        self.empty                = False
        self.p_ch                 = p_ch
        self.T_ch                 = T_ch
        self.ke                   = ke
        self.mm                   = mm
        self.Re                   = R_UNIV / mm
        self.throttle             = [1.]
        self.min_throttle         = min_throttle
        self.can_throttle         = min_throttle != 1.
        self.throttle_window      = throttle_window
        self.throttle_rate        = (min_throttle - 1.) / (throttle_window[1] - throttle_window[0])
        self.throttle_y_intercept = 1. - self.throttle_rate * throttle_window[0]
        
        # Throat pressure      [Pa]
        self.p_t                  = p_ch * (1 + (ke - 1)/2)**(-ke /(ke - 1))
        # Throat temperature   [K]
        self.T_t                  = T_ch*(1 / (1 + (ke - 1)/2))
        # Throat area          [m^2]
        self.A_t                  = (mdot / self.p_t) * np.sqrt(self.Re * self.T_t / ke)
        # Expansion ratio
        self.ex                   = (2 / (ke + 1))**(1 / (ke - 1)) * (p_ch / p_e)**(1/ke) / np.sqrt(
                                        (ke + 1) / (ke - 1) * (1 - (p_e / p_ch)**(((ke - 1) / ke))))
        # Exit area [m^2] 
        self.A_e                  = self.ex * self.A_t
        
        L                         = 0.1 # m from throat to exit
        # divergence angle
        alpha                     = np.arctan((np.sqrt(self.A_e) - np.sqrt(self.A_t))/(L * np.sqrt(math.pi)))
        # Thrust divergence losses correction,
        self.lam                  = 0.5 * (1 + np.cos(alpha))
        # exit velocity with 3% performance loss assumed from film cooling and whatever else
        self.Ve                   = 0.97 * self.lam * np.sqrt(2*self.ke / (self.ke - 1) * self.Re * self.T_ch *
                             (1 - (self.p_e/p_ch)**((self.ke - 1)/self.ke)))
    
    # linear throttle for simplicity
    def throttle_engine(self, drag):
        if self.can_throttle == False:
            return 1.
        elif drag <= self.throttle_window[0]:
            return 1.
        elif drag >= self.throttle_window[1]:
            return self.min_throttle
        else:
            return drag * self.throttle_rate + self.throttle_y_intercept
    
    # calculates engine thrust at a height
    def thrust(self, p_a, thrtl):
        return thrtl*self.mdot*self.Ve + (self.p_e - p_a)*self.A_e

# An RCS system is basically a bunch of little engines. We're not really using this yet.
class RCS(Engine):
    def __init__(self, mdot, p_e, p_ch, T_ch, ke, mm, throttle_window, min_throttle, is_RCS):
        self.up       = np.array([0,0,1]) # for convenience
        self.bang_on  = 0.00872665 # radians
        self.bang_off = 0.00872665 # for now, no hysteresis
        self.active   = False
        super().__init__(mdot, p_e, p_ch, T_ch, ke, mm, throttle_window, min_throttle, is_RCS)
        
    def controller(self, q, p_a, rkt_CoM):
        num_thrusters = 0
        up_body       = frame_rotation(q, self.up)
        deviation     = np.arccos(np.clip(np.dot(up_body, self.up), -1, 1))
        # kinda naive switching control
        if (self.tank.has_spare_gas(p_ch) and
                ((not self.active and deviation >= self.bang_on) or
                (self.active and deviation >= self.bang_off))):
            self.active   = True
            sign_vec      = -np.sign((np.cross(self.up, np.cross(self.up, up_body))))
            num_thrusters = sum([abs(e) for e in sign_vec])
            rcs_thrust    = self.thrust(p_a, 1.0)
            force         = rcs_thrust * sign_vec
            actuators     = self.out_rad*sign_vec + self.CoM - rkt_CoM
            torque        = np.cross(actuators, force)
        else:
            self.active   = False
            num_thrusters = 0
            force         = np.zeros(3)
            torque        = np.zeros(3)
        return (num_thrusters, force, torque)





###
# Structures
###

class Structure:
    def __init__(self, name, static=False):
        self.name   = name
        self.static = static
        self.parts  = []
        
    
    def add_part(self, name, material, length, radius, thickness, shape, 
                 height_coord, x_coord=0, y_coord=0, mass=0, prepend=False):
        coords = np.array([x_coord, y_coord, height_coord])
        if prepend:
            self.parts.insert(0, Component(name, material, length, radius, thickness, shape, mass))
            self.parts[0].center_of_mass(coords)
        else:
            self.parts.append(Component(name, material, length, radius, thickness, shape, mass))
            self.parts[-1].center_of_mass(coords)

            
    def add_component(self, component, height_coord, x_coord=0, y_coord=0, prepend=False):
        """Adds a component object to the structure, assuming it is on the axis"""
        coords = np.array([x_coord, y_coord, height_coord])
        if prepend:
            self.parts.insert(0, component)
            self.parts[0].center_of_mass(coords)
        else:
            self.parts.append(component)
            self.parts[-1].center_of_mass(coords)
            
        
    def add_engine(self, mdot, p_e, p_ch, T_ch, ke, Re, throttle_window, min_throttle, 
                   height_coord, x_coord=0, y_coord=0, is_RCS=False):
        coords = np.array([x_coord, y_coord, height_coord])
        if is_RCS:
            self.parts.append(RCS(mdot, p_e, p_ch, T_ch, ke, Re, throttle_window, min_throttle, is_RCS))
        else:
            self.parts.append(Engine(mdot, p_e, p_ch, T_ch, ke, Re, throttle_window, min_throttle, is_RCS))
        self.parts[-1].center_of_mass(coords)
    
    def add_structure(self, struct):
        self.parts.append(struct)
    
    def sum_parts(self):
        def parallel_axis(CoM_S, CoM_P, I, m):
            delta = CoM_P - CoM_S
            return I + m * (np.dot(delta, delta)*np.eye(3) - np.outer(delta, delta))
        
        for i, part in enumerate(self.parts):
            if hasattr(part, 'parts') and (not self.static or not hasattr(part, 'moment')):
                part.sum_parts()
                
        self.mass   = sum([p.mass for p in self.parts])
        self.CoM    = np.sum([p.CoM * p.mass for p in self.parts],
                          axis=0) / self.mass
        self.moment = np.sum([parallel_axis(self.CoM,
                                            p.CoM,
                                            p.moment,
                                            p.mass) for p in self.parts], axis=0)
        
    def read_out(self, description=[], indents=0):
        if indents==0:
            self.description = []
            description      = self.description
            self.description.append('\n' + self.name + "\tGLOW:" + str(self.GLOW) + "\tCurrent Mass:" + str(self.mass))
        for i, part in enumerate(self.parts):
            part_description = str('\n' + "  " * indents + str(i) + ".  " + part.name + 
                               "\t" + ("   " * indents) + "Mass: " + str(round(part.mass, 2)))
            description.append(part_description.expandtabs(32))
            if hasattr(part, 'parts'):
                part.read_out(description, indents + 1)

class Module(Structure):
    def __init__(self, name, material, length, radius, 
                 height_coord, x_coord=0, y_coord=0, 
                 shape='Shell', mono_thickness=0, wt_thickness=0):
        self.name   = name
        self.parts  = []
        self.length = length
        self.static = True
        self.material     = material
        self.p_0          = 101325 # 1 atm pressure
        
        coords = np.array([x_coord, y_coord, height_coord])
        self.coords = coords
        
        # airframe layup layer: material, distance from ID, thickness
        if material is ALUM: # for micromodules
            layers  = [[material, 0, wt_thickness*1.2]] # fudge factor
        elif shape == 'Shell': # for composite modules
            layers  = [[material, 0, INNER_CF_THK],
                      [NOMEX, INNER_CF_THK, NOMEX_THK],
                      [material, INNER_CF_THK + NOMEX_THK, OUTER_CF_THK]]
        else: # for conical part of nosecone
            layers  = [[material, 0, 0.00127]]
            
        self.thickness   = sum([layer[2] for layer in layers]) if material is not ALUM else mono_thickness
        
        for i in range(len(layers)):
            self.add_part(layers[i][0]['name'] + ' Layer ' + str(i+1), 
                          layers[i][0], length, radius + layers[i][1], layers[i][2], shape,
                         height_coord, x_coord, y_coord)

                
        if shape == 'Shell': 
            # so we don't add these to conical part of nosecone
            #self.add_part(coords + np.array([0,0,length]), 'CouplingRings', ALUM,
            #          0.056642, 0.14505, 0.00735, 'Shell')
            bulkhead_component = Component('Bulkhead', ALUM, 0.00406375, 0.2, 0.05, 'Shell', mass=0.5)
            self.add_part('Bulkhead', ALUM, 0.00406375, 0.2, 0.05, 'Shell', 
                          height_coord, x_coord, y_coord, mass=0.5)
            self.add_part('Bulkhead', ALUM, 0.00406375, 0.2, 0.05, 'Shell', 
                          height_coord+length, x_coord, y_coord, mass=0.5)

class Nosecone(Structure):
    def __init__(self, material, cyl_length, cone_length, radius, thickness, 
                 height_coord, x_coord=0, y_coord=0):
        self.name   = 'Nosecone'
        self.material     = material
        self.parts  = []
        self.length = cyl_length + cone_length
        self.thickness = thickness
        self.static = True
        self.p_0          = 101325 # 1 atm pressure
        
        coords = np.array([x_coord, y_coord, height_coord])
        self.coords = coords
        self.add_structure(Module('Cylinder', material, cyl_length, radius, height_coord, x_coord, y_coord))
        self.parts[-1].parts.pop() # throw away the upper bulkhead
        self.add_structure(Module('Cone', material, cone_length, radius, 
                                  height_coord+cyl_length, x_coord, y_coord, shape='Cone'))
    
class Tank(Structure):
    def __init__(self, material, in_radius, mono_thickness, wt_thickness, 
                 prop_mass, prop_matrl, p_0, height_coord, x_coord=0, y_coord=0, 
                 insulation=True, ullage=1.1, compressible=False):
        self.name         = prop_matrl['name']+' Tank'
        self.prop_matrl   = prop_matrl
        self.material     = material
        self.in_radius    = in_radius
        self.parts        = []
        self.static       = False
        fl_ht             = self.fluid_height(prop_mass)
        self.compressible = compressible
        length            = ullage * fl_ht
        self.length       = length
        self.volume       = length * np.pi * in_radius**2
        self.p_0          = p_0
        self.thickness    = mono_thickness
        
        coords = np.array([x_coord, y_coord, height_coord])
        self.coords       = coords
        self.add_part('Tank Walls', material, length, in_radius, wt_thickness, 'Shell', 
                     height_coord, x_coord, y_coord)
        if insulation:
            self.add_part('Insulation', CRYOGEL, length,
                          in_radius+SKIN_T, INSULTN_THKNS, 'Shell',
                         height_coord, x_coord, y_coord)
        self.add_part('Bottom Cap', material, wt_thickness*2, 0, in_radius, 'Shell', 
                     height_coord, x_coord, y_coord)
        self.add_part('Top Cap', material, wt_thickness*2, 0, in_radius, 'Shell',
                     height_coord+length, x_coord, y_coord)
        if not compressible: # so we don't add these to N2 tank
            #self.add_part(coords + np.array([0,0,length]), '2Rings&6Clamps', ALUM,
            #          0.056642, 0.14505, 0.00735, 'Shell')
            bulkhead_component = Component('Bulkhead', ALUM, 0.00406375, 0.2, 0.05, 'Shell', mass=0.5)
            self.add_part('Bulkhead', ALUM, 0.00406375, 0.2, 0.05, 'Shell', 
                          height_coord, x_coord, y_coord, mass=0.5)
            self.add_part('Bulkhead', ALUM, 0.00406375, 0.2, 0.05, 'Shell', 
                          height_coord+length, x_coord, y_coord, mass=0.5)
        self.dry_m = sum([part.mass for part in self.parts])
        self.add_part('Propellant', prop_matrl, fl_ht, 0, in_radius, 'Shell',
                     height_coord, x_coord, y_coord)

    def fluid_height(self, prop_mass):
        return prop_mass / (self.prop_matrl['rho'] * math.pi * self.in_radius**2)
        
    def drain(self, delta):
        prop     = self.parts.pop()
        new_prop = prop.mass - abs(delta)
        if self.compressible:
            self.prop_matrl['rho'] = new_prop / self.volume
            new_ht                 = self.length
        else:
            new_ht                 = self.fluid_height(new_prop)
        self.add_part('Propellant', self.prop_matrl, new_ht, 0, self.in_radius, 'Shell',
                     self.coords[2], self.coords[0], self.coords[1])

# jank wip
class N2_tank(Tank):
    def __init__(self, mdot, m0, T_tank, mm, tank_l, tank_r, height_coord, x_coord=0, y_coord=0):
        self.mdot     = mdot
        self.T_tank   = T_tank
        self.mm       = mm
        self.tank_V   = tank_l * tank_r**2 * np.pi
        self.tank_r   = tank_r
        self.tank_l   = tank_l
        self.gas_mass = m0
        density       = self.gas_mass / self.tank_V
        p_0           = density * T_tank * R_UNIV / mm
        
        super().__init__(CFIBER, tank_r, 0.003175, 0.003175, self.gas_mass,
                         Material('N2', density), p_0, height_coord, x_coord, y_coord,
                         insulation=False, ullage=1, compressible=True)
        
    # at some point, include function for keeping fuel tanks pressurized (at what pressure??)
    def ideal_mass(self, P, V, M, T):
        return P * V * M / (T * R_UNIV)
    def pressure(self):
        return self.parts[-1].mass * self.T_tank * R_UNIV / (self.tank_V * self.mm)
    # at some point, introduce check to make sure we don't need rest of gas for fuel tanks
    def has_spare_gas(self, p_ch):
        return self.parts[-1].mass > 2*self.mdot and self.pressure() > p_ch
        
class Rocket(Structure):
    def has_fuel(self):
        return self.lox_tank.parts[-1].mass > self.residual_o and self.ipa_tank.parts[-1].mass > self.residual_f




####
# Create a rocket from a function instead of from a config file
####

def create_rocket(mprop, mdot, p_e,
                  p_ch, T_ch, ke, mm,
                  throttle_window, min_throttle,
                  airfrm_in_rad, ipa_wt, of,
                  rcs_mdot, rcs_p_e, rcs_p_ch,
                  ballast, root, tip, sweep, span, thickness, con_nose_l,
                  tank_p_o, tank_p_f, rib_t, num_radl_dvsns):
    rocket    = Rocket('LV4')
    rocket.unadjusted_mprop = mprop
    fuel = Material('IPA/H20', ipa_wt/100* IPA['rho'] + (100 - ipa_wt)/100 * H20['rho'],
                             mu=0.00192, p_v=8840) # kg/m^3  Density of 64.8% IPA / 35.2% H20
    mdot_o, mdot_f = proportion(mdot, of)
    m_o, m_f  = proportion(mprop, of)
    coolant = m_f * 0.06
    burntime = mprop / mdot
    rocket.mdot_c = coolant / burntime
    m_f += coolant
    rocket.residual_o = 0.02 * m_o
    rocket.residual_f = 0.02 * m_f
    out_rad   = airfrm_in_rad + AIRFRAME_THICKNESS
    height    = 0
    
    rocket.rib_t = rib_t
    rocket.num_radl_dvsns = num_radl_dvsns
    cell_height     = TANK_OD * np.pi / num_radl_dvsns
    alpha           = (rib_t * RIB_DEPTH) / (SKIN_T * cell_height) # Web non-dimensional ratio, pg. 2.0.008
    beta            = (3 * alpha * (1 + DELTA)**2 + (1 + alpha) * (1 + alpha * DELTA**2))**0.5
    t_star          = SKIN_T * beta / (1 + alpha) # equivilent monocoque thickness, Eq. 2.5.3
    rocket.TANK_EQV_WT_T   = SKIN_T * (1 + 3 * alpha) # tank equivalent weight thickness, m
    rocket.TANK_MAX_P      = 0.666 * ALUM['Su'] * SKIN_T * (1 + alpha) / TANK_IN_RAD # Pa, prop tank max pressure
    

    def cpl_ring(size, name=None):
        """Generates a coupling ring. Pass object to create_rocket in the future
        size = the factor to multiply the height by
        see HALF_CPL_RING and THREE_QTR_CPL_RING in Display_Information
        """
        name = 'Coupling Rings {:.2} ht'.format(size) if name is None else name
        return Component(name, ALUM, size*COUPLING_RING_HT, AIRFRM_IN_RAD, CPLNG_RING_THK, 'Shell')
    
    def gen_fin(name):
        """Generates a fin with a given name. Pass object to create_rocket in the future"""
        return Fin(name, ALUM, root, tip, sweep, span, thickness)
        
    
    # Add parts and structures
    engine_subsystem = Structure('Engine Subsystem')
    engine_subsystem.add_structure(Module('Fin Can Module', CFIBER, FIN_CAN_L, airfrm_in_rad, 
                                          height+THRST_PLT+HALF_CPL_RING))
    engine_subsystem.parts[-1].add_part('Thrust Plate', ALUM, THRST_PLT, 0.0254, airfrm_in_rad - 0.0254, 'Shell',
                                       height)
    engine_subsystem.parts[-1].add_component(cpl_ring(0.5), height+THRST_PLT)
    engine_subsystem.parts[-1].add_part('Launch Button', 0.2051, 0.02, 0, 0, 'Point',
                                       height+THRST_PLT, x_coord=-out_rad)

    fin_set = Structure('Fins')
    fin_set.add_component(gen_fin('Front'), height+FIN_ROOT_HEIGHT, y_coord=out_rad)
    fin_set.add_component(gen_fin('Back'), height+FIN_ROOT_HEIGHT, y_coord=-out_rad)
    fin_set.add_component(gen_fin('Left'), height+FIN_ROOT_HEIGHT, x_coord=out_rad)
    fin_set.add_component(gen_fin('Right'), height+FIN_ROOT_HEIGHT, x_coord=-out_rad)
    rocket.fin = fin_set.parts[0]
    engine_subsystem.parts[-1].add_structure(fin_set)
    
    engine_subsystem.parts[-1].add_engine(mdot, p_e, p_ch, T_ch, ke, mm, throttle_window, min_throttle, height - ENG_CLEARANCE)
    rocket.engine    = engine_subsystem.parts[-1].parts[-1] # for convenience
    if USE_EFS:
        engine_subsystem.parts[-1].add_part('EMS', ALUM, 0.1016, 0.041, 0, 'Blob',
                                           L_ENGINE)
        engine_subsystem.parts[-1].add_part('EFS Plumbing', PLUMBING_M, L_FEED, 0, 0, 'Point',
                                           L_ENGINE+L_EMS)
    efs_height       = L_ENGINE + L_EMS + 0.5 * L_FEED
    fin_can_module   = engine_subsystem.parts[-1]
    
    rocket.add_structure(engine_subsystem.parts[-1])
    
    
    height += THRST_PLT + HALF_CPL_RING + FIN_CAN_L
    
    #rocket.add_structure(Module(np.array([0, 0, height+ THREE_QTR_CPL_RING]), 'Engine Passthru Module', ALUM, PAS_L, airfrm_in_rad, 'Shell', t_star, rocket.TANK_EQV_WT_T))
    rocket.add_structure(Module('Engine Passthru Module', ALUM, PAS_L, airfrm_in_rad, height+THREE_QTR_CPL_RING,
                                mono_thickness=AIRFRAME_THICKNESS, wt_thickness=AIRFRAME_THICKNESS))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    height += THREE_QTR_CPL_RING
    rocket.parts[-1].add_part('Launch Button', 0.2051, 0.02, 0, 0, 'Point',
                             height, x_coord=-out_rad)
    rocket.parts[-1].add_part('N2 to Eng Pipe', ALUM, N2_TO_ENG_L, N2_PIPE_IR, PIPE_THK, 'Shell', 
                              height, x_coord=-out_rad) # 1 in OD, 1/8th in thick
    rocket.parts[-1].add_part('IPA to Eng Pipe', ALUM, IPA_TO_ENG_L, FUEL_PIPE_IR, PIPE_THK, 'Shell', 
                              height, x_coord=out_rad) # 1 in OD, 1/8th in thick
    lug_1 = height + 0.5 * PAS_L 
    height += PAS_L
    
    engine_subsystem.add_structure(Tank(ALUM, TANK_IN_RAD, t_star, rocket.TANK_EQV_WT_T, m_o, LOX, tank_p_o, height+HALF_CPL_RING))
    engine_subsystem.parts[-1].add_component(cpl_ring(0.5), height, prepend=True)
    height += HALF_CPL_RING
    rocket.m_tank_o = engine_subsystem.parts[-1].dry_m
    rocket.l_o      = engine_subsystem.parts[-1].length
    rocket.lox_tank = engine_subsystem.parts[-1] # for convenience
    rocket.add_structure(engine_subsystem.parts[-1])
    surface_height_lox = height + rocket.lox_tank.parts[-1].length
    height += rocket.l_o
    
    #rocket.add_structure(Module(np.array([0, 0, height+HALF_CPL_RING]), 'Tank Passthru Module', ALUM, PAS_L, airfrm_in_rad, 'Shell', t_star, rocket.TANK_EQV_WT_T))
    rocket.add_structure(Module('1st Tank Passthru Module', ALUM, PAS_L, airfrm_in_rad, height+HALF_CPL_RING, 
                                mono_thickness=AIRFRAME_THICKNESS, wt_thickness=AIRFRAME_THICKNESS))
    rocket.parts[-1].add_component(cpl_ring(0.5), height)
    height += HALF_CPL_RING
    rocket.parts[-1].add_part('Plumbing', BETWEEN_TANKS_M, 0.02, 0, 0, 'Point',
                             height)
    rocket.parts[-1].add_part('Launch Button', 0.2051, 0.02, 0, 0, 'Point',
                             height, x_coord=-out_rad)
    rocket.parts[-1].add_part('N2 to LOX Pipe', ALUM, N2_TO_LOX_L, N2_PIPE_IR, PIPE_THK, 'Shell',
                             height, y_coord=out_rad) # 1 in OD, 1/8th in thick
    lug_2 = height + 0.5 * PAS_L
    height += PAS_L
    
    rocket.add_structure(Module('1st Propulsion MicroModule', CFIBER, MICRO_L, airfrm_in_rad,
                               height+THREE_QTR_CPL_RING))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    # three quarter, prepend false 
    height += THREE_QTR_CPL_RING + MICRO_L
    
    rocket.add_structure(Module('2nd Tank Passthru Module', ALUM, PAS_L, airfrm_in_rad, height+THREE_QTR_CPL_RING, 
                                mono_thickness=AIRFRAME_THICKNESS, wt_thickness=AIRFRAME_THICKNESS))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    # three quarter, prepend false
    height += THREE_QTR_CPL_RING + PAS_L
    
    # ADD THE LOX TANK
    engine_subsystem.add_structure(Tank(ALUM, TANK_IN_RAD, t_star, 
                                        rocket.TANK_EQV_WT_T, m_f, fuel, tank_p_f, height+HALF_CPL_RING))
    engine_subsystem.parts[-1].add_component(cpl_ring(0.5), height, prepend=True)
    height += HALF_CPL_RING
    rocket.m_tank_f = engine_subsystem.parts[-1].dry_m
    rocket.l_f      = engine_subsystem.parts[-1].length
    rocket.ipa_tank = engine_subsystem.parts[-1] # for convenience
    rocket.add_structure(engine_subsystem.parts[-1])
    surface_height_ipa = height + rocket.ipa_tank.parts[-1].length
    height += rocket.l_f
    
    # ADD THE EFS TO THE TANK
    rocket.delp_regen = REGEN_MULT / fuel['rho'] * (mdot_f/(REGEN_N * np.pi * 0.25 * REGEN_D**2))**2
    rocket.v_lfets_o, rocket.p_out_o, rocket.pow_o, rocket.rpm_o = power_req(p_ch, mdot_o, LOX, tank_p_o, FRIC_O, PLUMBING_L_O, DELP_INJ_O, surface_height_lox-efs_height)
    rocket.v_lfets_f, rocket.p_out_f, rocket.pow_f, rocket.rpm_f = power_req(p_ch, mdot_f, fuel, tank_p_f, FRIC_F, PLUMBING_L_F, DELP_INJ_F + rocket.delp_regen, surface_height_ipa-efs_height)
    total_pow = rocket.pow_o + rocket.pow_f
    efs_motors_m = total_pow / MOT_SPEC_POW 
    efs_bats_m   = 2 * total_pow / BAT_SPEC_POW
    if USE_EFS:
        fin_can_module.add_part('EFS 2 Motors/ESCs', efs_motors_m, L_FEED, 0, 0, 'Point', L_ENGINE+L_EMS)
        fin_can_module.add_part('EFS 4 Batteries', efs_bats_m, L_FEED, 0, 0, 'Point', L_ENGINE+L_EMS)
    
    
    rocket.eng_sys = engine_subsystem
    rocket.eng_sys.sum_parts()
    
    #rocket.add_structure(Module(np.array([0, 0, height+HALF_CPL_RING]), 'Tank Passthru Module', ALUM, PAS_L, airfrm_in_rad, 'Shell', t_star, rocket.TANK_EQV_WT_T))
    rocket.add_structure(Module('3rd Tank Passthru Module', ALUM, PAS_L, airfrm_in_rad, height+HALF_CPL_RING, 
                                mono_thickness=AIRFRAME_THICKNESS, wt_thickness=AIRFRAME_THICKNESS))
    rocket.parts[-1].add_component(cpl_ring(0.5), height)
    rocket.parts[-1].add_part('Plumbing', ABOVE_FUEL_M, 0.02, 0, 0, 'Point', height)
    rocket.parts[-1].add_part('N2 to IPA Pipe', ALUM, N2_TO_IPA_L, N2_PIPE_IR, PIPE_THK, 'Shell',
                             height, y_coord=-out_rad) # 1 in OD, 1/8th in thick
    height += PAS_L + HALF_CPL_RING
    
    rocket.add_structure(Module('2nd Propulsion MicroModule', CFIBER, MICRO_L, airfrm_in_rad,
                               height+THREE_QTR_CPL_RING))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    height += THREE_QTR_CPL_RING + MICRO_L
    
    rocket.add_structure(Module('4th Tank Passthru Module', ALUM, PAS_L, airfrm_in_rad, height+THREE_QTR_CPL_RING, 
                                mono_thickness=AIRFRAME_THICKNESS, wt_thickness=AIRFRAME_THICKNESS))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    height += THREE_QTR_CPL_RING + PAS_L
    
    rocket.add_structure(Module('AV Module', FIBERGLASS, AV_L, airfrm_in_rad, height+THREE_QTR_CPL_RING))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    height += THREE_QTR_CPL_RING
    rocket.parts[-1].add_part('AV/360', ALUM, 0.3, .0433, 0, 'Blob', height)
    height += AV_L
    
    #rocket.add_structure(Module(np.array([0, 0, height+ THREE_QTR_CPL_RING]), 'CAM/N2 Passthru Module', ALUM, CAM_L, airfrm_in_rad, 'Shell', t_star, rocket.TANK_EQV_WT_T))
    rocket.add_structure(Module('CAM/N2 Passthru Module', ALUM, CAM_L, airfrm_in_rad, height+THREE_QTR_CPL_RING,
                                mono_thickness=AIRFRAME_THICKNESS, wt_thickness=AIRFRAME_THICKNESS))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    height += CAM_L + THREE_QTR_CPL_RING
    rocket.parts[-1].add_part('Plumbing', UNDER_N2_M, 0.02, 0, 0, 'Point', height)
    rocket.parts[-1].add_part('Cameras', ALUM, 0.05, 0.05, 0, 'Blob', height)
    
    min_n2_mass = 1.5 * (n2_prop_reqs(rocket) + 1) # safety factor and assumed RCS needs
    n2_tank_vol = min_n2_mass * R_UNIV * N2_TEMP / (N2_MM * MAX_N2_TANK_P)
    n2_tank_l   = n2_tank_vol / (np.pi * N2_TANK_OR**2)
    if n2_tank_l < 0.45:
        n2_tank_l = 0.45
    # maximum tank length assuming standard module and need for 4 inch of room above or below
    elif n2_tank_l > N2_L - 0.1016:
        n2_tank_l = N2_L - 0.1016  
        
    rocket.add_structure(Module('N2 Module', CFIBER, N2_L, airfrm_in_rad, height+THREE_QTR_CPL_RING))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    height += THREE_QTR_CPL_RING
    rocket.parts[-1].add_structure(N2_tank(rcs_mdot, min_n2_mass,
                                          N2_TEMP, N2_MM, n2_tank_l, N2_TANK_OR, height))
    rocket.rcs_tank = rocket.parts[-1].parts[-1] # for convenience
    height += N2_L
    
    #rocket.add_structure(Module(np.array([0, 0, height+ THREE_QTR_CPL_RING]), 'RCS Module', ALUM, RCS_L, airfrm_in_rad, 'Shell', t_star, rocket.TANK_EQV_WT_T))
    rocket.add_structure(Module('RCS Module', ALUM, RCS_L, airfrm_in_rad, height+THREE_QTR_CPL_RING,
                                mono_thickness=AIRFRAME_THICKNESS, wt_thickness=AIRFRAME_THICKNESS))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    height += THREE_QTR_CPL_RING
    rocket.parts[-1].add_engine(rcs_mdot, rcs_p_e, rcs_p_ch, N2_TEMP, N2_KE, N2_MM, (200,100), 1, height, is_RCS=True)
    rocket.parts[-1].parts[-1].out_rad = out_rad
    rocket.parts[-1].parts[-1].tank = rocket.rcs_tank
    rocket.rcs_sys = rocket.parts[-1].parts[-1] # for convenience
    height += RCS_L
    
    #rocket.add_structure(Module(np.array([0, 0, height+HALF_CPL_RING]), 'ERS Module', ALUM, ERS_L, airfrm_in_rad, 'Shell', t_star, rocket.TANK_EQV_WT_T))
    rocket.add_structure(Module('ERS Module', ALUM, ERS_L, airfrm_in_rad, height+HALF_CPL_RING,
                                mono_thickness=AIRFRAME_THICKNESS, wt_thickness=AIRFRAME_THICKNESS))
    rocket.parts[-1].add_component(cpl_ring(0.5), height)
    height += HALF_CPL_RING
    rocket.parts[-1].add_part('ERS', ALUM, 0.1225, 0.07, 0, 'Blob', height)
    height += ERS_L
    
    # ADD NOSECONE
    rocket.add_structure(Nosecone(CFIBER, CYL_NOSE_L,  con_nose_l, airfrm_in_rad, 0.00762, height+THREE_QTR_CPL_RING))
    rocket.parts[-1].add_component(cpl_ring(0.75), height)
    height += THREE_QTR_CPL_RING
    rocket.parts[-1].add_part('Parachutes', ALUM, 0.0575, 0.057, 0, 'Blob', height+0.3)
    rocket.nose_l = con_nose_l + CYL_NOSE_L
    height += rocket.nose_l
    rocket.parts[-1].add_part('Nosetip', NOSETIP, 0.1, 0, 0, 'Point', height-0.1)
    rocket.parts[-1].add_part('Ballast', max(ballast, 0), 0.1, 0, 0, 'Point', height-0.2)
    
    # could calculate this based on tank length
    rocket.lug_separation   = lug_2 - lug_1 # m, between fore and aft launch lugs.
    rocket.tip_off_error    = False
    rocket.drag_perturb     = 1
    rocket.sum_parts()
    rocket.GLOW             = rocket.mass
    rocket.m_o, rocket.m_f  = m_o, m_f
    rocket.OF               = of
    rocket.inr_r            = airfrm_in_rad
    rocket.diameter         = TANK_OD
    rocket.out_rad          = out_rad
    rocket.length           = height
    rocket.frontal_area     = math.pi * out_rad**2
    rocket.eng_sys_dry_mass = rocket.m_tank_o + rocket.m_tank_f + engine_subsystem.parts[0].mass
    rocket.eng_sys_len      = FIN_CAN_L + PAS_L + rocket.l_o + PAS_L + rocket.l_f
    rocket.ballast          = max(ballast, 0)
    return rocket


