

# original credit to Julio, edits by Cory, clean up by Max


"""
Global varialbes used

A_PIPE
LFETS_PIPE_AREA
D_PIPE
G_N (gravity)
PUMP_EFF (pump efficiency)
"""

from system_definition import *

RtoD = EPSILON_PIPE / D_PIPE # Rougness to Diameter Ratio



def pipe_flow(mdot, material):
    """

    args:
        mdot: mass flow rate (kg/s)
        material: material object (fluid)

    returns:
        qdot: volumetric flow rate (m^3/s)
        v: bulk fluid velocity
        v_lfets: bulk velocity of LFETS pipe
        Re: reynolds number
    """
    qdot = mdot / material['rho'] # m^3/s Volumetric Flowrate
    v    = qdot / A_PIPE # m/s Fluid Velocity
    v_lfets = qdot / LFETS_PIPE_AREA
    Re   = material['rho'] * v * D_PIPE / material['mu'] # Reynolds Number
    return qdot, v, v_lfets, Re




def head_loss(fric, plum_len, v):
    """ calculates total head loss (major + minor) """
    h_LMajor = fric * plum_len * v**2 / (D_PIPE * 2 * G_N) # m Head Major Loss
    h_LMinor = K_L * v**2 / (2 * G_N) # m Head Minor Loss
    return h_LMajor + h_LMinor



"""
Find the power requirements using the bernouli energy equation

Assumptions

v_out = v_in: no leaks in system
z_out = z_in: changes in height are not significant
gamma = const: the specific weight of the fluid does not change

Final form

h = ((p_out - p_in)/gamma) + h_L
"""

def power_req(p_ch, mdot, material, p_in, fric, plum_len, delP_inj, static_head):
    """Pump power requirements
    """
    qdot, v, v_lfets, Re = pipe_flow(mdot, material)
    h_L         = head_loss(fric, plum_len, v) # Total Head Loss
    p_out       = 1.25 * (p_ch + delP_inj) # Required Outlet Pressure including injector and venturi losses
    delP        = p_out - p_in
    spec_dens   = G_N * material['rho']
    h_s         = delP / material['rho'] + h_L * G_N # Shaft Work Head times standard gravity
    if h_s <= 0:
        return v_lfets, p_out, 0, 0
    W           = mdot * h_s / PUMP_EFF # power requirement
    tau         = 2 # no less than 2, to avoid cavitation
    inlet_head  = p_in / spec_dens
    vap_p_head  = material['p_v'] / spec_dens
    npsh_a      = inlet_head + static_head - vap_p_head - h_L # available net pressure suction head
    npsh_r      = max(0, npsh_a / tau) # rotational. pardon the kludge...
    rpm         = U_SS * (npsh_r*3.281)**0.75 / (21.2 * np.sqrt(qdot*35.31)) # magic numbers are unit conversion weirdness
    #print(material['name'], ' SH: ',static_head, ' HL: ', h_L)
    return v_lfets, p_out, W * 1.111, rpm * 1.111 # FOS so that operating conditions are 80% of maximum



