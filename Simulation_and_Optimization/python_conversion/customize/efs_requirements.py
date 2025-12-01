

# original credit to Julio, edits by Cory, clean up by Max


import numpy as np

# if you really wanted, you could replace all 'accel' arguments with standard gravity constant



def pipe_flow(material, mdot, d_pipe):
    """Returns additional flow properties

    args:
        mdot: mass flow rate (kg/s)
        d_pipe: inner diameter of pipe (m)
        material: material object (fluid)

    returns:
        qdot: volumetric flow rate (m^3/s)
        v: bulk fluid velocity
        Re: reynolds number
    """
    qdot = mdot / material['rho'] # m^3/s Volumetric Flowrate
    a_pipe = (np.pi/4)*(d_pipe**2)
    v    = qdot / (a_pipe)  # m/s Fluid Velocity
    Re   = material['rho'] * v * d_pipe / material['mu'] # Reynolds Number
    return qdot, v, Re




def head_loss(fric, plum_len, v, d_pipe, k_l, accel):
    """ calculates total head loss (major + minor) 

    args:
        fric: skin friction
        plum_len: length of plumbing
        v: fluid velocity (probably bulk)
        d_pipe: pipe_diameter
        k_l: loss coefficient (original notes say for 90 deg flanged elbow, should really be for each)
        accel: acceleration (standard gravity at sea level and not moving)

    returns:
        major plus minor head loss
    """
    # double check these
    h_LMajor = fric * plum_len * v**2 / (d_pipe * 2 * accel) # m Head Major Loss
    h_LMinor = k_l * v**2 / (2 * accel) # m Head Minor Loss
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

def power_req(p_ch, mdot, material, p_in, fric, plum_len, delP_inj, static_head, d_pipe, pump_eff, k_l, u_ss, accel, fos=1.111):
    """Pump power requirements

    args:
        p_ch:
        mdot: mass flow rate
        material: fluid material properties
        p_in:
        fric: friction coefficient
        plub_len: plumbing length
        delP_inj: injector pressure
        static_head: preexisting pressure
        d_pipe: piping diameter
        k_l: loss coefficient for 90 deg flanged elbow
        u_ss: pump suction speed
        accel: acceleration (standard gravity at sea level and not moving)
        fos: factor of safety

    returns:
        v_lfets, p_out, W*fos, rpm*fos
    """
    qdot, v, Re = pipe_flow(material, mdot, d_pipe)
    qdot_lfets, v_lfets, Re_lfets = pipe_flow(material, mdot, d_pipe) # for LFETS

    h_L         = head_loss(fric, plum_len, v, d_pipe, k_l, accel) # Total Head Loss
    p_out       = 1.25 * (p_ch + delP_inj) # Required Outlet Pressure including injector and venturi losses
    delP        = p_out - p_in
    spec_dens   = accel * material['rho']
    h_s         = delP / material['rho'] + h_L * accel # Shaft Work Head times standard gravity
    if h_s <= 0:
        return v_lfets, p_out, 0, 0
    W           = mdot * h_s / pump_eff # power requirement
    tau         = 2 # no less than 2, to avoid cavitation
    inlet_head  = p_in / spec_dens
    vap_p_head  = material['p_v'] / spec_dens
    npsh_a      = inlet_head + static_head - vap_p_head - h_L # available net pressure suction head
    npsh_r      = max(0, npsh_a / tau) # rotational. pardon the kludge...
    rpm         = u_ss * (npsh_r*3.281)**0.75 / (21.2 * np.sqrt(qdot*35.31)) # magic numbers are unit conversion weirdness
    #print(material['name'], ' SH: ',static_head, ' HL: ', h_L)
    return v_lfets, p_out, W * fos, rpm * fos # FOS included




#if __name__ == "__main__":
#    EPSILON_PIPE       = 1.5 *10**(-6) # m Drawn Tubing Relative Roughness
#    R_to_D = EPSILON_PIPE / D_PIPE # roughness to diameter ratio
