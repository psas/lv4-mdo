
from system_definition import *

from structural_model import *
"""
This function is for calculating static loads and bending moments on the fly. We aren't yet doing anything with this information, but it will eventually be necessary to vet our structural designs for safety as we get closer to actual engineering work. I'm still not sure of the best way to explain what is happening here, but worked from examples given to me by A. Fung which are in the Google Drive.

It likely would be a good idea to account for elasticity and uncertainty by multiplying the loads by 2 or 3. I'm not sure whether there are other things I'm missing here (rotary inertia, shear effects?).

"""

def structural_analysis(rkt, env, state, gust):
    FoS = 2
    air, wind         = env.atmo(state[1], gust)
    p_a, rho, T_a, mu = air
    x, q, v, w, _        = state[1:]
    # this is inefficient, but we're only only calling this function when we need to
    force_body, torque_body, v0, dyn_press, Ma, alpha, CoP, fin_flutter, CN, CDax = env.aero_model.aero(state, rkt, air, wind)
    mult              = dyn_press * rkt.frontal_area # for convenience
    # Aerodynamic Normal Forces
    CN_nose           = mult * alpha * env.aero_model.C_N_alpha_nose(alpha)
    CN_body           = mult * alpha * env.aero_model.C_N_alpha_lift(alpha)
    CN_fin            = mult * alpha * env.aero_model.normal_force_coefficient(alpha, Ma)

    # Aerodynamic Axial Forces, unadjusted for angle of attack
    # friction, nose pressure, fin pressure, base drag, parasitic drags
    # friction is whole rocket (for simplicity, assume evenly distributed, though fins account for some)
    # pipes are split between two isogrid tanks
    # lugs are split between lowest two passthru micromodules
    axial_drags       = env.aero_model.C_D_0(env.aero_model.Reynold(v0), Ma)
    # adjustment multiplier from angle of attack
    aoa_cd_mult       = env.aero_model.C_D_axial(alpha, 1)
    # Aerodynamic Axial Forces
    CD_f, CD_p1, CD_p2, CD_b, CD_pipes, CD_lugs = [mult * aoa_cd_mult * drag for drag in axial_drags]

    # Thrust axial force
    thrust            = rkt.engine.thrust(air[0],
                            rkt.engine.throttle_engine(np.linalg.norm(force_body))) if not rkt.empty else 0

    # distribute friction force by length
    # why am i subtracting rings again?
    friction_forces   = CD_f * np.array([part.length / (rkt.length - 10 * COUPLING_RING_HT)
                                         for part in reversed(rkt.parts)])
    # distribute pipe and lug forces to passthru modules (assume evenly for simplicity)
    passthru_modules  = np.array([int('Passthru' in part.name) for part in reversed(rkt.parts)])
    para_drag_per_mod = (CD_pipes + CD_lugs) / sum(passthru_modules)
    parasitic_drags   = passthru_modules * para_drag_per_mod
    other_ax_forces   = np.array([CD_p1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, CD_p2 + CD_b - thrust])
    # array of all axial forces
    axial_forces      = friction_forces + parasitic_drags + other_ax_forces
    axial_forces     *= FoS

    # distribute body lift by length
    # why am i subtracting rings again?
    body_lifts        = CN_body * np.array([part.length / (rkt.length - NOSE_L - 10 * COUPLING_RING_HT)
                                            for part in reversed(rkt.parts)])
    body_lifts[0]     = 0 # nosecone already accounted for
    # other normal forces
    other_normals     = np.array([CN_nose, 0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, CN_fin])
    # array of all normal forces
    normal_forces     = body_lifts + other_normals
    normal_forces    *= FoS

    # Total forces
    total_normal      = sum(normal_forces)
    total_axial       = sum(axial_forces)

    # Total mass
    total_mass        = rkt.mass

    # load factors
    lat_load_factor   = -total_normal / total_mass
    ax_load_factor    = -total_axial / total_mass

    # inertial loads
    intl_lat_loads    = lat_load_factor * np.array([part.mass for part in reversed(rkt.parts)])
    intl_ax_loads     = ax_load_factor * np.array([part.mass for part in reversed(rkt.parts)])

    # combined loads
    cmbn_lat_loads    = normal_forces + intl_lat_loads
    cmbn_ax_loads     = axial_forces + intl_ax_loads

    # load at tops
    lat_load_at_top   = np.zeros(len(rkt.parts))
    for i, load_above in enumerate(lat_load_at_top):
        if i == len(rkt.parts)-1:
            pass
        else:
            lat_load_at_top[i + 1] = cmbn_lat_loads[i] + load_above

    ax_load_at_top    = np.zeros(len(rkt.parts))
    for i, load_above in enumerate(ax_load_at_top):
        if i == len(rkt.parts)-1:
            pass
        else:
            ax_load_at_top[i + 1]  = cmbn_ax_loads[i] + load_above

    bending_moments   = np.zeros(len(rkt.parts))
    for i, element in enumerate(reversed(rkt.parts)):
        if i == len(rkt.parts)-1:
            pass
        else:
            shear_at_top           = lat_load_at_top[i] * element.length
            force_at_middle        = cmbn_lat_loads[i] * (element.CoM[2] - element.coords[2])
            bending_moments[i+1]   = bending_moments[i] + shear_at_top + force_at_middle

    [axial_loading, bend_loading, pressure_loading, stress, hoop_stress,
    axial_stress_ratio, hoop_stress_ratio, axial_safe_margin, hoop_safe_margin] = stress_analysis(rkt, ax_load_at_top, bending_moments, p_a)
    #print('Axial Safety Margins',axial_safe_margin)
    #print('Hoop Safety Margins', hoop_safe_margin)
    return ax_load_at_top, lat_load_at_top, bending_moments

def stress_analysis(rkt, ax_load_at_top, bending_moments, amb_pressure):
    FoS = 1.25
    num = len(rkt.parts)
    axial_loading, bend_loading, pressure_loading = np.zeros(num), np.zeros(num), np.zeros(num)
    hoop_stress, stress, hoop_stress_ratio = np.zeros(num), np.zeros(num), np.zeros(num)
    axial_stress_ratio, axial_safe_margin, hoop_safe_margin = np.zeros(num), np.zeros(num), np.zeros(num)
    for i, element in enumerate(reversed(rkt.parts)):
        cross_section_area = 2 * np.pi * rkt.inr_r * element.thickness
        moment_of_area     = np.pi * rkt.inr_r**3 * element.thickness # approximation
        axial_loading[i]   = FoS * ax_load_at_top[i] / cross_section_area
        bend_loading[i]    = FoS * bending_moments[i] * rkt.out_rad / moment_of_area
        pressure_loading[i] = 0 #FoS * (element.p_0 - amb_pressure) * rkt.out_rad * 0.5 / element.thickness
        if type(element) is Tank:
            pressure_loading[i] = FoS * (element.p_0 - amb_pressure) * rkt.out_rad * 0.5 / element.thickness
            pressure_loading[i] += element.prop_matrl['rho'] * G_N * element.parts[-1].length
        stress[i]           = axial_loading[i] + bend_loading[i] - pressure_loading[i]
        hoop_stress[i]      = -2 * pressure_loading[i]
        axial_stress_ratio[i]     = 100 * abs(stress[i]) / element.material['Sy']
        hoop_stress_ratio[i]     = 100 * abs(hoop_stress[i]) / element.material['Sy']

        if stress[i] == 0:
            axial_safe_margin[i]  = 42069
        else:
            axial_safe_margin[i]  = (100/axial_stress_ratio[i] - 1) * 100
        if hoop_stress[i] == 0:
            hoop_safe_margin[i]  = 42069
        else:
            hoop_safe_margin[i]  = (100/hoop_stress_ratio[i] - 1) * 100

    return [axial_loading, bend_loading, pressure_loading, stress, hoop_stress,
           axial_stress_ratio, hoop_stress_ratio, axial_safe_margin, hoop_safe_margin]
