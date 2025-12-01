
# from .system_definition import R_UNIV

# # The purpose of this code is to determine the requirements for pressurant based on the specifications of downstream subsystems. Refer to Huzel and Huang for explanations.

def mass(p, v, t_g, z, mm, r_univ):
    """Mass of pressurant: Huzel and Huang Eq. 5-1

    args:
        p: pressure (Pa)
        v: volume (m^3)
        t_g: temperature (K)
        z: compressibility factor
        mm: molar mass (g/mol)
        r_univ: universal gas constant

    returns:
        mass (kg) of pressurant
    """
    return p * v * z * mm / (r_univ * t_g)



# # If we neglect heat and mass transfer, the ideal gas law (above) is sufficient. We next consider the case with heat transfer between pressurant and propellant, but not from tank walls.


def Q_pres_to_vap_prop(h, a, t, t_u, t_e):
    """Total heat transfer: Huzel and Huang Eq. 5-2

    args:
        h: heat transfer coefficient (J s^-1 m^-2 K^-1)
        a: area (m^2)
        t: duration (s)
        t_u: temperature of gas after expulsion (K)
        t_e: temperature of the propellant (K)

    returns:
        total heat transferred from pressurant gas to vaporized propellant (J)
    """
 
    return h * a * t * (t_u - t_e)


def vap_prop_mass(q, c_pl, h_v, c_pv, t_v, t_u, t_e):
    """Mass of vaporized propellant: Huzel and Huang Eq. 5-3

    args:
        q: heat transferred (J)
        c_pl: specific heat of liquid propellant (J kg^-1 K^-1)
        h_v: heat of vaporization of propellant (J/kg)
        c_pv: specific heat of propellant vapor (J kg^-1 K^-1)
        t_v: vaporization temperature of propellant (K)
        t_u: temperature of gas aver expulsion (K)
        t_e: temperature of propellant (K)

    returns:
        mass of vaporized propellant (kg)
    """
    return q / (c_pl * (t_v - t_e) + h_v + c_pv * (t_u - t_v))


def vap_prop_vol(m_v, z, mm, t_u, p, r_univ):
    """Volume of vaporized propellant: Huzel and Huang Eq. 5-4

    args:
        m_v: mass of vaporized propellant (kg)
        z: compressibility factor of mixture
        mm: molar mass of propellant vapor (g/mol)
        t_u: temperature of gas aver expulsion (K)
        p: pressure after expulsion
        t_e: temperature of propellant (K)
        r_univ: universal gas constant

    returns:
        volume of vaporized propellant (m^3)
    """
    return m_v * z * r_univ * t_u / (p * mm)


def gas_temp(m, c_pg, t_u, q):
    """Mean temperature of enginering pressured: Huzel and Huang Eq. 5-7

    args:
        m: mass of propellant (kg)
        c_pg: specific heat of pressurant (J kg^-1 K^-1)
        t_u: temperature of gass after expulsion (K)
        q: heat transferred (J)

    returns:
        mean temperature of entering pressured (K)
    """
    return q / (m * c_pg) + t_u

# eq 5-5, 5-6
# returns required mass and temperature of pressurant assuming heat transfer only to propellant
def pressurant_reqs_2(p_tank, v_tank, z_g, mm_g,
                     h, a, t, t_u, t_e,
                     c_pl, h_v, c_pv, t_v,
                     z_p, mm_p,
                     c_pg):
    """Required mass and temperature of propellant assuming heat transfer is only to propellant
    Huzel and Huang eq. 5-5, 5-6

    args:
        p_tank:
        v_tank:
        z_g: compressibility (gas?) 
        mm_g: molar mass (gas?)
        h: heat transfer coefficient
        a: area
        t_u: temperature after expulsion
        t_e: temperature on entry
        c_pl: heat transfer coefficient of [] liquid
        h_v: heat transfer coefficient of vapor
        c_pv: specific heat of pressurant vapor
        t_v: vaporization temperature of propellant
        z_p: compressibility of pressurant
        mm_p: molar mass pressurant
        c_pg: specific heat of pressurant gas
        r_univ: universal gas constant

    returns:
        m_g: mass of pressurant
        t_g: temperature of pressurant
    """
    q   = Q_pres_to_vap_prop(h, a, t, t_u, t_e)
    m_v = vap_prop_mass(q, c_pl, h_v, c_pv, t_v, t_u, t_e)
    v_v = vap_prop_vol(m_v, z_p, mm_p, t_u, p_tank, r_univ)
    m_g = mass(p_tank, v_tank - v_v, t_u, z_g, mm_g, r_univ)
    t_g = gas_temp(m_g, c_pg, t_u, q)
    return m_g, t_g


# include heat transfer
# eq 5-9, 5-10
# be wary of signs of heat transfer, make sure you pick correct sign for situation
# returns required mass and temperature of pressurant
# should just overload the previous function with additional optional arguments
def pressurant_reqs_3(p_tank, v_tank, z_g, mm_g,
                     h, a, t, t_u, t_e,
                     c_pl, h_v, c_pv, t_v,
                     z_p, mm_p,
                     c_pg,
                     q_g_tank, q_tank_p, r_univ):
    q   = Q_pres_to_vap_prop(h, a, t, t_u, t_e)
    m_v = vap_prop_mass(q + q_tank_p, c_pl, h_v, c_pv, t_v, t_u, t_e)
    v_v = vap_prop_vol(m_v, z_p, mm_p, t_u, p_tank)
    m_g = mass(p_tank, v_tank - v_v, t_u, z_g, mm_g, r_univ)
    t_g = gas_temp(m_g, c_pg, t_u, q + q_g_tank)
    return m_g, t_g



# 1st order approximation
def n2_prop_reqs(rkt, n2_temp, n2_z, n2_mm, r_univ):
    """n2 for lox and ipa in rocket

    args:
        rkt: rocket
        n2_temp: n2 gas temperature
        n2_z: n2 gas compressibility factor
        n2_mm: n2 gas molar mass

    returns:
        requirements for lox and ipa
    """
    m_g_lox = mass(rkt.lox_tank.p_0, rkt.lox_tank.volume, n2_temp, n2_z, n2_mm, r_univ)
    m_g_ipa = mass(rkt.ipa_tank.p_0, rkt.ipa_tank.volume, n2_temp, n2_z, n2_mm, r_univ)
    return m_g_lox + m_g_ipa

# n2 requirements with heat transfer
def n2_prop_reqs_detailed(sim, z_g, c_pg,
                   z_p_lox, mm_lox, c_pl_lox, h_v_lox, c_pv_lox, t_v_lox, h_lox, t_u_lox,
                   z_p_ipa, mm_ipa, c_pl_ipa, h_v_ipa, c_pv_ipa, t_v_ipa, h_ipa, t_u_ipa, 
                          r_univ):
    """Pressurant requirements for both lox and ipa"""

    # N2 requirements for LOX
    m_g_lox, t_g_lox = pressurant_reqs_3(p_tank=sim.LV4.lox_tank.p_0, 
                                         v_tank=sim.LV4.lox_tank.volume, 
                                         z_g=z_g, 
                                         mm_g=N2_MM,
                                         h=h_lox, 
                                         a=np.pi * sim.LV4.lox_tank.in_radius**2, 
                                         t=sim.t[sim.F_index], 
                                         t_u=t_u_lox, 
                                         t_e=90.18,
                                         c_pl=c_pl_lox, 
                                         h_v=h_v_lox, 
                                         c_pv=c_pv_lox, 
                                         t_v=t_v_lox,
                                         z_p=z_p_lox, 
                                         mm_p=mm_lox,
                                         c_pg=c_pg,
                                         q_g_tank=0, 
                                         q_tank_p=0, 
                                         r_univ=r_univ)

    # N2 requirements for IPA
    m_g_ipa, t_g_ipa = pressurant_reqs_3(p_tank=sim.LV4.ipa_tank.p_0, 
                                         v_tank=sim.LV4.ipa_tank.volume, 
                                         z_g=z_g, 
                                         mm_g=N2_MM,
                                         h=h_ipa, 
                                         a=np.pi * sim.LV4.ipa_tank.in_radius**2, 
                                         t=sim.t[sim.F_index], 
                                         t_u=t_u_ipa, 
                                         t_e=298.15,
                                         c_pl=c_pl_ipa, 
                                         h_v=h_v_ipa, 
                                         c_pv=c_pv_ipa, 
                                         t_v=t_v_ipa,
                                         z_p=z_p_ipa, 
                                         mm_p=mm_ipa,
                                         c_pg=c_pg,
                                         q_g_tank=0, 
                                         q_tank_p=0,
                                         r_univ=r_univ)

    # total requirements is the sum of the two
    return m_g_lox + m_g_ipa, t_g_lox, t_g_ipa



