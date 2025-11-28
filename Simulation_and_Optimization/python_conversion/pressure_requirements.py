

# # The purpose of this code is to determine the requirements for pressurant based on the specifications of downstream subsystems. Refer to Huzel and Huang for explanations.

# eq 5-1
# p is pressure (Pa), v is volume (m^3), t_g is temp (K), z is compressibility factor, mm is molar mass (g/mol)
# returns mass (kg). note these quantities are for pressurant.
def mass(p, v, t_g, z, mm):
    return p * v * z * mm / (R_UNIV * t_g)



# # If we neglect heat and mass transfer, the ideal gas law (above) is sufficient. We next consider the case with heat transfer between pressurant and propellant, but not from tank walls.


# eq 5-2
# h is heat transfer coefficient (J s^-1 m^-2 K^-1), a is area (m^2), t is duration (s),
# t_u is temp of gas after expulsion (K), t_e is temp of propellant (K)
# returns total heat transferred from pressurant gas to vaporized propellant (J)
def Q_pres_to_vap_prop(h, a, t, t_u, t_e):
    return h * a * t * (t_u - t_e)

# eq 5-3
# q is heat transferred (J), c_pl is specific heat of liquid propellant (J kg^-1 K^-1),
# h_v is heat of vaporization of propellant (J / kg), c_pv is specific heat of propellant vapor (J kg^-1 K^-1),
# t_v is vaporization temp of propellant (K), t_u is temp of gas after expulsion (K), t_e is temp of propellant (K)
# returns mass of vaporized propellant (kg)
def vap_prop_mass(q, c_pl, h_v, c_pv, t_v, t_u, t_e):
    return q / (c_pl * (t_v - t_e) + h_v + c_pv * (t_u - t_v))

# eq 5-4
# m_v is mass of vaporized propellant (kg), z is compressibility factor of mixture, mm is molar mass of propellant vapor (g/mol),
# t_u is temp of gas after expulsion (K), p is pressure after expulsion (Pa)
# returns volume of vaporized propellant (m^3)
def vap_prop_vol(m_v, z, mm, t_u, p):
    return m_v * z * R_UNIV * t_u / (p * mm)

# eq 5-7
# m is mass of propellant (kg), c_pg is specific heat of pressurant (J kg^-1 K^-1),
# t_u is temp of gas after expulsion (K), q is heat transferred (J)
# returns required mean temperature of entering pressured (K)
def gas_temp(m, c_pg, t_u, q):
    return q / (m * c_pg) + t_u

# eq 5-5, 5-6
# returns required mass and temperature of pressurant assuming heat transfer only to propellant
def pressurant_reqs_2(p_tank, v_tank, z_g, mm_g,
                     h, a, t, t_u, t_e,
                     c_pl, h_v, c_pv, t_v,
                     z_p, mm_p,
                     c_pg):
    q   = Q_pres_to_vap_prop(h, a, t, t_u, t_e)
    m_v = vap_prop_mass(q, c_pl, h_v, c_pv, t_v, t_u, t_e)
    v_v = vap_prop_vol(m_v, z_p, mm_p, t_u, p_tank)
    m_g = mass(p_tank, v_tank - v_v, t_u, z_g, mm_g)
    t_g = gas_temp(m_g, c_pg, t_u, q)
    return m_g, t_g


# include heat transfer
# eq 5-9, 5-10
# be wary of signs of heat transfer, make sure you pick correct sign for situation
# returns required mass and temperature of pressurant
def pressurant_reqs_3(p_tank, v_tank, z_g, mm_g,
                     h, a, t, t_u, t_e,
                     c_pl, h_v, c_pv, t_v,
                     z_p, mm_p,
                     c_pg,
                     q_g_tank, q_tank_p):
    q   = Q_pres_to_vap_prop(h, a, t, t_u, t_e)
    m_v = vap_prop_mass(q + q_tank_p, c_pl, h_v, c_pv, t_v, t_u, t_e)
    v_v = vap_prop_vol(m_v, z_p, mm_p, t_u, p_tank)
    m_g = mass(p_tank, v_tank - v_v, t_u, z_g, mm_g)
    t_g = gas_temp(m_g, c_pg, t_u, q + q_g_tank)
    return m_g, t_g





# 1st order approximation
def n2_prop_reqs(rkt):
    m_g_lox = mass(rkt.lox_tank.p_0, rkt.lox_tank.volume, N2_TEMP, 0.95, N2_MM)
    m_g_ipa = mass(rkt.ipa_tank.p_0, rkt.ipa_tank.volume, N2_TEMP, 0.95, N2_MM)
    return m_g_lox + m_g_ipa

# n2 requirements with heat transfer
def n2_prop_reqs_detailed(sim, z_g, c_pg,
                   z_p_lox, mm_lox, c_pl_lox, h_v_lox, c_pv_lox, t_v_lox, h_lox, t_u_lox,
                   z_p_ipa, mm_ipa, c_pl_ipa, h_v_ipa, c_pv_ipa, t_v_ipa, h_ipa, t_u_ipa):
    m_g_lox, t_g_lox = pressurant_reqs_3(sim.LV4.lox_tank.p_0, sim.LV4.lox_tank.volume, z_g, N2_MM,
                     h_lox, np.pi * sim.LV4.lox_tank.in_radius**2, sim.t[sim.F_index], t_u_lox, 90.18,
                     c_pl_lox, h_v_lox, c_pv_lox, t_v_lox,
                     z_p_lox, mm_lox,
                     c_pg,
                     0, 0)

    m_g_ipa, t_g_ipa = pressurant_reqs_3(sim.LV4.ipa_tank.p_0, sim.LV4.ipa_tank.volume, z_g, N2_MM,
                     h_ipa, np.pi * sim.LV4.ipa_tank.in_radius**2, sim.t[sim.F_index], t_u_ipa, 298.15,
                     c_pl_ipa, h_v_ipa, c_pv_ipa, t_v_ipa,
                     z_p_ipa, mm_ipa,
                     c_pg,
                     0, 0)
    return m_g_lox + m_g_ipa, t_g_lox, t_g_ipa
