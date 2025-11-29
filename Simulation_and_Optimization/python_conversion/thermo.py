
"""

The intention of this notebook is to assemble in one place all the functions required to do aerothermal analyses of the airframe in flight and thermal analyses of the engine.

The strategy chosen for the MDO is to create thermal circuits as in "Heat Transfer Modeling: An Inductive Approach" by George Sidebotham. This will be both accurate and fast enough for our system-level purposes.

The heat transfer model of the engine is yet to come.

Currently, I have begun implementing the functions from "Real-Time Aerodynamic Heating and Surface Temperature Calculations for Hypersonic Flight Simulation" by Quinn and Gong, and then I plan to look at "A Method for Calculating Transient Surface Temperatures and Surface Heating Rates for High-Speed Aircraft" by the same authors.

"""




import numpy as np
from scipy import interpolate

def stage_2(Delta, S, alpha_hat, Sc, Ste):
    factor = 2 * alpha_hat / (Delta - S)
    Delta_dot = factor * (3 + 2 * Sc) - 2 * Ste
    S_dot     = Ste - factor * Sc
    return Delta_dot, S_dot



Temp_array = [0,450,900,1350,1800,
              2250,2700,3150,3600,4050,
              4500,4950,5400,5850,6300,
              6750,7200,7650,8100,8550,
              9000,9450,9900,10350,10800,
              11250,11900,12150,12600,13050,
              13500,13750]

# lumped system analysis of heat transfer, assuming skin thin enough for uniform temperature
def Tw_dot(T_w, rho_w, c_p_w, thickness, h, H_ext, H_w, beta):
    radiation = -beta * T_w**4
    convect   = h * (H_ext - H_w)
    #print(H_ext , H_w)
    coeff     = rho_w * c_p_w * thickness
    return (radiation + convect) / coeff

# heat transfer coeff for 3d flow
def h_3d(K1, rho_mu_st, rho_w, mu_w, grad_U):
    return 0.94 * K1 * rho_mu_st**0.4 * (rho_w * mu_w)**0.1 * np.sqrt(grad_U)

# stagnation enthalpy
def H_st(H_inf, U_inf, cos_Lambda_square):
    return H_inf + U_inf**2 * cos_Lambda_square * 0.25 #/ 50103

# stagnation velocity gradient
def grad_U(R, M_cos_Lam_inf_square, P_inf, rho_inf):
    term1 = max(0,7 * M_cos_Lam_inf_square - 1)
    return 1/R * np.sqrt(term1 * (M_cos_Lam_inf_square + 5) / M_cos_Lam_inf_square * (P_inf/rho_inf))

# stagnation density times dynamic viscosity
def rho_mu_st(rho_inf, mu_inf, M_cos_Lam_inf_square, T_st, T_inf):
    return rho_inf * mu_inf * 6 * M_cos_Lam_inf_square / (M_cos_Lam_inf_square + 5) * (T_st / T_inf)**0.75

# dynamic viscosity
def mu(T):
    return 1.05e-7 * T**0.75

# static pressure
def P(P_inf, M_cos_Lam_inf_square):
    return max(0, P_inf * (7 * M_cos_Lam_inf_square - 1) / 6)

#def rho_w(P_w, T_w):
#    return P_w / (53.3 * T_w)

# T is Rankine, P in lb/ft^2, H in BTU/lbm


def flatten_nested_list(nested_list):
    """Effectively reshapes the nested list to an N-by-1 list"""
    flat_list = list()

    for the_list in nested_list:
        flat_list.extend(the_list)

    return flat_list


air_enthalpy_x = np.array([21160.0, 2116.0, 211.6, 21.16, 2.116, 0.212]*32)
air_enthalpy_y = np.array(flatten_nested_list([[val]*6 for val in Temp_array]))

air_enthalpy_data = [[0 for i in range(6)], [108 for i in range(6)], [217 for i in range(6)],
                                      [334 for i in range(6)], [451 for i in range(6)], [578 for i in range(6)],
                                      [704 for i in range(6)], [836,836,836,836,836,864], [968,968,971,978,1004,1091],
                                      [1111,1111,1131,1175,1280,1625], [1244,1264,1322,1486,1894,2771],
                                      [1379,1438,1595,1984,2682,2838], [1578,1706,2040,2641,2968,3038],
                                      [1788,2017,2456,3016,3170,3236], [2051,2481,2970,3279,3360,3467],
                                      [2352,2874,3346,3498,3637,3938],
                                      [2740,3321,3554,3724,3995,4844], [3103,3659,3841,4068,4644,6688],
                                      [3471,3807,4081,4514,5849,9305], [3782,4132,4441,5262,7890,12036],
                                      [4013,4400,4916,6457,10261,14452], [4392,4710,5730,8201,12573,15409],
                                      [4669,5147,6599,10260,14619,15951], [4970,5688,7952,12347,15723,16426],
                                      [5280,6435,9548,14252,16210,16825], [5712,7341,11369,15578,16681,17187],
                                      [6204,8521,13171,16259,17112,17603], [6772,9760,14510,16757,17520,18383],
                                      [7564,11155,15783,17095,17908,19265], [8498,12789,16600,17637,18527,20573],
                                      [9452,14267,17170,18163,19205,22421], [10623,15508,18181,18551,20172,24979]]

air_enthalpy_data = np.array(flatten_nested_list(air_enthalpy_data))

air_enthalpy = interpolate.LinearNDInterpolator(list(zip(air_enthalpy_x, air_enthalpy_y)), air_enthalpy_data)




K1 = interpolate.interp1d([0,1,5,10], [1,1,1.16,1.14])


def enthalpy(T, P):
    H_BTU = air_enthalpy(P * 0.02089, T * 1.8) # 1 Pa = 0.02089 lb/ft^2, 1 K = 1.8 R
    return H_BTU[0] * 2326 # 1 BTU/lbm = 2326 J/kg

def wall_temperature(T_w,
                     thickness, R, rho_w, c_p_w, beta,
                     T_inf, P_inf, rho_inf, U_inf, Ma_inf, mu_inf):
    lam      = 0 # leading-edge sweep angle
    c_Lam_sq = np.cos(lam)**2
    M_lam_sq = c_Lam_sq * Ma_inf**2
    
    P_w_st   = P(P_inf, M_lam_sq)
    
    H_inf    = enthalpy(T_inf, P_inf)
    H_ext    = H_st(H_inf, U_inf, c_Lam_sq)
    H_w      = enthalpy(T_w, P_w_st)
    
    H_st_arr = [abs(H_ext - enthalpy(T_i, P_w_st)) for T_i in Temp_array]
    H_st_ndx = np.argmin(H_st_arr)
    T_st     = Temp_array[H_st_ndx] # can make this more accurate if I really have to
    
    mu_inf   = mu(T_inf)
    rhomu_st = rho_mu_st(rho_inf, mu_inf, M_lam_sq, T_st, T_inf)
    gradU    = grad_U(R, M_lam_sq, P_inf, rho_inf)
    #rhow     = rho_w(P_w_st, T_w)
    mu_w     = mu(T_w)
    h        = h_3d(K1(Ma_inf), rhomu_st, rho_w, mu_w, gradU)
    #k = 235 # W /(m K)
    #biot = h * thickness / k
    #print(biot)
    
    Twdot    = Tw_dot(T_w, rho_w, c_p_w, thickness, h, H_ext, H_w, beta)
    
    #print(Twdot, rhow, h, H_ext, H_w)
    #K = 112.4 * 0.31832 #
    #qw = K * np.sqrt(P_inf/R) * (H_inf - H_w) -beta * T_w**4
    #Twdot1 = qw / (rho_w * c_p_w * thickness)
    #print(Twdot, Twdot1)
    #qw = 90 * np.sqrt(P_w_st/R) * ((H_ext - H_w)*1e-6)**1.17 -beta * T_w**4
    return Twdot



'''# heat transfer coefficient for cone
def h(C5, h_o, h_alpha):
    return C5 * (h_o + h_alpha)

# corrects 2d coefficients to conical flow
def C5(turbulent):
    return 1.15 if turbulent else 1.73

# flat plate at 0 angle of attack
def h_o(turbulent, C1, C3, rho_U_inf, mu_inf, x, T_inf, T_star):
    T_T   = T_inf / T_star
    if turbulent:
        return C1 * 0.0375 * rho_U_inf**0.8 * mu_inf**0.2 * x**(-0.2) * T_T **0.65
    else:
        return C3 * 0.421 * (rho_U_inf * mu_inf / x)**0.5 * T_T**0.125

# angle of attack and wedge/cone angles
def h_alpha(turbulent, A1, A2, rho_U_inf, x, delta, alpha):
    angle = delta + alpha # ignore lower/upper surface distinction
    if turbulent:
        return A1 * rho_U_inf**0.8 * x**(-0.2) * angle
    else:
        return A2 * (rho_U_inf / x)**0.5 * angle


rho_U_inf = rho_inf * U_inf'''
