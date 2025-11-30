

import numpy as np
import pandas as pd
from customize.system_definition import M_PROP, MDOT, P_E, THROTTLE_WINDOW, MIN_THROTTLE, RCS_MDOT, RCS_P_E, RCS_P_CH
from customize.system_definition import BALLAST, FIN_ROOT, FIN_TIP, FIN_SWEEP_ANGLE, FIN_SEMISPAN, FIN_THICKNESS
from customize.system_definition import CON_NOSE_L, LOX_TANK_P, IPA_TANK_P, RIB_T, NUM_RADL_DVSNS, AIRFRM_IN_RAD
from customize.system_definition import IPA_WT, OF, ENG_P_CH, ENG_T_CH, ENG_KE, ENG_MM
from customize.system_definition import AZ_PERTURB, EL_PERTURB
from trajectory_simulation import trajectory
from models.environment_model import Environment




amount = 20 # number of simulations per run
runs   = 5
# perturbations 'deg N', 'deg E', 'Launch Az', 'Launch El', 'Tip-Off',
            #   'Thrust Pitch', 'Thrust Yaw', 'mdot', 'Ve',
            #    'mass', 'drag', 'Wind'

def mc_init(num_sims):
    results = []
    perturbation_list = []

    for i in range(num_sims):
        perturbation = []
        perturbation.append(np.random.normal(0, 0.001)) # 111 m std dev
        perturbation.append(np.random.normal(0, 0.001)) # 111 m std dev
        perturbation.append(np.random.normal(AZ_PERTURB, 0.333333)) # degrees
        perturbation.append(np.random.normal(EL_PERTURB, 0.333333)) # degrees
        perturbation.append(np.random.rand() < 0.5)
        perturbation.append(np.random.normal(0, 0.0333333)) # degrees
        perturbation.append(np.random.normal(0, 0.0333333)) # degrees
        perturbation.append(np.random.normal(0, 0.02))
        perturbation.append(np.random.normal(0, 0.01))
        perturbation.append(np.random.exponential(2.5))
        perturbation.append(np.random.normal(0, 0.01))
        perturbation.append(np.random.rand() < 0.5)

        perturbation_list.append(perturbation)

    print('random perturbations obtained\n')

    fin_staging = False # temorary

    for i, perturbation in enumerate(perturbation_list):
        if i % 10 == 0: print('iterations:', i)

        sim = trajectory(fin_staging=False, stage_drop_ECEF=0, stage_root=0, stage_tip=0, stage_sweep=0, stage_span=0, stage_thickness=0, mass_red=0, m_prop=M_PROP, mdot=MDOT, 
                         p_e=P_E,
               throttle_window=THROTTLE_WINDOW, min_throttle=MIN_THROTTLE,
               rcs_mdot=RCS_MDOT, rcs_p_e=RCS_P_E, rcs_p_ch=RCS_P_CH,
               ballast=BALLAST, root=FIN_ROOT, tip=FIN_TIP, sweep=FIN_SWEEP_ANGLE, span=FIN_SEMISPAN, thickness=FIN_THICKNESS, con_nose_l=CON_NOSE_L,
               tank_p_o=LOX_TANK_P, tank_p_f=IPA_TANK_P, rib_t=RIB_T, num_radl_dvsns=NUM_RADL_DVSNS,
               airfrm_in_rad=AIRFRM_IN_RAD, ipa_wt=IPA_WT, of=OF, p_ch=ENG_P_CH, T_ch=ENG_T_CH, ke=ENG_KE, mm=ENG_MM,
                          perturbations=perturbation,
                          dt=0.025, adaptive=True, tol=0.005, descend=False, early_return=True, recovery=True)
        x, y, z = sim.raw_states[-1][0][1]
        perturbation.append(x)
        perturbation.append(y)
        perturbation.append(z)

        results.append(perturbation)
    print('done simulations!\n')

    return results




def lat_long(df):
    array = df.loc[:, ['x', 'y', 'z']]
    X = [x for x in array['x']]
    Y = [y for y in array['y']]
    Z = [z for z in array['z']]
    coords = [Environment(None, 17.7, 100, 100, 100).ECEF_to_geodetic([X[i], Y[i], Z[i]]) for i in range(len(array))]
    df['lat'] = [coord[0] for coord in coords]
    df['long'] = [coord[1] for coord in coords]
    df['height'] = [coord[2] for coord in coords]
    
def handle_data(results):
    # landing coordinates
    results_df = pd.DataFrame.from_records(results,
                                           columns=['deg N', 'deg E', 'Launch Az', 'Launch El', 'Tip-Off',
                                                    'Thrust Pitch', 'Thrust Yaw', 'mdot', 'Ve',
                                                    'mass', 'drag', 'Wind',
                                                    'x', 'y', 'z'])
    lat_long(results_df)
    print(results_df.describe())
    return results_df



for i in range(runs):
    results = mc_init(num_sims=amount)
    sim_data = handle_data(results)
    # make sure path exists first
    sim_data.to_csv(path_or_buf='./dispersion_sample_data/sim_data_' + str(i+4) + '.csv')



