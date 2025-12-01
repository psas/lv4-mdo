

# Nelder-Mead simplex search
# %matplotlib inline
# %run OpenRocket_Interface.ipynb
# %run Display_Information.ipynb
# %run Fin_Staging.ipynb
#%run System_Definition.ipynb
#%run Trajectory_Simulation.ipynb


from customize.system_definition import *
from trajectory_simulation import trajectory
from display_information import *
from openrocket_interface import *

# i'm sorry for global vars...
global allvectors, dbz, allobjval
dbz = 0 # arithmetic error tracker (from crossing boundary constraints)
allvectors = []               # array for all design vecs, global variable
allobjfun = []                # array for tracking objective function evaluations
fin_staging = False           # simulate fin staging or not

# Optimize engine thrust to < 6 kN, therefore set to ~ 5.5 kN as constraint is lenient
CONS_THRUST = 5800



####
# Functions of merit
####

# some of these are copies from pressure requirements


# all of our comparisons are ratios instead of subtractions because
# it's normalized, instead of dependent on magnitudes of variables and constraints

# minimize this, **2 makes it well behaved w.r.t. when var=cons
def objective(var, cons):
    return (var/cons)**2 / 2

def objective_additive(var, cons):
    return np.linalg.norm(var - cons)**2 / 2

# **2 because i like it more than abs(), but that also works
def exact(var, cons):
    return (var/cons - 1)**2 / 2

# this is your basic exterior penalty, either punishes for unfeasibility or is inactive
def exterior(var, cons, good_if_less_than=False):
    if good_if_less_than:
        return np.max([0, var/cons - 1])**2 / 2
    else:
        return np.max([0, -(var/cons - 1)])**2 / 2

# this barrier function restricts our objective function to the strictly feasible region
# make rockets great again, build that wall, etc, watch out for undefined operations
def barrier(var, cons, int_point=False, good_if_less_than=True):
    global dbz
    try: # just in case we accidentally leave feasible region
        if not int_point:
            if good_if_less_than:
                return -log(-(var/cons - 1))
            else:
                return -log(var/cons - 1)
        elif int_point:
            def interior(g): return 1/g # in case we don't like logarithms, which is a mistake
            if good_if_less_than:
                return -interior(var/cons - 1)
            else:
                return -interior(-(var/cons - 1))
    except:
        dbz += 1 # keep track of arithmetic errors, side effect
        return float('inf') # ordinarily, this is bad practice since it could confuse the optimizer
                            # however, since this is a barrier function not an ordinary penalty, i think it's fine




#############
# Optimization problem
#############


# this manages all our constraints
# penalty parameters: mu -> 0 and rho -> infinity
def penalty(sim, mu, rho):
    # barrier penalties have "less smooth" behavior
    b = [#barrier(#sim.alt[-1], CONS_ALT, int_point=False, good_if_less_than=False),
         #barrier(#sim.alt[-1], CONS_CEILING, int_point=False, good_if_less_than=True),
         #-log(sim.LV4.ballast),
         #barrier(sim.min_fin_flutter, 1.0, int_point=False, good_if_less_than=False)
        ]
    eq = []
    ext = [exterior(sim.alt[-1], CONS_ALT, good_if_less_than=False),
           exterior(sim.alt[-1], CONS_CEILING, good_if_less_than=True),
           #exterior(sim.min_fin_flutter, 1.0, good_if_less_than=False),
           exterior(sim.thrust[0], CONS_THRUST, good_if_less_than=True),
           exterior(sim.LV4.v_lfets_o, CONS_V_LFETS, good_if_less_than=True),
           exterior(sim.LV4.v_lfets_f, CONS_V_LFETS, good_if_less_than=True),
           #EFS removed
           #exterior(sim.LV4.pow_o, CONS_EFS, good_if_less_than=True),
           #exterior(sim.LV4.pow_f, CONS_EFS, good_if_less_than=True),
           # these dynamic constraints are reflective of the statistical correlation
           # between power and speed of motors commercially available. R^2 = 0.87, which is approximately 1.
           exterior(sim.LV4.rpm_o, 101003 - 6.96 * sim.LV4.pow_o, good_if_less_than=True),
           exterior(sim.LV4.rpm_f, 101003 - 6.96 * sim.LV4.pow_f, good_if_less_than=True),
           #exterior(sim.tip_off_aoa, CONS_AOA, good_if_less_than=True),
           #exterior(sim.launch_speed, CONS_LS, good_if_less_than=False),
           #exterior(sim.ld_ratio, CONS_LD, good_if_less_than=True),
           exterior(sim.TWR, CONS_TWR, good_if_less_than=False),
           exterior(sim.S_crit, CONS_S_CRIT, good_if_less_than=False),
           exterior(sim.max_g_force, CONS_ACCEL, good_if_less_than=True),
           #exterior(sim.min_stability, CONS_STBLTY, good_if_less_than=False),
           exterior(sim.impulse, CONS_IMPLS, good_if_less_than=True),
           exterior(sim.LV4.lox_tank.p_0, sim.LV4.TANK_MAX_P, good_if_less_than=True),
           exterior(sim.LV4.ipa_tank.p_0, sim.LV4.TANK_MAX_P, good_if_less_than=True),
           exterior(sim.LV4.lox_tank.p_0, CONS_TANK_MIN, good_if_less_than=False),
           exterior(sim.LV4.ipa_tank.p_0, CONS_TANK_MIN, good_if_less_than=False)
          ]
    return mu*sum(b) + rho*(sum(eq) + sum(ext))

# Pseudo-objective merit function
# x is array of design parameters, n is index of penalty and barrier functions
# print blocks are sanity checks so i'm not staring at a blank screen and can see what various tweaks actually do
def f(x, n=8):
    global allvectors, allobjfun
    ipa_wt, of_ratio, p_ch, Tc, MW, gamma, _ = propellant_optimizer(x[5])
    # get trajectory data
    sim = trajectory(fin_staging, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2],
               THROTTLE_WINDOW, MIN_THROTTLE,
               RCS_MDOT, RCS_P_E, RCS_P_CH,
               x[9], x[6], x[7], FIN_SWEEP_ANGLE, x[8], FIN_THICKNESS, CON_NOSE_L,
                x[3], x[4], RIB_T, NUM_RADL_DVSNS,
               AIRFRM_IN_RAD, ipa_wt, of_ratio, x[5], Tc, gamma, MW,
               [0, 0, AZ_PERTURB, EL_PERTURB, True, 0, 0, 0, 0, 0, 0, True],
                0.025, True, 0.005, False, False, False)

    obj_func = (1 * objective(sim.LV4.GLOW, CONS_MASS))
                # + 0.4 * (np.linalg.norm(sim.thrust[sim.F_index])/8000)**10)
                # base 11 constraint + 0.1 * objective(x[0]/x[1], 40)) # minimize GLOW
    # then, calculate penalization from trajectory performance
    pen_func = penalty(sim, MU_0 / (2**n), RHO_0 * (2**n)) # initial mu and rho selected for nice behavior
    # add objective and penalty functions
    merit_func = obj_func + pen_func
    if np.isnan(merit_func): merit_func = np.inf # crude hack
    allvectors.append(x) # maintains a list of every design, side effect
    allobjfun.append(log(merit_func)) # log plot is cleaner to look at, for later
    return merit_func

# we want to iterate our optimizer for theoretical "convergence" reasons (given some assumptions)
# n = number of sequential iterations
def iterate(f, x_0, n):
    x = x_0 # initial design vector
    global dbz
    designs = []
    for i in range(n):
        print("Iteration " + str(i+1) + ":")
        res = minimize(f, x, args=(i+3), method='nelder-mead', options={'disp': True, 'adaptive':True, 'xatol': 1, 'fatol': 0.1})
        x = res.x # feed optimal design vec into next iteration

        designs.append(res.x)   # we want to compare sequential objectives
                                # so we can stop when convergence criteria met
        alt = trajectory(fin_staging, stage_drop_ECEF, stage_root, stage_tip, stage_sweep, stage_span, stage_thickness, mass_red, x[0], x[1], x[2],
               THROTTLE_WINDOW, MIN_THROTTLE,
               RCS_MDOT, RCS_P_E, RCS_P_CH,
               BALLAST, FIN_ROOT, FIN_TIP, FIN_SWEEP_ANGLE, FIN_SEMISPAN, FIN_THICKNESS, CON_NOSE_L,
                x[3], x[4], RIB_T, NUM_RADL_DVSNS,
               AIRFRM_IN_RAD, IPA_WT, OF, ENG_P_CH, ENG_T_CH, ENG_KE, ENG_MM,
               [0, 0, AZ_PERTURB, EL_PERTURB, False, 0, 0, 0, 0, 0, 0, True],
                          0.025, True, 0.005, True, False, False).alt

        print("         Arithmetic errors (from violations of acceptable altitude window): "+str(dbz))
        print("Propellant mass (kg): "+str(x[0]))
        print("Mass flow rate (kg/s): "+str(x[1]))
        print("Exit pressure (Pa): "+str(x[2]))
        print("Peak Altitude (km): "+str(alt[-1]/1000))
        print('')
        dbz=0 # I only care about divisions by zero in each individual iteration, side effect
        if (i > 0) and (np.linalg.norm(designs[-1] - designs[-2]) < DELTA):
            print("Early termination! The Euclidean distance between the last two designs was < " + str(DELTA))
            break
    return x

def breed_rockets(f):
    res = differential_evolution(f, [(80, 200), (2.0, 4.5), (45000, 150000)],
                                 strategy='best1bin', popsize=80, mutation=(.1, .8), recombination=.05,
                                 updating='immediate', disp=True, atol=0.1, tol=0.1)
                                 #polish=True,workers=-1, disp=True, atol=0.1, tol=0.1)
    return res.x

def benchmark(x):
        return x[0]*x[1] - x[2]


################
# Optimization
###############



import os
import yaml

# Results, this is the big boi function
if __name__ == '__main__':

    # Get bonmin and ipopt paths from the path file
    #path_file_dir = os.path.dirname(os.path.realpath("__file__"))
    #path_file = open(os.path.join(path_file_dir, "bonmin_paths.yaml"))
    #paths = yaml.safe_load(path_file)

    init_array = [M_PROP, MDOT, P_E, LOX_TANK_P, IPA_TANK_P, ENG_P_CH]
    print("Began optimization at", datetime.now())
    #np.seterr(all='raise')

    # RBFOpt
    """
    bb = rbfopt.RbfoptUserBlackBox(len(init_array),
                                   np.array([75, 1.5, 86346/4, 805000/2.5, 805000/2.5, 344738*1.5]),
                                   np.array([375, 7, 86346*2.5, 2.068e6, 2.068e6, 4.137e6]), # design vector boundaries
                               np.array(['R']*len(init_array)), f)
    settings = rbfopt.RbfoptSettings(minlp_solver_path=paths["bonmin_path"],
                                     nlp_solver_path=paths["ipopt_path"],
                                    max_evaluations=800, eps_impr=1.0e-3)
    alg = rbfopt.RbfoptAlgorithm(settings, bb)
    val, x, itercount, evalcount, fast_evalcount = alg.optimize()
    """

    # iterative nelder-mead
    # feed initial design into iterative optimizer, get most (locally) feasible design
    '''test = trajectory(init_array[0], init_array[1], init_array[2],
               THROTTLE_WINDOW, MIN_THROTTLE,
               RCS_MDOT, RCS_P_E, RCS_P_CH, BALLAST, FIN_ROOT, FIN_TIP, FIN_SWEEP_ANGLE, FIN_SEMISPAN, FIN_THICKNESS, CON_NOSE_L,
               LOX_TANK_P, IPA_TANK_P, RIB_T, NUM_RADL_DVSNS,
               AIRFRM_IN_RAD, IPA_WT, OF, ENG_P_CH, ENG_T_CH, ENG_KE, ENG_MM,
               [0, 0, 0, 0, True, 0, 0, 0, 0, 0, 0, True],
                          0.05, True, 0.045, False, False)

    if (test.alt[-1] < CONS_ALT) or (test.alt[-1] > CONS_CEILING): # rudimentary error handling to save heartache
        raise Exception('Rocket apogee out of bounds! Apogee {:.3f} km'.format(test.alt[-1]/1000))
    x = iterate(f, init_array, ITERATIONS)'''


    # probe design space, darwin style. takes forever, literally.
    #res = breed_rockets(f)

    # simplicial homology

    res = shgo(f, bounds=[*zip([75,2,86346/1.5, 805000/2,805000/2, 344738*1.5, 0.8, 0.5, 0.15, 8],
                                [300,6,86346*1.5,2.068e6, 2.068e6, 4.137e6, 1.8, 0.9, 0.9, 40])],   # design vector boundaries
                   n=10, iters=1, sampling_method='simplicial',
                    #minimizer_kwargs={'method':'SLSQP', 'options':{'disp': True,'maxiter':250}},
               minimizer_kwargs={'method':'COBYLA', 'options':{'disp': True, 'adaptive':True, 'maxfev':500, 'xatol': 0.01, 'fatol': 0.05}},
                  options={'disp':True})
    x = res.x
    xl = res.xl

    print('Global min:', x)
    print('Local min:', xl)
    print('Function values:', res.funl)


    print("Ended optimization at", datetime.now())
    print("Function evaluations:",len(allvectors))
    print("Optimization done!")

    ipa_wt, of_ratio, p_ch, Tc, MW, gamma, _ = propellant_optimizer(x[5])

    if fin_staging:
        smaller_fin_sim = trajectory(fin_staging, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2],
                   THROTTLE_WINDOW, MIN_THROTTLE,
                   RCS_MDOT, RCS_P_E, RCS_P_CH,
                   x[9], x[6], x[7], FIN_SWEEP_ANGLE, x[8], 0.003175, CON_NOSE_L, #enter smaller fin parameters
                    x[3], x[4], RIB_T, NUM_RADL_DVSNS,
                   AIRFRM_IN_RAD, ipa_wt, of_ratio, x[5], Tc, gamma, MW,
                   [0, 0, AZ_PERTURB, EL_PERTURB, False, 0, 0, 0, 0, 0, 0, True],
                                     0.025, True, 0.005, True, False, False)

        drop_time = find_drop_time(smaller_fin_sim.stability_margin)
        time_diff = []
        for item in smaller_fin_sim.t:
            time_diff.append(abs(item - drop_time))
        min_index = time_diff.index(min(time_diff))

        def altitude_to_ECEF(alt):
            return (alt + 6.371e3)

        drop_ECEF = altitude_to_ECEF(smaller_fin_sim.alt[min_index])
    else:
        drop_ECEF = 0

    # get trajectory info from optimal design
    sim = trajectory(fin_staging, drop_ECEF, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2],
               THROTTLE_WINDOW, MIN_THROTTLE,
               RCS_MDOT, RCS_P_E, RCS_P_CH,
               x[9], x[6], x[7], FIN_SWEEP_ANGLE, x[8], FIN_THICKNESS, CON_NOSE_L, #enter larger fin parameters
               x[3], x[4], RIB_T, NUM_RADL_DVSNS,
               AIRFRM_IN_RAD, ipa_wt, of_ratio, x[5], Tc, gamma, MW,
               [0, 0, AZ_PERTURB, EL_PERTURB, True, 0, 0, 0, 0, 0, 0, False],
                          0.025, True, 0.005, False, False, False)

    if fin_staging:
        time_diff_2 = []
        for item in sim.t:
            time_diff_2.append(abs(item - drop_time))
        min_index_2 = time_diff_2.index(min(time_diff_2))
        print("Fin drop altitude: ", smaller_fin_sim.alt[min_index]/1000, "km")
        print("Mach number upon fin dropping: ", sim.Ma[min_index_2])

    textlist = print_results(sim, True)
    # draw pretty pictures of optimized trajectory
    rocket_plot(sim.t, sim.alt, sim.v, sim.a, sim.thrust,
                sim.dyn_press, sim.Ma, sim.stability_margin, sim.m, sim.p_a, sim.drag, sim.throttle, sim.fin_flutter, sim, True, None, None)
    Init_Stability(sim)
    # structural analysis
    structural_plot(sim.LV4)
    print(sim.env.ECEF_to_geodetic(sim.raw_states[-1][0][1]))
    print(sim.apogee)
    print(max([state[1][1][5][5] * 180 / np.pi
       for state in sim.raw_states[sim.LV4.tower_index:sim.LV4.F_index]]))
    print(sim.tip_off_aoa)
    # get/print info about our trajectory and rocket
    for line in textlist:
        print(line[1:])
    save_trajectory(sim)

    print('\nMaking an OpenRocket rocket and corresponding engine!')
    # create an openrocket file with matching engine for our design (and print/save trajectory data)
    make_engine(x[1], sim.m_prop, sim.thrust[0:sim.F_index + 1],
                sim.LV4.inr_r, AIRFRAME_THICKNESS, sim.LV4.l_o, sim.LV4.l_f, sim.LV4.m_tank_o, sim.LV4.m_tank_f,
                sim.t[sim.F_index], sim.LV4.engine.Ve/G_N,
                sim.LV4.eng_sys_dry_mass, sim.LV4.eng_sys_len, sim.openrocket_CoM,
                sim.LV4.ballast, sim.LV4.fin)

    # draw more pretty pictures, but of the optimizer guts
    design_grapher(allvectors, allobjfun)
