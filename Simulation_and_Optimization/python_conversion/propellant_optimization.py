
# current issues: optmization does not bind ipa properly, binds it after optimization


from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel
import rocketcea.cea_obj

import numpy as np
from scipy.optimize import minimize, shgo
import os
import contextlib



global PROPELLANT_SET
PROPELLANT_SET=False




# 
allvectors = []               # array for all design vecs, global variable
allobjfun = []                # array for tracking objective function evaluations


# Optimization Helper Functions

def objective(var, cons):
    return (var/cons)**2 / 2

def exterior(var, cons, good_if_less_than=False):
    if good_if_less_than:
        return np.max([0, var/cons - 1])**2 / 2
    else:
        return np.max([0, -(var/cons - 1)])**2 / 2

def exact(var, cons):
    return (var/cons - 1)**2 / 2

def proportion(amount, ratio):
    top = amount * ratio/(1 + ratio)
    bottom = amount * 1/(1 + ratio)
    return top, bottom





def get_propellant_properties(alc_wt, of_ratio, p_ch, exp_ratio, output=False):
    """ Gets the propellant properties

    args:
        acl_wt: alcholol weight ratio, 0% to 99%
        of_ratio: mixture ratios
        p_ch: chamber pressures
        exp_ratio: list of supersonic area ratios (may be called expansion ratios)
        output: prints the output

    returns:
        isp, IspVac, Cstar, Tc, MW, gamma, string
    """

    ipa_wt = min(alc_wt, 99)
    ipa_str = '''
    fuel C3H8O-2propanol C 3 H 8 O 1    wt%=''' + str(ipa_wt) + '''
    h,cal=-65133.     t(k)=298.15   rho,g/cc=0.786
    fuel water H 2 O 1  wt%=''' + str(100 - ipa_wt) + '''
    h,cal=-68308.  t(k)=298.15 rho,g/cc=0.9998
    '''

    add_new_fuel('LV4_Fuel', ipa_str)
    #add_new_fuel('LV4_Fuel', eth_str)

    if not output:
        PROPELLANT = CEA_Obj(oxName='LOX', fuelName='LV4_Fuel',
                        pressure_units='Pa', temperature_units='K', cstar_units='m/s',
                         density_units='kg/m^3', isp_units='sec', specific_heat_units='J/kg-K')
    else:
        PROPELLANT = rocketcea.cea_obj.CEA_Obj(oxName='LOX', fuelName='LV4_Fuel')
        p_ch *= 1.45e-4

    isp = PROPELLANT.get_Isp(Pc=p_ch, MR=of_ratio, eps=exp_ratio)
    IspVac, Cstar, Tc, MW, gamma = PROPELLANT.get_IvacCstrTc_ChmMwGam(Pc=p_ch, MR=of_ratio, eps=exp_ratio)
    string = PROPELLANT.get_full_cea_output(Pc=p_ch, MR=of_ratio, eps=exp_ratio) if output else ''
    return [isp, IspVac, Cstar, Tc, MW, gamma, string]

def prop_cost(x):
    """Cost function for propellant.

    Args:
        x: values of ipa_wt, of_ratio, p_ch, exp_ratio
            ipa_wt: ipa weight ratio
            of_ratio: mixture ratios
            p_ch: chamber pressures
            exp_ratio: supersonic area ratios (expansion ratios)

    Returns:
        merit: evaluation of propellant
    """

    # expand the current design vector
    ipa_wt, of_ratio, p_ch, exp_ratio = x
    # save the propellant properties
    isp, IspVac, Cstar, Tc, MW, gamma, string = get_propellant_properties(ipa_wt, of_ratio, p_ch, exp_ratio)
    # calcualte the merit of the propellant
    merit = -objective(IspVac, 250) + 100* exact(p_ch, 2413166) + 100 * exact(exp_ratio, 4.5495)
    return merit


def propellant_cost(x, chamber_pressure):
      ipa_wt, of_ratio, p_ch, exp_ratio = x
      isp, IspVac, Cstar, Tc, MW, gamma, string = get_propellant_properties(ipa_wt, of_ratio, p_ch, exp_ratio)
      merit = -objective(IspVac, 250) + 100 * exact(p_ch, chamber_pressure) + 100 * exact(exp_ratio, 4.5495)
      return merit


def custom_prop_cost(x, ref_x):
    """Customizable propellant cost function

    args:
        x: design vector
            ipa_wt: ipa weight ratio (99%)
            of_ratio: mixture ratio
            p_ch: chamber pressure
            exp_ratio: expansion ratio (supersonic area ratio)
        ref_x: reference values for design (dictionary)
            'IspVac': value for vacuum isp
            'p_ch': target chamber pressure
            'exp_ratio': target expansion ratio

    returns:
        merit: merit value
    """
    ipa_wt, of_ratio, p_ch, exp_ratio = x;
    isp, IspVac, Cstar, Tc, MW, gamma, string = get_propellant_properties(ipa_wt, of_ratio, p_ch, exp_ratio)
    merit = -objective(IspVac, ref_x['IspVac']) + 100 * exact(p_ch, ref_x['p_ch']) + 100*exact(exp_ratio, ref_x['exp_ratio'])
    return merit
    

def prop_run():
    # runs optimiztion and saves values
    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
        # x = prop_opt_min([64.8, 1.3, 2413166, 5.988])


        # # may replace other function below
        #res = minimize(custom_prop_cost, [64.8, 1.3, 2413166, 5.988], 
        #               args={'IspVac': 250, 'p_ch': 2413166, 'exp_ratio': 4.5495},
        #               bounds=[[0.1, 100], [0.1, 10], [0.1, 2413166*2], [3, 7]],
        #               method='nelder-mead', options={'adaptive':True})
        
        res = shgo(custom_prop_cost,
               args=[{'IspVac': 250, 'p_ch': 2413166, 'exp_ratio': 4.5495}],
               n=50, iters=2, sampling_method='sobol',
               bounds=[[0.1, 99], [0.1, 10], [0.1, 2413166*2], [3, 7]],
               minimizer_kwargs={'method':'Nelder-Mead', 'options':{'adaptive':True}})

        ipa_wt, of_ratio, p_ch, exp_ratio = res.x
        ipa_wt = min(ipa_wt, 99)
        isp, IspVac, Cstar, Tc, MW, gamma, _ = get_propellant_properties(ipa_wt, of_ratio, p_ch, exp_ratio)
    
    global PROPELLANT_SET
    PROPELLANT_SET = True
    return ipa_wt, of_ratio, p_ch, Tc, MW, gamma, _



def propellant_optimizer(chamber_pressure):
    with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
        res = minimize(custom_prop_cost, 
                       [64.8, 1.3, chamber_pressure, 4.5495], # initial design vector 
                       args={'IspVac': 250, 'p_ch': chamber_pressure, 'exp_ratio': 4.5495}, # normalization or target vector
                       method='nelder-mead', options={'adaptive':True})
        ipa_wt, of_ratio, p_ch, exp_ratio = res.x
        ipa_wt = min(ipa_wt, 99)
        isp, IspVac, Cstar, Tc, MW, gamma, _ = get_propellant_properties(ipa_wt, of_ratio, p_ch, exp_ratio)
    return ipa_wt, of_ratio, p_ch, Tc, MW, gamma, _



if __name__ == '__main__': 

    ipa_wt, of_ratio, p_ch, Tc, MW, gamma, prop_string = prop_run()
    print('Alcohol Wt %: \t', ipa_wt)
    print('OF ratio: \t', of_ratio)
    print('P_ch (Pa): \t', p_ch)
    # print('Expansion ratio: \t', exp_ratio)
    print()
    # print('Vacuum ISP (s): \t', IspVac)
    print('Chamber Temp (K): \t', Tc)
    print('Molar Wt (1/n): \t', MW)
    print('Spec Heat: \t', gamma)


    print()
    ipa_wt, of_ratio, p_ch, Tc, MW, gamma, prop_string = propellant_optimizer(2413166)
    print('Alcohol Wt %: \t', ipa_wt)
    print('OF ratio: \t', of_ratio)
    print('P_ch (Pa): \t', p_ch)
    # print('Expansion ratio: \t', exp_ratio)
    print()
    #print('Vacuum ISP (s): \t', IspVac)
    print('Chamber Temp (K): \t', Tc)
    print('Molar Wt (1/n): \t', MW)
    print('Spec Heat: \t', gamma)

       
