"""
This is adapted from a script created by natronics, 7deeptide, and jalouke to create engine files for OpenRocket. Its original purpose was to determine some crucial properties of the engine system based on initial assumptions and specifications a few years before the multidisciplinary design optimization for LV4 began.

It has since been reformed to take its parameters from the MDO model in order to utilize OpenRocket as a downstream verifier of potential designs. It should be noted that OpenRocket was not made with high-powered rockets in mind, and so one should not be too surprised if the final apogee differs between OpenRocket and our trajectory model by as much as 20 km.

Basically, this script creates an OpenRocket engine file, which is a black box with a time-dependant thrust curve and which accounts for its own mass properties and center of mass. T
"""


"""
Our optimization will spit out a list of relevant details, which is passed to this function so we can save everything potentially useful in one nice way.

Our rocket template is an xml file set to use one node to determine the diameter of each subsystem so we just have to climb that tree to change out the airframe radius and the two tank lengths and viola! Note, this effects just the rocket file, it is the next function that creates the engine file.

There is probably a nice way to use the imported XML library instead of manually inserting text, but I don't quite care enough since it's perfectly functional. This takes its parameters straight from the trajectory simulation of an optimized rocket. It doesn't do anything particularly exciting, just creates a thrust curve and calculates some odds and ends that OpenRocket likes. It's unfortunate that thrust in OpenRocket doesn't account for altitude, but the discrepancy in final apogee is usually no more than 2 or 3 km (our trajectory simulation seems to be conservative).

Be aware that I have not yet tested this in Windows or MACOSX, or outside of Ubuntu 18 for that matter.
"""

from customize.system_definition import *
from trajectory_simulation import trajectory
from display_information import *

# create a rocket file for our engine's dimensions and characteristics
def update_body(index, out_r, l_o, l_f, ballast, fin):
    unzip() # unpack template

    with open('rocket.ork', 'rb') as xml_file:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for child in root.iter():
            #set radius, this node propagates itself in openrocket
            for kid in child.iterfind('aftradius'):
                kid.text = str(out_r)
            for kid in child.iterfind("*[name='Tip Ballast']"):
                kid.find('mass').text = str(ballast) # set ballast mass
            for kid in child.iterfind("*[name='Fuel Tank']"):
                kid.find('length').text = str(l_f) # set fuel tank length
            for kid in child.iterfind("*[name='LOX Tank']"):
                kid.find('length').text = str(l_o) # set lox tank length
            for kid in child.iterfind("*[name='Trapezoidal fin set']"):
                kid.find('rootchord').text = str(fin.root) # set fin geometry
                kid.find('tipchord').text = str(fin.tip)
                kid.find('sweeplength').text = str(fin.sweep_length)
                kid.find('height').text = str(fin.semispan)
                kid.find('thickness').text = str(fin.thickness)

        tree.write('rocket.ork')

    zipit(index) # repack template
    print("New rocket generated from template!")

# when a father engine and a mother engine love each other very much...
# they make an openrocket engine
# that is approximately equivalent to the trajectory profile from the optimimization
def make_engine(mdot, prop_mass, Thrust,
                inr_radius, thickness,
                l_o, l_f, m_tank_o, m_tank_f,
                Burn_time, Isp, sys_mass, length, CoM,
                ballast, fin):

    index = str(get_index())

    radius = inr_radius + thickness
    eng_sys_dry_mass = sys_mass

    n = len(Thrust)
    peak = max(Thrust)
    average = float(sum(Thrust) / n)

    file_head = """<engine-database>
  <engine-list>
    <engine  mfg="PSAS" code="{code}" Type="Liquid" dia="{diameter}" len="{length}"
    initWt="{total_mass}" propWt="{M_prop}" delays="0" auto-calc-mass="0" auto-calc-cg="0"
    avgThrust="{a_thrust}" peakThrust="{p_thrust}" throatDia="0." exitDia="0." Itot="{impulse}"
    burn-time="{burn_time}" massFrac="{m_frac}" Isp="{Isp}" tDiv="10" tStep="-1." tFix="1"
    FDiv="10" FStep="-1." FFix="1" mDiv="10" mStep="-1." mFix="1" cgDiv="10"
    cgStep="-1." cgFix="1">
    <comments>Optimized engine</comments>
    <data>
""".format(**{'code': 'PSAS '+index,
              'diameter': radius*2*1000,
                  'length': length*1000,
                  'total_mass': (prop_mass[0] + eng_sys_dry_mass)*1000,
                  'M_prop': prop_mass[0]*1000,
                  'a_thrust': average,
                  'p_thrust': peak,
                  'burn_time': Burn_time,
                  'm_frac': eng_sys_dry_mass/(eng_sys_dry_mass + prop_mass[0]),
                  'impulse': average*Burn_time,
                  'Isp': Isp,
        })

    data = [] # this is going to be our thrust curve!
              # may be slightly inaccurate if/when altitudes vary between trajectory.py and openrocket
    resolution = float(Burn_time/n) # sec per step
    for i in range(n):
        data.append('     <eng-data  t="{t}" f="{thrust}" m="{mass}" cg="{cg}"/>\n'.format(**{
            't': i * resolution, # sec
            'thrust': Thrust[i],
            'mass': (eng_sys_dry_mass + prop_mass[i]) * 1000,
            'cg': (length - CoM[i]) * 1000
        }))

    file_tail = """
    </data>
  </engine>
</engine-list>
</engine-database>"""

    # we're gonna put the engine file in the default location
    prefix = "./"
    # # this code is old and needs to be updated
    #if 'linux' in _platform:
    #    home = os.path.expanduser("~")
    #    prefix =  os.path.join(home, '.openrocket/ThrustCurves/')
    #elif "darwin" in _platform:
    #    home = os.path.expanduser("~")
    #    prefix =  os.path.join(home, 'Library/Application Support/OpenRocket/ThrustCurves/')
    #elif "win" in _platform:
    #    home = os.getenv("APPDATA")
    #    prefix = os.path.join(home, "OpenRocket/ThrustCurves/")

    if not os.path.isdir(prefix):
        os.mkdirs(prefix) # in case openrocket not installed to prevent crash

    # now write the file, great job!
    with open(os.path.join(prefix, 'psas_motor_'+index+'.rse'), 'w') as eng:
        eng.write(file_head)
        for d in data:
            eng.write(d)
        eng.write(file_tail)

    update_body(index, radius, l_o, l_f, ballast, fin) # make a rocket to correspond with our new engine
    print("Reopen OpenRocket to run simulation.")



# This is for in case you want to use this program without going through the optimizer first.


if __name__ == '__main__':
    m_prop    = X0[0]
    mdot = X0[1]
    p_e  = X0[2]

    sim = trajectory(fin_staging=False, 
                          stage_drop_ECEF=0, 
                          stage_root=0, 
                          stage_tip=0, 
                          stage_sweep=0, 
                          stage_span=0, 
                          stage_thickness=0, 
                          mass_red=0, 
                          m_prop=M_PROP, 
                          mdot=6.0758, 
                          p_e=57563.9331, 
                          throttle_window=THROTTLE_WINDOW, 
                          min_throttle=MIN_THROTTLE, 
                          rcs_mdot=RCS_MDOT, # JUST FOR ROCKET 
                          rcs_p_e=RCS_P_E, # JUST FOR ROCKET 
                          rcs_p_ch=RCS_P_CH, # JUST FOR ROCKET
                          ballast=BALLAST, # BALLAST    fin optimization stuff
                          root=FIN_ROOT, # fin optimization stuff
                          tip=FIN_TIP, # fin optmization stuff
                          sweep=FIN_SWEEP_ANGLE, # fin sweep angle for stuff
                          span=FIN_SEMISPAN, # fin span for other stuff
                          thickness=FIN_THICKNESS, # fin thickness for other stuff
                          con_nose_l=CON_NOSE_L, # JUST FOR ROCKET
                          tank_p_o=LOX_TANK_P, # JUST FOR ROCKET
                          tank_p_f=IPA_TANK_P, # JUST FOR ROCKET 
                          rib_t=RIB_T, # JUST FOR ROCKET
                          num_radl_dvsns=NUM_RADL_DVSNS, # JUST FOR ROCKET
                          airfrm_in_rad=AIRFRM_IN_RAD, # JUST FOR ROCKET
                          ipa_wt=IPA_WT, # propellant optimization for cea output
                          of=OF, # propellant optimization for cea output
                          p_ch=ENG_P_CH, # propellant optimization for cea output
                          T_ch=ENG_T_CH, # JUST FOR ROCKET
                          ke=ENG_KE, # JUST FOR ROCKET
                          mm=ENG_MM, # JUST FOR ROCKET
                          perturbations=[0, 0, 0, 0, True, 0, 0, 0, 0, 0, 0, True],
                          dt=0.025, # could be DT 
                          adaptive=True, 
                          tol=0.005, 
                          descend=True, 
                          early_return=False, 
                          recovery=False)

    # get trajectory info from optimal design
    #sim = trajectory(False, 0, 0, 0, 0, 0, 0, 0, m_prop, mdot, p_e,
    #           THROTTLE_WINDOW, MIN_THROTTLE,
    #           RCS_MDOT, RCS_P_E, RCS_P_CH,
    #           0, FIN_ROOT, FIN_TIP, FIN_SWEEP_ANGLE, FIN_SEMISPAN, FIN_THICKNESS,
    #            LOX_TANK_P, IPA_TANK_P,
    #           AIRFRM_IN_RAD, OF, ENG_P_CH, ENG_T_CH, ENG_KE, ENG_MM,
    #           [0, 0, 0, 0, True, 0, 0, 0, 0, 0, 0, True],
    #                      DT, False, 0.045, False, False)

    print('finished sim')

    print(sim.LV4.launch_speed, "m/s")
    print(sim.thrust[0], "N")


    textlist = print_results(sim, True)
    # draw pretty pictures of optimized trajectory
    # what is our stability???????
    rocket_plot(sim.t, sim.alt, sim.v, sim.a, sim.thrust,
                sim.dyn_press, sim.Ma, sim.stability_margin, sim.m, sim.p_a, sim.drag, sim.throttle, sim.fin_flutter, sim, True, None, None)

    # get/print info about our trajectory and rocket
    for line in textlist:
        print(line)

    print('\nMaking an OpenRocket rocket and corresponding engine!')
    # create an openrocket file with matching engine for our design (and print/save trajectory data)
    make_engine(mdot, sim.m_prop, sim.thrust[0:sim.F_index + 1],
                sim.LV4.inr_r, AIRFRAME_THICKNESS, sim.LV4.l_o, sim.LV4.l_f, sim.LV4.m_tank_o, sim.LV4.m_tank_f,
                sim.t[sim.F_index], sim.LV4.engine.Ve/G_N,
                sim.LV4.eng_sys_dry_mass, sim.LV4.eng_sys_len, sim.openrocket_CoM,
                sim.LV4.ballast, sim.LV4.fin)

