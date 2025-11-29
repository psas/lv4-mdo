

import matplotlib.pyplot as plt

from customize.system_definition import *
import numpy as np
from datetime import datetime


def save_trajectory(sim):
    T_a     = np.array([state[1][1][4][2] for state in sim.raw_states[1:]])
    rel_vel = np.array([state[1][1][5][2] for state in sim.raw_states[1:]])
    alpha   = np.array([state[1][1][5][5] for state in sim.raw_states[1:]])
    arr     = np.array([sim.t, sim.alt, sim.v, sim.a, sim.thrust, sim.m,
                    rel_vel, alpha, sim.Ma, sim.dyn_press, sim.p_a, sim.rho, T_a]).T
    np.savetxt(RKT_PREFIX +'psas_rocket_'+str(get_index())+'_traj.csv', arr, delimiter=',',
                header='time (s),altitude (m),absolute velocity (m/s),acceleration (m/s^2),thrust (N),mass (kg),freestream velocity (m/s),angle of attack (rads),Mach,dynamic pressure (Pa),air pressure (Pa),air density (kg/m^3),air temperature (K)')






# this creates a list of strings for relevant data of trajectory
# arranging in sensible order is low priority. there might be more information to extract also.
# maybe get input from on high about most relevant info and how to organize it nicely
def print_results(sim, save):
    max_tabs = 42
    text_base = []
    text_base.append(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    text_base.append('\nDESIGN VECTOR')
    text_base.append('\n-----------------------------')
    text_base.append('\ndesign total propellant mass               = {:.4f} kg'.format(sim.m_prop[0]))
    text_base.append('\ndesign unadjusted propellant mass          = {:.4f} kg'.format(sim.LV4.unadjusted_mprop))
    text_base.append('\ndesign mass flow rate                      = {:.4f} kg/s'.format(sim.LV4.engine.mdot))
    text_base.append('\ndesign nozzle exit pressure                = {:.4f} Pa'.format(sim.LV4.engine.p_e))
    text_base.append('\ntotal tankage length (after adjustment)    = {:.4f} m'.format(sim.LV4.l_o + sim.LV4.l_f))
    text_base.append('\ndesign airframe diameter                   = {:.4f} m.'.format(sim.LV4.diameter))
    text_base.append('\ndesign airframe total length               = {:.4f} m.'.format(sim.LV4.length))
    text_base.append('\ndesign GLOW                                = {:.4f} kg'.format(sim.LV4.GLOW))
    text_base.append('\ndesign ballast mass                        = {:.4f} kg'.format(sim.LV4.ballast))
    text_base.append('\nconical part of nosecone length            = {:.4f} m'.format(sim.LV4.nose_l-CYL_NOSE_L))
    text_base.append('\ndesign fin root chord                      = {:.4f} m'.format(sim.LV4.fin.root))
    text_base.append('\ndesign fin tip chord                       = {:.4f} m'.format(sim.LV4.fin.tip))
    text_base.append('\ndesign fin sweep angle                     = {:.4f} deg'.format(np.degrees(sim.LV4.fin.sweep_angle)))
    text_base.append('\ndesign fin span                            = {:.4f} m'.format(sim.LV4.fin.semispan))
    text_base.append('\ndesign fin thickness                       = {:.4f} mm'.format(sim.LV4.fin.thickness*1000))

    text_base.append('\n')
    text_base.append('\nCONSTRAINTS')
    text_base.append('\n-----------------------------')
    text_base.append('\nLargest angle of attack (c.f. < {})         = {:.4f} deg'.format(
                                                                            CONS_AOA, sim.tip_off_aoa))
    text_base.append('\nL/D ratio (c.f. < {})                       = {:.4f}'.format(
                                                                            CONS_LD, sim.ld_ratio))
    text_base.append('\nfin flutter ratio (c.f. > {})                       = {:.4f}'.format(
                                                                            1.0, sim.min_fin_flutter))
    text_base.append('\nSommerfield criterion (c.f. pe/pa >= {})    = {:.4f}'.format(CONS_S_CRIT, sim.S_crit))
    text_base.append("\nmax acceleration (c.f. < {})                = {:.4f} gs".format(
                                                                                CONS_ACCEL, sim.max_g_force))
    text_base.append('\nTWR at lift off (c.f. > {})                 = {:.4f}'.format(CONS_TWR, np.linalg.norm(sim.thrust[0]) / (sim.LV4.GLOW * 9.807)))
    text_base.append('\nLowest stability margin caliber (c.f. > {}) = {:.4f}'.format(CONS_STBLTY, sim.min_stability))
    text_base.append('\nspeed when leaving launch rail (c.f. > {})  = {:.4f} m/s'.format(CONS_LS, sim.launch_speed))
    text_base.append('\naltitude at apogee (c.f. > {})              = {:.4f} km'.format(
                                                                                CONS_ALT/1000, sim.apogee/1000))
    text_base.append('\nLFETS LOX velocity (c.f. < {})    = {:.4f} m/s'.format(CONS_V_LFETS, sim.LV4.v_lfets_o))
    text_base.append('\nLFETS IPA velocity (c.f. < {})    = {:.4f} m/s'.format(CONS_V_LFETS, sim.LV4.v_lfets_f))
    text_base.append('\ndesign thrust (ground level) (c.f. < {})    = {:.4f} N'.format(CONS_THRUST,
                                                                                        np.linalg.norm(sim.thrust[0])))

    text_base.append('\n')
    text_base.append('\nADDITIONAL INFORMATION')
    text_base.append('\n-----------------------------')
    text_base.append("\nstarted at (lat, long, height): " + str(sim.env.ECEF_to_geodetic(sim.raw_states[0][0][1])))
    text_base.append("\nended at (lat, long, height): " + str(sim.env.ECEF_to_geodetic(sim.raw_states[-1][0][1])))
    text_base.append('\ndesign thrust (vacuum)                     = {:.2f} kN'.format(np.linalg.norm(sim.thrust[sim.F_index])/1000))
    text_base.append('\ndesign total dry mass                      = {:.4f} kg'.format(sim.LV4.GLOW - sim.m_prop[0]))
    text_base.append('\nmission time at landing                    = {:.4f} s'.format(sim.t[-1]))
    text_base.append('\nmission time at apogee                     = {:.4f} s'.format(sim.t[sim.ap_index]))
    text_base.append('\nmission time at burnout                    = {:.4f} s'.format(sim.t[sim.F_index]))
    text_base.append('\nmax dynamic pressure                       = {:.4f} Pa'.format(sim.maxq))
    text_base.append('\ndesign dV                                  = {:.4f} km/s'.format(sim.dV1))
    text_base.append('\nestimated minimum required dV              = {:.4f} km/s'.format(sqrt(2*G_N*sim.alt[sim.ap_index])/1000))

    # Creates a line with entry text with unites(e), value (v), and tab size (s)
    prep = lambda e, v, s=42: ('\n{}:\t= {:.4f}'.format(e, v)).expandtabs(s)

    text_base.append("\n")
    text_base.append("\nRCS SYSTEM DETAILS")
    text_base.append("\n-----------------------------")
    text_base.append(prep('RCS impulse budget (N s)', sim.rcs_impulse_budget))
    text_base.append(prep('avg RCS nozzle thrust (N)', sim.rcs_design_thrust))
    text_base.append(prep('N2 tank init pressure (Pa)', sim.LV4.rcs_tank.p_0))
    text_base.append(prep('N2 tank radius (m)', sim.LV4.rcs_tank.tank_r))
    text_base.append(prep('N2 tank length (m)', sim.LV4.rcs_tank.tank_l))
    text_base.append(prep('N2 actual mass (kg)', sim.LV4.rcs_tank.gas_mass))
    text_base.append(prep('N2 actual tank volume (m^3)', sim.LV4.rcs_tank.tank_V))
    text_base.append(prep('N2 desired mass (kg)', sim.design_n2_mass))
    text_base.append(prep('N2 desired tank volume (m^3)', sim.design_n2_tank_volume))


    text_base.append("\n")
    text_base.append("\nISOGRID TANK DETAILS")
    text_base.append("\n-----------------------------")
    text_base.append(prep("Isogrid rib thickness (cm)", sim.LV4.rib_t*100))
    text_base.append(prep("Isogrid # of radial cell divisions", sim.LV4.num_radl_dvsns))
    text_base.append(prep("Ox mass (kg)", sim.LV4.m_o))
    text_base.append(prep("Fuel mass (kg)", sim.LV4.m_f))
    text_base.append(prep("Ox tank length + ullage (m)", sim.LV4.l_o))
    text_base.append(prep("Fuel tank length + ullage (m)", sim.LV4.l_f))
    text_base.append(prep("Ox tank pressure (Pa)", sim.LV4.lox_tank.p_0))
    text_base.append(prep("Fuel tank pressure (Pa)", sim.LV4.ipa_tank.p_0))
    text_base.append(prep("Max tank pressure (Pa)", sim.LV4.TANK_MAX_P))


    text_base.append("\n")
    text_base.append("\nENGINE SYSTEM DETAILS")
    text_base.append("\n-----------------------------")
    mdot_o, mdot_f = proportion(sim.LV4.engine.mdot, sim.LV4.OF)
    text_base.append("\nOx flow:                                   = {:.4f} kg/s".format(mdot_o))
    text_base.append("\nFuel flow:                                 = {:.4f} kg/s".format(mdot_f))
    text_base.append("\nOF Ratio:                                  = {:.4f}".format(sim.LV4.OF))
    text_base.append('\ndesign chamber pressure                    = {:.4f} Pa'.format(sim.LV4.engine.p_ch))
    text_base.append('\nest. regen pressure drop                   = {:.4f} Pa'.format(sim.LV4.delp_regen))
    text_base.append('\ndesign expansion ratio                     = {:.4f}'.format(sim.LV4.engine.ex))
    text_base.append('\ndesign Exit area                           = {:.4f} cm^2'.format(sim.LV4.engine.A_e*10e3))
    text_base.append('\ndesign throat area                         = {:.4f} cm^2'.format(sim.LV4.engine.A_t*10e3))
    text_base.append('\ndesign Throat pressure                     = {:.4f} Pa'.format(sim.LV4.engine.p_t))
    text_base.append('\ndesign Throat temperature                  = {:.4f} K'.format(sim.LV4.engine.T_t))
    text_base.append('\ndesign Chamber temperature                 = {:.4f} K'.format(sim.LV4.engine.T_ch))
    text_base.append('\ndesign exit velocity                       = {:.4f} m/s'.format(sim.LV4.engine.Ve))
    text_base.append('\ndesign isp                                 = {:.4f} s'.format(sim.LV4.engine.Ve/G_N))
    text_base.append('\ndesign average impulse                     = {:.4f} N*s'.format(sim.impulse))
    text_base.append("\n")

    text_base.append("\nEFS SYSTEM DETAILS")
    text_base.append("\n-----------------------------")
    text_base.append("\nOutlet Pressure req. (Ox): " + str(sim.LV4.p_out_o) + " Pa")
    text_base.append("\nPower requirement (Ox): " + str(sim.LV4.pow_o) + " W")
    text_base.append("\nRotational Speed (Ox): " + str(sim.LV4.rpm_o) + " RPM")
    text_base.append("\nOutlet Pressure req. (Fuel): " + str(sim.LV4.p_out_f) + " Pa")
    text_base.append("\nPower requirement (Fuel): " + str(sim.LV4.pow_f) + " W")
    text_base.append("\nRotational Speed (Fuel): " + str(sim.LV4.rpm_f) + " RPM")

    sim.LV4.read_out()
    text_base.append('\n')
    text_base.append('\nPOST-FLIGHT MASS BUDGET')
    text_base.append('\n-----------------------------\n')
    for line in sim.LV4.description:
        text_base.append(line)

    def print_structural_analysis(title, loads):
        """Makes three plots: axial, lateral, bending moment

        title = the title of the plot
        loads = an array that contains:
                    - the axial load at the top of each component
                    - the laterial load at the top of each component
                    - the bending moment at each location
        """
        result_base = []
        axial_load, lateral_load, bending_moments = loads
        result_base.append('\n----\n'+title+'\n----')
        result_base.append('\nModule\t\tAxial \tLateral \tBending moments'.expandtabs(16))
        for i, element in enumerate(reversed(sim.LV4.parts)):
            loads=('%f\t%f\t%f'%(axial_load[i], lateral_load[i], bending_moments[i])).expandtabs(16)
            line=('\n  %s\t'%element.name).expandtabs(32) + loads
            result_base.append(line)
        return result_base

    text_base.append('\n\n\nSTRUCTURAL ANALYSIS')
    text_base.append('\n-----------------------------\n')
    text_base.append('\nAxial and Lateral loads occur at the top of the module\n')
    elevation_top = np.array([part.coords[2] + part.length for part in reversed(sim.LV4.parts)])

    text_base += print_structural_analysis('Loads at Launch', sim.LV4.loads_at_launch)
    text_base += print_structural_analysis('Loads at Max-Q', sim.LV4.loads_at_max_q)
    text_base += print_structural_analysis('Loads before Burnout', sim.LV4.loads_before_burnout)
    text_base += print_structural_analysis('Loads after Burnout', sim.LV4.loads_after_burnout)


    text_base.append('\n\n')
    text_base.append(sim.LV4.cea_output)

    if save:
        # create a file with all this info in it
        with open(RKT_PREFIX + 'psas_rocket_' + str(get_index()) + '_traj.txt', 'w') as traj:
            for line in text_base:
                traj.write(line)

    return text_base




# this creates a nice set of plots of our trajectory data and saves it to rocket_farm
def rocket_plot(t, alt, v, a, F, q, Ma, stability, m, p_a, D, throttle, fin_v, sim, save, start, end):
    # pylab.rcParams['figure.figsize'] = (20.0, 20.0)

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111, projection='3d')
    xyz = [state[0][1] for state in sim.raw_states[start:end]]
    x, y, z = [state[0] for state in xyz], [state[1] for state in xyz], [state[2] for state in xyz]
    ax.plot(x, y, z)
    ax.grid()
    plt.title('Trajectory')
    plt.show()

    #fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10) = plt.subplots(10, sharex=True)
    #fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax8, ax9) = plt.subplots(8, sharex=True)
    fig, (ax1, ax2, ax3, ax5, ax6, ax7, ax8) = plt.subplots(7, sharex=True)
    #for n in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10):
        #n.spines['top'].set_visible(False)
        #n.spines['right'].set_visible(False)
        #n.yaxis.set_ticks_position('left')
        #n.xaxis.set_ticks_position('bottom')
        #n.yaxis.labelpad = 40

    ax1.plot(t[start:end], alt[start:end]/1000, 'k')
    ax1.set_ylabel("Altitude (km)")
    ax1.yaxis.major.locator.set_params(nbins=10)
    ax1.set_title('LV4 Trajectory')

    ax2.plot(t[start:end], v[start:end], 'k')
    ax2.yaxis.major.locator.set_params(nbins=10)
    ax2.set_ylabel("Velocity (m/s)")

    ax3.plot(t[start:end], a[start:end]/G_N, 'k')
    ax3.yaxis.major.locator.set_params(nbins=10)
    ax3.set_ylabel("Acceleration/g_n")

    #ax4.plot(t[start:end], F[start:end]/1000, 'k')
    #ax4.yaxis.major.locator.set_params(nbins=10)
    #ax4.set_ylabel("Thrust (kN)")

    ax5.plot(t[start:end], q[start:end]/1000, 'k')
    ax5.yaxis.major.locator.set_params(nbins=10)
    ax5.set_ylabel("Dynamic Pressure (kPa)")

    ax6.plot(t[start:end], Ma[start:end], 'k')
    ax6.yaxis.major.locator.set_params(nbins=10)
    ax6.set_ylabel("Mach number")
    ax6.set_xlabel("t (s)")

    ax7.plot(t[start:end], stability[start:end], 'k')
    ax7.yaxis.major.locator.set_params(nbins=10)
    ax7.set_ylabel("Calibers")
    ax7.set_xlabel("t (s)")

    angle_attack = []
    for state in sim.raw_states[sim.LV4.tower_index:sim.LV4.F_index]:
        angle_attack.append(state[1][1][5][5]*180/np.pi)

    ax8.plot(t[sim.LV4.tower_index:sim.LV4.F_index], angle_attack)
    ax8.yaxis.major.locator.set_params(nbins=10)
    ax8.set_ylabel("Angle of Attack")
    ax8.set_xlabel("t (s)")

    #ax7.plot(t, np.array(m)*0.666*np.array(a), 'k')
    #ax7.plot(t[start:end], fin_v[start:end], 'k')
    #ax7.yaxis.major.locator.set_params(nbins=8)
    #ax7.set_ylabel("Fin Flutter Ratio")
    #ax7.set_ylabel("LOX Tank Axial Load")
    #ax7.set_xlabel("t (s)")

    #ax8.plot(t[start:end], D[start:end], 'k')
    #ax8.yaxis.major.locator.set_params(nbins=10)
    #ax8.set_ylabel("Drag (N)")

    #ax9.plot(t[start:end], p_a[start:end]/1000, 'k')
    #ax9.yaxis.major.locator.set_params(nbins=8)
    #ax9.set_ylabel("Air Pressure (kPa)")

    #ax10.plot(t[start:end], throttle[start:end], 'k')
    #ax10.yaxis.major.locator.set_params(nbins=5)
    #ax10.set_ylabel("Throttle (%)")

    # we save the nice figures we make and then display them
    if save:
        plt.savefig(RKT_PREFIX + 'psas_rocket_' + str(get_index()) + '_traj.svg')
    plt.show()

def Init_Stability(sim):
    off_rail = sim.t[sim.tower_index]
    plt.plot(sim.t[0:200], sim.stability_margin[0:200], 'k')
    plt.axvline(x=off_rail)
    plt.ylabel("Calibers")
    plt.xlabel("t (s)")
    plt.show()
    return None




# this function plots the data from the basic structural analysis. at some point, would be nice to put an image of the rocket in too for reference.
def structural_plot(rkt):
    elevation_top = np.array([part.coords[2] + part.length for part in reversed(rkt.parts)])

    def make_structural_plot(title, loads):
        """Makes three plots: axial, lateral, bending moment

        title = the title of the plot
        loads = an array that contains:
                    - the axial load at the top of each component
                    - the laterial load at the top of each component
                    - the bending moment at each location
        """
        axial_load_at_top, lateral_load_at_top, bending_moments = loads
        fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
        ax1.set_title(title)
        ax1.plot(elevation_top, axial_load_at_top, 'k')
        ax1.set_ylabel("Axial Load (N)")
        ax1.yaxis.major.locator.set_params(nbins=12)
        ax2.plot(elevation_top, lateral_load_at_top, 'k')
        ax2.set_ylabel("Lateral Load (N)")
        ax2.yaxis.major.locator.set_params(nbins=12)
        ax3.plot(elevation_top, bending_moments, 'k')
        ax3.set_ylabel("Bending Moment (N m)")
        ax3.yaxis.major.locator.set_params(nbins=12)
        ax3.set_xlabel("Elevation, Top (m)")
        ax3.xaxis.major.locator.set_params(nbins=16)
        plt.show()


    make_structural_plot('Loads at Tip-Off', rkt.loads_at_launch)
    make_structural_plot('Loads at Max-Q', rkt.loads_at_max_q)
    # make_structural_plot('Loads past Max-Q', rkt.loads_past_max_q)
    make_structural_plot('Loads before Burnout', rkt.loads_before_burnout)
    make_structural_plot('Loads after Burnout', rkt.loads_after_burnout)





# this creates some plots of the phase spaces of all our designs, doesn't save them
def phase_plot(m_prop, mdot_0, p_e):
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

    ax1.plot(m_prop, p_e)
    ax1.set_title('Design Vectors')
    ax1.yaxis.major.locator.set_params(nbins=6)
    ax1.set_xlabel("Propellant (kg)")
    ax1.set_ylabel("Exit pressure (kPa)")

    ax2.plot(mdot_0, p_e)
    ax2.yaxis.major.locator.set_params(nbins=6)
    ax2.set_xlabel("Mass flow rate (kg/s)")
    ax2.set_ylabel("Exit pressure (kPa)")

    ax3.plot(m_prop, mdot_0)
    ax3.yaxis.major.locator.set_params(nbins=6)
    ax3.set_ylabel("Mass flow rate (kg/s)")
    ax3.set_xlabel("Propellant (kg)")

    # we display the first diagram of projected 2d phase portraits
    plt.show()

    fig2 = plt.figure()
    ax = fig2.add_subplot(111, projection='3d')
    ax.set_title('Design Space Trajectory')
    ax.plot(m_prop, mdot_0, p_e)
    ax.set_xlabel("Propellant (kg)")
    ax.set_ylabel("Mass flow rate (kg/s)")
    ax.set_zlabel("Exit pressure (kPa)")

    # we display the interactive 3d phase portrait
    plt.show()
    # note, we're choosing to not automatically save these, but they can be saved from the interface
    fig3 = plt.figure()
    ax = fig3.add_subplot(111)
    ax.set_title('Design Values as % of optimum design value')
    ax.plot(m_prop/(m_prop[-1]), "r", label="m_prop")
    ax.plot(mdot_0/(mdot_0[-1]), "g", label="mdot_0")
    ax.plot(p_e/(p_e[-1]), "b", label="p_e")
    ax.legend()
    plt.show()

# this function makes some nice plots of the phase space of design vectors
def design_grapher(allvectors): # no arguments because I'm using dirty dirty global variables
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Objective + Constraint convergence')
    ax.plot(allobjfun, "r", label="Objective Function Evaluations (Logarithmic Scale)")
    ax.legend()
    plt.show()

    if len(allvectors[0]) == 2:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Design space')
        ax.plot([v[0] for v in allvectors], [v[1] for v in allvectors])
        ax.legend()
        plt.show()
    else:
        designplot = [[],[],[]]
        for i in range(0, len(allvectors)):
            for j in range(0, len(designplot)):
                designplot[j].append(allvectors[i][j])
        phase_plot(designplot[0], designplot[1], designplot[2])

