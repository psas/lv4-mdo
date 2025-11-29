




################################################
### optimization constraints
################################################



CONS_IMPLS   = 889600                    # maximum impulse, N s
CONS_AOA     = 8.                       # maximum angle of attack
CONS_MASS    = 450.                      # GLOW constraint, kg, somewhat arbitrary
CONS_LS      = 22.                       # min launch speed from 60' tower constraint, m/s
CONS_TWR     = 2.                        # TWR constraint
CONS_S_CRIT  = 0.35                      # Critical pressure ratio constraint
CONS_ACCEL   = 15.                       # Max acceleration constraint, g's
CONS_LD      = 25.                       # L/D ratio constraint, slightly arbitrary
CONS_ALT     = 105000. #* 0.5                   # Min altitude constraint, m
CONS_THRUST  = 10000                      # max ground-level thrust, N
CONS_CEILING = 150000.                   # base-11 maximum apogee requirement, m
CONS_STBLTY  = 2.0                       # minimum in flight stability margin caliber
#CONS_EFS     = 11000                      # maximum EFS pump power, W
CONS_EFS     = 0  # EFS is being removed
CONS_V_LFETS = 9.144                    # maximum fluid velocity in test stand
CONS_TANK_MIN = 689476 # Pa, minimum tank pressure (so FLIPS can use regulators)



