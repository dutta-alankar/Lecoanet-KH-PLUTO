#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENO3
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            11

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             EXPLICIT
#define  VISCOSITY                      EXPLICIT
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  DEL_RHO_BY_RHO0                0
#define  REYNOLDS                       1
#define  AMP                            2
#define  Y1                             3
#define  Y2                             4
#define  SIGMA                          5
#define  TANH_A                         6
#define  U_FLOW                         7
#define  RHO0                           8
#define  PRS0                           9
#define  LENGTH                         10

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH                    (g_inputParam[LENGTH])
#define  UNIT_DENSITY                   (g_inputParam[RHO0])
#define  UNIT_VELOCITY                  (g_inputParam[U_FLOW])
#define  MULTIPLE_LOG_FILES             YES

/* [End] user-defined constants (do not change this line) */
