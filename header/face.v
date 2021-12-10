face
{ lname = 6
}
***** Dims:
nspc integer +dim +params /1/ #
ngrd integer +dim +params /100/ # number of grid points
ndt  integer +dim +params  /1/ # var size of storage for time dependent variables

***** SpeciesHeader:
namespc(1:nspc) _character*6 /'D'/ +params #

***** GridHeader:
length  real  [m]  +params /0.01/ # wall_thickness
alpha real  +params  /1.15305056/ # Cell width scaling factor
grid_dx0 real  [m]  +params /1e-9/ # grid first cell length
grid_type character*256  +params /"A"/  #A: antisym (smallest cell at x=0) S: symmetric
grid_gen_mode character*256  +params /"alpha"/  #grid generation: [alpha]alpha=cell_scaling_factor|[seed]dx0=grid_dx0


***** SurfaceHeader:
left_surface_model(1:nspc) _character*1 +params /'S'/  #'B: Gamamaout=Kdes*cb^2 S: Gammaout=Kcs^2 N: no flux
right_surface_model(1:nspc) _character*1 +params /'S'/ #B: Gamamaout=Kdes*cb^2 S: Gammaout=Kcs^2 N: no flux

order_desorption_left(1:nspc) _real  +params /2.0/ # order of desorption cs^order or cb ^order
order_desorption_right(1:nspc) _real  +params /2.0/ # order of desorption cs^order or cb ^order
Eabs_l(1:nspc) _real [eV]  +params /0.1/ #Energy of absorption (vacuum->surface)
Edes_l(1:nspc) _real [eV]  +params /1.4/ # Energy of desorption (surface->vacuum)
Edes_lsat(1:nspc) _real [eV]  +params /0.7/ # Energy of desorption when surface is saturated (surface->vacuum)
Eb_l(1:nspc) _real [eV]  +params /2.0/ # Energy of bulk absortion (surface->bulk)
Eads_l(1:nspc) _real [eV]  +params /1.0/# Energy of adsorption (bulk->surface)

Eabs_r(1:nspc) _real [eV]  +params /0.1/ # Energy of absorption (vacuum->surface)
Edes_r(1:nspc) _real [eV]  +params /1.4/ # Energy of desorption (surface->vacuum)
Edes_rsat(1:nspc) _real [eV]  +params /0.7/ # Energy of desorption when surface is saturated (surface->vacuum)
Eb_r(1:nspc) _real [eV]  +params /2.0/ # Energy of bulk absortion (surface->bulk)
Eads_r(1:nspc) _real [eV]  +params /1.0/ # Energy of adsorption (bulk->surface)

nu (1:nspc) _real +params [s^-1] /1e13/ #Debye frequency for surface reaction rates  

dsrfl0(1:nspc) _real [m^-2]  +params /1e10/ #Initial left surface density of species
dsrfr0(1:nspc) _real [m^-2]  +params /1e10/ #Initial right surface density of species
dsrfm(1:nspc) _real [m^-2]  +params /1e19/ # Maximum surface density of species


***** VerboseHeader:
verbose_parser logical /FALSE/ #
verbose_input logical /FALSE/ #
verbose_init logical /FALSE/ #
verbose_step logical  /FALSE/ #
verbose_cap logical  /FALSE/ #
verbose_debug logical  /FALSE/ #
verbose_couple logical  /FALSE/ #
verbose_restore logical  /FALSE/ #
verbose_maths logical  /FALSE/ #
verbose_interface logical  /FALSE/ #
verbose_help logical  /FALSE/ #
verbose_version logical  /FALSE/ #
verbose_header logical  /FALSE/ #
verbose_surface logical  /FALSE/ #
verbose_compute logical  /FALSE/ #
verbose_update_dt logical  /FALSE/ #

***** BulkHeader:

  dens0(1:nspc) _real [m^-3]  +params /1e10/ # initial density of particles
  gxmax(1:nspc) _real [m]  +params /1.0/# 
  gsigm(1:nspc) _real [m]  +params /1.0/#
  gprof(1:nspc) _character*1  # Initial density profile is Gaussian (G), step(S), linear(L),peak(P),flat(F)'
  cdif0(1:nspc) _real [m^2s^-1]  +params /1e-7/ # pre-exponential factor of diffusion coefficients
  edif (1:nspc) _real [eV]  +params /0.4/ # activation energy of diffusion
   etr   (nspc) _real [eV]  +params /10.0/ # trapping energy 
   edtr  (nspc) _real [eV]  +params /0.0/ # detrapping energy
   densm(1:nspc) _real [m^-3]  +params /1e29/ # max denisty of particles
     nuth0(1:nspc,1:nspc) _real  +params /0.0/ #
     kbin0(1:nspc,1:nspc,1:nspc) _real  +params /0.0/ #
     eth  (1:nspc,1:nspc) _real  +params  /1000.0/ #
     ebin (1:nspc,1:nspc,1:nspc) _real  +params  /1000.0/  #
 cero      real  +params /0/ [m.s^-1] # erosion velocity
 cero_min real  +params /0/#
 cero_max real  +params /0/#
 gamero real  +params  /0/ #
 
***** MatTempHeader:
  solve_heat_eq logical /FALSE/ #
  
***** HeatFluxHeader:
  qflx_in  real  +params  /0 /# incoming heat flux from plasma
  rad      real  +params  /0/ #
 rad_min   real  +params  /0/ #
  rad_max  real  +params  /0/ #
 # t1        real  +params  /0/ #
 #  t2       real  +params  /0/ #
 # t3        real  +params  /0/ #
 # tpulse    #      real  +params  /0 /# period of plasma pulse
 
**** ParticleFluxHeader:
implantation_model(1:nspc) _character*(1) /'S'/ # G: gaussian S:Step E: ERFC
implantation_depth(1:nspc) _real [m]  +params /5e-9/ # implentation  depth
diagnostic_depth(1:nspc) _real  [m]  +params /5e-9/ # diagnostic depth
implantation_width(1:nspc) _real [m]  +params /5e-9/ # implentation width
enrg(1:nspc) _real [eV]  +params /0.0/ # Impact energy of ionized species
inflx(1:nspc)  _real  [m^-2 s^- 1] +params /1.e20/ # influx of particles (may differ from nominal particle flux in pulsed_plasma mode)
j_implantation_depth(1:nspc) _integer  #
j_diagnostic_depth(1:nspc) _integer  #
 gas_pressure(1:nspc) _real +params #
  gas_temp(1:nspc) _real +params #
  mass(1:nspc) _real+params  #
  
***** MaterialHeader:
lambda real  +params  /1e-10/ #
cvlm real  +params  /1.0/ #
csrf real  +params  /1.0/ #
clng real  +params  /1.0/ #
thcond   real  +params  /0/ #
rho      real  +params  /0/ #
cp       real  +params  /0/ #
rhocp     real  +params  /0/ #
emiss    real  +params  /0/ #

***** PhysConstsHeader: # Some physical and mathematical constants
ee real /1.602176462e-19/ #
eps0 real /8.854187817e-12/ #
amass real /1.66053886e-27/ #
pi real /3.14159265358979/ #
twopi real /6.28318530717959/ #
sqrt2 real /1.4142135623730951/ #
kb real /1.3806504e-23/ #
eekb real /1.160450595e+04/ #
sigma_sb real /5.670400e-08/ #
eps_machine real /1e-10/

***** TimeHeader:
 dt_face +input     real /1e-15/  #>@var current solver time step
 #dt_face_old      real /1e-15/ #>@var previous solver time step
 #dt_face_last     real /1e-15/ #>@var last solver time step with successful convergence
 #min_dt_face /1e-15/   real  #>@var current solver time step
 #max_dt_face /1e99/   real  #>@var current solver time step
 reduction_factor_dt_spc real  /0.1/ #
 reduction_factor_dt_heat real  /0.1/ #
 adjust_reduction_factor logical  /FALSE/ #
 #adjust_reduction_factor_string character*256  #
Nstep_increase_dt integer /10/ #
#end_time   real [s] # end time of simulations
time   real [s] # current time of simulations
#time_savevol   real [s] #
#time_savetime   real [s] #
#start_time  real [s] # start time  of simulation
cdt      real  # factor for solver time step (dt)
variable_timestep logical  #
#variable_timestep_string character*256  #


***** SolverHeader:
solver_eps real  /3.e-3/ #
solver_udspl real /9.e-1/#
solver_fdspl real /9.e0/#
solver_gdspl real /1.e-3/#
jac_eps real /1e-8/#
 solver_fstp real  /1.e-1/ #
 iter_solver_max integer  /150/ #
 iter_solver_max_first integer  /150/ #
 finalcheck logical  /TRUE/ #
 solver_step_count integer  /0/ #
 max_iter real  /1e9/ #
 nucut real /1d99/ #
 delta real /0.0/ #
 iter_solver integer /0 / # #of solver iterations at each time step
iteration integer /0/ #
 order_solver integer /1/ #Numerical order of solver: 1|2|5
 neq integer  #
 is_initialized logical /FALSE/#
enforce_error logical /TRUE/ #
compute_spc logical  /TRUE/ #
critical_reduction logical /FALSE/ #
loop_reduction_old logical /FALSE/ #
loop_reduction logical /FALSE/ #
counter_reduction integer  /0/ #
is_allocated logical /FALSE/ #
***** TimerHeader:

  tcpustart real #
  tcpufinish real #
  walltime_start real
  walltime_end real   #
  Nprint_run_info integer  /2 / # print info on current run every Nprint_run_info steps




***** Grid:
x(0:ngrd) _real  #
dx(0:ngrd) _real  #


***** BDFsolver:
# ** coefficients for BDF
#     --- 1st order BDF ---
a11 real /1.0/ #
a12 real /1.0/ #
#     --- 2nd order BDF ---
a21 real / 1.33333333/ #
a22 real /-0.33333333/ #
a23 real / 0.66666666/ #
#     --- 5th order BDF ---
a51 real / 1.92/ #
a52 real /-1.44/ #
a53 real / 0.64/ #
a54 real /-0.12/ #
a55 real / 0.48/ #



***** IOHeader:
max_ifile integer  /10000/ #
sfln_voldata integer /0/ #
sfln_srfdata integer /0/ #
sfln_heatdata integer /0/ #
normf real /0.0/ #
restart_filename character*256  #
read_input_file logical /TRUE/ #
restore_state_temp logical /TRUE/ #
input_filename character*256  #
logfile character*256  #
read_restart_file character*256  /'no'/ #read Restart file: yes|no|filename (yes:default "dsave.rst")'
read_state_file character*256  #
#steady_state_string character*256  #
steady_state logical /FALSE/ #
final_state_file character*256  #
casename character*256 /"run"/#
error_status integer /0/
iout integer /6/ #
   path_folder character*256 /"."/ # top folder where simulations files anf folders are written in
 dat_folder character*256 /"data"/ # top folder where vol,srf and heat data files are written in"

***** CapHeader:
active_cap_bulk logical  /FALSE/ #
active_cap_surface logical  /FALSE/ #




***** Surface:
Kabs_l(1:nspc) _real  #
Kdes_l(1:nspc) _real  #
Kb_l(1:nspc) _real  #
Kads_l(1:nspc) _real  #
Kabs_r(1:nspc) _real  #
Kdes_r(1:nspc) _real  #
Kb_r(1:nspc) _real  #
Kads_r(1:nspc) _real  #

K0abs_l(1:nspc) _real  #
K0des_l(1:nspc) _real  #
K0b_l(1:nspc) _real  #
K0ads_l(1:nspc) _real  #

K0abs_r(1:nspc) _real  #
K0des_r(1:nspc) _real  #
K0b_r(1:nspc) _real  #
K0ads_r(1:nspc) _real  #
# left surface at x(j=0)
dsrfl(0:ndt,1:nspc)    _real [m^-2] # density on left surface
 Gsrf_l(0:ndt,1:nspc)  _real  [m^-2s^-1]# net flux of species onto the left surface
 Gabs_l(0:ndt,1:nspc) _real  [m^-2s^-1]#
 Gdes_l(0:ndt,1:nspc) _real  [m^-2s^-1]#
 Gb_l(0:ndt,1:nspc) _real  [m^-2s^-1]#
 Gads_l(0:ndt,1:nspc) _real  [m^-2s^-1]#

#right surface at x(j=n+1)
 dsrfr(0:ndt,1:nspc)    _real  [m^-2]# density on right surface
 Gsrf_r(0:ndt,1:nspc)  _real  [m^-2s^-1]# net flux of species onto the right surface
 Gabs_r(0:ndt,1:nspc) _real  [m^-2s^-1]#
 Gdes_r(0:ndt,1:nspc) _real  [m^-2s^-1]#
 Gb_r(0:ndt,1:nspc) _real  [m^-2s^-1]#
 Gads_r(0:ndt,1:nspc) _real  [m^-2s^-1]#
left_surface_model_int(1:nspc) _integer  #
right_surface_model_int(1:nspc) _integer  #
  surf_model_B integer /999/ #
  surf_model_N integer /998/ #
  surf_model_S integer /997/ #
  








***** MatTemp:
 temp(0:ndt,0:ngrd)  _real /300.0/ [K] #


***** ReactionRatesHeader:
  nuth (0:ndt,0:ngrd,1:nspc,1:nspc) _real /0.0/  #
  kbin (0:ndt,0:ngrd,1:nspc,1:nspc,1:nspc) _real /0.0/  #


***** Bulk:
  dens(0:ndt,0:ngrd,1:nspc) _real [m^-3] # density of particles
  srs(0:ndt,0:ngrd,1:nspc) _real  #
  src_profile(0:ngrd,1:nspc) _real  #
  srb (0:ndt,0:ngrd,nspc,nspc) _real  #
  src (0:ndt,0:ngrd,nspc) _real  #
  cdif(0:ndt,0:ngrd,nspc) _real [m^2s^-1]  # coefficient of diffusion
  rct (0:ndt,0:ngrd,nspc) _real  # 
  ero_flx (0:ndt,0:ngrd,nspc)   _real [m^2s^-1] # erosion flux
  dif_flx (0:ndt,0:ngrd,nspc)   _real [m^2s^-1] # diffusion flux
  flx (0:ndt,0:ngrd,nspc) _real [m^2s^-1]  # total flux of particles
  rate_d (0:ndt,0:ngrd,nspc) _real  #
  

  j0 (1:nspc) _real /0.0/ #
   jout(0:ndt,1:nspc) _real  #
    rate_t(0:ndt,0:ngrd) _real  #
    qflx(0:ndt,0:ngrd) _real  # heat flux
    ero_qflx(0:ndt,0:ngrd) _real  #
     min_rate_surface real /1.0e-20/ #
  




***** test: # added by J.Guterl
print_hello() subroutine
initialize_wrapper() subroutine
initializec() subroutine
allocate() subroutine
do_step(status_step:logical) subroutine
compute_dt() real function
shift_array subroutine