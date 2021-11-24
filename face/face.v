face
{ lname = 6;
}
***** SpeciesHeader:
 nspc integer /1/ #
namespc(1:nspc) _character*6 /'D'/ #


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

***** TimeHeader:
 dt_face     real  #>@var current solver time step
 dt_face_old     real  #>@var previous solver time step
 dt_face_last     real  #>@var last solver time step with successful convergence
 min_dt_face    real  #>@var current solver time step
 max_dt_face    real  #>@var current solver time step
 dt0_face  real  # nominal time step
 reduction_factor_dt_spc real  /1.0/ #
 reduction_factor_dt_heat real  /1.0/ #
 adjust_reduction_factor logical  /FALSE/ #
 adjust_reduction_factor_string character*256  #
Nstep_increase_dt integer /10/ #
end_time   real  # end time of simulations
time   real  # current time of simulations
time_savevol   real  #
time_savetime   real  #
start_time  real [s] # start time  of simulation
cdt      real  # factor for solver time step (dt)
variable_timestep logical  #
variable_timestep_string character*256  #


 **** SolverHeader:
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
 status_step logical /TRUE/ #
enforce_error logical /FALSE/ #
compute_spc logical  /TRUE/ #
critical_reduction logical /FALSE/ #
loop_reduction_old logical /FALSE/ #
loop_reduction logical /FALSE/ #
counter_reduction integer  /0/ #

***** TimerHeader:

  tcpustart real #
  tcpufinish real #
  walltime_start real
  walltime_end real   #
  Nprint_run_info integer  /2 / # print info on current run every Nprint_run_info steps

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


***** GridHeader:
ngrd integer /100/ # number of grid points
length  real [m] # wall_thickness
x(0:ngrd) _real  #
dx(0:ngrd) _real  #
alpha real  /1.15305056/ # Cell width scaling factor
grid_dx0 real  [m] /1e-9/ # grid first cell length
grid_type character*256 /"A"/  #A: antisym (smallest cell at x=0) S: symmetric
grid_gen_mode character*256 /"alpha"/  #grid generation: [alpha]alpha=cell_scaling_factor|[seed]dx0=grid_dx0


**** BDFsolverHeader:
ndt  integer   /1/ # var size of storage for time dependent variables

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

***** DataDumpHeader:
 dump_space_dt real /0/ #
 dump_time_dt real /0/ #
 dump_restart_dt real /0/ #
 dump_space logical /FALSE/ #
 dump_time logical /FALSE/ #
 dump_restart logical /FALSE/ #
 dump_vol_append logical /FALSE/ #
 dump_srf_append logical /FALSE/ #
 dump_time_append logical /FALSE/ #
 dump_space_string character*256  #
 dump_time_string character*256  #
 dump_restart_string character*256  #
first_voldump logical /TRUE/ #
dump_vol_append_string character*256  #
dump_srf_append_string character*256  #


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
steady_state_string character*256  #
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
active_cap_surface_string character*256  #
active_cap_bulk_string character*256  #

***** onthefly_inventories:


***** inventories:
inv_Nnetbulk(1:nspc) _real
inv_Nnetsrf(1:nspc) _real
inv_Ntotbulk(1:nspc) _real
inv_Ntotsrf(1:nspc) _real
inv_Enetbulk(1:nspc) _real
inv_Etotbulk(1:nspc) _real


***** particle_balances:
      pb_Nnet(1:nspc) _real
      pb_Ninflux(1:nspc) _real
      pb_Noutflux(1:nspc) _real
      pb_p_net(1:nspc) _real /0.0/
      pb_p_max(1:nspc) _real /0.0/
      pb_f_lost(1:nspc) _real /0.0/
      pb_Noutflux_l(1:nspc) _real  /0.0/
      pb_Noutflux_r(1:nspc) _real /0.0/


***** energy_balances:
      eb_Enet real /0.0/
      eb_Ein real /0.0/
      eb_Eout real /0.0/
      eb_p_net real /0.0/
      eb_p_max real /0.0/
      eb_f_lost real /0.0/
      eb_Eout_l real /0.0/
      eb_Eout_r real /0.0/

***** outgassing_fluxes:
       out_Gdes(1:nspc) _real  /0.0/         #
       out_min_Gdes(1:nspc) _real /0.0/          #
       out_max_Gdes(1:nspc) _real /0.0/          #
       out_ave_Gdes(1:nspc) _real /0.0/          # ave deviation of Gdes over FACE run
       out_sig_Gdes(1:nspc) _real /0.0/          # sdt deviation of Gdes over FACE run
       out_Gpermeation(1:nspc)  _real  /0.0/         # =Gdes_r

***** wall_temperatures:
       wt_srf_temp_l real /0.0/ #
       wt_srf_temp_r real /0.0/          #
       wt_mean_temp real  /0.0/         #
       wt_max_temp real  /0.0/         #
       wt_min_temp real  /0.0/         #

***** InventoryHeader:
print_onthefly_inventory logical  /FALSE/ #
print_onthefly_inventory_string character*256  #
onthefly_int_dens(1:nspc) _real  #
onthefly_int_dsrf(1:nspc) _real  #
onthefly_net_int_dens(1:nspc) _real  #
onthefly_net_int_dsrf(1:nspc) _real  #
onthefly_int_des(1:nspc) _real #
onthefly_int_src(1:nspc) _real
trace_flux_sum_inflx(1:nspc) _real  #
trace_flux_sum_qflx(1:nspc) _real  #
trace_flux_sum_Q_l(1:nspc) _real  #
trace_flux_sum_Q_r(1:nspc) _real  #
trace_flux_sum_Gdes_l(1:nspc) _real  #
trace_flux_sig_Gdes_l(1:nspc) _real  #
trace_flux_min_Gdes_l(1:nspc) _real  #
trace_flux_max_Gdes_l(1:nspc) _real  #
trace_flux_sum_Gdes_r(1:nspc) _real  #
trace_flux_sig_Gdes_r(1:nspc) _real  #
trace_flux_min_Gdes_r(1:nspc) _real  #
trace_flux_max_Gdes_r(1:nspc) _real  #

***** SurfaceHeader:
order_desorption_left(1:nspc) _real /2.0/ # order of desorption cs^order or cb ^order
order_desorption_right(1:nspc) _real /2.0/ # order of desorption cs^order or cb ^order
Eabs_l(1:nspc) _real [eV] /0.1/ #Energy of absorption (vacuum->surface)
Edes_l(1:nspc) _real [eV] /1.4/ # Energy of desorption (surface->vacuum)
Edes_lsat(1:nspc) _real [eV] /0.7/ # Energy of desorption when surface is saturated (surface->vacuum)
Eb_l(1:nspc) _real [eV] /2.0/ # Energy of bulk absortion (surface->bulk)
Eads_l(1:nspc) _real [eV] /1.0/# Energy of adsorption (bulk->surface)

Eabs_r(1:nspc) _real [eV] /0.1/ # Energy of absorption (vacuum->surface)
Edes_r(1:nspc) _real [eV] /1.4/ # Energy of desorption (surface->vacuum)
Edes_rsat(1:nspc) _real [eV] /0.7/ # Energy of desorption when surface is saturated (surface->vacuum)
Eb_r(1:nspc) _real [eV] /2.0/ # Energy of bulk absortion (surface->bulk)
Eads_r(1:nspc) _real [eV] /1.0/ # Energy of adsorption (bulk->surface)

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
dsrfl(1:ndt,1:nspc)    _real [m^-2] # density on left surface
 Gsrf_l(1:ndt,1:nspc)  _real  [m^-2s^-1]# net flux of species onto the left surface
 Gabs_l(1:ndt,1:nspc) _real  [m^-2s^-1]#
 Gdes_l(1:ndt,1:nspc) _real  [m^-2s^-1]#
 Gb_l(1:ndt,1:nspc) _real  [m^-2s^-1]#
 Gads_l(1:ndt,1:nspc) _real  [m^-2s^-1]#

#right surface at x(j=n+1)
 dsrfr(1:ndt,1:nspc)    _real  [m^-2]# density on right surface
 Gsrf_r(1:ndt,1:nspc)  _real  [m^-2s^-1]# net flux of species onto the right surface
 Gabs_r(1:ndt,1:nspc) _real  [m^-2s^-1]#
 Gdes_r(1:ndt,1:nspc) _real  [m^-2s^-1]#
 Gb_r(1:ndt,1:nspc) _real  [m^-2s^-1]#
 Gads_r(1:ndt,1:nspc) _real  [m^-2s^-1]#

 left_surface_model_string(1:nspc) _character*6  #
 right_surface_model_string(1:nspc) _character*6  #
 left_surface_model(1:nspc) _integer  #
 right_surface_model(1:nspc) _integer  #
  surf_model_B integer /999/ #
  surf_model_N integer /998/ #
  surf_model_S integer /997/ #
  dsrfl0(1:nspc) _real [m^-2] /1e19/ #Initial left surface density of species
  dsrfr0(1:nspc) _real [m^-2] /1e19/ #Initial right surface density of species
  dsrfm(1:nspc) _real [m^-2] /1e19/ # Maximum surface density of species

**** ParticleFluxHeader:
implantation_model(1:nspc) _character*(1) /'S'/ # G: gaussian S:Step E: ERFC
implantation_depth(1:nspc) _real [m] /5e-9/ # implentation  depth
diagnostic_depth(1:nspc) _real  [m] /5e-9/ # diagnostic depth
j_implantation_depth(1:nspc) _integer  #
j_diagnostic_depth(1:nspc) _integer  #
implantation_width(1:nspc) _real [m] /5e-9/ # implentation width
enrg(1:nspc) _real [eV] /0.0/ # Impact energy of ionized species
inflx(1:nspc)  _real  [m^-2 s^- 1] /1.e20/ # influx of particles (may differ from nominal particle flux in pulsed_plasma mode)
Gamma_in_max(1:nspc)  _real [m^-2 s^- 1] /0.e20/ # max influx of particles
Gamma_in_base(1:nspc)  _real [m^-2 s^- 1] /1.e20/ # nominal influx of particles
Gamma_in_pulse_period(1:nspc)  _real [s] /1e99/ # max influx of particles
Gamma_in_pulse_duration(1:nspc)  _real [s] /1e99/ # max influx of particles
Gamma_in_pulse_starttime(1:nspc)  _real  # max influx of particles
Gamma_in_pulse_type(1:nspc)  _character*(6) /'N'/ # pulsed partice flux N: no S: sin R: rectangle E: ELM
Gamma_in_pulse(1:nspc)  _integer /998/  # nominal influx of particles
 Gamma_in_pulse_R integer /999/ #
 Gamma_in_pulse_N integer /998/ #
 Gamma_in_pulse_S integer /997/ #
 Gamma_in_pulse_B integer /996/ #
 gas_pressure(1:nspc) _real  #
  gas_temp(1:nspc) _real  #
  mass(1:nspc) _real  #

***** HeatFluxHeader:
T_pulse_type character*6 /'N'/ # nominal influx of particles/"
T_pulse  integer  # nominal influx of particles
 T_pulse_R integer /999/ #
 T_pulse_N integer /998/ #
 T_pulse_S integer /997/ #
 T_pulse_B integer /996/ #
T_pulse_max  real   [K]# max influx of particles
T_pulse_period  real   [s] # max influx of particles
T_pulse_duration  real   [s]# max influx of particles
T_pulse_starttime  real   [s]# max influx of particles
Q_in_base  real [W.m^-2]                      # nominal heat flux
Q_in_max  real   [W.m^-2]                    # Maximal external heat flux
Q_in_pulse_period  real  [s]                    #
Q_in_pulse_duration  real  [s]                    #
Q_in_pulse_starttime  real   [s]                   #
Q_pulse_type  character*6 /'N'/ # pulsed heat flux N: no S: sin R: rectangle E: ELM
Q_in_pulse  integer  #
 Q_in_pulse_R integer /999/ #
 Q_in_pulse_N integer /998/ #
 Q_in_pulse_S integer /997/ #
 Q_in_pulse_B integer /996/ #
 qform     real /0/ #
  qflx_in  real /0 / # incoming heat flux from plasma
  rad      real /0/ #
 rad_min   real /0/ #
  rad_max  real /0/ #
  t1        real /0/ #
   t2       real /0/ #
  t3        real /0/ #
  tpulse    #      real /0 / # period of plasma pulse


***** MaterialHeader:
lambda real /0.0/ #
cvlm real /1.0/ #
csrf real /1.0/ #
clng real /1.0/ #
thcond   real /0/ #
rho      real /0/ #
cp       real /0/ #
rhocp     real /0/ #
emiss    real /0/ #

***** MatTempHeader:
 temp(1:ndt,0:ngrd) [K] _real  #
 nramp integer #
 rtime(1:nramp) _real  #
 rtemp(1:nramp) _real  #
 temp_init real [K] /300.0/ #
 temp_final real [K] /300.0/ #
 dtemp real /0.0/ #
  tramp0 real /0/ #
  tramp1 real /0/ #
  framp_string character*256 /'none'/ #
  framp integer  #
  framp_none integer /999/ #
  framp_readfile integer /998/ #
  solve_heat_eq logical /FALSE/ #
  solve_heat_eq_string character*256  #

***** ReactionRatesHeader:
  nuth (1:ndt,0:ngrd,1:nspc,1:nspc) _real  #
  kbin (1:ndt,0:ngrd,1:nspc,1:nspc,1:nspc) _real  #
  nuth0(1:nspc,1:nspc) _real  #
  kbin0(1:nspc,1:nspc,1:nspc) _real  #
  eth  (1:nspc,1:nspc) _real /1.0/ #
  ebin (1:nspc,1:nspc,1:nspc) _real  #

***** BulkHeader:
  srs(1:ndt,0:ngrd,1:nspc) _real  #
  src_profile(0:ngrd,nspc) _real  #
  srb (ndt,0:ngrd,nspc,nspc) _real  #
  src (ndt,0:ngrd,nspc) _real  #
  cdif(ndt,0:ngrd,nspc) _real  #
  rct (ndt,0:ngrd,nspc) _real  #
  ero_flx (ndt,0:ngrd,nspc)   _real  # erosion flux
  dif_flx (ndt,0:ngrd,nspc)   _real  # erosion flux
  cdif0(nspc) _real  #
  edif (nspc) _real  #
  flx (ndt,0:ngrd,nspc) _real  #
  rate_d (ndt,0:ngrd,nspc) _real  #
  densm(1:nspc) _real  #
  etr   (nspc) _real  #
  edtr  (nspc) _real  #
  nu (1:nspc) _real  #
  j0 (1:nspc) _real  #
  dens(ndt,0:ngrd,nspc) _real  #
  dens0(1:nspc) _real  #
  gxmax(1:nspc) _real  #
  gsigm(1:nspc) _real  #
  gprof(1:nspc) _character*6  #
 jout(1:ndt,1:nspc) _real  #
 cero      real  # erosion velocity
 cero_min real  #
 cero_max real  #
 gamero real  #
 rate_t(ndt,0:ngrd) _real  #
 qflx(ndt,0:ngrd) _real  #
 ero_qflx(ndt,0:ngrd) _real  #
  min_rate_surface real /1.0e-20/ #



***** test: # added by J.Guterl
print_hello() subroutine
initialize_wrapper() subroutine
