solver
{ lname = 6
}



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
 dt_face     real /1e-15/  #>@var current solver time step
 dt_face_old      real /1e-15/ #>@var previous solver time step
 dt_face_last     real /1e-15/ #>@var last solver time step with successful convergence
 min_dt_face /1e-15/   real  #>@var current solver time step
 max_dt_face /1e99/   real  #>@var current solver time step
 dt0_face  real  # nominal time step
 reduction_factor_dt_spc real  /1.0/ #
 reduction_factor_dt_heat real  /1.0/ #
 adjust_reduction_factor logical  /FALSE/ #
 #adjust_reduction_factor_string character*256  #
Nstep_increase_dt integer /10/ #
end_time   real  # end time of simulations
time   real  # current time of simulations
time_savevol   real  #
time_savetime   real  #
start_time  real [s] # start time  of simulation
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
left_surface_model_int(1:nspc) _integer  #
right_surface_model_int(1:nspc) _integer  #
  surf_model_B integer /999/ #
  surf_model_N integer /998/ #
  surf_model_S integer /997/ #
  
**** ParticleFlux:
j_implantation_depth(1:nspc) _integer  #
j_diagnostic_depth(1:nspc) _integer  #






***** MatTemp:
 temp(1:ndt,0:ngrd)  _real /300.0/ [K] #


***** ReactionRatesHeader:
  nuth (1:ndt,0:ngrd,1:nspc,1:nspc) _real  #
  kbin (1:ndt,0:ngrd,1:nspc,1:nspc,1:nspc) _real  #


***** Bulk:
  dens(ndt,0:ngrd,nspc) _real [m^-3] # density of particles
  srs(1:ndt,0:ngrd,1:nspc) _real  #
  src_profile(0:ngrd,nspc) _real  #
  srb (ndt,0:ngrd,nspc,nspc) _real  #
  src (ndt,0:ngrd,nspc) _real  #
  cdif(ndt,0:ngrd,nspc) _real [m^2s^-1]  # coefficient of diffusion
  rct (ndt,0:ngrd,nspc) _real  # 
  ero_flx (ndt,0:ngrd,nspc)   _real [m^2s^-1] # erosion flux
  dif_flx (ndt,0:ngrd,nspc)   _real [m^2s^-1] # diffusion flux
  flx (ndt,0:ngrd,nspc) _real [m^2s^-1]  # total flux of particles
  rate_d (ndt,0:ngrd,nspc) _real  #
  

  nu (1:nspc) _real [s^-1] # 
  j0 (1:nspc) _real  #
   jout(1:ndt,1:nspc) _real  #
    rate_t(ndt,0:ngrd) _real  #
    qflx(ndt,0:ngrd) _real  # heat flux
    ero_qflx(ndt,0:ngrd) _real  #
     min_rate_surface real /1.0e-20/ #
  




***** test: # added by J.Guterl
print_hello() subroutine
initialize_wrapper() subroutine
initializec() subroutine
allocate() subroutine
do_step(status_step:logical) subroutine
update_dt() subroutine