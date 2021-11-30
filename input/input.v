input
{ lname = 6

}

***** SpeciesHeader:
nspc integer /1/ #
namespc(1:nspc) _character*6 /'D'/ #

***** GridHeader:
ngrd integer /100/ # number of grid points
length  real [m] /0.01/ # wall_thickness
alpha real  /1.15305056/ # Cell width scaling factor
grid_dx0 real  [m] /1e-9/ # grid first cell length
grid_type character*256 /"A"/  #A: antisym (smallest cell at x=0) S: symmetric
grid_gen_mode character*256 /"alpha"/  #grid generation: [alpha]alpha=cell_scaling_factor|[seed]dx0=grid_dx0

***** BDFsolverHeader:
ndt  integer   /1/ # var size of storage for time dependent variables

***** SurfaceHeader:
left_surface_model(1:nspc) _character*1 /'S'/  #'B: Gamamaout=Kdes*cb^2 S: Gammaout=Kcs^2 N: no flux
right_surface_model(1:nspc) _character*1 /'S'/ #B: Gamamaout=Kdes*cb^2 S: Gammaout=Kcs^2 N: no flux

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




dsrfl0(1:nspc) _real [m^-2] /1e10/ #Initial left surface density of species
dsrfr0(1:nspc) _real [m^-2] /1e10/ #Initial right surface density of species
dsrfm(1:nspc) _real [m^-2] /1e19/ # Maximum surface density of species


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

***** BulkHeader:

  dens0(1:nspc) _real [m^-3] /1e10/ # initial density of particles
  gxmax(1:nspc) _real [m] /1.0/# 
  gsigm(1:nspc) _real [m] /1.0/#
  gprof(1:nspc) _character*1  # Initial density profile is Gaussian (G), step(S), linear(L),peak(P),flat(F)'
  cdif0(1:nspc) _real [m^2s^-1] /1e-7/ # pre-exponential factor of diffusion coefficients
  edif (1:nspc) _real [eV] /0.4/ # activation energy of diffusion
   etr   (nspc) _real [eV] /10.0/ # trapping energy 
   edtr  (nspc) _real [eV] /0.0/ # detrapping energy
   densm(1:nspc) _real [m^-3] /1e29/ # max denisty of particles
     nuth0(1:nspc,1:nspc) _real  /0.0/ #
     kbin0(1:nspc,1:nspc,1:nspc) _real  /0.0/ #
     eth  (1:nspc,1:nspc) _real /1000.0/ #
     ebin (1:nspc,1:nspc,1:nspc) _real /1000.0/  #
 cero      real  /0/ [m.s^-1] # erosion velocity
 cero_min real  /0/#
 cero_max real  /0/#
 gamero real /0/ #
 
***** MatTempHeader:
  solve_heat_eq logical /FALSE/ #
  
***** HeatFluxHeader:
#T_pulse character*6 /'N'/ # nominal influx of particles/"
#T_pulse_int  integer  # nominal influx of particles
# T_pulse_R integer /999/ #
# T_pulse_N integer /998/ #
# T_pulse_S integer /997/ #
# T_pulse_B integer /996/ #
#T_pulse_max  real   [K]# max influx of particles
#T_pulse_period  real   [s] # max influx of particles
#T_pulse_duration  real   [s]# max influx of particles
#T_pulse_starttime  real   [s]# max influx of particles
#Q_in_base  real [W.m^-2]                      # nominal heat flux
#Q_in_max  real   [W.m^-2]                    # Maximal external heat flux
#Q_in_pulse_period  real  [s]                    #
#Q_in_pulse_duration  real  [s]                    #
#Q_in_pulse_starttime  real   [s]                   #
#Q_in_pulse  character*6 /'N'/ # pulsed heat flux N: no S: sin R: rectangle E: ELM
#Q_in_pulse_int  integer  #
# Q_in_pulse_R integer /999/ #
# Q_in_pulse_N integer /998/ #
# Q_in_pulse_S integer /997/ #
# Q_in_pulse_B integer /996/ #
# qform     real /0/ #
  qflx_in  real /0 / # incoming heat flux from plasma
  rad      real /0/ #
 rad_min   real /0/ #
  rad_max  real /0/ #
 # t1        real /0/ #
 #  t2       real /0/ #
 # t3        real /0/ #
 # tpulse    #      real /0 / # period of plasma pulse
 
**** ParticleFluxHeader:
implantation_model(1:nspc) _character*(1) /'S'/ # G: gaussian S:Step E: ERFC
implantation_depth(1:nspc) _real [m] /5e-9/ # implentation  depth
diagnostic_depth(1:nspc) _real  [m] /5e-9/ # diagnostic depth
implantation_width(1:nspc) _real [m] /5e-9/ # implentation width
enrg(1:nspc) _real [eV] /0.0/ # Impact energy of ionized species
inflx(1:nspc)  _real  [m^-2 s^- 1] /1.e20/ # influx of particles (may differ from nominal particle flux in pulsed_plasma mode)

 gas_pressure(1:nspc) _real  #
  gas_temp(1:nspc) _real  #
  mass(1:nspc) _real  #
  
***** MaterialHeader:
lambda real /1e-10/ #
cvlm real /1.0/ #
csrf real /1.0/ #
clng real /1.0/ #
thcond   real /0/ #
rho      real /0/ #
cp       real /0/ #
rhocp     real /0/ #
emiss    real /0/ #