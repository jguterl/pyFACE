
module header
  use HeatFluxHeader
  use TimerHeader
  use Dims
  use SpeciesHeader
  use GridHeader
  use Grid
  use ReactionRatesHeader
  use bulkHeader
  use bulk
  use SurfaceHeader
  use Surface
  use materialHeader
  use MatTempHeader
  use MatTemp
  use IOHeader
  use BDFsolver
  use PhysConstsHeader
  use CapHeader
  use VerboseHeader
  use TimeHeader

  use ParticleFluxHeader
  use SolverHeader
  integer,parameter:: DP= kind(1.d0) ! precision
  
end module header
