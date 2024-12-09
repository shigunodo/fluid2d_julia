module Fluid2d

mutable struct GenStructCoord
  # number of grids
  NI::Int64
  NJ::Int64
  NB::Int64
  # coordinate: size (NI,NJ)
  x::Array{Float64, 2}
  y::Array{Float64, 2}
  # metrices devided by Jacobian: size (NI-2*NB+2,NJ-2*NB+2)
  ixs::Array{Float64, 2}
  iys::Array{Float64, 2}
  jxs::Array{Float64, 2}
  jys::Array{Float64, 2}
  # inverse of Jacobian: size (NI-1,NJ-1)
  s::Array{Float64, 2}
  # dx for CFL condition: size (NI-2*NB,NJ-2*NB)
  dx::Array{Float64, 2}
  # constructor
  GenStructCoord(NI,NJ,NB) = new(NI,NJ,NB,zeros(NI,NJ),zeros(NI,NJ),zeros(NI-2*NB+2,NJ-2*NB+2),zeros(NI-2*NB+2,NJ-2*NB+2),zeros(NI-2*NB+2,NJ-2*NB+2),zeros(NI-2*NB+2,NJ-2*NB+2),zeros(NI-1,NJ-1),zeros(NI-2*NB,NJ-2*NB))
end

mutable struct BasicVarHD
  # number of grids
  NI::Int64
  NJ::Int64
  NB::Int64
  # density
  rho::Array{Float64, 2}
  # x-component of velocity
  u::Array{Float64, 2}
  # y-component of velocity
  v::Array{Float64, 2}
  # total energy per volume
  e::Array{Float64, 2}
  # constructor
  BasicVarHD(NI, NJ, NB) = new(NI, NJ, NB, zeros(NI,NJ), zeros(NI,NJ), zeros(NI,NJ), zeros(NI,NJ))
end


struct ArrayBuffer
  arr_q0::Array{Float64, 3}
  arr_q1::Array{Float64, 3}
  arr_q2::Array{Float64, 3}
  arr_fi::Array{Float64, 3}
  arr_fj::Array{Float64, 3}
  arr_rhs::Array{Float64, 3}
  function ArrayBuffer(basic::BasicVarHD)
    (; NI, NB, NJ, NB) = basic
    arr_q0 = zeros(NB,NI-2*NB,NJ-2*NB)
    arr_q1 = zeros(NB,NI-2*NB,NJ-2*NB)
    arr_q2 = zeros(NB,NI-2*NB,NJ-2*NB)
    arr_fi = zeros(NB,NI-2*NB+1,NJ-2*NB)
    arr_fj = zeros(NB,NI-2*NB,NJ-2*NB+1)
    arr_rhs = zeros(NB,NI-2*NB,NJ-2*NB)
    new(arr_q0, arr_q1, arr_q2, arr_fi, arr_fj, arr_rhs)
  end
end

struct EulerEq
  # no member variables
end

struct IdealEoS
  # specific heat
  gam::Float64
  # constructor
  IdealEoS(gam) = new(gam)
end

# structure representing compressible invicid ideal gas
# without turbulence scheme
mutable struct IdealGas
  # general structured coordinate
  coord::GenStructCoord
  # basic variables of hydrodynamics /wo turbulence
  basic::BasicVarHD
  # compressible Euler equation
  eq::EulerEq
  # EoS of ideal gas
  eos::IdealEoS
  # constructor
  IdealGas(NI,NJ,NB,gam) = new(GenStructCoord(NI,NJ,NB),BasicVarHD(NI,NJ,NB),EulerEq(),IdealEoS(gam))
end

# export GenStructCoord, BasicVarHD, IdealEoS, EulerEq

include("Fluid2d/Fnd.jl")
include("Fluid2d/Eos.jl")
include("Fluid2d/Coordinate.jl")
include("Fluid2d/BC.jl")
include("Fluid2d/IO.jl")
include("Fluid2d/Conserved.jl")
include("Fluid2d/FluxScheme.jl")
include("Fluid2d/Eq.jl")
include("Fluid2d/Marching.jl")

using .Marching: march_ssprk3, calc_cfl
using .IO: input, output
using .Coordinate: calc_metrices_dx
export IdealGas, march_ssprk3, calc_cfl, input, output, calc_metrices_dx

end # module Fluid2d
