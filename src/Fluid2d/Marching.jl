module Marching

import ..ArrayBuffer
import ..BasicVarHD
import ..GenStructCoord
import ..EulerEq
using ..Eos: calc_cs
using ..Conserved: calc_conserv, calc_basic
using ..Eq: calc_rhs
using ..BC: reflect_bc

# calc CFL condition
function calc_cfl(cfl_coeff, basic::BasicVarHD, coord::GenStructCoord, eos)

  NI = basic.NI
  NJ = basic.NJ
  NB = basic.NB

  nu = 1.0e+10
  for j = 1:NJ-2*NB
    for i = 1:NI-2*NB
      rho = basic.rho[i+NB,j+NB]
      u = basic.u[i+NB,j+NB]
      v = basic.v[i+NB,j+NB]
      e = basic.e[i+NB,j+NB]
      tmp = coord.dx[i,j] / (calc_cs(rho,u,v,e,eos) + sqrt(u*u+v*v))
      nu = min(nu, tmp)
    end
  end
  return cfl_coeff * nu
end

# marching dt with 3rd order SSP Runge-Kutta method
function march_ssprk3(arrbuff::ArrayBuffer, dt, bc_type, reconstruction, flux_scheme, basic::BasicVarHD, coord::GenStructCoord, eos, eq)

  NI = basic.NI
  NJ = basic.NJ
  NB = basic.NB

  # construct conservative var
  arr_q0 = arrbuff.arr_q0
  @inbounds for j = 1:NJ-2*NB
    for i = 1:NI-2*NB
      # inverse of Jacobian: averaging adjacent 4 values
      s_a = 0.25 * (coord.s[NB+i-1,NB+j-1] + coord.s[NB+i,NB+j-1] + coord.s[NB+i-1,NB+j] + coord.s[NB+i,NB+j])
      arr_q0[:,i,j] .= calc_conserv(basic.rho[NB+i,NB+j],basic.u[NB+i,NB+j],basic.v[NB+i,NB+j],basic.e[NB+i,NB+j],s_a)
    end
  end

  # 1st stage
  arr_q1 = arrbuff.arr_q1
  arr_q1 .= calc_rhs(arrbuff, reconstruction, flux_scheme, basic, coord, eos, eq)
  arr_q1 .= arr_q0 .+ dt .* arr_q1
  @inbounds for j = 1:NJ-2*NB
    for i = 1:NI-2*NB
      # inverse of Jacobian: averaging adjacent 4 values
      s_a = 0.25 * (coord.s[NB+i-1,NB+j-1] + coord.s[NB+i,NB+j-1] + coord.s[NB+i-1,NB+j] + coord.s[NB+i,NB+j])
      # update basic var
      basic.rho[NB+i,NB+j], basic.u[NB+i,NB+j], basic.v[NB+i,NB+j], basic.e[NB+i,NB+j] = calc_basic(@view(arr_q1[:,i,j]),s_a)
    end
  end
  reflect_bc(bc_type, basic, eos)

  # 2nd stage
  arr_q2 = arrbuff.arr_q2
  arr_q2 .= calc_rhs(arrbuff, reconstruction, flux_scheme, basic, coord, eos, eq)
  arr_q2 .= 0.75 .* arr_q0 .+ 0.25 .* (arr_q1 .+ dt .* arr_q2)
  @inbounds  for j = 1:NJ-2*NB
    for i = 1:NI-2*NB
      # inverse of Jacobian: averaging adjacent 4 values
      s_a = 0.25 * (coord.s[NB+i-1,NB+j-1] + coord.s[NB+i,NB+j-1] + coord.s[NB+i-1,NB+j] + coord.s[NB+i,NB+j])
      # update basic var
      basic.rho[NB+i,NB+j], basic.u[NB+i,NB+j], basic.v[NB+i,NB+j], basic.e[NB+i,NB+j] = calc_basic(@view(arr_q2[:,i,j]),s_a)
    end
  end
  reflect_bc(bc_type, basic, eos)

  # 3rd stage
  arr_q1 .= calc_rhs(arrbuff, reconstruction, flux_scheme, basic, coord, eos, eq)
  arr_q0 .= arr_q0 ./ 3.0 .+ 2.0 ./ 3.0 .* (arr_q2 .+ dt .* arr_q1)
  @inbounds for j = 1:NJ-2*NB
    for i = 1:NI-2*NB
      # inverse of Jacobian: averaging adjacent 4 values
      s_a = 0.25 * (coord.s[NB+i-1,NB+j-1] + coord.s[NB+i,NB+j-1] + coord.s[NB+i-1,NB+j] + coord.s[NB+i,NB+j])
      # update basic var
      basic.rho[NB+i,NB+j], basic.u[NB+i,NB+j], basic.v[NB+i,NB+j], basic.e[NB+i,NB+j] = calc_basic(@view(arr_q0[:,i,j]),s_a)
    end
  end
  reflect_bc(bc_type, basic, eos)

  return
end # function march_ssprk3

# export calc_cfl, march_ssprk3

end # module Marching