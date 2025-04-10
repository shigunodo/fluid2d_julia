module Eq

import ..ArrayBuffer
import ..BasicVarHD
import ..GenStructCoord
import ..EulerEq
using ..FluxScheme: reconst_by_basic_muscl, reconst_by_basic_mp5, calc_num_flux

# calc RHS with selected reconstruction/flux scheme
function calc_rhs(arrbuff::ArrayBuffer, reconstruction, flux_scheme, basic::BasicVarHD, coord::GenStructCoord, eos, eq::EulerEq)

  NI = basic.NI
  NJ = basic.NJ
  NB = basic.NB
  arr_fi = arrbuff.arr_fi
  arr_fj = arrbuff.arr_fj

  # i-direction
  # evaluating numerical flux at (i+0.5,j)
  @inbounds for j = 1:NJ-2*NB
    for i = 1:NI-2*NB+1

      # inverse of Jacobian is evaluated
      # by averaging values at (i+0.5,j+0.5) and (i+0.5,j-0.5)
      s_a = 0.5 * (coord.s[NB+i-1,NB+j-1] + coord.s[NB+i-1,NB+j])

      # metrices devided by Jacobian are evaluated
      # by averaging values at (i,j) and (i+1,j)
      ixs_a = 0.5 * (coord.ixs[i,j+1] + coord.ixs[i+1,j+1])
      iys_a = 0.5 * (coord.iys[i,j+1] + coord.iys[i+1,j+1])

      # reconstruction with selected scheme
      if reconstruction == "MUSCL_minmod_basic"
        rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r = reconst_by_basic_muscl(basic.rho[NB+i-2:NB+i+1,NB+j],basic.u[NB+i-2:NB+i+1,NB+j],basic.v[NB+i-2:NB+i+1,NB+j],basic.e[NB+i-2:NB+i+1,NB+j])
      elseif reconstruction == "MP5_basic"
        rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r = reconst_by_basic_mp5(
          @view(basic.rho[NB+i-3:NB+i+2,NB+j]),
          @view(basic.u[NB+i-3:NB+i+2,NB+j]),
          @view(basic.v[NB+i-3:NB+i+2,NB+j]),
          @view(basic.e[NB+i-3:NB+i+2,NB+j])
        )
      else
        throw(DomainError(reconstruction, "reconstruction scheme failed to be identified"))
      end

      # print(basic.rho[NB+i-2:NB+i+1,NB+j])

      # convective flux
      arr_fi[:,i,j] .= calc_num_flux(flux_scheme,rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r,ixs_a,iys_a,s_a,eos)

    end # for i
  end # for j

  # j-direction
  # evaluating numerical flux at (i,j+0.5)
  @inbounds for j = 1:NJ-2*NB+1
    for i = 1:NI-2*NB

      # inverse of Jacobian is evaluated
      # by averaging values at (i+0.5,j+0.5) and (i-0.5,j+0.5)
      s_a = 0.5 * (coord.s[NB+i-1,NB+j-1] + coord.s[NB+i,NB+j-1])

      # metrices devided by Jacobian are evaluated
      # by averaging values at (i,j) and (i,j+1)
      jxs_a = 0.5 * (coord.jxs[i+1,j] + coord.jxs[i+1,j+1])
      jys_a = 0.5 * (coord.jys[i+1,j] + coord.jys[i+1,j+1])

      # reconstruction with selected scheme
      if reconstruction == "MUSCL_minmod_basic"
        rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r = reconst_by_basic_muscl(basic.rho[NB+i,NB+j-2:NB+j+1],basic.u[NB+i,NB+j-2:NB+j+1],basic.v[NB+i,NB+j-2:NB+j+1],basic.e[NB+i,NB+j-2:NB+j+1])
      elseif reconstruction == "MP5_basic"
        rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r = reconst_by_basic_mp5(
        @view(basic.rho[NB+i,NB+j-3:NB+j+2]),
        @view(basic.u[NB+i,NB+j-3:NB+j+2]),
        @view(basic.v[NB+i,NB+j-3:NB+j+2]),
        @view(basic.e[NB+i,NB+j-3:NB+j+2]),
        )
      else
        throw(DomainError(reconstruction, "reconstruction scheme failed to be identified"))
      end

      # convective flux
      arr_fj[:,i,j] .= calc_num_flux(flux_scheme,rho_l,u_l,v_l,e_l,rho_r,u_r,v_r,e_r,jxs_a,jys_a,s_a,eos)

    end # for i
  end # for j

  # evaluating RHS
  fill!(arrbuff.arr_rhs, 0.0)
  arr_rhs = arrbuff.arr_rhs
  arr_rhs .+= @view arr_fi[:,1:NI-2*NB,:]
  arr_rhs .-= @view arr_fi[:,2:NI-2*NB+1,:]
  arr_rhs .+= @view arr_fj[:,:,1:NJ-2*NB]
  arr_rhs .-= @view arr_fj[:,:,2:NJ-2*NB+1]
  return arr_rhs
end # function calc_rhs

# export calc_rhs

end # module Eq