module Fnd
# reconstruction schemes are defined

#=
@inline function minmod(x, y)
  sx = sign(x)
  return sx * max(min(abs(x), sx*y), zero(x))
end
=#

function minmod(x, y)
  if x > 0
	  sx = one(x)
      sx * max(min(x, sx*y), zero(x))
  elseif x < 0
	  sx = -one(x)
	  sx * max(min(-x, sx*y), zero(x))
  else
	  zero(x)
  end
end

#=
@inline function minmod(x1,x2,x3,x4)
  s1 = sign(x1)
  return 0.5*(s1+sign(x2))*abs(0.5*(s1+sign(x3))*0.5*(s1+sign(x4)))*min(abs(x1),abs(x2),abs(x3),abs(x4))
end
=#

function minmod(x1,x2,x3,x4)
  s1 = sign(x1)
  return 0.5*0.5*0.5*
    (s1+sign(x2)) *
    abs((s1+sign(x3))*
      (s1+sign(x4))) *
     min(min(abs(x1),abs(x2)),min(abs(x3),abs(x4)))
end

# MUSCL-minmod
@inline function muscl_minmod(qm,q,qp,q2p)
  k = 1.0/3.0
  b = (3.0-k)/(1.0-k)
  # calc of L
  dp = qp - q
  dm = q - qm
  ddp = minmod(dp,b*dm)
  ddm = minmod(dm,b*dp)
  q_l = q + 0.25*((1.0-k)*ddm + (1.0+k)*ddp)
  # calc of R
  dp = q2p - qp
  dm = qp - q
  ddp = minmod(dp,b*dm)
  ddm = minmod(dm,b*dp)
  q_r = qp - 0.25*((1.0-k)*ddp + (1.0+k)*ddm)
  return q_l, q_r
end

@inline function median(x,y,z)
  return x + minmod(y-x,z-x)
end

@inline function mp5_sub(q2m,qm,q,qp,q2p)
  alp = 2.0
  q_l = (2.0 * q2m - 13.0 * qm + 47.0 * q + 27.0 * qp - 3.0 * q2p) / 60.0
  q_mp = q + minmod(qp - q, alp * (q - qm))
  if (q_l-q)*(q_l-q_mp) > 1.0e-10
    dm = q2m + q - 2.0 * qm
    d = qm + qp - 2.0 * q
    dp = q + q2p - 2.0 * qp
    d_mm = minmod(4.0 * dm - d, 4.0 * d - dm, dm, d)
    d_mp = minmod(4.0 * d - dp, 4.0 * dp - d, d, dp)
    q_ul = q + alp * (q - qm)
    q_md = 0.5 * (q + qp) - 0.5 * d_mp
    q_lc = q + 0.5 * (q - qm) + 4.0 / 3.0 * d_mm
    qmin = max(min(min(q,qp),min(q_md)),min(min(q,q_ul),min(q_lc)))
    qmax = min(max(max(q,qp),max(q_md)),max(max(q,q_ul),max(q_lc)))
    q_l = median(q_l,qmin,qmax)
  end
  return q_l
end

# MP5
@inline function mp5(q2m,qm,q,qp,q2p,q3p)
  q_l = mp5_sub(q2m,qm,q,qp,q2p)
  q_r = mp5_sub(q3p,q2p,qp,q,qm)
  return q_l, q_r
end

export muscl_minmod, mp5

end # module Fnd