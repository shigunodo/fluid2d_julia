using Fluid2d
using Dates

function iteration(tstep, dt_out, t, iter, fluid::IdealGas)
  qbuff = Fluid2d.ArrayBuffer(fluid.basic)
  while t < tstep * dt_out
    # CFL condition
    dt = calc_cfl(0.7, fluid.basic, fluid.coord, fluid.eos)
    # marching
    march_ssprk3(qbuff, dt, "periodical_in_i", "MP5_basic", "Roe_FDS", fluid.basic, fluid.coord, fluid.eos, fluid.eq)
    t = t + dt
    iter = iter + 1
  end
  return t, iter
end

function main()
  NI = 408
  NJ = 408
  NB = 4
  gamma = 1.4
  t_max = 3.0
  n_out = 100
  dir = joinpath(pkgdir(Fluid2d), "data") ####ワーキングディレクトリーを設定してください。####
  dir_o = dir     ####アウトプット用ディレクトリーを設定してください。####
  f_settings = joinpath(dir, "settings.dat") # 設定と log を出力するファイル
  f_coord = joinpath(dir, "coordinate.dat") # 計算格子の座標を取得するためのファイル

  start = Dates.now()

  # create fluid object
  fluid = IdealGas(NI,NJ,NB,gamma)

  # input coordinate
  input(f_coord, fluid.coord)
  # set metrices/Jacobian/dx for CFL
  calc_metrices_dx(fluid.coord)
  # input initial condition
  input(joinpath(dir_o, "b0000000.dat"), fluid.basic) # 初期条件が書かれたファイル
  t = 0.0
  # output settings
  output(f_settings, t_max, n_out, fluid.basic, fluid.eos, fluid.eq)

  dt_out = t_max / n_out
  # loop for output
  for tstep = 1:n_out
    # loop for marching
    iter = 0
    t, iter = iteration(tstep, dt_out, t, iter, fluid::IdealGas)
    # output
    now = Dates.now()
    cpu_time = floor(Int64, (now - start) / Dates.Millisecond(1000))
    rest_time = div(cpu_time * (n_out - tstep), tstep)
    output(dir_o, f_settings, tstep, t, iter, fluid.basic, cpu_time, rest_time)
  end

  print("Main program ended.")
  return
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end