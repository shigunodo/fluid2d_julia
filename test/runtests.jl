using Test

using Fluid2d
using Fluid2d.Fnd: minmod
using Dates

@testset "Fluid2d" begin
    NI = 18
    NJ = 18
    NB = 4
    gamma = 1.4
    dir = joinpath(pkgdir(Fluid2d), "data")
    dir_o = dir
    f_settings = joinpath(dir, "settings.dat")
    f_coord = joinpath(dir, "coordinate.dat")

    start = Dates.now()

    # create fluid object
    fluid = IdealGas(NI,NJ,NB,gamma)

    # input coordinate
    input(f_coord, fluid.coord)
    # set metrices/Jacobian/dx for CFL
    calc_metrices_dx(fluid.coord)
    # input initial condition
    input(joinpath(dir_o, "b0000000.dat"), fluid.basic)
    t = 0.0
    # output settings
    @test size(fluid.coord.x) == (NI, NJ)
end

function minmod_orig(x, y)
  return sign(x) * max(min(abs(x), sign(x)*y), zero(x))
end

function minmod_orig(x1,x2,x3,x4)
  s1 = sign(x1)
  return 0.5*(s1+sign(x2))*abs(0.5*(s1+sign(x3))*0.5*(s1+sign(x4)))*min(abs(x1),abs(x2),abs(x3),abs(x4))
end

@testset "minmod(x,y)" begin
    y = rand()
    x = zero(y)
    @test minmod(0,0) == minmod_orig(0,0)
    @test minmod(0.0, y) == minmod_orig(0.0, y)
    @test minmod(0.0, -y) == minmod_orig(0.0, -y)
    for _ in 1:100
        x = rand()
        y = rand()
        @test minmod(x, y) == minmod_orig(x, y)
    end
end

@testset "minmod(x1, x2, x3, x4)" begin
    y = rand()
    x = zero(y)
    @test minmod(0,0) == minmod_orig(0,0)
    @test minmod(0.0, y) == minmod_orig(0.0, y)
    @test minmod(0.0, -y) == minmod_orig(0.0, -y)
    for _ in 1:100
        x = rand()
        y = rand()
        @test minmod_orig(x, y) == minmod_orig(x, y)
    end
end