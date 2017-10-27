# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
#  Translation to Java and to MultiThreaded Code
#           Michael A. Frumkin
#           Mathew Schultz
#  Translation to Ruby
#           T. Nose
#           H. Tomari
# -------------------------------------------------------------------------
require "mda"
require "Runnable"
class SPBase < Runnable
  @@SPBaseConstInit = false
  def initialize(clss, np)
    @CLASS = 'S'
    @IMAX = 0
    @JMAX = 0
    @KMAX = 0
    @problem_size = 0
    @nx2 = 0
    @ny2 = 0
    @nz2 = 0
    @grid_points = [0, 0, 0]
    @niter_default = 0
    @dt_default = 0.0
    if not @@SPBaseConstInit then
      @@BMName = "SP"
      @@SPBaseConstInit = true
      @@tx1 = 0.0
      @@tx2 = 0.0
      @@tx3 = 0.0
      @@ty1 = 0.0
      @@ty2 = 0.0
      @@ty3 = 0.0
      @@tz1 = 0.0
      @@tz2 = 0.0
      @@tz3 = 0.0
      @@dx1 = 0.0
      @@dx2 = 0.0
      @@dx3 = 0.0
      @@dx4 = 0.0
      @@dx5 = 0.0
      @@dy1 = 0.0
      @@dy2 = 0.0
      @@dy3 = 0.0
      @@dy4 = 0.0
      @@dy5 = 0.0
      @@dz1 = 0.0
      @@dz2 = 0.0
      @@dz3 = 0.0
      @@dz4 = 0.0
      @@dz5 = 0.0
      @@dssp = 0.0
      @@dt = 0.0
      @@dxmax = 0.0
      @@dymax = 0.0
      @@dzmax = 0.0
      @@xxcon1 = 0.0
      @@xxcon2 = 0.0
      @@xxcon3 = 0.0
      @@xxcon4 = 0.0
      @@xxcon5 = 0.0
      @@dx1tx1 = 0.0
      @@dx2tx1 = 0.0
      @@dx3tx1 = 0.0
      @@dx4tx1 = 0.0
      @@dx5tx1 = 0.0
      @@yycon1 = 0.0
      @@yycon2 = 0.0
      @@yycon3 = 0.0
      @@yycon4 = 0.0
      @@yycon5 = 0.0
      @@dy1ty1 = 0.0
      @@dy2ty1 = 0.0
      @@dy3ty1 = 0.0
      @@dy4ty1 = 0.0
      @@dy5ty1 = 0.0
      @@zzcon1 = 0.0
      @@zzcon2 = 0.0
      @@zzcon3 = 0.0
      @@zzcon4 = 0.0
      @@zzcon5 = 0.0
      @@dz1tz1 = 0.0
      @@dz2tz1 = 0.0
      @@dz3tz1 = 0.0
      @@dz4tz1 = 0.0
      @@dz5tz1 = 0.0
      @@dnxm1 = 0.0
      @@dnym1 = 0.0
      @@dnzm1 = 0.0
      @@c1c2 = 0.0
      @@c1c5 = 0.0
      @@c3c4 = 0.0
      @@c1345 = 0.0
      @@conz1 = 0.0
      @@c1 = 0.0
      @@c2 = 0.0
      @@c3 = 0.0
      @@c4 = 0.0
      @@c5 = 0.0
      @@c4dssp = 0.0
      @@c5dssp = 0.0
      @@dtdssp = 0.0
      @@dttx1 = 0.0
      @@bt = 0.0
      @@dttx2 = 0.0
      @@dtty1 = 0.0
      @@dtty2 = 0.0
      @@dttz1 = 0.0
      @@dttz2 = 0.0
      @@c2dttx1 = 0.0
      @@c2dtty1 = 0.0
      @@c2dttz1 = 0.0
      @@comz1 = 0.0
      @@comz4 = 0.0
      @@comz5 = 0.0
      @@comz6 = 0.0
      @@c3c4tx3 = 0.0
      @@c3c4ty3 = 0.0
      @@c3c4tz3 = 0.0
      @@c2iv = 0.0
      @@con43 = 0.0
      @@con16 = 0.0
      @@t_total = 1
      @@t_rhsx = 2
      @@t_rhsy = 3
      @@t_rhsz = 4
      @@t_rhs = 5
      @@t_xsolve = 6
      @@t_ysolve = 7
      @@t_zsolve = 8
      @@t_rdis1 = 9
      @@t_rdis2 = 10
      @@t_txinvr = 11
      @@t_pinvr = 12
      @@t_ninvr = 13
      @@t_tzetar = 14
      @@t_add = 15
      @@t_last = 15
    end
    @ce = [2.0, 1.0, 2.0, 2.0, 5.0,
           0.0, 0.0, 2.0, 2.0, 4.0,
           0.0, 0.0, 0.0, 0.0, 3.0, 
           4.0, 0.0, 0.0, 0.0, 2.0,
           5.0, 1.0, 0.0, 0.0, 0.1,
           3.0, 2.0, 2.0, 2.0, 0.4,
           0.5, 3.0, 3.0, 3.0, 0.3,
           0.02, 0.01, 0.04, 0.03, 0.05,
           0.01, 0.03, 0.03, 0.05, 0.04,
           0.03, 0.02, 0.05, 0.04, 0.03, 
           0.5, 0.4, 0.3, 0.2, 0.1, 
           0.4, 0.3, 0.5, 0.1, 0.3, 
           0.3, 0.5, 0.4, 0.3, 0.2]

    @timeron = false
    @timer = NPB3_0_RUB::Timer.new()

    @CLASS = clss
    @num_threads = np
    case clss
    when 'S'
      @IMAX = @JMAX = @KMAX = @problem_size = @grid_points[0] = @grid_points[1] = @grid_points[2] = 12
      @dt_default = 0.015
      @niter_default = 100
    when 'W'
      @IMAX = @JMAX = @KMAX = @problem_size = @grid_points[0] = @grid_points[1] = @grid_points[2] = 36
      @dt_default = 0.0015
      @niter_default = 400
    when 'A'
      @IMAX = @JMAX = @KMAX = @problem_size = @grid_points[0] = @grid_points[1] = @grid_points[2] = 64
      @dt_default = 0.0015
      @niter_default = 400
    when 'B'
      @IMAX = @JMAX = @KMAX = @problem_size = @grid_points[0] = @grid_points[1] = @grid_points[2] = 102
      @dt_default = 0.001
      @niter_default = 400
    when 'C'
      @IMAX = @JMAX = @KMAX = @problem_size = @grid_points[0] = @grid_points[1] = @grid_points[2] = 162
      @dt_default = 0.00067
      @niter_default = 400
    end

    @isize1 = 5
    @jsize1 = 5 * (@IMAX + 1)
    @ksize1 = 5 * (@IMAX + 1) * (@JMAX + 1)
    @u = Array.new(5 * (@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)
    @rhs = Array.new(5 * (@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)
    @forcing = Array.new(5 * (@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)

    @jsize2 = (@IMAX + 1)
    @ksize2 = (@IMAX + 1) * (@JMAX + 1)
    @us = Array.new((@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)
    @vs = Array.new((@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)
    @ws = Array.new((@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)
    @qs = Array.new((@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)
    @rho_i = Array.new((@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)
    @speed = Array.new((@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)
    @square = Array.new((@IMAX + 1) * (@JMAX + 1) * @KMAX, 0.0)

    @jsize3 = @problem_size
    @ue = Array.new(@problem_size * 5, 0.0)
    @buf = Array.new(@problem_size * 5, 0.0)

    @jsize4 = 5

    @lhs = Array.new(5 * (@problem_size + 1), 0.0)
    @lhsp = Array.new(5 * (@problem_size + 1), 0.0)
    @lhsm = Array.new(5 * (@problem_size + 1), 0.0)

    @cv = Array.new(@problem_size, 0.0)
    @rhon = Array.new(@problem_size, 0.0)
    @rhos = Array.new(@problem_size, 0.0)
    @rhoq = Array.new(@problem_size, 0.0)
    @cuf = Array.new(@problem_size, 0.0)
    @q = Array.new(@problem_size, 0.0)
  end

  def set_interval(problem_size, interval)
    interval[0] = problem_size / @num_threads
    for i in 1..@num_threads-1 do
      interval[i] = interval[0]
    end
    remainder = problem_size % @num_threads; 
    for i in 0..remainder-1 do
      interval[i] += 1
    end
  end

  def set_partition(start, interval, prtn)
    prtn[0][0] = start
    if start == 0 then
      prtn[0][1] = interval[0] - 1
    else
      prtn[0][1] = interval[0]
    end
    for i in 1..interval.length-1 do
        prtn[i][0] = prtn[i - 1][1] + 1
        prtn[i][1] = prtn[i - 1][1] + interval[i]
    end
  end

  def dmax1(a, b)
    if a < b then
      return b
    else
      return a
    end
  end

  def dmax1ex(a, b, c, d)
    return dmax1(dmax1(a, b), dmax1(c, d))
  end

  def checksum(array, size, arrayname, stop)
    sum = 0; 
    for i in 0..size-1 do
      sum += array[i]
    end
    puts("array:" + arrayname + " checksum is: " + sum.to_s)
    if stop then
      exit(0)
    else
    end
  end

  def exact_solution(xi, eta, zeta, dtemp, offset)
    for m in 0..4 do
        dtemp[m + offset] = @ce[m + 0 * 5] + xi * (@ce[m + 1 * 5] + xi * (@ce[m + 4 * 5] + xi * (@ce[m + 7 * 5] + xi * @ce[m + 10 * 5]))) + eta * (@ce[m + 2 * 5] + eta * (@ce[m + 5 * 5] + eta * (@ce[m + 8 * 5] + eta * @ce[m + 11 * 5]))) + zeta * (@ce[m + 3 * 5] + zeta * (@ce[m + 6 * 5] + zeta * (@ce[m + 9 * 5] + zeta * @ce[m + 12 * 5])))
    end
  end

  def lhsinit(size)
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    i = 0; 
    while i <= size do
        for n in 0..4 do
            @lhs[n + i * @jsize4] = 0.0
            @lhsp[n + i * @jsize4] = 0.0
            @lhsm[n + i * @jsize4] = 0.0
        end
      i += size
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    i = 0; 
    while i <= size do
        @lhs[2 + i * @jsize4] = 1.0
        @lhsp[2 + i * @jsize4] = 1.0
        @lhsm[2 + i * @jsize4] = 1.0
      i += size
    end
  end

  def _initialize()
    
    pface = Array.new(5 * 3 * 2, 0.0); temp = Array.new(5, 0.0); 
    for k in 0..@grid_points[2]-1 do
        for j in 0..@grid_points[1]-1 do
            for i in 0..@grid_points[0]-1 do
                @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = 1.0
                @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = 0.0
                @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = 0.0
                @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = 0.0
                @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = 1.0
            end
        end
    end
    for k in 0..@grid_points[2]-1 do
        zeta = k * @@dnzm1
        for j in 0..@grid_points[1]-1 do
            eta = j * @@dnym1
            for i in 0..@grid_points[0]-1 do
                xi = i * @@dnxm1
                for ix in 0..1 do
                    pxi = ix
                    exact_solution(pxi, eta, zeta, pface, 0 + 0 * 5 + ix * 15)
                end
                for iy in 0..1 do
                    peta = iy
                    exact_solution(xi, peta, zeta, pface, 0 + 1 * 5 + iy * 15)
                end
                for iz in 0..1 do
                    pzeta = iz
                    exact_solution(xi, eta, pzeta, pface, 0 + 2 * 5 + iz * 15)
                end
                for m in 0..4 do
                    pxi = xi * pface[m + 0 * 5 + 1 * 15] + (1.0 - xi) * pface[m + 0 * 5 + 0 * 15]
                    peta = eta * pface[m + 1 * 5 + 1 * 15] + (1.0 - eta) * pface[m + 1 * 5 + 0 * 15]
                    pzeta = zeta * pface[m + 2 * 5 + 1 * 15] + (1.0 - zeta) * pface[m + 2 * 5 + 0 * 15]
                    @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = pxi + peta + pzeta - pxi * peta - pxi * pzeta - peta * pzeta + pxi * peta * pzeta
                end
            end
        end
    end
    xi = 0.0
    i = 0
    for k in 0..@grid_points[2]-1 do
        zeta = k * @@dnzm1
        for j in 0..@grid_points[1]-1 do
            eta = j * @@dnym1
            exact_solution(xi, eta, zeta, temp, 0)
            for m in 0..4 do
                @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = temp[m]
            end
        end
    end
    xi = 1.0
    i = @grid_points[0] - 1
    for k in 0..@grid_points[2]-1 do
        zeta = k * @@dnzm1
        for j in 0..@grid_points[1]-1 do
            eta = j * @@dnym1
            exact_solution(xi, eta, zeta, temp, 0)
            for m in 0..4 do
                @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = temp[m]
            end
        end
    end
    eta = 0.0
    j = 0
    for k in 0..@grid_points[2]-1 do
        zeta = k * @@dnzm1
        for i in 0..@grid_points[0]-1 do
            xi = i * @@dnxm1
            exact_solution(xi, eta, zeta, temp, 0)
            for m in 0..4 do
                @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = temp[m]
            end
        end
    end
    eta = 1.0
    j = @grid_points[1] - 1
    for k in 0..@grid_points[2]-1 do
        zeta = k * @@dnzm1
        for i in 0..@grid_points[0]-1 do
            xi = i * @@dnxm1
            exact_solution(xi, eta, zeta, temp, 0)
            for m in 0..4 do
                @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = temp[m]
            end
        end
    end
    zeta = 0.0
    k = 0
    for i in 0..@grid_points[0]-1 do
        xi = i * @@dnxm1
        for j in 0..@grid_points[1]-1 do
            eta = j * @@dnym1
            exact_solution(xi, eta, zeta, temp, 0)
            for m in 0..4 do
                @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = temp[m]
            end
        end
    end
    zeta = 1.0
    k = @grid_points[2] - 1
    for i in 0..@grid_points[0]-1 do
        xi = i * @@dnxm1
        for j in 0..@grid_points[1]-1 do
            eta = j * @@dnym1
            exact_solution(xi, eta, zeta, temp, 0)
            for m in 0..4 do
                @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = temp[m]
            end
        end
    end
  end

  def set_constants(ndid)
    @ce[0] = 2.0 * (1.0 + (ndid) * 0.01)
    @@c1 = 1.4
    @@c2 = 0.4
    @@c3 = 0.1
    @@c4 = 1.0
    @@c5 = 1.4
    @@bt = Math.sqrt(0.5)
    @@dnxm1 = 1.0 / (@grid_points[0] - 1)
    @@dnym1 = 1.0 / (@grid_points[1] - 1)
    @@dnzm1 = 1.0 / (@grid_points[2] - 1)
    @@c1c2 = @@c1 * @@c2
    @@c1c5 = @@c1 * @@c5
    @@c3c4 = @@c3 * @@c4
    @@c1345 = @@c1c5 * @@c3c4
    @@conz1 = (1.0 - @@c1c5)
    @@tx1 = 1.0 / (@@dnxm1 * @@dnxm1)
    @@tx2 = 1.0 / (2.0 * @@dnxm1)
    @@tx3 = 1.0 / @@dnxm1
    @@ty1 = 1.0 / (@@dnym1 * @@dnym1)
    @@ty2 = 1.0 / (2.0 * @@dnym1)
    @@ty3 = 1.0 / @@dnym1
    @@tz1 = 1.0 / (@@dnzm1 * @@dnzm1)
    @@tz2 = 1.0 / (2.0 * @@dnzm1)
    @@tz3 = 1.0 / @@dnzm1
    @@dx1 = 0.75
    @@dx2 = 0.75
    @@dx3 = 0.75
    @@dx4 = 0.75
    @@dx5 = 0.75
    @@dy1 = 0.75
    @@dy2 = 0.75
    @@dy3 = 0.75
    @@dy4 = 0.75
    @@dy5 = 0.75
    @@dz1 = 1.0
    @@dz2 = 1.0
    @@dz3 = 1.0
    @@dz4 = 1.0
    @@dz5 = 1.0
    @@dxmax = dmax1(@@dx3, @@dx4)
    @@dymax = dmax1(@@dy2, @@dy4)
    @@dzmax = dmax1(@@dz2, @@dz3)
    @@dssp = 0.25 * dmax1(@@dx1, dmax1(@@dy1, @@dz1))
    @@c4dssp = 4.0 * @@dssp
    @@c5dssp = 5.0 * @@dssp
    @@dttx1 = @@dt * @@tx1
    @@dttx2 = @@dt * @@tx2
    @@dtty1 = @@dt * @@ty1
    @@dtty2 = @@dt * @@ty2
    @@dttz1 = @@dt * @@tz1
    @@dttz2 = @@dt * @@tz2
    @@c2dttx1 = 2.0 * @@dttx1
    @@c2dtty1 = 2.0 * @@dtty1
    @@c2dttz1 = 2.0 * @@dttz1
    @@dtdssp = @@dt * @@dssp
    @@comz1 = @@dtdssp
    @@comz4 = 4.0 * @@dtdssp
    @@comz5 = 5.0 * @@dtdssp
    @@comz6 = 6.0 * @@dtdssp
    @@c3c4tx3 = @@c3c4 * @@tx3
    @@c3c4ty3 = @@c3c4 * @@ty3
    @@c3c4tz3 = @@c3c4 * @@tz3
    @@dx1tx1 = @@dx1 * @@tx1
    @@dx2tx1 = @@dx2 * @@tx1
    @@dx3tx1 = @@dx3 * @@tx1
    @@dx4tx1 = @@dx4 * @@tx1
    @@dx5tx1 = @@dx5 * @@tx1
    @@dy1ty1 = @@dy1 * @@ty1
    @@dy2ty1 = @@dy2 * @@ty1
    @@dy3ty1 = @@dy3 * @@ty1
    @@dy4ty1 = @@dy4 * @@ty1
    @@dy5ty1 = @@dy5 * @@ty1
    @@dz1tz1 = @@dz1 * @@tz1
    @@dz2tz1 = @@dz2 * @@tz1
    @@dz3tz1 = @@dz3 * @@tz1
    @@dz4tz1 = @@dz4 * @@tz1
    @@dz5tz1 = @@dz5 * @@tz1
    @@c2iv = 2.5
    @@con43 = 4.0 / 3.0
    @@con16 = 1.0 / 6.0
    @@xxcon1 = @@c3c4tx3 * @@con43 * @@tx3
    @@xxcon2 = @@c3c4tx3 * @@tx3
    @@xxcon3 = @@c3c4tx3 * @@conz1 * @@tx3
    @@xxcon4 = @@c3c4tx3 * @@con16 * @@tx3
    @@xxcon5 = @@c3c4tx3 * @@c1c5 * @@tx3
    @@yycon1 = @@c3c4ty3 * @@con43 * @@ty3
    @@yycon2 = @@c3c4ty3 * @@ty3
    @@yycon3 = @@c3c4ty3 * @@conz1 * @@ty3
    @@yycon4 = @@c3c4ty3 * @@con16 * @@ty3
    @@yycon5 = @@c3c4ty3 * @@c1c5 * @@ty3
    @@zzcon1 = @@c3c4tz3 * @@con43 * @@tz3
    @@zzcon2 = @@c3c4tz3 * @@tz3
    @@zzcon3 = @@c3c4tz3 * @@conz1 * @@tz3
    @@zzcon4 = @@c3c4tz3 * @@con16 * @@tz3
    @@zzcon5 = @@c3c4tz3 * @@c1c5 * @@tz3
  end

# *** public ***
  def BMName
    return @@BMName
  end
  def BMName=(val)
    @@BMName = val
  end

  def clss
    return @CLASS
  end

  attr_accessor :timeron

  attr_accessor :timer

  def t_total
    return @@t_total
  end
  def t_total=(val)
    @@t_total = val
  end

  def t_rhsx
    return @@t_rhsx
  end
  def t_rhsx=(val)
    @@t_rhsx = val
  end

  def t_rhsy
    return @@t_rhsy
  end
  def t_rhsy=(val)
    @@t_rhsy = val
  end

  def t_rhsz
    return @@t_rhsz
  end
  def t_rhsz=(val)
    @@t_rhsz = val
  end

  def t_rhs
    return @@t_rhs
  end
  def t_rhs=(val)
    @@t_rhs = val
  end

  def t_xsolve
    return @@t_xsolve
  end
  def t_xsolve=(val)
    @@t_xsolve = val
  end

  def t_ysolve
    return @@t_ysolve
  end
  def t_ysolve=(val)
    @@t_ysolve = val
  end

  def t_zsolve
    return @@t_zsolve
  end
  def t_zsolve=(val)
    @@t_zsolve = val
  end

  def t_rdis1
    return @@t_rdis1
  end
  def t_rdis1=(val)
    @@t_rdis1 = val
  end

  def t_rdis2
    return @@t_rdis2
  end
  def t_rdis2=(val)
    @@t_rdis2 = val
  end

  def t_txinvr
    return @@t_txinvr
  end
  def t_txinvr=(val)
    @@t_txinvr = val
  end

  def t_pinvr
    return @@t_pinvr
  end
  def t_pinvr=(val)
    @@t_pinvr = val
  end

  def t_ninvr
    return @@t_ninvr
  end
  def t_ninvr=(val)
    @@t_ninvr = val
  end

  def t_tzetar
    return @@t_tzetar
  end
  def t_tzetar=(val)
    @@t_tzetar = val
  end

  def t_add
    return @@t_add
  end
  def t_add=(val)
    @@t_add = val
  end

  def t_last
    return @@t_last
  end
  def t_last=(val)
    @@t_last = val
  end

# *** protected ***
  def imax
    return @IMAX
  end
  def imax=(val)
    @IMAX = val
  end
  def jmax
    return @JMAX
  end
  def jmax=(val)
    @JMAX = val
  end
  def kmax
    return @KMAX
  end
  def kmax=(val)
    @KMAX = val
  end

  attr_accessor :problem_size

  protected :problem_size
  attr_accessor :nx2

  protected :nx2
  attr_accessor :ny2

  protected :ny2
  attr_accessor :nz2

  protected :nz2
  attr_accessor :grid_points

  protected :grid_points
  attr_accessor :niter_default

  protected :niter_default
  attr_accessor :dt_default

  protected :dt_default
  attr_accessor :u

  protected :u
  attr_accessor :rhs

  protected :rhs
  attr_accessor :forcing

  protected :forcing
  attr_accessor :isize1

  protected :isize1
  attr_accessor :jsize1

  protected :jsize1
  attr_accessor :ksize1

  protected :ksize1
  attr_accessor :us

  protected :us
  attr_accessor :vs

  protected :vs
  attr_accessor :ws

  protected :ws
  attr_accessor :qs

  protected :qs
  attr_accessor :rho_i

  protected :rho_i
  attr_accessor :speed

  protected :speed
  attr_accessor :square

  protected :square
  attr_accessor :jsize2

  protected :jsize2
  attr_accessor :ksize2

  protected :ksize2
  attr_accessor :ue

  protected :ue
  attr_accessor :buf

  protected :buf
  attr_accessor :jsize3

  protected :jsize3
  attr_accessor :lhs

  protected :lhs
  attr_accessor :lhsp

  protected :lhsp
  attr_accessor :lhsm

  protected :lhsm
  attr_accessor :jsize4

  protected :jsize4
  attr_accessor :cv

  protected :cv
  attr_accessor :rhon

  protected :rhon
  attr_accessor :rhos

  protected :rhos
  attr_accessor :rhoq

  protected :rhoq
  attr_accessor :cuf

  protected :cuf
  attr_accessor :q

  protected :q
  def tx1
    return @@tx1
  end
  def tx1=(val)
    @@tx1 = val
  end

  protected :tx1
  def tx2
    return @@tx2
  end
  def tx2=(val)
    @@tx2 = val
  end

  protected :tx2
  def tx3
    return @@tx3
  end
  def tx3=(val)
    @@tx3 = val
  end

  protected :tx3
  def ty1
    return @@ty1
  end
  def ty1=(val)
    @@ty1 = val
  end

  protected :ty1
  def ty2
    return @@ty2
  end
  def ty2=(val)
    @@ty2 = val
  end

  protected :ty2
  def ty3
    return @@ty3
  end
  def ty3=(val)
    @@ty3 = val
  end

  protected :ty3
  def tz1
    return @@tz1
  end
  def tz1=(val)
    @@tz1 = val
  end

  protected :tz1
  def tz2
    return @@tz2
  end
  def tz2=(val)
    @@tz2 = val
  end

  protected :tz2
  def tz3
    return @@tz3
  end
  def tz3=(val)
    @@tz3 = val
  end

  protected :tz3
  def dx1
    return @@dx1
  end
  def dx1=(val)
    @@dx1 = val
  end

  protected :dx1
  def dx2
    return @@dx2
  end
  def dx2=(val)
    @@dx2 = val
  end

  protected :dx2
  def dx3
    return @@dx3
  end
  def dx3=(val)
    @@dx3 = val
  end

  protected :dx3
  def dx4
    return @@dx4
  end
  def dx4=(val)
    @@dx4 = val
  end

  protected :dx4
  def dx5
    return @@dx5
  end
  def dx5=(val)
    @@dx5 = val
  end

  protected :dx5
  def dy1
    return @@dy1
  end
  def dy1=(val)
    @@dy1 = val
  end

  protected :dy1
  def dy2
    return @@dy2
  end
  def dy2=(val)
    @@dy2 = val
  end

  protected :dy2
  def dy3
    return @@dy3
  end
  def dy3=(val)
    @@dy3 = val
  end

  protected :dy3
  def dy4
    return @@dy4
  end
  def dy4=(val)
    @@dy4 = val
  end

  protected :dy4
  def dy5
    return @@dy5
  end
  def dy5=(val)
    @@dy5 = val
  end

  protected :dy5
  def dz1
    return @@dz1
  end
  def dz1=(val)
    @@dz1 = val
  end

  protected :dz1
  def dz2
    return @@dz2
  end
  def dz2=(val)
    @@dz2 = val
  end

  protected :dz2
  def dz3
    return @@dz3
  end
  def dz3=(val)
    @@dz3 = val
  end

  protected :dz3
  def dz4
    return @@dz4
  end
  def dz4=(val)
    @@dz4 = val
  end

  protected :dz4
  def dz5
    return @@dz5
  end
  def dz5=(val)
    @@dz5 = val
  end

  protected :dz5
  def dssp
    return @@dssp
  end
  def dssp=(val)
    @@dssp = val
  end

  protected :dssp
  def dt
    return @@dt
  end
  def dt=(val)
    @@dt = val
  end

  protected :dt
  def dxmax
    return @@dxmax
  end
  def dxmax=(val)
    @@dxmax = val
  end

  protected :dxmax
  def dymax
    return @@dymax
  end
  def dymax=(val)
    @@dymax = val
  end

  protected :dymax
  def dzmax
    return @@dzmax
  end
  def dzmax=(val)
    @@dzmax = val
  end

  protected :dzmax
  def xxcon1
    return @@xxcon1
  end
  def xxcon1=(val)
    @@xxcon1 = val
  end

  protected :xxcon1
  def xxcon2
    return @@xxcon2
  end
  def xxcon2=(val)
    @@xxcon2 = val
  end

  protected :xxcon2
  def xxcon3
    return @@xxcon3
  end
  def xxcon3=(val)
    @@xxcon3 = val
  end

  protected :xxcon3
  def xxcon4
    return @@xxcon4
  end
  def xxcon4=(val)
    @@xxcon4 = val
  end

  protected :xxcon4
  def xxcon5
    return @@xxcon5
  end
  def xxcon5=(val)
    @@xxcon5 = val
  end

  protected :xxcon5
  def dx1tx1
    return @@dx1tx1
  end
  def dx1tx1=(val)
    @@dx1tx1 = val
  end

  protected :dx1tx1
  def dx2tx1
    return @@dx2tx1
  end
  def dx2tx1=(val)
    @@dx2tx1 = val
  end

  protected :dx2tx1
  def dx3tx1
    return @@dx3tx1
  end
  def dx3tx1=(val)
    @@dx3tx1 = val
  end

  protected :dx3tx1
  def dx4tx1
    return @@dx4tx1
  end
  def dx4tx1=(val)
    @@dx4tx1 = val
  end

  protected :dx4tx1
  def dx5tx1
    return @@dx5tx1
  end
  def dx5tx1=(val)
    @@dx5tx1 = val
  end

  protected :dx5tx1
  def yycon1
    return @@yycon1
  end
  def yycon1=(val)
    @@yycon1 = val
  end

  protected :yycon1
  def yycon2
    return @@yycon2
  end
  def yycon2=(val)
    @@yycon2 = val
  end

  protected :yycon2
  def yycon3
    return @@yycon3
  end
  def yycon3=(val)
    @@yycon3 = val
  end

  protected :yycon3
  def yycon4
    return @@yycon4
  end
  def yycon4=(val)
    @@yycon4 = val
  end

  protected :yycon4
  def yycon5
    return @@yycon5
  end
  def yycon5=(val)
    @@yycon5 = val
  end

  protected :yycon5
  def dy1ty1
    return @@dy1ty1
  end
  def dy1ty1=(val)
    @@dy1ty1 = val
  end

  protected :dy1ty1
  def dy2ty1
    return @@dy2ty1
  end
  def dy2ty1=(val)
    @@dy2ty1 = val
  end

  protected :dy2ty1
  def dy3ty1
    return @@dy3ty1
  end
  def dy3ty1=(val)
    @@dy3ty1 = val
  end

  protected :dy3ty1
  def dy4ty1
    return @@dy4ty1
  end
  def dy4ty1=(val)
    @@dy4ty1 = val
  end

  protected :dy4ty1
  def dy5ty1
    return @@dy5ty1
  end
  def dy5ty1=(val)
    @@dy5ty1 = val
  end

  protected :dy5ty1
  def zzcon1
    return @@zzcon1
  end
  def zzcon1=(val)
    @@zzcon1 = val
  end

  protected :zzcon1
  def zzcon2
    return @@zzcon2
  end
  def zzcon2=(val)
    @@zzcon2 = val
  end

  protected :zzcon2
  def zzcon3
    return @@zzcon3
  end
  def zzcon3=(val)
    @@zzcon3 = val
  end

  protected :zzcon3
  def zzcon4
    return @@zzcon4
  end
  def zzcon4=(val)
    @@zzcon4 = val
  end

  protected :zzcon4
  def zzcon5
    return @@zzcon5
  end
  def zzcon5=(val)
    @@zzcon5 = val
  end

  protected :zzcon5
  def dz1tz1
    return @@dz1tz1
  end
  def dz1tz1=(val)
    @@dz1tz1 = val
  end

  protected :dz1tz1
  def dz2tz1
    return @@dz2tz1
  end
  def dz2tz1=(val)
    @@dz2tz1 = val
  end

  protected :dz2tz1
  def dz3tz1
    return @@dz3tz1
  end
  def dz3tz1=(val)
    @@dz3tz1 = val
  end

  protected :dz3tz1
  def dz4tz1
    return @@dz4tz1
  end
  def dz4tz1=(val)
    @@dz4tz1 = val
  end

  protected :dz4tz1
  def dz5tz1
    return @@dz5tz1
  end
  def dz5tz1=(val)
    @@dz5tz1 = val
  end

  protected :dz5tz1
  def dnxm1
    return @@dnxm1
  end
  def dnxm1=(val)
    @@dnxm1 = val
  end

  protected :dnxm1
  def dnym1
    return @@dnym1
  end
  def dnym1=(val)
    @@dnym1 = val
  end

  protected :dnym1
  def dnzm1
    return @@dnzm1
  end
  def dnzm1=(val)
    @@dnzm1 = val
  end

  protected :dnzm1
  def c1c2
    return @@c1c2
  end
  def c1c2=(val)
    @@c1c2 = val
  end

  protected :c1c2
  def c1c5
    return @@c1c5
  end
  def c1c5=(val)
    @@c1c5 = val
  end

  protected :c1c5
  def c3c4
    return @@c3c4
  end
  def c3c4=(val)
    @@c3c4 = val
  end

  protected :c3c4
  def c1345
    return @@c1345
  end
  def c1345=(val)
    @@c1345 = val
  end

  protected :c1345
  def conz1
    return @@conz1
  end
  def conz1=(val)
    @@conz1 = val
  end

  protected :conz1
  def c1
    return @@c1
  end
  def c1=(val)
    @@c1 = val
  end

  protected :c1
  def c2
    return @@c2
  end
  def c2=(val)
    @@c2 = val
  end

  protected :c2
  def c3
    return @@c3
  end
  def c3=(val)
    @@c3 = val
  end

  protected :c3
  def c4
    return @@c4
  end
  def c4=(val)
    @@c4 = val
  end

  protected :c4
  def c5
    return @@c5
  end
  def c5=(val)
    @@c5 = val
  end

  protected :c5
  def c4dssp
    return @@c4dssp
  end
  def c4dssp=(val)
    @@c4dssp = val
  end

  protected :c4dssp
  def c5dssp
    return @@c5dssp
  end
  def c5dssp=(val)
    @@c5dssp = val
  end

  protected :c5dssp
  def dtdssp
    return @@dtdssp
  end
  def dtdssp=(val)
    @@dtdssp = val
  end

  protected :dtdssp
  def dttx1
    return @@dttx1
  end
  def dttx1=(val)
    @@dttx1 = val
  end

  protected :dttx1
  def bt
    return @@bt
  end
  def bt=(val)
    @@bt = val
  end

  protected :bt
  def dttx2
    return @@dttx2
  end
  def dttx2=(val)
    @@dttx2 = val
  end

  protected :dttx2
  def dtty1
    return @@dtty1
  end
  def dtty1=(val)
    @@dtty1 = val
  end

  protected :dtty1
  def dtty2
    return @@dtty2
  end
  def dtty2=(val)
    @@dtty2 = val
  end

  protected :dtty2
  def dttz1
    return @@dttz1
  end
  def dttz1=(val)
    @@dttz1 = val
  end

  protected :dttz1
  def dttz2
    return @@dttz2
  end
  def dttz2=(val)
    @@dttz2 = val
  end

  protected :dttz2
  def c2dttx1
    return @@c2dttx1
  end
  def c2dttx1=(val)
    @@c2dttx1 = val
  end

  protected :c2dttx1
  def c2dtty1
    return @@c2dtty1
  end
  def c2dtty1=(val)
    @@c2dtty1 = val
  end

  protected :c2dtty1
  def c2dttz1
    return @@c2dttz1
  end
  def c2dttz1=(val)
    @@c2dttz1 = val
  end

  protected :c2dttz1
  def comz1
    return @@comz1
  end
  def comz1=(val)
    @@comz1 = val
  end

  protected :comz1
  def comz4
    return @@comz4
  end
  def comz4=(val)
    @@comz4 = val
  end

  protected :comz4
  def comz5
    return @@comz5
  end
  def comz5=(val)
    @@comz5 = val
  end

  protected :comz5
  def comz6
    return @@comz6
  end
  def comz6=(val)
    @@comz6 = val
  end

  protected :comz6
  def c3c4tx3
    return @@c3c4tx3
  end
  def c3c4tx3=(val)
    @@c3c4tx3 = val
  end

  protected :c3c4tx3
  def c3c4ty3
    return @@c3c4ty3
  end
  def c3c4ty3=(val)
    @@c3c4ty3 = val
  end

  protected :c3c4ty3
  def c3c4tz3
    return @@c3c4tz3
  end
  def c3c4tz3=(val)
    @@c3c4tz3 = val
  end

  protected :c3c4tz3
  def c2iv
    return @@c2iv
  end
  def c2iv=(val)
    @@c2iv = val
  end

  protected :c2iv
  def con43
    return @@con43
  end
  def con43=(val)
    @@con43 = val
  end

  protected :con43
  def con16
    return @@con16
  end
  def con16=(val)
    @@con16 = val
  end
  protected :con16
  attr_accessor :ce
  protected :ce
  attr_accessor :master
  protected :master
  attr_accessor :num_threads
  protected :num_threads
  attr_accessor :rhscomputer
  protected :rhscomputer
  attr_accessor :txinverse
  protected :txinverse
  attr_accessor :xsolver
  protected :xsolver
  attr_accessor :ysolver
  protected :ysolver
  attr_accessor :zsolver
  protected :zsolver
  attr_accessor :rhsadder
  protected :rhsadder
end
