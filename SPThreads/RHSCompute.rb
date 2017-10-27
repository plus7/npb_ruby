# -*- coding: utf-8 -*-
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
class RHSCompute < SPBase
  def initialize(sp, low1, high1, low2, high2)
    super(sp.clss, sp.num_threads)
    @master =  nil 
    @done = true
    Init(sp)
    @lower_bound1 = low1
    @upper_bound1 = high1
    @lower_bound2 = low2
    @upper_bound2 = high2
    @state = 1
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = sp
  end

  def Init(sp)
    @IMAX = sp.imax
    @JMAX = sp.jmax
    @KMAX = sp.kmax
    @problem_size = sp.problem_size
    @nx2 = sp.nx2
    @ny2 = sp.ny2
    @nz2 = sp.nz2
    @grid_points = sp.grid_points
    @niter_default = sp.niter_default
    @dt_default = sp.dt_default
    @u = sp.u
    @rhs = sp.rhs
    @forcing = sp.forcing
    @isize1 = sp.isize1
    @jsize1 = sp.jsize1
    @ksize1 = sp.ksize1
    @us = sp.us
    @vs = sp.vs
    @ws = sp.ws
    @qs = sp.qs
    @rho_i = sp.rho_i
    @speed = sp.speed
    @square = sp.square
    @jsize2 = sp.jsize2
    @ksize2 = sp.ksize2
    @ue = sp.ue
    @buf = sp.buf
    @jsize3 = sp.jsize3
    @lhs = sp.lhs
    @lhsp = sp.lhsp
    @lhsm = sp.lhsm
    @jsize4 = sp.jsize4
    @cv = sp.cv
    @rhon = sp.rhon
    @rhos = sp.rhos
    @rhoq = sp.rhoq
    @cuf = sp.cuf
    @q = sp.q
    @ce = sp.ce
  end

  def run()
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    while true do
         self.synchronize do
            while @done == true do
                begin
                    cond.wait
                    @master.synchronize do
                        @master.cond.signal
                    end
                rescue InterruptedException => ie
                ensure
                end
            end
            step()
            @master.synchronize do
                @done = true
                @master.cond.signal
            end
        end
    end
  end

  def step()
    
    
    case @state
    when 1
      for k in @lower_bound1..@upper_bound1 do
        for j in 0..@grid_points[1] - 1 do
          for i in 0..@grid_points[0] - 1 do
            rho_inv = 1.0 / @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @rho_i[i + j * @jsize2 + k * @ksize2] = rho_inv
            @us[i + j * @jsize2 + k * @ksize2] = @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * rho_inv
            @vs[i + j * @jsize2 + k * @ksize2] = @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * rho_inv
            @ws[i + j * @jsize2 + k * @ksize2] = @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * rho_inv
            @square[i + j * @jsize2 + k * @ksize2] = 0.5 * (@u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]) * rho_inv
            @qs[i + j * @jsize2 + k * @ksize2] = @square[i + j * @jsize2 + k * @ksize2] * rho_inv
            aux = @@c1c2 * rho_inv * (@u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @square[i + j * @jsize2 + k * @ksize2])
            @speed[i + j * @jsize2 + k * @ksize2] = Math.sqrt(aux)
          end
        end
      end
    when 2
      for k in @lower_bound1..@upper_bound1 do
        for j in 0..@grid_points[1] - 1 do
          for i in 0..@grid_points[0] - 1 do
            for m in 0..4 do
              @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            end
          end
        end
      end
    when 3
      for k in @lower_bound2..@upper_bound2 do
        for j in 1..@ny2 do
          for i in 1..@nx2 do
            uijk = @us[i + j * @jsize2 + k * @ksize2]
            up1 = @us[i + 1 + j * @jsize2 + k * @ksize2]
            um1 = @us[i - 1 + j * @jsize2 + k * @ksize2]
            @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx1tx1 * (@u[0 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) - @@tx2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] -  @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx2tx1 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) + @@xxcon2 * @@con43 * (up1 - 2.0 * uijk + um1) - @@tx2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * up1 - @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * um1 + (@u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - @square[(i + 1) + j * @jsize2 + k * @ksize2] - @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + @square[(i - 1) + j * @jsize2 + k * @ksize2]) * @@c2)
            @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx3tx1 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) + @@xxcon2 * (@vs[i + 1 + j * @jsize2 + k * @ksize2] - 2.0 * @vs[i + j * @jsize2 + k * @ksize2] + @vs[i - 1 + j * @jsize2 + k * @ksize2]) - @@tx2 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * up1 - @u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * um1)
            @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx4tx1 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) + @@xxcon2 * (@ws[i + 1 + j * @jsize2 + k * @ksize2] - 2.0 * @ws[i + j * @jsize2 + k * @ksize2] + @ws[i - 1 + j * @jsize2 + k * @ksize2]) - @@tx2 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * up1 - @u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * um1)
            @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx5tx1 * (@u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) + @@xxcon3 * (@qs[i + 1 + j * @jsize2 + k * @ksize2] - 2.0 * @qs[i + j * @jsize2 + k * @ksize2] + @qs[i - 1 + j * @jsize2 + k * @ksize2]) + @@xxcon4 * (up1 * up1 - 2.0 * uijk * uijk + um1 * um1) + @@xxcon5 * (@u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + 1 + j * @jsize2 + k * @ksize2] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize2 + k * @ksize2] + @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i - 1 + j * @jsize2 + k * @ksize2]) - @@tx2 * ((@@c1 * @u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - @@c2 * @square[i + 1 + j * @jsize2 + k * @ksize2]) * up1 - (@@c1 * @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - @@c2 * @square[i - 1 + j * @jsize2 + k * @ksize2]) * um1)
            
          end
          i = 1
          for m in 0..4 do
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + (i + 2) * @isize1 + j * @jsize1 + k * @ksize1])
          end
          i = 2
          for m in 0..4 do
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + (i + 2) * @isize1 + j * @jsize1 + k * @ksize1])
          end
          for i in 3..@nx2 - 2 do
            for m in 0..4 do
              @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + (i - 2) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + (i + 2) * @isize1 + j * @jsize1 + k * @ksize1])
            end
          end
          i = @nx2 - 1
          for m in 0..4 do
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + (i - 2) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
          end
          i = @nx2
          for m in 0..4 do
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + (i - 2) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1])
          end
        end
      end
    when 4
      for k in @lower_bound2..@upper_bound2 do
        for j in 1..@ny2 do
          for i in 1..@nx2 do
            vijk = @vs[i + j * @jsize2 + k * @ksize2]
            vp1 = @vs[i + (j + 1) * @jsize2 + k * @ksize2]
            vm1 = @vs[i + (j - 1) * @jsize2 + k * @ksize2]
            @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy1ty1 * (@u[0 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) - @@ty2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1])
            @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy2ty1 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon2 * (@us[i + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @us[i + j * @jsize2 + k * @ksize2] + @us[i + (j - 1) * @jsize2 + k * @ksize2]) - @@ty2 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * vp1 - @u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * vm1)
            @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy3ty1 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon2 * @@con43 * (vp1 - 2.0 * vijk + vm1) - @@ty2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * vp1 - @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * vm1 + (@u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - @square[i + (j + 1) * @jsize2 + k * @ksize2] - @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + @square[i + (j - 1) * @jsize2 + k * @ksize2]) * @@c2)
            @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy4ty1 * (@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon2 * (@ws[i + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @ws[i + j * @jsize2 + k * @ksize2] + @ws[i + (j - 1) * @jsize2 + k * @ksize2]) - @@ty2 * (@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * vp1 - @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * vm1)
            @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy5ty1 * (@u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon3 * (@qs[i + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @qs[i + j * @jsize2 + k * @ksize2] + @qs[i + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon4 * (vp1 * vp1 - 2.0 * vijk * vijk + vm1 * vm1) + @@yycon5 * (@u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @rho_i[i + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize2 + k * @ksize2] + @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * @rho_i[i + (j - 1) * @jsize2 + k * @ksize2]) - @@ty2 * ((@@c1 * @u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - @@c2 * @square[i + (j + 1) * @jsize2 + k * @ksize2]) * vp1 - (@@c1 * @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - @@c2 * @square[i + (j - 1) * @jsize2 + k * @ksize2]) * vm1)
          end
        end
        j = 1
        for i in 1..@nx2 do
          for m in 0..4 do
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + (j + 2) * @jsize1 + k * @ksize1])
          end
        end
        j = 2
        for i in 1..@nx2 do
          for m in 0..4 do
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + (j + 2) * @jsize1 + k * @ksize1])
          end
        end
        for j in 3..@ny2 - 2 do
          for i in 1..@nx2 do
            for m in 0..4 do
              @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + (j - 2) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + (j + 2) * @jsize1 + k * @ksize1])
            end
          end
        end
        j = @ny2 - 1
        for i in 1..@nx2 do
          for m in 0..4 do
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + (j - 2) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
          end
        end
        j = @ny2
        for i in 1..@nx2 do
          for m in 0..4 do
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + (j - 2) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1])
          end
        end
      end
    when 5
      for j in @lower_bound2..@upper_bound2 do
          for k in 1..@ny2 do
              for i in 1..@nx2 do
                  wijk = @ws[i + j * @jsize2 + k * @ksize2]
                  wp1 = @ws[i + j * @jsize2 + (k + 1) * @ksize2]
                  wm1 = @ws[i + j * @jsize2 + (k - 1) * @ksize2]
                  @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz1tz1 * (@u[0 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) - @@tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
                  @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz2tz1 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon2 * (@us[i + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @us[i + j * @jsize2 + k * @ksize2] + @us[i + j * @jsize2 + (k - 1) * @ksize2]) - @@tz2 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * wp1 - @u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * wm1)
                  @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz3tz1 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon2 * (@vs[i + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @vs[i + j * @jsize2 + k * @ksize2] + @vs[i + j * @jsize2 + (k - 1) * @ksize2]) - @@tz2 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * wp1 - @u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * wm1)
                  @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz4tz1 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon2 * @@con43 * (wp1 - 2.0 * wijk + wm1) - @@tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * wp1 - @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * wm1 + (@u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - @square[i + j * @jsize2 + (k + 1) * @ksize2] - @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + @square[i + j * @jsize2 + (k - 1) * @ksize2]) * @@c2)
                  @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz5tz1 * (@u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon3 * (@qs[i + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @qs[i + j * @jsize2 + k * @ksize2] + @qs[i + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon4 * (wp1 * wp1 - 2.0 * wijk * wijk + wm1 * wm1) + @@zzcon5 * (@u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @rho_i[i + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize2 + k * @ksize2] + @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * @rho_i[i + j * @jsize2 + (k - 1) * @ksize2]) - @@tz2 * ((@@c1 * @u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - @@c2 * @square[i + j * @jsize2 + (k + 1) * @ksize2]) * wp1 - (@@c1 * @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] - @@c2 * @square[i + j * @jsize2 + (k - 1) * @ksize2]) * wm1)
              end
          end
          k = 1
          for i in 1..@nx2 do
              for m in 0..4 do
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + (k + 2) * @ksize1])
              end
          end
          k = 2
          for i in 1..@nx2 do
              for m in 0..4 do
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + (k + 2) * @ksize1])
              end
          end
          for k in 3..@nz2 - 2 do
              for i in 1..@nx2 do
                  for m in 0..4 do
                      @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + j * @jsize1 + (k - 2) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + (k + 2) * @ksize1])
                  end
              end
          end
          k = @nz2 - 1
          for i in 1..@nx2 do
              for m in 0..4 do
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + j * @jsize1 + (k - 2) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
              end
          end
          k = @nz2
          for i in 1..@nx2 do
              for m in 0..4 do
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + j * @jsize1 + (k - 2) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1])
              end
          end
      end
    when 6
      for k in @lower_bound2..@upper_bound2 do
          for j in 1..@ny2 do
              for i in 1..@nx2 do
                  for m in 0..4 do
                      @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] * @@dt
                  end
              end
          end
      end
    end
    @state += 1
    if @state == 7 then
      @state = 1
    else
    end
  end

# *** public ***
  attr_accessor :id

  attr_accessor :done
end
