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
class RHSCompute < BTBase
  def initialize(bt, low1, high1, low2, high2)
    super(bt.clss, bt.num_threads)
    @master =  nil 
    @done = true
    Init(bt)
    @lower_bound1 = low1
    @upper_bound1 = high1
    @lower_bound2 = low2
    @upper_bound2 = high2
    @state = 1
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = bt
  end

  def Init(bt)
    @IMAX = bt.imax
    @JMAX = bt.jmax
    @KMAX = bt.kmax
    @problem_size = bt.problem_size
    @grid_points = bt.grid_points
    @niter_default = bt.niter_default
    @dt_default = bt.dt_default
    @u = bt.u
    @rhs = bt.rhs
    @forcing = bt.forcing
    @cv = bt.cv
    @q = bt.q
    @cuf = bt.cuf
    @isize2 = bt.isize2
    @jsize2 = bt.jsize2
    @ksize2 = bt.ksize2
    @us = bt.us
    @vs = bt.vs
    @ws = bt.ws
    @qs = bt.qs
    @rho_i = bt.rho_i
    @square = bt.square
    @jsize1 = bt.jsize1
    @ksize1 = bt.ksize1
    @ue = bt.ue
    @buf = bt.buf
    @jsize3 = bt.jsize3
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
            @rho_inv = 1.0 / @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2]
            @rho_i[i + j * @jsize1 + k * @ksize1] = @rho_inv
            @us[i + j * @jsize1 + k * @ksize1] = @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_inv
            @vs[i + j * @jsize1 + k * @ksize1] = @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_inv
            @ws[i + j * @jsize1 + k * @ksize1] = @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_inv
            @square[i + j * @jsize1 + k * @ksize1] = 0.5 * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @rho_inv
            @qs[i + j * @jsize1 + k * @ksize1] = @square[i + j * @jsize1 + k * @ksize1] * @rho_inv
          end
        end
      end
      for k in @lower_bound1..@upper_bound1 do
        for j in 0..@grid_points[1] - 1 do
          for i in 0..@grid_points[0] - 1 do
            for m in 0..4 do
              @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2]
            end
          end
        end
      end
    when 2
      for k in @lower_bound2..@upper_bound2 do
        for j in 1..@grid_points[1] - 2 do
          for i in 1..@grid_points[0] - 2 do
            @uijk = @us[i + j * @jsize1 + k * @ksize1]
            @up1 = @us[(i + 1) + j * @jsize1 + k * @ksize1]
            @um1 = @us[(i - 1) + j * @jsize1 + k * @ksize1]
            @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx1tx1 * (@u[0 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[0 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) - @@tx2 * (@u[1 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - @u[1 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2])
            @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx2tx1 * (@u[1 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[1 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) + @@xxcon2 * @@con43 * (@up1 - 2.0 * @uijk + @um1) - @@tx2 * (@u[1 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] * @up1 - @u[1 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] * @um1 + (@u[4 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - @square[(i + 1) + j * @jsize1 + k * @ksize1] - @u[4 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + @square[(i - 1) + j * @jsize1 + k * @ksize1]) * @@c2)
            @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx3tx1 * (@u[2 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[2 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) + @@xxcon2 * (@vs[(i + 1) + j * @jsize1 + k * @ksize1] - 2.0 * @vs[i + j * @jsize1 + k * @ksize1] + @vs[(i - 1) + j * @jsize1 + k * @ksize1]) - @@tx2 * (@u[2 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] * @up1 - @u[2 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] * @um1)
            @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx4tx1 * (@u[3 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[3 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) + @@xxcon2 * (@ws[(i + 1) + j * @jsize1 + k * @ksize1] - 2.0 * @ws[i + j * @jsize1 + k * @ksize1] + @ws[(i - 1) + j * @jsize1 + k * @ksize1]) - @@tx2 * (@u[3 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] * @up1 - @u[3 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] * @um1)
            @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx5tx1 * (@u[4 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[4 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) + @@xxcon3 * (@qs[(i + 1) + j * @jsize1 + k * @ksize1] - 2.0 * @qs[i + j * @jsize1 + k * @ksize1] + @qs[(i - 1) + j * @jsize1 + k * @ksize1]) + @@xxcon4 * (@up1 * @up1 - 2.0 * @uijk * @uijk + @um1 * @um1) + @@xxcon5 * (@u[4 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[(i + 1) + j * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[i + j * @jsize1 + k * @ksize1] + @u[4 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[(i - 1) + j * @jsize1 + k * @ksize1]) - @@tx2 * ((@@c1 * @u[4 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - @@c2 * @square[(i + 1) + j * @jsize1 + k * @ksize1]) * @up1 - (@@c1 * @u[4 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] - @@c2 * @square[(i - 1) + j * @jsize1 + k * @ksize1]) * @um1)
          end
        end
        for j in 1..@grid_points[1] - 2 do
          i = 1
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] + @u[m + (i + 2) * @isize2 + j * @jsize2 + k * @ksize2])
          end
          i = 2
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @u[m + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] + @u[m + (i + 2) * @isize2 + j * @jsize2 + k * @ksize2])
          end
          for m in 0..4 do
            for i in 3..@grid_points[0] - 4 do
              @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + (i - 2) * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] + @u[m + (i + 2) * @isize2 + j * @jsize2 + k * @ksize2])
            end
          end
          i = @grid_points[0] - 3
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + (i - 2) * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2])
          end
          i = @grid_points[0] - 2
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + (i - 2) * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + 5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2])
          end
        end
      end
    when 3
      for k in @lower_bound2..@upper_bound2 do
        for j in 1..@grid_points[1] - 2 do
          for i in 1..@grid_points[0] - 2 do
            @vijk = @vs[i + j * @jsize1 + k * @ksize1]
            @vp1 = @vs[i + (j + 1) * @jsize1 + k * @ksize1]
            @vm1 = @vs[i + (j - 1) * @jsize1 + k * @ksize1]
            @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy1ty1 * (@u[0 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[0 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) - @@ty2 * (@u[2 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - @u[2 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2])
            @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy2ty1 * (@u[1 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[1 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon2 * (@us[i + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @us[i + j * @jsize1 + k * @ksize1] + @us[i + (j - 1) * @jsize1 + k * @ksize1]) - @@ty2 * (@u[1 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] * @vp1 - @u[1 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] * @vm1)
            @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy3ty1 * (@u[2 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[2 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon2 * @@con43 * (@vp1 - 2.0 * @vijk + @vm1) - @@ty2 * (@u[2 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] * @vp1 - @u[2 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] * @vm1 + (@u[4 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - @square[i + (j + 1) * @jsize1 + k * @ksize1] - @u[4 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + @square[i + (j - 1) * @jsize1 + k * @ksize1]) * @@c2)
            @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy4ty1 * (@u[3 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[3 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon2 * (@ws[i + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @ws[i + j * @jsize1 + k * @ksize1] + @ws[i + (j - 1) * @jsize1 + k * @ksize1]) - @@ty2 * (@u[3 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] * @vp1 - @u[3 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] * @vm1)
            @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy5ty1 * (@u[4 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[4 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon3 * (@qs[i + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @qs[i + j * @jsize1 + k * @ksize1] + @qs[i + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon4 * (@vp1 * @vp1 - 2.0 * @vijk * @vijk + @vm1 * @vm1) + @@yycon5 * (@u[4 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] * @rho_i[i + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[i + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] * @rho_i[i + (j - 1) * @jsize1 + k * @ksize1]) - @@ty2 * ((@@c1 * @u[4 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - @@c2 * @square[i + (j + 1) * @jsize1 + k * @ksize1]) * @vp1 - (@@c1 * @u[4 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] - @@c2 * @square[i + (j - 1) * @jsize1 + k * @ksize1]) * @vm1)
          end
        end
        for i in 1..@grid_points[0] - 2 do
          j = 1
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] + @u[m + i * @isize2 + (j + 2) * @jsize2 + k * @ksize2])
          end
          j = 2
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @u[m + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] + @u[m + i * @isize2 + (j + 2) * @jsize2 + k * @ksize2])
          end
        end
        for j in 3..@grid_points[1] - 4 do
          for i in 1..@grid_points[0] - 2 do
            for m in 0..4 do
              @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + (j - 2) * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] + @u[m + i * @isize2 + (j + 2) * @jsize2 + k * @ksize2])
            end
          end
        end
        for i in 1..@grid_points[0] - 2 do
          j = @grid_points[1] - 3
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + (j - 2) * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2])
          end
          j = @grid_points[1] - 2
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + (j - 2) * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + 5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2])
          end
        end
      end
    when 4
      for j in @lower_bound2..@upper_bound2 do
        for k in 1..@grid_points[1] - 2 do
          for i in 1..@grid_points[0] - 2 do
            @wijk = @ws[i + j * @jsize1 + k * @ksize1]
            @wp1 = @ws[i + j * @jsize1 + (k + 1) * @ksize1]
            @wm1 = @ws[i + j * @jsize1 + (k - 1) * @ksize1]
            @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz1tz1 * (@u[0 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[0 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) - @@tz2 * (@u[3 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - @u[3 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2])
            @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz2tz1 * (@u[1 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[1 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon2 * (@us[i + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @us[i + j * @jsize1 + k * @ksize1] + @us[i + j * @jsize1 + (k - 1) * @ksize1]) - @@tz2 * (@u[1 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] * @wp1 - @u[1 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] * @wm1)
            @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz3tz1 * (@u[2 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[2 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon2 * (@vs[i + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @vs[i + j * @jsize1 + k * @ksize1] + @vs[i + j * @jsize1 + (k - 1) * @ksize1]) - @@tz2 * (@u[2 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] * @wp1 - @u[2 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] * @wm1)
            @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz4tz1 * (@u[3 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[3 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon2 * @@con43 * (@wp1 - 2.0 * @wijk + @wm1) - @@tz2 * (@u[3 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] * @wp1 - @u[3 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] * @wm1 + (@u[4 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - @square[i + j * @jsize1 + (k + 1) * @ksize1] - @u[4 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + @square[i + j * @jsize1 + (k - 1) * @ksize1]) * @@c2)
            @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz5tz1 * (@u[4 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[4 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon3 * (@qs[i + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @qs[i + j * @jsize1 + k * @ksize1] + @qs[i + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon4 * (@wp1 * @wp1 - 2.0 * @wijk * @wijk + @wm1 * @wm1) + @@zzcon5 * (@u[4 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] * @rho_i[i + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[i + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] * @rho_i[i + j * @jsize1 + (k - 1) * @ksize1]) - @@tz2 * ((@@c1 * @u[4 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - @@c2 * @square[i + j * @jsize1 + (k + 1) * @ksize1]) * @wp1 - (@@c1 * @u[4 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] - @@c2 * @square[i + j * @jsize1 + (k - 1) * @ksize1]) * @wm1)
          end
        end
        for i in 1..@grid_points[0] - 2 do
          k = 1
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] + @u[m + i * @isize2 + j * @jsize2 + (k + 2) * @ksize2])
          end
          k = 2
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @u[m + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] + @u[m + i * @isize2 + j * @jsize2 + (k + 2) * @ksize2])
          end
        end
        for k in 3..@grid_points[2] - 4 do
          for i in 1..@grid_points[0] - 2 do
            for m in 0..4 do
              @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + j * @jsize2 + (k - 2) * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] + @u[m + i * @isize2 + j * @jsize2 + (k + 2) * @ksize2])
            end
          end
        end
        for i in 1..@grid_points[0] - 2 do
          k = @grid_points[2] - 3
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + j * @jsize2 + (k - 2) * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2])
          end
          k = @grid_points[2] - 2
          for m in 0..4 do
            @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + j * @jsize2 + (k - 2) * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + 5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2])
          end
        end
      end
    when 5
      for k in @lower_bound2..@upper_bound2 do
        for j in 1..@grid_points[1] - 2 do
          for i in 1..@grid_points[0] - 2 do
            for m in 0..4 do
              @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] * @@dt
            end
          end
        end
      end
    end
    @state += 1
    if @state == 6 then
      @state = 1
    else
    end
  end

  def print_arrays()
    rhs_density = 0; rhs_x_momentum = 0; rhs_y_momentum = 0; rhs_z_momentum = 0; rhs_energy = 0; 
    for i in 0..@grid_points[2]-1 do
      for j in 0..@grid_points[1]-1 do
        for k in 0..@grid_points[0]-1 do
          rhs_density += @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2]
          rhs_x_momentum += @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
          rhs_y_momentum += @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
          rhs_z_momentum += @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
          rhs_energy += @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2]
        end
      end
    end
    puts(" density: " + rhs_density.to_s)
    puts(" x_momentum: " + rhs_x_momentum.to_s)
    puts(" y_momentum: " + rhs_y_momentum.to_s)
    puts(" z_momentum: " + rhs_z_momentum.to_s)
    puts(" energy: " + rhs_energy.to_s)
  end

  # *** public ***
  attr_accessor :id

  attr_accessor :done
end
