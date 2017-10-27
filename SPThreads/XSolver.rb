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
class XSolver < SPBase
  def initialize(sp, low, high)
    super(sp.clss, sp.num_threads)
    @master =  nil 
    @done = true
    @state = 1
    Init(sp)
    @lower_bound = low
    @upper_bound = high
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = sp
    @lhs = Array.new(5 * (@problem_size + 1), 0.0)
    @lhsp = Array.new(5 * (@problem_size + 1), 0.0)
    @lhsm = Array.new(5 * (@problem_size + 1), 0.0)
    @cv = Array.new(@problem_size, 0.0)
    @rhon = Array.new(@problem_size, 0.0)
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
      for k in @lower_bound..@upper_bound do
          for j in 1..@ny2 do
              for i in 0..@grid_points[0] - 1 do
                  ru1 = @@c3c4 * @rho_i[i + j * @jsize2 + k * @ksize2]
                  @cv[i] = @us[i + j * @jsize2 + k * @ksize2]
                  @rhon[i] = dmax1ex(@@dx2 + @@con43 * ru1, @@dx5 + @@c1c5 * ru1, @@dxmax + ru1, @@dx1)
              end
              lhsinit(@grid_points[0] - 1)
              for i in 1..@nx2 do
                  @lhs[0 + i * @jsize4] = 0.0
                  @lhs[1 + i * @jsize4] = -@@dttx2 * @cv[(i - 1)] - @@dttx1 * @rhon[i - 1]
                  @lhs[2 + i * @jsize4] = 1.0 + @@c2dttx1 * @rhon[i]
                  @lhs[3 + i * @jsize4] = @@dttx2 * @cv[i + 1] - @@dttx1 * @rhon[i + 1]
                  @lhs[4 + i * @jsize4] = 0.0
              end
              i = 1
              @lhs[2 + i * @jsize4] = @lhs[2 + i * @jsize4] + @@comz5
              @lhs[3 + i * @jsize4] = @lhs[3 + i * @jsize4] - @@comz4
              @lhs[4 + i * @jsize4] = @lhs[4 + i * @jsize4] + @@comz1
              @lhs[1 + (i + 1) * @jsize4] = @lhs[1 + (i + 1) * @jsize4] - @@comz4
              @lhs[2 + (i + 1) * @jsize4] = @lhs[2 + (i + 1) * @jsize4] + @@comz6
              @lhs[3 + (i + 1) * @jsize4] = @lhs[3 + (i + 1) * @jsize4] - @@comz4
              @lhs[4 + (i + 1) * @jsize4] = @lhs[4 + (i + 1) * @jsize4] + @@comz1
              for i in 3..@grid_points[0] - 4 do
                  @lhs[0 + i * @jsize4] = @lhs[0 + i * @jsize4] + @@comz1
                  @lhs[1 + i * @jsize4] = @lhs[1 + i * @jsize4] - @@comz4
                  @lhs[2 + i * @jsize4] = @lhs[2 + i * @jsize4] + @@comz6
                  @lhs[3 + i * @jsize4] = @lhs[3 + i * @jsize4] - @@comz4
                  @lhs[4 + i * @jsize4] = @lhs[4 + i * @jsize4] + @@comz1
              end
              i = @grid_points[0] - 3
              @lhs[0 + i * @jsize4] = @lhs[0 + i * @jsize4] + @@comz1
              @lhs[1 + i * @jsize4] = @lhs[1 + i * @jsize4] - @@comz4
              @lhs[2 + i * @jsize4] = @lhs[2 + i * @jsize4] + @@comz6
              @lhs[3 + i * @jsize4] = @lhs[3 + i * @jsize4] - @@comz4
              @lhs[0 + (i + 1) * @jsize4] = @lhs[0 + (i + 1) * @jsize4] + @@comz1
              @lhs[1 + (i + 1) * @jsize4] = @lhs[1 + (i + 1) * @jsize4] - @@comz4
              @lhs[2 + (i + 1) * @jsize4] = @lhs[2 + (i + 1) * @jsize4] + @@comz5
              for i in 1..@nx2 do
                  @lhsp[0 + i * @jsize4] = @lhs[0 + i * @jsize4]
                  @lhsp[1 + i * @jsize4] = @lhs[1 + i * @jsize4] - @@dttx2 * @speed[(i - 1) + j * @jsize2 + k * @ksize2]
                  @lhsp[2 + i * @jsize4] = @lhs[2 + i * @jsize4]
                  @lhsp[3 + i * @jsize4] = @lhs[3 + i * @jsize4] + @@dttx2 * @speed[i + 1 + j * @jsize2 + k * @ksize2]
                  @lhsp[4 + i * @jsize4] = @lhs[4 + i * @jsize4]
                  @lhsm[0 + i * @jsize4] = @lhs[0 + i * @jsize4]
                  @lhsm[1 + i * @jsize4] = @lhs[1 + i * @jsize4] + @@dttx2 * @speed[i - 1 + j * @jsize2 + k * @ksize2]
                  @lhsm[2 + i * @jsize4] = @lhs[2 + i * @jsize4]
                  @lhsm[3 + i * @jsize4] = @lhs[3 + i * @jsize4] - @@dttx2 * @speed[i + 1 + j * @jsize2 + k * @ksize2]
                  @lhsm[4 + i * @jsize4] = @lhs[4 + i * @jsize4]
              end
              for i in 0..@grid_points[0] - 3 do
                  i1 = i + 1
                  i2 = i + 2
                  fac1 = 1.0 / @lhs[2 + i * @jsize4]
                  @lhs[3 + i * @jsize4] = fac1 * @lhs[3 + i * @jsize4]
                  @lhs[4 + i * @jsize4] = fac1 * @lhs[4 + i * @jsize4]
                  for m in 0..2 do
                      @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  end
                  @lhs[2 + i1 * @jsize4] = @lhs[2 + i1 * @jsize4] - @lhs[1 + i1 * @jsize4] * @lhs[3 + i * @jsize4]
                  @lhs[3 + i1 * @jsize4] = @lhs[3 + i1 * @jsize4] - @lhs[1 + i1 * @jsize4] * @lhs[4 + i * @jsize4]
                  for m in 0..2 do
                      @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  end
                  @lhs[1 + i2 * @jsize4] = @lhs[1 + i2 * @jsize4] - @lhs[0 + i2 * @jsize4] * @lhs[3 + i * @jsize4]
                  @lhs[2 + i2 * @jsize4] = @lhs[2 + i2 * @jsize4] - @lhs[0 + i2 * @jsize4] * @lhs[4 + i * @jsize4]
                  for m in 0..2 do
                      @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[0 + i2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  end
              end
              i = @grid_points[0] - 2
              i1 = @grid_points[0] - 1
              fac1 = 1.0 / @lhs[2 + i * @jsize4]
              @lhs[3 + i * @jsize4] = fac1 * @lhs[3 + i * @jsize4]
              @lhs[4 + i * @jsize4] = fac1 * @lhs[4 + i * @jsize4]
              for m in 0..2 do
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              end
              @lhs[2 + i1 * @jsize4] = @lhs[2 + i1 * @jsize4] - @lhs[1 + i1 * @jsize4] * @lhs[3 + i * @jsize4]
              @lhs[3 + i1 * @jsize4] = @lhs[3 + i1 * @jsize4] - @lhs[1 + i1 * @jsize4] * @lhs[4 + i * @jsize4]
              for m in 0..2 do
                  @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              end
              fac2 = 1.0 / @lhs[2 + i1 * @jsize4]
              for m in 0..2 do
                  @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = fac2 * @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1]
              end
              for i in 0..@grid_points[0] - 3 do
                  i1 = i + 1
                  i2 = i + 2
                  m = 3
                  fac1 = 1.0 / @lhsp[2 + i * @jsize4]
                  @lhsp[3 + i * @jsize4] = fac1 * @lhsp[3 + i * @jsize4]
                  @lhsp[4 + i * @jsize4] = fac1 * @lhsp[4 + i * @jsize4]
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  @lhsp[2 + i1 * @jsize4] = @lhsp[2 + i1 * @jsize4] - @lhsp[1 + i1 * @jsize4] * @lhsp[3 + i * @jsize4]
                  @lhsp[3 + i1 * @jsize4] = @lhsp[3 + i1 * @jsize4] - @lhsp[1 + i1 * @jsize4] * @lhsp[4 + i * @jsize4]
                  @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  @lhsp[1 + i2 * @jsize4] = @lhsp[1 + i2 * @jsize4] - @lhsp[0 + i2 * @jsize4] * @lhsp[3 + i * @jsize4]
                  @lhsp[2 + i2 * @jsize4] = @lhsp[2 + i2 * @jsize4] - @lhsp[0 + i2 * @jsize4] * @lhsp[4 + i * @jsize4]
                  @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[0 + i2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  m = 4
                  fac1 = 1.0 / @lhsm[2 + i * @jsize4]
                  @lhsm[3 + i * @jsize4] = fac1 * @lhsm[3 + i * @jsize4]
                  @lhsm[4 + i * @jsize4] = fac1 * @lhsm[4 + i * @jsize4]
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  @lhsm[2 + i1 * @jsize4] = @lhsm[2 + i1 * @jsize4] - @lhsm[1 + i1 * @jsize4] * @lhsm[3 + i * @jsize4]
                  @lhsm[3 + i1 * @jsize4] = @lhsm[3 + i1 * @jsize4] - @lhsm[1 + i1 * @jsize4] * @lhsm[4 + i * @jsize4]
                  @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  @lhsm[1 + i2 * @jsize4] = @lhsm[1 + i2 * @jsize4] - @lhsm[0 + i2 * @jsize4] * @lhsm[3 + i * @jsize4]
                  @lhsm[2 + i2 * @jsize4] = @lhsm[2 + i2 * @jsize4] - @lhsm[0 + i2 * @jsize4] * @lhsm[4 + i * @jsize4]
                  @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[0 + i2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              end
              i = @grid_points[0] - 2
              i1 = @grid_points[0] - 1
              m = 3
              fac1 = 1.0 / @lhsp[2 + i * @jsize4]
              @lhsp[3 + i * @jsize4] = fac1 * @lhsp[3 + i * @jsize4]
              @lhsp[4 + i * @jsize4] = fac1 * @lhsp[4 + i * @jsize4]
              @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              @lhsp[2 + i1 * @jsize4] = @lhsp[2 + i1 * @jsize4] - @lhsp[1 + i1 * @jsize4] * @lhsp[3 + i * @jsize4]
              @lhsp[3 + i1 * @jsize4] = @lhsp[3 + i1 * @jsize4] - @lhsp[1 + i1 * @jsize4] * @lhsp[4 + i * @jsize4]
              @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              m = 4
              fac1 = 1.0 / @lhsm[2 + i * @jsize4]
              @lhsm[3 + i * @jsize4] = fac1 * @lhsm[3 + i * @jsize4]
              @lhsm[4 + i * @jsize4] = fac1 * @lhsm[4 + i * @jsize4]
              @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              @lhsm[2 + i1 * @jsize4] = @lhsm[2 + i1 * @jsize4] - @lhsm[1 + i1 * @jsize4] * @lhsm[3 + i * @jsize4]
              @lhsm[3 + i1 * @jsize4] = @lhsm[3 + i1 * @jsize4] - @lhsm[1 + i1 * @jsize4] * @lhsm[4 + i * @jsize4]
              @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              @rhs[3 + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i1 * @isize1 + j * @jsize1 + k * @ksize1] / @lhsp[2 + i1 * @jsize4]
              @rhs[4 + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i1 * @isize1 + j * @jsize1 + k * @ksize1] / @lhsm[2 + i1 * @jsize4]
              i = @grid_points[0] - 2
              i1 = @grid_points[0] - 1
              for m in 0..2 do
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[3 + i * @jsize4] * @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1]
              end
              @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[3 + i * @jsize4] * @rhs[3 + i1 * @isize1 + j * @jsize1 + k * @ksize1]
              @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[3 + i * @jsize4] * @rhs[4 + i1 * @isize1 + j * @jsize1 + k * @ksize1]
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
              i = @grid_points[0] - 3
              while i >= 0 do
                  i1 = i + 1
                  i2 = i + 2
                  for m in 0..2 do
                      @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[3 + i * @jsize4] * @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[4 + i * @jsize4] * @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1]
                  end
                  @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[3 + i * @jsize4] * @rhs[3 + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[4 + i * @jsize4] * @rhs[3 + i2 * @isize1 + j * @jsize1 + k * @ksize1]
                  @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[3 + i * @jsize4] * @rhs[4 + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[4 + i * @jsize4] * @rhs[4 + i2 * @isize1 + j * @jsize1 + k * @ksize1]
                i -= 1
              end
          end
      end
    when 2
      for k in @lower_bound..@upper_bound do
          for j in 1..@ny2 do
              for i in 1..@nx2 do
                  r1 = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  r2 = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  r3 = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  r4 = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  r5 = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  t1 = @@bt * r3
                  t2 = 0.5 * (r4 + r5)
                  @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = -r2
                  @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = r1
                  @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @@bt * (r4 - r5)
                  @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = -t1 + t2
                  @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = t1 + t2
              end
          end
      end
    end
    @state += 1
    if @state == 3 then
      @state = 1
    else
    end
  end

  def lhsinit(size)
    
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    i = 0
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
    i = 0
    while i <= size do
        @lhs[2 + i * @jsize4] = 1.0
        @lhsp[2 + i * @jsize4] = 1.0
        @lhsm[2 + i * @jsize4] = 1.0
      i += size
    end
  end

# *** public ***

  attr_accessor :id

  attr_accessor :done

end
