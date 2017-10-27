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
class YSolver < SPBase
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
    @rhoq = Array.new(@problem_size, 0.0)
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
         self .synchronize do
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
          for i in 1..@grid_points[0] - 2 do
              for j in 0..@grid_points[1] - 1 do
                  ru1 = @@c3c4 * @rho_i[i + j * @jsize2 + k * @ksize2]
                  @cv[j] = @vs[i + j * @jsize2 + k * @ksize2]
                  @rhoq[j] = dmax1ex(@@dy3 + @@con43 * ru1, @@dy5 + @@c1c5 * ru1, @@dymax + ru1, @@dy1)
              end
              lhsinit(@grid_points[1] - 1)
              for j in 1..@grid_points[1] - 2 do
                  @lhs[0 + j * @jsize4] = 0.0
                  @lhs[1 + j * @jsize4] = -@@dtty2 * @cv[j - 1] - @@dtty1 * @rhoq[j - 1]
                  @lhs[2 + j * @jsize4] = 1.0 + @@c2dtty1 * @rhoq[j]
                  @lhs[3 + j * @jsize4] = @@dtty2 * @cv[j + 1] - @@dtty1 * @rhoq[j + 1]
                  @lhs[4 + j * @jsize4] = 0.0
              end
              j = 1
              @lhs[2 + j * @jsize4] = @lhs[2 + j * @jsize4] + @@comz5
              @lhs[3 + j * @jsize4] = @lhs[3 + j * @jsize4] - @@comz4
              @lhs[4 + j * @jsize4] = @lhs[4 + j * @jsize4] + @@comz1
              @lhs[1 + (j + 1) * @jsize4] = @lhs[1 + (j + 1) * @jsize4] - @@comz4
              @lhs[2 + (j + 1) * @jsize4] = @lhs[2 + (j + 1) * @jsize4] + @@comz6
              @lhs[3 + (j + 1) * @jsize4] = @lhs[3 + (j + 1) * @jsize4] - @@comz4
              @lhs[4 + (j + 1) * @jsize4] = @lhs[4 + (j + 1) * @jsize4] + @@comz1
              for j in 3..@grid_points[1] - 4 do
                  @lhs[0 + j * @jsize4] = @lhs[0 + j * @jsize4] + @@comz1
                  @lhs[1 + j * @jsize4] = @lhs[1 + j * @jsize4] - @@comz4
                  @lhs[2 + j * @jsize4] = @lhs[2 + j * @jsize4] + @@comz6
                  @lhs[3 + j * @jsize4] = @lhs[3 + j * @jsize4] - @@comz4
                  @lhs[4 + j * @jsize4] = @lhs[4 + j * @jsize4] + @@comz1
              end
              j = @grid_points[1] - 3
              @lhs[0 + j * @jsize4] = @lhs[0 + j * @jsize4] + @@comz1
              @lhs[1 + j * @jsize4] = @lhs[1 + j * @jsize4] - @@comz4
              @lhs[2 + j * @jsize4] = @lhs[2 + j * @jsize4] + @@comz6
              @lhs[3 + j * @jsize4] = @lhs[3 + j * @jsize4] - @@comz4
              @lhs[0 + (j + 1) * @jsize4] = @lhs[0 + (j + 1) * @jsize4] + @@comz1
              @lhs[1 + (j + 1) * @jsize4] = @lhs[1 + (j + 1) * @jsize4] - @@comz4
              @lhs[2 + (j + 1) * @jsize4] = @lhs[2 + (j + 1) * @jsize4] + @@comz5
              for j in 1..@grid_points[1] - 2 do
                  @lhsp[0 + j * @jsize4] = @lhs[0 + j * @jsize4]
                  @lhsp[1 + j * @jsize4] = @lhs[1 + j * @jsize4] - @@dtty2 * @speed[i + (j - 1) * @jsize2 + k * @ksize2]
                  @lhsp[2 + j * @jsize4] = @lhs[2 + j * @jsize4]
                  @lhsp[3 + j * @jsize4] = @lhs[3 + j * @jsize4] + @@dtty2 * @speed[i + (j + 1) * @jsize2 + k * @ksize2]
                  @lhsp[4 + j * @jsize4] = @lhs[4 + j * @jsize4]
                  @lhsm[0 + j * @jsize4] = @lhs[0 + j * @jsize4]
                  @lhsm[1 + j * @jsize4] = @lhs[1 + j * @jsize4] + @@dtty2 * @speed[i + (j - 1) * @jsize2 + k * @ksize2]
                  @lhsm[2 + j * @jsize4] = @lhs[2 + j * @jsize4]
                  @lhsm[3 + j * @jsize4] = @lhs[3 + j * @jsize4] - @@dtty2 * @speed[i + (j + 1) * @jsize2 + k * @ksize2]
                  @lhsm[4 + j * @jsize4] = @lhs[4 + j * @jsize4]
              end
              for j in 0..@grid_points[1] - 3 do
                  j1 = j + 1
                  j2 = j + 2
                  fac1 = 1.0 / @lhs[2 + j * @jsize4]
                  @lhs[3 + j * @jsize4] = fac1 * @lhs[3 + j * @jsize4]
                  @lhs[4 + j * @jsize4] = fac1 * @lhs[4 + j * @jsize4]
                  for m in 0..2 do
                      @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  end
                  @lhs[2 + j1 * @jsize4] = @lhs[2 + j1 * @jsize4] - @lhs[1 + j1 * @jsize4] * @lhs[3 + j * @jsize4]
                  @lhs[3 + j1 * @jsize4] = @lhs[3 + j1 * @jsize4] - @lhs[1 + j1 * @jsize4] * @lhs[4 + j * @jsize4]
                  for m in 0..2 do
                      @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhs[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  end
                  @lhs[1 + j2 * @jsize4] = @lhs[1 + j2 * @jsize4] - @lhs[0 + j2 * @jsize4] * @lhs[3 + j * @jsize4]
                  @lhs[2 + j2 * @jsize4] = @lhs[2 + j2 * @jsize4] - @lhs[0 + j2 * @jsize4] * @lhs[4 + j * @jsize4]
                  for m in 0..2 do
                      @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] - @lhs[0 + j2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  end
              end
              j = @grid_points[1] - 2
              j1 = @grid_points[1] - 1
              fac1 = 1.0 / @lhs[2 + j * @jsize4]
              @lhs[3 + j * @jsize4] = fac1 * @lhs[3 + j * @jsize4]
              @lhs[4 + j * @jsize4] = fac1 * @lhs[4 + j * @jsize4]
              for m in 0..2 do
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              end
              @lhs[2 + j1 * @jsize4] = @lhs[2 + j1 * @jsize4] - @lhs[1 + j1 * @jsize4] * @lhs[3 + j * @jsize4]
              @lhs[3 + j1 * @jsize4] = @lhs[3 + j1 * @jsize4] - @lhs[1 + j1 * @jsize4] * @lhs[4 + j * @jsize4]
              for m in 0..2 do
                  @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhs[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              end
              fac2 = 1.0 / @lhs[2 + j1 * @jsize4]
              for m in 0..2 do
                  @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = fac2 * @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1]
              end
              for j in 0..@grid_points[1] - 3 do
                  j1 = j + 1
                  j2 = j + 2
                  m = 3
                  fac1 = 1.0 / @lhsp[2 + j * @jsize4]
                  @lhsp[3 + j * @jsize4] = fac1 * @lhsp[3 + j * @jsize4]
                  @lhsp[4 + j * @jsize4] = fac1 * @lhsp[4 + j * @jsize4]
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  @lhsp[2 + j1 * @jsize4] = @lhsp[2 + j1 * @jsize4] - @lhsp[1 + j1 * @jsize4] * @lhsp[3 + j * @jsize4]
                  @lhsp[3 + j1 * @jsize4] = @lhsp[3 + j1 * @jsize4] - @lhsp[1 + j1 * @jsize4] * @lhsp[4 + j * @jsize4]
                  @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsp[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  @lhsp[1 + j2 * @jsize4] = @lhsp[1 + j2 * @jsize4] - @lhsp[0 + j2 * @jsize4] * @lhsp[3 + j * @jsize4]
                  @lhsp[2 + j2 * @jsize4] = @lhsp[2 + j2 * @jsize4] - @lhsp[0 + j2 * @jsize4] * @lhsp[4 + j * @jsize4]
                  @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] - @lhsp[0 + j2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  m = 4
                  fac1 = 1.0 / @lhsm[2 + j * @jsize4]
                  @lhsm[3 + j * @jsize4] = fac1 * @lhsm[3 + j * @jsize4]
                  @lhsm[4 + j * @jsize4] = fac1 * @lhsm[4 + j * @jsize4]
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  @lhsm[2 + j1 * @jsize4] = @lhsm[2 + j1 * @jsize4] - @lhsm[1 + j1 * @jsize4] * @lhsm[3 + j * @jsize4]
                  @lhsm[3 + j1 * @jsize4] = @lhsm[3 + j1 * @jsize4] - @lhsm[1 + j1 * @jsize4] * @lhsm[4 + j * @jsize4]
                  @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsm[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                  @lhsm[1 + j2 * @jsize4] = @lhsm[1 + j2 * @jsize4] - @lhsm[0 + j2 * @jsize4] * @lhsm[3 + j * @jsize4]
                  @lhsm[2 + j2 * @jsize4] = @lhsm[2 + j2 * @jsize4] - @lhsm[0 + j2 * @jsize4] * @lhsm[4 + j * @jsize4]
                  @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] - @lhsm[0 + j2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              end
              j = @grid_points[1] - 2
              j1 = @grid_points[1] - 1
              m = 3
              fac1 = 1.0 / @lhsp[2 + j * @jsize4]
              @lhsp[3 + j * @jsize4] = fac1 * @lhsp[3 + j * @jsize4]
              @lhsp[4 + j * @jsize4] = fac1 * @lhsp[4 + j * @jsize4]
              @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              @lhsp[2 + j1 * @jsize4] = @lhsp[2 + j1 * @jsize4] - @lhsp[1 + j1 * @jsize4] * @lhsp[3 + j * @jsize4]
              @lhsp[3 + j1 * @jsize4] = @lhsp[3 + j1 * @jsize4] - @lhsp[1 + j1 * @jsize4] * @lhsp[4 + j * @jsize4]
              @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsp[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              m = 4
              fac1 = 1.0 / @lhsm[2 + j * @jsize4]
              @lhsm[3 + j * @jsize4] = fac1 * @lhsm[3 + j * @jsize4]
              @lhsm[4 + j * @jsize4] = fac1 * @lhsm[4 + j * @jsize4]
              @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              @lhsm[2 + j1 * @jsize4] = @lhsm[2 + j1 * @jsize4] - @lhsm[1 + j1 * @jsize4] * @lhsm[3 + j * @jsize4]
              @lhsm[3 + j1 * @jsize4] = @lhsm[3 + j1 * @jsize4] - @lhsm[1 + j1 * @jsize4] * @lhsm[4 + j * @jsize4]
              @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsm[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
              @rhs[3 + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j1 * @jsize1 + k * @ksize1] / @lhsp[2 + j1 * @jsize4]
              @rhs[4 + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j1 * @jsize1 + k * @ksize1] / @lhsm[2 + j1 * @jsize4]
              j = @grid_points[1] - 2
              j1 = @grid_points[1] - 1
              for m in 0..2 do
                  @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[3 + j * @jsize4] * @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1]
              end
              @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[3 + j * @jsize4] * @rhs[3 + i * @isize1 + j1 * @jsize1 + k * @ksize1]
              @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[3 + j * @jsize4] * @rhs[4 + i * @isize1 + j1 * @jsize1 + k * @ksize1]
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
              j = @grid_points[1] - 3
              while j >= 0 do
                  j1 = j + 1
                  j2 = j + 2
                  for m in 0..2 do
                      @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[3 + j * @jsize4] * @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhs[4 + j * @jsize4] * @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1]
                  end
                  @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[3 + j * @jsize4] * @rhs[3 + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsp[4 + j * @jsize4] * @rhs[3 + i * @isize1 + j2 * @jsize1 + k * @ksize1]
                  @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[3 + j * @jsize4] * @rhs[4 + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsm[4 + j * @jsize4] * @rhs[4 + i * @isize1 + j2 * @jsize1 + k * @ksize1]
                j -= 1
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
                  t1 = @@bt * r1
                  t2 = 0.5 * (r4 + r5)
                  @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @@bt * (r4 - r5)
                  @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = -r3
                  @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = r2
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
