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
class ZSolver < SPBase
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
    @rhos = Array.new(@problem_size, 0.0)
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
    
    
    rtmp = Array.new(5 * (@KMAX + 1), 0.0); 
    case @state
    when 1
      for j in @lower_bound..@upper_bound do
          for i in 1..@nx2 do
              for k in 0..@nz2 + 1 do
                  ru1 = @@c3c4 * @rho_i[i + j * @jsize2 + k * @ksize2]
                  @cv[k] = @ws[i + j * @jsize2 + k * @ksize2]
                  @rhos[k] = dmax1ex(@@dz4 + @@con43 * ru1, @@dz5 + @@c1c5 * ru1, @@dzmax + ru1, @@dz1)
              end
              lhsinit(@grid_points[2] - 1)
              for k in 1..@nz2 do
                  @lhs[0 + k * @jsize4] = 0.0
                  @lhs[1 + k * @jsize4] = -@@dttz2 * @cv[k - 1] - @@dttz1 * @rhos[k - 1]
                  @lhs[2 + k * @jsize4] = 1.0 + @@c2dttz1 * @rhos[k]
                  @lhs[3 + k * @jsize4] = @@dttz2 * @cv[k + 1] - @@dttz1 * @rhos[k + 1]
                  @lhs[4 + k * @jsize4] = 0.0
              end
              k = 1
              @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz5
              @lhs[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@comz4
              @lhs[4 + k * @jsize4] = @lhs[4 + k * @jsize4] + @@comz1
              k = 2
              @lhs[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@comz4
              @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz6
              @lhs[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@comz4
              @lhs[4 + k * @jsize4] = @lhs[4 + k * @jsize4] + @@comz1
              for k in 3..@nz2 - 2 do
                  @lhs[0 + k * @jsize4] = @lhs[0 + k * @jsize4] + @@comz1
                  @lhs[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@comz4
                  @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz6
                  @lhs[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@comz4
                  @lhs[4 + k * @jsize4] = @lhs[4 + k * @jsize4] + @@comz1
              end
              k = @nz2 - 1
              @lhs[0 + k * @jsize4] = @lhs[0 + k * @jsize4] + @@comz1
              @lhs[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@comz4
              @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz6
              @lhs[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@comz4
              k = @nz2
              @lhs[0 + k * @jsize4] = @lhs[0 + k * @jsize4] + @@comz1
              @lhs[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@comz4
              @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz5
              for k in 1..@nz2 do
                  @lhsp[0 + k * @jsize4] = @lhs[0 + k * @jsize4]
                  @lhsp[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@dttz2 * @speed[i + j * @jsize2 + (k - 1) * @ksize2]
                  @lhsp[2 + k * @jsize4] = @lhs[2 + k * @jsize4]
                  @lhsp[3 + k * @jsize4] = @lhs[3 + k * @jsize4] + @@dttz2 * @speed[i + j * @jsize2 + (k + 1) * @ksize2]
                  @lhsp[4 + k * @jsize4] = @lhs[4 + k * @jsize4]
                  @lhsm[0 + k * @jsize4] = @lhs[0 + k * @jsize4]
                  @lhsm[1 + k * @jsize4] = @lhs[1 + k * @jsize4] + @@dttz2 * @speed[i + j * @jsize2 + (k - 1) * @ksize2]
                  @lhsm[2 + k * @jsize4] = @lhs[2 + k * @jsize4]
                  @lhsm[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@dttz2 * @speed[i + j * @jsize2 + (k + 1) * @ksize2]
                  @lhsm[4 + k * @jsize4] = @lhs[4 + k * @jsize4]
              end
              for k in 0..@nz2 + 1 do
                  rtmp[0 + k * 5] = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  rtmp[1 + k * 5] = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  rtmp[2 + k * 5] = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  rtmp[3 + k * 5] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  rtmp[4 + k * 5] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
              end
              for k in 0..@grid_points[2] - 3 do
                  k1 = k + 1
                  k2 = k + 2
                  fac1 = 1.0 / @lhs[2 + k * @jsize4]
                  @lhs[3 + k * @jsize4] = fac1 * @lhs[3 + k * @jsize4]
                  @lhs[4 + k * @jsize4] = fac1 * @lhs[4 + k * @jsize4]
                  for m in 0..2 do
                      rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
                  end
                  @lhs[2 + k1 * @jsize4] = @lhs[2 + k1 * @jsize4] - @lhs[1 + k1 * @jsize4] * @lhs[3 + k * @jsize4]
                  @lhs[3 + k1 * @jsize4] = @lhs[3 + k1 * @jsize4] - @lhs[1 + k1 * @jsize4] * @lhs[4 + k * @jsize4]
                  for m in 0..2 do
                      rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhs[1 + k1 * @jsize4] * rtmp[m + k * 5]
                  end
                  @lhs[1 + k2 * @jsize4] = @lhs[1 + k2 * @jsize4] - @lhs[0 + k2 * @jsize4] * @lhs[3 + k * @jsize4]
                  @lhs[2 + k2 * @jsize4] = @lhs[2 + k2 * @jsize4] - @lhs[0 + k2 * @jsize4] * @lhs[4 + k * @jsize4]
                  for m in 0..2 do
                      rtmp[m + k2 * 5] = rtmp[m + k2 * 5] - @lhs[0 + k2 * @jsize4] * rtmp[m + k * 5]
                  end
              end
              k = @grid_points[2] - 2
              k1 = @grid_points[2] - 1
              fac1 = 1.0 / @lhs[2 + k * @jsize4]
              @lhs[3 + k * @jsize4] = fac1 * @lhs[3 + k * @jsize4]
              @lhs[4 + k * @jsize4] = fac1 * @lhs[4 + k * @jsize4]
              for m in 0..2 do
                  rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
              end
              @lhs[2 + k1 * @jsize4] = @lhs[2 + k1 * @jsize4] - @lhs[1 + k1 * @jsize4] * @lhs[3 + k * @jsize4]
              @lhs[3 + k1 * @jsize4] = @lhs[3 + k1 * @jsize4] - @lhs[1 + k1 * @jsize4] * @lhs[4 + k * @jsize4]
              for m in 0..2 do
                  rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhs[1 + k1 * @jsize4] * rtmp[m + k * 5]
              end
              fac2 = 1.0 / @lhs[2 + k1 * @jsize4]
              for m in 0..2 do
                  rtmp[m + k1 * 5] = fac2 * rtmp[m + k1 * 5]
              end
              for k in 0..@grid_points[2] - 3 do
                  k1 = k + 1
                  k2 = k + 2
                  m = 3
                  fac1 = 1.0 / @lhsp[2 + k * @jsize4]
                  @lhsp[3 + k * @jsize4] = fac1 * @lhsp[3 + k * @jsize4]
                  @lhsp[4 + k * @jsize4] = fac1 * @lhsp[4 + k * @jsize4]
                  rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
                  @lhsp[2 + k1 * @jsize4] = @lhsp[2 + k1 * @jsize4] - @lhsp[1 + k1 * @jsize4] * @lhsp[3 + k * @jsize4]
                  @lhsp[3 + k1 * @jsize4] = @lhsp[3 + k1 * @jsize4] - @lhsp[1 + k1 * @jsize4] * @lhsp[4 + k * @jsize4]
                  rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhsp[1 + k1 * @jsize4] * rtmp[m + k * 5]
                  @lhsp[1 + k2 * @jsize4] = @lhsp[1 + k2 * @jsize4] - @lhsp[0 + k2 * @jsize4] * @lhsp[3 + k * @jsize4]
                  @lhsp[2 + k2 * @jsize4] = @lhsp[2 + k2 * @jsize4] - @lhsp[0 + k2 * @jsize4] * @lhsp[4 + k * @jsize4]
                  rtmp[m + k2 * 5] = rtmp[m + k2 * 5] - @lhsp[0 + k2 * @jsize4] * rtmp[m + k * 5]
                  m = 4
                  fac1 = 1.0 / @lhsm[2 + k * @jsize4]
                  @lhsm[3 + k * @jsize4] = fac1 * @lhsm[3 + k * @jsize4]
                  @lhsm[4 + k * @jsize4] = fac1 * @lhsm[4 + k * @jsize4]
                  rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
                  @lhsm[2 + k1 * @jsize4] = @lhsm[2 + k1 * @jsize4] - @lhsm[1 + k1 * @jsize4] * @lhsm[3 + k * @jsize4]
                  @lhsm[3 + k1 * @jsize4] = @lhsm[3 + k1 * @jsize4] - @lhsm[1 + k1 * @jsize4] * @lhsm[4 + k * @jsize4]
                  rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhsm[1 + k1 * @jsize4] * rtmp[m + k * 5]
                  @lhsm[1 + k2 * @jsize4] = @lhsm[1 + k2 * @jsize4] - @lhsm[0 + k2 * @jsize4] * @lhsm[3 + k * @jsize4]
                  @lhsm[2 + k2 * @jsize4] = @lhsm[2 + k2 * @jsize4] - @lhsm[0 + k2 * @jsize4] * @lhsm[4 + k * @jsize4]
                  rtmp[m + k2 * 5] = rtmp[m + k2 * 5] - @lhsm[0 + k2 * @jsize4] * rtmp[m + k * 5]
              end
              k = @grid_points[2] - 2
              k1 = @grid_points[2] - 1
              m = 3
              fac1 = 1.0 / @lhsp[2 + k * @jsize4]
              @lhsp[3 + k * @jsize4] = fac1 * @lhsp[3 + k * @jsize4]
              @lhsp[4 + k * @jsize4] = fac1 * @lhsp[4 + k * @jsize4]
              rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
              @lhsp[2 + k1 * @jsize4] = @lhsp[2 + k1 * @jsize4] - @lhsp[1 + k1 * @jsize4] * @lhsp[3 + k * @jsize4]
              @lhsp[3 + k1 * @jsize4] = @lhsp[3 + k1 * @jsize4] - @lhsp[1 + k1 * @jsize4] * @lhsp[4 + k * @jsize4]
              rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhsp[1 + k1 * @jsize4] * rtmp[m + k * 5]
              m = 4
              fac1 = 1.0 / @lhsm[2 + k * @jsize4]
              @lhsm[3 + k * @jsize4] = fac1 * @lhsm[3 + k * @jsize4]
              @lhsm[4 + k * @jsize4] = fac1 * @lhsm[4 + k * @jsize4]
              rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
              @lhsm[2 + k1 * @jsize4] = @lhsm[2 + k1 * @jsize4] - @lhsm[1 + k1 * @jsize4] * @lhsm[3 + k * @jsize4]
              @lhsm[3 + k1 * @jsize4] = @lhsm[3 + k1 * @jsize4] - @lhsm[1 + k1 * @jsize4] * @lhsm[4 + k * @jsize4]
              rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhsm[1 + k1 * @jsize4] * rtmp[m + k * 5]
              rtmp[3 + k1 * 5] = rtmp[3 + k1 * 5] / @lhsp[2 + k1 * @jsize4]
              rtmp[4 + k1 * 5] = rtmp[4 + k1 * 5] / @lhsm[2 + k1 * @jsize4]
              k = @grid_points[2] - 2
              k1 = @grid_points[2] - 1
              for m in 0..2 do
                  rtmp[m + k * 5] = rtmp[m + k * 5] - @lhs[3 + k * @jsize4] * rtmp[m + k1 * 5]
              end
              rtmp[3 + k * 5] = rtmp[3 + k * 5] - @lhsp[3 + k * @jsize4] * rtmp[3 + k1 * 5]
              rtmp[4 + k * 5] = rtmp[4 + k * 5] - @lhsm[3 + k * @jsize4] * rtmp[4 + k1 * 5]
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
              k = @grid_points[2] - 3
              while k >= 0 do
                  k1 = k + 1
                  k2 = k + 2
                  for m in 0..2 do
                      rtmp[m + k * 5] = rtmp[m + k * 5] - @lhs[3 + k * @jsize4] * rtmp[m + k1 * 5] - @lhs[4 + k * @jsize4] * rtmp[m + k2 * 5]
                  end
                  rtmp[3 + k * 5] = rtmp[3 + k * 5] - @lhsp[3 + k * @jsize4] * rtmp[3 + k1 * 5] - @lhsp[4 + k * @jsize4] * rtmp[3 + k2 * 5]
                  rtmp[4 + k * 5] = rtmp[4 + k * 5] - @lhsm[3 + k * @jsize4] * rtmp[4 + k1 * 5] - @lhsm[4 + k * @jsize4] * rtmp[4 + k2 * 5]
                k -= 1
              end
              for k in 0..@nz2 + 1 do
                  @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[0 + k * 5]
                  @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[1 + k * 5]
                  @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[2 + k * 5]
                  @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[3 + k * 5]
                  @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[4 + k * 5]
              end
          end
      end
    when 2
      for k in @lower_bound..@upper_bound do
          for j in 1..@ny2 do
              for i in 1..@nx2 do
                  xvel = @us[i + j * @jsize2 + k * @ksize2]
                  yvel = @vs[i + j * @jsize2 + k * @ksize2]
                  zvel = @ws[i + j * @jsize2 + k * @ksize2]
                  ac = @speed[i + j * @jsize2 + k * @ksize2]
                  ac2u = ac * ac
                  r1 = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  r2 = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  r3 = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  r4 = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  r5 = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  uzik1 = @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                  btuz = @@bt * uzik1
                  t1 = btuz / ac * (r4 + r5)
                  t2 = r3 + t1
                  t3 = btuz * (r4 - r5)
                  @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = t2
                  @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = -uzik1 * r2 + xvel * t2
                  @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = uzik1 * r1 + yvel * t2
                  @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = zvel * t2 + t3
                  @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = uzik1 * (-xvel * r2 + yvel * r1) + @qs[i + j * @jsize2 + k * @ksize2] * t2 + @@c2iv * ac2u * t1 + zvel * t3
              end
          end
      end
    end
;
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
