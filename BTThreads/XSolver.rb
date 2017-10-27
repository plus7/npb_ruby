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
class XSolver < BTBase
  def initialize(bt, low, high)
    super(bt.clss, bt.num_threads)
    @master =  nil 
    @done = true
    @fjac =  nil 
    @njac =  nil 
    @lhs =  nil 
    Init(bt)
    @lower_bound = low
    @upper_bound = high
    @fjac = Array.new(5 * 5 * (@problem_size + 1), 0.0)
    @njac = Array.new(5 * 5 * (@problem_size + 1), 0.0)
    @lhs = Array.new(5 * 5 * 3 * (@problem_size + 1), 0.0)
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
    
    isize = @grid_points[0] - 1
    for k in @lower_bound..@upper_bound do
        for j in 1..@grid_points[1] - 2 do
            for i in 0..isize do
                @tmp1 = @rho_i[i + j * @jsize1 + k * @ksize1]
                @tmp2 = @tmp1 * @tmp1
                @tmp3 = @tmp1 * @tmp2
                @fjac[0 + 0 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[0 + 1 * @@isize4 + i * @@jsize4] = 1.0
                @fjac[0 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[0 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[0 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[1 + 0 * @@isize4 + i * @@jsize4] = -(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@c2 * @qs[i + j * @jsize1 + k * @ksize1]
                @fjac[1 + 1 * @@isize4 + i * @@jsize4] = (2.0 - @@c2) * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] / @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2])
                @fjac[1 + 2 * @@isize4 + i * @@jsize4] = -@@c2 * (@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1)
                @fjac[1 + 3 * @@isize4 + i * @@jsize4] = -@@c2 * (@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1)
                @fjac[1 + 4 * @@isize4 + i * @@jsize4] = @@c2
                @fjac[2 + 0 * @@isize4 + i * @@jsize4] = -(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[2 + 1 * @@isize4 + i * @@jsize4] = @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 2 * @@isize4 + i * @@jsize4] = @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[2 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[3 + 0 * @@isize4 + i * @@jsize4] = -(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[3 + 1 * @@isize4 + i * @@jsize4] = @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[3 + 3 * @@isize4 + i * @@jsize4] = @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[4 + 0 * @@isize4 + i * @@jsize4] = (@@c2 * 2.0 * @square[i + j * @jsize1 + k * @ksize1] - @@c1 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2]) * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2)
                @fjac[4 + 1 * @@isize4 + i * @@jsize4] = @@c1 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1 - @@c2 * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2 + @qs[i + j * @jsize1 + k * @ksize1])
                @fjac[4 + 2 * @@isize4 + i * @@jsize4] = -@@c2 * (@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[4 + 3 * @@isize4 + i * @@jsize4] = -@@c2 * (@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[4 + 4 * @@isize4 + i * @@jsize4] = @@c1 * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1)
                @njac[0 + 0 * @@isize4 + i * @@jsize4] = 0.0
                @njac[0 + 1 * @@isize4 + i * @@jsize4] = 0.0
                @njac[0 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @njac[0 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @njac[0 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @njac[1 + 0 * @@isize4 + i * @@jsize4] = -@@con43 * @@c3c4 * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[1 + 1 * @@isize4 + i * @@jsize4] = @@con43 * @@c3c4 * @tmp1
                @njac[1 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @njac[1 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @njac[1 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @njac[2 + 0 * @@isize4 + i * @@jsize4] = -@@c3c4 * @tmp2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[2 + 1 * @@isize4 + i * @@jsize4] = 0.0
                @njac[2 + 2 * @@isize4 + i * @@jsize4] = @@c3c4 * @tmp1
                @njac[2 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @njac[2 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @njac[3 + 0 * @@isize4 + i * @@jsize4] = -@@c3c4 * @tmp2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[3 + 1 * @@isize4 + i * @@jsize4] = 0.0
                @njac[3 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @njac[3 + 3 * @@isize4 + i * @@jsize4] = @@c3c4 * @tmp1
                @njac[3 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @njac[4 + 0 * @@isize4 + i * @@jsize4] = -(@@con43 * @@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - (@@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - (@@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - @@c1345 * @tmp2 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 1 * @@isize4 + i * @@jsize4] = (@@con43 * @@c3c4 - @@c1345) * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 2 * @@isize4 + i * @@jsize4] = (@@c3c4 - @@c1345) * @tmp2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 3 * @@isize4 + i * @@jsize4] = (@@c3c4 - @@c1345) * @tmp2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 4 * @@isize4 + i * @@jsize4] = (@@c1345) * @tmp1
            end
            lhsinit(@lhs, isize)
            for i in 1..isize - 1 do
                @tmp1 = @@dt * @@tx1
                @tmp2 = @@dt * @@tx2
                @lhs[0 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx1
                @lhs[0 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 1 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 2 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 3 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 4 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 0 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx2
                @lhs[1 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 2 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 3 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 4 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 0 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 1 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx3
                @lhs[2 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 3 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 4 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 0 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 1 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 2 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx4
                @lhs[3 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 4 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 0 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 1 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 2 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 3 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx5
                @lhs[0 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[0 + 0 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx1
                @lhs[0 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 1 * @@isize4 + i * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 2 * @@isize4 + i * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 3 * @@isize4 + i * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 4 * @@isize4 + i * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 0 * @@isize4 + i * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[1 + 1 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx2
                @lhs[1 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 2 * @@isize4 + i * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 3 * @@isize4 + i * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 4 * @@isize4 + i * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 0 * @@isize4 + i * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 1 * @@isize4 + i * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[2 + 2 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx3
                @lhs[2 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 3 * @@isize4 + i * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 4 * @@isize4 + i * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 0 * @@isize4 + i * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 1 * @@isize4 + i * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 2 * @@isize4 + i * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[3 + 3 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx4
                @lhs[3 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 4 * @@isize4 + i * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 0 * @@isize4 + i * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 1 * @@isize4 + i * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 2 * @@isize4 + i * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 3 * @@isize4 + i * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[4 + 4 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx5
                @lhs[0 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx1
                @lhs[0 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 1 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 2 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 3 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 4 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 0 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx2
                @lhs[1 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 2 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 3 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 4 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 0 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 1 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx3
                @lhs[2 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 3 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 4 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 0 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 1 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 2 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx4
                @lhs[3 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 4 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 0 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 1 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 2 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 3 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx5
            end
            binvcrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + 0 * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + 0 * @@ksize4, @rhs, 0 + 0 * @isize2 + j * @jsize2 + k * @ksize2)
            for i in 1..isize - 1 do
                matvec_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4, @rhs, 0 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2, @rhs, 0 + i * @isize2 + j * @jsize2 + k * @ksize2)
                matmul_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + (i - 1) * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4)
                binvcrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4, @rhs, 0 + i * @isize2 + j * @jsize2 + k * @ksize2)
            end
            matvec_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + isize * @@ksize4, @rhs, 0 + (isize - 1) * @isize2 + j * @jsize2 + k * @ksize2, @rhs, 0 + isize * @isize2 + j * @jsize2 + k * @ksize2)
            matmul_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + isize * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + (isize - 1) * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + isize * @@ksize4)
            binvrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + isize * @@ksize4, @rhs, 0 + isize * @isize2 + j * @jsize2 + k * @ksize2)
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
            i = isize - 1
            while i >= 0 do
                for m in 0..@@BLOCK_SIZE - 1 do
                    for n in 0..@@BLOCK_SIZE - 1 do
                        @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] -= @lhs[m + n * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] * @rhs[n + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2]
                    end
                end
              i -= 1
            end
        end
    end
  end

  attr_accessor :id

  attr_accessor :done

end
