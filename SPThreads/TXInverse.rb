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
class TXInverse < SPBase
  def initialize(sp, low, high)
    super(sp.clss, sp.num_threads)
    @master =  nil 
    @done = true
    Init(sp)
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = sp
    @lower_bound = low
    @upper_bound = high
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
    
    
    for k in @lower_bound..@upper_bound do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                ru1 = @rho_i[i + j * @jsize2 + k * @ksize2]
                uu = @us[i + j * @jsize2 + k * @ksize2]
                vv = @vs[i + j * @jsize2 + k * @ksize2]
                ww = @ws[i + j * @jsize2 + k * @ksize2]
                ac = @speed[i + j * @jsize2 + k * @ksize2]
                ac2inv = 1.0 / (ac * ac)
                r1 = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r2 = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r3 = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r4 = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r5 = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                t1 = @@c2 * ac2inv * (@qs[i + j * @jsize2 + k * @ksize2] * r1 - uu * r2 - vv * r3 - ww * r4 + r5)
                t2 = @@bt * ru1 * (uu * r1 - r2)
                t3 = (@@bt * ru1 * ac) * t1
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = r1 - t1
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = -ru1 * (ww * r1 - r4)
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = ru1 * (vv * r1 - r3)
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = -t2 + t3
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = t2 + t3
            end
        end
    end
  end

# *** public ***
  attr_accessor :id

  attr_accessor :done

end
