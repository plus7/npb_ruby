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
class RHSAdder < SPBase
  def initialize(sp, low, high)
    super(sp.clss, sp.num_threads)
    @master =  nil 
    @done = true
    Init(sp)
    @lower_bound = low
    @upper_bound = high
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
    
    for k in @lower_bound..@upper_bound do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                for m in 0..4 do
                # puts sprintf("@u[%d]=%f",
                #              m + i * @isize1 + j * @jsize1 + k * @ksize1,
                #              @u[m + i * @isize1 + j * @jsize1 + k * @ksize1])

                    @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] + @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]

                # puts sprintf("@rhs[%d]=%f",
                #              m + i * @isize1 + j * @jsize1 + k * @ksize1,
                #              @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1])
                # puts sprintf("%d, %d, %d", k, j, i)
                end
            end
        end
    end
  end

# *** public ***
  attr_accessor :id

  attr_accessor :done

end
