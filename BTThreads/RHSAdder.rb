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
class RHSAdder < BTBase
  def initialize(bt, low, high)
    super(bt.clss, bt.num_threads)
    @master =  nil 
    @done = true
    Init(bt)
    @lower_bound = low
    @upper_bound = high
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
            for k in @lower_bound..@upper_bound do
                for j in 1..@grid_points[1] - 2 do
                    for i in 1..@grid_points[0] - 2 do
                        for m in 0..4 do
                            @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] += @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2]
                        end
                    end
                end
            end
            @master.synchronize do
                @done = true
                @master.cond.signal
            end
        end
    end
  end

# *** public ***
  attr_accessor :id

  attr_accessor :done
end
