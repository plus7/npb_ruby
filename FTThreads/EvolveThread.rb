# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
#  Translation to Java and to MultiThreaded Code
#           M. Frumkin
#           M. Schultz
#  Translation to Ruby
#           T. Nose
#           H. Tomari
# -------------------------------------------------------------------------

require "mda"
class EvolveThread < FTBase
  def initialize(ft, low1, high1)
    super(ft.clss, ft.num_threads, true)
    @kt = 0
    @done = true
    @@ap = (-4.0 * @@alpha * @@pi * @@pi)
    Init(ft)
    @lower_bound1 = low1
    @upper_bound1 = high1
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = ft
  end

  def Init(ft)
    @xtr = ft.xtr
    @xnt = ft.xnt
    @nx = ft.nx
    @ny = ft.ny
    @nz = ft.nz
    @ixtr = ft.isize3
    @jxtr = ft.jsize3
    @kxtr = ft.ksize3
    @ixnt = ft.isize4
    @jxnt = ft.jsize4
    @kxnt = ft.ksize4
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
    for i in @lower_bound1..@upper_bound1 do
        ii = i - (i / (@nx / 2)) * @nx; 
        ii2 = ii * ii; 
        for k in 0..@nz-1 do
            kk = k - (k / (@nz / 2)) * @nz; 
            ik2 = ii2 + kk * kk; 
            for j in 0..@ny-1 do
                jj = j - (j / (@ny / 2)) * @ny; 
                lexp = Math.exp((@@ap * (jj * jj + ik2)) * (@kt + 1)); 
                xntidx = j * @ixnt + k * @jxnt + i * @kxnt; 
                xtridx = j * @ixtr + i * @jxtr + k * @kxtr; 
                @xnt[@@REAL + xntidx] = lexp * @xtr[@@REAL + xtridx]
                @xnt[@@IMAG + xntidx] = lexp * @xtr[@@IMAG + xtridx]
            end
        end
    end
  end

# *** public ***

  attr_accessor :kt

  attr_accessor :id

  attr_accessor :done

end
