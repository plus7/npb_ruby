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
class Scale < LUBase
  def initialize(lu, low1, high1)
    super(lu.clss, lu.num_threads)
    @master =  nil 
    @rhscomputer =  nil 
    @scaler =  nil 
    @adder =  nil 
    @lowerjac =  nil 
    @upperjac =  nil 
    @done = true
    Init(lu)
    @lower_bound1 = low1
    @upper_bound1 = high1
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = lu
  end

  def Init(lu)
    @isiz1 = lu.isiz1
    @isiz2 = lu.isiz2
    @isiz3 = lu.isiz3
    @itmax_default = lu.itmax_default
    @dt_default = lu.dt_default
    @inorm_default = lu.inorm_default
    @u = lu.u
    @rsd = lu.rsd
    @frct = lu.frct
    @isize1 = lu.isize1
    @jsize1 = lu.jsize1
    @ksize1 = lu.ksize1
    @flux = lu.flux
    @isize2 = lu.isize2
    @qs = lu.qs
    @rho_i = lu.rho_i
    @jsize3 = lu.jsize3
    @ksize3 = lu.ksize3
    @a = lu.a
    @b = lu.b
    @c = lu.c
    @d = lu.d
    @isize4 = lu.isize4
    @jsize4 = lu.jsize4
    @ksize4 = lu.ksize4
    @nx = lu.nx
    @ny = lu.ny
    @nz = lu.nz
    @nx0 = lu.nx0
    @ny0 = lu.ny0
    @nz0 = lu.nz0
    @ist = lu.ist
    @iend = lu.iend
    @jst = lu.jst
    @jend = lu.jend
    @ii1 = lu.ii1
    @ii2 = lu.ii2
    @ji1 = lu.ji1
    @ji2 = lu.ji2
    @ki1 = lu.ki1
    @ki2 = lu.ki2
    @dxi = lu.dxi
    @deta = lu.deta
    @dzeta = lu.dzeta
    @tx1 = lu.tx1
    @tx2 = lu.tx2
    @tx3 = lu.tx3
    @ty1 = lu.ty1
    @ty2 = lu.ty2
    @ty3 = lu.ty3
    @tz1 = lu.tz1
    @tz2 = lu.tz1
    @tz3 = lu.tz3
    @dx1 = lu.dx1
    @dx2 = lu.dx2
    @dx3 = lu.dx3
    @dx4 = lu.dx4
    @dx5 = lu.dx5
    @dy1 = lu.dy1
    @dy2 = lu.dy2
    @dy3 = lu.dy3
    @dy4 = lu.dy4
    @dy5 = lu.dy5
    @dz1 = lu.dz1
    @dz2 = lu.dz2
    @dz3 = lu.dz3
    @dz4 = lu.dz4
    @dz5 = lu.dz5
    @dssp = lu.dssp
    @dt = lu.dt
    @omega = lu.omega
    @frc = lu.frc
    @ttotal = lu.ttotal
  end

  def run()
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    while true do
         self.synchronize do
            while @done == true do
                begin
                    cond.wait()
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
    
    for k in @lower_bound1..@upper_bound1 do
        for j in @jst - 1..@jend - 1 do
            for i in @ist - 1..@iend - 1 do
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @dt * @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
            end
        end
    end
  end

end
