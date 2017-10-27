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
class UpperJac < LUBase
  attr_accessor :neighbor
  attr_accessor :todo
  def initialize(lu, low1, high1)
    super(lu.clss, lu.num_threads)
    @master =  nil 
    @rhscomputer =  nil 
    @scaler =  nil 
    @adder =  nil 
    @lowerjac =  nil 
    @upperjac =  nil 
    @done = true
    @synch = Array.new(@isiz1 + 1, 0.0)
    @todo = 0
    Init(lu)
    @lower_bound1 = low1
    @upper_bound1 = high1
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @tv = Array.new(5 * (@isiz1 / 2 * 2 + 1) * @isiz2, 0.0)
    @master = lu
  end

  def Init(lu)
    @num_threads = lu.num_threads
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
    @tz2 = lu.tz2
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
    @udz = lu.c
    @udy = lu.b
    @udx = lu.a
    @d = lu.d
    @v = lu.rsd
    @ldmx = lu.isiz1
    @ldmy = lu.isiz2
    @ldmz = lu.isiz3
  end

  def run()
    
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    while true do
         self.synchronize do
            while @done == true do
                #begin
                    cond.wait
                #rescue InterruptedException => ie
                #ensure
                #end
            end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
            k = @nz - 2
            while k >= 1 do
                step(k)
                if @id > 0 then
                    while @todo <= 0 do
                        begin
                            cond.wait
                        rescue Exception => e
                        ensure
                        end
                        @master.synchronize do
                            @master.cond.signal
                        end
                    end
                else
                end
                step2(k)
                @todo -= 1
                if @id < @num_threads - 1 then
                  @neighbor[@id + 1].synchronize do
                      @neighbor[id + 1].todo += 1
                      @neighbor[@id + 1].cond.signal
                  end
                else
                end
              k -= 1
            end
            @done = true
            if @id == @num_threads - 1 then
              @master.synchronize do
                  @master.cond.signal
              end
            else
              @neighbor[@id + 1].synchronize do
                  @neighbor[@id + 1].cond.signal
              end
            end
        end
    end
  end

  def step(k)
    
    
    
    
    
    tmat = Array.new(5 * 5, 0.0); 
    tmp = 1.0 / (@omega * (2.0 - @omega))
    r43 = (4.0 / 3.0)
    c1345 = @c1 * @c3 * @c4 * @c5
    c34 = @c3 * @c4
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    j = @upper_bound1
    while j >= @lower_bound1 do
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
        i = @iend - 1
        while i >= @ist - 1 do
            tmp1 = @rho_i[i + j * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @d[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * (@tx1 * @dx1 + @ty1 * @dy1 + @tz1 * @dz1)
            @d[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (-@tx1 * r43 - @ty1 - @tz1) * (c34 * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 * r43 + @ty1 + @tz1) + @dt * 2.0 * (@tx1 * @dx2 + @ty1 * @dy2 + @tz1 * @dz2)
            @d[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (-@tx1 - @ty1 * r43 - @tz1) * (c34 * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 + @ty1 * r43 + @tz1) + @dt * 2.0 * (@tx1 * @dx3 + @ty1 * @dy3 + @tz1 * @dz3)
            @d[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (-@tx1 - @ty1 - @tz1 * r43) * (c34 * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 + @ty1 + @tz1 * r43) + @dt * 2.0 * (@tx1 * @dx4 + @ty1 * @dy4 + @tz1 * @dz4)
            @d[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * 2.0 * (((@tx1 * (r43 * c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (c34 - c1345)) * Math_dot_pow(@u[1 + i * @isize1 + j * @jsize1 + k * @ksize1], 2) + (@tx1 * (c34 - c1345) + @ty1 * (r43 * c34 - c1345) + @tz1 * (c34 - c1345)) * Math_dot_pow(@u[2 + i * @isize1 + j * @jsize1 + k * @ksize1], 2) + (@tx1 * (c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (r43 * c34 - c1345)) * Math_dot_pow(@u[3 + i * @isize1 + j * @jsize1 + k * @ksize1], 2)) * tmp3 + (@tx1 + @ty1 + @tz1) * c1345 * tmp2 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (@tx1 * (r43 * c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (c34 - c1345)) * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (@tx1 * (c34 - c1345) + @ty1 * (r43 * c34 - c1345) + @tz1 * (c34 - c1345)) * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (@tx1 * (c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (r43 * c34 - c1345)) * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * (@tx1 + @ty1 + @tz1) * c1345 * tmp1 + @dt * 2.0 * (@tx1 * @dx5 + @ty1 * @dy5 + @tz1 * @dz5)
            tmp1 = @rho_i[(i + 1) + j * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @a[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx1 * @dx1
            @a[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2
            @a[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-Math_dot_pow((@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1), 2) + @c2 * @qs[(i + 1) + j * @jsize3 + k * @ksize3] * tmp1) - @dt * @tx1 * (-r43 * c34 * tmp2 * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @a[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * ((2.0 - @c2) * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)) - @dt * @tx1 * (r43 * c34 * tmp1) - @dt * @tx1 * @dx2
            @a[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-@c2 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1))
            @a[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-@c2 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1))
            @a[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * @c2
            @a[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-(@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (-c34 * tmp2 * @u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @a[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)
            @a[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @dt * @tx1 * (c34 * tmp1) - @dt * @tx1 * @dx3
            @a[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-(@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (-c34 * tmp2 * @u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @a[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)
            @a[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @dt * @tx1 * (c34 * tmp1) - @dt * @tx1 * @dx4
            @a[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * ((@c2 * 2.0 * @qs[(i + 1) + j * @jsize3 + k * @ksize3] - @c1 * @u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp2)) - @dt * @tx1 * (-(r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - c1345 * tmp2 * @u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @a[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@c1 * (@u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @c2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp2 + @qs[(i + 1) + j * @jsize3 + k * @ksize3] * tmp1)) - @dt * @tx1 * (r43 * c34 - c1345) * tmp2 * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @a[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-@c2 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (c34 - c1345) * tmp2 * @u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @a[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-@c2 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (c34 - c1345) * tmp2 * @u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @a[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@c1 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)) - @dt * @tx1 * c1345 * tmp1 - @dt * @tx1 * @dx5
            tmp1 = @rho_i[i + (j + 1) * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @b[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty1 * @dy1
            @b[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2
            @b[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-(@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (-c34 * tmp2 * @u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            @b[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1) - @dt * @ty1 * (c34 * tmp1) - @dt * @ty1 * @dy2
            @b[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1)
            @b[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-Math_dot_pow((@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1), 2) + @c2 * (@qs[i + (j + 1) * @jsize3 + k * @ksize3] * tmp1)) - @dt * @ty1 * (-r43 * c34 * tmp2 * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            @b[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-@c2 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1))
            @b[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * ((2.0 - @c2) * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1)) - @dt * @ty1 * (r43 * c34 * tmp1) - @dt * @ty1 * @dy3
            @b[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-@c2 * (@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1))
            @b[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * @c2
            @b[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-(@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (-c34 * tmp2 * @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            @b[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1)
            @b[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1) - @dt * @ty1 * (c34 * tmp1) - @dt * @ty1 * @dy4
            @b[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * ((@c2 * 2.0 * @qs[i + (j + 1) * @jsize3 + k * @ksize3] - @c1 * @u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp2)) - @dt * @ty1 * (-(c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1], 2) - (r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1], 2) - c1345 * tmp2 * @u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            @b[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-@c2 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (c34 - c1345) * tmp2 * @u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]
            @b[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@c1 * (@u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1) - @c2 * (@qs[i + (j + 1) * @jsize3 + k * @ksize3] * tmp1 + @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp2)) - @dt * @ty1 * (r43 * c34 - c1345) * tmp2 * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]
            @b[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-@c2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (c34 - c1345) * tmp2 * @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]
            @b[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@c1 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1)) - @dt * @ty1 * c1345 * tmp1 - @dt * @ty1 * @dy5
            tmp1 = @rho_i[i + j * @jsize3 + (k + 1) * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @c[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz1 * @dz1
            @c[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2
            @c[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-(@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * tmp2) - @dt * @tz1 * (-c34 * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            @c[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1) - @dt * @tz1 * c34 * tmp1 - @dt * @tz1 * @dz2
            @c[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1)
            @c[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-(@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * tmp2) - @dt * @tz1 * (-c34 * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            @c[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1) - @dt * @tz1 * (c34 * tmp1) - @dt * @tz1 * @dz3
            @c[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1)
            @c[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-Math_dot_pow((@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1), 2) + @c2 * (@qs[i + j * @jsize3 + (k + 1) * @ksize3] * tmp1)) - @dt * @tz1 * (-r43 * c34 * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            @c[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-@c2 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1))
            @c[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-@c2 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1))
            @c[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (2.0 - @c2) * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1) - @dt * @tz1 * (r43 * c34 * tmp1) - @dt * @tz1 * @dz4
            @c[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * @c2
            @c[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * ((@c2 * 2.0 * @qs[i + j * @jsize3 + (k + 1) * @ksize3] - @c1 * @u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp2)) - @dt * @tz1 * (-(c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1], 2) - (r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1], 2) - c1345 * tmp2 * @u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            @c[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-@c2 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * tmp2) - @dt * @tz1 * (c34 - c1345) * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]
            @c[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-@c2 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * tmp2) - @dt * @tz1 * (c34 - c1345) * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]
            @c[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@c1 * (@u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1) - @c2 * (@qs[i + j * @jsize3 + (k + 1) * @ksize3] * tmp1 + @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp2)) - @dt * @tz1 * (r43 * c34 - c1345) * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]
            @c[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@c1 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1)) - @dt * @tz1 * c1345 * tmp1 - @dt * @tz1 * @dz5
          i -= 1
        end
      j -= 1
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    j = @upper_bound1
    while j >= @lower_bound1 do
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
        i = @iend - 1
        while i >= @ist - 1 do
            for m in 0..4 do
                @tv[m + i * @isize1 + j * @jsize1] = @omega * (@udz[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * @v[0 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @udz[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * @v[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @udz[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * @v[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @udz[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * @v[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @udz[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * @v[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            end
          i -= 1
        end
      j -= 1
    end
  end

  def step2(k)
    
    
    
    
    
    tmat = Array.new(5 * 5, 0.0); 
    tmp = 1.0 / (@omega * (2.0 - @omega))
    r43 = (4.0 / 3.0)
    c1345 = @c1 * @c3 * @c4 * @c5
    c34 = @c3 * @c4
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    j = @upper_bound1
    while j >= @lower_bound1 do
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
        i = @iend - 1
        while i >= @ist - 1 do
            for m in 0..4 do
                @tv[m + i * @isize1 + j * @jsize1] = @tv[m + i * @isize1 + j * @jsize1] + @omega * (@udy[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * @v[0 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @udx[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * @v[0 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @udy[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * @v[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @udx[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * @v[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @udy[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * @v[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @udx[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * @v[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @udy[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * @v[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @udx[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * @v[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @udy[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * @v[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @udx[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * @v[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            end
            for m in 0..4 do
                tmat[m + 0 * 5] = @d[m + 0 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 1 * 5] = @d[m + 1 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 2 * 5] = @d[m + 2 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 3 * 5] = @d[m + 3 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 4 * 5] = @d[m + 4 * @isize4 + i * @jsize4 + j * @ksize4]
            end
            tmp1 = 1.0 / tmat[0 + 0 * 5]
            tmp = tmp1 * tmat[1 + 0 * 5]
            tmat[1 + 1 * 5] = tmat[1 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[1 + 2 * 5] = tmat[1 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[1 + 3 * 5] = tmat[1 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[1 + 4 * 5] = tmat[1 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            @tv[1 + i * @isize1 + j * @jsize1] = @tv[1 + i * @isize1 + j * @jsize1] - @tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[2 + 0 * 5]
            tmat[2 + 1 * 5] = tmat[2 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[2 + 2 * 5] = tmat[2 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[2 + 3 * 5] = tmat[2 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[2 + 4 * 5] = tmat[2 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            @tv[2 + i * @isize1 + j * @jsize1] = @tv[2 + i * @isize1 + j * @jsize1] - @tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[3 + 0 * 5]
            tmat[3 + 1 * 5] = tmat[3 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[3 + 2 * 5] = tmat[3 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            @tv[3 + i * @isize1 + j * @jsize1] = @tv[3 + i * @isize1 + j * @jsize1] - @tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 0 * 5]
            tmat[4 + 1 * 5] = tmat[4 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[4 + 2 * 5] = tmat[4 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            @tv[4 + i * @isize1 + j * @jsize1] = @tv[4 + i * @isize1 + j * @jsize1] - @tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[1 + 1 * 5]
            tmp = tmp1 * tmat[2 + 1 * 5]
            tmat[2 + 2 * 5] = tmat[2 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[2 + 3 * 5] = tmat[2 + 3 * 5] - tmp * tmat[1 + 3 * 5]
            tmat[2 + 4 * 5] = tmat[2 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            @tv[2 + i * @isize1 + j * @jsize1] = @tv[2 + i * @isize1 + j * @jsize1] - @tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[3 + 1 * 5]
            tmat[3 + 2 * 5] = tmat[3 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[3 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            @tv[3 + i * @isize1 + j * @jsize1] = @tv[3 + i * @isize1 + j * @jsize1] - @tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 1 * 5]
            tmat[4 + 2 * 5] = tmat[4 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[1 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            @tv[4 + i * @isize1 + j * @jsize1] = @tv[4 + i * @isize1 + j * @jsize1] - @tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[2 + 2 * 5]
            tmp = tmp1 * tmat[3 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[2 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[2 + 4 * 5]
            @tv[3 + i * @isize1 + j * @jsize1] = @tv[3 + i * @isize1 + j * @jsize1] - @tv[2 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[2 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[2 + 4 * 5]
            @tv[4 + i * @isize1 + j * @jsize1] = @tv[4 + i * @isize1 + j * @jsize1] - @tv[2 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[3 + 3 * 5]
            tmp = tmp1 * tmat[4 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[3 + 4 * 5]
            @tv[4 + i * @isize1 + j * @jsize1] = @tv[4 + i * @isize1 + j * @jsize1] - @tv[3 + i * @isize1 + j * @jsize1] * tmp
            @tv[4 + i * @isize1 + j * @jsize1] = @tv[4 + i * @isize1 + j * @jsize1] / tmat[4 + 4 * 5]
            @tv[3 + i * @isize1 + j * @jsize1] = @tv[3 + i * @isize1 + j * @jsize1] - tmat[3 + 4 * 5] * @tv[4 + i * @isize1 + j * @jsize1]
            @tv[3 + i * @isize1 + j * @jsize1] = @tv[3 + i * @isize1 + j * @jsize1] / tmat[3 + 3 * 5]
            @tv[2 + i * @isize1 + j * @jsize1] = @tv[2 + i * @isize1 + j * @jsize1] - tmat[2 + 3 * 5] * @tv[3 + i * @isize1 + j * @jsize1] - tmat[2 + 4 * 5] * @tv[4 + i * @isize1 + j * @jsize1]
            @tv[2 + i * @isize1 + j * @jsize1] = @tv[2 + i * @isize1 + j * @jsize1] / tmat[2 + 2 * 5]
            @tv[1 + i * @isize1 + j * @jsize1] = @tv[1 + i * @isize1 + j * @jsize1] - tmat[1 + 2 * 5] * @tv[2 + i * @isize1 + j * @jsize1] - tmat[1 + 3 * 5] * @tv[3 + i * @isize1 + j * @jsize1] - tmat[1 + 4 * 5] * @tv[4 + i * @isize1 + j * @jsize1]
            @tv[1 + i * @isize1 + j * @jsize1] = @tv[1 + i * @isize1 + j * @jsize1] / tmat[1 + 1 * 5]
            @tv[0 + i * @isize1 + j * @jsize1] = @tv[0 + i * @isize1 + j * @jsize1] - tmat[0 + 1 * 5] * @tv[1 + i * @isize1 + j * @jsize1] - tmat[0 + 2 * 5] * @tv[2 + i * @isize1 + j * @jsize1] - tmat[0 + 3 * 5] * @tv[3 + i * @isize1 + j * @jsize1] - tmat[0 + 4 * 5] * @tv[4 + i * @isize1 + j * @jsize1]
            @tv[0 + i * @isize1 + j * @jsize1] = @tv[0 + i * @isize1 + j * @jsize1] / tmat[0 + 0 * 5]
            @v[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @v[0 + i * @isize1 + j * @jsize1 + k * @ksize1] - @tv[0 + i * @isize1 + j * @jsize1]
            @v[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @v[1 + i * @isize1 + j * @jsize1 + k * @ksize1] - @tv[1 + i * @isize1 + j * @jsize1]
            @v[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @v[2 + i * @isize1 + j * @jsize1 + k * @ksize1] - @tv[2 + i * @isize1 + j * @jsize1]
            @v[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @v[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @tv[3 + i * @isize1 + j * @jsize1]
            @v[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @v[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @tv[4 + i * @isize1 + j * @jsize1]
          i -= 1
        end
      j -= 1
    end
  end

end
