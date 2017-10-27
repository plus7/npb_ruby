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
class FFTThread < FTBase
  def initialize(ft, low1, high1, low2, high2)
    super(ft.clss, ft.num_threads, true)
    @done = true
    Init(ft)
    @state = 1
    @lb1 = low1
    @ub1 = high1
    @lb2 = low2
    @ub2 = high2
    @plane = Array.new(2 * (@maxdim + 1) * @maxdim, 0.0)
    @scr = Array.new(2 * (@maxdim + 1) * @maxdim, 0.0)
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = ft
  end

  def Init(ft)
    @maxdim = ft.maxdim
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
            @state += 1
            if @state == 4 then
              @state = 1
            else
            end
            @master.synchronize do
                @done = true
                @master.cond.signal
            end
        end
    end
  end

  def setVariables(sign1, tr, x1, exp11, exp21, exp31)
    @sign = sign1
    @x = x1
    @exp1 = exp11
    @exp2 = exp21
    @exp3 = exp31
    @n1 = @exp1.length >> 1
    @n2 = @exp2.length >> 1
    @n3 = @exp3.length >> 1
    if tr then
        @lower_bound1 = @lb2
        @upper_bound1 = @ub2
        @lower_bound2 = @lb1
        @upper_bound2 = @ub1
    else
        @lower_bound1 = @lb1
        @upper_bound1 = @ub1
        @lower_bound2 = @lb2
        @upper_bound2 = @ub2
    end
  end

  def step()
    log = ilog2(@n2); 
    isize3 = 2; jsize3 = isize3 * (@n1 + 1); ksize3 = jsize3 * @n2; 
    isize1 = 2; 
    jsize1 = 2 * (@n2 + 1); 
    case @state
    when 1
      for k in @lower_bound1..@upper_bound1 do
        swarztrauber(@sign, log, @n1, @n2, @x, k * ksize3, @n1, @exp2, @scr)
      end
    when 2
      log = ilog2(@n1)
      for k in @lower_bound1..@upper_bound1 do
          for j in 0..@n2-1 do
              for i in 0..@n1-1 do
                  @plane[@@REAL + j * isize1 + i * jsize1] = @x[@@REAL + i * isize3 + j * jsize3 + k * ksize3]
                  @plane[@@IMAG + j * isize1 + i * jsize1] = @x[@@IMAG + i * isize3 + j * jsize3 + k * ksize3]
              end
          end
          swarztrauber(@sign, log, @n2, @n1, @plane, 0, @n2, @exp1, @scr)
          for j in 0..@n2-1 do
              for i in 0..@n1-1 do
                  @x[@@REAL + i * isize3 + j * jsize3 + k * ksize3] = @plane[@@REAL + j * isize1 + i * jsize1]
                  @x[@@IMAG + i * isize3 + j * jsize3 + k * ksize3] = @plane[@@IMAG + j * isize1 + i * jsize1]
              end
          end
      end
    when 3
      log = ilog2(@n3)
      jsize1 = 2 * (@n1 + 1)
      for k in @lower_bound2..@upper_bound2 do
          for i in 0..@n3-1 do
              for j in 0..@n1-1 do
                  @plane[@@REAL + j * isize1 + i * jsize1] = @x[@@REAL + j * isize3 + k * jsize3 + i * ksize3]
                  @plane[@@IMAG + j * isize1 + i * jsize1] = @x[@@IMAG + j * isize3 + k * jsize3 + i * ksize3]
              end
          end
          swarztrauber(@sign, log, @n1, @n3, @plane, 0, @n1, @exp3, @scr)
          for i in 0..@n3-1 do
              for j in 0..@n1-1 do
                  @x[@@REAL + j * isize3 + k * jsize3 + i * ksize3] = @plane[@@REAL + j * isize1 + i * jsize1]
                  @x[@@IMAG + j * isize3 + k * jsize3 + i * ksize3] = @plane[@@IMAG + j * isize1 + i * jsize1]
              end
          end
      end
    end
  end

# *** public ***

  attr_accessor :id

  attr_accessor :done

end
