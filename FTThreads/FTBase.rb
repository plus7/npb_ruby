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
require "Runnable"
class FTBase < Runnable
  @@FTBaseStaticInit = false
  def initialize(clss, np, serial)
    @CLASS = 'S'
    @timeron = false
    @timer = NPB3_0_RUB::Timer.new()
    if not @@FTBaseStaticInit then
      @@FTBaseStaticInit = true
      @@BMName = "FT"
      @@REAL = 0
      @@IMAG = 1
      @@pi = Math::PI
      @@alpha = 0.000001
      @@fftblock_default = 4 * 4096
      @@fftblock = 0
    end
    @master =  nil 
    @CLASS = clss
    @num_threads = np
    case @CLASS
    when 'S'
      @nx = @ny = @nz = 64
      @niter_default = 6
    when 'W'
      @nx = @ny = 128
      @nz = 32
      @niter_default = 6
    when 'A'
      @nx = 256
      @ny = 256
      @nz = 128
      @niter_default = 6
    when 'B'
      @nx = 512
      @ny = @nz = 256
      @niter_default = 20
    when 'C'
      @nx = @ny = @nz = 512
      @niter_default = 20
    end
    @maxdim = max(@nx, max(@ny, @nx))
    if serial then
        @scr = Array.new(2 * (@maxdim + 1) * @maxdim, 0.0)
        @plane = Array.new(2 * (@maxdim + 1) * @maxdim, 0.0)
    else
    end
    @isize2 = 2
    @isize3 = 2
    @jsize3 = 2 * (@ny + 1)
    @ksize3 = 2 * (@ny + 1) * @nx
    @isize4 = 2
    @jsize4 = 2 * (@ny + 1)
    @ksize4 = 2 * (@ny + 1) * @nz
    @checksum = Array.new(2 * @niter_default, 0.0)
    @xtr = Array.new(2 * (@ny + 1) * @nx * @nz, 0.0)
    @xnt = Array.new(2 * (@ny + 1) * @nz * @nx, 0.0)
    @exp1 = Array.new(2 * @nx, 0.0)
    @exp2 = Array.new(2 * @ny, 0.0)
    @exp3 = Array.new(2 * @nz, 0.0)
  end

  def set_interval(threads, problem_size, interval)
    interval[0] = problem_size / threads
    for i in 1..threads-1 do
      interval[i] = interval[0]
    end
    remainder = problem_size % threads; 
    for i in 0..remainder-1 do
      interval[i] += 1
    end
  end

  def set_partition(start, interval, prt)
    prt[0][0] = start
    if start == 0 then
      prt[0][1] = interval[0] - 1
    else
      prt[0][1] = interval[0]
    end
    for i in 1..interval.length-1 do
        prt[i][0] = prt[i - 1][1] + 1
        prt[i][1] = prt[i - 1][1] + interval[i]
    end
  end

  def max(a, b)
    if a > b then
      return a
    else
      return b
    end
  end

  def CompExp(n, exponent)
    nu = n; 
    m = ilog2(n); 
    exponent[0] = m
    eps = 1.0E-16; 
    ku = 1; 
    ln = 1; 
    for j in 1..m do
        t = @@pi / ln; 
        for i in 0..ln - 1 do
            ti = i * t; 
            idx = (i + ku) * 2; 
            exponent[@@REAL + idx] = Math.cos(ti)
            exponent[@@IMAG + idx] = Math.sin(ti)
            if Math_dot_abs(exponent[@@REAL + idx]) < eps then
              exponent[@@REAL + idx] = 0
            else
            end
            if Math_dot_abs(exponent[@@IMAG + idx]) < eps then
              exponent[@@IMAG + idx] = 0
            else
            end
        end
        ku = ku + ln
        ln = 2 * ln
    end
  end

  def initial_conditions(u0, d1, d2, d3)
    tmp = Array.new(2 * @maxdim, 0.0); 
    ranStarts = Array.new(@maxdim, 0.0); 
    seed = 314159265; a = Math_dot_pow(5.0, 13); 
    start = seed; 
    rng = NPB3_0_RUB::Random.new(seed); 
    an = rng.ipow(a, 0); #ipow46
    rng.randlc(seed, an)
    an = rng.ipow(a, 2 * d1 * d2) #ipow46
    ranStarts[0] = start
    for k in 1..d3-1 do
        seed = rng.randlc(start, an)
        ranStarts[k] = start = seed
    end
    for k in 0..d3-1 do
        x0 = ranStarts[k]; 
        for j in 0..d1-1 do
            x0 = rng.vranlc(2 * d2, x0, a, tmp, 0)
            for i in 0..d2-1 do
                u0[@@REAL + j * @isize3 + i * @jsize3 + k * @ksize3] = tmp[@@REAL + i * 2]
                u0[@@IMAG + j * @isize3 + i * @jsize3 + k * @ksize3] = tmp[@@IMAG + i * 2]
            end
        end
    end
  end

  def ilog2(n)
    
    if n == 1 then
      return 0
    else
    end
    lg = 1
    nn = 2
    while nn < n do
        nn = nn * 2
        lg = lg + 1
    end
    return lg
  end

  def swarztrauber(is, m, len, n, x, xoffst, xd1, exponent, scr)
    j = 0; 
    
    
    isize1 = 2; jsize1 = 2 * (xd1 + 1); 
    u1 = Array.new(2, 0.0); x11 = Array.new(2, 0.0); x21 = Array.new(2, 0.0); 
    if @timeron then
      @timer.start(4)
    else
    end
    @@fftblock = @@fftblock_default / n
    if @@fftblock < 8 then
      @@fftblock = 8
    else
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    blockStart = 0
    while blockStart < len do
        blockEnd = blockStart + @@fftblock - 1
        if blockEnd >= len then
          blockEnd = len - 1
        else
        end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
        l = 1
        while l <= m do
            n1 = n / 2
            lk = Math_dot_pow(2, l - 1)
            li = Math_dot_pow(2, m - l)
            lj = 2 * lk
            ku = li
            for i in 0..li - 1 do
                i11 = i * lk
                i12 = i11 + n1
                i21 = i * lj
                i22 = i21 + lk
                u1[@@REAL] = exponent[@@REAL + (ku + i) * 2]
                if is >= 1 then
                    u1[@@IMAG] = exponent[@@IMAG + (ku + i) * 2]
                else
                    u1[@@IMAG] = -exponent[@@IMAG + (ku + i) * 2]
                end
                for k in 0..lk - 1 do
                    for j in blockStart..blockEnd do
                        x11[@@REAL] = x[@@REAL + j * 2 + (i11 + k) * jsize1 + xoffst]
                        x11[@@IMAG] = x[@@IMAG + j * 2 + (i11 + k) * jsize1 + xoffst]
                        x21[@@REAL] = x[@@REAL + j * 2 + (i12 + k) * jsize1 + xoffst]
                        x21[@@IMAG] = x[@@IMAG + j * 2 + (i12 + k) * jsize1 + xoffst]
                        scr[@@REAL + j * isize1 + (i21 + k) * jsize1] = x11[@@REAL] + x21[@@REAL]
                        scr[@@IMAG + j * isize1 + (i21 + k) * jsize1] = x11[@@IMAG] + x21[@@IMAG]
                        scr[@@REAL + j * 2 + (i22 + k) * jsize1] = u1[@@REAL] * (x11[@@REAL] - x21[@@REAL]) - u1[@@IMAG] * (x11[@@IMAG] - x21[@@IMAG])
                        scr[@@IMAG + j * 2 + (i22 + k) * jsize1] = u1[@@IMAG] * (x11[@@REAL] - x21[@@REAL]) + u1[@@REAL] * (x11[@@IMAG] - x21[@@IMAG])
                    end
                end
            end
            if l == m then
                for k in 0..n-1 do
                    for j in blockStart..blockEnd do
                        x[@@REAL + j * 2 + k * jsize1 + xoffst] = scr[@@REAL + j * isize1 + k * jsize1]
                        x[@@IMAG + j * 2 + k * jsize1 + xoffst] = scr[@@IMAG + j * isize1 + k * jsize1]
                    end
                end
            else
                n1 = n / 2
                lk = Math_dot_pow(2, l)
                li = Math_dot_pow(2, m - l - 1)
                lj = 2 * lk
                ku = li
                for i in 0..li - 1 do
                    i11 = i * lk
                    i12 = i11 + n1
                    i21 = i * lj
                    i22 = i21 + lk
                    u1[@@REAL] = exponent[@@REAL + (ku + i) * 2]
                    if is >= 1 then
                        u1[@@IMAG] = exponent[@@IMAG + (ku + i) * 2]
                    else
                        u1[@@IMAG] = -exponent[@@IMAG + (ku + i) * 2]
                    end
                    for k in 0..lk - 1 do
                        for j in blockStart..blockEnd do
                            x11[@@REAL] = scr[@@REAL + j * isize1 + (i11 + k) * jsize1]
                            x11[@@IMAG] = scr[@@IMAG + j * isize1 + (i11 + k) * jsize1]
                            x21[@@REAL] = scr[@@REAL + j * isize1 + (i12 + k) * jsize1]
                            x21[@@IMAG] = scr[@@IMAG + j * isize1 + (i12 + k) * jsize1]
                            x[@@REAL + j * 2 + (i21 + k) * jsize1 + xoffst] = x11[@@REAL] + x21[@@REAL]
                            x[@@IMAG + j * 2 + (i21 + k) * jsize1 + xoffst] = x11[@@IMAG] + x21[@@IMAG]
                            x[@@REAL + j * 2 + (i22 + k) * jsize1 + xoffst] = u1[@@REAL] * (x11[@@REAL] - x21[@@REAL]) - u1[@@IMAG] * (x11[@@IMAG] - x21[@@IMAG])
                            x[@@IMAG + j * 2 + (i22 + k) * jsize1 + xoffst] = u1[@@IMAG] * (x11[@@REAL] - x21[@@REAL]) + u1[@@REAL] * (x11[@@IMAG] - x21[@@IMAG])
                        end
                    end
                end
            end
          l += 2
        end
      blockStart += @@fftblock
    end
    if @timeron then
      @timer.stop(4)
    else
    end
  end

# *** public ***
  def BMName
    return @@BMName
  end
  def BMName=(val)
    @@BMName = val
  end

  def clss
    return @CLASS
  end

  attr_accessor :timeron

  attr_accessor :timer

# *** protected ***
  attr_accessor :nx

  protected :nx
  attr_accessor :ny

  protected :ny
  attr_accessor :nz

  protected :nz
  attr_accessor :maxdim

  protected :maxdim
  attr_accessor :niter_default

  protected :niter_default
  attr_accessor :scr

  protected :scr
  attr_accessor :plane

  protected :plane
  attr_accessor :isize2

  protected :isize2
  attr_accessor :isize3

  protected :isize3
  attr_accessor :jsize3

  protected :jsize3
  attr_accessor :ksize3

  protected :ksize3
  attr_accessor :isize4

  protected :isize4
  attr_accessor :jsize4

  protected :jsize4
  attr_accessor :ksize4

  protected :ksize4
  attr_accessor :checksum

  protected :checksum
  attr_accessor :xtr

  protected :xtr
  attr_accessor :xnt

  protected :xnt
  attr_accessor :exp1

  protected :exp1
  attr_accessor :exp2

  protected :exp2
  attr_accessor :exp3

  protected :exp3
  def real
    return @@REAL
  end
  def real=(val)
    @@REAL = val
  end

  protected :real
  def imag
    return @@IMAG
  end
  def imag=(val)
    @@IMAG = val
  end

  protected :imag
  def pi
    return @@pi
  end
  def pi=(val)
    @@pi = val
  end

  protected :pi
  def alpha
    return @@alpha
  end
  def alpha=(val)
    @@alpha = val
  end

  protected :alpha
  attr_accessor :master

  protected :master
  attr_accessor :num_threads

  protected :num_threads
  attr_accessor :doFFT

  protected :doFFT
  attr_accessor :doEvolve

  protected :doEvolve
  def fftblock_default
    return @@fftblock_default
  end
  def fftblock_default=(val)
    @@fftblock_default = val
  end

  protected :fftblock_default
  def fftblock
    return @@fftblock
  end
  def fftblock=(val)
    @@fftblock = val
  end

  protected :fftblock
end
