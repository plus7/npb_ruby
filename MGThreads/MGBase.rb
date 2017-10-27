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
class MGBase < Runnable
  @@MGBaseConstInit = false
  def initialize(clss, np, serial)
    @CLASS = 'S'
    @uoff = 0.0 # ad hoc fix
    if not @@MGBaseConstInit then
      @@MGBaseConstInit = true
      @@BMName = "MG"
      @@maxlevel = 11
      @@nit = 0
      @@T_total = 0
      @@T_init = 1
      @@T_bench = 2
      @@T_mg3P = 3
      @@T_psinv = 4
      @@T_resid = 5
      @@T_resid2 = 6
      @@T_rprj3 = 7
      @@T_interp = 8
      @@T_norm2 = 9
      @@T_last = 9
    end
    @timeron = false
    @timer = NPB3_0_RUB::Timer.new()
    @num_threads = 0
    @CLASS = clss
    @num_threads = np
    @nx = Array.new(@@maxlevel, 0.0)
    @ny = Array.new(@@maxlevel, 0.0)
    @nz = Array.new(@@maxlevel, 0.0)
    @ir = Array.new(@@maxlevel, 0.0)
    @m1 = Array.new(@@maxlevel, 0.0)
    @m2 = Array.new(@@maxlevel, 0.0)
    @m3 = Array.new(@@maxlevel, 0.0)
    case @CLASS
    when 'S'
      @nx_default = 32
      @ny_default = 32
      @nz_default = 32
      @nit_default = 4
      @lm = 5
      @lt_default = 5
      @ndim1 = 5
      @ndim2 = 5
      @ndim3 = 5
      @lt = @lt_default
      @@nit = @nit_default
      @nx[@lt - 1] = @nx_default
      @ny[@lt - 1] = @ny_default
      @nz[@lt - 1] = @nz_default
    when 'W'
      @nx_default = 64
      @ny_default = 64
      @nz_default = 64
      @nit_default = 40
      @lm = 6
      @lt_default = 6
      @ndim1 = 6
      @ndim2 = 6
      @ndim3 = 6
      @lt = @lt_default
      @@nit = @nit_default
      @nx[@lt - 1] = @nx_default
      @ny[@lt - 1] = @ny_default
      @nz[@lt - 1] = @nz_default
    when 'A'
      @nx_default = 256
      @ny_default = 256
      @nz_default = 256
      @nit_default = 4
      @lm = 8
      @lt_default = 8
      @ndim1 = 8
      @ndim2 = 8
      @ndim3 = 8
      @lt = @lt_default
      @@nit = @nit_default
      @nx[@lt - 1] = @nx_default
      @ny[@lt - 1] = @ny_default
      @nz[@lt - 1] = @nz_default
    when 'B'
      @nx_default = 256
      @ny_default = 256
      @nz_default = 256
      @nit_default = 20
      @lm = 8
      @lt_default = 8
      @ndim1 = 8
      @ndim2 = 8
      @ndim3 = 8
      @lt = @lt_default
      @@nit = @nit_default
      @nx[@lt - 1] = @nx_default
      @ny[@lt - 1] = @ny_default
      @nz[@lt - 1] = @nz_default
    when 'C'
      @nx_default = 512
      @ny_default = 512
      @nz_default = 512
      @nit_default = 20
      @lm = 9
      @lt_default = 9
      @ndim1 = 9
      @ndim2 = 9
      @ndim3 = 9
      @lt = @lt_default
      @@nit = @nit_default
      @nx[@lt - 1] = @nx_default
      @ny[@lt - 1] = @ny_default
      @nz[@lt - 1] = @nz_default
    end
    @nm = 2 + (1 << @lm)
    @nv = (2 + (1 << @ndim1)) * (2 + (1 << @ndim2)) * (2 + (1 << @ndim3))
    @nr = (8 * (@nv + @nm * @nm + 5 * @nm + 7 * @lm)) / 7
    @nm2 = 2 * @nm * @nm
    @r = Array.new(@nr, 0.0)
    @v = Array.new(@nv, 0.0)
    @u = Array.new(@nr, 0.0)
    @a = Array.new(4, 0.0)
    @c = Array.new(4, 0.0)
    @a[0] = -8.0 / 3.0
    @a[1] = 0.0
    @a[2] = 1.0 / 6.0
    @a[3] = 1.0 / 12.0
    if @CLASS == 'A' or @CLASS == 'S' or @CLASS == 'W' then
        @c[0] = -3.0 / 8.0
        @c[1] = +1.0 / 32.0
        @c[2] = -1.0 / 64.0
        @c[3] = 0.0
    else
        @c[0] = -3.0 / 17.0
        @c[1] = +1.0 / 33.0
        @c[2] = -1.0 / 61.0
        @c[3] = 0.0
    end
  end

  def checksum(arr, name, stop)
    csum = 0; 
    for i in 0..arr.length-1 do
      csum += arr[i]
    end
    puts(name + " checksum MG " + csum)
    if stop then
      System.exit(0)
    else
    end
  end

  def dmax1(a, b)
    if a < b then
      return b
    else
      return a
    end
  end

  def comm3(u, off, n1, n2, n3)
    
    for i3 in 1..n3 - 1-1 do
      for i2 in 1..n2 - 1-1 do
          u[off + n1 * (i2 + n2 * i3)] = u[off + n1 - 2 + n1 * (i2 + n2 * i3)]
          u[off + n1 - 1 + n1 * (i2 + n2 * i3)] = u[off + 1 + n1 * (i2 + n2 * i3)]
      end
    end
    for i3 in 1..n3 - 1-1 do
      for i1 in 0..n1-1 do
          u[off + i1 + n1 * n2 * i3] = u[off + i1 + n1 * (n2 - 2 + n2 * i3)]
          u[off + i1 + n1 * (n2 - 1 + n2 * i3)] = u[off + i1 + n1 * (1 + n2 * i3)]
      end
    end
    for i2 in 0..n2-1 do
      for i1 in 0..n1-1 do
          u[off + i1 + n1 * i2] = u[off + i1 + n1 * (i2 + n2 * (n3 - 2))]
          u[off + i1 + n1 * (i2 + n2 * (n3 - 1))] = u[off + i1 + n1 * (i2 + n2)]
      end
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

  def maxlevel
    return @@maxlevel
  end
  def maxlevel=(val)
    @@maxlevel = val
  end

  def nit
    return @@nit
  end
  def nit=(val)
    @@nit = val
  end

  attr_accessor :nx_default

  attr_accessor :ny_default

  attr_accessor :nz_default

  attr_accessor :nit_default

  attr_accessor :lm

  attr_accessor :lt_default

  attr_accessor :ndim1

  attr_accessor :ndim2

  attr_accessor :ndim3

  attr_accessor :nm

  attr_accessor :nv

  attr_accessor :nr

  attr_accessor :nm2

  attr_accessor :nx

  attr_accessor :ny

  attr_accessor :nz

  attr_accessor :ir

  attr_accessor :m1

  attr_accessor :m2

  attr_accessor :m3

  attr_accessor :lt

  attr_accessor :lb

  attr_accessor :u

  attr_accessor :v

  attr_accessor :r

  attr_accessor :a

  attr_accessor :c

  attr_accessor :zoff

  attr_accessor :zsize3

  attr_accessor :zsize2

  attr_accessor :zsize1

  attr_accessor :uoff

  attr_accessor :usize1

  attr_accessor :usize2

  attr_accessor :usize3

  attr_accessor :roff

  attr_accessor :rsize1

  attr_accessor :rsize2

  attr_accessor :rsize3

  attr_accessor :timeron

  attr_accessor :timer

  attr_accessor :num_threads

  attr_accessor :wstart

  attr_accessor :wend

  attr_accessor :interp

  attr_accessor :psinv

  attr_accessor :rprj

  attr_accessor :resid

  attr_accessor :master

# *** protected ***
  def T_total
    return @@T_total
  end
  def T_total=(val)
    @@T_total = val
  end

  protected :T_total
  def T_init
    return @@T_init
  end
  def T_init=(val)
    @@T_init = val
  end

  protected :T_init
  def T_bench
    return @@T_bench
  end
  def T_bench=(val)
    @@T_bench = val
  end

  protected :T_bench
  def T_mg3P
    return @@T_mg3P
  end
  def T_mg3P=(val)
    @@T_mg3P = val
  end

  protected :T_mg3P
  def T_psinv
    return @@T_psinv
  end
  def T_psinv=(val)
    @@T_psinv = val
  end

  protected :T_psinv
  def T_resid
    return @@T_resid
  end
  def T_resid=(val)
    @@T_resid = val
  end

  protected :T_resid
  def T_resid2
    return @@T_resid2
  end
  def T_resid2=(val)
    @@T_resid2 = val
  end

  protected :T_resid2
  def T_rprj3
    return @@T_rprj3
  end
  def T_rprj3=(val)
    @@T_rprj3 = val
  end

  protected :T_rprj3
  def T_interp
    return @@T_interp
  end
  def T_interp=(val)
    @@T_interp = val
  end

  protected :T_interp
  def T_norm2
    return @@T_norm2
  end
  def T_norm2=(val)
    @@T_norm2 = val
  end

  protected :T_norm2
  def T_last
    return @@T_last
  end
  def T_last=(val)
    @@T_last = val
  end

  protected :T_last
end
