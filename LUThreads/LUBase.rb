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
require "Runnable"
class LUBase < Runnable
  def clss
    return @CLASS
  end
  attr_accessor :id
  attr_accessor :num_threads
  attr_accessor :isiz1 
  attr_accessor :isiz2 
  attr_accessor :isiz3
  attr_accessor :itmax_default
  attr_accessor :dt_default
  attr_accessor :inorm_default
  attr_accessor :u
  attr_accessor :rsd
  attr_accessor :frct
  attr_accessor :isize1
  attr_accessor :jsize1
  attr_accessor :ksize1
  attr_accessor :flux
  attr_accessor :isize2
  attr_accessor :qs
  attr_accessor :rho_i
  attr_accessor :jsize3
  attr_accessor :ksize3
  attr_accessor :a
  attr_accessor :b
  attr_accessor :c
  attr_accessor :d
  attr_accessor :isize4
  attr_accessor :jsize4
  attr_accessor :ksize4
  attr_accessor :nx
  attr_accessor :ny
  attr_accessor :nz
  attr_accessor :nx0
  attr_accessor :ny0
  attr_accessor :nz0
  attr_accessor :ist
  attr_accessor :iend
  attr_accessor :jst
  attr_accessor :jend
  attr_accessor :ii1
  attr_accessor :ii2
  attr_accessor :ji1
  attr_accessor :ji2
  attr_accessor :ki1
  attr_accessor :ki2
  attr_accessor :dxi
  attr_accessor :deta
  attr_accessor :dzeta
  attr_accessor :tx1
  attr_accessor :tx2
  attr_accessor :tx3
  attr_accessor :ty1
  attr_accessor :ty2
  attr_accessor :ty3
  attr_accessor :tz1
  attr_accessor :tz2
  attr_accessor :tz3
  attr_accessor :dx1
  attr_accessor :dx2
  attr_accessor :dx3
  attr_accessor :dx4
  attr_accessor :dx5
  attr_accessor :dy1
  attr_accessor :dy2
  attr_accessor :dy3
  attr_accessor :dy4
  attr_accessor :dy5
  attr_accessor :dz1
  attr_accessor :dz2
  attr_accessor :dz3
  attr_accessor :dz4
  attr_accessor :dz5
  attr_accessor :dssp
  attr_accessor :dt
  attr_accessor :omega
  attr_accessor :frc
  attr_accessor :ttotal
  attr_accessor :done

  def initialize(cls, np)
    @BMName = "LU"
    @CLASS = 'S'
    @ipr_default = 1
    @omega_default = 1.2
    @tolrsd1_def = 0.00000001
    @tolrsd2_def = 0.00000001
    @tolrsd3_def = 0.00000001
    @tolrsd4_def = 0.00000001
    @tolrsd5_def = 0.00000001
    @c1 = 1.4
    @c2 = 0.4
    @c3 = 0.1
    @c4 = 1
    @c5 = 1.4
    @tolrsd = Array.new(5, 0.0)
    @rsdnm = Array.new(5, 0.0)
    @errnm = Array.new(5, 0.0)
    @ce = [2.0, 1.0, 2.0, 2.0, 5.0, 0.0, 0.0, 2.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 2.0, 5.0, 1.0, 0.0, 0.0, 0.1, 3.0, 2.0, 2.0, 2.0, 0.4, 0.5, 3.0, 3.0, 3.0, 0.3, 0.02, 0.01, 0.04, 0.03, 0.05, 0.01, 0.03, 0.03, 0.05, 0.04, 0.03, 0.02, 0.05, 0.04, 0.03, 0.5, 0.4, 0.3, 0.2, 0.1, 0.4, 0.3, 0.5, 0.1, 0.3, 0.3, 0.5, 0.4, 0.3, 0.2]
    @t_total = 1
    @t_rhsx = 2
    @t_rhsy = 3
    @t_rhsz = 4
    @t_rhs = 5
    @t_jacld = 6
    @t_blts = 7
    @t_jacu = 8
    @t_buts = 9
    @t_add = 10
    @t_l2norm = 11
    @t_last = 11
    @timer = NPB3_0_RUB::Timer.new()
    @CLASS = cls
    @num_threads = np
    case cls
    when 'S'
      @isiz1 = @isiz2 = @isiz3 = 12
      @itmax_default = @inorm_default = 50
      @dt_default = 0.5
    when 'W'
      @isiz1 = @isiz2 = @isiz3 = 33
      @itmax_default = @inorm_default = 300
      @dt_default = 0.0015
    when 'A'
      @isiz1 = @isiz2 = @isiz3 = 64
      @itmax_default = @inorm_default = 250
      @dt_default = 2
    when 'B'
      @isiz1 = @isiz2 = @isiz3 = 102
      @itmax_default = @inorm_default = 250
      @dt_default = 2
    when 'C'
      @isiz1 = @isiz2 = @isiz3 = 162
      @itmax_default = @inorm_default = 250
      @dt_default = 2
    end
    @u = Array.new(5 * (@isiz1 / 2 * 2 + 1) * (@isiz2 / 2 * 2 + 1) * @isiz3, 0.0)
    @rsd = Array.new(5 * (@isiz1 / 2 * 2 + 1) * (@isiz2 / 2 * 2 + 1) * @isiz3, 0.0)
    @frct = Array.new(5 * (@isiz1 / 2 * 2 + 1) * (@isiz2 / 2 * 2 + 1) * @isiz3, 0.0)
    @isize1 = 5
    @jsize1 = 5 * (@isiz1 / 2 * 2 + 1)
    @ksize1 = 5 * (@isiz1 / 2 * 2 + 1) * (@isiz2 / 2 * 2 + 1)
    @flux = Array.new(5 * @isiz1, 0.0)
    @isize2 = 5
    @qs = Array.new((@isiz1 / 2 * 2 + 1) * (@isiz2 / 2 * 2 + 1) * @isiz3, 0.0)
    @rho_i = Array.new((@isiz1 / 2 * 2 + 1) * (@isiz2 / 2 * 2 + 1) * @isiz3, 0.0)
    @jsize3 = (@isiz1 / 2 * 2 + 1)
    @ksize3 = (@isiz1 / 2 * 2 + 1) * (@isiz2 / 2 * 2 + 1)
    @a = Array.new(5 * 5 * (@isiz1 / 2 * 2 + 1) * (@isiz2), 0.0)
    @b = Array.new(5 * 5 * (@isiz1 / 2 * 2 + 1) * (@isiz2), 0.0)
    @c = Array.new(5 * 5 * (@isiz1 / 2 * 2 + 1) * (@isiz2), 0.0)
    @d = Array.new(5 * 5 * (@isiz1 / 2 * 2 + 1) * (@isiz2), 0.0)
    @isize4 = 5
    @jsize4 = 5 * 5
    @ksize4 = 5 * 5 * (@isiz1 / 2 * 2 + 1)
  end

  def checksum(array, size, arrayname, stop)
    sum = 0; 
    for i in 0..size-1 do
      sum += array[i]
    end
    puts("array:" + arrayname + " checksum is: " + sum)
    if stop then
      System.exit(0)
    else
    end
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

  def set_partition(start, interval, array)
    array[0][0] = start
    if start == 0 then
      array[0][1] = interval[0] - 1
    else
      array[0][1] = interval[0]
    end
    for i in 1..interval.length-1 do
        array[i][0] = array[i - 1][1] + 1
        array[i][1] = array[i - 1][1] + interval[i]
    end
  end

  def exact(i, j, k, u0)
    xi = (i - 1.0) / (@nx0 - 1); 
    eta = (j - 1.0) / (@ny0 - 1); 
    zeta = (k - 1.0) / (@nz - 1); 
    for m in 0..4 do
        u0[m] = @ce[m + 0 * 5] + (@ce[m + 1 * 5] + (@ce[m + 4 * 5] + (@ce[m + 7 * 5] + @ce[m + 10 * 5] * xi) * xi) * xi) * xi + (@ce[m + 2 * 5] + (@ce[m + 5 * 5] + (@ce[m + 8 * 5] + @ce[m + 11 * 5] * eta) * eta) * eta) * eta + (@ce[m + 3 * 5] + (@ce[m + 6 * 5] + (@ce[m + 9 * 5] + @ce[m + 12 * 5] * zeta) * zeta) * zeta) * zeta
    end
  end

  #def max(a, b)
  #  if a < b then
  #    return b
  #  else
  #    return a
  #  end
  #end

  def max(a, b, c=nil)
    if c != nil then
      return max(a, max(b, c))
    else
      if a < b then
        return b
      else
        return a
      end
    end
  end

end
