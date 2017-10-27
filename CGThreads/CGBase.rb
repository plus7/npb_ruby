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
class CGBase < Runnable
  attr_accessor :dmaster
  attr_accessor :rhomaster 
  attr_accessor :rnormmaster 
  attr_accessor :colidx
  attr_accessor :rowstr 
  attr_accessor :a
  attr_accessor :p
  attr_accessor :q
  attr_accessor :r
  attr_accessor :x
  attr_accessor :z

  def initialize(clss, np, serial)
    @BMName = "CG"
    @cgitmax = 25
    @CLASS = 'S'
    @t_init = 1
    @t_bench = 2
    @t_conj_grad = 3
    @t_last = 3
    @CLASS = clss
    @num_threads = np
    case @CLASS
    when 'S'
      @na = 1400
      @nonzer = 7
      @shift = 10
      @niter = 15
      @rcond = 0.1
      @zeta_verify_value = 8.5971775078648
    when 'W'
      @na = 7000
      @nonzer = 8
      @shift = 12
      @niter = 15
      @rcond = 0.1
      @zeta_verify_value = 10.362595087124
    when 'A'
      @na = 14000
      @nonzer = 11
      @shift = 20
      @niter = 15
      @rcond = 0.1
      @zeta_verify_value = 17.130235054029
    when 'B'
      @na = 75000
      @nonzer = 13
      @shift = 60
      @niter = 75
      @rcond = 0.1
      @zeta_verify_value = 22.712745482631
    when 'C'
      @na = 150000
      @nonzer = 15
      @shift = 110
      @niter = 75
      @rcond = 0.1
      @zeta_verify_value = 28.973605592845
    end
    @t_names = Array.new(@t_last + 1, 0.0)
    @timer = NPB3_0_RUB::Timer.new()
    @nz = (@na * (@nonzer + 1) * (@nonzer + 1) + @na * (@nonzer + 2))
    @colidx = Array.new(@nz + 1, 0)
    @rowstr = Array.new(@na + 2, 0)
    @iv = Array.new(2 * @na + 2, 0)
    @arow = Array.new(@nz + 1, 0)
    @acol = Array.new(@nz + 1, 0)
    @v = Array.new(@na + 2, 0.0)
    @aelt = Array.new(@nz + 1, 0.0)
    @a = Array.new(@nz + 1, 0.0)
    @p = Array.new(@na + 3, 0.0)
    @q = Array.new(@na + 3, 0.0)
    @r = Array.new(@na + 3, 0.0)
    @x = Array.new(@na + 3, 0.0)
    @z = Array.new(@na + 3, 0.0)
  end

end
