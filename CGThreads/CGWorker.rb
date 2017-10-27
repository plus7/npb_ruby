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
require "monitor"
class CGWorker < CGBase
  attr_accessor :id
  attr_accessor :taskOrder
  attr_accessor :alpha
  attr_accessor :beta
  attr_accessor :done

  def initialize(clss, np, cg, st, _end)
    super(clss, np, false)
    @BMName = "CG"
    @cgitmax = 25
    @CLASS = 'S'
    @t_init = 1
    @t_bench = 2
    @t_conj_grad = 3
    @t_last = 3
    @done = true
    Init(cg)
    @start1 = st
    @end1 = _end
    @done = true
    #setDaemon(true)
    #setPriority(Thread.MAX_PRIORITY)
    @master = cg
  end

  def Init(cg)
    @dmaster = cg.dmaster
    @rhomaster = cg.rhomaster
    @rnormmaster = cg.rnormmaster
    @colidx = cg.colidx
    @rowstr = cg.rowstr
    @a = cg.a
    @p = cg.p
    @q = cg.q
    @r = cg.r
    @x = cg.x
    @z = cg.z
  end

  def run()
    state = 0; 
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    while true do
         self.synchronize do
            while @done do
                begin
                    cond.wait
                    @master.synchronize do
                        @master.cond.signal
                    end
                rescue InterruptedException => ie
                ensure
                end
            end
            case @taskOrder
            when 0
              step0()
            when 1
              step1()
            when 2
              step2()
            when 3
              step3()
            when 4
              endWork()
            end
            @master.synchronize do
                @done = true
                @master.cond.signal
            end
        end
    end
  end

  def step0()
    for j in @start1..@end1 do
        sum = 0.0; 
        for k in @rowstr[j]..@rowstr[j + 1]-1 do
            sum = sum + @a[k] * @p[@colidx[k]]
        end
        @q[j] = sum
    end
    sum = 0.0; 
    for j in @start1..@end1 do
      sum += @p[j] * @q[j]
    end
    @dmaster[@id] = sum
  end

  def step1()
    for j in @start1..@end1 do
        @z[j] = @z[j] + @alpha * @p[j]
        @r[j] = @r[j] - @alpha * @q[j]
    end
    rho = 0.0; 
    for j in @start1..@end1 do
      rho += @r[j] * @r[j]
    end
    @rhomaster[@id] = rho
  end

  def step2()
    for j in @start1..@end1 do
      @p[j] = @r[j] + @beta * @p[j]
    end
  end

  def step3()
    rho = 0.0; 
    for j in @start1..@end1 do
        @q[j] = 0.0
        @z[j] = 0.0
        @r[j] = @x[j]
        @p[j] = @x[j]
        rho += @x[j] * @x[j]
    end
    @rhomaster[@id] = rho
  end

  def endWork()
    for j in @start1..@end1 do
        sum = 0.0; 
        for k in @rowstr[j]..@rowstr[j + 1] - 1 do
            sum += @a[k] * @z[@colidx[k]]
        end
        @r[j] = sum
    end
    sum = 0.0; 
    for j in @start1..@end1 do
      sum += (@x[j] - @r[j]) * (@x[j] - @r[j])
    end
    @rnormmaster[@id] = sum
  end

end
