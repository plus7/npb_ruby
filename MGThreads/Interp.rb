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
class Interp < MGBase
  def initialize(mg)
    @done = true
    @state = 0
    Init(mg)
    m = 535; 
    @z1 = Array.new(m, 0.0)
    @z2 = Array.new(m, 0.0)
    @z3 = Array.new(m, 0.0)
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = mg
  end

  def Init(mg)
    @num_threads = mg.num_threads
    @u = mg.u
  end

  def run()
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    while true do
         self .synchronize do
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
            GetWork()
            step()
            @master.synchronize do
                @done = true
                @master.cond.signal
            end
        end
    end
  end

  def step()
    if @work == 0 then
      return 
    else
    end
    
    if @n1 != 3 and @n2 != 3 and @n3 != 3 then
        for i3 in @start..@end do
            for i2 in 1..@mm2 - 1 do
                for i1 in 1..@mm1 do
                    @z1[i1 - 1] = @u[@zoff + i1 - 1 + @mm1 * (i2 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))]
                    @z2[i1 - 1] = @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * i3)] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))]
                    @z3[i1 - 1] = @u[@zoff + i1 - 1 + @mm1 * (i2 + @mm2 * i3)] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * i3)] + @z1[i1 - 1]
                end
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 2 + @n1 * (2 * i2 - 2 + @n2 * (2 * i3 - 2))] += @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))]
                    @u[@uoff + 2 * i1 - 1 + @n1 * (2 * i2 - 2 + @n2 * (2 * i3 - 2))] += 0.5 * (@u[@zoff + i1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))])
                end
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 2 + @n1 * (2 * i2 - 1 + @n2 * (2 * i3 - 2))] += 0.5 * @z1[i1 - 1]
                    @u[@uoff + 2 * i1 - 1 + @n1 * (2 * i2 - 1 + @n2 * (2 * i3 - 2))] += 0.25 * (@z1[i1 - 1] + @z1[i1])
                end
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 2 + @n1 * (2 * i2 - 2 + @n2 * (2 * i3 - 1))] += 0.5 * @z2[i1 - 1]
                    @u[@uoff + 2 * i1 - 1 + @n1 * (2 * i2 - 2 + @n2 * (2 * i3 - 1))] += 0.25 * (@z2[i1 - 1] + @z2[i1])
                end
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 2 + @n1 * (2 * i2 - 1 + @n2 * (2 * i3 - 1))] += 0.25 * @z3[i1 - 1]
                    @u[@uoff + 2 * i1 - 1 + @n1 * (2 * i2 - 1 + @n2 * (2 * i3 - 1))] += 0.125 * (@z3[i1 - 1] + @z3[i1])
                end
            end
        end
    else
        if @n1 == 3 then
            d1 = 2
            t1 = 1
        else
            d1 = 1
            t1 = 0
        end
        if @n2 == 3 then
            d2 = 2
            t2 = 1
        else
            d2 = 1
            t2 = 0
        end
        if @n3 == 3 then
            d3 = 2
            t3 = 1
        else
            d3 = 1
            t3 = 0
        end
        for i3 in @start..@end do
            for i2 in 1..@mm2 - 1 do
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 1 - d1 + @n1 * (2 * i2 - 1 - d2 + @n2 * (2 * i3 - 1 - d3))] += @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))]
                end
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 1 - t1 + @n1 * (2 * i2 - 1 - d2 + @n2 * (2 * i3 - 1 - d3))] += 0.5 * (@u[@zoff + i1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))])
                end
            end
            for i2 in 1..@mm2 - 1 do
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 1 - d1 + @n1 * (2 * i2 - 1 - t2 + @n2 * (2 * i3 - 1 - d3))] += 0.5 * (@u[@zoff + i1 - 1 + @mm1 * (i2 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))])
                end
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 1 - t1 + @n1 * (2 * i2 - 1 - t2 + @n2 * (2 * i3 - 1 - d3))] += 0.25 * (@u[@zoff + i1 + @mm1 * (i2 + @mm2 * (i3 - 1))] + @u[@zoff + i1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))])
                end
            end
        end
        for i3 in @start..@end do
            for i2 in 1..@mm2 - 1 do
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 1 - d1 + @n1 * (2 * i2 - 1 - d2 + @n2 * (2 * i3 - 1 - t3))] = 0.5 * (@u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * i3)] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))])
                end
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 1 - t1 + @n1 * (2 * i2 - 1 - d2 + @n2 * (2 * i3 - 1 - t3))] += 0.25 * (@u[@zoff + i1 + @mm1 * (i2 - 1 + @mm2 * i3)] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * i3)] + @u[@zoff + i1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))])
                end
            end
            for i2 in 1..@mm2 - 1 do
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 1 - d1 + @n1 * (2 * i2 - 1 - t2 + @n2 * (2 * i3 - 1 - t3))] += 0.25 * (@u[@zoff + i1 - 1 + @mm1 * (i2 + @mm2 * i3)] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * i3)] + @u[@zoff + i1 - 1 + @mm1 * (i2 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))])
                end
                for i1 in 1..@mm1 - 1 do
                    @u[@uoff + 2 * i1 - 1 - t1 + @n1 * (2 * i2 - 1 - t2 + @n2 * (2 * i3 - 1 - t3))] += 0.125 * (@u[@zoff + i1 + @mm1 * (i2 + @mm2 * i3)] + @u[@zoff + i1 + @mm1 * (i2 - 1 + @mm2 * i3)] + @u[@zoff + i1 - 1 + @mm1 * (i2 + @mm2 * i3)] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * i3)] + @u[@zoff + i1 + @mm1 * (i2 + @mm2 * (i3 - 1))] + @u[@zoff + i1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 + @mm2 * (i3 - 1))] + @u[@zoff + i1 - 1 + @mm1 * (i2 - 1 + @mm2 * (i3 - 1))])
                end
            end
        end
    end
  end

  def step2()
    if @work == 0 then
      return 
    else
    end
  end

  def GetWork()
    workpt = (@wend - @wstart) / @num_threads; 
    remainder = @wend - @wstart - workpt * @num_threads; 
    if workpt == 0 then
        if @id < @wend - @wstart then
            @work = 1
            @start = @end = @wstart + @id
        else
            @work = 0
        end
    else
        if @id < remainder then
            workpt += 1
            @start = @wstart + workpt * @id
            @end = @start + workpt - 1
            @work = workpt
        else
            @start = @wstart + remainder + workpt * @id
            @end = @start + workpt - 1
            @work = workpt
        end
    end
  end

# *** public ***

  attr_accessor :id

  attr_accessor :done

  attr_accessor :mm1

  attr_accessor :mm2

  attr_accessor :mm3

  attr_accessor :n1

  attr_accessor :n2

  attr_accessor :n3

  attr_accessor :zoff

  attr_accessor :uoff

end
