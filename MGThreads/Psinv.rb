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
class Psinv < MGBase
  def initialize(mg)
    super(mg.clss, mg.num_threads, true)
    @done = true
    @state = 0
    Init(mg)
    @r1 = Array.new(@nm + 1, 0.0)
    @r2 = Array.new(@nm + 1, 0.0)
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = mg
  end

  def Init(mg)
    @num_threads = mg.num_threads
    @r = mg.r
    @u = mg.u
    @c = mg.c
    @nm = mg.nm
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
    for i3 in @start..@end do
        for i2 in 1..@n2 - 1-1 do
            for i1 in 0..@n1-1 do
                @r1[i1] = @r[@roff + i1 + @n1 * (i2 - 1 + @n2 * i3)] + @r[@roff + i1 + @n1 * (i2 + 1 + @n2 * i3)] + @r[@roff + i1 + @n1 * (i2 + @n2 * (i3 - 1))] + @r[@roff + i1 + @n1 * (i2 + @n2 * (i3 + 1))]
                @r2[i1] = @r[@roff + i1 + @n1 * (i2 - 1 + @n2 * (i3 - 1))] + @r[@roff + i1 + @n1 * (i2 + 1 + @n2 * (i3 - 1))] + @r[@roff + i1 + @n1 * (i2 - 1 + @n2 * (i3 + 1))] + @r[@roff + i1 + @n1 * (i2 + 1 + @n2 * (i3 + 1))]
            end
            for i1 in 1..@n1 - 1-1 do
                @u[@uoff + i1 + @n1 * (i2 + @n2 * i3)] += @c[0] * @r[@roff + i1 + @n1 * (i2 + @n2 * i3)] + @c[1] * (@r[@roff + i1 - 1 + @n1 * (i2 + @n2 * i3)] + @r[@roff + i1 + 1 + @n1 * (i2 + @n2 * i3)] + @r1[i1]) + @c[2] * (@r2[i1] + @r1[i1 - 1] + @r1[i1 + 1])
            end
        end
    end
    for i3 in @start..@end do
      for i2 in 1..@n2 - 1-1 do
          @u[@uoff + @n1 * (i2 + @n2 * i3)] = @u[@uoff + @n1 - 2 + @n1 * (i2 + @n2 * i3)]
          @u[@uoff + @n1 - 1 + @n1 * (i2 + @n2 * i3)] = @u[@uoff + 1 + @n1 * (i2 + @n2 * i3)]
      end
    end
    for i3 in @start..@end do
      for i1 in 0..@n1-1 do
          @u[@uoff + i1 + @n1 * @n2 * i3] = @u[@uoff + i1 + @n1 * (@n2 - 2 + @n2 * i3)]
          @u[@uoff + i1 + @n1 * (@n2 - 1 + @n2 * i3)] = @u[@uoff + i1 + @n1 * (1 + @n2 * i3)]
      end
    end
  end

  def GetWork()
    workpt = (@wend - @wstart) / @num_threads; 
    remainder = @wend - @wstart - workpt * @num_threads; 
    if workpt == 0 then
        if @id <= @wend - @wstart then
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

  attr_accessor :n1

  attr_accessor :n2

  attr_accessor :n3

  attr_accessor :roff

  attr_accessor :uoff

end
