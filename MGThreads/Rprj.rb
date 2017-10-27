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
class Rprj < MGBase
  def initialize(mg)
    super(mg.clss, mg.num_threads, true)
    @done = true
    @state = 0
    Init(mg)
    @x1 = Array.new(@nm + 1, 0.0)
    @y1 = Array.new(@nm + 1, 0.0)
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
    @master = mg
  end

  def Init(mg)
    @num_threads = mg.num_threads
    @r = mg.r
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
    if @m1k == 3 then
      d1 = 2
    else
      d1 = 1
    end
    if @m2k == 3 then
      d2 = 2
    else
      d2 = 1
    end
    if @m3k == 3 then
      d3 = 2
    else
      d3 = 1
    end
    for j3 in @start..@end do
        i3 = 2 * j3 - d3 - 1
        for j2 in 2..@m2j - 1 do
            i2 = 2 * j2 - d2 - 1
            for j1 in 2..@m1j do
                i1 = 2 * j1 - d1 - 1
                @x1[i1 - 1] = @r[@roff + i1 - 1 + @m1k * (i2 - 1 + @m2k * i3)] + @r[@roff + i1 - 1 + @m1k * (i2 + 1 + @m2k * i3)] + @r[@roff + i1 - 1 + @m1k * (i2 + @m2k * (i3 - 1))] + @r[@roff + i1 - 1 + @m1k * (i2 + @m2k * (i3 + 1))]
                @y1[i1 - 1] = @r[@roff + i1 - 1 + @m1k * (i2 - 1 + @m2k * (i3 - 1))] + @r[@roff + i1 - 1 + @m1k * (i2 - 1 + @m2k * (i3 + 1))] + @r[@roff + i1 - 1 + @m1k * (i2 + 1 + @m2k * (i3 - 1))] + @r[@roff + i1 - 1 + @m1k * (i2 + 1 + @m2k * (i3 + 1))]
            end
            for j1 in 2..@m1j - 1 do
                i1 = 2 * j1 - d1 - 1
                y2 = @r[@roff + i1 + @m1k * (i2 - 1 + @m2k * (i3 - 1))] + @r[@roff + i1 + @m1k * (i2 - 1 + @m2k * (i3 + 1))] + @r[@roff + i1 + @m1k * (i2 + 1 + @m2k * (i3 - 1))] + @r[@roff + i1 + @m1k * (i2 + 1 + @m2k * (i3 + 1))]
                x2 = @r[@roff + i1 + @m1k * (i2 - 1 + @m2k * i3)] + @r[@roff + i1 + @m1k * (i2 + 1 + @m2k * i3)] + @r[@roff + i1 + @m1k * (i2 + @m2k * (i3 - 1))] + @r[@roff + i1 + @m1k * (i2 + @m2k * (i3 + 1))]
                @r[@zoff + j1 - 1 + @m1j * (j2 - 1 + @m2j * (j3 - 1))] = 0.5 * @r[@roff + i1 + @m1k * (i2 + @m2k * i3)] + 0.25 * (@r[@roff + i1 - 1 + @m1k * (i2 + @m2k * i3)] + @r[@roff + i1 + 1 + @m1k * (i2 + @m2k * i3)] + x2) + 0.125 * (@x1[i1 - 1] + @x1[i1 + 1] + y2) + 0.0625 * (@y1[i1 - 1] + @y1[i1 + 1])
            end
        end
    end
    for j3 in @start - 1..@end - 1 do
      for j2 in 1..@m2j - 1 do
          @r[@zoff + @m1j * (j2 + @m2j * j3)] = @r[@zoff + @m1j - 2 + @m1j * (j2 + @m2j * j3)]
          @r[@zoff + @m1j - 1 + @m1j * (j2 + @m2j * j3)] = @r[@zoff + 1 + @m1j * (j2 + @m2j * j3)]
      end
    end
    for j3 in @start - 1..@end - 1 do
      for j1 in 0..@m1j do
          @r[@zoff + j1 + @m1j * @m2j * j3] = @r[@zoff + j1 + @m1j * (@m2j - 2 + @m2j * j3)]
          @r[@zoff + j1 + @m1j * (@m2j - 1 + @m2j * j3)] = @r[@zoff + j1 + @m1j * (1 + @m2j * j3)]
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

  attr_accessor :m1k

  attr_accessor :m2k

  attr_accessor :m3k

  attr_accessor :m1j

  attr_accessor :m2j

  attr_accessor :m3j

  attr_accessor :zoff

  attr_accessor :roff

end
