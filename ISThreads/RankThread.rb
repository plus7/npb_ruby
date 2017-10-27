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
class RankThread < ISBase
  attr_accessor :rankthreads
  attr_accessor :done
  attr_accessor :local_hist
  #attr_accessor :iteration
  def self.iteration
    return @@iteration
  end
  def self.iteration=(n)
    @@iteration = n
  end

  def initialize(clss, is, id, s1, e1, s2, e2)
    super(clss, is.num_threads, false)
    @done = true
    @@iteration = 0
    Init(is)
    @master = is
    @id = id
    @start = s1
    @end = e1
    @rstart = s2
    @rend = e2
    @local_hist = Array.new(@MAX_KEY, 0.0)
    @state = 0
    #setPriority(Thread.MAX_PRIORITY)
    #setDaemon(true)
  end

  def Init(is)
    # @num_threads = is.num_threads
    @MAX_KEY = is.MAX_KEY
    @key_array = is.key_array
    @test_index_array = is.test_index_array
    @master_hist = is.master_hist
    @partial_verify_vals = is.partial_verify_vals
  end

  def run()
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    while true do
      self.synchronize do #synchronized: this
        while @done == true do
          cond.wait
          @master.synchronize do #synchronized: master
            @master.cond.signal# @master.notify()
          end
        end
        
        case @state
        when 0
          step1()
          @state = 1
        when 1
          step2()
          @state = 0
        end

        @master.synchronize do #synchronized: master
          @done = true
          @master.cond.signal # @master.notify()
        end
        end
    end
  end

  def step1()
    self.synchronize {

    @key_array[@@iteration] = @@iteration

    @key_array[@@iteration + @MAX_ITERATIONS] = @MAX_KEY - @@iteration

    for i in 0..@TEST_ARRAY_SIZE-1 do
        @partial_verify_vals[i] = @key_array[@test_index_array[i]]
    end

    for i in 0..@MAX_KEY-1 do
      @local_hist[i] = 0
    end

    for i in @start..@end do
      @local_hist[@key_array[i]] += 1
    end

    for i in 0..@MAX_KEY - 1-1 do
      @local_hist[i + 1] += @local_hist[i]
    end

    }
  end

  def step2()
    for i in @rstart..@rend do
        for j in 0..@num_threads-1 do
            @master_hist[i] += @rankthreads[j].local_hist[i]
        end
    end
  end

end
