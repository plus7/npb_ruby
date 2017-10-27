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
require 'monitor'
class ISBase #< Thread
  #include MonitorMixin
  def initialize(clss, np, serial)
    super()                    # for MonitorMixin
    @cond = nil
    @BMName = "IS"
    @CLASS = 'S'
    @MAX_ITERATIONS = 10
    @TEST_ARRAY_SIZE = 5
    @timeron = false
    @timer = NPB3_0_RUB::Timer.new()
    @S_test_index_array = [48427, 17148, 23627, 62548, 4431]
    @S_test_rank_array = [0, 18, 346, 64917, 65463]
    @W_test_index_array = [357773, 934767, 875723, 898999, 404505]
    @W_test_rank_array = [1249, 11698, 1039987, 1043896, 1048018]
    @A_test_index_array = [2112377, 662041, 5336171, 3642833, 4250760]
    @A_test_rank_array = [104, 17523, 123928, 8288932, 8388264]
    @B_test_index_array = [41869, 812306, 5102857, 18232239, 26860214]
    @B_test_rank_array = [33422937, 10244, 59149, 33135281, 99]
    @C_test_index_array = [44172927, 72999161, 74326391, 129606274, 21736814]
    @C_test_rank_array = [61147, 882988, 266290, 133997595, 133525895]
    @CLASS = clss
    @num_threads = np
    case @CLASS
    when 'S'
      @test_index_array = @S_test_index_array
      @test_rank_array = @S_test_rank_array
      @TOTAL_KEYS_LOG_2 = 16
      @MAX_KEY_LOG_2 = 11
      @NUM_BUCKETS_LOG_2 = 9
    when 'W'
      @test_index_array = @W_test_index_array
      @test_rank_array = @W_test_rank_array
      @TOTAL_KEYS_LOG_2 = 20
      @MAX_KEY_LOG_2 = 16
      @NUM_BUCKETS_LOG_2 = 10
    when 'A'
      @test_index_array = @A_test_index_array
      @test_rank_array = @A_test_rank_array
      @TOTAL_KEYS_LOG_2 = 23
      @MAX_KEY_LOG_2 = 19
      @NUM_BUCKETS_LOG_2 = 10
    when 'B'
      @test_index_array = @B_test_index_array
      @test_rank_array = @B_test_rank_array
      @TOTAL_KEYS_LOG_2 = 25
      @MAX_KEY_LOG_2 = 21
      @NUM_BUCKETS_LOG_2 = 10
    when 'C'
      @test_index_array = @C_test_index_array
      @test_rank_array = @C_test_rank_array
      @TOTAL_KEYS_LOG_2 = 27
      @MAX_KEY_LOG_2 = 23
      @NUM_BUCKETS_LOG_2 = 10
    end
    @TOTAL_KEYS = (1 << @TOTAL_KEYS_LOG_2)
    @MAX_KEY = (1 << @MAX_KEY_LOG_2)
    @NUM_BUCKETS = (1 << @NUM_BUCKETS_LOG_2)
    @NUM_KEYS = @TOTAL_KEYS
    @SIZE_OF_BUFFERS = @NUM_KEYS
    @key_array = Array.new(@SIZE_OF_BUFFERS, 0.0)
    @master_hist = Array.new(@MAX_KEY, 0.0)
    @partial_verify_vals = Array.new(@TEST_ARRAY_SIZE, 0.0)
    for i in 0..@MAX_KEY-1 do
      @master_hist[i] = 0
    end
  end

  def cond
    if @cond == nil
      @cond = self.new_cond
    end
    return @cond
  end

  def checksum(array, name, stop)
    check = 0; 
    for i in 0..array.length-1 do
      check += array[i]
    end
    puts(name + " checksum is " + check)
    if stop then
      exit(0)
    else
    end
  end

  attr_reader :num_threads
  attr_reader :MAX_KEY
  attr_reader :key_array
  attr_reader :test_index_array
  attr_reader :master_hist
  attr_reader :partial_verify_vals
  
  def start()
    @thread_obj = Thread.new(self){|r| r.run()}
  end

end
