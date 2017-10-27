# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
#  Authors: M. Yarrow
#  Translation to Java and to MultiThreaded Code
# 	    M. Frumkin
# 	    M. Schultz 
#  Translation to Ruby
#           T. Nose
#           H. Tomari
# -------------------------------------------------------------------------

require "mda"
require "ISThreads/ISBase"
require "ISThreads/RankThread"
require "BMInOut/BMArgs"
require "BMInOut/BMResults"
require 'Random'
require 'Timer'
require 'monitor'

class IS < ISBase
  def initialize(clss, np, ser)
    @bid = -1
    @serial = false
    @amult = 1220703125.0
    @passed_verification = 0
    super(clss, np, ser)
    @serial = ser
    @rng =  NPB3_0_RUB::Random.new()
  end

  def run()
    runBenchMark()
  end

  def runBenchMark()
    BMArgs.banner(@BMName, @CLASS, @serial, @num_threads)
    puts(" Size:  " + @TOTAL_KEYS.to_s + " Iterations:   " + @MAX_ITERATIONS.to_s)
    @timer = NPB3_0_RUB::Timer.new()
    @timer.resetTimer(0)
    initKeys(@amult)
    if @serial then
        rank(1)
    else
        setupThreads( self )
        RankThread.iteration = 1
        doSort()
        for i in 0..@MAX_KEY-1 do
            @master_hist[i] = 0
        end
        doSort()
        partial_verify(1)
    end
    #puts @passed_verification
    @passed_verification = 0
    if @CLASS != 'S' then
      puts("\n     iteration#")
    else
    end
    @timer.start(0)
    for it in 1..@MAX_ITERATIONS do
        if @CLASS != 'S' then
          puts("	  " + it.to_s)
        else
        end
        if @serial then
            rank(it)
        else
            RankThread.iteration = it
            doSort()
            for i in 0..@MAX_KEY-1 do
                @master_hist[i] = 0
            end
            doSort()
            partial_verify(it)
        end
    end
    @timer.stop(0)
    full_verify()
    verified = 0; 
    if @passed_verification == 5 * @MAX_ITERATIONS + 1 then
      verified = 1
    else
    end
    BMResults.printVerificationStatus(@CLASS, verified, @BMName)
    tm = @timer.readTimer(0); 
    res = BMResults.new(@BMName, @CLASS, @TOTAL_KEYS, 0, 0, @MAX_ITERATIONS, tm, getMOPS(tm, @MAX_ITERATIONS, @TOTAL_KEYS), "keys ranked", verified, @serial, @num_threads, @bid); 
    res.print()
  end

  def getMOPS(total_time, niter, num_keys)
    mops = 0.0; 
    if total_time > 0 then
        mops = niter + num_keys
        mops *= niter / (total_time * 1000000.0)
    else
    end
    return mops
  end

  def rank(iteration)
    @key_array[iteration] = iteration
    @key_array[iteration + @MAX_ITERATIONS] = @MAX_KEY - iteration
    for i in 0..@TEST_ARRAY_SIZE-1 do
        @partial_verify_vals[i] = @key_array[@test_index_array[i]]
    end
    for i in 0..@MAX_KEY-1 do
      @master_hist[i] = 0
    end
    for i in 0..@NUM_KEYS-1 do
      @master_hist[@key_array[i]] += 1 
    end
    for i in 0..@MAX_KEY - 1-1 do
        @master_hist[i + 1] += @master_hist[i]
    end
    partial_verify(iteration)
  end

  def partial_verify(iteration)
    for i in 0..@TEST_ARRAY_SIZE-1 do
        k = @partial_verify_vals[i]; 
        offset = iteration; 
        if 0 <= k and k <= @NUM_KEYS - 1 then
            case @CLASS
            when 'S'
              if i <= 2 then
                offset = iteration
              else
                offset = -iteration
              end
            when 'W'
              if i < 2 then
                offset = iteration - 2
              else
                offset = -iteration
              end
            when 'A'
              if i <= 2 then
                offset = iteration - 1
              else
                offset = -iteration + 1
              end
            when 'B'
              if i == 1 or i == 2 or i == 4 then
                offset = iteration
              else
                offset = -iteration
              end
            when 'C'
              if i <= 2 then
                offset = iteration
              else
                offset = -iteration
              end
            end
            if @master_hist[k - 1] != @test_rank_array[i] + offset then
                puts("Failed partial verification: " + "iteration" + iteration.to_s + ", test key " + i.to_s)
            else
              #puts @passed_verification
              @passed_verification += 1
            end
        else
        end
    end
  end

  def full_verify()
    key = 0; idx = 0; 
    for i in 0..@NUM_KEYS-1 do
        while idx == @master_hist[key] do
            key += 1
            if key >= @MAX_KEY or idx >= @NUM_KEYS then
              break
            else
            end
        end
        @key_array[idx] = key
        idx += 1
    end
    count = 0; 
    for i in 1..@NUM_KEYS-1 do
      if @key_array[i - 1] > @key_array[i] then
        count += 1
      else
      end
    end
    if count != 0 then
        puts("Full_verify: number of keys out of sort: " + count)
    else
      #puts @passed_verification
      @passed_verification += 1
    end
    return @passed_verification
  end

  def initKeys(a)
    
    k = (@MAX_KEY / 4).to_i; 
    for i in 0..@NUM_KEYS-1 do
        x = @rng.randlc(a)
        x += @rng.randlc(a)
        x += @rng.randlc(a)
        x += @rng.randlc(a)
        @key_array[i] = (x * k).to_i
    end
  end

  def doSort()
    synchronize {
    for m in 0..@num_threads-1 do
      @rankthreads[m].synchronize do# synchronized (rankthreads[m]) {
        @rankthreads[m].done = false
        @rankthreads[m].cond.signal
      end
    end
    for m in 0..@num_threads-1 do
      while not @rankthreads[m].done do
        self.cond.wait
        self.cond.broadcast
      end
    end
    }
  end

  def getTime()
    return @timer.readTimer(0)
  end

  def finalize()
    puts("IS: is about to be garbage collected")
    super.finalize()
  end

  def setupThreads(is)
    start = 0; _end = 0; remainder = @TOTAL_KEYS % @num_threads; offset = 0; 
    rstart = 0; rend = 0; rremainder = @MAX_KEY % @num_threads; roffset = 0; 
    @rankthreads = Array.new(@num_threads, 0.0)
    for i in 0..@num_threads-1 do
        start = i * (@TOTAL_KEYS / @num_threads) + offset
        _end = i * (@TOTAL_KEYS / @num_threads) + (@TOTAL_KEYS / @num_threads) - 1 + offset
        if remainder > 0 then
            remainder -= 1
            offset += 1
            _end += 1
        else
        end
        rstart = i * (@MAX_KEY / @num_threads) + roffset
        rend = i * (@MAX_KEY / @num_threads) + (@MAX_KEY / @num_threads) - 1 + roffset
        if rremainder > 0 then
            rremainder -= 1
            roffset += 1
            rend += 1
        else
        end
        @rankthreads[i] = RankThread.new(@CLASS, is, i, start, _end, rstart, rend)
        @rankthreads[i].extend(MonitorMixin)
        @rankthreads[i].start()
    end
    for i in 0..@num_threads-1 do
        @rankthreads[i].rankthreads = @rankthreads
    end
  end
end


  def main(argv)
    #is =  nil ; 
    BMArgs.parseCmdLineArgs(argv, @BMName)
    clss = BMArgs.clss; 
    np = BMArgs.num_threads; 
    serial = BMArgs.serial; 
    #begin
    is = IS.new(clss, np, serial)
    is.extend(MonitorMixin)
    #rescue OutOfMemoryError => e
    #    BMArgs.outOfMemoryMessage()
    #    System.exit(0)
    #ensure
    #end
    is.runBenchMark()
  end

main(ARGV)
