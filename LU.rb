# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
#  Authors: R. Van der Wijngaart
#           T. Harris
#           M. Yarrow
#  Translation to Java and MultiThreaded Code
#           M. Frumkin
#           M. Schultz
#  Translation to Ruby
#           T. Nose
#           H. Tomari
# -------------------------------------------------------------------------

require "mda"
require "mathutil"
require "LUThreads/LUBase"
require "LUThreads/RHSCompute"
require "LUThreads/Scale"
require "LUThreads/Adder"
require "LUThreads/LowerJac"
require "LUThreads/UpperJac"
require "BMInOut/BMResults"
require "BMInOut/BMArgs"
require "Random"
require "Timer"
require "monitor"
class LU < LUBase
  def initialize(clss, np, ser)
    @bid = -1
    @serial = false
    super(clss, np);
    @serial = ser
  end


  def run()
    runBenchMark()
  end

  def runBenchMark()
    BMArgs.banner(@BMName, @CLASS, @serial, @num_threads)
    numTimers = @t_last + 1; 
    t_names = Array.new(numTimers, 0.0); 
    trecs = Array.new(numTimers, 0.0); 
    setTimers(t_names)
    getInputPars()
    domain()
    setcoeff() #success
    if not @serial then
      setupThreads( self )
    else
    end
    setbv()
    #puts "u" 
    #puts @u #fail
    setiv()
    erhs()
    
    if @serial then
      tm = sssor()
    else
      tm = ssor()
    end
    error()
    pintgr()
    verified = verify(@rsdnm, @errnm, @frc); 
    @results = BMResults.new(@BMName, @CLASS, @nx0, @ny0, @nz0, @itmax, tm, getMFLOPS(@itmax, tm), "floating point", verified, @serial, @num_threads, @bid)
    @results.print()
    if @timeron then
      printTimers(t_names, trecs, tm)
    else
    end
  end

  def getMFLOPS(itmax, tm)
    mflops = 0.0; 
    if tm > 0 then
        mflops = 1984.77 * @nx0 * @ny0 * @nz0 - 10923.3 * Math_dot_pow((@nx0 + @ny0 + @nz0) / 3.0, 2) + 27770.9 * (@nx0 + @ny0 + @nz0) / 3.0 - 144010.0
        mflops *= itmax / (tm * 1000000.0)
    else
    end
    return mflops
  end

  def printTimers(t_names, trecs, tm)
    fmt = DecimalFormat.new("0.000"); 
    puts("  SECTION     Time (secs)")
    for i in 0..@t_last-1 do
      trecs[i] = @timer.readTimer(i)
    end
    if tm == 0.0 then
      tm = 1.0
    else
    end
    for i in 1..@t_last-1 do
        puts("  " + t_names[i] + ":" + fmt.format(trecs[i]) + "  (" + fmt.format(trecs[i] * 100.0/ tm) + "%)")
        if i == @t_rhs then
            t = trecs[@t_rhsx] + trecs[@t_rhsy] + trecs[@t_rhsz]; 
            puts("     " + "--> total " + "sub-rhs" + ":" + fmt.format(t) + "  (" + fmt.format(t * 100.0/ tm) + "%)")
            t = trecs[i] - t
            puts("     " + "--> total " + "rest-rhs" + ":" + fmt.format(t) + "  (" + fmt.format(t * 100.0/ tm) + "%)")
        else
        end
    end
  end

  def setTimers(t_names)
    #f1 = File.new("timer.flag"); 
    @timeron = false
    if File.exist?("timer.flag") then
        @timeron = true
        t_names[@t_total] = "total"
        t_names[@t_rhsx] = "rhsx"
        t_names[@t_rhsy] = "rhsy"
        t_names[@t_rhsz] = "rhsz"
        t_names[@t_rhs] = "rhs"
        t_names[@t_jacld] = "jacld"
        t_names[@t_blts] = "blts"
        t_names[@t_jacu] = "jacu"
        t_names[@t_buts] = "buts"
        t_names[@t_add] = "add"
        t_names[@t_l2norm] = "l2norm"
    else
    end
  end

  def blts(ldmx, ldmy, ldmz, nx, ny, nz, k, omega, v, tv, ldz, ldy, ldx, d, ist, iend, jst, jend, nx0, ny0)
    
    
    
    tmat = Array.new(5 * 5, 0.0); 
    for j in jst - 1..jend - 1 do
        for i in ist - 1..iend - 1 do
            for m in 0..4 do
                tv[m + i * @isize1 + j * @jsize1] = v[m + i * @isize1 + j * @jsize1 + k * @ksize1] - omega * (ldz[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * v[0 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + ldz[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * v[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + ldz[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * v[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + ldz[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * v[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + ldz[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * v[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
            end
        end
    end
    for j in jst - 1..jend - 1 do
        for i in ist - 1..iend - 1 do
            for m in 0..4 do
                tv[m + i * @isize1 + j * @jsize1] = tv[m + i * @isize1 + j * @jsize1] - omega * (ldy[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * v[0 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + ldx[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * v[0 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + ldy[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * v[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + ldx[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * v[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + ldy[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * v[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + ldx[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * v[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + ldy[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * v[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + ldx[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * v[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + ldy[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * v[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + ldx[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * v[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1])
            end
            for m in 0..4 do
                tmat[m + 0 * 5] = d[m + 0 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 1 * 5] = d[m + 1 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 2 * 5] = d[m + 2 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 3 * 5] = d[m + 3 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 4 * 5] = d[m + 4 * @isize4 + i * @jsize4 + j * @ksize4]
            end
            tmp1 = 1.0 / tmat[0 + 0 * 5]
            tmp = tmp1 * tmat[1 + 0 * 5]
            tmat[1 + 1 * 5] = tmat[1 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[1 + 2 * 5] = tmat[1 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[1 + 3 * 5] = tmat[1 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[1 + 4 * 5] = tmat[1 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            tv[1 + i * @isize1 + j * @jsize1] = tv[1 + i * @isize1 + j * @jsize1] - tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[2 + 0 * 5]
            tmat[2 + 1 * 5] = tmat[2 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[2 + 2 * 5] = tmat[2 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[2 + 3 * 5] = tmat[2 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[2 + 4 * 5] = tmat[2 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            tv[2 + i * @isize1 + j * @jsize1] = tv[2 + i * @isize1 + j * @jsize1] - tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[3 + 0 * 5]
            tmat[3 + 1 * 5] = tmat[3 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[3 + 2 * 5] = tmat[3 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] - tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 0 * 5]
            tmat[4 + 1 * 5] = tmat[4 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[4 + 2 * 5] = tmat[4 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] - tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[1 + 1 * 5]
            tmp = tmp1 * tmat[2 + 1 * 5]
            tmat[2 + 2 * 5] = tmat[2 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[2 + 3 * 5] = tmat[2 + 3 * 5] - tmp * tmat[1 + 3 * 5]
            tmat[2 + 4 * 5] = tmat[2 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            tv[2 + i * @isize1 + j * @jsize1] = tv[2 + i * @isize1 + j * @jsize1] - tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[3 + 1 * 5]
            tmat[3 + 2 * 5] = tmat[3 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[1 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] - tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 1 * 5]
            tmat[4 + 2 * 5] = tmat[4 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[1 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] - tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[2 + 2 * 5]
            tmp = tmp1 * tmat[3 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[2 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[2 + 4 * 5]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] - tv[2 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[2 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[2 + 4 * 5]
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] - tv[2 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[3 + 3 * 5]
            tmp = tmp1 * tmat[4 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[3 + 4 * 5]
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] - tv[3 + i * @isize1 + j * @jsize1] * tmp
            v[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = tv[4 + i * @isize1 + j * @jsize1] / tmat[4 + 4 * 5]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] - tmat[3 + 4 * 5] * v[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
            v[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = tv[3 + i * @isize1 + j * @jsize1] / tmat[3 + 3 * 5]
            tv[2 + i * @isize1 + j * @jsize1] = tv[2 + i * @isize1 + j * @jsize1] - tmat[2 + 3 * 5] * v[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - tmat[2 + 4 * 5] * v[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
            v[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = tv[2 + i * @isize1 + j * @jsize1] / tmat[2 + 2 * 5]
            tv[1 + i * @isize1 + j * @jsize1] = tv[1 + i * @isize1 + j * @jsize1] - tmat[1 + 2 * 5] * v[2 + i * @isize1 + j * @jsize1 + k * @ksize1] - tmat[1 + 3 * 5] * v[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - tmat[1 + 4 * 5] * v[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
            v[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = tv[1 + i * @isize1 + j * @jsize1] / tmat[1 + 1 * 5]
            tv[0 + i * @isize1 + j * @jsize1] = tv[0 + i * @isize1 + j * @jsize1] - tmat[0 + 1 * 5] * v[1 + i * @isize1 + j * @jsize1 + k * @ksize1] - tmat[0 + 2 * 5] * v[2 + i * @isize1 + j * @jsize1 + k * @ksize1] - tmat[0 + 3 * 5] * v[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - tmat[0 + 4 * 5] * v[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
            v[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = tv[0 + i * @isize1 + j * @jsize1] / tmat[0 + 0 * 5]
        end
    end
  end

  def buts(ldmx, ldmy, ldmz, nx, ny, nz, k, omega, v, tv, d, udx, udy, udz, ist, iend, jst, jend, nx0, ny0)
    
    
    tmat = Array.new(5 * 5, 0.0); 
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    j = jend - 1
    while j >= jst - 1 do
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
        i = iend - 1
        while i >= ist - 1 do
            for m in 0..4 do
                tv[m + i * @isize1 + j * @jsize1] = omega * (udz[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * v[0 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + udz[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * v[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + udz[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * v[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + udz[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * v[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + udz[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * v[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            end
          i -= 1
        end
      j -= 1
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    j = jend - 1
    while j >= jst - 1 do
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
        i = iend - 1
        while i >= ist - 1 do
            for m in 0..4 do
                tv[m + i * @isize1 + j * @jsize1] = tv[m + i * @isize1 + j * @jsize1] + omega * (udy[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * v[0 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + udx[m + 0 * @isize4 + i * @jsize4 + j * @ksize4] * v[0 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + udy[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * v[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + udx[m + 1 * @isize4 + i * @jsize4 + j * @ksize4] * v[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + udy[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * v[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + udx[m + 2 * @isize4 + i * @jsize4 + j * @ksize4] * v[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + udy[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * v[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + udx[m + 3 * @isize4 + i * @jsize4 + j * @ksize4] * v[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + udy[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * v[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + udx[m + 4 * @isize4 + i * @jsize4 + j * @ksize4] * v[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            end
            for m in 0..4 do
                tmat[m + 0 * 5] = d[m + 0 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 1 * 5] = d[m + 1 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 2 * 5] = d[m + 2 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 3 * 5] = d[m + 3 * @isize4 + i * @jsize4 + j * @ksize4]
                tmat[m + 4 * 5] = d[m + 4 * @isize4 + i * @jsize4 + j * @ksize4]
            end
            tmp1 = 1.0 / tmat[0 + 0 * 5]
            tmp = tmp1 * tmat[1 + 0 * 5]
            tmat[1 + 1 * 5] = tmat[1 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[1 + 2 * 5] = tmat[1 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[1 + 3 * 5] = tmat[1 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[1 + 4 * 5] = tmat[1 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            tv[1 + i * @isize1 + j * @jsize1] = tv[1 + i * @isize1 + j * @jsize1] - tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[2 + 0 * 5]
            tmat[2 + 1 * 5] = tmat[2 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[2 + 2 * 5] = tmat[2 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[2 + 3 * 5] = tmat[2 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[2 + 4 * 5] = tmat[2 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            tv[2 + i * @isize1 + j * @jsize1] = tv[2 + i * @isize1 + j * @jsize1] - tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[3 + 0 * 5]
            tmat[3 + 1 * 5] = tmat[3 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[3 + 2 * 5] = tmat[3 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] - tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 0 * 5]
            tmat[4 + 1 * 5] = tmat[4 + 1 * 5] - tmp * tmat[0 + 1 * 5]
            tmat[4 + 2 * 5] = tmat[4 + 2 * 5] - tmp * tmat[0 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[0 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[0 + 4 * 5]
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] - tv[0 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[1 + 1 * 5]
            tmp = tmp1 * tmat[2 + 1 * 5]
            tmat[2 + 2 * 5] = tmat[2 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[2 + 3 * 5] = tmat[2 + 3 * 5] - tmp * tmat[1 + 3 * 5]
            tmat[2 + 4 * 5] = tmat[2 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            tv[2 + i * @isize1 + j * @jsize1] = tv[2 + i * @isize1 + j * @jsize1] - tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[3 + 1 * 5]
            tmat[3 + 2 * 5] = tmat[3 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[3 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] - tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 1 * 5]
            tmat[4 + 2 * 5] = tmat[4 + 2 * 5] - tmp * tmat[1 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[1 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[1 + 4 * 5]
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] - tv[1 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[2 + 2 * 5]
            tmp = tmp1 * tmat[3 + 2 * 5]
            tmat[3 + 3 * 5] = tmat[3 + 3 * 5] - tmp * tmat[2 + 3 * 5]
            tmat[3 + 4 * 5] = tmat[3 + 4 * 5] - tmp * tmat[2 + 4 * 5]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] - tv[2 + i * @isize1 + j * @jsize1] * tmp
            tmp = tmp1 * tmat[4 + 2 * 5]
            tmat[4 + 3 * 5] = tmat[4 + 3 * 5] - tmp * tmat[2 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[2 + 4 * 5]
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] - tv[2 + i * @isize1 + j * @jsize1] * tmp
            tmp1 = 1.0 / tmat[3 + 3 * 5]
            tmp = tmp1 * tmat[4 + 3 * 5]
            tmat[4 + 4 * 5] = tmat[4 + 4 * 5] - tmp * tmat[3 + 4 * 5]
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] - tv[3 + i * @isize1 + j * @jsize1] * tmp
            tv[4 + i * @isize1 + j * @jsize1] = tv[4 + i * @isize1 + j * @jsize1] / tmat[4 + 4 * 5]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] - tmat[3 + 4 * 5] * tv[4 + i * @isize1 + j * @jsize1]
            tv[3 + i * @isize1 + j * @jsize1] = tv[3 + i * @isize1 + j * @jsize1] / tmat[3 + 3 * 5]
            tv[2 + i * @isize1 + j * @jsize1] = tv[2 + i * @isize1 + j * @jsize1] - tmat[2 + 3 * 5] * tv[3 + i * @isize1 + j * @jsize1] - tmat[2 + 4 * 5] * tv[4 + i * @isize1 + j * @jsize1]
            tv[2 + i * @isize1 + j * @jsize1] = tv[2 + i * @isize1 + j * @jsize1] / tmat[2 + 2 * 5]
            tv[1 + i * @isize1 + j * @jsize1] = tv[1 + i * @isize1 + j * @jsize1] - tmat[1 + 2 * 5] * tv[2 + i * @isize1 + j * @jsize1] - tmat[1 + 3 * 5] * tv[3 + i * @isize1 + j * @jsize1] - tmat[1 + 4 * 5] * tv[4 + i * @isize1 + j * @jsize1]
            tv[1 + i * @isize1 + j * @jsize1] = tv[1 + i * @isize1 + j * @jsize1] / tmat[1 + 1 * 5]
            tv[0 + i * @isize1 + j * @jsize1] = tv[0 + i * @isize1 + j * @jsize1] - tmat[0 + 1 * 5] * tv[1 + i * @isize1 + j * @jsize1] - tmat[0 + 2 * 5] * tv[2 + i * @isize1 + j * @jsize1] - tmat[0 + 3 * 5] * tv[3 + i * @isize1 + j * @jsize1] - tmat[0 + 4 * 5] * tv[4 + i * @isize1 + j * @jsize1]
            tv[0 + i * @isize1 + j * @jsize1] = tv[0 + i * @isize1 + j * @jsize1] / tmat[0 + 0 * 5]
            v[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = v[0 + i * @isize1 + j * @jsize1 + k * @ksize1] - tv[0 + i * @isize1 + j * @jsize1]
            v[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = v[1 + i * @isize1 + j * @jsize1 + k * @ksize1] - tv[1 + i * @isize1 + j * @jsize1]
            v[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = v[2 + i * @isize1 + j * @jsize1 + k * @ksize1] - tv[2 + i * @isize1 + j * @jsize1]
            v[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = v[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - tv[3 + i * @isize1 + j * @jsize1]
            v[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = v[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - tv[4 + i * @isize1 + j * @jsize1]
          i -= 1
        end
      j -= 1
    end
  end

  def domain()
    @nx = @nx0
    @ny = @ny0
    @nz = @nz0
    if (@nx < 4) or (@ny < 4) or (@nz < 4) then
        puts("     " + "SUBDOMAIN SIZE IS TOO SMALL - ")
        puts("     " + "ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS")
        puts("     " + "SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL")
        puts("     " + "TO 4 THEY ARE CURRENTLY " + @nx + " " + @ny + " " + @nz)
        System.exit(0)
    else
    end
    if (@nx > @isiz1) or (@ny > @isiz2) or (@nz > @isiz3) then
        puts("     " + "SUBDOMAIN SIZE IS TOO LARGE - ")
        puts("     " + "ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS")
        puts("     " + "SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL")
        puts("     " + "TO ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY." + " THEY ARE")
        puts("     " + " CURRENTLY " + @nx + " " + @ny + " " + @nz)
        System.exit(0)
    else
    end
    @ist = 2
    @iend = @nx - 1
    @jst = 2
    @jend = @ny - 1
    @ii1 = 2
    @ii2 = @nx0 - 1
    @ji1 = 2
    @ji2 = @ny0 - 2
    @ki1 = 3
    @ki2 = @nz0 - 1
  end

  def erhs()
    
    
    
    
    
    
    
    
    
    
    
    for k in 0..@nz - 1 do
        for j in 0..@ny - 1 do
            for i in 0..@nx - 1 do
                for m in 0..4 do
                    @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] = 0.0
                end
            end
        end
    end
    for k in 0..@nz - 1 do
        zeta = k / (@nz - 1.0)
        for j in 0..@ny - 1 do
            eta = j / (@ny0 - 1.0)
            for i in 0..@nx - 1 do
                xi = i / (@nx0 - 1.0)
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @ce[m + 0 * 5] + (@ce[m + 1 * 5] + (@ce[m + 4 * 5] + (@ce[m + 7 * 5] + @ce[m + 10 * 5] * xi) * xi) * xi) * xi + (@ce[m + 2 * 5] + (@ce[m + 5 * 5] + (@ce[m + 8 * 5] + @ce[m + 11 * 5] * eta) * eta) * eta) * eta + (@ce[m + 3 * 5] + (@ce[m + 6 * 5] + (@ce[m + 9 * 5] + @ce[m + 12 * 5] * zeta) * zeta) * zeta) * zeta
                end
            end
        end
    end
    for k in 1..@nz - 2 do
        for j in @jst - 1..@jend - 1 do
            for i in 0..@nx - 1 do
                @flux[0 + i * @isize2] = @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u21 = @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                q = 0.50 * (@rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1]) / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                @flux[1 + i * @isize2] = @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * u21 + @c2 * (@rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - q)
                @flux[2 + i * @isize2] = @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * u21
                @flux[3 + i * @isize2] = @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * u21
                @flux[4 + i * @isize2] = (@c1 * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @c2 * q) * u21
            end
            for i in @ist - 1..@iend - 1 do
                for m in 0..4 do
                    @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @tx2 * (@flux[m + (i + 1) * @isize2] - @flux[m + (i - 1) * @isize2])
                end
            end
            for i in @ist - 1..@nx - 1 do
                tmp = 1.0 / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u21i = tmp * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u31i = tmp * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u41i = tmp * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u51i = tmp * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                tmp = 1.0 / @rsd[0 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                u21im1 = tmp * @rsd[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                u31im1 = tmp * @rsd[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                u41im1 = tmp * @rsd[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                u51im1 = tmp * @rsd[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                @flux[1 + i * @isize2] = (4.0 / 3.0) * @tx3 * (u21i - u21im1)
                @flux[2 + i * @isize2] = @tx3 * (u31i - u31im1)
                @flux[3 + i * @isize2] = @tx3 * (u41i - u41im1)
                @flux[4 + i * @isize2] = 0.50 * (1.0 - @c1 * @c5) * @tx3 * ((Math_dot_pow(u21i, 2) + Math_dot_pow(u31i, 2) + Math_dot_pow(u41i, 2)) - (Math_dot_pow(u21im1, 2) + Math_dot_pow(u31im1, 2) + Math_dot_pow(u41im1, 2))) + (1.0 / 6.0) * @tx3 * (Math_dot_pow(u21i, 2) - Math_dot_pow(u21im1, 2)) + @c1 * @c5 * @tx3 * (u51i - u51im1)
            end
            for i in @ist - 1..@iend - 1 do
                @frct[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @dx1 * @tx1 * (@rsd[0 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[0 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @frct[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tx3 * @c3 * @c4 * (@flux[1 + (i + 1) * @isize2] - @flux[1 + i * @isize2]) + @dx2 * @tx1 * (@rsd[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @frct[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tx3 * @c3 * @c4 * (@flux[2 + (i + 1) * @isize2] - @flux[2 + i * @isize2]) + @dx3 * @tx1 * (@rsd[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @frct[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tx3 * @c3 * @c4 * (@flux[3 + (i + 1) * @isize2] - @flux[3 + i * @isize2]) + @dx4 * @tx1 * (@rsd[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @frct[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tx3 * @c3 * @c4 * (@flux[4 + (i + 1) * @isize2] - @flux[4 + i * @isize2]) + @dx5 * @tx1 * (@rsd[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            end
            for m in 0..4 do
                @frct[m + 1 * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + 1 * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (+5.0 * @rsd[m + 1 * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + 2 * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[m + 3 * @isize1 + j * @jsize1 + k * @ksize1])
                @frct[m + 2 * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + 2 * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (-4.0 * @rsd[m + 1 * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @rsd[m + 2 * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + 3 * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[m + 4 * @isize1 + j * @jsize1 + k * @ksize1])
            end
            for i in 3..@nx - 4 do
                for m in 0..4 do
                    @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@rsd[m + (i - 2) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[m + (i + 2) * @isize1 + j * @jsize1 + k * @ksize1])
                end
            end
            for m in 0..4 do
                @frct[m + (@nx - 3) * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + (@nx - 3) * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@rsd[m + (@nx - 5) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + (@nx - 4) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @rsd[m + (@nx - 3) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + (@nx - 2) * @isize1 + j * @jsize1 + k * @ksize1])
                @frct[m + (@nx - 2) * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + (@nx - 2) * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@rsd[m + (@nx - 4) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + (@nx - 3) * @isize1 + j * @jsize1 + k * @ksize1] + 5.0 * @rsd[m + (@nx - 2) * @isize1 + j * @jsize1 + k * @ksize1])
            end
        end
    end
    for k in 1..@nz - 2 do
        for i in @ist - 1..@iend - 1 do
            for j in 0..@ny - 1 do
                @flux[0 + j * @isize2] = @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u31 = @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                q = 0.50 * (@rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1]) / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                @flux[1 + j * @isize2] = @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * u31
                @flux[2 + j * @isize2] = @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * u31 + @c2 * (@rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - q)
                @flux[3 + j * @isize2] = @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * u31
                @flux[4 + j * @isize2] = (@c1 * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @c2 * q) * u31
            end
            for j in @jst - 1..@jend - 1 do
                for m in 0..4 do
                    @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @ty2 * (@flux[m + (j + 1) * @isize2] - @flux[m + (j - 1) * @isize2])
                end
            end
            for j in @jst - 1..@ny - 1 do
                tmp = 1.0 / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u21j = tmp * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u31j = tmp * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u41j = tmp * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u51j = tmp * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                tmp = 1.0 / @rsd[0 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                u21jm1 = tmp * @rsd[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                u31jm1 = tmp * @rsd[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                u41jm1 = tmp * @rsd[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                u51jm1 = tmp * @rsd[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                @flux[1 + j * @isize2] = @ty3 * (u21j - u21jm1)
                @flux[2 + j * @isize2] = (4.0 / 3.0) * @ty3 * (u31j - u31jm1)
                @flux[3 + j * @isize2] = @ty3 * (u41j - u41jm1)
                @flux[4 + j * @isize2] = 0.50 * (1.0 - @c1 * @c5) * @ty3 * ((Math_dot_pow(u21j, 2) + Math_dot_pow(u31j, 2) + Math_dot_pow(u41j, 2)) - (Math_dot_pow(u21jm1, 2) + Math_dot_pow(u31jm1, 2) + Math_dot_pow(u41jm1, 2))) + (1.0 / 6.0) * @ty3 * (Math_dot_pow(u31j, 2) - Math_dot_pow(u31jm1, 2)) + @c1 * @c5 * @ty3 * (u51j - u51jm1)
            end
            for j in @jst - 1..@jend - 1 do
                @frct[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @dy1 * @ty1 * (@rsd[0 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[0 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
                @frct[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @ty3 * @c3 * @c4 * (@flux[1 + (j + 1) * @isize2] - @flux[1 + j * @isize2]) + @dy2 * @ty1 * (@rsd[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
                @frct[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @ty3 * @c3 * @c4 * (@flux[2 + (j + 1) * @isize2] - @flux[2 + j * @isize2]) + @dy3 * @ty1 * (@rsd[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
                @frct[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @ty3 * @c3 * @c4 * (@flux[3 + (j + 1) * @isize2] - @flux[3 + j * @isize2]) + @dy4 * @ty1 * (@rsd[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
                @frct[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @ty3 * @c3 * @c4 * (@flux[4 + (j + 1) * @isize2] - @flux[4 + j * @isize2]) + @dy5 * @ty1 * (@rsd[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            end
            for m in 0..4 do
                @frct[m + i * @isize1 + 1 * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + 1 * @jsize1 + k * @ksize1] - @dssp * (+5.0 * @rsd[m + i * @isize1 + 1 * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + i * @isize1 + 2 * @jsize1 + k * @ksize1] + @rsd[m + i * @isize1 + 3 * @jsize1 + k * @ksize1])
                @frct[m + i * @isize1 + 2 * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + 2 * @jsize1 + k * @ksize1] - @dssp * (-4.0 * @rsd[m + i * @isize1 + 1 * @jsize1 + k * @ksize1] + 6.0 * @rsd[m + i * @isize1 + 2 * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + i * @isize1 + 3 * @jsize1 + k * @ksize1] + @rsd[m + i * @isize1 + 4 * @jsize1 + k * @ksize1])
            end
            for j in 3..@ny - 4 do
                for m in 0..4 do
                    @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@rsd[m + i * @isize1 + (j - 2) * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 6.0 * @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @rsd[m + i * @isize1 + (j + 2) * @jsize1 + k * @ksize1])
                end
            end
            for m in 0..4 do
                @frct[m + i * @isize1 + (@ny - 3) * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + (@ny - 3) * @jsize1 + k * @ksize1] - @dssp * (@rsd[m + i * @isize1 + (@ny - 5) * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + i * @isize1 + (@ny - 4) * @jsize1 + k * @ksize1] + 6.0 * @rsd[m + i * @isize1 + (@ny - 3) * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + i * @isize1 + (@ny - 2) * @jsize1 + k * @ksize1])
                @frct[m + i * @isize1 + (@ny - 2) * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + (@ny - 2) * @jsize1 + k * @ksize1] - @dssp * (@rsd[m + i * @isize1 + (@ny - 4) * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + i * @isize1 + (@ny - 3) * @jsize1 + k * @ksize1] + 5.0 * @rsd[m + i * @isize1 + (@ny - 2) * @jsize1 + k * @ksize1])
            end
        end
    end
    for j in @jst - 1..@jend - 1 do
        for i in @ist - 1..@iend - 1 do
            for k in 0..@nz - 1 do
                @flux[0 + k * @isize2] = @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u41 = @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                q = 0.50 * (@rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1]) / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                @flux[1 + k * @isize2] = @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * u41
                @flux[2 + k * @isize2] = @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * u41
                @flux[3 + k * @isize2] = @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * u41 + @c2 * (@rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - q)
                @flux[4 + k * @isize2] = (@c1 * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @c2 * q) * u41
            end
            for k in 1..@nz - 2 do
                for m in 0..4 do
                    @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @tz2 * (@flux[m + (k + 1) * @isize2] - @flux[m + (k - 1) * @isize2])
                end
            end
            for k in 1..@nz - 1 do
                tmp = 1.0 / @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u21k = tmp * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u31k = tmp * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u41k = tmp * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u51k = tmp * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                tmp = 1.0 / @rsd[0 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                u21km1 = tmp * @rsd[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                u31km1 = tmp * @rsd[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                u41km1 = tmp * @rsd[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                u51km1 = tmp * @rsd[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                @flux[1 + k * @isize2] = @tz3 * (u21k - u21km1)
                @flux[2 + k * @isize2] = @tz3 * (u31k - u31km1)
                @flux[3 + k * @isize2] = (4.0 / 3.0) * @tz3 * (u41k - u41km1)
                @flux[4 + k * @isize2] = 0.50 * (1.0 - @c1 * @c5) * @tz3 * ((Math_dot_pow(u21k, 2) + Math_dot_pow(u31k, 2) + Math_dot_pow(u41k, 2)) - (Math_dot_pow(u21km1, 2) + Math_dot_pow(u31km1, 2) + Math_dot_pow(u41km1, 2))) + (1.0 / 6.0) * @tz3 * (Math_dot_pow(u41k, 2) - Math_dot_pow(u41km1, 2)) + @c1 * @c5 * @tz3 * (u51k - u51km1)
            end
            for k in 1..@nz - 2 do
                @frct[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @dz1 * @tz1 * (@rsd[0 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[0 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
                @frct[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tz3 * @c3 * @c4 * (@flux[1 + (k + 1) * @isize2] - @flux[1 + k * @isize2]) + @dz2 * @tz1 * (@rsd[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
                @frct[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tz3 * @c3 * @c4 * (@flux[2 + (k + 1) * @isize2] - @flux[2 + k * @isize2]) + @dz3 * @tz1 * (@rsd[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
                @frct[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tz3 * @c3 * @c4 * (@flux[3 + (k + 1) * @isize2] - @flux[3 + k * @isize2]) + @dz4 * @tz1 * (@rsd[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
                @frct[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tz3 * @c3 * @c4 * (@flux[4 + (k + 1) * @isize2] - @flux[4 + k * @isize2]) + @dz5 * @tz1 * (@rsd[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @rsd[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
            end
            for m in 0..4 do
                @frct[m + i * @isize1 + j * @jsize1 + 1 * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + 1 * @ksize1] - @dssp * (+5.0 * @rsd[m + i * @isize1 + j * @jsize1 + 1 * @ksize1] - 4.0 * @rsd[m + i * @isize1 + j * @jsize1 + 2 * @ksize1] + @rsd[m + i * @isize1 + j * @jsize1 + 3 * @ksize1])
                @frct[m + i * @isize1 + j * @jsize1 + 2 * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + 2 * @ksize1] - @dssp * (-4.0 * @rsd[m + i * @isize1 + j * @jsize1 + 1 * @ksize1] + 6.0 * @rsd[m + i * @isize1 + j * @jsize1 + 2 * @ksize1] - 4.0 * @rsd[m + i * @isize1 + j * @jsize1 + 3 * @ksize1] + @rsd[m + i * @isize1 + j * @jsize1 + 4 * @ksize1])
            end
            for k in 3..@nz - 4 do
                for m in 0..4 do
                    @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@rsd[m + i * @isize1 + j * @jsize1 + (k - 2) * @ksize1] - 4.0 * @rsd[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 6.0 * @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @rsd[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @rsd[m + i * @isize1 + j * @jsize1 + (k + 2) * @ksize1])
                end
            end
            for m in 0..4 do
                @frct[m + i * @isize1 + j * @jsize1 + (@nz - 3) * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + (@nz - 3) * @ksize1] - @dssp * (@rsd[m + i * @isize1 + j * @jsize1 + (@nz - 5) * @ksize1] - 4.0 * @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 4) * @ksize1] + 6.0 * @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 3) * @ksize1] - 4.0 * @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 2) * @ksize1])
                @frct[m + i * @isize1 + j * @jsize1 + (@nz - 2) * @ksize1] = @frct[m + i * @isize1 + j * @jsize1 + (@nz - 2) * @ksize1] - @dssp * (@rsd[m + i * @isize1 + j * @jsize1 + (@nz - 4) * @ksize1] - 4.0 * @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 3) * @ksize1] + 5.0 * @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 2) * @ksize1])
            end
        end
    end
  end

  def error()
    
    
    u000ijk = Array.new(5, 0.0); 
    for m in 0..4 do
        @errnm[m] = 0.0
    end
    for k in 1..@nz - 2 do
        for j in @jst - 1..@jend - 1 do
            for i in @ist - 1..@iend - 1 do
                exact(i + 1, j + 1, k + 1, u000ijk)
                for m in 0..4 do
                    tmp = (u000ijk[m] - @u[m + i * @isize1 + j * @jsize1 + k * @ksize1])
                    @errnm[m] = @errnm[m] + Math_dot_pow(tmp, 2)
                end
            end
        end
    end
    for m in 0..4 do
        @errnm[m] = Math.sqrt(@errnm[m] / ((@nx0 - 2.0) * (@ny0 - 2.0) * (@nz0 - 2.0)))
    end
  end

  def jacld(k)
    
    
    
    
    
    r43 = (4.0 / 3.0)
    c1345 = @c1 * @c3 * @c4 * @c5
    c34 = @c3 * @c4
    for j in @jst - 1..@jend - 1 do
        for i in @ist - 1..@iend - 1 do
            tmp1 = @rho_i[i + j * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @d[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * (@tx1 * @dx1 + @ty1 * @dy1 + @tz1 * @dz1)
            @d[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * 2.0 * (@tx1 * r43 + @ty1 + @tz1) * c34 * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 * r43 + @ty1 + @tz1) + @dt * 2.0 * (@tx1 * @dx2 + @ty1 * @dy2 + @tz1 * @dz2)
            @d[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * 2.0 * (@tx1 + @ty1 * r43 + @tz1) * c34 * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 + @ty1 * r43 + @tz1) + @dt * 2.0 * (@tx1 * @dx3 + @ty1 * @dy3 + @tz1 * @dz3)
            @d[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * 2.0 * (@tx1 + @ty1 + @tz1 * r43) * c34 * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 + @ty1 + @tz1 * r43) + @dt * 2.0 * (@tx1 * @dx4 + @ty1 * @dy4 + @tz1 * @dz4)
            @d[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * 2.0 * (((@tx1 * (r43 * c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (c34 - c1345)) * (Math_dot_pow(@u[1 + i * @isize1 + j * @jsize1 + k * @ksize1], 2)) + (@tx1 * (c34 - c1345) + @ty1 * (r43 * c34 - c1345) + @tz1 * (c34 - c1345)) * (Math_dot_pow(@u[2 + i * @isize1 + j * @jsize1 + k * @ksize1], 2)) + (@tx1 * (c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (r43 * c34 - c1345)) * (Math_dot_pow(@u[3 + i * @isize1 + j * @jsize1 + k * @ksize1], 2))) * tmp3 + (@tx1 + @ty1 + @tz1) * c1345 * tmp2 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * (@tx1 * (r43 * c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (c34 - c1345))
            @d[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * (@tx1 * (c34 - c1345) + @ty1 * (r43 * c34 - c1345) + @tz1 * (c34 - c1345))
            @d[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * (@tx1 * (c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (r43 * c34 - c1345))
            @d[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * (@tx1 + @ty1 + @tz1) * c1345 * tmp1 + @dt * 2.0 * (@tx1 * @dx5 + @ty1 * @dy5 + @tz1 * @dz5)
            tmp1 = @rho_i[i + j * @jsize3 + (k - 1) * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @a[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz1 * @dz1
            @a[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2
            @a[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (-(@u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) * tmp2) - @dt * @tz1 * (-c34 * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
            @a[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1) - @dt * @tz1 * c34 * tmp1 - @dt * @tz1 * @dz2
            @a[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (@u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1)
            @a[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (-(@u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) * tmp2) - @dt * @tz1 * (-c34 * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
            @a[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1) - @dt * @tz1 * (c34 * tmp1) - @dt * @tz1 * @dz3
            @a[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (@u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1)
            @a[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (-Math_dot_pow((@u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1), 2) + @c2 * @qs[i + j * @jsize3 + (k - 1) * @ksize3] * tmp1) - @dt * @tz1 * (-r43 * c34 * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
            @a[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (-@c2 * (@u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1))
            @a[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (-@c2 * (@u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1))
            @a[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (2.0 - @c2) * (@u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1) - @dt * @tz1 * (r43 * c34 * tmp1) - @dt * @tz1 * @dz4
            @a[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * @c2
            @a[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * ((@c2 * 2.0 * @qs[i + j * @jsize3 + (k - 1) * @ksize3] - @c1 * @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp2) - @dt * @tz1 * (-(c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1], 2) - (r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1], 2) - c1345 * tmp2 * @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
            @a[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (-@c2 * (@u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) * tmp2) - @dt * @tz1 * (c34 - c1345) * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
            @a[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (-@c2 * (@u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) * tmp2) - @dt * @tz1 * (c34 - c1345) * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
            @a[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (@c1 * (@u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1) - @c2 * (@qs[i + j * @jsize3 + (k - 1) * @ksize3] * tmp1 + @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp2)) - @dt * @tz1 * (r43 * c34 - c1345) * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
            @a[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz2 * (@c1 * (@u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * tmp1)) - @dt * @tz1 * c1345 * tmp1 - @dt * @tz1 * @dz5
            tmp1 = @rho_i[i + (j - 1) * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @b[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty1 * @dy1
            @b[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2
            @b[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (-(@u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (-c34 * tmp2 * @u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1])
            @b[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1) - @dt * @ty1 * (c34 * tmp1) - @dt * @ty1 * @dy2
            @b[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (@u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1)
            @b[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (-Math_dot_pow((@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1), 2) + @c2 * (@qs[i + (j - 1) * @jsize3 + k * @ksize3] * tmp1)) - @dt * @ty1 * (-r43 * c34 * tmp2 * @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1])
            @b[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (-@c2 * (@u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1))
            @b[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * ((2.0 - @c2) * (@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1)) - @dt * @ty1 * (r43 * c34 * tmp1) - @dt * @ty1 * @dy3
            @b[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (-@c2 * (@u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1))
            @b[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * @c2
            @b[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (-(@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (-c34 * tmp2 * @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1])
            @b[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (@u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1)
            @b[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1) - @dt * @ty1 * (c34 * tmp1) - @dt * @ty1 * @dy4
            @b[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * ((@c2 * 2.0 * @qs[i + (j - 1) * @jsize3 + k * @ksize3] - @c1 * @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) * (@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp2)) - @dt * @ty1 * (-(c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1], 2) - (r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1], 2) - c1345 * tmp2 * @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1])
            @b[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (-@c2 * (@u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (c34 - c1345) * tmp2 * @u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
            @b[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (@c1 * (@u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1) - @c2 * (@qs[i + (j - 1) * @jsize3 + k * @ksize3] * tmp1 + @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp2)) - @dt * @ty1 * (r43 * c34 - c1345) * tmp2 * @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
            @b[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (-@c2 * (@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (c34 - c1345) * tmp2 * @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
            @b[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty2 * (@c1 * (@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * tmp1)) - @dt * @ty1 * c1345 * tmp1 - @dt * @ty1 * @dy5
            tmp1 = @rho_i[(i - 1) + j * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @c[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx1 * @dx1
            @c[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2
            @c[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (-Math_dot_pow((@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1), 2) + @c2 * @qs[(i - 1) + j * @jsize3 + k * @ksize3] * tmp1) - @dt * @tx1 * (-r43 * c34 * tmp2 * @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @c[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * ((2.0 - @c2) * (@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)) - @dt * @tx1 * (r43 * c34 * tmp1) - @dt * @tx1 * @dx2
            @c[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (-@c2 * (@u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1))
            @c[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (-@c2 * (@u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1))
            @c[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * @c2
            @c[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (-(@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (-c34 * tmp2 * @u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @c[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (@u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)
            @c[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @dt * @tx1 * (c34 * tmp1) - @dt * @tx1 * @dx3
            @c[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (-(@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (-c34 * tmp2 * @u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @c[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (@u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)
            @c[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @dt * @tx1 * (c34 * tmp1) - @dt * @tx1 * @dx4
            @c[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * ((@c2 * 2.0 * @qs[(i - 1) + j * @jsize3 + k * @ksize3] - @c1 * @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) * @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp2) - @dt * @tx1 * (-(r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - c1345 * tmp2 * @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @c[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (@c1 * (@u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @c2 * (@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp2 + @qs[(i - 1) + j * @jsize3 + k * @ksize3] * tmp1)) - @dt * @tx1 * (r43 * c34 - c1345) * tmp2 * @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @c[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (-@c2 * (@u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (c34 - c1345) * tmp2 * @u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @c[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (-@c2 * (@u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (c34 - c1345) * tmp2 * @u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @c[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx2 * (@c1 * (@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)) - @dt * @tx1 * c1345 * tmp1 - @dt * @tx1 * @dx5
        end
    end
  end

  def jacu(k)
    
    
    
    
    
    r43 = (4.0 / 3.0)
    c1345 = @c1 * @c3 * @c4 * @c5
    c34 = @c3 * @c4
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    j = @jend - 1
    while j >= @jst - 1 do
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
        i = @iend - 1
        while i >= @ist - 1 do
            tmp1 = @rho_i[i + j * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @d[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * (@tx1 * @dx1 + @ty1 * @dy1 + @tz1 * @dz1)
            @d[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (-@tx1 * r43 - @ty1 - @tz1) * (c34 * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 * r43 + @ty1 + @tz1) + @dt * 2.0 * (@tx1 * @dx2 + @ty1 * @dy2 + @tz1 * @dz2)
            @d[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (-@tx1 - @ty1 * r43 - @tz1) * (c34 * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 + @ty1 * r43 + @tz1) + @dt * 2.0 * (@tx1 * @dx3 + @ty1 * @dy3 + @tz1 * @dz3)
            @d[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (-@tx1 - @ty1 - @tz1 * r43) * (c34 * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * c34 * tmp1 * (@tx1 + @ty1 + @tz1 * r43) + @dt * 2.0 * (@tx1 * @dx4 + @ty1 * @dy4 + @tz1 * @dz4)
            @d[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @d[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * 2.0 * (((@tx1 * (r43 * c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (c34 - c1345)) * Math_dot_pow(@u[1 + i * @isize1 + j * @jsize1 + k * @ksize1], 2) + (@tx1 * (c34 - c1345) + @ty1 * (r43 * c34 - c1345) + @tz1 * (c34 - c1345)) * Math_dot_pow(@u[2 + i * @isize1 + j * @jsize1 + k * @ksize1], 2) + (@tx1 * (c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (r43 * c34 - c1345)) * Math_dot_pow(@u[3 + i * @isize1 + j * @jsize1 + k * @ksize1], 2)) * tmp3 + (@tx1 + @ty1 + @tz1) * c1345 * tmp2 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1])
            @d[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (@tx1 * (r43 * c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (c34 - c1345)) * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (@tx1 * (c34 - c1345) + @ty1 * (r43 * c34 - c1345) + @tz1 * (c34 - c1345)) * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * 2.0 * (@tx1 * (c34 - c1345) + @ty1 * (c34 - c1345) + @tz1 * (r43 * c34 - c1345)) * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
            @d[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 1.0 + @dt * 2.0 * (@tx1 + @ty1 + @tz1) * c1345 * tmp1 + @dt * 2.0 * (@tx1 * @dx5 + @ty1 * @dy5 + @tz1 * @dz5)
            tmp1 = @rho_i[(i + 1) + j * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @a[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tx1 * @dx1
            @a[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2
            @a[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-Math_dot_pow((@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1), 2) + @c2 * @qs[(i + 1) + j * @jsize3 + k * @ksize3] * tmp1) - @dt * @tx1 * (-r43 * c34 * tmp2 * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @a[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * ((2.0 - @c2) * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)) - @dt * @tx1 * (r43 * c34 * tmp1) - @dt * @tx1 * @dx2
            @a[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-@c2 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1))
            @a[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-@c2 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1))
            @a[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * @c2
            @a[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-(@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (-c34 * tmp2 * @u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @a[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)
            @a[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @dt * @tx1 * (c34 * tmp1) - @dt * @tx1 * @dx3
            @a[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-(@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (-c34 * tmp2 * @u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @a[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)
            @a[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @dt * @tx1 * (c34 * tmp1) - @dt * @tx1 * @dx4
            @a[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @a[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * ((@c2 * 2.0 * @qs[(i + 1) + j * @jsize3 + k * @ksize3] - @c1 * @u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp2)) - @dt * @tx1 * (-(r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) - c1345 * tmp2 * @u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            @a[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@c1 * (@u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1) - @c2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp2 + @qs[(i + 1) + j * @jsize3 + k * @ksize3] * tmp1)) - @dt * @tx1 * (r43 * c34 - c1345) * tmp2 * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @a[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-@c2 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (c34 - c1345) * tmp2 * @u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @a[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (-@c2 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]) * tmp2) - @dt * @tx1 * (c34 - c1345) * tmp2 * @u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1]
            @a[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tx2 * (@c1 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * tmp1)) - @dt * @tx1 * c1345 * tmp1 - @dt * @tx1 * @dx5
            tmp1 = @rho_i[i + (j + 1) * @jsize3 + k * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @b[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @ty1 * @dy1
            @b[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2
            @b[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-(@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (-c34 * tmp2 * @u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            @b[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1) - @dt * @ty1 * (c34 * tmp1) - @dt * @ty1 * @dy2
            @b[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1)
            @b[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-Math_dot_pow((@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1), 2) + @c2 * (@qs[i + (j + 1) * @jsize3 + k * @ksize3] * tmp1)) - @dt * @ty1 * (-r43 * c34 * tmp2 * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            @b[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-@c2 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1))
            @b[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * ((2.0 - @c2) * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1)) - @dt * @ty1 * (r43 * c34 * tmp1) - @dt * @ty1 * @dy3
            @b[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-@c2 * (@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1))
            @b[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * @c2
            @b[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-(@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (-c34 * tmp2 * @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            @b[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1)
            @b[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1) - @dt * @ty1 * (c34 * tmp1) - @dt * @ty1 * @dy4
            @b[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @b[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * ((@c2 * 2.0 * @qs[i + (j + 1) * @jsize3 + k * @ksize3] - @c1 * @u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp2)) - @dt * @ty1 * (-(c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1], 2) - (r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1], 2) - c1345 * tmp2 * @u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            @b[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-@c2 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (c34 - c1345) * tmp2 * @u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]
            @b[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@c1 * (@u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1) - @c2 * (@qs[i + (j + 1) * @jsize3 + k * @ksize3] * tmp1 + @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp2)) - @dt * @ty1 * (r43 * c34 - c1345) * tmp2 * @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]
            @b[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (-@c2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]) * tmp2) - @dt * @ty1 * (c34 - c1345) * tmp2 * @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1]
            @b[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @ty2 * (@c1 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * tmp1)) - @dt * @ty1 * c1345 * tmp1 - @dt * @ty1 * @dy5
            tmp1 = @rho_i[i + j * @jsize3 + (k + 1) * @ksize3]
            tmp2 = tmp1 * tmp1
            tmp3 = tmp1 * tmp2
            @c[0 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = -@dt * @tz1 * @dz1
            @c[0 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[0 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[0 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2
            @c[0 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[1 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-(@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * tmp2) - @dt * @tz1 * (-c34 * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            @c[1 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1) - @dt * @tz1 * c34 * tmp1 - @dt * @tz1 * @dz2
            @c[1 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[1 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1)
            @c[1 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[2 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-(@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * tmp2) - @dt * @tz1 * (-c34 * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            @c[2 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[2 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1) - @dt * @tz1 * (c34 * tmp1) - @dt * @tz1 * @dz3
            @c[2 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1)
            @c[2 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
            @c[3 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-Math_dot_pow((@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1), 2) + @c2 * (@qs[i + j * @jsize3 + (k + 1) * @ksize3] * tmp1)) - @dt * @tz1 * (-r43 * c34 * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            @c[3 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-@c2 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1))
            @c[3 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-@c2 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1))
            @c[3 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (2.0 - @c2) * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1) - @dt * @tz1 * (r43 * c34 * tmp1) - @dt * @tz1 * @dz4
            @c[3 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * @c2
            @c[4 + 0 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * ((@c2 * 2.0 * @qs[i + j * @jsize3 + (k + 1) * @ksize3] - @c1 * @u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp2)) - @dt * @tz1 * (-(c34 - c1345) * tmp3 * Math_dot_pow(@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1], 2) - (c34 - c1345) * tmp3 * Math_dot_pow(@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1], 2) - (r43 * c34 - c1345) * tmp3 * Math_dot_pow(@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1], 2) - c1345 * tmp2 * @u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            @c[4 + 1 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-@c2 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * tmp2) - @dt * @tz1 * (c34 - c1345) * tmp2 * @u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]
            @c[4 + 2 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (-@c2 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]) * tmp2) - @dt * @tz1 * (c34 - c1345) * tmp2 * @u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]
            @c[4 + 3 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@c1 * (@u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1) - @c2 * (@qs[i + j * @jsize3 + (k + 1) * @ksize3] * tmp1 + @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp2)) - @dt * @tz1 * (r43 * c34 - c1345) * tmp2 * @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1]
            @c[4 + 4 * @isize4 + i * @jsize4 + j * @ksize4] = @dt * @tz2 * (@c1 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * tmp1)) - @dt * @tz1 * c1345 * tmp1 - @dt * @tz1 * @dz5
          i -= 1
        end
      j -= 1
    end
  end

  def l2norm(ldx, ldy, ldz, nx0, ny0, nz0, ist, iend, jst, jend, v, sum)
    
    for m in 0..4 do
        sum[m] = 0.0
    end
    for k in 1..nz0 - 2 do
        for j in jst - 1..jend - 1 do
            for i in ist - 1..iend - 1 do
                for m in 0..4 do
                    sum[m] = sum[m] + v[m + i * @isize1 + j * @jsize1 + k * @ksize1] * v[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
            end
        end
    end
    for m in 0..4 do
        sum[m] = Math.sqrt(sum[m] / ((nx0 - 2.0) * (ny0 - 2.0) * (nz0 - 2.0)))
    end
  end

  def pintgr()
    
    
    
    phi1 = Array.new((@isiz2 + 2) * (@isiz3 + 2), 0.0); phi2 = Array.new((@isiz2 + 2) * (@isiz3 + 2), 0.0); 
    
    isize5 = (@isiz2 + 2); 
    ibeg = @ii1
    ifin = @ii2
    jbeg = @ji1
    jfin = @ji2
    ifin1 = ifin - 1
    jfin1 = jfin - 1
    for i in 0..@isiz2 + 1 do
        for k in 0..@isiz3 + 1 do
            phi1[i + k * isize5] = 0.0
            phi2[i + k * isize5] = 0.0
        end
    end
    for j in jbeg - 1..jfin - 1 do
        for i in ibeg - 1..ifin - 1 do
            k = @ki1 - 1
            phi1[i + j * isize5] = @c2 * (@u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - 0.50 * (Math_dot_pow(@u[1 + i * @isize1 + j * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[2 + i * @isize1 + j * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[3 + i * @isize1 + j * @jsize1 + k * @ksize1], 2)) / @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1])
            k = @ki2 - 1
            phi2[i + j * isize5] = @c2 * (@u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - 0.50 * (Math_dot_pow(@u[1 + i * @isize1 + j * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[2 + i * @isize1 + j * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[3 + i * @isize1 + j * @jsize1 + k * @ksize1], 2)) / @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1])
        end
    end
    frc1 = 0.0
    for j in jbeg - 1..jfin1 - 1 do
        for i in ibeg - 1..ifin1 - 1 do
            frc1 = frc1 + (phi1[i + j * isize5] + phi1[(i + 1) + j * isize5] + phi1[i + (j + 1) * isize5] + phi1[(i + 1) + (j + 1) * isize5] + phi2[i + j * isize5] + phi2[(i + 1) + j * isize5] + phi2[i + (j + 1) * isize5] + phi2[(i + 1) + (j + 1) * isize5])
        end
    end
    frc1 = @dxi * @deta * frc1
    for i in 0..@isiz2 + 1 do
        for k in 0..@isiz3 + 1 do
            phi1[i + k * isize5] = 0.0
            phi2[i + k * isize5] = 0.0
        end
    end
    if jbeg == @ji1 then
        for k in @ki1 - 1..@ki2 - 1 do
            for i in ibeg - 1..ifin - 1 do
                phi1[i + k * isize5] = @c2 * (@u[4 + i * @isize1 + (jbeg - 1) * @jsize1 + k * @ksize1] - 0.50 * (Math_dot_pow(@u[1 + i * @isize1 + (jbeg - 1) * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[2 + i * @isize1 + (jbeg - 1) * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[3 + i * @isize1 + (jbeg - 1) * @jsize1 + k * @ksize1], 2)) / @u[0 + i * @isize1 + (jbeg - 1) * @jsize1 + k * @ksize1])
            end
        end
    else
    end
    if jfin == @ji2 then
        for k in @ki1 - 1..@ki2 - 1 do
            for i in ibeg - 1..ifin - 1 do
                phi2[i + k * isize5] = @c2 * (@u[4 + i * @isize1 + (jfin - 1) * @jsize1 + k * @ksize1] - 0.50 * (Math_dot_pow(@u[1 + i * @isize1 + (jfin - 1) * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[2 + i * @isize1 + (jfin - 1) * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[3 + i * @isize1 + (jfin - 1) * @jsize1 + k * @ksize1], 2)) / @u[0 + i * @isize1 + (jfin - 1) * @jsize1 + k * @ksize1])
            end
        end
    else
    end
    frc2 = 0.0
    for k in @ki1 - 1..@ki2 - 2 do
        for i in ibeg - 1..ifin1 - 1 do
            frc2 = frc2 + (phi1[i + k * isize5] + phi1[(i + 1) + k * isize5] + phi1[i + (k + 1) * isize5] + phi1[(i + 1) + (k + 1) * isize5] + phi2[i + k * isize5] + phi2[(i + 1) + k * isize5] + phi2[i + (k + 1) * isize5] + phi2[(i + 1) + (k + 1) * isize5])
        end
    end
    frc2 = @dxi * @dzeta * frc2
    for i in 0..@isiz2 + 1 do
        for k in 0..@isiz3 + 1 do
            phi1[i + k * isize5] = 0.0
            phi2[i + k * isize5] = 0.0
        end
    end
    if ibeg == @ii1 then
        for k in @ki1 - 1..@ki2 - 1 do
            for j in jbeg - 1..jfin - 1 do
                phi1[j + k * isize5] = @c2 * (@u[4 + (ibeg - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 0.50 * (Math_dot_pow(@u[1 + (ibeg - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[2 + (ibeg - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[3 + (ibeg - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2)) / @u[0 + (ibeg - 1) * @isize1 + j * @jsize1 + k * @ksize1])
            end
        end
    else
    end
    if ifin == @ii2 then
        for k in @ki1 - 1..@ki2 - 1 do
            for j in jbeg - 1..jfin - 1 do
                phi2[j + k * isize5] = @c2 * (@u[4 + (ifin - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 0.50 * (Math_dot_pow(@u[1 + (ifin - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[2 + (ifin - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2) + Math_dot_pow(@u[3 + (ifin - 1) * @isize1 + j * @jsize1 + k * @ksize1], 2)) / @u[0 + (ifin - 1) * @isize1 + j * @jsize1 + k * @ksize1])
            end
        end
    else
    end
    frc3 = 0.0
    for k in @ki1 - 1..@ki2 - 2 do
        for j in jbeg - 1..jfin1 - 1 do
            frc3 = frc3 + (phi1[j + k * isize5] + phi1[(j + 1) + k * isize5] + phi1[j + (k + 1) * isize5] + phi1[(j + 1) + (k + 1) * isize5] + phi2[j + k * isize5] + phi2[(j + 1) + k * isize5] + phi2[j + (k + 1) * isize5] + phi2[(j + 1) + (k + 1) * isize5])
        end
    end
    frc3 = @deta * @dzeta * frc3
    @frc = 0.25 * (frc1 + frc2 + frc3)
  end

  def getInputPars()
    #f2 = File.new("inputlu.data"); 
    if File.exist?("inputlu.data") then
         
        begin
            datafile = open("inputlu.data")
            datafile.binmode
            puts("Reading from input file inputlu.data")
            tmp = datafile.read(4)
            @ipr = tmp.unpack("i*")
            tmp = datafile.read(4)
            @inorm = tmp.unpack("i*")
            tmp = datafile.read(4)
            @itmax = tmp.unpack("i*")
            tmp = datafile.read(8)
            @dt = tmp.unpack("d*")
            tmp = datafile.read(8)
            @omega = tmp.unpack("d*")
            tmp = datafile.read(8)
            @tolrsd[0] = tmp.unpack("d*")
            tmp = datafile.read(8)
            @tolrsd[1] = tmp.unpack("d*")
            tmp = datafile.read(8)
            @tolrsd[2] = tmp.unpack("d*")
            tmp = datafile.read(8)
            @tolrsd[3] = tmp.unpack("d*")
            tmp = datafile.read(8)
            @tolrsd[4] = tmp.unpack("d*")
            tmp = datafile.read(4)
            @nx0 = tmp.unpack("i*")
            tmp = datafile.read(4)
            @ny0 = tmp.unpack("i*")
            tmp = datafile.read(4)
            @nz0 = tmp.unpack("i*")
            datafile.close
        rescue Exception => e
            STDERR.puts("exception caught!")
        ensure
        end
    else
        @ipr = @ipr_default
        @inorm = @inorm_default
        @itmax = @itmax_default
        @dt = @dt_default
        @omega = @omega_default
        @tolrsd[0] = @tolrsd1_def
        @tolrsd[1] = @tolrsd2_def
        @tolrsd[2] = @tolrsd3_def
        @tolrsd[3] = @tolrsd4_def
        @tolrsd[4] = @tolrsd5_def
        @nx0 = @isiz1
        @ny0 = @isiz2
        @nz0 = @isiz3
    end
    if (@nx0 < 4) or (@ny0 < 4) or (@nz0 < 4) then
        puts("     PROBLEM SIZE IS TOO SMALL - ")
        puts("     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5")
        System.exit(0)
    else
    end
    if (@nx0 > @isiz1) or (@ny0 > @isiz2) or (@nz0 > @isiz3) then
        puts("     PROBLEM SIZE IS TOO LARGE - ")
        puts("     NX, NY AND NZ SHOULD BE EQUAL TO")
        puts("     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY")
        System.exit(0)
    else
    end
    puts("LU: Iterations=" + @itmax.to_s + " dt=" + @dt.to_s)
  end

  def rhs()
    
    
    
    
    
    
    
    
    
    
    for k in 0..@nz - 1 do
        for j in 0..@ny - 1 do
            for i in 0..@nx - 1 do
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = -@frct[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
                tmp = 1.0 / @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                @rho_i[i + j * @jsize3 + k * @ksize3] = tmp
                @qs[i + j * @jsize3 + k * @ksize3] = 0.50 * (@u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]) * tmp
            end
        end
    end
    if @timeron then
      @timer.start(@t_rhsx)
    else
    end
    for k in 1..@nz - 2 do
        for j in @jst - 1..@jend - 1 do
            for i in 0..@nx - 1 do
                @flux[0 + i * @isize2] = @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u21 = @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize3 + k * @ksize3]
                q = @qs[i + j * @jsize3 + k * @ksize3]
                @flux[1 + i * @isize2] = @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * u21 + @c2 * (@u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - q)
                @flux[2 + i * @isize2] = @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * u21
                @flux[3 + i * @isize2] = @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * u21
                @flux[4 + i * @isize2] = (@c1 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @c2 * q) * u21
            end
            for i in @ist - 1..@iend - 1 do
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @tx2 * (@flux[m + (i + 1) * @isize2] - @flux[m + (i - 1) * @isize2])
                end
            end
            for i in @ist - 1..@nx - 1 do
                tmp = @rho_i[i + j * @jsize3 + k * @ksize3]
                u21i = tmp * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u31i = tmp * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u41i = tmp * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u51i = tmp * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                tmp = @rho_i[(i - 1) + j * @jsize3 + k * @ksize3]
                u21im1 = tmp * @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                u31im1 = tmp * @u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                u41im1 = tmp * @u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                u51im1 = tmp * @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]
                @flux[1 + i * @isize2] = (4.0 / 3.0) * @tx3 * (u21i - u21im1)
                @flux[2 + i * @isize2] = @tx3 * (u31i - u31im1)
                @flux[3 + i * @isize2] = @tx3 * (u41i - u41im1)
                @flux[4 + i * @isize2] = 0.50 * (1.0 - @c1 * @c5) * @tx3 * ((Math_dot_pow(u21i, 2) + Math_dot_pow(u31i, 2) + Math_dot_pow(u41i, 2)) - (Math_dot_pow(u21im1, 2) + Math_dot_pow(u31im1, 2) + Math_dot_pow(u41im1, 2))) + (1.0 / 6.0) * @tx3 * (Math_dot_pow(u21i, 2) - Math_dot_pow(u21im1, 2)) + @c1 * @c5 * @tx3 * (u51i - u51im1)
            end
            for i in @ist - 1..@iend - 1 do
                @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @dx1 * @tx1 * (@u[0 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tx3 * @c3 * @c4 * (@flux[1 + (i + 1) * @isize2] - @flux[1 + i * @isize2]) + @dx2 * @tx1 * (@u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tx3 * @c3 * @c4 * (@flux[2 + (i + 1) * @isize2] - @flux[2 + i * @isize2]) + @dx3 * @tx1 * (@u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tx3 * @c3 * @c4 * (@flux[3 + (i + 1) * @isize2] - @flux[3 + i * @isize2]) + @dx4 * @tx1 * (@u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tx3 * @c3 * @c4 * (@flux[4 + (i + 1) * @isize2] - @flux[4 + i * @isize2]) + @dx5 * @tx1 * (@u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            end
            for m in 0..4 do
                @rsd[m + 1 * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + 1 * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (+5.0 * @u[m + 1 * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + 2 * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + 3 * @isize1 + j * @jsize1 + k * @ksize1])
                @rsd[m + 2 * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + 2 * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (-4.0 * @u[m + 1 * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + 2 * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + 3 * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + 4 * @isize1 + j * @jsize1 + k * @ksize1])
            end
            for i in 3..@nx - 4 do
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@u[m + (i - 2) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + (i + 2) * @isize1 + j * @jsize1 + k * @ksize1])
                end
            end
            for m in 0..4 do
                @rsd[m + (@nx - 3) * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + (@nx - 3) * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@u[m + (@nx - 5) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (@nx - 4) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + (@nx - 3) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (@nx - 2) * @isize1 + j * @jsize1 + k * @ksize1])
                @rsd[m + (@nx - 2) * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + (@nx - 2) * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@u[m + (@nx - 4) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (@nx - 3) * @isize1 + j * @jsize1 + k * @ksize1] + 5.0 * @u[m + (@nx - 2) * @isize1 + j * @jsize1 + k * @ksize1])
            end
        end
    end
    if @timeron then
      @timer.stop(@t_rhsx)
    else
    end
    if @timeron then
      @timer.start(@t_rhsy)
    else
    end
    for k in 1..@nz - 2 do
        for i in @ist - 1..@iend - 1 do
            for j in 0..@ny - 1 do
                @flux[0 + j * @isize2] = @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u31 = @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize3 + k * @ksize3]
                q = @qs[i + j * @jsize3 + k * @ksize3]
                @flux[1 + j * @isize2] = @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * u31
                @flux[2 + j * @isize2] = @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * u31 + @c2 * (@u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - q)
                @flux[3 + j * @isize2] = @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * u31
                @flux[4 + j * @isize2] = (@c1 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @c2 * q) * u31
            end
            for j in @jst - 1..@jend - 1 do
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @ty2 * (@flux[m + (j + 1) * @isize2] - @flux[m + (j - 1) * @isize2])
                end
            end
            for j in @jst - 1..@ny - 1 do
                tmp = @rho_i[i + j * @jsize3 + k * @ksize3]
                u21j = tmp * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u31j = tmp * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u41j = tmp * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u51j = tmp * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                tmp = @rho_i[i + (j - 1) * @jsize3 + k * @ksize3]
                u21jm1 = tmp * @u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                u31jm1 = tmp * @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                u41jm1 = tmp * @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                u51jm1 = tmp * @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]
                @flux[1 + j * @isize2] = @ty3 * (u21j - u21jm1)
                @flux[2 + j * @isize2] = (4.0 / 3.0) * @ty3 * (u31j - u31jm1)
                @flux[3 + j * @isize2] = @ty3 * (u41j - u41jm1)
                @flux[4 + j * @isize2] = 0.50 * (1.0 - @c1 * @c5) * @ty3 * ((Math_dot_pow(u21j, 2) + Math_dot_pow(u31j, 2) + Math_dot_pow(u41j, 2)) - (Math_dot_pow(u21jm1, 2) + Math_dot_pow(u31jm1, 2) + Math_dot_pow(u41jm1, 2))) + (1.0 / 6.0) * @ty3 * (Math_dot_pow(u31j, 2) - Math_dot_pow(u31jm1, 2)) + @c1 * @c5 * @ty3 * (u51j - u51jm1)
            end
            for j in @jst - 1..@jend - 1 do
                @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @dy1 * @ty1 * (@u[0 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
                @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @ty3 * @c3 * @c4 * (@flux[1 + (j + 1) * @isize2] - @flux[1 + j * @isize2]) + @dy2 * @ty1 * (@u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
                @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @ty3 * @c3 * @c4 * (@flux[2 + (j + 1) * @isize2] - @flux[2 + j * @isize2]) + @dy3 * @ty1 * (@u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
                @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @ty3 * @c3 * @c4 * (@flux[3 + (j + 1) * @isize2] - @flux[3 + j * @isize2]) + @dy4 * @ty1 * (@u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
                @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @ty3 * @c3 * @c4 * (@flux[4 + (j + 1) * @isize2] - @flux[4 + j * @isize2]) + @dy5 * @ty1 * (@u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            end
            for m in 0..4 do
                @rsd[m + i * @isize1 + 1 * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + 1 * @jsize1 + k * @ksize1] - @dssp * (+5.0 * @u[m + i * @isize1 + 1 * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + 2 * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + 3 * @jsize1 + k * @ksize1])
                @rsd[m + i * @isize1 + 2 * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + 2 * @jsize1 + k * @ksize1] - @dssp * (-4.0 * @u[m + i * @isize1 + 1 * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + 2 * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + 3 * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + 4 * @jsize1 + k * @ksize1])
            end
            for j in 3..@ny - 4 do
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@u[m + i * @isize1 + (j - 2) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + (j + 2) * @jsize1 + k * @ksize1])
                end
            end
            for m in 0..4 do
                @rsd[m + i * @isize1 + (@ny - 3) * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + (@ny - 3) * @jsize1 + k * @ksize1] - @dssp * (@u[m + i * @isize1 + (@ny - 5) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (@ny - 4) * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + (@ny - 3) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (@ny - 2) * @jsize1 + k * @ksize1])
                @rsd[m + i * @isize1 + (@ny - 2) * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + (@ny - 2) * @jsize1 + k * @ksize1] - @dssp * (@u[m + i * @isize1 + (@ny - 4) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (@ny - 3) * @jsize1 + k * @ksize1] + 5.0 * @u[m + i * @isize1 + (@ny - 2) * @jsize1 + k * @ksize1])
            end
        end
    end
    if @timeron then
      @timer.stop(@t_rhsy)
    else
    end
    if @timeron then
      @timer.start(@t_rhsz)
    else
    end
    for j in @jst - 1..@jend - 1 do
        for i in @ist - 1..@iend - 1 do
            for k in 0..@nz - 1 do
                @flux[0 + k * @isize2] = @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u41 = @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize3 + k * @ksize3]
                q = @qs[i + j * @jsize3 + k * @ksize3]
                @flux[1 + k * @isize2] = @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * u41
                @flux[2 + k * @isize2] = @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * u41
                @flux[3 + k * @isize2] = @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * u41 + @c2 * (@u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - q)
                @flux[4 + k * @isize2] = (@c1 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @c2 * q) * u41
            end
            for k in 1..@nz - 2 do
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @tz2 * (@flux[m + (k + 1) * @isize2] - @flux[m + (k - 1) * @isize2])
                end
            end
            for k in 1..@nz - 1 do
                tmp = @rho_i[i + j * @jsize3 + k * @ksize3]
                u21k = tmp * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u31k = tmp * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u41k = tmp * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                u51k = tmp * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                tmp = @rho_i[i + j * @jsize3 + (k - 1) * @ksize3]
                u21km1 = tmp * @u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                u31km1 = tmp * @u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                u41km1 = tmp * @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                u51km1 = tmp * @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]
                @flux[1 + k * @isize2] = @tz3 * (u21k - u21km1)
                @flux[2 + k * @isize2] = @tz3 * (u31k - u31km1)
                @flux[3 + k * @isize2] = (4.0 / 3.0) * @tz3 * (u41k - u41km1)
                @flux[4 + k * @isize2] = 0.50 * (1.0 - @c1 * @c5) * @tz3 * ((Math_dot_pow(u21k, 2) + Math_dot_pow(u31k, 2) + Math_dot_pow(u41k, 2)) - (Math_dot_pow(u21km1, 2) + Math_dot_pow(u31km1, 2) + Math_dot_pow(u41km1, 2))) + (1.0 / 6.0) * @tz3 * (Math_dot_pow(u41k, 2) - Math_dot_pow(u41km1, 2)) + @c1 * @c5 * @tz3 * (u51k - u51km1)
            end
            for k in 1..@nz - 2 do
                @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @dz1 * @tz1 * (@u[0 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
                @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tz3 * @c3 * @c4 * (@flux[1 + (k + 1) * @isize2] - @flux[1 + k * @isize2]) + @dz2 * @tz1 * (@u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
                @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tz3 * @c3 * @c4 * (@flux[2 + (k + 1) * @isize2] - @flux[2 + k * @isize2]) + @dz3 * @tz1 * (@u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
                @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tz3 * @c3 * @c4 * (@flux[3 + (k + 1) * @isize2] - @flux[3 + k * @isize2]) + @dz4 * @tz1 * (@u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
                @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @tz3 * @c3 * @c4 * (@flux[4 + (k + 1) * @isize2] - @flux[4 + k * @isize2]) + @dz5 * @tz1 * (@u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            end
            for m in 0..4 do
                @rsd[m + i * @isize1 + j * @jsize1 + 1 * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + 1 * @ksize1] - @dssp * (+5.0 * @u[m + i * @isize1 + j * @jsize1 + 1 * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + 2 * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + 3 * @ksize1])
                @rsd[m + i * @isize1 + j * @jsize1 + 2 * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + 2 * @ksize1] - @dssp * (-4.0 * @u[m + i * @isize1 + j * @jsize1 + 1 * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + 2 * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + 3 * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + 4 * @ksize1])
            end
            for k in 3..@nz - 4 do
                for m in 0..4 do
                    @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @dssp * (@u[m + i * @isize1 + j * @jsize1 + (k - 2) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + (k + 2) * @ksize1])
                end
            end
            for m in 0..4 do
                @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 3) * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 3) * @ksize1] - @dssp * (@u[m + i * @isize1 + j * @jsize1 + (@nz - 5) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (@nz - 4) * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + (@nz - 3) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (@nz - 2) * @ksize1])
                @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 2) * @ksize1] = @rsd[m + i * @isize1 + j * @jsize1 + (@nz - 2) * @ksize1] - @dssp * (@u[m + i * @isize1 + j * @jsize1 + (@nz - 4) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (@nz - 3) * @ksize1] + 5.0 * @u[m + i * @isize1 + j * @jsize1 + (@nz - 2) * @ksize1])
            end
        end
    end
    if @timeron then
      @timer.stop(@t_rhsz)
    else
    end
  end

  def setcoeff()
    @dxi = 1.0 / (@nx0 - 1)
    @deta = 1.0 / (@ny0 - 1)
    @dzeta = 1.0 / (@nz0 - 1)
    @tx1 = 1.0 / (@dxi * @dxi)
    @tx2 = 1.0 / (2.0 * @dxi)
    @tx3 = 1.0 / @dxi
    @ty1 = 1.0 / (@deta * @deta)
    @ty2 = 1.0 / (2.0 * @deta)
    @ty3 = 1.0 / @deta
    @tz1 = 1.0 / (@dzeta * @dzeta)
    @tz2 = 1.0 / (2.0 * @dzeta)
    @tz3 = 1.0 / @dzeta
    @dx1 = 0.75
    @dx2 = @dx1
    @dx3 = @dx1
    @dx4 = @dx1
    @dx5 = @dx1
    @dy1 = 0.75
    @dy2 = @dy1
    @dy3 = @dy1
    @dy4 = @dy1
    @dy5 = @dy1
    @dz1 = 1.00
    @dz2 = @dz1
    @dz3 = @dz1
    @dz4 = @dz1
    @dz5 = @dz1
    @dssp = (max(@dx1, @dy1, @dz1)) / 4.0
  end

  def setbv()
    
    temp1 = Array.new(5, 0.0); temp2 = Array.new(5, 0.0); 
    for j in 0..@ny - 1 do
        for i in 0..@nx - 1 do
            exact(i + 1, j + 1, 1, temp1)
            exact(i + 1, j + 1, @nz, temp2)
            for m in 0..4 do
                @u[m + i * @isize1 + j * @jsize1 + 0 * @ksize1] = temp1[m]
                @u[m + i * @isize1 + j * @jsize1 + (@nz - 1) * @ksize1] = temp2[m]
            end
        end
    end
    for k in 0..@nz - 1 do
        for i in 0..@nx - 1 do
            exact(i + 1, 1, k + 1, temp1)
            exact(i + 1, @ny, k + 1, temp2)
            for m in 0..4 do
                @u[m + i * @isize1 + 0 * @jsize1 + k * @ksize1] = temp1[m]
                @u[m + i * @isize1 + (@ny - 1) * @jsize1 + k * @ksize1] = temp2[m]
            end
        end
    end
    for k in 0..@nz - 1 do
        for j in 0..@ny - 1 do
            exact(1, j + 1, k + 1, temp1)
            exact(@nx, j + 1, k + 1, temp2)
            for m in 0..4 do
                @u[m + 0 * @isize1 + j * @jsize1 + k * @ksize1] = temp1[m]
                @u[m + (@nx - 1) * @isize1 + j * @jsize1 + k * @ksize1] = temp2[m]
            end
        end
    end
  end

  def setiv()
    
    
    
    ue_1jk = Array.new(5, 0.0); ue_nx0jk = Array.new(5, 0.0); ue_i1k = Array.new(5, 0.0); ue_iny0k = Array.new(5, 0.0); ue_ij1 = Array.new(5, 0.0); ue_ijnz = Array.new(5, 0.0); 
    for k in 1..@nz - 2 do
        zeta = k / (@nz - 1.0)
        for j in 1..@ny - 2 do
            eta = j / (@ny0 - 1.0)
            for i in 1..@nx - 2 do
                xi = i / (@nx0 - 1.0)
                exact(1, j + 1, k + 1, ue_1jk)
                exact(@nx0, j + 1, k + 1, ue_nx0jk)
                exact(i + 1, 1, k + 1, ue_i1k)
                exact(i + 1, @ny0, k + 1, ue_iny0k)
                exact(i + 1, j + 1, 1, ue_ij1)
                exact(i + 1, j + 1, @nz, ue_ijnz)
                for m in 0..4 do
                    pxi = (1.0 - xi) * ue_1jk[m] + xi * ue_nx0jk[m]
                    peta = (1.0 - eta) * ue_i1k[m] + eta * ue_iny0k[m]
                    pzeta = (1.0 - zeta) * ue_ij1[m] + zeta * ue_ijnz[m]
                    @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] = pxi + peta + pzeta - pxi * peta - peta * pzeta - pzeta * pxi + pxi * peta * pzeta
                end
            end
        end
    end
  end

  def ssor()
    
    
    delunm = Array.new(5, 0.0); tv = Array.new(5 * (@isiz1 / 2 * 2 + 1) * @isiz2, 0.0); 
    tmp = 1.0 / (@omega * (2.0 - @omega)); 
    for j in 0..@isiz2 - 1 do
        for i in 0..@isiz1 - 1 do
            for n in 0..4 do
                for m in 0..4 do
                    @a[m + n * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
                    @b[m + n * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
                    @c[m + n * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
                    @d[m + n * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
                end
            end
        end
    end
    doRHSiteration()
    doRHSiteration()
    doRHSiteration()
    doRHSiteration()
    l2norm(@isiz1, @isiz2, @isiz3, @nx0, @ny0, @nz0, @ist, @iend, @jst, @jend, @rsd, @rsdnm)
    @timer.resetAllTimers()
    @timer.start(1)
    for istep in 1..@itmax do
        if true then #istep % 20 == 0 or istep == @itmax or istep == 1 then
            puts(" Time step " + istep.to_s)
        else
        end
      self.synchronize do
            for m in 0..@num_threads-1 do
              @scaler[m].synchronize do
                  @scaler[m].done = false
                  @scaler[m].cond.signal
              end
            end
            for m in 0..@num_threads-1 do
              while not @scaler[m].done do
                  #begin
                      cond.wait
                  #rescue InterruptedException => e
                  #ensure
                  #end
                  cond.broadcast
              end
            end
        end
        if @timeron then
          @timer.start(@t_jacld)
        else
        end
      self.synchronize do
            for m in 0..@num_threads-1 do
              @lowerjac[m].synchronize do
                  @lowerjac[m].done = false
                  @lowerjac[m].cond.signal
                  cond.broadcast
              end
            end
            while not @lowerjac[num_threads - 1].done do
                #begin
                    cond.wait
                #rescue InterruptedException => e
                #ensure
                #end
                cond.broadcast
            end
        end
        if @timeron then
          @timer.stop(@t_jacld)
        else
        end
        if @timeron then
          @timer.start(@t_jacu)
        else
        end
      self.synchronize do
            for m in 0..@num_threads-1 do
              @upperjac[m].synchronize do
                  @upperjac[m].done = false
                  @upperjac[m].cond.signal
                  cond.broadcast
              end
            end
            while not @upperjac[num_threads - 1].done do
                #begin
                    cond.wait
                #rescue InterruptedException => e
                #ensure
                #end
                cond.broadcast
            end
        end
        if @timeron then
          @timer.stop(@t_jacu)
        else
        end
        if @timeron then
          @timer.start(@t_add)
        else
        end

  self.synchronize do

            for m in 0..@num_threads-1 do
              @adder[m].synchronize do
                  @adder[m].done = false
                  @adder[m].cond.signal
              end
            end
            for m in 0..@num_threads-1 do
              while not @adder[m].done do
                  #begin
                      cond.wait
                  #rescue InterruptedException => e
                  #ensure
                  #end
                  cond.broadcast
              end
            end
        end
        if @timeron then
          @timer.stop(@t_add)
        else
        end
        if istep % @inorm == 0 then
            if @timeron then
              @timer.start(@t_l2norm)
            else
            end
            l2norm(@isiz1, @isiz2, @isiz3, @nx0, @ny0, @nz0, @ist, @iend, @jst, @jend, @rsd, delunm)
            if @timeron then
              @timer.stop(@t_l2norm)
            else
            end
        else
        end
        if @timeron then
          @timer.start(@t_rhs)
        else
        end
        doRHSiteration()
        if @timeron then
          @timer.start(@t_rhsx)
        else
        end
        doRHSiteration()
        if @timeron then
          @timer.stop(@t_rhsx)
        else
        end
        if @timeron then
          @timer.start(@t_rhsy)
        else
        end
        doRHSiteration()
        if @timeron then
          @timer.stop(@t_rhsy)
        else
        end
        if @timeron then
          @timer.start(@t_rhsz)
        else
        end
        doRHSiteration()
        if @timeron then
          @timer.stop(@t_rhsz)
        else
        end
        if @timeron then
          @timer.stop(@t_rhs)
        else
        end
        if istep % @inorm == 0 or istep == @itmax then
            if @timeron then
              @timer.start(@t_l2norm)
            else
            end
            l2norm(@isiz1, @isiz2, @isiz3, @nx0, @ny0, @nz0, @ist, @iend, @jst, @jend, @rsd, @rsdnm)
            if @timeron then
              @timer.stop(@t_l2norm)
            else
            end
        else
        end
        if (@rsdnm[0] < @tolrsd[0]) and (@rsdnm[1] < @tolrsd[1]) and (@rsdnm[2] < @tolrsd[2]) and (@rsdnm[3] < @tolrsd[3]) and (@rsdnm[4] < @tolrsd[4]) then
            @timer.stop(1)
            return @timer.readTimer(1)
        else
        end
    end
    @timer.stop(1)
    return @timer.readTimer(1)
  end

  def doRHSiteration()
    synchronize {
    for m in 0..@num_threads-1 do
      @rhscomputer[m].synchronize do
          @rhscomputer[m].done = false
          @rhscomputer[m].cond.signal
      end
    end
    for m in 0..@num_threads-1 do
      while not @rhscomputer[m].done do
          begin
              cond.wait
          rescue InterruptedException => e
          ensure
          end
          cond.broadcast
      end
    end
    }
  end

  def sssor()
    
    
    
    delunm = Array.new(5, 0.0); tv = Array.new(5 * (@isiz1 / 2 * 2 + 1) * @isiz2, 0.0); 
    tmp = 1.0 / (@omega * (2.0 - @omega))
    for j in 0..@isiz2 - 1 do
        for i in 0..@isiz1 - 1 do
            for n in 0..4 do
                for m in 0..4 do
                    @a[m + n * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
                    @b[m + n * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
                    @c[m + n * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
                    @d[m + n * @isize4 + i * @jsize4 + j * @ksize4] = 0.0
                end
            end
        end
    end
    @timer.resetAllTimers()
    rhs()
    l2norm(@isiz1, @isiz2, @isiz3, @nx0, @ny0, @nz0, @ist, @iend, @jst, @jend, @rsd, @rsdnm)
    @timer.resetAllTimers()
    @timer.start(1)
    for istep in 1..@itmax do
        if true then #istep % 20 == 0 or istep == @itmax or istep == 1 then
            puts(" Time step " + istep.to_s)
        else
        end
        if @timeron then
          @timer.start(@t_rhs)
        else
        end
        for k in 1..@nz - 2 do
            for j in @jst - 1..@jend - 1 do
                for i in @ist - 1..@iend - 1 do
                    for m in 0..4 do
                        @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @dt * @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                    end
                end
            end
        end
        if @timeron then
          @timer.stop(@t_rhs)
        else
        end
        for k in 1..@nz - 2 do
            if @timeron then
              @timer.start(@t_jacld)
            else
            end
            jacld(k)
            if @timeron then
              @timer.stop(@t_jacld)
            else
            end
            if @timeron then
              @timer.start(@t_blts)
            else
            end
            blts(@isiz1, @isiz2, @isiz3, @nx, @ny, @nz, k, @omega, @rsd, tv, @a, @b, @c, @d, @ist, @iend, @jst, @jend, @nx0, @ny0)
            if @timeron then
              @timer.stop(@t_blts)
            else
            end
        end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
        k = @nz - 2
        while k >= 1 do
            if @timeron then
              @timer.start(@t_jacu)
            else
            end
            jacu(k)
            if @timeron then
              @timer.stop(@t_jacu)
            else
            end
            if @timeron then
              @timer.start(@t_buts)
            else
            end
            buts(@isiz1, @isiz2, @isiz3, @nx, @ny, @nz, k, @omega, @rsd, tv, @d, @a, @b, @c, @ist, @iend, @jst, @jend, @nx0, @ny0)
            if @timeron then
              @timer.stop(@t_buts)
            else
            end
          k -= 1
        end
        if @timeron then
          @timer.start(@t_add)
        else
        end
        for k in 1..@nz - 2 do
            for j in @jst - 1..@jend - 1 do
                for i in @ist - 1..@iend - 1 do
                    for m in 0..4 do
                        @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] += +tmp * @rsd[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                    end
                end
            end
        end
        if @timeron then
          @timer.stop(@t_add)
        else
        end
        if istep % @inorm == 0 then
            if @timeron then
              @timer.start(@t_l2norm)
            else
            end
            l2norm(@isiz1, @isiz2, @isiz3, @nx0, @ny0, @nz0, @ist, @iend, @jst, @jend, @rsd, delunm)
            if @timeron then
              @timer.stop(@t_l2norm)
            else
            end
        else
        end
        if @timeron then
          @timer.start(@t_rhs)
        else
        end
        rhs()
        if @timeron then
          @timer.stop(@t_rhs)
        else
        end
        if istep % @inorm == 0 or istep == @itmax then
            if @timeron then
              @timer.start(@t_l2norm)
            else
            end
            l2norm(@isiz1, @isiz2, @isiz3, @nx0, @ny0, @nz0, @ist, @iend, @jst, @jend, @rsd, @rsdnm)
            if @timeron then
              @timer.stop(@t_l2norm)
            else
            end
        else
        end
        if (@rsdnm[0] < @tolrsd[0]) and (@rsdnm[1] < @tolrsd[1]) and (@rsdnm[2] < @tolrsd[2]) and (@rsdnm[3] < @tolrsd[3]) and (@rsdnm[4] < @tolrsd[4]) then
            @timer.stop(1)
            return @timer.readTimer(1)
        else
        end
    end
    @timer.stop(1)
    return @timer.readTimer(1)
  end

  def verify(xcr, xce, xci)
    xcrref = Array.new(5, 0.0); xceref = Array.new(5, 0.0); xciref = 0; xcrdif = Array.new(5, 0.0); xcedif = Array.new(5, 0.0); xcidif = 0; dtref = 0; 
    
    verified = -1; 
    clss = 'U'; 
    for m in 0..4 do
        xcrref[m] = 1.0
        xceref[m] = 1.0
    end
    xciref = 1.0
    if (@nx0 == 12) and (@ny0 == 12) and (@nz0 == 12) and (@itmax == 50) then
        clss = 'S'
        dtref = 0.5
        xcrref[0] = 1.6196343210976702E-2
        xcrref[1] = 2.1976745164821318E-3
        xcrref[2] = 1.5179927653399185E-3
        xcrref[3] = 1.5029584435994323E-3
        xcrref[4] = 3.4264073155896461E-2
        xceref[0] = 6.4223319957960924E-4
        xceref[1] = 8.4144342047347926E-5
        xceref[2] = 5.8588269616485186E-5
        xceref[3] = 5.8474222595157350E-5
        xceref[4] = 1.3103347914111294E-3
        xciref = 7.8418928865937083
    else
      if (@nx0 == 33) and (@ny0 == 33) and (@nz0 == 33) and (@itmax == 300) then
          clss = 'W'
          dtref = 1.5E-3
          xcrref[0] = 0.1236511638192E+02
          xcrref[1] = 0.1317228477799E+01
          xcrref[2] = 0.2550120713095E+01
          xcrref[3] = 0.2326187750252E+01
          xcrref[4] = 0.2826799444189E+02
          xceref[0] = 0.4867877144216
          xceref[1] = 0.5064652880982E-1
          xceref[2] = 0.9281818101960E-1
          xceref[3] = 0.8570126542733E-1
          xceref[4] = 0.1084277417792E+01
          xciref = 0.1161399311023E+02
      else
        if (@nx0 == 64) and (@ny0 == 64) and (@nz0 == 64) and (@itmax == 250) then
            clss = 'A'
            dtref = 2.0
            xcrref[0] = 7.7902107606689367E+02
            xcrref[1] = 6.3402765259692870E+01
            xcrref[2] = 1.9499249727292479E+02
            xcrref[3] = 1.7845301160418537E+02
            xcrref[4] = 1.8384760349464247E+03
            xceref[0] = 2.9964085685471943E+01
            xceref[1] = 2.8194576365003349
            xceref[2] = 7.3473412698774742
            xceref[3] = 6.7139225687777051
            xceref[4] = 7.0715315688392578E+01
            xciref = 2.6030925604886277E+01
        else
          if (@nx0 == 102) and (@ny0 == 102) and (@nz0 == 102) and (@itmax == 250) then
              clss = 'B'
              dtref = 2.0
              xcrref[0] = 3.5532672969982736E+03
              xcrref[1] = 2.6214750795310692E+02
              xcrref[2] = 8.8333721850952190E+02
              xcrref[3] = 7.7812774739425265E+02
              xcrref[4] = 7.3087969592545314E+03
              xceref[0] = 1.1401176380212709E+02
              xceref[1] = 8.1098963655421574
              xceref[2] = 2.8480597317698308E+01
              xceref[3] = 2.5905394567832939E+01
              xceref[4] = 2.6054907504857413E+02
              xciref = 4.7887162703308227E+01
          else
            if (@nx0 == 162) and (@ny0 == 162) and (@nz0 == 162) and (@itmax == 250) then
                clss = 'C'
                dtref = 2.0
                xcrref[0] = 1.03766980323537846E+04
                xcrref[1] = 8.92212458801008552E+02
                xcrref[2] = 2.56238814582660871E+03
                xcrref[3] = 2.19194343857831427E+03
                xcrref[4] = 1.78078057261061185E+04
                xceref[0] = 2.15986399716949279E+02
                xceref[1] = 1.55789559239863600E+01
                xceref[2] = 5.41318863077207766E+01
                xceref[3] = 4.82262643154045421E+01
                xceref[4] = 4.55902910043250358E+02
                xciref = 6.66404553572181300E+01
            else
            end
          end
        end
      end
    end
    for m in 0..4 do
        xcrdif[m] = Math_dot_abs((xcr[m] - xcrref[m]) / xcrref[m])
        xcedif[m] = Math_dot_abs((xce[m] - xceref[m]) / xceref[m])
    end
    xcidif = Math_dot_abs((xci - xciref) / xciref)
    epsilon = 1.0E-8; 
    if clss != 'U' then
        puts(" Verification being performed for class " + clss)
        puts(" Accuracy setting for epsilon = " + epsilon.to_s)
        if Math_dot_abs(@dt - dtref) <= epsilon then
            if verified == -1 then
              verified = 1
            else
            end
        else
            verified = 0
            clss = 'U'
            puts(" DT= " + @dt.to_s + " does not match the reference value of " + dtref.to_s)
        end
    else
        puts(" Unknown class")
        verified = -1
    end
    if clss != 'U' then
        puts(" Comparison of RMS-norms of residual")
    else
        puts(" RMS-norms of residual")
        verified = -1
    end
    verified = BMResults.printComparisonStatusArr(clss, verified, epsilon, xcr, xcrref, xcrdif)
    if clss != 'U' then
        puts(" Comparison of RMS-norms of solution error")
    else
        puts(" RMS-norms of solution error")
    end
    verified = BMResults.printComparisonStatusArr(clss, verified, epsilon, xce, xceref, xcedif)
    if clss != 'U' then
        puts(" Comparison of surface integral")
    else
        puts(" Surface integral")
    end
    verified = BMResults.printComparisonStatus(clss, verified, epsilon, xci, xciref, xcidif)
    BMResults.printVerificationStatus(clss, verified, @BMName)
    return verified
  end

  def checksum(array, size, arrayname, stop)
    sum = 0; 
    for i in 0..size-1 do
        sum += array[i]
    end
    puts("array:" + arrayname.to_s + " checksum is: " + sum.to_s)
    if stop then
      exit(0)
    else
    end
  end

  def checkSum(arr)
    csum = 0.0; 
    for k in 0..@nz - 1 do
        for j in 0..@ny - 1 do
            for i in 0..@nx - 1 do
                for m in 0..4 do
                    offset = m + i * @isize1 + j * @jsize1 + k * @ksize1; 
                    csum += (arr[offset] * arr[offset]) / (@nx * @ny * @nz * 5.0)
                end
            end
        end
    end
    return csum
  end

  def getTime()
    return @timer.readTimer(1)
  end

  def finalize()
    puts("LU: is about to be garbage collected")
    super.finalize()
  end

  def setupThreads(lu)
    @master = lu
    if @num_threads > @isiz1 - 2 then
      @num_threads = @isiz1 - 2
    else
    end
    interval1 = Array.new(@num_threads, 0.0); 
    interval2 = Array.new(@num_threads, 0.0); 
    set_interval(@num_threads, @isiz1, interval1)
    set_interval(@num_threads, @isiz1 - 2, interval2)
    partition1 = MultiDimArray(interval1.length, 2); 
    partition2 = MultiDimArray(interval2.length, 2); 
    set_partition(0, interval1, partition1)
    set_partition(1, interval2, partition2)
    @rhscomputer = Array.new(@num_threads, 0.0)
    @scaler = Array.new(@num_threads, 0.0)
    @adder = Array.new(@num_threads, 0.0)
    @lowerjac = Array.new(@num_threads, 0.0)
    @upperjac = Array.new(@num_threads, 0.0)
    for ii in 0..@num_threads-1 do
        @rhscomputer[ii] = RHSCompute.new(lu, partition1[ii][0], partition1[ii][1], partition2[ii][0], partition2[ii][1])
      @rhscomputer[ii].extend(MonitorMixin)
        @rhscomputer[ii].id = ii
        @rhscomputer[ii].start()
        @scaler[ii] = Scale.new(lu, partition2[ii][0], partition2[ii][1])
        @scaler[ii].extend(MonitorMixin)
        @scaler[ii].id = ii
        @scaler[ii].start()
        @adder[ii] = Adder.new(lu, partition2[ii][0], partition2[ii][1])
        @adder[ii].extend(MonitorMixin)
        @adder[ii].id = ii
        @adder[ii].start()
        @lowerjac[ii] = LowerJac.new(lu, partition2[ii][0], partition2[ii][1])
        @lowerjac[ii].extend(MonitorMixin)
        @lowerjac[ii].id = ii
        @lowerjac[ii].neighbor = @lowerjac
        @lowerjac[ii].start()
        @upperjac[ii] = UpperJac.new(lu, partition2[num_threads - ii - 1][0], partition2[num_threads - ii - 1][1])
        @upperjac[ii].extend(MonitorMixin)
        @upperjac[ii].id = ii
        @upperjac[ii].neighbor = @upperjac
        @upperjac[ii].start()
    end
  end

end

  def main(argv)
    lu =  nil ; 
    BMArgs.parseCmdLineArgs(argv, @BMName)
    clss = BMArgs.clss; 
    np = BMArgs.num_threads; 
    serial = BMArgs.serial; 
    #begin
    lu = LU.new(clss, np, serial)
    lu.extend(MonitorMixin)
    #rescue OutOfMemoryError => e
    #BMArgs.outOfMemoryMessage()
    #exit(0)
    #ensure
    #end
    lu.runBenchMark()
  end

main(ARGV)
