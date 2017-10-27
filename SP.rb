# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
#  Authors: R. Van der Wijngaart
# 	    T. Harris
# 	    M. Yarrow
#  Modified for PBN (Programming Baseline for NPB):
# 	    H. Jin
#  Translation to Java and MultiThreaded Code
# 	    M. Frumkin
# 	    M. Schultz
#  Translation to Ruby
#           T. Nose
#           H. Tomari
# -------------------------------------------------------------------------

require "mda"
require "mathutil"
require "SPThreads/SPBase"
require "SPThreads/RHSCompute"
require "SPThreads/XSolver"
require "SPThreads/YSolver"
require "SPThreads/ZSolver"
require "SPThreads/RHSAdder"
require "SPThreads/TXInverse"
require "BMInOut/BMResults"
require "BMInOut/BMArgs"
require "Timer"
require "Random"
require "Runnable"
require "monitor"

class SP < SPBase
  def initialize(clss, np, ser)
    @master =  nil 
    @bid = -1
    @serial = false
    super(clss, np)
    @serial = ser
  end

  def run()
    runBenchMark()
  end

  def runBenchMark()
    BMArgs.banner(@@BMName, @CLASS, @serial, @num_threads)
    numTimers = @@t_last + 1; 
    t_names = Array.new(numTimers, 0.0); 
    trecs = Array.new(numTimers, 0.0); 
    setTimers(t_names)
    niter = getInputPars(); 
    set_constants(0)
    _initialize()
    exact_rhs()
    if not @serial then
      setupThreads( self )
    else
    end
    if @serial then
      adi_serial()
    else
      adi()
    end
    _initialize()
    @timer.resetAllTimers()
    @timer.start(@@t_total)
    for step in 1..niter do
        if step % 20 == 0 or step == 1 or step == niter then
            puts("Time step " + step.to_s)
        else
        end
        if @serial then
          adi_serial()
        else
          adi()
        end
    end
    @timer.stop(1)
    verified = verify(niter); 
    time = @timer.readTimer(@@t_total); 
    @results = BMResults.new(@@BMName, @CLASS, @grid_points[0], @grid_points[1], @grid_points[2], niter, time, getMFLOPS(time, niter), "floating point", verified, @serial, @num_threads, @bid)
    @results.print()
    if @timeron then
      printTimers(t_names, trecs, time)
    else
    end
  end

  def getMFLOPS(total_time, niter)
    mflops = 0.0; 
    if total_time > 0 then
        n3 = @grid_points[0] * @grid_points[1] * @grid_points[2]; 
        t = (@grid_points[0] + @grid_points[1] + @grid_points[2]) / 3.0; 
        mflops = 881.174 * n3 - 4683.91 * t * t + 11484.5 * t - 19272.4
        mflops *= niter / (total_time * 1000000.0)
    else
    end
    return mflops
  end

  def adi_serial()
    if @timeron then
      @timer.start(@@t_rhs)
    else
    end
    compute_rhs()
    if @timeron then
      @timer.stop(@@t_rhs)
    else
    end
    if @timeron then
      @timer.start(@@t_txinvr)
    else
    end
    txinvr()
    if @timeron then
      @timer.stop(@@t_txinvr)
    else
    end
    x_solve()
    y_solve()
    z_solve()
    if @timeron then
      @timer.start(@@t_add)
    else
    end
    add()
    if @timeron then
      @timer.stop(@@t_add)
    else
    end
  end

  def adi()
    if @timeron then
      @timer.start(@@t_rhs)
    else
    end
    doRHS()
    doRHS()
    if @timeron then
      @timer.start(@@t_rhsx)
    else
    end
    doRHS()
    if @timeron then
      @timer.stop(@@t_rhsx)
    else
    end
    if @timeron then
      @timer.start(@@t_rhsy)
    else
    end
    doRHS()
    if @timeron then
      @timer.stop(@@t_rhsy)
    else
    end
    if @timeron then
      @timer.start(@@t_rhsz)
    else
    end
    doRHS()
    if @timeron then
      @timer.stop(@@t_rhsz)
    else
    end
    doRHS()
    if @timeron then
      @timer.stop(@@t_rhs)
    else
    end
    if @timeron then
      @timer.start(@@t_txinvr)
    else
    end
     self.synchronize do
        for m in 0..@num_threads-1 do
          @txinverse[m].synchronize do
              @txinverse[m].done = false
              @txinverse[m].cond.signal
          end
        end
        for m in 0..@num_threads-1 do
          while not @txinverse[m].done do
              begin
                  cond.wait
              rescue InterruptedException => e
              ensure
              end
              cond.broadcast
          end
        end
    end
    if @timeron then
      @timer.stop(@@t_txinvr)
    else
    end
    if @timeron then
      @timer.start(@@t_xsolve)
    else
    end
    doXsolve()
    if @timeron then
      @timer.stop(@@t_xsolve)
    else
    end
    if @timeron then
      @timer.start(@@t_ninvr)
    else
    end
    doXsolve()
    if @timeron then
      @timer.stop(@@t_ninvr)
    else
    end
    if @timeron then
      @timer.start(@@t_ysolve)
    else
    end
    doYsolve()
    if @timeron then
      @timer.stop(@@t_ysolve)
    else
    end
    if @timeron then
      @timer.start(@@t_pinvr)
    else
    end
    doYsolve()
    if @timeron then
      @timer.stop(@@t_pinvr)
    else
    end
    if @timeron then
      @timer.start(@@t_zsolve)
    else
    end
    doZsolve()
    if @timeron then
      @timer.stop(@@t_zsolve)
    else
    end
    if @timeron then
      @timer.start(@@t_tzetar)
    else
    end
    doZsolve()
    if @timeron then
      @timer.stop(@@t_tzetar)
    else
    end
    if @timeron then
      @timer.start(@@t_add)
    else
    end

     self .synchronize do
        for m in 0..@num_threads-1 do
          @rhsadder[m].synchronize do
              @rhsadder[m].done = false
              @rhsadder[m].cond.signal
          end
        end
        for m in 0..@num_threads-1 do
          while not @rhsadder[m].done do
              begin
                  cond.wait
              rescue InterruptedException => e
              ensure
              end
              cond.broadcast
          end
        end
    end
    if @timeron then
      @timer.stop(@@t_add)
    else
    end
  end

  def doRHS()
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

  def doXsolve()
    synchronize {
    
    for m in 0..@num_threads-1 do
      @xsolver[m].synchronize do
          @xsolver[m].done = false
          @xsolver[m].cond.signal
      end
    end
    for m in 0..@num_threads-1 do
      while not @xsolver[m].done do
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

  def doYsolve()
    synchronize {
    
    for m in 0..@num_threads-1 do
      @ysolver[m].synchronize do
          @ysolver[m].done = false
          @ysolver[m].cond.signal
      end
    end
    for m in 0..@num_threads-1 do
      while not @ysolver[m].done do
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

  def doZsolve()
    synchronize {
    
    for m in 0..@num_threads-1 do
      @zsolver[m].synchronize do
          @zsolver[m].done = false
          @zsolver[m].cond.signal
      end
    end
    for m in 0..@num_threads-1 do
      while not @zsolver[m].done do
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

  def getInputPars()
    niter = 0; 
    #f2 = File.new("inputsp.data"); 
    if File.exist?("inputsp.data") then
        begin
            fis = FileInputStream.new(f2); 
            datafile = DataInputStream.new(fis); 
            puts("Reading from input file inputsp.data")
            niter = datafile.readInt()
            @@dt = datafile.readDouble()
            @grid_points[0] = datafile.readInt()
            @grid_points[1] = datafile.readInt()
            @grid_points[2] = datafile.readInt()
            fis.close()
        rescue Exception => e
            STDERR.puts("exception caught!")
        ensure
        end
    else
        puts("No input file inputsp.data," + "Using compiled defaults")
        niter = @niter_default
        @@dt = @dt_default
        @grid_points[0] = @problem_size
        @grid_points[1] = @problem_size
        @grid_points[2] = @problem_size
    end
    puts("Size: " + @grid_points[0].to_s + " X " + @grid_points[1].to_s + " X " + @grid_points[2].to_s)
    if (@grid_points[0] > @IMAX) or (@grid_points[1] > @JMAX) or (@grid_points[2] > @KMAX) then
        puts("Problem size too big for array")
        exit(0)
    else
    end
    puts("Iterations: " + niter.to_s + " dt: " + @@dt.to_s)
    @nx2 = @grid_points[0] - 2
    @ny2 = @grid_points[1] - 2
    @nz2 = @grid_points[2] - 2
    return niter
  end

  def setTimers(t_names)
    #f1 = File.new("timer.flag"); 
    @timeron = false
    if File.exist?("timer.flag") then
        @timeron = true
        t_names[@@t_total] = "total"
        t_names[@@t_rhsx] = "rhsx"
        t_names[@@t_rhsy] = "rhsy"
        t_names[@@t_rhsz] = "rhsz"
        t_names[@@t_rhs] = "rhs"
        t_names[@@t_xsolve] = "xsolve"
        t_names[@@t_ysolve] = "ysolve"
        t_names[@@t_zsolve] = "zsolve"
        t_names[@@t_rdis1] = "redist1"
        t_names[@@t_rdis2] = "redist2"
        t_names[@@t_tzetar] = "tzetar"
        t_names[@@t_ninvr] = "ninvr"
        t_names[@@t_pinvr] = "pinvr"
        t_names[@@t_txinvr] = "txinvr"
        t_names[@@t_add] = "add"
    else
    end
  end

  def printTimers(t_names, trecs, tmax)
    #fmt = DecimalFormat.new("0.000"); 
    
    puts("  SECTION   Time (secs)")
    for i in 1..@@t_last do
        trecs[i] = @timer.readTimer(i)
    end
    if tmax == 0.0 then
      tmax = 1.0
    else
    end
    for i in 1..@@t_last-1 do
        puts(t_names[i] + ":" + (trecs[i]).to_s + "  (" + (trecs[i] * 100 / tmax).to_s + "%)")
        if i == @@t_rhs then
            t = trecs[@@t_rhsx] + trecs[@@t_rhsy] + trecs[@@t_rhsz]
            puts("    --> total " + "sub-rhs" + ":" + (t).to_s + "  (" + (t * 100.0 / tmax).to_s + "%)")
            t = trecs[@@t_rhs] - t
            puts("    --> total " + "rest-rhs" + ":" + (t).to_s + "  (" + (t * 100.0 / tmax).to_s + "%)")
        else
          if i == @@t_zsolve then
              t = trecs[@@t_zsolve] - trecs[@@t_rdis1] - trecs[@@t_rdis2]
              puts("    --> total " + "sub-zsol" + ":" + (t).to_s + "  (" + (t * 100.0 / tmax).to_s + "%)")
          else
            if i == @@t_rdis2 then
                t = trecs[@@t_rdis1] + trecs[@@t_rdis2]
                puts("    --> total " + "redist" + ":" + (t).to_s + "  (" + (t * 100.0 / tmax).to_s + "%)")
            else
            end
          end
        end
    end
  end

  def add()
    
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                for m in 0..4 do
                    @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] += @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
            end
        end
    end
  end

  def error_norm(rms)
    
    u_exact = Array.new(5, 0.0); 
    for m in 0..4 do
        rms[m] = 0.0
    end
    for k in 0..@grid_points[2] - 1 do
        zeta = k * @@dnzm1
        for j in 0..@grid_points[1] - 1 do
            eta = j * @@dnym1
            for i in 0..@grid_points[0] - 1 do
                xi = i * @@dnxm1
                exact_solution(xi, eta, zeta, u_exact, 0)
                for m in 0..4 do
                    add = @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - u_exact[m]
                    rms[m] = rms[m] + add * add
                end
            end
        end
    end
    for m in 0..4 do
        for d in 0..2 do
            rms[m] = rms[m] / (@grid_points[d] - 2)
        end
        rms[m] = Math.sqrt(rms[m])
    end
  end

  def rhs_norm(rms)
    
    
    for m in 0..4 do
        rms[m] = 0.0
    end
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                for m in 0..4 do
                    add = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                    rms[m] = rms[m] + add * add
                end
            end
        end
    end
    for m in 0..4 do
        for d in 0..2 do
            rms[m] = rms[m] / (@grid_points[d] - 2)
        end
        rms[m] = Math.sqrt(rms[m])
    end
  end

  def exact_rhs()
    dtemp = Array.new(5, 0.0); 
    
    for k in 0..@grid_points[2] - 1 do
        for j in 0..@grid_points[1] - 1 do
            for i in 0..@grid_points[0] - 1 do
                for m in 0..4 do
                    @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = 0.0
                end
            end
        end
    end
    for k in 1..@grid_points[2] - 2 do
        zeta = k * @@dnzm1
        for j in 1..@grid_points[1] - 2 do
            eta = j * @@dnym1
            for i in 0..@grid_points[0] - 1 do
                xi = i * @@dnxm1
                exact_solution(xi, eta, zeta, dtemp, 0)
                for m in 0..4 do
                    @ue[i + m * @jsize3] = dtemp[m]
                end
                dtpp = 1.0 / dtemp[0]
                for m in 1..4 do
                    @buf[i + m * @jsize3] = dtpp * dtemp[m]
                end
                @cuf[i] = @buf[i + 1 * @jsize3] * @buf[i + 1 * @jsize3]
                @buf[i + 0 * @jsize3] = @cuf[i] + @buf[i + 2 * @jsize3] * @buf[i + 2 * @jsize3] + @buf[i + 3 * @jsize3] * @buf[i + 3 * @jsize3]
                @q[i] = 0.5 * (@buf[i + 1 * @jsize3] * @ue[i + 1 * @jsize3] + @buf[i + 2 * @jsize3] * @ue[i + 2 * @jsize3] + @buf[i + 3 * @jsize3] * @ue[i + 3 * @jsize3])
            end
            for i in 1..@grid_points[0] - 2 do
                im1 = i - 1
                ip1 = i + 1
                @forcing[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[0 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tx2 * (@ue[ip1 + 1 * @jsize3] - @ue[im1 + 1 * @jsize3]) + @@dx1tx1 * (@ue[ip1 + 0 * @jsize3] - 2.0 * @ue[i + 0 * @jsize3] + @ue[im1 + 0 * @jsize3])
                @forcing[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[1 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tx2 * ((@ue[ip1 + 1 * @jsize3] * @buf[ip1 + 1 * @jsize3] + @@c2 * (@ue[ip1 + 4 * @jsize3] - @q[ip1])) - (@ue[im1 + 1 * @jsize3] * @buf[im1 + 1 * @jsize3] + @@c2 * (@ue[im1 + 4 * @jsize3] - @q[im1]))) + @@xxcon1 * (@buf[ip1 + 1 * @jsize3] - 2.0 * @buf[i + 1 * @jsize3] + @buf[im1 + 1 * @jsize3]) + @@dx2tx1 * (@ue[ip1 + 1 * @jsize3] - 2.0 * @ue[i + 1 * @jsize3] + @ue[im1 + 1 * @jsize3])
                @forcing[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[2 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tx2 * (@ue[ip1 + 2 * @jsize3] * @buf[ip1 + 1 * @jsize3] - @ue[im1 + 2 * @jsize3] * @buf[im1 + 1 * @jsize3]) + @@xxcon2 * (@buf[ip1 + 2 * @jsize3] - 2.0 * @buf[i + 2 * @jsize3] + @buf[im1 + 2 * @jsize3]) + @@dx3tx1 * (@ue[ip1 + 2 * @jsize3] - 2.0 * @ue[i + 2 * @jsize3] + @ue[im1 + 2 * @jsize3])
                @forcing[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tx2 * (@ue[ip1 + 3 * @jsize3] * @buf[ip1 + 1 * @jsize3] - @ue[im1 + 3 * @jsize3] * @buf[im1 + 1 * @jsize3]) + @@xxcon2 * (@buf[ip1 + 3 * @jsize3] - 2.0 * @buf[i + 3 * @jsize3] + @buf[im1 + 3 * @jsize3]) + @@dx4tx1 * (@ue[ip1 + 3 * @jsize3] - 2.0 * @ue[i + 3 * @jsize3] + @ue[im1 + 3 * @jsize3])
                @forcing[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tx2 * (@buf[ip1 + 1 * @jsize3] * (@@c1 * @ue[ip1 + 4 * @jsize3] - @@c2 * @q[ip1]) - @buf[im1 + 1 * @jsize3] * (@@c1 * @ue[im1 + 4 * @jsize3] - @@c2 * @q[im1])) + 0.5 * @@xxcon3 * (@buf[ip1 + 0 * @jsize3] - 2.0 * @buf[i + 0 * @jsize3] + @buf[im1 + 0 * @jsize3]) + @@xxcon4 * (@cuf[ip1] - 2.0 * @cuf[i] + @cuf[im1]) + @@xxcon5 * (@buf[ip1 + 4 * @jsize3] - 2.0 * @buf[i + 4 * @jsize3] + @buf[im1 + 4 * @jsize3]) + @@dx5tx1 * (@ue[ip1 + 4 * @jsize3] - 2.0 * @ue[i + 4 * @jsize3] + @ue[im1 + 4 * @jsize3])
            end
            for m in 0..4 do
                i = 1
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @ue[i + m * @jsize3] - 4.0 * @ue[i + 1 + m * @jsize3] + @ue[i + 2 + m * @jsize3])
                i = 2
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @ue[i - 1 + m * @jsize3] + 6.0 * @ue[i + m * @jsize3] - 4.0 * @ue[i + 1 + m * @jsize3] + @ue[i + 2 + m * @jsize3])
            end
            for m in 0..4 do
                for i in 3..@grid_points[0] - 4 do
                    @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[i - 2 + m * @jsize3] - 4.0 * @ue[i - 1 + m * @jsize3] + 6.0 * @ue[i + m * @jsize3] - 4.0 * @ue[i + 1 + m * @jsize3] + @ue[i + 2 + m * @jsize3])
                end
            end
            for m in 0..4 do
                i = @grid_points[0] - 3
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[i - 2 + m * @jsize3] - 4.0 * @ue[i - 1 + m * @jsize3] + 6.0 * @ue[i + m * @jsize3] - 4.0 * @ue[i + 1 + m * @jsize3])
                i = @grid_points[0] - 2
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[i - 2 + m * @jsize3] - 4.0 * @ue[i - 1 + m * @jsize3] + 5.0 * @ue[i + m * @jsize3])
            end
        end
    end
    for k in 1..@grid_points[2] - 2 do
        zeta = k * @@dnzm1
        for i in 1..@grid_points[0] - 2 do
            xi = i * @@dnxm1
            for j in 0..@grid_points[1] - 1 do
                eta = j * @@dnym1
                exact_solution(xi, eta, zeta, dtemp, 0)
                for m in 0..4 do
                    @ue[j + m * @jsize3] = dtemp[m]
                end
                dtpp = 1.0 / dtemp[0]
                for m in 1..4 do
                    @buf[j + m * @jsize3] = dtpp * dtemp[m]
                end
                @cuf[j] = @buf[j + 2 * @jsize3] * @buf[j + 2 * @jsize3]
                @buf[j + 0 * @jsize3] = @cuf[j] + @buf[j + 1 * @jsize3] * @buf[j + 1 * @jsize3] + @buf[j + 3 * @jsize3] * @buf[j + 3 * @jsize3]
                @q[j] = 0.5 * (@buf[j + 1 * @jsize3] * @ue[j + 1 * @jsize3] + @buf[j + 2 * @jsize3] * @ue[j + 2 * @jsize3] + @buf[j + 3 * @jsize3] * @ue[j + 3 * @jsize3])
            end
            for j in 1..@grid_points[1] - 2 do
                jm1 = j - 1
                jp1 = j + 1
                @forcing[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[0 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@ty2 * (@ue[jp1 + 2 * @jsize3] - @ue[jm1 + 2 * @jsize3]) + @@dy1ty1 * (@ue[jp1 + 0 * @jsize3] - 2.0 * @ue[j + 0 * @jsize3] + @ue[jm1 + 0 * @jsize3])
                @forcing[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[1 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@ty2 * (@ue[jp1 + 1 * @jsize3] * @buf[jp1 + 2 * @jsize3] - @ue[jm1 + 1 * @jsize3] * @buf[jm1 + 2 * @jsize3]) + @@yycon2 * (@buf[jp1 + 1 * @jsize3] - 2.0 * @buf[j + 1 * @jsize3] + @buf[jm1 + 1 * @jsize3]) + @@dy2ty1 * (@ue[jp1 + 1 * @jsize3] - 2.0 * @ue[j + 1 * @jsize3] + @ue[jm1 + 1 * @jsize3])
                @forcing[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[2 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@ty2 * ((@ue[jp1 + 2 * @jsize3] * @buf[jp1 + 2 * @jsize3] + @@c2 * (@ue[jp1 + 4 * @jsize3] - @q[jp1])) - (@ue[jm1 + 2 * @jsize3] * @buf[jm1 + 2 * @jsize3] + @@c2 * (@ue[jm1 + 4 * @jsize3] - @q[jm1]))) + @@yycon1 * (@buf[jp1 + 2 * @jsize3] - 2.0 * @buf[j + 2 * @jsize3] + @buf[jm1 + 2 * @jsize3]) + @@dy3ty1 * (@ue[jp1 + 2 * @jsize3] - 2.0 * @ue[j + 2 * @jsize3] + @ue[jm1 + 2 * @jsize3])
                @forcing[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@ty2 * (@ue[jp1 + 3 * @jsize3] * @buf[jp1 + 2 * @jsize3] - @ue[jm1 + 3 * @jsize3] * @buf[jm1 + 2 * @jsize3]) + @@yycon2 * (@buf[jp1 + 3 * @jsize3] - 2.0 * @buf[j + 3 * @jsize3] + @buf[jm1 + 3 * @jsize3]) + @@dy4ty1 * (@ue[jp1 + 3 * @jsize3] - 2.0 * @ue[j + 3 * @jsize3] + @ue[jm1 + 3 * @jsize3])
                @forcing[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@ty2 * (@buf[jp1 + 2 * @jsize3] * (@@c1 * @ue[jp1 + 4 * @jsize3] - @@c2 * @q[jp1]) - @buf[jm1 + 2 * @jsize3] * (@@c1 * @ue[jm1 + 4 * @jsize3] - @@c2 * @q[jm1])) + 0.5 * @@yycon3 * (@buf[jp1 + 0 * @jsize3] - 2.0 * @buf[j + 0 * @jsize3] + @buf[jm1 + 0 * @jsize3]) + @@yycon4 * (@cuf[jp1] - 2.0 * @cuf[j] + @cuf[jm1]) + @@yycon5 * (@buf[jp1 + 4 * @jsize3] - 2.0 * @buf[j + 4 * @jsize3] + @buf[jm1 + 4 * @jsize3]) + @@dy5ty1 * (@ue[jp1 + 4 * @jsize3] - 2.0 * @ue[j + 4 * @jsize3] + @ue[jm1 + 4 * @jsize3])
            end
            for m in 0..4 do
                j = 1
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @ue[j + m * @jsize3] - 4.0 * @ue[j + 1 + m * @jsize3] + @ue[j + 2 + m * @jsize3])
                j = 2
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @ue[j - 1 + m * @jsize3] + 6.0 * @ue[j + m * @jsize3] - 4.0 * @ue[j + 1 + m * @jsize3] + @ue[j + 2 + m * @jsize3])
            end
            for m in 0..4 do
                for j in 3..@grid_points[1] - 4 do
                    @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[j - 2 + m * @jsize3] - 4.0 * @ue[j - 1 + m * @jsize3] + 6.0 * @ue[j + m * @jsize3] - 4.0 * @ue[j + 1 + m * @jsize3] + @ue[j + 2 + m * @jsize3])
                end
            end
            for m in 0..4 do
                j = @grid_points[1] - 3
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[j - 2 + m * @jsize3] - 4.0 * @ue[j - 1 + m * @jsize3] + 6.0 * @ue[j + m * @jsize3] - 4.0 * @ue[j + 1 + m * @jsize3])
                j = @grid_points[1] - 2
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[j - 2 + m * @jsize3] - 4.0 * @ue[j - 1 + m * @jsize3] + 5.0 * @ue[j + m * @jsize3])
            end
        end
    end
    for j in 1..@grid_points[1] - 2 do
        eta = j * @@dnym1
        for i in 1..@grid_points[0] - 2 do
            xi = i * @@dnxm1
            for k in 0..@grid_points[2] - 1 do
                zeta = k * @@dnzm1
                exact_solution(xi, eta, zeta, dtemp, 0)
                for m in 0..4 do
                    @ue[k + m * @jsize3] = dtemp[m]
                end
                dtpp = 1.0 / dtemp[0]
                for m in 1..4 do
                    @buf[k + m * @jsize3] = dtpp * dtemp[m]
                end
                @cuf[k] = @buf[k + 3 * @jsize3] * @buf[k + 3 * @jsize3]
                @buf[k + 0 * @jsize3] = @cuf[k] + @buf[k + 1 * @jsize3] * @buf[k + 1 * @jsize3] + @buf[k + 2 * @jsize3] * @buf[k + 2 * @jsize3]
                @q[k] = 0.5 * (@buf[k + 1 * @jsize3] * @ue[k + 1 * @jsize3] + @buf[k + 2 * @jsize3] * @ue[k + 2 * @jsize3] + @buf[k + 3 * @jsize3] * @ue[k + 3 * @jsize3])
            end
            for k in 1..@grid_points[2] - 2 do
                km1 = k - 1
                kp1 = k + 1
                @forcing[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[0 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tz2 * (@ue[kp1 + 3 * @jsize3] - @ue[km1 + 3 * @jsize3]) + @@dz1tz1 * (@ue[kp1 + 0 * @jsize3] - 2.0 * @ue[k + 0 * @jsize3] + @ue[km1 + 0 * @jsize3])
                @forcing[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[1 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tz2 * (@ue[kp1 + 1 * @jsize3] * @buf[kp1 + 3 * @jsize3] - @ue[km1 + 1 * @jsize3] * @buf[km1 + 3 * @jsize3]) + @@zzcon2 * (@buf[kp1 + 1 * @jsize3] - 2.0 * @buf[k + 1 * @jsize3] + @buf[km1 + 1 * @jsize3]) + @@dz2tz1 * (@ue[kp1 + 1 * @jsize3] - 2.0 * @ue[k + 1 * @jsize3] + @ue[km1 + 1 * @jsize3])
                @forcing[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[2 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tz2 * (@ue[kp1 + 2 * @jsize3] * @buf[kp1 + 3 * @jsize3] - @ue[km1 + 2 * @jsize3] * @buf[km1 + 3 * @jsize3]) + @@zzcon2 * (@buf[kp1 + 2 * @jsize3] - 2.0 * @buf[k + 2 * @jsize3] + @buf[km1 + 2 * @jsize3]) + @@dz3tz1 * (@ue[kp1 + 2 * @jsize3] - 2.0 * @ue[k + 2 * @jsize3] + @ue[km1 + 2 * @jsize3])
                @forcing[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tz2 * ((@ue[kp1 + 3 * @jsize3] * @buf[kp1 + 3 * @jsize3] + @@c2 * (@ue[kp1 + 4 * @jsize3] - @q[kp1])) - (@ue[km1 + 3 * @jsize3] * @buf[km1 + 3 * @jsize3] + @@c2 * (@ue[km1 + 4 * @jsize3] - @q[km1]))) + @@zzcon1 * (@buf[kp1 + 3 * @jsize3] - 2.0 * @buf[k + 3 * @jsize3] + @buf[km1 + 3 * @jsize3]) + @@dz4tz1 * (@ue[kp1 + 3 * @jsize3] - 2.0 * @ue[k + 3 * @jsize3] + @ue[km1 + 3 * @jsize3])
                @forcing[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @@tz2 * (@buf[kp1 + 3 * @jsize3] * (@@c1 * @ue[kp1 + 4 * @jsize3] - @@c2 * @q[kp1]) - @buf[km1 + 3 * @jsize3] * (@@c1 * @ue[km1 + 4 * @jsize3] - @@c2 * @q[km1])) + 0.5 * @@zzcon3 * (@buf[kp1 + 0 * @jsize3] - 2.0 * @buf[k + 0 * @jsize3] + @buf[km1 + 0 * @jsize3]) + @@zzcon4 * (@cuf[kp1] - 2.0 * @cuf[k] + @cuf[km1]) + @@zzcon5 * (@buf[kp1 + 4 * @jsize3] - 2.0 * @buf[k + 4 * @jsize3] + @buf[km1 + 4 * @jsize3]) + @@dz5tz1 * (@ue[kp1 + 4 * @jsize3] - 2.0 * @ue[k + 4 * @jsize3] + @ue[km1 + 4 * @jsize3])
            end
            for m in 0..4 do
                k = 1
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @ue[k + m * @jsize3] - 4.0 * @ue[k + 1 + m * @jsize3] + @ue[k + 2 + m * @jsize3])
                k = 2
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @ue[k - 1 + m * @jsize3] + 6.0 * @ue[k + m * @jsize3] - 4.0 * @ue[k + 1 + m * @jsize3] + @ue[k + 2 + m * @jsize3])
            end
            for m in 0..4 do
                for k in 3..@grid_points[2] - 4 do
                    @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[k - 2 + m * @jsize3] - 4.0 * @ue[k - 1 + m * @jsize3] + 6.0 * @ue[k + m * @jsize3] - 4.0 * @ue[k + 1 + m * @jsize3] + @ue[k + 2 + m * @jsize3])
                end
            end
            for m in 0..4 do
                k = @grid_points[2] - 3
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[k - 2 + m * @jsize3] - 4.0 * @ue[k - 1 + m * @jsize3] + 6.0 * @ue[k + m * @jsize3] - 4.0 * @ue[k + 1 + m * @jsize3])
                k = @grid_points[2] - 2
                @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@ue[k - 2 + m * @jsize3] - 4.0 * @ue[k - 1 + m * @jsize3] + 5.0 * @ue[k + m * @jsize3])
            end
        end
    end
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                for m in 0..4 do
                    @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1] = -1. * @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
            end
        end
    end
  end

  def ninvr()
    
    
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                r1 = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r2 = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r3 = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r4 = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r5 = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                t1 = @@bt * r3
                t2 = 0.5 * (r4 + r5)
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = -r2
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = r1
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @@bt * (r4 - r5)
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = -t1 + t2
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = t1 + t2
            end
        end
    end
  end

  def pinvr()
    
    
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                r1 = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r2 = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r3 = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r4 = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r5 = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                t1 = @@bt * r1
                t2 = 0.5 * (r4 + r5)
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @@bt * (r4 - r5)
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = -r3
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = r2
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = -t1 + t2
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = t1 + t2
            end
        end
    end
  end

  def compute_rhs()
    
    
    for k in 0..@grid_points[2] - 1 do
        for j in 0..@grid_points[1] - 1 do
            for i in 0..@grid_points[0] - 1 do
                rho_inv = 1.0 / @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                @rho_i[i + j * @jsize2 + k * @ksize2] = rho_inv
                @us[i + j * @jsize2 + k * @ksize2] = @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * rho_inv
                @vs[i + j * @jsize2 + k * @ksize2] = @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * rho_inv
                @ws[i + j * @jsize2 + k * @ksize2] = @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * rho_inv
                @square[i + j * @jsize2 + k * @ksize2] = 0.5 * (@u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1]) * rho_inv
                @qs[i + j * @jsize2 + k * @ksize2] = @square[i + j * @jsize2 + k * @ksize2] * rho_inv
                aux = @@c1c2 * rho_inv * (@u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @square[i + j * @jsize2 + k * @ksize2])
                @speed[i + j * @jsize2 + k * @ksize2] = Math.sqrt(aux)
            end
        end
    end
    for k in 0..@grid_points[2] - 1 do
        for j in 0..@grid_points[1] - 1 do
            for i in 0..@grid_points[0] - 1 do
                for m in 0..4 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @forcing[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
            end
        end
    end
    if @timeron then
      @timer.start(@@t_rhsx)
    else
    end
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                uijk = @us[i + j * @jsize2 + k * @ksize2]
                up1 = @us[i + 1 + j * @jsize2 + k * @ksize2]
                um1 = @us[i - 1 + j * @jsize2 + k * @ksize2]
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx1tx1 * (@u[0 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) - @@tx2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1])
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx2tx1 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) + @@xxcon2 * @@con43 * (up1 - 2.0 * uijk + um1) - @@tx2 * (@u[1 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * up1 - @u[1 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * um1 + (@u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - @square[(i + 1) + j * @jsize2 + k * @ksize2] - @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + @square[(i - 1) + j * @jsize2 + k * @ksize2]) * @@c2)
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx3tx1 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) + @@xxcon2 * (@vs[i + 1 + j * @jsize2 + k * @ksize2] - 2.0 * @vs[i + j * @jsize2 + k * @ksize2] + @vs[i - 1 + j * @jsize2 + k * @ksize2]) - @@tx2 * (@u[2 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * up1 - @u[2 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * um1)
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx4tx1 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) + @@xxcon2 * (@ws[i + 1 + j * @jsize2 + k * @ksize2] - 2.0 * @ws[i + j * @jsize2 + k * @ksize2] + @ws[i - 1 + j * @jsize2 + k * @ksize2]) - @@tx2 * (@u[3 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * up1 - @u[3 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * um1)
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dx5tx1 * (@u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1]) + @@xxcon3 * (@qs[i + 1 + j * @jsize2 + k * @ksize2] - 2.0 * @qs[i + j * @jsize2 + k * @ksize2] + @qs[i - 1 + j * @jsize2 + k * @ksize2]) + @@xxcon4 * (up1 * up1 - 2.0 * uijk * uijk + um1 * um1) + @@xxcon5 * (@u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + 1 + j * @jsize2 + k * @ksize2] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize2 + k * @ksize2] + @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i - 1 + j * @jsize2 + k * @ksize2]) - @@tx2 * ((@@c1 * @u[4 + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] - @@c2 * @square[i + 1 + j * @jsize2 + k * @ksize2]) * up1 - (@@c1 * @u[4 + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] - @@c2 * @square[i - 1 + j * @jsize2 + k * @ksize2]) * um1)
            end
            i = 1
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + (i + 2) * @isize1 + j * @jsize1 + k * @ksize1])
            end
            i = 2
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + (i + 2) * @isize1 + j * @jsize1 + k * @ksize1])
            end
            for i in 3..@nx2 - 2 do
                for m in 0..4 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + (i - 2) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1] + @u[m + (i + 2) * @isize1 + j * @jsize1 + k * @ksize1])
                end
            end
            i = @nx2 - 1
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + (i - 2) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i + 1) * @isize1 + j * @jsize1 + k * @ksize1])
            end
            i = @nx2
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + (i - 2) * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + (i - 1) * @isize1 + j * @jsize1 + k * @ksize1] + 5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1])
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_rhsx)
    else
    end
    if @timeron then
      @timer.start(@@t_rhsy)
    else
    end
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                vijk = @vs[i + j * @jsize2 + k * @ksize2]
                vp1 = @vs[i + (j + 1) * @jsize2 + k * @ksize2]
                vm1 = @vs[i + (j - 1) * @jsize2 + k * @ksize2]
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy1ty1 * (@u[0 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) - @@ty2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1])
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy2ty1 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon2 * (@us[i + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @us[i + j * @jsize2 + k * @ksize2] + @us[i + (j - 1) * @jsize2 + k * @ksize2]) - @@ty2 * (@u[1 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * vp1 - @u[1 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * vm1)
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy3ty1 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon2 * @@con43 * (vp1 - 2.0 * vijk + vm1) - @@ty2 * (@u[2 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * vp1 - @u[2 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * vm1 + (@u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - @square[i + (j + 1) * @jsize2 + k * @ksize2] - @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + @square[i + (j - 1) * @jsize2 + k * @ksize2]) * @@c2)
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy4ty1 * (@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon2 * (@ws[i + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @ws[i + j * @jsize2 + k * @ksize2] + @ws[i + (j - 1) * @jsize2 + k * @ksize2]) - @@ty2 * (@u[3 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * vp1 - @u[3 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * vm1)
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dy5ty1 * (@u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon3 * (@qs[i + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @qs[i + j * @jsize2 + k * @ksize2] + @qs[i + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon4 * (vp1 * vp1 - 2.0 * vijk * vijk + vm1 * vm1) + @@yycon5 * (@u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] * @rho_i[i + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize2 + k * @ksize2] + @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] * @rho_i[i + (j - 1) * @jsize2 + k * @ksize2]) - @@ty2 * ((@@c1 * @u[4 + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] - @@c2 * @square[i + (j + 1) * @jsize2 + k * @ksize2]) * vp1 - (@@c1 * @u[4 + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] - @@c2 * @square[i + (j - 1) * @jsize2 + k * @ksize2]) * vm1)
            end
        end
        j = 1
        for i in 1..@nx2 do
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + (j + 2) * @jsize1 + k * @ksize1])
            end
        end
        j = 2
        for i in 1..@nx2 do
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + (j + 2) * @jsize1 + k * @ksize1])
            end
        end
        for j in 3..@ny2 - 2 do
            for i in 1..@nx2 do
                for m in 0..4 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + (j - 2) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1] + @u[m + i * @isize1 + (j + 2) * @jsize1 + k * @ksize1])
                end
            end
        end
        j = @ny2 - 1
        for i in 1..@nx2 do
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + (j - 2) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j + 1) * @jsize1 + k * @ksize1])
            end
        end
        j = @ny2
        for i in 1..@nx2 do
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + (j - 2) * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + (j - 1) * @jsize1 + k * @ksize1] + 5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1])
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_rhsy)
    else
    end
    if @timeron then
      @timer.start(@@t_rhsz)
    else
    end
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                wijk = @ws[i + j * @jsize2 + k * @ksize2]
                wp1 = @ws[i + j * @jsize2 + (k + 1) * @ksize2]
                wm1 = @ws[i + j * @jsize2 + (k - 1) * @ksize2]
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz1tz1 * (@u[0 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[0 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) - @@tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1])
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz2tz1 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[1 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon2 * (@us[i + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @us[i + j * @jsize2 + k * @ksize2] + @us[i + j * @jsize2 + (k - 1) * @ksize2]) - @@tz2 * (@u[1 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * wp1 - @u[1 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * wm1)
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz3tz1 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[2 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon2 * (@vs[i + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @vs[i + j * @jsize2 + k * @ksize2] + @vs[i + j * @jsize2 + (k - 1) * @ksize2]) - @@tz2 * (@u[2 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * wp1 - @u[2 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * wm1)
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz4tz1 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[3 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon2 * @@con43 * (wp1 - 2.0 * wijk + wm1) - @@tz2 * (@u[3 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * wp1 - @u[3 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * wm1 + (@u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - @square[i + j * @jsize2 + (k + 1) * @ksize2] - @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + @square[i + j * @jsize2 + (k - 1) * @ksize2]) * @@c2)
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @@dz5tz1 * (@u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon3 * (@qs[i + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @qs[i + j * @jsize2 + k * @ksize2] + @qs[i + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon4 * (wp1 * wp1 - 2.0 * wijk * wijk + wm1 * wm1) + @@zzcon5 * (@u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] * @rho_i[i + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[4 + i * @isize1 + j * @jsize1 + k * @ksize1] * @rho_i[i + j * @jsize2 + k * @ksize2] + @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] * @rho_i[i + j * @jsize2 + (k - 1) * @ksize2]) - @@tz2 * ((@@c1 * @u[4 + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] - @@c2 * @square[i + j * @jsize2 + (k + 1) * @ksize2]) * wp1 - (@@c1 * @u[4 + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] - @@c2 * @square[i + j * @jsize2 + (k - 1) * @ksize2]) * wm1)
            end
        end
    end
    k = 1
    for j in 1..@ny2 do
        for i in 1..@nx2 do
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + (k + 2) * @ksize1])
            end
        end
    end
    k = 2
    for j in 1..@ny2 do
        for i in 1..@nx2 do
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (-4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + (k + 2) * @ksize1])
            end
        end
    end
    for k in 3..@nz2 - 2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                for m in 0..4 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + j * @jsize1 + (k - 2) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1] + @u[m + i * @isize1 + j * @jsize1 + (k + 2) * @ksize1])
                end
            end
        end
    end
    k = @nz2 - 1
    for j in 1..@ny2 do
        for i in 1..@nx2 do
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + j * @jsize1 + (k - 2) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 6.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k + 1) * @ksize1])
            end
        end
    end
    k = @nz2
    for j in 1..@ny2 do
        for i in 1..@nx2 do
            for m in 0..4 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @@dssp * (@u[m + i * @isize1 + j * @jsize1 + (k - 2) * @ksize1] - 4.0 * @u[m + i * @isize1 + j * @jsize1 + (k - 1) * @ksize1] + 5.0 * @u[m + i * @isize1 + j * @jsize1 + k * @ksize1])
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_rhsz)
    else
    end
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                for m in 0..4 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] * @@dt
                end
            end
        end
    end
  end

  def txinvr()
    
    
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                ru1 = @rho_i[i + j * @jsize2 + k * @ksize2]
                uu = @us[i + j * @jsize2 + k * @ksize2]
                vv = @vs[i + j * @jsize2 + k * @ksize2]
                ww = @ws[i + j * @jsize2 + k * @ksize2]
                ac = @speed[i + j * @jsize2 + k * @ksize2]
                ac2inv = 1.0 / (ac * ac)
                r1 = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r2 = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r3 = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r4 = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r5 = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                t1 = @@c2 * ac2inv * (@qs[i + j * @jsize2 + k * @ksize2] * r1 - uu * r2 - vv * r3 - ww * r4 + r5)
                t2 = @@bt * ru1 * (uu * r1 - r2)
                t3 = (@@bt * ru1 * ac) * t1
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = r1 - t1
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = -ru1 * (ww * r1 - r4)
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = ru1 * (vv * r1 - r3)
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = -t2 + t3
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = t2 + t3
            end
        end
    end
  end

  def tzetar()
    
    
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 1..@nx2 do
                xvel = @us[i + j * @jsize2 + k * @ksize2]
                yvel = @vs[i + j * @jsize2 + k * @ksize2]
                zvel = @ws[i + j * @jsize2 + k * @ksize2]
                ac = @speed[i + j * @jsize2 + k * @ksize2]
                ac2u = ac * ac
                r1 = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r2 = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r3 = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r4 = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                r5 = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
                uzik1 = @u[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                btuz = @@bt * uzik1
                t1 = btuz / ac * (r4 + r5)
                t2 = r3 + t1
                t3 = btuz * (r4 - r5)
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = t2
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = -uzik1 * r2 + xvel * t2
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = uzik1 * r1 + yvel * t2
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = zvel * t2 + t3
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = uzik1 * (-xvel * r2 + yvel * r1) + @qs[i + j * @jsize2 + k * @ksize2] * t2 + @@c2iv * ac2u * t1 + zvel * t3
            end
        end
    end
  end

  def verify(no_time_steps)
    xcrref = Array.new(5, 0.0); xceref = Array.new(5, 0.0); xcrdif = Array.new(5, 0.0); xcedif = Array.new(5, 0.0); xce = Array.new(5, 0.0); xcr = Array.new(5, 0.0); dtref = 0; 
    
    verified = -1; 
    clss = 'U'; 
    error_norm(xce)
    compute_rhs()
    rhs_norm(xcr)
    for m in 0..4 do
      xcr[m] = xcr[m] / @@dt
    end
    for m in 0..4 do
        xcrref[m] = 1.0
        xceref[m] = 1.0
    end
    if @grid_points[0] == 12 and @grid_points[1] == 12 and @grid_points[2] == 12 and no_time_steps == 100 then
        clss = 'S'
        dtref = 0.015
        xcrref[0] = 2.7470315451339479E-2
        xcrref[1] = 1.0360746705285417E-2
        xcrref[2] = 1.6235745065095532E-2
        xcrref[3] = 1.5840557224455615E-2
        xcrref[4] = 3.4849040609362460E-2
        xceref[0] = 2.7289258557377227E-5
        xceref[1] = 1.0364446640837285E-5
        xceref[2] = 1.6154798287166471E-5
        xceref[3] = 1.5750704994480102E-5
        xceref[4] = 3.4177666183390531E-5
    else
      if (@grid_points[0] == 36) and (@grid_points[1] == 36) and (@grid_points[2] == 36) and (no_time_steps == 400) then
          clss = 'W'
          dtref = 0.0015
          xcrref[0] = 0.1893253733584E-2
          xcrref[1] = 0.1717075447775E-3
          xcrref[2] = 0.2778153350936E-3
          xcrref[3] = 0.2887475409984E-3
          xcrref[4] = 0.3143611161242E-2
          xceref[0] = 0.7542088599534E-4
          xceref[1] = 0.6512852253086E-5
          xceref[2] = 0.1049092285688E-4
          xceref[3] = 0.1128838671535E-4
          xceref[4] = 0.1212845639773E-3
      else
        if (@grid_points[0] == 64) and (@grid_points[1] == 64) and (@grid_points[2] == 64) and (no_time_steps == 400) then
            clss = 'A'
            dtref = 0.0015
            xcrref[0] = 2.4799822399300195
            xcrref[1] = 1.1276337964368832
            xcrref[2] = 1.5028977888770491
            xcrref[3] = 1.4217816211695179
            xcrref[4] = 2.1292113035138280
            xceref[0] = 1.0900140297820550E-4
            xceref[1] = 3.7343951769282091E-5
            xceref[2] = 5.0092785406541633E-5
            xceref[3] = 4.7671093939528255E-5
            xceref[4] = 1.3621613399213001E-4
        else
          if (@grid_points[0] == 102) and (@grid_points[1] == 102) and (@grid_points[2] == 102) and (no_time_steps == 400) then
              clss = 'B'
              dtref = 0.001
              xcrref[0] = 0.6903293579998E+02
              xcrref[1] = 0.3095134488084E+02
              xcrref[2] = 0.4103336647017E+02
              xcrref[3] = 0.3864769009604E+02
              xcrref[4] = 0.5643482272596E+02
              xceref[0] = 0.9810006190188E-02
              xceref[1] = 0.1022827905670E-02
              xceref[2] = 0.1720597911692E-02
              xceref[3] = 0.1694479428231E-02
              xceref[4] = 0.1847456263981E-01
          else
            if (@grid_points[0] == 162) and (@grid_points[1] == 162) and (@grid_points[2] == 162) and (no_time_steps == 400) then
                clss = 'C'
                dtref = 0.00067
                xcrref[0] = 0.5881691581829E+03
                xcrref[1] = 0.2454417603569E+03
                xcrref[2] = 0.3293829191851E+03
                xcrref[3] = 0.3081924971891E+03
                xcrref[4] = 0.4597223799176E+03
                xceref[0] = 0.2598120500183
                xceref[1] = 0.2590888922315E-01
                xceref[2] = 0.5132886416320E-01
                xceref[3] = 0.4806073419454E-01
                xceref[4] = 0.5483377491301
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
    epsilon = 1.0E-8; 
    if clss != 'U' then
        puts(" Verification being performed for class " + clss)
        puts(" Accuracy setting for epsilon = " + epsilon.to_s)
        if Math_dot_abs(@@dt - dtref) <= epsilon then
            if verified == -1 then
              verified = 1
            else
            end
        else
            verified = 0
            clss = 'U'
            puts("DT does not match the reference value of " + dtref.to_s)
        end
        puts(" Comparison of RMS-norms of residual")
    else
        puts(" Unknown CLASS")
        puts(" RMS-norms of residual")
    end
    verified = BMResults.printComparisonStatusArr(clss, verified, epsilon, xcr, xcrref, xcrdif)
    if clss != 'U' then
        puts(" Comparison of RMS-norms of solution error")
    else
        puts(" RMS-norms of solution error")
    end
    verified = BMResults.printComparisonStatusArr(clss, verified, epsilon, xce, xceref, xcedif)
    BMResults.printVerificationStatus(clss, verified, @@BMName)
    return verified
  end

  def x_solve()
    
    
    if @timeron then
      @timer.start(@@t_xsolve)
    else
    end
    for k in 1..@nz2 do
        for j in 1..@ny2 do
            for i in 0..@grid_points[0] - 1 do
                ru1 = @@c3c4 * @rho_i[i + j * @jsize2 + k * @ksize2]
                @cv[i] = @us[i + j * @jsize2 + k * @ksize2]
                @rhon[i] = dmax1ex(@@dx2 + @@con43 * ru1, @@dx5 + @@c1c5 * ru1, @@dxmax + ru1, @@dx1)
            end
            lhsinit(@grid_points[0] - 1)
            for i in 1..@nx2 do
                @lhs[0 + i * @jsize4] = 0.0
                @lhs[1 + i * @jsize4] = -@@dttx2 * @cv[(i - 1)] - @@dttx1 * @rhon[i - 1]
                @lhs[2 + i * @jsize4] = 1.0 + @@c2dttx1 * @rhon[i]
                @lhs[3 + i * @jsize4] = @@dttx2 * @cv[i + 1] - @@dttx1 * @rhon[i + 1]
                @lhs[4 + i * @jsize4] = 0.0
            end
            i = 1
            @lhs[2 + i * @jsize4] = @lhs[2 + i * @jsize4] + @@comz5
            @lhs[3 + i * @jsize4] = @lhs[3 + i * @jsize4] - @@comz4
            @lhs[4 + i * @jsize4] = @lhs[4 + i * @jsize4] + @@comz1
            @lhs[1 + (i + 1) * @jsize4] = @lhs[1 + (i + 1) * @jsize4] - @@comz4
            @lhs[2 + (i + 1) * @jsize4] = @lhs[2 + (i + 1) * @jsize4] + @@comz6
            @lhs[3 + (i + 1) * @jsize4] = @lhs[3 + (i + 1) * @jsize4] - @@comz4
            @lhs[4 + (i + 1) * @jsize4] = @lhs[4 + (i + 1) * @jsize4] + @@comz1
            for i in 3..@grid_points[0] - 4 do
                @lhs[0 + i * @jsize4] = @lhs[0 + i * @jsize4] + @@comz1
                @lhs[1 + i * @jsize4] = @lhs[1 + i * @jsize4] - @@comz4
                @lhs[2 + i * @jsize4] = @lhs[2 + i * @jsize4] + @@comz6
                @lhs[3 + i * @jsize4] = @lhs[3 + i * @jsize4] - @@comz4
                @lhs[4 + i * @jsize4] = @lhs[4 + i * @jsize4] + @@comz1
            end
            i = @grid_points[0] - 3
            @lhs[0 + i * @jsize4] = @lhs[0 + i * @jsize4] + @@comz1
            @lhs[1 + i * @jsize4] = @lhs[1 + i * @jsize4] - @@comz4
            @lhs[2 + i * @jsize4] = @lhs[2 + i * @jsize4] + @@comz6
            @lhs[3 + i * @jsize4] = @lhs[3 + i * @jsize4] - @@comz4
            @lhs[0 + (i + 1) * @jsize4] = @lhs[0 + (i + 1) * @jsize4] + @@comz1
            @lhs[1 + (i + 1) * @jsize4] = @lhs[1 + (i + 1) * @jsize4] - @@comz4
            @lhs[2 + (i + 1) * @jsize4] = @lhs[2 + (i + 1) * @jsize4] + @@comz5
            for i in 1..@nx2 do
                @lhsp[0 + i * @jsize4] = @lhs[0 + i * @jsize4]
                @lhsp[1 + i * @jsize4] = @lhs[1 + i * @jsize4] - @@dttx2 * @speed[(i - 1) + j * @jsize2 + k * @ksize2]
                @lhsp[2 + i * @jsize4] = @lhs[2 + i * @jsize4]
                @lhsp[3 + i * @jsize4] = @lhs[3 + i * @jsize4] + @@dttx2 * @speed[i + 1 + j * @jsize2 + k * @ksize2]
                @lhsp[4 + i * @jsize4] = @lhs[4 + i * @jsize4]
                @lhsm[0 + i * @jsize4] = @lhs[0 + i * @jsize4]
                @lhsm[1 + i * @jsize4] = @lhs[1 + i * @jsize4] + @@dttx2 * @speed[i - 1 + j * @jsize2 + k * @ksize2]
                @lhsm[2 + i * @jsize4] = @lhs[2 + i * @jsize4]
                @lhsm[3 + i * @jsize4] = @lhs[3 + i * @jsize4] - @@dttx2 * @speed[i + 1 + j * @jsize2 + k * @ksize2]
                @lhsm[4 + i * @jsize4] = @lhs[4 + i * @jsize4]
            end
            for i in 0..@grid_points[0] - 3 do
                i1 = i + 1
                i2 = i + 2
                fac1 = 1. / @lhs[2 + i * @jsize4]
                @lhs[3 + i * @jsize4] = fac1 * @lhs[3 + i * @jsize4]
                @lhs[4 + i * @jsize4] = fac1 * @lhs[4 + i * @jsize4]
                for m in 0..2 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
                @lhs[2 + i1 * @jsize4] = @lhs[2 + i1 * @jsize4] - @lhs[1 + i1 * @jsize4] * @lhs[3 + i * @jsize4]
                @lhs[3 + i1 * @jsize4] = @lhs[3 + i1 * @jsize4] - @lhs[1 + i1 * @jsize4] * @lhs[4 + i * @jsize4]
                for m in 0..2 do
                    @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
                @lhs[1 + i2 * @jsize4] = @lhs[1 + i2 * @jsize4] - @lhs[0 + i2 * @jsize4] * @lhs[3 + i * @jsize4]
                @lhs[2 + i2 * @jsize4] = @lhs[2 + i2 * @jsize4] - @lhs[0 + i2 * @jsize4] * @lhs[4 + i * @jsize4]
                for m in 0..2 do
                    @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[0 + i2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
            end
            i = @grid_points[0] - 2
            i1 = @grid_points[0] - 1
            fac1 = 1. / @lhs[2 + i * @jsize4]
            @lhs[3 + i * @jsize4] = fac1 * @lhs[3 + i * @jsize4]
            @lhs[4 + i * @jsize4] = fac1 * @lhs[4 + i * @jsize4]
            for m in 0..2 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            end
            @lhs[2 + i1 * @jsize4] = @lhs[2 + i1 * @jsize4] - @lhs[1 + i1 * @jsize4] * @lhs[3 + i * @jsize4]
            @lhs[3 + i1 * @jsize4] = @lhs[3 + i1 * @jsize4] - @lhs[1 + i1 * @jsize4] * @lhs[4 + i * @jsize4]
            for m in 0..2 do
                @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            end
            fac2 = 1. / @lhs[2 + i1 * @jsize4]
            for m in 0..2 do
                @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = fac2 * @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1]
            end
            for i in 0..@grid_points[0] - 3 do
                i1 = i + 1
                i2 = i + 2
                m = 3
                fac1 = 1. / @lhsp[2 + i * @jsize4]
                @lhsp[3 + i * @jsize4] = fac1 * @lhsp[3 + i * @jsize4]
                @lhsp[4 + i * @jsize4] = fac1 * @lhsp[4 + i * @jsize4]
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                @lhsp[2 + i1 * @jsize4] = @lhsp[2 + i1 * @jsize4] - @lhsp[1 + i1 * @jsize4] * @lhsp[3 + i * @jsize4]
                @lhsp[3 + i1 * @jsize4] = @lhsp[3 + i1 * @jsize4] - @lhsp[1 + i1 * @jsize4] * @lhsp[4 + i * @jsize4]
                @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                @lhsp[1 + i2 * @jsize4] = @lhsp[1 + i2 * @jsize4] - @lhsp[0 + i2 * @jsize4] * @lhsp[3 + i * @jsize4]
                @lhsp[2 + i2 * @jsize4] = @lhsp[2 + i2 * @jsize4] - @lhsp[0 + i2 * @jsize4] * @lhsp[4 + i * @jsize4]
                @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[0 + i2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                m = 4
                fac1 = 1. / @lhsm[2 + i * @jsize4]
                @lhsm[3 + i * @jsize4] = fac1 * @lhsm[3 + i * @jsize4]
                @lhsm[4 + i * @jsize4] = fac1 * @lhsm[4 + i * @jsize4]
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                @lhsm[2 + i1 * @jsize4] = @lhsm[2 + i1 * @jsize4] - @lhsm[1 + i1 * @jsize4] * @lhsm[3 + i * @jsize4]
                @lhsm[3 + i1 * @jsize4] = @lhsm[3 + i1 * @jsize4] - @lhsm[1 + i1 * @jsize4] * @lhsm[4 + i * @jsize4]
                @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                @lhsm[1 + i2 * @jsize4] = @lhsm[1 + i2 * @jsize4] - @lhsm[0 + i2 * @jsize4] * @lhsm[3 + i * @jsize4]
                @lhsm[2 + i2 * @jsize4] = @lhsm[2 + i2 * @jsize4] - @lhsm[0 + i2 * @jsize4] * @lhsm[4 + i * @jsize4]
                @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[0 + i2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            end
            i = @grid_points[0] - 2
            i1 = @grid_points[0] - 1
            m = 3
            fac1 = 1. / @lhsp[2 + i * @jsize4]
            @lhsp[3 + i * @jsize4] = fac1 * @lhsp[3 + i * @jsize4]
            @lhsp[4 + i * @jsize4] = fac1 * @lhsp[4 + i * @jsize4]
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            @lhsp[2 + i1 * @jsize4] = @lhsp[2 + i1 * @jsize4] - @lhsp[1 + i1 * @jsize4] * @lhsp[3 + i * @jsize4]
            @lhsp[3 + i1 * @jsize4] = @lhsp[3 + i1 * @jsize4] - @lhsp[1 + i1 * @jsize4] * @lhsp[4 + i * @jsize4]
            @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            m = 4
            fac1 = 1. / @lhsm[2 + i * @jsize4]
            @lhsm[3 + i * @jsize4] = fac1 * @lhsm[3 + i * @jsize4]
            @lhsm[4 + i * @jsize4] = fac1 * @lhsm[4 + i * @jsize4]
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            @lhsm[2 + i1 * @jsize4] = @lhsm[2 + i1 * @jsize4] - @lhsm[1 + i1 * @jsize4] * @lhsm[3 + i * @jsize4]
            @lhsm[3 + i1 * @jsize4] = @lhsm[3 + i1 * @jsize4] - @lhsm[1 + i1 * @jsize4] * @lhsm[4 + i * @jsize4]
            @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[1 + i1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            @rhs[3 + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i1 * @isize1 + j * @jsize1 + k * @ksize1] / @lhsp[2 + i1 * @jsize4]
            @rhs[4 + i1 * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i1 * @isize1 + j * @jsize1 + k * @ksize1] / @lhsm[2 + i1 * @jsize4]
            i = @grid_points[0] - 2
            i1 = @grid_points[0] - 1
            for m in 0..2 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[3 + i * @jsize4] * @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1]
            end
            @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[3 + i * @jsize4] * @rhs[3 + i1 * @isize1 + j * @jsize1 + k * @ksize1]
            @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[3 + i * @jsize4] * @rhs[4 + i1 * @isize1 + j * @jsize1 + k * @ksize1]
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
            i = @grid_points[0] - 3
            while i >= 0 do
                i1 = i + 1
                i2 = i + 2
                for m in 0..2 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[3 + i * @jsize4] * @rhs[m + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[4 + i * @jsize4] * @rhs[m + i2 * @isize1 + j * @jsize1 + k * @ksize1]
                end
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[3 + i * @jsize4] * @rhs[3 + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[4 + i * @jsize4] * @rhs[3 + i2 * @isize1 + j * @jsize1 + k * @ksize1]
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[3 + i * @jsize4] * @rhs[4 + i1 * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[4 + i * @jsize4] * @rhs[4 + i2 * @isize1 + j * @jsize1 + k * @ksize1]
              i -= 1
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_xsolve)
    else
    end
    if @timeron then
      @timer.start(@@t_ninvr)
    else
    end
    ninvr()
    if @timeron then
      @timer.stop(@@t_ninvr)
    else
    end
  end

  def y_solve()
    
    
    if @timeron then
      @timer.start(@@t_ysolve)
    else
    end
    for k in 1..@grid_points[2] - 2 do
        for i in 1..@grid_points[0] - 2 do
            for j in 0..@grid_points[1] - 1 do
                ru1 = @@c3c4 * @rho_i[i + j * @jsize2 + k * @ksize2]
                @cv[j] = @vs[i + j * @jsize2 + k * @ksize2]
                @rhoq[j] = dmax1ex(@@dy3 + @@con43 * ru1, @@dy5 + @@c1c5 * ru1, @@dymax + ru1, @@dy1)
            end
            lhsinit(@grid_points[1] - 1)
            for j in 1..@grid_points[1] - 2 do
                @lhs[0 + j * @jsize4] = 0.0
                @lhs[1 + j * @jsize4] = -@@dtty2 * @cv[j - 1] - @@dtty1 * @rhoq[j - 1]
                @lhs[2 + j * @jsize4] = 1.0 + @@c2dtty1 * @rhoq[j]
                @lhs[3 + j * @jsize4] = @@dtty2 * @cv[j + 1] - @@dtty1 * @rhoq[j + 1]
                @lhs[4 + j * @jsize4] = 0.0
            end
            j = 1
            @lhs[2 + j * @jsize4] = @lhs[2 + j * @jsize4] + @@comz5
            @lhs[3 + j * @jsize4] = @lhs[3 + j * @jsize4] - @@comz4
            @lhs[4 + j * @jsize4] = @lhs[4 + j * @jsize4] + @@comz1
            @lhs[1 + (j + 1) * @jsize4] = @lhs[1 + (j + 1) * @jsize4] - @@comz4
            @lhs[2 + (j + 1) * @jsize4] = @lhs[2 + (j + 1) * @jsize4] + @@comz6
            @lhs[3 + (j + 1) * @jsize4] = @lhs[3 + (j + 1) * @jsize4] - @@comz4
            @lhs[4 + (j + 1) * @jsize4] = @lhs[4 + (j + 1) * @jsize4] + @@comz1
            for j in 3..@grid_points[1] - 4 do
                @lhs[0 + j * @jsize4] = @lhs[0 + j * @jsize4] + @@comz1
                @lhs[1 + j * @jsize4] = @lhs[1 + j * @jsize4] - @@comz4
                @lhs[2 + j * @jsize4] = @lhs[2 + j * @jsize4] + @@comz6
                @lhs[3 + j * @jsize4] = @lhs[3 + j * @jsize4] - @@comz4
                @lhs[4 + j * @jsize4] = @lhs[4 + j * @jsize4] + @@comz1
            end
            j = @grid_points[1] - 3
            @lhs[0 + j * @jsize4] = @lhs[0 + j * @jsize4] + @@comz1
            @lhs[1 + j * @jsize4] = @lhs[1 + j * @jsize4] - @@comz4
            @lhs[2 + j * @jsize4] = @lhs[2 + j * @jsize4] + @@comz6
            @lhs[3 + j * @jsize4] = @lhs[3 + j * @jsize4] - @@comz4
            @lhs[0 + (j + 1) * @jsize4] = @lhs[0 + (j + 1) * @jsize4] + @@comz1
            @lhs[1 + (j + 1) * @jsize4] = @lhs[1 + (j + 1) * @jsize4] - @@comz4
            @lhs[2 + (j + 1) * @jsize4] = @lhs[2 + (j + 1) * @jsize4] + @@comz5
            for j in 1..@grid_points[1] - 2 do
                @lhsp[0 + j * @jsize4] = @lhs[0 + j * @jsize4]
                @lhsp[1 + j * @jsize4] = @lhs[1 + j * @jsize4] - @@dtty2 * @speed[i + (j - 1) * @jsize2 + k * @ksize2]
                @lhsp[2 + j * @jsize4] = @lhs[2 + j * @jsize4]
                @lhsp[3 + j * @jsize4] = @lhs[3 + j * @jsize4] + @@dtty2 * @speed[i + (j + 1) * @jsize2 + k * @ksize2]
                @lhsp[4 + j * @jsize4] = @lhs[4 + j * @jsize4]
                @lhsm[0 + j * @jsize4] = @lhs[0 + j * @jsize4]
                @lhsm[1 + j * @jsize4] = @lhs[1 + j * @jsize4] + @@dtty2 * @speed[i + (j - 1) * @jsize2 + k * @ksize2]
                @lhsm[2 + j * @jsize4] = @lhs[2 + j * @jsize4]
                @lhsm[3 + j * @jsize4] = @lhs[3 + j * @jsize4] - @@dtty2 * @speed[i + (j + 1) * @jsize2 + k * @ksize2]
                @lhsm[4 + j * @jsize4] = @lhs[4 + j * @jsize4]
            end
            for j in 0..@grid_points[1] - 3 do
                j1 = j + 1
                j2 = j + 2
                fac1 = 1. / @lhs[2 + j * @jsize4]
                @lhs[3 + j * @jsize4] = fac1 * @lhs[3 + j * @jsize4]
                @lhs[4 + j * @jsize4] = fac1 * @lhs[4 + j * @jsize4]
                for m in 0..2 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
                @lhs[2 + j1 * @jsize4] = @lhs[2 + j1 * @jsize4] - @lhs[1 + j1 * @jsize4] * @lhs[3 + j * @jsize4]
                @lhs[3 + j1 * @jsize4] = @lhs[3 + j1 * @jsize4] - @lhs[1 + j1 * @jsize4] * @lhs[4 + j * @jsize4]
                for m in 0..2 do
                    @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhs[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
                @lhs[1 + j2 * @jsize4] = @lhs[1 + j2 * @jsize4] - @lhs[0 + j2 * @jsize4] * @lhs[3 + j * @jsize4]
                @lhs[2 + j2 * @jsize4] = @lhs[2 + j2 * @jsize4] - @lhs[0 + j2 * @jsize4] * @lhs[4 + j * @jsize4]
                for m in 0..2 do
                    @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] - @lhs[0 + j2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                end
            end
            j = @grid_points[1] - 2
            j1 = @grid_points[1] - 1
            fac1 = 1. / @lhs[2 + j * @jsize4]
            @lhs[3 + j * @jsize4] = fac1 * @lhs[3 + j * @jsize4]
            @lhs[4 + j * @jsize4] = fac1 * @lhs[4 + j * @jsize4]
            for m in 0..2 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            end
            @lhs[2 + j1 * @jsize4] = @lhs[2 + j1 * @jsize4] - @lhs[1 + j1 * @jsize4] * @lhs[3 + j * @jsize4]
            @lhs[3 + j1 * @jsize4] = @lhs[3 + j1 * @jsize4] - @lhs[1 + j1 * @jsize4] * @lhs[4 + j * @jsize4]
            for m in 0..2 do
                @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhs[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            end
            fac2 = 1. / @lhs[2 + j1 * @jsize4]
            for m in 0..2 do
                @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = fac2 * @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1]
            end
            for j in 0..@grid_points[1] - 3 do
                j1 = j + 1
                j2 = j + 2
                m = 3
                fac1 = 1. / @lhsp[2 + j * @jsize4]
                @lhsp[3 + j * @jsize4] = fac1 * @lhsp[3 + j * @jsize4]
                @lhsp[4 + j * @jsize4] = fac1 * @lhsp[4 + j * @jsize4]
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                @lhsp[2 + j1 * @jsize4] = @lhsp[2 + j1 * @jsize4] - @lhsp[1 + j1 * @jsize4] * @lhsp[3 + j * @jsize4]
                @lhsp[3 + j1 * @jsize4] = @lhsp[3 + j1 * @jsize4] - @lhsp[1 + j1 * @jsize4] * @lhsp[4 + j * @jsize4]
                @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsp[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                @lhsp[1 + j2 * @jsize4] = @lhsp[1 + j2 * @jsize4] - @lhsp[0 + j2 * @jsize4] * @lhsp[3 + j * @jsize4]
                @lhsp[2 + j2 * @jsize4] = @lhsp[2 + j2 * @jsize4] - @lhsp[0 + j2 * @jsize4] * @lhsp[4 + j * @jsize4]
                @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] - @lhsp[0 + j2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                m = 4
                fac1 = 1. / @lhsm[2 + j * @jsize4]
                @lhsm[3 + j * @jsize4] = fac1 * @lhsm[3 + j * @jsize4]
                @lhsm[4 + j * @jsize4] = fac1 * @lhsm[4 + j * @jsize4]
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                @lhsm[2 + j1 * @jsize4] = @lhsm[2 + j1 * @jsize4] - @lhsm[1 + j1 * @jsize4] * @lhsm[3 + j * @jsize4]
                @lhsm[3 + j1 * @jsize4] = @lhsm[3 + j1 * @jsize4] - @lhsm[1 + j1 * @jsize4] * @lhsm[4 + j * @jsize4]
                @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsm[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
                @lhsm[1 + j2 * @jsize4] = @lhsm[1 + j2 * @jsize4] - @lhsm[0 + j2 * @jsize4] * @lhsm[3 + j * @jsize4]
                @lhsm[2 + j2 * @jsize4] = @lhsm[2 + j2 * @jsize4] - @lhsm[0 + j2 * @jsize4] * @lhsm[4 + j * @jsize4]
                @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1] - @lhsm[0 + j2 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            end
            j = @grid_points[1] - 2
            j1 = @grid_points[1] - 1
            m = 3
            fac1 = 1. / @lhsp[2 + j * @jsize4]
            @lhsp[3 + j * @jsize4] = fac1 * @lhsp[3 + j * @jsize4]
            @lhsp[4 + j * @jsize4] = fac1 * @lhsp[4 + j * @jsize4]
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            @lhsp[2 + j1 * @jsize4] = @lhsp[2 + j1 * @jsize4] - @lhsp[1 + j1 * @jsize4] * @lhsp[3 + j * @jsize4]
            @lhsp[3 + j1 * @jsize4] = @lhsp[3 + j1 * @jsize4] - @lhsp[1 + j1 * @jsize4] * @lhsp[4 + j * @jsize4]
            @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsp[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            m = 4
            fac1 = 1. / @lhsm[2 + j * @jsize4]
            @lhsm[3 + j * @jsize4] = fac1 * @lhsm[3 + j * @jsize4]
            @lhsm[4 + j * @jsize4] = fac1 * @lhsm[4 + j * @jsize4]
            @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = fac1 * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            @lhsm[2 + j1 * @jsize4] = @lhsm[2 + j1 * @jsize4] - @lhsm[1 + j1 * @jsize4] * @lhsm[3 + j * @jsize4]
            @lhsm[3 + j1 * @jsize4] = @lhsm[3 + j1 * @jsize4] - @lhsm[1 + j1 * @jsize4] * @lhsm[4 + j * @jsize4]
            @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsm[1 + j1 * @jsize4] * @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1]
            @rhs[3 + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j1 * @jsize1 + k * @ksize1] / @lhsp[2 + j1 * @jsize4]
            @rhs[4 + i * @isize1 + j1 * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j1 * @jsize1 + k * @ksize1] / @lhsm[2 + j1 * @jsize4]
            j = @grid_points[1] - 2
            j1 = @grid_points[1] - 1
            for m in 0..2 do
                @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[3 + j * @jsize4] * @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1]
            end
            @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[3 + j * @jsize4] * @rhs[3 + i * @isize1 + j1 * @jsize1 + k * @ksize1]
            @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[3 + j * @jsize4] * @rhs[4 + i * @isize1 + j1 * @jsize1 + k * @ksize1]
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
            j = @grid_points[1] - 3
            while j >= 0 do
                j1 = j + 1
                j2 = j + 2
                for m in 0..2 do
                    @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[m + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhs[3 + j * @jsize4] * @rhs[m + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhs[4 + j * @jsize4] * @rhs[m + i * @isize1 + j2 * @jsize1 + k * @ksize1]
                end
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsp[3 + j * @jsize4] * @rhs[3 + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsp[4 + j * @jsize4] * @rhs[3 + i * @isize1 + j2 * @jsize1 + k * @ksize1]
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] - @lhsm[3 + j * @jsize4] * @rhs[4 + i * @isize1 + j1 * @jsize1 + k * @ksize1] - @lhsm[4 + j * @jsize4] * @rhs[4 + i * @isize1 + j2 * @jsize1 + k * @ksize1]
              j -= 1
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_ysolve)
    else
    end
    if @timeron then
      @timer.start(@@t_pinvr)
    else
    end
    pinvr()
    if @timeron then
      @timer.stop(@@t_pinvr)
    else
    end
  end

  def z_solve()
    
    rtmp = Array.new(5 * (@KMAX + 1), 0.0); 
    if @timeron then
      @timer.start(@@t_zsolve)
    else
    end
    for j in 1..@ny2 do
        for i in 1..@nx2 do
            for k in 0..@nz2 + 1 do
                ru1 = @@c3c4 * @rho_i[i + j * @jsize2 + k * @ksize2]
                @cv[k] = @ws[i + j * @jsize2 + k * @ksize2]
                @rhos[k] = dmax1ex(@@dz4 + @@con43 * ru1, @@dz5 + @@c1c5 * ru1, @@dzmax + ru1, @@dz1)
            end
            lhsinit(@grid_points[2] - 1)
            for k in 1..@nz2 do
                @lhs[0 + k * @jsize4] = 0.0
                @lhs[1 + k * @jsize4] = -@@dttz2 * @cv[k - 1] - @@dttz1 * @rhos[k - 1]
                @lhs[2 + k * @jsize4] = 1.0 + @@c2dttz1 * @rhos[k]
                @lhs[3 + k * @jsize4] = @@dttz2 * @cv[k + 1] - @@dttz1 * @rhos[k + 1]
                @lhs[4 + k * @jsize4] = 0.0
            end
            k = 1
            @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz5
            @lhs[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@comz4
            @lhs[4 + k * @jsize4] = @lhs[4 + k * @jsize4] + @@comz1
            k = 2
            @lhs[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@comz4
            @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz6
            @lhs[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@comz4
            @lhs[4 + k * @jsize4] = @lhs[4 + k * @jsize4] + @@comz1
            for k in 3..@nz2 - 2 do
                @lhs[0 + k * @jsize4] = @lhs[0 + k * @jsize4] + @@comz1
                @lhs[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@comz4
                @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz6
                @lhs[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@comz4
                @lhs[4 + k * @jsize4] = @lhs[4 + k * @jsize4] + @@comz1
            end
            k = @nz2 - 1
            @lhs[0 + k * @jsize4] = @lhs[0 + k * @jsize4] + @@comz1
            @lhs[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@comz4
            @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz6
            @lhs[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@comz4
            k = @nz2
            @lhs[0 + k * @jsize4] = @lhs[0 + k * @jsize4] + @@comz1
            @lhs[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@comz4
            @lhs[2 + k * @jsize4] = @lhs[2 + k * @jsize4] + @@comz5
            for k in 1..@nz2 do
                @lhsp[0 + k * @jsize4] = @lhs[0 + k * @jsize4]
                @lhsp[1 + k * @jsize4] = @lhs[1 + k * @jsize4] - @@dttz2 * @speed[i + j * @jsize2 + (k - 1) * @ksize2]
                @lhsp[2 + k * @jsize4] = @lhs[2 + k * @jsize4]
                @lhsp[3 + k * @jsize4] = @lhs[3 + k * @jsize4] + @@dttz2 * @speed[i + j * @jsize2 + (k + 1) * @ksize2]
                @lhsp[4 + k * @jsize4] = @lhs[4 + k * @jsize4]
                @lhsm[0 + k * @jsize4] = @lhs[0 + k * @jsize4]
                @lhsm[1 + k * @jsize4] = @lhs[1 + k * @jsize4] + @@dttz2 * @speed[i + j * @jsize2 + (k - 1) * @ksize2]
                @lhsm[2 + k * @jsize4] = @lhs[2 + k * @jsize4]
                @lhsm[3 + k * @jsize4] = @lhs[3 + k * @jsize4] - @@dttz2 * @speed[i + j * @jsize2 + (k + 1) * @ksize2]
                @lhsm[4 + k * @jsize4] = @lhs[4 + k * @jsize4]
            end
            for k in 0..@nz2 + 1 do
                rtmp[0 + k * 5] = @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1]
                rtmp[1 + k * 5] = @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1]
                rtmp[2 + k * 5] = @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1]
                rtmp[3 + k * 5] = @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1]
                rtmp[4 + k * 5] = @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1]
            end
            for k in 0..@grid_points[2] - 3 do
                k1 = k + 1
                k2 = k + 2
                fac1 = 1. / @lhs[2 + k * @jsize4]
                @lhs[3 + k * @jsize4] = fac1 * @lhs[3 + k * @jsize4]
                @lhs[4 + k * @jsize4] = fac1 * @lhs[4 + k * @jsize4]
                for m in 0..2 do
                    rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
                end
                @lhs[2 + k1 * @jsize4] = @lhs[2 + k1 * @jsize4] - @lhs[1 + k1 * @jsize4] * @lhs[3 + k * @jsize4]
                @lhs[3 + k1 * @jsize4] = @lhs[3 + k1 * @jsize4] - @lhs[1 + k1 * @jsize4] * @lhs[4 + k * @jsize4]
                for m in 0..2 do
                    rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhs[1 + k1 * @jsize4] * rtmp[m + k * 5]
                end
                @lhs[1 + k2 * @jsize4] = @lhs[1 + k2 * @jsize4] - @lhs[0 + k2 * @jsize4] * @lhs[3 + k * @jsize4]
                @lhs[2 + k2 * @jsize4] = @lhs[2 + k2 * @jsize4] - @lhs[0 + k2 * @jsize4] * @lhs[4 + k * @jsize4]
                for m in 0..2 do
                    rtmp[m + k2 * 5] = rtmp[m + k2 * 5] - @lhs[0 + k2 * @jsize4] * rtmp[m + k * 5]
                end
            end
            k = @grid_points[2] - 2
            k1 = @grid_points[2] - 1
            fac1 = 1. / @lhs[2 + k * @jsize4]
            @lhs[3 + k * @jsize4] = fac1 * @lhs[3 + k * @jsize4]
            @lhs[4 + k * @jsize4] = fac1 * @lhs[4 + k * @jsize4]
            for m in 0..2 do
                rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
            end
            @lhs[2 + k1 * @jsize4] = @lhs[2 + k1 * @jsize4] - @lhs[1 + k1 * @jsize4] * @lhs[3 + k * @jsize4]
            @lhs[3 + k1 * @jsize4] = @lhs[3 + k1 * @jsize4] - @lhs[1 + k1 * @jsize4] * @lhs[4 + k * @jsize4]
            for m in 0..2 do
                rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhs[1 + k1 * @jsize4] * rtmp[m + k * 5]
            end
            fac2 = 1. / @lhs[2 + k1 * @jsize4]
            for m in 0..2 do
                rtmp[m + k1 * 5] = fac2 * rtmp[m + k1 * 5]
            end
            for k in 0..@grid_points[2] - 3 do
                k1 = k + 1
                k2 = k + 2
                m = 3
                fac1 = 1. / @lhsp[2 + k * @jsize4]
                @lhsp[3 + k * @jsize4] = fac1 * @lhsp[3 + k * @jsize4]
                @lhsp[4 + k * @jsize4] = fac1 * @lhsp[4 + k * @jsize4]
                rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
                @lhsp[2 + k1 * @jsize4] = @lhsp[2 + k1 * @jsize4] - @lhsp[1 + k1 * @jsize4] * @lhsp[3 + k * @jsize4]
                @lhsp[3 + k1 * @jsize4] = @lhsp[3 + k1 * @jsize4] - @lhsp[1 + k1 * @jsize4] * @lhsp[4 + k * @jsize4]
                rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhsp[1 + k1 * @jsize4] * rtmp[m + k * 5]
                @lhsp[1 + k2 * @jsize4] = @lhsp[1 + k2 * @jsize4] - @lhsp[0 + k2 * @jsize4] * @lhsp[3 + k * @jsize4]
                @lhsp[2 + k2 * @jsize4] = @lhsp[2 + k2 * @jsize4] - @lhsp[0 + k2 * @jsize4] * @lhsp[4 + k * @jsize4]
                rtmp[m + k2 * 5] = rtmp[m + k2 * 5] - @lhsp[0 + k2 * @jsize4] * rtmp[m + k * 5]
                m = 4
                fac1 = 1. / @lhsm[2 + k * @jsize4]
                @lhsm[3 + k * @jsize4] = fac1 * @lhsm[3 + k * @jsize4]
                @lhsm[4 + k * @jsize4] = fac1 * @lhsm[4 + k * @jsize4]
                rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
                @lhsm[2 + k1 * @jsize4] = @lhsm[2 + k1 * @jsize4] - @lhsm[1 + k1 * @jsize4] * @lhsm[3 + k * @jsize4]
                @lhsm[3 + k1 * @jsize4] = @lhsm[3 + k1 * @jsize4] - @lhsm[1 + k1 * @jsize4] * @lhsm[4 + k * @jsize4]
                rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhsm[1 + k1 * @jsize4] * rtmp[m + k * 5]
                @lhsm[1 + k2 * @jsize4] = @lhsm[1 + k2 * @jsize4] - @lhsm[0 + k2 * @jsize4] * @lhsm[3 + k * @jsize4]
                @lhsm[2 + k2 * @jsize4] = @lhsm[2 + k2 * @jsize4] - @lhsm[0 + k2 * @jsize4] * @lhsm[4 + k * @jsize4]
                rtmp[m + k2 * 5] = rtmp[m + k2 * 5] - @lhsm[0 + k2 * @jsize4] * rtmp[m + k * 5]
            end
            k = @grid_points[2] - 2
            k1 = @grid_points[2] - 1
            m = 3
            fac1 = 1. / @lhsp[2 + k * @jsize4]
            @lhsp[3 + k * @jsize4] = fac1 * @lhsp[3 + k * @jsize4]
            @lhsp[4 + k * @jsize4] = fac1 * @lhsp[4 + k * @jsize4]
            rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
            @lhsp[2 + k1 * @jsize4] = @lhsp[2 + k1 * @jsize4] - @lhsp[1 + k1 * @jsize4] * @lhsp[3 + k * @jsize4]
            @lhsp[3 + k1 * @jsize4] = @lhsp[3 + k1 * @jsize4] - @lhsp[1 + k1 * @jsize4] * @lhsp[4 + k * @jsize4]
            rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhsp[1 + k1 * @jsize4] * rtmp[m + k * 5]
            m = 4
            fac1 = 1.0 / @lhsm[2 + k * @jsize4]
            @lhsm[3 + k * @jsize4] = fac1 * @lhsm[3 + k * @jsize4]
            @lhsm[4 + k * @jsize4] = fac1 * @lhsm[4 + k * @jsize4]
            rtmp[m + k * 5] = fac1 * rtmp[m + k * 5]
            @lhsm[2 + k1 * @jsize4] = @lhsm[2 + k1 * @jsize4] - @lhsm[1 + k1 * @jsize4] * @lhsm[3 + k * @jsize4]
            @lhsm[3 + k1 * @jsize4] = @lhsm[3 + k1 * @jsize4] - @lhsm[1 + k1 * @jsize4] * @lhsm[4 + k * @jsize4]
            rtmp[m + k1 * 5] = rtmp[m + k1 * 5] - @lhsm[1 + k1 * @jsize4] * rtmp[m + k * 5]
            rtmp[3 + k1 * 5] = rtmp[3 + k1 * 5] / @lhsp[2 + k1 * @jsize4]
            rtmp[4 + k1 * 5] = rtmp[4 + k1 * 5] / @lhsm[2 + k1 * @jsize4]
            k = @grid_points[2] - 2
            k1 = @grid_points[2] - 1
            for m in 0..2 do
                rtmp[m + k * 5] = rtmp[m + k * 5] - @lhs[3 + k * @jsize4] * rtmp[m + k1 * 5]
            end
            rtmp[3 + k * 5] = rtmp[3 + k * 5] - @lhsp[3 + k * @jsize4] * rtmp[3 + k1 * 5]
            rtmp[4 + k * 5] = rtmp[4 + k * 5] - @lhsm[3 + k * @jsize4] * rtmp[4 + k1 * 5]
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
            k = @grid_points[2] - 3
            while k >= 0 do
                k1 = k + 1
                k2 = k + 2
                for m in 0..2 do
                    rtmp[m + k * 5] = rtmp[m + k * 5] - @lhs[3 + k * @jsize4] * rtmp[m + k1 * 5] - @lhs[4 + k * @jsize4] * rtmp[m + k2 * 5]
                end
                rtmp[3 + k * 5] = rtmp[3 + k * 5] - @lhsp[3 + k * @jsize4] * rtmp[3 + k1 * 5] - @lhsp[4 + k * @jsize4] * rtmp[3 + k2 * 5]
                rtmp[4 + k * 5] = rtmp[4 + k * 5] - @lhsm[3 + k * @jsize4] * rtmp[4 + k1 * 5] - @lhsm[4 + k * @jsize4] * rtmp[4 + k2 * 5]
              k -= 1
            end
            for k in 0..@nz2 + 1 do
                @rhs[0 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[0 + k * 5]
                @rhs[1 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[1 + k * 5]
                @rhs[2 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[2 + k * 5]
                @rhs[3 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[3 + k * 5]
                @rhs[4 + i * @isize1 + j * @jsize1 + k * @ksize1] = rtmp[4 + k * 5]
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_zsolve)
    else
    end
    if @timeron then
      @timer.start(@@t_tzetar)
    else
    end
    tzetar()
    if @timeron then
      @timer.stop(@@t_tzetar)
    else
    end
  end

  def checkSum(arr)
    csum = 0.0; 
    for k in 0..@grid_points[2] - 1 do
        for j in 0..@grid_points[1] - 1 do
            for i in 0..@grid_points[0] - 1 do
                for m in 0..4 do
                    offset = m + i * @isize1 + j * @jsize1 + k * @ksize1; 
                    csum += (arr[offset] * arr[offset]) / (@grid_points[2] * @grid_points[1] * @grid_points[0] * 5)
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

  def setupThreads(sp)
    @master = sp
    if @num_threads > @problem_size - 2 then
      @num_threads = @problem_size - 2
    else
    end
    interval1 = Array.new(@num_threads, 0.0); 
    interval2 = Array.new(@num_threads, 0.0); 
    set_interval(@problem_size, interval1)
    set_interval(@problem_size - 2, interval2)
    partition1 = MultiDimArray(interval1.length, 2); 
    partition2 = MultiDimArray(interval2.length, 2); 
    set_partition(0, interval1, partition1)
    set_partition(1, interval2, partition2)
    @rhscomputer = Array.new(@num_threads, 0.0)
    @txinverse = Array.new(@num_threads, 0.0)
    @xsolver = Array.new(@num_threads, 0.0)
    @ysolver = Array.new(@num_threads, 0.0)
    @zsolver = Array.new(@num_threads, 0.0)
    @rhsadder = Array.new(@num_threads, 0.0)
    for ii in 0..@num_threads-1 do
      @rhscomputer[ii] = RHSCompute.new(sp, partition1[ii][0], partition1[ii][1], partition2[ii][0], partition2[ii][1])
      @rhscomputer[ii].extend(MonitorMixin)
      @rhscomputer[ii].id = ii
      @rhscomputer[ii].start()

      @xsolver[ii] = XSolver.new(sp, partition2[ii][0], partition2[ii][1])
      @xsolver[ii].extend(MonitorMixin)
      @xsolver[ii].id = ii
      @xsolver[ii].start()

      @txinverse[ii] = TXInverse.new(sp, partition2[ii][0], partition2[ii][1])
      @txinverse[ii].extend(MonitorMixin)
      @txinverse[ii].id = ii
      @txinverse[ii].start()

      @ysolver[ii] = YSolver.new(sp, partition2[ii][0], partition2[ii][1])
      @ysolver[ii].extend(MonitorMixin)
      @ysolver[ii].id = ii
      @ysolver[ii].start()

      @zsolver[ii] = ZSolver.new(sp, partition2[ii][0], partition2[ii][1])
      @zsolver[ii].extend(MonitorMixin)
      @zsolver[ii].id = ii
      @zsolver[ii].start()

      @rhsadder[ii] = RHSAdder.new(sp, partition2[ii][0], partition2[ii][1])
      @rhsadder[ii].extend(MonitorMixin)
      @rhsadder[ii].id = ii
      @rhsadder[ii].start()
    end
  end

# *** public ***

  attr_accessor :bid

  attr_accessor :results

  attr_accessor :serial

end



  def main(argv)
    sp =  nil 
    BMArgs.parseCmdLineArgs(argv, "SP")
    clss = BMArgs.clss; 
    np = BMArgs.num_threads; 
    serial = BMArgs.serial; 
    #begin
    sp = SP.new(clss, np, serial)
    sp.extend(MonitorMixin)
    #rescue OutOfMemoryError => e
    #    BMArgs.outOfMemoryMessage()
    #    System.exit(0)
    #ensure
    #end
    sp.runBenchMark()
  end

main(ARGV)
