# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
#  Authors: R. Van der Wijngaart
#	    T. Harris
#	    M. Yarrow
#  Modified for PBN (Programming Baseline for NPB):
#	    H. Jin
#  Translation to Java and to MultiThreaded Code
#	    M. Frumkin
#	    M. Schultz
#  Translation to Ruby
#           T. Nose
#           H. Tomari
# -------------------------------------------------------------------------

require "mda"
require "mathutil"
require "BTThreads/BTBase"
require "BTThreads/RHSCompute"
require "BTThreads/XSolver"
require "BTThreads/YSolver"
require "BTThreads/ZSolver"
require "BTThreads/RHSAdder"

require "BMInOut/BMResults"
require "BMInOut/BMArgs"
require "Timer"
require "Random"
require "Runnable"
require "monitor"

class BT < BTBase
  def initialize(clss, threads, ser)
    @master =  nil 
    @bid = -1
    @serial = true
    super(clss, threads)
    @serial = ser
    @fjac = Array.new(5 * 5 * (@problem_size + 1), 0.0)
    @njac = Array.new(5 * 5 * (@problem_size + 1), 0.0)
    @lhs = Array.new(5 * 5 * 3 * (@problem_size + 1), 0.0)
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
    set_constants()
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
    @timer.stop(@@t_total)
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
        navg = (@grid_points[0] + @grid_points[1] + @grid_points[2]) / 3.0; 
        mflops = 3478.8 * n3 - 17655.7 * Math_dot_pow(navg, 2) + 28023.7 * navg
        mflops *= niter / (total_time * 1000000.0)
    else
    end
    return mflops
  end

  def adi_serial()
    compute_rhs()
    x_solve()
    y_solve()
    z_solve()
    add()
  end

  def adi()
    if @timeron then
      @timer.start(@@t_rhs)
    else
    end
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
      @timer.start(@@t_xsolve)
    else
    end
     self .synchronize do
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
    end
    if @timeron then
      @timer.stop(@@t_xsolve)
    else
    end
    if @timeron then
      @timer.start(@@t_ysolve)
    else
    end
     self .synchronize do
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
    end
    if @timeron then
      @timer.stop(@@t_ysolve)
    else
    end
    if @timeron then
      @timer.start(@@t_zsolve)
    else
    end
     self .synchronize do
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
    end
    if @timeron then
      @timer.stop(@@t_zsolve)
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

  def printTimers(t_names, trecs, tmax)
    #fmt = DecimalFormat.new("0.000"); 
    
    puts("SECTION  Time           (secs)")
    for i in 1..@@t_last do
      trecs[i] = @timer.readTimer(i)
    end
    if tmax == 0.0 then
      tmax = 1.0
    else
    end
    for i in 1..@@t_last do
        puts(t_names[i] + ":" + (trecs[i]).to_s + ":" + "  (" + (trecs[i] * 100 / tmax).to_s + "%)")
        if i == @@t_rhs then
            t = trecs[@@t_rhsx] + trecs[@@t_rhsy] + trecs[@@t_rhsz]
            print("    --> total ")
            print("sub-rhs ")
            print((t))
            print("  (")
            print((t * 100 / tmax))
            puts("%)")
            t = trecs[@@t_rhs] - trecs[@@t_rhsx] + trecs[@@t_rhsy] + trecs[@@t_rhsz]
            print("    --> total ")
            print("rest-rhs ")
            print((t))
            print("  (")
            print((t * 100 / tmax))
            puts("%)")
        else
          if i == @@t_zsolve then
              t = trecs[@@t_zsolve] - trecs[@@t_rdis1] - trecs[@@t_rdis2]
              print("    --> total ")
              print("sub-zsol ")
              print((t))
              print("  ")
              puts((t * 100 / tmax))
              puts()
          else
            if i == @@t_rdis2 then
                t = trecs[@@t_rdis1] + trecs[@@t_rdis2]
                print("    --> total ")
                print("redist ")
                print((t))
                print("  ")
                puts((t * 100 / tmax))
            else
            end
          end
        end
    end
  end

  def getInputPars()
    niter = 0; 
    #f2 = File.new("inputbt.data"); 
    if File.exist?("inputbt.data") then
        begin
            fis = FileInputStream.new(f2); 
            datafile = DataInputStream.new(fis); 
            puts("Reading from input file inputbt.data")
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
        puts("No input file inputbt.data, Using compiled defaults")
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
    return niter
  end

  def setTimers(t_names)
    #f1 = File.new("timer.flag"); 
    @timeron = false
    if File.exist?("timer.flag") then
        @timeron = true
        t_names[@@t_total] = String.new("total    ")
        t_names[@@t_rhsx] = String.new("rhsx     ")
        t_names[@@t_rhsy] = String.new("rhsy     ")
        t_names[@@t_rhsz] = String.new("rhsz     ")
        t_names[@@t_rhs] = String.new("rhs      ")
        t_names[@@t_xsolve] = String.new("xsolve   ")
        t_names[@@t_ysolve] = String.new("ysolve   ")
        t_names[@@t_zsolve] = String.new("zsolve   ")
        t_names[@@t_rdis1] = String.new("redist1  ")
        t_names[@@t_rdis2] = String.new("redist2  ")
        t_names[@@t_add] = String.new("add      ")
    else
    end
  end

  def rhs_norm(rms, rmsoffst)
    
    
    for m in 0..rms.length-1 do
      rms[m + rmsoffst] = 0.0
    end
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                for m in 0..rms.length-1 do
                    add = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2]
                    rms[m] += add * add
                end
            end
        end
    end
    for m in 0..rms.length-1 do
        for d in 0..2 do
            rms[m] /= @grid_points[d] - 2
        end
        rms[m] = Math.sqrt(rms[m + rmsoffst])
    end
  end

  def error_norm(rms, rmsoffst)
    
    u_exact = Array.new(5, 0.0); 
    for m in 0..rms.length-1 do
      rms[m + rmsoffst] = 0.0
    end
    for k in 0..@grid_points[2] - 1 do
        zeta = k * @@dnzm1
        for j in 0..@grid_points[1] - 1 do
            eta = j * @@dnym1
            for i in 0..@grid_points[0] - 1 do
                xi = i * @@dnxm1
                exact_solution(xi, eta, zeta, u_exact, 0)
                for m in 0..rms.length-1 do
                    add = @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - u_exact[m]
                    rms[m] += add * add
                end
            end
        end
    end
    for m in 0..rms.length-1 do
        for d in 0..2 do
            rms[m] /= @grid_points[d] - 2
        end
        rms[m] = Math.sqrt(rms[m])
    end
  end

  def verify(no_time_steps)
    xcrref = Array.new(5, 0.0); xceref = Array.new(5, 0.0); xcrdif = Array.new(5, 0.0); xcedif = Array.new(5, 0.0); xce = Array.new(5, 0.0); xcr = Array.new(5, 0.0); dtref = 0; 
    
    verified = -1; 
    clss = 'U'; 
    error_norm(xce, 0)
    compute_rhs()
    rhs_norm(xcr, 0)
    for m in 0..xcr.length-1 do
      xcr[m] = xcr[m] / @@dt
    end
    for m in 1..xcrref.length-1 do
        xcrref[m] = 1.0
        xceref[m] = 1.0
    end
    if @grid_points[0] == 12 and @grid_points[1] == 12 and @grid_points[2] == 12 and no_time_steps == 60 then
        clss = 'S'
        dtref = 0.01
        xcrref[0] = 1.7034283709541311E-01
        xcrref[1] = 1.2975252070034097E-02
        xcrref[2] = 3.2527926989486055E-02
        xcrref[3] = 2.6436421275166801E-02
        xcrref[4] = 1.9211784131744430E-01
        xceref[0] = 4.9976913345811579E-04
        xceref[1] = 4.5195666782961927E-05
        xceref[2] = 7.3973765172921357E-05
        xceref[3] = 7.3821238632439731E-05
        xceref[4] = 8.9269630987491446E-04
    else
      if (@grid_points[0] == 24) and (@grid_points[1] == 24) and (@grid_points[2] == 24) and (no_time_steps == 200) then
          clss = 'W'
          dtref = 0.0008
          xcrref[0] = 0.1125590409344E+03
          xcrref[1] = 0.1180007595731E+02
          xcrref[2] = 0.2710329767846E+02
          xcrref[3] = 0.2469174937669E+02
          xcrref[4] = 0.2638427874317E+03
          xceref[0] = 0.4419655736008E+01
          xceref[1] = 0.4638531260002
          xceref[2] = 0.1011551749967E+01
          xceref[3] = 0.9235878729944
          xceref[4] = 0.1018045837718E+02
      else
        if (@grid_points[0] == 64) and (@grid_points[1] == 64) and (@grid_points[2] == 64) and (no_time_steps == 200) then
            clss = 'A'
            dtref = 0.0008
            xcrref[0] = 1.0806346714637264E+02
            xcrref[1] = 1.1319730901220813E+01
            xcrref[2] = 2.5974354511582465E+01
            xcrref[3] = 2.3665622544678910E+01
            xcrref[4] = 2.5278963211748344E+02
            xceref[0] = 4.2348416040525025
            xceref[1] = 4.4390282496995698E-01
            xceref[2] = 9.6692480136345650E-01
            xceref[3] = 8.8302063039765474E-01
            xceref[4] = 9.7379901770829278
        else
          if (@grid_points[0] == 102) and (@grid_points[1] == 102) and (@grid_points[2] == 102) and (no_time_steps == 200) then
              clss = 'B'
              dtref = 0.0003
              xcrref[0] = 1.4233597229287254E+03
              xcrref[1] = 9.9330522590150238E+01
              xcrref[2] = 3.5646025644535285E+02
              xcrref[3] = 3.2485447959084092E+02
              xcrref[4] = 3.2707541254659363E+03
              xceref[0] = 5.2969847140936856E+01
              xceref[1] = 4.4632896115670668
              xceref[2] = 1.3122573342210174E+01
              xceref[3] = 1.2006925323559144E+01
              xceref[4] = 1.2459576151035986E+02
          else
            if (@grid_points[0] == 162) and (@grid_points[1] == 162) and (@grid_points[2] == 162) and (no_time_steps == 200) then
                clss = 'C'
                dtref = 0.0001
                xcrref[0] = 0.62398116551764615E+04
                xcrref[1] = 0.50793239190423964E+03
                xcrref[2] = 0.15423530093013596E+04
                xcrref[3] = 0.13302387929291190E+04
                xcrref[4] = 0.11604087428436455E+05
                xceref[0] = 0.16462008369091265E+03
                xceref[1] = 0.11497107903824313E+02
                xceref[2] = 0.41207446207461508E+02
                xceref[3] = 0.37087651059694167E+02
                xceref[4] = 0.36211053051841265E+03
            else
            end
          end
        end
      end
    end
    for m in 0..xcr.length-1 do
        xcrdif[m] = Math_dot_abs((xcr[m] - xcrref[m]) / xcrref[m])
        xcedif[m] = Math_dot_abs((xce[m] - xceref[m]) / xceref[m])
    end
    epsilon = 1.0 * Math_dot_pow(0.1, 8); 
    if clss != 'U' then
        puts("Verification being performed for class " + clss)
        puts("accuracy setting for epsilon = " + epsilon.to_s)
        if Math_dot_abs(@@dt - dtref) <= epsilon then
            verified = 1
        else
            verified = 0
            clss = 'U'
            puts("DT does not match the reference value of " + dtref.to_s)
        end
    else
        puts("Unknown class")
    end
    if clss != 'U' then
      puts("Comparison of RMS-norms of residual")
    else
      puts("RMS-norms of residual")
    end
    verified = BMResults.printComparisonStatusArr(clss, verified, epsilon, xcr, xcrref, xcrdif)
    if clss != 'U' then
        puts("Comparison of RMS-norms of solution error")
    else
        puts("RMS-norms of solution error")
    end
    verified = BMResults.printComparisonStatusArr(clss, verified, epsilon, xce, xceref, xcedif)
    BMResults.printVerificationStatus(clss, verified, @@BMName)
    return verified
  end

  def add()
    
    if @timeron then
      @timer.start(@@t_add)
    else
    end
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                for m in 0..4 do
                    @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] += @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2]
                end
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_add)
    else
    end
  end

  def exact_rhs()
    dtemp = Array.new(5, 0.0); 
    
    for k in 0..@grid_points[2] - 1 do
        for j in 0..@grid_points[1] - 1 do
            for i in 0..@grid_points[0] - 1 do
                for m in 0..4 do
                    @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = 0.0
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
                @forcing[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[0 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tx2 * (@ue[ip1 + 1 * @jsize3] - @ue[im1 + 1 * @jsize3]) + @@dx1tx1 * (@ue[ip1 + 0 * @jsize3] - 2.0 * @ue[i + 0 * @jsize3] + @ue[im1 + 0 * @jsize3])
                @forcing[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[1 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tx2 * ((@ue[ip1 + 1 * @jsize3] * @buf[ip1 + 1 * @jsize3] + @@c2 * (@ue[ip1 + 4 * @jsize3] - @q[ip1])) - (@ue[im1 + 1 * @jsize3] * @buf[im1 + 1 * @jsize3] + @@c2 * (@ue[im1 + 4 * @jsize3] - @q[im1]))) + @@xxcon1 * (@buf[ip1 + 1 * @jsize3] - 2.0 * @buf[i + 1 * @jsize3] + @buf[im1 + 1 * @jsize3]) + @@dx2tx1 * (@ue[ip1 + 1 * @jsize3] - 2.0 * @ue[i + 1 * @jsize3] + @ue[im1 + 1 * @jsize3])
                @forcing[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[2 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tx2 * (@ue[ip1 + 2 * @jsize3] * @buf[ip1 + 1 * @jsize3] - @ue[im1 + 2 * @jsize3] * @buf[im1 + 1 * @jsize3]) + @@xxcon2 * (@buf[ip1 + 2 * @jsize3] - 2.0 * @buf[i + 2 * @jsize3] + @buf[im1 + 2 * @jsize3]) + @@dx3tx1 * (@ue[ip1 + 2 * @jsize3] - 2.0 * @ue[i + 2 * @jsize3] + @ue[im1 + 2 * @jsize3])
                @forcing[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[3 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tx2 * (@ue[ip1 + 3 * @jsize3] * @buf[ip1 + 1 * @jsize3] - @ue[im1 + 3 * @jsize3] * @buf[im1 + 1 * @jsize3]) + @@xxcon2 * (@buf[ip1 + 3 * @jsize3] - 2.0 * @buf[i + 3 * @jsize3] + @buf[im1 + 3 * @jsize3]) + @@dx4tx1 * (@ue[ip1 + 3 * @jsize3] - 2.0 * @ue[i + 3 * @jsize3] + @ue[im1 + 3 * @jsize3])
                @forcing[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[4 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tx2 * (@buf[ip1 + 1 * @jsize3] * (@@c1 * @ue[ip1 + 4 * @jsize3] - @@c2 * @q[ip1]) - @buf[im1 + 1 * @jsize3] * (@@c1 * @ue[im1 + 4 * @jsize3] - @@c2 * @q[im1])) + 0.5 * @@xxcon3 * (@buf[ip1 + 0 * @jsize3] - 2.0 * @buf[i + 0 * @jsize3] + @buf[im1 + 0 * @jsize3]) + @@xxcon4 * (@cuf[ip1] - 2.0 * @cuf[i] + @cuf[im1]) + @@xxcon5 * (@buf[ip1 + 4 * @jsize3] - 2.0 * @buf[i + 4 * @jsize3] + @buf[im1 + 4 * @jsize3]) + @@dx5tx1 * (@ue[ip1 + 4 * @jsize3] - 2.0 * @ue[i + 4 * @jsize3] + @ue[im1 + 4 * @jsize3])
            end
            for m in 0..4 do
                i = 1
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @ue[i + m * @jsize3] - 4.0 * @ue[i + 1 + m * @jsize3] + @ue[i + 2 + m * @jsize3])
                i = 2
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @ue[i - 1 + m * @jsize3] + 6.0 * @ue[i + m * @jsize3] - 4.0 * @ue[i + 1 + m * @jsize3] + @ue[i + 2 + m * @jsize3])
            end
            for m in 0..4 do
                for i in 3..@grid_points[0] - 4 do
                    @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[i - 2 + m * @jsize3] - 4.0 * @ue[i - 1 + m * @jsize3] + 6.0 * @ue[i + m * @jsize3] - 4.0 * @ue[i + 1 + m * @jsize3] + @ue[i + 2 + m * @jsize3])
                end
            end
            for m in 0..4 do
                i = @grid_points[0] - 3
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[i - 2 + m * @jsize3] - 4.0 * @ue[i - 1 + m * @jsize3] + 6.0 * @ue[i + m * @jsize3] - 4.0 * @ue[i + 1 + m * @jsize3])
                i = @grid_points[0] - 2
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[i - 2 + m * @jsize3] - 4.0 * @ue[i - 1 + m * @jsize3] + 5.0 * @ue[i + m * @jsize3])
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
                @forcing[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[0 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@ty2 * (@ue[jp1 + 2 * @jsize3] - @ue[jm1 + 2 * @jsize3]) + @@dy1ty1 * (@ue[jp1 + 0 * @jsize3] - 2.0 * @ue[j + 0 * @jsize3] + @ue[jm1 + 0 * @jsize3])
                @forcing[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[1 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@ty2 * (@ue[jp1 + 1 * @jsize3] * @buf[jp1 + 2 * @jsize3] - @ue[jm1 + 1 * @jsize3] * @buf[jm1 + 2 * @jsize3]) + @@yycon2 * (@buf[jp1 + 1 * @jsize3] - 2.0 * @buf[j + 1 * @jsize3] + @buf[jm1 + 1 * @jsize3]) + @@dy2ty1 * (@ue[jp1 + 1 * @jsize3] - 2.0 * @ue[j + 1 * @jsize3] + @ue[jm1 + 1 * @jsize3])
                @forcing[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[2 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@ty2 * ((@ue[jp1 + 2 * @jsize3] * @buf[jp1 + 2 * @jsize3] + @@c2 * (@ue[jp1 + 4 * @jsize3] - @q[jp1])) - (@ue[jm1 + 2 * @jsize3] * @buf[jm1 + 2 * @jsize3] + @@c2 * (@ue[jm1 + 4 * @jsize3] - @q[jm1]))) + @@yycon1 * (@buf[jp1 + 2 * @jsize3] - 2.0 * @buf[j + 2 * @jsize3] + @buf[jm1 + 2 * @jsize3]) + @@dy3ty1 * (@ue[jp1 + 2 * @jsize3] - 2.0 * @ue[j + 2 * @jsize3] + @ue[jm1 + 2 * @jsize3])
                @forcing[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[3 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@ty2 * (@ue[jp1 + 3 * @jsize3] * @buf[jp1 + 2 * @jsize3] - @ue[jm1 + 3 * @jsize3] * @buf[jm1 + 2 * @jsize3]) + @@yycon2 * (@buf[jp1 + 3 * @jsize3] - 2.0 * @buf[j + 3 * @jsize3] + @buf[jm1 + 3 * @jsize3]) + @@dy4ty1 * (@ue[jp1 + 3 * @jsize3] - 2.0 * @ue[j + 3 * @jsize3] + @ue[jm1 + 3 * @jsize3])
                @forcing[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[4 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@ty2 * (@buf[jp1 + 2 * @jsize3] * (@@c1 * @ue[jp1 + 4 * @jsize3] - @@c2 * @q[jp1]) - @buf[jm1 + 2 * @jsize3] * (@@c1 * @ue[jm1 + 4 * @jsize3] - @@c2 * @q[jm1])) + 0.5 * @@yycon3 * (@buf[jp1 + 0 * @jsize3] - 2.0 * @buf[j + 0 * @jsize3] + @buf[jm1 + 0 * @jsize3]) + @@yycon4 * (@cuf[jp1] - 2.0 * @cuf[j] + @cuf[jm1]) + @@yycon5 * (@buf[jp1 + 4 * @jsize3] - 2.0 * @buf[j + 4 * @jsize3] + @buf[jm1 + 4 * @jsize3]) + @@dy5ty1 * (@ue[jp1 + 4 * @jsize3] - 2.0 * @ue[j + 4 * @jsize3] + @ue[jm1 + 4 * @jsize3])
            end
            for m in 0..4 do
                j = 1
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @ue[j + m * @jsize3] - 4.0 * @ue[j + 1 + m * @jsize3] + @ue[j + 2 + m * @jsize3])
                j = 2
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @ue[j - 1 + m * @jsize3] + 6.0 * @ue[j + m * @jsize3] - 4.0 * @ue[j + 1 + m * @jsize3] + @ue[j + 2 + m * @jsize3])
            end
            for m in 0..4 do
                for j in 3..@grid_points[1] - 4 do
                    @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[j - 2 + m * @jsize3] - 4.0 * @ue[j - 1 + m * @jsize3] + 6.0 * @ue[j + m * @jsize3] - 4.0 * @ue[j + 1 + m * @jsize3] + @ue[j + 2 + m * @jsize3])
                end
            end
            for m in 0..4 do
                j = @grid_points[1] - 3
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[j - 2 + m * @jsize3] - 4.0 * @ue[j - 1 + m * @jsize3] + 6.0 * @ue[j + m * @jsize3] - 4.0 * @ue[j + 1 + m * @jsize3])
                j = @grid_points[1] - 2
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[j - 2 + m * @jsize3] - 4.0 * @ue[j - 1 + m * @jsize3] + 5.0 * @ue[j + m * @jsize3])
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
                @forcing[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[0 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tz2 * (@ue[kp1 + 3 * @jsize3] - @ue[km1 + 3 * @jsize3]) + @@dz1tz1 * (@ue[kp1 + 0 * @jsize3] - 2.0 * @ue[k + 0 * @jsize3] + @ue[km1 + 0 * @jsize3])
                @forcing[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[1 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tz2 * (@ue[kp1 + 1 * @jsize3] * @buf[kp1 + 3 * @jsize3] - @ue[km1 + 1 * @jsize3] * @buf[km1 + 3 * @jsize3]) + @@zzcon2 * (@buf[kp1 + 1 * @jsize3] - 2.0 * @buf[k + 1 * @jsize3] + @buf[km1 + 1 * @jsize3]) + @@dz2tz1 * (@ue[kp1 + 1 * @jsize3] - 2.0 * @ue[k + 1 * @jsize3] + @ue[km1 + 1 * @jsize3])
                @forcing[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[2 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tz2 * (@ue[kp1 + 2 * @jsize3] * @buf[kp1 + 3 * @jsize3] - @ue[km1 + 2 * @jsize3] * @buf[km1 + 3 * @jsize3]) + @@zzcon2 * (@buf[kp1 + 2 * @jsize3] - 2.0 * @buf[k + 2 * @jsize3] + @buf[km1 + 2 * @jsize3]) + @@dz3tz1 * (@ue[kp1 + 2 * @jsize3] - 2.0 * @ue[k + 2 * @jsize3] + @ue[km1 + 2 * @jsize3])
                @forcing[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[3 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tz2 * ((@ue[kp1 + 3 * @jsize3] * @buf[kp1 + 3 * @jsize3] + @@c2 * (@ue[kp1 + 4 * @jsize3] - @q[kp1])) - (@ue[km1 + 3 * @jsize3] * @buf[km1 + 3 * @jsize3] + @@c2 * (@ue[km1 + 4 * @jsize3] - @q[km1]))) + @@zzcon1 * (@buf[kp1 + 3 * @jsize3] - 2.0 * @buf[k + 3 * @jsize3] + @buf[km1 + 3 * @jsize3]) + @@dz4tz1 * (@ue[kp1 + 3 * @jsize3] - 2.0 * @ue[k + 3 * @jsize3] + @ue[km1 + 3 * @jsize3])
                @forcing[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[4 + i * @isize2 + j * @jsize2 + k * @ksize2] - @@tz2 * (@buf[kp1 + 3 * @jsize3] * (@@c1 * @ue[kp1 + 4 * @jsize3] - @@c2 * @q[kp1]) - @buf[km1 + 3 * @jsize3] * (@@c1 * @ue[km1 + 4 * @jsize3] - @@c2 * @q[km1])) + 0.5 * @@zzcon3 * (@buf[kp1 + 0 * @jsize3] - 2.0 * @buf[k + 0 * @jsize3] + @buf[km1 + 0 * @jsize3]) + @@zzcon4 * (@cuf[kp1] - 2.0 * @cuf[k] + @cuf[km1]) + @@zzcon5 * (@buf[kp1 + 4 * @jsize3] - 2.0 * @buf[k + 4 * @jsize3] + @buf[km1 + 4 * @jsize3]) + @@dz5tz1 * (@ue[kp1 + 4 * @jsize3] - 2.0 * @ue[k + 4 * @jsize3] + @ue[km1 + 4 * @jsize3])
            end
            for m in 0..4 do
                k = 1
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @ue[k + m * @jsize3] - 4.0 * @ue[k + 1 + m * @jsize3] + @ue[k + 2 + m * @jsize3])
                k = 2
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @ue[k - 1 + m * @jsize3] + 6.0 * @ue[k + m * @jsize3] - 4.0 * @ue[k + 1 + m * @jsize3] + @ue[k + 2 + m * @jsize3])
            end
            for m in 0..4 do
                for k in 3..@grid_points[2] - 4 do
                    @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[k - 2 + m * @jsize3] - 4.0 * @ue[k - 1 + m * @jsize3] + 6.0 * @ue[k + m * @jsize3] - 4.0 * @ue[k + 1 + m * @jsize3] + @ue[k + 2 + m * @jsize3])
                end
            end
            for m in 0..4 do
                k = @grid_points[2] - 3
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[k - 2 + m * @jsize3] - 4.0 * @ue[k - 1 + m * @jsize3] + 6.0 * @ue[k + m * @jsize3] - 4.0 * @ue[k + 1 + m * @jsize3])
                k = @grid_points[2] - 2
                @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@ue[k - 2 + m * @jsize3] - 4.0 * @ue[k - 1 + m * @jsize3] + 5.0 * @ue[k + m * @jsize3])
            end
        end
    end
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                for m in 0..4 do
                    @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2] = -1.0 * @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2]
                end
            end
        end
    end
  end

  def x_solve()
    
    if @timeron then
      @timer.start(@@t_xsolve)
    else
    end
    isize = @grid_points[0] - 1
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 0..isize do
                @tmp1 = @rho_i[i + j * @jsize1 + k * @ksize1]
                @tmp2 = @tmp1 * @tmp1
                @tmp3 = @tmp1 * @tmp2
                @fjac[0 + 0 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[0 + 1 * @@isize4 + i * @@jsize4] = 1.0
                @fjac[0 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[0 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[0 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[1 + 0 * @@isize4 + i * @@jsize4] = -(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@c2 * @qs[i + j * @jsize1 + k * @ksize1]
                @fjac[1 + 1 * @@isize4 + i * @@jsize4] = (2.0 - @@c2) * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] / @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2])
                @fjac[1 + 2 * @@isize4 + i * @@jsize4] = -@@c2 * (@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1)
                @fjac[1 + 3 * @@isize4 + i * @@jsize4] = -@@c2 * (@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1)
                @fjac[1 + 4 * @@isize4 + i * @@jsize4] = @@c2
                @fjac[2 + 0 * @@isize4 + i * @@jsize4] = -(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[2 + 1 * @@isize4 + i * @@jsize4] = @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 2 * @@isize4 + i * @@jsize4] = @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[2 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[3 + 0 * @@isize4 + i * @@jsize4] = -(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[3 + 1 * @@isize4 + i * @@jsize4] = @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[3 + 3 * @@isize4 + i * @@jsize4] = @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @fjac[4 + 0 * @@isize4 + i * @@jsize4] = (@@c2 * 2.0 * @square[i + j * @jsize1 + k * @ksize1] - @@c1 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2]) * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2)
                @fjac[4 + 1 * @@isize4 + i * @@jsize4] = @@c1 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1 - @@c2 * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2 + @qs[i + j * @jsize1 + k * @ksize1])
                @fjac[4 + 2 * @@isize4 + i * @@jsize4] = -@@c2 * (@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[4 + 3 * @@isize4 + i * @@jsize4] = -@@c2 * (@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[4 + 4 * @@isize4 + i * @@jsize4] = @@c1 * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1)
                @njac[0 + 0 * @@isize4 + i * @@jsize4] = 0.0
                @njac[0 + 1 * @@isize4 + i * @@jsize4] = 0.0
                @njac[0 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @njac[0 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @njac[0 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @njac[1 + 0 * @@isize4 + i * @@jsize4] = -@@con43 * @@c3c4 * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[1 + 1 * @@isize4 + i * @@jsize4] = @@con43 * @@c3c4 * @tmp1
                @njac[1 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @njac[1 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @njac[1 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @njac[2 + 0 * @@isize4 + i * @@jsize4] = -@@c3c4 * @tmp2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[2 + 1 * @@isize4 + i * @@jsize4] = 0.0
                @njac[2 + 2 * @@isize4 + i * @@jsize4] = @@c3c4 * @tmp1
                @njac[2 + 3 * @@isize4 + i * @@jsize4] = 0.0
                @njac[2 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @njac[3 + 0 * @@isize4 + i * @@jsize4] = -@@c3c4 * @tmp2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[3 + 1 * @@isize4 + i * @@jsize4] = 0.0
                @njac[3 + 2 * @@isize4 + i * @@jsize4] = 0.0
                @njac[3 + 3 * @@isize4 + i * @@jsize4] = @@c3c4 * @tmp1
                @njac[3 + 4 * @@isize4 + i * @@jsize4] = 0.0
                @njac[4 + 0 * @@isize4 + i * @@jsize4] = -(@@con43 * @@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - (@@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - (@@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - @@c1345 * @tmp2 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 1 * @@isize4 + i * @@jsize4] = (@@con43 * @@c3c4 - @@c1345) * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 2 * @@isize4 + i * @@jsize4] = (@@c3c4 - @@c1345) * @tmp2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 3 * @@isize4 + i * @@jsize4] = (@@c3c4 - @@c1345) * @tmp2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 4 * @@isize4 + i * @@jsize4] = (@@c1345) * @tmp1
            end
            lhsinit(@lhs, isize)
            for i in 1..isize - 1 do
                @tmp1 = @@dt * @@tx1
                @tmp2 = @@dt * @@tx2
                @lhs[0 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx1
                @lhs[0 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 1 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 2 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 3 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[0 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[0 + 4 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 0 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx2
                @lhs[1 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 2 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 3 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[1 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[1 + 4 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 0 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 1 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx3
                @lhs[2 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 3 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[2 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[2 + 4 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 0 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 1 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 2 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx4
                @lhs[3 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[3 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[3 + 4 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 0 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 0 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 1 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 1 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 2 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 2 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 3 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 3 * @@isize4 + (i - 1) * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4] = -@tmp2 * @fjac[4 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @njac[4 + 4 * @@isize4 + (i - 1) * @@jsize4] - @tmp1 * @@dx5
                @lhs[0 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[0 + 0 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx1
                @lhs[0 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 1 * @@isize4 + i * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 2 * @@isize4 + i * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 3 * @@isize4 + i * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 4 * @@isize4 + i * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 0 * @@isize4 + i * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[1 + 1 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx2
                @lhs[1 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 2 * @@isize4 + i * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 3 * @@isize4 + i * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 4 * @@isize4 + i * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 0 * @@isize4 + i * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 1 * @@isize4 + i * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[2 + 2 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx3
                @lhs[2 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 3 * @@isize4 + i * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 4 * @@isize4 + i * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 0 * @@isize4 + i * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 1 * @@isize4 + i * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 2 * @@isize4 + i * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[3 + 3 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx4
                @lhs[3 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 4 * @@isize4 + i * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 0 * @@isize4 + i * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 1 * @@isize4 + i * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 2 * @@isize4 + i * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 3 * @@isize4 + i * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[4 + 4 * @@isize4 + i * @@jsize4] + @tmp1 * 2.0 * @@dx5
                @lhs[0 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx1
                @lhs[0 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 1 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 2 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 3 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[0 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[0 + 4 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 0 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx2
                @lhs[1 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 2 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 3 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[1 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[1 + 4 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 0 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 1 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx3
                @lhs[2 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 3 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[2 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[2 + 4 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 0 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 1 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 2 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx4
                @lhs[3 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[3 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[3 + 4 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 0 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 0 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 1 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 1 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 2 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 2 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 3 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 3 * @@isize4 + (i + 1) * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] = @tmp2 * @fjac[4 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @njac[4 + 4 * @@isize4 + (i + 1) * @@jsize4] - @tmp1 * @@dx5
            end
            binvcrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + 0 * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + 0 * @@ksize4, @rhs, 0 + 0 * @isize2 + j * @jsize2 + k * @ksize2)
            for i in 1..isize - 1 do
                matvec_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4, @rhs, 0 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2, @rhs, 0 + i * @isize2 + j * @jsize2 + k * @ksize2)
                matmul_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + i * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + (i - 1) * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4)
                binvcrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + i * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + i * @@ksize4, @rhs, 0 + i * @isize2 + j * @jsize2 + k * @ksize2)
            end
            matvec_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + isize * @@ksize4, @rhs, 0 + (isize - 1) * @isize2 + j * @jsize2 + k * @ksize2, @rhs, 0 + isize * @isize2 + j * @jsize2 + k * @ksize2)
            matmul_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + isize * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + (isize - 1) * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + isize * @@ksize4)
            binvrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + isize * @@ksize4, @rhs, 0 + isize * @isize2 + j * @jsize2 + k * @ksize2)
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
            i = isize - 1
            while i >= 0 do
                for m in 0..@@BLOCK_SIZE - 1 do
                    for n in 0..@@BLOCK_SIZE - 1 do
                        @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @lhs[m + n * @@isize4 + @@cc * @@jsize4 + i * @@ksize4] * @rhs[n + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2]
                    end
                end
              i -= 1
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_xsolve)
    else
    end
  end

  def compute_rhs()
    
    
    if @timeron then
      @timer.start(@@t_rhs)
    else
    end
    for k in 0..@grid_points[2] - 1 do
        for j in 0..@grid_points[1] - 1 do
            for i in 0..@grid_points[0] - 1 do
                rho_inv = 1.0 / @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @rho_i[i + j * @jsize1 + k * @ksize1] = rho_inv
                @us[i + j * @jsize1 + k * @ksize1] = @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * rho_inv
                @vs[i + j * @jsize1 + k * @ksize1] = @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * rho_inv
                @ws[i + j * @jsize1 + k * @ksize1] = @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * rho_inv
                @square[i + j * @jsize1 + k * @ksize1] = 0.5 * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * rho_inv
                @qs[i + j * @jsize1 + k * @ksize1] = @square[i + j * @jsize1 + k * @ksize1] * rho_inv
            end
        end
    end
    for k in 0..@grid_points[2] - 1 do
        for j in 0..@grid_points[1] - 1 do
            for i in 0..@grid_points[0] - 1 do
                for m in 0..4 do
                    @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @forcing[m + i * @isize2 + j * @jsize2 + k * @ksize2]
                end
            end
        end
    end
    if @timeron then
      @timer.start(@@t_rhsx)
    else
    end
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                uijk = @us[i + j * @jsize1 + k * @ksize1]
                up1 = @us[(i + 1) + j * @jsize1 + k * @ksize1]
                um1 = @us[(i - 1) + j * @jsize1 + k * @ksize1]
                @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx1tx1 * (@u[0 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[0 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) - @@tx2 * (@u[1 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - @u[1 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2])
                @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx2tx1 * (@u[1 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[1 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) + @@xxcon2 * @@con43 * (up1 - 2.0 * uijk + um1) - @@tx2 * (@u[1 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] * up1 - @u[1 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] * um1 + (@u[4 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - @square[(i + 1) + j * @jsize1 + k * @ksize1] - @u[4 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + @square[(i - 1) + j * @jsize1 + k * @ksize1]) * @@c2)
                @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx3tx1 * (@u[2 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[2 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) + @@xxcon2 * (@vs[(i + 1) + j * @jsize1 + k * @ksize1] - 2.0 * @vs[i + j * @jsize1 + k * @ksize1] + @vs[(i - 1) + j * @jsize1 + k * @ksize1]) - @@tx2 * (@u[2 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] * up1 - @u[2 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] * um1)
                @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx4tx1 * (@u[3 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[3 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) + @@xxcon2 * (@ws[(i + 1) + j * @jsize1 + k * @ksize1] - 2.0 * @ws[i + j * @jsize1 + k * @ksize1] + @ws[(i - 1) + j * @jsize1 + k * @ksize1]) - @@tx2 * (@u[3 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] * up1 - @u[3 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] * um1)
                @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dx5tx1 * (@u[4 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[4 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2]) + @@xxcon3 * (@qs[(i + 1) + j * @jsize1 + k * @ksize1] - 2.0 * @qs[i + j * @jsize1 + k * @ksize1] + @qs[(i - 1) + j * @jsize1 + k * @ksize1]) + @@xxcon4 * (up1 * up1 - 2.0 * uijk * uijk + um1 * um1) + @@xxcon5 * (@u[4 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[(i + 1) + j * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[i + j * @jsize1 + k * @ksize1] + @u[4 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[(i - 1) + j * @jsize1 + k * @ksize1]) - @@tx2 * ((@@c1 * @u[4 + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] - @@c2 * @square[(i + 1) + j * @jsize1 + k * @ksize1]) * up1 - (@@c1 * @u[4 + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] - @@c2 * @square[(i - 1) + j * @jsize1 + k * @ksize1]) * um1)
            end
        end
        for j in 1..@grid_points[1] - 2 do
            i = 1
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] + @u[m + (i + 2) * @isize2 + j * @jsize2 + k * @ksize2])
            end
            i = 2
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @u[m + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] + @u[m + (i + 2) * @isize2 + j * @jsize2 + k * @ksize2])
            end
            for m in 0..4 do
                for i in 3..@grid_points[0] - 4 do
                    @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + (i - 2) * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2] + @u[m + (i + 2) * @isize2 + j * @jsize2 + k * @ksize2])
                end
            end
            i = @grid_points[0] - 3
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + (i - 2) * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i + 1) * @isize2 + j * @jsize2 + k * @ksize2])
            end
            i = @grid_points[0] - 2
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + (i - 2) * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + (i - 1) * @isize2 + j * @jsize2 + k * @ksize2] + 5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2])
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
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                vijk = @vs[i + j * @jsize1 + k * @ksize1]
                vp1 = @vs[i + (j + 1) * @jsize1 + k * @ksize1]
                vm1 = @vs[i + (j - 1) * @jsize1 + k * @ksize1]
                @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy1ty1 * (@u[0 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[0 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) - @@ty2 * (@u[2 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - @u[2 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2])
                @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy2ty1 * (@u[1 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[1 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon2 * (@us[i + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @us[i + j * @jsize1 + k * @ksize1] + @us[i + (j - 1) * @jsize1 + k * @ksize1]) - @@ty2 * (@u[1 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] * vp1 - @u[1 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] * vm1)
                @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy3ty1 * (@u[2 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[2 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon2 * @@con43 * (vp1 - 2.0 * vijk + vm1) - @@ty2 * (@u[2 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] * vp1 - @u[2 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] * vm1 + (@u[4 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - @square[i + (j + 1) * @jsize1 + k * @ksize1] - @u[4 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + @square[i + (j - 1) * @jsize1 + k * @ksize1]) * @@c2)
                @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy4ty1 * (@u[3 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[3 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon2 * (@ws[i + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @ws[i + j * @jsize1 + k * @ksize1] + @ws[i + (j - 1) * @jsize1 + k * @ksize1]) - @@ty2 * (@u[3 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] * vp1 - @u[3 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] * vm1)
                @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dy5ty1 * (@u[4 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[4 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2]) + @@yycon3 * (@qs[i + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @qs[i + j * @jsize1 + k * @ksize1] + @qs[i + (j - 1) * @jsize1 + k * @ksize1]) + @@yycon4 * (vp1 * vp1 - 2.0 * vijk * vijk + vm1 * vm1) + @@yycon5 * (@u[4 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] * @rho_i[i + (j + 1) * @jsize1 + k * @ksize1] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[i + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] * @rho_i[i + (j - 1) * @jsize1 + k * @ksize1]) - @@ty2 * ((@@c1 * @u[4 + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] - @@c2 * @square[i + (j + 1) * @jsize1 + k * @ksize1]) * vp1 - (@@c1 * @u[4 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] - @@c2 * @square[i + (j - 1) * @jsize1 + k * @ksize1]) * vm1)
            end
        end
        for i in 1..@grid_points[0] - 2 do
            j = 1
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] + @u[m + i * @isize2 + (j + 2) * @jsize2 + k * @ksize2])
            end
            j = 2
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @u[m + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] + @u[m + i * @isize2 + (j + 2) * @jsize2 + k * @ksize2])
            end
        end
        for j in 3..@grid_points[1] - 4 do
            for i in 1..@grid_points[0] - 2 do
                for m in 0..4 do
                    @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + (j - 2) * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2] + @u[m + i * @isize2 + (j + 2) * @jsize2 + k * @ksize2])
                end
            end
        end
        for i in 1..@grid_points[0] - 2 do
            j = @grid_points[1] - 3
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + (j - 2) * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2])
            end
            j = @grid_points[1] - 2
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + (j - 2) * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2] + 5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2])
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
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                wijk = @ws[i + j * @jsize1 + k * @ksize1]
                wp1 = @ws[i + j * @jsize1 + (k + 1) * @ksize1]
                wm1 = @ws[i + j * @jsize1 + (k - 1) * @ksize1]
                @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz1tz1 * (@u[0 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[0 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) - @@tz2 * (@u[3 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - @u[3 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2])
                @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz2tz1 * (@u[1 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[1 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon2 * (@us[i + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @us[i + j * @jsize1 + k * @ksize1] + @us[i + j * @jsize1 + (k - 1) * @ksize1]) - @@tz2 * (@u[1 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] * wp1 - @u[1 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] * wm1)
                @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz3tz1 * (@u[2 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[2 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon2 * (@vs[i + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @vs[i + j * @jsize1 + k * @ksize1] + @vs[i + j * @jsize1 + (k - 1) * @ksize1]) - @@tz2 * (@u[2 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] * wp1 - @u[2 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] * wm1)
                @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz4tz1 * (@u[3 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[3 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon2 * @@con43 * (wp1 - 2.0 * wijk + wm1) - @@tz2 * (@u[3 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] * wp1 - @u[3 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] * wm1 + (@u[4 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - @square[i + j * @jsize1 + (k + 1) * @ksize1] - @u[4 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + @square[i + j * @jsize1 + (k - 1) * @ksize1]) * @@c2)
                @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @@dz5tz1 * (@u[4 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] + @u[4 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2]) + @@zzcon3 * (@qs[i + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @qs[i + j * @jsize1 + k * @ksize1] + @qs[i + j * @jsize1 + (k - 1) * @ksize1]) + @@zzcon4 * (wp1 * wp1 - 2.0 * wijk * wijk + wm1 * wm1) + @@zzcon5 * (@u[4 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] * @rho_i[i + j * @jsize1 + (k + 1) * @ksize1] - 2.0 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @rho_i[i + j * @jsize1 + k * @ksize1] + @u[4 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] * @rho_i[i + j * @jsize1 + (k - 1) * @ksize1]) - @@tz2 * ((@@c1 * @u[4 + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] - @@c2 * @square[i + j * @jsize1 + (k + 1) * @ksize1]) * wp1 - (@@c1 * @u[4 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] - @@c2 * @square[i + j * @jsize1 + (k - 1) * @ksize1]) * wm1)
            end
        end
    end
    for j in 1..@grid_points[1] - 2 do
        for i in 1..@grid_points[0] - 2 do
            k = 1
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] + @u[m + i * @isize2 + j * @jsize2 + (k + 2) * @ksize2])
            end
            k = 2
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (-4.0 * @u[m + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] + @u[m + i * @isize2 + j * @jsize2 + (k + 2) * @ksize2])
            end
        end
    end
    for k in 3..@grid_points[2] - 4 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                for m in 0..4 do
                    @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + j * @jsize2 + (k - 2) * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2] + @u[m + i * @isize2 + j * @jsize2 + (k + 2) * @ksize2])
                end
            end
        end
    end
    for j in 1..@grid_points[1] - 2 do
        for i in 1..@grid_points[0] - 2 do
            k = @grid_points[2] - 3
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + j * @jsize2 + (k - 2) * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + 6.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2])
            end
            k = @grid_points[2] - 2
            for m in 0..4 do
                @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @@dssp * (@u[m + i * @isize2 + j * @jsize2 + (k - 2) * @ksize2] - 4.0 * @u[m + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2] + 5.0 * @u[m + i * @isize2 + j * @jsize2 + k * @ksize2])
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_rhsz)
    else
    end
    for k in 1..@grid_points[2] - 2 do
        for j in 1..@grid_points[1] - 2 do
            for i in 1..@grid_points[0] - 2 do
                for m in 0..4 do
                    @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] * @@dt
                end
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_rhs)
    else
    end
  end

  def print_lhs()
    count1 = 0; count2 = 0; count3 = 0; 
    for i in 0..5-1 do
        for j in 0..5-1 do
            for m in 0..@problem_size + 1-1 do
                count2 += @njac[i + j * @@isize4 + m * @@jsize4]
                count3 += @fjac[i + j * @@isize4 + m * @@jsize4]
                for k in 0..3-1 do
                    count1 += @lhs[i + j * @@isize4 + k * @@jsize4 + m * @@ksize4]
                end
            end
        end
    end
    print("lhs checksum is: ")
    puts(count1)
    print("fjac checksum is: ")
    puts(count3)
    print("njac checksum is: ")
    puts(count2)
  end

  def y_solve()
    
    if @timeron then
      @timer.start(@@t_ysolve)
    else
    end
    jsize = @grid_points[1] - 1
    for k in 1..@grid_points[2] - 2 do
        for i in 1..@grid_points[0] - 2 do
            for j in 0..jsize do
                @tmp1 = @rho_i[i + j * @jsize1 + k * @ksize1]
                @tmp2 = @tmp1 * @tmp1
                @tmp3 = @tmp1 * @tmp2
                @fjac[0 + 0 * @@isize4 + j * @@jsize4] = 0.0
                @fjac[0 + 1 * @@isize4 + j * @@jsize4] = 0.0
                @fjac[0 + 2 * @@isize4 + j * @@jsize4] = 1.0
                @fjac[0 + 3 * @@isize4 + j * @@jsize4] = 0.0
                @fjac[0 + 4 * @@isize4 + j * @@jsize4] = 0.0
                @fjac[1 + 0 * @@isize4 + j * @@jsize4] = -(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[1 + 1 * @@isize4 + j * @@jsize4] = @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[1 + 2 * @@isize4 + j * @@jsize4] = @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[1 + 3 * @@isize4 + j * @@jsize4] = 0.0
                @fjac[1 + 4 * @@isize4 + j * @@jsize4] = 0.0
                @fjac[2 + 0 * @@isize4 + j * @@jsize4] = -(@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2) + @@c2 * @qs[i + j * @jsize1 + k * @ksize1]
                @fjac[2 + 1 * @@isize4 + j * @@jsize4] = -@@c2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 2 * @@isize4 + j * @@jsize4] = (2.0 - @@c2) * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 3 * @@isize4 + j * @@jsize4] = -@@c2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 4 * @@isize4 + j * @@jsize4] = @@c2
                @fjac[3 + 0 * @@isize4 + j * @@jsize4] = -(@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[3 + 1 * @@isize4 + j * @@jsize4] = 0.0
                @fjac[3 + 2 * @@isize4 + j * @@jsize4] = @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 3 * @@isize4 + j * @@jsize4] = @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 4 * @@isize4 + j * @@jsize4] = 0.0
                @fjac[4 + 0 * @@isize4 + j * @@jsize4] = (@@c2 * 2.0 * @square[i + j * @jsize1 + k * @ksize1] - @@c1 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2
                @fjac[4 + 1 * @@isize4 + j * @@jsize4] = -@@c2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2
                @fjac[4 + 2 * @@isize4 + j * @@jsize4] = @@c1 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1 - @@c2 * (@qs[i + j * @jsize1 + k * @ksize1] + @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2)
                @fjac[4 + 3 * @@isize4 + j * @@jsize4] = -@@c2 * (@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[4 + 4 * @@isize4 + j * @@jsize4] = @@c1 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @njac[0 + 0 * @@isize4 + j * @@jsize4] = 0.0
                @njac[0 + 1 * @@isize4 + j * @@jsize4] = 0.0
                @njac[0 + 2 * @@isize4 + j * @@jsize4] = 0.0
                @njac[0 + 3 * @@isize4 + j * @@jsize4] = 0.0
                @njac[0 + 4 * @@isize4 + j * @@jsize4] = 0.0
                @njac[1 + 0 * @@isize4 + j * @@jsize4] = -@@c3c4 * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[1 + 1 * @@isize4 + j * @@jsize4] = @@c3c4 * @tmp1
                @njac[1 + 2 * @@isize4 + j * @@jsize4] = 0.0
                @njac[1 + 3 * @@isize4 + j * @@jsize4] = 0.0
                @njac[1 + 4 * @@isize4 + j * @@jsize4] = 0.0
                @njac[2 + 0 * @@isize4 + j * @@jsize4] = -@@con43 * @@c3c4 * @tmp2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[2 + 1 * @@isize4 + j * @@jsize4] = 0.0
                @njac[2 + 2 * @@isize4 + j * @@jsize4] = @@con43 * @@c3c4 * @tmp1
                @njac[2 + 3 * @@isize4 + j * @@jsize4] = 0.0
                @njac[2 + 4 * @@isize4 + j * @@jsize4] = 0.0
                @njac[3 + 0 * @@isize4 + j * @@jsize4] = -@@c3c4 * @tmp2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[3 + 1 * @@isize4 + j * @@jsize4] = 0.0
                @njac[3 + 2 * @@isize4 + j * @@jsize4] = 0.0
                @njac[3 + 3 * @@isize4 + j * @@jsize4] = @@c3c4 * @tmp1
                @njac[3 + 4 * @@isize4 + j * @@jsize4] = 0.0
                @njac[4 + 0 * @@isize4 + j * @@jsize4] = -(@@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - (@@con43 * @@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - (@@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - @@c1345 * @tmp2 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 1 * @@isize4 + j * @@jsize4] = (@@c3c4 - @@c1345) * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 2 * @@isize4 + j * @@jsize4] = (@@con43 * @@c3c4 - @@c1345) * @tmp2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 3 * @@isize4 + j * @@jsize4] = (@@c3c4 - @@c1345) * @tmp2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 4 * @@isize4 + j * @@jsize4] = (@@c1345) * @tmp1
            end
            lhsinit(@lhs, jsize)
            for j in 1..jsize - 1 do
                @tmp1 = @@dt * @@ty1
                @tmp2 = @@dt * @@ty2
                @lhs[0 + 0 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[0 + 0 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[0 + 0 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @@dy1
                @lhs[0 + 1 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[0 + 1 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[0 + 1 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[0 + 2 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[0 + 2 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[0 + 3 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[0 + 3 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[0 + 4 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[0 + 4 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[1 + 0 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[1 + 0 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[1 + 1 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[1 + 1 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @@dy2
                @lhs[1 + 2 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[1 + 2 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[1 + 2 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[1 + 3 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[1 + 3 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[1 + 4 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[1 + 4 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[2 + 0 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[2 + 0 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[2 + 1 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[2 + 1 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[2 + 2 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[2 + 2 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @@dy3
                @lhs[2 + 3 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[2 + 3 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[2 + 3 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[2 + 4 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[2 + 4 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[3 + 0 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[3 + 0 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[3 + 1 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[3 + 1 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[3 + 2 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[3 + 2 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[3 + 3 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[3 + 3 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @@dy4
                @lhs[3 + 4 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[3 + 4 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[3 + 4 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[4 + 0 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[4 + 0 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[4 + 1 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[4 + 1 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[4 + 2 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[4 + 2 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[4 + 3 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[4 + 3 * @@isize4 + (j - 1) * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4] = -@tmp2 * @fjac[4 + 4 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @njac[4 + 4 * @@isize4 + (j - 1) * @@jsize4] - @tmp1 * @@dy5
                @lhs[0 + 0 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[0 + 0 * @@isize4 + j * @@jsize4] + @tmp1 * 2.0 * @@dy1
                @lhs[0 + 1 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 1 * @@isize4 + j * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 2 * @@isize4 + j * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 3 * @@isize4 + j * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 4 * @@isize4 + j * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 0 * @@isize4 + j * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[1 + 1 * @@isize4 + j * @@jsize4] + @tmp1 * 2.0 * @@dy2
                @lhs[1 + 2 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 2 * @@isize4 + j * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 3 * @@isize4 + j * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 4 * @@isize4 + j * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 0 * @@isize4 + j * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 1 * @@isize4 + j * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[2 + 2 * @@isize4 + j * @@jsize4] + @tmp1 * 2.0 * @@dy3
                @lhs[2 + 3 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 3 * @@isize4 + j * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 4 * @@isize4 + j * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 0 * @@isize4 + j * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 1 * @@isize4 + j * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 2 * @@isize4 + j * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[3 + 3 * @@isize4 + j * @@jsize4] + @tmp1 * 2.0 * @@dy4
                @lhs[3 + 4 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 4 * @@isize4 + j * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 0 * @@isize4 + j * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 1 * @@isize4 + j * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 2 * @@isize4 + j * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 3 * @@isize4 + j * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[4 + 4 * @@isize4 + j * @@jsize4] + @tmp1 * 2.0 * @@dy5
                @lhs[0 + 0 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[0 + 0 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[0 + 0 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @@dy1
                @lhs[0 + 1 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[0 + 1 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[0 + 1 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[0 + 2 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[0 + 2 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[0 + 3 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[0 + 3 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[0 + 4 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[0 + 4 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[1 + 0 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[1 + 0 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[1 + 1 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[1 + 1 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @@dy2
                @lhs[1 + 2 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[1 + 2 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[1 + 2 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[1 + 3 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[1 + 3 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[1 + 4 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[1 + 4 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[2 + 0 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[2 + 0 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[2 + 1 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[2 + 1 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[2 + 2 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[2 + 2 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @@dy3
                @lhs[2 + 3 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[2 + 3 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[2 + 3 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[2 + 4 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[2 + 4 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[3 + 0 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[3 + 0 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[3 + 1 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[3 + 1 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[3 + 2 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[3 + 2 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[3 + 3 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[3 + 3 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @@dy4
                @lhs[3 + 4 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[3 + 4 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[3 + 4 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[4 + 0 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[4 + 0 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[4 + 1 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[4 + 1 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[4 + 2 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[4 + 2 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[4 + 3 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[4 + 3 * @@isize4 + (j + 1) * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] = @tmp2 * @fjac[4 + 4 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @njac[4 + 4 * @@isize4 + (j + 1) * @@jsize4] - @tmp1 * @@dy5
            end
            binvcrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + 0 * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + 0 * @@ksize4, @rhs, 0 + i * @isize2 + 0 * @jsize2 + k * @ksize2)
            for j in 1..jsize - 1 do
                matvec_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4, @rhs, 0 + i * @isize2 + (j - 1) * @jsize2 + k * @ksize2, @rhs, 0 + i * @isize2 + j * @jsize2 + k * @ksize2)
                matmul_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + j * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + (j - 1) * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4)
                binvcrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + j * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + j * @@ksize4, @rhs, 0 + i * @isize2 + j * @jsize2 + k * @ksize2)
            end
            matvec_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + jsize * @@ksize4, @rhs, 0 + i * @isize2 + (jsize - 1) * @jsize2 + k * @ksize2, @rhs, 0 + i * @isize2 + jsize * @jsize2 + k * @ksize2)
            matmul_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + jsize * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + (jsize - 1) * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + jsize * @@ksize4)
            binvrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + jsize * @@ksize4, @rhs, 0 + i * @isize2 + jsize * @jsize2 + k * @ksize2)
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
            j = jsize - 1
            while j >= 0 do
                for m in 0..@@BLOCK_SIZE - 1 do
                    for n in 0..@@BLOCK_SIZE - 1 do
                        @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] = @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] - @lhs[m + n * @@isize4 + @@cc * @@jsize4 + j * @@ksize4] * @rhs[n + i * @isize2 + (j + 1) * @jsize2 + k * @ksize2]
                    end
                end
              j -= 1
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_ysolve)
    else
    end
  end

  def z_solve()
    
    if @timeron then
      @timer.start(@@t_zsolve)
    else
    end
    ksize = @grid_points[2] - 1
    for j in 1..@grid_points[1] - 2 do
        for i in 1..@grid_points[0] - 2 do
            for k in 0..ksize do
                @tmp1 = 1.0 / @u[0 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @tmp2 = @tmp1 * @tmp1
                @tmp3 = @tmp1 * @tmp2
                @fjac[0 + 0 * @@isize4 + k * @@jsize4] = 0.0
                @fjac[0 + 1 * @@isize4 + k * @@jsize4] = 0.0
                @fjac[0 + 2 * @@isize4 + k * @@jsize4] = 0.0
                @fjac[0 + 3 * @@isize4 + k * @@jsize4] = 1.0
                @fjac[0 + 4 * @@isize4 + k * @@jsize4] = 0.0
                @fjac[1 + 0 * @@isize4 + k * @@jsize4] = -(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[1 + 1 * @@isize4 + k * @@jsize4] = @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[1 + 2 * @@isize4 + k * @@jsize4] = 0.0
                @fjac[1 + 3 * @@isize4 + k * @@jsize4] = @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[1 + 4 * @@isize4 + k * @@jsize4] = 0.0
                @fjac[2 + 0 * @@isize4 + k * @@jsize4] = -(@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[2 + 1 * @@isize4 + k * @@jsize4] = 0.0
                @fjac[2 + 2 * @@isize4 + k * @@jsize4] = @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 3 * @@isize4 + k * @@jsize4] = @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[2 + 4 * @@isize4 + k * @@jsize4] = 0.0
                @fjac[3 + 0 * @@isize4 + k * @@jsize4] = -(@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2) + @@c2 * @qs[i + j * @jsize1 + k * @ksize1]
                @fjac[3 + 1 * @@isize4 + k * @@jsize4] = -@@c2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 2 * @@isize4 + k * @@jsize4] = -@@c2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 3 * @@isize4 + k * @@jsize4] = (2.0 - @@c2) * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @fjac[3 + 4 * @@isize4 + k * @@jsize4] = @@c2
                @fjac[4 + 0 * @@isize4 + k * @@jsize4] = (@@c2 * 2.0 * @square[i + j * @jsize1 + k * @ksize1] - @@c1 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2
                @fjac[4 + 1 * @@isize4 + k * @@jsize4] = -@@c2 * (@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[4 + 2 * @@isize4 + k * @@jsize4] = -@@c2 * (@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]) * @tmp2
                @fjac[4 + 3 * @@isize4 + k * @@jsize4] = @@c1 * (@u[4 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1) - @@c2 * (@qs[i + j * @jsize1 + k * @ksize1] + @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp2)
                @fjac[4 + 4 * @@isize4 + k * @@jsize4] = @@c1 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2] * @tmp1
                @njac[0 + 0 * @@isize4 + k * @@jsize4] = 0.0
                @njac[0 + 1 * @@isize4 + k * @@jsize4] = 0.0
                @njac[0 + 2 * @@isize4 + k * @@jsize4] = 0.0
                @njac[0 + 3 * @@isize4 + k * @@jsize4] = 0.0
                @njac[0 + 4 * @@isize4 + k * @@jsize4] = 0.0
                @njac[1 + 0 * @@isize4 + k * @@jsize4] = -@@c3c4 * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[1 + 1 * @@isize4 + k * @@jsize4] = @@c3c4 * @tmp1
                @njac[1 + 2 * @@isize4 + k * @@jsize4] = 0.0
                @njac[1 + 3 * @@isize4 + k * @@jsize4] = 0.0
                @njac[1 + 4 * @@isize4 + k * @@jsize4] = 0.0
                @njac[2 + 0 * @@isize4 + k * @@jsize4] = -@@c3c4 * @tmp2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[2 + 1 * @@isize4 + k * @@jsize4] = 0.0
                @njac[2 + 2 * @@isize4 + k * @@jsize4] = @@c3c4 * @tmp1
                @njac[2 + 3 * @@isize4 + k * @@jsize4] = 0.0
                @njac[2 + 4 * @@isize4 + k * @@jsize4] = 0.0
                @njac[3 + 0 * @@isize4 + k * @@jsize4] = -@@con43 * @@c3c4 * @tmp2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[3 + 1 * @@isize4 + k * @@jsize4] = 0.0
                @njac[3 + 2 * @@isize4 + k * @@jsize4] = 0.0
                @njac[3 + 3 * @@isize4 + k * @@jsize4] = @@con43 * @@c3 * @@c4 * @tmp1
                @njac[3 + 4 * @@isize4 + k * @@jsize4] = 0.0
                @njac[4 + 0 * @@isize4 + k * @@jsize4] = -(@@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[1 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - (@@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[2 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - (@@con43 * @@c3c4 - @@c1345) * @tmp3 * (Math_dot_pow(@u[3 + i * @isize2 + j * @jsize2 + k * @ksize2], 2)) - @@c1345 * @tmp2 * @u[4 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 1 * @@isize4 + k * @@jsize4] = (@@c3c4 - @@c1345) * @tmp2 * @u[1 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 2 * @@isize4 + k * @@jsize4] = (@@c3c4 - @@c1345) * @tmp2 * @u[2 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 3 * @@isize4 + k * @@jsize4] = (@@con43 * @@c3c4 - @@c1345) * @tmp2 * @u[3 + i * @isize2 + j * @jsize2 + k * @ksize2]
                @njac[4 + 4 * @@isize4 + k * @@jsize4] = (@@c1345) * @tmp1
            end
            lhsinit(@lhs, ksize)
            for k in 1..ksize - 1 do
                @tmp1 = @@dt * @@tz1
                @tmp2 = @@dt * @@tz2
                @lhs[0 + 0 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[0 + 0 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[0 + 0 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @@dz1
                @lhs[0 + 1 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[0 + 1 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[0 + 1 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[0 + 2 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[0 + 2 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[0 + 3 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[0 + 3 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[0 + 4 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[0 + 4 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[1 + 0 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[1 + 0 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[1 + 1 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[1 + 1 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @@dz2
                @lhs[1 + 2 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[1 + 2 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[1 + 2 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[1 + 3 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[1 + 3 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[1 + 4 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[1 + 4 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[2 + 0 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[2 + 0 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[2 + 1 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[2 + 1 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[2 + 2 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[2 + 2 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @@dz3
                @lhs[2 + 3 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[2 + 3 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[2 + 3 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[2 + 4 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[2 + 4 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[3 + 0 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[3 + 0 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[3 + 1 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[3 + 1 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[3 + 2 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[3 + 2 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[3 + 3 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[3 + 3 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @@dz4
                @lhs[3 + 4 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[3 + 4 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[3 + 4 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[4 + 0 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[4 + 0 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[4 + 1 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[4 + 1 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[4 + 2 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[4 + 2 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[4 + 3 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[4 + 3 * @@isize4 + (k - 1) * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4] = -@tmp2 * @fjac[4 + 4 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @njac[4 + 4 * @@isize4 + (k - 1) * @@jsize4] - @tmp1 * @@dz5
                @lhs[0 + 0 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[0 + 0 * @@isize4 + k * @@jsize4] + @tmp1 * 2.0 * @@dz1
                @lhs[0 + 1 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 1 * @@isize4 + k * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 2 * @@isize4 + k * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 3 * @@isize4 + k * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[0 + 4 * @@isize4 + k * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 0 * @@isize4 + k * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[1 + 1 * @@isize4 + k * @@jsize4] + @tmp1 * 2.0 * @@dz2
                @lhs[1 + 2 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 2 * @@isize4 + k * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 3 * @@isize4 + k * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[1 + 4 * @@isize4 + k * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 0 * @@isize4 + k * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 1 * @@isize4 + k * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[2 + 2 * @@isize4 + k * @@jsize4] + @tmp1 * 2.0 * @@dz3
                @lhs[2 + 3 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 3 * @@isize4 + k * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[2 + 4 * @@isize4 + k * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 0 * @@isize4 + k * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 1 * @@isize4 + k * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 2 * @@isize4 + k * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[3 + 3 * @@isize4 + k * @@jsize4] + @tmp1 * 2.0 * @@dz4
                @lhs[3 + 4 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[3 + 4 * @@isize4 + k * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 0 * @@isize4 + k * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 1 * @@isize4 + k * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 2 * @@isize4 + k * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = @tmp1 * 2.0 * @njac[4 + 3 * @@isize4 + k * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4] = 1.0 + @tmp1 * 2.0 * @njac[4 + 4 * @@isize4 + k * @@jsize4] + @tmp1 * 2.0 * @@dz5
                @lhs[0 + 0 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[0 + 0 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[0 + 0 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @@dz1
                @lhs[0 + 1 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[0 + 1 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[0 + 1 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[0 + 2 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[0 + 2 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[0 + 2 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[0 + 3 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[0 + 3 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[0 + 3 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[0 + 4 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[0 + 4 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[0 + 4 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[1 + 0 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[1 + 0 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[1 + 0 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[1 + 1 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[1 + 1 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[1 + 1 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @@dz2
                @lhs[1 + 2 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[1 + 2 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[1 + 2 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[1 + 3 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[1 + 3 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[1 + 3 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[1 + 4 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[1 + 4 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[1 + 4 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[2 + 0 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[2 + 0 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[2 + 0 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[2 + 1 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[2 + 1 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[2 + 1 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[2 + 2 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[2 + 2 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[2 + 2 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @@dz3
                @lhs[2 + 3 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[2 + 3 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[2 + 3 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[2 + 4 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[2 + 4 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[2 + 4 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[3 + 0 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[3 + 0 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[3 + 0 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[3 + 1 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[3 + 1 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[3 + 1 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[3 + 2 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[3 + 2 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[3 + 2 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[3 + 3 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[3 + 3 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[3 + 3 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @@dz4
                @lhs[3 + 4 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[3 + 4 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[3 + 4 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[4 + 0 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[4 + 0 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[4 + 0 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[4 + 1 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[4 + 1 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[4 + 1 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[4 + 2 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[4 + 2 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[4 + 2 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[4 + 3 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[4 + 3 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[4 + 3 * @@isize4 + (k + 1) * @@jsize4]
                @lhs[4 + 4 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] = @tmp2 * @fjac[4 + 4 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @njac[4 + 4 * @@isize4 + (k + 1) * @@jsize4] - @tmp1 * @@dz5
            end
            binvcrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + 0 * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + 0 * @@ksize4, @rhs, 0 + i * @isize2 + j * @jsize2 + 0 * @ksize2)
            for k in 1..ksize - 1 do
                matvec_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4, @rhs, 0 + i * @isize2 + j * @jsize2 + (k - 1) * @ksize2, @rhs, 0 + i * @isize2 + j * @jsize2 + k * @ksize2)
                matmul_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + k * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + (k - 1) * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4)
                binvcrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + k * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + k * @@ksize4, @rhs, 0 + i * @isize2 + j * @jsize2 + k * @ksize2)
            end
            matvec_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + ksize * @@ksize4, @rhs, 0 + i * @isize2 + j * @jsize2 + (ksize - 1) * @ksize2, @rhs, 0 + i * @isize2 + j * @jsize2 + ksize * @ksize2)
            matmul_sub(@lhs, 0 + 0 * @@isize4 + @@aa * @@jsize4 + ksize * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@cc * @@jsize4 + (ksize - 1) * @@ksize4, @lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + ksize * @@ksize4)
            binvrhs(@lhs, 0 + 0 * @@isize4 + @@bb * @@jsize4 + ksize * @@ksize4, @rhs, 0 + i * @isize2 + j * @jsize2 + ksize * @ksize2)
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
            k = ksize - 1
            while k >= 0 do
                for m in 0..@@BLOCK_SIZE - 1 do
                    for n in 0..@@BLOCK_SIZE - 1 do
                        @rhs[m + i * @isize2 + j * @jsize2 + k * @ksize2] += -@lhs[m + n * @@isize4 + @@cc * @@jsize4 + k * @@ksize4] * @rhs[n + i * @isize2 + j * @jsize2 + (k + 1) * @ksize2]
                    end
                end
              k -= 1
            end
        end
    end
    if @timeron then
      @timer.stop(@@t_zsolve)
    else
    end
  end

  def checkSum(arr)
    csum = 0.0; 
    for k in 0..@grid_points[2] - 1 do
        for j in 0..@grid_points[1] - 1 do
            for i in 0..@grid_points[0] - 1 do
                for m in 0..4 do
                    offset = m + i * @isize2 + j * @jsize2 + k * @ksize2; 
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
    puts("BT: is about to be garbage collected")
    super.finalize()
  end

  def setupThreads(bt)
    @master = bt
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
    @xsolver = Array.new(@num_threads, 0.0)
    @ysolver = Array.new(@num_threads, 0.0)
    @zsolver = Array.new(@num_threads, 0.0)
    @rhsadder = Array.new(@num_threads, 0.0)
    for ii in 0..@num_threads-1 do
      @rhscomputer[ii] = RHSCompute.new(bt, partition1[ii][0], partition1[ii][1], partition2[ii][0], partition2[ii][1])
      @rhscomputer[ii].extend(MonitorMixin)
      @rhscomputer[ii].id = ii
      @rhscomputer[ii].start()
      @xsolver[ii] = XSolver.new(bt, partition2[ii][0], partition2[ii][1])
      @xsolver[ii].extend(MonitorMixin)
      @xsolver[ii].id = ii
      @xsolver[ii].start()
      @ysolver[ii] = YSolver.new(bt, partition2[ii][0], partition2[ii][1])
      @ysolver[ii].extend(MonitorMixin)
      @ysolver[ii].id = ii
      @ysolver[ii].start()
      @zsolver[ii] = ZSolver.new(bt, partition2[ii][0], partition2[ii][1])
      @zsolver[ii].extend(MonitorMixin)
      @zsolver[ii].id = ii
      @zsolver[ii].start()
      @rhsadder[ii] = RHSAdder.new(bt, partition2[ii][0], partition2[ii][1])
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
    bt =  nil ; 
    BMArgs.parseCmdLineArgs(argv, "BT")
    clss = BMArgs.clss
    np = BMArgs.num_threads
    serial = BMArgs.serial 
    #begin
    bt = BT.new(clss, np, serial)
    bt.extend(MonitorMixin)
    #rescue OutOfMemoryError => e
    #    BMArgs.outOfMemoryMessage()
    #    System.exit(0)
    #ensure
    #end
    bt.runBenchMark()
  end

main(ARGV)
