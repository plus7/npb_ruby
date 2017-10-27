# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
#  Authors: D. Bailey
#	    W. Saphir
#  Translation to Java and MultiThreaded Code
# 	    M. Frumkin
#	    M. Schultz
#  Translation to Ruby
#           T. Nose
#           H. Tomari
# -------------------------------------------------------------------------

require "mda"
require "mathutil"
require "FTThreads/FTBase"
require "FTThreads/FFTThread"
require "FTThreads/EvolveThread"
require "BMInOut/BMResults"
require "BMInOut/BMArgs"
require "Timer"
require "Random"
require "Runnable"
require "monitor"
class FT < FTBase
  def initialize(clss, np, ser)
    @bid = -1
    @serial = true
    @done = false
    super(clss, np, ser)
    @serial = ser
  end

  def run()
    runBenchMark()
  end

  def runBenchMark()
    BMArgs.banner(@@BMName, @CLASS, @serial, @num_threads)
    puts(" Size = " + @nx.to_s + " X " + @ny.to_s + " X " + @nz.to_s + " niter = " + @niter_default.to_s)
    setTimers()
    @timer.resetAllTimers()
    if @serial then
      appft_serial()
    else
      appft()
    end
    if @timeron then
      @timer.start(14)
    else
    end
    verified = verify(4, @nx, @ny, @nz, @niter_default, @checksum); 
    if @timeron then
      @timer.stop(14)
    else
    end
    @timer.stop(1)
    time = @timer.readTimer(1); 
    @results = BMResults.new(@@BMName, @CLASS, @nx, @ny, @nz, @niter_default, time, getMFLOPS(time, @nx, @ny, @nz), "floating point", verified, @serial, @num_threads, @bid)
    @results.print()
    if @timeron then
      printTimers()
    else
    end
    @done = true
  end

  def appft_serial()
    if @timeron then
      @timer.start(2)
    else
    end
    initial_conditions(@xtr, @ny, @nx, @nz)
    CompExp(@nx, @exp1)
    CompExp(@ny, @exp2)
    CompExp(@nz, @exp3)
    fftXYZ(1, @xtr, @exp2, @exp1, @exp3, @ny, @nx, @nz)
    if @timeron then
      @timer.stop(2)
    else
    end
    @timer.start(1)
    if @timeron then
      @timer.start(12)
    else
    end
    initial_conditions(@xtr, @ny, @nx, @nz)
    if @timeron then
      @timer.stop(12)
    else
    end
    if @timeron then
      @timer.start(15)
    else
    end
    fftXYZ(1, @xtr, @exp2, @exp1, @exp3, @ny, @nx, @nz)
    if @timeron then
      @timer.stop(15)
    else
    end
    ap = (-4.0 * @@alpha * Math_dot_pow(@@pi, 2)); 
    n12 = @nx / 2; 
    n22 = @ny / 2; 
    n32 = @nz / 2; 
    for it in 0..@niter_default-1 do
        if @timeron then
          @timer.start(11)
        else
        end
        for i in 0..@nx-1 do
            ii = i - ((i) / n12) * @nx; 
            ii2 = ii * ii; 
            for k in 0..@nz-1 do
                kk = k - ((k) / n32) * @nz; 
                ik2 = ii2 + kk * kk; 
                for j in 0..@ny-1 do
                    jj = j - ((j) / n22) * @ny; 
                    @xnt[@@REAL + j * @isize4 + k * @jsize4 + i * @ksize4] = @xtr[@@REAL + j * @isize3 + i * @jsize3 + k * @ksize3] * Math.exp((ap * (jj * jj + ik2)) * (it + 1))
                    @xnt[@@IMAG + j * @isize4 + k * @jsize4 + i * @ksize4] = @xtr[@@IMAG + j * @isize3 + i * @jsize3 + k * @ksize3] * Math.exp((ap * (jj * jj + ik2)) * (it + 1))
                end
            end
        end
        if @timeron then
          @timer.stop(11)
        else
        end
        if @timeron then
          @timer.start(15)
        else
        end
        fftXYZ(-1, @xnt, @exp2, @exp3, @exp1, @ny, @nz, @nx)
        if @timeron then
          @timer.stop(15)
        else
        end
        if @timeron then
          @timer.start(10)
        else
        end
        calculateChecksum(@checksum, @@REAL + it * @isize2, it, @xnt, @ny, @nz, @nx)
        if @timeron then
          @timer.stop(10)
        else
        end
    end
  end

  def appft()
    if @timeron then
      @timer.start(2)
    else
    end
    initial_conditions(@xtr, @ny, @nx, @nz)
    CompExp(@nx, @exp1)
    CompExp(@ny, @exp2)
    CompExp(@nz, @exp3)
    setupThreads( self )
    for m in 0..@num_threads-1 do
      @doFFT[m].synchronize do
          @doFFT[m].setVariables(1, false, @xtr, @exp2, @exp1, @exp3)
      end
    end
    doFFT()
    doFFT()
    doFFT()
    if @timeron then
      @timer.stop(2)
    else
    end
    @timer.start(1)
    if @timeron then
      @timer.start(12)
    else
    end
    initial_conditions(@xtr, @ny, @nx, @nz)
    if @timeron then
      @timer.stop(12)
    else
    end
    if @timeron then
      @timer.start(15)
    else
    end
    for m in 0..@num_threads-1 do
      @doFFT[m].synchronize do
          @doFFT[m].setVariables(1, false, @xtr, @exp2, @exp1, @exp3)
      end
    end
    doFFT()
    doFFT()
    doFFT()
    if @timeron then
      @timer.stop(15)
    else
    end
    for it in 0..@niter_default-1 do
        if @timeron then
          @timer.start(11)
        else
        end
        doEvolve(it)
        if @timeron then
          @timer.stop(11)
        else
        end
        if @timeron then
          @timer.start(15)
        else
        end
        for m in 0..@num_threads-1 do
          @doFFT[m].synchronize do
              @doFFT[m].setVariables(-1, true, @xnt, @exp2, @exp3, @exp1)
          end
        end
        if @timeron then
          @timer.start(3)
        else
        end
        if @timeron then
          @timer.start(7)
        else
        end
        doFFT()
        if @timeron then
          @timer.stop(7)
        else
        end
        if @timeron then
          @timer.start(8)
        else
        end
        doFFT()
        if @timeron then
          @timer.stop(8)
        else
        end
        if @timeron then
          @timer.start(9)
        else
        end
        doFFT()
        if @timeron then
          @timer.stop(9)
        else
        end
        if @timeron then
          @timer.stop(3)
        else
        end
        if @timeron then
          @timer.stop(15)
        else
        end
        if @timeron then
          @timer.start(10)
        else
        end
        calculateChecksum(@checksum, @@REAL + it * @isize2, it, @xnt, @ny, @nz, @nx)
        if @timeron then
          @timer.stop(10)
        else
        end
    end
  end

  def doFFT()
    synchronize {
    
    for m in 0..@num_threads-1 do
      @doFFT[m].synchronize do
          @doFFT[m].done = false
          @doFFT[m].cond.signal
      end
    end
    for m in 0..@num_threads-1 do
      while not @doFFT[m].done do
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

  def doEvolve(it)
    synchronize {
    
    for m in 0..@num_threads-1 do
      @doEvolve[m].synchronize do
          @doEvolve[m].done = false
          @doEvolve[m].kt = it
          @doEvolve[m].cond.signal
      end
    end
    for m in 0..@num_threads-1 do
      while not @doEvolve[m].done do
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

  def setTimers()
    #f1 = File.new("timer.flag"); 
    @timeron = false
    if File.exist?("timer.flag") then
      @timeron = true
    else
    end
  end

  def printTimers()
    #fmt = DecimalFormat.new("0.000"); 
    puts("  SECTION   Time (secs)")
    puts("FT time =		      " + (@timer.readTimer(1)).to_s)
    puts("WarmUp time =		      " + (@timer.readTimer(2)).to_s)
    puts("ffXYZ body time =	      " + (@timer.readTimer(3)))
    puts("Swarztrauber body time =      " + (@timer.readTimer(4).to_s))
    puts("Redistribution time =	      " + (@timer.readTimer(5)).to_s)
    puts("Transposition time =	      " + (@timer.readTimer(6)).to_s)
    puts("X time =		      " + (@timer.readTimer(7)).to_s)
    puts("Y time =		      " + (@timer.readTimer(8)).to_s)
    puts("Z time =		      " + (@timer.readTimer(9)).to_s)
    puts("calculateChecksum =	      " + (@timer.readTimer(10)).to_s)
    puts("evolve =		      " + (@timer.readTimer(11)).to_s)
    puts("compute_initial_conditions =  " + (@timer.readTimer(12)).to_s)
    puts("twiddle =		      " + (@timer.readTimer(13)).to_s)
    puts("verify =		      " + (@timer.readTimer(14)).to_s)
    puts("fftXYZ =		      " + (@timer.readTimer(15)).to_s)
  end

  def getMFLOPS(total_time, nx, ny, nz)
    mflops = 0.0; 
    ntotal = nx * ny * nz; 
    if total_time > 0 then
        mflops = 14.8157 + 7.19641 * Math.log(ntotal) + (5.23518 + 7.21113 * Math.log(ntotal)) * @niter_default
        mflops *= ntotal / (total_time * 1000000.0)
    else
    end
    return mflops
  end

  def calculateChecksum(csum, csmffst, iterN, u, d1, d2, d3)
    
    isize3 = 2; jsize3 = isize3 * (d1 + 1); ksize3 = jsize3 * d2; 
    csum[@@REAL + csmffst] = 0.0
    csum[@@IMAG + csmffst] = 0.0
    csumr = 0.0; csumi = 0.0; 
    for i in 1..1024 do
        ii = (1 * i) % d3
        ji = (3 * i) % d1
        ki = (5 * i) % d2
        csumr += u[@@REAL + ji * isize3 + ki * jsize3 + ii * ksize3]
        csumi += u[@@IMAG + ji * isize3 + ki * jsize3 + ii * ksize3]
    end
    csum[@@REAL + csmffst] = csumr / (d1 * d2 * d3)
    csum[@@IMAG + csmffst] = csumi / (d1 * d2 * d3)
  end

  def fftXYZ(sign, x, exp1, exp2, exp3, n1, n2, n3)
    i = 0; j = 0; 
    isize3 = 2; 
    jsize3 = isize3 * (n1 + 1)
    ksize3 = jsize3 * n2
    isize1 = 2; jsize1 = 2 * (n2 + 1); 
    if @timeron then
      @timer.start(3)
    else
    end
    log = ilog2(n2)
    if @timeron then
      @timer.start(7)
    else
    end
    for k in 0..n3-1 do
      swarztrauber(sign, log, n1, n2, x, k * ksize3, n1, exp2, @scr)
    end
    if @timeron then
      @timer.stop(7)
    else
    end
    log = ilog2(n1)
    if @timeron then
      @timer.start(8)
    else
    end
    for k in 0..n3-1 do
        for j in 0..n2-1 do
            for i in 0..n1-1 do
                @plane[@@REAL + j * isize1 + i * jsize1] = x[@@REAL + i * isize3 + j * jsize3 + k * ksize3]
                @plane[@@IMAG + j * isize1 + i * jsize1] = x[@@IMAG + i * isize3 + j * jsize3 + k * ksize3]
            end
        end
        swarztrauber(sign, log, n2, n1, @plane, 0, n2, exp1, @scr)
        for j in 0..n2-1 do
            for i in 0..n1-1 do
                x[@@REAL + i * isize3 + j * jsize3 + k * ksize3] = @plane[@@REAL + j * isize1 + i * jsize1]
                x[@@IMAG + i * isize3 + j * jsize3 + k * ksize3] = @plane[@@IMAG + j * isize1 + i * jsize1]
            end
        end
    end
    if @timeron then
      @timer.stop(8)
    else
    end
    log = ilog2(n3)
    if @timeron then
      @timer.start(9)
    else
    end
    jsize1 = 2 * (n1 + 1)
    for k in 0..n2-1 do
        for i in 0..n3-1 do
            for j in 0..n1-1 do
                @plane[@@REAL + j * isize1 + i * jsize1] = x[@@REAL + j * isize3 + k * jsize3 + i * ksize3]
                @plane[@@IMAG + j * isize1 + i * jsize1] = x[@@IMAG + j * isize3 + k * jsize3 + i * ksize3]
            end
        end
        swarztrauber(sign, log, n1, n3, @plane, 0, n1, exp3, @scr)
        for i in 0..n3-1 do
            for j in 0..n1-1 do
                x[@@REAL + j * isize3 + k * jsize3 + i * ksize3] = @plane[@@REAL + j * isize1 + i * jsize1]
                x[@@IMAG + j * isize3 + k * jsize3 + i * ksize3] = @plane[@@IMAG + j * isize1 + i * jsize1]
            end
        end
    end
    if @timeron then
      @timer.stop(9)
    else
    end
    if @timeron then
      @timer.stop(3)
    else
    end
  end

  def verify(ires, n1, n2, n3, nt, cksum)
    verified = -1; 
    temp = Array.new(@niter_default, 0.0); 
    cexpd = Array.new(2 * 21, 0.0); 
    if (n1 == 64) and (n2 == 64) and (n3 == 64) and (nt == 6) then
        cexpd[@@REAL + 0 * 2] = 554.6087004964
        cexpd[@@REAL + 1 * 2] = 554.6385409189
        cexpd[@@REAL + 2 * 2] = 554.6148406171
        cexpd[@@REAL + 3 * 2] = 554.5423607415
        cexpd[@@REAL + 4 * 2] = 554.4255039624
        cexpd[@@REAL + 5 * 2] = 554.2683411902
        cexpd[@@IMAG + 0 * 2] = 484.5363331978
        cexpd[@@IMAG + 1 * 2] = 486.5304269511
        cexpd[@@IMAG + 2 * 2] = 488.3910722336
        cexpd[@@IMAG + 3 * 2] = 490.1273169046
        cexpd[@@IMAG + 4 * 2] = 491.7475857993
        cexpd[@@IMAG + 5 * 2] = 493.2597244941
    else
      if (n1 == 128) and (n2 == 128) and (n3 == 32) and (nt == 6) then
          cexpd[@@REAL + 0 * 2] = 567.3612178944
          cexpd[@@REAL + 1 * 2] = 563.1436885271
          cexpd[@@REAL + 2 * 2] = 559.4024089970
          cexpd[@@REAL + 3 * 2] = 556.0698047020
          cexpd[@@REAL + 4 * 2] = 553.0898991250
          cexpd[@@REAL + 5 * 2] = 550.4159734538
          cexpd[@@IMAG + 0 * 2] = 529.3246849175
          cexpd[@@IMAG + 1 * 2] = 528.2149986629
          cexpd[@@IMAG + 2 * 2] = 527.0996558037
          cexpd[@@IMAG + 3 * 2] = 526.0027904925
          cexpd[@@IMAG + 4 * 2] = 524.9400845633
          cexpd[@@IMAG + 5 * 2] = 523.9212247086
      else
        if (n1 == 256) and (n2 == 256) and (n3 == 128) and (nt == 6) then
            cexpd[@@REAL + 0 * 2] = 504.6735008193
            cexpd[@@REAL + 1 * 2] = 505.9412319734
            cexpd[@@REAL + 2 * 2] = 506.9376896287
            cexpd[@@REAL + 3 * 2] = 507.7892868474
            cexpd[@@REAL + 4 * 2] = 508.5233095391
            cexpd[@@REAL + 5 * 2] = 509.1487099959
            cexpd[@@IMAG + 0 * 2] = 511.4047905510
            cexpd[@@IMAG + 1 * 2] = 509.8809666433
            cexpd[@@IMAG + 2 * 2] = 509.8144042213
            cexpd[@@IMAG + 3 * 2] = 510.1336130759
            cexpd[@@IMAG + 4 * 2] = 510.4914655194
            cexpd[@@IMAG + 5 * 2] = 510.7917842803
        else
          if (n1 == 512) and (n2 == 256) and (n3 == 256) and (nt == 20) then
              cexpd[@@REAL + 0 * 2] = 517.7643571579
              cexpd[@@REAL + 1 * 2] = 515.4521291263
              cexpd[@@REAL + 2 * 2] = 514.6409228649
              cexpd[@@REAL + 3 * 2] = 514.2378756213
              cexpd[@@REAL + 4 * 2] = 513.9626667737
              cexpd[@@REAL + 5 * 2] = 513.7423460082
              cexpd[@@REAL + 6 * 2] = 513.5547056878
              cexpd[@@REAL + 7 * 2] = 513.3910925466
              cexpd[@@REAL + 8 * 2] = 513.2470705390
              cexpd[@@REAL + 9 * 2] = 513.1197729984
              cexpd[@@REAL + 10 * 2] = 513.0070319283
              cexpd[@@REAL + 11 * 2] = 512.9070537032
              cexpd[@@REAL + 12 * 2] = 512.8182883502
              cexpd[@@REAL + 13 * 2] = 512.7393733383
              cexpd[@@REAL + 14 * 2] = 512.6691062020
              cexpd[@@REAL + 15 * 2] = 512.6064276004
              cexpd[@@REAL + 16 * 2] = 512.5504076570
              cexpd[@@REAL + 17 * 2] = 512.5002331720
              cexpd[@@REAL + 18 * 2] = 512.4551951846
              cexpd[@@REAL + 19 * 2] = 512.4146770029
              cexpd[@@IMAG + 0 * 2] = 507.7803458597
              cexpd[@@IMAG + 1 * 2] = 508.8249431599
              cexpd[@@IMAG + 2 * 2] = 509.6208912659
              cexpd[@@IMAG + 3 * 2] = 510.1023387619
              cexpd[@@IMAG + 4 * 2] = 510.3976610617
              cexpd[@@IMAG + 5 * 2] = 510.5948019802
              cexpd[@@IMAG + 6 * 2] = 510.7404165783
              cexpd[@@IMAG + 7 * 2] = 510.8576573661
              cexpd[@@IMAG + 8 * 2] = 510.9577278523
              cexpd[@@IMAG + 9 * 2] = 511.0460304483
              cexpd[@@IMAG + 10 * 2] = 511.1252433800
              cexpd[@@IMAG + 11 * 2] = 511.1968077718
              cexpd[@@IMAG + 12 * 2] = 511.2616233064
              cexpd[@@IMAG + 13 * 2] = 511.3203605551
              cexpd[@@IMAG + 14 * 2] = 511.3735928093
              cexpd[@@IMAG + 15 * 2] = 511.4218460548
              cexpd[@@IMAG + 16 * 2] = 511.4656139760
              cexpd[@@IMAG + 17 * 2] = 511.5053595966
              cexpd[@@IMAG + 18 * 2] = 511.5415130407
              cexpd[@@IMAG + 19 * 2] = 511.5744692211
          else
            if (n1 == 512) and (n2 == 512) and (n3 == 512) and (nt == 20) then
                cexpd[@@REAL + 0 * 2] = 519.5078707457
                cexpd[@@REAL + 1 * 2] = 515.5422171134
                cexpd[@@REAL + 2 * 2] = 514.4678022222
                cexpd[@@REAL + 3 * 2] = 514.0150594328
                cexpd[@@REAL + 4 * 2] = 513.7550426810
                cexpd[@@REAL + 5 * 2] = 513.5811056728
                cexpd[@@REAL + 6 * 2] = 513.4569343165
                cexpd[@@REAL + 7 * 2] = 513.3651975661
                cexpd[@@REAL + 8 * 2] = 513.2955192805
                cexpd[@@REAL + 9 * 2] = 513.2410471738
                cexpd[@@REAL + 10 * 2] = 513.1971141679
                cexpd[@@REAL + 11 * 2] = 513.1605205716
                cexpd[@@REAL + 12 * 2] = 513.1290734194
                cexpd[@@REAL + 13 * 2] = 513.1012720314
                cexpd[@@REAL + 14 * 2] = 513.0760908195
                cexpd[@@REAL + 15 * 2] = 513.0528295923
                cexpd[@@REAL + 16 * 2] = 513.0310107773
                cexpd[@@REAL + 17 * 2] = 513.0103090133
                cexpd[@@REAL + 18 * 2] = 512.9905029333
                cexpd[@@REAL + 19 * 2] = 512.9714421109
                cexpd[@@IMAG + 0 * 2] = 514.9019699238
                cexpd[@@IMAG + 1 * 2] = 512.7578201997
                cexpd[@@IMAG + 2 * 2] = 512.2251847514
                cexpd[@@IMAG + 3 * 2] = 512.1090289018
                cexpd[@@IMAG + 4 * 2] = 512.1143685824
                cexpd[@@IMAG + 5 * 2] = 512.1496764568
                cexpd[@@IMAG + 6 * 2] = 512.1870921893
                cexpd[@@IMAG + 7 * 2] = 512.2193250322
                cexpd[@@IMAG + 8 * 2] = 512.2454735794
                cexpd[@@IMAG + 9 * 2] = 512.2663649603
                cexpd[@@IMAG + 10 * 2] = 512.2830879827
                cexpd[@@IMAG + 11 * 2] = 512.2965869718
                cexpd[@@IMAG + 12 * 2] = 512.3075927445
                cexpd[@@IMAG + 13 * 2] = 512.3166486553
                cexpd[@@IMAG + 14 * 2] = 512.3241541685
                cexpd[@@IMAG + 15 * 2] = 512.3304037599
                cexpd[@@IMAG + 16 * 2] = 512.3356167976
                cexpd[@@IMAG + 17 * 2] = 512.3399592211
                cexpd[@@IMAG + 18 * 2] = 512.3435588985
                cexpd[@@IMAG + 19 * 2] = 512.3465164008
            else
            end
          end
        end
      end
    end
    epsilon = 1.0E-12; 
    if nt <= 0 then
    else
        for it in 0..nt-1 do
            csumr = (cksum[@@REAL + it * 2] - cexpd[@@REAL + it * 2]) / cexpd[@@REAL + it * 2]; 
            csumi = (cksum[@@IMAG + it * 2] - cexpd[@@IMAG + it * 2]) / cexpd[@@IMAG + it * 2]; 
            if Math_dot_abs(csumr) <= epsilon or Math_dot_abs(csumi) <= epsilon then
                if verified == -1 then
                  verified = 1
                else
                end
            else
                verified = 0
            end
        end
    end
    BMResults.printVerificationStatus(@CLASS, verified, @@BMName)
    return verified
  end

  def getTime()
    return @timer.readTimer(1)
  end

  def isDone()
    return @done
  end

  def finalize()
    puts("FT: is about to be garbage collected")
    super.finalize()
  end

  def setupThreads(ft)
    @master = ft
    if @num_threads > @nz then
      @num_threads = @nz
    else
    end
    if @num_threads > @nx then
      @num_threads = @nx
    else
    end
    interval1 = Array.new(@num_threads, 0.0); 
    interval2 = Array.new(@num_threads, 0.0); 
    set_interval(@num_threads, @nz, interval1)
    set_interval(@num_threads, @nx, interval2)
    partition1 = MultiDimArray(interval1.length, 2); 
    partition2 = MultiDimArray(interval2.length, 2); 
    partition1 = MultiDimArray(interval1.length, 2)
    partition2 = MultiDimArray(interval2.length, 2)
    set_partition(0, interval1, partition1)
    set_partition(0, interval2, partition2)
    @doFFT = Array.new(@num_threads, 0.0)
    @doEvolve = Array.new(@num_threads, 0.0)
    for ii in 0..@num_threads-1 do
      @doFFT[ii] = FFTThread.new(ft, partition1[ii][0], partition1[ii][1], partition2[ii][0], partition2[ii][1])
      @doFFT[ii].extend(MonitorMixin)
      @doFFT[ii].id = ii
      @doFFT[ii].start()
      @doEvolve[ii] = EvolveThread.new(ft, partition2[ii][0], partition2[ii][1])
      @doEvolve[ii].extend(MonitorMixin)
      @doEvolve[ii].id = ii
      @doEvolve[ii].start()
    end
  end

# *** public ***

  attr_accessor :bid

  attr_accessor :results

  attr_accessor :serial

end

  def main(argv)
    ft =  nil ; 
    BMArgs.parseCmdLineArgs(argv, "FT")
    clss = BMArgs.clss; 
    np = BMArgs.num_threads; 
    serial = BMArgs.serial; 
    #begin
    ft = FT.new(clss, np, serial)
    ft.extend(MonitorMixin)
    #rescue OutOfMemoryError => e
    #    BMArgs.outOfMemoryMessage()
    #    System.exit(0)
    #ensure
    #end
    ft.runBenchMark()
  end

main(ARGV)
