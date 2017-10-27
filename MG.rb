# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
# Authors: E. Barszcz                                                     
#          P. Frederickson					         
#          A. Woo					                 
#          M. Yarrow					                 
# Translation to Java and MultiThreaded Code				  
#          M. Frumkin							 
#          M. Schultz
#  Translation to Ruby
#          T. Nose
#          H. Tomari						 
# -------------------------------------------------------------------------

require "mda"
require "mathutil"
require "MGThreads/MGBase"
require "MGThreads/Interp"
require "MGThreads/Psinv"
require "MGThreads/Resid"
require "MGThreads/Rprj"
require "BMInOut/BMResults"
require "BMInOut/BMArgs"
require "Timer"
require "Random"
require "Runnable"
require "monitor"

class MG < MGBase
  def initialize(clss, np, ser)
    @bid = -1
    @serial = false
    @timeron = false
    super(clss, np, ser)
    @serial = ser
  end

  def run()
    runBenchMark()
  end

  def runBenchMark()
    BMArgs.banner(@@BMName, @CLASS, @serial, @num_threads)
    niter = getInputPars(); 
    @nsizes = Array.new(3, 0.0)
    setup(@nsizes)
    @n1 = @nsizes[0]
    @n2 = @nsizes[1]
    @n3 = @nsizes[2]
    setTimers()
    @timer.resetAllTimers()
    @timer.start(@@T_init)
    zero3(@u, 0, @n1, @n2, @n3)
    zran3(@v, @n1, @n2, @n3, @nx[@lt - 1], @ny[@lt - 1])
    if not @serial then
      setupThreads( self )
    else
    end
    if @serial then
      resid(@u, @v, @r, 0, @n1, @n2, @n3)
    else
      residMaster(@u, @v, @r, 0, @n1, @n2, @n3)
    end
    if @serial then
        mg3P(@u, @v, @r, @n1, @n2, @n3)
        resid(@u, @v, @r, 0, @n1, @n2, @n3)
    else
        mg3Pmaster(@u, @v, @r, @n1, @n2, @n3)
        residMaster(@u, @v, @r, 0, @n1, @n2, @n3)
    end
    zero3(@u, 0, @n1, @n2, @n3)
    zran3(@v, @n1, @n2, @n3, @nx[@lt - 1], @ny[@lt - 1])
    @timer.stop(@@T_init)
    @timer.start(@@T_bench)
    if @timeron then
      @timer.start(@@T_resid2)
    else
    end
    if @serial then
      resid(@u, @v, @r, 0, @n1, @n2, @n3)
    else
      residMaster(@u, @v, @r, 0, @n1, @n2, @n3)
    end
    if @timeron then
      @timer.stop(@@T_resid2)
    else
    end
    for it in 1..@@nit do
        if @timeron then
          @timer.start(@@T_mg3P)
        else
        end
        if @serial then
          mg3P(@u, @v, @r, @n1, @n2, @n3)
        else
          mg3Pmaster(@u, @v, @r, @n1, @n2, @n3)
        end
        if @timeron then
          @timer.stop(@@T_mg3P)
        else
        end
        if @timeron then
          @timer.start(@@T_resid2)
        else
        end
        if @serial then
          resid(@u, @v, @r, 0, @n1, @n2, @n3)
        else
          residMaster(@u, @v, @r, 0, @n1, @n2, @n3)
        end
        if @timeron then
          @timer.stop(@@T_resid2)
        else
        end
    end
    @timer.stop(@@T_bench)
    tinit = @timer.readTimer(@@T_init); 
    puts(" Initialization time: " + tinit.to_s + " seconds")
    @rnm2 = norm2u3(@r, @n1, @n2, @n3, @rnmu, @nx[@lt - 1], @ny[@lt - 1], @nz[@lt - 1])
    @verified = verify(@rnm2)
    tm = @timer.readTimer(@@T_bench); 
    @results = BMResults.new("MG", @CLASS, @nx[@lt - 1], @ny[@lt - 1], @nz[@lt - 1], @@nit, tm, getMFLOPS(tm, @@nit), "floating point", @verified, @serial, @num_threads, @bid)
    @results.print()
    if @timeron then
      printTimers()
    else
    end
  end

  def verify(rnm2)
    verify_value = 0.0; 
    @epsilon = 1.0E-8
    if @CLASS != 'U' then
        if @CLASS == 'S' then
            verify_value = 0.530770700573E-4
        else
          if @CLASS == 'W' then
              verify_value = 0.250391406439E-17
          else
            if @CLASS == 'A' then
                verify_value = 0.2433365309E-5
            else
              if @CLASS == 'B' then
                  verify_value = 0.180056440132E-5
              else
                if @CLASS == 'C' then
                    verify_value = 0.570674826298E-6
                else
                end
              end
            end
          end
        end
        puts(" L2 Norm is " + rnm2.to_s)
        if Math_dot_abs(rnm2 - verify_value) < @epsilon then
            @verified = 1
            puts(" Deviation is   " + (rnm2 - verify_value).to_s)
        else
            @verified = 0
            puts(" The correct L2 Norm is " + verify_value.to_s)
        end
    else
        @verified = -1
    end
    BMResults.printVerificationStatus(@CLASS, @verified, @@BMName)
    return @verified
  end

  def getMFLOPS(tm, niter)
    mflops = 0.0; 
    if tm > 0.0 then
        mflops = 58.0 * @n1 * @n2 * @n3
        mflops *= niter / (tm * 1000000.0)
    else
    end
    return mflops
  end

  def getInputPars()
    lnx = 32; lny = 32; lnz = 32; 
    #f2 = File.new("mg.input"); 
    if File.exist?("mg.input") then
        puts("Reading from input file mg.input")
        begin
            fis = FileInputStream.new(f2); 
            datafile = DataInputStream.new(fis); 
            @lt = datafile.readInt()
            if @lt > @@maxlevel then
                puts("lt=" + @lt + " Maximum allowable=" + @@maxlevel)
                System.exit(0)
            else
            end
            lnx = datafile.readInt()
            lny = datafile.readInt()
            lnz = datafile.readInt()
            @@nit = datafile.readInt()
            fis.close()
        rescue Exception => e
            STDERR.puts("Error reading from file mg.input")
        ensure
        end
        if lnx != lny or lnx != lnz then
            @CLASS = 'U'
        else
          if lnx == 32 and @@nit == 4 then
              @CLASS = 'S'
          else
            if lnx == 64 and @@nit == 40 then
                @CLASS = 'W'
            else
              if lnx == 256 and @@nit == 20 then
                  @CLASS = 'B'
              else
                if lnx == 512 and @@nit == 20 then
                    @CLASS = 'C'
                else
                  if lnx == 256 and @@nit == 4 then
                      @CLASS = 'A'
                  else
                      @CLASS = 'U'
                  end
                end
              end
            end
          end
        end
    else
        puts(" No input file mg.input, Using compiled defaults")
    end
    puts(" Size:  " + @nx[@lt - 1].to_s + "x" + @ny[@lt - 1].to_s + "x" + @nz[@lt - 1].to_s + " Iterations:   " + @@nit.to_s)
    return @@nit
  end

  def setTimers()
    #f1 = File.new("timer.flag"); 
    if File.exist?("timer.flag") then
        @timeron = true
        @t_names = Array.new(16, 0.0)
        @t_names[@@T_init] = "init"
        @t_names[@@T_bench] = "benchmark"
        @t_names[@@T_mg3P] = "mg3P"
        @t_names[@@T_psinv] = "psinv"
        @t_names[@@T_resid] = "resid"
        @t_names[@@T_rprj3] = "rprj3"
        @t_names[@@T_interp] = "interp"
        @t_names[@@T_norm2] = "norm2"
    else
    end
  end

  def printTimers()
    #fmt = DecimalFormat.new("0.000"); 
    puts("  SECTION   Time (secs)")
    tmax = @timer.readTimer(@@T_bench); 
    if tmax == 0.0 then
      tmax = 1.0
    else
    end
    for i in @@T_bench..@@T_last do
        t = @timer.readTimer(i); 
        if i == @@T_resid2 then
            t = @timer.readTimer(@@T_resid) - t
            puts("	  --> total mg-resid " + (t).to_s + " (" + (t * 100.0 / tmax).to_s + "%)")
        else
            puts("	" + @t_names[i].to_s + "  " + (t).to_s + " (" + (t * 100.0 / tmax).to_s + "%)")
        end
    end
  end

  def setup(nsizes)
    
    
    
    size1 = 3; size2 = 10; 
    mi = Array.new(size1 * size2, 0.0); 
    ng = Array.new(size1 * size2, 0.0); 
    
    @lb = 1
    ng[(@lt - 1) * size1] = @nx[@lt - 1]
    ng[1 + (@lt - 1) * size1] = @ny[@lt - 1]
    ng[2 + (@lt - 1) * size1] = @nz[@lt - 1]
    for ax in 0..size1-1 do
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
      k = @lt - 2
      while k >= 0 do
        ng[ax + k * size1] = ng[ax + (k + 1) * size1] / 2
        k -= 1
      end
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    k = @lt - 2
    while k >= 0 do
        @nx[k] = ng[k * size1]
        @ny[k] = ng[1 + k * size1]
        @nz[k] = ng[2 + k * size1]
      k -= 1
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    k = @lt - 1
    while k >= 0 do
        for ax in 0..size1-1 do
            mi[ax + k * size1] = 2 + ng[ax + k * size1]
        end
        @m1[k] = mi[k * size1]
        @m2[k] = mi[1 + k * size1]
        @m3[k] = mi[2 + k * size1]
      k -= 1
    end
    k = @lt - 1
    @is1 = 2 + ng[k * size1] - ng[k * size1]
    @ie1 = 1 + ng[k * size1]
    @n1 = nsizes[0] = 3 + @ie1 - @is1
    @is2 = 2 + ng[1 + k * size1] - ng[1 + k * size1]
    @ie2 = 1 + ng[1 + k * size1]
    @n2 = nsizes[1] = 3 + @ie2 - @is2
    @is3 = 2 + ng[2 + k * size1] - ng[2 + k * size1]
    @ie3 = 1 + ng[2 + k * size1]
    @n3 = nsizes[2] = 3 + @ie3 - @is3
    @ir[@lt - 1] = 0
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    j = @lt - 2
    while j >= 0 do
        @ir[j] = @ir[j + 1] + @m1[j + 1] * @m2[j + 1] * @m3[j + 1]
      j -= 1
    end
  end

  def zero3(z, off, n1, n2, n3)
    
    for i3 in 0..n3-1 do
      for i2 in 0..n2-1 do
        for i1 in 0..n1-1 do
          z[off + i1 + n1 * (i2 + n2 * i3)] = 0.0
        end
      end
    end
  end

  def zran3(z, n1, n2, n3, nx, ny)
    
    mm = 10; 
    
    ten = Array.new(mm * 2, 0.0); 
    
    
    j1 = Array.new(mm * 2, 0.0); j2 = Array.new(mm * 2, 0.0); j3 = Array.new(mm * 2, 0.0); 
    jg = Array.new(4 * mm * 2, 0.0); jg_temp = Array.new(4, 0.0); 
    zero3(z, 0, n1, n2, n3)
    i = @is1 - 2 + nx * (@is2 - 2 + ny * (@is3 - 2))
    d1 = @ie1 - @is1 + 1
    e1 = @ie1 - @is1 + 2
    e2 = @ie2 - @is2 + 2
    e3 = @ie3 - @is3 + 2
    seed = 314159265.0; a = Math_dot_pow(5.0, 13); 
    rng = NPB3_0_RUB::Random.new(); 
    a1 = rng.power(a, nx)
    a2 = rng.power(a, nx * ny)
    ai = rng.power(a, i)
    x0 = rng.randlc(seed, ai)
    for i3 in 2..e3 do
        x1 = x0
        for i2 in 2..e2 do
            xx = x1
            rng.vranlc(d1, xx, a, z, (1 + n1 * (i2 - 1 + n2 * (i3 - 1))))
            x1 = rng.randlc(x1, a1)
        end
        x0 = rng.randlc(x0, a2)
    end
    for i in 0..mm-1 do
        ten[i + mm] = 0.0
        j1[i + mm] = 0
        j2[i + mm] = 0
        j3[i + mm] = 0
        ten[i] = 1.0
        j1[i] = 0
        j2[i] = 0
        j3[i] = 0
    end
    for i3 in 1..n3 - 1-1 do
        for i2 in 1..n2 - 1-1 do
            for i1 in 1..n1 - 1-1 do
                if z[i1 + n1 * (i2 + n2 * i3)] > ten[mm] then
                    ten[mm] = z[i1 + n1 * (i2 + n2 * i3)]
                    j1[mm] = i1
                    j2[mm] = i2
                    j3[mm] = i3
                    bubble(ten, j1, j2, j3, mm, 1)
                else
                end
                if z[i1 + n1 * (i2 + n2 * i3)] < ten[0] then
                    ten[0] = z[i1 + n1 * (i2 + n2 * i3)]
                    j1[0] = i1
                    j2[0] = i2
                    j3[0] = i3
                    bubble(ten, j1, j2, j3, mm, 0)
                else
                end
            end
        end
    end
    i1 = mm
    i0 = mm
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    i = mm - 1
    while i >= 0 do
        best = z[j1[i1 - 1 + mm] + n1 * (j2[i1 - 1 + mm] + n2 * (j3[i1 - 1 + mm]))]
        if best == z[j1[i1 - 1 + mm] + n1 * (j2[i1 - 1 + mm] + n2 * (j3[i1 - 1 + mm]))] then
            jg[4 * (i + mm)] = 0
            jg[1 + 4 * (i + mm)] = @is1 - 2 + j1[i1 - 1 + mm]
            jg[2 + 4 * (i + mm)] = @is2 - 2 + j2[i1 - 1 + mm]
            jg[3 + 4 * (i + mm)] = @is3 - 2 + j3[i1 - 1 + mm]
            i1 = i1 - 1
        else
            jg[4 * (i + mm)] = 0
            jg[1 + 4 * (i + mm)] = 0
            jg[2 + 4 * (i + mm)] = 0
            jg[3 + 4 * (i + mm)] = 0
        end
        ten[i + mm] = best
        best = z[j1[i0 - 1] + n1 * (j2[i0 - 1] + n2 * (j3[i0 - 1]))]
        if best == z[j1[i0 - 1] + n1 * (j2[i0 - 1] + n2 * (j3[i0 - 1]))] then
            jg[4 * i] = 0
            jg[1 + 4 * i] = @is1 - 2 + j1[i0 - 1]
            jg[2 + 4 * i] = @is2 - 2 + j2[i0 - 1]
            jg[3 + 4 * i] = @is3 - 2 + j3[i0 - 1]
            i0 = i0 - 1
        else
            jg[4 * i] = 0
            jg[1 + 4 * i] = 0
            jg[2 + 4 * i] = 0
            jg[3 + 4 * i] = 0
        end
        ten[i] = best
      i -= 1
    end
    m1 = i1 + 1
    m0 = i0 + 1
    for i3 in 0..n3-1 do
      for i2 in 0..n2-1 do
        for i1 in 0..n1-1 do
          z[i1 + n1 * (i2 + n2 * i3)] = 0.0
        end
      end
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    i = mm
    while i >= m0 do
      z[j1[i - 1] + n1 * (j2[i - 1] + n2 * (j3[i - 1]))] = -1.0
      i -= 1
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    i = mm
    while i >= m1 do
      z[j1[i - 1 + mm] + n1 * (j2[i - 1 + mm] + n2 * (j3[i - 1 + mm]))] = 1.0
      i -= 1
    end
    comm3(z, 0, n1, n2, n3)
  end

  def norm2u3(r, n1, n2, n3, rnmu, nx, ny, nz)
    if @timeron then
      @timer.start(@@T_norm2)
    else
    end
    rnmu = 0.0
    rnm2 = 0.0; 
    for i3 in 1..n3 - 1-1 do
      for i2 in 1..n2 - 1-1 do
        for i1 in 1..n1 - 1-1 do
            rnm2 += r[i1 + n1 * (i2 + n2 * i3)] * r[i1 + n1 * (i2 + n2 * i3)]
            a = Math_dot_abs(r[i1 + n1 * (i2 + n2 * i3)]); 
            rnmu = dmax1(rnmu, a)
        end
      end
    end
    rnm2 = Math.sqrt(rnm2 / (nx * ny * nz))
    if @timeron then
      @timer.stop(@@T_norm2)
    else
    end
    return rnm2
  end

  def TestNorm(r, n1, n2, n3)
    rnm2 = 0.0; 
    for i3 in 1..n3 - 1-1 do
      for i2 in 1..n2 - 1-1 do
        for i1 in 1..n1 - 1-1 do
            rnm2 += r[i1 + n1 * (i2 + n2 * i3)] * r[i1 + n1 * (i2 + n2 * i3)]
        end
      end
    end
    rnm2 = Math.sqrt(rnm2 / (n1 * n2 * n3))
    puts("*****TestNorm  " + rnm2)
    return rnm2
  end

  def bubble(ten, j1, j2, j3, m, ind)
    
    j_temp = 0; 
    if ind == 1 then
        for i in 0..m - 1-1 do
            if ten[i + m * ind] > ten[i + 1 + m * ind] then
                temp = ten[i + 1 + m * ind]
                ten[i + 1 + m * ind] = ten[i + m * ind]
                ten[i + m * ind] = temp
                j_temp = j1[i + 1 + m * ind]
                j1[i + 1 + m * ind] = j1[i + m * ind]
                j1[i + m * ind] = j_temp
                j_temp = j2[i + 1 + m * ind]
                j2[i + 1 + m * ind] = j2[i + m * ind]
                j2[i + m * ind] = j_temp
                j_temp = j3[i + 1 + m * ind]
                j3[i + 1 + m * ind] = j3[i + m * ind]
                j3[i + m * ind] = j_temp
            else
                return 
            end
        end
    else
        for i in 0..m - 1-1 do
            if ten[i + m * ind] < ten[i + 1 + m * ind] then
                temp = ten[i + 1 + m * ind]
                ten[i + 1 + m * ind] = ten[i + m * ind]
                ten[i + m * ind] = temp
                j_temp = j1[i + 1 + m * ind]
                j1[i + 1 + m * ind] = j1[i + m * ind]
                j1[i + m * ind] = j_temp
                j_temp = j2[i + 1 + m * ind]
                j2[i + 1 + m * ind] = j2[i + m * ind]
                j2[i + m * ind] = j_temp
                j_temp = j3[i + 1 + m * ind]
                j3[i + 1 + m * ind] = j3[i + m * ind]
                j3[i + m * ind] = j_temp
            else
                return 
            end
        end
    end
  end

  def resid(u, v, r, off, n1, n2, n3)
    
    u1 = Array.new(@nm + 1, 0.0); 
    u2 = Array.new(@nm + 1, 0.0); 
    if @timeron then
      @timer.start(@@T_resid)
    else
    end
    for i3 in 1..n3 - 1-1 do
      for i2 in 1..n2 - 1-1 do
          for i1 in 0..n1-1 do
              u1[i1] = u[off + i1 + n1 * (i2 - 1 + n3 * i3)] + u[off + i1 + n1 * (i2 + 1 + n3 * i3)] + u[off + i1 + n1 * (i2 + n3 * (i3 - 1))] + u[off + i1 + n1 * (i2 + n3 * (i3 + 1))]
              u2[i1] = u[off + i1 + n1 * (i2 - 1 + n3 * (i3 - 1))] + u[off + i1 + n1 * (i2 + 1 + n3 * (i3 - 1))] + u[off + i1 + n1 * (i2 - 1 + n3 * (i3 + 1))] + u[off + i1 + n1 * (i2 + 1 + n3 * (i3 + 1))]
          end
          for i1 in 1..n1 - 1-1 do
              r[off + i1 + n1 * (i2 + n3 * i3)] = v[off + i1 + n1 * (i2 + n3 * i3)] - @a[0] * u[off + i1 + n1 * (i2 + n3 * i3)] - @a[2] * (u2[i1] + u1[i1 - 1] + u1[i1 + 1]) - @a[3] * (u2[i1 - 1] + u2[i1 + 1])
          end
      end
    end
    comm3(r, off, n1, n2, n3)
    if @timeron then
      @timer.stop(@@T_resid)
    else
    end
  end

  def mg3P(u, v, r, n1, n2, n3)
    
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    k = @lt - 1
    while k >= @lb do
        j = k - 1
        rprj3(r, @ir[k], @m1[k], @m2[k], @m3[k], @ir[j], @m1[j], @m2[j], @m3[j])
      k -= 1
    end
    k = @lb - 1
    zero3(u, @ir[k], @m1[k], @m2[k], @m3[k])
    psinv(r, @ir[k], u, @ir[k], @m1[k], @m2[k], @m3[k])
    for k in @lb..@lt - 1-1 do
        j = k - 1
        zero3(u, @ir[k], @m1[k], @m2[k], @m3[k])
        interp(u, @ir[j], @m1[j], @m2[j], @m3[j], @ir[k], @m1[k], @m2[k], @m3[k])
        resid(u, r, r, @ir[k], @m1[k], @m2[k], @m3[k])
        psinv(r, @ir[k], u, @ir[k], @m1[k], @m2[k], @m3[k])
    end
    j = @lt - 2
    k = @lt - 1
    interp(u, @ir[j], @m1[j], @m2[j], @m3[j], 0, n1, n2, n3)
    resid(u, v, r, 0, n1, n2, n3)
    psinv(r, 0, u, 0, n1, n2, n3)
  end

  def mg3Pmaster(u, v, r, n1, n2, n3)
    
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    k = @lt - 1
    while k >= @lb do
        j = k - 1
        rprj3Master(r, @ir[k], @m1[k], @m2[k], @m3[k], @ir[j], @m1[j], @m2[j], @m3[j])
      k -= 1
    end
    k = @lb - 1
    zero3(u, @ir[k], @m1[k], @m2[k], @m3[k])
    psinvMaster(r, @ir[k], u, @ir[k], @m1[k], @m2[k], @m3[k])
    for k in @lb..@lt - 1-1 do
        j = k - 1
        zero3(u, @ir[k], @m1[k], @m2[k], @m3[k])
        interpMaster(u, @ir[j], @m1[j], @m2[j], @m3[j], @ir[k], @m1[k], @m2[k], @m3[k])
        residMaster(u, r, r, @ir[k], @m1[k], @m2[k], @m3[k])
        psinvMaster(r, @ir[k], u, @ir[k], @m1[k], @m2[k], @m3[k])
    end
    j = @lt - 2
    k = @lt - 1
    interpMaster(u, @ir[j], @m1[j], @m2[j], @m3[j], 0, n1, n2, n3)
    residMaster(u, v, r, 0, n1, n2, n3)
    psinvMaster(r, 0, u, 0, n1, n2, n3)
  end

  def rprj3(r, roff, m1k, m2k, m3k, soff, m1j, m2j, m3j)
    
    
    x1 = Array.new(@nm + 1, 0.0); y1 = Array.new(@nm + 1, 0.0); 
    if @timeron then
      @timer.start(@@T_rprj3)
    else
    end
    if m1k == 3 then
        d1 = 2
    else
        d1 = 1
    end
    if m2k == 3 then
        d2 = 2
    else
        d2 = 1
    end
    if m3k == 3 then
        d3 = 2
    else
        d3 = 1
    end
    for j3 in 2..m3j - 1 do
        i3 = 2 * j3 - d3 - 1
        for j2 in 2..m2j - 1 do
            i2 = 2 * j2 - d2 - 1
            for j1 in 2..m1j do
                i1 = 2 * j1 - d1 - 1
                x1[i1 - 1] = r[roff + i1 - 1 + m1k * (i2 - 1 + m2k * i3)] + r[roff + i1 - 1 + m1k * (i2 + 1 + m2k * i3)] + r[roff + i1 - 1 + m1k * (i2 + m2k * (i3 - 1))] + r[roff + i1 - 1 + m1k * (i2 + m2k * (i3 + 1))]
                y1[i1 - 1] = r[roff + i1 - 1 + m1k * (i2 - 1 + m2k * (i3 - 1))] + r[roff + i1 - 1 + m1k * (i2 - 1 + m2k * (i3 + 1))] + r[roff + i1 - 1 + m1k * (i2 + 1 + m2k * (i3 - 1))] + r[roff + i1 - 1 + m1k * (i2 + 1 + m2k * (i3 + 1))]
            end
            for j1 in 2..m1j - 1 do
                i1 = 2 * j1 - d1 - 1
                y2 = r[roff + i1 + m1k * (i2 - 1 + m2k * (i3 - 1))] + r[roff + i1 + m1k * (i2 - 1 + m2k * (i3 + 1))] + r[roff + i1 + m1k * (i2 + 1 + m2k * (i3 - 1))] + r[roff + i1 + m1k * (i2 + 1 + m2k * (i3 + 1))]
                x2 = r[roff + i1 + m1k * (i2 - 1 + m2k * i3)] + r[roff + i1 + m1k * (i2 + 1 + m2k * i3)] + r[roff + i1 + m1k * (i2 + m2k * (i3 - 1))] + r[roff + i1 + m1k * (i2 + m2k * (i3 + 1))]
                r[soff + j1 - 1 + m1j * (j2 - 1 + m2j * (j3 - 1))] = 0.5 * r[roff + i1 + m1k * (i2 + m2k * i3)] + 0.25 * (r[roff + i1 - 1 + m1k * (i2 + m2k * i3)] + r[roff + i1 + 1 + m1k * (i2 + m2k * i3)] + x2) + 0.125 * (x1[i1 - 1] + x1[i1 + 1] + y2) + 0.0625 * (y1[i1 - 1] + y1[i1 + 1])
            end
        end
    end
    comm3(r, soff, m1j, m2j, m3j)
    if @timeron then
      @timer.stop(@@T_rprj3)
    else
    end
  end

  def interp(u, zoff, mm1, mm2, mm3, uoff, n1, n2, n3)
    
    m = 535; 
    z1 = Array.new(m, 0.0); z2 = Array.new(m, 0.0); z3 = Array.new(m, 0.0); 
    if @timeron then
      @timer.start(@@T_interp)
    else
    end
    if n1 != 3 and n2 != 3 and n3 != 3 then
        for i3 in 1..mm3 - 1 do
            for i2 in 1..mm2 - 1 do
                for i1 in 1..mm1 do
                    z1[i1 - 1] = u[zoff + i1 - 1 + mm1 * (i2 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))]
                    z2[i1 - 1] = u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * i3)] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))]
                    z3[i1 - 1] = u[zoff + i1 - 1 + mm1 * (i2 + mm2 * i3)] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * i3)] + z1[i1 - 1]
                end
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 2 + n1 * (2 * i2 - 2 + n2 * (2 * i3 - 2))] += u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))]
                    u[uoff + 2 * i1 - 1 + n1 * (2 * i2 - 2 + n2 * (2 * i3 - 2))] += 0.5 * (u[zoff + i1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))])
                end
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 2 + n1 * (2 * i2 - 1 + n2 * (2 * i3 - 2))] += 0.5 * z1[i1 - 1]
                    u[uoff + 2 * i1 - 1 + n1 * (2 * i2 - 1 + n2 * (2 * i3 - 2))] += 0.25 * (z1[i1 - 1] + z1[i1])
                end
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 2 + n1 * (2 * i2 - 2 + n2 * (2 * i3 - 1))] += 0.5 * z2[i1 - 1]
                    u[uoff + 2 * i1 - 1 + n1 * (2 * i2 - 2 + n2 * (2 * i3 - 1))] += 0.25 * (z2[i1 - 1] + z2[i1])
                end
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 2 + n1 * (2 * i2 - 1 + n2 * (2 * i3 - 1))] += 0.25 * z3[i1 - 1]
                    u[uoff + 2 * i1 - 1 + n1 * (2 * i2 - 1 + n2 * (2 * i3 - 1))] += 0.125 * (z3[i1 - 1] + z3[i1])
                end
            end
        end
    else
        if n1 == 3 then
            d1 = 2
            t1 = 1
        else
            d1 = 1
            t1 = 0
        end
        if n2 == 3 then
            d2 = 2
            t2 = 1
        else
            d2 = 1
            t2 = 0
        end
        if n3 == 3 then
            d3 = 2
            t3 = 1
        else
            d3 = 1
            t3 = 0
        end
        for i3 in 1..mm3 - 1 do
            for i2 in 1..mm2 - 1 do
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 1 - d1 + n1 * (2 * i2 - 1 - d2 + n2 * (2 * i3 - 1 - d3))] += u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))]
                end
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 1 - t1 + n1 * (2 * i2 - 1 - d2 + n2 * (2 * i3 - 1 - d3))] += 0.5 * (u[zoff + i1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))])
                end
            end
            for i2 in 1..mm2 - 1 do
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 1 - d1 + n1 * (2 * i2 - 1 - t2 + n2 * (2 * i3 - 1 - d3))] += 0.5 * (u[zoff + i1 - 1 + mm1 * (i2 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))])
                end
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 1 - t1 + n1 * (2 * i2 - 1 - t2 + n2 * (2 * i3 - 1 - d3))] += 0.25 * (u[zoff + i1 + mm1 * (i2 + mm2 * (i3 - 1))] + u[zoff + i1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))])
                end
            end
        end
        for i3 in 1..mm3 - 1 do
            for i2 in 1..mm2 - 1 do
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 1 - d1 + n1 * (2 * i2 - 1 - d2 + n2 * (2 * i3 - 1 - t3))] = 0.5 * (u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * i3)] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))])
                end
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 1 - t1 + n1 * (2 * i2 - 1 - d2 + n2 * (2 * i3 - 1 - t3))] += 0.25 * (u[zoff + i1 + mm1 * (i2 - 1 + mm2 * i3)] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * i3)] + u[zoff + i1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))])
                end
            end
            for i2 in 1..mm2 - 1 do
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 1 - d1 + n1 * (2 * i2 - 1 - t2 + n2 * (2 * i3 - 1 - t3))] += 0.25 * (u[zoff + i1 - 1 + mm1 * (i2 + mm2 * i3)] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * i3)] + u[zoff + i1 - 1 + mm1 * (i2 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))])
                end
                for i1 in 1..mm1 - 1 do
                    u[uoff + 2 * i1 - 1 - t1 + n1 * (2 * i2 - 1 - t2 + n2 * (2 * i3 - 1 - t3))] += 0.125 * (u[zoff + i1 + mm1 * (i2 + mm2 * i3)] + u[zoff + i1 + mm1 * (i2 - 1 + mm2 * i3)] + u[zoff + i1 - 1 + mm1 * (i2 + mm2 * i3)] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * i3)] + u[zoff + i1 + mm1 * (i2 + mm2 * (i3 - 1))] + u[zoff + i1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 + mm2 * (i3 - 1))] + u[zoff + i1 - 1 + mm1 * (i2 - 1 + mm2 * (i3 - 1))])
                end
            end
        end
    end
    if @timeron then
      @timer.stop(@@T_interp)
    else
    end
  end

  def psinv(r, roff, u, uoff, n1, n2, n3)
    
    r1 = Array.new(@nm + 1, 0.0); r2 = Array.new(@nm + 1, 0.0); 
    if @timeron then
      @timer.start(@@T_psinv)
    else
    end
    for i3 in 1..n3 - 1-1 do
        for i2 in 1..n2 - 1-1 do
            for i1 in 0..n1-1 do
                r1[i1] = r[roff + i1 + n1 * (i2 - 1 + n2 * i3)] + r[roff + i1 + n1 * (i2 + 1 + n2 * i3)] + r[roff + i1 + n1 * (i2 + n2 * (i3 - 1))] + r[roff + i1 + n1 * (i2 + n2 * (i3 + 1))]
                r2[i1] = r[roff + i1 + n1 * (i2 - 1 + n2 * (i3 - 1))] + r[roff + i1 + n1 * (i2 + 1 + n2 * (i3 - 1))] + r[roff + i1 + n1 * (i2 - 1 + n2 * (i3 + 1))] + r[roff + i1 + n1 * (i2 + 1 + n2 * (i3 + 1))]
            end
            for i1 in 1..n1 - 1-1 do
                u[uoff + i1 + n1 * (i2 + n2 * i3)] += @c[0] * r[roff + i1 + n1 * (i2 + n2 * i3)] + @c[1] * (r[roff + i1 - 1 + n1 * (i2 + n2 * i3)] + r[roff + i1 + 1 + n1 * (i2 + n2 * i3)] + r1[i1]) + @c[2] * (r2[i1] + r1[i1 - 1] + r1[i1 + 1])
            end
        end
    end
    comm3(u, uoff, n1, n2, n3)
    if @timeron then
      @timer.stop(@@T_psinv)
    else
    end
  end

  def residMaster(u, v, r, off, n1, n2, n3)
    if @timeron then
      @timer.start(@@T_resid)
    else
    end
    if @num_threads == 1 then
      resid(u, v, r, off, n1, n2, n3)
    else
        visr = false; 
        if v == r then
          visr = true
        else
        end
         self.synchronize do
            for l in 0..@num_threads-1 do
              @resid[l].synchronize do
                  @resid[l].done = false
                  @resid[l].visr = visr
                  @resid[l].wstart = 1
                  @resid[l].wend = n3
                  @resid[l].n1 = n1
                  @resid[l].n2 = n2
                  @resid[l].n3 = n3
                  @resid[l].off = off
                  @resid[l].cond.signal
              end
            end
            for l in 0..@num_threads-1 do
              while not @resid[l].done do
                  begin
                      cond.wait
                  rescue InterruptedException => e
                  ensure
                  end
                  cond.broadcast
              end
            end
        end
        comm3(r, off, n1, n2, n3)
    end
    if @timeron then
      @timer.stop(@@T_resid)
    else
    end
  end

  def psinvMaster(r, roffl, u, uoffl, n1, n2, n3)
    if @timeron then
      @timer.start(@@T_psinv)
    else
    end
    if @num_threads == 1 then
      psinv(r, roffl, u, uoffl, n1, n2, n3)
    else
         self.synchronize do
            for l in 0..@num_threads-1 do
              @psinv[l].synchronize do
                  @psinv[l].done = false
                  @psinv[l].wstart = 1
                  @psinv[l].wend = n3
                  @psinv[l].n1 = n1
                  @psinv[l].n2 = n2
                  @psinv[l].n3 = n3
                  @psinv[l].roff = roffl
                  @psinv[l].uoff = uoffl
                  @psinv[l].cond.signal
              end
            end
            for l in 0..@num_threads-1 do
              while not @psinv[l].done do
                  begin
                      cond.wait
                  rescue InterruptedException => e
                  ensure
                  end
                  cond.broadcast
              end
            end
        end
        comm3(u, uoffl, n1, n2, n3)
    end
    if @timeron then
      @timer.stop(@@T_psinv)
    else
    end
  end

  def interpMaster(u, zoffl, mm1, mm2, mm3, uoffl, n1, n2, n3)
    if @timeron then
      @timer.start(@@T_interp)
    else
    end
    if @num_threads == 1 then
      interp(u, zoffl, mm1, mm2, mm3, uoffl, n1, n2, n3)
    else
         self.synchronize do
            for l in 0..@num_threads-1 do
              @interp[l].synchronize do
                  @interp[l].done = false
                  @interp[l].wstart = 1
                  @interp[l].wend = mm3
                  @interp[l].mm1 = mm1
                  @interp[l].mm2 = mm2
                  @interp[l].mm3 = mm3
                  @interp[l].n1 = n1
                  @interp[l].n2 = n2
                  @interp[l].n3 = n3
                  @interp[l].zoff = zoffl
                  @interp[l].uoff = uoffl
                  @interp[l].cond.signal
              end
            end
            for l in 0..@num_threads-1 do
              while not @interp[l].done do
                  begin
                      cond.wait
                  rescue InterruptedException => e
                  ensure
                  end
                  cond.broadcast
              end
            end
        end
    end
    if @timeron then
      @timer.stop(@@T_interp)
    else
    end
  end

  def rprj3Master(r, roffl, m1k, m2k, m3k, soffl, m1j, m2j, m3j)
    if @timeron then
      @timer.start(@@T_rprj3)
    else
    end
    if @num_threads == 1 then
      rprj3(r, roffl, m1k, m2k, m3k, soffl, m1j, m2j, m3j)
    else
         self .synchronize do
            for l in 0..@num_threads-1 do
              @rprj[l].synchronize do
                  @rprj[l].done = false
                  @rprj[l].wstart = 2
                  @rprj[l].wend = m3j
                  @rprj[l].m1k = m1k
                  @rprj[l].m2k = m2k
                  @rprj[l].m3k = m3k
                  @rprj[l].m1j = m1j
                  @rprj[l].m2j = m2j
                  @rprj[l].m3j = m3j
                  @rprj[l].roff = roffl
                  @rprj[l].zoff = soffl
                  @rprj[l].cond.signal
              end
            end
            for l in 0..@num_threads-1 do
                while not @rprj[l].done do
                    begin
                        cond.wait
                    rescue InterruptedException => e
                    ensure
                    end
                    cond.broadcast
                end
            end
        end
        comm3(r, soffl, m1j, m2j, m3j)
    end
    if @timeron then
      @timer.stop(@@T_rprj3)
    else
    end
  end

  def getTime()
    return @timer.readTimer(@@T_bench)
  end

  def finalize()
    puts("MG: is about to be garbage collected")
    super.finalize()
  end

# *** public ***

  attr_accessor :bid

  attr_accessor :results

  attr_accessor :serial

  attr_accessor :timeron

  attr_accessor :rnm2

  attr_accessor :rnmu

  attr_accessor :epsilon

  attr_accessor :n1

  attr_accessor :n2

  attr_accessor :n3

  attr_accessor :nn

  def setupThreads(mg)
    @interp = Array.new(@num_threads, 0.0)
    @psinv = Array.new(@num_threads, 0.0)
    @rprj = Array.new(@num_threads, 0.0)
    @resid = Array.new(@num_threads, 0.0)
    for i in 0..@num_threads-1 do
      @interp[i] = Interp.new(mg)
      @interp[i].extend(MonitorMixin)
      @interp[i].id = i
      @interp[i].start()

      @psinv[i] = Psinv.new(mg)
      @psinv[i].extend(MonitorMixin)
      @psinv[i].id = i
      @psinv[i].start()

      @rprj[i] = Rprj.new(mg)
      @rprj[i].extend(MonitorMixin)
      @rprj[i].id = i
      @rprj[i].start()

      @resid[i] = Resid.new(mg)
      @resid[i].extend(MonitorMixin)
      @resid[i].id = i
      @resid[i].start()
    end
  end
end

  def main(argv)
    mg =  nil ; 
    BMArgs.parseCmdLineArgs(argv, "MG")
    clss = BMArgs.clss; 
    np = BMArgs.num_threads; 
    serial = BMArgs.serial; 
    #begin
    mg = MG.new(clss, np, serial)
    mg.extend(MonitorMixin)
    #rescue OutOfMemoryError => e
    #BMArgs.outOfMemoryMessage()
    #exit(0)
    #ensure
    #end
    mg.runBenchMark()
  end

main(ARGV)
