# -------------------------------------------------------------------------
#  NAS Parallel Benchmarks 3.0 Ruby version
# -------------------------------------------------------------------------
#  Authors: M. Yarrow
#           C. Kuszmaul
#  Translation to Java and to MultiThreaded Code
# 	    M. Frumkin
# 	    M. Schultz 
#  Translation to Ruby
#           T. Nose
#           H. Tomari
# -------------------------------------------------------------------------

require "mda"
require "mathutil"
require "CGThreads/CGBase"
require "CGThreads/CGWorker"
require "Runnable"
require "BMInOut/BMResults"
require "BMInOut/BMArgs"
require "Random"
require "Timer"
require "monitor"
class CG < CGBase
  def initialize(clss, np, ser)
    @bid = -1
    @serial = true
    @amult = 1220703125.0
    super(clss, np, ser)
    @serial = ser
    @rng = NPB3_0_RUB::Random.new()
  end

  def run()
    runBenchMark()
  end

  def runBenchMark()
    
    
    rnorm = 0; 
    BMArgs.banner(@BMName, @CLASS, @serial, @num_threads)
    puts(" Size: " + @na.to_s + " Iterations: " + @niter.to_s)
    @timer.resetAllTimers()
    setTimers()
    @timer.start(@t_init)
    @firstrow = 1
    @lastrow = @na
    @firstcol = 1
    @lastcol = @na
    @naa = @na
    @nzz = @nz
    zeta = @rng.randlc(@amult)
    if not @serial then
      setupThreads( self )
    else
    end
    makea(@naa, @nzz, @a, @colidx, @rowstr, @nonzer, @rcond, @arow, @acol, @aelt, @v, @iv, @shift)
    for j in 1..@lastrow - @firstrow + 1 do
        for k in @rowstr[j]..@rowstr[j + 1] - 1 do
            @colidx[k] = @colidx[k] - @firstcol + 1
        end
    end
    for i in 1..@na + 1 do
      @x[i] = 1.0
    end
    zeta = 0.0
    rnorm = conj_grad(@colidx, @rowstr, @x, @z, @a, @p, @q, @r, rnorm)
    tnorm1 = 0.0; 
    tnorm2 = 0.0; 
    for j in 1..@lastcol - @firstcol + 1 do
        tnorm1 += @x[j] * @z[j]
        tnorm2 += @z[j] * @z[j]
    end
    tnorm2 = 1.0 / Math.sqrt(tnorm2)
    for j in 1..@lastcol - @firstcol + 1 do
        @x[j] = tnorm1 * @z[j]
    end
    for i in 1..@na + 1 do
      @x[i] = 1.0
    end
    zeta = 0.0
    @timer.stop(@t_init)
    @timer.start(@t_bench)
    for it in 1..@niter do
        if @timeron then
          @timer.start(@t_conj_grad)
        else
        end
        if @serial then
            rnorm = conj_grad(@colidx, @rowstr, @x, @z, @a, @p, @q, @r, rnorm)
        else
            executeTask(3)
            rho = 0.0; 
            for m in 0..@num_threads-1 do
              rho += @rhomaster[m]
            end
            for ii in 0..@cgitmax-1 do
                executeTask(0)
                dcff = 0.0; 
                for m in 0..@num_threads-1 do
                  dcff += @dmaster[m]
                end
                @alpha = rho / dcff
                rho0 = rho; 
                executeTask(1)
                rho = 0.0
                for m in 0..@num_threads-1 do
                  rho += @rhomaster[m]
                end
                @beta = rho / rho0
                executeTask(2)
            end
            executeTask(4)
            rnorm = 0.0
            for m in 0..@num_threads-1 do
              rnorm += @rnormmaster[m]
            end
            rnorm = Math.sqrt(rnorm)
        end
        if @timeron then
          @timer.stop(@t_conj_grad)
        else
        end
        tnorm1 = 0.0
        tnorm2 = 0.0
        for j in 1..@lastcol - @firstcol + 1 do
            tnorm1 += @x[j] * @z[j]
            tnorm2 += @z[j] * @z[j]
        end
        tnorm2 = 1.0 / Math.sqrt(tnorm2)
        zeta = @shift + 1.0 / tnorm1
        puts("    " + it.to_s + "       " + rnorm.to_s + " " + zeta.to_s)
        for j in 1..@lastcol - @firstcol + 1 do
            @x[j] = tnorm2 * @z[j]
        end
    end
    @timer.stop(@t_bench)
    verified = verify(zeta); 
    time = @timer.readTimer(@t_bench); 
    @results = BMResults.new(@BMName, @CLASS, @na, 0, 0, @niter, time, getMFLOPS(time, @niter), "floating point", verified, @serial, @num_threads, @bid)
    @results.print()
    if @timeron then
      PrintTimers()
    else
    end
  end

  def setTimers()
    #fp = File.new("timer.flag"); 
    if File.exist?("timer.flag") then
        @timeron = true
        @t_names[@t_init] = ("init")
        @t_names[@t_bench] = ("benchmark")
        @t_names[@t_conj_grad] = ("conjugate gradient")
    else
        @timeron = false
    end
  end

  def getMFLOPS(total_time, niter)
    mflops = 0.0; 
    if total_time != 0.0 then
        mflops = (2 * @na) * (3.0 + (@nonzer * (@nonzer + 1)) + 25.0 * (5.0 + (@nonzer * (@nonzer + 1))) + 3.0)
        mflops *= niter / (total_time * 1000000.0)
    else
    end
    return mflops
  end

  def verify(zeta)
    verified = 0; 
    epsilon = 1.0E-10; 
    if @CLASS != 'U' then
        puts(" Zeta is   " + zeta.to_s)
        if Math_dot_abs(zeta - @zeta_verify_value) <= epsilon then
            verified = 1
            puts(" Deviation is   " + (zeta - @zeta_verify_value).to_s)
        else
            verified = 0
            puts(" The correct zeta is " + @zeta_verify_value.to_s)
        end
    else
        verified = -1
    end
    BMResults.printVerificationStatus(@CLASS, verified, @BMName)
    return verified
  end

  def makea(n, nz, a, colidx, rowstr, nonzer, rcond, arow, acol, aelt, v, iv, shift)
    
    
    size = 1.0
    ratio = Math_dot_pow(rcond, (1.0 / n))
    nnza = 0
    for i in 1..n do
        colidx[n + i] = 0
    end
    for iouter in 1..n do
        nzv = nonzer
        sprnvc(n, nzv, v, iv, colidx, 0, colidx, n)
        nzv = vecset(n, v, iv, nzv, iouter, 0.5)
        for ivelt in 1..nzv do
            jcol = iv[ivelt]
            if jcol >= @firstcol and jcol <= @lastcol then
                scale = size * v[ivelt]
                for ivelt1 in 1..nzv do
                    irow = iv[ivelt1]
                    if irow >= @firstrow and irow <= @lastrow then
                        nnza = nnza + 1
                        if nnza > nz then
                            puts("Space for matrix elements exceeded in makea")
                            puts("nnza, nzmax = " + nnza.to_s + ", " + nz.to_s)
                            puts(" iouter = " + iouter.to_s)
                            exit(0)
                        else
                        end
                        acol[nnza] = jcol
                        arow[nnza] = irow
                        aelt[nnza] = v[ivelt1] * scale
                    else
                    end
                end
            else
            end
        end
        size = size * ratio
    end
    for i in @firstrow..@lastrow do
        if i >= @firstcol and i <= @lastcol then
            iouter = n + i; 
            nnza = nnza + 1
            if nnza > nz then
                puts("Space for matrix elements exceeded in makea")
                puts("nnza, nzmax = " + nnza.to_s + ", " + nz.to_s)
                puts(" iouter = " + iouter.to_s)
                exit(0)
            else
            end
            acol[nnza] = i
            arow[nnza] = i
            aelt[nnza] = rcond - shift
        else
        end
    end
    sparse(a, colidx, rowstr, n, arow, acol, aelt, v, iv, 0, iv, n, nnza)
    return 
  end

  def sprnvc(n, nz, v, iv, nzloc, nzloc_offst, mark, mark_offst)
    nzrow = 0; nzv = 0; 
    nn1 = 1; 
    while true do
        nn1 = 2 * nn1
        if nn1 >= n then
          break
        else
        end
    end
    while true do
        if nzv >= nz then
            for ii in 1..nzrow do
                idx = nzloc[ii + nzloc_offst]
                mark[idx + mark_offst] = 0
            end
            return 
        else
        end
        vecelt = @rng.randlc(@amult); 
        vecloc = @rng.randlc(@amult); 
        idx = ((vecloc * nn1) + 1).to_i
        if idx > n then
          next
        else
        end
        if mark[idx + mark_offst] == 0 then
            mark[idx + mark_offst] = 1
            nzrow = nzrow + 1
            nzloc[nzrow + nzloc_offst] = idx
            nzv = nzv + 1
            v[nzv] = vecelt
            iv[nzv] = idx
        else
        end
    end
  end

  def vecset(n, v, iv, nzv, ival, val)
    set = false; 
    for k in 1..nzv do
        if iv[k] == ival then
            v[k] = val
            set = true
        else
        end
    end
    if not set then
        nzv = nzv + 1
        v[nzv] = val
        iv[nzv] = ival
    else
    end
    return nzv
  end

  def sparse(a, colidx, rowstr, n, arow, acol, aelt, x, mark, mark_offst, nzloc, nzloc_offst, nnza)
    
    
    
    nrows = @lastrow - @firstrow + 1
    for j in 1..n do
        rowstr[j] = 0
        mark[j + mark_offst] = 0
    end
    rowstr[n + 1] = 0
    for nza in 1..nnza do
        j = (arow[nza] - @firstrow + 1) + 1
        rowstr[j] = rowstr[j] + 1
    end
    rowstr[1] = 1
    for j in 2..nrows + 1 do
        rowstr[j] = rowstr[j] + rowstr[j - 1]
    end
    for nza in 1..nnza do
        j = arow[nza] - @firstrow + 1
        k = rowstr[j]
        a[k] = aelt[nza]
        colidx[k] = acol[nza]
        rowstr[j] = rowstr[j] + 1
    end
# this loop cannot be expressed as a for-loop in Ruby.
# translated into while-do loop.
    j = nrows - 1
    while j >= 0 do
        rowstr[j + 1] = rowstr[j]
      j -= 1
    end
    rowstr[1] = 1
    nza = 0
    for i in 1..n do
        x[i] = 0.0
        mark[i + mark_offst] = 0
    end
    jajp1 = rowstr[1]
    for j in 1..nrows do
        nzrow = 0
        for k in jajp1..(rowstr[j + 1] - 1) do
            i = colidx[k]
            if false then
                puts("k:" + k.to_s)
                puts("x[i]:" + (x[i]).to_s)
                puts("a[k]:" + (a[k]).to_s)
                puts("i:" + i.to_s)
                puts("mark_offst:" + mark_offst.to_s)
                exit(-1)
            else
            end
            x[i] = x[i] + a[k]
            if (mark[i + mark_offst] == 0) and (x[i] != 0) then
                mark[i + mark_offst] = 1
                nzrow = nzrow + 1
                nzloc[nzrow + nzloc_offst] = i
            else
            end
        end
        for k in 1..nzrow do
            i = nzloc[k + nzloc_offst]
            mark[i + mark_offst] = 0
            xi = x[i]
            x[i] = 0
            if xi != 0 then
                nza = nza + 1
                a[nza] = xi
                colidx[nza] = i
            else
            end
        end
        jajp1 = rowstr[j + 1]
        rowstr[j + 1] = nza + rowstr[1]
    end
    return 
  end

  def conj_grad(colidx, rowstr, x, z, a, p, q, r, rnorm)
    
    
    
    for j in 1..@naa + 1 do
        q[j] = 0.0
        z[j] = 0.0
        r[j] = x[j]
        p[j] = r[j]
    end
    rho = 0.0
    for j in 1..@lastcol - @firstcol + 1 do
      rho += r[j] * r[j]
    end
    for cgit in 1..@cgitmax do
        for j in 1..@lastrow - @firstrow + 1 do
            sum = 0.0
            for k in rowstr[j]..rowstr[j + 1] - 1 do
                sum += a[k] * p[colidx[k]]
            end
            q[j] = sum
        end
        d = 0.0
        for j in 1..@lastcol - @firstcol + 1 do
          d += p[j] * q[j]
        end
        @alpha = rho / d
        for j in 1..@lastcol - @firstcol + 1 do
            z[j] = z[j] + @alpha * p[j]
            r[j] = r[j] - @alpha * q[j]
        end
        rho0 = rho
        rho = 0.0
        for j in 1..@lastcol - @firstcol + 1 do
          rho += r[j] * r[j]
        end
        @beta = rho / rho0
        for j in 1..@lastcol - @firstcol + 1 do
            p[j] = r[j] + @beta * p[j]
        end
    end
    for j in 1..@lastrow - @firstrow + 1 do
        sum = 0.0
        for k in rowstr[j]..rowstr[j + 1] - 1 do
            sum += a[k] * z[colidx[k]]
        end
        r[j] = sum
    end
    sum = 0.0
    for j in 1..@lastcol - @firstcol + 1 do
      sum += (x[j] - r[j]) * (x[j] - r[j])
    end
    return Math.sqrt(sum)
  end

  def endWork()
    
    for j in 1..@lastrow - @firstrow + 1 do
        sum = 0.0
        for k in @rowstr[j]..@rowstr[j + 1] - 1 do
            sum += @a[k] * @z[@colidx[k]]
        end
        @r[j] = sum
    end
    sum = 0.0
    for j in 1..@lastcol - @firstcol + 1 do
      sum += (@x[j] - @r[j]) * (@x[j] - @r[j])
    end
    return Math.sqrt(sum)
  end

  def PrintTimers()
    #fmt = DecimalFormat.new("0.000"); 
    puts("  SECTION   Time (secs)")
    ttot = @timer.readTimer(@t_bench); 
    if ttot == 0.0 then
      ttot = 1.0
    else
    end
    for i in 1..@t_last do
        tm = @timer.readTimer(i); 
        if i == @t_init then
            puts("  " + @t_names[i] + ":" + (tm).to_s)
        else
            puts("  " + @t_names[i] + ":" + (tm).to_s + "  (" + (tm * 100.0 / ttot).to_s + "%)")
            if i == @t_conj_grad then
                tm = ttot - tm
                puts("    --> total rest :" + (tm).to_s + "  (" + (tm * 100.0 / ttot).to_s + "%)")
            else
            end
        end
    end
  end

  def executeTask(orderNum)
    synchronize {
    for m in 0..@num_threads-1 do
        @worker[m].synchronize do
            @worker[m].taskOrder = orderNum
            @worker[m].done = false
            @worker[m].alpha = @alpha
            @worker[m].beta = @beta
            @worker[m].cond.signal #notify()
        end
    end
    for m in 0..@num_threads-1 do
        while not @worker[m].done do
            #begin
          self.cond.wait
            #rescue InterruptedException => e
            #ensure
            #end
          self.cond.broadcast #notifyAll()
        end
    end
    }
  end

  def getTime()
    return @timer.readTimer(@t_bench)
  end

  def finalize()
    puts("CG: is about to be garbage collected")
    super.finalize()
  end

  def setupThreads(cg)
    @worker = Array.new(@num_threads, 0.0)
    @master = cg
    div = @na / @num_threads; 
    rem = @na % @num_threads; 
    start = 1; _end = 0; 
    @dmaster = Array.new(@num_threads, 0.0)
    @rhomaster = Array.new(@num_threads, 0.0)
    @rnormmaster = Array.new(@num_threads, 0.0)
    for i in 0..@num_threads-1 do
      _end += div
      if rem != 0 then
        rem -= 1
        _end += 1
      else
      end
      @worker[i] = CGWorker.new(@CLASS, @num_threads, cg, start, _end)
      @worker[i].id = i
      @worker[i].extend(MonitorMixin)
      start = _end + 1
    end
    for i in 0..@num_threads-1 do
      @worker[i].start()
    end
  end

end


  def main(argv)
    #cg =  nil ; 
    BMArgs.parseCmdLineArgs(argv, @BMName)
    clss = BMArgs.clss; 
    np = BMArgs.num_threads; 
    serial = BMArgs.serial; 
    #begin
    cg = CG.new(clss, np, serial)
    cg.extend(MonitorMixin)
    #rescue OutOfMemoryError => e
    #    BMArgs.outOfMemoryMessage()
    #    System.exit(0)
    #ensure
    #end
    cg.runBenchMark()
  end

main(ARGV)
