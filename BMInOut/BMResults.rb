# -*- coding: utf-8 -*-
class BMResults
  def initialize(bid)
    @pid = bid
    @clss = "S"
    @optype = "floating point"
  end
  def initialize(bname, clss,
                 bn1, 
                 bn2,
                 bn3,
                 bniter,
                 btime,
                 bmops,
                 boptype,
                 passed_verification,
                 bserial,
                 num_threads,
                 bid)
      @pid=bid
      @name=bname
      @clss=clss
      @n1=bn1
      @n2=bn2
      @n3=bn3
      @niter=bniter
      @time=btime
      @mops=bmops
      @optype=boptype
      @verified=passed_verification
      @serial=bserial
      @numthreads=num_threads
  end

  def print()
    content = <<EOS
***** NAS Parallel Benchmarks Ruby version (NPB3_0_RUB) #{@name} *****
* Class             = #{@clss}
* Size              = #{@n1}X#{@n2}X#{@n3}
* Iterations        = #{@niter}
* Time in seconds   = #{@time}
* ACCTime           = #{@acctime}
* Mops total        = #{@mops}
* Operation type    = #{@optype}
* Verification      = #{if @verified==1 then "Successful" else "Failed" end}
* 
* Please send all errors/feedbacks to:
* ahya365@gmail.com / https://twitter.com/plus7
***************************************************************
EOS
    lines = content.split("\n")
    width = lines[0].length
    for i in 1..lines.length-1 do
      spaces = width - lines[i].length - 1
      str = ""
      spaces.times{ str = str + " " }
      lines[i] = lines[i] + str + "*"
    end
    puts lines.join("\n")
  end

  def getFromFile(filename)
    @verified = -1
    f = open(filename)
    #TODO: エラー処理
    #エラーが出たら0を返すらしい
    begin
      f.each_line {|line|
        if line.index("Time in seconds =")>=0 then
          keyword="Time in seconds ="
          idx1=line.index(keyword)
          idx1+=keyword.length()
          dbl = line[idx1..-1].to_f
          acctime=dbl
          wctime=dbl
        elsif line.index("Verification    =")>=0 then
          @verified=0
          if ((line.index("successful")>=0 and line.index("successful")<0) or (line.index("SUCCESSFUL")>=0 and line.index("UNSUCCESSFUL")<0)) then
            @verified=1
          end
        elsif line.index("Mop/s total     =")>=0 then
          keyword="Mop/s total     ="
          idx1=line.index(keyword)
          idx1+=keyword.length()
          @mops = line[idx1..-1].to_f
        end
      }
    ensure
     f.close
    end
#    print();
    return 1
  end
  def self.printVerificationStatus(clss, verified, bmname)
    if clss == "U" or verified == -1 then
      verified = -1
      puts " Problem size unknown"
      puts bmname+"."+clss+": Verification Not Performed"
    elsif verified == 1 then
      puts bmname+"."+clss+": Verification Successful"
    else
      puts bmname+"."+clss+": Verification Failed"
    end
  end

  # RubyはオーバーロードできないのでArr接尾辞をつけた
  def self.printComparisonStatusArr( clss, verified, epsilon, 
                             xcr_arr, xcrref_arr, xcrdif_arr)
    for m, xcr, xcrref, xcrdif in (0..xcr_arr.length-1).zip(xcr_arr, xcrref_arr, xcrdif_arr)
      if clss == "U" then
        puts " " + xcr.to_s
      else
        if xcrdif <= epsilon then
          if verified == -1 then verified = 1 end
        else
          verified = 0
          puts "FAILURE: "
        end
        puts m.to_s + ". " + xcr.to_s + " " + xcrref.to_s + " " + xcrdif.to_s
      end
    end
    return verified
  end

  def self.printComparisonStatus( clss, verified, epsilon,
                             xcr, xcrref, xcrdif)
    if clss == "U" then
      puts " " + xcr.to_s
    else
      if xcrdif <= epsilon then
    	if verified == -1 then verified = 1 end
      else
    	verified = 0
    	puts "FAILURE: "
      end
      puts xcr.to_s + " " + xcrref.to_s + " " + xcrdif.to_s
    end
    return verified
  end

end
