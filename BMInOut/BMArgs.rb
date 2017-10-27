# -*- coding: utf-8 -*-
class BMArgs
  @@CLASS = 'U'
  @@num_threads=4
  @@serial = true
  def initialize()
    @@CLASS='U'
    @@num_threads=4
    @@serial=true
  end
  def self.num_threads
    return @@num_threads
  end
  def self.clss
    return @@CLASS
  end
  def self.serial
    return @@serial
  end
  def self.parseCmdLineArgs(argv, bmname)
    for arg in argv do
      if arg=="SERIAL" or
         arg=="serial" or
         arg=="-serial" or
         arg=="-SERIAL" then
        @@serial = true
      elsif arg.start_with?("class=") or
          arg.start_with?("CLASS=") or
          arg.start_with?("-class") or
          arg.start_with?("-CLASS") then
        if arg.length()>6 then @@CLASS = arg[6].chr.upcase() end
        if @@CLASS != "A" and
           @@CLASS != "B" and
           @@CLASS != "C" and
           @@CLASS != "S" and
           @@CLASS != "W" then
          puts "classes allowed are A,B,C,W and S."
          commandLineError(bmname)
        end
      elsif arg.start_with?("np=") or
          arg.start_with?("NP=") or
          arg.start_with?("-NP") or
          arg.start_with?("-np") then
        if arg.length()>3 then 
          @@num_threads = arg[3..-1].to_i 
          @@serial = false
        end #TODO: 厳密なエラーチェック
      end
    end
  end
  def self.commandLineError(bmname)
    puts "synopsis: java "+bmname+" CLASS=[ABCWS] -serial [-NPnnn]"
    puts "[ABCWS] is the size class \n" +
         "-serial specifies the serial version and\n" +
         "-NP specifies number of threads where nnn "+
         "is an integer"
    exit(1)
  end
  def self.outOfMemoryMessage()
    puts "The java maximum heap size is "+
         "to small to run this benchmark class"
    puts "To allocate more memory, use the -mxn option" +
         " where n is the number of bytes to be allocated"
  end
  def self.banner(bmname, clss, serial, np)
    puts " NAS Parallel Benchmarks Ruby version (NPB3_0_RUB)"
    if serial then puts " Serial Version "+bmname+"."+clss
    else puts " Multithreaded Version "+bmname+"."+clss+" np="+np.to_s end
  end
end

