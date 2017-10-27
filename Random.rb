#/*
#!-------------------------------------------------------------------------!
#!                                                                         !
#!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.0         !
#!                                                                         !
#!                       J A V A         V E R S I O N                     !
#!                                                                         !
#!                              R A N D O M                                !
#!                                                                         !
#!-------------------------------------------------------------------------!
#!                                                                         !
#!    This benchmark is a serial version of the NPB3_0_JAV Random number   !
#!    generating code.                                                     !
#!                                                                         !
#!    Permission to use, copy, distribute and modify this software         !
#!    for any purpose with or without fee is hereby granted.  We           !
#!    request, however, that all derived work reference the NAS            !
#!    Parallel Benchmarks 3.0. This software is provided "as is"           !
#!    without express or implied warranty.                                 !
#!                                                                         !
#!    Information on NPB 3.0, including the Technical Report NAS-02-008    !
#!    "Implementation of the NAS Parallel Benchmarks in Java",             !
#!    original specifications, source code, results and information        !
#!    on how to submit new results, is available at:                       !
#!                                                                         !
#!           http://www.nas.nasa.gov/Software/NPB/                         !
#!                                                                         !
#!    Send comments or suggestions to  npb@nas.nasa.gov                    !
#!                                                                         !
#!          NAS Parallel Benchmarks Group                                  !
#!          NASA Ames Research Center                                      !
#!          Mail Stop: T27A-1                                              !
#!          Moffett Field, CA   94035-1000                                 !
#!                                                                         !
#!          E-mail:  npb@nas.nasa.gov                                      !
#!          Fax:     (650) 604-3957                                        !
#!                                                                         !
#!-------------------------------------------------------------------------!
#!     Translation to Java and to MultiThreaded Code:                      !
#!     Michael A. Frumkin                                                  !
#!-------------------------------------------------------------------------!
#*/
module NPB3_0_RUB
  class Random
    #//Random Number Multiplier
    AMULT = 1220703125.0 #//=Math.pow(5.0,13);
    #//constants
    D2M46=0.5**46
    I246M1=2**46-1
    def initialize(sd=nil)
      #//default seed
      @tran = 314159265.0   #//First 9 digits of PI
      if sd then
        @seed=sd
      end
    end
    #//Random number generator with an internal/external seed
    # use randlc(a) or randlc(x,a) as in NPB3_0_JAV.Random
    def randlc(x,a=nil)
      if not a then
        a=x
        x=@tran
        updatetran=true
      end
      r23=0.5**23
      r46=r23*r23
      t23=2.0**23
      t46=t23*t23
      #//---------------------------------------------------------------------
      #//   Break A into two parts such that A = 2^23 * A1 + A2.
      #//---------------------------------------------------------------------
      t1=r23*a
      a1=t1.to_i()
      a2=a-t23*a1
      #//---------------------------------------------------------------------
      #//   Break X into two parts such that X = 2^23 * X1 + X2, compute
      #//   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
      #//   X = 2^23 * Z + A2 * X2  (mod 2^46).
      #//---------------------------------------------------------------------
      t1=r23*x
      x1=t1.to_i()
      x2=x-t23*x1
      t1=a1*x2+a2*x1
      t2=(r23*t1).to_i()
      z=t1-t23*t2
      t3=t23*z+a2*x2
      t4=(r46*t3).to_i()
      x=t3-t46*t4
      if updatetran then
        @tran=x
        r46*@tran
      else
        x
      end
    end
    def vranlc(n,x,a,y,offset)
      n.times { |i|
        x=(x.to_i()*a.to_i()) & I246M1.to_i()
        y[offset+i]=x*D2M46
      }
      x
    end
    def ipow(a,exponent)
      #//---------------------------------------------------------------------
      #// Use
      #//   a^n = a^(n/2)*a^(n/2) if n even else
      #//   a^n = a*a^(n-1)       if n odd
      #//---------------------------------------------------------------------
      if exponent != 0 then
        q=a
        r=1
        n=exponent
        while(n>1)
          n2=n/2
          if(n2*2==n)
            @seed=randlc(q,q)
            q=@seed
            n=n2
          else
            @seed=randlc(r,q)
            r=seed
            n=n-1
          end
        end
        @seed=randlc(r,q)
      end
      @seed
    end
    def power(a,n)
      #//c---------------------------------------------------------------------
      #//c     power  raises an integer, disguised as a double
      #//c     precision real, to an integer power
      #//c---------------------------------------------------------------------
      pow=1.0
      nj=n
      aj=a
      while nj!=0
        if nj%2==1 then
          @seed=randlc(pow,aj)
          pow=@seed
        end
        ajj=aj
        @seed=randlc(aj,ajj)
        aj=@seed
        nj=nj/2
      end
      pow
    end
  end
end
