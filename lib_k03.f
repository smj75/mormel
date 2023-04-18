C
C MELTING PARAMETERISATION OF KATZ ET AL., 2003, G-CUBED
C
      function XK03(ip)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
      double precision frac(ixzmax,ixzmax),nfrac(ixzmax,ixzmax)
C MODULE FOR USE BY ROOT FINDER WHEN WATER IS PRESENT
      external fkatz
C PARAMETERS FROM KATZ ET AL. TABLE 2
      parameter (a1=1085.7, a2=132.9, a3=-5.1)
      parameter (b1=1475.0, b2=80.0, b3=-3.2)
      parameter (c1=1780.0, c2=45.0, c3=-2.0)
      parameter (r1=0.50, r2=0.08)
      parameter (beta1=1.5, beta2=1.5)
C PARAMETERS FOR ROOT FINDING
      parameter (fmax=1.)
      parameter (fmin=0.)
      parameter (facc=1.e-4)
C OTHER PARAMETERS
      parameter(ierr=-999)
C PRESSURE, TEMPERATURE, WATER CONTENT, MODAL CPX IN HERE
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C      common /deltem/x,xsat,dtem
      common /deltem/dtem
C COMMUNICATION WITH FKATZ
      common /tkatz/td,tsol,tlliq,tliq,fcpx,tcpx,rcpx
C      write(6,*)'input:',p,t,xh2o,cpxm
C
C INITIALIZE OUTPUT VALUE
C
      pkatz=0.0
C
C SOLIDUS (EQUATION 4)
C
      tsol=a1+a2*p+a3*p*p
C      write(6,*)'tsol,p:',tsol,p,t
C
C BELOW SOLIDUS?-MOVED TO A NEW LOCATION TO ALLOW
C PLOTTING OF 
C
C      if (xh2o.eq.0. .and. t.le.tsol) then
C         pkatz=0.
C         return
C      endif
C
C LHERZOLITE LIQUIDUS AND LIQUIDUS (EQUATIONS 5 & 10)
C
      tlliq=b1+b2*p+b3*p*p
      tliq=c1+c2*p+c3*p*p
C      write(6,*)'tlliq,tliq:',tlliq,tliq
C
C MELT FRACTION AT WHICH CPX IS EXHAUSTED (EQUATIONS 6 & 7)
C
      rcpx=r1+r2*p
      fcpx=cpxm/rcpx
C      write(6,*)'rcpx,fcpx:',rcpx,fcpx
C
C WRITE OUT DATA FOR KATZ FIGURE 1
C
C      write(9,*)p,tsol,tlliq,tliq,fcpx
C
C BELOW SOLIDUS?
C
C      if (xh2o.eq.0. .and. t.le.tsol) then
C         write(6,*)'t less than tsol:',t-tsol
C         pkatz=0.
C         return
C      endif
C
C IF WATER CONTENT ZERO, MELT FRACTION CAN BE CALUCLATED DIRECTLY
C
      if (xh2o.eq.0.) then
        if (t.le.tsol) then
c          write(6,*)'t less than tsol:',t-tsol
          pkatz=0.
		      XK03=pkatz
		      return
        else
C
C MELT FRACTION WHEN CPX AVAILABLE FOR MELTING (EQUATIONS 2 & 3)
C
          td=(t-tsol)/(tlliq-tsol)
          pkatz=td**beta1
c         write(6,*)'results:',t,td,pkatz
C
C MELT FRACTION WHEN OPX IS ALSO MELTED (EQUATIONS 8-10)
C
          if (pkatz.gt.fcpx) then
C            write(6,*)'Dry OPX melting'
            tcpx=fcpx**(1./beta1)*(tlliq-tsol)+tsol
C CHANGING TSOL IN THE NEXT EQN FOR TCPX
            td=(t-tcpx)/(tliq-tcpx)
            pkatz=fcpx+(1.-fcpx)*td**beta2
C            write(6,*)'tcpx,td,pkatz:',tcpx,td,pkatz
          endif
        endif
C
C WRITE OUT DATA FOR KATZ FIGURE 2 PAGE 4
C
C         write(10,*)t,pkatz
C
C IF WATER PRESENT, NEED TO FIND MELT FRACTION USING ROOT FINDER
C SINCE IT OCCURS ON BOTH SIDES OF EQUATION 19
C
      elseif (xh2o.gt.0.) then
C
C CHECK THERE IS A ROOT BETWEEN MINIMUM AND MAXIMUM MELT FRACTIONS
C
C         write(6,*)'Checking for a root'
         f1=fkatz(fmin)
         f2=fkatz(fmax)
C         write(6,*)'fmin,fmax:',f1,f2
         if(f1*f2.ge.0.) then
C         write(6,*)'No Melting f1f2 positive'
         endif
         if (f1*f2.lt.0.) then
C            write(6,*)'f1*f2.lt.0?: ',f1*f2
C
C YES, FIND ROOT BY BISECTION
C
C            write(6,*)'Finding Root....'
            pkatz=rtbis(fkatz,fmin,fmax,facc)
C            write(6,*)'Found a Root....'
C            write(6,*)'pkatz:',pkatz
         else
C
C NO, NO MELTING
C
            pkatz=0.
         endif
      endif
C
C SET RETURN VALUE FOR THE FUNCTION
C
      XK03=pkatz
C      write(6,*)'Melt Fraction: ',pkatz 
      return
      end
C
C------------------------------------------------------------------
C 
C GIVEN ESTIMATE OF MELT FRACTION, CALCULATE MELT FRACTION WHEN WATER IS PRESENT
C RETURN DIFFERENCE BETWEEN STARTING AND CALULATED MELT FRACTION
C WHEN STARTING VALUE IS CORRECT, DIFFERENCE IS ZERO
C
      function FKATZ(fin)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
      double precision frac(ixzmax,ixzmax),nfrac(ixzmax,ixzmax)
C
C PARAMETERS FROM KATZ ET AL. TABLE 2
C
      parameter (beta1=1.5, beta2=1.5)
      parameter (pk=43.0, gamma=0.75, dh2o=0.01)
      parameter (chi1=12.0, chi2=1.0, plam=0.6)
C
C COMMUNICATION WITH PKATZ
C
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C      common /deltem/x,xsat,dtem
      common /deltem/dtem
      common /tkatz/td,tsol,tlliq,tliq,fcpx,tcpx,rcpx
      common /fig3/xsol
C
C SATURATION CONCENTRATION OF WATER (EQUATION 17)
C
      xsat=chi1*p**plam+chi2*p
C      write(6,*)'p,xsat:',p,xsat
C 
C ASIDE EQUATION JUST TO WRITE OUT WATER SATURATED SOLIDUS
C WHICH IS NOT DEPENDENT ON CALCULATIONG THE MELT FRACTION
C
C      xsol2=tsol-(pk*xsat**gamma)
C      write(11,*)p,xsol2
C      write(6,*)'Written fig3 data'
C
C CONCENTRATION OF WATER (EQUATION 18)
C I THINK THERE IS A PROBLEM HERE:
C X=XSAT SHOULD BE IF X>XSAT, NOT XH20!
C
C      if (xh2o.gt.xsat) then
C         x=xsat
C      else
C         x=xh2o/(dh2o+fin*(1.-dh2o))
C      endif
C WATER CONCENTRATION FOR EQUILLIBRIUM MELTING
       x=xh2o/(dh2o+fin*(1.-dh2o))
C WATER CONCENTRATION FOR FRACTIONAL MELTING
C       x=(xh2o/dh2o)*(1.-fin)**((1./dh2o)-1.)
C       write(6,*)'x,fin:',x,fin
       if(x.gt.xsat) then
       x=xsat
       endif
C
C CORRECTION TO SOLIDUS AND LIQUIDUS FOR HYDROUS MELTING (EQUATION 16)
C
      dtem=pk*x**gamma
C      write(6,*)'xsat,x,dt:',xsat,x,dtem
C
C WRITE OUT DATA FOR KATZ FIGURE 3 PAGE 5
C
      xsol=tsol-dtem
C      write(11,*)p,tsol,xsol
C      write(6,*)'tsol,xsol',tsol,xsol,t
C
C IS TEMPERATURE ABOVE SOLIDUS?
C
C      if (t.le.tsol-dtem) then
       if (t.le.xsol) then
C         write(6,*)'Hydrous temperature less than solidus'
C         write(6,*)'No Melting!'
C         write(6,*)'t:',t,'xsol:',xsol
C BECAUSE FOUT WOULD BE ZERO, SO FKATZ=FIN-FOUT
C WOULD GIVE THE SAME RESULT
         fkatz=fin
         return
      endif
C
C MELT FRACTION WHEN CPX AVAILABLE FOR MELTING (EQUATION 19)
C
      td=(t-(tsol-dtem))/(tlliq-tsol)
C      write(6,*)'t,tlliq,td:',t,tlliq,td
      fout=td**beta1
C      write(6,*)'fout:',fout
C
C MELT FRACTION WHEN OPX IS ALSO MELTED (EQUATIONS 8-10 MODIFIED)
C
      if (fout.gt.fcpx) then
c         write(6,*)'wet opx melting:',fout
         tcpx=fcpx**(1./beta1)*(tlliq-tsol)+tsol-dtem
C TLIQ ALSO IS EFFECTED BY THE WATER CONTENTS, THEREFORE -DTEM
C ONLY REASON TLLIQ AND TSOL DON'T HAVE ABOVE IS BECAUSE THEY
C WOULD CANCEL EACH OTHER OUT WITH THE MINUS SIGN
         td=(t-tcpx)/(tliq-dtem-tcpx)
         fout=fcpx+(1.-fcpx)*td**beta2
C         write(6,*)'fcpx,tcpx,tliq,td,fout:',fcpx,tcpx,tliq,td,fout
      endif
C
C DIFFERENCE BETWEEN INPUT AND OUTPUT MELT FRACTION
C
      fkatz=fin-fout
C      write(6,*)'fkatz:',fin,fout,fkatz
      return
      end
C
C------------------------------------------------------------------
C 
C FUNCTION TO CALCULATE SOLIDUS TEMPERATURE
C
      function SK03(p)
      implicit double precision(a-h,o-z)
C
C PARAMETERS FROM KATZ ET AL. TABLE 2
C
      parameter (a1=1085.7, a2=132.9, a3=-5.1)
C
      sk03=a1+a2*p+a3*p*p
C
      return
      end
C
C------------------------------------------------------------------
C 
C CALCULATE VARIOUS DERIVATIVES NEEDED TO CALCULATE INSTANTANEOUS
C MELT PRODUCTION RATES
C
      function DK03(ip)
      implicit double precision(a-h,o-z)
C
C COMMUNICATION OF VELOCITIES FROM BW94 PROGRAM FOR NUMERICAL
C SOLUTION TO DF/Dt
C
      parameter (ixzmax=2001)
      parameter (fmode=0)
      double precision vx(ixzmax,ixzmax)
      double precision vz(ixzmax,ixzmax)
      double precision frac(ixzmax,ixzmax)
      double precision dpda(ixzmax,ixzmax)
      double precision dtrda(ixzmax,ixzmax)
      double precision dfda(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      common /grid/xmax,zmax,dx,dz
      common /velo/vx,vz,uhalf,vmax
      common /melfrc/f1x,f2x,f1z,f2z,fx
C
C COMMUNICATION OF DERIVATIVES, WHERE
C  A=TIME
C  P=PRESSURE
C  T=DIMENSIONLESS TEMPERATURE
C  TL=LIQUIDUS TEMPERATURE
C  TR=REAL TEMPERATURE (BEFORE CORRECTION FOR LATENT HEAT OF MELTING) 
C  TS=SOLIDUS TEMPERATURE
C  X=MELT FRACTION
C TIME DERIVATIVES ARE INPUT (I.E. DP/DA, DT/DA, DTR/DA)
C REMAINING DERIVATTIVES ARE OUTPUT
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
C
C PRESSURE, TEMPERATURE, WATER CONTENT, MODAL CPX INPUT HERE
C NB WATER CONTENT AND MODAL CPX NOT USED IN THIS PARAMETERISATION
C
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C      common /deltem/x,xsat,dtem
      common /deltem/dtem
      common /intder/dpdac,dtrdac,f,dfdac
C
C DIMENSIONLESS, SOLIDUS AND LIQUIDUS TEMPERATURE IN HERE
C
      common /tkatz/td,tsol,tlliq,tliq,fcpx,tcpx,rcpx
C
C CONSTANTS IN PARAMETERISATIONS FOR SOLIDUS, LIQUIDUS & MELT FRACTION
C PARAMETERS FROM KATZ ET AL. TABLE 2
C
      parameter (a1=1085.7, a2=132.9, a3=-5.1)
      parameter (b1=1475.0, b2=80.0, b3=-3.2)
      parameter (c1=1780.0, c2=45.0, c3=-2.0)
      parameter (r1=0.50, r2=0.08)
      parameter (beta1=1.5, beta2=1.5)
      parameter (pk=43.0, gamma=0.75, dh2o=0.01)
      parameter (chi1=12.0, chi2=1.0, plam=0.6)
C
C FIND DIMENSIONLESS, SOLIDUS AND LIQUIDUS TEMPERATURES AND 
C MELT FRACTION
C
C      dk03=XK03(ip)
      dk03=f
C      write(6,*)'f:',dk03
C
C IF BELOW SOLIDUS, DX/Dt MUST BE ZERO
C
      if (dk03.le.0.) then
        dfdac=0.
C        write(6,*)'No melting so return'
        return
      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C RECALCULATE VALUES ORIGINALLY CALCULATED IN XK03 FOR USE IN
C CALCULATING THE DERIVATIVES
C VALUES REQUIRED FOR:
C ANHYDROUS: TSOL,TLLIQ,TLIQ,TCPX,RCPX,FCPX
C ADDITIONAL FOR HYDROUS: DTEM,X,XSAT
C
C TSOL:
      tsol=a1+a2*p+a3*p*p
C TLLIQ:
      tlliq=b1+b2*p+b3*p*p
C TLIQ:
      tliq=c1+c2*p+c3*p*p
C RCPX:
      rcpx=r1+r2*p
C FCPX:
      fcpx=cpxm/rcpx
C TCPX:
      tcpx=fcpx**(1./beta1)*(tlliq-tsol)+tsol
C ADDITIONAL CALCULATIONS FOR HYDROUS CALCULATIONS:
C
      x=0.0
      if (xh2o.gt.0.) then
C XSAT:
        xsat=chi1*p**plam+chi2*p
C X:
C WATER CONCENTRATION FOR EQUILLIBRIUM MELTING
        x=xh2o/(dh2o+f*(1.-dh2o))
C WATER CONCENTRATION FOR FRACTIONAL MELTING
C       x=(xh2o/dh2o)*(1.-f)**((1./dh2o)-1.)
C
        if(x.gt.xsat) then
          x=xsat
        endif
C DTEM:
        dtem=pk*x**gamma
C TCPX FOR HYDROUS:
        tcpx=fcpx**(1./beta1)*(tlliq-tsol)+tsol-dtem
      endif
C      write(6,*)'xsat,x,dt:',xsat,x,dtem
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C TESTING: NUMERICAL SOLUTION FOR DF/DT
C DF/Da=vx*df/dx+vz*df/dz, WHERE d IS DELTA
C CALCULATE GRADIENTS FROM PREVIOUS NODE TO THE ONE AFTER THE 
C NODE IN CONSIDERATION, NOTE NEED TO CHANGE THE VX AND VZ FOR
C ONES CALCULATED IN THE INTEG SUBROUTINE
C
C      if (fmode.eq.1) then
CC          write(6,*)'Calculating numerical solution for DF/Dt'
C          do ix=2,2
C             do iz=30,30
CC                write(6,*)vx(ix,iz),f2x,f1x,dx
CC                write(6,*)vz(ix,iz),f2z,f1z,dz
C                dfdx=vx(ix,iz)*(f2x-f1x)/(2.*dx)
C                dfdz=vz(ix,iz)*(f1z-f2z)/(2.*dz)
C                dfdac=dfdx+dfdz
C                write(6,*)'Numerical Melt Production Rate: ',dfda
C             enddo
C          enddo
C 
C ANALYTICAL SOLUTION FOR DF/Dt, FROM DIFFERENTIATING KATZ
C EQUATIONS AS GIVEN
C
      if (fmode.eq.0) then
C         write(6,*)'Calculating analytical solution for DF/Dt'
C
C DF/Dt FOR ANHYDROUS MELTING, BELOW FCPX-OUT
C        
         if (dk03.le.fcpx .and. xh2o.eq.0.) then
C           write(6,*)'Anhydrous below Fcpx'
C
C SET DIMENSIONLESS TEMPERAUTRE TDIM=(T-TSOL)/(TLLIQ-TSOL)
C
           tdim=(t-tsol)/(tlliq-tsol)
C           write(6,*)'t,tsol,tlliq,tdim: ',t,tsol,tlliq,tdim
C
C CALCULATE dT'/dp AND dT'/dT, THE CHANGE IN DIMENSIONLESS
C TEMPERATURE WRT PRESSURE AND TEMPERATURE
C
           dtdp=((-a2-2.*a3*p)*(tlliq-tsol)-(t-tsol)*
     &          (b2+2.*b3*p-a2-2.*a3*p))/((tlliq-tsol)**2.)
           dtdtr=1./(tlliq-tsol)
C           write(6,*)'dtdp: ',dtdp,'dtdtr: ',dtdtr
C           write(6,*)'dpdac: ',dpdac,'dtrdac: ',dtrdac
C
C CALCULATE DX/Dt
C THE CHANGE IN MELT FRACTION WRT DIMENSIONLESS TEMPERATURE
C AND THE CHANGE IN MELT FRACTION WRT TIME
C
           dfdt=beta1*(tdim**(beta1-1.))
C           write(6,*)'dfdt: ',dfdt
           dfdac=dfdt*(dpdac*dtdp+dtrdac*dtdtr)
C           write(6,*)'Melt Rate Lib_K03: ',dfdac
C
C DF/Dt FOR ANHYDROUS MELTING, ABOVE FCPX-OUT
C
         else if (dk03.gt.fcpx .and. xh2o.eq.0.) then
C           write(6,*)'Anhydrous above Fcpx'
C
C SET DIMENSIONLESS TEMPERATURE TDIM=(T-TCPX)/(TLIQ-TCPX)
C TCPX IS CALCULATED IN XK03
C
C           tdim=(t-tcpx)/(tliq-tcpx)
C
C CALCULATE DERIVATIVES OF FCPX-OUT WRT PRESSURE AND TEMPERATURE
C
C           dfcdp=-1.*r2*cpxm/(rcpx**2.)
C           dfcdtr=0.
C
C CALCULATE THE CHANGE IN DIMENSIONLESS TEMPERATURE WRT PRESSURE
C
C           fq=t-tcpx
C           gq=tliq-tcpx
C           fqdash=(1./beta1)*(fcpx**(1./beta1-1.))
C     &            *(r2*cpxm/(rcpx**2.))*(tlliq-tsol)
C     &            -(fcpx**(1./beta1))*(b2+2.*b3*p-a2-2.*a3*p)
C     &            -a2-2.*a3*p
C           gqdash=c2+2.*c3*p+fqdash
C           dtdp=(fqdash*gq-fq*gqdash)/(gq**2.)
C
C CALCULATE THE CHANGE IN DIMENSIONLESS TEMPERATURE WRT TEMPERATURE
C
C           dtdtr=1./(tliq-tcpx)**2.
C
C CALCULATE DX/DT' THE CHANGE IN MELT FRACTION WRT DIMENSIONLESS
C TEMPERATURE 
C
C           dfdt=beta2*tdim**(beta2-1.)
C
C CALCULATE THE CHANGE IN MELT FRACTION WRT PRESSURE
C
C           dfdp=dfcdp*(1.-tdim**beta2)+dfdt*dtdp*(1.-fcpx)
C
C CALCULATE THE CHNAGE IN MELT FRACTION WRT TEMPERATURE
C
C           dfdtr=dfcdtr*(1.-tdim**beta2)+dfdt*dtdtr*(1.-fcpx)
C
C CALCULATE THE CHANGE IN MELT FRACTION WRT TIME
C
C           write(6,*)'Melt Rate calculation: ',dpdac,dfdp,dtrdac,dfdtr
C           dfdac=dpdac*dfdp+dtrdac*dfdtr
C           write(6,*)'Melt Rate: ',dfdac
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NEW NEATER EQUATIONS FOR ANHYDROUS ABOVE FCPX-OUT CCCCCCCCCCCCCCCCCCC
C
C DEFINING A DIMENSIONLESS TEMPERATURE TDIM
C           write(6,*)'Anhydrous above fcpx!'
C           stop
C
           tdim=(t-tcpx)/(tliq-tcpx)
C
C CHANGE IN FCPX W.R.T. PRESSURE (DFCDP)
C
           dfcdp=(-1.)*r2*cpxm/(rcpx**2.)
C
C CHANGE IN DIMENSIONLESS TEMPERATURE W.R.T. TEMPERATURE (DTDTR)
C
           dtdtr=1./(tliq-tcpx)
C
C CHANGE IN TCPX W.R.T. PRESSURE (DTCDP)
C
           dtcdp=(1./beta1)*(fcpx**(1./beta1-1.))*dfcdp*(tlliq-tsol)
     &           +(b2+2.*b3*p-a2-2.*a3*p)*(fcpx**(1./beta1))
     &           +a2+2.*a3*p
C
C CHANGE IN DIMENSIONLESS TEMPERATURE W.R.T. PRESSURE (DTDP)
C
           dtdp=((-1.)*dtcdp*(tliq-tcpx)-(c2+2.*c3*p-dtcdp)*(t-tcpx))
     &          /((tliq-tcpx)**2.)
C
C CHANGE IN MELT FRACTION W.R.T. PRESSURE (DFDP)
C WRONG!!!EQUATION.
C           dfdp=dfcdp*(1.-(tdim**beta2))
C     &          +beta2*(tdim**(beta2-1.))*dtcdp*(1.-fcpx)
           dfdp=dfcdp*(1.-(tdim**beta2))
     &          +beta2*(tdim**(beta2-1.))*dtdp*(1.-fcpx)
C
C CHANGE IN MELT FRACTION W.R.T. TEMPERATURE (DFDTR)
C
           dfdtr=beta2*(tdim**(beta2-1.))*dtdtr*(1.-fcpx)
C
C CHANGE IN MELT FRACTION WITH RESPECT TO TIME (DFDAC)
C 
           dfdac=dpdac*dfdp+dtrdac*dfdtr
C           write(6,*)'Melt Rate: ',dfdac
C
C
C DF/Dt FOR HYDROUS MELTING, BELOW FCPX-OUT
C
         else if (dk03.le.fcpx .and. xh2o.gt.0.) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C SET DIMENSIONLESS TEMPERATURE TDIM=(T-TSOL+DTEM)/(TLLIQ-TSOL)
C DTEM CALCULATED IN XKO3
C
C           tdim=(t-(tsol-dtem))/(tlliq-tsol)  
C
C CALCULATE THE CHANGE IN MELT FRACTION WRT PRESSURE
C SUFFIXES OF A-E FOR EASE OF CALCULATION
C
C           dfdpa=beta1*tdim**(beta1-1.)/(tlliq-tsol)
C           dfdpb=beta1*tdim**(beta1-1.)*(-a2-2.*a3*p)/(tlliq-tsol)**2.
C           dfdpc=pk*(xh2o/(dh2o+dk03*(1.-dh20)))**gamma
C           dfdpd=pk*gamma*(xh2o/(dh2o+dk03*(1.-dh20)))**(gamma-1.)
C           dfdpe=(1.-dh2o)*xh2o/(dh2o+dk03*(1.-dh2o))**2.
C COMBINING ALL THE PARTS TOGETHER
C           dfdp=(dfdpa*(-a2-2.*a3*p)-dfdpb*(t-tsol+dfdpc)
C     &          *(b2+2.*b3*p-a2-2.*a3*p))/(1.+dfdpa*dfdpd*dfdpe)
C
C CALCULATE THE CHANGE IN MELT FRACTION WRT TEMPERATURE
C USING SOME EXPRESSION CALCULATED AS A-E ABOVE
C
C           dfdtr=dfdpa/(1.+dfdpa*dfdpd*dfdpe)
C
C CALCULATE THE CHANGE IN MELT FRACTION WRT TIME
C     
C           dfdac=dpdac*dfdp+dtrdac*dfdtr
C           write(6,*)'Melt Rate: ',dfdac         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NEW NEATER AND ULTIMATELY CORRECT EQUATIONS FOR THE HYDROUS
C PARAMETERISATION OF MELTING BELOW FCPX
C
C SET DIMENSIONLESS TEMPERATURE TDIM=(TREAL-(TSOL-DTEM))/(TLLIQ-TSOL)
C
           tdim=(t-(tsol-dtem))/(tlliq-tsol)
C
C SET THE VARIOUS TEMPERATURE DERIVATIVES W.R.T. PRESSURE
C
           dtsolp=a2+2.*a3*p
           dtlliq=b2+2.*b3*p
           dtliqp=c2+2.*c3*p
C
C CALCULATE THE CHANGE IN MELT FRACTION W.R.T. PRESSURE
C
C DFDP FOR BELOW CPX, EQUILIBRIUM MELTING
           dfdp=((-dtsolp)*(tlliq-tsol)
     &           -(dtlliq-dtsolp)*(t-(tsol-dtem)))
     &          /((((tlliq-tsol)**2.)/(beta1*tdim**(beta1-1.)))
     &           +((pk*gamma*x**(gamma-1.)*(1.-dh2o)*xh2o*(tlliq-tsol))
     &          /(dh2o+dk03*(1.-dh2o))**2.))
C DFDP FOR BELOW CPX, FRACTIONAL MELTING
C           dfdp=((-dtsolp)*(tlliq-tsol)
C     &           -(dtlliq-dtsolp)*(t-(tsol-dtem)))
C     &          /((((tlliq-tsol)**2.)/(beta1*tdim**(beta1-1.)))
C     &           +((pk*gamma*x**(gamma-1.)*(xh2o/dh2o)*(1./dh2o)
C     &          *((1.-dk03)**((1./dh2o)-2.))*(tlliq-tsol))))
C
C CALCULATE THE CHANGE IN MELT FRACTION W.R.T. REAL TEMPERATURE
C
C DFDT FOR BELOW CPX, EQUILLIBRIUM MELTING
           dfdtr=1./((tlliq-tsol)/(beta1*tdim**(beta1-1.))
     &              +(pk*gamma*x**(gamma-1.)*(1.-dh2o)*xh2o)
     &              /(dh2o+dk03*(1.-dh2o))**2.)
C DFDT FOR BELOW CPX, FRACTIONAL MELTING
C           dfdtr=1./((tlliq-tsol)/(beta1*tdim**(beta1-1.))
C     &              +(pk*gamma*x**(gamma-1.)*(xh2o/dh2o)*(1./dh2o)
C     &              *((1.-dk03)**((1./dh2o)-2.))))
C    
C           write(6,*)'t,p,tlliq,tsol,dtlliq,dtsol,dtem,tdim,dk03:',
C     &         t,p,tlliq,tsol,dtlliq,dtsol,dtem,tdim,dk03
C
C CALCULATE THE CHANGE IN MELT FRACTION WITH RESPECT TO TIME
C
           dfdac=dpdac*dfdp+dtrdac*dfdtr
C           write(6,*)'Melt Rate: ',dfdac
C           write(6,*)'dpdac,dfdp,dtrdac,dfdtr:',dpdac,dfdp,dtrdac,dfdtr
C
C DF/Dt FOR HYDROUS MELTING, ABOVE FCPX-OUT
C
         else if (dk03.gt.fcpx .and. xh2o.gt.0.) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SET DIMENSIONLESS TEMPERATURE TDIM=(T-TCPX)/(TLIQ-TCPX)
C TCPX IS CALCULATED IN XK03
C
C           tdim=(t-tcpx)/(tliq-tcpx)
C 
C CALCULATE THE CHANGE IN MELT FRACTION WRT PRESSURE
C SUFFIXES A- FOR EASE OF CALCULATION
C
C           dfdpa=beta2*tdim**(beta2-1.)*(1.-fcpx)/(tliq-tcpx)
C           dfdpb=beta2*tdim**(beta2-1.)*(1.-fcpx)/(tliq-tcpx)**2.
C           dfdpc=fcpx**(1./beta1-1.)
C           dfdpd=1./beta1*fcpx**(1./beta1-1.)
C           dfdpe=b2+2.*b3*p-a2-2.*a3*p
C           dfdpf=r2*cpxm/rcpx**2.
C           dfdpg=pk*gamma*(xh2o/(dh2o+dk03*(1.-dh20)))**(gamma-1.)
C           dfdph=(1.-dh2o)*xh2o/(dh2o+dk03*(1.-dh2o))**2.
CC COMBINING ALL THE PARTS TOGETHER
C           dfdp=(-dfdpf*(1.-tdim**beta2)
C     &          +dfdpa*(dfdpd*dfdpf*(tlliq-tsol)
C     &          -dfdpc*dfdpe-a2-2.*a3*p)
C     &          -dfdpb*(c2+2.*c3*p+dfdpd*dfdpf*(tlliq-tsol)
C     &          -dfdpc*dfdpe-a2-2.*a3*p)*(t-tcpx))
C     &          /(1.-dfdpa*dfdpg*dfdph+dfdpb*(t-tcpx)*dfdpg*dfdph)
C
C CALCULATE THE CHANGE IN MELT FRACTION WRT TEMPERATURE
C USING SOME OF THE EXPRESSION CALCULATED ABOVE
C
C           dfdtr=dfdpa/(1.+dfdpb*dfdpg*dfdph*(t-tcpx)-1./(tliq-tcpx))
C
C CALCULATE THE CHANGE IN MELT FRACTION WRT TIME
C
C           dfdac=dpdac*dfdp+dtrdac*dfdtr
C           write(6,*)'Melt Rate: ',dfdac
C           dfda(ix,iz)=0.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NEW CORRECT EQNS FOR F GREATER THAN FCPX
C           write(6,*)'Hydrous above fcpx'
C           stop
C
C SET TDIM=(T-TCPX)/(TLIQ-DTEM-TCPX)
C
           tdim=(t-tcpx)/(tliq-dtem-tcpx)
C
C SET THE VARIOUS TEMPERATURE DERIVATIVES W.R.T. PRESSURE
C
           dtsolp=a2+2.*a3*p
           dtlliq=b2+2.*b3*p
           dtliqp=c2+2.*c3*p
C
C CALCULATE THE CHANGE IN MELT FRACTION W.R.T. PRESSURE
C
C LET AA, BB, CC, DD EQUAL VARIOUS COMPLICATED PARTS OF THE
C EQUATION
C
           aa=(1./beta1)*fcpx**(1./beta1-1.)*(((-r2)*cpxm)
     &        /(rcpx**2.))
     &        *(tlliq-tsol)
     &        +(dtlliq-dtsolp)*fcpx**(1./beta1-1.)+dtsolp
           bb=((1.-tdim**beta2)*(((-r2)*cpxm)/(rcpx**2.))
     &        *((tliq-dtem-tcpx)**2.))
     &        /((1.-fcpx)*beta2*tdim**(beta2-1.))
           cc=((tliq-dtem-tcpx)**2.)/((1.-fcpx)*beta2*tdim**(beta2-1.))
C FOR FOPX, EQUILLIBRIUM MELTING
           dd=(pk*gamma*x**(gamma-1.)*(1.-dh2o)*xh2o)
     &        /((dh2o+dk03*(1.-dh2o))**2.)
C FOR FOPX, FRACTIONAL MELTING
C           dd=pk*gamma*x**(gamma-1.)*(xh2o/dh2o)
C     &        *(1./dh2o)*(1.-dk03)**((1./dh2o)-2.)
C NOW PUT TOGETHER TO COMPLETE DFDP
           dfdp=(bb+aa*(dtem+tcpx-tliq)-(dtliqp-aa)*(t-tcpx))
     &          /(cc+dd*(tliq-dtem-tcpx))
C
C CALCULATE THE CHANGE IN MELT FRACTION W.R.T. TEMPERATURE
C
C LET DD EQUAL PART OF DDTEM/DP
C
C           dd=(pk*gamma*x**(gamma-1.)*(1.-dh2o)*xh2o)
C     &        /((dh2o+dk03*(1.-dh2o))**2.)
C NOW COMPLETE DFDTR
           dfdtr=(tliq-dtem-tcpx)
     &           /(((tliq-dtem-tcpx)**2.)
     &           /(beta2*tdim**(beta2-1.)*(1.-fcpx)) 
     &            -dd*(tliq-dtem-tcpx))
C
C CALCULATE THE CHANGE IN MELT FRACTION W.R.T. TIME
           dfdac=dpdac*dfdp+dtrdac*dfdtr
C           write(6,*)'dpdac,dfdp,dtrdac,dfdtr:',dpdac,dfdp,dtrdac,dfdtr
C     &                                        ,tdim
C           write(6,*)'Melt Rate: ',dfdac
C          
C
         endif
      endif
      return
      end
C
C
C------------------------------------------------------------------
C
C FIND ROOT USING BISECTION METHOD FROM NUMERICAL RECIPES
C
      FUNCTION rtbis(func,x1,x2,xacc)
      implicit double precision(a-h,o-z)
      INTEGER JMAX
      double precision rtbis,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      double precision dx,f,fmid,xmid
C      common /mf/p,t,xh2o,cpxm
C      common /tkatz/td,tsol,tlliq,tliq,fcpx
C      common /fig3/xsol
C      write(6,*)'Using rtbis'
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.)then 
            write(6,*)'rtbis-no root','f,fmid: ',f,fmid 
            stop
      endif
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis=xmid
C        write(6,*)'Absolute of dx','xacc',abs(dx),xacc
        if(abs(dx).lt.xacc .or. fmid.eq.0.) then
C         write(6,*) 'rtbis succcessful'
         return
C        write(11,*)t,tsol,xsol
C        write(6,*)'Written data for figure 3!'
C        return
        endif
11    continue
      write(6,*)'too many bisections in rtbis' 
      stop
      END
