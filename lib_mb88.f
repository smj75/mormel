C RACHEL WALTERS 27TH JUNE 2006
C MELTING PARAMETERISATION OF MCKENZIE & BICKLE, 1988,
C JOURNAL OF PETROOGY, V29, P625-679.
C
C------------------------------------------------------------------
C 
C FUNCTION TO RETURN MELT FRACTION X GIVEN PRESSURE AND TEMPERATURE
C
      function XMB88(ip)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
      double precision frac(ixzmax,ixzmax),nfrac(ixzmax,ixzmax)
C
C PRESSURE, TEMPERATURE, WATER CONTENT, MODAL CPX INPUT HERE
C NB WATER CONTENT AND MODAL CPX NOT USED HERE
C
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C
C DIMENSIONLESS, SOLIDUS AND LIQUIDUS TEMPERATURE OUTPUT HERE
C
      common /tmb88/td,tsol,tliq
C
C CONSTANTS IN PARAMETERISATIONS FOR SOLIDUS, LIQUIDUS & MELT FRACTION
C
      data as,bs,cs,t0s/4.968e-4,1.2e-2,136.,1100./
      data al,bl,cl,dl/1736.2,4.343,180.,2.2169/
      data ax,bx,cx,dx/0.5,0.25,0.4256,2.988/
C
C SOLIDUS (EQUATION 18)
C
      tsol=smb88(p)
C
C IF BELOW SOLIDUS, MELT FRACTION IN ZERO
C
      if (t.le.tsol) then
         xmb88=0.
         return
      endif
C
C LIQUIDUS (EQUATION 19)
C
      tliq=al+bl*p+cl*atan2(p,dl)
C      write(6,*)'lines:',tsol,tliq
C
C DIMENSIONLESS TEMPERATURE (EQUATION 20)
C
      td=(t-0.5*(tsol+tliq))/(tliq-tsol)
C
C MELT FRACTION (EQUATION 21)
C
      xmb88=ax+td+(td**2.-bx)*(cx+dx*td)
C         write(6,*)'xmb88:',tsol,t,tliq,td,xmb88,p
      return
      end
C
C------------------------------------------------------------------
C 
C CALCULATE VARIOUS DERIVATIVES NEEDED TO CALCULATE INSTANTANEOUS
C MELT PRODUCTION RATES
C
      function DMB88(ip)
      implicit double precision(a-h,o-z)
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
C
      parameter(ixzmax=2001)
      double precision frac(ixzmax,ixzmax),dpda(ixzmax,ixzmax),
     &     dtrda(ixzmax,ixzmax),dfda(ixzmax,ixzmax),
     &     dtuda(2,ixzmax,ixzmax),nfrac(ixzmax,ixzmax),
     &     dtsda(ixzmax,ixzmax)      
C      common /derivs/dpda,dtrda,dtdp,dtdtl,dtdtr,dtdts,dtldp,
C     &     dtrda,dtsda,dtsdp,dxdt
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
      common /dermb/dtsdac,dtdac,dtuda,dtsda
C
C PRESSURE, TEMPERATURE, WATER CONTENT, MODAL CPX INPUT HERE
C NB WATER CONTENT AND MODAL CPX NOT USED IN THIS PARAMETERISATION
C
      common /mf/p,t,xh2o,cpxm,frac,nfrac
      common /intder/dpdac,dtrdac,f,dfdac
C
C DIMENSIONLESS, SOLIDUS AND LIQUIDUS TEMPERATURE IN HERE
C
      common /tmb88/td,tsol,tliq
C
C CONSTANTS IN PARAMETERISATIONS FOR SOLIDUS, LIQUIDUS & MELT FRACTION
C
      data as,bs,cs,t0s/4.968e-4,1.2e-2,136.05,1100./
      data al,bl,cl,dl/1736.2,4.343,180.,2.2169/
      data ax,bx,cx,dx/0.5,0.253,0.4256,2.988/
C
C FIND DIMENSIONLESS, SOLIDUS AND LIQUIDUS TEMPERATURES AND MELT 
C FRACTION
C
C NEED TO RE RUN XMB88 SO THAT THE CORRECT TSOL, TLIQ AND TD 
C IS CALCULATED
C
      dmb88=XMB88(ip)
      if (f.eq.0.) then
        dmb88=0.
      endif
C      dmb88=f
C      write(6,*)'xmb88: ',dmb88
C      write(6,*)'f: ',f
C
C IF BELOW SOLIDUS, DX/Dt  MUST BE ZERO
C
      if (dmb88.le.0.) then
        dfdac=0.
C        return
        goto 2
      endif

C      write(6,*)'xmb88: ',dmb88
C      write(6,*)'f: ',f
C
C CHANGE IN MELT FRACTION WRT DIMENSIONLESS TEMPERATURE 
C (DIFFERENTIATING EQ 21)
C
      dfdt=1.+2.*cx*td+3.*dx*td**2.-bx*dx
C
C CHANGE IN LIQUIDUS TEMPEARTURE WRT PRESSURE (DIFFERENTIATING EQ 19)
C
      dtldp=bl+(cl/dl)/(1.+(p/dl)**2.)
C
C CHANGE IN DIMENSIONLESS TEMPERATURE WRT LIQUIDUS TEMPERATURE FOR CONSTANT
C SOLIDUS TEMPERATURE (DIFFERENTIATING EQ 20)
C
      deltsq=(tliq-tsol)**2.
      dtdtl=(tsol-t)/deltsq
C
C CHANGE IN SOLIDUS TEMPEARTURE WRT PRESSURE (DIFFERENTIATING EQ 18)
C
      dpdts=1./cs+as*bs*exp(bs*(tsol-t0s))
      dtsdp=1./dpdts
C
C CHANGE IN DIMENSIONLESS TEMPERATURE WRT SOLIDUS TEMPERATURE FOR CONSTANT
C LIQUIDUS TEMPERATURE (DIFFERENTIATING EQ 20)
C
      dtdts=(t-tliq)/deltsq
C
C CHANGE IN DIMENSTIONLESS TEMPERATURE WITH RESPECT TO PRESSURE AT
C CONSTANT TEMPERATURE (FROM THE EXPRESSIONS ABOVE)
C
      dtdp=dtldp*dtdtl+dtsdp*dtdts
C
C CHANGE IN DIMENSIONLESS TEMPERATURE WRT TEMPERATURE AT CONSTANT
C PRESSURE (DIFFERENTIATING EQ 20)
C
      dtdtr=1./(tliq-tsol)
C
C CHANGE IN SOLIDUS TEMPERATURE WITH RESPECT TO TIME
C
      dtsdac=dtsdp*dpdac
C
C CHANGE IN MELT FRACTION WRT TIME DX/Dt
C
C      write(6,*)'dtdp: ',dtdp,'dtdtr: ',dtdtr
C      write(6,*)'dpdac: ',dpdac,'dtrdac: ',dtrdac
C      write(6,*)'dfdt: ',dfdt,'dtsdac:',dtsdac
      dfdac=dfdt*(dpdac*dtdp+dtrdac*dtdtr)
C      write(6,*)'Melt Rate: ',dfdac
C
 2    continue
      return
      end
C
C------------------------------------------------------------------
C 
C FUNCTION TO CALCULATE SOLIDUS TEMPERATURE
C
      function SMB88(p)
      implicit double precision(a-h,o-z)
C
C CONSTANTS IN PARAMETERISATIONS FOR SOLIDUS
C
      data as,bs,cs,t0s/4.968e-4,1.2e-2,136.05,1100./
C
C SOLIDUS (EQUATION 18).  EXPRESSION HAS P ON BOTH SIDES SO IS SOLVED BY ITERATION.
C
C      write(6,*)'Entering SMB88'
      tsol=cs*p+t0s
C     write(6,*)'smb88:p,tsol',p,tsol
      n=0
 10   n=n+1
      tsmt0s=tsol-t0s
      pd=tsmt0s/cs+as*exp(bs*tsmt0s)-p
      dpodts=1./cs+as*bs*exp(bs*tsmt0s)
      dt=(-pd)/dpodts
      tsol=tsol+dt
C      write(6,*)tsol
      err=abs(dt/tsol)
      if ( (err.gt.1.e-6) .and. (n.lt.10) ) goto 10
      smb88=tsol
C      write(6,*)'Final Ts',tsol
      return
      end
C 
C------------------------------------------------------------------
C
