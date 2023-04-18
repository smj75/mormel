C
C ROUTINES USED TO CALCULATE FIRSTLY FICTIVE ELEMENTS
C AND THEN REAL ELEMENTS
C
C FIRST SUBROUTINE CALLED BY WJ.F TO CHOOSE WHICH ELEMENTS
C ARE TO BE CALCULATED.
C
      subroutine COMPSW(n1,re1,ntp1)
C
C REAL DECLARATION
C 
      implicit double precision(a-h,o-z)
C 
C 
      parameter (iemax=20)
      character*2 elabel(iemax)
      character*100 tzfile(iemax)
      common /comp/tzfile,elabel,gdepth,icomp,isourc,nel
C
C NO COMPOSITION CALCULATION
C
      if (icomp.eq.0) then
        return
C
C FICTIVE ELEMENTS
C
      elseif (icomp.eq.1) then
c        write(6,*)'Calculating Fictive Elements'
        call FICOMP(n1,re1)
C
C SPIEGELMAN TRACE ELEMENTS
C
      elseif (icomp.eq.2) then
        write(6,*)'Calculating Spiegelman trace elements, variable D'
        write(6,*)'Routine not completed: stopping'
        stop
        call SPCOMP(n1,re1)
C
C REAL ELEMENTS USING LOOK-UP TABLES EXPORTED FROM INVMEL
C
      elseif (icomp.eq.3) then
        call CAINT(n1,ntp1)
C
C INVALID MELTING MODEL
C
      else
        write(6,*)'Invalid option for icomp'
        write(6,*)'Choose 0,1,2'
        write(6,*)'STOP PROGRAM'
        stop
      endif
C
      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C FICOMP - FICTIVE COMPOSITION
C
      subroutine FICOMP(n1,re1)
C
C REAL DECLARATION
C 
      implicit double precision(a-h,o-z)

C PARAMETER SETTING FOR ARRAY SIZES
C MAX NUMBER OF GRID NODES IN ONE DIRECTION
C
      parameter (ixzmax=2001)
C
C PARAMETER SETTING FOR FICTIVE ELEMENTS
C BULK D
C
      parameter (del1=5.,del2=1.,del3=0.5,del4=0.1,del5=0.01,
     &           del6=0.001)
C
C DEPTH FOR SPINEL-GARNET TRANSITION
C
      
C
C PARAMTER SETTING FOR ORIGINAL CONCENTRATIONS
C SET TO 1 FOR ALL FICTIVE ELEMENTS TO START WITH
C
      parameter (c01=1.)
C
C GRID
C
      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
C
C MELTING
C
      double precision frac(ixzmax,ixzmax),nfrac(ixzmax,ixzmax)
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C
C RESIDUE
C RES - ARRAYS REQUIRED TO CALCULATE THE NEW RESIDUE STRUCTURE
C REST - STARTING RESIDUE (ZERO TO START WITH)
C
      double precision res(3,ixzmax,ixzmax),rest(ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C PARAMETERS
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &        texpf,texps
C
C TIMING
C
      common /calc/tend,icalc
C
C COMMON FOR MELT RATE
C
      double precision dpda(ixzmax,ixzmax)
      double precision dtrda(ixzmax,ixzmax)
      double precision dfda(ixzmax,ixzmax)
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
C
C COMMON FOR FICTIVE COMPOSITIONS: 
C FCI = FICTIVE COMPOSITION - INSTANTANEOUS
C FCA = ACCUMMULATED INSTANTANEOUS COMPOSITION
C SCI = SOLID COMPOSITION
C
      double precision fci(6,ixzmax,ixzmax),fca(6,ixzmax,ixzmax),
     &     tcomp(6),carea(6),sci(6,ixzmax,ixzmax)
      common /fakecp/fci,fca,sci,tcomp,carea
C
C COMMUNICATION OF R (INTEGRATED MELT PRODUCTION RATE)
C
      common /intgra/tarea,tcrust
C
C FIRST CALCULATE THE SOLID COMPOSITION
C
        do ix=1,nx
          do iz=1,nz
            sci(1,ix,iz)=c01*(1.-frac(ix,iz))**((1./del1)-1.)
            sci(2,ix,iz)=c01*(1.-frac(ix,iz))**((1./del2)-1.)
            sci(3,ix,iz)=c01*(1.-frac(ix,iz))**((1./del3)-1.)
            sci(4,ix,iz)=c01*(1.-frac(ix,iz))**((1./del4)-1.)
            sci(5,ix,iz)=c01*(1.-frac(ix,iz))**((1./del5)-1.)
            sci(6,ix,iz)=c01*(1.-frac(ix,iz))**((1./del6)-1.)
          enddo
        enddo
C
C CALCULATE THE INSTANTANEOUS MELT COMPOSITION
C EQUATION FOR INSTANTANEOUS MELT CONCENTRATION: 
C CL/C0 = 1/D(1-X)**(1/D-1)
C
      do ix=1,nx
       do iz=1,nz
          if(nfrac(ix,iz).le.0.) then
C          write(6,*) 'NO MELTING SO NO COMPOSITION! DERR'
          fci(1,ix,iz)=0.
          fci(2,ix,iz)=0.
          fci(3,ix,iz)=0.
          fci(4,ix,iz)=0.
          fci(5,ix,iz)=0.
          fci(6,ix,iz)=0.
          goto 12
          endif
C IF STATEMENTS TO MAKE CL=O IF THE SOURCE HAS RUN OUT
C TOLERANCE SET TO 1.E-10
          if(sci(1,ix,iz).lt. 1.e-10) then
          fci(1,ix,iz)=0.
          goto 12
          endif
          if(sci(2,ix,iz).lt. 1.e-10) then
          fci(2,ix,iz)=0.
          goto 12
          endif
          if(sci(3,ix,iz).lt. 1.e-10) then
          fci(3,ix,iz)=0.
          goto 12
          endif
          if(sci(4,ix,iz).lt. 1.e-10) then
          fci(4,ix,iz)=0.
          goto 12
          endif
          if(sci(5,ix,iz).lt. 1.e-10) then
          fci(5,ix,iz)=0.
          goto 12
          endif
          if(sci(6,ix,iz).lt. 1.e-10) then
          fci(6,ix,iz)=0.
          goto 12
          endif
          fci(1,ix,iz)=c01*(1.-frac(ix,iz))**(1./del1)
     &                 /(del1*(1.-frac(ix,iz)))
          fci(2,ix,iz)=c01*(1.-frac(ix,iz))**(1./del2)
     &                 /(del2*(1.-frac(ix,iz)))
          fci(3,ix,iz)=c01*(1.-frac(ix,iz))**(1./del3)
     &                 /(del3*(1.-frac(ix,iz)))
          fci(4,ix,iz)=c01*(1.-frac(ix,iz))**(1./del4)
     &                 /(del4*(1.-frac(ix,iz)))
          fci(5,ix,iz)=c01*(1.-frac(ix,iz))**(1./del5)
     &                 /(del5*(1.-frac(ix,iz)))
          fci(6,ix,iz)=c01*(1.-frac(ix,iz))**(1./del6)
     &                 /(del6*(1.-frac(ix,iz)))
 12       continue
       enddo
      enddo
C
C CALCULATED THE ACCUMULATED COMPOSITION BY INTEGRATING
C THE CL CURVE AND AVERGAING OVE RHTE MELT FRACTION IN
C THE VERTICAL DIRECTION
C
      do ix=1,nx
         do iz=1,nz
            if(nfrac(ix,iz).le.0.) then
             fca(1,ix,iz)=0.
             fca(2,ix,iz)=0.
             fca(3,ix,iz)=0.
             fca(4,ix,iz)=0.
             fca(5,ix,iz)=0.
             fca(6,ix,iz)=0.
C             write(6,*)'nfrac is zero'
             goto 14
            endif
            if(iz.eq.nz) then
             fca(1,ix,iz)=fci(1,ix,iz)
             fca(2,ix,iz)=fci(2,ix,iz)
             fca(3,ix,iz)=fci(3,ix,iz)
             fca(4,ix,iz)=fci(4,ix,iz)
             fca(5,ix,iz)=fci(5,ix,iz)
             fca(6,ix,iz)=fci(6,ix,iz)
C             write(6,*)'at bottom:',z(iz)*1.e-3
             goto 14
            endif
C INTEGRATION FOR THE FIRST ELEMENT:
C DELF = DIFFERENCE IN F WHICH GOES ON THE BOTTOM AND IS
C SAME FOR ALL ELEMENTS
C 
          delf=frac(ix,iz)-frac(ix,iz+1)
C          write(6,*)'delf=',delf
C          if(delf.le.0.) then
C            write(6,*)z(iz)*1.e-3,z(iz+1)*1.e-3,'nfrac1:',nfrac(ix,iz),
C     &                nfrac(ix,iz+1),dfda(ix,iz),dfda(ix,iz+1)
C          endif
C INTEGRATION OF CL, FOR F1,ELEMENT 1
          cli1a=(-1.)*(((1.-frac(ix,iz+1))**(1./del1))/(1./del1))
C          write(6,*)'cli1a=',cli1a
C INTEGRATION OF CL FOR F2, ELEMENT 1
          cli1b=(-1.)*(((1.-frac(ix,iz))**(1./del1))/(1./del1))
C          write(6,*)'cli1b=',cli1b
C TOTAL FOR TOP INTEGRAL,ELEMENT 1:
          cltop=(c01/del1)*(cli1b-cli1a)
C          write(6,*)'cltop=',cltop
C INTEGRATION OF ELEMENT 1 COMPLETE:
          if(delf.eq.0.) then
            fca(1,ix,iz)=cltop
C          elseif (delf.lt.0.) then
C            fca(1,ix,iz)=cltop/(delf*-1.)
          else
          fca(1,ix,iz)=cltop/delf
          endif
C          write(6,*)'fci:',fci(1,ix,iz),'fca:',fca(1,ix,iz)
C
C ELEMENT 2
C 
          cli2a=(-1.)*(((1.-frac(ix,iz+1))**(1./del2))/(1./del2)) 
          cli2b=(-1.)*(((1.-frac(ix,iz))**(1./del2))/(1./del2))
          cltop=(c01/del2)*(cli2b-cli2a)
          if(delf.eq.0.) then
            fca(2,ix,iz)=cltop
C          elseif (delf.lt.0.) then
C            fca(2,ix,iz)=cltop/(delf*-1.)
          else
          fca(2,ix,iz)=cltop/delf
          endif
C          fca(2,ix,iz)=cltop/delf
C 
C ELEMENT 3
C 
          cli3a=(-1.)*(((1.-frac(ix,iz+1))**(1./del3))/(1./del3)) 
          cli3b=(-1.)*(((1.-frac(ix,iz))**(1./del3))/(1./del3))
          cltop=(c01/del3)*(cli3b-cli3a)
          if(delf.eq.0.) then
            fca(3,ix,iz)=cltop
C          elseif (delf.lt.0.) then
C            fca(3,ix,iz)=cltop/(delf*-1.)
          else
          fca(3,ix,iz)=cltop/delf
          endif
C          fca(3,ix,iz)=cltop/delf
C
C ELEMENT 4
C
          cli4a=(-1.)*(((1.-frac(ix,iz+1))**(1./del4))/(1./del4)) 
          cli4b=(-1.)*(((1.-frac(ix,iz))**(1./del4))/(1./del4))
          cltop=(c01/del4)*(cli4b-cli4a)
          if(delf.eq.0.) then
            fca(4,ix,iz)=cltop
C          elseif (delf.lt.0.) then
C            fca(4,ix,iz)=cltop/(delf*-1.)
          else
          fca(4,ix,iz)=cltop/delf
          endif
C          fca(4,ix,iz)=cltop/delf
C
C ELEMENT 5
C
          cli5a=(-1.)*(((1.-frac(ix,iz+1))**(1./del5))/(1./del5)) 
          cli5b=(-1.)*(((1.-frac(ix,iz))**(1./del5))/(1./del5))
          cltop=(c01/del5)*(cli5b-cli5a)
          if(delf.eq.0.) then
            fca(5,ix,iz)=cltop
C          elseif (delf.lt.0.) then
C            fca(5,ix,iz)=cltop/(delf*-1.)
          else
          fca(5,ix,iz)=cltop/delf
          endif
C          fca(5,ix,iz)=cltop/delf
C
C ELEMENT 6
C
          cli6a=(-1.)*(((1.-frac(ix,iz+1))**(1./del6))/(1./del6)) 
          cli6b=(-1.)*(((1.-frac(ix,iz))**(1./del6))/(1./del6))
          cltop=(c01/del6)*(cli6b-cli6a)
          if(delf.eq.0.) then
            fca(6,ix,iz)=cltop
C          elseif (delf.lt.0.) then
C            fca(6,ix,iz)=cltop/(delf*-1.)
          else
          fca(6,ix,iz)=cltop/delf
          endif
C          fca(6,ix,iz)=cltop/delf
C
C
 14       continue
         enddo
       enddo
C
C
C SET AREAS TO ZERO FOR NEXT CALCULATION
C
      do ic=1,6
         carea(ic)=0.
      enddo
C
C
C INTEGRATE THE INSTANTANEOUS MELT COMPOSITION TO FIND
C FCA - THE AVERAGE MELT COMPOSITION USING THE FOLLOWING
C EQUATION: FCA = 1/R (INTEGRAL(CL*DX/Dt)) 
C
      carea1 = 0.0
      do ic=1,6
C        do ix=1,nx
C         do iz=1,nz-1
        do iz=1,nz-1
          do ix=1,nx
C
C TRY AND PUT FRAC AS NFRAC AND SEE WHAT HAPPENS: 
C DOESN'T MAKE MUCH DIFFERENCE TO FINAL ANSWER
C PLUS THIKNK IT MIGHT MEAN THIS ANYWAY
C         fr1=dfda(ix,iz)/((1.-frac(ix,iz))**2.)
C         fr2=dfda(ix,iz+1)/((1.-frac(ix,iz+1))**2.)
C CHANGED FCI TO FCA TO TRY AND REGISTER THE TINIEST MELT
C FRACTIONS IN THE CALUCLATION WHERE THEY ARE MISSED AT THE
C EDGES OF THE MELTING REGION
C
C         write(6,*)ic,fca(ic,ix,iz),res(re1,ix,iz),dfda(ix,iz),
C     &   nfrac(ix,iz)
C
            if (dfda(ix,iz).le.0.) then
              cr1=0.
            else
              cr1=fca(ic,ix,iz)*dfda(ix,iz)/((res(re1,ix,iz)
     &        -nfrac(ix,iz))**2.)
            endif
            if (dfda(ix,iz+1).le.0.) then
              cr2=0.
            else
              cr2=fca(ic,ix,iz+1)*dfda(ix,iz+1)/((res(re1,ix,iz+1)
     &         -nfrac(ix,iz+1))**2.)
            endif
            carea2=(cr1+cr2)*dz/2.
C            write(6,*)'carea2',carea2
c            if (ix.eq.1) then
c              carea1=carea2
c              goto 13
c            elseif(ix.gt.1) then
c              carea(ic)=carea(ic)+(carea1+carea2)*dx/2.
cC              write(6,*)'Totalarea:',carea(ic)
c              carea1=carea2
c            endif
            if (ix.gt.1)
     &        carea(ic)=carea(ic)+(carea1+carea2)*dx/2.
C              write(6,*)'Totalarea:',carea(ic)
            carea1=carea2
          enddo
        enddo
      enddo
C
C CALCULATE FINAL COMPOSITION BY DIVIDING INTEGRAL PRESUMABLY INCLUDE
C THE DENSITY CONTRAST TOO - makes no difference as they cancel out??
C
       do ic=1,6
           if(tarea.eq.0.) then
           tcomp(ic)=0.
           else
C          tcomp(ic)=rhos/rhom*carea(ic)/tarea
           tcomp(ic)=carea(ic)/tarea
           endif
       enddo
C
C END SUBROUTINE
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SPCOMP - SPIEGELMAN COMPOSITION
C 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
C      subroutine SPCOMP(n1,re1)
CC
CC REAL DECLARATION
CC 
C      implicit double precision*8(a-h,o-z)
C
CC PARAMETER SETTING FOR ARRAY SIZES
CC MAX NUMBER OF GRID NODES IN ONE DIRECTION
CC
C      parameter (ixzmax=2001)
CC
CC PARAMETER SETTING FOR RARE EARTH ELEMENTS, SET UP dlag (D La Garnet)
CC
CC BULK D GARNET
CC
C      parameter (dlag=0.0057, dceg=0.0105, dprg=0.0192, dndg=0.0285,
C     &           dsmg=0.0515, deug=0.0705, dgdg=0.0962, dtbg=0.1341, 
C     &           ddyg=0.1811, dhog=0.2487, derg=0.3170, dtmg=0.4634,
C     &           dybg=0.6140, dlug=0.8295)
CC
CC BULK D SPINEL
CC
C      parameter (dlas=0.0065, dces=0.0112, dprs=0.0171, dnds=0.0237,
C     &           dsms=0.0298, deus=0.0358, dgds=0.0356, dtbs=0.0374, 
C     &           ddys=0.0404, dhos=0.0395, ders=0.0396, dtms=0.0416,
C     &           dybs=0.0432, dlus=0.0464)
CC
CC 
CC COMMON FOR GDEPTH (COMPOSITION)
c      common /comp/tzfile,elabel,gdepth,icomp,isourc,nel
CC
CC PARAMTER SETTING FOR ORIGINAL CONCENTRATIONS
CC McKENZIE AND O'NIONS
C      if(isourc=0) then
C      parameter (c0la=0.206, c0ce=0.722, c0pr=0.143, c0nd=0.815, 
C     &           c0sm=0.299, c0eu=0.115, c0gd=0.419, c0tb=0.077, 
C     &           c0dy=0.525, c0ho=0.120, c0er=0.347, c0tm=0.054,
C     &           c0yb=0.347, c0lu=0.057)   
C      else
C      write(6,*)'Currently invalid option, change isourc to 0'
C      stop
C      endif
CC
CC GRID
CC
C      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
C      common /grid/xmax,zmax,dx,dz
C      common /grido/pres,x,z,nx,nz
CC
CC MELTING
CC
C      double precision frac(ixzmax,ixzmax),nfrac(ixzmax,ixzmax)
C      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
C      common /mf/p,t,xh2o,cpxm,frac,nfrac
CC
CC RESIDUE
CC RES - ARRAYS REQUIRED TO CALCULATE THE NEW RESIDUE STRUCTURE
CC REST - STARTING RESIDUE (ZERO TO START WITH)
C      double precision res(3,ixzmax,ixzmax),rest(ixzmax,ixzmax)
C      common /resdue/rest,res,resc
CC
CC PARAMETERS
CC
C      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
C     &        texpf,texps
CC
CC TIMING
CC
C      common /calc/tend,icalc
CC
CC COMMON FOR MELT RATE
CC
C      double precision dpda(ixzmax,ixzmax),dtrda(ixzmax,ixzmax),dfda(ixzmax,ixzmax)
C      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
C     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
CC
CC COMMON FOR REE COMPOSITIONS: 
CC LAI = La INSTANTANEOUS
CC LAA = La ACCUMMULATED INSTANTANEOUS COMPOSITION
CC LAS = La SOLID COMPOSITION
CC LAMELT = AVERAGE LA COMPOSITION IN MELT
CC LAINT = LA INTEGRAL
CC
C      double precision lai(ixzmax,ixzmax),laa(ixzmax,ixzmax),
C     &     lamelt,laint,las(ixzmax,ixzmax),cgla(3,ixzmax,izmax),
C     &     cglas(ixzmax,ixzmax)
C      common /lanth/lai,laa,las,lamelt,laint,cgla,cglas
CC
CC COMMUNICATION OF R (INTEGRATED MELT PRODUCTION RATE)
C      common /intgra/tarea,tcrust
C
CC
CC FIRST CALCULATE THE SOLID COMPOSITION
CC CGLA IS AN ARRA
CC
C        do ix=1,nx
C          do iz=1,nz
C           if(z(iz).gt.gdepth) then
C            las(ix,iz)=cgla(com1,ix,iz)*c0la*
C     &                 (1.-frac(ix,iz))**((1./dlag)-1.)
C           else
C            las(ix,iz)=cgla(com1,ix,iz)*c0la*
C     &                 (1.-frac(ix,iz))**((1./dlas)-1.)
C          enddo
C        enddo
CC
CC CALCULATE THE INSTANTANEOUS MELT COMPOSITION
CC EQUATION FOR INSTANTANEOUS MELT CONCENTRATION: 
CC CL/C0 = 1/D(1-X)**(1/D-1)
CC
C      do ix=1,nx
C       do iz=1,nz
C          if(nfrac(ix,iz).le.0.) then
CC          write(6,*) 'NO MELTING SO NO COMPOSITION! DERR'
C          lai(ix,iz)=0.
C          goto 12
C          endif
CC IF STATEMENTS TO MAKE CL=O IF THE SOURCE HAS RUN OUT
CC TOLERANCE SET TO 1.E-10
C          if(las(ix,iz).lt. 1.e-10) then
C          lai(ix,iz)=0.
C          goto 12
C          endif
C          if(z(iz).gt.gdepth) then
C          lai(ix,iz)=mfla(com1,ix,iz)*c0la*(1.-frac(ix,iz))**(1./dlag)
C     &                 /(dlag*(1.-frac(ix,iz)))
C          else
C          lai(ix,iz)=mfla(com1,ix,iz)*c0la*(1.-frac(ix,iz))**(1./dlas)
C     &                 /(dlas*(1.-frac(ix,iz)))
C 12       continue
C       enddo
C      enddo
CC
CC CALCULATED THE ACCUMULATED COMPOSITION BY INTEGRATING
CC THE CL CURVE AND AVERGAING OVE RHTE MELT FRACTION IN
CC THE VERTICAL DIRECTION
CC
C      do ix=1,nx
C         do iz=1,nz
C            if(nfrac(ix,iz).le.0.) then
C             laa(ix,iz)=0.
CC             write(6,*)'nfrac is zero'
C             goto 14
C            endif
C            if(iz.eq.nz) then
C             laa(ix,iz)=lai(ix,iz)
CC             write(6,*)'at bottom:',z(iz)*1.e-3
C             goto 14
C            endif
CC INTEGRATION FOR THE FIRST ELEMENT:
CC DELF = DIFFERENCE IN F WHICH GOES ON THE BOTTOM AND IS
CC SAME FOR ALL ELEMENTS
CC 
C          delf=frac(ix,iz)-frac(ix,iz+1)
CC          write(6,*)'delf=',delf
C          if(delf.le.0.) then
C            write(6,*)z(iz)*1.e-3,z(iz+1)*1.e-3,'nfrac1:',nfrac(ix,iz),
C     &                nfrac(ix,iz+1),dfda(ix,iz),dfda(ix,iz+1)
C          endif
CC INTEGRATION OF CL, FOR F1,ELEMENT 1
C          if(z(iz).gt.gdepth) then
C          cli1a=-1.*(mfla(com1,ix,iz)*c0la/dalg)
C     &          *(((1.-frac(ix,iz+1))**(1./dlag))/(1./dlag))
CC          write(6,*)'cli1a=',cli1a
C          else
C          cli1a=-1.*(mfla(com1,ix,iz)*c0la/dals)
C     &          *(((1.-frac(ix,iz+1))**(1./dlas))/(1./dlas))
C          endif
CC INTEGRATION OF CL FOR F2, ELEMENT 1
C          if(z(iz).gt.gdepth) then
C          cli1b=-1.*(mfla(com1,ix,iz)*c0la/dalg)
C     &          *(((1.-frac(ix,iz))**(1./dlag))/(1./dlag))
CC          write(6,*)'cli1b=',cli1b
C          else
C          cli1b=-1.*(mfla(com1,ix,iz)*c0la/dals)
C     &          *(((1.-frac(ix,iz))**(1./dlas))/(1./dlas))
CC TOTAL FOR TOP INTEGRAL,ELEMENT 1:
C          cltop=cli1b-cli1a
CC          write(6,*)'cltop=',cltop
CC INTEGRATION OF ELEMENT 1 COMPLETE:
C          laa(ix,iz)=cltop/delf
CC          write(6,*)'La-inst:',lai(ix,iz),'La-acc:',laa(ix,iz)
CC
CC
CC
C 14       continue
C         enddo
C       enddo
CC
CC
CC SET AREAS TO ZERO FOR NEXT CALCULATION
CC
C      do ic=1,6
C         carea(ic)=0.
C      enddo
CC
CC
CC INTEGRATE THE INSTANTANEOUS MELT COMPOSITION TO FIND
CC FCA - THE AVERAGE MELT COMPOSITION USING THE FOLLOWING
CC EQUATION: FCA = 1/R (INTEGRAL(CL*DX/Dt)) 
CC
C      do ic=1,6
C      do ix=1,nx
C         do iz=1,nz-1
CC TRY AND PUT FRAC AS NFRAC AND SEE WHAT HAPPENS: 
CC DOESN'T MAKE MUCH DIFFERENCE TO FINAL ANSWER
CC PLUS THIKNK IT MIGHT MEAN THIS ANYWAY
CC         fr1=dfda(ix,iz)/((1.-frac(ix,iz))**2.)
CC         fr2=dfda(ix,iz+1)/((1.-frac(ix,iz+1))**2.)
CC CHANGED FCI TO FCA TO TRY AND REGISTER THE TINIEST MELT
CC FRACTIONS IN THE CALUCLATION WHERE THEY ARE MISSED AT THE
CC EDGES OF THE MELTING REGION
CC
C         cr1=fca(ic,ix,iz)*dfda(ix,iz)/((res(re1,ix,iz)
C     &        -nfrac(ix,iz))**2.)
C         cr2=fca(ic,ix,iz+1)*dfda(ix,iz+1)/((res(re1,ix,iz+1)
C     &         -nfrac(ix,iz+1))**2.)
C         carea2=(cr1+cr2)*dz/2.
C         if(ix.eq.1) then
C         carea1=carea2
C         goto 13
C         elseif(ix.gt.1) then
C         carea(ic)=carea(ic)+(carea1+carea2)*dx/2.
CC         write(6,*)'Totalarea:',carea(ic)
C         carea1=carea2
C         endif
C 13      continue
C         enddo
C      enddo
C      enddo
CC
CC CALCULATE FINAL COMPOSITION BY DIVIDING INTEGRAL PRESUMABLY INCLUDE
CC THE DENSITY CONTRAST TOO - makes no difference as they cancel out??
CC
C       do ic=1,6
C          tcomp(ic)=rhos/rhom*carea(ic)/tarea
C       enddo
CC
CC END SUBROUTINE
CC
C      return
C      end
C
C--------------------------------------------------------------
C
C CALCULATE AVERAGE MELT PRODUCTION RATE 
C AND AVERAGE MELTING RATE AS A FUNCTION OF DEPTH
C BY INTEGRATION OVER MELTING REGION
C
      subroutine FZ(n1,re1)
C
C DOUBLE PRECISION DECLARATION
C 
      implicit double precision(a-h,o-z)
C
C PARAMETER SETTING FOR ARRAY SIZES
C MAX NUMBER OF GRID NODES IN ONE DIRECTION
C
      parameter (ixzmax=2001)
C
C GRID
C
      double precision pres(ixzmax)
      double precision x(ixzmax)
      double precision z(ixzmax)
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
C
C MELTING
C
      double precision frac(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C
C RESIDUE
C RES - ARRAYS REQUIRED TO CALCULATE THE NEW RESIDUE STRUCTURE
C REST - STARTING RESIDUE (ZERO TO START WITH)
C
      double precision res(3,ixzmax,ixzmax)
      double precision rest(ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C PARAMETERS
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &        texpf,texps
C
C TIMING
C
      common /calc/tend,icalc
C
C COMMON FOR MELT RATE
C
      double precision dpda(ixzmax,ixzmax)
      double precision dtrda(ixzmax,ixzmax)
      double precision dfda(ixzmax,ixzmax)
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr

      common /countr/itym,ntime,iendpr
C
C COMMUNICATION OF R (INTEGRATED MELT PRODUCTION RATE)
C
      common /intgra/tarea,tcrust
C
C COMMUNICATION OF DEPTH DEPENDENT X
C
      double precision fidarr(ixzmax)
      double precision fiarea(ixzmax)
      double precision firate(ixzmax)
      double precision ftdarr(ixzmax)
      double precision ftarea(ixzmax)
      double precision ftrate(ixzmax)
      common /fidepth/fidarr,fiarea,firate
      common /ftdepth/ftdarr,ftarea,ftrate
C
C       /
C       |            dF
C       |   F(x , z) -- dA
C       |            dt
C       / A
C  C =  -------------------
C            /
C            |   dF
C            |   -- dA
C            |   dt
C            / A
C
C A: CALCULATE MEAN F FOR INSTANTANEOUS F.
C
C 1. CALCULATE THE TOP INTEGRAL F*DFDA
C
      do iz=1,nz
        fiara1=0.
        do ix=1,nx-1
          if (nfrac(ix,iz).le.0. .and. dfda(ix,iz).gt.0.) then
            write(6,*)'ERROR: Melting rate not zero when finstant 0'
            stop
          endif
          fir1=0.
          fir1=nfrac(ix,iz)*dfda(ix,iz)
          fir2=nfrac(ix+1,iz)*dfda(ix+1,iz)
          fiara1=fiara1+((fir1+fir2)*dx/2.)
        enddo
        fiarea(iz)=fiara1
      enddo
C
C 2. CALCULATE THE BOTTOM INTEGRAL OF MELT RATE
C
      do iz=1,nz
         firat1=0.
         do ix=1,nx-1
            fira1=dfda(ix,iz)
            fira2=dfda(ix+1,iz)
            firat1=firat1+((fira1+fira2)*dx/2.)
         enddo
         firate(iz)=firat1
      enddo
C
C 3. CALCULATE THE TOTAL ITEGRAL
C
      do iz=1,nz
        if(firate(iz).le.0.) then
        fidarr(iz)=0.
        else
        fidarr(iz)=fiarea(iz)/firate(iz)
        endif
      enddo
C
C B: CALCULATE MEAN F FOR TOTAL F, WHEN MELTING RATE IS NOT ZERO
C 1. CALCULATE THE TOP INTEGRAL F*DFDA
C
      do iz=1,nz
         ftara1=0.
         do ix=1,nx-1
            ftr1=0.
            if (nfrac(ix,iz).le.0.) then
              ftr1=0.
            else
              ftr1=frac(ix,iz)*dfda(ix,iz)
            endif
            ftr2=0.
            if (nfrac(ix+1,iz).le.0.) then
              ftr2=0.
            else
              ftr2=frac(ix+1,iz)*dfda(ix+1,iz)
            endif
C            write(6,*)ix,iz
C            write(6,*)'frac ix',frac(ix,iz),
C     &                'frac ix+1',frac(ix+1,iz)
C            write(6,*)'nfrac ix',nfrac(ix,iz),
C     &                'nfrac ix+1',nfrac(ix+1,iz)
C            if(itym.eq.4 .and. iz.gt.10) then
C            stop
C            endif
            ftara1=ftara1+((ftr1+ftr2)*dx/2.)
         enddo
         ftarea(iz)=ftara1
      enddo
C
C 2. CALCULATE THE BOTTOM INTEGRAL OF MELT RATE
C
      do iz=1,nz
        ftrat1=0.
        do ix=1,nx-1
          ftra1=dfda(ix,iz)
          ftra2=dfda(ix+1,iz)
          ftrat1=ftrat1+((ftra1+ftra2)*dx/2.)
        enddo
        ftrate(iz)=ftrat1
      enddo
C
C 3. CALCULATE THE TOTAL ITEGRAL
C
      do iz=1,nz
        if (ftrate(iz).eq.0.) then
          ftdarr(iz)=0.
        else
          ftdarr(iz)=ftarea(iz)/ftrate(iz)
        endif
      enddo
C
C
      return
      end
C
C--------------------------------------------------------------
C
C CALCULATE POINT-AND-DEPTH AVERAGE MELTING COMPOSITION
C BY INTEGRATING OVER MELTING REGION AND USING LOOK-UP
C TABLE OF INSTANTANEOUS MELT COMPOSITIONS
C
      subroutine CAINT(n1,ntp1)
C
C DOUBLE PRECISION DECLARATION
C 
      implicit double precision(a-h,o-z)
C
C PARAMETER SETTING FOR ARRAY SIZES
C MAX NUMBER OF GRID NODES IN ONE DIRECTION
C MAX NUMBER OF ELEMENTS IN COMPOSITION CALCULATION
C
      parameter (ixzmax=2001)
      parameter (iemax=20)
C
C GRID
C
      double precision pres(ixzmax)
      double precision x(ixzmax)
      double precision z(ixzmax)
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
C
C PARAMETERS
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &        texpf,texps
C
C TIMING
C
      common /calc/tend,icalc
      common /countr/itym,ntime,iendpr
C
C MELT RATE
C
      double precision dpda(ixzmax,ixzmax)
      double precision dtrda(ixzmax,ixzmax)
      double precision dfda(ixzmax,ixzmax)
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
C
C COMMUNICATION OF R (INTEGRATED MELT PRODUCTION RATE)
C
      common /intgra/tarea,tcrust
C
C TEMPERATURE FIELD
C
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      common /temo/tem,tems
C
C POTENTIAL TEMPERATURE
C  TP, POTENTIAL TEMPERATURE
C  TPS, STARTING POT TEMP
C
      double precision tp(3,ixzmax,ixzmax)
      double precision tps(ixzmax,ixzmax)
      common /pottem/tp,tps
C
C INSTANTANEOUS COMPOSITION DATA
C
      character*2 elabel(iemax)
      character*100 tzfile(iemax)
      double precision c(iemax,ixzmax,ixzmax)
      double precision ctp(ixzmax)
      double precision cz(ixzmax)
      common /cinst/c,cz,ctp,ncz,nctp
      common /comp/tzfile,elabel,gdepth,icomp,isourc,nel
C
C STOREAGE FOR POINT AND DEPTH AVERAGES
C
      double precision cpda(iemax)
      common /cav/cpda
C
C EQUATION TO BE SOLVED
C
C        /
C        |               dF
C        |   c( T_p, z ) -- dA
C        |               dt
C        / A
C  C =  -----------------------
C             /
C             |   dF
C             |   -- dA
C             |   dt
C             / A
C
C
C INSTANTANEOUS MELT COMPOSITIONS c(T_p,z) ARE
C INTERPOLATED FORM THE LOOK-UP TABLE IN 
C FUNCTION CINTRP
C
C POTENTIAL TEMPERATURE IS ESTIMATED FROM THE
C UNCORRECTED TEMPERATURE USING THE ADIABATIC GRADIENT
C EQUIVALENT TO ASSUMING THAT CONDUCTION IS NEGLIGABLE
C
C DEPTH IN THE LOOK-UP TABLE IS IN KM +VE DOWNWARDS
C POTENTIAL TEMPERATURE IN THE LOOK-UP TABLE IS IN C
C
C
C
C FOR EVERY ELEMENT FOR WHICH INSTANTANEOUS COMPS ARE
C SPECIFIED
C
      do iel=1,nel
C
C LOOP OVER DEPTH
C
        cfxz=0.0
        fxz=0.0
        cfx0=0.0
        fx0=0.0
        do iz=1,nz
C
C DEPTH
C
          zc=(-z(iz))*1.0e-3
C
C LOOP OVER X
C
          cfx1=0.0
          fx1=0.0
          dfda0=0.0
          cinst0=0.0
          do ix=1,nx
C
C MELTING RATE
C
            dfda1=max(dfda(ix,iz),0.0)
C
C POTENTIAL TEMPERATURE
C
            tpc1=tp(ntp1,ix,iz)
c            write(6,*)'CAINT: tp = ',tpc1
c            tpc1=1300.0
C
C GUARD AGAINST EXTRAPOLATION
C
            if (dfda1.gt.0.0) then
              if (zc.lt.cz(1) .or. zc.gt.cz(ncz)) then
                write(6,*)'ERROR: Extrapolation in CAINT'
                write(6,*)'Min & Max z: ',cz(1),cz(ncz)
                write(6,*)'Required z: ',zc
                return
              endif
              if (tpc1.lt.ctp(1) .or. tpc1.gt.ctp(nctp)) then
                write(6,*)'ERROR: Extrapolation in CAINT'
                write(6,*)'x,z: ',x(ix),z(iz)
                write(6,*)'Min & Max Tp: ',ctp(1),ctp(nctp)
                write(6,*)'Required Tp: ',tpc1
                return
              endif
            endif
C
C INSTANTANEOUS MELT COMPOSITION
C
            cinst1=CINTRP(tpc1,zc,iel,itym)
c            write(6,*)x(ix)*1.0e-3,z(iz)*1.0e-3,cinst1
C
C INTEGRATE OVER X
C
            if (ix.gt.1) then
              cfx1=cfx1+0.5*dx*(cinst0*dfda0+cinst1*dfda1)
              fx1=fx1+0.5*dx*(dfda0+dfda1)
            endif
            cinst0=cinst1
            dfda0=dfda1
          enddo
C
C INTEGRATE OVER Z
C
          if (iz.gt.1) then
            cfxz=cfxz+0.5*dz*(cfx0+cfx1)
            fxz=fxz+0.5*dz*(fx0+fx1)
          endif
          cfx0=cfx1
          fx0=fx1
        enddo
C
C POINT AND DEPTH AVERAGE USING EQUATION ABOVE
C
        if (fxz.le.0.0) then
          cpda(iel)=0.0
        else
          cpda(iel)=cfxz/fxz
        endif
C
C END OF LOOP FOR EACH ELEMENT
C
 110    format(a2,1x,'P&D av: ',f10.4,5x,'Melt rate: ',f10.4)
c        write(6,110)elabel(iel),cfxz,fxz     
      enddo
      return
      end
C
C--------------------------------------------------------------
C
CC 
CC SUBROUTINE DRAFT FOR TRACKING 2-LAYER SOURCE COMPOSITION ROUNDS
CC AROUND THE GRID
CC
CC
CC
C      subroutine COMFVI(re0,re1)
C 1. CALCULATES RESIDUE "STRUCTURE" USING ADVECTION ONLY
CC PARTIAL DIFFERENTIAL EQUATION:
CC dR/dt = -VXdR/dX -VZdR/dZ
CC
CC 2. 
CC
CC REAL DECLARATION
C      implicit double precision(a-h,o-z)
CC
CC PARAMETER SETTING FOR ARRAY SIZES
CC MAX NUMBER OF GRID NODES IN ONE DIRECTION
CC
C      parameter (ixzmax=2001)
CC
CC ARRAY FOR TRIDIAGONAL MATRIX
CC
C      double precision a(ixzmax),b(ixzmax),c(ixzmax),r(ixzmax),u(ixzmax)
CC
CC GRID
CC
C      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
C      common /grid/xmax,zmax,dx,dz
C      common /grido/pres,x,z,nx,nz
CC
CC MELTING
CC
C C     common /melt/ds,xh2o,cpxm,iparam,imelt,itemp
CC
CC RESIDUE
CC RES - ARRAYS REQUIRED TO CALCULATE THE NEW RESIDUE STRUCTURE
CC REST - STARTING RESIDUE (ZERO TO START WITH)
C      double precision res(3,ixzmax,ixzmax),rest(ixzmax,ixzmax)
C      common /resdue/rest,res,resc
CC
CC VELOCITY GRID ON ORIGINAL NODES
CC
C      double precision vx(ixzmax,ixzmax),vz(ixzmax)
C      common /velo/vx,vz,uhalf,vmax
C      common /vel/ufull,alpha,viss,crust,ialpha
CC
CC VELOCITY ON STAGGERED GRID
CC
C      double precision vxs(ixzmax,ixzmax),vzs(ixzmax,ixzmax)
C      common /velstg/vxs,vzs
CC
CC PARAMETERS
CC
C      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
C     &        texpf,texps
CC
CC TIMING - FOR CALCULATING AOP
CC
C      double precision time(ixzmax)
C      common /timing/time,it
C      common /calc/tend,icalc
CC
CC TOLERANCES FOR CHECKING CONVERGENCE OF RESDIUE STRUCTURE
CC RTOLMX = MAXIMUM AND RTOLME = MAX MEAN FRACTION DIFFERENCE
CC ACROSS THE GRID BETWEEN VERTICAL AND HORIZONTAL SWEEPS
CC 0.0001 IS JUST A GUESS ????
CC
C      parameter (rtolmx=0.00001, rtolme=0.00001)
CC
CC MAXIMUM NUMBER OF ITERATIONS FOR RESIDUE CONVERGENCE
CC
C      parameter (maxit=10000)
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC END OF SET UP
CC START OF CALCULATIONS
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC IF STEADY-STATE CALCULATIONTHEN NO NEED TO CALCULATE
CC THE RESIDUE 
C      if(icalc.eq.0) then
C      goto 104
C      endif
CC
CC ARRAY COUNTERS 
CC REO IS THE RESIDUE VALUES AT THE BEGINNING
CC RE1,RE2 ARE USED IN CONVERGENCE LOOP
CC OUTPUT RESIDUE IS SET TO R1 AT THE END
CC
C      re2=6-re0-re1 
CC
CC SET RE0 RESIDUE PROFILE AS STARTING RESIDUE PROFILE
CC TESTED: ARRAYS SET UP FINE
C      open(31,file='res_test_starting')
C      do ix=1,nx
C        do iz=1,nz
C          res(re0,ix,iz)=rest(ix,iz)
C          res(re1,ix,iz)=rest(ix,iz)
C          res(re2,ix,iz)=rest(ix,iz)
C          write(31,*)ix,iz,rest(ix,iz),res(re0,ix,iz),res(re1,ix,iz),
C     &               res(re2,ix,iz)
C        enddo
C      enddo
C      close(31)
CC 
CC COUNTERS FOR VERTICAL AND HORIZONTAL PASS
CC
C      ivpass=0
C      ihpass=0
CC Q:DO THESE NEED TO SEPERATE NAMES TO THE TEMPERATURE ONES?
CC
CC
CC LOOP HERE FOR ITERATION WITHIN THIS PASS TO ALLOW CONVERGENCE
CC
C 101  continue
CC
CC **MAKE UPWRADS VERTICAL PASSES**
CC
C      ivpass=ivpass+1
C      do ix=1,nx
CC       do ix=1,1
CC
CC BASAL BOUNDARY CONDITION FIXED AT RESIDUE = 1, THE PROGRAMME 
CC DOES NOT ALLOW ANY MELTING TO OCCUR AT THE BASAL NODES
CC
C      b(1)=1.
C      c(1)=0.
C      r(1)=res(re0,ix,nz)
CC      write(6,*)'BOTTOM:',r(1)
CC
CC VERTICAL LOOP TO LOAD ARRAYS FOR TRIDIAGONAL MATRIX ALGORITHM.
CC COUNTER IZUP COUNTS UP FROM BOTTOM TO TOP OF CALCULATION BOX
CC THIS IS THE DIRECTION USED IN TRIDAG SINCE IT IS IN THE DIRECTION
CC OF ADVECTION OF THE RESIDUE
CC COUNTER IZDOWN COUNTS UP FROM TOP TO BOTTOM OF CALCULATION BOX
CC THIS IS THE DIRECTION USED FOR THE RESIDUE AND VELOCITY ARRAY
CC
C      do izdown=2,nz-1
C        izup=nz-izdown+1
CC
CC CALCULATE NEIGHBOUR COEFFICIENTS
CC
C      call FVRS2D(anbe,anbw,anbn,anbs,ix,nx,izdown,nz,re1)
CC
CC NEIGHBOUR COEFFICIENTS FOR HORIZONTALLY ADJACENT BOUNDARY POINTS
CC AXIAL BOUNDARY: RESIDUE FIELD IS ASSUMED SYMMETRICAL AS TEMP
CC FIELD IS SYMMETRICAL
CC OFF-AXIAL BOUNDARY: ASSUME VERY LARGE PECLET NUMBER, SO ANBE=0
CC
C      if(ix.eq.1) anbw=anbe
C      if(ix.eq.nx) anbe=0.
CC
CC LOADING TRIDIAGONAL MATRIX FOR UPWARD SWEEP
CC A_S = TRIDAG A()
CC A_N = TRIDAG C()
CC A_P = TRIDAG P()
CC THESE COEFFICIENTS BECOME NEGATIVE WHEN REARRANGED TO
CC BE INPUT INTO THE NUMERICAL RECIPES TRIDAG ROUTINE
CC
C      a(izup)=-anbs
C      c(izup)=-anbn
CC
CC SOURCE TERM IS ZERO, SET SC AND SP EQUAL TO ZERO
CC
C      sc=0.
C      sp=0.
C      if(sp.gt.0.) then
C        write(6,*)'Source term S_P is +ve'
C      endif
CC IF STATEMENT IN PREPARATION A POSSIBLE CHANGE IN
CC THE SOURCE TERM AS THE PROGRAM IS DEVELOPED FURTHER
CC
C      a0p=rhos*dx*dz/(time(it)-time(it-1))
CC
CC CONSTANT TERM B (EQ. 5.57F) ALSO INCLUDES EAST AND WEST
CC NEIGHBOUR COEFFICIENTS WHICH ARE CONSTANT IN THE UPWARD
CC SWEEP. B IS R() IN TRIDAG, RIGHT HAND SIDE OF THE ARRAY
CC 
C      if(ix.gt.1) then
C         ranbw=res(re1,ix-1,izdown)*anbw
C      else
C         ranbw=res(re1,ix+1,izdown)*anbe
C      endif
C      if(ix.lt.nx) then
C         ranbe=res(re1,ix+1,izdown)*anbe
C      else
CC
CC 0 OR 1??
C         ranbe=0.
C      endif
CC B TERM THE TERMS RANBE AND RANBW ARE REARRANGED TO BE INCLUDED
CC IN THE B TERM FROM THE MAIN DISCRETIZATION EQUATION AS THEY 
CC REMAIN CONSTANT IN THE UPWARD SWEEP
C      r(izup)=sc*dx*dz+a0p*res(re0,ix,izdown)+ranbw+ranbe
CC
CC CENTRAL POINT COEFFICIENT A_P INTO ARRAY B()
CC
C      b(izup)=-a(izup)-c(izup)+anbe+anbw+a0p-sp*dxdz
CC
CC END OF VERTICAL LOADING LOOP
CC
C      enddo
CC
CC TOP BOUNDARY CONDITION
CC NO MELTING SHOULD OCCUR IN THE CRUST SO HOPEFULLY IF 
CC FIX BOUNDARY AT 1
CC
C      a(nz)=0.
C      b(nz)=1.
C      r(nz)=res(re0,ix,1)
C      write(6,*)'TOP:',r(nz)
C
C FIND VERTICAL RESIDUE PROFILE USING TRIDIAGONAL MATRIX
C ALGORITHM
C
C      call TRIDAG(a,b,c,r,u,nz)
C MAY NEED SEPERATE TRIDAG ROUTINE
C
C LOAD RESULTS INTO NEW RESIDUE ARRAY
C
C      do izdown=1,nz
C         izup=nz-izdown+1
C           res(re2,ix,izdown)=u(izup)
C           if(ix.eq.1) then
C           test=z(izdown)/z(nz)*res(re0,ix,nz)
C           write(6,*)x(ix)/1000.,z(izdown)/1000.,
C     &      res(re0,ix,izdown),res(re2,ix,izdown),test
C           endif
C      enddo
C
C END OF VERTICAL PASSES
C
C      enddo
C 
C CHECK FOR CONVERGENCE IN THE RESIDUE FIELD
C MAXIMUM AND RMS DIFFERENCES IN RESIDUE BETWEEN
C TWO PASSES
C
C      call RCNVGE(res,rdfmxv,rdfmev,nx,nz,re1,re2)
C       write(6,*)'Max and Mean Difference RESV:',rdfmxv,rdfmev,
C     &        'V AND H LOOPS',ihpass,ivpass
C
C MAKE FINAL RESIDUE FIELD AFTER VERTICAL PASSES
C THE STARTING RESIDUE FIELD FOR HORIZONTAL PASSES
C SHOULD BE ABLE TO USE OTHER SWAP.
C      write(6,*)'Before Swap: re1:',re1,' re2:',re2
C      call SWAP2(re1,re2)
CC      write(6,*)'After Swap: re1:',re1,' re2:',re2
CC
CC **MAKE HORIZONTAL PASSES IN DIRECTION OF PLATE SPREADING**
CC
C 102  continue
C      ihpass=ihpass+1
C      do izdown=2,nz-1
CC       do izdown=nz-1,nz-1
C         izup=nz-izdown+1
CC
CC HORIZONTAL LOOP TO LOAD ARRAYS FOR TRIDIAGONAL MATRIX
CC
C      do ix=1,nx
CC
CC CALCULATE NEIGHBOUR COEFFICIENTS
CC
C      call FVRS2D(anbe,anbw,anbn,anbs,ix,nx,izdown,nz,re1)
CC
CC THE OFF AXIS BOUNDARY - PERHAPS THIS SHOULD BE MADE EQUAL
CC TO THE END NODES IF ALL VELOCITIES ARE HORIZONTAL, OR
CC PERHAPS IT DOESN'T MATTER BECAUSE EVERYTHING IS FLOWING OUT
CC
CC AXIAL BOUNDARY: SYMMETRICAL AS TEMPERATURE FIELD IS SYMMETRICAL
CC
C      if(ix.eq.1) anbe=anbe*2
CC
CC LOADING TRIDIAGONAL MATRIX
CC WESTERN COEFFICIENT - A_W - A()
CC EASTERN COEFFICIENT - A_E - C()
CC COEFFICIENTS BECOME NEGATIVE BECAUSE OF NUMERICAL RECIPES
CC
C      a(ix)=-anbw
C      c(ix)=-anbe
CC
CC SOURCE TERM IS ZERO, SET SC AND SP EQUAL TO ZERO
CC
C      sc=0.
C      sp=0.
C      if(sp.gt.0.) then
C        write(6,*)'Source term S_P is +ve'
C      endif
CC IF STATEMENT IN PREPARATION A POSSIBLE CHANGE IN
CC THE SOURCE TERM AS THE PROGRAM IS DEVELOPED FURTHER
CC
C      a0p=rhos*dx*dz/(time(it)-time(it-1))
CC
CC CONSTANT TERM B IN PATANKHAR - R() IN TRIDAG
CC INCLUDES N AND S COEFFICIENTS
CC
C      r(ix)=sc*dx*dz+a0p*res(re0,ix,izdown)
C     &        +anbn*res(re1,ix,izdown-1)+anbs*res(re1,ix,izdown+1)
CC
CC CENTRAL POINT COEFFICIENT - A_P - B()
CC
C      b(ix)=-a(ix)-c(ix)+anbn+anbs+a0p-sp*dx*dz
CC
CC END OF HORIZONTAL LOADING LOOP
CC
C      enddo
CC 
CC FIND NEW HORIZONTAL RESIDUE PROFILE
CC
CC       write(6,*)'Calling Tridag'
C      call TRIDAG(a,b,c,r,u,nx)
CC       write (6,*)'Return from Tridag'
CC
CC LOAD RESULTS INTO NEW ARRAY
CC
C      do ix=1,nx
C         res(re2,ix,izdown)=u(ix)
C      enddo
CC
CC END OF HORIZONTAL PASSES
CC
C      enddo
CC
CC CHECK FOR CONVERGENCE IN THE RESIDUE FIELD
CC MAXIMUM DIFFERENCES AND RMS DIFFERENCES IN RESIDUE
CC BETWEEN TWO FIELDS
CC
C      call RCNVGE(res,rdfmxh,rdfmeh,nx,nz,re1,re2)
CC       write(6,*)'Max and Mean Difference RESH:',rdfmxh,rdfmea,
CC     &        'V AND H LOOPS',ihpass,ivpass
CC
CC MAKE FINAL RESIDUE FIELD AFETR HORIZONTAL PASSES
CC INTO THE STARTING FIELD FOR VERTICAL PASSES
CC
C      call SWAP2(re1,re2)
C
CC HAS RESIDUE FIELD CONVERGED? 
C IF SO CONTINUES BUT REPORT IF NOT WITHIN SPECIFIED TOLERANCE
C
C 103  continue
C      rdfmx=abs(rdfmxv-rdfmxh)
C      rdfme=abs(rdfmev-rdfmeh)
C      if(rdfmx.lt.1.e-10 .and. rdfme.lt.1.e-10) then
C         write(6,*)'Residue calculation has converged'
C         write(6,*)'***Outside specified tolerance***'
CC WHY IS THIS OUTSIDE THE SPECIFIED TOLERANCE?
C         if(rdfmxh.gt.rtolmx) then
C            write(6,*)'Max resdiue change:',rdfmxh,
C     &                '; tolerated, ',rtolmx
C         endif
C         if(rdfmeh.gt.rtolme) then
C            write(6,*)'Mean resdiue change:',rdfmeh,
C     &                '; tolerated, ',rtolme
C         endif
C         write(6,*)'in',ihpass,' iterations'
C         goto 104
C      endif
CC
CC MAXIMUM NUMBER OF ITERATIONS REACHED WITHOUT CONVERGENCE
CC
C      if (ihpass.gt.maxit) then
C        write(6,*)'Residue calculation in RESFVI:'
C        write(6,*)'Max. iterations exceeded without convergence'
C        write(6,*)'(',ihpass,' iterations )'
C        write(6,*)'Max. residue difference V & H:',rdfmxv,rdfmxh
C        write(6,*)'Max. difference tolerated:',rtolmx
C        write(6,*)'Mean. residue difference V & H:',rdfmev,rdfmeh
C        write(6,*)'Max. difference tolerated:',rtolme
C        write(6,*)'STOPPING'
C        stop
C      endif
CC
CC HAS RESIDUE FIELD CONVERGED WITHIN SPECIFIED TOLERANCES
CC IF NOT LOOP BACK FOR ANOTHER VERTICAL AND HORIZONTAL PASS
CC
C      if(rdfmxv.gt.rtolmx .or. rdfmxh.gt.rtolmx .or.
C     &   rdfmev.gt.rtolme .or. rdfmeh.gt.rtolme) goto 101
CC
CC RESIDUE FIELD CONVERGED
CC  
C 104  continue
CC
C WRITE OUT RESIDUE FIELD
CC
CC
CC END SUBROUTINE
CC
C      return
C      end
CC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C SUBROUTINE TO CALCULATE NEIGHBOUR COEFFICIENTS FOR RESIDUE
C
C      subroutine FVRS2D(anbe,anbw,anbn,anbs,ix,nx,iz,nz,re1)
CC
CC DIMENSIONS OF THE FIELD NX AND NZ ARE USED TO CHECK WHETHER THE
CC POINT IS ON A BOUNDARY. COEFFICIENTS ARE NOT CALCULATED OUTSIDE THE
CC BOUNDARY
CC OUTPUT: NEIGHBOUR COEFFICIENTS ANBE, ANBW, ANBN, ANBS.
CC
C      implicit double precision(a-h,o-z)
CC
CC MAXIMUM SIZE OF GRID
CC
C      parameter (ixzmax=2001)
CC
CC RESIDUE GRID
CC
C      double precision rest(ixzmax,ixzmax),res(3,ixzmax,ixzmax)
C      common /resdue/rest,res,resc
CC
CC VELOCITY GRID STAGGERED WITH RESPECT TO TEMPERATURE GRID
CC NOTE: VELOCITY IS STAGGERED TO THE EAST AND SOUTH SIDE OF 
CC EACH NODE 
CC
C      double precision vxs(ixzmax,ixzmax),vzs(ixzmax,ixzmax)
C      common /velstg/vxs,vzs
CC
CC GRID
CC
C      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
C      common /grid/xmax,zmax,dx,dz
CC      common /grido/pres,x,z,nx,nz
CC REMOVED ABOVE COMMON TO ALLOW NX TO TRANSFER IN FUNCTION
CC
CC PARAMETERS
CC
C      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
C     &       texpf,texps
CC
CC NEIGHBOUR COEFFICIENTS CACULATED USING ADVECTION ONLY
CC GAMMA THE DIFFUSION PARAMETER IS ZERO, SO THE CONDUCTANCES
CC D BECOME ZERO. THE PECLET NUMBER IS VERY LARGE (DOMINATED BY
CC CONVECTION BUT THE POWER LAW AND PECLET NUMBER IS IRRELEVANT
CC BECAUSE THE EQUATION FOR THE NEIGHBOUR COEFFICIENT MULTIPLIES
C THE A(P) FUNCTION BY D WHICH IS ZERO.
CC
CC F. IS THE MASS FLOW RATE, R IS THE DENSITY
C
C EASTERN NEIGHBOUR COEFFICIENT
C
C      if (ix.lt.nx) then
C         re=rhos
C         fe=re*vxs(ix,iz)*dz
C         anbe=max(-fe,0.)
C      else
C         anbe=0.
C      endif
CC
CC WESTERN NEIGHBOUR COEFFICIENT
C
C      if(ix.gt.1) then
C        rw=rhos
C        fw=rw*vxs(ix-1,iz)*dz
C        anbw=max(fw,0.)
C      else
C        anbw=0.
C      endif
CC
CC NORTHERN NEIGHBOUR COEFFICIENT
CC
C      if(iz.gt.1) then
C        rn=rhos
C        fn=rn*vzs(ix,iz-1)*dx
C        anbn=max(-fn,0.)
C      else
C        anbn=0.
C      endif
CC
CC SOUTHERN NEIGHBOUR COEFFICIENT
CC
C      if(iz.lt.nz) then
CC        rs=rhos
C        fs=rs*vzs(ix,iz)*dx
C        anbs=max(fs,0.)
C      else
C        anbs=0.
C      endif
CC
CC END OF SUBROUTINE
CC
C      return
C      end
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
C      subroutine RCNVGE(res,rdfmax,rdfmea,nx,nz,re0,re1)
CC            call RCNVGE(res,rdfmxv,rdfmev,nx,nz,re1,re2)
CC            call RCNVGE(res,rdfmxh,rdfmeh,nx,nz,re1,re2)
CC
CC TEST FOR CONVERGENCE BETWEEN STARTING AND FINISHING
CC RESIDUE FIELDS. RETURNS:
CC MAXIMUM RESIDUE DIFF, RDFMAX
CC MEAN RESIDUE DIFF, RDFMEA
CC
CC REAL DECLARATION
CC
C      implicit double precision(a-h,o-z)
CC
CC MAXIMUM SIZE OF GRID
CC
C      parameter (ixzmax=2001)
CC 
CC RESIDUE GRID
CC
C      double precision res(3,ixzmax,ixzmax)
CC
CC COMPARE RESIDUE GRIDS
CC
C      rdfmax=0.
C      rdfmea=0.
C      do ix=1,nx
C       do iz=1,nz
C         dif=abs(res(re1,ix,iz)-res(re0,ix,iz))
C         if(dif.gt.rdfmax) rdfmax=dif
C         rdfmea=rdfmea+dif
C       enddo
C      enddo
C      rdfmea=rdfmea/float(nx*nz)
CC
CC END SUBROUTINE
CC
C      return
C      end
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
CC SUBROUTINE TO SWAP RE1 AND RE2 OVER
C      subroutine SWAP2(re1,re2)
C      implicit double precision(a-h,o-z)
C      rtmp=re1
C      re1=re2
C      re2=rtmp
C      return
C      end
C
C-----------------------------------------------------------
C
      function CINTRP(tp0,z0,iel,itym)
C
C BILINEAR INTERPOLATION
C RETURNS COMPOSITION AT GIVEN POTENTIAL TEMPERATURE AND DEPTH
C STORED IN LOOK-UP TABLE C( EL, TP, Z)
C
      implicit double precision(a-h,o-z)
      parameter (iemax=20)
      parameter (ixzmax=2001)
      character*2 elabel(iemax)
      character*100 tzfile(iemax)
      common /comp/tzfile,elabel,gdepth,icomp,isourc,nel
      double precision c(iemax,ixzmax,ixzmax),ctp(ixzmax),
     &                 cz(ixzmax)
      common /cinst/c,cz,ctp,ncz,nctp
C
C JLOTP AND JLOZ LOCATE THE GRID SQUARE FOR INTERPOLATION
C AND ARE SAVED TO SAVE TIME WHEN SEARCHING NEXT TIME
C
      save jlotp,jloz
      if (itym.le.1) then
        jlotp=1
        jloz=1
      endif
C
C FIND GRID SQUARE FOR INTERPOLATION 
C
      call HUNT(ctp,nctp,tp0,jlotp)
      call HUNT(cz,ncz,z0,jloz)
c      write(6,*)jlotp,jloz
C
C BILINEAR INTERPOLATION COEFFICIENTS (NUM. REC. 3.6.4)
C
      ftp=(tp0-ctp(jlotp))/(ctp(jlotp+1)-ctp(jlotp))
      fz=(z0-cz(jloz))/(cz(jloz+1)-cz(jloz))
      ftpc=1.0-ftp
      fzc=1.0-fz
c      write(6,*)ftp,fz,ftpc,fzc
C
C BILINEAR INTERPOLATION (NUM. REC. 3.6.3 & 3.6.5)
C
      cinst=ftpc*fzc*c(iel,jlotp,jloz)
     &    +ftp*fzc*c(iel,jlotp+1,jloz)
     &    +ftp*fz*c(iel,jlotp+1,jloz+1)
     &    +ftpc*fz*c(iel,jlotp,jloz+1)
c      write(6,*)c(iel,jlotp,jloz),
c     &    c(iel,jlotp+1,jloz),
c     &    c(iel,jlotp+1,jloz+1),
c     &    c(iel,jlotp,jloz+1)
C
C ENSURE INSTANTANEOUS COMPOSITION IS NOT NEGATIVE
C
      if (cinst.gt.0.0) then
        CINTRP=cinst
      else 
        CINTRP=0.0
      endif
C
      return
      end
C
C-----------------------------------------------------------
C
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      DOUBLE PRECISION x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
