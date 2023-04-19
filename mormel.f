C RACHEL WALTERS JUNE 27TH 2006
C STEVE JONES MAY 2005
C PROGRAM TO CALCULATE THICKNESS OF OCEANIC CRUST
C BASED ON BOWN AND WHITE, EARTH AND PLANETARY SCIENCE LETTERS, 1994.
C AND JOURNAL OF GEOPHYSICAL RESEARCH, 1995
C
C LENGTH IN M
C MASS IN KG
C TIME IN YR
C PRESSURE IN GPA
C TEMPERATURE IN DEGREES CELCIUS
C
C THE FOLLOWING ROUTINES ARE HERE
C	SETP
C	READP
C READC
C	SETGRD
C	TYMING
C	TILOOP
C 	SETVEL
C	SETRES
C	SETTEM
C	POPEN
C	FINISH
C	INTRES
C
C
      program WJ
      implicit double precision(a-h,o-z)
C
C MAXIMUM ARRAY SIZES
C
      parameter (ixzmax=2001)
      parameter (itimax=2001)
      parameter (iamax=10)
      parameter (iemax=20)
C
C TEMPERATURE INPUT AND OUTPUT FILES
C
      character*50 pfileh,pfilev,tfilei,tfileo
      common /tem/ptem,tfilei,tfileo,item
C
      common /countr/itym,ntime,iendpr
C
C AVERAGE MELT PRODUCTION RATE & CRUSTAL THICKNESS
C
      common /intgra/tarea,tcrust
C 
C PLATE SPREADING HISTORY AND CORNER FLOW GEOMETRY
C 
      double precision ufull(iamax)
      double precision utime(iamax)
      common /vel/ufull,utime,alpha,viss,crust,nu,ialpha
C
      common /struct/dlith,djump,iriftr,ijump
C
C TIMING
C
      double precision tyme(itimax)
      double precision ufulla(itimax)
      double precision ufullt(itimax)
      double precision age(itimax)
      double precision rtyme(itimax)
      common /uftym/ufullt,ufulla,age,tyme,rtyme
C
C
C
      common /temcor/temc
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C
C GRID COORDINATES
C
      double precision pres(ixzmax)
      double precision x(ixzmax)
      double precision z(ixzmax)
      common /grido/pres,x,z,nx,nz
C
C 
C
      double precision temc(ixzmax,ixzmax)
      double precision frac(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
C
C POTENTIAL TEMPERATURE
C  TP, POTENTIAL TEMPERATURE
C  TPS, STARTING POT TEMP
C
      double precision tp(3,ixzmax,ixzmax)
      double precision tps(ixzmax,ixzmax)
      common /pottem/tp,tps
C
C RESIDUE
C
      double precision rest(ixzmax,ixzmax)
      double precision res(3,ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C COMPOSITION
C
      character*2 elabel(iemax)
      character*100 tzfile(iemax)
      double precision fci(6,ixzmax,ixzmax)
      double precision fca(6,ixzmax,ixzmax)
      double precision tcomp(6)
      double precision carea(6)
      double precision sci(6,ixzmax,ixzmax)
      double precision c(iemax,ixzmax,ixzmax)
      double precision ctp(ixzmax)
      double precision cz(ixzmax)
      double precision cpda(iemax)
      common /comp/tzfile,elabel,gdepth,icomp,isourc,nel
      common /cinst/c,cz,ctp,ncz,nctp
      common /cav/cpda
      common /fakecp/fci,fca,sci,tcomp,carea
C
C COMMON FOR MELT RATE
C
      double precision dpda(ixzmax,ixzmax),dtrda(ixzmax,ixzmax),
     &     dfda(ixzmax,ixzmax)
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
C
C COMMUNICATION OF DEPTH DEPENDENT X
C FI = INSTANTANEOUS MELT FRACTION
C FT = TOTAL MELT FRACTION (WHERE FI>0)
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
C
C OPEN FILES FOR QUANTITIES TO BE WRITTEN OUT AT EACH TIME STEP
C
      open(52,file='WJ_crust_out')
      open(53,file='WJ_comp_out')
      open(55,file='WJ_meanf_out')
C
C PARAMETERS WHICH REMAIN CONSTANT
C
      write(6,*)
      call SETP
C
C READ IN PARAMETERS SPECIFIED BY USER
C
      call READP
C
C READ IN INSTANTANEOUS COMPOSITION DATA
C
      if (icomp.eq.3) call READC
C
C SET UP GRID IN X AND Z
C
      call SETGRD
C
C SET UP RESIDUE GRID
C
      if (item.ne.1) call SETRES
C
C SET UP POTENTIAL TEMPERATURE GRID
C
      call SETTP
C
C HEADER FOR OUTPUT FILE
C
 130  format('Col ',i2,2x,a40)
      write(52,130) 1,'time step number'
      write(52,130) 2,'time (Myr)'
      write(52,130) 3,'spreading rate (km/Myr)'
      write(52,130) 4,'melt thickness (km)'
      write(52,130) 5,'mean melt production rate (km^2/yr?)'
      write(52,130) 6,'depth to maximum Tp (km)'
      do iel=1,nel
        write(52,130) 6+iel, elabel(iel)//' concentration ratio'      
      enddo
 131  format(23(a2,1x))
      write(52,131)'EL','2','3','4','5','6',(elabel(iel), iel=1,nel )
 132  format('DATA starts here')
      write(52,132)
C
C LOOP FOR EACH TIME STEP
C
 6    continue
C
C CALL THE TIMING SUBROUTINE TO CHANGE SWITCHES AND CONTROL PROGRAM
C
      call TYMING
C
C LOOP FOR MAKING VELOCITIES COMPATABLE WITH CRUSTAL THICKNESS
C
 5    continue
C
C CALCULATE VELOCITES AT EACH NODE FROM CORNER FLOW MODEL
C
      call SETVEL
C
C SET UP STARTING THERMAL STRUCTURE
C EITHER READ IN OR CALCULATED FROM HALF SPACE COOLING MODEL
C
      call SETTEM
C
C CHANGE BASAL TEMPERATURE BOUNDARY CONDITION
C TEST TIME-DEPENDENT LOWER BOUNDARY LAYER
C
      call TEMBAS(temmax,tpmax) 
C
C SET UP 1D OUTPUT PROFILES
C
cc      call POPEN
C
C CALCULATE (VELOCITY) TEMPERATURE AND MELTING
C
      call CALVTM(n1,re1,ntp0,ntp1)
C
C DETERMINE WHETHER CRUST ITERATION IS NECESSARY
C
      errcst=0.1*1.e3
C      dcrust=tcrust-crust
C      if (abs(dcrust).gt.errcst) then
C        write(6,*)'Crust iteration, dcrust:',dcrust*1.e-3
         crust=tcrust
C        goto 5
C      else
C        if (itym.eq.1) then
C          if (itype.eq.1 .or. itype.eq.3) then
C            write(52,*)0,0.,0.,0.
C          endif
C       endif
C
C WRITE TO SCREEN 
C
      call AZMAX(tp,1,ntp1,aval,zval)
 101  format('Time: ',f7.2,1x,
     &       'Crust:',f6.3,1x,
     &       'Max basal Tp: ',f6.1,1x,
     &       'Max Tp:',f8.2,1x,
     &       'at depth: ',f6.1)
c 102  format('Time: ',f7.2,5x,'max Tp:',f8.2,5x,'at depth: ',f6.1)
      write(6,101) tyme(itym)*1.e-6,tcrust*1.e-3,tpmax,
     &             aval,zval*1.e-3
c      write(6,102) tyme(itym)*1.e-6,aval,zval*1.e-3
c      write(6,*)'Time: ',tyme(itym)*1.e-6,
c     &          '  Crust: ',tcrust*1.e-3,
c     &          '  Max basal temp.: ',temmax
C
C FILE: WJ_CRUST_OUT
C      Col 1:  step number
C      Col 2:  time (Myr)
C      Col 3:  spreading rate (km/Myr)
C      Col 4:  melt thickness (km)
C      Col 5:  mean melt production rate (m^2/yr?)
C      Col 6:  depth to maximum Tp (km)
C      Col 7:  concentration of element 1
C        ...
C      Col 7+NEL:  concentration of element NEL
C
 103  format(i4,1x,f7.2,1x,f6.1,1x,f6.2,1x,f6.3,1x,f6.1,20(1x,f7.3))
      write(52,103)itym,tyme(itym)*1.e-6,ufullt(itym)*1.e3,
     &           tcrust*1.e-3,tarea,zval*1.0e-3,
     &           (cpda(iel), iel=1,nel )
C
C FILE: WJ_TIMED_OUT
C
C      do ix=1,nx
C         do iz=1,nz
C         write(51,*)itym,tyme(itym)*1.e-6,
C      &              x(ix)*1.e-3,z(iz)*1.e-3,frac(ix,iz),
C      &              temc(ix,iz),rest(ix,iz)
cc         write(51,*)itym,tyme(itym)*1.e-6,
cc     & x(ix)*1.e-3,z(iz)*1.e-3,frac(ix,iz),temc(ix,iz),nfrac(ix,iz)
C         enddo
C      enddo
C      endif
C
C WRITE OUT FINAL DATA AND CLOSE FILES
C
cc      call FINISH(n1)
C
C CALL SUBROUTINE TO CHOOSE WHICH COMPOSITION TO GO WITH
C AND THEN CALCULATE COMPOSITION
C
      call COMPSW(n1,re1,ntp1)
      if (icomp.eq.3) then
 112  format('Time: ',f7.2,5x,a2,f6.2)
        do iel=1,nel
          write(6,112) tyme(itym)*1.e-6,elabel(iel),cpda(iel)
        enddo
c        stop
      endif
C
C INTEGRATION TO GET MELT FRACTION VERSUS DEPTH
C FILE 53: WJ_COMPC_OUT
C FILE 55: WRITE OUT THE F DEPENDENT ON Z ARRAY FOR SFI PROPOSAL
C A SIMPLE INTEGRATION FROM SUBROUTINE FZ
C
c      call FZ(n1,re1)
c      do ic=1,6
c        write(53,*)itym,tyme(itym),ic,tcomp(ic)
cC        write(6,*)'Element:',ic,'Integral:',tcomp(ic)
c      enddo
c      do iz=1,nz
c        write(55,*)itym,tyme(itym),z(iz)*1.e-3,
c     &            fidarr(iz),ftdarr(iz)
c      enddo
C
C NEED TO CALL AN END OF TIME LOOP SUBROUTINE
C 
      call TILOOP(n1,re1)
C
C CLAUSE TO FINISH PROGRAM OR LOOP THE TIME
C
      if (iendpr.eq.0) then
        goto 6
      elseif (iendpr.eq.1) then
C
C STATE WHERE DATA IS
c      write(52,*)itym,tyme(itym)*1.e-6,ufull*1.e3,tcrust*1.e-3
C
      write(6,*) 'Time series is in file "WJ_crust_out"'
      write(6,*) '  Col 1:  step number'
      write(6,*) '  Col 2:  time (Myr)'
      write(6,*) '  Col 3:  spreading rate (km/Myr)'
      write(6,*) '  Col 4:  melt thickness (km)'
      write(6,*) '  Col 5:  mean melt production rate (km^2/yr?)'
      write(6,*) '  Col 6:  depth to maximum Tp (km)'
      write(6,*) '  Col 7:  concentration of element 1'      
      write(6,*) '      ...'      
      write(6,*) '  Col 7+nel:  concentration of element nel'      
C 
C WRITE OUT CORRECTED TEMPERATURE STRUCTURE TEMC(IX,IZ)
C
        open(19,file='WJ_tem_res_out')
C5       format (4(f10.2,1x))
        do iz=1,nz
          do ix=1,nx
            write(19,*)x(ix)*1.e-3,z(iz)*1e-3,
     &                 temc(ix,iz),rest(ix,iz),
     &                 tp(ntp1,ix,iz),dfda(ix,iz)
          enddo
        enddo
        close(19)
        write(6,*) 'Written final spatial data to "WJ_tem_res_out"'
        write(6,*) '  Col 1:  x (km)'
        write(6,*) '  Col 2:  z (km)'
        write(6,*) '  Col 3:  temperature (deg C)'
        write(6,*) '  Col 4:  unmelted fraction'
        write(6,*) '  Col 5:  potential temperature'
        write(6,*) '  Col 6:  overall melting rate (km^2/yr?)'
C
c        open(54,file='WJ_melt_out')
c        do ix=1,nx
c          do iz=1,nz
c            write(54,*)x(ix)*1.e-3,z(iz)*1.e-3,dfda(ix,iz)
c          enddo
c        enddo
c        write(6,*)
c     & 'Written final melting rates to "WJ_melt_out"'
c        close(54)

c        close(51)
        close(52)
        close(53)
        close(55)
        write(6,*)'Program Finished!'
      endif
C
 7    end
C
C----------------------------------------------------------
C
      subroutine SETP
      implicit double precision(a-h,o-z)
C
C VALUES FOR CONSTANTS USED THROUGHOUT THE PROGRAM
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C
C WRITE
C
      write(6,*)
      write(6,*)'SETP: Fixed parameters (recompile to change)'
C
C PI
C
      pi=3.141592654
C
C RADIANS PER DEGREE
C
      rad=pi/180.
C
C SECONDS PER YEAR
C
      sy=60.*60.*24.*365.25
      write(6,*)'Seconds per year:',sy
C
C ACCELERATION DUE TO GRAVITY (M/YR^2)
C
      g=9.81
      write(6,*)'Gravity (m/s^2):',g
      g=g*sy*sy
      write(6,*)'Gravity (m/yr^2):',g
C
C THERMAL DIFFUSIVITY (M^2/YR)
C BOWN & WHITE USE 0.804E-6
C SPIEGELMAN & MCKENZIE USE 1.E-6
C
      tdiff=0.804e-6
      write(6,*)'Thermal diffusion (m^2/s):',tdiff
      tdiff=tdiff*sy
      write(6,*)'Thermal diffusion (m^2/yr):',tdiff
C
C SPECIFIC HEAT AT CONSTANT PRESSURE (J/KG/K OR M^2/S^2/K)
C
C      cp=1.2e3
      cp=1.0e3
      write(6,*)'Specific heat at constant P (J/kg/K):',cp
      cp=cp*sy*sy
      write(6,*)'Specific heat at constant P (m^2/yr^2/K):',cp
C
C THERMAL EXPANSION COEFFICIENT OF SOLID MANTLE (/K)
C FROM BOWN & WHITE 1995
C
      texps=3.2e-5
      write(6,*)'Thermal expansion of matrix (/K):',texps
C
C THERMAL EXPANSION COEFFICIENT OF MAGMA (/K)
C FROM BOWN & WHITE 1995
C
      texpf=6.8e-5
      write(6,*)'Thermal expansion of melt (/K):',texpf
C
C DENSITY OF MANTLE MATRIX
C
      rhos=3.3e3
      write(6,*)'Density of matrix (kg/m^3):',rhos
C
C DENSITY OF MELT
C
      rhom=2.8e3
      write(6,*)'Density of melt (kg/m^3):',rhom
C
C DIFFERENCE BETWEEN MATRIX AND MELT DENSITIES (KG/M^3)
C
      drho=rhos-rhom
      write(6,*)'Melt-matrix density difference (kg/m^3):',drho
C
C CONSTANT THERMAL CONDUCTIVITY SPECIFIC HEAT ALREADY GIVEN
C CALCULATED FROM DIFFUSION, DENSITY AND 
C SET TO BE COMPATIBLE WITH BOWN & WHITE (1994)
C
      tcon=tdiff*cp*rhos
      write(6,*)'Thermal conductivity (W/m/K):',tcon/sy/sy/sy
      write(6,*)'Thermal conductivity (kg.m/yr^3/K):',tcon
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine READP
C
C READ IN PARAMETERS REQUIRED FOR MODELLING
C GRID: PARAMETERS TO SET UP GRID (X,Z) TO SOLVE PDE
C VEL: VELOCITY FIELD
C TEM: TEMPERATURE BOUNDARY CONDITIONS
C MELT: MELTING CALCULATION
C OUT: TIME INTERVAL FOR OUTPUT
C
      implicit double precision(a-h,o-z)
      parameter (iamax=10)
      parameter (iemax=20)
      parameter (ixzmax=2001)
      logical found
      character*50 pfileh,pfilev,tfilei,tfileo
      character*2 elabel(iemax)
      character*100 tzfile(iemax)
      double precision c(iemax,ixzmax,ixzmax),ctp(ixzmax),
     &                 cz(ixzmax)
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
      common /grid/xmax,zmax,dx,dz
C 
C PLATE SPREADING HISTORY AND CORNER FLOW GEOMETRY
C 
      double precision ufull(iamax)
      double precision utime(iamax)
      common /vel/ufull,utime,alpha,viss,crust,nu,ialpha
C
      common /struct/dlith,djump,iriftr,ijump
      common /tem/ptem,tfilei,tfileo,item
      common /melt/ds,xh2o,cpxm,iparam,imelt,itemp
C      common /calc/itype,timfin,timstp,tend,icalc
      common/calc/tend,icalc
      common /prfout/xprf,wprf,zprf,tprf,pfileh,pfilev
      common /countr/itym,ntime,iendpr
      common /timset/timfin,timstp,itype
c      common /comp/icomp,gdepth,isourc
      common /comp/tzfile,elabel,gdepth,icomp,isourc,nel
      common /cinst/c,cz,ctp,ncz,nctp
			common /minusa/aplumex,aplumey,aplumer,
     &               abslith,absasth,aplvflx
			common /minusb/bpuldur,bpulst,bpultem,ibpul
      namelist /grid/xmax,zmax,dx
      namelist /vel/ufull,utime,alpha,viss,crust,ialpha
      namelist /struct/dlith,djump,iriftr,ijump
      namelist /tem/ptem,tfilei,tfileo,item
      namelist /melt/ds,xh2o,cpxm,iparam,imelt,itemp
      namelist /calc/itype,timfin,timstp,tend,icalc
      namelist /timset/itype,timfin,timstp
      namelist /prfout/pfilev,pfileh,tprf,wprf,xprf,zprf
      namelist /comp/icomp,gdepth,isourc,nel,nctp,ncz,
     &               elabel,tzfile
			namelist /minusa/aplumex,aplumey,aplumer,
     &                 abslith,absasth,aplvflx
			namelist /minusb/ibpul,bpuldur,bpulst,bpultem
C
      inquire (file='wj_parameters',exist=found)
      if (.not.found) then
        write(6,*) 'File "wj_parameters" not found.'
        stop
      endif
      write(6,*)
      write(6,*)'READP: Parameters from input file'
      open(4,file='wj_parameters')
C
      read(4,grid)
      write(6,grid)
      write(52,grid)
C
      read(4,vel)
      write(6,vel)
      write(52,vel)
      do ia=2,iamax
        if (utime(ia).le.utime(ia-1)) then
          nu=ia-1
          goto 10
         endif
      enddo
C
 10   read (4,struct)
      write(6,struct)
      write(52,struct)
C
      read(4,tem)
      write(6,tem)
      write(52,tem)
C
      read(4,melt)
      write(6,melt)
      write(52,melt)
C
      read(4,calc)
      write(6,calc)
      write(52,calc)
C
      read(4,timset)
      write(6,timset)
      write(52,timset)
C
      read(4,comp)
      write(6,comp)
      write(52,comp)
C
      read(4,prfout)
      write(6,prfout)
      write(52,prfout)
C
      read(4,minusa)
      write(6,minusa)
      write(52,minusa)
C
      read(4,minusb)
      write(6,minusb)
      write(52,minusb)
C
      close(4)
C
C CHANGE FROM INPUT TO WORKING UNITS
C DISTANCE KM -> M
C
      crust=crust*1.e3
      dx=dx*1.e3
      do ia=1,nu
        ufull(ia)=ufull(ia)*1.e-3
      enddo
      xmax=xmax*1.e3
      zmax=zmax*1.e3
      wprf=wprf*1.e3
      xprf=xprf*1.e3
      zprf=zprf*1.e3
      gdepth=gdepth*1.e3
      dlith=dlith*1.e3
      djump=djump*1.e3
      aplumex = aplumex*1.0e3 
      aplumey = aplumey*1.0e3
      abslith = abslith*1.0e3
      absasth = absasth*1.0e3
      aplvflx = aplvflx*1.0e9      
C 
C DEGREES -> RADIANS
C
      alpha=alpha*rad
C
C TIME S -> YR
C
      viss=viss*sy
C
C TIME MYR -> YR
C
      tend=tend*1.e6
      tprf=tprf*1.e6
      timfin=timfin*1.e6
      timstp=timstp*1.e6
      bpuldur = bpuldur*1.0e6
      bpulst = bpulst*1.0e6
      do ia=1,nu
        utime(ia)=utime(ia)*1.0e6
      enddo
C
C SET ITYM TO ZERO AS IS START OF PROGRAM
C
      itym=1
C
      return
      end
C
C----------------------------------------------------------
C
C READ IN INSTANTANEOUS MELT COMPOSITIONS
C
      subroutine READC
      implicit double precision(a-h,o-z)
C
C READS INFORMATION FROM THE NEL FILES WHOSE NAMES ARE
C STORED IN TZFILE.  EACH FILE CONTAINS INSTANTANEOUS
C MELT COMPOSITION RATIOS FOR ONE ELEMENT.  THE ELEMENT
C ATOMIC SYMBOLS ARE STORED IN ELABEL.
C EACH FILE MUST CONTAIN ASCII TEXT IN 3 SPACE-SEPARATED COLUMNS: 
C   POTENTIAL TEMPERATURE, DEPTH, COMPOSITION RATIO
C THE POTENTIAL TEMERATURE AND DEPTH VALUES MUST BE THE
C SAME IN ALL OF THE FILES; THIS IS CHECKED FOR AND THE
C PROGRAM STOPS IF THE TP, Z VALUES ARE NOT THE SAME
C THIS INFORMATION IS USED LATER TO ESTIMATE THE AVERAGE
C MELT COMPOSITION GENERATED OVER THE WHOLE MELTING REGION
C AS A FUNCTION OF TIME.  
C 
      parameter (iemax=20)
      parameter (ixzmax=2001)
      logical found
      character*2 elabel(iemax)
      character*100 tzfile(iemax)
      double precision c(iemax,ixzmax,ixzmax),ctp(ixzmax),
     &                 cz(ixzmax)
      common /comp/tzfile,elabel,gdepth,icomp,isourc,nel
      common /cinst/c,cz,ctp,ncz,nctp
C
C LOOP OVER EACH ELEMENT
C
      write(6,*)'READC: Reading in compositional data'
      do iel=1,nel
C
C CHECK THAT FILE EXISTS AND OPEN IT
C
        inquire (file=tzfile(iel),exist=found)
        if (.not.found) then
            write(6,*) 'Melt composition file not found: ',
     &                 tzfile(iel)
            stop
        endif
        open (20,file=tzfile(iel))
        write(6,*)
     & 'Reading in instantaneous compositions from file ',
     &  tzfile(iel)
C
C IF FIRST ELEMENT, COUNT AND STORE TEMPERATURE AND DEPTH
C VALUES AS WE READ IN
C
        if (iel.eq.1) then
          do icz=ncz,1,-1
            do ictp=1,nctp
              read(20,*,end=10) ctp(ictp), cz(icz), 
     &                          c(iel,ictp,icz)
c              write(6,*)ictp,icz,ctp(ictp),cz(icz),
c     &                  c(iel,ictp,icz)
            enddo
          enddo
C
C FINISHED READING IN FIRST FILE
C CHECK THAT THE TP AND Z INCREASE MONOTONICALLY
C
          close (20)
          do icz=2,ncz
            if (cz(icz).le.cz(icz-1)) then
              write(6,*)'Depths do not increase monotonically'
              stop
            endif
          enddo
          do ictp=2,nctp
            if (ctp(ictp).le.cz(ictp-1)) then
              write(6,*)'Temperatures dot no increase monotonically'
              stop
            endif
          enddo
C
C CHECK THAT SUBSEQUENT FILES HAVE THE SAME TP AND Z VALUES
C
        else
          do icz=ncz,1,-1
            do ictp=1,nctp
              read(20,*)tptmp,ztmp,c(iel,ictp,icz)
              tpdif=abs(tptmp-ctp(ictp))
              zdif=abs(ztmp-cz(icz))
c              write(6,*)ictp,icz,ctp(ictp),cz(icz),
c     &                  c(iel,ictp,icz)
              if (tpdif.gt.1.e-4 .or. zdif.gt.1.e-4) then
                write(6,*)'Input file not correctly ordered',ictp,icz
                write(6,*)'  Expected ',ctp(ictp),cz(icz)
                write(6,*)'  Read ',tptmp,ztmp
                stop
              endif
            enddo
          enddo
          close (20)
C
C END OF LOOP OVER ELEMENTS
C
        endif
      enddo
C
      return
C
C END OF FILE REACHED UNEXPECTEDLY
C
 10   close (20)
      write(6,*)"ERROR: End of file reached unexpectedly"
      write(6,*)"       Check values of NCTP and NCZ"
      stop
C
      end
C
C----------------------------------------------------------
C
C SET UP GRID IN X AND Z
C
      subroutine SETGRD
      implicit double precision(a-h,o-z)
C
C MAXIMUM NUMBER OF NODES IN BOTH X AND Z DIRECTIONS
C
      parameter (ixzmax=2001)
C MAXIMUM NUMBER OF GRID CHANGES IN EACH DIRECTION
C     parameter (igmax=100)
      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
C      double precision xgrid(igmax),zgrid(igmax),dxt(igmax),dzt(igmax),
C     & dxe(ixzmax,ixzmax),dxw(ixzmax,ixzmax),dzn(ixzmax,ixzmax),
C     & dzs(ixzmax,ixzmax),nxt(igmax),nzt(igmax)
C
C INPUT PARAMETERS SPECIFING GRID DIMENSIONS AND MESH SIZE
C
      common /grid/xmax,zmax,dx,dz
C     common /grid/xmax,zmax,xgrid,zgrid,dxt,dzt,nxt,nzt
C
C COMMUNITATION OF GRID WITHIN PROGRAM
C
      common /grido/pres,x,z,nx,nz
C     common /grido/pres,x,z,nx,nz,dxe,dxw,dzn,dzs
C
C COMMUNICATION OF FIXED PARAMETERS
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C
C WRITE TO SCREEN
C
      write(6,*)
      write(6,*)'SETGRD: Set up grid geometries'
C
C SET UP VARIABLE GRID MESH
C 1. READ IN GRID SPECIFICATIONS FROM FILE
C      inquire (file='wj_grid_specs',exist=found)
C      if (.not.found) then
C        write(6,*) 'File "wj_grid_specs" not found.'
C	stop
C      endif
C      write(6,*)'READ in grid specifications'
C      open(7,file='wj_grid_specs')
C       igrid=1
C XGRID IS THE X VALUE WHERE THE GRID SPACING CHANGES AND DXT 
C IS THE GRID SPACING BETWEEN THE PREVIOUS XGRID AND THE PRESENT SO 
C IF XGRID(1)=10 THEN DXT(1) COVERS 0 TO 10KM
C 8     read(7,*,end=9)xgrid(igrid),dxt(igrid),zgrid(igrid),dzt(igrid)
C       if(dxt(igrid) .leq. 0. .or. dxz(igrid) .leq. 0.) then
C        write(6,*)'Negative or zero values given for dx and dz'
C        stop
C       endif
C       if(xgrid(igrid) .lt. 0. .or. zgrid(igrid) .lt. 0.) then
C        write(6,*)'Negative values given for xgrid and zgrid'
C        stop
C       endif
C CONVERT VALUES TO THE RIGHT UNITS
C       xgrid(igrid)=xgrid(igrid)*1.e3
C       dxt(igrid)=dxt(igrid)*1.e3
C       zgrid(igrid)=zgrid(igrid)*1.e3
C       dzt(igrid)=dzt(igrid)*1.e3
C IGRID IS THE NUMBER OF THE GRID CHANGE SET UP
C       igrid=igrid+1
C       goto 8
C NGRID IS THE MAXIMUM NUMBER OF CHANGES
C 9     ngrid=igrid-1
C       if(xgrid(ngrid).neq.xmax .or. zgrid(ngrid).neq.zmax)then
C       write(6,*)'Grid specifications not complete, xgrid and zgrid
C    & last entries must have the same value as xmax and zmax'
C       endif        
C      close(7)
C 2.SET UP GRID IN THE X DIRECTION
C      do igrid=1,ngrid
C         if(igrid.eq.1) then
C           nxt(igrid)=int(xgrid(igrid)/dxt(igrid))+1
C           do ix=1,nxt(igrid)
C              x(ix)=float(ix-1)*dxt(igrid)
C           enddo
C         else
C           nxt(igrid)=nxt(igrid-1)+int(xgrid(igrid)/dxt(igrid))
C          do ix=nxt(igrid-1)+1,nxt(igrid)
C             x(ix)=x(ix-1)+float(ix-nxt(igrid-1))*dxt(igrid)
C          enddo
C         endif
C      enddo
C SET FINAL TOTAL OF X NODES
C      nx=nxt(ngrid)
C
C 3.SET UP GRID IN THE Z DIRECTION NEGATIVELY INCREASING DOWNWARDS
C      do igrid=1,ngrid
C         if(igrid.eq.1) then
C           nzt(igrid)=int(zgrid(igrid)/dzt(igrid))+1
C           do iz=1,nzt(igrid)
C              z(iz)=-float(iz-1)*dzt(igrid)
C
C WRITE OUT THE ARRAYS TO CHECK
C
C
C
C
C
C
C GRID MESH SIZE ASSUMED SAME IN X AND Z
C
      dz=dx
C
C HOW MANY X AND Z NODES?
C
      nx=int(xmax/dx)+1
      if (nx.gt.ixzmax) then
         write(6,*) 'IXZMAX too small for profile in X direction'
         stop
      endif
      nz=int(zmax/dz)+1
      if (nz.gt.ixzmax) then
         write(6,*) 'IXZMAX too small for profile in Z direction'
         stop
      endif
C
C COORDINATES OF NODES IN X AND Z
C X SPACING ASSUMED SAME AS Z SPACING
C ALSO PRESSURE IN ** GPA ** I.E. TIME UNITS S, NOT YR
C
      do ix=1,nx
         x(ix)=float(ix-1)*dx
      enddo
      do iz=1,nz
         z(iz)=(-1.0)*float(iz-1)*dz
         pres(iz)=(-rhos)*g*z(iz)*1.e-9/sy/sy
      enddo
C
      write(6,*)'X nodes,',nx,';  Z nodes,',nz
      write(6,*)'Maximum pressure (GPa):',pres(nz)
      return
      end
C
C----------------------------------------------------------
C      
      subroutine TYMING
      implicit double precision(a-h,o-z)
C
C SUBROUTINE CONTROLS THE TIME LOOP FOR TIME-DEPENDENT EQUATIONS
C
C COMMONS REQUIRED TO COMMUNICATE WITH REST OF PROGRAM
C CALC - TEND,ICALC
C TIME OUT - NTIME,TYME(ITIME)
C
C MAXIMUM ARRAY SIZES
C
      parameter (iamax=10)
      parameter (itimax=2001)
      parameter (ixzmax=2001)
C
C
C
      logical found
      double precision pres(ixzmax)
      double precision x(ixzmax)
      double precision z(ixzmax)
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
C
c      common /calc/itype,timfin,timstp,tend,icalc
      common /calc/tend,icalc
C
      common /timset/timfin,timstp,itype
      common /countr/itym,ntime,iendpr
C 
C PLATE SPREADING HISTORY AND CORNER FLOW GEOMETRY
C 
      double precision ufull(iamax)
      double precision utime(iamax)
      common /vel/ufull,utime,alpha,viss,crust,nu,ialpha
C
      common /struct/dlith,djump,iriftr,ijump
C
C SPREADING RATE HISTORY
C
      double precision tyme(itimax)
      double precision ufulla(itimax)
      double precision ufullt(itimax)
      double precision age(itimax)
      double precision rtyme(itimax)
      common /uftym/ufullt,ufulla,age,tyme,rtyme
C
C
C
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
C
C FOR A STEADY-STATE CALCULATION ITYPE = 0
C
      if(itype.eq.0) then
        icalc=0
        itemp=1
        write(6,*)'icalc:',icalc,ufull
C
C FOR A TIME-DEPENDENT CALCULATION
C
      elseif(itype.eq.1) then
C
C CLAUSE: IF YOU WANT TO RUN ONE TIME DEP CALCULATION,
C EITHER SET TIMSTP=0. OR TIMSTP=TIMFIN
C
        if(itym.eq.1) then
          if(timstp.le.0.) then
            write(6,*)'Time-step must be greater than zero'
            stop
          endif
          if(timstp.gt.timfin) then
            timstp=timfin
          endif
          if(timstp.eq.timfin) then
			tend=timfin
            ntime=1
          else
            ntime=int(timfin/timstp)
          endif
C
C CHECK THERE IS ENOUGH SPACE IN THE ARRAY
C
          if (ntime.gt.itimax) then
            write(6,*)
            write(6,*)'ERROR: Not enough space in time array'
            write(6,*)'itimax is ',itimax,' but ',ntime,' is needed'
            stop
          endif
C
C SET UP A 1-D TIME ARRAY TO COUNT THE TIME
C
          do itime=1,ntime
            tyme(itime)=float(itime)*timstp
          enddo
C             
          if(tyme(ntime).lt.timfin) then
            ntime=ntime+1
            tyme(ntime)=timfin
          endif
C
C ENDIF FOR TIME ARRAY
C
C      endif
C
C ENDIF FOR ITYM=1
C
        endif
C
C SET TEND, IE THE NEXT JUMP FOR THE PROGRAM TO CALCULATE
C
        if(itym.eq.1) then
          icalc=1
          itemp=0
          tend=tyme(itym)
        else
          icalc=1
          itemp=0
          tend=tyme(itym)-tyme(itym-1)
        endif
C      
C TIME-DEPENDENT CALCULATION WITH STARTING STEADY-STATE TEMPERATURE
C FIELD   
C
      elseif(itype.eq.2) then
C
C SET UP ARRAYS FOR SPREADING RATE IN ITIME, ONLY FIRST TIME ROUND
C 
        if(itym.eq.1) then
          call READUJ(nivaru)
C
C SET UP TIME ARRAY TO INCLUDE TIME STEPS AND SPREADING RATE CHANGES
C
          if(timstp.le.0.) then
            write(6,*)'**INCORRECT VALUE FOR TIME STEP, 
     &                 MUST BE GREATER THAN ZERO**'
            stop
          endif
C
          if (timstp.gt.timfin) then
      	    timstp=timfin
          endif
          if (timstp.eq.timfin) then
            ntime=2
            do itimea=1,ntime
              tyme(itimea)=float(itimea-1)*timstp
            enddo
          else
            ntime=int(timfin/timstp)+1
C
C CHECK THERE IS ENOUGH SPACE IN THE ARRAY
C
            if (ntime.gt.itimax) then
              write(6,*)
              write(6,*)'ERROR: Not enough space in time array'
              write(6,*)'itimax is ',itimax,' but ',ntime,' is needed'
              stop
            endif
C
C SET UP A 1-D TIME ARRAY TO COUNT THE TIME
C
            do itimea=1,ntime
              rtyme(itimea)=float(itimea-1)*timstp
            enddo
C             
            if (rtyme(ntime).lt.timfin) then
              ntime=ntime+1
              rtyme(ntime)=timfin
            endif
          endif
          do itimea=1,ntime
            write(6,*)rtyme(itimea),itimea,timstp,timfin
          enddo
C
C NOW COMBINE ARRAYS OF TIME AND SPREADING RATE TO SET UP ALL THE
C TIMESTEPS REQUIRED
C
      ivaru=1
      itimea=1
      maxtim=ntime+nivaru
      do itime=1,maxtim
         if(itime.eq.1) then
           if(age(ivaru).eq.0.)then
              ufullt(itime)=ufulla(ivaru)
              tyme(itime)=rtyme(itimea)
              itimea=itimea+1
              ivaru=ivaru+1
           else
              ufullt(itime)=ufull(1)
              tyme(itime)=rtyme(itimea)
              itimea=itimea+1
            endif
          else
            if (ivaru.gt.nivaru) then
              tyme(itime)=rtyme(itimea)
              ufullt(itime)=ufulla(ivaru-1)
              itimea=itimea+1
C              write(6,*)'Endof spreading file'
              goto 9
            endif
            if(rtyme(itimea).lt.age(ivaru))then
              tyme(itime)=rtyme(itimea)
C
C CONSTANT UFULL
C
              ufullt(itime)=ufulla(ivaru-1)
              itimea=itimea+1
            elseif(rtyme(itimea).eq.age(ivaru))then
              ufullt(itime)=ufullt(itime-1)
              tyme(itime)=age(ivaru)
              itimea=itimea+1
              ivaru=ivaru+1
            else
              ufullt(itime)=ufullt(itime-1)
              tyme(itime)=age(ivaru)
              ivaru=ivaru+1
              ntime=ntime+1
              if (ntime.gt.itimax) then
                write(6,*)'Too many steps, increase itimax in program!'
                stop
              endif
            endif
          endif
 9        continue
        enddo
C
C WRITE OUT SPREADING ARRAY, AND TYME ARRAY TO CHECK
C
c      do itime=1,ntime
c         write(6,*)'Spreading Rate:',ufullt(itime),'Time:',tyme(itime)
c      enddo
C
C
C ENDIF FOR THE FIRST TIME ROUND
C
      endif
C
C FIRST CALCULATION FOR THIS ITYPE=2 IS STEADY-STATE
C
      if(itym.eq.1) then
         ufull(1)=ufullt(itym)
         icalc=0 
         itemp=1 
         write(6,*)'Doing steady-state calculation!'
         write(6,*)'icalc',icalc
      else
         ufull(1)=ufullt(itym)
         tend=tyme(itym)-tyme(itym-1)
         icalc=1
         itemp=0
      endif
C
C FOR A TIME DEPENDENT CALCULATION WITH VARIABLE SPREADING RATE
C NOT STARTING FROM A STEADY-STATE
C 
      elseif(itype.eq.3) then
      	if (itym.eq.1) then
C          call READUW(nivaru)
          call READUJ(nivaru)
C
C SET UP TIME ARRAY TO INCLUDE TIME STEPS AND SPREADING RATE CHANGES
C
      if(timstp.le.0.) then
      write(6,*)'**INCORRECT VALUE FOR TIME STEP, 
     & MUST BE GREATER THAN ZERO**'
      stop
      endif
      if(timstp.gt.timfin) then
      timstp=timfin
      endif
      if(timstp.eq.timfin) then
         ntime=1
      else
         ntime=int(timfin/timstp)
C
C CHECK THERE IS ENOUGH SPACE IN THE ARRAY
C
          if (ntime.gt.itimax) then
            write(6,*)
            write(6,*)'ERROR: Not enough space in time array'
            write(6,*)'itimax is ',itimax,' but ',ntime,' is needed'
            stop
          endif
C
C SET UP A 1-D TIME ARRAY TO COUNT THE TIME
C
         do itimea=1,ntime
            rtyme(itimea)=float(itimea)*timstp
         enddo
C             
         if(rtyme(ntime).lt.timfin) then
           ntime=ntime+1
           rtyme(ntime)=timfin
         endif
      endif
C
C NOW COMBINE ARRAYS OF TIME AND SPREADING RATE TO SET UP ALL THE
C TIMESTEPS REQUIRED
C
      ivaru=1
      itimea=1
      maxtim=ntime+nivaru
      do itime=1,maxtim
         if(itime.eq.1) then
           if(age(ivaru).eq.0.)then
              ufullt(itime)=ufulla(ivaru)
              tyme(itime)=rtyme(itimea)
              itimea=itimea+1
              ivaru=ivaru+1
C          write(6,*)'Starting Rate'
C          write(6,*)'Spreading Rate:',ufullt(itime),'Time:',tyme(itime)
           else
              ufullt(itime)=ufull(1)
              tyme(itime)=rtyme(itimea)
              itimea=itimea+1
C          write(6,*)'Starting Rate'
C          write(6,*)'Spreading Rate:',ufullt(itime),'Time:',tyme(itime)
           endif
        else
C
C IF ITIME IS GREATER THAN 1, ITIME IS COUNTER
C
          if(ivaru.gt.nivaru) then
            tyme(itime)=rtyme(itimea)
            ufullt(itime)=ufulla(ivaru-1)
            itimea=itimea+1
C          write(6,*)'End of spreading rate changes'
C          write(6,*)'Spreading Rate:',ufullt(itime),'Time:',tyme(itime)
C            write(6,*)'Endof spreading file'
            goto 8
          endif
          if(rtyme(itimea).lt.age(ivaru))then
            tyme(itime)=rtyme(itimea)
C
C CONSTANT UFULL
C
            ufullt(itime)=ufulla(ivaru-1)
            itimea=itimea+1
C          write(6,*)'ivaru:',ivaru
C          write(6,*)'Spreading rate:',ufulla(ivaru-1),ufulla(ivaru)
C          write(6,*)'Spreading Rate:',ufullt(itime),'Time:',tyme(itime)
          elseif(rtyme(itimea).eq.age(ivaru))then
            ufullt(itime)=ufullt(itime-1)
            tyme(itime)=age(ivaru)
            itimea=itimea+1
            ivaru=ivaru+1
C          write(6,*)'Spreading rate change?',ivaru,ufulla(ivaru-1)
C          write(6,*)'Spreading Rate:',ufullt(itime),'Time:',tyme(itime)
          else
            write(6,*)'Calling last if on array change.'
            ufullt(itime)=ufullt(itime-1)
            tyme(itime)=age(ivaru)
            ivaru=ivaru+1
            ntime=ntime+1
C          write(6,*)'Spreading Rate:',ufullt(itime),'Time:',tyme(itime)
            if(ntime.gt.itimax) then
               write(6,*)'Too many steps, increase itimax in program!'
               stop
            endif
          endif
        endif
 8     continue
      enddo
C
C WRITE OUT SPREADING ARRAY, AND TYME ARRAY TO CHECK
C
c      do itime=1,ntime
c         write(6,*)'Spreading Rate:',ufullt(itime),'Time:',tyme(itime)
c      enddo
C
C
C ENDIF FOR THE FIRST TIME ROUND
C
      endif
C
C SET TEND, IE THE NEXT JUMP FOR THE PROGRAM TO CALCULATE
C
       if(itym.eq.1) then
         icalc=1
         itemp=0
         tend=tyme(itym)
         ufull(1)=ufullt(itym)
C         write(6,*)'set itemp = 0'
       else
         icalc=1
         itemp=0
C         write(6,*)'set itemp = 0'
         tend=tyme(itym)-tyme(itym-1)
         ufull(1)=ufullt(itym)
       endif
C
C CLAUSE TO PICK UP AN INCORRECT OPTION FOR ITYPE
C
      else
      write(6,*)'**INCORRECT OPTION FOR ITYPE:0-3**',itype
      stop
      endif
C      
C RETURN TO PROGRAM
      return
      end
C
C----------------------------------------------------------
C
      subroutine READUW(nivaru)
      implicit double precision(a-h,o-z)
C
C READ IN SPREADING RATE HISTORY: RACHEL WALTERS' VERSION
C
      logical found
C
C MAXIMUM ARRAY SIZES
C
      parameter (iamax=10)
      parameter (itimax=2001)
C
C SPREADING RATE HISTORY
C
      double precision tyme(itimax)
      double precision ufulla(itimax)
      double precision ufullt(itimax)
      double precision age(itimax)
      double precision rtyme(itimax)
      common /uftym/ufullt,ufulla,age,tyme,rtyme
C
C FIND AND OPEN INPUT FILE
C
c      write(6,*)
c      write(6,*) 'Entering READUW'
     	inquire (file='wj_spreading_rates',exist=found)
     	if (.not.found) then
       	write(6,*) 'File "wj_spreading_rates" not found.'
	     	stop
      endif
c      write(6,*)'READ in spreading array'
      open(5,file='wj_spreading_rates')
C
C UFULLA AND AGE ARE THE READ IN ARRAYS OF FULL SPREADING
C RATES AND TIMES
C
      ivaru=1
 4    read(5,*,end=7) ufulla(ivaru), age(ivaru)
      ivaru=ivaru+1
      goto 4
 7    nivaru=ivaru-1
      close(5)
C
C CONVERT UNITS
C (WRITE OUT THE ARRAYS TO CHECK)
C
      do ivaru=1,nivaru
        ufulla(ivaru)=ufulla(ivaru)*1.e-3
        age(ivaru)=age(ivaru)*1.e6
c        write(6,*)ivaru,' Spreading Rate:', ufulla(ivaru),
c     &            '  Age:',age(ivaru)
      enddo
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine READUJ(nivaru)
      implicit double precision(a-h,o-z)
C
C READ IN SPREADING RATE HISTORY: STEVE JONES' VERSION
C
C MAXIMUM ARRAY SIZES
C
      parameter (iamax=10)
      parameter (itimax=2001)
C
C SPREADING RATE HISTORY
C
      double precision tyme(itimax)
      double precision ufulla(itimax)
      double precision ufullt(itimax)
      double precision age(itimax)
      double precision rtyme(itimax)
      common /uftym/ufullt,ufulla,age,tyme,rtyme
C 
C PLATE SPREADING HISTORY AND CORNER FLOW GEOMETRY
C 
      double precision ufull(iamax)
      double precision utime(iamax)
      common /vel/ufull,utime,alpha,viss,crust,nu,ialpha
C
c      write(6,*)
c      write(6,*) 'Entering READUJ'
C
C INPUT DATA IS IN UFULL AND UTIME
C UFULLA AND AGE ARE THE READ IN ARRAYS OF FULL SPREADING
C RATES AND TIMES
C
      if (nu.eq.1) then
        nivaru=2
        ufulla(1)=ufull(1)
		    age(1)=0.0
        ufulla(2)=ufull(1)
		    age(2)=1000.0e6
C
      else
        nivaru=nu
        do ivaru=1,nivaru
          ufulla(ivaru)=ufull(ivaru)
		  		age(ivaru)=utime(ivaru)
        enddo
      endif
C
C (WRITE OUT THE ARRAYS TO CHECK)
C
c      do ivaru=1,nivaru
c        write(6,*)ivaru,' Spreading Rate:', ufulla(ivaru),
c     &            '  Age:',age(ivaru)
c      enddo
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine TILOOP(n1,re1)
      implicit double precision(a-h,o-z)
C 
      parameter (iamax=10)
      parameter (ixzmax=2001)
      parameter (itimax=2001)
      double precision tem(3,ixzmax,ixzmax),tems(ixzmax,ixzmax),
     &                 temc(ixzmax,ixzmax)
      double precision x(ixzmax),z(ixzmax),pres(ixzmax)
C
C TIMING
C
      double precision tyme(itimax)
      double precision ufulla(itimax)
      double precision ufullt(itimax)
      double precision age(itimax)
      double precision rtyme(itimax)
      common /uftym/ufullt,ufulla,age,tyme,rtyme
C
C
C
      common /temo/tem,tems
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
      common /countr/itym,ntime,iendpr
      common /timset/timfin,timstp,itype
c      common /calc/itype,timfin,timstp,tend,icalc
      common /calc/tend,icalc
C 
C PLATE SPREADING HISTORY AND CORNER FLOW GEOMETRY
C 
      double precision ufull(iamax)
      double precision utime(iamax)
      common /vel/ufull,utime,alpha,viss,crust,nu,ialpha
C
      common /struct/dlith,djump,iriftr,ijump
      common /temcor/temc
      common /temtes/itest
C
C RESIDUE
C
      double precision rest(ixzmax,ixzmax),res(3,ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C MELT FRACTION
C
      double precision frac(ixzmax,ixzmax),nfrac(ixzmax,ixzmax)
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C

C FOR 1 STEADY-STATE CALCULATION ONLY:
C 
      if(itype.eq.0) then
      iendpr=1
C      
C      
C FOR TIME DEPENDENT CALCULATION ONLY:
C NFRAC IS NEW MELT FRACTION, THAT IS CORRECTED FOR RESIDUE
C SET RE0 RESIDUE PROFILE AS STARTING RESIDUE PROFILE
C
      elseif(itype.eq.1 .or. itype.eq.2 .or. itype.eq.3) then
         if(itym.eq.ntime) then
         iendpr=1
           do ix=1,nx
             do iz=1,nz
               rest(ix,iz)=res(re1,ix,iz)-nfrac(ix,iz)
C               write(6,*)ix,iz,res(re1,ix,iz),rest(ix,iz)
c               res(1,ix,iz)=rest(ix,iz)
c               res(2,ix,iz)=rest(ix,iz)
c               res(3,ix,iz)=rest(ix,iz)
             enddo
           enddo
C
C 11/12/2008
C CHANGE: ADDITIONAL SWITCH TO SHIFT TEMPERATURE FIELD HORIZONTALLY FOR
C NEW RIDGE JUMP.
C CRUST STAYS SAME AS BEFORE AS IF REACHED STEADY STATE SHOULD BE JUMPING
C INTO SAME THICKNESS OF CRUST.
C
         elseif(iriftr.eq.1 .and. itym.eq.ijump) then
          write(6,*)'Ridge Shift:',djump*1.e-3
C FIND NODE IX JUMPING TO, CALL JX
           jx=0
           do ix=1,nx
             if (x(ix).eq.djump) jx=ix-1
           enddo
C DEFINE NEW STARTING TEMEPRATURE ARRAY AS SHIFTED
C CORRECTED TEMEPRATURE FIELD
          do iz=1,nz
           do ix=1,nx
            if(ix.le.nx-jx) then
            tems(ix,iz)=temc(ix+jx,iz)
            rest(ix,iz)=res(re1,ix,iz)-nfrac(ix,iz)
C            rest(ix,iz)=rest(ix+jx,iz)
C CAN'T DO IN SAME DO LOOP AS IX+JX NOT CALCULATED YET
            endif
           enddo
          enddo
		  do iz=1,nz
               do ix=1,nx
		  rest(ix,iz)=rest(ix+jx,iz)
		   enddo
	      enddo
C now either reset the grid to be smaller horizontally, or
c define a temperature array for missing stuff.
C RESET THE GRID TO BE SMALLER.
          xmax=xmax-djump
          jnx=nx-jx
          nx=int(xmax/dx)+1
C          write(6,*)'new nx',nx,'check same',jnx
          itym=itym+1
          iendpr=0
         else
         itym=itym+1
C         n1=itest
         do ix=1,nx
          do iz=1,nz
           tems(ix,iz)=temc(ix,iz)
           rest(ix,iz)=res(re1,ix,iz)-nfrac(ix,iz)
C           tems(ix,iz)=tem(n1,ix,iz)
          enddo
         enddo
         iendpr=0
      endif
      endif
C
C COULD PUT IN CODE TO WRITE THE STUFF YOU WANT OUT AT
C EVERY TIME STEP
C
      return
      end
C
C
C----------------------------------------------------------
C
      subroutine SETVEL
      implicit double precision(a-h,o-z)
C
C CALCULATE VELOCITES AT EACH NODE FROM CORNER FLOW MODEL
C HORIZONTAL COMPONENT IN VX(X,Z)
C VERTICAL COMPONENT IN VZ(X,Z)
C
      parameter (iamax=10)
      parameter (ixzmax=2001)
      double precision pres(ixzmax),vx(ixzmax,ixzmax),
     &     vz(ixzmax,ixzmax),x(ixzmax),z(ixzmax),dpda(ixzmax,ixzmax),
     &     dtrda(ixzmax,ixzmax),dfda(ixzmax,ixzmax)
C
C VELOCITY GRID STAGGERED WITH RESPECT TO TEMPERATURE GRID
C
      double precision vxs(ixzmax,ixzmax),vzs(ixzmax,ixzmax)
      common /velstg/vxs,vzs
C
C INPUT PARAMETERS SPECIFING VELOCITY FIELD
C
      double precision ufull(iamax)
      double precision utime(iamax)
      common /vel/ufull,utime,alpha,viss,crust,nu,ialpha
C
      common /struct/dlith,djump,iriftr,ijump
C
C COMMUNICATION OF VELOCITY FIELD WITHIN PROGRAM
C
      common /velo/vx,vz,uhalf,vmax
C
C OTHER PARAMETERS 
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
      common /countr/itym,ntime,iendpr
C
C COMMUNICATION OF DERIVATIVES THROUGHOUT THE PROGRAM AND
C VARIOUS PARAMETERISATIONS
C
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
C
C COMMUNICATION OF CORNER FLOW CONSTANTS
C
      common /corner/a,b
C
C WRITE TO SCREEN
C
      if (itym.eq.1) then
        write(6,*)
        write(6,*)'SETVEL: Starting velocity field'
      endif
C
C HALF SPREADING RATE
C
      uhalf=ufull(1)/2.
C
C IF IALPHA=1, ESTIMATE WEDGE ANGLE USING SCALING ARGUMENT
C IN APPENDIX B OF SPIEGELMAN & MCKENZIE, EPSL, 1987.
C THE VALUE OF PHI0 MAKES LITTLE DIFFERENCE TO THE
C CONCLUSIONS OF BOWN & WHITE 1994 - VARYING PHI0 FROM
C 0-0.03 CHANGES THE WEDGE ANGLE BY ABOUT A DEGREE.
C
      phi0=0.0
      if (ialpha.eq.1) then
         if (uhalf.eq.0.) then
           write(6,*)'Cannot use UFULL=0. with IAPHA=1'
           write(6,*)'(results in division by 0.)'
           stop
         endif
         sl=sqrt(viss*uhalf/(1.-phi0)/drho/g)
         write(6,*)'Melt extraction length scale (km):',sl*1.e-3
         ain=0.
 1       b=2./(pi-2.*ain-sin(2.*ain))
C         write(6,*)'b',b
         slr=sl*sqrt(2.*b)
C         write(6,*)'slr',slr
         dlr=2.*sqrt(pi*tdiff*slr/uhalf)
C         write(6,*)'dlr',dlr
         alpha=atan2(dlr,slr)
C         write(6,*)'alpha',alpha
         if (abs(alpha-ain).gt.1.e-3) then
            ain=alpha
            goto 1
         endif
         write(6,*)'Wedge angle from half space cooling:',alpha/rad
         write(6,*)'Wedge angle from half space cooling(r):',alpha
      else if (ialpha.eq.0) then
         if (itym.eq.1) write(6,*)'Wedge angle input:',alpha/rad
      else
         if (itym.eq.1) write(6,*)'Incorrect option for IALPHA:',ialpha
      endif
C
C CONSTANTS FOR CORNER FLOW SOLUTION
C
      b=2.*uhalf/(pi-2.*alpha-sin(2.*alpha))
      a=b*sin(alpha)*sin(alpha)
      vmax=(-2.0)*uhalf*(sin(alpha)*sin(alpha)-1.)/
     &     (pi-2.*alpha-sin(2.*alpha))
      if (itym.eq.1) write(6,*)'Corner flow constants A & B:',a,b
C
C CHECKED AS NEGATIVE WITH ORIGINAL EQUATION, MUST HAVE BEEN
C IN A DIFFERENT COORDINATE SYSTEM WHERE Z INCREASES POSITIVELY
C DOWNWARDS, MEANING THAT ANY UPWELLING VELOCITY WOULD BE NEGATIVE
C SINCE OUR COORDINATE SYSTEM HAS Z INCREASING UPWARDS, THE 
C UPWELLING RATE SHOULD BE POSITIVE, THEREFORE THE EQUATION FOR
C VMAX AND VZ HAVE BEEN MULTIPLIED BY -1.
C
      if (itym.eq.1) 
     &  write(6,*)'Maximum upwelling rate (km/Myr):',vmax*1.e3
C
C DETERMINE VELOCITY COMPONENTS IN HORIZONTAL AND VERTICAL DIRECTIONS
C ON NON-STAGGERED GRID, 
C I.E. VELOCITY NODES ARE SAME AS TEMPERATURE NODES
C EQUATIONS FROM JULL'S PHD DISSERTATION
C NOTE Z-AXIS OPPOSITE TO BOWN & WHITE
C LOOP OVER X AND Z
C OUTPUT VELOCITY FIELD TO bw94_vel_out
C
C      crust=0.
cc      open(7,file='bw94_vel_out')
cc      write(7,3)uhalf*1.e3,alpha/rad
 3    format ('uhalf:',f5.1,' angle:',f5.1)
 2    format (4(f20.2,1x))
      do ix=1,nx
         do iz=1,nz
C
C VELOCITIES AT DEPTHS BELOW BASE OF CRUST
C
            if (z(iz).lt.-crust) then
C
C DISTANCE IN X DIRECTION ON NON-STAGGERED GRID
C
               xtmp=x(ix)
C
C DEPTH BELOW BASE OF CRUST ON NON-STAGGERED GRID
C N.B. Z0 IS +VE
C ALSO X (DISTANCE) SQUARED TIMES Z (DEPTH) SQUARED
C
               z0=-z(iz)-crust
               xspzs=xtmp*xtmp+z0*z0
C
C HORIZONTAL VELOCITY ON NON-STAGGERED GRID
C
               vx(ix,iz)=b*(atan2(xtmp,z0)-xtmp*z0/xspzs)
               if (vx(ix,iz).gt.uhalf) vx(ix,iz)=uhalf
C
C VERTICAL VELOCITY ON NON-STAGGERED GRID
C VZ HAS BEEN MULTIPLIED BY MINUS ONE SINCE IT WOULD BE NEGATIVE
C IN THE ORIGINAL COORDINATE SYSTEM THAT THESE EQNS WERE WRITTEN
C IN, OUR HAS DEPTH INCREASING UPWARDS (ALL DEPTHS ARE -VE) SO
C THE VERTICAL VELOCITY SHOULD BE POSITIVE. 
C
               vz(ix,iz)=-a+b*z0*z0/xspzs
               if (vz(ix,iz).lt.0.) vz(ix,iz)=0.
C
C VELOCITY VALUES IN THE CRUST.
C
            else
               vx(ix,iz)=uhalf
               vz(ix,iz)=0.
            endif
C
C WRITE TO FILE FOR PLOTTING
C
cc            write(7,2)x(ix)*1.e-3,z(iz)*1.e-3,
cc     &           vx(ix,iz)*1.e3,vz(ix,iz)*1.e3
         enddo
      enddo
cc      close (7)
cc      write(6,*)'Written velocity field to "bw94_vel_out" (km, km/Myr)'
C
C NOW CALUCALTE VELOCITY FIELD ON STAGGERED GRID
C USED IN PATANKAR'S FINITE VOLUME METHOD (1980)
C THE STAGGERED GRID IS REQUIRED TO STABLIZE CALCULATION OF AN
C ARBITRARY VELOCITY FIELD 
C
C DISTANCE IN X DIRECTION ON STAGGERED GRID
C
      do ix=1,nx
         xtmp=x(ix)
         xtmps=xtmp+0.5*dx
C 
C DEPTH ON STAGGERED GRID
C
         do iz=1,nz
            zs=z(iz)-0.5*dz
C
C VELOCITIES AT DEPTHS BELOW BASE OF CRUST
C
            if (zs.lt.-crust) then
C
C DEPTH BELOW BASE OF CRUST ON STAGGERED GRID
C N.B. Z0 IS +VE
C
               z0=-z(iz)-crust
               z0s=-zs-crust
C
C STAGGERED HORIZONTAL VELOCITY GRID
C HORIZONTAL DISTANCES ARE MIDWAY BETWEEN NODES USED FOR TEMPERATURE CALC
C VERTICAL DEPTHS ARE SAME AS NODES USED FOR TEMPERATURE CALC
C SEE PATANKAR (1980) FIGURE 6.6
C HORIZONTAL VELOCITY VXS(IX,IZ) AT DISTANCE X=X(IX)+DX/2 AND DEPTH Z=Z(IZ)
C
               xspzs=xtmps*xtmps+z0*z0
               vxs(ix,iz)=b*(atan2(xtmps,z0)-xtmps*z0/xspzs)
               if (vxs(ix,iz).gt.uhalf) vxs(ix,iz)=uhalf
C
C STAGGERED VERTICAL VELOCITY GRID
C HORIZONTAL NODES ARE SAME AS NODES USED FOR TEMPERATURE CALC
C VERTICAL NODES ARE MIDWAY BETWEEN NODES USED FOR TEMPERATURE CALC
C SEE PATANKAR (1980) FIGURE 6.6
C VERTICAL VELOCITY VZS(IX,IZ) AT DISTANCE X=X(IX) AND DEPTH Z=Z(IZ)-DZ/2
C
               xspzs=xtmp*xtmp+z0s*z0s
               vzs(ix,iz)=-a+b*z0s*z0s/xspzs
               if (vzs(ix,iz).lt.0.) vzs(ix,iz)=0.
C
C VELOCITY VALUES IN THE CRUST
C
            else
               vzs(ix,iz)=0.
            endif
            if (z(iz).ge.-crust) vxs(ix,iz)=uhalf
C
C END OF STAGGERED SET-UP GRID LOOPS
C
         enddo
      enddo
C
C OPEN OUTPUT FILE
C
cc      open(8,file='bw94_vel_stag_out')
 5    format (3(f10.2,1x))
cc      write(8,3)uhalf*1.e3,alpha/rad
C
C WRITE HORIZONTAL STRIPS
C
cc      do iz=1,nz
cc        do ix=1,nx
cc            write(8,5)(x(ix)+dx/2.)*1.e-3,z(iz)*1.e-3,
cc     &           vxs(ix,iz)*1.e3
cc        enddo
cc      enddo
C
C WRITE HORIIZONTAL STRIPS
C
cc      do ix=1,nx
cc        do iz=nz,1,-1
ccc            write(6,*)ix,'X',iz,'Z',vz(ix,iz)*1.e3,
ccc     &                       vzs(ix,iz)*1.e3,vz(ix,iz+1)*1.e3
cc            write(8,5)x(ix)*1.e-3,(z(iz)-dz/2.)*1.e-3,
cc     &           vzs(ix,iz)*1.e3
cc        enddo
cc      enddo
C
C CLOSE FILE
C
cc      close (8)
cc      write(6,*)
cc     & 'Written staggered velocity field to "bw94_vel_stag_out"'
C
C *** THIS WILL HAVE TO BE CALCULATED EVERY TIME STEP EVENTUALLY ***
C CALCULATE DP/Dt FROM THE VELOCITY FIELD, FOR USE IN CALCULATING
C THE RATE OF CHANGE IN MELT FRACTION, NEGATIVE BECAUSE EVERY
C PARTICLE IS UPWELLING. P=RHO*G*(-Z), DIFFERENTIATE TO GET
C DP/DZ=-RHO*G, HOWEVER, I THINK YOU IGNORE THE -VE Z AS THIS IS
C A FUNCTION OF THE COORDINATE SYSTEM AND NOT OF THE ORIGINAL
C EQUATION. 
C UNITS ARE ALSO CORRECT TO GET GPA/YR
C
cc      open (15,file='bw94_dpda_out')
 4    format (4(f20.15,1x))
      do ix=1,nx
         do iz=1,nz
            dpda(ix,iz)=(-vz(ix,iz))*rhos*g*1.e-9/sy/sy
cc            write(15,4)x(ix)*1.e-3,z(iz)*1.e-3,
cc     &              vz(ix,iz)*1.e3,dpda(ix,iz)
         enddo
      enddo
cc      close(15)
cc      write(6,*)'Calculated pressure-time derivatives DP/Da'
cc      write(6,*)'Written pressure derivatives to "bw94_dpda_out"'
C
      return
      end
Cpres(iz)=-rhos*g*z(iz)*1.e-9/sy/sy
C
C----------------------------------------------------------
C
      subroutine SETRES
C
C SET STARTING RESIDUE GRID AS NO MELTING HAS OCCURRED
C SO ALL NODES HAVE 1 RESIDUE
C
C REAL DECLARATION
C
      implicit double precision(a-h,o-z)
C
C PARAMETER SETTING FOR NO. OF GRID NODES IN A DIRECTION
C
      parameter (ixzmax=2001)
C
C GRID
C
      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
C
C RESIDUE
C RES - ARRAYS REQUIRED TO CALCULATE THE NEW RESIDUE STRUCTURE
C REST - STARTING RESIDUE (ZERO TO START WITH)
C
      double precision res(3,ixzmax,ixzmax),rest(ixzmax,ixzmax)
      common /resdue/rest,res,resc
C 
C ITYM COUNTER
C
      common /countr/itym,ntime,iendpr
C
C TEMPORARY COMMON
C
      double precision p(ixzmax,ixzmax)
      common /gauss/p
C
C SET ALL NODES EQUAL TO 1 IN THE FIRST TIME STEP
C
      if (itym.eq.1) then
        write(6,*)'SETTING RESIDUE FIELD'
        x0=100.
        z0=-50.
        open(91,file='starting_residue')
        do ix=1,nx
          do iz=1,nz
C           p(ix,iz)=1./(2.*pi*1.)*exp(-(abs(x(ix)*1.e-3-x0)**2.+ 
C     &     abs(z(iz)*1.e-3-z0)**2.)/(2.*1.*1.))
            p(ix,iz)=exp((-1.0)*(abs(x(ix)*1.e-3-x0)**2.+ 
     &      abs(z(iz)*1.e-3-z0)**2.)/500.)
C           write(6,*)'x,z,p',x(ix),z(iz),p(ix,iz),xd1,zd1
            if (p(ix,iz).lt.0.00001) then
              p(ix,iz)=0.
            endif
C
C CODE LINE FOR GAUSSIAN TEST CIRCLE
C           rest(ix,iz)=1.-p(ix,iz)
C CODE LINE FOR RUNNING PROGRAM PROPERLY
            rest(ix,iz)=1.
           write(91,*)x(ix)*1.e-3,z(iz)*1.e-3,rest(ix,iz)
        enddo
      enddo
      close(91)
      write(6,*)'END OF SETTING RESIDUE - STOP'
      endif
C
C REST IS RESET AT THE END OF THE MELTING CALCULATION AND
C SAVED FOR THE NEXT TIMESTEP.
C
C END SUBROUTINE
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine SETTP
C
C SET STARTING POTENTIAL TEMPERATURE GRID
C TO BE USED TO ESTIMATE COMPOSITION
C ALL NODES BEGIN AT STARTING POTENTIAL TEMPERATURE
C
C REAL DECLARATION
C
      implicit double precision(a-h,o-z)
C
C MAX NUMBER OF GRID NODES
C
      parameter (ixzmax=2001)
C
C GRID
C
      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
C
C POTENTIAL TEMPERATURE
C  TP, POTENTIAL TEMPERATURE
C  TPS, STARTING POT TEMP
C
      double precision tp(3,ixzmax,ixzmax)
      double precision tps(ixzmax,ixzmax)
      common /pottem/tp,tps
C
C STARTING POTENTIAL TEMP
C
      character*50 tfilei,tfileo
      common /tem/ptem,tfilei,tfileo,item
C
C INITIALIZE POTENTIAL TEMPERATURE GRID
C
      write(6,*)'SETTP: Initializing potential temperature grid'
      do ix=1,nx
        do iz=1,nz
          tps(ix,iz)=ptem
          tp(1,ix,iz)=ptem
          tp(2,ix,iz)=ptem
          tp(3,ix,iz)=ptem
        enddo
      enddo
C
      return
      end
C
C
C----------------------------------------------------------
C
      subroutine SETTEM
C
C SET STARTING TEMPERATURE FIELD
C  ITEM=0, HALFSPACE COOLING MODEL
C  ITEM=1, READ IN PRE-EXISTING FILE
C  ITEM=2, EQUILIBRIUM CONDUCTIVE GEOTHERM (NO ADIABATIC GRADIENT)
C  ITEM=3, ZERO EVERYWHERE
C 
C STARTING TEMPERATURE FIELD STORED IN TEMS(X,Z)
C ALSO COPIED INTO FIELD USED FOR CALCULATION, TEM(X,Z)
C
      implicit double precision(a-h,o-z)
      parameter (iamax=10)
      parameter (ixzmax=2001)
      double precision pres(ixzmax),
     &   vx(ixzmax,ixzmax),vz(ixzmax,ixzmax),x(ixzmax),z(ixzmax)
      character*50 tfilei,tfileo
c      character*3 cownt
      logical found
C
C INPUT PARAMETERS SPECIFING TEMPERATURE FIELD
C
      common /tem/ptem,tfilei,tfileo,item
C
C COMMUNICATION OF TEMPERATURE FIELD WITHIN PROGRAM
C
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      common /temo/tem,tems
C
C OTHER PARAMETERS 
C
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C 
C PLATE SPREADING HISTORY AND CORNER FLOW GEOMETRY
C 
      double precision ufull(iamax)
      double precision utime(iamax)
      common /vel/ufull,utime,alpha,viss,crust,nu,ialpha
C
      common /struct/dlith,djump,iriftr,ijump
      common /velo/vx,vz,uhalf,vmax
      common /countr/itym,ntime,iendpr
C
C RESIDUE
C
      double precision rest(ixzmax,ixzmax),res(3,ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C WRITE TO SCREEN
C
cc      write(6,*)
cc      write(6,*)'SETTEM: Starting temperature field'
C
C FIRST TIME SUBROUTINE IS CALLED
C
      if (itym.eq.1) then
C
C SURFACE BOUNDARY CONDITION, TEMPERATURE ZERO
C BOTTOM BOUNDARY CONDITION FROM POTENTIAL TEMP
C
        t0=0.
        t1=(ptem+273.)/exp(g*texps*z(nz)/cp)-273.
        gadiab=g*texps/cp
        write(6,*)'Surface and Basal temperatures (oC)',t0,t1
        write(6,*)'Adiabatic gradients:',(t0+273.)*gadiab,
     &                                 (t1+273.)*gadiab
        do ix=1,nx
          tems(ix,1)=t0
          tem(1,ix,1)=t0
          tem(2,ix,1)=t0
          tem(3,ix,1)=t0
          tems(ix,nz)=t1
          tem(1,ix,nz)=t1
          tem(2,ix,nz)=t1
          tem(3,ix,nz)=t1
        enddo
C
C ITEM=1, READ IN EXISTING TEMPERATURE PROFILE IF REQUIRED
C AND THE ACCOMPANYING STARTING RESIDUE DEPLETION VALUES
C NB THE X AND Z VALUES MUST MATCH THE PRESENT GRID
C
        if (item.eq.1) then         
          inquire (file=tfilei,exist=found)
          if (.not.found) then
            write(6,*) 'Temperature/Residue file not found: ',tfilei
            stop
          endif
          open (20,file=tfilei)
          write(6,*)
     & 'Read in starting temperature and residue depletion from file: ',
     & tfilei
          do iz=1,nz
            do ix=1,nx
              read(20,*)xtmp,ztmp,tem(1,ix,iz),rest(ix,iz)
              xdif=abs(xtmp-x(ix)/1000.)
              zdif=abs(ztmp-z(iz)/1000.)
c              write(6,*)'xdif,zdif:',xtmp,ztmp,tem(1,ix,iz)
              if (xdif.gt.1.e-4 .or. zdif.gt.1.e-4) then
                write(6,*)'Input file not correctly ordered',ix,iz
                write(6,*)'  Expected:',x(ix)/1000.,z(iz)/1000.
                write(6,*)'  Read:',xtmp,ztmp
                stop
              endif
              tems(ix,iz)=tem(1,ix,iz)
              tem(2,ix,iz)=tem(1,ix,iz)
              tem(3,ix,iz)=tem(1,ix,iz)
              res(1,ixzmax,ixzmax)=rest(ix,iz)
              res(2,ixzmax,ixzmax)=rest(ix,iz)
              res(3,ixzmax,ixzmax)=rest(ix,iz)
            enddo
          enddo
C
C ITEM=0, SET TEMPERATURE STRUCTURE FROM HALFSPACE COOLING MODEL
C
        else if (item.eq.0) then
          do ix=1,nx
            age=x(ix)/uhalf
            sl=2.*sqrt(tdiff*age)
            do iz=2,nz-1
              zd=(-z(iz))/sl
C               if (ix.eq.2 .and. iz.eq.100) then
C               write(6,*)'zd stuff',x(ix),z(iz),zd
C               endif
              tem(1,ix,iz)=t1*erf(zd)
              tem(2,ix,iz)=tem(1,ix,iz)
              tem(3,ix,iz)=tem(1,ix,iz)
              tems(ix,iz)=tem(1,ix,iz)
            enddo
          enddo
          write(6,*)'Starting temperature field: Cooling halfspace'
C
C ITEM=2, SET EQUILIBRIUM CONDUCTIVE GEOTHERM
C 11/12/2008
C CHANGE: SET CONDUCTIVE GEOTHERM TO BASE OF LITHOSPHERE AS
C SET BY PARAMETER FILE, THEN ADIABATIC GRADIENT BEYONG THAT
C
C OLD CODE WITHOUT LITHOSPHERE: JUST UNCOMMENT TO RE-IMPLEMENT
C      else if (item.eq.2) then
C         do iz=2,nz-1
C            t2=t1/z(nz)*z(iz)
Cc            write(6,*)t1,z(iz),z(nz),t2
C            do ix=1,nx
C               tem(1,ix,iz)=t2
C               tem(2,ix,iz)=tem(1,ix,iz)
C               tem(3,ix,iz)=tem(1,ix,iz)
C               tems(ix,iz)=tem(1,ix,iz)
CC               write(6,*)'tems, ix, iz:',ix,iz,tems(ix,iz)
C            enddo
C         enddo
C      write(6,*)'Starting temperature field: Conductive geotherm'
C OLD CODE END
C 
        else if (item.eq.2) then
C TEMPERATURE AT THE BASE OF THE LITHOSPHERE
          tbase=(ptem+273.)/exp(g*texps*(-dlith)/cp)-273.
          do iz=2,nz-1
            if (z(iz).ge.-dlith) then
              t2=tbase/(-dlith)*z(iz)
c              write(6,*)t1,z(iz),z(nz),t2
              do ix=1,nx
                tem(1,ix,iz)=t2
                tem(2,ix,iz)=tem(1,ix,iz)
                tem(3,ix,iz)=tem(1,ix,iz)
                tems(ix,iz)=tem(1,ix,iz)
C                write(6,*)'tems, ix, iz:',ix,iz,tems(ix,iz)
              enddo
            else
              t2=(ptem+273.)/exp(g*texps*z(iz)/cp)-273.
              do ix=1,nx
                tem(1,ix,iz)=t2
                tem(2,ix,iz)=tem(1,ix,iz)
                tem(3,ix,iz)=tem(1,ix,iz)
                tems(ix,iz)=tem(1,ix,iz)
C                write(6,*)'tems, ix, iz:',ix,iz,tems(ix,iz)
              enddo
            endif
          enddo
          write(6,*)'Starting temperature field: Conductive geotherm'


C
C ITEM=3, STARTING TEMPERATURE ZERO
C
        else if (item.eq.3) then
          do iz=2,nz-1
            do ix=1,nx
               tem(1,ix,iz)=0.
               tem(2,ix,iz)=tem(1,ix,iz)
               tem(3,ix,iz)=tem(1,ix,iz)
               tems(ix,iz)=tem(1,ix,iz)
            enddo
          enddo
          write(6,*)'Starting temperature field: Zero'
C
C UNKNOWN ITEM VALUE
C
        else
          write(6,*)'Incorrect option for ITEM:',item
          stop
        endif
C
C TIME IF BLOCK, IF IT IS NOT THE FIRST CALCULATION IN A TIME
C DEPENDENT LOOP THEN
C
      elseif (itym.ne.1) then
        do iz=1,nz
          do ix=1,nx
            tem(1,ix,iz)=tems(ix,iz)
            tem(2,ix,iz)=tem(1,ix,iz)
            tem(3,ix,iz)=tem(1,ix,iz)
          enddo
        enddo
      endif
C
C WRITE OUT STARTING TEMPERATURE STRUCTURE:
C
C      write(cownt,'(i3.3)')itym
C      open(15,file='tems_out_'//cownt) 
C      do ix=1,nx
C         do iz=1,nz
C         write(15,*)ix,iz,tems(ix,iz)
C         enddo
C      enddo
C      close(15)
      if(itym.eq.1) then
      open(15,file='tems_out')
       do ix=1,nx
         do iz=1,nz
         write(15,*)x(ix)*1.e-3,z(iz)*1.e-3,tems(ix,iz)
         enddo
      enddo
      close(15)
      endif
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine POPEN
      implicit double precision(a-h,o-z)
c      implicit integer(i-n)
      parameter (ixzmax=2001)
C
C DECIDE WHETHER 1D TEMPERATURE PROFILES ARE TO BE WRITTEN
C IF SO, OPEN 1D TEMPERATURE PROFILE FILES
C
C OUTPUT VERTICAL TEMPERATURE PROFILE
C
      common /prfout/xprf,wprf,zprf,tprf,pfileh,pfilev
      character*50 pfileh,pfilev
      common /iprf/iprfv1,iprfv2,iprfh
C
C COMMUNICATION OF TEMPERATURE FIELD WITHIN PROGRAM
C
      double precision tem(3,ixzmax,ixzmax),tems(ixzmax,ixzmax)
      common /temo/tem,tems
C
C TEMPERATURE GRID
C
      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
      common /grido/pres,x,z,nx,nz
      common /grid/xmax,zmax,dx,dz
C
C WRITE TO SCREEN
C
      write(6,*)
      write(6,*)'POPEN: Set up 1D profiles'
C
C IS 1D GEOTHERM TO BE WRITTEN OUT?
C YES, IF AN OUTPUT FILE NAME IS GIVEN.
C
      jlen=index(pfilev(1:),' ')
      if (jlen.ne.0) then
        write(6,*)'1D vertical profile to be written to file ',pfilev
C
C DETERMINE HORIZONTAL NODES VERTICAL PROFILE TO BE CALCULATED BETWEEN
C (AVERAGE OVER REGION X=XPRF TO X=XPRF+WPRF)
C
        iprfv1=nint(xprf/dx)+1
        iprfv2=nint((xprf+wprf)/dx)+1
        write(6,*)' between',x(iprfv1)/1000.,' km (node',iprfv1,
     &           ') and',x(iprfv2)/1000.,' km (node',iprfv2,')'
        open(10,file=pfilev)
C
C WRITE OUT STARTING TEMPERATURE 1D GEOTHERM
C
 1      format ('> ',f8.2) 
        write(10,1) 0.
C        geo1d=0.
        do iz=1,nz
          geo1d=0.
          do ix=iprfv1,iprfv2
            geo1d=geo1d+tems(ix,iz)
          enddo
          geo1d=geo1d/float(iprfv2-iprfv1+1)
C          write (6,*)z(iz)/1000.,geo1d
          write(10,*)z(iz)/1000.,geo1d
        enddo
C
C IF NO VERTICAL PROFILE TO BE WRITTEN
C
      else
        iprfv1=0
        iprfv2=0
      endif
C
C IS 1D HORIZONTAL PROFILE TO BE WRITTEN OUT?
C YES, IF AN OUTPUT FILE NAME IS GIVEN
C
      jlen=index(pfileh(1:),' ')
      if (jlen.ne.0) then
C
C DETERMINE VERTICAL NODE PROFILE IS TO BE CALCULATED AT
C AND OPEN OUTPUT FILE
C
        write(6,*)'1D horizontal profile to be written to file ',pfileh
        iprfh=nint(zprf/dz)+1
        write(6,*)' at depth',(-z(iprfh))/1000.,' km (node',iprfh,')'
        open(11,file=pfileh)
C
C WRITE OUT STARTING TEMPERATURE 1D HORIZONTAL PROFILE
C
        write(11,1) 0.
        do ix=1,nx
          write(11,1)x(iz)/1000.,tems(ix,iprfh)
        enddo
C
C CLOSE HORIZONTAL PROFILE FILE
C
        close (11)
C
C IF NO HORIZONTAL PROFILE TO BE WRITTEN
C
      else
        iprfh=0
      endif
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine FINISH(n1)
C
C WRITE OUT REMAINING DATA
C CLOSE FILES
C
      implicit double precision(a-h,o-z)
C
C MAXIMUM SIZE OF GRIDS
C
      parameter (ixzmax=2001)
C
C TEMPERATURE GRID
C
      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
      common /grido/pres,x,z,nx,nz
      common /grid/xmax,zmax,dx,dz
C
C TEMPERATURE FIELD
C
      double precision tem(3,ixzmax,ixzmax),tems(ixzmax,ixzmax)
      common /temo/tem,tems
C
C TIME
C
      common /timing/time,it
      double precision time(ixzmax)
C
C WRITE TO SCREEN
C
c      write(6,*)
c      write(6,*)'FINISH: Final output'
C
C CLOSE 1D GEOTHERM FILE
C
 11   format ('>')
      write(10,11) 
      close (10)
c      write(6,*)'Closed 1D profile file'
C
C WRITE OUT STARTING AND FINAL TEMPERATURE STRUCTURES
C NOTE TEM(N1,IX,IZ) IS NOT THE CORRECTED TEMPERATURE
C STRUCTURE, NEEDS TO BE TEMC(IX,IZ)
C
      open(8,file='bw94_tem_out')
 22   format (4(f10.2,1x))
      do ix=1,nx
        do iz=1,nz
           write(8,22)x(ix)*1.e-3,z(iz)*1.e-3,tems(ix,iz),tem(n1,ix,iz)
        enddo
      enddo
      close (8)
c      write(6,*)'Written temperature field to "bw94_tem_out" (km, oC)'
C
      return
      end
C
C------------------------------------------------------------------
C
C      subroutine INTRES
C
CC TEMPRORARY SUBROUTINE TO INTEGRATE UNDER THE RESIDUE
CC TO TEST FOR CONSISTENT AREA
CC
CC GRID PARAMETER
CC
C      parameter (ixzmax=1000)
CC
CC GRID
CC
C      double precision pres(ixzmax,ixzmax),x(ixzmax),z(ixzmax)
C      common /grid/xmax,zmax,dx,dz
C      common /grido/pres,x,z,nx,nz
CC
CC RESIDUE COMMON
CC
C      double precision res(3,ixzmax,ixzmax),rest(ixzmax,ixzmax)
C      common /resdue/rest,res,resc
CC
CC TESTING COMMON
CC
C      common /gauss/area,volume
CC
CC
CC AREA - SET UP IS THAT IF
C      do ix=1,nx
C        do iz=1,nz
C           g1=
C
C------------------------------------------------------------------
C
      subroutine AZMAX(array,ix,m,amax,zmax)
C 
C RETURN DEPTH OF MAXIMUM VALUE OF ARRAY 
C AT X LOCATION IX AND TIME LOCATION M
C
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
C
C GRID COORDINATES
C
      double precision pres(ixzmax)
      double precision x(ixzmax)
      double precision z(ixzmax)
      common /grido/pres,x,z,nx,nz
C
C ARRAY TO BE ANALYZED HERE
C
      double precision array(3,ixzmax,ixzmax)
C
      amax = 1.0d-30
      zmax = 1.0d-30
      do iz=1,nz
        if (array(m,ix,iz).gt.amax) then
          amax = array(m,ix,iz)
          zmax = z(iz)
        endif
c        write(6,*) array(m,ix,iz),z(iz),amax,zmax
      enddo
      return
      end
C
C------------------------------------------------------------------
C
