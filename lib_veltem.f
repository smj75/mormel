C
C ROUTINES FOR CALCULATING VELOCITY AND TEMPERATURE FIELDS
C STEPHEN JONES SUMMER 2005, AUTUMN 2006
C 
C CONTAINS FOLLOWING ROUTINES
C	CALVTM
C	VELCAL
C	TEMFVI
C	FVNB2D
C	TTCON
C	APOW
C	TEMPDE
C	SWAP
C	TCNVGE
C	TRIDAG
C
C----------------------------------------------------------
C
      subroutine CALVTM(n1,re1,ntp0,ntp1)
C
C VELOCITY, TEMPERATURE AND MELTING CALCULATION
C 1. VELOCITY FIELD CAN BE FIXED (E.G. FROM A CORNER
C    FLOW MODEL, OR IS CALCULATED FROM THE TEMPERATURE
C    FIELD AND A VISCOSITY-TEMPERATURE RELATIONSHIP
C    USING SIMPLER ALGORITHM OF PATANKAR (1980).
C 2. TEMPERATURE FIELD CALCULATED EITHER USING A FINITE
C    VOLUME IMPLICIT SCHEME OR AN EXPLICIT SCHEME.
C 3. MELTING CALCULATIONS USE THE PARAMETERIZATIONS OF
C    MCKENZIE & BICKLE (1988) AND KATZ ET AL.
C    UPDATED BY RACHEL.
C
      implicit double precision(a-h,o-z)
C
C MAXIMUM SIZE OF GRIDS
C
      parameter (ixzmax=2001)
      parameter (iemax=20)
C
C TEMPERATURE GRID
C
      double precision pres(ixzmax)
      double precision x(ixzmax)
      double precision z(ixzmax)
      common /grido/pres,x,z,nx,nz
      common /grid/xmax,zmax,dx,dz

      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      double precision temc(ixzmax,ixzmax)
      common /temo/tem,tems
C
C COMMONS REQUIRED FOR CONTINUITY TO ALLOW A NO MELTING CALCULATION
C HOPEFULLY NOT AFFECTING THE NORMAL RUNNING OF THE PROGRAM
C
      common /temcor/temc
      common /vel/ufull,alpha,viss,crust,ialpha
      common /intgra/tarea,tcrust
      common /countr/itym,ntime,iendpr
C
C MELTING
C
      common /melt/ds,xh2o,cpxm,iparam,imelt,itemp
C
C RESIDUE
C
      double precision res(3,ixzmax,ixzmax)
      double precision rest(ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C POTENTIAL TEMPERATURE
C  TP, POTENTIAL TEMPERATURE
C  TPS, STARTING POT TEMP
C
      double precision tp(3,ixzmax,ixzmax)
      double precision tps(ixzmax,ixzmax)
      common /pottem/tp,tps
C
C COMPOSITION
C
      character*2 elabel(iemax)
      character*100 tzfile(iemax)
      common /comp/tzfile,elabel,gdepth,icomp,isourc,nel
C
C TIME
C
      common /timing/time,it
c      common /calc/itype,timfin,timstp,tend,icalc
      common /calc/tend,icalc
      double precision time(ixzmax)
C
C TOLERANCES FOR CONVERGENCE OF RESIDUE & POTENTIAL
C TEMPERATURE FIELDS
C
      parameter (rtolmx=0.00001, rtolme=0.00001)
      parameter (ttolmx=0.01, ttolme=0.01)
C
C WRITE TO SCREEN
C
cc      write(6,*)
cc      write(6,*)'CALVTM: Coupled vel/temp/melting calculation'
C
C SET IVEL HERE FOR NOW
C
      ivel=0
C
C START TIME 
C
      it=1
      time(it)=0. 
C
C INITIATE LABELS FOR BEFORE AND AFTER TEMPERATURES
C
      n0=1
      n1=2
C
C LABELS FOR BEFORE AND AFTER POTENTIAL TEMPERATURES
C
      if (itym.eq.1) then
        ntp0=1
        ntp1=2
      else
        call SWAP(ntp0,ntp1)
      endif
c      write(6,*)'ntp0, ntp1',ntp0,ntp1
C
C INITIATE LABELS FOR BEFORE AND AFTER RESIDUES
C
      re0=1
      re1=2
C
C LOOP FOR TIMESTEP
C
 100  continue
C
C TIME AT END OF THIS TIMESTEP
C
      it=it+1
      time(it)=tend
cc      write(6,*)'Aiming for time',time(it)*1.e-6,' Myr'
C
C IF VELOCITY FIELD IS FIXED (INDEPENDENT OF TEMPERATURE)
C PROCEED DIRECTLY TO TEMPERATURE CALCULATION
C
      if (ivel.eq.0) goto 200
C
C UPDATE VELOCOITY FIELD
C
      write(6,*)'About to call VELCAL'
      call VELCAL(n1)
C
C ITERATE TO FIND VELOCITY FIELD?
C
C
C CALCULATE RESIDUE PROFILE
C FIRST SET RE0 RESIDUE PROFILE AS STARTING RESIDUE PROFILE
C SUBROUTINE NEEDS TO BE EXECUTED BEFORE CALLING THE
C MELTING ROUTINE BECAUSE THE INTEGRATION IS LINKED INTO
C THE MELTING SUBROUTINE AND WE NEED TO FIX THE MELT
C BEFORE INTEGRATING
C
 200  do ix=1,nx
        do iz=1,nz
          res(1,ix,iz)=rest(ix,iz)
          res(2,ix,iz)=rest(ix,iz)
          res(3,ix,iz)=rest(ix,iz)
        enddo
      enddo
      m0=nint(re0)
      m1=nint(re1)
      call ADVFVI(res, rtolmx, rtolme, m0, m1)
      re0=m0
      re1=m1
C
C ADVECT FIELD OF ORIGINAL POTENTIAL TEMPERATURES ROUND GRID
C
      if (icomp.eq.3) then
        call ADVFVI(tp, ttolmx, ttolme, ntp0, ntp1)
c        write(6,*)'ntp0, ntp1',ntp0,ntp1
      endif
C
C CALCULATE TEMPERATURE PROFILE
C
 300  call TEMFVI(n0,n1)
C
C MELTING CALCULATION
C
cc      write(6,*)'imelt:',imelt
C SWITCH TO COPE WITH DYING RIDGE TO HAVE NO MELTING
C      if (itym.eq.11) imelt=0
C      if (imelt.ne.0) then
CR ADDED RE1 TO MELT
      call MEL(n1,re1)
C      else 
C FOR CONTINUITY WHEN NO MELTING IS CALCULATED
C       do ix=1,nx
C         do iz=1,nz
C         temc(ix,iz)=tem(n1,ix,iz)
C         enddo
C       enddo
C      tcrust=crust
C      endif
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine VELCAL(n1)
C
C UPDATE PREVIOUSLY GUESSED VELOCITY FIELD
C TO BE COMPATIBLE WITH NEW TEMPERATURE FIELD
C COMPRISES STEPS 2-6 OF SIMPLER ALGORITHM
C OF PATANKAR (1980)
C
      implicit double precision(a-h,o-z)
C
C MAXIMUM SIZE OF GRID
C
      parameter (ixzmax=2001)
C
C ARRAYS FOR TRIDIAGONAL SYSTEM OF SIMULTANEOUS EQUATIONS
C USED TO SOLVE IMPLICIT DISCRETIZATION EQUATIONS
C
c      double precision a(ixzmax)
c      double precision b(ixzmax)
c      double precision c(ixzmax)
c      double precision r(ixzmax)
c      double precision u(ixzmax)
C
C CALCULATE COEFFICIENTS FROM MOMENTUM EQUATION
C

C
C CALCULATE PSEUDOVELOCITIES
C

C
C COEFFICIENTS FOR PRESSURE EQUATION
C

C
C SOLVE PRESSURE EQUATION TO GET GUESSED PRESSURE FIELD, P*
C

C
C SOLVE MOMENTUM EQUATIONS TO OBTAIN GUESSED VELOCITY FIELDS U*, V*
C

C
C CALCULATE MASS SOURCE AND SOLVE PRESSURE CORRECTION P'
C

C
C CORRECT VELOCITY FIELD
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine TEMFVI(n0,n1)
C
C (1) CALCULATES TEMPERATURE STRUCTURE FROM ADVECTION-DIFFUSION
C PARTIAL DIFFERENTIAL EQUATION:
C DT/Dt = K (D2T/DX2 + D2T/DZ2) - VX DT/DX - VZ (DT/DZ + H).
C USES FINITE VOLUME IMPLICIT DISCRETIZATION METHOD IN PATANKAR (1980)
C
C (2) CALLS ROUTINE FOR MELTING CALUCLATIONS.
C IF ICALC=0, MELTING CALCULATIONS CARRIED OUT ONLY AT END
C OF TEMPERATURE CALCULATIONS.
C IF ICALC=1, MELTING CALCULATIONS CARRIED OUT AT EACH
C TIMESTEP.
C
C THE FOLLOWING CHECKS HAVE BEEN MADE
C   1. WHEN UHALF=0 AND NO ADABATIC GRADIENT, LINEAR CONDUCTIVE
C      GEOTHERM IS RETRIEVED AFTER A TIME LONG IN COMPARISON WITH
C      THE THERMAL TIME CONSTANT
C   2. STARTING FROM ZERO TEMPERATURE FIELD WHEN UHALF IS BIG
C      BASAL TEMPERATURE PROPAGATES TO SURFACE AT UPWELLING SPEED
C
      implicit double precision(a-h,o-z)
C
C MAXIMUM SIZE OF GRID
C
      parameter (ixzmax=2001)
C
C ARRAYS FOR TRIDIAGONAL SYSTEM OF SIMULTANEOUS EQUATIONS
C USED TO SOLVE IMPLICIT DISCRETIZATION EQUATIONS
C
      double precision a(ixzmax)
      double precision b(ixzmax)
      double precision c(ixzmax)
      double precision r(ixzmax)
      double precision u(ixzmax)
C
      double precision pres(ixzmax)
      double precision dpda(ixzmax,ixzmax)
      double precision dtrda(ixzmax,ixzmax)
      double precision dfda(ixzmax,ixzmax)
C
C TEMPERATURE GRID
C
      double precision x(ixzmax)
      double precision z(ixzmax)
      common /grid/xmax,zmax,dx,dz
C
C VELOCITY ON SAME GRID AS TEMPERATURES
C
      double precision vx(ixzmax,ixzmax)
      double precision vz(ixzmax,ixzmax)
      common /velo/vx,vz,uhalf,vmax
C
C VELOCITY GRID STAGGERED WITH RESPECT TO TEMPERATURE GRID
C
      double precision vxs(ixzmax,ixzmax)
      double precision vzs(ixzmax,ixzmax)
      common /velstg/vxs,vzs
C
C OUTPUT 1D GEOTHERM
C
      common /iprf/iprfv1,iprfv2,iprfh
C
C INPUT AND OUTPUT FILES
C
      character*30 tfilei,tfileo
C
C TIME
C
      common /timing/time,it
c      common /calc/itype,timfin,timstp,tend,icalc
      common /calc/tend,icalc
      double precision time(ixzmax)
C
C COMMUNICATION WITHIN PROGRAM
C
      common /grido/pres,x,z,nx,nz
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C
C TEMPERATURE FIELD
C
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      common /temo/tem,tems
      common /tem/ptem,tfilei,tfileo,item
C
C ARE THESE USED ANY MORE HERE?
C
      common /melt/ds,xh2o,cpxm,iparam,imelt,itemp
      common /vel/ufull,alpha,viss,crust,ialpha
C
C COMMON FOR MELT FRACTION FOR DF/Dt
C
      common /melfrc/f1x,f2x,f1z,f2z,fx
C
C COMMUNICATION OF DERIVATIVES THROUGHOUT PROGRAM
C
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
C
C TOLERANCE FOR CHECKING CONVERGENCE OF TEMPERATURE FIELD
C  TOLMX & TOLME, MAXIMUM AND MAXIMUM MEAN TEMPERATURE DIFFERENCE
C  BETWEEN TEMPERATURE FIELD AFTER VERTICAL SWEEPS AND THAT AFTER
C  HORIZONTAL SWEEPS.
C
      parameter (tolmx=0.1, tolme=0.1)
C
C MAXIMUM NUMBER OF ITERATIONS FOR TEMPERATURE CONVERGENCE
C
      parameter (maxit=10000)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C END OF SETUP
C START OF CALCULATIONS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C CONSTANT FOR CALCULATING ADIABATIC GRADIENT
C
      hcon=g*texps/cp
C
C TEMPERATURE ARRAY LABELS FOR 
C NO IS THE TEMPERATURE AT THE BEGINNING OF THIS TIMESTEP
C    REMAINS THE SAME THROUGHOUT THIS ROUTINE
C N1,N2 ARE USED IN THE CONVERGENCE LOOP IN THIS ROUTINE
C    OUTPUT TEMPERATURE PROFILE IS SET TO N1 AT THE END
C
      n2=6-n0-n1
C
C COUNTERS FOR VERTICAL AND HORIZONTAL PASSES
C
      ivpass=0
      ihpass=0
C
C LOOP HERE FOR ITERATION WITHIN THIS PASS TO ACCOUNT FOR
C TEMPERATURE DEPENDENT TERMS (ADIABATIC GRADIENT, DENSITY, CP)
C
 101  continue
C
C (^) MAKE UPWARDS VERTICAL PASSES (^)
C
      ivpass=ivpass+1
      do ix=1,nx
c      do ix=1,1
C
C BASAL BOUNDARY CONDITION
C TEMPERATURE FIXED AT MANTLE POTENTIAL TEMPERATURE
C
        b(1)=1.
        c(1)=0.
        r(1)=tem(n0,ix,nz)
C
C VERTICAL LOOP TO LOAD ARRAYS FOR TRIDIAGONAL MATRIX ALGORITHM.
C COUNTER IZUP COUNTS UP FROM BOTTOM TO TOP OF CALCULATION BOX
C THIS IS THE DIRECTION USED IN TRIDAG SINCE IT IS THE DIRECTION
C OF ADVECTION OF HEAT.
C COUNTER IZDOWN COUNTS UP FROM TOP TO BOTTOM OF CALCULATION BOX
C THIS IS THE DIRECTION USED FOR THE TEMPERATURE AND VELOCITY ARRAYS
C
        do izdown=2,nz-1
          izup=nz-izdown+1          
C
C CALCULATE NEIGHBOUR COEFFICIENTS
C
          call FVNB2D(anbe,anbw,anbn,anbs,ix,nx,izdown,nz,n1)
C
C NEIGHBOUR POINT COEFFICIENTS FOR HORIZONTALLY ADJACENT BOUNDARY POINTS 
C AXIAL BOUNDARY: TEMPERATURE FIELD IS SYMMETRICAL, SO ANBW=ANBE
C OFF-AXIS BOUNDARY: ASSUME VERY LARGE PECLET NUMBER, SO ANBE=0.
C
          if (ix.eq.1) anbw=anbe
          if (ix.eq.nx) anbe=0.
C
C FOR UPWARD SWEEP, 
C  SOUTHERN POINT (COEFFICIENT A_S) -> TRIDAG A()
C  MIDDLE POINT (COEFFICIENT A_P) -> TRIDAG B()
C  NORTHERN POINT (COEFFICIENT A_N) -> TRIDAG C().
C NOTE THESE COEFFICIENTS BECOME -VE IN COMPARISON WITH 
C PATANKAR WHEN REARRANGED TO BE INPUT TO THE NUMERICAL 
C RECIPES ROUTINE TRIDAG
C
          a(izup)=-anbs
          c(izup)=-anbn
C
C ADIABATIC TEMPERATURE GRADIENT GOES INTO THE SOURCE TERM S
C IN PATANKAR. WHEN LINEARIZING THE SOURCE TERM, THE CONSTANT
C PART SC MUST BE +VE AND THE LINEAR PART SP MUST BE -VE.  
C BUT THE ADIABATIC TEMPERATURE DECREASES DURING UPWELLING, IE -VE.
C THEREFORE, HANDLE THE TERM AS DESCRIBED IN PATANKAR 7.2.2.
C V IS THE VELOCITY AT THE TEMPERATURE NODE (REMEBER THE VELOCITY
C GRID IS STAGGERED)
C
          sc=0.
          v=0.5*(vzs(ix,izdown-1)+vzs(ix,izdown))
          h=hcon*(tem(n1,ix,izdown)+273.)
          if (tem(n1,ix,izdown).eq.0) then
            sp=0.
          elseif (tem(n1,ix,izdown).gt.0) then
            sp=(-rhos)*v*h/(tem(n1,ix,izdown)+273.)
          else
            write(6,*)'Source term S_P is +ve'
            stop 
          endif
c        sp=0.
C
C STEADY-STATE OR TIME-DEPENDENT CALCULATION (EQ 5.57E)
C
          if (icalc.eq.0) then
            a0p=0.
          else
            a0p=rhos*dx*dz/(time(it)-time(it-1))
	  endif
C
C CONSTANT TERM B IN PATANKAR (EQ 5.57F)
C ALSO INCLUDES WEST AND EAST NEIGHBOUR POINTS, 
C WHICH ARE CONSTANT IN THE UPWARD SWEEP
C FORMS RIGHT HAND SIDE ARRAY R() OF TRIDAG
C
          if (ix.gt.1) then
            tanbw=tem(n1,ix-1,izdown)*anbw
          else
            tanbw=tem(n1,ix+1,izdown)*anbe
          endif
          if (ix.lt.nx) then
            tanbe=tem(n1,ix+1,izdown)*anbe
          else
            tanbe=0.
          endif
          r(izup)=sc*dx*dz+a0p*tem(n0,ix,izdown)+tanbw+tanbe
C
C CENTRAL POINT COEFFICIENTS (A_P, EQ 5.57G) GO INTO 
C ARRAY B() FOR TRIDAG
C
          b(izup)=-a(izup)-c(izup)+anbe+anbw+a0p-sp*dx*dz
C
C END OF VERTICAL LOADING LOOP
C
        enddo
C
C TOP BOUNDARY CONDITION
C TEMPERATURE FIXED AT EARTH SURFACE TEMPERATURE
C
        a(nz)=0.
        b(nz)=1.
        r(nz)=tem(n0,ix,1)
C
C FIND VERTICAL TEMPERATURE PROFILE USING TRIDIAGONAL MATRIX ALGORITHM
C        
        call TRIDAG(a,b,c,r,u,nz)
C
C LOAD RESULTS INTO TEMPERATURE ARRAY
C
        do izdown=1,nz
          izup=nz-izdown+1          
            tem(n2,ix,izdown)=u(izup)
c            if (ix.eq.1) then
c              test=z(izdown)/z(nz)*tem(n0,ix,nz)
c              write(6,*)x(ix)/1000.,z(izdown)/1000.,
c     &         tem(n0,ix,izdown),tem(n2,ix,izdown),test
c          endif
        enddo
C
C END OF VERTICAL PASSES
C
      enddo
C
C CHECK FOR CONVERGENCE IN THE TEMPERATURE FIELD
C MAXIMUM AND RMS DIFFERENCES IN TEMPERATURE BETWEEN TWO PASSES
C
      call TCNVGE(tem,difmxv,difmev,nx,nz,n1,n2)
c      write(6,*)'Max and Mean Difference',difmxv,difmev,
c     &    ' V and H loops:',ivpass,ihpass
C
C MAKE FINAL TEMPERATURE FIELD AFTER VERTICAL PASSES
C STARTING TEMPERATURE FIELD FOR HORIZONTAL PASSES
C
      call SWAP(n1,n2)
c      if (ivpass.eq.1) goto 101
c      goto 104
C
C (>) MAKE HORIZONTAL PASSES IN DIRECTION OF PLATE SPREADING (>)
C
 102  continue
      ihpass=ihpass+1
      do izdown=2,nz-1
c      do izdown=nz-1,nz-1
        izup=nz-izdown+1          
C
C HORIZONTAL LOOP TO LOAD ARRAYS FOR TRIDIAGONAL MATRIX ALGORITHM
C
        do ix=1,nx
C
C CALCULATE NEIGHBOUR COEFFICIENTS
C
          call FVNB2D(anbe,anbw,anbn,anbs,ix,nx,izdown,nz,n1)
C
C OFF-AXIS BOUNDARY IS AN OUTFLOW BOUNDARY: NO BOUNDARY 
C INFORMATION REQUIRED IF WE ASSUME THAT THE PECLET NUMBER IS 
C LARGE, SINCE THE VALUE OF THE EASTERN COEFFICIENT AE IS THEN ZERO.
C THIS PRACTICE IS EQUIVALENT TO THAT OF BOWN & WHITE (1994) WHO
C ASSUMED D2T/DX2=0 WHICH MAKES CONDUCTION ZERO AT THE BOUNDARY
C IMPLYING A LARGE PECLET NUMBER.
C NOTE THAT THE PECLET NUMBER IS NOT ACTUALLY VERY LARGE (>10)
C UNLESS THE FULL PLATE SPREADING RATE IS OVER ABOUT 50 KM/MYR.
C NO NEED TO DO ANYTHING TO IMPLEMENT THIS CONDITION SINCE
C THE EASTERN NEIGHBOUR COEFFICIENT AT THE EASTERN BOUNDARY
C IS NOT DEFINED (EQUIVALENT TO ZERO) IN SUBROUTINE TRIDAG.
C
C AXIAL BOUNDARY CONDITION: DT/DX=0, NO HEAT FLUX ACROSS THE BOUNDARY.
C SINCE TEMPERATURE FIELD IS SYMMETRICAL ABOUT THE BOUNDARY, THE
C WESTERN NEIGHBOUR COEFFICIENT (ANBW) MUST BE THE SAME AS THE
C EASTERN COEFFICIENT (ANBE).  THEREFORE MULTIPLY THE EASTERN
C COEFFICIENT BY 2 TO ACCOUNT FOR BOTH WESTERN & EASTERN BOUNDARIES.
C THEN THERE IS NO NEED TO CHANGE THE DX OR DZ TERMS SINCE WE
C ARE WORKING WITH A NORMAL CONTROL VOLUME, NOT A HALF VOLUME.
C
          if (ix.eq.1) anbe=anbe*2.
C
C
C NEIGHBOUR POINT COEFFICIENTS, EQS 5.57D AND 5.57C. 
C FOR SWEEP AWAY FROM AXIS, 
C  WESTERN POINT (COEFFICIENT A_W) -> TRIDAG A()
C  EASTERN POINT (COEFFICIENT A_E) -> TRIDAG C().
C NOTE THESE COEFFICIENTS BECOME -VE IN COMPARISON WITH 
C PATANKAR WHEN REARRANGED TO BE INPUT TO THE NUMERICAL 
C RECIPES ROUTINE TRIDAG
C
          a(ix)=-anbw
          c(ix)=-anbe
C
C ADIABATIC GRADIENT ACCOUNTED FOR IN SOURCE TERM AS ABOVE
C
          sc=0.
          v=0.5*(vzs(ix,izdown-1)+vzs(ix,izdown))
          h=hcon*(tem(n1,ix,izdown)+273.)
          if (tem(n1,ix,izdown).eq.0) then
            sp=0.
          elseif (tem(n1,ix,izdown).gt.0) then
            sp=(-rhos)*v*h/(tem(n1,ix,izdown)+273.)
          else
            write(6,*)'Source term S_P is +ve'
            stop 
          endif
C
C STEADY-STATE OR TIME-DEPENDENT CALCULATION (EQ 5.57E)
C
          if (icalc.eq.0) then
            a0p=0.
          else
            a0p=rhos*dx*dz/(time(it)-time(it-1))
	        endif
C
C CONSTANT TERM B IN PATANKAR (EQ 5.57F)
C ALSO INCLUDES NORTH AND SOUTH NEIGHBOUR POINTS, 
C WHICH ARE CONSTANT IN THE HORIZONTAL SWEEP
C FORMS RIGHT HAND SIDE ARRAY R() OF TRIDAG
C
          r(ix)=sc*dx*dz+a0p*tem(n0,ix,izdown)
     &         +anbn*tem(n1,ix,izdown-1)+anbs*tem(n1,ix,izdown+1)
C
C CENTRAL POINT COEFFICIENTS (A_P, EQ 5.57G) GO INTO 
C ARRAY B() FOR TRIDAG
C
          b(ix)=-a(ix)-c(ix)+anbn+anbs+a0p-sp*dx*dz
C
C END OF HORIZONTAL LOADING LOOP
C
        enddo
C
C FIND NEW VERTICAL TEMPERATURE PROFILE USING TRIDIAGONAL
C MATRIX ALGORITHM
C        
        call TRIDAG(a,b,c,r,u,nx)
C
C LOAD RESULTS INTO TEMPERATURE ARRAY
C
        do ix=1,nx
          tem(n2,ix,izdown)=u(ix)
c          if (izdown.eq.5) then
c            test=(1.-x(ix)/x(nx))*tem(n0,1,izdown)
c            write(6,*)x(ix)/1000.,z(izdown)/1000.,
c     &              tem(n0,ix,izdown),tem(n2,ix,izdown),test
c          endif
        enddo
C
C END OF HORIZONTAL PASSES
C
      enddo
C
C CHECK FOR CONVERGENCE IN THE TEMPERATURE FIELD
C MAXIMUM AND RMS DIFFERENCES IN TEMPERATURE BETWEEN TWO PASSES
C
      call TCNVGE(tem,difmxh,difmeh,nx,nz,n1,n2)
c      write(6,*)'Max and Mean Difference',difmxh,difmeh,
c     &    ' V and H loops:',ihpass,ihpass
C
C MAKE FINAL TEMPERATURE FIELD AFTER HORIZONTAL PASSES
C STARTING TEMPERATURE FIELD FOR VERTICAL PASSES
C
      call SWAP(n1,n2)
C
C HAS TEMPERATURE FIELD CONVERGED?
C IF SO CONTINUE, BUT REPORT IF NOT WITHIN SPECIFIED TOLERANCE
C
 103  continue
      difmx=abs(difmxv-difmxh)
      difme=abs(difmev-difmeh)
      if (difmx.lt.1.e-10 .and. difme.lt.1.e-10) then
        write(6,*)'Temperature calculation in TEMFVI converged'
        write(6,*)'*** outside specified tolerance ***'
        if (difmxh.gt.tolmx) then
          write(6,*)'Max. temp. change:  actual,',difmxh,
     &              ';  tolerated,',tolmx
        endif
        if (difmeh.gt.tolme) then
          write(6,*)'Mean temp. change:  actual,',difmeh,
     &              ';  tolerated,',tolme
        endif
        write(6,*)'in',ihpass,' iterations'
        goto 104
      endif
C
C MAXIMUM NUMBER OF ITERATIONS REACHED WITHOUT CONVERGENCE
C
      if (ihpass.gt.maxit) then
        write(6,*)'Temperature calculation in TEMFVI:'
        write(6,*)'Max. iterations exceeded without convergence'
        write(6,*)'(',ihpass,' iterations )'
        write(6,*)'Max. temp. difference V & H:',difmxv,difmxh
        write(6,*)'Max. difference tolerated:',tolmx
        write(6,*)'Mean. temp. difference V & H:',difmev,difmeh
        write(6,*)'Mean. difference tolerated:',tolme
        write(6,*)'STOPPING'
        stop
      endif
C
C HAS TEMPERATURE FIELD CONVERGED WITHIN SPECIFIED TOLERANCES?
C IF NOT, LOOP BACK FOR ANOTHER HORIZONTAL AND VERTICAL PASS.
C
      if (difmxv.gt.tolmx .or. difmxh.gt.tolmx .or.
     &    difmev.gt.tolme .or. difmeh.gt.tolme) goto 101
C
C TEMPERATURE FIELD CONVERGED
C
 104  continue
C
C WRITE OUT 1D VERTICAL PROFILE
C  COL 1, DEPTH
C  COL 2, TEMPERATURE
C  COL 3, HORIZONTAL VELOCITY
C  COL 4, VERTICAL VELOCITY 
C  COL 5, PRESSURE (BAR)
C
      rhop=3000.0
 1    format ('> ',f8.2) 
      if (iprfv1.ne.0) then
        write(10,1) time(it)/1.e6
        do iz=1,nz
          geo1d=0.
          do ix=iprfv1,iprfv2
            geo1d=geo1d+tem(n1,ix,iz)
          enddo
          geo1d=geo1d/float(iprfv2-iprfv1+1)
c          velu=(vxz(ix,iz)+vxz(ix,iz-1))/2.
c          velv=
          pres1d=(-rhop)*g*z(iz)*1.e-5/sy/sy
          write(10,*)z(iz)/1000.,geo1d,pres1d
        enddo
      endif
C
C WRITE OUT HORIZONTAL PROFILE
C
      if (iprfh.ne.0) then
        write(11,1) time(it)/1.e6
        do ix=1,nx
c          velu=(vxz(ix,iz)+vxz(ix,iz-1))/2.
c          velv=
          pres1d=(-rhop)*g*z(iz)*1.e-5/sy/sy
          write(11,*)x(ix)/1000.,tem(n1,ix,iprfh),pres1d
        enddo
      endif
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine FVNB2D(anbe,anbw,anbn,anbs,ix,nx,iz,nz,n1)
C
C CALCULATE NEIGHBOUR POINT WEIGHTING COEFFICIENTS
C FOR THE FINITE VOLUME DISCRETIZATION SCHEME OF
C PATANKAR (1980).
C THE NEIGHBOURS SURROUND THE POINT (IX,IZ).
C THE DIMENSIONS OF THE FIELD (NX,NZ) ARE USED TO CHECK WHETHER
C THE POINT IS ON A BOUNDARY; COEFFICIENTS ARE NOT CALCULATED
C FOR THE REGION OUTSIDE THE BOUNDARY. 
C OUTPUT: NEIGHBOUR COEFFICIENTS ANBE, ANBW, ANBN, ANBS.
C
      implicit double precision(a-h,o-z)
C
C MAXIMUM SIZE OF GRID
C
      parameter (ixzmax=2001)
C
C TEMPERATURE GRID
C
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      common /temo/tem,tems
C
C VELOCITY GRID STAGGERED WITH RESPECT TO TEMPERATURE GRID
C (INPUT, NOT ALTERED HERE)
C
      double precision vxs(ixzmax,ixzmax)
      double precision vzs(ixzmax,ixzmax)
      common /velstg/vxs,vzs
C
C FIXED PHYSICAL PARAMETERS AND GRID SPACING
C
      common /grid/xmax,zmax,dx,dz
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C
C THE NEIGHBOUR COEFFICIENTS ARE CALCULATED AS FOLLOWS.
C 1. G?, RATIO OF THERMAL CONDUCTIVIY AND SPECIFIC HEAT CAPACITY
C    (THIS IS THE DIFFUSION PARAMETER GAMMA IN PATANKAR 1980)
C 2. D?, CONDUCTANCES FROM EQ 5.58, ASSUMING THAT THE FINITE 
C    VOLUME BOUNDARIES ARE MIDWAY BETWEEN THE NODES SO 
C    DX=DELTAX AND DZ=DELTAZ
C 3. R?, DENSITIES AT FINITE VOLUME BOUNDARIES
C 4. F?, FLOW RATES FROM EQ 5.52
C    NO NEED TO INTERPOLATE BECAUSE OF STAGGERED VELOCITY GRID
C 5. P?, PECLET NUMBERS (ADVECTIVE/DIFFUSIVE TRANSPORT) FROM EQ 5.59
C 6. A?, POWER-LAW DIFFERENCING SCHEME FUNCTION A(|P|), FROM EQ 5.60
C 7. ANB?, NEIGHBOUR POINT COEFFICIENTS, EQS 5.57A-D. 
C
C TEMPERATURE
C
      temp=tem(n1,ix,iz)
C
C EASTERN NEIGHBOUR COEFFICIENT
C
      if (ix.lt.nx) then          
        ge=TTCON(temp)/cp
        de=ge*dz/dx
        re=rhos
        fe=re*vxs(ix,iz)*dz
c        fe=re*0.01*dz
        pe=fe/de
        ae=APOW(pe)
        anbe=de*ae+max(-fe,0.)
      else
        anbe=0.
      endif
C
C WESTERN NEIGHBOUR COEFFICIENT
C
      if (ix.gt.1) then
        gw=TTCON(temp)/cp
        dw=gw*dz/dx
        rw=rhos
        fw=rw*vxs(ix-1,iz)*dz
c        fw=rw*0.01*dz
        pw=fw/dw
        aw=APOW(pw)
        anbw=dw*aw+max(fw,0.)
      else
        anbw=0.
      endif 
C
C NORTHERN NEIHBOUR COEFFICIENT
C
      if (iz.gt.1) then
        gn=TTCON(temp)/cp
        dn=gn*dx/dz
        rn=rhos
        fn=rn*vzs(ix,iz-1)*dx
        pn=fn/dn
        an=APOW(pn)
        anbn=dn*an+max(-fn,0.)
      else
        anbn=0.
      endif
C
C SOUTHERN NEIGHBOUR COEFFICIENT
C
      if (iz.lt.nz) then
        gs=TTCON(temp)/cp
        ds=gs*dx/dz
        rs=rhos
        fs=rs*vzs(ix,iz)*dx
        ps=fs/ds
        as=APOW(ps)
        anbs=ds*as+max(fs,0.)
      else
        anbs=0.
      endif
C
      return
      end
C
C----------------------------------------------------------
C
      function TTCON(tem)
C
C RETURNS THERMAL CONDUCTIVITY (CAN BE A FUNCTION OF TEMPERATURE)
C  ITCON=1, CONSTANT CALCULATED FROM DIFFUSIVITY AND SPECIFIC HEAT
C  ITCON=2, FUNCTION OF TEMPERATURE, (NOT YET IMPLEMENTED)
C AT PRESENT ITCON IS HARD-WIRED HERE
C
      implicit double precision(a-h,o-z)
      parameter (itcon=1)
C
C IMPORT CONSTANT VALUE OF TCON CALCULATED EARLIER
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C
C CONSTANTS FOR PARAMETERIZATION IN MCKENZIE ET AL., EPSL, 2005
C
      data b /5.3/
      data c /0.0015/
      double precision d(4)
      data d /1.753e-2,-1.0365e-4,2.2451e-7,-3.4071e-11/
C
C CONSTANT VALUE OF TCON
C
      if (itcon.eq.1) then
        ttcon=tcon
C
C TCON DEPENDS ON TEMPERATURE (MCKENZIE ET AL., 2005, EQ 4)
C
      elseif (itcon.eq.2) then
        ttcon=b/(1.+c*tem)
        do i=1,4
          ttcon=ttcon+d(i)*(tem+273.)**float(i-1)
        enddo
        ttcon=ttcon*sy*sy*sy
C
C INVALID OPTION
C
      else
        write(6,*)'Invalid option for Thermal Conductivity, ITCON:',
     &            itcon
        stop
      endif
      return
      end
C
C----------------------------------------------------------
C
      function APOW(p)
C
C CALCULATES WEIGHTING COEFFICIENT A FROM INPUT PECLET
C NUMBER P ACCORDING TO POWER LAW SCHEME OF PATANKAR (1980)
C SECTION 5.2.6 AND TABLE 5.2
C CHOICE OF SHCEME IS HARDWIRED AS PARAMETER 1A
C  IA=2, UPWIND DIFFERENCING
C  IA=5, POWER LAW DIFFERENCING
C
      implicit double precision(a-h,o-z)
      parameter (ia=5)
      if (ia.eq.2) then
        apow=1.
      else if (ia.eq.5) then
        p1=1.-0.1*abs(p)
        p2=p1*p1
        p5=p2*p2*p1
        apow=max(0.,p5)
      endif
      return
      end
C
C----------------------------------------------------------
C
      subroutine TEMPDE
      implicit double precision(a-h,o-z)
C
C (1) CALCULATES TEMPERATURE STRUCTURE FROM ADVECTION-DIFFUSION
C PARTIAL DIFFERENTIAL EQUATION:
C DT/Dt = K (D2T/DX2 + D2T/DZ2) - VX DT/DX - VZ (DT/DZ + H).
C ALTERNATING-DIRECTION IMPLICIT METHOD FOR DIFFUSIVE PORTION.
C UPWIND DIFFERENCING FOR ADVECTIVE PORTION.
C CALCUULATION PROCEEDS UNTIL STEADY STATE TEMPERATURE HAS
C BEEN REACHED OR UNTIL TIME=TOUT, WHICHEVER IS SOONER.
C
C (2) CALLS ROUTINE FOR MELTING CALUCLATIONS.
C IF ICALC=0, MELTING CALCULATIONS CARRIED OUT ONLY AT END
C OF TEMPERATURE CALCULATIONS.
C IF ICALC=1, MELTING CALCULATIONS CARRIED OUT AT EACH
C TIMESTEP.
C
      parameter (ixzmax=2001)
      double precision pres(ixzmax)
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      double precision vx(ixzmax,ixzmax)
      double precision vz(ixzmax,ixzmax)
      double precision x(ixzmax)
      double precision z(ixzmax)
      double precision dpda(ixzmax,ixzmax)
      double precision dtrda(ixzmax,ixzmax)
      double precision dfda(ixzmax,ixzmax)
      character*30 tfilei,tfileo
C
C COMMUNICATION WITHIN PROGRAM
C
c      common /calc/itype,timfin,timstp,tend,icalc
      common /calc/tend,icalc
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
      common /tem/ptem,tfilei,tfileo,item
      common /temo/tem,tems
      common /vel/ufull,alpha,viss,crust,ialpha
      common /velo/vx,vz,uhalf,vmax
      common /melt/ds,xh2o,cpxm,iparam,imelt,itemp
C COMMON FOR MELT FRACTION FOR DF/Dt
      common /melfrc/f1x,f2x,f1z,f2z,fx
C COMMUNICATION OF DERIVATIVES THROUGHOUT PROGRAM
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
C
C DETERMINE TIMESCALE
C IMPLICIT METHOD FOR DIFFUSIVE PORTION UNCONDITIONALLY STABLE
C USE COURANT CRITERION FOR ADVECTIVE PORTION
C EACH TIMESTEP OF PERIOD DT IS DIVIDED INTO 4 EQUAL SUBSTEPS 
C OF PERIOD DT4, ONE FOR EACH PORTION OF THE EQUATION
C
      smaxx=abs(uhalf)
      smaxz=abs(vmax)
      smax=max(smaxx,smaxz)
      dtc=dx/sqrt(2.)/smax
      write(6,*)'Maximum speed (km/Myr):',smax*1.e3
      write(6,*)'Courant timestep (Myr):',dtc*1.e-6
      dtd=dx*dx/2./tdiff  /5.
      write(6,*)'Diffusion timestep (Myr):',dtd*1.e-6
      dt=min(dtc,dtd)
      write(6,*)'Timestep used (Myr):',dt*1.e-6
C
C Constants For The Diffusive And Advective Terms
C And Adiabatic Temperature Gradient
C
      Cdiff=Tdiff*Dt/Dx/Dx
      Cadv=Dt/Dx
C      Ch=-G*Texps/Cp
      Ch=G*Texps/Cp
C
C OVER-RELAXATION FACTOR USEFUL IF MELTING NOT TO BE CALCULATED
C AT EACH TIMESTEP
C
      if (icalc.eq.0) then
         for=1.5
      else
         for=1.
      endif
C
C SET TIME COUNTERS TO ZERO
C
      t=0.
      it=0
C
C FINITE DIFFERENCING TIME LOOP BEGINS HERE
C
      n0=1
      n1=2
C 
 1    continue
C CHANGED T FOR TIM TO TRY AND DISTINGUISH BETWEEN TIME AND TEMP
      if (icalc.ne.0 .and. (t+dt).gt.tend) dt=tend-t
      t=t+dt
      it=it+1
      dtemmx=0.
C
C SPATIAL LOOPS
C
      do ix=1,nx
         do iz=2,nz-1
C
C      do ix=1,1
C         do iz=10,10
C
C HORIZONTAL AND VERTICAL DIFFUSION TERMS
C
            if (ix.eq.1 .or. ix.eq.nx) then
               xd=0.
            else
               xd=(tem(n0,ix+1,iz)-2.*tem(n0,ix,iz)+tem(n0,ix-1,iz))
            endif
            zd=(tem(n0,ix,iz+1)-2.*tem(n0,ix,iz)+tem(n0,ix,iz-1))
C
C UPWIND DIFFERENCING FOR HORIZONTAL ADVECTION TERM
C
            if (ix.eq.1) then
               xa=0.
            else
               xa=tem(n0,ix,iz)-tem(n0,ix-1,iz)
            endif
            za=tem(n0,ix,iz)-tem(n0,ix,iz+1)
C
C ADIABATIC TEMPERATURE GRADIENT
C
            h=ch*(tem(n0,ix,iz)+273.)
C
C TEMPERATURE CHANGE (STORE MAXIMUM)
C
C            write(6,*)ix,iz,cdiff,cadv,tem(n0,ix,iz),
C     &       tem(n0,ix,iz+1),tem(n0,ix+1,iz),tem(n0,ix,iz-1),
C     &       tem(n0,ix-1,iz),xd,zd,xa,za,h
            dtem=cdiff*(xd+zd)-cadv*(vx(ix,iz)*xa+vz(ix,iz)*(za+h))
C            dtem=cdiff*(xd+zd)-cadv*(vx(ix,iz)*xa+vz(ix,iz)*za)
C     &           -dt*vz(ix,iz)*h
C            write(6,*)'cdiff',cdiff,'xd',xd,'zd',zd,'cadv',
C     &                cadv,'vx',vx(ix,iz),'vz',vz(ix,iz),
C     &                'xa',xa,'za',za,'h',h,'dtem',dtem
C WHEN THE MAXIMUM TEMPERATURE DIFFERENCE OVER THE WHOLE GRID IS LESS
C THAN 1, THEN ALL THE DTEMS ARE LESS THAN 1.
            if (abs(dtem).gt.dtemmx) dtemmx=abs(dtem)
C
C UPDATE TEMPERATURE
C
C            tem(n1,ix,iz)=tem(n0,ix,iz)+for*dtem
            tem(n1,ix,iz)=tem(n0,ix,iz)+dtem
         enddo
      enddo
C
C IF AIMING FOR STEADY STATE TEMPERATURE STRUCTURE
C
      write(6,*)it,' Time (Myr):',t*1.e-6,',    Max dTem:',dtemmx
      if (icalc.eq.0) then
         if (dtemmx.gt.0.1) then
            call SWAP(n0,n1)
            goto 1
         else
          write(6,*)'Steady-state temperature structure achieved after',
     &           t*1.e-6,' Myr'
            write(6,*)'Calculating melt volume...'
            call MEL(n1,re1)
         endif
C
C IF AIMING FOR PARTICULAR TIME
C
      else if (icalc.eq.1) then
         call MEL(n1,re1)
         call SWAP(n0,n1)
         if (t.lt.tend) goto 1
      endif
C
C WRITE OUT STARTING AND FINAL TEMPERATURE STRCTUCTURES
C
      open(7,file='bw94_tem_out')
 2    format (4(f10.2,1x))
      do ix=1,nx
         do iz=1,nz
            write(7,2)x(ix)*1.e-3,z(iz)*1.e-3,tems(ix,iz),tem(n1,ix,iz)
         enddo
      enddo
      close (7)
      write(6,*)'Written temperature field to "bw94_tem_out" (km, oC)'
C
C
C WRITE OUT 1D GEOTHERM AT RIDGE
C
      open(8,file='bw94_1d_geotherm_out')
      open(9,file='bw94_adiabat_ptpath_out')
C
C DRIDGE IS THE WIDTH OF THE ZONE OF THE AVERAGING
C XOFFS IS THE START OF THE AVERAGED REGION AND XOFFF
C IS THE FINISH (IN M), NRIDGS IS THE STARTING NODE
C AND NRIDGF IS THE FINISHING NODE, GEO1D IS THE GEOTHERM
C
      rhop=3000.0
      xoffs=0.0
      xofff=5000.0
      nridgs=xoffs/dx+1
      nridgf=xofff/dx+1
      do iz=1,nz 
        geo1d=0. 
        do ix=nridgs,nridgf
          geo1d=geo1d+tem(n1,ix,iz) 
        enddo
        geo1d=geo1d/(nridgf-nridgs+1)
        z(iz)=(iz-1)*(-dz)
        write(8,*)z(iz)/1000.,geo1d
      enddo 
      close (8)
C
C CONVERT PRESSURE AND OUTPUT FROM THE BASE OF THE CRUST
C NCRUST IS THE NODE AT THE BASE OF THE CRUST, PRES1D IS
C THE PRESSURE (IN BARS, AS REQUIRED FOR THE ADIABAT PROGRAM)
C  
      ncrust=crust/dz+1
      do iz=ncrust,nz 
        geo1d=0. 
        do ix=nridgs,nridgf
           geo1d=geo1d+tem(n1,ix,iz) 
        enddo
        geo1d=geo1d/(nridgf-nridgs+1)
        z(iz)=(iz-1)*(-dz)
        pres1d=(-rhop)*g*z(iz)*1.e-5/sy/sy
        write(9,*)pres1d,geo1d
      enddo    
C      
      close (9)
C
C
      return
      end
C
C----------------------------------------------------------
C
      subroutine SWAP(n0,n1)
      implicit double precision(a-h,o-z)
      ntmp=n1
      n1=n0
      n0=ntmp
      return
      end
C
C----------------------------------------------------------
C
      subroutine TCNVGE(tem,difmax,difmea,nx,nz,n0,n1)
C
C TEST FOR CONVERGENCE BETWEEN STARTING AND FINISHING 
C TEMPERATURE FIELDS.  RETURNS:
C   MAXIMUM TEMPERATURE DIFFERENCE, DIFMAX
C   MEAN TEMPERATURE DIFFERENCE, DIFMEA
C
      implicit double precision(a-h,o-z)
C
C MAXIMUM SIZE OF GRID
C
      parameter (ixzmax=2001)
      double precision tem(3,ixzmax,ixzmax)
C
C COMPARE TEMPERATURE GRIDS
C
      difmax=0.
      difmea=0.
      do ix=1,nx
        do iz=1,nz
          dif=abs(tem(n1,ix,iz)-tem(n0,ix,iz))
          if (dif.gt.difmax) difmax=dif
          difmea=difmea+dif
        enddo
      enddo
      difmea=difmea/float(nx*nz)
C
      return
      end
C
C----------------------------------------------------------
C
C SOLVE SIMULTANEOUS EQUATIONS IN TRIDIAGONAL MATRIX
C ROUTINE FROM NUMERICALL RECIPES
C
      SUBROUTINE tridag(a,b,c,r,u,n)
      implicit double precision(a-h,o-z)
      INTEGER n,NMAX
      double precision a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=2000)
      INTEGER j
      double precision bet,gam(NMAX)
      if(b(1).eq.0.) then
      write(6,*) 'tridag: rewrite equations'
      stop
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
        write(6,*) 'tridag failed'
        stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
C
C------------------------------------------------------------------
C

