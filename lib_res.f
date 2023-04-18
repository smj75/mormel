C 
C SUBROUTINES FOR TRACKING RESIDUE AROUND THE GRID
C RACHEL WALTERS 2007-9
C MODIFIED STEPHEN JONES SPRING 2011
C 
C CONTAINS FOLLOWING ROUTINES
C	RESFVI
C	FVRS2D
C	RCNVGE
C	SWAP2
C
      subroutine ADVFVI(array,tolmx,tolme,m0,m1)
C
C IMPLICIT FINITE VOLUME STEP USING ADVECTION ONLY
C PARTIAL DIFFERENTIAL EQUATION:
C   dR         dR       dR
C   --  =  -Vx --  - Vz --
C   dt         dx       dz
C
C INPUTS:
C   ARRAY CONTAINS FIELD TO BE TRANSPORTED AROUND GRID
C   TOLMX
C   TOLME ARE THE MAXIMUM AND MEAN VARIABILITIES ALLOWED 
C         BEFORE SOLUTION IS ASSUMED TO HAVE CONVERGED
C   N0 INDEX OF STARTING FIELD
C   N1 INDEX OF FINISHING FIELD
C
      implicit double precision(a-h,o-z)
C
C MAX NUMBER OF GRID NODES
C
      parameter (ixzmax=2001)
C
C ARRAY CONTAINING FIELD TO BE TRANSPORTED ROUND GRID
C
      double precision array(3,ixzmax,ixzmax)
C
C GRID
C
      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
C
C ARRAY FOR TRIDIAGONAL MATRIX
C
      double precision a(ixzmax)
      double precision b(ixzmax)
      double precision c(ixzmax)
      double precision r(ixzmax)
      double precision u(ixzmax)
C
C MELTING
C
      common /melt/ds,xh2o,cpxm,iparam,imelt,itemp
C
C VELOCITY GRID ON ORIGINAL NODES
C
      double precision vx(ixzmax,ixzmax),vz(ixzmax)
      common /velo/vx,vz,uhalf,vmax
      common /vel/ufull,alpha,viss,crust,ialpha
C
C VELOCITY ON STAGGERED GRID
C
      double precision vxs(ixzmax,ixzmax),vzs(ixzmax,ixzmax)
      common /velstg/vxs,vzs
C
C PARAMETERS
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &        texpf,texps
C
C TIMING - FOR CALCULATING AOP
C
      double precision time(ixzmax)
      common /timing/time,it
      common /calc/tend,icalc
C
C MAXIMUM NUMBER OF ITERATIONS FOR CONVERGENCE
C
      parameter (maxit=10000)
C
C TINY VALUE USED TO ASSESS WHETHER CALCULATION HAS CONVERGED
C
      parameter (tiny=1.0d-20)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C END OF SET UP
C START OF CALCULATIONS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C IF STEADY-STATE CALCULATION THEN NO NEED TO CALCULATE
C THE RESIDUE 
C
      if (icalc.eq.0) return
C
C ARRAY COUNTERS 
C REO IS THE RESIDUE VALUES AT THE BEGINNING
C RE1,RE2 ARE USED IN CONVERGENCE LOOP
C OUTPUT RESIDUE IS SET TO R1 AT THE END
C
      m2=6-m0-m1 
C 
C COUNTERS FOR VERTICAL AND HORIZONTAL PASS
C
      ivpass=0
      ihpass=0
C
C LOOP HERE FOR ITERATION WITHIN THIS PASS TO ALLOW CONVERGENCE
C
 101  continue
C
C **MAKE UPWRADS VERTICAL PASSES**
C
      ivpass=ivpass+1
      do ix=1,nx
C       do ix=1,1
C
C BASAL BOUNDARY CONDITION FIXED AT RESIDUE = 1, THE PROGRAMME 
C DOES NOT ALLOW ANY MELTING TO OCCUR AT THE BASAL NODES
C
        b(1)=1.
        c(1)=0.
        r(1)=array(m0,ix,nz)
C        write(6,*)'BOTTOM:',r(1)
C
C VERTICAL LOOP TO LOAD ARRAYS FOR TRIDIAGONAL MATRIX ALGORITHM.
C COUNTER IZUP COUNTS UP FROM BOTTOM TO TOP OF CALCULATION BOX
C THIS IS THE DIRECTION USED IN TRIDAG SINCE IT IS IN THE DIRECTION
C OF ADVECTION OF THE RESIDUE
C COUNTER IZDOWN COUNTS UP FROM TOP TO BOTTOM OF CALCULATION BOX
C THIS IS THE DIRECTION USED FOR THE RESIDUE AND VELOCITY ARRAY
C
        do izdown=2,nz-1
          izup=nz-izdown+1
C
C CALCULATE NEIGHBOUR COEFFICIENTS
C
          call FVRS2D(anbe,anbw,anbn,anbs,ix,nx,izdown,nz)
C
C NEIGHBOUR COEFFICIENTS FOR HORIZONTALLY ADJACENT BOUNDARY POINTS
C AXIAL BOUNDARY: RESIDUE FIELD IS ASSUMED SYMMETRICAL AS TEMP
C FIELD IS SYMMETRICAL
C OFF-AXIAL BOUNDARY: ASSUME VERY LARGE PECLET NUMBER, SO ANBE=0
C
          if (ix.eq.1) anbw=anbe
          if (ix.eq.nx) anbe=0.
C
C LOADING TRIDIAGONAL MATRIX FOR UPWARD SWEEP
C A_S = TRIDAG A()
C A_N = TRIDAG C()
C A_P = TRIDAG P()
C THESE COEFFICIENTS BECOME NEGATIVE WHEN REARRANGED TO
C BE INPUT INTO THE NUMERICAL RECIPES TRIDAG ROUTINE
C
          a(izup)=-anbs
          c(izup)=-anbn
C
C SOURCE TERM IS ZERO, SET SC AND SP EQUAL TO ZERO
C
          sc=0.
          sp=0.
          if(sp.gt.0.) then
            write(6,*)'Source term S_P is +ve'
          endif
C
C IF STATEMENT IN PREPARATION A POSSIBLE CHANGE IN
C THE SOURCE TERM AS THE PROGRAM IS DEVELOPED FURTHER
C
          a0p=rhos*dx*dz/(time(it)-time(it-1))
C
C CONSTANT TERM B (EQ. 5.57F) ALSO INCLUDES EAST AND WEST
C NEIGHBOUR COEFFICIENTS WHICH ARE CONSTANT IN THE UPWARD
C SWEEP. B IS R() IN TRIDAG, RIGHT HAND SIDE OF THE ARRAY
C 
          if (ix.gt.1) then
            ranbw=array(m1,ix-1,izdown)*anbw
          else
            ranbw=array(m1,ix+1,izdown)*anbe
          endif
          if (ix.lt.nx) then
            ranbe=array(m1,ix+1,izdown)*anbe
          else
C 0 OR 1??
            ranbe=0.
          endif
C
C B TERM THE TERMS RANBE AND RANBW ARE REARRANGED TO BE INCLUDED
C IN THE B TERM FROM THE MAIN DISCRETIZATION EQUATION AS THEY 
C REMAIN CONSTANT IN THE UPWARD SWEEP
C
          r(izup)=sc*dx*dz+a0p*array(m0,ix,izdown)+ranbw+ranbe
C
C CENTRAL POINT COEFFICIENT A_P INTO ARRAY B()
C
          b(izup)=-a(izup)-c(izup)+anbe+anbw+a0p-sp*dxdz
C
C END OF VERTICAL LOADING LOOP
C
        enddo
C
C TOP BOUNDARY CONDITION
C NO MELTING SHOULD OCCUR IN THE CRUST SO HOPEFULLY IF 
C FIX BOUNDARY AT 1
C
        a(nz)=0.
        b(nz)=1.
        r(nz)=array(m0,ix,1)
C
C FIND VERTICAL PROFILE USING TRIDIAGONAL MATRIX ALGORITHM
C
        call TRIDAG(a,b,c,r,u,nz)
C
C LOAD RESULTS INTO NEW RESIDUE ARRAY
C
        do izdown=1,nz
          izup=nz-izdown+1
          array(m2,ix,izdown)=u(izup)
        enddo
C
C END OF VERTICAL PASSES
C
      enddo
C 
C CHECK FOR CONVERGENCE
C MAXIMUM AND RMS DIFFERENCES IN FIELD BETWEEN TWO PASSES
C
      call RCNVGE(array,rdfmxv,rdfmev,nx,nz,m1,m2)
C
C PUT FINAL FIELD AFTER VERTICAL PASSES
C INTO THE STARTING FIELD FOR HORIZONTAL PASSES
C
      call SWAP(m1,m2)
C
C
C **MAKE HORIZONTAL PASSES IN DIRECTION OF PLATE SPREADING**
C
 102  continue
      ihpass=ihpass+1
      do izdown=2,nz-1
C       do izdown=nz-1,nz-1
         izup=nz-izdown+1
C
C HORIZONTAL LOOP TO LOAD ARRAYS FOR TRIDIAGONAL MATRIX
C
      do ix=1,nx
C
C CALCULATE NEIGHBOUR COEFFICIENTS
C
      call FVRS2D(anbe,anbw,anbn,anbs,ix,nx,izdown,nz)
C
C THE OFF AXIS BOUNDARY - PERHAPS THIS SHOULD BE MADE EQUAL
C TO THE END NODES IF ALL VELOCITIES ARE HORIZONTAL, OR
C PERHAPS IT DOESN'T MATTER BECAUSE EVERYTHING IS FLOWING OUT
C
C AXIAL BOUNDARY: SYMMETRICAL AS TEMPERATURE FIELD IS SYMMETRICAL
C
      if(ix.eq.1) anbe=anbe*2
C
C LOADING TRIDIAGONAL MATRIX
C WESTERN COEFFICIENT - A_W - A()
C EASTERN COEFFICIENT - A_E - C()
C COEFFICIENTS BECOME NEGATIVE BECAUSE OF NUMERICAL RECIPES
C
      a(ix)=-anbw
      c(ix)=-anbe
C
C SOURCE TERM IS ZERO, SET SC AND SP EQUAL TO ZERO
C
      sc=0.
      sp=0.
      if(sp.gt.0.) then
        write(6,*)'Source term S_P is +ve'
      endif
C
C IF STATEMENT IN PREPARATION A POSSIBLE CHANGE IN
C THE SOURCE TERM AS THE PROGRAM IS DEVELOPED FURTHER
C
      a0p=rhos*dx*dz/(time(it)-time(it-1))
C
C CONSTANT TERM B IN PATANKHAR - R() IN TRIDAG
C INCLUDES N AND S COEFFICIENTS
C
      r(ix)=sc*dx*dz+a0p*array(m0,ix,izdown)
     &        +anbn*array(m1,ix,izdown-1)+anbs*array(m1,ix,izdown+1)
C
C CENTRAL POINT COEFFICIENT - A_P - B()
C
      b(ix)=-a(ix)-c(ix)+anbn+anbs+a0p-sp*dx*dz
C
C END OF HORIZONTAL LOADING LOOP
C
      enddo
C 
C FIND NEW HORIZONTAL RESIDUE PROFILE
C
      call TRIDAG(a,b,c,r,u,nx)
C
C LOAD RESULTS INTO NEW ARRAY
C
      do ix=1,nx
         array(m2,ix,izdown)=u(ix)
      enddo
C
C END OF HORIZONTAL PASSES
C
      enddo
C
C CHECK FOR CONVERGENCE IN THE FIELD
C
      call RCNVGE(array,rdfmxh,rdfmeh,nx,nz,m1,m2)
C
C PUT FINAL FIELD AFTER HORIZONTAL PASSES
C INTO THE STARTING FIELD FOR VERTICAL PASSES
C
      call SWAP(m1,m2)
C
C HAS CALCULATION CONVERGED? 
C IF SO CONTINUES BUT REPORT IF NOT WITHIN SPECIFIED TOLERANCE
C
 103  continue
      rdfmx=abs(rdfmxv-rdfmxh)
      rdfme=abs(rdfmev-rdfmeh)
      if(rdfmx.lt.tiny .and. rdfme.lt.tiny) then
         write(6,*)'ADVFVI: Calculation has converged'
         write(6,*)'***outside specified tolerance***'
c         if(rdfmxh.gt.tolmx) then
            write(6,*)'Max change:',rdfmxh,
     &                '; tolerated, ',tolmx
c         endif
c         if(rdfmeh.gt.tolme) then
            write(6,*)'Mean change:',rdfmeh,
     &                '; tolerated, ',tolme
c         endif
         write(6,*)'in',ihpass,' iterations'
         goto 104
      endif
C
C MAXIMUM NUMBER OF ITERATIONS REACHED WITHOUT CONVERGENCE
C
      if (ihpass.gt.maxit) then
        write(6,*)'In subroutine ADVFVI:'
        write(6,*)'Max. iterations exceeded without convergence'
        write(6,*)'(',ihpass,' iterations )'
        write(6,*)'Max. difference V & H:',rdfmxv,rdfmxh
        write(6,*)'Max. difference tolerated:',tolmx
        write(6,*)'Mean difference V & H:',rdfmev,rdfmeh
        write(6,*)'Mean difference tolerated:',tolme
        write(6,*)'STOPPING'
        stop
      endif
C
C HAS RESIDUE FIELD CONVERGED WITHIN SPECIFIED TOLERANCES
C IF NOT LOOP BACK FOR ANOTHER VERTICAL AND HORIZONTAL PASS
C
      if (rdfmxv.gt.tolmx .or. rdfmxh.gt.tolmx .or.
     &   rdfmev.gt.tolme .or. rdfmeh.gt.tolme) goto 101
C
C RESIDUE FIELD CONVERGED
C  
 104  continue
C
C END SUBROUTINE
C
      return
      end
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C SUBROUTINE TO CALCULATE NEIGHBOUR COEFFICIENTS FOR RESIDUE
C
      subroutine FVRS2D(anbe,anbw,anbn,anbs,ix,nx,iz,nz)
C
C DIMENSIONS OF THE FIELD NX AND NZ ARE USED TO CHECK WHETHER THE
C POINT IS ON A BOUNDARY. COEFFICIENTS ARE NOT CALCULATED OUTSIDE THE
C BOUNDARY
C OUTPUT: NEIGHBOUR COEFFICIENTS ANBE, ANBW, ANBN, ANBS.
C
      implicit double precision(a-h,o-z)
C
C MAXIMUM SIZE OF GRID
C
      parameter (ixzmax=2001)
C
C RESIDUE GRID
C
c      double precision rest(ixzmax,ixzmax),res(3,ixzmax,ixzmax)
c      common /resdue/rest,res,resc
C
C VELOCITY GRID STAGGERED WITH RESPECT TO TEMPERATURE GRID
C NOTE: VELOCITY IS STAGGERED TO THE EAST AND SOUTH SIDE OF 
C EACH NODE 
C
      double precision vxs(ixzmax,ixzmax),vzs(ixzmax,ixzmax)
      common /velstg/vxs,vzs
C
C GRID
C
c      double precision pres(ixzmax),x(ixzmax),z(ixzmax)
      common /grid/xmax,zmax,dx,dz
C      common /grido/pres,x,z,nx,nz
C REMOVED ABOVE COMMON TO ALLOW NX TO TRANSFER IN FUNCTION
C
C PARAMETERS
C
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &       texpf,texps
C
C NEIGHBOUR COEFFICIENTS CACULATED USING ADVECTION ONLY
C GAMMA THE DIFFUSION PARAMETER IS ZERO, SO THE CONDUCTANCES
C D BECOME ZERO. THE PECLET NUMBER IS VERY LARGE (DOMINATED BY
C CONVECTION BUT THE POWER LAW AND PECLET NUMBER IS IRRELEVANT
C BECAUSE THE EQUATION FOR THE NEIGHBOUR COEFFICIENT MULTIPLIES
C THE A(P) FUNCTION BY D WHICH IS ZERO.
C
C F. IS THE MASS FLOW RATE, R IS THE DENSITY
C
C EASTERN NEIGHBOUR COEFFICIENT
C
      if (ix.lt.nx) then
         re=rhos
         fe=re*vxs(ix,iz)*dz
         anbe=max(-fe,0.)
      else
         anbe=0.
      endif
C
C WESTERN NEIGHBOUR COEFFICIENT
C
      if(ix.gt.1) then
        rw=rhos
        fw=rw*vxs(ix-1,iz)*dz
        anbw=max(fw,0.)
      else
        anbw=0.
      endif
C
C NORTHERN NEIGHBOUR COEFFICIENT
C
      if(iz.gt.1) then
        rn=rhos
        fn=rn*vzs(ix,iz-1)*dx
        anbn=max(-fn,0.)
      else
        anbn=0.
      endif
C
C SOUTHERN NEIGHBOUR COEFFICIENT
C
      if(iz.lt.nz) then
        rs=rhos
        fs=rs*vzs(ix,iz)*dx
        anbs=max(fs,0.)
      else
        anbs=0.
      endif
C
C END OF SUBROUTINE
C
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine RCNVGE(array,rdfmax,rdfmea,nx,nz,m0,m1)
C
C TEST FOR CONVERGENCE BETWEEN STARTING AND FINISHING
C RESIDUE FIELDS. RETURNS:
C MAXIMUM RESIDUE DIFF, RDFMAX
C MEAN RESIDUE DIFF, RDFMEA
C
C REAL DECLARATION
C
      implicit double precision(a-h,o-z)
C
C MAXIMUM SIZE OF GRID
C
      parameter (ixzmax=2001)
C 
C RESIDUE GRID
C
      double precision array(3,ixzmax,ixzmax)
C
C COMPARE RESIDUE GRIDS
C
      rdfmax=0.
      rdfmea=0.
      do ix=1,nx
       do iz=1,nz
         dif=abs(array(m1,ix,iz)-array(m0,ix,iz))
         if(dif.gt.rdfmax) rdfmax=dif
         rdfmea=rdfmea+dif
       enddo
      enddo
      rdfmea=rdfmea/float(nx*nz)
C
C END SUBROUTINE
C
      return
      end
C
C----------------------------------------------------------
C
