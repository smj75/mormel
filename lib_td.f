C
C THE FOLLOWING SUBROUTINES ARE HERE
C	TEMBAS
C	POISEUILLE
C
C
C CALCULATE TEMPERATURE FIELD ALONG BOTTOM BOUNDARY OF MODEL
C ASSUMING THAT THE CORNER FLOW MODEL IS SITTING ON TOP OF
C AN ASTHENOSPHERE CHANNEL WITH RADIALLY SYMMETRICAL OUTFLOW
C WITHIN A PLUME HEAD
C
      subroutine TEMBAS(temmax,tpmax)
      implicit double precision(a-h,o-z)
C
C MAXIMUM ARRAY SIZES
C
      parameter (iamax=10)
      parameter (itimax=2001)
      parameter (ixzmax=2001)
C
      double precision pres(ixzmax)
      double precision x(ixzmax)
      double precision z(ixzmax)
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      character*30 tfilei,tfileo
C
C COMMUNICATION WITHIN PROGRAM
C
      common /calc/tend,icalc
      common /corner/a,b
      common /countr/itym,ntime,iendpr
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
      common /temo/tem,tems
      common /tem/ptem,tfilei,tfileo,item
			common /minusa/aplumex,aplumey,aplumer,
     &               abslith,absasth,aplvflx
			common /minusb/bpuldur,bpulst,bpultem,ibpul
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
C PLATE SPREADING HISTORY AND CORNER FLOW GEOMETRY
C 
      double precision ufull(iamax)
      double precision utime(iamax)
      common /vel/ufull,utime,alpha,viss,crust,nu,ialpha
C
C POTENTIAL TEMPERATURE
C  TP, POTENTIAL TEMPERATURE
C  TPS, STARTING POT TEMP
C
      double precision tp(3,ixzmax,ixzmax)
      double precision tps(ixzmax,ixzmax)
      common /pottem/tp,tps
C
C TIME
C
c      return
      time = tyme(itym)
C
C DEFINE VSR MODEL PARAMETERS: get_user_parameters ( int argc, char **argv )
C
c      bpuldur = 0.5 * 1.0e6
c      bpulst = 2. * 1.0e6
c      bpulst = 200. * 1.0e6
c      bpultem = 55.
c      aplumex = 150. * 1.0e3 
c      aplumey = 0. * 1.0e3
c      abslith = 100. * 1.0e3
c      absasth = 200. * 1.0e3
c      aplvflx = 32.e6 * (1.0e9 / 1.0e6)
c      ryslice = 950. * 1.0e3

c      write(6,*)"TIME: ",time,tyme(itym)
      temmax = -999.0
      tpmax = -999.0
C
C IF RADIAL PLUME COORDINATE IS PROVIDED, 
C CALCULATE ALONG RIDGE DISTANE
C
      if ( aplumer .ne. 0.0 ) then
        aplumey = sqrt(aplumer**2 - aplumex**2)
			endif
			ryslice = sqrt(aplumex**2 + aplumey**2)
C
C CALCULATE VARIOUS MODEL INPUTS: plume_velocity_field ()
C
      TWOPI = 2.0*pi
      asthorig = (-0.5) * (absasth + abslith)
      asthhh = 0.5 * (absasth - abslith)
      avaflx = aplvflx / 2.0 / asthhh
      pulaflx = 3.0 * avaflx / TWOPI
C
C CALCULATE VARIOUS MODEL INPUTS: analytic_radial ()
C
      totait = bpultem * avaflx * bpuldur
      cnorm = totait / avaflx / bpuldur
      ckt = pulaflx * time
      cktrel = pulaflx * ( time - bpuldur )
      cqt = avaflx * time
      cqtrel = avaflx * ( time - bpuldur )
      tdelsq = 2.0 * bpuldur * bpuldur
      widbase = abslith / tan(alpha)
c      write(6,*)"avaflx",avaflx
c      write(6,*)"totait:",totait
c      write(6,*)"cnorm:",cnorm
c      write(6,*)"ckt, cktrel:",ckt,cktrel
c      write(6,*)"cqt, cqtrel:",cqt,cqtrel
c      write(6,*)"tdelsq:",tdelsq
c      write(6,*)"widbase:",widbase
C
C FIND INTERSECTION OF HORIZONTAL STREAMLINE AND 
C OUTER EDGE OF TRIANGULAR MELTING REGION
C
      if ( aplumex .ge. widbase ) then
        xedge = widbase
      elseif ( aplumex .le. -widbase ) then
        xedge = -widbase
      else
        xedge = aplumex
      endif
c      write(6,*)"xedge:",xedge
C
C HORIZONTAL LOOP
C
      do 1 ix=1,nx
        bastp = ptem
        dtp = 0.0
C
C IBPUL=0: NO PULSE
C
				if (ibpul.eq.0) then
          dtp = 0.0
c					t1 = tembkgr
C
C IBPUL=10: SQUARE WAVE
C
				else if (ibpul.eq.10) then
					n = int((time-bpulst)/bpuldur/2.0)
          tmid = (n+0.5)*2.0*bpuldur
					if ( (time-bpulst-tmid) .lt. 0.0 ) then
            dtp = 0.5*bpultem
c						t1 = tembkgr + 0.5*bpultem
					else
            dtp = (-0.5)*bpultem
c						t1 = tembkgr - 0.5*bpultem
					endif				
C
C IBPUL=1: SINGLE SQUARE
C
				else if (ibpul.eq.1) then
					if ( time.ge.bpulst .and. time.le.(bpulst+bpuldur) ) then
            dtp = bpultem
c						t1 = tembkgr + bpultem
					else
            dtp = 0.0
c						t1 = tembkgr
					endif
C
C IBPUL=20: SINE WAVE
C
				else if (ibpul.eq.20) then
          dtp = bpultem*sin(TWOPI*(time-bpulst)/bpuldur)
c					t1 = tembkgr +
c     &      bpultem*sin(TWOPI*(time-bpulst)/bpuldur)
C
C IBPUL=2: SINGLE COSINE
C
				else if (ibpul.eq.2) then
					if ( time.ge.bpulst .and. time.le.(bpulst+bpuldur) ) then
            dtp = bpultem*(0.5-0.5*cos((time-bpulst)/TWOPI*bpuldur))
c						t1 = tembkgr + 
c     &        bpultem*(0.5-0.5*cos((time-bpulst)/TWOPI*bpuldur))
					else
            dtp = 0.0
c						t1 = tembkgr
					endif
C
C IBPUL=3: SINGLE GAUSSIAN
C
				else if (ibpul.eq.3) then
					if ( time.gt.(bpulst-5.0*bpuldur) .and.
     &         time.lt.(bpulst+5.0*bpuldur) ) then
            dtp = bpultem*exp((-1.0)*((time-bpulst)/bpuldur)**2)
c						t1 = tembkgr + 
c     &        bpultem*exp((-1.0)*((time-bpulst)/bpuldur)**2)
					else
            dtp = 0.0
c						t1 = tembkgr
					endif
C
C IBPUL=300: TEST BY SIMPLY WRITING A GAUSSIAN CURVE
C
				else if (ibpul.eq.300) then
	        xpos=1.0-x(ix)/widbase
	        if (xpos.lt.0.0) xpos=0.0
  	      dtp = bpultem*xpos*exp((-1.0)*((time-bpulst)/bpuldur)**2)
c  	      t1 = tembkgr+bpultem*xpos
c     &  	     *exp((-1.0)*((time-bpulst)/bpuldur)**2)
C
C INTEGRATE TO FIND AVERAGE VERTICAL UPWELLING RATE
C
				else if (ibpul.eq.5) then
        	if ( x(ix) .lt. widbase ) then
          	x1 = x(ix)
        	else
         	  x1 = widbase
        	endif
        	x1mx0 = x1 - xedge
        	if ( x1mx0 .ne. 0.0 ) then
          	x1ol = x1 / abslith
          	x0ol = xedge / abslith
          avupwell = ( a*x1mx0 + b*abslith*( atan(x0ol) - atan(x1ol) ) 	)
     &             / x1mx0
        	else
         	 	avupwell = 0.0
        	endif
c        	if ( x(ix) .gt. widbase ) avupwell = 0.0
C
C READ TEMPERATURE FROM MODEL OF RUDGE ET AL., EPSL, 2008
C
       		pirsq = PI * ( ( (x(ix)-aplumex)**2.0 +
     &                   (ryslice-aplumey)**2.0) )
        	zstar = asthhh + avupwell*time
        	tstar = time - bpulst
        	fzh = 3.0 / 2.0 * ( 1.0 - ( (zstar/asthhh)**2.0) )
        	if (fzh.gt.0.) then
          	dtp = cnorm*
     &      	    exp((pirsq/avaflx/fzh - tstar)**2.0 /(-tdelsq) )
c          	t1 = tembkgr + cnorm*
c     &      	    exp((pirsq/avaflx/fzh - tstar)**2.0 /(-tdelsq) )
c         	 write(6,*)ix,zstar,fzh,tstar,t1
        	else
            dtp = 0.0
c         	 t1 = tembkgr
c         	 write(6,*)ix,t1
        	endif


				endif
C
C BASAL POTENTIAL TEMPERATURE
C
        bastp = ptem + dtp
C
C BASAL REAL TEMPERATURE
C
        bastem = (bastp+273.)/exp(g*texps*z(nz)/cp)-273.
C
C LOAD BASAL TEMEPERATURE INTO TEMPERATURE ARRAYS
C
c      	 if (ix.eq.1) write(6,*)ix,t1
        tems(ix,nz)=bastem
        tem(1,ix,nz)=bastem
        tem(2,ix,nz)=bastem
        tem(3,ix,nz)=bastem
C
C LOAD BASAL TEMPERATURE INTO POTENTIAL TEMPERATURE ARRAYS
C

        tps(ix,nz)=bastp
        tp(1,ix,nz)=bastp
        tp(2,ix,nz)=bastp
        tp(3,ix,nz)=bastp    
C
C STORE MAXIMUM BASAL TEMPERATURE
C
        if ( tem(1,ix,nz) .gt. temmax ) temmax = tem(1,ix,nz)
        if ( tp(1,ix,nz) .gt. tpmax ) tpmax = tp(1,ix,nz)
C
C WRITE OUT TO TEST
C
c    2   format(i3,1x,f10.3,1x,f9.4,1x,f12.2,1x,f16.5)
c      	 if (ix.eq.1) write(6,*)ix,t1
c        write(6,2) ix, x(ix)*1.0e-3, avupwell*1.0e3, zstar*1.0e-3,
c     &			   tem(n0,ix,nz)

C
C END OF HORIZONTAL LOOP
C
    1 enddo
c	  write (6,*) "TEMBAS: Max. basal temp. ",temmax
C
      return
      end
C
