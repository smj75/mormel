C
C ROUTINES FOR MELTING CALCULATIONS
C STEPHEN JONES SUMMER 2005
C RACHEL WALTERS 2006...
C 
C CONTAINS THE FOLLOWING ROUTINES
C	MEL
C	INTEG
C	FRCMEL
C	SOL
C	TLHF
C	DTDA
C	XDS
C       DFXDS
C	RTBIS2,3,4
C
C----------------------------------------------------------
C
CR PUT RE1 IN MEL NAME
      subroutine MEL(n1,re1)
      implicit double precision(a-h,o-z)
C
C DETERMINE TOTAL AMOUNT OF MELT PRODUCED EITHER IN STEADY STATE
C OR IN THIS TIMESTEP
C
      parameter (ixzmax=2001)
      double precision pres(ixzmax)
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      double precision vx(ixzmax,ixzmax)
      double precision vz(ixzmax,ixzmax)
      double precision x(ixzmax),z(ixzmax)
      double precision frac(ixzmax,ixzmax)
      double precision temc(ixzmax,ixzmax)
      double precision dpda(ixzmax,ixzmax)
      double precision dtrda(ixzmax,ixzmax)
      double precision dfda(ixzmax,ixzmax)
      double precision dtuda(2,ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      double precision ftemp(ixzmax,ixzmax)
      double precision dtsda(ixzmax,ixzmax)
      double precision dtida(ixzmax,ixzmax)
      character*30 tfilei,tfileo
C
C COMMUNICATION WITHIN PROGRAM
C
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
      common /mf/p,t,xh2o,cpxm,frac,nfrac
      common /tem/ptem,tfilei,tfileo,item
      common /temo/tem,tems
C TEMCOR=CORRECTED TEMPERATURE
      common /temcor/temc
      common /vel/ufull,alpha,viss,crust,ialpha
      common /velo/vx,vz,uhalf,vmax
C COMMON FOR MELT FRACTION FOR DF/Dt NUMERICAL SOLUTION
      common /melfrc/f1x,f2x,f1z,f2z,fx
C COMMUNICATION OF DERIVATIVES THROUGHOUT THE PROGRAM
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
      common /intder/dpdac,dtrdac,f,dfdac
      common /dermb/dtsdac,dtudac,dtuda,dtsda
      common /intgra/tarea,tcrust
      common /ftempy/ftemp
      common /deltem/dtem
      common /temtes/itest
C
C COMMUNICATION WITH FUNCTION XDS
C ACCURACY OF MELT FRACTION WHEN FOUND BY ROOT BISECTION 
C
      external XDS
      common /tinit/ti
      parameter (xacc=1.e-6)
C
C TEST FUNCTION TO SEE IF THERE IS ENOUGH MELTING MORE
C THAN RESIDUE
C
      external XDTEST
C
C COMMUNICATION WITH A NEW FUNCTION DFXDS, TO CALCULATE
C TIME DEPENDENT CHANGE IN REAL TEMPERATURE
C
      external DFXDS
      common /derkat/dtidac,dtida,fi
      parameter (dfacc=1.e-6)
C
C COMMON - RESIDUE
C
      double precision rest(ixzmax,ixzmax)
      double precision res(3,ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C XH2O AND CPXM HAVE TO BE COPIED SINCE THEY ARE IN TWO COMMON BLOCKS
C
      xh2o=xh2o1
      cpxm=cpxm1
      itc=n1
      irec=re1
      t1=-999.0
C
C 
C SET INTEGRATION TOTAL EQUAL TO ZERO
C      tarea=0.
C
C IF NO MELTING CALCULATION IS SELECTED THEN TESTING TEMP
C EVOLUTION
C
      if(imelt.eq.0) then
       do ix=1,nx
         do iz=1,nz
         temc(ix,iz)=tem(n1,ix,iz)
         enddo
       enddo
       tcrust=crust
       goto 2
      else
       continue
      endif

C LOOPS OVER VERTICAL AND HORIZONTAL
C
      do iz=1,nz
         do ix=1,nx
            frac(ix,iz)=0.
            nfrac(ix,iz)=0.
C
C LOOPS FOR NUMERICAL SOLUTION TO DF/Dt
C      do iz=29,31
C         do ix=1,3
C
C ANY MELTING?
C
            p=pres(iz)
            tmp=tem(n1,ix,iz)
            tsol=SOL(iparam)
C            write(6,*)'p',p,'tmp',tmp,'tsol',tsol
            if (tmp.le.tsol .and. xh2o.eq.0.) then
C            if (tmp.le.tsol) then
C            write(6,*)'Temperature below solidus'
C
C IF TEMPERATURE DOES NOT NEED TO BE CORRECTED FOR LATENT HEAT 
C OF MELTING THEN TRANSFER TO CORRECTED TEMP STRUCTURE
             frac(ix,iz)=1.-res(re1,ix,iz)
C             frac(ix,iz)=0.
             nfrac(ix,iz)=0.
             temc(ix,iz)=tem(n1,ix,iz)
C             rest(ix,iz)=res(re1,ix,iz)
C NOT NEEDED HERE BECAUSE DON'T WANT TO CYCLE THE NEW RESIDUE
C TILL CRUST HAS BEEN ITERATED.           
            goto 1
            endif
C             
C WRITE OUT P AND TSOL FOR KATZ FIGURE 1 (ABOVE)
C
ccc            write(6,*)'MELTING'
C
C CORRECT THERMAL STRUCTURE FOR LATENT HEAT OF MELTING USING
C PARAMETERISATION OF WATSON AND MCKENZIE (1992)
C BUT ONLY IF STEADY-STATE CALCULATION
C
            if (imelt.eq.1 .and. xh2o.eq.0.) then
               t=TLHF(tmp,p,tsol)
               t1=t
C               write(6,*)'t in,t out,tsol:',tmp,t,tsol
            elseif(imelt.eq.1 .and. xh2o.gt.0.) then
                write(6,*)'Cannot do hydrous McKenzie & Bickle'
                stop
C
C USE METHOD OF BOWN & WHITE 1995
C IS THIS INTENDED TO BE A TIME-DEPENDENT CALCULATION? IS THIS
C JUST THE SAME AS THE W&M PARAMETERISATION
C

C
C USE METHOD OF KATZ ET AL. 2003
C
            elseif (imelt.eq.3 .and. xh2o.eq.0.) then
               ti=tmp
               resc=res(re1,ix,iz)
C TEST TO SEE IF THERE IS ENOUGH MELTING TO STILL MELT WHEN RESIDUE
C IS TAKEN INTO ACCOUNT. MAX MELT FRACTION AT TMP, TEMP IS REDUCED
C AND CORRECTED FOR LH WHEN MELTING OCCURS.
C
               mtest=XDTEST(tmp,itest)
               if(itest.eq.1) then
                  frac(ix,iz)=1.-resc
                  nfrac(ix,iz)=0.
                  temc(ix,iz)=tem(n1,ix,iz)
                  goto 1
               endif
C               write(6,*)'Temp: ',ti,'Tsol: ',tsol, 'p:',p
               t=RTBIS2(XDS,tsol,ti,xacc)
C               write(6,*)'Found solution to rtbis2'              
               t1=t
C
C KATZ HYDROUS PARAMETERISATION: TEMPERATURE MUST BE ABOVE THE
C LOWEST POSSIBLE TEMPERATURE FOR THE SATURATED SOLIDUS
C BUT THERE MAY STILL NOT BE ANY MELTING
C XSOLB IS THE BOUNDARY TEMPERATURE FOR THE HYDROUS SOLIDII
C
            elseif (imelt.eq.3 .and. xh2o.gt.0.) then
               xsolb=900.
C               xsolb=tsol
               resc=res(re1,ix,iz)
               if(tmp.le.xsolb) then
               temc(ix,iz)=tem(n1,ix,iz)
               frac(ix,iz)=1.-resc
               nfrac(ix,iz)=0.
C               rest(ix,iz)=res(re1,ix,iz)
               goto 1
               endif
               mtest=XDTEST(tmp,itest)
               if(itest.eq.1) then
                  frac(ix,iz)=1.-resc
                  nfrac(ix,iz)=0.
                  temc(ix,iz)=tem(n1,ix,iz)
                  goto 1
               endif
               ti=tmp
C               write(6,*)'start ti:',ti
C               t=RTBIS3(XDS,tsol,ti,xacc)
C TSOL IS LOWERED WHEN THE HYDROUS CALCULATION IS IMPLEMENTED
C SO USING TSOL AS THE LOWER BOUND WOULD NOT ACCOUNT FOR ALL
C MELTING.
               t=RTBIS3(XDS,xsolb,ti,xacc)
               if(t.eq.0.) then
               t=tem(n1,ix,iz)
               frac(ix,iz)=1.-resc
               nfrac(ix,iz)=0.
               temc(ix,iz)=tem(n1,ix,iz)
               goto 1
               endif
               t1=t
C               write(6,*)'ti:',ti,'tout:',t1
C            elseif (imelt.eq.0) then
C            temc(ix,iz)=tem(n1,ix,iz)
            else
            write(6,*)'Incorrect options for imelt and xh2o'
            endif
C IF NO MELTING CALCULATION IS SELECTED THEN TESTING TEMP
C EVOLUTION
C      if(imelt.eq.0) then
C       tcrust=crust
C       goto 2
C      else
C       continue
C      endif
C
C TEMPERATURE FIELD HAS NOW BEEN CORRECTED FOR LATENT HEAT OF
C MELTING, SAVE AS AN ARRAY IN IX AND IZ
C
            temc(ix,iz)=t1 
C
C CALCULATE THE MELT FRACTION USING CHOSEN PARAMETERISATION
C
C            write(6,*)'ix,iz:',ix,iz
            frac(ix,iz)=frcmel(iparam)
C            write(6,*)'Melt Fraction: ',frac(ix,iz)
C            if (frac(ix,iz).gt.0. .and. vz(ix,iz).eq.0.) then
C           write(6,*)'Melting either in crust/outside 
C     &                      wedge Nodes:',ix,iz
C            stop
C            endif
C
C CLAUSE TO STOP PROGRAM RUNNING IF MELTING OCCURS AT THE BOTTOM
C BOUNDARY I.E. IZ=NZ, THE GRID NEEDS TO BE EXTENDED!
C
C            if(iz.eq.nz .and. nfrac(ix,iz).gt.0.) then
C            write(6,*)'MELTING AT BASE - EXTEND GRID'
C            stop
C            endif
 1        continue 
        enddo
      enddo
C
C
C
C ONCE MELT FRACTION HAS REACHED A MAXIMUM, ANY FURTHER NODES
C ZERO. WHEN THE MELT FRACTION STARTS TO DECREASE, IT IS
C ASSUMED THAT ALL THE MELT IS EXTRACTED IMMEDIATELY WITH
C NO FURTHER REACTION WITH THE MATRIX AND SOLIDIFIES TO
C FORM THE OCEANIC CRUST
C
CR      do ix=1,nx
CR         do iz=1,nz
CR            ftemp(ix,iz)=0.
CR            if(iz.eq.nz) then
CR               ftemp(ix,iz)=0.
CR            else         
CR               ftemp(ix,iz)=frac(ix,iz)-frac(ix,iz+1)
CR            endif
CR         enddo
CR      enddo
C
C WRITE OUT THE MELT FRACTION SO THAT CAN CHECK EXACTLY WHAT DEPTH
C MELTING STARTS AND FINISHES, AND WHAT THE MAXIMUM MELT FRACTION IS.
C      
CR      open(27,file='bw94_meltf_out')
CR      do ix=1,nx
CR         do iz=1,nz
CR         if (ftemp(ix,iz).lt.0.) then
CR         frac(ix,iz)=0.
C ABOVE THE MAX MELT FRACTION, ASSUME THAT NO FURTHER MELTING OCCURS
C AND THEREFORE NO CORRECTION OF TEMPERATURE IS REQUIRED 
C IF YOU DON'T CORRECT THE TEMP, THEN END UP WITH A FUNNY INCREASE IN
C THE TEMP AS YOU DECREASE IN DEPTH
C THIS IS NOT QUITE ACCURATE, NEED TO CONSIDER WHERE THE PACKET IS COMING
C FROM EITHER VERTICAL OR HORIZONTAL??
C          if(vz(ix,iz).ge.vx(ix,iz))then
CR          temc(ix,iz)=tem(n1,ix,iz)-(tem(n1,ix,iz+1)-temc(ix,iz+1))
C          else
C          temc(ix,iz)=tem(n1,ix,iz)-(tem(n1,ix-1,iz)-temc(ix-1,iz))
C          endif
CR         endif
CR         write(27,*)x(ix)*1.e-3,z(iz)*1.e-3,frac(ix,iz),
CR     & tem(n1,ix,iz),temc(ix,iz)
CR         enddo
CR      enddo
CR      close(27)    
C
C CALCULATE ACTUAL MELT FRACTION AT SOME POINT BY CONSIDERING
C THE RESIDUE THE ABOVE CALCULATION OF WHERE THE MELT FRACTION
C DECREASES IS NOT REQUIRED IF THE RESIDUE IS CONSIDERED, BUT 
C MAY NEED A BETTER ESTIMATE/RECALCULATION OF THE LATENT HEAT
C OF MELTING
C
C FIND CORRECTED MELT FRACTION (AT THE MOMENT KATZ ONLY)
C
      do ix=1,nx
        do iz=1,nz
C          if(nfrac(ix,iz).eq.0.) then
C          goto 15
C          endif
          nfrac(ix,iz)=frac(ix,iz)-(1.-res(re1,ix,iz))
C          write(6,*)'nfrac:',nfrac(ix,iz),frac(ix,iz), res(re1,ix,iz)
C          write(6,*)'re1',re1
C CLAUSE: IF THE MELT FRACTION IS SMALLER THAN WHAT HAS
C ALREADY BEEN REMOVED FROM RESIDUE THEN NO MELTING OCCURS!
          if (nfrac(ix,iz).le. 0.) then
          nfrac(ix,iz)=0.
          frac(ix,iz)=1.-res(re1,ix,iz)
C CHANGED TO SEE IF IT IS RESIDUE AND FRAC GETTING MIXED UP FOR
C INTEGRAL. IT DOESN'T!
C           frac(ix,iz)=0.
          endif
C MOVED FROM FRAC(IX,IZ) JUST ABOVE
          if(iz.eq.nz .and. nfrac(ix,iz).gt.0.) then
             write(6,*)'MELTING AT BASE - EXTEND GRID'
             stop
          endif
15         continue
        enddo
      enddo
C
C 
C WRITE OUT CORRECTED TEMPERATURE STRUCTURE TEMC(IX,IZ)
C
cc      open(19,file='bw94_temc_out')
ccC5    format (3(f10.2,1x))
cc      do iz=1,nx
cc         do ix=1,nz
cc            write(19,*)x(ix)*1.e-3,z(iz)*1e-3,temc(ix,iz)
cc         enddo
cc      enddo
cc      close(19)
cc      write(6,*)'Written corrected temperature field to "bw94_temc"'
C
C WRITING TEMP FIELD OUT TO CHECK CURVATURE CALCULATIONS
C
C 5    format (4(f10.2,1x))
C      open(20,file='bw94_temps_out')
C      do iz=1,nz
C         do ix=1,nx
C         write(20,5)x(ix)*1.e-3,z(iz)*1.e-3,tem(n1,ix,iz),temc(ix,iz)
C         enddo
C      enddo
C      close(20)
C
C--------------------------------------------------------------------
C CALCULATE DT/Dt FROM THE CORRECTED TEMPERATURE FIELD ARRAY TEMC
C EXPLAIN THE EQUATIONS ETC.
C
C
C      write(6,*)'calculating DT/Da'
C
C CALCULATE DTDA FOR MB88 PARAMETERISATION
C
      if(iparam.eq.1) then
       if(itemp.eq.1) then
C
C CALCULATE THE TWO PARTS OF THE ADVECTION DIFFUSION EQUATION FOR
C THE UNCORRECTED TEMPERATURE FIELD AS THEY SHOULD BE EQUAL
C
       open(24,file='wj_dtuda_out')
        do ix=1,nx
          do iz=1,nz
C CALCULATE THE ADIABATIC TEMPERATURE GRADIENT
            ch=g*texps/cp
            h=ch*(tem(n1,ix,iz)+273.)
            if (ix.eq.1) then
             dtudx=(tem(n1,ix+1,iz)-2.*tem(n1,ix,iz)
     &              +tem(n1,ix+1,iz))/(dx**2.)
            else if (ix.eq.nx) then
             dtudx=0.
C BECAUSE THIS BOUNDARY IS NOT SYMMETRICAL
            else
             dtudx=(tem(n1,ix+1,iz)-2.*tem(n1,ix,iz)
     &              +tem(n1,ix-1,iz))/(dx**2.)
            endif
            if(iz.eq.1 .or. iz.eq.nz) then
            dtudz=0.
            else
            dtudz=(tem(n1,ix,iz+1)-2.*tem(n1,ix,iz)
     &            +tem(n1,ix,iz-1))/(dz**2.)
            endif
            dtuda(1,ix,iz)=tdiff*(dtudx+dtudz)-vz(ix,iz)*h
            write(24,*)ix,iz,dtudx,dtudz,vz(ix,iz),h,dtuda(1,ix,iz)
C 2          continue
          enddo
        enddo
C        write(6,*)'Written dtuda out'
       close(24) 
C
C CALCULATE V.GRADT OF UNCORRECTED TEMPERATURE TO 
C CHECK IF IS SAME AS
C
       open(23,file='wj_dtudacheck_out')
         do ix=1,nx
            do iz=1,nz
              if(ix.eq.1) then
              xadv=(((tem(n1,ix+1,iz)+tem(n1,ix,iz))/2.)
     &                 -tem(n1,ix,iz))/(0.5*dx)
C              xadv=0.
              else
              xadv=(tem(n1,ix,iz)-tem(n1,ix-1,iz))/dx
C              xadv=(((tem(n1,ix,iz)+tem(n1,ix+1,iz))/2.)
C     &              -((tem(n1,ix,iz)+tem(n1,ix-1,iz))/2.))/dx
              endif
C              if(iz.eq.1) then
C              zadv=(tem(n1,ix,iz)-((tem(n1,ix,iz+1)+
C     &                 tem(n1,ix,iz))/2.))/(dz/2.)
              if(iz.eq.nz) then
C              zadv=(((tem(n1,ix,iz-1)+tem(n1,ix,iz))/2.)
C     &                 -tem(n1,ix,iz))/(dz/2.)
C              zadv=tem(n1,ix,iz-1)-tem(n1,ix,iz)
              zadv=0.
              else
              zadv=(tem(n1,ix,iz)-tem(n1,ix,iz+1))/dz
C              zadv2=(((tem(n1,ix,iz)+tem(n1,ix,iz-1))/2.)
C     &              -((tem(n1,ix,iz)+tem(n1,ix,iz+1))/2.))/dz
              endif
              dtuda(2,ix,iz)=(zadv*vz(ix,iz)+xadv*vx(ix,iz))
              write(23,*)ix,iz,dtuda(1,ix,iz),dtuda(2,ix,iz),
     &                   dtuda(1,ix,iz)-dtuda(2,ix,iz)
            enddo
         enddo
       close(23)
C
       open(26,file='wj_dtrda1_out')
         do ix=1,nx
            do iz=1,nz
              dtrda2=0.
              if(ix.eq.1) then
C              xadv=temc(ix+1,iz)-temc(ix,iz)
              xadv=(((temc(ix+1,iz)+temc(ix,iz))/2.)
     &                 -temc(ix,iz))/(0.5*dx)
C              xadv=0.
              else
              xadv=(temc(ix,iz)-temc(ix-1,iz))/dx
              endif
              if(iz.eq.nz) then
C              zadv=temc(ix,iz-1)-temc(ix,iz)
C
C IF THE MELTING STARTS AT THE BASE OF THE GRID
C MAKE THE GRID BIGGER SO THAT THE IZ=NZ IS IRRELEVANT
C BECAUSE NO MELTING IS HAPPENING THERE.
C THIS SHOULD ALWAYS BE THE CASE WITH NORMAL POT TEMP
C BUT IF YOU PUT IN 1500 THEN MELTING STARTS BEFORE 100KM
C
              zadv=0.
              else
              zadv=(temc(ix,iz)-temc(ix,iz+1))/dz
              endif
              dtrda(ix,iz)=(zadv*vz(ix,iz)+xadv*vx(ix,iz))
              write(26,*)ix,iz,dtrda(ix,iz)
            enddo
         enddo
       close(26)
C
C---------------------------------------------------------
C
C CALCULATE DT/Dt USING BOWN & WHITE 1995
C
       else if(itemp.eq.0) then 
C
C CALCULATE DTUDA, USING THE DIFFUSION PART OF THE EQN
C
       open(24,file='wj_dtuda_out')
        do ix=1,nx
          do iz=1,nz
C CALCULATE THE ADIABATIC TEMPERATURE GRADIENT
            ch=g*texps/cp
            h=ch*(tem(n1,ix,iz)+273.)
            if (ix.eq.1) then
             dtudx=(tem(n1,ix+1,iz)-2.*tem(n1,ix,iz)
     &              +tem(n1,ix+1,iz))/(dx**2.)
C             dtudx=0.
C MAKING DTUDX 0 MAKES NO DIFFERENCE AT ALL, WHY IS THIS?
            else if (ix.eq.nx) then
             dtudx=0.
C BECAUSE THIS BOUNDARY IS NOT SYMMETRICAL
            else
             dtudx=(tem(n1,ix+1,iz)-2.*tem(n1,ix,iz)
     &              +tem(n1,ix-1,iz))/(dx**2.)
            endif
            if(iz.eq.1 .or. iz.eq.nz) then
            dtudz=0.
            else
            dtudz=(tem(n1,ix,iz+1)-2.*tem(n1,ix,iz)
     &            +tem(n1,ix,iz-1))/(dz**2.)
            endif
            dtuda(1,ix,iz)=tdiff*(dtudx+dtudz)-vz(ix,iz)*h
C            write(24,*)ix,iz,dtudx,dtudz,vz(ix,iz),h,dtuda(1,ix,iz)
            write(24,*)ix,iz,dtuda(1,ix,iz)
          enddo
        enddo
C       write(6,*)'Written dtuda out'
       close(24)
C
C TEST (ONLY FOR STEADY STATE) TO SEE IF CALCULATING USING V.GRADT
C INSTEAD OF THE
C DIFFUSION PART OF THE EQUATION WILL MAKE A DIFFERENCE TO THE FINAL
C OUTCOME
C
C RESULT: MAKES THE FINAL CRUST SLIGHTLY THICKER THAN  USING THE
C DIFFUSION PART OF THE EQUATION. THE TWO WAYS DON'T PRODUCE EXACTLY
C THE SAME DTUDA AS CHECKED PREVIOUSLY, MUST BE TO DO WITH THE ACCURACY
C OF THE STEADY-STATE CALCULATION.
C      open(24,file='wj_dtuda_out')
C         do ix=1,nx
C            do iz=1,nz
C              if(ix.eq.1) then
C              xadv=(((tem(n1,ix+1,iz)+tem(n1,ix,iz))/2.)
C     &                 -tem(n1,ix,iz))/(0.5*dx)
CC              xadv=0.
C              else
C              xadv=(tem(n1,ix,iz)-tem(n1,ix-1,iz))/dx
CC              xadv=(((tem(n1,ix,iz)+tem(n1,ix+1,iz))/2.)
CC     &              -((tem(n1,ix,iz)+tem(n1,ix-1,iz))/2.))/dx
C              endif
CC              if(iz.eq.1) then
CC              zadv=(tem(n1,ix,iz)-((tem(n1,ix,iz+1)+
CC     &                 tem(n1,ix,iz))/2.))/(dz/2.)
C              if(iz.eq.nz) then
CC              zadv=(((tem(n1,ix,iz-1)+tem(n1,ix,iz))/2.)
CC     &                 -tem(n1,ix,iz))/(dz/2.)
CC              zadv=tem(n1,ix,iz-1)-tem(n1,ix,iz)
C              zadv=0.
C              else
C              zadv=(tem(n1,ix,iz)-tem(n1,ix,iz+1))/dz
CC              zadv2=(((tem(n1,ix,iz)+tem(n1,ix,iz-1))/2.)
CC     &              -((tem(n1,ix,iz)+tem(n1,ix,iz+1))/2.))/dz
C              endif
C              dtuda(1,ix,iz)=(zadv*vz(ix,iz)+xadv*vx(ix,iz))
C              write(24,*)ix,iz,dtuda(1,ix,iz)                   
C            enddo
C         enddo
C       close(24)
C USE DTDA FUNCTION TO CALCULATE DTRDA USING BOWN AND WHITE METHOD
C
        do ix=1,nx
          do iz=1,nz
          f=frac(ix,iz)
          p=pres(iz)
          t=temc(ix,iz)
          dtudac=dtuda(1,ix,iz)
          dpdac=dpda(ix,iz)
C          write(6,*)'dtudac:',dtudac,'dpdac:',dpdac
          if(f.le.0.) then
C
C IS THIS CORRECT? IT DOESN'T MATTER BECAUSE YOU ARE ONLY INTEGRATING
C OVER THE AREA WHERE THERE IS MELTING, BUT FOR TIME DEPENDENT
C CALCULATIONS NEED TO CHECK THAT THERE IS NO FEEDBACK FROM THIS
C CALCULATION. THINK IT SHOULD JUST BE THE TEMPERATURE FIELD AND
C MELTING THAT IS FED BACK
C
           dtrda(ix,iz)=dtudac
           goto 5
          else
CCsmj CHANGED FTEMP TO FTEMP(IX,IZ)
           ftemp(ix,iz)=dmb88(iparam)
           dtsda(ix,iz)=dtsdac
           dtrdac=dtda(itemp) 
           dtrda(ix,iz)=dtrdac
          endif
 5        continue
          enddo
        enddo
      endif
C
C IF PARAM=2, USE KATZ DERIVATIVES
C
      elseif(iparam.eq.2) then    
       if(itemp.eq.1) then
C
C JUST NEED TO DELETE THIS NEXT BIT ONCE THE
C HYDROUS DERIVATIVES HAVE BEEN FIXED
C 
C        if(xh2o.gt.0.) then
C        write(6,*)'Cannot do hydrous derivatives yet, due to
C     &             unsaved values'
C         stop
C        endif
       write(6,*)'Use steady-state calc vgradt'
        open(26,file='wj_dtrda1_out')
         do ix=1,nx
            do iz=1,nz
              dtrda2=0.
              if(ix.eq.1) then
C              xadv=temc(ix+1,iz)-temc(ix,iz)
              xadv=(((temc(ix+1,iz)+temc(ix,iz))/2.)
     &                 -temc(ix,iz))/(0.5*dx)
C              xadv=0.
              else
              xadv=(temc(ix,iz)-temc(ix-1,iz))/dx
              endif
              if(iz.eq.nz) then
C              zadv=temc(ix,iz-1)-temc(ix,iz)
C
C IF THE MELTING STARTS AT THE BASE OF THE GRID
C MAKE THE GRID BIGGER SO THAT THE IZ=NZ IS IRRELEVANT
C BECAUSE NO MELTING IS HAPPENING THERE.
C THIS SHOULD ALWAYS BE THE CASE WITH NORMAL POT TEMP
C BUT IF YOU PUT IN 1500 THEN MELTING STARTS BEFORE 100KM
C
              zadv=0.
              else
              zadv=(temc(ix,iz)-temc(ix,iz+1))/dz
              endif
              dtrda(ix,iz)=(zadv*vz(ix,iz)+xadv*vx(ix,iz))
              write(26,*)ix,iz,dtrda(ix,iz)
            enddo
         enddo
        close(26)
C
C USE THE TIME DEPENDENT CALCULATION FOR CHANGE IN TEMP
C WITH TIME NEED TO DO NEW DIFFERENTIATING FOR THIS.
C
       elseif(itemp.eq.0) then
C
C CALCULATE THE CHANGE IN INTIAL TEMPERATURE W.R.T TO TIME
C USING K(D2T/DX2+D2T/DZ2)-VZH
C
C
cc        open(24,file='wj_dtida_out')
        do ix=1,nx
          do iz=1,nz
C CALCULATE THE ADIABATIC TEMPERATURE GRADIENT
            ch=g*texps/cp
            h=ch*(tem(n1,ix,iz)+273.)
            if (ix.eq.1) then
             dtidx=(tem(n1,ix+1,iz)-2.*tem(n1,ix,iz)
     &              +tem(n1,ix+1,iz))/(dx**2.)
            else if (ix.eq.nx) then
             dtidx=0.
C BECAUSE THIS BOUNDARY IS NOT SYMMETRICAL
            else
             dtidx=(tem(n1,ix+1,iz)-2.*tem(n1,ix,iz)
     &              +tem(n1,ix-1,iz))/(dx**2.)
            endif
            if(iz.eq.1 .or. iz.eq.nz) then
            dtidz=0.
            else
            dtidz=(tem(n1,ix,iz+1)-2.*tem(n1,ix,iz)
     &            +tem(n1,ix,iz-1))/(dz**2.)
            endif
            dtida(ix,iz)=tdiff*(dtidx+dtidz)-vz(ix,iz)*h
cc            write(24,*)ix,iz,dtida(ix,iz),vz(ix,iz),dtidx,dtidz,
cc     &                 dx,dz,tem(n1,ix,iz)
          enddo
        enddo
C       write(6,*)'Written dtida out crust katz:',crust*1.e-3
cc       close(24)
C
C USE RTBIS TO FIND THE VALUE OF THE CHANGE IN REAL TEMPERATURE
C W.R.T. TIME, WHICH IS INTERDEPENDENT WITH THE MELT PRODUCTIVITY
C
      do ix=1,nx
       do iz=1,nz
C CAREFUL THAT YOU HAVE EVERYTHING NEEDED TO RUN DK03 FOR EACH NODE!!!
        f=frac(ix,iz)
        fi=nfrac(ix,iz)
C		write(6,*)fi,f,res(re1,ix,iz)
C        write(6,*)'nfrac:',nfrac(ix,iz),fi
        p=pres(iz)
        t=temc(ix,iz)
        dpdac=dpda(ix,iz)
        dtidac=dtida(ix,iz)
        if(fi.le.0.) then
C        write(6,*)'fi less than 0.',fi
        dtrda(ix,iz)=dtidac
        goto 6
        else
C CHECK FOR A ROOT BY RUNNING DK03, IF THE
        dtr1=1
        dtr2=-1
C IF RETURN IS NO ROOT THEN EXPAND DTR2.
C        write(6,*)'ix,iz:',ix,iz
C        write(6,*)'Calling DFXDS'
        dtrdac=rtbis4(DFXDS,dtr1,dtr2,dfacc)
C        dtrdac=
        dtrda(ix,iz)=dtrdac
        endif
 6      continue
C        write(6,*)fi,dtrdac,dtidac
       enddo
      enddo
C
       endif
      endif
C
C-----------------------------------------------------------
C
C CALL THE INTEGRATION SUBROUTINE TO CALCULATE THE TOTAL MELT
C PRODUCTION RATE
C
      call INTEG(re1)
C
C WRITE OUT THE MELT RATE TO A FILE TO CHECK POSITIVES AND NEGATIVES
C
C      open(7,file='bw94_meltrate_out')
C      do ix=1,nx
C         do iz=1,nz
C            write(7,*)dfda(ix,iz)
C         enddo
C      enddo
C      close(7)
C
C CALCULATE THE TOTAL CRUSTAL THICKNESS P.443, BOWN AND WHITE 94
C
      tcrust=2.*tarea/(ufull)
cc      write(6,*)'Total Crust(km): ',tcrust*1.e-3
C
C END OF CRUSTAL ITERATION LOOP, ERRCST IS THE ERROR REQUIRED
C FOR THE CRUST DIFFERENCE BETWEEN INPUT AND OUTPUT TO BE LESS THAN
C OR EQUAL TO. DCRUST IS THE DIFFERENCE BETWEEN THE TWO CRUSTAL VALUES
C SHOULD RETURN TO SETVEL SUBROUTINE IN WJ.F
C      errcst=0.1.
C      dcrust=tcrust-crust
C      if(abs(dcrust).gt.errcst) then
C      write(6,*)'Crust iteration'
C      crust=tcrust
C      else
C      write(6,*)'Final Crust(km):',tcrust
C      endif
C
 2    continue
      return
      end
C
C---------------------------------------------------------------------
C
      subroutine INTEG(re1)
      implicit double precision(a-h,o-z)
C
C CALCULATE THE TOTAL VOLUMETRIC MELT PRODUCTION RATE FROM
C INSTANTANEOUS MELT PRODUCTION RATE CALCULATED IN EACH OF THE
C PARAMETERISATIONS AND THE MELT FRACTION
C
      parameter (ixzmax=2001)
      character*30 tfilei,tfileo
      double precision pres(ixzmax)
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      double precision vx(ixzmax,ixzmax)
      double precision vz(ixzmax,ixzmax)
      double precision x(ixzmax),z(ixzmax)
      double precision frac(ixzmax,ixzmax)
      double precision temc(ixzmax,ixzmax)
      double precision dpda(ixzmax,ixzmax)
      double precision dtrda(ixzmax,ixzmax)
      double precision dfda(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      double precision dtsda(ixzmax,ixzmax)
      double precision dtuda(2,ixzmax,ixzmax)
C
C COMMUNICATION WITHIN PROGRAM
C
      common /grid/xmax,zmax,dx,dz
      common /grido/pres,x,z,nx,nz
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
      common /mf/p,t,xh2o,cpxm,frac,nfrac
      common /intder/dpdac,dtrdac,f,dfdac
      common /tem/ptem,tfilei,tfileo,item
      common /temo/tem,tems
C TEMCOR=CORRECTED TEMPERATURE
      common /temcor/temc
      common /vel/ufull,alpha,viss,crust,ialpha
      common /velo/vx,vz,uhalf,vmax
C COMMON FOR MELT FRACTION FOR DF/Dt NUMERICAL SOLUTION
      common /melfrc/f1x,f2x,f1z,f2z,fx
C COMMUNICATION OF DERIVATIVES THROUGHOUT THE PROGRAM
      common /derivs/dpda,dtrda,dfda,dtdp,dtdtr,dfdt,dtldp,
     &               dtdtl,dtsdp,dtdts,dfcdp,dfcdtr
      common /dermb/dtsdac,dtudac,dtuda,dtsda
C
C COMMUNICATION OF MELTING RATES 
      common /intgra/tarea,tcrust
      common /deltem/dtem
C      common /deltem/x,xsat,dtem
C
C COMMON - RESIDUE
C
      double precision rest(ixzmax,ixzmax)
      double precision res(3,ixzmax,ixzmax)
      common /resdue/rest,res,resc
C	  
C COMMON (INTERNAL TO INTEG SR)
C
	  double precision xarea(ixzmax)
	  common /xcrust/xarea
C
C 
C SET TOTAL AREA TO ZERO FOR START OF INTEGRATION
C
      tarea=0.
C 
C START SPATIAL LOOPS IN IZ THEN IX, THIS SHOULD DEFINITELY BE
C CORRECTED TEMPERATURE BECAUSE IN MEL SUBROUTINE, TEMPERATURE
C IS CORRECTED FIRST BEFORE MELT FRACTION IS CALCULATED
C
      do iz=1,nz
         do ix=1,nx
            t=temc(ix,iz)
            p=pres(iz)
            dpdac=dpda(ix,iz)
C            write(6,*)'DPDA:',dpdac
            dtrdac=dtrda(ix,iz)
            f=frac(ix,iz)
			if(nfrac(ix,iz).le.0.) then
			 dfda(ix,iz)=0.
			 dfdac=0.
			else
			
C            write(6,*)'BW94_f: ',f
C
C CALCULATE THE VARIOUS DERIVATIVES AND LOCAL INSTANTANEOUS MELT
C PRODUCTION RATE WITH CHOSEN MELTING PARAMETERISATION, FTMP IS
C A TEMPORARY VARIABLE TO ALLOW THE FUNCTION TO RUN, SINCE THE MELT
C FRACTION HAS ALREADY BEEN CALCULATED USING FRCMEL(IP). ALSO,
C BECAUSE THIS IS IN A SEPARATE LOOP TO THE MELT FRACTION CALCULATION
C THE MELT FRACTION FUNCTION WILL NEED TO BE RE-RUN TO CALCULATE 
C VALUES REQUIRED IN THE DERIVATIVE FUNCTION, WHICH ARE SPECIFIC TO
C EACH NODE
C
             if (iparam.eq.1) then
C              write(6,*)'Using mb88, derivatives'
              ftmp=dmb88(iparam)
             elseif (iparam.eq.2) then
C              write(6,*)'Using katz, derivatives'
              ftmp=dk03(iparam)
             else
              write(6,*)'Derivative calculation not implemented'
              stop
             endif
C
C CONDITION THAT IF THE MELTING RATE IS NEGATIVE, MEANS THAT
C THE NODE IS COOLING, SO I'LL EXCLUDE IT FROM THE INTEGRATION
C ALTHOUGH NOT SURE IF THIS IS EXACTLY RIGHT.
C IN BW94, P.443: INTEGRATE OVER THE REGION WHERE dX/dT IS
C GREATER THAN ZERO.
C IF NFRAC IS ZERO THEN NO FURTHER MELTING IS ACHIEVED SO THEREFORE
C THE MELTING RATE MUST BE ZERO.
            if (dfdac.le.0. .or. nfrac(ix,iz).le.0.) then
            dfdac=0.
            nfrac(ix,iz)=0.
            frac(ix,iz)=1.-res(re1,ix,iz)
            endif
C
C
C SAVE MELT RATE IN AN ARRAY
C
            dfda(ix,iz)=dfdac
			endif
C
C WRITE TO SCREEN MELT RATE
C
C            write(6,*)'Melt Rate BW94: ',dfda(ix,iz)
C
C END DERIVATIVES LOOPS
C
         enddo
      enddo 
C
C INTEGRATE OVER THE TOTAL AREA TO GET THE TOTAL VOLUMETRIC 
C MELT PRODUCTION RATE (BOWN & WHITE 94, EQN:4)
C
C FR IS THE VOLUMETRIC PRODUCTION RATE (OR THE START
C OF CALCULATING IT). 
C R=(RHOS/RHOC)*INTEGRAL(DFDA/(1-F**2)
C
C            fr=dfda(ix,iz)/((1.-frac(ix,iz))**2.)
C            if(ix.eq.1) then
C             fp=fr
C              goto 1
C            else if(ix.gt.1) then
C              tarea=tarea+fp+fr
C              write(6,*)'doing tarea calculation'
C              fp=fr
C            endif    
C 1         continue    
C         enddo
C      enddo
C
C INTEGRATE OVER THE TOTAL AREA TO GET THE TOTAL VOLUMETRIC 
C PRODUCTION RATE (BOWN & WHITE 94, EQN 4)
C
C
C FOR EACH X-NODE NEED TO CALCULATE THE AREA OF THE TRAPEZOID
C BETWEEN TWO ADJACENT Z-NODES AND THE AREA OF THE TRAPEZOID 
C BETWEEN THE TWO PREVIOUS Z-NODES
C      do ix=1,nx
C         do iz=1,nz-1
C         fr1=dfda(ix,iz)/((1.-frac(ix,iz))**2.)
C         if (fr1.eq.0.) then
C         area2=0.
C         else
C         fr2=dfda(ix,iz+1)/((1.-frac(ix,iz+1))**2.)
C         area2=(fr1+fr2)*dz/2.
C         endif
C         if(ix.eq.1) then
C         area1=area2
C         goto 1
C         elseif(ix.gt.1) then
C         tarea=tarea+(area1+area2)*dx/2.
C         write(6,*)'Totalarea:',tarea
C         area1=area2
C         endif
C 1       continue
C         enddo
C      enddo
C      do ix=1,nx
C         do iz=1,nz-1
C       do iz=1,nz-1
C          do ix=1,nx
C 1/1/09
C INTEGRATE IN X DIRECTION FIRST
        do iz=1,nz
		  xtarea=0.
		do ix=1,nx-1
C TRY AND PUT FRAC AS NFRAC AND SEE WHAT HAPPENS: 
C DOESN'T MAKE MUCH DIFFERENCE TO FINAL ANSWER
C PLUS THIKNK IT MIGHT MEAN THIS ANYWAY
C 1/1/09 TRY REMOVING EXTRACTION CORRECTION TO SEE IF CAN FIX
C HYDROUS VS ANHYDROUS CRUST AT HIGHER TEMPERATURES....
C         fr1=dfda(ix,iz)/((1.-frac(ix,iz))**2.)
C         fr2=dfda(ix,iz+1)/((1.-frac(ix,iz+1))**2.)
C         fr1=dfda(ix,iz)/((1.-nfrac(ix,iz))**2.)
C         fr2=dfda(ix,iz+1)/((1.-nfrac(ix,iz+1))**2.
C         fr1=dfda(ix,iz)/((res(re1,ix,iz)-nfrac(ix,iz))**2.)
C         fr2=dfda(ix,iz+1)/((res(re1,ix,iz+1)-nfrac(ix,iz+1))**2.)
C 1/1/09 CHANGE INTEGRATION TO SEE IF RID OF HYDROUS CRUST
C LESS THAN ANHYDROUS AT HIGH T
         fr1=dfda(ix,iz)/((res(re1,ix,iz)-nfrac(ix,iz))**2.)
         fr2=dfda(ix+1,iz)/((res(re1,ix+1,iz)-nfrac(ix+1,iz))**2.)
C         fr1=dfda(ix,iz)
C         fr2=dfda(ix,iz+1)
C         area2=(fr1+fr2)*dx/2.
		 xtarea=xtarea+(fr1+fr2)*dx/2.
C         if(ix.eq.1) then
C         area1=area2
C         goto 1
C         elseif(ix.gt.1) then
C 1/1/09         tarea=tarea+(area1+area2)*dx/2.
C         write(6,*)'Totalarea:',tarea
C         area1=area2
C         endif
 1       continue
         enddo
		 xarea(iz)=xtarea
      enddo
C INTEGRATE IN Z DIRECTION
	  do iz=1,nz-1
	     zfr1=xarea(iz)
		 zfr2=xarea(iz+1)
		 tarea=tarea+(zfr1+zfr2)*dz/2.
      enddo
C
C CALCULATE THE FINAL TOTAL...
C RHOS: DENSITY OF MANTLE MATRIX, RHOM: DENSITY OF CRUST
C
       tarea=tarea*rhos/rhom
C      tarea=tarea*dx*rhos/2./rhom
cc       write(6,*)'Total integrated melt production rate', tarea
C
C WRITE OUT THE MELT FRACTION VS DEPTH FOR JUST UNDERNEATH THE RIDGE 
C AXIS, SO TRY NODE IX=1 TO START BUT MIGHT NEED 2.
C
C      open(17,file='bw94_figure7_data')
C      do ix=20,20
C         do iz=1,nz
C         write(17,*)z(iz)*1.e-3,frac(ix,iz)
C         enddo
C      enddo
C      close(17)
C
      return
      end
C
C
C----------------------------------------------------------
C
      function FRCMEL(ip)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
      double precision frac(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
C
C COMMON STATEMENTS
C
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C
C RESIDUE
C
      double precision res(3,ixzmax,ixzmax)
      double precision rest(ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C FUNCTION TO FIND MELT FRACTION FROM VARIOUS PARAMETERISATIONS
C SIGNALLED BY IPARAM
C IPARAM=1, MCKENZIE & BICKLE, 1988
C IPARAM=2, KATZ ET AL., G-CUBED, 4(9), 2003
      if (ip.eq.1) then
C         write(6,*)'Using mb88, frcmel'
         frcmel=xmb88(ip)
      elseif (ip.eq.2) then
         frcmel=xk03(ip)
C 
C         write(6,*) 'frcmel:',frcmel
      else
         write(6,*)'Melting parameterisation not implemented:',iparam
         stop
      endif
C      write(6,*)'frcmel:',frcmel
      return
      end
C
C----------------------------------------------------------
C
      function SOL(ip)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
      double precision frac(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      common /mf/p,t,xh2o,cpxm,frac,nfrac
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
C
C FUNCTION TO FIND MELT FRACTION FROM VARIOUS PARAMETERISATIONS
C SIGNALLED BY IPARAM
C IPARAM=1, MCKENZIE & BICKLE, 1988
C IPARAM=2, KATZ ET AL., G-CUBED, 4(9), 2003
C
      if (ip.eq.1) then
         sol=SMB88(p)
      elseif (ip.eq.2) then
         sol=sk03(p)
      else
         write(6,*)'Melting parameterisation not implemented:',iparam
         stop
      endif
C      write(6,*)'sol:',sol
      return
      end
C
C****************************************************************
C
      function TLHF(tlh,p,ts)
      implicit double precision(a-h,o-z)
C
C CALCULATES THE REAL TEMPERATURE OF A PARTICLE IN THE MELTING
C REGION, ACCOUNTING FOR LOSS OF LATENT HEAT OF FUSION.
C INPUTS:
C  T IS THE TEMPERATURE IN THE ABSENCE OF MELTING
C  P IS THE PRESSURE
C OUTPUT:
C  TLHF IS THE REAL TEMPERATURE
C USES THE PARAMETERISATION DESCRIBED IN WATSON & MCKENZIE (1991).
C CONSTANTS A00-A33 MUST BE DETERMINED BY FITTING THEIR MODEL TO
C ADIABATIC DECOMPRESSION PATHS; A DIFFERENT SET IS NEEDED FOR
C EACH VALUE OF THE ENTROPY OF FUSION.
C
C THE FOLLOWING CONSTANTS FOR ENTROPY OF FUSION OF 400 J/KG/oC
C FROM BOWN DISSERTATION PAGE 56
C
      data a00,a10,a01,a02,a03,a11
     & /0.39653,0.00019,-0.03839,-0.00270,-0.00002,-0.00003/
C
C THE FOLLLOWING CONSTANTS FOR ENTROPY OF FUSION OF 400 J/KG/oC
C FROM BOWN & WHITE 1995 P.18014
C      data a00,a10,a01,a02,a03,a11,a12,a13
C     & /0.398188,0.000122,-0.038345,-0.004680,0.000554,0.000014,
C     & 0.000003,-0.000004/
C
C THE FOLLOWING CONSTANTS FOR ENTROPY OF FUSION OF 250 J/KG/oC
C FROM BOWN DISSERTATION PAGE 56
C
C      data a00,a10,a01,a02,a03,a11
C     & /0.50307,0.00036,-0.03472,-0.00467,-0.00005,-0.00005/
C
      delt=tlh-ts
C FOR USE WITH BOWN DISSERTATION CONSTANTS
      b=a00+a01*p+a02*p**2.+a03*p**3.+a10*delt+a11*delt*p
C FOR USE WITH BOWN & WHITE 1995 CONSTANTS
C      b=a00+a01*p+a02*p**2.+a03*p**3.+a10*delt+a11*delt*p
C     &  +a12*delt*p**2.+a13*delt*p**3.
      tlhf=ts+delt*b
      return
      end
C
C****************************************************************
C
C FUNCTION TO CALCULATE THE REAL TEMPERATURE CHANGE WITH TIME
C FOR STEADY-STATE OR TIME-DEPENDENT CALCULATIONS
C
      function DTDA(itempr)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001) 
      double precision sa(4,4)
      double precision frac(ixzmax,ixzmax)
      double precision dtuda(2,ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      double precision tem(3,ixzmax,ixzmax)
      double precision tems(ixzmax,ixzmax)
      double precision dtsda(ixzmax,ixzmax)
C
CC COMMONS REQUIRED FOR STEADY-STATE
CC
C      common /grid/xmax,zmax,dx,dz
C      common /grido/pres,x,z,nx,nz
C      common /velo/vx,vz,uhalf,vmax
C      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
C      common /temcor/temc
C
C FURTHER COMMONS REQUIRED FOR TIME DEPENDENT
C
      common /temo/tem,tems
      common /dermb/dtsdac,dtudac,dtuda,dtsda
      common /intder/dpdac,dtrdac,f,dfdac
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp      
      common /mf/p,t,xh2o,cpxm,frac,nfrac
C
C THE FOLLOWING CONSTANTS FOR ENTROPY OF FUSION OF 400 J/KG/oC
C FROM BOWN DISSERTATION PAGE 56
C
C      sa(1,1)=0.39653
C      sa(1,2)=-0.03839
C      sa(1,3)=-0.00270
C      sa(1,4)=-0.00002
C      sa(2,1)=0.00019
C      sa(2,2)=-0.00003
C      sa(2,3)=0.
C      sa(2,4)=0.
C
C      efft=t-tsol
C      sum=0.
C      do m=1,2
C         mm=m-1
C         do n=1,4
C            nn=n-1
C            b=sa(m,n)*p**(nn-1)*efft**mm
C            c=p*(mm+1)*(dutemodt-dtsdt)+efft*nn*dpdt
C            sum=sum+b*c
C     write(6,*)'b,c,sum',b,c,sum
C         enddo
C      enddo
C      dtda=dtsda+sum
C      return
C      end
C
C IF ITEMP=0 THEN USE TIME DEPENDENT CALCULATION, METHOD IN BOWN
C AND WHITE 1995
C THE FOLLOWING CONSTANTS FOR ENTROPY OF FUSION OF 400 J/KG/oC
C FROM BOWN DISSERTATION PAGE 56
C
      sa(1,1)=0.39653
      sa(1,2)=-0.03839
      sa(1,3)=-0.00270
      sa(1,4)=-0.00002
      sa(2,1)=0.00019
      sa(2,2)=-0.00003
      sa(2,3)=0.
      sa(2,4)=0.
C 
C THE FOLLOWING CONSTANTS FOR ENTROPY OF FUSION OF 400 J/KG/oC
C FROM BOWN AND WHITE 1995 P. 18014
C      sa(1,1)=0.398188
C      sa(1,2)=-0.038345
C      sa(1,3)=-0.004680
C      sa(1,4)=0.000554
C      sa(2,1)=0.000122
C      sa(2,2)=0.000014
C      sa(2,3)=0.000003
C      sa(2,4)=-0.000004
C
C NEED TO CALCULATE TSOL AGAIN
C
      tsol=SOL(iparam)
      efft=t-tsol
C TSOL PART WORKS
      sum=0.
      do m=1,2
         mm=m-1
         do n=1,4
            nn=n-1
            b=sa(m,n)*p**(nn-1)*efft**mm
            c=p*(mm+1)*(dtudac-dtsdac)+efft*nn*dpdac
            sum=sum+b*c
         enddo
      enddo
      dtda=dtsdac+sum
      return
      end
CC
CC
CC****************************************************************
CC
      function XDS(tem)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
      double precision frac(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      common /mf/p,t,xh2o,cpxm,frac,nfrac
      common /tinit/ti
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C
C RESIDUE COMMUNICATION
C
      double precision rest(ixzmax,ixzmax)
      double precision res(3,ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C FIND MELT FRACTION WHILST CONSERVING ENTHALPY DURING MELTING
C BY ITERATION USING METHOD DESCRIBED IN KATZ ET AL., 2003
C
      t=tem
C      x1=frcmel(iparam)
C RESIDUE INCLUDED IN X1 SO AS TO AFFECT THE TEMPERATURE IN X2
C EQUATION
      x1=frcmel(iparam)-(1.-resc)
C       write(6,*)'x1,fracmel(iparam),resc:',x1,frcmel(iparam),resc
C      write(6,*)'frcmel:',x1
      x2=cp/sy/sy*(ti/t-1.)/ds
C      write(6,*)'ti,t:',ti,t
C      write(6,*)'x2:',x2
      xds=x1-x2
      return
      end
CC****************************************************************
CC
      function XDTEST(tem,itest)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
      double precision frac(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      common /mf/p,t,xh2o,cpxm,frac,nfrac
      common /tinit/ti
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C
C RESIDUE COMMUNICATION
C
      double precision rest(ixzmax,ixzmax)
      double precision res(3,ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C FIND MELT FRACTION WHILST CONSERVING ENTHALPY DURING MELTING
C BY ITERATION USING METHOD DESCRIBED IN KATZ ET AL., 2003
C
      t=tem
C      x1=frcmel(iparam)
C RESIDUE INCLUDED IN X1 SO AS TO AFFECT THE TEMPERATURE IN X2
C EQUATION
      x1=frcmel(iparam)-(1.-resc)
      if(x1.le.0.) then
        itest=1
        xdtest=1.
      else
        itest=0
        xdtest=0.
      endif
        

      return
      end
CC
CC
CC****************************************************************
CC
      function DFXDS(dtrc)
      implicit double precision(a-h,o-z)
      parameter (ixzmax=2001)
      double precision dtida(ixzmax,ixzmax)
      double precision frac(ixzmax,ixzmax)
      double precision nfrac(ixzmax,ixzmax)
      common /derkat/dtidac,dtida,fi
      common /melt/ds,xh2o1,cpxm1,iparam,imelt,itemp
      common /mf/p,t,xh2o,cpxm,frac,nfrac
      common /intder/dpdac,dtrdac,f,dfdac
      common /param/cp,drho,g,pi,rad,sy,rhom,rhos,tcon,tdiff,
     &     texpf,texps
C
C RESIDUE COMMUNICATION
C
      double precision rest(ixzmax,ixzmax)
      double precision res(3,ixzmax,ixzmax)
      common /resdue/rest,res,resc
C
C
C FIND THE MELT PRODUCTION RATE WHILST CONSERVING ENTROPY
C
      dtrdac=dtrc
C      write(6,*)'Using DFXDS dk03'
C      write(6,*)'melt fraction:',f
C      write(6,*)'melt fraction FI:',fi
      dfx1=dk03(ip)
      dfx1=dfdac
C      write(6,*)'dfdac',dfdac
C CHANGED F FOR FI, SINCE THE SMALLER MELTFRACTION IS
C THE ONE FOR WHICH THE IN THE ORIGINAL XDS FUNCTION
C
      dfx2=(cp/sy/sy*dtidac-dtrdac*(cp/sy/sy+fi*ds))/(ds*t)
C      write(6,*)'dtidac:',dtidac
      dfxds=dfx1-dfx2
C      write(6,*)'dfx1,dfx2,dfxds,dtidac:',dfx1,dfx2,dfxds,dtidac
      return
      end
C
C----------------------------------------------------------
C
C FIND ROOT USING BISECTION METHOD FROM NUMERICAL RECIPES
C
      FUNCTION rtbis2(func,x1,x2,xacc)
      implicit double precision(a-h,o-z)
      INTEGER JMAX
      double precision rtbis2,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      double precision dx,f,fmid,xmid
C      write(6,*)'Using rtbis2'
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.)then 
            write(6,*)'rtbis2-no root'
            stop
      endif
      if(f.lt.0.)then
        rtbis2=x1
        dx=x2-x1
      else
        rtbis2=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis2+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis2=xmid
C        write(6,*)'Absolute of dx','xacc',abs(dx),xacc
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      write(6,*)'too many bisections in rtbis2' 
      stop
      END
C
C-------------------------------------------------------------
C FIND ROOT USING BISECTION METHOD FROM NUMERICAL RECIPES
C
      FUNCTION rtbis3(func,x1,x2,xacc)
      implicit double precision(a-h,o-z)
      INTEGER JMAX
      double precision rtbis3,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      double precision dx,f,fmid,xmid
C      write(6,*)'Using rtbis3'
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.)then 
            rtbis3=0.
            return
C            write(6,*)'rtbis2-no root'
C            stop
      endif
      if(f.lt.0.)then
        rtbis3=x1
        dx=x2-x1
      else
        rtbis3=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis3+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis3=xmid
C        write(6,*)'Absolute of dx','xacc',abs(dx),xacc
        if(abs(dx).lt.xacc .or. fmid.eq.0.) then
C        write(6,*)'Found solution to RTBIS3'
        return
        endif 
11    continue
      write(6,*)'too many bisections in rtbis3' 
      stop
      END
C-------------------------------------------------------------
C
C FIND ROOT USING BISECTION METHOD FROM NUMERICAL RECIPES
C
      FUNCTION rtbis4(func,x1,x2,dfacc)
      implicit double precision(a-h,o-z)
      INTEGER JMAX
      double precision rtbis4,x1,x2,dfacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      double precision dx,f,fmid,xmid
C      write(6,*)'Using rtbis4'
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.)then 
            write(6,*)'f,fmid:',f,fmid
            write(6,*)'rtbis4-no root'
            stop
      endif
      if(f.lt.0.)then
        rtbis4=x1
        dx=x2-x1
      else
        rtbis4=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis4+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis4=xmid
C        write(6,*)'Absolute of dx','xacc',abs(dx),xacc
        if(abs(dx).lt.dfacc .or. fmid.eq.0.) return
11    continue
      write(6,*)'too many bisections in rtbis4' 
      stop
      END
C
C-------------------------------------------------------------
