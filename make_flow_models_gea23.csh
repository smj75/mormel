#!/bin/csh

# RUN SET OF MORMEL MODELS WITH VARYING WATER CONTENT
# FOR GERNON, JONES ET AL KIMBERLITE PAPER
# SCRIPT ADAPTED FROM make_ce0806_models.csh THAT
# PRODUCES SET OF 4 VSR MODELS USED IN CE0806 PAPER

# THIS SCRIPT DRIVES "run_steady_state_auto_gea23.csh"
# WHICH GENERATES THE MORMEL INPUT PARAMETER FILE
# FOR INDIVIDUAL MODEL RUNS

set MAKE_MELTING_MODELS=1
set MAKE_COMPOSITION_GRIDS=0
set RUN_VSR_COMPOSITION_MODELS=0

# POTENTIAL TEMPERATURE IS FIXED IN THIS STUDY
# WE ARE INTERESTED IN NORMAL SITUTATION - 
# PRESENCE OF A HOT PLUME WOULD GENERATE SIILL MORE MELT

set T_PS="1300"

# WEDGE ANGLES - FIXED
# UNIMPORTANT IN ITSELF IN THIS CASE
# BECAUSE WE ARE ONLY INTERESTED IN THE SINGLE
# MELTING COLUMN BENEATH THE AXIS OF UPWELLING. 
# IN PRACTICE, THE VALUE IS CHOSEN TOGETHER
# WITH UFULL IN THE MORMEL INPUT FILE TO 
# ENSURE THE REQUIRED VALUE FOR UPWELLING AT
# THE CENTRE OF THE MODEL.   
# IN THIS CASE, ANGLE=80 DEG, UFULL=230 KM/MYR 
# GIVES REQUIRED UPWELLING RATE OF 30 KM/MYR

set ANGLES="80"

# TIME STEP

set TIME_STEP=0.02

# BULK WATER CONTENTS
# A MAIN REQUIREMENT IN THIS STUDY IS TO TEST
# HOW MUCH WATER IS REQUIRED TO GENERATE MELTING
# WITH NORMAL TEMPERATURE MANTLE BENEATH VERY THICK
# LITHOSPHERE

set XH2OS="0.0 0.05 0.1 0.15 0.2"
#set XH2OS="0.15"

#
# MAKE STEADY-STATE MELTING MODELS USED TO CALCULATE THE COMPOSITION GRIDS
#

if ( $MAKE_MELTING_MODELS == 1 ) then
	
	set ANGLE=$ANGLES
	set STORE=.
	set END_TIME=1.0
	set TP=$T_PS

	foreach XH2O ($XH2OS)
		set MODEL=${TP}_${ANGLE}_${TIME_STEP}_${XH2O}
		echo "Calculating steady-state melting model $MODEL..."

# MAKE DIRECTORY

		set DIR=${STORE}/${TP}_${ANGLE}
		set DIR=${STORE}/${TP}_${ANGLE}_${TIME_STEP}
		set DIR=${STORE}/${TP}_${ANGLE}_${TIME_STEP}_${XH2O}
		if ( -d $DIR ) then 
    		echo "$DIR exists already" 
    	else
     	 mkdir $DIR
     	 echo "Made directory $DIR" 
    	endif
    
#

		set RUN_MODEL = 1
		while ($RUN_MODEL == 1)

# CRUSTAL THICKNESS

			if ( -f $DIR/WJ_crust_steady ) then 
    			set TMP=`tail -1 $DIR/WJ_crust_steady`
    			set CRUST=$TMP[4]
    		else
    			set CRUST=7.78
    		endif
    		echo "Crustal thickness: $CRUST" 

# PARAMETERS FILE

			date
			source ./run_steady_state_auto_gea23.csh
			date

# REPEAT WITH CORRECT CRUSTAL THICKNESS?

    		set TMP=`tail -1 $DIR/WJ_crust_steady`
    		set RUN_MODEL=`awk -v c0=$CRUST -v c1=$TMP[4] 'BEGIN{dif=c0-c1; dif2=dif*dif; if (dif2>0) dif=sqrt(dif2); if (dif<0.2) {print 0} else {print 1}}'`
    		set CRUST=$TMP[4]

    		set RUN_MODEL=0

		end
	end
endif



