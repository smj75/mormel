# mormel
**Time-dependent mid-ocean ridge melting calculator.**

MORMEL calculates the degree and composition of decompressional mantle melting when the geometry of mantle motion and/or the temperature of the mantle source can vary through time.  The code was designed to calculate mid-ocean ridge melting, and it has also been used to estimate melt producitivity in other decompressional melting situations.

The following studies have used MORMEL.  The code version in this repository in Jea14.

- [W10:](http://www.tara.tcd.ie/handle/2262/78244)  Walters RL.  Geochemical signature of rift relocations at Iceland.  PhD Thesis, University of Dublin, Trinity College.  *Describes code development in Chapter 6 and some applications to modelling rift relocations at Iceland in Chapter 7.  This original version had a simplified method for calculating trace element compositions that is superceded by the versions below.*

- [Wea13:](http://doi.org/10.1016/j.epsl.2013.06.040)  Walters RL, Jones SM, Maclennan J.  Renewed melting at the abandoned Húnafloí Rift, northern Iceland, caused by plume pulsing.  Earth and Planetary Science Letters 377–378 (2013) 227–238, 10.1016/j.epsl.2013.06.040.  *The first peer-reviewed publication to use MORMEL.  Model description is in §5.3 of the main paper and §B of the supplementary material.  The calculations estimate the amount, timing and composition of magma produced during reactivation of an abandoned spreading axis.  The compositional calculation works in combination with Dan McKenzie's INVMEL code.  Note that the version of MORMEL in this repository is Jea14, which has an improved and more convenient composition calculation.*  

- [Jea14:](http://doi.org/10.1016/j.epsl.2013.09.029) Jones SM, Murton BJ, Fitton JG, White NJ, Maclennan J, Walters RL.  A joint geochemical-geophysical record of time-dependent mantle convection south of Iceland.  Earth and Planetary Science Letters 386 (2014) 86–97, 10.1016/j.epsl.2013.09.029.  *The MORMEL code version in this repository.  Augments the Wea13 version to provide a more self-consistent estimate of the trace element composition that is delivered entirely by MORMEL.  Model description is in §6 of the main paper and §B of the supplementary material.  The calculations estimate the amount, timing and composition of magma produced when pulses of hotter and cooler mantle travel beneath a mid-ocean ridge.*  

- [Gea23:](http://doi.org/10.21203/rs.3.rs-986686/v1)  Gernon TM, Jones SM, Brune S, Hincks TK, Merdith AS, Palmer MR, Schumacher JC, Primiceri RM, Field M, Griffin WL, O’Reilly SY, Keir D, Spencer CJ, Glerum A.  Kimberlite ascent by rift-driven disruption of cratonic mantle keels.  Accepted in principle  by Nature (April 2023).  *The Jea14 code was used to make Figure 3c.  The calculations estimate the magma productivity and magma water content associated with delamination of basal cratonic lithosphere.*

 

Make using "make". For examples of operation, see the shell and GMT (Generic Mapping Tools) scripts below that were used to create some key figures in the references above.  

## *mormel.f*, *lib_comp.f*, *lib_k03.f*, *lib_mb88.f*, *lib_mel.f*, *lib_res.f*, *lib_comp.f*, *lib_td.f*, *lib_veltem.f*

Source code for MORMEL.

## *wj_parameters*

Example input file for MORMEL.  The comments in this file are not intended as a comprehensive set of instructions, and are not guaranteed to be accurate because they might be relict from earlier versions of the code.  To understand the functionality and be sure that the code is doing what you want, you should check the source code files for yourself and carefully inspect the outputs.  You can also refer to the scripts below that were used to generate the some of the figures in the peer-reviewed papers.

## *deep_hydrous_gea23.gmt*, *make_flow_models_gea23.csh*, *run_steady_state_auto_gea23.csh*

Scripts used to calculate and plot Figure 3c of Gernon, Jones et al. (2023).  The calculations estimate the magma productivity and magma water content associated with delamination of basal cratonic lithosphere.  See the paper itself and the comments in the scripts for more detail.  


