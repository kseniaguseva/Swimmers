This code produces the first part of results of the paper Guseva & Feudel
(2023), "Advantages of run-reverse motility pattern for bacteria in dynamic
fluid environments".

It is written in Python, with additional functions imported from C++ using
Boost.Python. To run the code make sure to have installed Python 3.11 and
Boost.Python.

# Main components of each model:
-------------------------------

## Module_homogeneous:
Models microbial dynamics in homogeneous conditions (no
food sources). 

1. **randomflow_mod.py**: Models a random flow field. For a detailed
description of the implementation of the random flow field and its
properties please see H. Sigurgeirsson and A. M. Stuart, Physics of
Fluids 14, 4352 (2002) or J. García-Ojalvo, J. M. Sancho, and
L. Ramírez-Piscina, Physical Review A 46, 4670 (1992).

2. **swimmers_mod.py**: Implements swimmer motility. Swimmers are
advected and rotated by the flow field, while moving with a certain
constant speed. They are also able to randomly switching direction
after certain periods of time. 

3. **fluid.cc**: functions written in C++, to improve the speed of
computations.

## Module_source:
Models microbial dynamics in the presence of passively advected food sources, we
are interested in the interplay between: flow, swimming, reorientation (motility
modes) and chemotaxis of bacteria.

1. **randomflow_mod.py**: same as above.

2. **swimmers_mod.py**: same as above.

3. **chemotaxis_mod.py**: models extension and contraction of swimming runs
   depending on the changes in concentration along running trajectory.

4. **foodrad_mod**: models the passively advected source of nutrients.

5. **fluid.cc**: functions written in C++, to improve the speed of
computations.


## Additional files to link C++ and Python

**numpy_bind**, **demangle.cc** and **demangle.hh** Functions for moving arrays
between Numpy and C++, requires Boost.Python.


# Run files for different scenarios:
-----------------------------------

The files generate the results of Guseva & Feudel (2023). Create a file
directory "data" and "Fig".

1. **space_swimmers.py**: Generates the data points for Fig.3 (a).
Run using three parameters as input:

    pyhton ./space_swimmers.py us tau alpha

where:
	us  -- velocity of the bacteria relative to the flow velocity (set to 0.5)
	tau -- running time of bacteria between direction changes
	alpha -- shape of the bacteria (set to 0.98 in the figure)
	
the script set ups simulations for all three different swimmer types.

other parameters set to:
      number of swimmers  -- 50000  (N_bac, line 33)
      grid size  -- 250
      dt -- 0.01 s

the output is set after 1 min.
      	       	   
2. **align.py**: Generates the data poins for Fig.3 (b, c, d).
Run using four parameters as input:
    	  pyhton ./align_us.py us tau alpha type
	 
for the definition of parameters see item 1 and also:

       	type -- type of motility pattern (can be "NN" -- simple swimmers,
       	     	"RT" -- swimmers doing run and tumble, "RR" -- swimmers with run
	       	reverse motility)
        us -- swimming velocity of bacteria in relation to the fluid velocity
        tau -- time between reorientation events
        alpha -- shape of bacteria (from 0 to 1; 0 -- spherical, 1 -- ellongated)


        other parameters set to:
        number of swimmers  -- 1000  (N_bac)
      	grid size  -- 250
	    transient time set to -- 1 min (60 s)          
	    output time -- at the end 


to produce the data for the Fig. 3 use:
   (b) run --  submit_us.sh
   (c) run --  submit_tau.sh
   (d) run --  submit_alpha.sh
   
3. **swimmers_fig_motility.py** and **swimmers_fig_shape.py**, produce the Fig.4
   (a, b), please consult the manuscript to set up the parameters. The output is
   saved in the folder Fig. Please note that a series of figures is generated
   after a transient of 40s (see lines 81 and 82 to set up the output).
   
4. **dist_source.py** used to produce the data for Fig. 5, please consult the manuscript 
for parameters. The output is saved in the folder "data".
   

   


