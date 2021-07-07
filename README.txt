purpose: generate pure hexahedral mesh for random pebble bed in a cylinder container

input: random pebble bed in a cylinder container
output: pure quadratic hexahedral mesh in Nek5000 format

strategy: 
1. use the voronoi cell method to divide domain to into voronoi cells. Each voronoi cell contains one pebble.
2. the polygon faces of voronoi cell will be divided into quads. 
3. quads will be projected to pebble surface, while generating hexahedral elements


usage:

1. obtain random pebble coordinates from experiment or DEM calculation.
use the pebbles.dat file to store pebble coordinates.
pebbles.dat file also contains infomation about the cylinder.

2. use top_ghost.dat and bot_ghost.dat to store the pebbles just above and below pebbles in pebbles.dat.
pebbles in top_ghost.dat and bot_ghost.dat will be used to construct ghost pebbles.

To use top_ghost.dat and bot_ghost.dat is optional but prefered.
because use top_ghost.dat and bot_ghost.dat could reduce the max aspect ratio of elements.

NOTICE: the final mesh will be scaled to ensure pebble diameter is 1. 
But user should provide dimensional values in pebbles.dat, top_ghost.dat and bot_ghost.dat. 

3. run voro_to_hex_full_parallel.py to generate hfiles

if_top_bot_ghost_pebble_from_files is True if top_ghost.dat and bot_ghost.dat are provided.
otherwise if_top_bot_ghost_pebble_from_files should be False.

nprocs is the number of parallel processors you want to use.

In this step, you should watch if code report non-right-hand elements.
if non-right-hand elements are reported, then reduce the max value tol array from 0.12 to 0.1.

if unlucky, several iterations of tunning may needed in this step.

4. run hfiles_to_rea.py to gather parallel data and convert to rea files
nprocs should be consistent  with the value in voro_to_hex_full_parallel.py.
newReaFile is the exported rea file.

5. run Nek5000 tool reatore2 and gencon
reatore2 will convert rea mesh into re2 mesh.
gecon will generate co2 file from re2 file.
if this step is an iteration to fix negative and low jacobian elements, user should NOT rerun gencon.

6. run Nek5000/NekRS case. 
if negative jacobian elements are reported, store these element numbers in the nje array in hfiles_to_rea.py and redo step 4.

7. put the following code in userchk() to dump element with low scaled-jacobian.

ccc
     if (istep.eq.0) then	
c print out minimum scale-jac element number
c follow subroutine  mesh_metrics
   
      do ie = 1,nelt
         dratio = vlmin(JACM1(1,1,1,ie),nxyz)/
     $            vlmax(JACM1(1,1,1,ie),nxyz)
        if (dratio.le.7e-3) then
         eg = lglel(ie)
         write(6,*) 'please linearize element: ',eg
        endif
      enddo   

      endif
ccc

Store these element numbers in the lje array in hfiles_to_rea.py and redo step 4.

Several iterations may needed from step4 to 7 to ensure there is no negative jacobian element and minimum scaled jacobian is acceptable. 
NOTICE rerun gencon is not needed !!!!

