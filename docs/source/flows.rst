Flows
=====

The flows structure (or function) take care of handling the tasks.


Tasks
~~~~~

Tasks contain interdependent works.

YamboQPBSETasks
---------------

This task run a GW and Bethe-Salpeter calculation.

YamboQPTask
-----------


In addition to the Yambo-related tasks, yambopy has also pw-related tasks to perform ca
non-selfconsistent calculations, band structure calculations and cell optimization. 

PwRelaxTasks
------------

The relaxation task performs three concatenated calculations. First, the atomic relaxation is performed. The second calculation reads the
new atomic positions and it performs a cell relaxation. The third and last calculation is the self-consistent calculations of the density
with the optimized cell parameters and atomic positions.

atomic relaxation >> cell relaxation >> self-consistent calculation



PwNscfTasks
-----------




PwBandsTasks
------------
