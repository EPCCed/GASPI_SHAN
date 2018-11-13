# miniGhost: stencil computations with boundary exchange.

This version of the MiniGhost benchmark from Sandia has been modified.  
It is used to demonstrate a transition from flat MPI code to GASPI-SHAN.  
(GASPI Shared Notifications). GASPI-SHAN directly uses shared memory for  
both data layout and solver data. Node local sends here are replaced with  
notifications in shared memory. Communication to remote nodes is replaced  
with notified GASPI communication.







