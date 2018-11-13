# SHAN - (SHA)red (N)otification library.
The SHAN library makes use of both shared memory communication and the
GASPI communication library in order to provide an efficient interface
for a persistant neighborhood communication (e.g ghost cell exchanges).  

SHAN provides interfaces for

- memory allocation.  
  As the SHAN library directly uses shared memory, solver data which
  is being exchanged, has to exist in shared memory too.
  Applications which want to use SHAN need to replace existing memory allocation with
  a call to 'shan_alloc_shared' with 'SHAN_DATA' as allocation type.
  The pointer to the allocated memory for solver data can be accessed with
  'shan_get_shared_ptr' (see SHAN_comm.h).  
  
- neighborhood initialization  
  The SHAN lib establishs a persistant communication
  between neighbors with the call to 'shan_comm_init_comm'. 

- persistant communication
  SHAN uses a flexible type concept, where type information is published
  in shared memory. Local neighbors can access that information
  and directly convert remote data types into local data types without the intermediate
  step of writing or reading a linear buffer. 
  SHAN pre-allocates memory for the requested communication in the form of a so called 'SHAN type'.
  This includes buffer space for data types, and the communication ressources for 
  intra (via shared windows) and inter-node (via GASPI) messages.    
  A SHAN type can be seen as an independent persistant communication handle
  for both shared and remote memory. 
  
- setting communication meta data  
  SHAN users need to provide  for the numbers of elements, the element size
  and the element offsets both for receiving and sending. 
  Pointers to this meta data can be accessed via 'shan_comm_type_offset'
  The length of the actual message and the offsets can be 
  adjusted dynamically by changing these meta data values.

- writing of data  
  Node local communication will use reading rather than writing.
  As both type information and data is visble across the node,
  the SHAN lib direcly can access that data. A write then merely flags
  that data is available for reading.

- receiving/waiting and testing for data.  
  the SHAN library will directly convert types in the receive step
  if data is shared node locally or unpack if data is
  being sent from other nodes.

- waiting for sends.
  As there is no sending of data node-locally (but rather a shared memory notification)
  waiting for send requests actually is replaced by the wait for 'all other ranks have
  read the data'. For remote messages SHAN makes use of double buffering.
  Local buffers here can be reused if the remote data has arrived.
  This implicitly provided validity of remote buffers however is only valid for bidirectional communication.

  