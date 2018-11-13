/*
    Copyright (c) T-Systems SfR, C.Simmendinger <christian.simmendinger@t-systems.com>, 2018

    This file is part of SHAN.

    SHAN is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
	    the Free Software Foundation, either version 3 of the License, or
	        (at your option) any later version.

    SHAN is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
	    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	        GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
        along with SHAN.  If not, see <https://www.gnu.org/licenses/>.
	
*/


#include <mpi.h>
#include <unistd.h>
#include <stdio.h>
#include "SHAN_comm.h"

#define DATA_KEY 4710

void local_exchange_comm(int iProc
			 , int target
			 , int sz
			 , void *sbuf
			 , void *rbuf
			 , MPI_Comm MPI_COMM_ALL
			 )
{
  /* mutual exchange of offset array */
  if (target > iProc) /* first send */
    {
      MPI_Send(sbuf
	       , sz
	       , MPI_CHAR
	       , target
	       , DATA_KEY
	       , MPI_COMM_ALL
	       );
      MPI_Recv(rbuf
	       , sz
	       , MPI_CHAR
	       , target
	       , DATA_KEY
	       , MPI_COMM_ALL
	       , MPI_STATUS_IGNORE
	       );
    }
  else  /* first receive */
    {
      MPI_Recv(rbuf
	       , sz
	       , MPI_CHAR
	       , target
	       , DATA_KEY
	       , MPI_COMM_ALL
	       , MPI_STATUS_IGNORE
	       );
      MPI_Send(sbuf
	       , sz
	       , MPI_CHAR
	       , target
	       , DATA_KEY
	       , MPI_COMM_ALL
	       );
    }
}

