
C Copyright (c) 2010,2011,2012 RENCI.
C All rights reserved. This program and the accompanying materials
C MAY BE available under the terms of the RENCI Open Source License
C UNC at Chapel Hill which accompanies this distribution, and is available at
C http://www.renci.org/resources/open-source-software-license

C
C wrap IO calls to SIFS from C++
C MO integrals

       SUBROUTINE FTOCXXOPEN(iunit, filename, sizename, istat )
       CHARACTER*256 filename 
       INTEGER sizename
       INTEGER istat
       CHARACTER*256 FNAME
       INTEGER iunit
   
       FNAME(1:sizename) = filename 
       istat = 0
       OPEN(iunit,FILE=FNAME,STATUS='OLD',FORM='UNFORMATTED',ERR=999)
       RETURN
 999   istat = 1
       RETURN 
       END

C Formatted Open

       SUBROUTINE FTOCXXFORMATTEDOPEN(iunit, filename, sizename, istat )
       CHARACTER*256 filename
       INTEGER sizename
       INTEGER istat
       CHARACTER*256 FNAME
       INTEGER iunit
       istat = 0
       FNAME(1:sizename) = filename
       WRITE(6,*)"Open for MOs"
       WRITE(6,*) sizename
       WRITE(6,*) FNAME(1:sizename)
       WRITE(6,*) iunit
       OPEN(iunit,FILE=FNAME, STATUS='OLD', FORM='FORMATTED',ERR=999)
       WRITE(6,*) "Finished with MO fto open"
       RETURN
 999   istat = 1 
       RETURN
       END

C
C Formatted Open NEW : fail if file exists

       SUBROUTINE FTOCXXFORMATTEDNEWOPEN(iunit,filename, sizename,istat)
       CHARACTER*256 filename
       INTEGER sizename
       INTEGER istat
       CHARACTER*256 FNAME
       INTEGER iunit
       istat = 0
       FNAME(1:sizename) = filename
       WRITE(6,*)"Open for NOs"
       WRITE(6,*) sizename
       WRITE(6,*) FNAME(1:sizename)
       WRITE(6,*) iunit
       OPEN(iunit,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED',ERR=999)
       RETURN
 999   istat = 1
       RETURN
       END

C
       SUBROUTINE FTOCXXCLOSE(iunit)
       INTEGER iunit

       CLOSE(iunit)
       RETURN
       END

