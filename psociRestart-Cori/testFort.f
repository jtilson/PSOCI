      PROGRAM READINTS

      CHARACTER*256 A,B,C
      CHARACTER*64 FILENAME
      INTEGER*4 ntitle,nsym,nbft,ninfo,nenrgy,nmap,ierr

c hardcode in some arrays
      CHARACTER*80 title(20)
      INTEGER*4 nbpsy(4)
      CHARACTER*4 slabel(4)
      INTEGER*4 info(5)
      CHARACTER*8 bfnlab(86)
      INTEGER*4 ietype(1)
      REAL*8 energy(1)
      INTEGER*4 imtype(0)
      INTEGER*4 map(86,0)

      FILENAME='./INTEGRALS/RUO32BIT/moints'
      LUNIT=2

      OPEN(LUNIT,FILE=FILENAME,STATUS='OLD', FORM='UNFORMATTED')
      
      call sifrh1(LUNIT,ntitle,nsym,nbft,ninfo,nenrgy,nmap,ierr)

c
c -> This works fine
c
      WRITE(6,*)" Write parameters "
      WRITE(6,*)'ntitle,nsym,nbft,ninfo,nenrgy,nmap,ierr',
     +  ntitle,nsym,nbft,ninfo,nenrgy,nmap,ierr
c
c-> Try second header read
c
      call sifrh2(LUNIT,ntitle,nsym,nbft,ninfo,nenrgy,nmap,
     &  title, nbpsy,   slabel,  info,    bfnlab,
     &  ietype,  energy,  imtype,  map, ierr )

      WRITE(6,*)' AFTER sifrh2 '
      WRITE(6,*)' TWO TITLES are ', title(1)
      WRITE(6,*)' 2ND TITLE is ', title(2)
      WRITE(6,*)' 3 TITLE is ', title(3)
      WRITE(6,*)' 4 TITLE is ', title(4)
      WRITE(6,*)' 5 TITLE is ', title(5)
      WRITE(6,*)' 6 TITLE is ', title(6)
      WRITE(6,*)' 7 TITLE is ', title(7)

      WRITE(6,*)' NBPSY 1', nbpsy(1)
      WRITE(6,*)' NBPSY 2', nbpsy(2)
      WRITE(6,*)' NBPSY 3', nbpsy(3)
      WRITE(6,*)' NBPSY 4', nbpsy(4)

      WRITE(6,*)' SYMLABEL 1 ', slabel(1) 
      WRITE(6,*)' SYMLABEL 2 ', slabel(2) 
      WRITE(6,*)' SYMLABEL 3 ', slabel(3) 
      WRITE(6,*)' SYMLABEL 4 ', slabel(4) 

      WRITE(6,*)' INFO 1', info(1)
      WRITE(6,*)' INFO 2', info(2)
      WRITE(6,*)' INFO 3', info(3)
      WRITE(6,*)' INFO 4', info(4)
      WRITE(6,*)' INFO 5', info(5)

      WRITE(6,*)' bfnlab 1', bfnlab(1)
      WRITE(6,*)' bfnlab 2', bfnlab(2)
      WRITE(6,*)' bfnlab 3', bfnlab(3)
      WRITE(6,*)' bfnlab 4', bfnlab(4)
      WRITE(6,*)' bfnlab 86', bfnlab(86)

      WRITE(6,*)' map 1', map(1,1)
      WRITE(6,*)' map 2', map(2,1)
      WRITE(6,*)' map 3', map(3,1)
      WRITE(6,*)' map 4', map(4,1)
      WRITE(6,*)' map 86', map(86,1)

      WRITE(6,*)' IETYPE ', ietype(1)
      WRITE(6,*)' ENERGY terms are', nenrgy, energy(1)

c
c-> test reading MO coefficient file for CIDEN code
c
      END     

