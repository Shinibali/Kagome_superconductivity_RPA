!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Function: collecting bare susceptibility data, because f90 will be faster than matlab on cluster
!! Shinibali, NOV 26,2016
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM collectdata 
IMPLICIT NONE 

 INTEGER                 :: i,nwtktq=7,nsite,orbitals,nR,nobm,nq,nk,totq,Intdim,nelem
 real*8,allocatable      :: a(:),chire(:),chiim(:)
 Complex*16,allocatable  :: chio(:,:),rearrchio(:,:,:)
 CHARACTER(LEN=5)        :: str 
 CHARACTER(LEN=10000)    :: outfile
 CHARACTER(*), PARAMETER :: fin='./input_files/',fout='./bareresults/',fw='NsiteNorbNkxNqxNr.bin'
 CHARACTER(*), PARAMETER :: fchiore='chio_re.bin',fchioim='chio_im.bin'
 CHARACTER(*), PARAMETER :: fnx = 'nx.bin', fny = 'ny.bin'

allocate(a(nwtktq))
open(unit=10, file=fin//fw,status='old',access='direct',recl=nwtktq*8)
read(10,rec=1)a
nsite = int(a(1))
orbitals = int(a(2))
nk = int(a(3))
nq = int(a(4))
totq=nq**2
nR = int(a(7))
 close(10)
nobm = 1
Intdim = ( orbitals * nsite )**2
nelem = Intdim*Intdim*nR*nR*nobm

open(unit=8,file = fin//fnx ,status='old',access='direct',recl=totq*nk*8)
 close(8,status='delete')
open(unit=8,file = fin//fny ,status='old',access='direct',recl=totq*nk*8)
 close(8,status='delete')

allocate( chire(nelem), chiim(nelem), chio(nelem,totq) )

DO i = 1,totq
  
  write(str,"(I5.1)") i
  
  outfile ='Re_chio_qq'//trim(adjustl(str))//'.bin'
  open(unit=8,file = fout//trim(outfile) ,status='old',access='direct',recl=nelem*8)
  read(8,rec=1)chire
  close(8,status='delete')
    
  outfile ='Im_chio_qq'//trim(adjustl(str))//'.bin'
  open(unit=8,file = fout//trim(outfile) ,status='old',access='direct',recl=nelem*8)
  read(8,rec=1)chiim
  close(8,status='delete')

  chio(:,i) = cmplx(chire,chiim)
  
ENDDO  
 rearrchio = reshape( chio,(/ nelem,nq,nq /) )

  open(unit=8,file = fout//fchiore ,status='replace',form='unformatted',access='stream')
  write(8)real(rearrchio)
  close(8)
  
  open(unit=8,file = fout//fchioim ,status='replace',form='unformatted',access='stream')
  write(8)aimag(rearrchio)
  close(8)

END PROGRAM
