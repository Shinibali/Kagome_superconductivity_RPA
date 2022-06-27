!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Shinibali, Jan 30,2019
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM susceptibility 
use mymodule

IMPLICIT NONE 

 INTEGER                 :: INP,opt1,opt2,I,nsite,orbitals,nwan,qloop,qq,uu,bq
 INTEGER                 :: nqgrid,nkgrid,nR,Intdim,nobm=1,nwtktq=7
 integer,allocatable     :: iq(:),nx(:,:),ny(:,:),Hamind(:,:),chiind(:,:,:,:)
 real*8,allocatable      :: a(:),bosonfq(:),Ek(:,:,:),fermi(:,:,:),Rvec(:,:),kvec(:,:,:)
 real*8                  :: pi=4.d0*atan(1.d0),MatsuT
 Complex*16,allocatable  :: pmfqbare(:,:,:,:,:),ak(:,:,:,:)
 complex*16              :: jj=cmplx(0.d0,1.d0)
 CHARACTER(LEN=10)       :: str1,str2
 CHARACTER(LEN=10000)    :: arg,outfile
 CHARACTER(*), PARAMETER :: fin='./input_files/', fout='./bareresults/'
 CHARACTER(*), PARAMETER :: fw='NsiteNorbNkxNqxNr.bin'

INP = IARGC()
call GETARG(1,arg)
read(arg,*) opt1
call GETARG(2,arg)
read(arg,*) opt2
iq =(/ (I,I=opt1,opt2) /)

allocate(a(nwtktq) )
open(unit=10, file=fin//fw,status='old',access='direct',recl=nwtktq*8)
read(10,rec=1)a
nsite = int(a(1))
orbitals = int(a(2))
nkgrid = int(a(3))
nqgrid = int(a(4))
MatsuT = a(5)
nR = int(a(7))
 close(10)
deallocate(a)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nwan = nsite * orbitals !! length of the Hamiltonian along the row/col
allocate( nx(nkgrid,nqgrid**2), ny(nkgrid,nqgrid**2), Rvec(nR,2), kvec(2,nkgrid,nkgrid))
allocate( ak(nwan,nwan,nkgrid,nkgrid), Hamind(orbitals,nsite), chiind(orbitals,nsite,orbitals,nsite) )
allocate( Ek(nwan,nkgrid,nkgrid), fermi(nwan,nkgrid,nkgrid), bosonfq(nobm) )
CALL readinput_lindsus(orbitals,nsite,nwan,nobm,nkgrid,nqgrid,nx,ny, &
                        & Ek,ak,fermi,bosonfq,Hamind,chiind,nR,Rvec,kvec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Intdim = ( orbitals * nsite )**2
allocate( pmfqbare(Intdim,Intdim,nR,nR,nobm))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do qloop = lbound(iq,1),ubound(iq,1)
  qq = iq(qloop)

  CALL lindhardsus(orbitals,nsite,nwan,nobm,nkgrid,qq,pmfqbare,nx,ny, &
                    & ak,Ek,fermi,bosonfq,MatsuT,Hamind,chiind,nR,Rvec,kvec)

  write(str1,"(I10.1)") qq

  outfile ='Re_chio_qq'//trim(adjustl(str1))//'.bin'
  open(unit=8,file = fout//trim(outfile) ,status='replace',form='unformatted',access='stream')
  write(8)real(pmfqbare)
  close(8)
  
  outfile ='Im_chio_qq'//trim(adjustl(str1))//'.bin'
  open(unit=8,file = fout//trim(outfile) ,status='replace',form='unformatted',access='stream')
  write(8)aimag(pmfqbare)
  close(8)

enddo  !! qloop

END PROGRAM
