!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Shinibali, Jan 30,2019
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM Vvswm_eval 
use mymodule

IMPLICIT NONE 

 INTEGER                 :: nsite,orbitals,nR,nwan,bq,qq
 INTEGER                 :: nqgrid,nkgrid,nobm,nwtktq=7,Intdim
 integer,allocatable     :: nx(:,:),ny(:,:),Hamind(:,:),chiind(:,:,:,:)
 real*8,allocatable      :: Rvec(:,:), kvec(:,:,:)
 real*8,allocatable      :: a(:),bosonfq(:),Ek(:,:,:),fermi(:,:,:)
 real*8                  :: pi=4.d0*atan(1.d0),MatsuT
 Complex*16,allocatable  :: pmfqbare(:,:,:,:,:),ak(:,:,:,:)
 Complex*16,allocatable  :: bare(:,:,:,:,:,:)
 complex*16              :: jj=cmplx(0.d0,1.d0)
 CHARACTER(LEN=100)      :: arg,outfile
 CHARACTER(*), PARAMETER :: fin='./input_files/',fout='./bareresults/'
 CHARACTER(*), PARAMETER :: fw='NsiteNorbNkxNqxNr.bin',finp='totq_nkx.bin'

allocate(a(nwtktq) )
open(unit=10, file=fin//fw,status='old',access='direct',recl=nwtktq*8)
read(10,rec=1)a
nsite = int(a(1))
orbitals = int(a(2))
MatsuT = a(5)
nR = int(a(7))
 close(10)
deallocate(a)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(a(2))
open(unit=10, file=fin//finp,status='old',access='direct',recl=2*8)
read(10,rec=1)a
nqgrid = int(a(1))     !! This is for the new series of q-points along kpath
nkgrid = int(a(2))     !! This is for the new k-grid
 close(10)
deallocate(a)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nobm=1                  !! Don't change this,it is for the original 0 fq suceptibility calculation
nwan = nsite * orbitals !! length of the Hamiltonian along the row/col
allocate( nx(nkgrid,nqgrid), ny(nkgrid,nqgrid), Rvec(nR,2), kvec(2,nkgrid,nkgrid))
allocate( ak(nwan,nwan,nkgrid,nkgrid), Hamind(orbitals,nsite), chiind(orbitals,nsite,orbitals,nsite) )
allocate( Ek(nwan,nkgrid,nkgrid), fermi(nwan,nkgrid,nkgrid), bosonfq(nobm) )
CALL readinput_lindsus_qpath(orbitals,nsite,nwan,nobm,nkgrid,nqgrid,nx,ny, &
                        & Ek,ak,fermi,bosonfq,Hamind,chiind,nR,Rvec,kvec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Intdim = ( orbitals * nsite )**2
allocate( pmfqbare(Intdim,Intdim,nR,nR,nobm), bare(Intdim,Intdim,nR,nR,nobm,nqgrid) )
do qq = 1,nqgrid
   CALL lindhardsus(orbitals,nsite,nwan,nobm,nkgrid,qq,pmfqbare,nx,ny, &
                    & ak,Ek,fermi,bosonfq,MatsuT,Hamind,chiind,nR,Rvec,kvec)
   bare(:,:,:,:,:,qq) = pmfqbare 
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

outfile = 're_kpath_chi.bin'
open(unit=8,file = fout//trim(outfile) ,status='replace',form='unformatted',access='stream')
write(8)real(bare)
 close(8)
outfile = 'im_kpath_chi.bin'
open(unit=8,file = fout//trim(outfile) ,status='replace',form='unformatted',access='stream')
write(8)aimag(bare)
 close(8)

END PROGRAM