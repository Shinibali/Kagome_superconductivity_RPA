SUBROUTINE readinput_lindsus(orbitals,nsite,nwan,nobm,nkx,nqx,nx,ny, &
                             & Ek,ak,fermi,bosonfq,Hamind,chiind,nR,Rvec,kvec)

 IMPLICIT NONE

 INTEGER                 :: orbitals,nsite,nwan,nkx,nqx,nR,totq,dimen,totk,nobm
 integer,allocatable     :: nx(:,:),ny(:,:), Hamind(:,:), chiind(:,:,:,:)
 real*8,allocatable      :: a(:),Rvec(:,:),kvec(:,:,:),Ek(:,:,:),fermi(:,:,:)
 real*8,allocatable      :: akre(:,:,:,:),akim(:,:,:,:),bosonfq(:)
 Complex*16,allocatable  :: ak(:,:,:,:)
 CHARACTER(*), PARAMETER :: fin="./input_files/",fout="./bareresults/",fRvec='Rvec.bin',fkv='kvec.bin'
 CHARACTER(*), PARAMETER :: fek = 'ek.bin', ffe = 'fermi.bin', fakre = 'akre.bin', fakim = 'akim.bin'
 CHARACTER(*), PARAMETER :: fnx = 'nx.bin', fny = 'ny.bin', fHi = 'Hamind.bin', fIi = 'chiind.bin'
 logical                 :: file_exists

totq = nqx**2
totk = (nkx)**2 !!Sum over all K-points in the grid
dimen = nwan*nwan

bosonfq(nobm) = 10.d0 ** (-5)

allocate( a((nkx)*totq) )
INQUIRE(FILE=fin//fnx, EXIST=file_exists)
if ( file_exists ) then
  open(unit=10, file=fin//fnx,status='old',access='direct',recl=(nkx)*totq*8)
  read(10,rec=1)a
  nx = reshape( int(a),(/ nkx,totq /) )
  close(10)
  open(unit=10, file=fin//fny,status='old',access='direct',recl=(nkx)*totq*8)
  read(10,rec=1)a
  ny = reshape( int(a),(/ nkx,totq /) )
  close(10)
else
  a=0
  nx = reshape( int(a),(/ nkx,totq /) )
  ny = reshape( int(a),(/ nkx,totq /) )
end if
deallocate(a)

allocate( a(nwan*nkx*nkx) )
open(unit=10, file=fin//fek,status='old',access='direct',recl=nwan*nkx*nkx*8)
read(10,rec=1)a
Ek = reshape( a,(/ nwan,nkx,nkx /) )
 close(10)
open(unit=10, file=fin//ffe,status='old',access='direct',recl=nwan*nkx*nkx*8)
read(10,rec=1)a
fermi = reshape( a,(/ nwan,nkx,nkx /) )
 close(10)
deallocate(a)

allocate( a(nwan*nkx*nkx*nwan))
open(unit=10, file=fin//fakre,status='old',access='direct',recl=nwan*nkx*nkx*nwan*8)
read(10,rec=1)a
akre = reshape( a,(/ nwan,nwan,nkx,nkx /) )
 close(10)
open(unit=10, file=fin//fakim,status='old',access='direct',recl=nwan*nkx*nkx*nwan*8)
read(10,rec=1)a
akim = reshape( a,(/ nwan,nwan,nkx,nkx /) )
 close(10)
ak = cmplx(akre,akim)
deallocate(a)

allocate( a(nwan))
open(unit=10, file=fin//fHi,status='old',access='direct',recl=nwan*8)
read(10,rec=1)a
Hamind = reshape( int(a),(/ orbitals,nsite /) )
  close(10)
deallocate(a)

allocate( a(orbitals**2 * nsite**2))
open(unit=10, file=fin//fIi,status='old',access='direct',recl=orbitals**2 * nsite**2 *8)
read(10,rec=1)a
chiind = reshape( int(a),(/ orbitals,nsite,orbitals,nsite /) )
  close(10)
deallocate(a)

allocate( a(nR*2) )
open(unit=10, file=fin//fRvec,status='old',access='direct',recl=nR*2*8)
read(10,rec=1)a
Rvec = reshape( a,(/ nR,2 /) )
 close(10)
deallocate(a)

allocate( a(2*nkx*nkx) )
open(unit=10, file=fin//fkv,status='old',access='direct',recl=2*nkx*nkx*8)
read(10,rec=1)a
kvec = reshape( a,(/ 2,nkx,nkx /) )
 close(10)

 RETURN
END subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE readinput_lindsus_qpath(orbitals,nsite,nwan,nobm,nkx,nqx,nx,ny, &
                             & Ek,ak,fermi,bosonfq,Hamind,chiind,nR,Rvec,kvec)

 IMPLICIT NONE

 INTEGER                 :: orbitals,nsite,nwan,nkx,nqx,nR,totq,dimen,totk,nobm
 integer,allocatable     :: nx(:,:),ny(:,:), Hamind(:,:), chiind(:,:,:,:)
 real*8,allocatable      :: a(:),Rvec(:,:),kvec(:,:,:),Ek(:,:,:),fermi(:,:,:)
 real*8,allocatable      :: akre(:,:,:,:),akim(:,:,:,:),bosonfq(:)
 Complex*16,allocatable  :: ak(:,:,:,:)
 CHARACTER(*), PARAMETER :: fin="./input_files/",fout="./bareresults/",fRvec='Rvec.bin'
 CHARACTER(*), PARAMETER :: fkv='kvec_wm.bin', fek = 'ek_wm.bin', ffe = 'fermi_wm.bin'
 CHARACTER(*), PARAMETER :: fakre = 'akre_wm.bin', fakim = 'akim_wm.bin',fHi = 'Hamind.bin'
 CHARACTER(*), PARAMETER :: fnx = 'nx_wm.bin', fny = 'ny_wm.bin', fIi = 'chiind.bin'
 logical                 :: file_exists

totq = nqx
totk = (nkx)**2 !!Sum over all K-points in the grid
dimen = nwan*nwan

bosonfq(1) = 0.0

allocate( a((nkx)*totq) )
INQUIRE(FILE=fin//fnx, EXIST=file_exists)
if ( file_exists ) then
  open(unit=10, file=fin//fnx,status='old',access='direct',recl=(nkx)*totq*8)
  read(10,rec=1)a
  nx = reshape( int(a),(/ nkx,totq /) )
  close(10,status='delete')
  open(unit=10, file=fin//fny,status='old',access='direct',recl=(nkx)*totq*8)
  read(10,rec=1)a
  ny = reshape( int(a),(/ nkx,totq /) )
  close(10,status='delete')
else
  a=0
  nx = reshape( int(a),(/ nkx,totq /) )
  ny = reshape( int(a),(/ nkx,totq /) )
end if
deallocate(a)

allocate( a(nwan*nkx*nkx) )
open(unit=10, file=fin//fek,status='old',access='direct',recl=nwan*nkx*nkx*8)
read(10,rec=1)a
Ek = reshape( a,(/ nwan,nkx,nkx /) )
 close(10,status='delete')
open(unit=10, file=fin//ffe,status='old',access='direct',recl=nwan*nkx*nkx*8)
read(10,rec=1)a
fermi = reshape( a,(/ nwan,nkx,nkx /) )
 close(10,status='delete')
deallocate(a)

allocate( a(nwan*nkx*nkx*nwan))
open(unit=10, file=fin//fakre,status='old',access='direct',recl=nwan*nkx*nkx*nwan*8)
read(10,rec=1)a
akre = reshape( a,(/ nwan,nwan,nkx,nkx /) )
 close(10,status='delete')
open(unit=10, file=fin//fakim,status='old',access='direct',recl=nwan*nkx*nkx*nwan*8)
read(10,rec=1)a
akim = reshape( a,(/ nwan,nwan,nkx,nkx /) )
 close(10,status='delete')
ak = cmplx(akre,akim)
deallocate(a)

allocate( a(nwan))
open(unit=10, file=fin//fHi,status='old',access='direct',recl=nwan*8)
read(10,rec=1)a
Hamind = reshape( int(a),(/ orbitals,nsite /) )
  close(10)
deallocate(a)

allocate( a(orbitals**2 * nsite**2))
open(unit=10, file=fin//fIi,status='old',access='direct',recl=orbitals**2 * nsite**2 *8)
read(10,rec=1)a
chiind = reshape( int(a),(/ orbitals,nsite,orbitals,nsite /) )
  close(10)
deallocate(a)

allocate( a(nR*2) )
open(unit=10, file=fin//fRvec,status='old',access='direct',recl=nR*2*8)
read(10,rec=1)a
Rvec = reshape( a,(/ nR,2 /) )
 close(10)
deallocate(a)

allocate( a(2*nkx*nkx) )
open(unit=10, file=fin//fkv,status='old',access='direct',recl=2*nkx*nkx*8)
read(10,rec=1)a
kvec = reshape( a,(/ 2,nkx,nkx /) )
 close(10,status='delete')

 RETURN
END subroutine 
