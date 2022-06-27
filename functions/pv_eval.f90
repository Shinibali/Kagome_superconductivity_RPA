!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Shinibali, NOV 26,2016
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM pveval
use mymodule
 
IMPLICIT NONE 

 INTEGER                 :: INP,Uind,nwtktq=7,nsite,orbitals,nR,nwan,dimen,totk
 INTEGER                 :: mu1,mu2,mu3,mu4,n1,n2,k,kp,sp1,sp2,sp3,sp4,l,lp
 INTEGER                 :: Rblock_size
 integer,allocatable     :: bandind(:)
 real*8,allocatable      :: a(:),slp(:)
 Complex*16,allocatable  :: peigjk(:,:,:),meigjk(:,:,:),v_kkp(:,:,:,:),spin_vkkp(:,:,:,:,:,:)
 Complex*16,allocatable  :: gsij(:,:),sum0(:,:),gn1n2(:,:,:,:,:,:)
 complex*16              :: prod1,prod2
 Complex*16,allocatable  :: Gamma(:,:,:,:), Pauli_Gamma(:,:,:)
 character(:),allocatable:: fdid
 CHARACTER(LEN=10)       :: hhmmss,str
 CHARACTER(LEN=10000)    :: arg,outfile
 CHARACTER(*), PARAMETER :: fin='./input_files/',fout='./bareresults/'
 CHARACTER(*), PARAMETER :: fw='NsiteNorbNkxNqxNr.bin',ftotk='totk_onFS.bin',fics='interpchio_size.bin'

INP = IARGC()
call GETARG(1,arg)
fdid = trim(adjustl(arg))
call GETARG(2,arg)
read(arg,*) Uind

allocate(a(nwtktq) )
open(unit=10, file=fin//fw,status='old',access='direct',recl=nwtktq*8)
read(10,rec=1)a
nsite = int(a(1))
orbitals = int(a(2))
nR = int(a(7))
 close(10)
nwan = nsite * orbitals
dimen = nwan*nwan
Rblock_size = 2*2 * nsite**2 * orbitals**2
deallocate(a)

allocate(a(1))
open(unit=10, file=fin//ftotk,status='old',access='direct',recl=1*8)
read(10,rec=1)a
totk = int(a(1))
 close(10)
deallocate(a)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate( bandind(totk), peigjk(nwan,nwan,totk), meigjk(nwan,nwan,totk) )
allocate( v_kkp(Rblock_size,Rblock_size,totk,totk) )
CALL readinput_pairvertex(nR,nwan,totk,peigjk,meigjk,fdid,bandind,v_kkp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate( sum0(4,4),gn1n2(2,2,2,2,totk,totk) )

do kp = 1,totk
  n2 = bandind(kp)
  do k = 1,totk
     n1 = bandind(k)
     spin_vkkp = reshape( v_kkp(:,:,k,kp), (/ nwan,nwan,2*2,nwan,nwan,2*2 /) )

     sum0 = (0.d0,0.d0)
     do mu4 = 1,nwan
       do mu3 = 1,nwan
         do mu2 = 1,nwan
           do mu1 = 1,nwan
                                
              prod1 = conjg( peigjk( mu1,n1,k ) ) * conjg( meigjk( mu3,n1,k ) )
              prod2 = peigjk( mu4,n2,kp ) * meigjk( mu2,n2,kp )
              sum0 = sum0 + prod1 * spin_vkkp(mu1,mu2,:,mu3,mu4,:) * prod2
              
           enddo  !! mu1
         enddo  !! mu2
       enddo  !! mu3
     enddo  !! mu4
     gn1n2 (:,:,:,:,k,kp) = reshape( sum0 , (/ 2,2,2,2 /) )
  enddo  !! k
enddo  !! kp

allocate( Gamma(totk,totk,4,4), Pauli_Gamma(2,2,4) )
slp = (/ -1.0, 1.0, -1.0, 1.0 /)  !! x,z = +1; 0,y = -1
CALL pauli( Pauli_Gamma )
do lp = 1,4
  do l = 1,4
    do sp4 = 1,2
      do sp3 = 1,2
        do sp2 = 1,2
          do sp1 = 1,2
            Gamma(:,:,l,lp) = Gamma(:,:,l,lp) + Pauli_Gamma(sp3,sp1,l) * & 
                                              & gn1n2(sp1,sp2,sp3,sp4,:,:) * &
                                              & Pauli_Gamma(sp4,sp2,lp) * slp(lp)
          enddo !! sp1
        enddo !! sp2
      enddo !! sp3
    enddo !! sp4
  enddo !! l
enddo !! lp

write(str,"(I10.1)") Uind

outfile ='re_gamma_U'//trim(adjustl(str))//'.bin'
open(unit=8,file = fout//trim(outfile) ,status='replace',form='unformatted',access='stream')
write(8)real(Gamma)
 close(8)
outfile ='im_gamma_U'//trim(adjustl(str))//'.bin'
open(unit=8,file = fout//trim(outfile) ,status='replace',form='unformatted',access='stream')
write(8)aimag(Gamma)
 close(8)

!! CALL time()
!! write(*,*)'Finished evaluating Gamma(k,kp,l,lp) evaluation in pv_eval.f90'

END PROGRAM