SUBROUTINE readinput_pairvertex(nR,nwan,totk,peigjk,meigjk,fdid,bandind,v_kkp)

 IMPLICIT NONE

 INTEGER                 :: nR,nwan,dimen,totk,bandind(:),Rblock_size
 real*8,allocatable      :: a(:),eigre(:,:,:),eigim(:,:,:),re(:,:,:,:),im(:,:,:,:)
 Complex*16              :: peigjk(:,:,:),meigjk(:,:,:), v_kkp(:,:,:,:)
 character(:),allocatable:: fdid
 CHARACTER(*), PARAMETER :: fin='./input_files/',fout='./bareresults/'
 CHARACTER(*), PARAMETER :: fb='bandind.bin',fvre='re_vkkp.bin',fvim='im_vkkp.bin'
 CHARACTER(*), PARAMETER :: fpre='peig_re.bin',fpim='peig_im.bin',fmre='meig_re.bin',fmim='meig_im.bin'

fdid = trim(adjustl(fdid))
dimen = nwan*nwan

allocate( a(nwan*nwan*totk), eigre(nwan,nwan,totk), eigim(nwan,nwan,totk) )
open(unit=10, file=fin//fpre,status='old',access='direct',recl=nwan*nwan*totk*8)
read(10,rec=1)a
eigre = reshape( a,(/ nwan,nwan,totk /) )
 close(10)
open(unit=10, file=fin//fpim,status='old',access='direct',recl=nwan*nwan*totk*8)
read(10,rec=1)a
eigim = reshape( a,(/ nwan,nwan,totk /) )
 close(10)
peigjk = cmplx(eigre,eigim)

open(unit=10, file=fin//fmre,status='old',access='direct',recl=nwan*nwan*totk*8)
read(10,rec=1)a
eigre = reshape( a,(/ nwan,nwan,totk /) )
 close(10)
open(unit=10, file=fin//fmim,status='old',access='direct',recl=nwan*nwan*totk*8)
read(10,rec=1)a
eigim = reshape( a,(/ nwan,nwan,totk /) )
 close(10)
meigjk = cmplx(eigre,eigim)
deallocate(a)

allocate(a(totk))
open(unit=10, file=fin//fb,status='old',access='direct',recl=totk*8)
read(10,rec=1)a
bandind = int(a)
 close(10)
deallocate(a)

Rblock_size = 2*2 * nwan**2

allocate( a(Rblock_size**2 * totk**2 ))
allocate( re(Rblock_size,Rblock_size,totk,totk) )
allocate( im(Rblock_size,Rblock_size,totk,totk) )
open(unit=10, file=fdid//fvre,status='old',access='direct',recl=Rblock_size**2 * totk**2*8)
read(10,rec=1)a
re = reshape( a,(/ Rblock_size,Rblock_size,totk,totk /) )
 close(10)
open(unit=10, file=fdid//fvim,status='old',access='direct',recl=Rblock_size**2 * totk**2*8)
read(10,rec=1)a
im = reshape( a,(/ Rblock_size,Rblock_size,totk,totk /) )
 close(10)
v_kkp = cmplx(re,im)

 RETURN
END subroutine 