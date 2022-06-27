SUBROUTINE lindhardsus(orbitals,nsite,nwan,nobm,nkx,qq,pmfqbare, &
                       & nx,ny,ak,Ek,fermi,bosonfq,MatsuT,Hamind,chiind,nR,Rvec,kvec)

IMPLICIT NONE 

 INTEGER                 :: orbitals,nsite,nwan,nR,nkx,totk,k,l,bq,nobm,nxi,nyi
 INTEGER                 :: m1,m2,m3,m4,mu1,mu2,mu3,mu4,n1,n2,mn1,mn2,ir1,ir2
 integer                 :: row1,col1,row2,col2,qq,nx(:,:),ny(:,:),Hamind(:,:),chiind(:,:,:,:)
 real*8                  :: Ek(:,:,:),Rvec(:,:),kvec(:,:,:),fermi(:,:,:),bosonfq(:),beta,MatsuT
 real*8                  :: fermifac,den,den1,fermifac1,tol=10.d0**(-5),r1mr2(2),kxky(2)
 complex*16              :: psus,nsus,num,num1,jj=cmplx(0.d0,1.d0),ipbf,ak(:,:,:,:)
 Complex*16              :: pmfqbare(:,:,:,:,:),latfac

totk = (nkx)**2 !!Sum over all K-points in the grid
beta = 1./dble( 8.617 * 10.d0**(-5) * MatsuT )

do bq = 1,nobm
ipbf = 0.01*jj

do ir2 = 1,nR
do ir1 = 1,nR

r1mr2(1) = Rvec(ir1,1) - Rvec(ir2,1)
r1mr2(2) = Rvec(ir1,2) - Rvec(ir2,2)

do m4 = 1,nsite
 do mu4 = 1,orbitals
  do m3 = 1,nsite
   do mu3 = 1,orbitals
    do m2 = 1,nsite
     do mu2 = 1,orbitals
      do m1 = 1,nsite
       do mu1 = 1,orbitals

        psus = cmplx(0.d0,0.d0)   !! intialize chi_0 to zero

        do k = 1,nkx
          do l = 1,nkx  
           nyi = ny(l,qq)
           nxi = nx(k,qq)
  
           kxky(1:2) = kvec(1:2, l,k) 
           latfac = exp( jj * ( kxky(1) * r1mr2(1) + kxky(2) * r1mr2(2) ) )

              do mn1 = 1,nsite
              do n1 = 1,orbitals

              col1 = Hamind( n1,mn1 )
              row1 = Hamind( mu1,m1 )
              row2 = Hamind( mu4,m4 )
              num = conjg( ak(row1,col1,nyi,nxi) ) * ak(row2,col1,nyi,nxi)
              fermifac = fermi(col1,nyi,nxi)
              den = Ek(col1,nyi,nxi)

              do mn2 = 1,nsite
              do n2 = 1,orbitals

                col2 = Hamind( n2,mn2 )
                row1 = Hamind( mu2,m2 )
                row2 = Hamind( mu3,m3 )
                num1 = num * ak(row1,col2,l,k) * conjg( ak(row2,col2,l,k) )  
                den1 = den - Ek(col2,l,k) 

                if (abs(den1).lt.tol) then 
                  fermifac1 = - beta/4.0 * ( 1.0/cosh( den*beta/2.0 ) )**2
                  psus = psus + fermifac1 * latfac * num1
                else
                  fermifac1 = fermifac - fermi(col2,l,k)
                  psus = psus + fermifac1/(den1) * latfac * num1
                endif

              enddo  !! n2
              enddo  !! mn2
            enddo  !! n1
            enddo  !! mn1
            
          enddo  !! l       
        enddo  !! k

        row1 = chiind( mu1,m1,mu2,m2 )
        col1 = chiind( mu3,m3,mu4,m4 )
        pmfqbare(row1,col1,ir1,ir2,bq) = psus

      enddo  !! mu1
      enddo  !! m1
    enddo  !! mu2
    enddo  !! m2
  enddo  !! mu3
  enddo  !! m3
enddo  !! mu4
enddo  !! m4
enddo  !! ir1
enddo  !! ir2
enddo  !! bq

pmfqbare = - pmfqbare / dble(totk)

 RETURN
END subroutine 