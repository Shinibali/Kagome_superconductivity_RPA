SUBROUTINE matinv( C, invC )

 IMPLICIT NONE 

 INTEGER                 :: m,info, error
 integer,allocatable     :: IPIV(:)
 Complex*16              :: C(:,:), invC(:,:)
 Complex*16,allocatable  :: WORK(:)

invC = C
m = size(C,1)

allocate(WORK(m),IPIV(m),stat=error)
if (error.ne.0)then
  print *,"error:not enough memory"
  stop
end if

call ZGETRF(m,m,invC,m,IPIV,info)
if(info .ne. 0) then
  write(*,*)"failed LU decomposition"
end if

call ZGETRI(m,invC,m,IPIV,WORK,m,info)
if(info .ne. 0) then
  write(*,*)"failed inverting matrix"
  stop
end if

RETURN
END subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE multmat( A,B,C,prod )
 IMPLICIT NONE

 Complex*16              :: A(:,:), B(:,:), C(:,:), prod(:,:)

prod = matmul(A,B)
prod = matmul(prod,C) !! This gives prod = A * B * C

RETURN
END subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE pauli( Pauli_Gamma )
 IMPLICIT NONE

 Complex*16              :: jj=cmplx(0.0,1.0),zeros=cmplx(0.0,0.0),Pauli_Gamma(:,:,:)
 Complex*16              :: pauli_0(2,2),pauli_x(2,2),pauli_y(2,2),pauli_z(2,2)

pauli_0 = reshape( (/ 1.0, 0.0, 0.0, 1.0 /),  (/ 2,2 /) )
pauli_x = reshape( (/ 0.0, 1.0, 1.0, 0.0 /),  (/ 2,2 /) )
pauli_y = reshape( (/ zeros, jj, -jj, zeros  /),  (/ 2,2 /) )
pauli_z = reshape( (/ 1.0, 0.0, 0.0, -1.0 /), (/ 2,2 /) )

Pauli_Gamma(:,:,1) = matmul(pauli_0,pauli_y) * jj/sqrt(2.0)
Pauli_Gamma(:,:,2) = matmul(pauli_x,pauli_y) * jj/sqrt(2.0)
Pauli_Gamma(:,:,3) = matmul(pauli_y,pauli_y) * jj/sqrt(2.0)
Pauli_Gamma(:,:,4) = matmul(pauli_z,pauli_y) * jj/sqrt(2.0)

RETURN
END subroutine
