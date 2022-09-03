! Description: This code generates mean P_k. Parameters k and n can be varied.



  implicit none
  
  integer,parameter::n1=10000000,rea=10000,kup=10
  integer::i,j,irea,in,seed=-2853349,n,ai,k=1
  double precision::w(1:n1)=0.d0,z,ylist(1:rea,1:kup)=0.d0,y(1:kup)=0.d0,sm=0.d0,ytmp=0.d0
  double precision::y2(1:5)=0.d0,y2t(1:5)=0.d0,sav=0.d0,wmax=0.d0,wmaxav=0.d0
  double precision::ran2,kp=3.d0

  !!!!! parameters
  double precision::a=1.4d0

  
     n=4*10**4
  do ai=1,1
		!a=.1d0*ai

     ylist=0.d0
      w=0.d0
     y=0.d0

 do irea=1,rea

     
  z=ran2(seed)
  w(1)=(-log(z))**(-1.d0/a)
  do i=1,n-1
     z=ran2(seed)
     w(i+1)=w(i)*(1.d0-(w(i)**a)*log(z))**(-1.d0/a)
  enddo



     sm=0.d0
     do j=1,n
        sm=sm+w(j)
     enddo

     do i=1,6
        k=2.d0**(i*1.d0)
ytmp=0.d0
        do j=1,n
        ytmp=ytmp+(w(j)/sm)**(k)
        enddo
y(i)=y(i)+ytmp/(1.d0*rea)
     enddo


!print*,ytmp,wmax
  enddo

!print*,a,y(1),y(2)!,y(3)
enddo



!y2=y2*((1.d0*n)**(a-1.d0))



do i=1,6
print*,i,y(i)
enddo

!print*,rea
end program

FUNCTION ran2(idum)
!    gerador de números aleatórios baseado na combinacao
!    de dois geradores lineares congruenciais tendo um periodo
!    maior do que 2x10^18. A saida e' "baralhada" 
!     -> valor inicial de idum=-1  (Numerical recipes 2a. edicao)
IMPLICIT NONE
INTEGER, PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,&
     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,&
     ir2=3791,ntab=32,ndiv=1+imm1/ntab
DOUBLE PRECISION , PARAMETER ::   am=1.d0/im1,eps=1.d-14,rnmx=1.d0-eps
DOUBLE PRECISION :: ran2
INTEGER, DIMENSION(ntab) :: iv
INTEGER :: idum,idum2,j,k,iy

save iv,iy,idum2
data idum2/123456789/,iv /ntab*0/,iy /0/
      
if(idum.le.0) then
  idum=max(-idum,1)
  idum2=idum
  do j=ntab+8,1,-1
    k=idum/iq1
    idum=ia1*(idum-k*iq1)-ir1*k
    if(idum.lt.0) idum=idum+im1
    if(j.le.ntab) iv(j)=idum
   end do
   iy=iv(1)
endif

k=idum/iq1
idum=ia1*(idum-k*iq1)-ir1*k
if(idum.lt.0) idum=idum+im1
k=idum2/iq2
idum2=ia2*(idum2-k*iq2)-ir2*k

if(idum2.lt.0) idum2=idum2+im2

j=1+iy/ndiv
iy=iv(j)-idum2
iv(j)=idum
if (iy.lt.1)iy=iy+imm1
ran2=min(am*iy,rnmx)

END FUNCTION

  
     
  
