! Description: This code generates the distribtuion of P_2. 

  implicit none
  
  integer,parameter::n1=10000000,rea=2**15,kup=3,nu=200
  integer::i,j,irea,in,seed=-3254359,n,ai,bindex=1
  double precision::w(1:n1)=0.d0,z,ylist(1:rea,1:kup)=0.d0,y(1:kup)=0.d0,sm=0.d0,ytmp=0.d0,wmax=0.d0
  double precision::y2(1:5)=0.d0,y2t(1:5)=0.d0,sav=0.d0,bsize=0.d0,xy=0.d0,pdist(0:nu)=0.d0
  double precision::ran2,k=2.d0

  !!!!! parameters
  double precision::a=.5d0
   n=5000
  bsize=(1.d0/((1.d0*n)**(1.d0-1.d0)))/(1.d0*nu)
    

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

wmax=(w(1)/sm)**2

ytmp=0.d0
 do j=1,n
        ytmp=ytmp+(w(j)/sm)**(k)
     enddo
xy=wmax/bsize
bindex=int(xy)
!print*,bindex
pdist(bindex)=pdist(bindex)+1.d0
do i=1,5
  y2(i)=y2(i)+ytmp**i/(1.d0*rea)
enddo

  enddo

!print*,a,y(1),y(2)!,y(3)
enddo

sav=gamma(1.d0-1.d0/a)
y2t=0.d0
do i=1,5
y2t(i)=(a/sav**a)*gamma(i*k-a)*gamma(a)/gamma(k*i)
enddo


y2=y2*((1.d0*n)**(a-1.d0))

do i=1,5
!print*,i,y2(i),y2t(i)
enddo

!print*,'n',n,'rea',rea,'a',a,'k',k

pdist=pdist/(1.d0*rea*bsize)
do i=0,nu
print*,i*bsize,pdist(i)
enddo

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

  
     
  
