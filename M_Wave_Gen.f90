 Module M_Wave_Gen

use M_General, only: nx,ny,nz,pi,x,y,z,u,v,w,phi,landa,gz,dt,tpstar,lx,ly,Lz,Froude


implicit none


real(8),dimension ( 0:1,1:ny-1,1:nz-1) :: uWa
real(8),dimension (-1:0,1:ny-1,0:nz+1) :: wWa
real(8),dimension(:),Allocatable       :: rndRealArr
real(8),dimension(:),Allocatable       :: xgage,ygage,zgage
real(8),dimension(1:ny-1)              :: ZGageGen
real(8)                                :: hait,wavegen,UCurrent,HH
real(8)                                :: H13,T1,OmegaMin,OmegaMax
integer                                :: Sample,GageN


contains


subroutine Wave_Gen_Constant_Ini()
implicit none 

include "Par_Constant_Wave.txt"

return

end subroutine 


subroutine Regularwave ()

implicit none 

real(8)  :: HH1,period,kwave,wwave,LRamp
integer  :: k,j





kwave=2*pi/(landa)
period=1.0/(  1/(2*pi)* sqrt(kwave*abs(gz)*tanh(kwave*hait)) )
wwave=2*pi/(period)


if ( ((tpstar-1)*dt).lt.(0.1*period) ) then 
HH1=((tpstar-1)*dt)/(0.1*period)*HH
else 
HH1=HH
end if

print*,"wave height=",HH1,"tpstar=",tpstar

!do k=1, nz-1 
!print*, z(k)
!end do 


     ZGageGen(:)=0 
    do j=1,ny-1 
      do k=1,nz-1
       if ((0.5*( phi(1,j,k+1)+phi(0,j,k+1) ).ge.0 )) then 
          !! may be more precsiness can be done above !!      
             Zgagegen(j)=z(k)   +  &
          &   ( z(k+1)-z(k) )/ ( 0.5*( phi(1,j,k+1)+phi(0,j,k+1) ) &
          &     -0.5*( phi(1,j,k)+phi(0,j,k) ) ) *( -0.5*( phi(1,j,k)+phi(0,j,k) ) ) 
        
         exit 
        end if 
       end do 
     end do 
    

UWa=0
do k=1,nz-1
do j=1,ny-1 
       

 uWa(1,j,k)=UCurrent+&
 & HH1/2* wwave *cos(-wwave*(tpstar*dt))*cosH(kwave* ( -abs(z(k)-zgagegen(j))+hait) )/sinH(kwave*hait)!+ &
    !& 3/4*((pi*HH1)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))* &
    !& cosH( 2*kwave* ((z(k)-zgage(1))+hait) )/( (sinH(kwave*hait))**4 )
 
 uWa(0,j,k)=UCurrent+&
 & HH1/2* wwave*cos(kwave*(0.5*(x(0)+x(-1)))-wwave*(tpstar*dt))*cosH(kwave* ( -abs(z(k)-zgagegen(j))+hait) )/sinH(kwave*hait) !+ &
    !& 3/4*((pi*HH1)**2)/(period*landa)*cos(2*kwave*(0.5*(x(0)+x(-1)))-2*wwave*(tpstar*dt))* &
    !& cosH( 2*kwave* ((z(k)-zgage(1))+hait) )/( (sinH(kwave*hait))**4 )
                        
      end do
end do 




WWa=0
 do k=0,nz+1
 do j=1,ny-1   

 wWa(0,j,k)=HH1/2*wwave*  &
& sin(kwave*x(0)-wwave*(dble(tpstar)*dt))*sinH(kwave* ( -abs(z(k)-zgagegen(j))+hait) )/sinH(kwave*hait) !+ &
     !& 3/4*((pi*HH1)**2)/(period*landa)*sin(2*kwave*x(0)-2*wwave*(tpstar*dt))* &
     !& sinH( 2*kwave* ((z(k)-zgage(1))+hait) )/( (sinH(kwave*hait))**4 )

 wWa(-1,j,k)=HH1/2*wwave* &
& sin(kwave*x(-1)-wwave*(dble(tpstar)*dt))*sinH(kwave* ( -abs(z(k)-zgagegen(j))+hait) )/sinH(kwave*hait) !+ &
     !& 3/4*((pi*HH1)**2)/(period*landa)*sin(2*kwave*x(-1)-2*wwave*(tpstar*dt))* &
     !& sinH( 2*kwave* ((z(k)-zgage(1))+hait) )/( (sinH(kwave*hait))**4 )
 
 
 end do   
 end do 



!! ramp for air velocity to reach zero within half of the distance from top
do j=1,ny-1
 
 LRamp=0.5*(Lz-ZGageGen(j))

  do k=1,nz-1
       
   if (z(k).gt.ZGageGen(j)) then 

    uWa(1,j,k)=uWa(1,j,k)* 1.0/LRamp*max( (ZGageGen(j)+LRamp-z(k)),0.0 )  
    uWa(0,j,k)=uWa(0,j,k)* 1.0/LRamp*max( (ZGageGen(j)+LRamp-z(k)),0.0 )  
 
   end if 
  end do 
end do 



 do j=1,ny-1 
  
  LRamp=0.5*(Lz-ZGageGen(j))
  do k=0,nz+1

   if (z(k).gt.ZGageGen(j)) then 

     wWa(0,j,k) =wWa(0,j,k) * 1.0/LRamp*max( (ZGageGen(j)+LRamp-z(k)),0.0 )
     wWa(-1,j,k)=wWa(-1,j,k)* 1.0/LRamp*max( (ZGageGen(j)+LRamp-z(k)),0.0 )
   end if 
  end do 
 end do 




return

end subroutine 
 








subroutine Wavegage()

 IMPLICIT NONE
 

   integer igage,jgage,i,j,k,kk



 

   do kk=1, size(xgage)    
   igage=1000

     do  i=0,nx-1
 
       if (x(i).le.xgage(kk).AND.x(i+1).ge.xgage(kk)) then 
         igage=i
         exit 
       end if 
     end do

     do  j=1,ny-1
 
       if (y(j).le.ygage(kk).AND.y(j+1).ge.ygage(kk)) then 
         jgage=j
         exit 
       end if 
     end do
 

     zgage(kk)=0
     do k=1,nz-1
       if (0.5*( phi(igage,jGage,k)+phi(igage-1,jGage,k) ).le.0.AND.0.5*( phi(igage,jGage,k+1)+phi(igage-1,jGage,k+1) ).ge.0 ) then 
         !! may be more precsiness can be done above !!      
         zgage(kk)=z(k)+ &
         &   ( z(k+1)-z(k) )/ ( 0.5*( phi(igage,jGage,k+1)+phi(igage-1,jGage,k+1) ) &
         &     -0.5*( phi(igage,jGage,k)+phi(igage-1,jGage,k) ) ) *( -0.5*( phi(igage,jGage,k)+phi(igage-1,jGage,k) ) ) 
         exit 
       end if 
     end do 

     if (igage.eq.1000) then 
       print *, "error in wave gage"
     end if 
    
   end do 




             
return        
  end subroutine  


subroutine RandomGenerator()



!!local !! 
integer                            :: seedSize
integer, dimension(:), allocatable :: seed
integer, dimension(8)              :: dtVals

 call DATE_AND_TIME(VALUES=dtVals)

write (*,*) 'Clock Values=', dtVals
call RANDOM_SEED(SIZE=seedSize)
if(seedSize .gt. 8) then
     write (*,*) 'ERROR: Seed size too large to init with DATE_AND_TIME return'
     stop
  end if

  ! We know the seed size, so allocate space for it and query to see what it is
  allocate(seed(seedSize)) 
  call RANDOM_SEED(GET=seed)
  
  
  write (*,*) 'Old Seed: ', seed

call RANDOM_SEED(PUT=dtVals((9-seedSize):8))

  ! Re-query the seed to make sure it worked.
  call RANDOM_SEED(GET=seed)
  write (*,*) 'New Seed: ', seed
  ! Get an array of random numbers
  rndRealArr = 0

Allocate ( rndRealArr(sample) )

  call RANDOM_NUMBER(rndRealArr)
  write (*,*) 'Random Array: ', rndRealArr

return 
end subroutine  




Subroutine Random()


implicit none 






!!! local

real(8),dimension(Sample) :: Omega,SOmega,Am,Landa2,period2,KWave2,first
integer i,j,k,size1(sample)

  
!OPEN(1150,file='PiersonSpectra.plt')
!write(1150,*) 'variables="Frequency(rad/Sec)","Spectrum"'



print*, "Random wave"
do i=1,Sample
Omega(i)=OmegaMin+real(i)/real(Sample)*(OmegaMax-OmegaMin)
SOmega(i)=(H13**2)*T1*0.11/(2*pi)*(Omega(i)*T1/(2*pi))**(-5)*exp(-0.44*(Omega(i)*T1/(2*pi))**-4)
Am(i)=sqrt(2*SOmega(i)*(OmegaMax-OmegaMin)/real(Sample))
KWave2(i)=(Omega(i)**2)/abs(gz) !! i depth assumption !!
period2(i)=2*pi/Omega(i)
landa2(i)=2*pi/KWave2(i)
!write(1150,*) Omega(i),SOmega(i)
end do 

!UWa=0.0 ; WWa=0.0

!do k=1,nz-1


! if (phi(0,ny/2,k).lt.0 ) then

!   do i=1,Sample
!   UWa(1,k)=UWa(1,k)+ &
!   & Am(i) * Omega(i) *cos(-Omega(i)*(tpstar*dt)+2*pi*rndRealArr(i))*cosH(KWave2(i)* &
!   & ( (z(k)-zgage(1))+hait) )/sinH(KWave2(i)*hait)  !+ &
!  ! & 3/4*((pi*Am(i))**2)/(period2(i)*landa2(i))*cos(-2*Omega(i)*(tpstar*dt))*cosH( 2*kwave2(i)* ((z(k)-zgage(1)) & 
!  ! & +hait) )/( (sinH(KWave2(i)*hait))**4 )

   
!  end do

! else
 
!  do i=1,Sample
!  UWa(1,k)=UWa(1,k)+ &
!  & Am(i)* Omega(i) *cos(-Omega(i)*(tpstar*dt)+2*pi*rndRealArr(i))*cosH(KWave2(i)* &
!  & (-(z(k)-zgage(1))+hait) )/sinH(KWave2(i)*hait) !+ &
! ! & 3/4*((pi*Am(i))**2)/(period2(i)*landa2(i))*cos(-2*Omega(i)*(tpstar*dt))*cosH( 2*kwave2(i)* (-(z(k)-zgage(1)) & 
! ! & +hait) )/( (sinH(KWave2(i)*hait))**4 )

  
!  end do

! end if 
! Uwa(0,k)=Uwa(1,k)    !! approximation

!end do 


return 
end subroutine 






end module 





