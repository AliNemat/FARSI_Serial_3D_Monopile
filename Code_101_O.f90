
program Code_104
  use M_General
  use M_Wave_Gen
  use M_Solid
  use M_Mesh
  use M_Tether
  use M_Tower
  use M_Add_Force
   !use M_Graphic

  implicit none
  logical,parameter         :: START=.True.  

  !real(8),dimension (1:nx-1,1:ny-1,1:nz-1,8)        ::   Cof
  !real(8),dimension (1:nx-1,1:ny-1,1:nz-1,Iterate)  ::  presv

  real(8),dimension (1:nx,0:ny,0:nz)          ::   advectu
  real(8),dimension (0:nx,1:ny,0:nz)          ::   advectv
  real(8),dimension (0:nx,0:ny,1:nz)          ::   advectw

  real(8),dimension (2:nx-1,1:ny-1,1:nz-1)    ::   Tx   
  real(8),dimension (1:nx-1,2:ny-1,1:nz-1)    ::   Ty
  real(8),dimension (1:nx-1,1:ny-1,2:nz-1)    ::   Tz

  real(8),dimension (1:nx-1,1:ny-1,1:nz-1)    ::   divv
  real(8),dimension (0:nx,0:ny,0:nz)          ::   ro,p,miuv
  real(8),dimension (-1:nx+1,-1:ny+1,-1:nz+1) ::   phiold2,phiold

  real(8),dimension (0:nx+1,-1:ny+1,-1:nz+1)  ::   uold2,uold,uext
  real(8),dimension (-1:nx+1,0:ny+1,-1:nz+1)  ::   vold2,vold,vext
  real(8),dimension (-1:nx+1,-1:ny+1,0:nz+1)  ::   wold2,wold,wext
  
  real(8),dimension (3)                       ::   F_Wind,M_Wind

  real(8)                                     ::   div,divmax,maxp,eps,sumphi,gap
  real(8)                                     ::   EndTime,StartTime
  integer                                     ::   i,j,k,tp,numimp,torder 
  
  real(8)                                     ::   HV
  
  CHARACTER(len=5)                            :: NNUMBER 
  CHARACTER(len=15)                           :: FileName




350 format (3(1x,e15.7))
111 format (A110)
112 format (A200)
122 format (11(1x,e15.7)) 




   if (START) then

     call General_Constant_Ini()
     call General_Dynamic_Ini()
     Call meshgenerator()

     eps=2.5*( y(int(ny/2))-y(int(ny/2)-1) )    !  !2*hz(70)
     gap=0.3*( y(int(ny/2))-y(int(ny/2)-1) )   !1.2*hy(36)

     call Fluid_Dynamic_Ini(p)

     call Solid_Constant_Ini()
     call Solid_Finder_Dynamic_Ini(gap)
     call proprety(ro,rodrop,roair,rosolid,rocon,gap,eps)
     call Solid_Fluid_Dynamic_Ini(ro)
 
     call Tower_Constant_Ini()
     call Tower_Dynamic_Ini()

     call Tether_Constant_Ini()
     call Tether_Dynamic_Ini()
   
     call Wave_Gen_Constant_Ini()
 
     
   else !! Continue the Run 

     call General_Constant_Ini()
     call General_Dynamic_read()
     Call meshgenerator()

     eps=2.5*( y(int(ny/2))-y(int(ny/2)-1) )   
     gap=0.3*( y(int(ny/2))-y(int(ny/2)-1) )  
 
     call Fluid_Dynamic_Read(p)

     call Solid_Constant_Ini()
     call Solid_Dynamic_Read()
 
     call Tower_Constant_Ini()
     call Tower_Dynamic_Read()   !! mass of the tower is required for the pre-tension of tethers

     call Tether_Constant_Ini() 
     call Tether_Dynamic_read() 
   
     call Wave_Gen_Constant_Ini()


     
   end if 
   

   if (test) then 

     OPEN(195,file='grideval2DSolid.plt')
     write(195,*) 'zone i=',nx-1,' k=',nz-1
     do k=1,nz-1   
       do i=1,nx-1
         write(195,350) x(i),z(k),Ib_Solid(i,int(ny/2),k)
       end do 
     end do
     call flush (195)
    
   end if


        
   if (START) then
    
     OPEN(45,file='TLPlocation.plt')
     OPEN(55,file='waveanglemass.plt')
     OPEN(75,file='PlatformForce.plt')
     OPEN(1005,file='pressure2.plt')
     OPEN(1015,file='WaveGage.plt')
     OPEN(1025,file='TLPlocation2.plt')



     write(45,111) 'variables="t","xbar","ybar","zbar","ubar","vbar","wbar","omegax","omegay","omegaz"'  
     write(55,111) 'variables="Time(Sec)","TOTAL mass","Pitch angle"," Roll angle","Yaw angle","waveh","anacy","anacz"'
     write(75,112) 'variables="Time(Sec)","Drag_Int","Drag_Dif","Lift_Int","Lift_Dif","MDrag_Int(N.m)","MLift_Int(N.m)","M_Dif(N.m)","Drag_P(N.m)","MDrag_P(N.m)"'     
     write(1005,*) 'variables="x","y","z","pressure"'
     write(1015,112) 'variables="t","Gage1","Gage2","Gage3","Gage4","Gage5","Gage6","Gage7",&
                &        "Gage1b","Gage2b","Gage3b","Gage4b","Gage5b","Gage6b","Gage7b"'
     write(1025,111) 'variables="t","xbar","ybar","zbar","ubar","vbar","wbar","omegax","omegay","omegaz"'  

     write(unit=NNUMBER, fmt='(I5.5)') 0+10000 
     FileName='Result'//NNUMBER//'.dat'

     OPEN(unit=0+10000,file=FileName)   
     write(0+10000,111) 'variables="x","y","z","u","v","w","phi","ro","p","ib","div"' 
     write(0+10000,*) 'ZONE T="',0,'" STRANDID=',500000, 'SOLUTIONTIME=',0,'i=',nx-1,' j=',ny-1,' k=',nz-1 ,'DATAPACKING=POINT'
 
     do k=1,nz-1 ; 
       do j=1,ny-1 
         do i=1,nx-1
           write(0+10000,122) x(i),y(j),z(k),0.5*(u(i,j,k)+u(i+1,j,k)),0.5*(v(i,j,k)+v(i,j+1,k)), &
           &               0.5*(w(i,j,k)+w(i,j,k+1)),phi(i,j,k),ro(i,j,k),p(i,j,k),Ib_Solid(i,j,k),0.0
         end do 
       end do 
     end do
  
     !call solidmeshgen ( 0,dt,z(1))
                   
     call flush (0+10000)
     close(0+10000)
   else  !! continue
         
     OPEN(45,file='TLPlocation.plt',POSITION='APPEND',status='old')
     OPEN(55,file='waveanglemass.plt',POSITION='APPEND',status='old')
     OPEN(75,file='PlatformForce.plt',POSITION='APPEND',status='old')
     OPEN(1005,file='pressure2.plt',POSITION='APPEND',status='old')
     OPEN(1015,file='WaveGage.plt',POSITION='APPEND',status='old')
     OPEN(1025,file='TLPlocation2.plt',POSITION='APPEND',status='old')

   end if 






   call CPU_TIME(StartTime)
   print*,"Start time is", StartTime

   !do tp=1,tstep !time step
   do tp=tstepS,tstepE !time step      
      
     print*,"time step=",tp
     call flush (6)

     uold(:,:,:)=u(:,:,:)             
     vold(:,:,:)=v(:,:,:)             
     wold(:,:,:)=w(:,:,:)             
     phiold(:,:,:)=phi(:,:,:)             
     
     call Solid_Sec_Sav
     call Solid_Fluid_Sec_Sav
     call Tether_Sec_Sav
     Call Tower_Sec_Sav
       
     Force_Sav(:)=Force(:)
     Mom_Sav(:)=Mom(:)
      do torder=1,2
     
         !! Three Equations in Time !! 1- Solid motion, 2- Level set equation and 3- Navier Stokes equations
         uold2(:,:,:)=u(:,:,:)             
         vold2(:,:,:)=v(:,:,:)             
         wold2(:,:,:)=w(:,:,:)             
         phiold2(:,:,:)=phi(:,:,:)   
          
         call Solid_It_Sav
         call Solid_Fluid_It_Sav
         call Tether_It_Sav
         call Tower_It_Sav
       
         tpstar=tp-1+torder
 
         numimp=1 ; alphaconv=100 
 
         do while (numimp.lt.Iterate) 

           call proprety(ro,rodrop,roair,rosolid,rocon,gap,eps)
           call proprety(miuv,miudrop,miuair,miusolid,miucon,gap,eps)
           call advection(advectu,advectv,advectw)
           call viscoset(miuv,ro,Tx,Ty,Tz)
           !call AddForce(ro)   
           fbox(:,:,:)=0
           fboy(:,:,:)=0
           fboz(:,:,:)=0
   
           Do i=2,nx-1 
             Do j=1, ny-1 
               do k=1, nz-1
                 U(i,j,k)=Uold2(i,j,k)+(advectu(i,j,k)+Tx(i,j,k)+gx +fbox(i,j,k))*dt
               end do 
             end do 
           end do  

           do j=2,ny-1 
             do i=1,nx-1 
               do k=1, nz-1
                 V(i,j,k)=Vold2(i,j,k)+( advectv(i,j,k)+Ty(i,j,k)+gy  +fboy(i,j,k))*dt
               end do 
             end do 
           end do  

           do j=1,ny-1 
             do i=1,nx-1 
               do k=2, nz-1
                 W(i,j,k)=Wold2(i,j,k)+( advectW(i,j,k)+Tz(i,j,k)+gz  +fboz(i,j,k))*dt  
               end do 
             end do 
           end do  


           call BoundCond_UVW_Poisson() !! The velocities are defined and are not based on the 

           !  call pCofent(ro,dt,Cof)
           !p(:,:,:)=presv(:,:,:,numimp)
           if (isolver.eq.1) then
 
             ! call SORpoison (beta,pdif,p,maxSOR,COF,test)
             call poisson (ro,beta,pdif,p,dt)

           elseif(isolver.eq.2) then

             !  call BICGSTABpoison(toler,print_resid,nonzero_x,mxmatvec,p,Cof)

           else 

             !   call hyprepoison(maxIteration,maxError,p,cof)

           end if 
           !presv(:,:,:,numimp)=p(:,:,:)
     
                    

           do i=2,nx-1
             do j=1,ny-1 
               do k=1, nz-1
                 U(i,j,k)=U(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )   )*( p(i,j,k)-p(i-1,j,k) )*dt/( x(i)-x(i-1)  )
               end do 
             end do 
           end do 

           do i=1,nx-1 
             do j=2,ny-1 
               do k=1, nz-1
                 V(i,j,k)=V(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )   )*( p(i,j,k)-p(i,j-1,k) )*dt/( y(j)-y(j-1)  )
               end do 
             end do 
           end do 

           do i=1,nx-1 
             do j=1,ny-1 
               do k=2, nz-1
                 W(i,j,k)=W(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  )   )*( p(i,j,k)-p(i,j,k-1) )*dt/( z(k)-z(k-1)  )
               end do 
             end do 
           end do
 
           call boundarycond()    !! it is only for divergence in the next line other wise it could be merged with after solid if the solid is not on the boundaries

           divmax=0
           do i=1,nx-1 
             do j=1,ny-1 
               do k=1,nz-1

                 div= (U(i+1,j,k)-U(i,j,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
                 &(v(i,j+1,k)-v(i,j,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+&
                 &(W(i,j,k+1)-W(i,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) 
    
                divv(i,j,k)= (U(i+1,j,k)-U(i,j,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
                &(v(i,j+1,k)-v(i,j,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+&
                &(W(i,j,k+1)-W(i,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )  
 
                 if(abs(div).gt.divmax)then
                   divmax=abs(div)
                 end if 

               end do 
             end do 
           end do 
           print*,'DIVERGENCE before correction=',DIVMAX
           if (torder==2) then 
             call Solid_load_Ini(ro) 
           end if 

           call Solid_Fluid(ro)

           if (torder==2) then 
             call Solid_load_Mod(ro)
             call Solid_load_Final(ro,p,Advectu,Advectv,AdvectW,Tx,Ty,Tz,tpstar) 
           end if 

           call Solid_Dynamic_Update(Gap)
           call Tether_Dynamic_Update
           call Tether_FM           

           call Tower_Dynamic_Update
           call Tower_FM

           call boundarycond()  !! since the  solid has intersection with boundary  //for monopile activated

           F_Wind(:)=0 ;  M_Wind(:)=0
           F_Tow(1)=0 

           Force(:)=F_Tow(:)   +F_Wind(:)  +FTether(:) 
           Mom  (:)=M_Tow_Cg(:)+M_Wind(:)  +MomTether(:) 
          

           print*, "External force in x,y and z",Force (1),Force(2),Force(3)
           print*, "External moment in x,y and z",Mom(1),Mom(2),Mom(3)          
          divmax=0
          do i=1,nx-1 ;
            do j=1,ny-1 
              do k=1,nz-1
                div= (U(i+1,j,k)-U(i,j,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
                &(v(i,j+1,k)-v(i,j,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+&
                &(W(i,j,k+1)-W(i,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) 
 

 
                if(abs(div).gt.divmax)then
                  divmax=abs(div)
                end if 

              end do 
            end do 
          end do 


          call levelset(phiold2,dt) 
          !call contactangle(uext,vext,wext,gap,test) 
          !call levelset2(uext,vext,wext,0.25*dt,test)
          call reinitialize(eps,dt)



          print*,'DIVERGENCE Iterate=',DIVMAX
          print *,"numimp=",numimp,"alphaconverge=",alphaconv

          numimp=numimp+1



        end do !! implicite!!


      end do   !! torder !!

      u=0.5d0*(u+uold)             
      v=0.5d0*(v+vold)    
      w=0.5d0*(w+wold)          
      phi=0.5d0*(phi+phiold) 
      Force(:)=0.5*(Force(:)+Force_Sav(:))
      Mom(:)  =0.5*(Mom(:)  +Mom_Sav(:)  )

      call Solid_Sec_Cor            
      call Solid_Fluid_Sec_Cor
      call Tether_Sec_Cor  
      call Tower_Sec_Cor      
 

!! Post Processing !




!!pp1
write(45,2010) tp*dt,xbar,ybar,zbar,ubar,vbar,wbar,omegax,omegay,omegaz
call flush (45)

write(1025,2010) tp*dt,point_Cyl_1(1,1),point_Cyl_1(1,2),point_Cyl_1(1,3),ubar,vbar,wbar,omegax,omegay,omegaz
call flush (1025)


2010  format (10(1x,e15.7))

!!pp2
divmax=0
do i=1,nx-1 ;do j=1,ny-1 ;do k=1,nz-1
div= (U(i+1,j,k)-U(i,j,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
    &(v(i,j+1,k)-v(i,j,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+&
    &(W(i,j,k+1)-W(i,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) 

  !divv(i,j,k)= (U(i+1,j,k)-U(i,j,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
  !  &(v(i,j+1,k)-v(i,j,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+&
  !  &(W(i,j,k+1)-W(i,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )   

    
if(abs(div).gt.divmax)then
divmax=abs(div)
else
end if 
end do ; end do ;end do 

print*,'DIVERGENCE  solid =',DIVMAX



!!pp3
sumphi=0 
do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 
if (phi(i,j,k).lt.0) then 
sumphi=sumphi+1
end if 
end do ; end do  ; end do 

!write(55,2011) tp*dt,sumphi,ATAN( ( px(1,1)-px(2,2) )/sqrt(len**2-( px(1,1)-px(2,2) )**2))*180/pi, &             !! it became corrected !!
!                         &  ATAN( ( py(1,1)-py(2,2) )/sqrt(len**2-( py(1,1)-py(2,2) )**2))*180/pi, &
!                         &  ASIN( ( py(6,4)-py(5,4) )/(2*r2) )*180/pi,zfree,zfreeg


print*,"total mass=", sumphi
write(55,2011) tp*dt,sumphi,ATAN( (point_Cyl_1(1,1)-point_Cyl_1(2,1))/(point_Cyl_1(1,3)-point_Cyl_1(2,3)) )*180/pi, &             !! it became corrected (2) !!
                         &  ATAN( (point_Cyl_1(1,2)-point_Cyl_1(2,2))/(point_Cyl_1(1,3)-point_Cyl_1(2,3)) )*180/pi, & 
                         &  0.0,zgage(1),anacy,anacz
                         
2011  format (8(1x,e15.7))
call flush(55)



call Tether_Plot(tp)
call Tower_Plot(tp)



write(1015,2020) tp*dt,zgage(1),zgage(2),zgage(3),zgage(4),zgage(5),zgage(6),zgage(7), &
&                      zgage(8),zgage(9),zgage(10),zgage(11),zgage(12),zgage(13),zgage(14)
2020  format (15(1x,e15.7))
call flush (1015)


     !!pp5   
     countPlot=countPlot+1
     if (countPlot.eq.plot)then 

       countPlot=0
       print*, "Write data"
       write(unit=NNUMBER, fmt='(I5.5)') tp+10000 
       FileName='Result'//NNUMBER//'.dat'

       OPEN(unit=tp+10000,file=FileName)
       write(tp+10000,111) 'variables="x","y","z","u","v","w","phi","ro","p","ib","Div"'  
       write(tp+10000,*) 'ZONE T="',tp,'" STRANDID=',tp, 'SOLUTIONTIME=',tp*dt,'i=',nx-1,' j=',ny-1,' k=',nz-1 ,'DATAPACKING=POINT'
       !!write(35,*) 'zone i=',nx-1,'j=',ny-1 ,' k=',nz-1
       do k=1,nz-1 
         do j=1,ny-1 
           do i=1,nx-1
             write(tp+10000,122) x(i),y(j),z(k),0.5*(u(i,j,k)+u(i+1,j,k)),0.5*(v(i,j,k)+v(i,j+1,k)), &
             &            0.5*(w(i,j,k)+w(i,j,k+1)),phi(i,j,k),ro(i,j,k),p(i,j,k),Ib_Solid(i,j,k),abs(divv(i,j,k))

           end do 
         end do
       end do

       !call solidmeshgen ( tp,dt,z(1))
       call flush(tp+10000)
       close(tp+10000)

     end if



     call CPU_TIME(EndTime)
     print*,"Time elapsed after last backup",EndTime-StartTime,"seconds"
     if ((EndTime-StartTime).gt. (SavMin*60.0) ) then

       call CPU_TIME(StartTime)
       call General_Dynamic_Write(tp)
       call Fluid_Dynamic_Write(p)
       call Solid_Dynamic_Write()             
       call Tower_Dynamic_Write()
       call Tether_Dynamic_Write() 
      
     end if 

 


   end do !time step

write(*,*)'end'
read (*,*)

    end program 




!!!!!!!!!!!!!!!!!!!!!!! subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine boundarycond()

use M_Wave_Gen

implicit none 

!! Local !!
real(8) sumww,SumT
real(8) AAA
integer i,j,k



SumWW=0 
SumT=0




!!Wavegen 1=no wave, 2=periodic wave, 3=periodic wave+current, 4=random wave


 call boundarycondphi()

 call Wavegage()






 if (wavegen.eq.4) then
   call RandomGenerator() 
   call Random()

 else if (wavegen.eq.2 ) then 
   call regularWave()
 end if 

   
if (wavegen.eq.1) then 

   do j=1, ny-1 ; do k=1, nz-1  
     U(1,j,k)=0       !U(2,j,k)      
     U(0,j,k)=-U(2,j,k)
   end do ; end do
else

 

   do j=1 ,ny-1 ; do k=1,nz-1
     u(1,j,k)=Uwa(1,j,k)
     u(0,j,k)=Uwa(0,j,k)
      sumww=sumww+u(1,j,k)*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )
  
     ! sumww=sumww+u(1,j,k)*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )
 end do ; end do 
end if                         
 

 
    
 if (wavegen.eq.3) then !! wave + current

  do j=1 ,ny-1 ; do k=1,nz-1 
    u(nx,j,k)=u(nx-1,j,k)
    u(nx+1,j,k)=u(nx-2,j,k)
  end do ; end do 

else 

  do j=1 ,ny-1 ; do k=1,nz-1 
   u(nx,j,k)=0 
   u(nx+1,j,k)=-u(nx-1,j,k)
  end do ; end do 


end if 
 
  
print*," boundary summation=",sumww


Do i=0, nx+1 ;do k=1 , nz-1
U(i,0,k)=U(i,1,k)          
U(i,-1,k)=U(i,2,k)        
U(i,ny,k)=U(i,ny-1,k)
U(i,ny+1,k)=U(i,ny-2,k)
end do ; end do
do i=0,nx+1 ;do j=0,ny+1
U(i,j,-1)=-U(i,j,2)
U(i,j,0)=-U(i,j,1)
U(i,j,nz)=-U(i,j,nz-1)
U(i,j,nz+1)=-U(i,j,nz-2)
end do ; end do 






Do i=1, nx-1 ; do k=1, nz-1
V(i,1,k)=0
V(i,0,k)=-V(i,2,k)
V(i,ny,k)=0
V(i,ny+1,k)=-V(i,ny-1,k)
end do ; end do 



do j=0,ny+1 ;  do k=1 , nz-1
   
   V(0,j,k)=V(1,j,k)
   V(-1,j,k)=V(2,j,k)
   V(nx,j,k)=V(nx-1,j,k)
   V(nx+1,j,k)=V(nx-2,j,k)
    
end do ; end do 

    
do i=0,nx+1 ;do j=0,ny+1       
V(i,j,-1)=V(i,j,2)
V(i,j,0) =V(i,j,1)

V(i,j,nz)=V(i,j,nz-1)
V(i,j,nz+1)=V(i,j,nz-2)
end do ; end do





if (wavegen.eq.1.OR.wavegen.eq.3) then
 
   do j=1 ,ny-1 ; do i=1,nx-1
   W(i,j,1)=0
   w(i,j,0)=-w(i,j,2)
   W(i,j,nz)=0
   W(i,j,nz+1)=-W(i,j,nz-1)
   end do ; end do
   
else if (wavegen.eq.2.OR.wavegen.eq.4) then

   
!   do j=1 ,ny-1 ; do i=1,nx-1
!   SumT=SumT+( 2*W(i,j,nz-1)-W(i,j,nz-2) )*( 0.5*(x(i+1)+x(i))-0.5*(x(i)+x(i-1)) )*( 0.5*(y(j+1)+y(j))-0.5*(y(j)+y(j-1)) )
!   end do ; end do


!   print*,"SumT=",sumT,"SumWW=",sumww 

   !if (sumT.eq.0) then 

   do j=1 ,ny-1 ; do i=1,nx-1
    
    W(i,j,1)=0
    w(i,j,0)=-w(i,j,2)
    AAA=( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )*( 0.5*(y(j+1)+y(j))-0.5*(y(j)+y(j-1)) ) 
    W(i,j,nz)=sumww/( real( nx-1 )*real( ny-1)*AAA )
    W(i,j,nz+1)=W(i,j,nz)
   end do ; end do 
  ! else 
    ! do j=1 ,ny-1 ; do i=1,nx-1
    !  W(i,j,1)=0
    !  w(i,j,0)=-w(i,j,2)
    !  W(i,j,nz)=sumww/SumT*( 2*W(i,j,nz-1)-W(i,j,nz-2) ) 
    !  W(i,j,nz+1)=W(i,j,nz)
   
      ! W(i,j,nz)=W(i,j,nz-1)
      ! W(i,j,nz+1)=W(i,j,nz-2)
     !end do ; end do 
   ! end if 

end if 
  

if (wavegen.eq.1) then 


   do k=0,nz+1 ; do j=1,ny-1 
     w(0,j,k)=w(1,j,k)
     w(-1,j,k)=w(2,j,k)
     w(nx,j,k)=w(nx-1,j,k)
     w(nx+1,j,k)=w(nx-2,j,k)
   end do ; end do
print*, "Velocity is not defined in gravity direction"
   
else 

   do j=1,ny-1 ; do k=0,nz+1
   
     w(0,j,k)=WWa(0,j,k)
     w(-1,j,k)=WWa(-1,j,k)
     w(nx,j,k)=w(nx-1,j,k)
     w(nx+1,j,k)=w(nx-2,j,k)
   end do ; end do 
end if 
  

do i=0,nx+1 ; do k=0,nz+1
w(i,0,k)=w(i,1,k)
w(i,-1,k)=w(i,2,k)
w(i,ny,k)=w(i,ny-1,k)
w(i,ny+1,k)=w(i,ny-2,k)
end do ; end do


return 
end 

Subroutine BoundCond_UVW_Poisson()

 use M_Wave_Gen

  implicit none 

  !! Local !!
  real(8) sumww,SumT
  real(8) AAA
  integer i,j,k



   SumWW=0 
   SumT=0

   !!Wavegen 1=no wave, 2=periodic wave, 3=periodic wave+current, 4=random wave


   if (wavegen.eq.4) then
     call RandomGenerator() 
     call Random()
   else if (wavegen.eq.2 ) then 
     call regularWave()
   else
     print*, "Error Boundary Velocity Poisson"
   end if 


   do j=1 ,ny-1 
     do k=1,nz-1
       u(1,j,k)=Uwa(1,j,k)
       sumww=sumww+u(1,j,k)*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )
       ! sumww=sumww+u(1,j,k)*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )
     end do ;
   end do 
 
   do j=1 ,ny-1 
     do k=1,nz-1 
       u(nx,j,k)=0 
     end do 
    end do 


   do i=1, nx-1 
     do k=1, nz-1
       V(i,1,k)=0
       V(i,ny,k)=0
     end do 
   end do 

   
   do j=1 ,ny-1 
     do i=1,nx-1    
       W(i,j,1)=0
       AAA=( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )*( 0.5*(y(j+1)+y(j))-0.5*(y(j)+y(j-1)) ) 
       W(i,j,nz)=sumww/( real( nx-1 )*real( ny-1)*AAA )
     end do 
   end do 
    



 return 
end 



    subroutine Fluid_Dynamic_Ini(p)


      use M_General,              only: nx,ny,nz,x,y,z,u,v,w,phi,zfree
      implicit none

      real(8),dimension (0:nx,0:ny,0:nz)          ::   p
      integer i,j,k
 


      
      do i=1,nx-1 
        do j=1,ny-1 
          do k=1,nz-1
            phi(i,j,k)=(z(k)-zfree) 
          end do 
        end do 
      end do 

  
      U=0.0 
      V=0.0 
      w=0.0 
      p=0.0

      !!!!!!!!!!!!!!!!!! for free surface testing with Wu et al. (2001)

      !phi=0
      !do k=1,nz-1 ; do i=1,nx-1
      !
      !phi(i,3,k)=10000
      !
     
     !do ipr=1,nx-1
     !phiges=sqrt(   ( x(i)-x(ipr) )**2  +  (  z(k)-( 1.0+0.05*cos(pi*(x(ipr)) ) )  )**2   )
     !if (phiges.lt.phi(i,3,k) ) then 
     !phi(i,3,k)=phiges
     !end if 
     !end do 
     !
     !if ( z(k).lt.( 1.0+0.05*cos(pi*(x(i)) ) ) ) then 
     !phi(i,3,k)=-phi(i,3,k)
     !end if 
     !
     !end do ; end do 
     !
     !do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 
     !phi(i,j,k)=phi(i,3,k)
     !end do ; end do ; end do 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end testing !!!!!!!!!!!

     print*, "Fluid Dynamic Initial  condition finished"


     return 
   end subroutine 

 subroutine advection(advectu,advectv,advectw)
      
 use M_General, only: nx,ny,nz,x,y,z,u,v,w
 implicit none 
 integer i,j,k
 
 real(8) ux,uy,uz,vx,vy,vz,wx,wy,wz
 real(8),dimension (1:nx,0:ny,0:nz)       ::   dpux,dmux,dpuy,dmuy,dpuz,dmuz,advectu
 real(8),dimension (0:nx,1:ny,0:nz)       ::   dpvx,dmvx,dpvy,dmvy,dpvz,dmvz,advectv
 real(8),dimension (0:nx,0:ny,1:nz)       ::   dpwx,dmwx,dpwy,dmwy,dpwz,dmwz,advectw

 do i=1,nx ; do j=0,ny ; do k=0,nz
 Dpux(i,j,k)=( u(i+1,j,k)-u(i,j,k)   )/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )
 Dmux(i,j,k)=( u(i,j,k)  -u(i-1,j,k) )/( 0.5*(x(i)  +x(i-1))-0.5*(x(i-1)+x(i-2)) )
 Dpuy(i,j,k)=( u(i,j+1,k)-u(i,j,k)   )/( y(j+1)-y(j) )
 Dmuy(i,j,k)=( u(i,j,k)  -u(i,j-1,k) )/( y(j)-y(j-1) )
 Dpuz(i,j,k)=( u(i,j,k+1)-u(i,j,k)   )/( z(k+1)-z(k) )
 Dmuz(i,j,k)=( u(i,j,k)  -u(i,j,k-1) )/( z(k)-z(k-1) )
 end do ;end do ; end do 

 do i=2,nx-1 ; do j=1,ny-1 ; do k=1,nz-1     !!A!!
 
 
 !!U1    
 if (u(i,j,k).gt.0.0) then

    if (  abs(  Dmux(i,j,k)-Dmux(i-1,j,k) ).lt.abs(  Dpux(i,j,k)-Dpux(i-1,j,k) )   ) then
 ux=Dmux(i,j,k)+  0.5*(  Dmux(i,j,k)-Dmux(i-1,j,k)  ) 
    else
 ux=Dmux(i,j,k)+  0.5*(  Dpux(i,j,k)-Dpux(i-1,j,k)  )
    end if 
    
 else

    if (  abs(  Dmux(i+1,j,k)-Dmux(i,j,k) ).lt.abs(  Dpux(i+1,j,k)-Dpux(i,j,k) )   ) then
 ux=Dpux(i,j,k)-  0.5*(  Dmux(i+1,j,k)-Dmux(i,j,k)  ) 
    else
 ux=Dpux(i,j,k)-  0.5*(  Dpux(i+1,j,k)-Dpux(i,j,k)  )
    end if 

 end if 


!!U2
 if (  (0.25*( V(i,j,k)+V(i,j+1,k)+v(i-1,j,k)+v(i-1,j+1,k) )).gt.0.0) then

    if (  abs(  DmuY(i,j,k)-DmuY(i,j-1,k) ).lt.abs(  DpuY(i,j,k)-DpuY(i,j-1,k) )   ) then
 uY=DmuY(i,j,k)+  0.5*(  DmuY(i,j,k)-DmuY(i,j-1,k)  ) 
    else
 uY=DmuY(i,j,k)+  0.5*(  DpuY(i,j,k)-DpuY(i,j-1,k)  )
    end if 
    
 else

    if (  abs(  DmuY(i,j+1,k)-DmuY(i,j,k) ).lt.abs(  DpuY(i,j+1,k)-DpuY(i,j,k) )   ) then
 uY=DpuY(i,j,k)-  0.5*(  DmuY(i,j+1,k)-DmuY(i,j,k)  ) 
    else
 uY=DpuY(i,j,k)-  0.5*(  DpuY(i,j+1,k)-DpuY(i,j,k)  )
    end if 

 end if 
 
 
 !!!U3
  if (  (0.25*( W(i,j,k)+W(i-1,j,k)+W(i,j,k+1)+W(i-1,j,k+1) )).gt.0.0) then
  
    if (  abs(  DmuZ(i,j,k)-DmuZ(i,j,k-1) ).lt.abs(  DpuZ(i,j,k)-DpuZ(i,j,k-1) )   ) then
 uZ=DmuZ(i,j,k)+  0.5*(  DmuZ(i,j,k)-DmuZ(i,j,k-1)  ) 
    else
 uZ=DmuZ(i,j,k)+  0.5*(  DpuZ(i,j,k)-DpuZ(i,j,k-1)  )
    end if
    
 else 
    
    if (  abs(  Dmuz(i,j,k+1)-Dmuz(i,j,k) ).lt.abs(  DpuZ(i,j,k+1)-DpuZ(i,j,k) )   ) then
 uZ=DpuZ(i,j,k)-  0.5*(  DmuZ(i,j,k+1)-DmuZ(i,j,k)  ) 
    else
 uZ=DpuZ(i,j,k)-  0.5*(  DpuZ(i,j,k+1)-DpuZ(i,j,k)  )
    end if 
    
 end if 
    
 
 advectu(i,j,k)=-(u(i,j,k)*uX+0.25*( V(i,j,k)+V(i,j+1,k)+v(i-1,j,k)+v(i-1,j+1,k) )*uY+ &
 &                0.25*( W(i,j,k)+W(i-1,j,k)+W(i,j,k+1)+W(i-1,j,k+1) )*uZ)

 end do ; end do; end do 
  
 do i=0,nx ; do j=1,ny ; do k=0,nz
 Dpvx(i,j,k)=( v(i+1,j,k)-v(i,j,k)   )/( x(i+1)-x(i) )
 Dmvx(i,j,k)=( v(i,j,k)  -v(i-1,j,k) )/( x(i)-x(i-1) )
 Dpvy(i,j,k)=( v(i,j+1,k)-v(i,j,k)   )/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )
 Dmvy(i,j,k)=( v(i,j,k)  -v(i,j-1,k) )/( 0.5*(y(j)  +y(j-1))-0.5*(y(j-1)+y(j-2)) )
 Dpvz(i,j,k)=( v(i,j,k+1)-v(i,j,k)   )/( z(k+1)-z(k) )
 Dmvz(i,j,k)=( v(i,j,k)  -v(i,j,k-1) )/( z(k)-z(k-1) )
 end do ;end do ; end do 

 do i=1,nx-1 ; do j=2,ny-1  ; do k=1,nz-1 !!A!!
 
 
!!V1
 if (  (0.25*( u(i,j,k)+u(i+1,j,k)+u(i,j-1,k)+u(i+1,j-1,k) )).gt.0.0) then

    if (  abs(  Dmvx(i,j,k)-Dmvx(i-1,j,k) ).lt.abs(  Dpvx(i,j,k)-Dpvx(i-1,j,k) )   ) then
 vx=Dmvx(i,j,k)+  0.5*(  Dmvx(i,j,k)-Dmvx(i-1,j,k)  ) 
    else
 vx=Dmvx(i,j,k)+  0.5*(  Dpvx(i,j,k)-Dpvx(i-1,j,k)  )
    end if 
    
 else

    if (  abs(  Dmvx(i+1,j,k)-Dmvx(i,j,k) ).lt.abs(  Dpvx(i+1,j,k)-Dpvx(i,j,k) )   ) then
 vx=Dpvx(i,j,k)-  0.5*(  Dmvx(i+1,j,k)-Dmvx(i,j,k)  ) 
    else
 vx=Dpvx(i,j,k)-  0.5*(  Dpvx(i+1,j,k)-Dpvx(i,j,k)  )
    end if 

 end if 



!!V2
 if ( V(i,j,k).gt.0.0) then

    if (  abs(  DmvY(i,j,k)-DmvY(i,j-1,k) ).lt.abs(  DpvY(i,j,k)-DpvY(i,j-1,k) )   ) then
 vY=DmvY(i,j,k)+  0.5*(  DmvY(i,j,k)-DmvY(i,j-1,k)  ) 
    else
 vY=DmvY(i,j,k)+  0.5*(  DpvY(i,j,k)-DpvY(i,j-1,k)  )
    end if 
    
 else

    if (  abs(  DmvY(i,j+1,k)-DmvY(i,j,k) ).lt.abs(  DpvY(i,j+1,k)-DpvY(i,j,k) )   ) then
 vY=DpvY(i,j,k)-  0.5*(  DmvY(i,j+1,k)-DmvY(i,j,k)  ) 
    else
 vY=DpvY(i,j,k)-  0.5*(  DpvY(i,j+1,k)-DpvY(i,j,k)  )
    end if 

 end if 
 
 !!V3
 if (  (0.25*( W(i,j,k)+W(i,j-1,k)+W(i,j,k+1)+W(i,j-1,k+1) )).gt.0.0) then

    if (  abs(  DmvZ(i,j,k)-DmvZ(i,j,k-1) ).lt.abs(  DpvZ(i,j,k)-DpvZ(i,j,k-1) )   ) then
 vZ=DmvZ(i,j,k)+  0.5*(  DmvZ(i,j,k)-DmvZ(i,j,k-1)  ) 
    else
 vZ=DmvZ(i,j,k)+  0.5*(  DpvZ(i,j,k)-DpvZ(i,j,k-1)  )
    end if 
    
 else

    if (  abs(  DmvZ(i,j,k+1)-DmvZ(i,j,k) ).lt.abs(  DpvZ(i,j,k+1)-DpvZ(i,j,k) )   ) then
 vZ=DpvZ(i,j,k)-  0.5*(  DmvZ(i,j,k+1)-DmvZ(i,j,k)  ) 
    else
 vZ=DpvZ(i,j,k)-  0.5*(  DpvZ(i,j,k+1)-DpvZ(i,j,k)  )
    end if 

 end if 

 advectv(i,j,k)=-(0.25*(  u(i,j,k)+u(i+1,j,k)+u(i,j-1,k)+u(i+1,j-1,k)  )*vX+ V(i,j,k)*vY+ &
 &                        0.25*( W(i,j,k)+W(i,j,k+1)+W(i,j-1,k)+W(i,j-1,k+1) )*VZ)

 end do ; end do ; end do 
 
 
 
 do i=0,nx ; do j=0,ny ; do k=1,nz
 DpWx(i,j,k)=( W(i+1,j,k)-W(i,j,k)   )/( x(i+1)-x(i) )
 DmWx(i,j,k)=( W(i,j,k)  -W(i-1,j,k) )/( x(i)-x(i-1) )
 DpWy(i,j,k)=( W(i,j+1,k)-W(i,j,k)   )/( y(j+1)-y(j) )
 DmWy(i,j,k)=( W(i,j,k)  -W(i,j-1,k) )/( y(j)-y(j-1) )
 DpWz(i,j,k)=( W(i,j,k+1)-W(i,j,k)   )/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )
 DmWz(i,j,k)=( W(i,j,k)  -W(i,j,k-1) )/( 0.5*(z(k)  +z(k-1))-0.5*(z(k-1)+z(k-2)) )
 end do ;end do ; end do 


 do i=1,nx-1 ; do j=1,ny-1  ; do k=2,nz-1 !!A!!
 
 
!!W1
 if (  (0.25*( u(i,j,k)+u(i+1,j,k)+u(i,j,k-1)+u(i+1,j,k-1) )).gt.0.0) then

    if (  abs(  DmWx(i,j,k)-DmWx(i-1,j,k) ).lt.abs(  DpWx(i,j,k)-DpWx(i-1,j,k) )   ) then
 Wx=DmWx(i,j,k)+  0.5*(  DmWx(i,j,k)-DmWx(i-1,j,k)  ) 
    else
 Wx=DmWx(i,j,k)+  0.5*(  DpWx(i,j,k)-DpWx(i-1,j,k)  )
    end if 
    
 else

    if (  abs(  DmWx(i+1,j,k)-DmWx(i,j,k) ).lt.abs(  DpWx(i+1,j,k)-DpWx(i,j,k) )   ) then
 Wx=DpWx(i,j,k)-  0.5*(  DmWx(i+1,j,k)-DmWx(i,j,k)  ) 
    else
 Wx=DpWx(i,j,k)-  0.5*(  DpWx(i+1,j,k)-DpWx(i,j,k)  )
    end if 

 end if 



!!W2
 if (  (0.25*( V(i,j,k)+V(i,j+1,k)+v(i,j,k-1)+v(i,j+1,k-1) )).gt.0.0) then

    if (  abs(  DmWY(i,j,k)-DmWY(i,j-1,k) ).lt.abs(  DpWY(i,j,k)-DpWY(i,j-1,k) )   ) then
 WY=DmWY(i,j,k)+  0.5*(  DmWY(i,j,k)-DmWY(i,j-1,k)  ) 
    else
 WY=DmWY(i,j,k)+  0.5*(  DpWY(i,j,k)-DpWY(i,j-1,k)  )
    end if 
    
 else

    if (  abs(  DmWY(i,j+1,k)-DmWY(i,j,k) ).lt.abs(  DpWY(i,j+1,k)-DpWY(i,j,k) )   ) then
 WY=DpWY(i,j,k)-  0.5*(  DmWY(i,j+1,k)-DmWY(i,j,k)  ) 
    else
 WY=DpWY(i,j,k)-  0.5*(  DpWY(i,j+1,k)-DpWY(i,j,k)  )
    end if 

 end if 
 
 !!!W3
  if ( W(i,j,k).gt.0.0) then

    if (  abs(  DmWZ(i,j,k)-DmWZ(i,j,k-1) ).lt.abs(  DpWZ(i,j,k)-DpWZ(i,j,k-1) )   ) then
 WZ=DmWZ(i,j,k)+  0.5*(  DmWZ(i,j,k)-DmWZ(i,j,k-1)  ) 
    else
 WZ=DmWZ(i,j,k)+  0.5*(  DpWZ(i,j,k)-DpWZ(i,j,k-1)  )
    end if 
    
 else

    if (  abs(  DmWZ(i,j,k+1)-DmWZ(i,j,k) ).lt.abs(  DpWZ(i,j,k+1)-DpWZ(i,j,k) )   ) then
 WZ=DpWZ(i,j,k)-  0.5*(  DmWZ(i,j,k+1)-DmWZ(i,j,k)  ) 
    else
 WZ=DpWZ(i,j,k)-  0.5*(  DpWZ(i,j,k+1)-DpWZ(i,j,k)  )
    end if 

 end if 

 advectW(i,j,k)=-(0.25*( u(i,j,k)+u(i+1,j,k)+u(i,j,k-1)+u(i+1,j,k-1) )*WX+ &
 &                0.25*( V(i,j,k)+V(i,j+1,k)+v(i,j,k-1)+v(i,j+1,k-1) )*WY+W(i,j,k)*WZ)
 

 end do ; end do ; end do 


 return 

 end 


 subroutine viscoset(miuv,ro,Tx,Ty,Tz)
 use M_General, only: nx,ny,nz,x,y,z,u,v,w
 implicit none 
 integer i,j,k
 real(8) Txxr,Txxl,Tyxd,Tyxu,Tzxf,Tzxb,Tyyu,Tyyd,Txyr,Txyl,Tzyf,Tzyb,Tzzf,Tzzb,Tyzd,Tyzu,Txzr,Txzl
 real(8) Tx(2:nx-1,1:ny-1,1:nz-1),Ty(1:nx-1,2:ny-1,1:nz-1),Tz(1:nx-1,1:ny-1,2:nz-1)
 real(8),dimension (0:nx,0:ny,0:nz)       ::   ro,miuv


 Do i=2,nx-1 ; Do j=1, ny-1 ; do k=1, nz-1


 
 Txxr=2*miuv(i,j,k)  *( U(i+1,j,K)-U(i,j,K)   )/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )
 
 Txxl=2*miuv(i-1,j,k)*( U(i,j,K)  -U(i-1,j,K) )/( 0.5*(x(i)  +x(i-1))-0.5*(x(i-1)+x(i-2)) )

 Tyxu=0.25*( miuv(i,j,k)+miuv(i,j+1,k)+miuv(i-1,j,k)+miuv(i-1,j+1,k) )*&
   &(      ( U(i,j+1,K)-U(i,j,K)     )/( Y(j+1)-y(j) )+ ( V(i,j+1,K)-V(i-1,j+1,K) )/( x(i)-x(i-1) )     )
                 
 Tyxd=0.25*( miuv(i,j-1,k)+miuv(i,j,k)+miuv(i-1,j-1,k)+miuv(i-1,j,k) )*&
   &(      ( U(i,j,K)-U(i,j-1,K)     )/( y(j)-y(j-1) )+ ( V(i,j,K)  -V(i-1,j,K)   )/( x(i)-x(i-1) )     )

   
 Tzxf=0.25*( miuv(i,j,k)+miuv(i,j,k+1)+miuv(i-1,j,k)+miuv(i-1,j,k+1) )*&
   &(      ( U(i,j,K+1)-U(i,j,K)     )/( z(k+1)-z(k) )+ ( W(i,j,K+1)-W(i-1,j,K+1) )/( x(i)-x(i-1) )     )
  
 
 Tzxb=0.25*( miuv(i,j,k-1)+miuv(i,j,k)+miuv(i-1,j,k-1)+miuv(i-1,j,k) )*&
   &(      ( U(i,j,K)-U(i,j,K-1)     )/( z(k)-z(k-1) )+ ( W(i,j,K)  -W(i-1,j,K)   )/( x(i)-x(i-1) )     )
   
   

 Tx(i,j,k)= (  (Txxr-Txxl)/( x(i)-x(i-1)) + &
  &            (Tyxu-Tyxd)/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+ &
  &            (Tzxf-Tzxb)/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )  )/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )
 
   
end do ; end do ; end do  

Do j=2,ny-1 ; Do i=1,nx-1 ; do k=1, nz-1



Tyyu=2*miuv(i,j,k)   *( V(i,j+1,k)-V(i,j,k)   )/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )
 
 Tyyd=2*miuv(i,j-1,k)*( V(i,j,k)  -V(i,j-1,k) )/( 0.5*(y(j)  +y(j-1))-0.5*(y(j-1)+y(j-2)) )
 
 
 Txyr=0.25*( miuv(i,j,k)+miuv(i,j-1,k)+miuv(i+1,j,k)+miuv(i+1,j-1,k) )*&
  & (     ( V(i+1,j,k)-V(i,j,k)     )/( x(i+1)-x(i) )+ ( U(i+1,j,k)  -U(i+1,j-1,k) )/( y(j)-y(j-1)  )   )
  
 Txyl=0.25*( miuv(i-1,j,k)+miuv(i-1,j-1,k)+miuv(i,j,k)+miuv(i,j-1,k) )*&
  & (     ( V(i,j,k)-V(i-1,j,k)     )/( x(i)-x(i-1) )+ ( U(i,j,k)    -U(i,j-1,k)   )/( y(j)-y(j-1)  )   )
  
 Tzyf=0.25*( miuv(i,j,k)+miuv(i,j-1,k)+miuv(i,j,k+1)+miuv(i,j-1,k+1) )*&
 &  (     ( v(i,j,k+1)-V(i,j,k)     )/( z(k+1)-z(k) )+ (W(i,j,k+1)  -W(i,j-1,k+1)  )/( y(j)-y(j-1)  )   )
 
 Tzyb=0.25*( miuv(i,j,k-1)+miuv(i,j-1,k-1)+miuv(i,j,k)+miuv(i,j-1,k) )*&                         
 &  (     ( v(i,j,k)-V(i,j,k-1)     )/( z(k)-z(k-1) )+ (W(i,j,k)    -W(i,j-1,k)    )/( y(j)-y(j-1)  )   )  
  

 Ty(i,j,k)= (  (Tyyu-Tyyd)/( y(j)-y(j-1) ) + &
            &  (Txyr-Txyl)/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+ &
            &  (Tzyf-Tzyb)/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )  )/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )
 

   
end do ; end do ; end do 

Do j=1,ny-1 ; Do i=1,nx-1 ; do k=2, nz-1


Tzzf=2*miuv(i,j,k)  * ( w(i,j,k+1)-w(i,j,k)   )/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)   +z(k-1)) )
Tzzb=2*miuv(i,j,k-1)* ( w(i,j,k)  -w(i,j,k-1) )/( 0.5*(z(k)  +z(k-1))-0.5*(z(k-1) +z(k-2)) )

Txzr=0.25*( miuv(i,j,k)+miuv(i+1,j,k)+miuv(i,j,k-1)+miuv(i+1,j,k-1) )*&
   &(      ( U(i+1,j,K)-U(i+1,j,K-1)     )/( z(k)-z(k-1) )+ ( W(i+1,j,K)-W(i,j,K) )/( x(i+1)-x(i) )     )
  
Txzl=0.25*( miuv(i-1,j,k)+miuv(i,j,k)+miuv(i-1,j,k-1)+miuv(i,j,k-1) )*&
   &(      ( U(i,j,K)  -U(i,j,K-1)       )/( z(k)-z(k-1) )+ ( W(i,j,K)-W(i-1,j,K) )/( x(i)-x(i-1) )     )
   
Tyzu=0.25*( miuv(i,j,k)+miuv(i,j+1,k)+miuv(i,j,k-1)+miuv(i,j+1,k-1) )*&
   &(      ( V(i,j+1,K) -V(i,j+1,K-1)    )/( z(k)-z(k-1) )+ ( W(i,j+1,K)-W(i,j,K) )/( y(j+1)-y(j) )     )

Tyzd=0.25*( miuv(i,j-1,k)+miuv(i,j,k)+miuv(i,j-1,k-1)+miuv(i,j,k-1) )*&
   &(      ( V(i,j,K)   -V(i,j,K-1)      )/( z(k)-z(k-1) )+ ( W(i,j,K)-W(i,j-1,K) )/( y(j)-y(j-1) )     )
   
Tz(i,j,k)= (   (Tzzf-Tzzb)/(z(k)-z(k-1)) + &
           &   (Txzr-Txzl)/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+ &
           &   (Tyzu-Tyzd)/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )  )/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  )  

end do ; end do ; end do 

return 
end 



subroutine levelset(phiold2,dt) 
use M_General, only:nx,ny,nz,x,y,z,u,v,w,phi

implicit none 

real(8), DIMENSION (0:nx,0:ny,0:nz)       ::Dpphix,Dmphix,Dpphiy,Dmphiy,Dpphiz,Dmphiz
real(8), DIMENSION (1:nx-1,1:ny-1,1:nz-1) ::phix,phiy,phiz,Lphin
real(8),dimension (-1:nx+1,-1:ny+1,-1:nz+1)  ::   phiold2
real(8) Lphis,dt
integer i,j,k,kk


!phiold2(1:nx-1,1:ny-1,1:nz-1)=phi(1:nx-1,1:ny-1,1:nz-1)


do kk=1,2  !!prediction correction method!!


do i=0,nx ; do j=0,ny  ; do  k=0,nz                               
Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i)  )
Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1)  )
Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j)  )
Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1)  )
Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k)  )
Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1)  )
end do ;end do ; end do 


do i=1,nx-1 ; do j=1,ny-1  ; do k=1,nz-1   !!A!!

if (0.5*( u(i,j,k)+u(i+1,j,k) ).gt.0.0) then


  if (  abs(  Dmphix(i,j,k)-Dmphix(i-1,j,k) ).lt.abs(  Dpphix(i,j,k)-Dpphix(i-1,j,k) )   ) then
phix(i,j,k)=Dmphix(i,j,k)+  0.5*(  Dmphix(i,j,k)-Dmphix(i-1,j,k)  ) 
   else
phix(i,j,k)=Dmphix(i,j,k)+  0.5*(  Dpphix(i,j,k)-Dpphix(i-1,j,k)  )
   end if 
    
else

   if (  abs(  Dmphix(i+1,j,k)-Dmphix(i,j,k) ).lt.abs(  Dpphix(i+1,j,k)-Dpphix(i,j,k) )   ) then
phix(i,j,k)=Dpphix(i,j,k)-  0.5*(  Dmphix(i+1,j,k)-Dmphix(i,j,k)  ) 
   else
phix(i,j,k)=Dpphix(i,j,k)-  0.5*(  Dpphix(i+1,j,k)-Dpphix(i,j,k)  )
   end if 

end if 



if (0.5*( V(i,j,k)+V(i,j+1,k) ).gt.0.0) then


   if (  abs(  DmphiY(i,j,k)-DmphiY(i,j-1,k) ).lt.abs(  DpphiY(i,j,k)-DpphiY(i,j-1,k) )   ) then
phiY(i,j,k)=DmphiY(i,j,k)+  0.5*(  DmphiY(i,j,k)-DmphiY(i,j-1,k)  ) 
   else
phiY(i,j,k)=DmphiY(i,j,k)+  0.5*(  DpphiY(i,j,k)-DpphiY(i,j-1,k)  )
   end if 
    
else

  if (  abs(  DmphiY(i,j+1,k)-DmphiY(i,j,k) ).lt.abs(  DpphiY(i,j+1,k)-DpphiY(i,j,k) )   ) then
phiY(i,j,k)=DpphiY(i,j,k)-  0.5*(  DmphiY(i,j+1,k)-DmphiY(i,j,k)  ) 
   else
phiY(i,j,k)=DpphiY(i,j,k)-  0.5*(  DpphiY(i,j+1,k)-DpphiY(i,j,k)  )
   end if 


end if 

if (0.5*( W(i,j,k)+W(i,j,k+1) ).gt.0.0) then

   if (  abs(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1) ).lt.abs(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1) )   ) then
phiZ(i,j,k)=DmphiZ(i,j,k)+  0.5*(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1)  ) 
   else
phiZ(i,j,k)=DmphiZ(i,j,k)+  0.5*(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1)  )
   end if 
    
else

   if (  abs(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k) ).lt.abs(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k) )   ) then
phiZ(i,j,k)=DpphiZ(i,j,k)-  0.5*(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k)  ) 
   else
phiZ(i,j,k)=DpphiZ(i,j,k)-  0.5*(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k)  )
   end if 

end if 


     
end do ;end do ; end do   !!A!!

if (kk.eq.1) then

do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1
Lphin(i,j,k)=  (-0.5)*( u(i,j,k)+u(i+1,j,k) )* phix(i,j,k)+&
           &   (-0.5)*( V(i,j,k)+V(i,j+1,k) )* phiy(i,j,k)+&  
           &   (-0.5)*( W(i,j,k)+W(i,j,k+1) )* phiz(i,j,k)  
phi(i,j,k)=phiold2(i,j,k)+dt*( (-0.5)*( u(i,j,k)+u(i+1,j,k) )* phix(i,j,k)+&
                     &     (-0.5)*( V(i,j,k)+V(i,j+1,k) )* phiy(i,j,k)+&
                     &     (-0.5)*( W(i,j,k)+W(i,j,k+1) )* phiz(i,j,k) )     
end do ;end do ; end do 

                  
end if                        

   if (kk.eq.2) then
   
   do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1
!phi(i,j)=phi(i,j)+dt*( -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
 !                    & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j)  )
 
Lphis=   (-0.5)*( u(i,j,k)+u(i+1,j,k) )* phix(i,j,k)+&
       & (-0.5)*( V(i,j,k)+V(i,j+1,k) )* phiy(i,j,k)+& 
       & (-0.5)*( W(i,j,k)+W(i,j,k+1) )* phiz(i,j,k)       
phi(i,j,k)=phiold2(i,j,k)+0.5*dt*(Lphis+Lphin(i,j,k) )                       
   end do ;end do ; end do 

   end if 


end do  !!prediction corection method!!                      

call boundarycondphi()

RETURN
END 


subroutine reinitialize(eps,dt)
use M_General, only: nx,ny,nz,pi,x,y,z,phi

implicit none
real(8) eps
real(8) phixm,phiym,phizm,phixp,phiyp,phizp,phixr,phiyr,phizr,sphi,s,dtau
real(8), Dimension (0:nx,0:ny,0:nz)       ::  Dpphix,Dmphix,Dpphiy,Dmphiy,Dpphiz,Dmphiz
real(8), Dimension (1:nx-1,1:ny-1,1:nz-1) ::  phiold,phin,lphin
real(8) sgn
real(8) lphis,bb,avgh,dt
integer i,j,k,kk,m


phiold(1:nx-1,1:ny-1,1:nz-1)=phi(1:nx-1,1:ny-1,1:nz-1)
avgh=1000 !! large number!! 
do i=-1,nx ; do j=-1, ny ; do k=-1,nz
bb=min ( (x(i+1)-x(i)),(y(j+1)-y(j)),(z(k+1)-z(k)) )
if (bb.lt.avgh) then
avgh=bb
end if 
end do ; end do ; end do 

dtau=0.5*avgh

do kk=1,3 !? 

do m=1,2   !!runge kutta!!

do i=0,nx ; do j=0,ny ; do k=0,nz
Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i) )
Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1) )
Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j) )
Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1) )
Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k) )
Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1) )
end do ;end do ; end do


do i=1,nx-1 ; do j=1,ny-1  ; do k=1,nz-1 !!A!!  !!A!!


!sphi=phiold(i,j)/(  sqrt(phiold(i,j)*phiold(i,j)+h*h)  )
!sphi=sgn(phiold)




if (phiold(i,j,k).gt.eps) then
sphi=1.0
else if (phiold(i,j,k).lt.-eps) then
sphi=-1.0
else
sphi=phiold(i,j,k)/eps -(1/pi)*sin(pi*phiold(i,j,k)/eps) 
end if 

 

   if (  abs(  Dmphix(i,j,k)-Dmphix(i-1,j,k) ).lt.abs(  Dpphix(i,j,k)-Dpphix(i-1,j,k) )   ) then
phixm=Dmphix(i,j,k)+  0.5*(  Dmphix(i,j,k)-Dmphix(i-1,j,k)  ) 
   else
phixm=Dmphix(i,j,k)+  0.5*(  Dpphix(i,j,k)-Dpphix(i-1,j,k)  )
   end if 
    

    if (  abs(  Dmphix(i+1,j,k)-Dmphix(i,j,k) ).lt.abs(  Dpphix(i+1,j,k)-Dpphix(i,j,k) )   ) then
phixp=Dpphix(i,j,k)-  0.5*(  Dmphix(i+1,j,k)-Dmphix(i,j,k)  ) 
   else
phixp=Dpphix(i,j,k)-  0.5*(  Dpphix(i+1,j,k)-Dpphix(i,j,k)  )
   end if 



    if (  abs(  DmphiY(i,j,k)-DmphiY(i,j-1,k) ).lt.abs(  DpphiY(i,j,k)-DpphiY(i,j-1,k) )   ) then
phiYm=DmphiY(i,j,k)+  0.5*(  DmphiY(i,j,k)-DmphiY(i,j-1,k)  ) 
   else
phiYm=DmphiY(i,j,k)+  0.5*(  DpphiY(i,j,k)-DpphiY(i,j-1,k)  )
   end if 
    


   if (  abs(  DmphiY(i,j+1,k)-DmphiY(i,j,k) ).lt.abs(  DpphiY(i,j+1,k)-DpphiY(i,j,k) )   ) then
phiYp=DpphiY(i,j,k)-  0.5*(  DmphiY(i,j+1,k)-DmphiY(i,j,k)  ) 
   else
phiYp=DpphiY(i,j,k)-  0.5*(  DpphiY(i,j+1,k)-DpphiY(i,j,k)  )
   end if
  
  
   
   if (  abs(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1) ).lt.abs(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1) )   ) then
phiZm=DmphiZ(i,j,k)+  0.5*(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1)  ) 
   else
phiZm=DmphiZ(i,j,k)+  0.5*(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1)  )
   end if 
    

   if (  abs(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k) ).lt.abs(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k) )   ) then
phiZp=DpphiZ(i,j,k)-  0.5*(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k)  ) 
   else
phiZp=DpphiZ(i,j,k)-  0.5*(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k)  )
   end if 

   
        if (sphi*phixp.ge.0.0.AND.sphi*phixm.ge.0.0) then
   phixr=phixm
   else if (sphi*phixp.le.0.0.AND.sphi*phixm.le.0.0) then 
   phixr=phixp
   else if (sphi*phixp.gt.0.0.AND.sphi*phixm.lt.0.0) then
   phixr=0.0
   else if (sphi*phixp.lt.0.0.AND.sphi*phixm.gt.0.0) then
   s=sphi*( abs(phixp)-abs(phixm) )/(phixp-phixm)
                      if (s.gt.0.0) then
                      phixr=phixm
                   else
                      phixr=phixp 
                   end if
   end if
   
   
        if (sphi*phiyp.ge.0.0.AND.sphi*phiym.ge.0.0) then
   phiyr=phiym
   else if (sphi*phiyp.le.0.0.AND.sphi*phiym.le.0.0) then 
   phiyr=phiyp
   else if (sphi*phiyp.gt.0.0.AND.sphi*phiym.lt.0.0) then
   phiyr=0.0
   else if (sphi*phiyp.lt.0.0.AND.sphi*phiym.gt.0.0) then
   s=sphi*( abs(phiyp)-abs(phiym) )/(phiyp-phiym)
                      if (s.gt.0.0) then
                      phiyr=phiym
                   else
                      phiyr=phiyp 
                  end if
    end if
    
        if (sphi*phizp.ge.0.0.AND.sphi*phizm.ge.0.0) then
   phizr=phizm
   else if (sphi*phizp.le.0.0.AND.sphi*phizm.le.0.0) then 
   phizr=phizp
   else if (sphi*phizp.gt.0.0.AND.sphi*phizm.lt.0.0) then
   phizr=0.0
   else if (sphi*phizp.lt.0.0.AND.sphi*phizm.gt.0.0) then
   s=sphi*( abs(phizp)-abs(phizm) )/(phizp-phizm)
                      if (s.gt.0.0) then
                      phizr=phizm
                   else
                      phizr=phizp 
                   end if
   end if
   
   
   if (m.eq.1) then
   lphin(i,j,k)=sphi*(  1.0-sqrt(phiyr*phiyr+phixr*phixr+phizr*phizr)  )
   phin(i,j,k)=phi(i,j,k)
   phi(i,j,k)=phi(i,j,k)+dtau*sphi*(  1.0-sqrt(phiyr*phiyr+phixr*phixr+phizr*phizr)  )
   
   end if
   
   if (m.eq.2) then
   lphis   =sphi*(  1.0-sqrt(phiyr*phiyr+phixr*phixr+phizr*phizr)  )
   phi(i,j,k)=phin(i,j,k)+0.5*dtau*(  lphis+lphin(i,j,k)  )                   
   end if 

end do ; end do ; end do

end do 

 call boundarycondphi()

end do !!fictious time step!!

return 
END



subroutine boundarycondphi()
use M_General, only: nx,ny,nz,phi
implicit none 
integer i,j,k



do j=1,ny-1 ; do k=1,nz-1
phi(0,j,k)=2*phi(1,j,k)-phi(2,j,k)
phi(-1,j,k)=3*phi(1,j,k)-2*phi(2,j,k)

phi(nx,j,k)=2*phi(nx-1,j,k)-phi(nx-2,j,k)
phi(nx+1,j,k)=3*phi(nx-1,j,k)-2*phi(nx-2,j,k)
end do ; end do

do i=-1,nx+1 ; do k=1,nz-1
phi(i,0,k)=2*phi(i,1,k)-phi(i,2,k)
phi(i,-1,k)=3*phi(i,1,k)-2*phi(i,2,k)

phi(i,ny,k)=2*phi(i,ny-1,k)-phi(i,ny-2,k)
phi(i,ny+1,k)=3*phi(i,ny-1,k)-2*phi(i,ny-2,k)
end do ; end do 

do i=-1, nx+1 ; do j=-1 , ny+1 
phi(i,j,0)=2*phi(i,j,1)-phi(i,j,2)
phi(i,j,-1)=3*phi(i,j,1)-2*phi(i,j,2)

phi(i,j,nz)=2*phi(i,j,nz-1)-phi(i,j,nz-2)
phi(i,j,nz+1)=3*phi(i,j,nz-1)-2*phi(i,j,nz-2)
end do ; end do

return 
end  



   
real(8) function HV(phi,eps)
IMPLICIT  NONE

                 real(8) phi,pi,eps


                 pi=3.1415926540
                 if (phi.gt.eps) then
                   HV=1.0
                 else if (phi.lt.-eps) then
                   HV=0.0
                 else
                   HV=0.50*(  1.0+ phi/eps +(1/pi)*sin(pi*phi/eps)  )
                 end if 

return
end 





subroutine pCofent(ro,dt,Cof)
 use M_General, only: nx,ny,nz,x,y,z,u,v,w
 implicit none 
 
 
real(8),dimension(1:nx-1,1:ny-1,1:nz-1,8) ,intent(out)     ::   Cof
real(8),dimension(0:nx,0:ny,0:nz)         ,intent(in)      ::   ro
real(8),                                   intent(in)      ::   dt


!! local variablele 
real(8),dimension (1:nx-1,1:ny-1,1:nz-1)    ::   Apx,Amx,Apy,Amy,Apz,Amz,AP,Q,p
integer i,j,k


 do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1

 if (i==nx-1) then
 Apx(i,j,k)=0
 else
 !Apx(i,j,k)=1/(hx(i)*hx(i+1))/(  0.5*(ro(i,j,k)+ro(i+1,j,k))  )
 Apx(i,j,k)=1/( (x(i+1)-x(i))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i+1,j,k))  )
 end if 
 
 if (i==1)then 
 Amx(i,j,k)=0       
 else 
 !Amx(i,j,k)=1/(hx(i)*hx(i-1))/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )
 Amx(i,j,k)=1/( (x(i)-x(i-1))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )
 end if 
 
 if (j==ny-1) then
 Apy(i,j,k)=0
 else
 !Apy(i,j,k)=1/(hy(j)*hy(j+1))/(  0.5*(ro(i,j,k)+ro(i,j+1,k))  )
 Apy(i,j,k)=1/( (y(j+1)-y(j))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j+1,k))  )

 end if 
 if (j==1) then
 Amy(i,j,k)=0
 else
 !Amy(i,j,k)=1/(hy(j)*hy(j-1))/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )
 Amy(i,j,k)=1/( (y(j)-y(j-1))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )
 
 end if 
 if (k==nz-1) then
 Apz(i,j,k)=0
 else
 !Apz(i,j,k)=1/(hz(k)*hz(k+1))/(  0.5*(ro(i,j,k)+ro(i,j,k+1))  )
 Apz(i,j,k)=1/( (z(k+1)-z(k))*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j,k+1))  )  
 end if 
 
 if (k==1) then
 Amz(i,j,k)=0   
 else
 !Amz(i,j,k)=1/(hz(k)*hz(k-1))/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  ) 
 Amz(i,j,k)=1/( (z(k)-z(k-1))*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  )
 end if 

 AP(i,j,k)=-( Apx(i,j,k)+Amx(i,j,k)+Apy(i,j,k)+Amy(i,j,k)+Apz(i,j,k)+Amz(i,j,k) ) 
 
 
Q(I,J,K)=(  (U(I+1,J,k)-U(I,J,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
          &  (V(I,J+1,k)-V(I,J,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+& 
          &  (W(I,j,k+1)-W(I,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )   )/dt

 end do ; end do ; end do
 
 
 Cof(:,:,:,1)=Amx(:,:,:)
 Cof(:,:,:,2)=Apx(:,:,:)
 Cof(:,:,:,3)=Amy(:,:,:)
 Cof(:,:,:,4)=Apy(:,:,:)
 Cof(:,:,:,5)=Amz(:,:,:)
 Cof(:,:,:,6)=Apz(:,:,:)
 Cof(:,:,:,7)=Ap (:,:,:)
 Cof(:,:,:,8)=Q  (:,:,:)
 
          

   return 

   end 
   
   
  subroutine SORpoison(beta,pdif,p,maxSOR,COF,test)
 use M_General,only: nx,ny,nz

 implicit none 
 
 integer  n
 parameter ( n=(nx-1)*(ny-1)*(nz-1) ) 

 real(8)                                  ,intent(in)      ::  pdif,beta
 real(8),dimension(1:nx-1,1:ny-1,1:nz-1)  ,intent(inout)   ::  p
 real(8),dimension(1:nx-1,1:ny-1,1:nz-1,8),intent(in)      ::  Cof
 integer                                  ,intent(inout)   ::  maxSOR
 Logical, intent(in)                                       ::  test
 !! local variables !!
 real(8),dimension(-nx*ny:n+nx*ny)          ::   pp
 real(8),dimension (n)                      ::   rhs,xx
 real(8),dimension(7,n)                     ::   AA
 real(8) maxp
 integer kk,i,j,k,number2
 

if (test) then 
maxSOR=1000
end if 

  kk=0
  do k=1,nz-1 
  do j=1,ny-1 
  do i=1,nx-1 
  kk=kk+1
  AA(1,kk)=Cof(i,j,k,1)
  AA(2,kk)=Cof(i,j,k,2)
  AA(3,kk)=Cof(i,j,k,3)
  AA(4,kk)=Cof(i,j,k,4)
  AA(5,kk)=Cof(i,j,k,5)
  AA(6,kk)=Cof(i,j,k,6)
  AA(7,kk)=Cof(i,j,k,7)
  rhs(kk) =Cof(i,j,k,8)
  xx(i)= p(i,j,k)
  end do ; end do ; end do
  
   print *, "SOR in" 
    
  pp=0
  do i=1,n
   pp(i)=xx(i)
  end do
  
  
  number2=0 ; maxp=100.0
    
    do while (maxp.gt.pdif.AND.number2.lt.maxSOR) 
       maxp=0 ; number2=number2+1
 
     do i=1,n 

    xx(i)=pp(i)
    pp(i)=beta*( ( AA(1,i)*pp(i-1)+AA(2,i)*pp(i+1)+AA(3,i)*pp(i-(nx-1))+AA(4,i)*pp(i+(nx-1)) &
    &             +AA(5,i)*pp(i-((nx-1)*(ny-1))) &
    &             +AA(6,i)*pp(i+((nx-1)*(ny-1)))-rhs(i) )/(-AA(7,i)) )+(1-beta)*pp(i)
         
           if (abs( xx(i)-pp(i) ).gt.maxp) then 
           maxp=abs( xx(i)-pp(i) )
           end if 
 
     end do 
 
    end do  !!while !
            
    print *,"toler=",maxp,"number of iteration=",number2
    
  
  print*, "SOR out "  
    
           
    kk=0     
    do k=1,nz-1 
    do j=1,ny-1 
    do i=1,nx-1 
    kk=kk+1
    p(i,j,k)=real(xx(kk))
    end do
    end do 
    end do 
     
            

   return 

   end 



    
    
    
    subroutine hyprepoison(maxIteration,maxError,p,cof)
    use M_General, only: nx,ny,nz
!    use module_poisson
   
  implicit none 

  
 real(8),dimension(1:nx-1,1:ny-1,1:nz-1)   ,intent(inout)  ::  p
 real(8),dimension(1:nx-1,1:ny-1,1:nz-1,8) ,intent(in)     ::  Cof
 integer,intent(inout) :: maxIteration
 real(8),intent(inout) :: maxError
 real(8),dimension(1:nz-1,1:nx-1,1:ny-1)  :: Pnew
  integer :: is,ie,js,je,ks,ke,numIteration
  integer :: i,j,k 
  real(8), allocatable, dimension(:,:,:,:) ::   A
  

  is=1 ; ie=nz-1 
  js=1 ; je=nx-1
  ks=1 ; ke=ny-1
 
allocate( A(is:ie,js:je,ks:ke,8) )
  do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 
  A(k,i,j,1)= cof(i,j,k,5)
  A(k,i,j,2)= cof(i,j,k,6)
  A(k,i,j,3)= cof(i,j,k,1)
  A(k,i,j,4)= cof(i,j,k,2)
  A(k,i,j,5)= cof(i,j,k,3)
  A(k,i,j,6)= cof(i,j,k,4)
  A(k,i,j,7)=-cof(i,j,k,7)
  A(k,i,j,8)=-cof(i,j,k,8)
 end do ; end do ; end do 


 print*, "hypre in"
 call flush(6)
 do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 
 Pnew(k,i,j)=p(i,j,k)
 end do ; end do ; end do  
 ! call solve(A,Pnew,maxError,maxIteration,numIteration,is,ie,js,je,ks,ke)
 do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 
 p(i,j,k)=Pnew(k,i,j)
 end do ; end do ; end do 

write(*,1300)numIteration,maxval(p)
  1300 FORMAT('Solved in ',I3, ' iterations, Pmax=',e15.6e2)
  
 print*, "hypre out"
call flush (6) 
return 
   end 




subroutine contactangle(uext,vext,wext,gap,test)


use M_General,          only : nx,ny,nz,x,y,z,phi
use M_solid,            only :Ib_Solid
implicit none

real(8),intent(in)                                       ::   gap
Logical,intent(in)                                       ::   test
real(8),intent(out),dimension (0:nx+1,-1:ny+1,-1:nz+1)   ::   uext
real(8),intent(out),dimension (-1:nx+1,0:ny+1,-1:nz+1)   ::   vext
real(8),intent(out),dimension (-1:nx+1,-1:ny+1,0:nz+1)   ::   wext



!! local variable !! 
real(8),dimension (0:nx,0:ny,0:nz)           ::    divext,divextx,divexty,divextz
real(8),dimension (1:nx-1,1:ny-1,1:nz-1)     ::    Iphi
real(8),dimension (0:nx,0:ny,0:nz)           ::    Ibf
integer i,j,k
!!!!

Ibf=0
do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 
Ibf(i,j,k)= 0.508*(Ib_Solid(i,j,k))+ &
0.026*( Ib_Solid(i,j,k-1)+Ib_Solid(i,j,k+1)+Ib_Solid(i,j-1,k)+Ib_Solid(i,j+1,k)+Ib_Solid(i-1,j,k)+Ib_Solid(i+1,j,k) )+ &
0.018*( Ib_Solid(i,j-1,k-1)+Ib_Solid(i,j-1,k+1)+Ib_Solid(i,j+1,k-1)+Ib_Solid(i,j+1,k+1)  + &
        Ib_Solid(i-1,j,k-1)+Ib_Solid(i-1,j,k+1)+Ib_Solid(i+1,j,k-1)+Ib_Solid(i+1,j,k+1)  + &
        Ib_Solid(i-1,j-1,k)+Ib_Solid(i-1,j+1,k)+Ib_Solid(i+1,j-1,k)+Ib_Solid(i+1,j+1,k) )+ &
0.015*( Ib_Solid(i-1,j-1,k-1)+Ib_Solid(i-1,j-1,k+1)+ &
        Ib_Solid(i+1,j+1,k-1)+Ib_Solid(i+1,j+1,k+1)+ &
        Ib_Solid(i-1,j+1,k-1)+Ib_Solid(i-1,j+1,k+1)+ &
        Ib_Solid(i+1,j-1,k-1)+Ib_Solid(i+1,j-1,k+1) )

end do ; end do ; end do 


            
uext=0 ; vext=0 ; wext=0 ; divext=0 ; divextx=0 ; divexty=0 ; uext=0 ; vext=0 
do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 


Iphi(i,j,k)=(0.5d0+0.5d0*( phi(i,j,k)**3 + 1.5d0* gap*gap *phi(i,j,k) )  / &
 &          ( phi(i,j,k)*phi(i,j,k)+gap*gap )**1.5d0)
!if ( Ibf(i,j)    .ne.0.0.AND.Ibf(i+1,j)  .ne.0.0.AND.Ibf(i-1,j)  .ne.0.0.AND.Ibf(i,j-1)  .ne.0.0.AND.Ibf(i,j+1).ne.0.0.AND.& 
!&    Ibf(i+1,j+1).ne.0.0.AND.Ibf(i+1,j-1).ne.0.0.AND.Ibf(i-1,j+1).ne.0.0.And.Ibf(i-1,j-1).ne.0.0.AND. &
!&    Ibf(i,j)    .ne.1.0.AND.Ibf(i+1,j)  .ne.1.0.AND.Ibf(i-1,j)  .ne.1.0.AND.Ibf(i,j-1)  .ne.1.0.AND.Ibf(i,j+1).ne.1.0.AND.& 
!&    Ibf(i+1,j+1).ne.1.0.AND.Ibf(i+1,j-1).ne.1.0.AND.Ibf(i-1,j+1).ne.1.0.And.Ibf(i-1,j-1).ne.1.0) then

if (       Ibf(i,j,k)   .ne.0.0  .AND. Ibf(i,j,k)  .ne.1.0 &
 &    .AND.Ibf(i-1,j,k) .ne.0.0  .AND. Ibf(i-1,j,k).ne.1.0 &
 &    .AND.Ibf(i,j-1,k) .ne.0.0  .AND. Ibf(i,j-1,k).ne.1.0 &
 &    .AND.Ibf(i,j,k-1) .ne.0.0  .AND. Ibf(i,j,k-1).ne.1.0 &
 &    .AND.Iphi(i,j,k) .gt.0.001.AND. Iphi(i,j,k).lt.0.999  ) then 

 

divextx(i,j,k)=sqrt ( &
  ( ( Ibf(i,j,k)-Ibf(i-1,j,k) )/( x(i)-x(i-1) ) )**2 + &
  
  ( ( 0.25*(Ibf(i,j,k)+Ibf(i-1,j,k)  +Ibf(i,j+1,k)+Ibf(i-1,j+1,k)) &
     -0.25*(Ibf(i,j,k)+Ibf(i-1,j,k)  +Ibf(i,j-1,k)+Ibf(i-1,j-1,k))   )/( 0.5*(y(j+1)-y(j-1)) ) )**2  &      
+ ( ( 0.25*(Ibf(i,j,k)+Ibf(i-1,j,k)  +Ibf(i,j,k+1)+Ibf(i-1,j,k+1)) &
     -0.25*(Ibf(i,j,k)+Ibf(i-1,j,k)  +Ibf(i,j,k-1)+Ibf(i-1,j,k-1))   )/( 0.5*(z(k+1)-z(k-1)) ) )**2  )

divexty(i,j,k)=sqrt ( &
  ( ( Ibf(i,j,k)-Ibf(i,j-1,k))/(y(j)-y(j-1)) )**2  +  &
  
  ( ( 0.25*( Ibf(i,j,k)+Ibf(i,j-1,k)  +Ibf(i+1,j,k)+Ibf(i+1,j-1,k)) &
     -0.25*( Ibf(i,j,k)+Ibf(i,j-1,k)  +Ibf(i-1,j,k)+Ibf(i-1,j-1,k))   )/( 0.5*(x(i+1)-x(i-1)) ) )**2 &
+ ( ( 0.25*( Ibf(i,j,k)+Ibf(i,j-1,k)  +Ibf(i,j,k+1)+Ibf(i,j-1,k+1)) &
     -0.25*( Ibf(i,j,k)+Ibf(i,j-1,k)  +Ibf(i,j,k-1)+Ibf(i,j-1,k-1))   )/( 0.5*(z(k+1)-z(k-1)) ) )**2  )
    
divextz(i,j,k)=sqrt ( &
  ( ( Ibf(i,j,k)-Ibf(i,j,k-1))/(z(k)-z(k-1)) )**2   +  &
 
  ( (  0.25*(  Ibf(i,j,k)+Ibf(i,j,k-1) +Ibf(i+1,j,k)+Ibf(i+1,j,k-1)) &
      -0.25*(  Ibf(i,j,k)+Ibf(i,j,k-1) +Ibf(i-1,j,k)+Ibf(i-1,j,k-1))   )/( 0.5*(x(i+1)-x(i-1)) ) )**2  &   
+ ( (  0.25*(  Ibf(i,j,k)+Ibf(i,j,k-1) +Ibf(i,j+1,k)+Ibf(i,j+1,k-1)) &
      -0.25*(  Ibf(i,j,k)+Ibf(i,j,k-1) +Ibf(i,j-1,k)+Ibf(i,j-1,k-1))   )/( 0.5*(y(j+1)-y(j-1)) ) )**2   ) 


    

 

 
!if (divextx(i,j).gt.0.0001.AND.divexty(i,j).gt.0.0001) then 
uext(i,j,k)=8.0*( (0.5-0.5*(Ibf(i,j,k)+Ibf(i-1,j,k)))**3 )*(Ibf(i,j,k)-Ibf(i-1,j,k))/( x(i)-x(i-1) )/( divextx(i,j,k)+0.001**2)
vext(i,j,k)=8.0*( (0.5-0.5*(Ibf(i,j,k)+Ibf(i,j-1,k)))**3 )*(Ibf(i,j,k)-Ibf(i,j-1,k))/( y(j)-y(j-1) )/( divexty(i,j,k)+0.001**2) 
wext(i,j,k)=8.0*( (0.5-0.5*(Ibf(i,j,k)+Ibf(i,j,k-1)))**3 )*(Ibf(i,j,k)-Ibf(i,j,k-1))/( z(k)-z(k-1) )/( divextz(i,j,k)+0.001**2) 
else 
uext(i,j,k)=0
vext(i,j,k)=0
wext(i,j,k)=0
end if 

end do ; end do ; end do 

if (test) then 

OPEN(275,file='contactangle.plt')
write(275,*) 'variables="x","y","z","uext","vext","wext","divextx","divexty","divextz","ibf","ib"'
   
write(275,*) 'zone i=',nx-1,' j=',ny-1, 'k=',nz-1 
do k=1,nz-1 ; do j=1,ny-1 ;  Do i=1,nx-1

write(275,550) x(i),y(j),z(k),0.5*(uext(i,j,k)+uext(i+1,j,k)), &
&      0.5*(vext(i,j,k)+vext(i,j+1,k)),0.5*(wext(i,j,k)+wext(i,j,k+1)), &
&      0.5*(divextx(i,j,k)+divextx(i+1,j,k)),0.5*(divexty(i,j,k)+divexty(i,j+1,k)), &
&      0.5*(divextz(i,j,k)+divextz(i,j,k+1)),Ibf(i,j,k),Ib_Solid(i,j,k)
end do ; end do ; end do 
550 format (11(1x,e15.7)) 

end if 



 return 
 end 
 
 
 subroutine levelset2(uext,vext,wext,dt,test) 
use M_General,          only: nx,ny,nz,x,y,z,phi
use M_Solid,            only: Ib_Solid  !! only for plotting here

implicit none 

real(8), DIMENSION (0:nx,0:ny,0:nz)         ::Dpphix,Dmphix,Dpphiy,Dmphiy,Dpphiz,Dmphiz
real(8), DIMENSION (1:nx-1,1:ny-1,1:nz-1)   ::phix,phiy,phiz,Lphin,phiold
real(8) Lphis,dt
integer i,j,k,kk,LL
real(8),dimension (0:nx+1,-1:ny+1,-1:nz+1)  ::   uext
real(8),dimension (-1:nx+1,0:ny+1,-1:nz+1)  ::   vext
real(8),dimension (-1:nx+1,-1:ny+1,0:nz+1)  ::   wext
logical test


do LL=1,4 

phiold(1:nx-1,1:ny-1,1:nz-1)=phi(1:nx-1,1:ny-1,1:nz-1)


do kk=1,2  !!prediction correction method!!


do i=0,nx ; do j=0,ny  ; do  k=0,nz                               
Dpphix(i,j,k)=( phi(i+1,j,k)-phi(i,j,k)   )/( x(i+1)-x(i)  )
Dmphix(i,j,k)=( phi(i,j,k)  -phi(i-1,j,k) )/( x(i)-x(i-1)  )
Dpphiy(i,j,k)=( phi(i,j+1,k)-phi(i,j,k)   )/( y(j+1)-y(j)  )
Dmphiy(i,j,k)=( phi(i,j,k)  -phi(i,j-1,k) )/( y(j)-y(j-1)  )
Dpphiz(i,j,k)=( phi(i,j,k+1)-phi(i,j,k)   )/( z(k+1)-z(k)  )
Dmphiz(i,j,k)=( phi(i,j,k)  -phi(i,j,k-1) )/( z(k)-z(k-1)  )
end do ;end do ; end do 


do i=1,nx-1 ; do j=1,ny-1  ; do k=1,nz-1   !!A!!

if (0.5*( uext(i,j,k)+uext(i+1,j,k) ).gt.0.0) then


  if (  abs(  Dmphix(i,j,k)-Dmphix(i-1,j,k) ).lt.abs(  Dpphix(i,j,k)-Dpphix(i-1,j,k) )   ) then
phix(i,j,k)=Dmphix(i,j,k)+  0.5*(  Dmphix(i,j,k)-Dmphix(i-1,j,k)  ) 
   else
phix(i,j,k)=Dmphix(i,j,k)+  0.5*(  Dpphix(i,j,k)-Dpphix(i-1,j,k)  )
   end if 
    
else

   if (  abs(  Dmphix(i+1,j,k)-Dmphix(i,j,k) ).lt.abs(  Dpphix(i+1,j,k)-Dpphix(i,j,k) )   ) then
phix(i,j,k)=Dpphix(i,j,k)-  0.5*(  Dmphix(i+1,j,k)-Dmphix(i,j,k)  ) 
   else
phix(i,j,k)=Dpphix(i,j,k)-  0.5*(  Dpphix(i+1,j,k)-Dpphix(i,j,k)  )
   end if 

end if 



if (0.5*( Vext(i,j,k)+Vext(i,j+1,k) ).gt.0.0) then


   if (  abs(  DmphiY(i,j,k)-DmphiY(i,j-1,k) ).lt.abs(  DpphiY(i,j,k)-DpphiY(i,j-1,k) )   ) then
phiY(i,j,k)=DmphiY(i,j,k)+  0.5*(  DmphiY(i,j,k)-DmphiY(i,j-1,k)  ) 
   else
phiY(i,j,k)=DmphiY(i,j,k)+  0.5*(  DpphiY(i,j,k)-DpphiY(i,j-1,k)  )
   end if 
    
else

  if (  abs(  DmphiY(i,j+1,k)-DmphiY(i,j,k) ).lt.abs(  DpphiY(i,j+1,k)-DpphiY(i,j,k) )   ) then
phiY(i,j,k)=DpphiY(i,j,k)-  0.5*(  DmphiY(i,j+1,k)-DmphiY(i,j,k)  ) 
   else
phiY(i,j,k)=DpphiY(i,j,k)-  0.5*(  DpphiY(i,j+1,k)-DpphiY(i,j,k)  )
   end if 


end if 

if (0.5*( Wext(i,j,k)+Wext(i,j,k+1) ).gt.0.0) then

   if (  abs(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1) ).lt.abs(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1) )   ) then
phiZ(i,j,k)=DmphiZ(i,j,k)+  0.5*(  DmphiZ(i,j,k)-DmphiZ(i,j,k-1)  ) 
   else
phiZ(i,j,k)=DmphiZ(i,j,k)+  0.5*(  DpphiZ(i,j,k)-DpphiZ(i,j,k-1)  )
   end if 
    
else

   if (  abs(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k) ).lt.abs(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k) )   ) then
phiZ(i,j,k)=DpphiZ(i,j,k)-  0.5*(  DmphiZ(i,j,k+1)-DmphiZ(i,j,k)  ) 
   else
phiZ(i,j,k)=DpphiZ(i,j,k)-  0.5*(  DpphiZ(i,j,k+1)-DpphiZ(i,j,k)  )
   end if 

end if 


     
end do ;end do ; end do   !!A!!

if (kk.eq.1) then

do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1
Lphin(i,j,k)=  (-0.5)*( uext(i,j,k)+uext(i+1,j,k) )* phix(i,j,k)+&
           &   (-0.5)*( Vext(i,j,k)+Vext(i,j+1,k) )* phiy(i,j,k)+&  
           &   (-0.5)*( Wext(i,j,k)+Wext(i,j,k+1) )* phiz(i,j,k)  
phi(i,j,k)=phi(i,j,k)+dt*( (-0.5)*( uext(i,j,k)+uext(i+1,j,k) )* phix(i,j,k)+&
                     &     (-0.5)*( Vext(i,j,k)+Vext(i,j+1,k) )* phiy(i,j,k)+&
                     &     (-0.5)*( Wext(i,j,k)+Wext(i,j,k+1) )* phiz(i,j,k) )     
end do ;end do ; end do 

                  
end if                        

   if (kk.eq.2) then
   
   do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1

 
Lphis=   (-0.5)*( uext(i,j,k)+uext(i+1,j,k) )* phix(i,j,k)+&
       & (-0.5)*( Vext(i,j,k)+Vext(i,j+1,k) )* phiy(i,j,k)+& 
       & (-0.5)*( Wext(i,j,k)+Wext(i,j,k+1) )* phiz(i,j,k)       
phi(i,j,k)=phiold(i,j,k)+0.5*dt*(Lphis+Lphin(i,j,k) )                       
   end do ;end do ; end do 

   end if 


end do  !!prediction corection method!!                      

 call boundarycondphi()

if (test) then 
 
OPEN(135,file='Densitydis.plt')
write(135,*) 'variables="x","y","z","uext","vext","wext","phi","ib"'
   write(135,*) 'zone i=',nx-1,' j=',ny-1 ,' k=',nz-1
do k=1,nz-1 ; do j=1,ny-1 ;  Do i=1,nx-1

write(135,560) x(i),y(j),z(k),uext(i,j,k),vext(i,j,k),wext(i,j,k),phi(i,j,k),Ib_Solid(i,j,k)
end do ; end do ; end do 
560 format (8(1x,e15.7))
end if



end do 

 

RETURN
END 


  subroutine Proprety(pro,prodrop,proair,prosolid,procon,gap,eps)
    use M_General,             only: nx,ny,nz,phi
    use M_Solid,               only: Ib_Solid,Ic_Solid
    implicit none 

    integer i,j,k 
    double precision,intent(in)                       ::   eps 
    dOUBLE PRECISION gapN,gap,prodrop,prosolid,proair,procon
    double precision,dimension (0:nx,0:ny,0:nz)       ::   Iphi,pro

    real(8)                                     ::   HV

    gapN=3.33*gap
    do i=0,nx
      do j=0,ny
        do k=0 ,nz

          !Iphi(i,j,k)=(0.5d0+0.5d0*( phi(i,j,k)**3 + 1.5d0*GapN*GapN*phi(i,j,k) )  / ( phi(i,j,k)*phi(i,j,k)+GapN*GapN)**1.5d0)

          !if (Iphi(i,j,k).lt.0.001) then
          !  pro(i,j,k)=prodrop
          !else if (Iphi(i,j,k).gt.0.999) then
          !  pro(i,j,k)=proair
          !else
          !  pro(i,j,k)=prodrop+Iphi(i,j,k)*(proair-prodrop)
          !end if
          pro(i,j,k)=proair+HV( -phi(i,j,k),eps )*(prodrop-proair)
          pro(i,j,k)=pro(i,j,k)+Ib_Solid(i,j,k)*( prosolid-pro(i,j,k) )
          pro(i,j,k)=pro(i,j,k)+Ic_Solid(i,j,k)*( procon )

        end do
      end do
    end do
 


    return 

  end subroutine    

subroutine poisson (ro,beta,pdif,p,dt)
 use M_General,           only: nx,ny,nz,x,y,z,u,v,w
 use M_Solid,             only: xbar,ybar,zbar,ubar,vbar,wbar,omegax,omegay,omegaz,Ib_Solid
 implicit none 
 
 

 real(8),dimension (1:nx-1,1:ny-1,1:nz-1)    ::   Apx,Amx,Apy,Amy,Apz,Amz,AP,Q,QF
 real(8),dimension (0:nx,0:ny,0:nz)          ::   ro,p,pold,Ibpx,Ibmx,Ibpy,Ibmy,Ibpz,Ibmz
 real(8) pdif,beta,maxp,dt,sump,meanp  
 integer i,j,k,number2



Ibpx=0
Ibmx=0
Ibpy=0
Ibmy=0
Ibpz=0
Ibmz=0

!do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 
!Ibmx(i,j,k)= 0.5*(Ib(i,j,k)+Ib(i-1,j,k))
!Ibpx(i,j,k)= 0.5*(Ib(i,j,k)+Ib(i+1,j,k))
!Ibmy(i,j,k)= 0.5*(Ib(i,j,k)+Ib(i,j-1,k))
!Ibpy(i,j,k)= 0.5*(Ib(i,j,k)+Ib(i,j+1,k))
!Ibpz(i,j,k)= 0.5*(Ib(i,j,k)+Ib(i,j,k+1))
!Ibmz(i,j,k)= 0.5*(Ib(i,j,k)+Ib(i,j,k-1))
!end do ; end do ; end do




!do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 

!QF(i,j,k)=(0.5*(u(i,j,k)+u(i+1,j,k))-( ubar + omegay*(z(k)-zbar)  -omegaz*(y(j)-ybar) ))&
!         *(-Ibpx(i,j,k)+Ibmx(i,j,k))/( 0.5*(x(i+1)-x(i-1)) )+                           &
!&         (0.5*(v(i,j,k)+v(i,j+1,k))-( vbar - omegax*(z(k)-zbar)  +omegaz*(x(i)-xbar) ))&
!&        *(-Ibpy(i,j,k)+Ibmy(i,j,k))/( 0.5*(y(j+1)-y(j-1)) )+                           &
!&         (0.5*(w(i,j,k)+w(i,j,k+1))-( wbar + omegax*(y(j)-ybar)  -omegay*(x(i)-xbar) ))&
!&        *(-Ibpz(i,j,k)+Ibmz(i,j,k))/( 0.5*(z(k+1)-z(k-1)) ) 

!end do ; end do ; end do
 do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1

 if (i==nx-1) then
 Apx(i,j,k)=0
 else
 
 !Apx(i,j,k)=1/( (x(i+1)-x(i))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i+1,j,k))  )
 Apx(i,j,k)=(1-Ibpx(i,j,k))/( (x(i+1)-x(i))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i+1,j,k))  )
 end if 
 
 if (i==1)then 
 Amx(i,j,k)=0       
 else 

 Amx(i,j,k)=(1-Ibmx(i,j,k))/( (x(i)-x(i-1))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )
 end if 
 
 if (j==ny-1) then
 Apy(i,j,k)=0
 else
 
 Apy(i,j,k)=(1-Ibpy(i,j,k))/( (y(j+1)-y(j))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j+1,k))  )

 end if 
 if (j==1) then
 Amy(i,j,k)=0
 else
 
 Amy(i,j,k)=(1-Ibmy(i,j,k))/( (y(j)-y(j-1))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )
 
 end if 
 if (k==nz-1) then
 Apz(i,j,k)=0
 else

 Apz(i,j,k)=(1-Ibpz(i,j,k))/( (z(k+1)-z(k))*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j,k+1))  )  
 end if 
 
 if (k==1) then
 Amz(i,j,k)=0   
 else
 
 Amz(i,j,k)=(1-Ibmz(i,j,k))/( (z(k)-z(k-1))*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  )
 end if 

 AP(i,j,k)=-(   Apx(i,j,k)+ &
              & Amx(i,j,k)+ &
              & Apy(i,j,k)+ &
              & Amy(i,j,k)+ &
              & Apz(i,j,k)+ &
              & Amz(i,j,k) ) 

 !AP(i,j,k)=-(   (1-Ib(i,j,k))/( (x(i+1)-x(i))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i+1,j,k))  )  + &
 !             & (1-Ib(i,j,k))/( (x(i)-x(i-1))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )  + &
 !             & (1-Ib(i,j,k))/( (y(j+1)-y(j))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j+1,k))  )  + &
 !             & (1-Ib(i,j,k))/( (y(j)-y(j-1))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j-1,k))  )  + &
 !             & (1-Ib(i,j,k))/( (z(k+1)-z(k))*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j,k+1))  )  + &
 !             & (1-Ib(i,j,k))/( (z(k)-z(k-1))*( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) ) )/(  0.5*(ro(i,j,k)+ro(i,j,k-1))  )   ) 
 
  
!Q(I,J,K)=(1-Ib(i,j,k))*(   (U(I+1,J,k)-U(I,J,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
!          &  (V(I,J+1,k)-V(I,J,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+& 
!          &  (W(I,j,k+1)-W(I,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )   )/dt + QF(i,j,k)/dt

Q(I,J,K)=(   (U(I+1,J,k)-U(I,J,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
          &  (V(I,J+1,k)-V(I,J,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+& 
          &  (W(I,j,k+1)-W(I,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )   )/dt 


!Q(I,J,K)=(   ((1-Ibpx(i,j,k))*U(I+1,J,k)-(1-Ibmx(i,j,k))*U(I,J,k))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
!          &  ((1-Ibpy(i,j,k))*V(I,J+1,k)-(1-Ibmy(i,j,k))*V(I,J,k))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )+& 
!          &  ((1-Ibpz(i,j,k))*W(I,j,k+1)-(1-Ibmz(i,j,k))*W(I,j,k))/( 0.5*(z(k+1)+z(k)  )-0.5*(z(k)  +z(k-1)) )   )/dt

 end do ; end do ; end do
 
 
 
number2=0
maxp=1000.0
print*, "I am in poisson"
           do while (maxp.gt.pdif.AND.number2.lt.2500) 
 maxp=0 
 number2=number2+1
 
 do i=1,nx-1 ; do j=1,ny-1 ;do k=1, nz-1

!if ( Ap(i,j,k).eq.0 ) then 
!p(i,j,k)=0 
!else    
 

   pold(i,j,k)=p(i,j,k)
   p(i,j,k)=beta*( ( Apx(i,j,k)*P(I+1,j,k)+Amx(i,j,k)*P(I-1,j,k)+Apy(i,j,k)*P(I,j+1,k) &
   &                +Amy(i,j,k)*P(I,j-1,k)+Apz(i,j,k)*P(I,j,k+1)+Amz(i,j,k)*P(I,j,k-1)-Q(I,j,k)  )/(-AP(i,j,k))  ) &
         & +(1-beta)*p(i,j,k)
   if (abs( pold(i,j,k)-p(i,j,k) ).gt.maxp) then 
      maxp=abs( pold(i,j,k)-p(i,j,k) )
   end if

!end if    
 
! sump=sump+p(i,j,k) 
 end do ; end do ; end do 

!meanp=sump/(real(nx-1)*real(ny-1)*real(nz-1))
!p(1:nx-1,1:ny-1,1:nz-1)=p(1:nx-1,1:ny-1,1:nz-1)-meanp
!sump=0
            end do  !!while !
                        
 
   print *,maxp,number2 
 


             

   return 

   end 




subroutine Fluid_Dynamic_Write(p)

   use M_General, only: u,v,w,p,phi,nx,ny,nz,fileplace
   implicit none
   real(8),dimension (0:nx,0:ny,0:nz         ) ::   P
   real(8),dimension (-1:nx+1,-1:ny+1,-1:nz+1) ::   USav,VSav,WSav,PSav
   integer i,j,k


501  format (5(1x,e23.15))

   USav=0
   VSav=0
   WSav=0 
   PSav=0
   !! x ,y ,z and phi are in the range that we want to store the data!! 

   do k=0,nz 
     do j=0,ny
       do i=0,nx        
         PSav(i,j,k)=p(i,j,k)
       end do
     end do
   end do

   do k=-1,nz+1 
     do j=-1,ny+1
       do i=0,nx+1        
         USav(i,j,k)=U(i,j,k)
       end do
     end do
   end do

   do k=-1,nz+1 
     do j=0,ny+1
       do i=-1,nx+1        
         VSav(i,j,k)=V(i,j,k)
       end do
     end do
   end do

   do k=0,nz+1 
     do j=-1,ny+1
       do i=-1,nx+1        
         WSav(i,j,k)=W(i,j,k)
       end do
     end do
   end do

   OPEN(unit=3005,file=fileplace//"Fluid_Dynamic_Save.dat",STATUS='REPLACE')   
 
   do k=-1,nz+1 
     do j=-1,ny+1 
       do i=-1,nx+1
         write(3005,501) USav(i,j,k),VSav(i,j,k),WSav(i,j,k),phi(i,j,k),PSav(i,j,k)
       end do 
     end do
   end do  
                   
   call flush (3005)
   close(3005)

   return
end subroutine Fluid_Dynamic_Write 



Subroutine Fluid_Dynamic_Read(p)
   use M_General, only: u,v,w,p,phi,nx,ny,nz,fileplace
   

   implicit none

   real(8),dimension (0:nx,0:ny,0:nz         ) ::   P
   real(8),dimension (-1:nx+1,-1:ny+1,-1:nz+1) ::   USav,VSav,WSav,PSav
   integer i,j,k

501  format (5(1x,e23.15))

   OPEN(unit=3005,file=fileplace//"Fluid_Dynamic_Save.dat")   
 
   do k=-1,nz+1 
     do j=-1,ny+1 
       do i=-1,nx+1
         Read(3005,501) USav(i,j,k),VSav(i,j,k),WSav(i,j,k),phi(i,j,k),PSav(i,j,k)
       end do 
     end do
   end do  
                   
   call flush (3005)  !! ?? 
   close(3005)

   do k=0,nz 
     do j=0,ny
       do i=0,nx        
         P(i,j,k)=PSav(i,j,k)
       end do
     end do
   end do

   do k=-1,nz+1 
     do j=-1,ny+1
       do i=0,nx+1        
         U(i,j,k)=USav(i,j,k)
       end do
     end do
   end do

   do k=-1,nz+1 
     do j=0,ny+1
       do i=-1,nx+1        
         V(i,j,k)=VSav(i,j,k)
       end do
     end do
   end do

   do k=0,nz+1 
     do j=-1,ny+1
       do i=-1,nx+1        
         W(i,j,k)=WSav(i,j,k)
       end do
     end do
   end do

return
end subroutine Fluid_Dynamic_Read

















     







