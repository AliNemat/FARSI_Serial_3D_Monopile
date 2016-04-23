Module M_SolidFinder_Num
  use M_Math,                 only:  HFF
  use M_General,              only:  x,y,z,nx,ny,nz
  use M_Mesh,                 only:  hx,hy,hz
  implicit none
  parameter smooth=1
  parameter LinN=4
  integer         ,dimension(1:nx-1,1:nz-1)        :: ccL       
  double precision,dimension(1:nx-1,1:nz-1)        :: cc,CI,CF,CFS
  double precision,dimension(LinN)                 :: mm,mm1,mm2,alpha


  contains




    subroutine Solid_Box_Mark(pointO,Width,Ib_SolidFindB,npointsO)
      implicit none 
      double precision,intent(in),dimension(npointsO,3)        :: pointO
      real(8)         ,intent(out),dimension(0:nx,0:ny,0:nz)   :: Ib_SolidFindB
      integer         ,intent(in)                              :: npointsO
      real(8)         ,intent(in)                              :: Width

      double precision,dimension(npointsO+1,2)                 :: point
      integer                                                  :: i
                                     
      do i=1,npointsO
        point(i,1)=pointO(i,1)
        point(i,2)=pointO(i,3)
      end do 

      point(npointSO+1,1)=pointO(1,1)
      point(npointSO+1,2)=pointO(1,3)

      call BorderFinder(npointsO,point)
      call insideFinder(npointsO,point)
      call FractionFinder(npointsO,point)
      call WidthModel(npointsO,pointO,Ib_SolidFindB,width)
      ! call FractionSmoother()

      !CFSOut(:,:)=CFS(:,:)

      return
    end subroutine



     

    subroutine BorderFinder(npointsO,point)
      implicit none

      double precision,intent(in),dimension(npointsO+1,2)     :: point
      integer         ,intent(in)                             :: npointsO
      double precision,dimension(2)                           :: PointMax,PointMin
      integer kk,i,j,k
       real  sumc

  
      do kk=1,npointSO

        if (point(kk+1,1).ne.point(kk,1)) then 

          mm(kk)=(point(kk+1,2)-point(kk,2))/(point(kk+1,1)-point(kk,1))

          mm1(kk)=-mm(kk)/sqrt(1+mm(kk)*mm(kk))
          mm2(kk)=  1.0/sqrt(1+mm(kk)*mm(kk))
       
          alpha(kk)=(point(kk,2)-mm(kk)*point(kk,1))/sqrt(1+mm(kk)*mm(kk))

        else

          mm1(kk)=1
          mm2(kk)=0
          alpha(kk)=point(kk,1)

        end if

      end do 
      !print*, "mm1",mm1(:)
      !print*, "mm2",mm2(:)
      !print*, "alpha",alpha(:)

      PointMax=maxval (point,dim=1)
      PointMin=minval (point,dim=1)
!      print*, pointMin(1),pointMin(2)
!      print*,Pointmax(1),pointmax(2)
 !     print*,"point values"
 !     do i=1,5
 !      do j=1,2
 !        print*,point(i,j),i,j
 !      end do 
 !     end do 
 !     print*, "These are point values"
       
 !      pointMin(1)=400
 !      pointMin(2)=100

 !      pointMax(1)=600
 !       pointMax(2)=130
      

      ccl(:,:)=0
      do i=1, nx-1
        if      ( x(i).gt. ( pointMax(1)+0.2*(pointMax(1)-pointMin(1)))  ) then 
          cc(i,:)=0
        else if ( x(i).lt. ( pointMin(1)-0.2*(pointMax(1)-pointMin(1)))  ) then
          cc(i,:)=0
        else
          do k=1,nz-1
            if      ( Z(k).gt. ( pointMax(2)+0.2*(pointMax(2)-pointMin(2)))  ) then 
              cc(i,k)=0
            else if ( Z(k).lt. ( pointMin(2)-0.2*(pointMax(2)-pointMin(2)))  ) then
              cc(i,k)=0
            else 
              do kk=1,npointsO

                if (mm2(kk).eq.0) then
                  !   print*, " I am here1 " 
                  if ( (x(i)-0.5*hx(i)).lt.alpha(kk).AND.(x(i)+0.5*hx(i)).gt.Alpha(kk) ) then
                    if ( ( (x(i)-point(kk,1  ))**2 + (z(k)-point(kk,2  ))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                         ( (x(i)-point(kk+1,1))**2 + (z(k)-point(kk+1,2))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ) ) then

                           !print*, " I am here2 ",i,k,x(i),z(k) 
                      cc(i,k)=1.0
                      ccl(i,k)=kk
                      exit
                    end if  
                  end if 
                else
                  if (      max(  (Alpha(kk)-mm1(kk)*(x(i)-0.5*hx(i)))/mm2(kk),(Alpha(kk)-mm1(kk)*(x(i)+0.5*hx(i)))/mm2(kk) ).lt. (z(k)-0.5*hz(k)) ) then
                    cc(i,k)=0
                  else if ( min ( (Alpha(kk)-mm1(kk)*(x(i)-0.5*hx(i)))/mm2(kk),(Alpha(kk)-mm1(kk)*(x(i)+0.5*hx(i)))/mm2(kk) ).gt. (z(k)+0.5*hz(k)) ) then 
                    cc(i,k)=0
                  else
                 
                    if ( ( (x(i)-point(kk,1  ))**2 + (z(k)-point(kk,2  ))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                         ( (x(i)-point(kk+1,1))**2 + (z(k)-point(kk+1,2))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ) ) then

                           !print*, " I am here3 ",i,k,x(i),z(k) 
                      cc(i,k)=1.0
                      ccl(i,k)=kk
                      exit
                    end if  
                  end if 
                end if
              end do
            end if 
          end do
        end if
      end do

!sumc=0
!do i=1,nx-1 ; do k=1,nz-1

!sumc=cc(i,k)+sumc
!end do 
!end do 
!print*, sumc, "sumc"


      OPEN(2350,file='BorderFinder.plt')

      write(2350,*) 'zone i=',nx-1,' k=',nz-1
      do k=1,nz-1 ;  Do i=1,nx-1

      write(2350,3500) x(i),z(k),cc(i,k)
      end do ; end do
      3500  format (3(1x,e15.7))
      call flush (2350)


      return

    end subroutine


subroutine InsideFinder(npointsO,point)
implicit none

double precision,intent(in),dimension(npointsO+1,2)     :: point
integer         ,intent(in)                            :: npointsO
double precision                                       :: ZCross
double precision,dimension(2)                          :: PointMax,PointMin
integer i,k,kk




    PointMax=maxval (point,dim=1)
    PointMin=minval (point,dim=1)

    CI(:,:)=0
  do i=1,nx-1 
    if      ( x(i).gt. ( pointMax(1)+0.2*(pointMax(1)-pointMin(1)))  ) then 

      CI(i,:)=0

    else if ( x(i).lt. ( pointMin(1)-0.2*(pointMax(1)-pointMin(1)))  ) then

      CI(i,:)=0
    else

      do k=1,nz-1
        if      ( z(k).gt. ( pointMax(2)+0.2*(pointMax(2)-pointMin(2)))  ) then 
          CI(i,k)=0
        else if ( z(k).lt. ( pointMin(2)-0.2*(pointMax(2)-pointMin(2)))  ) then
          CI(i,k)=0
        else 
          do kk=1,npointsO
            !print*, " I come inside"            
            if (mm2(kk).ne.0) then
              ZCross=(Alpha(kk)-mm1(kk)*x(i))/mm2(kk)

              if (((ZCross-point(kk  ,2))**2 + (X(i)-point(kk  ,1))**2).le.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
              &   ((ZCross-point(kk+1,2))**2 + (X(i)-point(kk+1,1))**2).le.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
              &   ((Alpha(kk)-mm1(kk)*x(i))/mm2(kk)).gt.z(k)                                                                                      ) then
                    
                cI(i,k)=cI(i,k)+1
             !   print*, "I am inside"
              end if  
            end if
          
          end do
       
          if (cI(i,k).eq.1.AND.cc(i,k).eq.0) then
            CI(i,k)=1.0
          else 
            CI(i,k)=0
          end if
 
        end if 
      end do

    end if
  end do  
         
OPEN(2351,file='InsideFinder.plt')

write(2351,*) 'zone i=',nx-1,' k=',nz-1
do k=1,nz-1 ;  Do i=1,nx-1

write(2351,3501) x(i),z(k),cI(i,k),cc(i,k)
end do ; end do
3501  format (4(1x,e15.7))
call flush (2351)







 

return
end subroutine 


subroutine FractionFinder(npointsO,point)
implicit none
double precision,intent(in),dimension(npointsO+1,2)      :: point
double precision,dimension(npointsO)                    :: mm1N,mm2N,alphaN
integer         ,intent(in)                             :: npointsO
integer kk,i,k


  CF(:,:)=0 
  do i=1,nx-1
    do k=1,nz-1
      if (cc(i,k).ne.0) then
        kk=ccl(i,k) 
        if       ( (point(kk+1,1)-point(kk,1)).gt.0.AND.( point(kk+1,2)-point(kk,2) ).gt.0 ) then
          mm1N(kk)= mm1(kk)
          mm2N(kk)=-mm2(kk)
          AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)-0.5*hx(i)) -mm2(kk)*(z(k)+0.5*hz(k))
          !print*, "hx,hz",hx(i),hz(k)
          !print*,"old value",mm1(1),mm2(1), Alpha(1)
          !print*,"new value",mm1N(1),mm2N(1), AlphaN(1)
          !print*, "x(i),z(k)", x(i),z(k) 
 
          !print*,  "I am here 1"
          if (AlphaN(kk).le.0.0) then
            mm1N(kk)=-mm1N(kk)
            mm2N(kk)=-mm2N(kk)
            AlphaN(kk)=-AlphaN(kk)
          end if 
        else if  ( (point(kk+1,1)-point(kk,1)).le.0.AND.( point(kk+1,2)-point(kk,2) ).ge.0 ) then
          mm1N(kk)= mm1(kk)
          mm2N(kk)= mm2(kk)
          AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)-0.5*hx(i)) -mm2(kk)*(z(k)-0.5*hz(k))
        ! print *, " I am here 2"
          if (AlphaN(kk).le.0.0) then
            mm1N(kk)=-mm1N(kk)
            mm2N(kk)=-mm2N(kk)
            AlphaN(kk)=-AlphaN(kk)
          end if 
        else if ( (point(kk+1,1)-point(kk,1)).le.0.AND.( point(kk+1,2)-point(kk,2) ).le.0 ) then
        ! print*, "I am here 3 "
          mm1N(kk)=-mm1(kk)
          mm2N(kk)= mm2(kk)
          AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)+0.5*hx(i)) -mm2(kk)*(z(k)-0.5*hz(k))
          if (AlphaN(kk).le.0.0) then
            mm1N(kk)=-mm1N(kk)
            mm2N(kk)=-mm2N(kk)
            AlphaN(kk)=-AlphaN(kk)
          end if 
        else
          !print*,  "I am here 4"
          mm1N(kk)=-mm1(kk)
          mm2N(kk)=-mm2(kk)
          AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)+0.5*hx(i)) -mm2(kk)*(z(k)+0.5*hz(k))
          if (AlphaN(kk).le.0.0) then
            mm1N(kk)=-mm1N(kk)
            mm2N(kk)=-mm2N(kk)
            AlphaN(kk)=-AlphaN(kk)
          end if 
        end if
        
        if     (mm2N(kk).eq.0) then
            CF(i,k)=(AlphaN(kk)/mm1N(kk))/hx(i)
         !print*," I am here 5 "
        else if(mm1N(kk).eq.0) then
            CF(i,k)=(AlphaN(kk)/mm2N(kk))/hz(k)
        !print *,  "I am here 6" 
        else
          ! print*, "I am here 7"
            CF(i,k)=(AlphaN(kk)**2)/(2*mm1N(kk)*mm2N(kk))/(hx(i)*hz(k))*(1 &
                                                           & -HFF( AlphaN(kk)-mm1N(kk)*hx(i) )*( (AlphaN(kk)-mm1N(kk)*hx(i))/AlphaN(kk) )**2 &
                                                           & -HFF( AlphaN(kk)-mm2N(kk)*hz(k) )*( (AlphaN(kk)-mm2N(kk)*hz(k))/AlphaN(kk) )**2  )
                                                        
          !print*,"CF",CF(i,k) 
        end if
      
      end if
      if (cI(i,k).eq.1.0) then
         CF(i,k)=1.0
      end if 


    end do
  end do 
  


  !if (Test) then   

OPEN(2352,file='FractionFinder.plt') 

write(2352,*) 'zone i=',nx-1,' k=',nz-1 

 do k=1,nz-1 ;  Do i=1,nx-1

    write(2352,3502) x(i),z(k),cI(i,k),cc(i,k),CF(i,k),dble (ccl(i,k))
  end do ; end do
  3502  format (6(1x,e15.7))
  call flush (2352)
!end if 

return

end subroutine 


    subroutine WidthModel(npointsO,pointO,Ib_SolidFindB,width)
      implicit none

      double precision,intent(in),dimension(npointsO,3)             :: pointO
      integer         ,intent(in)                                   :: npointsO
      real(8),intent(out),dimension(0:nx,0:ny,0:nz)                 :: Ib_SolidFindB
      real(8),intent(in)                                            :: Width
 
      real(8),            dimension(1:nx-1,1:ny-1,1:nz-1)           :: Ib_SolidFind 
      integer                                                       :: i,j,k,JCr1,JCr2

      do i=1,nx-1
        do k=1,nz-1 
          
          do j=1,ny-1 
            Ib_SolidFind(i,j,k)=CF(i,k)
          end do 
        
        end do 
      end do 
    
      do j=1,ny-1
        if ( y(j).lt.pointO(1,2) ) then 
        else
          JCr1=j
          exit
        end if 
      end do

      do j=jcr1+1,ny-1
        if ( y(j).lt.(pointO(1,2)+width)) then 
        else
          JCr2=j
          exit
        end if 
      end do

       
      print*,'critical values are',Jcr1,Jcr2
 
      !!! constructing right sideO      if ( (pointO(1,2)-y(jcr1)).gt.(0.5*hy(jcr1)) ) then
      if ( (y(jcr1)-pointO(1,2)).gt.(0.5*hy(jcr1)) ) then
        Ib_SolidFind(:,1:jcr1-2,:)=0
        IB_SolidFind(:,jcr1-1,:)=CF(:,:)*(y(jcr1-1)+0.5*hy(jcr1-1)-pointO(1,2))/hy(jcr1-1)
!        print*,'I am here 3d 1'
      else

        Ib_SolidFind(:,1:jcr1-1,:)=0
        IB_SolidFind(:,jcr1,:)=CF(:,:) *(y(jcr1)+0.5*hy(jcr1)-pointO(1,2))/hy(jcr1)

 !       print*,'I am here 3d 2'
      endif 

      !! constructing left side
      if ( (y(jcr2)-(pointO(1,2)+width)).gt. (0.5*hy(jcr2)) ) then

        Ib_SolidFind(:,jcr2:ny-1,:)=0
        IB_SolidFind(:,jcr2-1,:)=CF(:,:) *(pointO(1,2)+width-(y(jcr2-1)-0.5*hy(jcr2-1)) )/hy(jcr2-1)

 !       print*,'I am here 3d 3'
      else

        Ib_SolidFind(:,jcr2+1:ny-1,:)=0
        IB_SolidFind(:,jcr2,:)=CF(:,:) *(pointO(1,2)+width-(y(jcr2)-0.5*hy(jcr2)) )/hy(jcr2)

  !      print*,'I am here 3d 4'
      endif 
      !! the previous values for the rest is okey

      Ib_SolidFindB(:,:,:)=0.0
      do i=1,nx-1
        do j=1,ny-1 
          do k=1,nz-1 
            IB_SolidFindB(i,j,k)=Ib_SolidFind(i,j,k)
          end do 
        end do 
      end do 

      return
    end subroutine 
     




Subroutine FractionSmoother()
implicit none
integer i,k,kk

 CFS(:,:)=CF(:,:)  !! for boundary points !!


  if (smooth.ne.0) then
 
    do kk=1,smooth
      do i=2,nx-2 ; do k=2,nz-2 
        CFS(i,k)= 0.25*(CF(i,k))+ &
        &   0.125*(CF(i+1,k)+CF(i-1,k)+CF(i,k-1)+CF(i,k+1)) + &
        &  0.0625*(CF(i+1,k+1)+CF(i+1,k-1)+CF(i-1,k+1)+CF(i-1,k-1))
      end do 
      end do 
    end do 
  end if 
  
return 
end subroutine  




 
          
      

end  module 


  

 
         
          
     
        
        
      
          


  

  
  
 


  
 


