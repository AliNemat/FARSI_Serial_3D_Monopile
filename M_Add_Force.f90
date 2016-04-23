Module M_Add_Force
  use M_General,          only: x,y,z,nx,ny,nz,test,Force,Mom
  use M_Solid,            only: Ib_Solid,XBar,YBar,ZBar
  use M_Math,             only: MTRXSOL

  real(8),dimension (1:nx-1,1:ny-1,1:nz-1)              ::   fbox,fboy,fboz

  contains

    Subroutine AddForce(ro) 

      implicit none 


      real(8),intent(in),dimension (0:nx,0:ny,0:nz)         ::   ro

      !! local variables !!
      real(8)                                               ::   MassPlat,PElement
      real(8),dimension(3)                                  ::   SumM,VAccPlat,VAnAcPlat
      real(8),dimension(3)                                  ::   rTmp
      real(8),dimension(3,3)                                ::   IPlat

      integer i,j,k 


      MassPlat=0  
      IPlat=0   
      do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1
        Pelement= 0.125*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*(z(k+1)-z(k-1))
        MassPlat=  MassPlat+                                      PElement*( Ib_Solid(i,j,k)*ro(i,j,k) )
        IPlat(3,3)=IPlat(3,3)+ ( (x(i)-xbar)**2+ (y(j)-ybar)**2 )*pelement*( Ib_Solid(i,j,k)*ro(i,j,k) ) 
        IPlat(2,2)=IPlat(2,2)+ ( (x(i)-xbar)**2+ (z(k)-zbar)**2 )*pelement*( Ib_Solid(i,j,k)*ro(i,j,k) ) 
        IPlat(1,1)=IPlat(1,1)+ ( (y(j)-ybar)**2+ (z(k)-zbar)**2 )*pelement*( Ib_Solid(i,j,k)*ro(i,j,k) ) 
        IPlat(1,2)=IPlat(1,2)+ ( (y(j)-ybar)   * (x(i)-xbar)    )*pelement*( Ib_Solid(i,j,k)*ro(i,j,k) ) 
        IPlat(1,3)=IPlat(1,3)+ ( (x(i)-xbar)   * (z(k)-zbar)    )*pelement*( Ib_Solid(i,j,k)*ro(i,j,k) )        !! mistake solved!!
        IPlat(2,3)=IPlat(2,3)+ ( (y(j)-ybar)   * (z(k)-zbar)    )*pelement*( Ib_Solid(i,j,k)*ro(i,j,k) ) 
      end do ; end do ; end do
      IPlat(1,2)=-IPlat(1,2)
      IPlat(1,3)=-IPlat(1,3)
      IPlat(2,3)=-IPlat(2,3)
       
      IPlat(2,1)=IPlat(1,2)
      IPlat(3,1)=IPlat(1,3)
      IPlat(3,2)=IPlat(2,3)

  
      print*,'IPlat(2,2)',IPlat(2,2)
      print*,'IPlat(3,3)',IPlat(3,3)
      print*,'Mom',mom(2)
      VAccPlat(:)=Force(:)/MassPlat           
      VAnAcPlat(:)=MTRXSOL(Mom,IPlat)
     
      print*,'Matrix Clac',VAnAcplat
   
      FBox=0 ; Fboy=0 ; FBoz=0 ;
      do i=1,nx-1 ; do j=1,ny-1 ; do k=1,nz-1 
        rTmp(1)=x(i)-XBar
        rTmp(2)=y(j)-YBar
        rTmp(3)=Z(k)-Zbar
        FBoX(i,j,k)= 0.5* (Ib_Solid(i-1,j,k)+Ib_Solid(i,j,k)) *(VAccPlat(1)+VAnAcPlat(2)*rTmp(3)-VAnAcPlat(3)*rTmp(2))
        FBoY(i,j,k)= 0.5* (Ib_Solid(i,j-1,k)+Ib_Solid(i,j,k)) *(VAccPlat(2)-VAnAcPlat(1)*rTmp(3)+VAnAcPlat(3)*rTmp(1))
        FBoZ(i,j,k)= 0.5* (Ib_Solid(i,j,k-1)+Ib_Solid(i,j,k)) *(VAccPlat(3)+VAnAcPlat(1)*rTmp(2)-VAnAcPlat(2)*rTmp(1))     
      end do ; end do ; end do
      
      return   
    end subroutine 

 
end module 




