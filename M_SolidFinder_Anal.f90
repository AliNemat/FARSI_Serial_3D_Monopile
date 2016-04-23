Module M_SolidFinder_Anal
  use M_Math,                 only:  LENGTHF
  use M_General,              only:  x,y,z,nx,ny,nz 
  implicit none


  contains


    subroutine Solid_Sphe_Mark(point_Sphe,r_Sphe,gap,Ib_Sphe,num_Sphe)

      implicit none
      integer num_Sphe
      real(8),dimension (num_Sphe,3)        :: point_Sphe
      real(8),dimension (0:nx,0:ny,0:nz)    :: Ib_Sphe
      real(8)                               :: r_Sphe,gap
      real(8) d3,x1,y1,z1
      integer i,j,k


      x1=point_Sphe(1,1)
      y1=point_Sphe(1,2)
      z1=point_Sphe(1,3)
      
        
      do i=0,nx 
        do j=0,ny 
          do k=0,nz

            d3=sqrt( (x(i)-x1)**2 + (y(j)-y1)**2 + (z(k)-z1)**2 )
            Ib_Sphe(i,j,k)=0.5+0.5*( (r_Sphe-d3)**3 + 1.50* gap*gap *(r_Sphe-d3) )  / ( (r_Sphe-d3)*(r_Sphe-d3)+gap*gap )**1.50 

            if (Ib_Sphe(i,j,k).le.0.001) then 
              Ib_Sphe(i,j,k)=0
            else if (Ib_Sphe(i,j,k).ge.0.999) then 
              Ib_Sphe(i,j,k)=1.0 
            end if


          end do 
        end do 
      end do 
 

      return 

    End subroutine 


    subroutine Solid_Cyl_Mark(point_Cyl,r_Cyl,gap,Ib_Cyl,num_Cyl)

      implicit none

      integer num_Cyl
      real(8),dimension (num_Cyl,3)        :: point_Cyl     
      real(8) x1,x2,y1,y2,z1,z2,xbar1,ybar1,zbar1
      real(8), dimension(0:nx,0:ny,0:nz)   :: Ib1,Ib2,Ib_Cyl
      real(8) ax,ay,az,bx,by,bz,r_Cyl,tt,cx,cy,cz,d1,d2,gap,len_Cyl
      integer i,j,k

      x1=point_Cyl(1,1)
      y1=point_Cyl(1,2)
      z1=point_Cyl(1,3)

      x2=point_Cyl(2,1)
      y2=point_Cyl(2,2)
      z2=point_Cyl(2,3)
      
      len_Cyl=LENGTHF(X1-X2,Y1-Y2,Z1-Z2)

      xbar1=0.5*(x1+x2) ; ybar1=0.5*(y1+y2) ; zbar1=0.5*(z1+z2)  

      do i=0,nx 
        do j=0,ny 
          do k=0,nz

            ax=x(i)-x1 ; ay=y(j)-y1 ; az=z(k)-z1
            bx=x(i)-x2 ; by=y(j)-y2 ; bz=z(k)-z2
 
            d1=sqrt( (ay*bz-az*by)**2 + (ax*bz-az*bx)**2 + (ax*by-ay*bx)**2 )/ sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
 
            tt= ax*( x2-x1 ) + ay*( y2-y1 ) +  az*( z2-z1 ) 
            tt=tt/( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
 
            cx=x1+ ( x2-x1 )*tt
            cy=y1+ ( y2-y1 )*tt
            cz=z1+ ( z2-z1 )*tt

            Ib1(i,j,k)=0.50+0.50*( (r_Cyl-d1)**3 + 1.50* gap*gap *(r_Cyl-d1) )  / ( (r_Cyl-d1)*(r_Cyl-d1)+gap*gap )**1.50  

            d2=sqrt( (cx-xbar1)**2+ (cy-ybar1)**2+(cz-zbar1)**2 )

            Ib2(i,j,k)=0.50+0.50*( (0.5*len_Cyl-d2)**3 + 1.50* gap*gap *(0.5*len_Cyl-d2) )  / ( (0.5*len_Cyl-d2)*(0.5*len_Cyl-d2)+gap*gap )**1.50  
            Ib_Cyl(i,j,k)=Ib1(i,j,k)*Ib2(i,j,k)


 
            if (Ib_Cyl(i,j,k).le.0.001) then 
              Ib_Cyl(i,j,k)=0
            else if (Ib_Cyl(i,j,k).ge.0.999) then 
              Ib_Cyl(i,j,k)=1.0 
            end if

          end do 
        end do 
      end do 
 

      return 

    End subroutine 


end module












