module M_Mesh
  use M_General,           only:x,y,z,nx,ny,nz,lx,ly,lz,landa,ChV,Bp

  real          :: hx(0:nx-1),hy(0:ny-1),hz(0:nz-1) 


  contains 
    subroutine meshgenerator()
      implicit none

      integer :: meshgen
      real    :: XE(0:nx-1),YE(0:ny-1),ZE(0:nz-1)
      real(8) :: Lx1,Lx2,Lx3,Lx4,Lz1,Lz2,Lz3,Ly1,Ly2,Ly3
      real(8) :: hxx1,hxx2,hxx3,hxx4,hyy1,hyy2,hyy3,hzz1,hzz2,hzz3
      real(8) :: Amesh12,Amesh14,Amesh21,Amesh23,Amesh31,Amesh32,Amesh33,betam12,betam14,betam21,betam23
      real(8) :: betam31,betam32,betam33,Dmesh12,Dmesh14,Dmesh21,Dmesh23,Dmesh31,Dmesh32,Dmesh33
      integer :: i,j,k,nx1,nx2,nx3,ny1,ny2,nz1,nz2


include 'Par_Constant_Mesh.txt'
if (meshgen.eq.1) then 
    x(1)=0.5*Lx/(nx-1)
    x(-1)=-3*x(1) ; x(0)=-x(1)
    do i=1,nx+1
      x(i)=x(i-1)+Lx/(nx-1)
    end do 
    hx(:)=Lx/(nx-1)


    y(1)=0.5*Ly/(ny-1)
    y(-1)=-3*y(1) ; y(0)=-y(1)
    do j=1,ny+1
      y(j)=y(j-1)+Ly/(ny-1)
    end do 
    hy(:)=Ly/(ny-1)



    z(1)=0.5*Lz/(nz-1)
    z(-1)=-3*z(1) ; z(0)=-z(1)
    do k=1,nz+1
      z(k)=z(k-1)+Lz/(nz-1)
    end do
    hz(:)=Lz/(nz-1)   



else

    nx2=nx2+nx1
    nx3=nx3+nx2
    nz2=nz2+nz1    
    ny2=ny2+ny1
    
    hxx1=1.0/dble(nx1-1) ; hxx2=1.0/dble(nx2-nx1) ;  hxx3=1.0/dble(nx3-nx2) ; hxx4=1.0/dble(nx-nx3)
    hyy1=1.0/dble(ny1-1) ; hyy2=1.0/dble(ny2-ny1) ;  hyy3=1.0/dble(ny-ny2)  
    hzz1=1.0/dble(nz1-1) ; hzz2=1.0/dble(nz2-nz1) ;  hzz3=1.0/dble(nz-nz2)

    Amesh12=1/(2*betam12)*log(  (  1+( exp(betam12)-1 )*Dmesh12/Lx2 )/(  1+( exp(-betam12)-1 )*Dmesh12 /Lx2   )    ) 
    Amesh14=1/(2*betam14)*log(  (  1+( exp(betam14)-1 )*Dmesh14/Lx4 )/(  1+( exp(-betam14)-1 )*Dmesh14 /Lx4   )    ) 


    Amesh21=1/(2*betam21)*log(  (  1+( exp(betam21)-1 )*Dmesh21/Ly1 )/(  1+( exp(-betam21)-1 )*Dmesh21 /Ly1   )    ) 
    Amesh23=1/(2*betam23)*log(  (  1+( exp(betam23)-1 )*Dmesh23/Ly3 )/(  1+( exp(-betam23)-1 )*Dmesh23 /Ly3   )    )

    Amesh31=1/(2*betam31)*log(  (  1+( exp(betam31)-1 )*Dmesh31/Lz1 )/(  1+( exp(-betam31)-1 )*Dmesh31 /Lz1   )    ) 
    Amesh32=1/(2*betam32)*log(  (  1+( exp(betam32)-1 )*Dmesh32/Lz2 )/(  1+( exp(-betam32)-1 )*Dmesh32 /Lz2   )    )  
    Amesh33=1/(2*betam33)*log(  (  1+( exp(betam33)-1 )*Dmesh33/Lz3 )/(  1+( exp(-betam33)-1 )*Dmesh33 /Lz3   )    )  
     
     
     do i=0,nx1-1
       XE(i)=hxx1*Lx1*dble(i)
     end do 
     
     do i=nx1,nx2-1  !nx-1
       xE(i)=XE(nx1-1)+Dmesh12* (  1+ (  sinh ( betam12*(dble(i-(nx1-1))*hxx2-Amesh12) )  )/sinh(betam12*Amesh12)  )
     end do 
     
     do i=nx2,nx3-1
       xE(i)=xE(nx2-1)+hxx3*Lx3*dble( i-(nx2-1) )
     end do 

     do i=nx3,nx-1
       xE(i)=xE(nx3-1)+Dmesh14* (  1+ (  sinh ( betam14*(dble(i-(nx3-1))*hxx4-Amesh14) )  )/sinh(betam14*Amesh14)  )
     end do 
     
     
     YE(0)=0
     do j=1,ny1-1 
       yE(j)=          Dmesh21* (  1+ (  sinh ( betam21*(dble(j)        *hyy1-Amesh21) )  )/sinh(betam21*Amesh21)  )
     end do           
     
     do j=ny1,ny2-1
       yE(j)=yE(ny1-1)+hyy2*Ly2*dble(j-(ny1-1))
     end do 
         
     do j=ny2,ny-1
       yE(j)=yE(ny2-1)+Dmesh23* (  1+ (  sinh ( betam23*(dble(j-(ny2-1))*hyy3-Amesh23) )  )/sinh(betam23*Amesh23)  )
     end do 
     

     ZE(0)=0
     do k=1,nz1-1
       zE(k)=         Dmesh31* (  1+ (  sinh ( betam31*(dble(k)         *hzz1-Amesh31) )  )/sinh(betam31*Amesh31)  )
     end do 
     
     do k=nz1,nz2-1
       zE(k)=zE(nz1-1)+Dmesh32* (  1+ (  sinh ( betam32*(dble(k-(nz1-1))*hzz2-Amesh32) )  )/sinh(betam32*Amesh32)  )
     end do 
     
     do k=nz2,nz-1
       zE(k)=zE(nz2-1)+Dmesh33* (  1+ (  sinh ( betam33*(dble(k-(nz2-1))*hzz3-Amesh33) )  )/sinh(betam33*Amesh33)  )
     end do 
  

     do i=1,nx-1
       x(i)=XE(i-1)+0.5*(XE(i)-XE(i-1))
       hx(i)=XE(i)-XE(i-1)
     end do
     x(-1)=-3*x(1) ; x(0)=-x(1)
     x(nx)=x(nx-1)+ ( x(nx-1)-x(nx-2) ) ; x(nx+1)=x(nx)+ ( x(nx)-x(nx-1) )

     do j=1,ny-1
       y(j)=YE(j-1)+0.5*(YE(j)-YE(j-1))
       hy(j)=YE(j)-YE(j-1)
     end do
     y(-1)=-3*y(1) ;y(0)=-y(1)
     y(ny)=y(ny-1)+ ( y(ny-1)-y(ny-2) ) ; y(ny+1)=y(ny)+ ( y(ny)-y(ny-1) )

     do k=1,nz-1
       Z(k)=ZE(k-1)+0.5*(ZE(k)-ZE(k-1))
       hz(k)=ZE(k)-ZE(k-1)
     end do
     z(-1)=-3*z(1) ;z(0)=-z(1)
     z(nz)=z(nz-1)+ ( z(nz-1)-z(nz-2) ) ; z(nz+1)=z(nz)+ ( z(nz)-z(nz-1) )
  
     

 end if 


OPEN(205,file='meshx.plt')
OPEN(215,file='meshy.plt')
OPEN(225,file='meshz.plt')
OPEN(235,file='grideval2DXZ.plt')

do i=1,nx-1
write(205,*) i,x(i)
end do 
call flush (205)
do j=1,ny-1
write(215,*) j,y(j)
end do 
call flush (215)
do k=1,nz-1
write(225,*) k,z(k)
end do 
call flush (225)

write(235,*) 'zone i=',nx-1,' k=',nz-1
do k=1,nz-1 ;  Do i=1,nx-1

write(235,351) x(i),z(k)
end do ; end do
351 format (3(1x,e15.7))
call flush (235)


 return 
  end subroutine 
  
  


     
end module 






