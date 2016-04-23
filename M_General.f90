module M_General

  implicit none 
  save
  integer     , parameter                           :: nx=217, ny=101, nz=151  !nx=139, ny=67, nz=101  !nx=211, ny=60, nz=71 !nx=191,ny=60  ,nz=91 !  ny=60  nz=91
  real(8)     , parameter                           :: Froude=(1.0/12.0)   
  real(8)     , parameter                           :: pi=3.141592653
  real(8)     , parameter                           :: Landa=Froude*75.4
  CHARACTER(*), parameter                           :: fileplace = "/home/ali/Research/2015/April/5th/Backup/"
  logical     ,parameter                            :: Test=.False.

  real(8),dimension (0:nx+1,-1:ny+1,-1:nz+1)   ::   u
  real(8),dimension (-1:nx+1,0:ny+1,-1:nz+1)   ::   v
  real(8),dimension (-1:nx+1,-1:ny+1,0:nz+1)   ::   w
  real(8),dimension (-1:nx+1,-1:ny+1,-1:nz+1)  ::   phi 
  real(8),dimension (1:1,1:3)                  ::   Bp
  real(8),dimension (-1:nx+1)                  ::   x
  real(8),dimension (-1:ny+1)                  ::   y
  real(8),dimension (-1:nz+1)                  ::   z
  real(8),dimension (1:3)                      ::   Chv
  real(8),dimension (3)                        ::   Force,Mom,Force_Sav,Mom_sav

  real(8)                                      ::   lx,ly,lz,dt,gx,gy,gz,zfree
  real(8)                                      ::   beta,maxSOR,pdif,toler,maxError
  real(8)                                      ::   miudrop,miuair,rodrop,roair

  integer                                      ::   tpstar,countPlot,plot,tstepS,tstepE,SavMin
  integer                                      ::   isolver,Iterate,mxmatvec,maxIteration
  logical                                      ::   print_resid,nonzero_x


  contains
    Subroutine General_Dynamic_Ini()
      implicit none

      Include "Par_Dynamic_General.txt"
      return
    end subroutine 
   

    Subroutine General_Constant_Ini()
      implicit none
   
      Include "Par_Constant_General.txt"    
      return 
    end subroutine 


    subroutine General_Dynamic_Write(tp)   
      implicit none
      integer,intent(in) :: tp
 
      OPEN(unit=3035,file=fileplace//"General_Dynamic_Save.dat",STATUS='REPLACE') 

      Write(3035,*) tp+1
      Write(3035,*) CountPlot
      Write(3035,511) Force(1),Force(2),Force(3)
      Write(3035,511) Mom(1),Mom(2),Mom(3) 

      511   format (3(1x,e23.15))
 
      call flush (3035)  
 
      close(3035)

      return
    end subroutine General_Dynamic_Write

    subroutine General_Dynamic_read()
      implicit none


      OPEN(unit=3035,file=fileplace//"General_Dynamic_Save.dat") 

      read(3035,*) tsteps
      read(3035,*) CountPlot
      read(3035,511) Force(1),Force(2),Force(3)
      read(3035,511) Mom(1),Mom(2),Mom(3) 

      511   format (3(1x,e23.15))
     
      close(3035)

      return
    end subroutine General_Dynamic_read


end module
