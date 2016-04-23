Module M_Tether

  use M_General,              only: pi,Ly,Bp,dt,RoDrop,gz,ZFree,fileplace              
  use M_Solid,                only: M_Cyl,tetang,R_Cyl_1,Solid_Bot_Ini,Xbar,YBar,ZBar,omegax,omegay,omegaz,ubar,vbar,wbar
  use M_Tower,                only: MNac,MTow,MRot
  use M_Math,                 only: LENGTHF,CROSS
  integer,parameter                               ::   tetn=3
  real(8),dimension(1:tetn,3)                     ::   point_Teth,point_Teth_old,point_Teth_old2
  real(8),dimension(1:tetn)                       ::   Ox,Oy,Oz,LZero,F_Teth_Amp    
  real(8),dimension(3)                            ::   FTether,MomTether
  real(8)                                         ::   KTeth,ETeth,Beta_Teth,DTethO,DTethI,TickTeth,leg   
    
    
    
  contains


    subroutine Tether_Constant_Ini()
      implicit none 

      real(8)              :: ResBuoy
      integer :: i 
      
      Include "Par_Constant_Tether.txt"
       
       
      Ox(1)=Bp(1,1)+(r_Cyl_1+leg)*cos(0.0)       ; Oy(1)=Bp(1,2)+(r_Cyl_1+leg)*sin(0.0)        
      Oz(1)=0.0

      Ox(2)=Bp(1,1)+(r_Cyl_1+leg)*cos(2*pi/3)    ; Oy(2)=Bp(1,2)+(r_Cyl_1+leg)*sin(2*pi/3)  
      Oz(2)=0.0

      Ox(3)=Bp(1,1)+(r_Cyl_1+leg)*cos(4*pi/3)    ; Oy(3)=Bp(1,2)+(r_Cyl_1+leg)*sin(4*pi/3)
      Oz(3)=0.0
       
      DTethO=DTethI+2*TickTeth

      KTeth=ETeth*pi/4*(DTethO**2-DTethI**2)/Solid_Bot_Ini
      ResBuoy=(ZFree-Solid_Bot_Ini)*(pi*r_Cyl_1**2)*RoDrop-(M_Cyl+MTow+MNac+MRot)

      do i=1,tetn
        Lzero(i)=Solid_Bot_Ini-ResBuoy*abs(gz)/(real(tetn)*KTeth)   !! inital position no angle 
      end do                       
      
      
      print*,'Reserver Buoyancy=',ResBuoy/1000,'Tonne'
      print*,"neutral lenght of tethers are:",Lzero(:)
      print*,"Tether spring coefficient=",KTeth/1000000,'MNewton'
      print*,"Tethers increased length in the statitcal stable position are=",ResBuoy*abs(gz)/(real(tetn)*KTeth)

      return 
    end subroutine 


    subroutine Tether_Dynamic_Ini()
      implicit none
      integer :: i 
      Include "Par_Dynamic_Tether.txt"

      do i=1,tetn
        print*,"Tether",i, "length in the first/current  time step configuration is:", &
        &  LENGTHF( ox(i)-point_Teth(i,1),oy(i)-point_Teth(i,2),oz(1)-point_Teth(i,3) ) 
      end do

 
      OPEN(65,file='Tether_Force.plt')
      write(65,111) 'variables="Time(sec)","Tether_1","Tether_2","Tether_3"'
      111 format (A110)

      print*,"Initial tethers values are ended"

      return
    end subroutine



    subroutine Tether_FM()
      implicit none     
      real(8),dimension (tetn)                                 :: Lnew,L_old2
      real(8),dimension(3)                                     :: F_Teth_Tmp,M_Teth_Tmp,r_Tmp
      integer i
     
    
      FTether(:)=0
      MomTether(:)=0
      do i=1,tetn   !! number of tethers !!
        
        L_old2(i)    =LENGTHF( ox(i)-point_Teth_old2(i,1),oy(i)-point_Teth_old2(i,2),oz(i)-point_Teth_old2(i,3) )
        Lnew(i)      =LENGTHF( ox(i)-point_Teth(i,1)     ,oy(i)-point_Teth(i,2)     ,oz(i)-point_Teth(i,3)      )
        F_Teth_Amp(i)=KTeth* ( Lnew(i)-Lzero(i) )+ Beta_Teth*KTeth*(Lnew(i)-L_old2(i))/dt
        !FSTether(1)=KTeth* ( max(Lzero(i),Lnew(i))-Lzero(i) )* ( Ox(i)-Px(i,4) )/Lnew(i)
        !FSTether(2)=KTeth* ( max(Lzero(i),Lnew(i))-Lzero(i) )* ( Oy(i)-Py(i,4) )/Lnew(i)
        !FSTether(3)=KTeth* ( max(Lzero(i),Lnew(i))-Lzero(i) )* ( Oz(i)-Pz(i,4) )/Lnew(i)
    
        F_Teth_Tmp(1)=F_Teth_Amp(i) * ( Ox(i)-point_Teth(i,1) )/Lnew(i)
        F_Teth_Tmp(2)=F_Teth_Amp(i) * ( Oy(i)-point_Teth(i,2) )/Lnew(i)
        F_Teth_Tmp(3)=F_Teth_Amp(i) * ( Oz(i)-point_Teth(i,3) )/Lnew(i)
      
        r_Tmp(1)=point_Teth(i,1)-XBar
        r_Tmp(2)=point_Teth(i,2)-YBar 
        r_Tmp(3)=point_Teth(i,3)-ZBar
        M_Teth_Tmp(:)=CROSS(r_Tmp,F_Teth_Tmp)
        FTether(:)  =  FTether(:)+  F_Teth_Tmp(:)
        MomTether(:)=MomTether(:)+  M_Teth_Tmp(:)
      end do 

      MomTether(1)=0
      MomTether(3)=0

      PRINT*,'tether moment',MomTether(2)

      return
    end subroutine 



    subroutine Tether_Sec_Sav()
      implicit none 

      point_Teth_old(:,:)=point_Teth(:,:)

      return           
    end subroutine
 
    subroutine Tether_It_Sav()
      implicit none   
 
      point_Teth_old2(:,:)=point_Teth(:,:)

      return           
    end subroutine



    subroutine Tether_Sec_Cor()
      implicit none 

      point_Teth(:,:)=0.5*(point_Teth(:,:)+point_Teth_old(:,:))

      return    
    end subroutine 



    subroutine Tether_Dynamic_Update()
      implicit none 

      call Tether_Dynamic_Update_Posi()
         
      return 
    end subroutine 


    subroutine Tether_Dynamic_Update_posi()
      implicit none 
 
      real(8),dimension(tetn,3)                       ::   point_Teth_n
      integer                                         ::   i

      do i=1,tetn
        point_Teth_n(i,1)=point_Teth_old2(i,1) +dt*( ubar + omegay*(point_Teth(i,3)-zbar)  -omegaz*(point_Teth(i,2)-ybar) )
        point_Teth_n(i,2)=point_Teth_old2(i,2) +dt*( vbar - omegax*(point_Teth(i,3)-zbar)  +omegaz*(point_Teth(i,1)-xbar) )
        point_Teth_n(i,3)=point_Teth_old2(i,3) +dt*( wbar + omegax*(point_Teth(i,2)-ybar)  -omegay*(point_Teth(i,1)-xbar) )
      end do 

      point_Teth(:,:)=point_Teth_n(:,:)
      
      return      
    end subroutine 

    
    subroutine Tether_Plot(tp)
      implicit none

      integer,intent(in) :: tp

      write(65,2012)  (tp+1)*dt,F_Teth_Amp(1),F_Teth_Amp(2),F_Teth_Amp(3)

      2012  format (4(1x,e15.7))
      call flush (65)
      
      return
    end subroutine  


    Subroutine Tether_Dynamic_Write()   
      implicit none 
      integer i
  
511   format (3(1x,e23.15))

      OPEN(unit=3045,file=fileplace//"Teth_Save.dat",STATUS='REPLACE')   
     
      do i=1,tetn
        write (3045,511) point_Teth(i,1),point_Teth(i,2),point_Teth(i,3)
      end do     
   

      close(3045)

      return
    end subroutine 



    Subroutine Tether_Dynamic_Read()   
      implicit none 
      integer i
  
      511   format (3(1x,e23.15))

      OPEN(unit=3045,file=fileplace//"Teth_Save.dat")   
     
      do i=1,tetn
        read (3045,511) point_Teth(i,1),point_Teth(i,2),point_Teth(i,3)
      end do     
   
      close(3045)



      OPEN(65,file='Tether_Force.plt',POSITION='APPEND',status='old')

      return
    end subroutine 


end module 






