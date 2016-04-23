 Module M_Tower

   use M_General,             only: pi,gz,dt,Froude,fileplace
   use M_Solid,               only: Point_Cyl_1,asolidx,AnAcY,xbar,ybar,zbar,omegay
   use M_math,                only: CROSS,LENGTHF

   real(8),dimension(3)   :: F_Tow,M_Tow,M_Tow_cg
   real(8),save           :: MRot,MNac,MTow,DiB,DiT,Di,TickB,TickT,Tick,TowEH,ETow,TurbH
   real(8),save           :: OffLen,DYawRotor,DYawNacelle,RSpinV,Le,ICTow,ITow,IRotY,IRotX,RRot,PlCgH
   real(8)                :: Gama_sav,Gama(3),CorTowSp,Beta_Tow,SolidTeta,K_Tow

   contains

     subroutine Tower_Constant_Ini()
       implicit none 

       include "Par_Constant_Tower.txt"


       Di=0.5*(DiB+DiT)
       Tick=0.5*(TickT+TickB)
       ICTow=pi/64.0*((Di+2*Tick)**4-Di**4)
       ITow=1.0/12.0*MTow*(  3*( (Di/2)**2 +((Di+2*Tick)/2)**2 )+TowEH**2  )
       IRotY=1.0/4.0*MRot*RRot**2

       return
      end subroutine 

      subroutine Tower_Dynamic_Ini()
        implicit none 

        PlCgH=LENGTHF(point_Cyl_1(1,1)-xbar,point_Cyl_1(1,2)-ybar,point_Cyl_1(1,3)-zbar)  !Towerb+Len-ZBar
        !! This Gama intial position is true if we do not wish independent inital tilting for the tower 
        Gama(:)=ATAN( (point_Cyl_1(1,1)-point_Cyl_1(2,1))/(point_Cyl_1(1,3)-point_Cyl_1(2,3)) )
        print*, "Zbar=",ZBar,"PlCgH=",PlCgH



        OPEN(2025,file='TowerMotion.plt')
        write(2025,2015) 'variables="t","Gama","Teta","aSolidX","AnAcy","mom","mom Cg","Forcex","Forcez","Mom damper"'  
        2015  format (A150)

        return
      end subroutine Tower_Dynamic_Ini

      subroutine Tower_Sec_Sav()

        implicit none  
        Gama_sav= Gama(1)
        
        return
      end subroutine

      subroutine Tower_It_Sav()

        implicit none  
        !! order is importent here!!
        Gama(3)=Gama(2)
        Gama(2)=Gama(1)
        
        return
      end subroutine

      subroutine Tower_Dynamic_Update()
        implicit none    

        K_Tow=CorTowSp*3*ETow*ICTow/(TurbH)  
        SolidTeta     =ATAN( (point_Cyl_1(1,1)     -point_Cyl_1(2,1))     /(point_Cyl_1(1,3)      -point_Cyl_1(2,3))     )

        !! difference of theta and theta dot are for two different time steps, however it should not highly change the results.
        Gama(1)=2*Gama(2)-Gama(3)                                                               &  
                 & +( dt*dt/(IRotY+ITow+(MRot+MNac)*TurbH**2) )                               &
                 &*(                                                                          &
                 & -  MTow             *TowEH               *aSolidX                          &
                 & - (MRot+MNac)       *TurbH               *aSolidX                          & 
                 & - (MTow)     * PlCgH*TowEH               *AnAcY                            &
                 & - (MRot+MNac)* PlCgH*TurbH               *AnAcY                            &
                                                                       
                 & +   MTow            *TowEH               *abs(gz)*Gama(2)                  &
                 & +  (MRot+MNac)      *TurbH               *abs(gz)*Gama(2)                  &
                  
                 & -           K_Tow*(  Gama(2)            -SolidTeta)                        &
                 & -  Beta_Tow*K_Tow*( (Gama(2)-Gama(3))/dt-omegay   )                        &
                 & )

        return
      end subroutine



      subroutine Tower_FM()
        implicit none
      
        real(8)              ,dimension(3)   :: R_Tmp


        F_Tow(1)= -(MNac+MRot)      *(TurbH*(Gama(1)-2*Gama(2)+Gama(3))/(dt**2)) &
            &     -MTow             *(TowEH*(Gama(1)-2*Gama(2)+Gama(3))/(dt**2)) &
            &     -(MNac+MRot+MTow) *(ASolidX+PlCgH*AnAcY) 
        F_Tow(2)=0
        F_Tow(3)=-(MTow+MRot+MNac)*abs(gz)   !! assuming the centifugal acceleration which is a nonlinear term is negligible
    
        M_Tow(1)=0
        M_Tow(2)=K_Tow*(Gama(1)-SolidTeta)+Beta_Tow*K_Tow*( (Gama(1)-Gama(2))/dt-omegay ) 
        M_Tow(3)=0

        R_Tmp(1)=Point_Cyl_1(1,1)-xbar 
        R_Tmp(2)=Point_Cyl_1(1,2)-ybar
        R_Tmp(3)=Point_Cyl_1(1,3)-zbar
   
        M_Tow_cg(:)=M_Tow(:)+CROSS(R_Tmp,F_Tow)

      return
     end subroutine


     subroutine Tower_Sec_Cor
       implicit none
       
       Gama(1)=0.5*(Gama(1)+Gama_Sav)
       Gama(2)=Gama_Sav
     end subroutine 
    

     subroutine Tower_Plot(tp)   
       implicit none

       integer,intent(in) :: tp
       !! it can be one time step shifting in force and angle due to prediction correction step.

       write(2025,2014)  tp*dt,Gama(1)*180/pi,SolidTeta*180/pi,asolidx,AnAcy*180/pi, &
                       & M_Tow(2),M_Tow_cg(2),F_Tow(1),F_Tow(3),Beta_Tow*K_Tow*( (Gama(1)-Gama(2))/dt-omegay )
                     
       2014  format (10(1x,e15.7))

     end subroutine   


     Subroutine Tower_Dynamic_Write()

   
       implicit none
   
       OPEN(unit=3025,file=fileplace//"Tower_Dynamic_Save.dat",STATUS='REPLACE') 

       522  format (1x,e23.15)

       write(3025,522) gama(1)
       write(3025,522) gama(2)
       write(3025,522) gama(3)
       write(3025,522) PlCgH

       call flush (3025)
       close(3025)

       return
     end subroutine Tower_Dynamic_Write


     Subroutine Tower_Dynamic_Read()

       implicit none

       OPEN(unit=3025,file=fileplace//"Tower_Dynamic_Save.dat")
       522  format (1x,e23.15)

       read(3025,522) gama(1)
       read(3025,522) gama(2)
       read(3025,522) gama(3)
       read(3025,522) PlCgH

       close(3025)

       OPEN(2025,file='TowerMotion.plt',POSITION='APPEND',STATUS='OLD')

       return
     end subroutine Tower_Dynamic_Read



  end module  


