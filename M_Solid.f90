

Module M_Solid

  use M_General,     only: dt,nx,ny,nz,pi,ChV,Bp,x,y,z,Froude,fileplace
  use M_SolidFinder_Anal
  use M_SolidFinder_Num



  implicit none
  !! 2 cylinder
  real(8),dimension(2,3)                      ::   point_Cyl_1,point_Cyl_2
  real(8),dimension(2,3)                      ::   point_Cyl_old_1,point_Cyl_old_2
  real(8),dimension(2,3)                      ::   point_Cyl_old2_1,point_Cyl_old2_2

  real(8)                                     ::   Len_Cyl_1,r_Cyl_1,Len_Cyl_2,r_Cyl_2,Per_Con
  real(8),dimension (0:nx,0:ny,0:nz)          ::   Ib_Cyl_1,Ic_Cyl_2

  !! 2 sphere 
  real(8),dimension(1,3)                      ::   point_Sphe_1,point_Sphe_2
  real(8),dimension(1,3)                      ::   point_Sphe_old_1,point_Sphe_old_2
  real(8),dimension(1,3)                      ::   point_Sphe_old2_1,point_Sphe_old2_2

  real(8)                                     ::   r_sphe_1,r_sphe_2
  real(8),dimension (0:nx,0:ny,0:nz)          ::   Ib_Sphe_1,Ib_Sphe_2

  !! 2 Box 
  real(8),dimension(4,3)                      ::   point_Box_1,point_box_2
  real(8),dimension(4,3)                      ::   point_Box_old_1,point_box_old_2
  real(8),dimension(4,3)                      ::   point_Box_old2_1,point_box_old2_2

  real(8)                                     ::   W_Box_1,W_Box_2,Len_Box_1,Len_Box_2,h_Box_1,h_Box_2
  real(8),dimension (0:nx,0:ny,0:nz)          ::   Ib_Box_1,Ib_Box_2
 


  real(8),dimension (0:nx,0:ny,0:nz)          ::   Ib_Solid,Ic_Solid,Ib_Solid_old,Ic_Solid_old
  Real(8)                                     ::   rocon,rosolid,miusolid,miucon,tetang,M_Cyl,Solid_Bot_Ini
  

  real(8)                                     ::   Tmass,xbar,ybar,zbar,ubar,vbar,wbar,omegax,omegay,omegaz
  real(8)                                     ::   asolidx,asolidy,asolidz,anacx,anacy,anacz
  real(8)                                     ::   SumU,SumV,SumW,SumIx,SumIy,SumIz

  real(8)                                     ::   ubarold2,vbarold2,wbarold2,omegaxold2,omegayold2,omegazold2,accx,accy,accz
  real(8)                                     ::   anacxgh,anacygh,anaczgh
  real(8)                                     ::   xbarold,ybarold,zbarold,ubarold,vbarold,wbarold,omegaxold,omegayold,omegazold
  real(8)                                     ::   asolidxold,asolidyold,asolidzold,anacxold,anacyold,anaczold,alphaconv
  real(8)                                     ::   sumuold,sumvold,sumwold,sumIxold,SumIyold,SumIzold


  real(8),dimension(3)                        ::   load_Ref
  real(8)                                     ::   SumU_l,SumV_l,SumW_l,SumIx_l,SumIy_l,SumIz_l



  contains


    Subroutine Solid_Constant_Ini()
      implicit none

      Include "Par_Constant_Solid.txt"


      M_Cyl= pi*(r_Cyl_1**2)*rosolid*len_Cyl_1+ &
           & pi*(r_Cyl_2**2)*rocon  *len_Cyl_2
      print *, 'Hull mass=',M_Cyl/1000,'Tonnes Theoritical'

    end subroutine 


    Subroutine Solid_Finder_Dynamic_Ini(gap)
      implicit none
      real(8)                                 :: gap

      Ib_Cyl_1(:,:,:)=0
      Include "Par_Dynamic_Cyl_1.txt"
      call Solid_Cyl_Mark(point_Cyl_1,r_Cyl_1,gap,Ib_Cyl_1,2)

      Ic_Cyl_2(:,:,:)=0
      !Include "Par_Dynamic_Cyl_2.txt"
      !call Solid_Cyl_Mark(point_Cyl_2,r_Cyl_2,gap,Ic_Cyl_2,2)

      Ib_Sphe_1(:,:,:)=0
!      Include "Par_Dynamic_Sphere_1.txt"
!      call Solid_Sphe_Mark(point_Sphe_1,r_Sphe_1,gap,Ib_Sphe_1,1)

      Ib_Sphe_2(:,:,:)=0
!      Include "Par_Dynamic_Sphere_2.txt"
!      call Solid_Sphe_Mark(point_Sphe_2,r_Sphe_2,gap,Ib_Sphe_2,1)

      Ib_Box_1(:,:,:)=0
!      Include "Par_Dynamic_Box_1.txt"
!      call Solid_Box_Mark(point_Box_1,w_Box_1,Ib_Box_1,4)

      Ib_Box_2(:,:,:)=0
!      Include "Par_Dynamic_Box_2.txt"
!      call Solid_Box_Mark(point_Box_2,W_Box_2,Ib_Box_2,4) 


      Ib_solid(:,:,:)=max(IB_Cyl_1(:,:,:),Ib_Sphe_1(:,:,:),Ib_Sphe_2(:,:,:),Ib_Box_1(:,:,:),Ib_Box_2(:,:,:))
      Ic_Solid(:,:,:)=Ic_Cyl_2(:,:,:)
       
      

      
    end subroutine


    subroutine Solid_Fluid_Dynamic_Ini(ro) 
      implicit none
      real(8),intent (in),dimension (0:nx,0:ny,0:nz)          ::   ro


      call Cg_Mass_SolidFinder(ro)   !!it is more logical to update it after updateing the marker function
 

      sumu=0    ; sumv=0    ; sumw=0   !! Linear momentum
      sumIx=0   ; SumIy=0   ; SumIz=0  !! Angular momentum 
      ubar=0    ; vbar=0    ; wbar=0
      omegax=0  ; omegay=0  ; omegaz=0  
      asolidx=0 ; asolidy=0 ; asolidz=0
      anacx=0   ; anacy=0   ; anacz=0

    return
    end subroutine 


    subroutine Solid_Fluid(ro)
      use M_General,              only: u,v,w  
                            
      implicit none
      integer i,j,k,tp
      real(8),dimension(0:nx,0:ny,0:nz),intent(in)    :: ro 
      real(8) Ixy,Ixz,Iyz,Izz,Iyy,Ixx
      real(8) AAA,BBB,CCC,detmat,dm
      real(8),dimension (1:3,1:3)       :: Im

 
      call Cg_Mass_SolidFinder(ro)   !!it is more logical to update it after updateing the marker function
 
      sumu =0 ;sumv =0 ;sumw =0 
      Izz =0 ;Iyy  =0 ;Ixx  =0
      Ixy =0 ;Ixz  =0 ;Iyz  =0 
      sumIx=0 ;sumIy=0 ;sumIz=0 !!!sumIx=Hx , SumIy=Hy, SumIz=Hz !!
      do i=1,nx-1 
        do j=1,ny-1 
          do k=1,nz-1

            dm=0.125*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*(z(k+1)-z(k-1))*( Ib_Solid(i,j,k)*ro(i,j,k) )
    
            AAA=0.5*( u(i,j,k)+u(i+1,j,k)  )*dm 
            BBB=0.5*( v(i,j,k)+v(i,j+1,k)  )*dm  
            CCC=0.5*( w(i,j,k)+w(i,j,k+1)  )*dm  
            sumu=sumu+AAA 
            sumv=sumv+BBB 
            sumw=sumw+CCC

            sumIx=sumIx+ (y(j)-ybar)*CCC -(z(k)-zbar)*BBB
            sumIy=sumIy- (x(i)-xbar)*CCC +(z(k)-zbar)*AAA              !! mistake solve!!
            sumIz=sumIz+ (x(i)-xbar)*BBB -(y(j)-ybar)*AAA   

 
            Izz=Izz+ ( (x(i)-xbar)**2+ (y(j)-ybar)**2 )*dm 
            Iyy=Iyy+ ( (x(i)-xbar)**2+ (z(k)-zbar)**2 )*dm 
            Ixx=Ixx+ ( (y(j)-ybar)**2+ (z(k)-zbar)**2 )*dm 
            Ixy=Ixy+ ( (y(j)-ybar)   * (x(i)-xbar)    )*dm 
            Ixz=Ixz+ ( (x(i)-xbar)   * (z(k)-zbar)    )*dm        !! mistake solved!!
            Iyz=Iyz+ ( (y(j)-ybar)   * (z(k)-zbar)    )*dm

          end do 
        end do 
      end do 


      detmat=Ixx*(Izz*Iyy-Iyz*Iyz)-Ixy*(Izz*Ixy+Iyz*Ixz)-Ixz*(Iyz*Ixy+Iyy*Ixz)

      Im(1,1)=Izz*Iyy-Iyz*Iyz ;  Im(1,2)=Izz*Ixy+Iyz*Ixz ; Im(1,3)=Iyz*Ixy+Iyy*Ixz 
      Im(2,1)=Im(1,2)         ;  Im(2,2)=Izz*Ixx-Ixz*Ixz ; Im(2,3)=Iyz*Ixx+Ixy*Ixz
      Im(3,1)=Im(1,3)         ;  Im(3,2)=Im(2,3)         ; Im(3,3)=Iyy*Ixx-Ixy*Ixy

      omegax=0 !( Im(1,1)*sumIx +Im(1,2)*sumIy +Im(1,3)*sumIz )/detmat
      omegay=0 !( Im(2,1)*sumIx +Im(2,2)*sumIy +Im(2,3)*sumIz )/detmat
      omegaz=0 !( Im(3,1)*sumIx +Im(3,2)*sumIy +Im(3,3)*sumIz )/detmat

      ubar  =0 !sumu/Tmass 
      vbar  =0 ! sumv/Tmass
      wbar  =0 ! sumw/Tmass
                                                                     

      
      do i=2,nx-1 
        do j=1,ny-1 
          do k=1, nz-1

            u(i,j,k)=u(i,j,k)+0.50*( Ib_Solid(i-1,j,k)+Ib_Solid(i,j,k) )* &
            &                ( ( ubar + omegay*(z(k)-zbar)  -omegaz*(y(j)-ybar) )-u(i,j,k) ) 

          end do 
        end do 
      end do 

      do i=1,nx-1 
        do j=2,ny-1 
          do k=1, nz-1

            v(i,j,k)=v(i,j,k)+0.50*( Ib_Solid(i,j-1,k)+Ib_Solid(i,j,k) )* &
            &                ( ( vbar - omegax*(z(k)-zbar)  +omegaz*(x(i)-xbar) )-v(i,j,k) )

          end do 
        end do 
      end do 

      do i=1,nx-1 
        do j=1,ny-1 
          do k=2, nz-1

            w(i,j,k)=w(i,j,k)+0.50*( Ib_Solid(i,j,k-1)+Ib_Solid(i,j,k) )* &
            &                ( ( wbar + omegax*(y(j)-ybar)  -omegay*(x(i)-xbar) )-w(i,j,k) )

          end do 
        end do 
      end do
     

      asolidx=(ubar-ubarold2)/dt     
      asolidy=(vbar-vbarold2)/dt         
      asolidz=(wbar-wbarold2)/dt

      anacxgh=anacx                      
      anacygh=anacy                      
      anaczgh=anacz
       
      anacx=(omegax-omegaxold2 )/dt      
      anacy=(omegay-omegayold2 )/dt      
      anacz=(omegaz-omegazold2 )/dt

      alphaconv=abs(anacy-anacygh)
      !print*,"alpha convergence=",alphaconv
 
      anacx=0.02*anacx+ 0.98*anacxgh  
      anacy=0.02*anacy+ 0.98*anacygh 
      anacz=0.02*anacz+ 0.98*anaczgh     !! it can be used in if condition !!

      print*,"I(2,2)",Iyy 
 

      return
    End subroutine 





    subroutine Solid_Sec_Sav()
      implicit none 

      
      point_Cyl_old_1(:,:)=point_Cyl_1(:,:)      
      point_Cyl_old_2(:,:)=point_Cyl_2(:,:)
     
      point_Sphe_old_1(:,:)=point_Sphe_1(:,:)      
      point_Sphe_old_2(:,:)=point_Sphe_2(:,:)
      
      point_Box_old_1(:,:)=point_Box_1(:,:)
      point_Box_old_2(:,:)=point_Box_2(:,:)

      Ib_Solid_old(:,:,:)=Ib_Solid(:,:,:)
      Ic_Solid_old(:,:,:)=Ic_Solid(:,:,:)

      return

    end subroutine


    subroutine Solid_Fluid_Sec_Sav()
      implicit none 

      xbarold=xbar       
      ybarold=ybar       
      zbarold=zbar   

      ubarold=ubar       
      vbarold=vbar       
      wbarold=wbar

      omegaxold=omegax   
      omegayold=omegay   
      omegazold=omegaz 

      asolidxold=asolidx 
      asolidyold=asolidy 
      asolidzold=asolidz

      anacxold=anacx     
      anacyold=anacy     
      anaczold=anacz 

      sumuold=sumu  
      sumvold=sumv   
      sumwold=sumw 
 
      sumIxold=sumIx  
      SumIyold=SumIy    
      SumIzold=SumIz  

      return
    end subroutine





    subroutine Solid_It_Sav()
      implicit none 


      point_Cyl_old2_1(:,:)=point_Cyl_1(:,:)
      point_Cyl_old2_2(:,:)=point_Cyl_2(:,:)

      point_Sphe_old2_1(:,:)=point_Sphe_1(:,:)
      point_Sphe_old2_2(:,:)=point_Sphe_2(:,:)

      point_Box_old2_1(:,:)=point_Box_1(:,:)
      point_Box_old2_2(:,:)=point_Box_2(:,:)

      return
    end subroutine


subroutine Solid_Fluid_It_Sav()
implicit none 

      ubarold2=ubar     ; vbarold2=vbar    ; wbarold2=wbar 
      omegaxold2=omegax ;omegayold2=omegay ; omegazold2=omegaz 

end subroutine






    subroutine Solid_Sec_Cor()
      implicit none 

      point_Cyl_1(:,:)=0.5*(point_Cyl_1(:,:)+point_Cyl_old_1(:,:))      
      point_Cyl_2(:,:)=0.5*(point_Cyl_2(:,:)+point_Cyl_old_2(:,:))
           
      point_Sphe_1(:,:)=0.5*(point_Sphe_1(:,:)+point_Sphe_old_1(:,:))     
      point_Sphe_2(:,:)=0.5*(point_Sphe_2(:,:)+point_Sphe_old_2(:,:))
      
      point_Box_1(:,:)=0.5*(point_Box_1(:,:)+point_Box_old_1(:,:))
      point_Box_2(:,:)=0.5*(point_Box_2(:,:)+point_Box_old_2(:,:))

      Ib_solid(:,:,:)=0.5*(Ib_Solid(:,:,:)+Ib_Solid_old(:,:,:)) 
      Ic_solid(:,:,:)=0.5*(Ic_Solid(:,:,:)+Ic_Solid_old(:,:,:))

      return
    end subroutine


    subroutine Solid_Fluid_Sec_Cor()
      implicit none 

      xbar=0.5d0*(xbarold+xbar)          ; ybar=0.5d0*(ybarold+ybar)          ;     zbar=0.5d0*(zbarold+zbar)   
      ubar=0.5d0*(ubarold+ubar)          ; vbar=0.5d0*(vbarold+vbar)          ;     wbar=0.5d0*(wbarold+wbar)      
      omegax=0.5d0*(omegaxold+omegax)    ; omegay=0.5d0*(omegayold+omegay)    ;     omegaz=0.5d0*(omegazold+omegaz)
      asolidx=0.5d0*(asolidxold+asolidx) ; asolidy=0.5d0*(asolidyold+asolidy) ;     asolidz=0.5d0*(asolidzold+asolidz) 
      anacx=0.5d0*(anacxold+anacx)       ; anacy=0.5d0*(anacyold+anacy)       ;     anacz=0.5d0*(anaczold+anacz)    
      sumu=0.5d0*(sumuold+sumu)          ; sumv=0.5d0*(sumvold+sumv)          ;     sumw=0.5d0*(sumwold+sumw) 
      SumIx=0.5d0*(sumIxold+sumIx)       ; SumIy=0.5d0*(SumIyold+SumIy)       ;     SumIz=0.5d0*(SumIzold+SumIz) 

      return
    end subroutine



    Subroutine  Solid_Dynamic_Update(gap)
      implicit none 
      real(8)                        :: gap

      call Solid_Cyl_Mark(point_Cyl_1,r_Cyl_1,gap,Ib_Cyl_1,2)
      call Solid_Update_Posi(point_Cyl_1,point_Cyl_old2_1,2)

      !call Solid_Cyl_Mark(point_Cyl_2,r_Cyl_2,gap,Ic_Cyl_2,2)
      !call Solid_Update_Posi(point_Cyl_2,point_Cyl_old2_2,2)

      !call Solid_Sphe_Mark(point_Sphe_1,r_Sphe_1,gap,Ib_Sphe_1,1)
      !call Solid_Update_Posi(point_Sphe_1,point_Sphe_old2_1,1)

      !call Solid_Sphe_Mark(point_Sphe_2,r_Sphe_2,gap,Ib_Sphe_2,1)
      !call Solid_Update_Posi(point_Sphe_2,point_Sphe_old2_2,1)

      !call Solid_Box_Mark(point_Box_1,w_Box_1,Ib_Box_1,4)
      !call Solid_Update_Posi(point_Box_1,point_Box_old2_1,4)

      !call Solid_Box_Mark(point_Box_2,w_Box_2,Ib_Box_2,4)
      !call Solid_Update_Posi(point_Box_2,point_Box_old2_2,4)

 

      Ib_solid(:,:,:)=max(IB_Cyl_1(:,:,:),Ib_Sphe_1(:,:,:),Ib_Sphe_2(:,:,:),Ib_Box_1(:,:,:),Ib_Box_2(:,:,:))
      Ic_Solid(:,:,:)=Ic_Cyl_2(:,:,:)


      return
    end subroutine


    subroutine Solid_Update_Posi(point_Gen,point_Gen_old2,num)
      implicit none 
      integer,intent(in)                            ::   num
      real(8),intent(in), dimension(num,3)          ::   point_Gen_old2
      real(8),intent(out),dimension(num,3)          ::   point_Gen

      real(8)            ,dimension(num,3)          ::   point_Gen_n
      integer                                       ::   i

      do i=1,num
        point_Gen_n(i,1)=point_Gen_old2(i,1) +dt*( ubar + omegay*(point_Gen(i,3)-zbar)  -omegaz*(point_Gen(i,2)-ybar) )
        point_Gen_n(i,2)=point_Gen_old2(i,2) +dt*( vbar - omegax*(point_Gen(i,3)-zbar)  +omegaz*(point_Gen(i,1)-xbar) )
        point_Gen_n(i,3)=point_Gen_old2(i,3) +dt*( wbar + omegax*(point_Gen(i,2)-ybar)  -omegay*(point_Gen(i,1)-xbar) )
      end do 

      point_Gen(:,:)=point_Gen_n(:,:)
      
      return
    end subroutine 
 
    subroutine Cg_Mass_SolidFinder(ro)
      implicit none 
 
      real(8),intent (in),dimension (0:nx,0:ny,0:nz)          ::   ro

      real(8) dm,Txm,Tym,Tzm
      integer i,j,k

      Tmass=0
      Txm=0
      Tym=0
      Tzm=0
      do i=1,nx-1 
        do j=1,ny-1 
          do k=1,nz-1 

            dm=0.125*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*(z(k+1)-z(k-1))*( Ib_Solid(i,j,k)*ro(i,j,k) )
            Tmass=Tmass+dm
            Txm=Txm+x(i)*dm
            Tym=Tym+y(j)*dm
            Tzm=Tzm+z(k)*dm

          end do 
        end do 
      end do 

      xbar=Txm/Tmass
      ybar=Tym/Tmass 
      zbar=Tzm/Tmass

      print*, "Platfrom cg in x,y,z",xbar,ybar,zbar
      print*,"mass numerical=",Tmass
      
      return
    end subroutine 

    subroutine Solid_load_Ini(ro)
      use M_General,              only: u,v,w  
      implicit none
 
      real(8),dimension(0:nx,0:ny,0:nz),intent(in)    ::   ro 

      real(8)                                         ::   AAA,BBB,CCC,dm
      integer                                         ::   i,j,k

      sumu_l =0 ;sumv_l=0 ;sumw_l =0 
      sumIx_l=0 ;sumIy_l=0 ;sumIz_l=0 !!!sumIx=Hx , SumIy=Hy, SumIz=Hz !!
      do i=1,nx-1 
        do j=1,ny-1 
          do k=1,nz-1

            dm=0.125*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*(z(k+1)-z(k-1))*( Ib_Solid(i,j,k)*ro(i,j,k) )
    
            AAA=0.5*( u(i,j,k)+u(i+1,j,k)  )*dm 
            BBB=0.5*( v(i,j,k)+v(i,j+1,k)  )*dm  
            CCC=0.5*( w(i,j,k)+w(i,j,k+1)  )*dm  
            sumu_l=sumu_l+AAA 
            sumv_l=sumv_l+BBB 
            sumw_l=sumw_l+CCC

            sumIx_l=sumIx_l+ (y(j)-Load_Ref(2))*CCC -(z(k)-Load_Ref(3))*BBB
            sumIy_l=sumIy_l- (x(i)-Load_Ref(1))*CCC +(z(k)-Load_Ref(3))*AAA              !! mistake solve!!
            sumIz_l=sumIz_l+ (x(i)-Load_Ref(1))*BBB -(y(j)-Load_Ref(2))*AAA   

          end do 
        end do 
      end do
 
      return 
    end subroutine 


    subroutine Solid_load_Mod(ro)
      use M_General,              only: u,v,w  
      implicit none

      real(8),dimension(0:nx,0:ny,0:nz),intent(in)    ::   ro
 
      real(8)                                         ::   SumU_lE,SumV_lE,SumW_lE,SumIx_lE,SumIy_lE,SumIz_lE
      real(8)                                         ::   AAA,BBB,CCC,dm
      integer                                         ::   i,j,k


      sumu_lE =0 ;sumv_lE =0 ;sumw_lE =0 
      sumIx_lE=0 ;sumIy_lE=0 ;sumIz_lE=0 !!!sumIx=Hx , SumIy=Hy, SumIz=Hz !!
      do i=1,nx-1 
        do j=1,ny-1 
          do k=1,nz-1


            dm=0.125*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*(z(k+1)-z(k-1))*( Ib_Solid(i,j,k)*ro(i,j,k) )
    
            AAA=0.5*( u(i,j,k)+u(i+1,j,k)  )*dm 
            BBB=0.5*( v(i,j,k)+v(i,j+1,k)  )*dm  
            CCC=0.5*( w(i,j,k)+w(i,j,k+1)  )*dm  
            sumu_lE=sumu_lE+AAA 
            sumv_lE=sumv_lE+BBB 
            sumw_lE=sumw_lE+CCC

            sumIx_lE=sumIx_lE+ (y(j)-Load_Ref(2))*CCC -(z(k)-Load_Ref(3))*BBB
            sumIy_lE=sumIy_lE- (x(i)-Load_Ref(1))*CCC +(z(k)-Load_Ref(3))*AAA              !! mistake solve!!
            sumIz_lE=sumIz_lE+ (x(i)-Load_Ref(1))*BBB -(y(j)-Load_Ref(2))*AAA   

            
          end do 
        end do 
      end do 

      
      sumu_l =sumu_l -sumu_lE 
      sumv_l =sumv_l -sumv_lE 
      sumw_l =sumw_l -sumw_lE 
      sumIx_l=sumIx_l-sumIx_lE    
      sumIy_l=sumIy_l-sumIy_lE    
      sumIz_l=sumIz_l-sumIz_lE    

      
      print*,sumu_lE/sumu_l*100,"% error drag" 
      print*,sumIy_lE/sumIy_l*100,"% error moment"


      return
    end subroutine 



    subroutine Solid_load_Final(ro,p,Advectu,Advectv,AdvectW,Tx,Ty,Tz,tpstar)
      use M_General,              only: gz  
      implicit none
      real(8),intent(in),dimension (1:nx,0:ny,0:nz)          ::   advectu
      real(8),intent(in),dimension (0:nx,1:ny,0:nz)          ::   advectv
      real(8),intent(in),dimension (0:nx,0:ny,1:nz)          ::   advectw

      real(8),intent(in),dimension (2:nx-1,1:ny-1,1:nz-1)    ::   Tx   
      real(8),intent(in),dimension (1:nx-1,2:ny-1,1:nz-1)    ::   Ty
      real(8),intent(in),dimension (1:nx-1,1:ny-1,2:nz-1)    ::   Tz
      real(8),intent(in),dimension (0:nx,0:ny,0:nz)          ::   ro,p
      integer,intent(in)                                     ::   tpstar

      real(8)                                                ::   Drag,Lift,MDrag,MLift,Drag_P,MDrag_P
      integer                                                ::   i,j,k


      Drag=0 ; Lift=0 ; MDrag=0 ; MLift=0
      Drag_P=0 ;  MDrag_P=0 
      do i=1,nx-1 
        do j=1,ny-1 
          do k=1,nz-1


            Drag=Drag+ 0.5*(ro(i,j,k)+ro(i-1,j,k))*0.5* (Ib_Solid(i,j,k)+Ib_Solid(i-1,j,k))* &
            & (    x(i)-x(i-1))*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )*( 0.5*(z(k)+z(k+1))-0.5*(z(k)+z(k-1)) )*  &
            & (    +advectu(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )   )*( p(i,j,k)-p(i-1,j,k) )/(x(i)-x(i-1))+Tx(i,j,k) )

            Drag_p=Drag_p+ 0.5* (Ib_Solid(i,j,k)+Ib_Solid(i-1,j,k))* &
            &     ( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )*( 0.5*(z(k)+z(k+1))-0.5*(z(k)+z(k-1)) )* (    -( p(i,j,k)-p(i-1,j,k) ) )

            Lift=Lift+ 0.5*(ro(i,j,k)+ro(i,j,k-1))*0.5*(Ib_Solid(i,j,k)+Ib_Solid(i,j,k-1))* &
            & (   0.5*(x(i)+x(i+1))-0.5*(x(i)+x(i-1)) )*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )*(z(k)-z(k-1))*  &
            & (   +advectw(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i,j,k-1)))   )*( p(i,j,k)-p(i,j,k-1) )/(z(k)-z(k-1))+Tz(i,j,k)+gz )

            MDrag_P=MDrag_P+(Z(k)-Load_Ref(3))*(0.5* (Ib_Solid(i,j,k)+Ib_Solid(i-1,j,k))* &
            &     ( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )*( 0.5*(z(k)+z(k+1))-0.5*(z(k)+z(k-1)) )* (    -( p(i,j,k)-p(i-1,j,k) ) ) )

            MDrag=MDrag+(Z(k)-Load_Ref(3))*(0.5*(ro(i,j,k)+ro(i-1,j,k))*0.5* (Ib_Solid(i,j,k)+Ib_Solid(i-1,j,k))* &
            & (x(i)-x(i-1))*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )*( 0.5*(z(k)+z(k+1))-0.5*(z(k)+z(k-1)) )*  &
            & ( +advectu(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i-1,j,k))  )   )*( p(i,j,k)-p(i-1,j,k) )/(x(i)-x(i-1))+Tx(i,j,k) ) )

            MLift=MLift+(x(i)-Load_Ref(1))*(0.5*(ro(i,j,k)+ro(i,j,k-1))*0.5*(Ib_Solid(i,j,k)+Ib_Solid(i,j,k-1))* &
            & (   0.5*(x(i)+x(i+1))-0.5*(x(i)+x(i-1)) )*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )*(z(k)-z(k-1))*  &
            & (   +advectw(i,j,k)-(  1/(  0.5*(ro(i,j,k)+ro(i,j,k-1)))   )*( p(i,j,k)-p(i,j,k-1) )/(z(k)-z(k-1))+Tz(i,j,k)+gz ) )


          end do 
        end do 
      end do
      

      write(75,2013)  tpstar*dt,Drag,(sumu_l-0)/dt,Lift,(sumw_l-0)/dt,MDrag,MLift,(sumIy_l-0)/dt,Drag_P,MDrag_P
      2013  format (10(1x,e15.7))
      call flush (75)


      return 
    end subroutine 





Subroutine Solid_Dynamic_Write()
   
  
   implicit none 
   integer i,j,k

502  format (2(1x,e23.15))  
511  format (3(1x,e23.15))
522  format (1x,e23.15)

   OPEN(unit=3015,file=fileplace//"Platform_Dynamic_Save.dat",STATUS='REPLACE')   

   do k=0,nz 
     do j=0,ny
       do i=0,nx
         write(3015,502) Ib_Solid(i,j,k),Ic_Solid(i,j,k)
       end do 
     end do
   end do 
   
     do i=1,2
       write (3015,511) point_Cyl_1(i,1),point_Cyl_1(i,2),point_Cyl_1(i,3)
       write (3015,511) point_Cyl_2(i,1),point_Cyl_2(i,2),point_Cyl_2(i,3)
     end do
     
     do i=1,1
       write (3015,511) point_Sphe_1(i,1),point_Sphe_1(i,2),point_Sphe_1(i,3)
       write (3015,511) point_Sphe_2(i,1),point_Sphe_2(i,3),point_Sphe_2(i,3)
     end do 

     do i=1,4
       write (3015,511) point_Box_1(i,1),point_Box_1(i,2),point_Box_1(i,3)
       write (3015,511) point_Box_2(i,1),point_Box_2(i,2),point_Box_2(i,3)
     end do
 



   write(3015, 522) Xbar
   write(3015, 522) Ybar
   write(3015, 522) Zbar
   write(3015, 522) Ubar
   write(3015, 522) Vbar
   write(3015, 522) Wbar
   write(3015, 522) OmegaX
   write(3015, 522) OmegaY
   write(3015, 522) OmegaZ
   write(3015, 522) aSolidX
   write(3015, 522) aSolidY
   write(3015, 522) aSolidZ
   write(3015, 522) anacX
   write(3015, 522) anacY
   write(3015, 522) anacZ
   write(3015, 522) sumU
   write(3015, 522) SumV
   write(3015, 522) SumW
   write(3015, 522) SumIx
   write(3015, 522) SumIy
   write(3015, 522) SumIz
    
   call flush (3015)   
   close(3015)

   return
end subroutine Solid_Dynamic_Write


Subroutine Solid_Dynamic_Read()
   
   implicit none 
   integer i,j,k

502  format (2(1x,e23.15))  
511  format (3(1x,e23.15))
522  format (1x,e23.15)


   OPEN(unit=3015,file=fileplace//"Platform_Dynamic_Save.dat")   

   do k=0,nz 
     do j=0,ny
       do i=0,nx
         read(3015,502) Ib_Solid(i,j,k),Ic_Solid(i,j,k)
       end do 
     end do
   end do 
   
   
     do i=1,2
       read (3015,511) point_Cyl_1(i,1),point_Cyl_1(i,2),point_Cyl_1(i,3)
       read (3015,511) point_Cyl_2(i,1),point_Cyl_2(i,2),point_Cyl_2(i,3)
     end do
     
     do i=1,1
       read (3015,511) point_Sphe_1(i,1),point_Sphe_1(i,2),point_Sphe_1(i,3)
       read (3015,511) point_Sphe_2(i,1),point_Sphe_2(i,3),point_Sphe_2(i,3)
     end do 

     do i=1,4
       read (3015,511) point_Box_1(i,1),point_Box_1(i,2),point_Box_1(i,3)
       read (3015,511) point_Box_2(i,1),point_Box_2(i,2),point_Box_2(i,3)
     end do
   




   read(3015, 522) Xbar
   read(3015, 522) Ybar
   read(3015, 522) Zbar
   read(3015, 522) Ubar
   read(3015, 522) Vbar
   read(3015, 522) Wbar
   read(3015, 522) OmegaX
   read(3015, 522) OmegaY
   read(3015, 522) OmegaZ
   read(3015, 522) aSolidX
   read(3015, 522) aSolidY
   read(3015, 522) aSolidZ
   read(3015, 522) anacX
   read(3015, 522) anacY
   read(3015, 522) anacZ
   read(3015, 522) sumU
   read(3015, 522) SumV
   read(3015, 522) SumW
   read(3015, 522) SumIx
   read(3015, 522) SumIy
   read(3015, 522) SumIz

   close(3015)

return
end subroutine Solid_Dynamic_Read

 

    


end module 
