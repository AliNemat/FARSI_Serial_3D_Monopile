module M_Graphic
use M_General, only: r2
use M_Platform_Constant 
use M_Platform_Dynamic                      
use M_Tether, only: Ox,Oy,Oz,DTethI
use M_Tower,  only: RSpinV,TurbH,DiT,DiB,RRot
use M_Math,   only: ROTATIONX,ROTATIONY,ROTATIONZ,LENGTHF

contains 

subroutine solidmeshgen (tp,dt,gapz)
     

implicit none 

integer total,np,Normal,refpoint

parameter (total=32,np=360,Normal=7, refpoint=2) 

    real(8),dimension   (1:total)     :: ppx,ppy,ppz,R,ap,bp,cp,a2p,b2p,c2p     
    real(8),dimension   (1:20000)     :: x,y,z
    real(8),dimension   (1:Normal)    :: aprim,bprim,cprim,a2prim,b2prim,c2prim,nnx,nny,nnz,dp,xp,yp,zp   
   
    real(8),dimension (1:refpoint)    :: pxref,pyref,pzref
    
    
    real(8) dt,len2,gapz,hh,gammma,teta1,teta2
    
    real(8) RHub,rnac,blade,Rnacbig,lenhnac,gapnac,Rblade,RSpoke
 
    
    integer i,j,k,tp,num,numc(1:total),L,M
    
  
    write(tp+10000,*) 'ZONE T="',-tp,'" STRANDID=',tp+100000, 'SOLUTIONTIME=',tp*dt, 'DATAPACKING=POINT' 
    write(tp+10000,*) 'NODES=',Total*(np+1),'ELEMENTS=',2*total*np,'ZONETYPE=FETRIANGLE'  
    
    
    
    
    Rspoke=Froude*1.5 
    lenhnac=Froude*6  !! half of nacelle length
    Rnac=Froude*2.5 !! 2  !! smaller side  of the naclle radius.Other side is 3 times larger
    RHub=Froude*0.1*Rnac !30      !63  !0.63
    gapnac=Froude*2.75 ; Rnacbig=Rnac
    Rblade=Froude*1.7
    blade=Froude*2.0*Rblade !0.008 !! tickness of modeling baldes !! 
   
   
    
    
     
    
    gammma=RspinV*360/60*tp*dt    !! in degree !!
    teta1=120.0
    teta2=240.0
  
    R(2)=0.5*DiB 
    R(1)=0.5*DiT
    R(3:4)=r2
    R(5:12)=0.5*DTethI   !! very 
   
 
    R(13:14)=Rnacbig
    R(15)=Rnac
    R(16)=R(13)
    R(17)=Rnac
    R(18)=RHub
    R(19:26)=Rspoke
    R(27)=Rblade ; R(28)=0.7/1.7*Rblade
    R(29)=Rblade ; R(30)=0.7/1.7*Rblade
    R(31)=Rblade ; R(32)=0.7/1.7*Rblade
    
     
     
!     
!     !! i'=a'i+b'j+c'k
!     !! j'=a"i+b"j+c"k
!     !! general point on the circle  (x',y')=rcos(@)i'+rsin(@)j'=>(x,y,z)=rcos(@)(a'i+b'j+c'k)+rsin(@)(a"i+b"j+c"k)
!     !! the cordinate system is the same just the center of coordinate is different  which we change here in a do loop!!! and then turing around the cirlce !! 

! 
     
     do i=1,4 
    
               if (i.eq.1) then
                        j=1 ; k=1 ; L=2 ; M=2  
                 else if (i.eq.2) then
                         j=5 ; k=4 ; L=6 ; M=4  
                 else if (i.eq.3) then
                         j=1 ; k=4 ; L=4 ; M=4
                else 
                         j=2 ; k=4 ; L=3 ; M=4
               end if 
     
     nnx(i)=(px(j,k)-px(L,M))/LENGTHF(px(j,k)-px(L,M),py(j,k)-py(L,M),pz(j,k)-pz(L,M))
     nny(i)=(py(j,k)-py(L,M))/LENGTHF(px(j,k)-px(L,M),py(j,k)-py(L,M),pz(j,k)-pz(L,M))
     nnz(i)=(pz(j,k)-pz(L,M))/LENGTHF(px(j,k)-px(L,M),py(j,k)-py(L,M),pz(j,k)-pz(L,M))
     dp(i)=px(j,k)*nnx(i)+py(j,k)*nny(i)+pz(j,k)*nnz(i)
   
     if (nny(i).ne.0.0   ) then 
          xp(i)=px(j,k)+10 ; zp(i)=pz(j,k)+12
          yp(i)=(dp(1)-nnx(i)*xp(i)-nnz(i)*zp(i))/nny(i)
     else if  (nnz(i).ne.0.0 ) then  
          xp(i)=px(j,k)+10 ; yp(i)=py(j,k)+12 
          Zp(i)=(dp(i)-nnx(i)*xp(i)-nny(i)*yp(i))/nnz(i)
     else
          xp(i)=px(j,k) ; yp(i)=py(j,k)+10 ; zp(i)=pz(j,k)+12   
     end if   
     
     aprim(i)=(xp(i)-px(j,k)) /LENGTHF( xp(i)-px(j,k), yp(i)-py(j,k), zp(i)-pz(j,k))  
     bprim(i)=(yp(i)-py(j,k)) /LENGTHF( xp(i)-px(j,k), yp(i)-py(j,k), zp(i)-pz(j,k))  
     cprim(i)=(zp(i)-pz(j,k)) /LENGTHF( xp(i)-px(j,k), yp(i)-py(j,k), zp(i)-pz(j,k))     
     a2prim(i)= (nny(i)*cprim(i)-nnz(i)*bprim(i))
     b2prim(i)=-(nnx(i)*cprim(i)-nnz(i)*aprim(i))
     c2prim(i)= (nnx(i)*bprim(i)-nny(i)*aprim(i))
     
     end do 
     
     
     
   
     ap(1:12)=aprim(1)   ;    bp (1:12)=bprim(1) ;    cp (1:12)=cprim(1)
     a2p(1:12)=a2prim(1) ;   b2p(1:12)=b2prim(1) ;    c2p(1:12)=c2prim(1)
     
     ap (13:18)=aprim(2)   ;    bp (13:18)=bprim(2) ;   cp (13:18)=cprim(2)
     a2p(13:18)=a2prim(2) ;    b2p(13:18)=b2prim(2) ;    c2p(13:18)=c2prim(2)
     
     ap (19:20)=aprim(3)   ;    bp (19:20)=bprim(3) ;   cp (19:20)=cprim(3)
     a2p(19:20)=a2prim(3) ;    b2p(19:20)=b2prim(3)  ;   c2p(19:20)=c2prim(3)
     
     ap(21:24)=aprim(4)    ;    bp (21:24)=bprim(4) ;   cp (21:24)=cprim(4)
     a2p(21:24)=a2prim(4) ;    b2p(21:24)=b2prim(4)  ;   c2p(21:24)=c2prim(4)
     
     ap(25:26)=aprim(3)    ;    bp (25:26)=bprim(3) ;   cp (25:26)=cprim(3)
     a2p(25:26)=a2prim(3) ;     b2p(25:26)=b2prim(3)  ;   c2p(25:26)=c2prim(3)
     
    
     
     
     
    ppz(3)=pz(1,1)
    ppy(3)=py(1,1)
    ppx(3)=px(1,1)
    
    ppz(4)=pz(2,2)
    ppy(4)=py(2,2)
    ppx(4)=px(2,2) 
    
    ppz(2)=pz(1,1)
    ppy(2)=py(1,1)
    ppx(2)=px(1,1)
    
    ppz(1)=pz(2,2)+ ( pz(1,1)-pz(2,2) )*(TurbH+len)/len
    ppy(1)=py(2,2)+ ( py(1,1)-py(2,2) )*(TurbH+len)/len  
    ppx(1)=px(2,2)+ ( px(1,1)-px(2,2) )*(TurbH+len)/len
  
   
     ppx(13)=ppx(1)+Rnacbig* nnx(1)-gapnac*nnx(2)  !! slightly goes up in tower direction and goes back in nacelle direction 
     ppy(13)=ppy(1)+Rnacbig* nny(1)-gapnac*nny(2)  !! Yaw is bring to account 
     ppz(13)=ppz(1)+Rnacbig* nnz(1)-gapnac*nnz(2)
    
     ppx(14)=ppx(13)+1.5*Lenhnac*nnx(2)
     ppy(14)=ppy(13)+1.5*Lenhnac*nny(2)
     ppz(14)=ppz(13)+1.5*Lenhnac*nnz(2)
     
     ppx(15)=ppx(13)-0.5*Lenhnac*nnx(2)
     ppy(15)=ppy(13)-0.5*Lenhnac*nny(2)
     ppz(15)=ppz(13)-0.5*Lenhnac*nnz(2)
    
     
     ppx(16)=ppx(13)
     ppy(16)=ppy(13)
     ppz(16)=ppz(13)
    
     ppx(17)=ppx(14)
     ppy(17)=ppy(14)
     ppz(17)=ppz(14)
     
     ppx(18)=ppx(17)+blade*nnx(2)
     ppy(18)=ppy(17)+blade*nny(2)
     ppz(18)=ppz(17)+blade*nnz(2)
     
     ppx(19)=px(1,4)+leg *(px(2,2)-px(1,4))/LENGTHF(px(2,2)-px(1,4),py(2,2)-py(1,4),pz(2,2)-pz(1,4))
     ppy(19)=py(1,4)+leg *(py(2,2)-py(1,4))/LENGTHF(px(2,2)-px(1,4),py(2,2)-py(1,4),pz(2,2)-pz(1,4))
     ppz(19)=pz(1,4)+leg *(pz(2,2)-pz(1,4))/LENGTHF(px(2,2)-px(1,4),py(2,2)-py(1,4),pz(2,2)-pz(1,4))
     
     ppx(20)=px(1,4)
     ppy(20)=py(1,4)
     ppz(20)=pz(1,4)
     
     
     ppx(21)=px(2,4)+leg *(px(2,2)-px(2,4))/LENGTHF(px(2,2)-px(2,4),py(2,2)-py(2,4),pz(2,2)-pz(2,4))
     ppy(21)=py(2,4)+leg *(py(2,2)-py(2,4))/LENGTHF(px(2,2)-px(2,4),py(2,2)-py(2,4),pz(2,2)-pz(2,4))
     ppz(21)=pz(2,4)+leg *(pz(2,2)-pz(2,4))/LENGTHF(px(2,2)-px(2,4),py(2,2)-py(2,4),pz(2,2)-pz(2,4))
     
     ppx(22)=px(2,4)
     ppy(22)=py(2,4)
     ppz(22)=pz(2,4)
     
       
     ppx(23)=px(3,4)+leg *(px(2,2)-px(3,4))/LENGTHF(px(2,2)-px(3,4),py(2,2)-py(3,4),pz(2,2)-pz(3,4))
     ppy(23)=py(3,4)+leg *(py(2,2)-py(3,4))/LENGTHF(px(2,2)-px(3,4),py(2,2)-py(3,4),pz(2,2)-pz(3,4))
     ppz(23)=pz(3,4)+leg *(pz(2,2)-pz(3,4))/LENGTHF(px(2,2)-px(3,4),py(2,2)-py(3,4),pz(2,2)-pz(3,4))
     
     ppx(24)=px(3,4)
     ppy(24)=py(3,4)
     ppz(24)=pz(3,4)
     
     
     !! fourth spoke and tether 
     
     !ppx(25)=px(4,4)+leg *(px(2,2)-px(4,4))/LENGTHF(px(2,2)-px(4,4),py(2,2)-py(4,4),pz(2,2)-pz(4,4))
     !ppy(25)=py(4,4)+leg *(py(2,2)-py(4,4))/LENGTHF(px(2,2)-px(4,4),py(2,2)-py(4,4),pz(2,2)-pz(4,4))
     !ppz(25)=pz(4,4)+leg *(pz(2,2)-pz(4,4))/LENGTHF(px(2,2)-px(4,4),py(2,2)-py(4,4),pz(2,2)-pz(4,4))
     
     !ppx(26)=px(4,4)
     !ppy(26)=py(4,4)
     !ppz(26)=pz(4,4)
      
      !! put it on top of the third one 
     ppx(25)=px(3,4)+leg *(px(2,2)-px(3,4))/LENGTHF(px(2,2)-px(3,4),py(2,2)-py(3,4),pz(2,2)-pz(3,4))
     ppy(25)=py(3,4)+leg *(py(2,2)-py(3,4))/LENGTHF(px(2,2)-px(3,4),py(2,2)-py(3,4),pz(2,2)-pz(3,4))
     ppz(25)=pz(3,4)+leg *(pz(2,2)-pz(3,4))/LENGTHF(px(2,2)-px(3,4),py(2,2)-py(3,4),pz(2,2)-pz(3,4))
     
     ppx(26)=px(3,4)
     ppy(26)=py(3,4)
     ppz(26)=pz(3,4)
     
     !!!!!!!!!!!!!!!!!!!!!! blade section !!!!!!!!!!!!!!!!!!!!!
     pxref(1)=ppx(14)+Rnacbig *nnx(1)-Rblade*nnx(2)
     pyref(1)=ppy(14)+Rnacbig *nny(1)-Rblade*nny(2)
     pzref(1)=ppz(14)+Rnacbig *nnz(1)-Rblade*nnz(2)
     
     pxref(2)=pxref(1)+RRot*nnx(1)
     pyref(2)=pyref(1)+RRot*nny(1)
     pzref(2)=pzref(1)+RRot*nnz(1)
     
     
     
     
     ppx(27)=ROTATIOnx( pxref(1),pyref(1),pzref(1),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),gammma)  !!blade root !!
     ppy(27)=ROTATIOny( pxref(1),pyref(1),pzref(1),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),gammma)
     ppz(27)=ROTATIONZ( pxref(1),pyref(1),pzref(1),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),gammma)
     
     ppx(28)=ROTATIONX( pxref(2),pyref(2),pzref(2),nnx(2),nnY(2),nnz(2),ppx(14),ppy(14),ppz(14),gammma)       !! blade tip!!
     ppy(28)=ROTATIOny( pxref(2),pyref(2),pzref(2),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),gammma)
     ppz(28)=ROTATIONZ( pxref(2),pyref(2),pzref(2),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),gammma)
     
     ppx(29)=ROTATIONX( ppx(27), ppy(27), ppz(27),nnx(2), nnY(2),nnz(2),ppx(14),ppy(14),ppz(14),teta1)
     ppY(29)=ROTATIONY( ppx(27), ppy(27), ppz(27),nnx(2), nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta1)
     ppZ(29)=ROTATIONZ( ppx(27), ppy(27), ppz(27),nnx(2), nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta1)
     
     
     
     ppx(30)=ROTATIONX( ppx(28), ppy(28), ppz(28),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta1)
     ppY(30)=ROTATIONY( ppx(28), ppy(28), ppz(28),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta1)
     ppZ(30)=ROTATIONZ( ppx(28), ppy(28), ppz(28),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta1)
    
     ppx(31)=ROTATIONX( ppx(27), ppy(27), ppz(27),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta2)
     ppY(31)=ROTATIONY( ppx(27), ppy(27), ppz(27),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta2)
     ppZ(31)=ROTATIONZ( ppx(27), ppy(27), ppz(27),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta2)
     
     ppx(32)=ROTATIONX( ppx(28), ppy(28), ppz(28),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta2)
     ppY(32)=ROTATIONY( ppx(28), ppy(28), ppz(28),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta2)
     ppZ(32)=ROTATIONZ( ppx(28), ppy(28), ppz(28),nnx(2),nny(2),nnz(2),ppx(14),ppy(14),ppz(14),teta2)
    
     
     !!!!!!!!!!!!!!!!! end of blade section !!!!!!!!!!!!!!!!!
     
     
    do i=5,11,2
    ppx(i)=ppx(15+i)
    ppy(i)=ppy(15+i)
    ppz(i)=ppz(15+i)
    end do
    
     ppx(6)=ox(1)
     ppy(6)=oy(1)
     ppz(6)=oz(1)+gapz
     
    
    
     ppx(8)=ox(2)
     ppy(8)=oy(2)
     ppz(8)=oz(2)+ gapz
   
     ppx(10)=ox(3)
     ppy(10)=oy(3)
     ppz(10)=oz(3)+gapz
    
     !! fourth tether 
     !ppx(12)=ox(4)
     !ppy(12)=oy(4)
     !ppz(12)=oz(4)+gapz
     
     !! put it on top of the third one 
     ppx(12)=ox(3)
     ppy(12)=oy(3)
     ppz(12)=oz(3)+gapz
     !!!!!!!!!!! end of defining point !!!


     
     !!!!!!!!!! blade section !!!!!!!!!!!!!!
    
     do i=5,7 
    
               if (i.eq.5) then
                        j=27  
                 else if (i.eq.6) then
                        j=29 
                 else if (i.eq.7) then
                        j=31
                
               end if 
     
     nnx(i)=(ppx(j+1)-ppx(j))/LENGTHF(ppx(j+1)-ppx(j),ppy(j+1)-ppy(j),ppz(j+1)-ppz(j))
     nny(i)=(ppy(j+1)-ppy(j))/LENGTHF(ppx(j+1)-ppx(j),ppy(j+1)-ppy(j),ppz(j+1)-ppz(j))
     nnz(i)=(ppz(j+1)-ppz(j))/LENGTHF(ppx(j+1)-ppx(j),ppy(j+1)-ppy(j),ppz(j+1)-ppz(j))
     dp(i)=ppx(j+1)*nnx(i)+ppy(j+1)*nny(i)+ppz(j+1)*nnz(i)
   
     if (nny(i).ne.0.0   ) then 
          xp(i)=ppx(j+1)+10 ; zp(i)=ppz(j+1)+12
          yp(i)=(dp(1)-nnx(i)*xp(i)-nnz(i)*zp(i))/nny(i)
     else if  (nnz(i).ne.0.0 ) then  
          xp(i)=ppx(j+1)+10 ; yp(i)=ppy(j+1)+12 
          Zp(i)=(dp(i)-nnx(i)*xp(i)-nny(i)*yp(i))/nnz(i)
     else
          xp(i)=ppx(j+1) ; yp(i)=ppy(j+1)+10 ; zp(i)=ppz(j+1)+12   
     end if   
     
     aprim(i)=(xp(i)-ppx(j+1)) /LENGTHF( xp(i)-ppx(j+1), yp(i)-ppy(j+1), zp(i)-ppz(j+1))  
     bprim(i)=(yp(i)-ppy(j+1)) /LENGTHF( xp(i)-ppx(j+1), yp(i)-ppy(j+1), zp(i)-ppz(j+1))  
     cprim(i)=(zp(i)-ppz(j+1)) /LENGTHF( xp(i)-ppx(j+1), yp(i)-ppy(j+1), zp(i)-ppz(j+1))    
     a2prim(i)= (nny(i)*cprim(i)-nnz(i)*bprim(i))
     b2prim(i)=-(nnx(i)*cprim(i)-nnz(i)*aprim(i))
     c2prim(i)= (nnx(i)*bprim(i)-nny(i)*aprim(i))
     
     end do 
     
     ap(27:28)=aprim(5)    ;  bp (27:28)=bprim(5)   ;   cp (27:28)=cprim(5)
    a2p(27:28)=a2prim(5)   ;  b2p(27:28)=b2prim(5)  ;   c2p(27:28)=c2prim(5)
    
    ap(29:30)=aprim(6)    ;  bp (29:30)=bprim(6)   ;   cp (29:30)=cprim(6)
    a2p(29:30)=a2prim(6)   ;  b2p(29:30)=b2prim(6)  ;   c2p(29:30)=c2prim(6)
    
    ap(31:32)=aprim(7)    ;  bp (31:32)=bprim(7)   ;   cp (31:32)=cprim(7)
    a2p(31:32)=a2prim(7)   ;  b2p(31:32)=b2prim(7)  ;   c2p(31:32)=c2prim(7)
    
    !!!!!!!!!!!!!!!!!!!!!! end of blade section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     do i=1,total  
    numc(i)=i-1
    end do
    x=0 ; y=0 ; z=0 
    
   do k=1,total
     x(numc(k)*(np+1)+1)=ppx(k)
     y(numc(k)*(np+1)+1)=ppy(k)
     z(numc(k)*(np+1)+1)=ppz(k)
       do i=1,np
        x(numc(k)*(np+1)+i+1)=ppx(k)+ R(k)*cos(pi*real(i)/180)*ap(k) + R(k)*sin(pi*real(i)/180)*a2p(k) 
        y(numc(k)*(np+1)+i+1)=ppy(k)+ R(k)*cos(pi*real(i)/180)*bp(k) + R(k)*sin(pi*real(i)/180)*b2p(k) 
        z(numc(k)*(np+1)+i+1)=ppz(k)+ R(k)*cos(pi*real(i)/180)*cp(k) + R(k)*sin(pi*real(i)/180)*c2p(k)  
       end do 
   end do
       
  
    do i=1,total*(np+1)
    write (tp+10000,*) x(i),y(i),z(i),1,1,1,1,1,1,1,1
    end do
  
  
    num=0 
    do k=1,total
      do i=1,np-1
      write (tp+10000,*) numc(k)*(np+1)+1,numc(k)*(np+1)+i+1,numc(k)*(np+1)+i+2 !! other elements!!
      num=num+1
      end do
      write(tp+10000,*)  numc(k)*(np+1)+1,numc(k)*(np+1)+np+1,numc(k)*(np+1)+2  !! one element !!
      num=num+1
   end do 
    
   
   
    do k=1,total-1,2 
      do i=2,np
        write(tp+10000,*) i    +numc(k)*(np+1),i+numc(k+1)*(np+1),i+1+numc(k  )*(np+1)
        write(tp+10000,*) i+ 1 +numc(k)*(np+1),i+numc(k+1)*(np+1),i+1+numc(k+1)*(np+1)
        num=num+2
      end do 
       write(tp+10000,*)  np+1 +numc(k)*(np+1),(np+1)+numc(k+1)*(np+1),2+numc(k  )*(np+1)
       write(tp+10000,*)  2    +numc(k)*(np+1),(np+1)+numc(k+1)*(np+1),2+numc(k+1)*(np+1)
       num=num+2
    end do 
   
   call flush (tp+10000)
   
    
    
    return  
    end subroutine 




end module 

