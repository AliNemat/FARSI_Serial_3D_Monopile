 Module M_Math
real(8)   :: pi=3.14159265359

contains 



FUNCTION CROSS(a,b)
   IMPLICIT  NONE
   Real(8),dimension (1:3)                     ::  CROSS
   real(8),dimension (1:3)                     ::   a,b
              
    CROSS(1)= a(2)*b(3)-a(3)*b(2)
    CROSS(2)=-a(1)*b(3)+a(3)*b(1)
    CROSS(3)= a(1)*b(2)-a(2)*b(1)
                           
              
 END FUNCTION CROSS  



FUNCTION MTRXSOL(FMtrx,KMtrx)
   IMPLICIT  NONE
   Real(8),dimension (3)                       ::  MtrxSol,FMtrx
   real(8),dimension (3,3)                     ::  Im,KMtrx
   real(8)                                     ::  DetMat

!! This is for a Symmertic 3 by 3 matrix !!

DetMat=  KMtrx(1,1)*KMtrx(2,2)*KMtrx(3,3) &
   &    +KMtrx(1,2)*KMtrx(2,3)*KMtrx(3,1) &
   &    +KMtrx(1,3)*KMtrx(2,1)*KMtrx(3,2) &
   &    -KMtrx(1,3)*KMtrx(2,2)*KMtrx(3,1) &
   &    -KMtrx(1,1)*KMtrx(2,3)*KMtrx(3,2) &
   &    -KMtrx(1,2)*KMtrx(2,1)*KMtrx(3,3)  


Im(1,1)= (KMtrx(2,2)*KMtrx(3,3)-KMtrx(2,3)*KMtrx(3,2))/DetMat
Im(1,2)=-(KMtrx(2,1)*KMtrx(3,3)-KMtrx(2,3)*KMtrx(3,1))/DetMat
Im(1,3)= (KMtrx(2,1)*KMtrx(3,2)-KMtrx(2,2)*KMtrx(3,1))/DetMat
Im(2,2)= (KMtrx(1,1)*KMtrx(3,3)-KMtrx(1,3)*KMtrx(3,1))/DetMat
Im(2,3)=-(KMtrx(1,1)*KMtrx(3,2)-KMtrx(1,2)*KMtrx(3,1))/DetMat
Im(3,3)= (KMtrx(1,1)*KMtrx(2,2)-KMtrx(1,2)*KMtrx(2,1))/DetMat

Im(2,1)=Im(1,2)          
Im(3,1)=Im(1,3)         
Im(3,2)=Im(2,3)         

MtrxSol(1)= Im(1,1)*FMtrx(1) +Im(1,2)*FMtrx(2) +Im(1,3)*FMtrx(3) 
MtrxSol(2)= Im(2,1)*FMtrx(1) +Im(2,2)*FMtrx(2) +Im(2,3)*FMtrx(3) 
MtrxSol(3)= Im(3,1)*FMtrx(1) +Im(3,2)*FMtrx(2) +Im(3,3)*FMtrx(3) 





              
                              
              
 END FUNCTION MTRXSOL  

real(8) function ROTATIONX(x1,y1,z1,u1,v1,w1,a1,b1,c1,teta1)
        implicit none 
                 real(8),intent(in) :: x1,y1,z1,u1,v1,w1,teta1,a1,b1,c1
              
                 ROTATIONX=(a1*(v1**2+w1**2)-u1*(b1*v1+c1*w1-u1*x1-v1*y1-w1*z1))*(1-cos(teta1*pi/180)) +&
                 &         x1*cos(teta1*pi/180)         +&
                 &         (-c1*v1+b1*w1-w1*y1+v1*z1)*sin(teta1*pi/180)
                
        return 
        end function ROTATIONX
        
real(8) function ROTATIONY(x1,y1,z1,u1,v1,w1,a1,b1,c1,teta1)
        implicit none 
                 real(8),intent(in) :: x1,y1,z1,u1,v1,w1,teta1,a1,b1,c1
                 
                 ROTATIONY=(b1*(u1**2+w1**2)-v1*(a1*u1+c1*w1-u1*x1-v1*y1-w1*z1))*(1-cos(teta1*pi/180)) +&
                 &         y1*cos(teta1*pi/180)         +&
                 &         (c1*u1-a1*w1+w1*x1-u1*z1)*sin(teta1*pi/180)
         return         
       end function ROTATIONY
       
real(8) function ROTATIONZ(x1,y1,z1,u1,v1,w1,a1,b1,c1,teta1)
        implicit none 
              real(8),intent(in) :: x1,y1,z1,u1,v1,w1,teta1,a1,b1,c1
              
 ROTATIONZ=(c1*((u1**2)+(v1**2))-w1*(a1*u1+b1*v1-u1*x1-v1*y1-w1*z1))*(1-cos(teta1*pi/180)) +&
                 &         z1*cos(teta1*pi/180)         +&
                 &         (-b1*u1+a1*v1-v1*x1+u1*y1)*sin(teta1*pi/180)
                           
             return  
        end function ROTATIONZ       

real(8) function LENGTHF(a1,b1,c1)
          implicit none 
          real(8),intent(in) :: a1,b1,c1
            LENGTHF=sqrt(a1**2+b1**2+c1**2)
          return 
         end function 

DOUBLE PRECISION FUNCTION HFF(a)

double precision :: a 
  if (a.gt.0.0) then 
    HFF=1.0
  else 
    HFF=0 
  end if 

end function HFF



SUBROUTINE M33INV (A, AINV, OK_FLAG)


!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************


      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV


    subroutine RotationInv_YCoord(theta,I_m,I_mRot)

    !****************************************************************************
    ! RotationInv_y-             Compute the second moment of inerta of the rotated body.


    ! Ttheta  = scalar input  angle in radian that the body is rotated
    ! I_m     = 3*3    input  mass moment of inertia of the unrotated body- fixed coordinate system 
    ! I_mRot  = 3*3    output mass moment of inertia of the rotated body - fixed coordinate system
    !****************************************************************************
 
      implicit none
      double precision,dimension(3,3),intent(in)      :: I_m
      double precision               ,intent(in)      :: Theta
      double precision,dimension(3,3),intent(out)     :: I_mRot
    
      double precision,dimension(3,3)                 :: R_InvRot

      !! Invese of the standard roation matrix   
      R_InvRot(1,1)=cos(theta)      
      R_InvRot(1,2)=0      
      R_InvRot(1,3)=-sin(theta)      
      R_InvRot(2,1)=0    
      R_InvRot(2,2)=1.0      
      R_InvRot(2,3)=0    
      R_InvRot(3,1)=sin(theta)      
      R_InvRot(3,2)=0      
      R_InvRot(3,3)=cos(theta)      

      I_mRot=matmul( matmul(R_InvRot,I_m),transpose(R_InvRot) )
      return

    end subroutine 
end module 


