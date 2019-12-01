PROGRAM main
IMPLICIT NONE
REAL*8::LE0,Lcm,dLMPi,dLMG,LMN,SE0,Scm,dSMPi,dSMG,SMN,XE0,Xcm,dXMPi,dXMG,XMN
INTEGER::K
open(unit=11,file='SET1.txt')
open(unit=12,file='SET2.txt')
open(unit=13,file='SET3.txt')


write(11,1)"  ","EB0","ECM","DMBPI","DMBG","MB"

K=1
call Corrections(K,LE0,Lcm,dLMPi,dLMG,LMN,SE0,Scm,dSMPi,dSMG,SMN,XE0,Xcm,dXMPi,dXMG,XMN)
write(11,2) LE0,Lcm,dLMPi,dLMG,LMN
write(11,3) SE0,Scm,dSMPi,dSMG,SMN
write(11,4) XE0,Xcm,dXMPi,dXMG,XMN


write(12,1)"  ","EB0","ECM","DMBPI","DMBG","MB"

K=2
call Corrections(K,LE0,Lcm,dLMPi,dLMG,LMN,SE0,Scm,dSMPi,dSMG,SMN,XE0,Xcm,dXMPi,dXMG,XMN)
write(12,2) LE0,Lcm,dLMPi,dLMG,LMN
write(12,3) SE0,Scm,dSMPi,dSMG,SMN
write(12,4) XE0,Xcm,dXMPi,dXMG,XMN



write(13,1)"  ","EB0","ECM","DMBPI","DMBG","MB"

K=3
call Corrections(K,LE0,Lcm,dLMPi,dLMG,LMN,SE0,Scm,dSMPi,dSMG,SMN,XE0,Xcm,dXMPi,dXMG,XMN)
write(13,2) LE0,Lcm,dLMPi,dLMG,LMN
write(13,3) SE0,Scm,dSMPi,dSMG,SMN
write(13,4) XE0,Xcm,dXMPi,dXMG,XMN



1 format(/,A1,5X,A2,5X,A3,5X,A4,5X,A5,5X,A6,/)
2 format("L",5x,f8.3, 5x,f7.3, 5x,f7.3, 5x,f7.3, 5x,f8.3)
3 format("S",5x,f8.3, 5x,f7.3, 5x,f7.3, 5x,f7.3, 5x,f8.3)
4 format("X",5x,f8.3, 5x,f7.3, 5x,f7.3, 5x,f7.3, 5x,f8.3)


end 

!-------------------------------------------------------------------


SUBROUTINE Corrections(K,LE0,Lcm,dLMPi,dLMG,LMN,SE0,Scm,dSMPi,dSMG,SMN,XE0,Xcm,dXMPi,dXMG,XMN)
    IMPLICIT NONE
    INTEGER::K
    REAL*8::set(6), set1(6), set2(6), set3(6)
    REAL*8::dM,dE,PI, hc, F_pi, piM, Gs
    REAL*8:: U,C,UA, CA,CV, UMP, UEP, UM, UE, UW, UR, CMP, CEP, CM, CE, CW, CR
    INTEGER::MAX=200,i    !二分法
    INTEGER::n            !pi介子修正
    REAL*8::tol=1.0D-10   !二分法
    REAL*8::a,b,x,DUU,DCC !二分法
    REAL*8::aa,ab         !pi介子修正  
    REAL*8::RUU,RUC,RCC,ALU,ALC, AL !胶子修正
    
  
    REAL*8::FNNPI2
    REAL*8::LE0,Lcm,dLMPi,dLMG,LMN
    REAL*8::SE0,Scm,dSMPi,dSMG,SMN
    REAL*8::XE0,Xcm,dXMPi,dXMG,XMN
    REAL*8,EXTERNAL::U1,C1,F
    COMMON /Meson_Mean_Field/ dM,dE                                              
    PARAMETER (PI = 3.141592653589793d0, hc = 197.3269631D0)
    parameter (F_pi = 93.0, piM = 140.0, AL=0.58D0)
    !          g_sigma
    PARAMETER (Gs = 0.0)
    
  !          Mass    a[MeV]^3        V_0 
    DATA set1/ 250,  0.579450D0,   -24.286601D0,    &
               1300, 0.118172D0,    284.58724D0/ 
    
    
    DATA set2/  300,  0.534296D0,   -62.257187D0,     &
                1350,  0.117312D0,    239.53994D0/
    
    
    DATA set3/  350,  0.495596D0,   -102.041575D0,    &
                1400,  0.116036D0,    193.67265D0/

    ! External sigma meson field, correction on mass
    dM = 0.0
    ! External rho & omega meson field correction on energy
    dE = 0.0
    
  
    
   
    SELECTCASE (K)
    CASE(1)
        set = set1
    CASE(2)
        set = set2
    CASE(3)
        set = set3
    END SELECT
   
    

!--------------------------------------        
UA=set(2) 
UMP = set(1) + set(3)/2.0D0
   
a=0d0
b=1000d0

DO i=1,MAX !二分法求根
x=(a+b)/2D0 
  IF (U1(x, UA, UMP)==0) EXIT
    IF(U1(a, UA, UMP)*U1(x, UA, UMP)<0d0) THEN
    b=x
    ELSE
    a=x
    END IF
  DUU=abs(U1(a, UA, UMP)-U1(b, UA, UMP))
  IF (dUU<tol) EXIT
end do
U=x !U夸克质量

    UEP = U
    UW = UEP + UMP
    UE =UEP + set(3)/2.0
    UM = set(1)
    UR = (set(2)*hc**3*UW)**(-1.0/4.0)
    
!--------------------------------------     
     
CA = set(5)
CMP = set(4) + set(6)/2.0D0
    
a=1000d0
b=3000d0

DO i=1,MAX !二分法求根
x=(a+b)/2D0 
  IF (C1(x, CA, CMP)==0) EXIT
    IF(C1(a, CA, CMP)*C1(x, CA, CMP)<0d0) THEN
    b=x
    ELSE
    a=x
    END IF
  DCC=abs(C1(a, CA, CMP)-C1(b, CA, CMP))
  IF (DCC<tol) EXIT
end do
C=x    !C夸克质量

    CEP = C
    CW = CEP + CMP
    CE =CEP + set(6)/2.0
    CM =set(4)
    CR = (CA*hc**3*CW)**(-1.0/4.0)
    
!--------------------------------------     
 
    LE0 = 2*UE+CE  !lambda超子零阶能量
    SE0 = LE0  !sigma超子零阶能量
    XE0 = 2*CE+UE  !xi超子零阶能量
   
 
!--------------------------------------

     !======================!
     !   C.M. correction    !
     !======================! 

 Lcm=1/(2*UM+CM)*(12*UM/((3*UEP+UMP)*UR**2)+6*CM/((3*CEP+CMP)*CR**2))&
    +3/(2*UM+CM)*(4*UA*hc**3*UM*(UEP+UMP)*UR**2/(3*UEP+UMP)+2*CA*hc**3*CM*(CEP+CMP)*CR**2/(3*CEP+CMP))&
    -3/(2*(2*UM+CM)**2)*(UA*hc**3*UM**2*(11*UEP+UMP)*UR**2/(3*UEP+UMP)+CA*hc**3*CM**2*(11*CEP+CMP)*CR**2/(2*(3*CEP+CMP)))&
    -1/(2*(2*UM+CM)**2)*(UA*hc**3*UM**2*(UEP+11*UMP)*UR**2/(3*UEP+UMP)+CA*hc**3*CM**2*(CEP+11*CMP)*CR**2/(2*(3*CEP+CMP)))&
    -1/(2*(2*UM+CM)**2)*(UA*hc**3*UM**2*(UEP+3*UMP)*(11*UEP+UMP)*UR**2/(3*UEP+UMP)**2+CA*hc**3*CM**2*(UEP+3*UMP)*(11*CEP+CMP)*CR**2/((3*UEP+UMP)*(3*CEP+CMP))&
    +UA*hc**3*UM**2*(CEP+3*CMP)*(11*UEP+UMP)*UR**2/((3*CEP+CMP)*(3*UEP+UMP)))
 Scm= Lcm
 Xcm= 1/(2*CM+UM)*(12*CM/((3*CEP+CMP)*CR**2)+6*UM/((3*UEP+UMP)*UR**2))&
    +3/(2*CM+UM)*(4*CA*hc**3*CM*(CEP+CMP)*CR**2/(3*CEP+CMP)+2*UA*hc**3*UM*(UEP+UMP)*UR**2/(3*UEP+UMP))&
    -3/(2*(2*CM+UM)**2)*(CA*hc**3*CM**2*(11*CEP+CMP)*CR**2/(3*CEP+CMP)+UA*hc**3*UM**2*(11*UEP+UMP)*UR**2/(2*(3*UEP+UMP)))&
    -1/(2*(2*CM+UM)**2)*(CA*hc**3*CM**2*(CEP+11*CMP)*CR**2/(3*CEP+CMP)+UA*hc**3*UM**2*(UEP+11*UMP)*UR**2/(2*(3*UEP+UMP)))&
    -1/(2*(2*CM+UM)**2)*(CA*hc**3*CM**2*(CEP+3*CMP)*(11*CEP+CMP)*CR**2/(3*CEP+CMP)**2+UA*hc**3*UM**2*(CEP+3*CMP)*(11*UEP+UMP)*UR**2/((3*CEP+CMP)*(3*UEP+UMP))&
    +CA*hc**3*CM**2*(UEP+3*UMP)*(11*CEP+CMP)*CR**2/((3*UEP+UMP)*(3*CEP+CMP)))

 
     !======================!
     ! Pi meson correction  !
     !======================!
     

aa=0d0
ab=2000d0
n=100000

call solve(F,C,aa,ab,n,UW,UEP,UMP,UR)
FNNPI2=(25*piM**2)/(1296D0*pi*F_pi**2)*(5*UEP+7*UMP)**2/(3D0*UEP+UMP)**2

dLMPi=(-108.0/25.0D0)*C*FNNPI2
dSMPi=(-12.0/5.0D0)*C*FNNPI2
dXMPi=(-27.0/25.0D0)*C*FNNPI2

 
 
     !======================!
     ! Gluon correction     !
     !======================!
 
 
RUU=sqrt(6/(UEP**2-UMP**2))
RUC=sqrt(3*(1/(UEP**2-UMP**2)+1/(CEP**2-CMP**2)))
RCC=sqrt(6/(CEP**2-CMP**2))
ALU=1/((UEP+UMP)*(3*UEP+UMP))
ALC=1/((CEP+CMP)*(3*CEP+CMP))
dLMG=AL*(-3)*256/(9*sqrt(pi))*1/RUU**3*1/(3*UEP+UMP)**2& 
   +AL*16/(3*sqrt(pi)*RUU)*(1-2*ALU/RUU**2+3*ALU**2/RUU**4)&
   +AL*(-2)*16/(3*sqrt(pi)*RUC)*(1-(ALU+ALC)/RUC**2+3*ALU*ALC/RUC**4)&
   +AL*16/(3*sqrt(pi)*RCC)*(1-2*ALC/RCC**2+3*ALC**2/RCC**4)
   
   
dSMG=AL*256/(9*sqrt(pi))*1/RUU**3*1/(3*UEP+UMP)**2&
    +AL*(-4)*256/(9*sqrt(pi))*1/RUC**3*1/((3*UEP+UMP)*(3*CEP+CMP))&
    +AL*16/(3*sqrt(pi)*RUU)*(1-2*ALU/RUU**2+3*ALU**2/RUU**4)&
    +AL*(-2)*16/(3*sqrt(pi)*RUC)*(1-(ALU+ALC)/RUC**2+3*ALU*ALC/RUC**4)&
    +AL*16/(3*sqrt(pi)*RCC)*(1-2*ALC/RCC**2+3*ALC**2/RCC**4)
    
dXMG=AL*256/(9*sqrt(pi))*1/RCC**3*1/(3*CEP+CMP)**2&
    +AL*(-4)*256/(9*sqrt(pi))*1/RUC**3*1/((3*UEP+UMP)*(3*CEP+CMP))&
    +AL*16/(3*sqrt(pi)*RUU)*(1-2*ALU/RUU**2+3*ALU**2/RUU**4)&
    +AL*(-2)*16/(3*sqrt(pi)*RUC)*(1-(ALU+ALC)/RUC**2+3*ALU*ALC/RUC**4)&
    +AL*16/(3*sqrt(pi)*RCC)*(1-2*ALC/RCC**2+3*ALC**2/RCC**4)

 
 LMN=LE0-Lcm+dLMPi+dLMG
 SMN=SE0-Scm+dSMPi+dSMG
 XMN=XE0-Xcm+dXMPi+dXMG
 

 
 END  
     
   
   
!--------------------------------------
   
FUNCTION U1(U, UA, UMP)
IMPLICIT NONE
REAL*8::U1
REAL*8::U, UA, UMP 
REAL*8,PARAMETER::hc = 197.3269631d0
U1=(U-UMP)*sqrt((U+UMP)/(UA*hc**3))-3
END FUNCTION U1 
   
!--------------------------------------   
   
FUNCTION C1(C, CA, CMP)
IMPLICIT NONE
REAL*8::C1
REAL*8::C, CA, CMP
REAL*8,PARAMETER::hc = 197.3269631d0
C1=(C-CMP)*sqrt((C+CMP)/(CA*hc**3))-3
END FUNCTION C1    
   
!--------------------------------------    
  
FUNCTION F(k,UW,UEP,UMP,UR)
IMPLICIT NONE
REAL*8::F,k,UW,UEP,UMP,UR
REAL*8,PARAMETER::pi=3.141592653589793d0
REAL*8,PARAMETER::piM=140d0
F=(1/(pi*piM**2))*(k**4/(k**2+piM**2))*(1-(3*k**2)/(2*UW*(5*UEP+7*UMP)))**2*exp(-k**2*UR**2/2)
END FUNCTION F 

!-------------------------------------- 

SUBROUTINE solve(F,C,aa,ab,n,UW,UEP,UMP,UR)
IMPLICIT NONE
REAL*8::C,aa,ab,UW,UEP,UMP,UR
REAL*8::h,F1,F2,F3,F4,t1,t2
INTEGER::n,k
REAL*8,external::F
C=0d0
h=(ab-aa)/n/2d0

F1=F(aa,UW,UEP,UMP,UR)
F2=F(ab,UW,UEP,UMP,UR)

C=F1+F2
!k=0 情况
F1=F(aa+h,UW,UEP,UMP,UR)
C=C+4d0*F1

DO k=1,n-1
t1=aa+(2d0*k+1)*h
t2=aa+2d0*k*h
F3=F(t1,UW,UEP,UMP,UR)
F4=F(t2,UW,UEP,UMP,UR)
C=C+F3*4d0+F4*2d0  
END DO

C=C*h/3d0

END SUBROUTINE solve