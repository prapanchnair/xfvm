subroutine Roe_interior(UL,UR,enormal,gamma,F_interface)
  
  implicit none
  
  
  integer::i
  real, dimension(4):: UL,UR,FL,FR,F_interface
  real :: r0,p0,c0,u0,v0,gamma,H0,q0,Vn0
  real :: rL,rR,pL,pR,u_L,u_R,v_L,v_R,eL,eR,HL,HR,qL,qR,Vn_L,Vn_R
  real :: delp,delVn,delu,delv,delr
  real,dimension(2) :: enormal
  real,dimension(4) :: F1,F23,F4,Term
  
  !  UL(1)=1.0; UR(1)=0.125
  !  UL(2)=0.75; UR(2)=0.0
  !  UL(3)=0.0; UR(3)=0.0
  !  UL(4)=2.78; UR(4)=0.285
  
  !  gamma =1.4
  !  enormal(1)=0.707;enormal(2)=0.707
  
  !find FL and FR
  !===============
  
  rL=UL(1)
  u_L = UL(2)/UL(1)
  v_L = UL(3)/UL(1)
  eL = UL(4)
  qL = (u_L**2 + v_L**2)**0.5
  
  PL = (gamma-1)*(UL(4)-(0.5*rL*qL**2))
!write(*,*) "PL",PL
  Vn_L= (u_L*enormal(1)) + (v_L*enormal(2))

  HL = (eL+pL)/rL

  !write(*,*) Vn_L  
  !---------------------

  rR=UR(1)
  u_R = UR(2)/UR(1)
  v_R = UR(3)/UR(1)
  eR = UR(4) 
  qR = (u_R**2 + v_R**2)**0.5

  PR = (gamma-1)*(UR(4)-(0.5*rR*qR**2))
  
  Vn_R = (u_R*enormal(1)) + (v_R*enormal(2))

  HR = (eR+pR)/rR

  !write(*,*) Vn_R
  !---------------------

  !FL
  !---
  
  FL(1) = Vn_L*rL
  FL(2) = Vn_L*UL(2) + (PL*enormal(1))
  FL(3) = Vn_L*UL(3) + (PL*enormal(2))
  FL(4) = Vn_L*(eL+PL)
  
  !FR
  !----
  
  FR(1) = Vn_R*rR
  FR(2) = Vn_R*UR(2) + (PR*enormal(1))
  FR(3) = Vn_R*UR(3) + (PR*enormal(2))
  FR(4) = Vn_R*(eR+PR)
  
  !write(*,*) FL
  !write(*,*) FR

  !find A(Qr-Ql)=Term
  !------------
  !Averaged quantities
  !-------------------
  
  r0 = sqrt(rL*rR)

  u0 = ((sqrt(rL)*u_L)+(sqrt(rR)*u_R))/(sqrt(rR)+sqrt(rL))

  v0 = ((sqrt(rL)*v_L)+(sqrt(rR)*v_R))/(sqrt(rR)+sqrt(rL))
  
  h0 = ((sqrt(rL)*HL)+(sqrt(rR)*HR))/(sqrt(rR)+sqrt(rL))
  
  q0 = sqrt(u0**2 + v0**2)

  c0 = ((gamma-1)*(h0 - (0.5*q0**2)))**0.5
    
  Vn0 = u0*enormal(1) + v0*enormal(2)
  !---------------------------------

  delp= pR-pL;delVn=Vn_R - Vn_L
  delr=rR-rL;delu=u_R-u_L
  delv=v_R-v_L

  !F1 
  !----

  F1(1) =abs(Vn0 - c0)*(0.5/c0**2)*(delp - (r0*c0*delVn))*(1.0)
  F1(2) =abs(Vn0 - c0)*(0.5/c0**2)*(delp - (r0*c0*delVn))*(u0- (c0*enormal(1)))
  F1(3) =abs(Vn0 - c0)*(0.5/c0**2)*(delp - (r0*c0*delVn))*(v0- (c0*enormal(2)))
  F1(4) =abs(Vn0 - c0)*(0.5/c0**2)*(delp - (r0*c0*delVn))*(h0- (c0*Vn0))



!F23
  !----

  F23(1) =abs(Vn0)*( (delr-(delp/(c0**2)))*(1.0) + 0 )
  F23(2) =abs(Vn0)*( ((delr-(delp/(c0**2)))*(u0)) + (r0*(delu - (enormal(1)*delVn))) )
  F23(3) =abs(Vn0)*( ((delr-(delp/(c0**2)))*(v0)) + (r0*(delv - (enormal(2)*delVn))) )
  F23(4) =abs(Vn0)*( ((delr-(delp/(c0**2)))*(0.5*q0**2)) + (r0*((u0*delu)+(v0*delv) - (Vn0*delVn))) )



  !F4
  !----

  F4(1) =abs(Vn0 + c0)*(0.5/c0**2)*(delp + (r0*c0*delVn))*(1.0)
  F4(2) =abs(Vn0 + c0)*(0.5/c0**2)*(delp + (r0*c0*delVn))*(u0+ (c0*enormal(1)))
  F4(3) =abs(Vn0 + c0)*(0.5/c0**2)*(delp + (r0*c0*delVn))*(v0+ (c0*enormal(2)))
  F4(4) =abs(Vn0 + c0)*(0.5/c0**2)*(delp + (r0*c0*delVn))*(h0+ (c0*Vn0))


  !Term
  !----

  Term(1) = (F1(1)) + (F23(1)) + (F4(1))
  Term(2) = (F1(2)) + (F23(2)) + (F4(2))
  Term(3) = (F1(3)) + (F23(3)) + (F4(3))
  Term(4) = (F1(4)) + (F23(4)) + (F4(4))
  
  do i=1,4
     
     F_interface(i) = 0.5*(FL(i)+FR(i)-Term(i))   
     
  end do


!write(*,*) " u is negative",u_R
!pause

!write(*,*) F_interface
!write(*,*) "c0" , c0


end subroutine Roe_interior

