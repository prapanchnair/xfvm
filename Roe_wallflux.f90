subroutine Roe_wallflux(UL,enormal,gamma,F_interface)

  implicit none


  integer::i
  real, dimension(4):: UL,FL,FR,F_interface
  real :: gamma
  real :: rL,rR,pL,pR,u_L,u_R,v_L,v_R,eL,eR,HL,HR,qL,qR,Vn_L,Vn_R
  
  real,dimension(2) :: enormal
  

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

  Vn_L= (u_L*enormal(1)) + (v_L*enormal(2))

  HL = (eL+pL)/rL

  !write(*,*) Vn_L  
  !---------------------

  !write(*,*) Vn_R
  !---------------------




  F_interface(1) = 0.0
  F_interface(2) = PL*enormal(1)
  F_interface(3) = PL*enormal(2)
  F_interface(4) = 0.0




  !write(*,*) F_interface

end subroutine Roe_wallflux

