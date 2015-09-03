subroutine timestep(ncells,gamma,U,A,CFL,dt)
  
  integer :: ncells
  
  real, dimension(4,ncells):: U
  
  real, dimension(ncells):: A
  
  real :: CFL
  
  real :: dt
  real, dimension(ncells) :: ddt
  
  real :: mint, pi

  real :: uu, vv, p, r, gamma,c,qq,L
  
  pi= 4.0*atan(1.0)
  
! write(*,*) "pi", pi
  do i=1,ncells
!write(*,*) A(i)     
     uu = U(2,i)/U(1,i)
     vv = U(3,i)/U(1,i)
     
     r  = U(1,i)
     
     qq = uu**2 + vv**2
     
     p  = (gamma -1.0)*(U(4,i) - (0.5*r*qq))

     c  = sqrt(gamma*p/r)
     
     L  = sqrt(A(i)/pi)
!write(*,*) "A(i) and pi",sqrt( A(i)/ pi)
 !write(*,*) "L",L    
     ddt(i) =  L / (sqrt(qq) + c)
     
!write(*,*) "ddddt" , ddt(i)     
  end do
  
  ! now find the min of dt
  dt = minval(ddt)
 ! write(*,*) "dt",dt
  dt=CFL*dt
 ! write(*,*) "after cfl dt" , dt
end subroutine timestep
