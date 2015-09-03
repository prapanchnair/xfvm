subroutine areacalc(ncells,nnodes, nodexy,cell2node, A, area)

! subroutine to calculate the cell areas corresponding to cell num
! ncells- num of cells
!nnodes - num of nodes
! nodexy- xy coordinates of nodes
! A - area vector(output)
! area - total area of the domain(output)

integer :: ncells, nnodes, i 

integer, dimension(3,ncells):: cell2node
integer:: AAAA
real, dimension(2, nnodes):: nodexy

real, dimension(ncells)::A,A1

real :: area, dx21,dy21,dx31,dy31

area = 0.0
do i=1,ncells

  dx21 = nodexy(1,cell2node(2,i)) - nodexy(1,cell2node(1,i));
  dy21 = nodexy(2,cell2node(2,i)) -nodexy(2,cell2node(1,i));
  
  dx31 = nodexy(1,cell2node(3,i)) -nodexy(1,cell2node(1,i));
  dy31 = nodexy(2,cell2node(3,i)) -nodexy(2,cell2node(1,i));

  A1(i) =( 0.5*(dx21*dy31 - dy21*dx31));

  A(i)=abs(A1(i))
area = area+A1(i)
!write(*,*) A(i)


end do

write(*,*) "total area =", area
!pause
read(*,*) AAAA
end subroutine areacalc
