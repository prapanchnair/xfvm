subroutine edgenormals(nnodes, nedges, nodexy, edge2cell, enormal,elength )

  ! ncells - number of cells

  ! nedges - number of edges

  ! nodexy - the coordinates of the nodes

  ! enormal - output of nedges size

  ! ed2cell - edge cell connectivity 


  ! normals calculated point from right cell to left cell

  integer :: i, ncells, nedges

  integer, dimension(4,nedges) :: edge2cell

  real, dimension(2,nnodes) :: nodexy

  real, dimension(2, nedges) :: enormal

  real, dimension(nedges) :: elength


  do i=1, nedges

     enormal(1,i) = nodexy(2, edge2cell(1,i))- nodexy(2, edge2cell(2,i))

     enormal(2,i) = nodexy(1, edge2cell(2,i))- nodexy(1, edge2cell(1,i))

     elength(i) = sqrt(enormal(1,i)**2.0 + enormal(2,i)**2.0)

     enormal(1,i) = enormal(1,i)/elength(i)

     enormal(2,i) = enormal(2,i)/elength(i)

  end do

  !do i=1,nedges

  !  write(*,'(100G14.5)') enormal(1,i), enormal(2,i)

  !end do


end subroutine edgenormals
