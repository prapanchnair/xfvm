program areainfo

implicit none

integer:: ndim,nzones,nbtypes,nnodes,nedges,ncells,maxppf,maxfpc
integer::i,j,ppf,p1,p2,n1,n2,n21,n22,testp1,testp2
real,allocatable,dimension(:,:)::nodexy, enormal
real,allocatable,dimension(:)::A
real:: dx21,dx31,dy21,dy31,area
integer, allocatable,dimension(:,:)::edge2cell, cell2node


open(11,file='squaresimple-COBALT.inp',form='formatted')

read(11,*) ndim , nzones, nbtypes

read(11,*) nnodes, nedges, ncells, maxppf, maxfpc

allocate(nodexy(2,nnodes))

allocate(edge2cell(4,nedges))

allocate(cell2node(3,ncells))

allocate(A(ncells))

allocate(enormal(2,nedges))

do i=1,nnodes
   
   read(11,*) nodexy(1,i),nodexy(2,i)

end do

do i=1, nedges
 
      read(11,*) ppf, edge2cell(1,i),edge2cell(2,i),&
           
           edge2cell(3,i), edge2cell(4,i)

end do


call edgearrange(nedges,edge2cell)
do i=1,nedges

   write(*,'(5I7)') (edge2cell(j,i),j=1,4)

end do

write(*,*) "=========================="

! Now write the cell number info in an array called cell2node
! go thro each column of edge2cell nedge times.. 
! the third column would correspond to the cell number-cell
! go to the cell'th column of cell2 node fill first 2 columns with the node numbers..
! use modulo operation
! check 1st vlue of the pair with that in cell2node, where it matches, put first! value there and the second value in the next colum - by a modulo 

!befor i call cell node, i have to arrange edge2cell
!so that all the boundary edges are on the top 
!of the list


call cellnode(ncells, nedges, edge2cell, cell2node)

!calculating areas of the cells.
! we have cell2node and nodexy
!=========================================
call areacalc(ncells, nnodes, nodexy, cell2node, A , area)
!===========================================
call edgenormals(nnodes, nedges, nodexy, edge2cell, enormal)
!===========================================

















close(11)

deallocate(nodexy,edge2cell)

end program





