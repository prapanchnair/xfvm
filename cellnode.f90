subroutine cellnode(ncells, nedges, edge2cell, cell2node)

! this is a subroutine that generates the array of 
! cell number to the 3 node numbers that make the triangular cell

! ncells - total number of cells

! nedges - total number of edges

! edge2cell- the arrage containing 4 rows
!            1st and 2nd are the starting and ending nodes

!            3 and 4th are the cell numbers to left and right
!            of the corresponding edge when moving from 1 to 2

! cell 2 node had 3 rows each representing the nodes making the triangle
! the nodes are in anticlockwise direction

! declarations

integer:: ncells, nedges

integer :: i, j, p1, p2, n1, n2, n21, n22, testp1, testp2

integer, dimension(4, nedges) :: edge2cell

integer, dimension(3, ncells) :: cell2node

cell2node=0

do i=1,nedges

   p1= edge2cell(3,i)
   
   p2= edge2cell(4,i)

   n1= edge2cell(1,i)

   n2= edge2cell(2,i)
testp1=0
testp2=0

   do j=1,3
 
      if (testp1==0) then
         if (n1==cell2node(j,p1).or. cell2node(j,p1)==0) then

            cell2node(j,p1)=n1
            cell2node(mod(j,3)+1,p1)=n2
            testp1=1


         else if(n1==cell2node(mod(j,3)+1,p1)) then

            cell2node(mod(j,3)+1,p1)=n1
            cell2node(mod(j,3)+2,p1)=n2
            testp1=1
         else 
            cell2node(mod(j,3)+2,p1)=n1
            cell2node((mod((mod(j,3)+2),3)+1),p1)=n2
            
            testp1=1
         endif
      end if
      !if (i==3) write(*,*)cell2node(1,3)
   end do



   if (p2>0) then
      
      n22=n1
      n21=n2


      do j=1,3
         
         if (testp2==0) then
            if (n21==cell2node(j,p2) .or. cell2node(j,p2)==0) then
               
               cell2node(j,p2)=n21
               cell2node(mod(j,3)+1,p2)=n22
               testp2 =1
               
            else if(n21==cell2node(mod(j,3)+1,p2))  then
               
               cell2node(mod(j,3)+1,p2)=n21
               cell2node(mod(j,3)+2,p2)=n22
               testp2 =1
               
            else
               
               cell2node(mod(j,3)+2,p2)=n21
               cell2node((mod(mod(j,3)+2,3)+1),p2)=n22
               testp2=1
               
            end if
         end if
      end do
      
   endif
   !write(*,'(5I3)') (cell2node(j,3),j=1,3)
   
end do

!do i=1,ncells
   
!      write(*,'(5I7)') (cell2node(j,i),j=1,3)
   
!end do


end subroutine cellnode
