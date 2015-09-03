subroutine edgearrange(nedges,edge2cell,nbedges)

  integer :: nedges,i,j,k,nbedges

  integer,dimension(4,nedges)::edge2cell,dummy

  integer,dimension(4,1)::swap

  do j=1,nedges
     do i=1,nedges

        if(edge2cell(4,i)>=0 .AND. edge2cell(4,i+1)<0) then

           !selection sort here
           do k=1,4
              swap(k,1)= edge2cell(k,i)
           end do
           do k=1,4
              edge2cell(k,i)=edge2cell(k,i+1)
           end do
           do k=1,4
              edge2cell(k,i+1)=swap(k,1)
           end do

        endif

     end do

  end do

  !Find the number of boundary edges

  nbedges=0

  do i=1,nedges

     if(edge2cell(4,i)<0) then
        
        nbedges=nbedges+1
     
     end if

  end do

end subroutine edgearrange
