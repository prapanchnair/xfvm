program euler2D

  implicit none

  !Declarations
  !=======================

  !Geometric Variables, filenames,time
  integer:: ndim,nzones,nbtypes,nnodes,nedges,ncells,maxppf,maxfpc
  integer::i,j,ppf,p1,p2,n1,n2,n21,n22,testp1,testp2,nbedges,maxiter,iter
  integer:: iL,iR
  real,allocatable,dimension(:,:)::nodexy, enormal
  real,allocatable,dimension(:)::A,elength,density,U_vel,V_vel,pressure
  real:: dx21,dx31,dy21,dy31,area,dtmin,Rtol,maxR,dt
  integer, allocatable,dimension(:,:)::edge2cell, cell2node
  character*30 :: gridfile
  integer:: AAAA,ios
  !Reference flow variables
  real::CFL,P_ref,den_ref, u_ref,v_ref,gamma,q_ref,Mx_ref,My_ref,c_ref,M_ref

  !Boundary variables
  real, dimension(4):: U_inf, U_init
  real ::  M_infx, M_infy, M_inf, M_init,add,rms

  !Other Flow Variables
  real,allocatable,dimension(:,:)::U,R,Unew,err
  real,dimension(4)::F_interface
  !Reading inputs
  !=======================

  open(9,file='fvminput',form='formatted')

  read(9,'(1X)')
  read(9,'(1X)')
  read(9,'(1X)')
  read(9,*) P_ref , den_ref, u_ref, v_ref
  read(9,'(1X)') 
  read(9,*) gridfile
  read(9,'(1X)') 
  read(9,*) CFL
  read(9,'(1X)')
  read(9,*) gamma
  write(*,*) "gamma" ,gamma
read(*,*) AAAA  
!pause

  !Reading grids
  !=======================

  open(11,file=gridfile,form='formatted')

  read(11,*) ndim , nzones, nbtypes
  read(11,*) nnodes, nedges, ncells, maxppf, maxfpc


  allocate(nodexy(2,nnodes))
  allocate(edge2cell(4,nedges))
  allocate(cell2node(3,ncells))
  allocate(A(ncells))
  allocate(enormal(2,nedges))
  !  allocate(dt(ncells))
  allocate(elength(nedges))
  do i=1,nnodes

     read(11,*) nodexy(1,i),nodexy(2,i)

  end do

  do i=1, nedges

     read(11,*) ppf, edge2cell(1,i),edge2cell(2,i),&
          
          edge2cell(3,i), edge2cell(4,i)

  end do

  !Arrange the edges with boundary edges on top
  !=========================================
  call edgearrange(nedges,edge2cell,nbedges)

  write(*,*) "=========================="
  write(*,*) "Edge to Cell info"

  do i=1,nedges

     write(*,'(5I7)') (edge2cell(j,i),j=1,4)

  end do




  write(*,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

  ! Now write the cell number info in an array called cell2node
  ! go thro each column of edge2cell nedge times.. 
  ! the third column would correspond to the cell number-cell
  ! go to the cell'th column of cell2 node fill first 2 columns with the node numbers..
  ! use modulo operation
  ! check 1st vlue of the pair with that in cell2node, where it matches, put first
  ! value there and the second value in the next colum - by a modulo 

  !befor i call cell node, i have to arrange edge2cell
  !so that all the boundary edges are on the top 
  !of the list


  call cellnode(ncells, nedges, edge2cell, cell2node)

  !calculating areas of the cells.
  ! we have cell2node and nodexy
  !-------------------------------------------
  call areacalc(ncells, nnodes, nodexy, cell2node, A , area)
  !-------------------------------------------
  call edgenormals(nnodes, nedges, nodexy, edge2cell, enormal,elength)
  !-------------------------------------------
  ! write(*,*) "Number of bedges =", nbedges



  ! have some reference values calculated
  !=====================================

  q_ref = sqrt(u_ref**2 + v_ref**2)
  c_ref = sqrt(gamma*P_ref/den_ref)
  write(*,*) "speed of sound", c_ref

  Mx_ref = u_ref/ c_ref

  My_ref = v_ref/ c_ref

  M_ref = q_ref/c_ref



  !Set the freestream values. _inf
  !=======================

  !  M_inf=2.0 !This is just temporary

  U_inf(1) = 1.0
  U_inf(2) = Mx_ref
  U_inf(3) = My_ref
  U_inf(4) = (1.0/(gamma*(gamma-1))) + (0.5*(M_ref**2))
  write(*,*) U_inf

  !Initialize the flow 
  !=======================

  allocate(U(4,ncells))
  allocate(Unew(4,ncells))
  allocate(R(4,ncells))
  allocate(err(4,ncells))

  U=1.0

  M_init = 0.0

  U_init(1) = 1.0
  U_init(2) = M_init
  U_init(3) = 0.0
  U_init(4) = (1.0/(gamma*(gamma-1))) + (0.5*(M_init**2))!put init instead of ref

  do i=1,ncells

     U(1,i) =  U_init(1)*U(1,i)
     U(2,i) =  U_init(2)*U(2,i)
     U(3,i) =  U_init(3)*U(3,i)
     U(4,i) =  U_init(4)*U(4,i)

  end do
  !write(*,*)" U "
  !write(*,*) U
  !pause
read(*,*) AAAA  
  call timestep(ncells,gamma,U,A,CFL,dt)

  ! write(*,*) "dt", dt


  !Start the time loop 
  !=======================

  do iter= 1,15000!while(iter<maxiter .and. maxR>Rtol)


     R = 0.0


     !For interior edges, calculate fluxes
     ! and hence find the residuals for
     !cells associated with these edges

     do i=nbedges+1,nedges

        iL = edge2cell(3,i)
        iR = edge2cell(4,i)
        F_interface = 0.0

        call Roe_interior(U(1:4,iL),U(1:4,iR),enormal(1:2,i),gamma,F_interface)

        R(1:4,iL)=R(1:4,iL) + F_interface(1:4)*elength(i)
        R(1:4,iR)=R(1:4,iR) - F_interface(1:4)*elength(i)
        !      write(*,*) "finterface", F_interface
     end do
     !Start here

     !For boundary edges calculate the fluxes
     !calssify the boundaries here. 
     ! There should be 2 types of boundaries.
     !boundary - wall has -4 as a flag

     do i=1,nbedges

        if( edge2cell(4,i) == -4 ) then !wall

           iL =  edge2cell(3,i)

           call Roe_wallflux(U(1:4,iL),enormal(1:2,i),gamma,F_interface)

        else  !farfield

           iL = edge2cell(3,i)

           call Roe_interior(U(1:4,iL),U_inf(1:4),enormal(1:2,i),gamma,F_interface)


        end if
        !       write(*,*) "fboundary", f_interface


        do j=1,4    
           R(j,iL)=R(j,iL) + (F_interface(j)*elength(i))
        end do

        !write(*,*) "R"
        !write(*,*) R
        !pause
     end do

     !Timestep

     !Forward Euler Step
     !=================

     do i=1,ncells

        Unew(1,i)= U(1,i) - (dt*R(1,i)/A(i))

        Unew(2,i)= U(2,i) - (dt*R(2,i)/A(i))

        Unew(3,i)= U(3,i) - (dt*R(3,i)/A(i))

        Unew(4,i)= U(4,i) - (dt*R(4,i)/A(i)) 

     end do


     call timestep(ncells,gamma,Unew, A, CFL, dt)

     !write(*,*) "dt" , dt
     !boundary - farfiel had -1 as a flag

     !create the grid again . with 2 kinds of boundaries. wall and inf

     ! for inf call the same Roe_interior.

     ! for wall call another modified version of Roe_interior

     !First time step should be really small

     !use euler time stepping and find U

     !find the primitive variables from that

     !use approach 2. for time step
     !get the elengths.
     !and find time for each step. 
     !find the minimum time step and use in the loop

     !Then calculate the time step

     err = Unew-U
     add=0.0
     do i=1,ncells
        do j=1,4
           add= add+(err(j,i)**2)
        end do
     end do
     rms=sqrt(add)

     write(*,*) "error",iter, rms
     U = Unew



     !End of time loop
     !=======================


     ! write out the grid file again

     ! and the solutin too.in cobalt format.

     !thats it

     !write(*,*) iter
  end do

  ! write(*,*) U(1,1:ncells)
  ! write(*,*) U(2,1:ncells)
  !write(*,*) U(3,1:ncells)
  ! write(*,*) U(4,1:ncells)

  allocate(density(ncells))
  allocate(u_vel(ncells))
  allocate(v_vel(ncells))
  allocate(pressure(ncells))

  do i=1,ncells

     density(i)= U(1,i)
     u_vel(i)= U(2,i)/U(1,i)
     V_vel(i)= U(3,i)/U(1,i)
     pressure(i)= (gamma-1)*(U(4,i)-(0.5*density(i)*(u_vel(i)**2+v_vel(i)**2)))
  end do

  open(33,file="output3.dat",form='formatted',status='replace',iostat=ios)
  if (ios /=0 ) then
	write (*,'(a)') 'Tecplot Write-Open - Fatal error!'
  end if
  open(34,file="density",form='formatted')

  write(33,*) "title = sample"
  write(33,*) "variables = x,y,d,u,v,p"
  write(33,*) "zone n=",nnodes,",e=",ncells,",DATAPACKING=BLOCK,ZONETYPE=FETRIANGLE,varlocation=([3,4,5,6]=cellcentered)"
  WRITE(33,'(3E14.4)') (nodexy(1,i),i=1,nnodes)
  write(33,'(3E14.4)') (nodexy(2,i),i=1,nnodes)
  write(33,'(3E14.4)') (density(i),i=1,ncells)
  write(33,'(3E14.4)') (u_vel(i),i=1,ncells)
  write(33,'(3E14.4)') (v_vel(i),i=1,ncells)
  write(33,'(3E14.4)') (pressure(i),i=1,ncells)
  do i=1,ncells
     write(33,*) (cell2node(j,i),j=1,3)
  end do
  close(33)
  close(34)
  close(11)

  deallocate(nodexy,edge2cell)

end program euler2D





