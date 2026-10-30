module utilities
!-----------------------------------------------------------------------
! Gas properties, run settings, command-line options, and helpers shared
! by the solver: state conversions, time step, residual, CFL ramp,
! solution output and small matrix routines.
!
! Cell states are stored as U(4, ncell+nbface): conserved variables
! (rho, rho*u, rho*v, rho*E) for each real cell, followed by one ghost
! cell per boundary face (see the mesh module).
!-----------------------------------------------------------------------
use mesh
use gradients
implicit none
real*8, parameter :: pi = acos(-1.0d0)
real*8, parameter :: gam = 1.4d0
real*8, parameter :: gamm = 1.4d0-1.0d0
real*8, parameter :: rair = 287.05d0
real*8, parameter :: Pr = 0.7d0
real*8, parameter :: cv = (5.0d0/2.0d0)*Rair
!Implicit solver: kmax+1 point-Jacobi sweeps per iteration
!(10 sweeps needed for stability at CFL ~50-100)
integer, parameter :: kmax = 9

!Boundary-condition types
integer, parameter :: BC_INFLOW = 1, BC_OUTFLOW = 2, BC_SLIPWALL = 3
integer, parameter :: BC_MIRROR = 4, BC_ADIABATIC = 5

!Freestream (fs = Mach, p, rho), run settings
real*8 :: fs(3),CFL,Id(4,4)
real*8 :: Kq,mu,lamu
integer :: iblend,impl,ivisc,iter,isave,cflramp
integer :: switch
!BC type of each mesh marker (filled from the input file)
integer, allocatable :: mbc(:)

!File names, set from the command line (see parse_args)
character(len=256) :: input_file = 'input.inp'
character(len=256) :: grid_file = ''
!--mesh-check: read and check an SU2 mesh, write VTK files, stop
logical :: mesh_check = .false.
!--gradient-test: check the cell gradients on the mesh and stop
logical :: grad_test = .false.
contains
   subroutine parse_args()
   !Reads command-line options:
   !  --input FILE / -i FILE / --input=FILE   input file (default input.inp)
   !  --grid FILE  / -g FILE / --grid=FILE    SU2 mesh file (required)
   !  --mesh-check                            check an SU2 mesh and stop
   !  --gradient-test                         check the cell gradients and stop
   !  --help / -h
   implicit none
   integer :: i,nargs,eq
   character(len=256) :: arg,key,val
   logical :: have_val
    nargs = command_argument_count()
    i = 1
    do while(i <= nargs)
      call get_command_argument(i,arg)
      eq = index(arg,'=')
      if(eq > 0) then
        key = arg(1:eq-1)
        val = arg(eq+1:)
        have_val = .true.
      else
        key = arg
        val = ''
        have_val = .false.
      endif
      select case(trim(key))
      case('--input','-i','--grid','-g')
        if(.NOT. have_val) then
          if(i == nargs) then
            print*,'Missing file name after ',trim(key)
            call print_usage()
            stop 1
          endif
          i = i+1
          call get_command_argument(i,val)
        endif
        if(trim(key) == '--input' .OR. trim(key) == '-i') then
          input_file = val
        else
          grid_file = val
        endif
      case('--mesh-check')
        mesh_check = .true.
      case('--gradient-test')
        grad_test = .true.
      case('--help','-h')
        call print_usage()
        stop
      case default
        print*,'Unknown argument: ',trim(arg)
        call print_usage()
        stop 1
      end select
      i = i+1
    enddo
   end subroutine

   subroutine print_usage()
   implicit none
    print '(A)','Usage: scuffs --grid MESH.su2 [--input FILE] [--mesh-check]'
    print '(A)','  -g, --grid FILE    SU2 mesh file (required)'
    print '(A)','  -i, --input FILE   input file (default: input.inp)'
    print '(A)','      --mesh-check   check the mesh, write mesh.vtk and'
    print '(A)','                     mesh_boundary.vtk, then stop'
    print '(A)','      --gradient-test  check the cell gradients on the mesh, then stop'
    print '(A)','  -h, --help         show this message'
   end subroutine

   subroutine check_file(fname)
   !Stops with a clear message if a file does not exist
   implicit none
   character(len=*), intent(in) :: fname
   logical :: ex
    inquire(file=trim(fname),exist=ex)
    if(.NOT. ex) then
      print*,'File not found: ',trim(fname)
      stop 1
    endif
   end subroutine

   subroutine read_input(fname)
   !Input file (list-directed; anything after the values on a line is a
   !comment):
   !  ivisc impl cfl_ramp
   !  p_inf rho_inf M_inf
   !  iterations save_interval
   !  CFL
   !  nbc                       number of boundary lines that follow
   !  marker_name bc_type       (nbc lines)
   !  initial_condition
   !The mesh must already be loaded, so marker names can be checked.
   implicit none
   character(len=*), intent(in) :: fname
   integer, parameter :: un = 10
   integer :: nbc,k,m,btype
   character(len=64) :: bname
    call check_file(fname)
    open(unit=un,file=trim(fname),status='old',action='read')
    read(un,*) ivisc, impl, cflramp
    read(un,*) fs(2), fs(3), fs(1)
    read(un,*) iter, isave
    read(un,*) CFL
    read(un,*) nbc
    allocate(mbc(nmark))
    mbc = 0
    do k = 1,nbc
      read(un,*) bname, btype
      m = marker_id(bname)
      if(m == 0) then
        print*,'Input error: marker "',trim(bname),'" is not in the mesh. Mesh markers:'
        do m = 1,nmark
          print*,'   ',trim(markname(m))
        enddo
        stop 1
      endif
      if(btype < 1 .OR. btype > 5) then
        print*,'Input error: BC type for marker "',trim(bname),'" must be 1-5'
        stop 1
      endif
      mbc(m) = btype
    enddo
    read(un,*) switch
    close(un)
    do m = 1,nmark
      if(mbc(m) == 0) then
        print*,'Input error: no boundary condition given for marker "',trim(markname(m)),'"'
        stop 1
      endif
    enddo
    if(isave < 1) isave = iter
   end subroutine

   subroutine init_gas()
   !Derived gas properties and a summary of the run
   implicit none
   integer :: k,m
   real*8 :: tt
   character(len=15) :: bcn(5)
    bcn = (/'inflow        ','outflow       ','slip wall     ','mirror        ', &
            'adiabatic wall'/)
    print*,'SCUFFS: Shock-capturing Compressible Unstructured Finite-volume Flow Solver'
    print*,'    1st order inviscid fluxes (MSW)'
    if(ivisc == 0) then
      print*,'    inviscid flow'
    else
      print*,'    viscous flow'
    endif
    if(impl == 0) then
      print*,'    1st order explicit time integration'
    else
      print*,'    1st order implicit time integration'
    endif
    print*,'    Boundary conditions:'
    do m = 1,nmark
      print '(8X,A20,A)',markname(m),trim(bcn(mbc(m)))
    enddo
    print*,'-----------------------------'
    print*,'Freestream flow'
    print*,'   Mach #=',fs(1)
    print*,'   Rho =',fs(3)
    print*,'   P =',fs(2)
    print*,'   U =',fs(1)*(sqrt(gam*fs(2)/fs(3)))

    Id(:,:) = 0.0d0
    do k = 1,4
      Id(k,k) = 1.0d0
    enddo
    tt = fs(2)/(Rair*fs(3))
    mu = (1.458d-6)*(tt**1.5d0)/(tt+110.4d0)
    lamu = -(2.0d0/3.0d0)*mu
    Kq = gam*mu*cv/Pr
    print*,'Gas Properties:'
    print*,'   Temp: ',tt
    print*,'   Dynamic Viscosity: ',mu
    print*,'   Thermal conductivity: ',Kq
    print*,'-----------------------------'
    print*,'Running for',iter,' iterations'
   end subroutine

   subroutine cfl_ramp(itr,cfl)
   !CFL schedule for implicit steady-state runs
   implicit none
   integer, intent(in) :: itr
   real*8, intent(inout) :: cfl
   integer, parameter :: cflnum = 9
   integer, parameter :: cfl_iter(cflnum) = (/1,20,30,40,500,700,1000,1250,1500/)
   real*8, parameter :: cfl_list(cflnum) = &
        (/0.4d0,0.8d0,1.6d0,3.2d0,6.4d0,12.8d0,25.6d0,50.0d0,100.0d0/)
   integer :: i
    do i = 1,cflnum
      if(itr == cfl_iter(i)) cfl = cfl_list(i)
    enddo
   end subroutine

   subroutine write_solution_vtk(U,fname)
   !Cell-centred rho, u, v, p, Mach, |grad rho| and velocity on the mesh
   !(legacy VTK). The ghost cells of U must be set (set_ghosts), since the
   !density gradient uses them.
   implicit none
   real*8, intent(in) :: U(:,:)
   character(len=*), intent(in) :: fname
   integer, parameter :: un = 100
   integer :: c,k
   real*8 :: V(4),Uc(4)
   real*8, allocatable :: P(:,:),Q(:,:),grad(:,:,:)
    !primitive variables of every cell (ghosts included)
    allocate(Q(4,ncell+nbface),P(6,ncell),grad(2,4,ncell))
    do c = 1,ncell+nbface
      Uc = U(:,c)
      call dfromp(Uc,V,-2)
      Q(:,c) = V
    enddo
    call lsq_gradients(Q,grad)
    do c = 1,ncell
      P(1:4,c) = Q(:,c)
      P(5,c) = sqrt(Q(2,c)**2 + Q(3,c)**2)/sqrt(gam*Q(4,c)/Q(1,c))
      P(6,c) = sqrt(grad(1,1,c)**2 + grad(2,1,c)**2)
    enddo

    open(unit=un,file=trim(fname),status='replace',action='write')
    write(un,'(A)') '# vtk DataFile Version 3.0'
    write(un,'(A)') 'SCUFFS solution'
    write(un,'(A)') 'ASCII'
    write(un,'(A)') 'DATASET UNSTRUCTURED_GRID'
    write(un,'(A,I0,A)') 'POINTS ',nnode,' double'
    do k = 1,nnode
      write(un,'(3ES24.15)') xy(1,k),xy(2,k),0.0d0
    enddo
    write(un,'(A,I0,1X,I0)') 'CELLS ',ncell,sum(nvert)+ncell
    do c = 1,ncell
      write(un,'(I0,4(1X,I0))') nvert(c),(cnode(k,c)-1,k=1,nvert(c))
    enddo
    write(un,'(A,I0)') 'CELL_TYPES ',ncell
    do c = 1,ncell
      if(nvert(c) == 3) then
        write(un,'(I0)') 5
      else
        write(un,'(I0)') 9
      endif
    enddo
    write(un,'(A,I0)') 'CELL_DATA ',ncell
    call scalar('rho',P(1,:))
    call scalar('u',P(2,:))
    call scalar('v',P(3,:))
    call scalar('p',P(4,:))
    call scalar('mach',P(5,:))
    call scalar('grad_rho',P(6,:))
    write(un,'(A)') 'VECTORS velocity double'
    do c = 1,ncell
      write(un,'(3ES24.15)') P(2,c),P(3,c),0.0d0
    enddo
    close(un)
    deallocate(P,Q,grad)
   contains
    subroutine scalar(name,vals)
    character(len=*), intent(in) :: name
    real*8, intent(in) :: vals(:)
    integer :: i
     write(un,'(A,A,A)') 'SCALARS ',name,' double 1'
     write(un,'(A)') 'LOOKUP_TABLE default'
     do i = 1,size(vals)
       write(un,'(ES24.15)') vals(i)
     enddo
    end subroutine
   end subroutine

   subroutine mat_mult(A,B,C,dims)
   !C = A*B for small explicit-shape matrices, dims = (rows A, cols A,
   !rows B, cols B)
   implicit none
   integer, intent(in) :: dims(:)
   real*8, intent(in) :: A(dims(1),dims(2)),B(dims(3),dims(4))
   real*8, intent(inout) :: C(dims(1),dims(4))
   integer :: i,j,k
    C(:,:) = 0.0d0
    do j = 1,dims(4)
    do k = 1,dims(3)
    do i = 1,dims(1)
      C(i,j) = C(i,j) + A(i,k)*B(k,j)
    enddo
    enddo
    enddo
   end subroutine

   subroutine dfromp(U,V,sgn)
   !Conversions between conserved U and primitive V:
   !  sgn =  2: V = (rho,u,v,p) -> U      sgn = -2: U -> V = (rho,u,v,p)
   !  sgn =  1: V = (rho,u,v,T) -> U      sgn = -1: U -> V = (rho,u,v,T)
   implicit none
   real*8, intent(inout) :: U(:)
   real*8, intent(inout) :: V(:)
   integer, intent(in) :: sgn
    if(sgn == 1) then
      U(1) = V(1)
      U(2) = V(1)*V(2)
      U(3) = V(1)*V(3)
      U(4) = V(1)*cv*V(4) + 0.5d0*V(1)*(V(2)*V(2) + V(3)*V(3))
    elseif(sgn == -1) then
      V(1) = U(1)
      V(2) = U(2)/U(1)
      V(3) = U(3)/U(1)
      V(4) = (U(4)/U(1) -  0.5d0*(V(2)*V(2) + V(3)*V(3)))/cv
    elseif(sgn == 2) then
      U(1) = V(1)
      U(2) = V(1)*V(2)
      U(3) = V(1)*V(3)
      U(4) = V(4)/gamm + (0.5d0)*(V(1))*(V(2)*V(2) + V(3)*V(3))
    elseif(sgn == -2) then
      V(1) = U(1)
      V(2) = U(2)/U(1)
      V(3) = U(3)/U(1)
      V(4) = gamm*(U(4) - 0.5d0*V(1)*(V(2)*V(2) + V(3)*V(3)))
    else
      print*,'not a valid option for dfromp'
      stop 1
    endif
   end subroutine

   subroutine calc_res(U1,U2,res)
   !RMS change of the four conserved variables over the real cells, each
   !normalised by its freestream scale so they are comparable
   implicit none
   real*8, intent(in) :: U1(:,:),U2(:,:)
   real*8, intent(out) :: res
   integer :: c,k
   real*8 :: uref(4),ainf,uinf
    ainf = sqrt(gam*fs(2)/fs(3))
    uinf = fs(1)*ainf
    uref(1) = fs(3)
    uref(2) = fs(3)*max(uinf,ainf)
    uref(3) = uref(2)
    uref(4) = fs(2)/gamm + 0.5d0*fs(3)*uinf*uinf
    res = 0.0d0
    do c = 1,ncell
      do k = 1,4
        res = res + ((U2(k,c)-U1(k,c))/uref(k))**2
      enddo
      if(isnan(res)) then
        print*,'NaN in cell',c,' at (x,y) =',cctr(1,c),cctr(2,c)
        return
      endif
    enddo
    res = sqrt(res/dble(4*ncell))
   end subroutine

   subroutine get_dt(U,dt)
   !Global time step: dt = CFL * min over cells of
   !   area / (1/2 * sum over faces of (|u.n| + a) * length)
   implicit none
   real*8, intent(in) :: U(:,:)
   real*8, intent(out) :: dt
   integer :: c,k,f
   real*8 :: V(4),Uc(4),a,lam
    dt = huge(1.0d0)
    do c = 1,ncell
      Uc = U(:,c)
      call dfromp(Uc,V,-2)
      a = sqrt(gam*V(4)/V(1))
      lam = 0.0d0
      do k = 1,ncface(c)
        f = cface(k,c)
        lam = lam + (abs(V(2)*fnrm(1,f) + V(3)*fnrm(2,f)) + a)*flen(f)
      enddo
      dt = min(dt, CFL*cvol(c)/(0.5d0*lam))
    enddo
   end subroutine

   subroutine inverse_4x4(A,B)
   implicit none
   real*8, intent(inout) :: A(:,:), B(:,:)
   real*8 :: detA
   !Calculating Determinant of A
   detA = A(1,1)*A(2,2)*A(3,3)*A(4,4) + A(1,1)*A(2,3)*A(3,4)*A(4,2) + A(1,1)*A(2,4)*A(3,2)*A(4,3)
   detA = detA + A(1,2)*A(2,1)*A(3,4)*A(4,3) + A(1,2)*A(2,3)*A(3,1)*A(4,4) + A(1,2)*A(2,4)*A(3,3)*A(4,1)
   detA = detA + A(1,3)*A(2,1)*A(3,2)*A(4,4) + A(1,3)*A(2,2)*A(3,4)*A(4,1) + A(1,3)*A(2,4)*A(3,1)*A(4,2)
   detA = detA + A(1,4)*A(2,1)*A(3,3)*A(4,2) + A(1,4)*A(2,2)*A(3,1)*A(4,3) + A(1,4)*A(2,3)*A(3,2)*A(4,1)
   detA = detA - A(1,1)*A(2,2)*A(3,4)*A(4,3) - A(1,1)*A(2,3)*A(3,2)*A(4,4) - A(1,1)*A(2,4)*A(3,3)*A(4,2)
   detA = detA - A(1,2)*A(2,1)*A(3,3)*A(4,4) - A(1,2)*A(2,3)*A(3,4)*A(4,1) - A(1,2)*A(2,4)*A(3,1)*A(4,3)
   detA = detA - A(1,3)*A(2,1)*A(3,4)*A(4,2) - A(1,3)*A(2,2)*A(3,1)*A(4,4) - A(1,3)*A(2,4)*A(3,2)*A(4,1)
   detA = detA - A(1,4)*A(2,1)*A(3,2)*A(4,3) - A(1,4)*A(2,2)*A(3,3)*A(4,1) - A(1,4)*A(2,3)*A(3,1)*A(4,2)
   
   !First row of B
   B(1,1) = A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) 
   B(1,1) = B(1,1) - A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4) - A(2,4)*A(3,3)*A(4,2)
   B(1,2) = A(1,2)*A(3,4)*A(4,3) + A(1,3)*A(3,2)*A(4,4) + A(1,4)*A(3,3)*A(4,2)
   B(1,2) = B(1,2) - A(1,2)*A(3,3)*A(4,4) - A(1,3)*A(3,4)*A(4,2) - A(1,4)*A(3,2)*A(4,3)
   B(1,3) = A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3)
   B(1,3) = B(1,3) - A(1,2)*A(2,4)*A(4,3) - A(1,3)*A(2,2)*A(4,4) - A(1,4)*A(2,3)*A(4,2)
   B(1,4) = A(1,2)*A(2,4)*A(3,3) + A(1,3)*A(2,2)*A(3,4) + A(1,4)*A(2,3)*A(3,2)
   B(1,4) = B(1,4) - A(1,2)*A(2,3)*A(3,4) - A(1,3)*A(2,4)*A(3,2) - A(1,4)*A(2,2)*A(3,3)

   !Second row of B
   B(2,1) = A(2,1)*A(3,4)*A(4,3) + A(2,3)*A(3,1)*A(4,4) + A(2,4)*A(3,3)*A(4,1) 
   B(2,1) = B(2,1) - A(2,1)*A(3,3)*A(4,4) - A(2,3)*A(3,4)*A(4,1) - A(2,4)*A(3,1)*A(4,3)
   B(2,2) = A(1,1)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,3) 
   B(2,2) = B(2,2) - A(1,1)*A(3,4)*A(4,3) - A(1,3)*A(3,1)*A(4,4) - A(1,4)*A(3,3)*A(4,1)
   B(2,3) = A(1,1)*A(2,4)*A(4,3) + A(1,3)*A(2,1)*A(4,4) + A(1,4)*A(2,3)*A(4,1) 
   B(2,3) = B(2,3) - A(1,1)*A(2,3)*A(4,4) - A(1,3)*A(2,4)*A(4,1) - A(1,4)*A(2,1)*A(4,3)
   B(2,4) = A(1,1)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,3) 
   B(2,4) = B(2,4) - A(1,1)*A(2,4)*A(3,3) - A(1,3)*A(2,1)*A(3,4) - A(1,4)*A(2,3)*A(3,1)

   !Third row of B
   B(3,1) = A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) 
   B(3,1) = B(3,1) - A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4) - A(2,4)*A(3,2)*A(4,1)
   B(3,2) = A(1,1)*A(3,4)*A(4,2) + A(1,2)*A(3,1)*A(4,4) + A(1,4)*A(3,2)*A(4,1) 
   B(3,2) = B(3,2) - A(1,1)*A(3,2)*A(4,4) - A(1,2)*A(3,4)*A(4,1) - A(1,4)*A(3,1)*A(4,2)
   B(3,3) = A(1,1)*A(2,2)*A(4,4) + A(1,2)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,2) 
   B(3,3) = B(3,3) - A(1,1)*A(2,4)*A(4,2) - A(1,2)*A(2,1)*A(4,4) - A(1,4)*A(2,2)*A(4,1)
   B(3,4) = A(1,1)*A(2,4)*A(3,2) + A(1,2)*A(2,1)*A(3,4) + A(1,4)*A(2,2)*A(3,1) 
   B(3,4) = B(3,4) - A(1,1)*A(2,2)*A(3,4) - A(1,2)*A(2,4)*A(3,1) - A(1,4)*A(2,1)*A(3,2)
   !Fourth row of B
   B(4,1) = A(2,1)*A(3,3)*A(4,2) + A(2,2)*A(3,1)*A(4,3) + A(2,3)*A(3,2)*A(4,1) 
   B(4,1) = B(4,1) - A(2,1)*A(3,2)*A(4,3) - A(2,2)*A(3,3)*A(4,1) - A(2,3)*A(3,1)*A(4,2)
   B(4,2) = A(1,1)*A(3,2)*A(4,3) + A(1,2)*A(3,3)*A(4,1) + A(1,3)*A(3,1)*A(4,2) 
   B(4,2) = B(4,2) - A(1,1)*A(3,3)*A(4,2) - A(1,2)*A(3,1)*A(4,3) - A(1,3)*A(3,2)*A(4,1)
   B(4,3) = A(1,1)*A(2,3)*A(4,2) + A(1,2)*A(2,1)*A(4,3) + A(1,3)*A(2,2)*A(4,1) 
   B(4,3) = B(4,3) - A(1,1)*A(2,2)*A(4,3) - A(1,2)*A(2,3)*A(4,1) - A(1,3)*A(2,1)*A(4,2)
   B(4,4) = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) 
   B(4,4) = B(4,4) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) - A(1,3)*A(2,2)*A(3,1)

   !Multiplying by 1/det(A)
   B = (1.0d0/detA)*B(:,:)
   end subroutine

end module
