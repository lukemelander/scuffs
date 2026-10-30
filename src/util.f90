module utilities
implicit none
real*8, parameter :: pi = acos(-1.0d0)
real*8, parameter :: gam = 1.4d0
real*8, parameter :: gamm = 1.4d0-1.0d0
real*8, parameter :: rair = 287.05d0
real*8, parameter :: Pr = 0.7d0
real*8, parameter :: cv = (5.0d0/2.0d0)*Rair
!Implicit solver: kmax+1 point-Jacobi sweeps per iteration
!(10 sweeps needed for stability at CFL ~50-100)
integer*8, parameter :: kmax = 9
!First wall face (cell index i) of the flat plate, stype = 4
integer*8, parameter :: iplate = 11
real*8, allocatable :: nxy(:,:,:), cxy(:,:,:), df(:,:),dxy(:,:),ndf(:,:)
real*8 :: fs(3),CFL,Id(4,4)
real*8 :: Kq, mu,lamu
integer*8, allocatable :: ncf(:,:,:), nfc(:,:,:),bclist(:)
integer*8 :: stype,iblend,impl,ivisc,nx,ny,bcs(4),iter,isave,fnum
integer :: switch
!File names, set from the command line (see parse_args)
character(len=256) :: input_file = 'input.inp'
character(len=256) :: grid_file = ''
!--mesh-check: read and check an SU2 mesh, write VTK files, stop
logical :: mesh_check = .false.
contains
   subroutine parse_args()
   !Reads command-line options:
   !  --input FILE / -i FILE / --input=FILE   input file (default input.inp)
   !  --grid FILE  / -g FILE / --grid=FILE    grid file (default: from the
   !                                          case type in the input file)
   !  --mesh-check                          check an SU2 mesh and stop
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
    print '(A)','Usage: scuffs [--input FILE] [--grid FILE]'
    print '(A)','  -i, --input FILE   input file (default: input.inp)'
    print '(A)','  -g, --grid FILE    grid file (default: picked from the case type'
    print '(A)','                     in the input file, e.g. 4 -> plate.dat)'
    print '(A)','      --mesh-check   read the SU2 mesh given by --grid, check it,'
    print '(A)','                     write mesh.vtk and mesh_boundary.vtk, then stop'
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

   subroutine init()
   implicit none
   character(LEN=256) :: fname
   character(LEN=10), parameter :: fname0='block.dat'
   character(LEN=10), parameter :: fname1='ramp.dat'
   character(LEN=10), parameter :: fname2='rampb.dat'
   character(LEN=10), parameter :: fname3='cyl.dat'
   character(LEN=10), parameter :: fname4='plate.dat'
   character(LEN=15),dimension(4) :: bcn
   character(LEN=17) :: dum3
   character(LEN=7) :: dum1
   character(LEN=2) :: dum2
   integer*8 :: i,j,k 
   real*8 :: tt
      print*,'SCUFFS: Structured Compressible Upwind Flow Finite-volume Solver'
      print*,'    1st order inviscid fluxes'

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
      fname = ''
      if(stype == 0) then
       fname = fname0
       print*,'    Block'
      elseif(stype == 1) then
       fname = fname1
       print*,'    Channel with Ramp'
      elseif(stype == 2) then
       fname = fname2
       print*,'    Channel with Ramp (shock aligned grid)'
      elseif(stype == 3) then
       fname = fname3
       print*,'    Cylinder'
      elseif(stype == 4) then
       fname = fname4
       print*,'    Flat Plate'
      endif

      do k = 1,4
       if(bcs(k) == 1) then
        bcn(k) = 'inflow'
       elseif(bcs(k) == 2) then
        bcn(k) = 'outflow'
       elseif(bcs(k) == 3) then
        bcn(k) = 'inviscid wall'
       elseif(bcs(k) == 4) then
        bcn(k) = 'mirror'
       elseif(bcs(k) == 5) then
        bcn(k) = 'adiabatic wall'
       endif
      enddo

      !A grid given on the command line overrides the default for stype
      if(len_trim(grid_file) > 0) fname = grid_file
      if(len_trim(fname) == 0) then
        print*,'No grid file: give one with --grid or use case type 0-4'
        stop 1
      endif
      call check_file(fname)
      print*,'    Grid file: ',trim(fname)
      open(unit=20,file=trim(fname),status='old',action='read')
      read(20,*)
      read(20,*) dum1,dum2,nx,dum2,ny
       allocate(nxy(nx,ny,2))
      read(20,*) ((nxy(i,j,1),i=1,nx),j=1,ny)
      read(20,*) ((nxy(i,j,2),i=1,nx),j=1,ny)
      close(20)
      print*, 'Read in grid'
      print*,'    Nx =',nx
      print*,'    Ny =',ny
      print*,'    BCs:'
      print*,'      i=1 ,',bcn(2)
      print*,'      i=Nx,',bcn(1)
      print*,'      j=1 ,',bcn(4)
      print*,'      j=Ny,',bcn(3)
 
      print*,'-----------------------------'
      print*,'Freestream flow'
      print*,'   Mach #=',fs(1)
      print*,'   Rho =',fs(3)
      print*,'   P =',fs(2)
      print*,'   U =',fs(1)*(sqrt(gam*fs(2)/fs(3)))
 
      !Deriving some useful quantities
      Id(:,:) = 0.0d0
      do k = 1,4
        Id(k,k) = 1.0d0
      enddo
      Tt = fs(2)/(Rair*fs(3))
      mu = (1.458d-6)*(tt**1.5d0)/(tt+110.4d0)
      lamu = -(2.0d0/3.0d0)*mu
      Kq = gam*mu*cv/Pr
      print*,'Gas Properties:'
      print*,'   Temp: ',Tt
      print*,'   Dynamic Viscosity: ',mu
      print*,'   Thermal conductivity: ',Kq
      print*,'-----------------------------'
      print*,'Running for',iter,' iterations' 
  end subroutine

   subroutine cfl_ramp(itr,dt,cfl)
   implicit none
   integer*8, intent(in) :: itr
   real*8, intent(inout) :: cfl,dt
   integer*8 :: cflnum,i,j,k
   real*8, allocatable :: cfl_list(:)
   integer*8, allocatable :: cfl_iter(:)
   cflnum = 9
   allocate(cfl_iter(cflnum),cfl_list(cflnum))
   cfl_iter(1) = 1
   cfl_list(1) = 0.4d0  
   cfl_iter(2) = 20
   cfl_list(2) = 0.8d0
   cfl_iter(3) = 30
   cfl_list(3) = 1.6d0
   cfl_iter(4) = 40
   cfl_list(4) = 3.2d0
   cfl_iter(5) = 500
   cfl_list(5) = 6.4d0
   cfl_iter(6) = 700
   cfl_list(6) = 12.8d0
   cfl_iter(7) = 1000
   cfl_list(7) = 25.6d0   
   cfl_iter(8) = 1250
   cfl_list(8) = 50.0d0
   cfl_iter(9) = 1500
   cfl_list(9) = 100.0d0
    do i = 1,cflnum
     if(itr == cfl_iter(i)) then
       cfl = cfl_list(i)
     endif
    enddo
   deallocate(cfl_iter,cfl_list)
   end subroutine

   subroutine build_metrics()
   implicit none
   integer i,j,k,fi
   real*8 :: a1,a2
    !Building cell center array
    do i = 2,nx
    do j = 2,ny
      cxy(i,j,1) = (nxy(i-1,j-1,1)+nxy(i,j-1,1)+nxy(i-1,j,1)+nxy(i,j,1))/4.0d0
      cxy(i,j,2) = (nxy(i-1,j-1,2)+nxy(i,j-1,2)+nxy(i-1,j,2)+nxy(i,j,2))/4.0d0
    enddo
    enddo

    do i = 2,nx
     cxy(i,1,:) = -cxy(i,2,:) + (nxy(i,1,:) + nxy(i-1,1,:))
     cxy(i,ny+1,:) = -cxy(i,ny,:) + (nxy(i,ny,:) + nxy(i-1,ny,:))
    enddo

    do j = 2,ny
     cxy(1,j,:) = -cxy(2,j,:) + (nxy(1,j,:)+nxy(1,j-1,:))
     cxy(nx+1,j,:) = -cxy(nx,j,:) + (nxy(nx,j,:)+nxy(nx,j-1,:))
    enddo

    do i = 2,nx
    do j = 2,ny

       !Evaluating dxy
       call calc_area(nxy(i,j,:),nxy(i-1,j-1,:),nxy(i,j-1,:),a1)
       call calc_area(nxy(i,j,:),nxy(i-1,j-1,:),nxy(i-1,j,:),a2)

       dxy(i,j) = a1 + a2

    enddo
    enddo
    !Getting dNu and dEta for each face
    fi = 1
    allocate(ndf(fnum,3),df(fnum,5))
    do j = 1,ny-1
    do i = 1,nx
      !Surface Normals
      ndf(fi,1) = nxy(i,j+1,2) - nxy(i,j,2) 
      ndf(fi,2) = -(nxy(i,j+1,1) - nxy(i,j,1))
      ndf(fi,3) = sqrt(ndf(fi,1)*ndf(fi,1) + ndf(fi,2)*ndf(fi,2))
      ndf(fi,1:2) = ndf(fi,1:2)/ndf(fi,3)
      
      !Surface Derivatives (dx/de,dx/dn,dy/de,dy,dn)
      df(fi,2) = nxy(i,j+1,1) - nxy(i,j,1)
      df(fi,4) = nxy(i,j+1,2) - nxy(i,j,2)
      if(i==1)then
       df(fi,1) = 0.5d0*((nxy(i+1,j,1)-nxy(i,j,1)) + (nxy(i+1,j+1,1)-nxy(i,j+1,1)))
       df(fi,3) = 0.5d0*((nxy(i+1,j,2)-nxy(i,j,2)) + (nxy(i+1,j+1,2)-nxy(i,j+1,2)))
      elseif(i==nx) then
       df(fi,1) = 0.5d0*((nxy(i,j,1)-nxy(i-1,j,1)) + (nxy(i,j+1,1)-nxy(i-1,j+1,1)))
       df(fi,3) = 0.5d0*((nxy(i,j,2)-nxy(i-1,j,2)) + (nxy(i,j+1,2)-nxy(i-1,j+1,2)))
      else
       df(fi,1) = 0.5d0*((nxy(i+1,j,1)-nxy(i-1,j,1))/2.0d0 + (nxy(i+1,j+1,1)-nxy(i-1,j+1,1))/2.0d0)
       df(fi,3) = 0.5d0*((nxy(i+1,j,2)-nxy(i-1,j,2))/2.0d0 + (nxy(i+1,j+1,2)-nxy(i-1,j+1,2))/2.0d0)
      endif
      df(fi,5) = df(fi,1)*df(fi,4) - df(fi,2)*df(fi,3)
      fi = fi+1
    enddo
    enddo

    do i = 1,nx-1
    do j = 1,ny
      !Surface normals
      ndf(fi,1) = -(nxy(i+1,j,2) - nxy(i,j,2))
      ndf(fi,2) = (nxy(i+1,j,1) - nxy(i,j,1))
      ndf(fi,3) = sqrt(ndf(fi,1)*ndf(fi,1) + ndf(fi,2)*ndf(fi,2))
      ndf(fi,1:2) = ndf(fi,1:2)/ndf(fi,3)

      !Surface Derivatives
      df(fi,1) = nxy(i+1,j,1) - nxy(i,j,1)
      df(fi,3) = nxy(i+1,j,2) - nxy(i,j,2)
      if(j==1)then
       df(fi,2) = 0.5d0*((nxy(i,j+1,1)-nxy(i,j,1)) + (nxy(i+1,j+1,1)-nxy(i+1,j,1)))
       df(fi,4) = 0.5d0*((nxy(i,j+1,2)-nxy(i,j,2)) + (nxy(i+1,j+1,2)-nxy(i+1,j,2)))
      elseif(j==ny) then
       df(fi,2) = 0.5d0*((nxy(i,j,1)-nxy(i,j-1,1)) + (nxy(i+1,j,1)-nxy(i+1,j-1,1)))
       df(fi,4) = 0.5d0*((nxy(i,j,2)-nxy(i,j-1,2)) + (nxy(i+1,j,2)-nxy(i+1,j-1,2))) 
      else
       df(fi,2) = 0.5d0*((nxy(i,j+1,1)-nxy(i,j-1,1))/2.0d0 + (nxy(i+1,j+1,1)-nxy(i+1,j-1,1))/2.0d0)
       df(fi,4) = 0.5d0*((nxy(i,j+1,2)-nxy(i,j-1,2))/2.0d0 + (nxy(i+1,j+1,2)-nxy(i+1,j-1,2))/2.0d0)
      endif
      df(fi,5) = df(fi,1)*df(fi,4) - df(fi,2)*df(fi,3)
      fi = fi+1
    enddo
    enddo
   
   open(unit=20,file='checkdf.dat')
   do i = 2,nx
    write(20,*) df(ncf(i,2,3),:)
   enddo
   write(20,*) 'break'
      do i = 2,nx
    write(20,*) df(ncf(i,3,3),:)
   enddo
   close(20)
   end subroutine
 
   subroutine calc_area(p1,p2,p3,a)
   implicit none
   real*8, intent(in) :: p1(:),p2(:),p3(:)
   real*8,intent(inout) ::a
   real*8 :: v1(2),v2(2)
   v1 = p1(:)-p2(:)
   v2 = p3(:)-p2(:)
   a = abs(0.5d0*(v1(1)*v2(2) - v1(2)*v2(1)))

   end subroutine
   

   subroutine write_vars(U,fname)
   implicit none
   integer*8 :: i,j,k
   real*8, intent(in) :: U(:,:,:)
   real*8, allocatable :: m(:,:),rho(:,:),uu(:,:),uv(:,:),p(:,:)
   character(LEN=8), intent(in) :: fname
   allocate(m(nx+1,ny+1),rho(nx+1,ny+1),uu(nx+1,ny+1),uv(nx+1,ny+1),p(nx+1,ny+1))
   do i = 1,nx+1
   do j = 1,ny+1
   rho(i,j) = U(i,j,1)
   uu(i,j) = U(i,j,2)/U(i,j,1)
   uv(i,j) = U(i,j,3)/U(i,j,1)
   P(i,j) =  (gam-1.0d0)*U(i,j,4) - 0.5d0*U(i,j,1)*(gam-1.0d0)*((U(i,j,2)/U(i,j,1))**2 + (U(i,j,3)/U(i,j,1))**2)
   m(i,j) = sqrt(uu(i,j)**2.0d0 + uv(i,j)**2.0d0)/sqrt(gam*P(i,j)/rho(i,j))
   enddo
   enddo 
   
   open(unit=100,file=fname,action='write')
   !Tecplot ASCII header (also read by ParaView's Tecplot reader):
   !space-separated variable names, cell-centered range for vars 3-7
   write(100,'(A)') 'VARIABLES = "X" "Y" "rho" "u" "v" "p" "m"'
   write(100,'(A,I0,A,I0,A)') 'ZONE I=',nx,', J=',ny, &
           ', DATAPACKING=BLOCK, VARLOCATION=([3-7]=CELLCENTERED)'
   do j =1,ny
      write(100,*) (nxy(i,j,1),i=1,nx)
   enddo
   do j=1,ny
      write(100,*) (nxy(i,j,2),i=1,nx)
   enddo
   
   do j=2,ny
      write(100,*) (rho(i,j),i=2,nx)
   enddo
   do j=2,ny
      write(100,*) (uu(i,j), i=2,nx)
   enddo
   do j=2,ny
      write(100,*) (uv(i,j), i=2,nx)
   enddo
   do j=2,ny
      write(100,*) (p(i,j), i=2,nx)
   enddo
   do j=2,ny
      write(100,*) (m(i,j), i=2,nx)
   enddo

   close(100)
   end subroutine

   subroutine mat_mult(A,B,C,dims)
      implicit none
      integer*8 :: i,j,k
      integer*8, intent(in) :: dims(:)
      real*8, intent(in) :: A(dims(1),dims(2)),B(dims(3),dims(4))
      real*8, intent(inout) :: C(dims(1),dims(4))
      real*8 :: sum1
      real*8, allocatable :: AA(:,:), BB(:,:)
       C(:,:) = 0.0d0
       do i = 1,dims(1)
       do j = 1,dims(4)
        do k = 1,dims(3)
          C(i,j) = C(i,j) + A(i,k)*B(k,j)
        enddo
       enddo
       enddo

   end subroutine

   subroutine dfromp(U,V,sgn)
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
      elseif(sgn== -2) then
        V(1) = U(1)
        V(2) = U(2)/U(1)
        V(3) = U(3)/U(1)
        V(4) = gamm*(U(4) - 0.5d0*V(1)*(V(2)*V(2) + V(3)*V(3)))
      else
         print*,'not a valid option for dfromp'
         stop
      endif
   end subroutine

   subroutine calc_res(U1,U2,res)
     implicit none
     real*8,intent(in) :: u1(:,:,:),U2(:,:,:)
     real*8,intent(inout) :: res
     integer*8 :: i,j,k
     real*8 :: uref(4),ainf,uinf
     !RMS change of all four conserved variables over the interior cells,
     !each normalised by its freestream scale so they are comparable
     ainf = sqrt(gam*fs(2)/fs(3))
     uinf = fs(1)*ainf
     uref(1) = fs(3)
     uref(2) = fs(3)*max(uinf,ainf)
     uref(3) = uref(2)
     uref(4) = fs(2)/gamm + 0.5d0*fs(3)*uinf*uinf
     res = 0.0d0
     do i = 2,nx
     do j = 2,ny
      do k = 1,4
       res = res + ((U2(i,j,k)-U1(i,j,k))/uref(k))**2
      enddo
      if(isnan(res) .eqv. .true.) then
         print*,'NaN at i=',i,'j=',j
         return
      endif
     enddo
     enddo
     
     res = sqrt(res/dble(4*(nx-1)*(ny-1)))
   end subroutine
   
   subroutine get_dt(U,dt,CFL)
   implicit none
   real*8, intent(inout) :: dt,CFL,U(:,:,:)
   real*8 :: a,dx,d,V(4),q(2),Sx(2),Sy(2),vol,t1,t2,Sx2,Sy2,dtnew
   real*8 :: lamv
   integer*8 :: i,j,k
   !Start from a large value so dt is the true CFL limit each call
   dt = huge(1.0d0)
   do i =2,nx
   do j =2,ny
    
    call dfromp(U(i,j,:),V,-2)
    a = sqrt(gam*V(4)/V(1))
    
    q(1) = V(2)
    q(2) = V(3)

    Sx(1) = nxy(i,j,2) - nxy(i,j-1,2)
    Sx(2) = -(nxy(i,j,1) - nxy(i,j-1,1))
    
    Sy(1) = nxy(i,j,2) - nxy(i-1,j,2)
    Sy(2) = -(nxy(i,j,1) - nxy(i-1,j,1))

    vol = dxy(i,j)

    t1 = abs(q(1)*Sx(1) + q(2)*Sx(2))
    t2 = abs(q(1)*Sy(1) + q(2)*Sy(2))
    Sx2 = Sx(1)*Sx(1) + Sx(2)*Sx(2)
    Sy2 = Sy(1)*Sy(1) + Sy(2)*Sy(2)
    !Viscous spectral radius (thin-layer: j-direction only). Needed for
    !stability of the explicit scheme near the wall, where the grid is
    !very fine in j. The implicit scheme treats viscous terms implicitly.
    lamv = 0.0d0
    if(ivisc == 1 .AND. impl == 0) then
      lamv = max(4.0d0/3.0d0, gam/Pr)*(mu/V(1))*Sy2/vol
    endif
    dtnew = CFL*(vol/(t1 + t2 + a*sqrt(Sx2 + Sy2) + 2.0d0*lamv))
   
    if(dt .GT. dtnew) then
      dt = dtnew
    endif
   

   enddo
   enddo
   end subroutine

   subroutine inverse_4x4(A,B)
   implicit none
   real*8, intent(inout) :: A(:,:), B(:,:)
   integer*8 :: i,j,k
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
