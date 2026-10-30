program main
!-----------------------------------------------------------------------
! SCUFFS: Shock-capturing Compressible Unstructured Finite-volume Flow Solver
! 2D, cell-centred, first order, on unstructured (triangle/quad) meshes.
!
!   scuffs --grid MESH.su2 [--input input.inp]
!   scuffs --grid MESH.su2 --mesh-check
!-----------------------------------------------------------------------
   use mesh
   use gradients
   use utilities
   use solver
   implicit none
   integer :: it
   real*8 :: dt,ttot,res
   real*8, allocatable :: U1(:,:),U2(:,:)

   !Command-line options
   call parse_args()
   if(len_trim(grid_file) == 0) then
     print*,'A mesh file is required: --grid FILE.su2'
     call print_usage()
     stop 1
   endif
   call check_file(grid_file)

   !Mesh and the least-squares gradient matrices
   call mesh_setup(trim(grid_file))
   call setup_gradients()

   !Gradient test mode: check the gradients on this mesh and stop
   if(grad_test) then
     call gradient_test()
     stop
   endif

   !Mesh check mode: write VTK files of the mesh and stop
   if(mesh_check) then
     call write_mesh_vtk('mesh.vtk')
     call write_boundary_vtk('mesh_boundary.vtk')
     stop
   endif

   !Input file (needs the mesh's marker names)
   print*,'Input file: ',trim(input_file)
   call read_input(input_file)
   if(ivisc == 1) then
     print*,'Viscous flow is not available yet on unstructured meshes (roadmap step A6).'
     print*,'Set ivisc = 0 in the input file.'
     stop 1
   endif
   call init_gas()

   !Initial condition. U holds the real cells followed by one ghost cell
   !per boundary face.
   allocate(U1(4,ncell+nbface),U2(4,ncell+nbface))
   iblend = 1
   call init_flow(U1)
   call set_ghosts(U1)
   U2 = U1
   ttot = 0.0d0
   res = 0.0d0

   print*,'...starting iterations'
   open(unit=99,file='residual.dat',action='write',status='replace')

   do it = 1,iter
     if(cflramp == 1) call cfl_ramp(it,CFL)
     !Wall and mirror ghosts from the current interior states
     call set_ghosts(U1)
     call get_dt(U1,dt)
     if(impl == 0) then
       call solver_exp(dt,U1,U2)
     else
       call solver_impl(dt,U1,U2)
     endif
     call calc_res(U1,U2,res)
     ttot = ttot + dt

     !Stop if the run diverged
     if(isnan(res)) exit
     print*,it,ttot,CFL,dt,res
     write(99,*) it,res

     if(modulo(it,isave) == 0) then
       call set_ghosts(U2)
       call write_solution_vtk(U2,'data.vtk')
       print*,'saved to data.vtk'
     endif
     !Real cells and ghosts (outflow ghost states carry over)
     U1 = U2
   enddo
   close(99)

   !Final solution (unless the run diverged)
   if(.NOT. isnan(res)) then
     call set_ghosts(U1)
     call write_solution_vtk(U1,'data.vtk')
     print*,'final solution saved to data.vtk'
   endif
   print*,'...done'
end program
