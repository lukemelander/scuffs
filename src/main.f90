program main
   use utilities
   use solver
   implicit none
   integer*8 :: i,j,k
   real*8 :: a,ttot
   !Read in variables
   integer*8 :: ivers, nb,cn,fi,cflramp
   real*8 :: th, c
   real*8 :: dt,Fc(4),Uc(2,4),V(4),np(2,2)
   real*8 :: dym,dxm,res,aa,aoa,w,dk(4)
   real*8, allocatable :: U1(:,:,:),U2(:,:,:)
   
   !Command-line options (--input, --grid)
   call parse_args()

   !Reading in input and grid files
   call check_file(input_file)
   print*,'Input file: ',trim(input_file)
   open(unit=10,file=trim(input_file),status='old',action='read')
     read(10,*) ivisc, impl, cflramp
     read(10,*) fs(2),fs(3),fs(1)
     read(10,*) stype, iter, isave
     read(10,*) CFL
     read(10,*) bcs(1),bcs(2),bcs(3),bcs(4)
     read(10,*) switch
   close(10)

   call init()
   !Computing grid metrics and connectivity
      print*,'Computing grid metrics and connectivity'
       allocate(cxy(nx+1,ny+1,2))                    !Cell center xy coords array
       allocate(dxy(nx+1,ny+1))    !df = face derivatives, dxy = volume
       fnum = nx*(ny-1)+ny*(nx-1)
       allocate(nfc(fnum,2,2),bclist(fnum))
       fi = 1
       do j = 2,ny
       do i = 1,nx
        nfc(fi,1,1) = i
        nfc(fi,1,2) = j
        nfc(fi,2,1) = i+1
        nfc(fi,2,2) = j
        if(i==1) then
          bclist(fi) = bcs(2)
        elseif(i==nx) then
          bclist(fi) = bcs(1)
        else
          bclist(fi) = 0
        endif
        fi = fi+1
       enddo
       enddo
       do i = 2,nx
       do j = 1,ny
        nfc(fi,1,1) = i
        nfc(fi,1,2) = j
        nfc(fi,2,1) = i
        nfc(fi,2,2) = j+1
        if(j==1) then
          if(stype == 4) then
            if(i .GE. iplate) then
              bclist(fi) = bcs(4)
            else
              bclist(fi) = 3
            endif
          else 
           bclist(fi) = bcs(4)
          endif 
        elseif(j==ny) then
          bclist(fi) = bcs(3)
        else
          bclist(fi) = 0
        endif
        fi = fi+1
       enddo
       enddo

       allocate(ncf(nx+1,ny+1,4))
       fi = 1
       do j = 2,ny
       do i = 2,nx
        ncf(i,j,1) = fi
        ncf(i,j,2) = fi+1
        fi = fi+1
       enddo
       fi=fi+1
       enddo
       do i = 2,nx
       do j = 2,ny
         ncf(i,j,3) = fi
         ncf(i,j,4) = fi+1
         fi = fi+1
       enddo
       fi = fi+1
       enddo
       print*,'   faces: ',fnum

      call build_metrics()
      print*,'done'
     
   !Step 3: Init flow, call main solver subroutine (using read in case # from input file)
      print*,'...Initializing flow'
       allocate(U1(nx+1,ny+1,4),U2(nx+1,ny+1,4))
       dt = 1.0d0
       ttot = 0.0d0
       iblend = 1
       call init_flow(U1) 
       call init_flow(U2)
       call get_dt(U1,dt,CFL)
      print*,'done'

      print*,'...starting iterations'
       open(unit=99,file='residual.dat',action='write')

       !Main solver loop       
       do i = 1,iter
         if(cflramp==1) then
          call cfl_ramp(i,dt,CFL)
         endif
         call get_dt(U1,dt,CFL)
         if(impl == 0) then
          call solver_exp(dt,U1,U2)
         elseif(impl == 1) then
          call solver_impl(dt,U1,U2)
         endif   
         call calc_res(U1,U2,res)
         ttot = ttot+dt

         !Checking for nan or finished
         if(isnan(res) .EQV. .TRUE.) exit
         print*,i,ttot,cfl,dt,res
         write(99,*) i,res

         !Writing out soln
         if(modulo(i,isave) == 0) then
            call write_vars(U2,'data.dat')
            print*,"saved to data.dat"
         endif
         U1(1:nx+1,1:ny+1,:) = U2(1:nx+1,1:ny+1,:)
       enddo
       close(99)

       !Write the final solution (unless the run diverged)
       if(.NOT. isnan(res)) then
          call write_vars(U1,'data.dat')
          print*,"final solution saved to data.dat"
       endif
      print*,'...done'
end program
