module gradients
!-----------------------------------------------------------------------
! Cell gradients on the unstructured mesh.
!
! Weighted least squares (main method). For cell i with neighbours j
! (the cells across its faces; the ghost cell for a boundary face):
!
!    d_ij = x_j - x_i,   dphi_ij = phi_j - phi_i,   w_ij = 1/|d_ij|^2
!    minimise  sum_j w_ij (grad_i . d_ij - dphi_ij)^2
!
! which gives the 2x2 system G grad_i = r with
!
!    G = sum_j w d d^T = | a b |      r = sum_j w dphi d
!                        | b c |
!
! G depends only on the geometry, so its inverse (three numbers per cell)
! is computed once in setup_gradients. The result is exact for linear
! fields on any mesh, and first-order accurate in general.
!
! Both cells of a face see the same w*dphi*d (d and dphi both change sign),
! so G and r are assembled in one loop over faces.
!
! Green-Gauss (cross-check only):
!    grad_i = (1/A_i) sum_faces phi_f n_f L_f,   phi_f = (phi_L + phi_R)/2
! exact for linear fields only on uniform/parallelogram meshes.
!
! Field arrays are phi(nv, ncell+nbface): nv variables for every real cell
! followed by the ghost cells. Gradients are grad(2, nv, ncell).
!-----------------------------------------------------------------------
use mesh
implicit none
private

real*8, allocatable :: lsq(:,:)     !(3,ncell): inverse of G as (g11, g12, g22)

public :: setup_gradients, lsq_gradients, green_gauss_gradients, gradient_test

contains

   subroutine setup_gradients()
   !Least-squares matrices G for every cell, stored inverted
   implicit none
   integer :: f,cl,cr,c
   real*8 :: d(2),w,det
   real*8, allocatable :: G(:,:)
    if(allocated(lsq)) deallocate(lsq)
    allocate(lsq(3,ncell),G(3,ncell))
    G(:,:) = 0.0d0
    do f = 1,nface
      cl = fcell(1,f)
      cr = fcell(2,f)
      d = cctr(:,cr) - cctr(:,cl)
      w = 1.0d0/(d(1)*d(1) + d(2)*d(2))
      G(1,cl) = G(1,cl) + w*d(1)*d(1)
      G(2,cl) = G(2,cl) + w*d(1)*d(2)
      G(3,cl) = G(3,cl) + w*d(2)*d(2)
      if(cr <= ncell) then
        G(1,cr) = G(1,cr) + w*d(1)*d(1)
        G(2,cr) = G(2,cr) + w*d(1)*d(2)
        G(3,cr) = G(3,cr) + w*d(2)*d(2)
      endif
    enddo
    do c = 1,ncell
      det = G(1,c)*G(3,c) - G(2,c)*G(2,c)
      if(det <= 1.0d-12*G(1,c)*G(3,c)) then
        print*,'Least-squares gradient matrix is singular in cell',c, &
               ' at (x,y) =',cctr(1,c),cctr(2,c)
        print*,'(its neighbours all lie on a line: check the mesh)'
        stop 1
      endif
      lsq(1,c) =  G(3,c)/det
      lsq(2,c) = -G(2,c)/det
      lsq(3,c) =  G(1,c)/det
    enddo
    deallocate(G)
   end subroutine

   subroutine lsq_gradients(phi,grad)
   !Weighted least-squares gradients of every variable in phi
   implicit none
   real*8, intent(in) :: phi(:,:)
   real*8, intent(out) :: grad(:,:,:)
   integer :: f,cl,cr,c,k,nv
   real*8 :: d(2),w,t
   real*8, allocatable :: r(:,:,:)
    nv = size(phi,1)
    allocate(r(2,nv,ncell))
    r(:,:,:) = 0.0d0
    do f = 1,nface
      cl = fcell(1,f)
      cr = fcell(2,f)
      d = cctr(:,cr) - cctr(:,cl)
      w = 1.0d0/(d(1)*d(1) + d(2)*d(2))
      do k = 1,nv
        t = w*(phi(k,cr) - phi(k,cl))
        r(1,k,cl) = r(1,k,cl) + t*d(1)
        r(2,k,cl) = r(2,k,cl) + t*d(2)
        if(cr <= ncell) then
          r(1,k,cr) = r(1,k,cr) + t*d(1)
          r(2,k,cr) = r(2,k,cr) + t*d(2)
        endif
      enddo
    enddo
    do c = 1,ncell
      do k = 1,nv
        grad(1,k,c) = lsq(1,c)*r(1,k,c) + lsq(2,c)*r(2,k,c)
        grad(2,k,c) = lsq(2,c)*r(1,k,c) + lsq(3,c)*r(2,k,c)
      enddo
    enddo
    deallocate(r)
   end subroutine

   subroutine green_gauss_gradients(phi,grad)
   !Green-Gauss gradients with face values = average of the two cells
   implicit none
   real*8, intent(in) :: phi(:,:)
   real*8, intent(out) :: grad(:,:,:)
   integer :: f,cl,cr,c,k,nv
   real*8 :: t
    nv = size(phi,1)
    grad(:,:,:) = 0.0d0
    do f = 1,nface
      cl = fcell(1,f)
      cr = fcell(2,f)
      do k = 1,nv
        t = 0.5d0*(phi(k,cl) + phi(k,cr))*flen(f)
        grad(:,k,cl) = grad(:,k,cl) + t*fnrm(:,f)
        if(cr <= ncell) grad(:,k,cr) = grad(:,k,cr) - t*fnrm(:,f)
      enddo
    enddo
    do c = 1,ncell
      grad(:,:,c) = grad(:,:,c)/cvol(c)
    enddo
   end subroutine

   subroutine gradient_test()
   !Gradients of analytic fields (set at cell and ghost centres) compared
   !with the exact gradients at the cell centres:
   !   constant  phi = 5                 (both methods: zero)
   !   linear    phi = 1 + 2x - 3y       (least squares: exact on any mesh)
   !   quadratic phi = x^2 + xy + y^2    (first-order error in general)
   implicit none
   integer, parameter :: nv = 3
   real*8, allocatable :: phi(:,:),gls(:,:,:),ggg(:,:,:),ex(:,:,:)
   logical, allocatable :: bcell(:)
   integer :: c,f,k
   real*8 :: x,y,h
   character(len=10) :: names(nv)
    names = (/'constant  ','linear    ','quadratic '/)
    allocate(phi(nv,ncell+nbface),gls(2,nv,ncell),ggg(2,nv,ncell),ex(2,nv,ncell))
    allocate(bcell(ncell))
    do c = 1,ncell+nbface
      x = cctr(1,c)
      y = cctr(2,c)
      phi(1,c) = 5.0d0
      phi(2,c) = 1.0d0 + 2.0d0*x - 3.0d0*y
      phi(3,c) = x*x + x*y + y*y
      if(c <= ncell) then
        ex(:,1,c) = (/0.0d0,0.0d0/)
        ex(:,2,c) = (/2.0d0,-3.0d0/)
        ex(:,3,c) = (/2.0d0*x + y, x + 2.0d0*y/)
      endif
    enddo
    bcell = .false.
    do f = nintface+1,nface
      bcell(fcell(1,f)) = .true.
    enddo

    call lsq_gradients(phi,gls)
    call green_gauss_gradients(phi,ggg)

    !typical cell size, for reference
    h = sqrt(sum(cvol)/dble(ncell))
    print*,'Gradient test: max |grad - exact| over all cells and over boundary cells'
    print '(A,ES10.3)','    typical cell size h = ',h
    print '(4X,A10,2X,A22,2X,A22)','field','least squares','Green-Gauss'
    print '(4X,A10,2X,2A11,2X,2A11)','','all','boundary','all','boundary'
    do k = 1,nv
      print '(4X,A10,2X,2ES11.3,2X,2ES11.3)',names(k), &
        maxerr(gls,k,.false.),maxerr(gls,k,.true.), &
        maxerr(ggg,k,.false.),maxerr(ggg,k,.true.)
    enddo
    if(maxerr(gls,2,.false.) < 1.0d-10*maxval(abs(ex(:,2,:)))) then
      print*,'Least squares is exact for the linear field: PASS'
    else
      print*,'Least squares is NOT exact for the linear field: FAIL'
    endif
    deallocate(phi,gls,ggg,ex,bcell)
   contains
    real*8 function maxerr(g,k,bnd_only)
    real*8, intent(in) :: g(:,:,:)
    integer, intent(in) :: k
    logical, intent(in) :: bnd_only
    integer :: i
     maxerr = 0.0d0
     do i = 1,ncell
       if(bnd_only .AND. .NOT. bcell(i)) cycle
       maxerr = max(maxerr, sqrt((g(1,k,i)-ex(1,k,i))**2 + (g(2,k,i)-ex(2,k,i))**2))
     enddo
    end function
   end subroutine

end module gradients
