module solver
!-----------------------------------------------------------------------
! First-order finite-volume solver on the unstructured mesh.
!
! Fluxes are computed in one loop over faces. Each face f has a left cell
! fcell(1,f) and a right cell fcell(2,f); the unit normal fnrm(:,f) points
! from left to right. For a boundary face the right cell is the ghost cell
! of that face, so boundary conditions only ever set the ghost (side 2).
!
! Flux convention: flux(:,f) is the flux through face f in the direction
! of its normal, already multiplied by the face length. Cell residuals are
!   R(:,left)  += flux(:,f)      R(:,right) -= flux(:,f)
! and dU/dt = -R/area.
!
! Flux Jacobians: A(:,:,1,f) = dF/dU(left), A(:,:,2,f) = dF/dU(right).
!-----------------------------------------------------------------------
use utilities
implicit none

contains

   !--------------------------------------------------------------------
   ! Initial condition
   !--------------------------------------------------------------------
   subroutine init_flow(U)
   !Initial state of every cell, including the ghost cells (inflow ghosts
   !keep this state for the whole run)
   !  switch = 1: uniform freestream
   !  switch = 2: shock tube, gas at rest; cells below the middle of the
   !              mesh's y range at 2x pressure and density
   !  switch = 3: freestream with a 1% bump in the cells at the middle of
   !              the mesh's x range
   implicit none
   real*8, intent(inout) :: U(:,:)
   real*8 :: a,ufs,V(4),xmid,ymid
   integer :: c
    a = sqrt(gam*fs(2)/fs(3))
    ufs = fs(1)*a
    xmid = 0.5d0*(minval(xy(1,:)) + maxval(xy(1,:)))
    ymid = 0.5d0*(minval(xy(2,:)) + maxval(xy(2,:)))
    do c = 1,ncell+nbface
      select case(switch)
      case(2)
        if(cctr(2,c) < ymid) then
          V = (/2.0d0*fs(3),0.0d0,0.0d0,2.0d0*fs(2)/)
        else
          V = (/fs(3),0.0d0,0.0d0,fs(2)/)
        endif
      case(3)
        V = (/fs(3),ufs,0.0d0,fs(2)/)
        if(c <= ncell) then
          if(abs(cctr(1,c)-xmid) <= 0.5d0*sqrt(cvol(c))) then
            V = (/1.01d0*fs(3),1.01d0*ufs,0.0d0,1.01d0*fs(2)/)
          endif
        endif
      case default
        V = (/fs(3),ufs,0.0d0,fs(2)/)
      end select
      call dfromp(U(:,c),V,2)
    enddo
   end subroutine

   !--------------------------------------------------------------------
   ! Face fluxes
   !--------------------------------------------------------------------
   subroutine calc_fluxes(dt,U1,U2,flux,A)
   !Fluxes and flux Jacobians for every face. Outflow ghost states are
   !updated here and stored in U2 so they carry over to the next iteration.
   implicit none
   real*8, intent(in) :: dt
   real*8, intent(in) :: U1(:,:)
   real*8, intent(inout) :: U2(:,:)
   real*8, intent(out) :: flux(:,:),A(:,:,:,:)
   integer :: f,cl,cr,bc
   real*8 :: Up(2,4),np(2,2),nf(3)
    A(:,:,:,:) = 0.0d0
    do f = 1,nface
      cl = fcell(1,f)
      cr = fcell(2,f)
      Up(1,:) = U1(:,cl)
      Up(2,:) = U1(:,cr)
      nf = (/fnrm(1,f),fnrm(2,f),flen(f)/)
      if(f <= nintface) then
        call calc_MSWexp(flux(:,f),Up,nf,A(:,:,:,f))
      else
        bc = mbc(fmark(f))
        np(1,:) = (/cctr(1,cl),cctr(1,cr)/)
        np(2,:) = (/cctr(2,cl),cctr(2,cr)/)
        if(bc == BC_SLIPWALL .OR. bc == BC_ADIABATIC) then
          !Walls: pressure flux only (the adiabatic wall's inviscid part
          !is a slip wall); its Jacobian goes on the cell's diagonal
          call exp_bcs(dt,np,Up,flux(:,f),BC_SLIPWALL,nf)
          call wall_jac(Up(1,:),nf,A(:,:,1,f))
        else
          call exp_bcs(dt,np,Up,flux(:,f),bc,nf)
          call calc_MSWexp(flux(:,f),Up,nf,A(:,:,:,f))
          if(bc == BC_OUTFLOW) U2(:,cr) = Up(2,:)
        endif
      endif
    enddo
   end subroutine

   subroutine assemble_residual(flux,Rc)
   !Cell residuals from face fluxes
   implicit none
   real*8, intent(in) :: flux(:,:)
   real*8, intent(out) :: Rc(:,:)
   integer :: f,cl,cr
    Rc(:,:) = 0.0d0
    do f = 1,nface
      cl = fcell(1,f)
      cr = fcell(2,f)
      Rc(:,cl) = Rc(:,cl) + flux(:,f)
      if(cr <= ncell) Rc(:,cr) = Rc(:,cr) - flux(:,f)
    enddo
   end subroutine

   !--------------------------------------------------------------------
   ! Time integration
   !--------------------------------------------------------------------
   subroutine solver_exp(dt,U1,U2)
   !Explicit (forward Euler) step
   implicit none
   real*8, intent(in) :: dt
   real*8, intent(in) :: U1(:,:)
   real*8, intent(inout) :: U2(:,:)
   real*8, allocatable :: flux(:,:),A(:,:,:,:),Rc(:,:)
   integer :: c
    allocate(flux(4,nface),A(4,4,2,nface),Rc(4,ncell))
    call calc_fluxes(dt,U1,U2,flux,A)
    call assemble_residual(flux,Rc)
    do c = 1,ncell
      U2(:,c) = U1(:,c) - (dt/cvol(c))*Rc(:,c)
    enddo
    deallocate(flux,A,Rc)
   end subroutine

   subroutine solver_impl(dt,U1,U2)
   !Implicit (backward Euler) step, linearised and solved with kmax+1
   !point-Jacobi sweeps:
   !   (I + dt/V * D_c) dU_c = -dt/V * (R_c + sum over neighbours O_cn dU_n)
   !D_c (diagonal block) and O_cn (neighbour blocks) are built from the face
   !Jacobians. Wall Jacobians sit on the diagonal; other boundary faces see
   !their ghost-cell increment, set by impl_bcs from the interior increment
   !of the previous sweep.
   implicit none
   real*8, intent(in) :: dt
   real*8, intent(in) :: U1(:,:)
   real*8, intent(inout) :: U2(:,:)
   real*8, allocatable :: flux(:,:),A(:,:,:,:),Rc(:,:),Minv(:,:,:)
   real*8, allocatable :: dUk1(:,:),dUk2(:,:),B(:,:)
   real*8 :: M(4,4),Up(2,4),dUp(2,4),np(2,2),nf(3)
   integer :: c,f,k,cl,cr,bc
    allocate(flux(4,nface),A(4,4,2,nface),Rc(4,ncell),Minv(4,4,ncell))
    allocate(dUk1(4,ncell+nbface),dUk2(4,ncell),B(4,ncell))
    call calc_fluxes(dt,U1,U2,flux,A)
    call assemble_residual(flux,Rc)

    !Diagonal blocks: +A(left) on the left cell, -A(right) on the right
    !cell (boundary faces only have a left cell)
    Minv(:,:,:) = 0.0d0
    do f = 1,nface
      cl = fcell(1,f)
      cr = fcell(2,f)
      Minv(:,:,cl) = Minv(:,:,cl) + A(:,:,1,f)
      if(cr <= ncell) Minv(:,:,cr) = Minv(:,:,cr) - A(:,:,2,f)
    enddo
    do c = 1,ncell
      M = Id + (dt/cvol(c))*Minv(:,:,c)
      call inverse_4x4(M,Minv(:,:,c))
    enddo

    !Point-Jacobi sweeps
    dUk1(:,:) = 0.0d0
    do k = 1,kmax+1
      B(:,:) = 0.0d0
      if(k > 1) then
        !Ghost increments from the interior increments
        do f = nintface+1,nface
          bc = mbc(fmark(f))
          if(bc == BC_SLIPWALL .OR. bc == BC_ADIABATIC) cycle
          cl = fcell(1,f)
          cr = fcell(2,f)
          Up(1,:) = U1(:,cl)
          Up(2,:) = U1(:,cr)
          dUp(1,:) = dUk1(:,cl)
          dUp(2,:) = dUk1(:,cr)
          np(1,:) = (/cctr(1,cl),cctr(1,cr)/)
          np(2,:) = (/cctr(2,cl),cctr(2,cr)/)
          nf = (/fnrm(1,f),fnrm(2,f),flen(f)/)
          call impl_bcs(dt,np,Up,dUp,bc,nf)
          dUk1(:,cr) = dUp(2,:)
        enddo
        !Neighbour terms
        do f = 1,nintface
          cl = fcell(1,f)
          cr = fcell(2,f)
          B(:,cl) = B(:,cl) + matmul(A(:,:,2,f),dUk1(:,cr))
          B(:,cr) = B(:,cr) - matmul(A(:,:,1,f),dUk1(:,cl))
        enddo
        do f = nintface+1,nface
          bc = mbc(fmark(f))
          if(bc == BC_SLIPWALL .OR. bc == BC_ADIABATIC) cycle
          cl = fcell(1,f)
          B(:,cl) = B(:,cl) + matmul(A(:,:,2,f),dUk1(:,fcell(2,f)))
        enddo
      endif
      do c = 1,ncell
        dUk2(:,c) = matmul(Minv(:,:,c), -(dt/cvol(c))*(Rc(:,c) + B(:,c)))
      enddo
      dUk1(:,1:ncell) = dUk2(:,:)
    enddo

    do c = 1,ncell
      U2(:,c) = U1(:,c) + dUk1(:,c)
    enddo
    deallocate(flux,A,Rc,Minv,dUk1,dUk2,B)
   end subroutine

   subroutine set_ghosts(U)
   !Ghost-cell states from the boundary conditions, so every ghost holds a
   !sensible state (used as a neighbour by the cell gradients):
   !  inflow:    freestream (unchanged)
   !  outflow:   saved characteristic state (unchanged)
   !  slip wall: interior state with the normal velocity reflected
   !  mirror:    copy of the interior
   !  adiabatic: interior state with both velocity components reflected
   implicit none
   real*8, intent(inout) :: U(:,:)
   integer :: f,cl,g,bc
   real*8 :: V(4),Ui(4),un
    do f = nintface+1,nface
      bc = mbc(fmark(f))
      cl = fcell(1,f)
      g = fcell(2,f)
      select case(bc)
      case(BC_SLIPWALL)
        Ui = U(:,cl)
        call dfromp(Ui,V,-2)
        un = V(2)*fnrm(1,f) + V(3)*fnrm(2,f)
        V(2) = V(2) - 2.0d0*un*fnrm(1,f)
        V(3) = V(3) - 2.0d0*un*fnrm(2,f)
        call dfromp(Ui,V,2)
        U(:,g) = Ui
      case(BC_MIRROR)
        U(:,g) = U(:,cl)
      case(BC_ADIABATIC)
        U(:,g) = (/U(1,cl),-U(2,cl),-U(3,cl),U(4,cl)/)
      end select
    enddo
   end subroutine

   !--------------------------------------------------------------------
   ! Boundary conditions. U(1,:) is the interior cell, U(2,:) the ghost;
   ! nf = (n_x, n_y, length) with n pointing out of the domain.
   !--------------------------------------------------------------------
   subroutine exp_bcs(dt,np,U,F,bc,nf)
   !Sets the ghost state U(2,:) (or, for a slip wall, the face flux F)
   implicit none
   integer, intent(in) :: bc
   real*8, intent(in) :: dt,np(2,2),nf(3)
   real*8, intent(inout) :: U(:,:),F(:)
   real*8 :: V(4),V1(4),V2(4),a,dxave,ch(3),r1,r2,r3,r4
   real*8 :: nxt,nyt,un1,ut1,un2,ut2,dpc,dun,un,ut,Ui(4)
    select case(bc)
    case(BC_INFLOW)
      !Supersonic inflow: the ghost keeps its freestream state

    case(BC_OUTFLOW)
      !Characteristic outflow
      Ui = U(1,:)
      call dfromp(Ui,V1,-2)
      Ui = U(2,:)
      call dfromp(Ui,V2,-2)
      !Distance between the two cell centres
      dxave = sqrt((np(1,2)-np(1,1))**2 + (np(2,2)-np(2,1))**2)
      a = sqrt(gam*V1(4)/V1(1))
      !Velocities in the (normal, tangential) frame of the face
      nxt = nf(1)
      nyt = nf(2)
      un1 =  V1(2)*nxt + V1(3)*nyt
      ut1 = -V1(2)*nyt + V1(3)*nxt
      un2 =  V2(2)*nxt + V2(3)*nyt
      ut2 = -V2(2)*nyt + V2(3)*nxt
      !Characteristic CFL numbers. Only outgoing waves (>0) are updated
      !from the interior, and c/(1+c) (implicit in the ghost cell) keeps
      !the update bounded at large CFL
      ch(1) = max(un1*dt/dxave,0.0d0)
      ch(2) = max((un1+a)*dt/dxave,0.0d0)
      ch(3) = max((un1-a)*dt/dxave,0.0d0)
      ch = ch/(1.0d0+ch)
      !Ghost-cell update: w_g(new) = w_g - c*(w_g - w_interior)
      r1 = -ch(1)*(V2(1)-V1(1) - (V2(4)-V1(4))/(a*a))
      r2 = -ch(2)*(V2(4)-V1(4) + V1(1)*a*(un2-un1))
      r3 = -ch(1)*(ut2-ut1)
      r4 = -ch(3)*(V2(4)-V1(4) - V1(1)*a*(un2-un1))
      dpc = (r2+r4)/2.0d0
      dun = (r2-r4)/(2.0d0*V1(1)*a)
      un  = un2 + dun
      ut  = ut2 + r3
      V(1) = V2(1) + r1 + dpc/(a*a)
      V(2) = un*nxt - ut*nyt
      V(3) = un*nyt + ut*nxt
      V(4) = V2(4) + dpc
      call dfromp(Ui,V,2)
      U(2,:) = Ui

    case(BC_SLIPWALL)
      !Inviscid wall / symmetry plane: pressure flux only
      Ui = U(1,:)
      call dfromp(Ui,V,-2)
      F(:) = (/0.0d0,nf(1)*V(4)*nf(3),nf(2)*V(4)*nf(3),0.0d0/)

    case(BC_MIRROR)
      !Ghost = interior
      U(2,:) = U(1,:)

    case(BC_ADIABATIC)
      !No-slip adiabatic wall: velocity reflected, same density and energy
      U(2,:) = (/U(1,1),-U(1,2),-U(1,3),U(1,4)/)
    end select
   end subroutine

   subroutine impl_bcs(dt,np,U,dUk,bc,nf)
   !Ghost-cell increment dUk(2,:) from the interior increment dUk(1,:)
   !(walls are handled through wall_jac on the diagonal instead)
   implicit none
   integer, intent(in) :: bc
   real*8, intent(in) :: dt,np(2,2),nf(3)
   real*8, intent(inout) :: U(:,:),dUk(:,:)
   real*8 :: V1(4),V2(4),a,dxave,nxt,nyt,un1,vsq,gm1,Ui(4)
   real*8 :: C1(4,4),C1n(4,4),Rp(4,4),lam(4,4),S(4,4)
   real*8 :: q1(4,4),q2(4,4),q3(4,4),q4(4,4),q5(4,4),D1(4,4),D1inv(4,4)
   real*8 :: dV(4),RHSn(4),RHS(4),D2(4)
   integer :: dims1(4),dims2(4)
    dims1 = (/4,4,4,4/)
    dims2 = (/4,4,4,1/)
    select case(bc)
    case(BC_INFLOW)
      !Ghost state is fixed: no increment

    case(BC_OUTFLOW)
      Ui = U(1,:)
      call dfromp(Ui,V1,-2)
      Ui = U(2,:)
      call dfromp(Ui,V2,-2)
      dxave = sqrt((np(1,2)-np(1,1))**2 + (np(2,2)-np(2,1))**2)
      a = sqrt(gam*V1(4)/V1(1))
      nxt = nf(1)
      nyt = nf(2)
      un1 = V1(2)*nxt + V1(3)*nyt
      !Rotation of primitive variables into the (normal, tangential) frame
      Rp(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
      Rp(2,:) = (/0.0d0,nxt,nyt,0.0d0/)
      Rp(3,:) = (/0.0d0,-nyt,nxt,0.0d0/)
      Rp(4,:) = (/0.0d0,0.0d0,0.0d0,1.0d0/)
      C1n(1,:) = (/1.0d0,0.0d0,0.0d0,-1.0d0/(a*a)/)
      C1n(2,:) = (/0.0d0,V1(1)*a,0.0d0,1.0d0/)
      C1n(3,:) = (/0.0d0,0.0d0,1.0d0,0.0d0/)
      C1n(4,:) = (/0.0d0,-V1(1)*a,0.0d0,1.0d0/)
      call mat_mult(C1n,Rp,C1,dims1)
      !Only outgoing characteristics (lambda > 0) come from the interior
      lam(:,:) = 0.0d0
      lam(1,1) = max(un1,0.0d0)
      lam(2,2) = max(un1+a,0.0d0)
      lam(3,3) = max(un1,0.0d0)
      lam(4,4) = max(un1-a,0.0d0)
      gm1 = gam-1.0d0
      vsq = 0.5d0*(V1(2)**2 + V1(3)**2)
      S(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
      S(2,:) = (/-V1(2)/V1(1),1.0d0/V1(1),0.0d0,0.0d0/)
      S(3,:) = (/-V1(3)/V1(1),0.0d0,1.0d0/V1(1),0.0d0/)
      S(4,:) = (/vsq*gm1,-V1(2)*gm1,-V1(3)*gm1,gm1/)
      q1 = (dt/dxave)*lam
      call mat_mult(q1,C1,q2,dims1)
      dV = V2 - V1
      call mat_mult(q2,dV,RHSn,dims2)
      !contribution from the interior cell's increment
      call mat_mult(q2,S,q3,dims1)
      D2 = matmul(q3,dUk(1,:))
      !solve for the ghost cell's increment
      q4 = Id + q1
      call mat_mult(q4,C1,q5,dims1)
      call mat_mult(q5,S,D1,dims1)
      call inverse_4x4(D1,D1inv)
      RHS = D2 - RHSn
      dUk(2,:) = matmul(D1inv,RHS)

    case(BC_MIRROR)
      dUk(2,:) = dUk(1,:)
    end select
   end subroutine

   subroutine wall_jac(Ui,nf,W)
   !Jacobian of the slip-wall flux (0, n_x*p*L, n_y*p*L, 0) with respect
   !to the conserved variables of the interior cell
   implicit none
   real*8, intent(in) :: Ui(:),nf(3)
   real*8, intent(inout) :: W(:,:)
   real*8 :: Ut(4),V(4),dpdU(4)
    Ut = Ui
    call dfromp(Ut,V,-2)
    dpdU = gamm*(/0.5d0*(V(2)*V(2)+V(3)*V(3)),-V(2),-V(3),1.0d0/)
    W(:,:) = 0.0d0
    W(2,:) = nf(1)*nf(3)*dpdU
    W(3,:) = nf(2)*nf(3)*dpdU
   end subroutine

   !--------------------------------------------------------------------
   ! Modified Steger-Warming flux
   !--------------------------------------------------------------------
   subroutine blend_vars(VL,VR,Vave,Vp,Vm)
   !Pressure-gradient weighted blending of the left/right states
   implicit none
   real*8, intent(in) :: VL(:),VR(:)
   real*8, intent(inout) :: Vave(:),Vp(:),Vm(:)
   real*8 :: pg,den,AL,AR,wt
    !minimum pressure of the two sides
    den = min(VL(4),VR(4))
    !pressure gradient across the face and weight
    pg = (VR(4)-VL(4))/den
    wt = 1.0d0/(1.0d0+pg**2)
    AL = sqrt(gam*VL(4)/VL(1))
    AR = sqrt(gam*VR(4)/VR(1))
    Vp(1) = Vave(1)
    Vp(2) = wt*((VL(2)+VR(2))/2.0d0)+(1.0d0-wt)*VL(2)
    Vp(3) = wt*((VL(3)+VR(3))/2.0d0)+(1.0d0-wt)*VL(3)
    Vp(4) = Vave(4)
    Vp(5) = wt*((AL+AR)/2.0d0)+(1.0d0-wt)*AL
    Vm(1) = Vave(1)
    Vm(2) = wt*((VL(2)+VR(2))/2.0d0)+(1.0d0-wt)*VR(2)
    Vm(3) = wt*((VL(3)+VR(3))/2.0d0)+(1.0d0-wt)*VR(3)
    Vm(4) = Vave(4)
    Vm(5) = wt*((AL+AR)/2.0d0)+(1.0d0-wt)*AR
   end subroutine

   subroutine calc_MSWexp(F,U,nf,Af)
   !MSW flux F = A+ U(1) + A- U(2) through a face with normal nf (n_x, n_y,
   !length). Af(:,:,1) = A+, Af(:,:,2) = A-.
   implicit none
   real*8, intent(in) :: nf(3)
   real*8, intent(inout) :: U(:,:),F(:),Af(:,:,:)
   integer :: i,k,dims(4),ipass,npass
   real*8 :: VL(4),VR(4),Vave(4),Uave(4),Vp(5),Vm(5),Ut(4)
   real*8 :: a,a2,vsq,gm1,eps,da,cosa,sina,sgn,up
   real*8 :: S(4,4),Ca(4,4),Ra(4,4)
   real*8 :: Sinv(4,4),Cainv(4,4),Rainv(4,4)
   real*8 :: Lambda(4,4,2)
   real*8 :: m1(4,4),m2(4,4),m3(4,4),m4(4,4),m5(4,4)
   real*8 :: F1(4),F2(4)
    Ut = U(1,:)
    call dfromp(Ut,VL,-2)
    Ut = U(2,:)
    call dfromp(Ut,VR,-2)
    Uave(:) = (U(1,:) + U(2,:))/2.0d0
    call dfromp(Uave,Vave,-2)
    gm1 = gam-1.0d0
    call blend_vars(VL,VR,Vave,Vp,Vm)
    cosa = nf(1)
    sina = nf(2)
    da = nf(3)

    !Pass 1: MSW split (pressure-weighted blending) -> face flux.
    !Pass 2 (implicit only): Jacobians for the implicit operator from an
    !unblended split (A+ at the left state, A- at the right state, i.e.
    !Steger-Warming-like). The blended MSW Jacobians are not consistent
    !with the flux at shocks, which makes the implicit scheme unstable at
    !high CFL; this split is more dissipative.
    npass = 1
    if(impl == 1) npass = 2
    do ipass = 1,npass
      if(ipass == 2) then
        Vp(2) = VL(2)
        Vp(3) = VL(3)
        Vp(5) = sqrt(gam*VL(4)/VL(1))
        Vm(2) = VR(2)
        Vm(3) = VR(3)
        Vm(5) = sqrt(gam*VR(4)/VR(1))
      endif
      Lambda(:,:,:) = 0.0d0
      do k = 1,2
        if(k == 1) then
          if(iblend == 1) then
            Vave(2) = Vp(2)
            Vave(3) = Vp(3)
            a = Vp(5)
          else
            a = sqrt(gam*Vave(4)/Vave(1))
          endif
          sgn = 1.0d0
        else
          if(iblend == 1) then
            Vave(2) = Vm(2)
            Vave(3) = Vm(3)
            a = Vm(5)
          else
            a = sqrt(gam*Vave(4)/Vave(1))
          endif
          sgn = -1.0d0
        endif

        !Kinetic energy from the (possibly blended) velocities so that S
        !and Sinv stay exact inverses
        vsq = 0.5d0*(Vave(2)**2 + Vave(3)**2)
        a2 = a*a
        eps = 0.001d0*a
        Sinv(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
        Sinv(2,:) = (/Vave(2),Vave(1),0.0d0,0.0d0/)
        Sinv(3,:) = (/Vave(3),0.0d0,Vave(1),0.0d0/)
        Sinv(4,:) = (/vsq,Vave(1)*Vave(2),Vave(1)*Vave(3),1.0d0/gm1/)
        S(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
        S(2,:) = (/-Vave(2)/Vave(1),1.0d0/Vave(1),0.0d0,0.0d0/)
        S(3,:) = (/-Vave(3)/Vave(1),0.0d0,1.0d0/Vave(1),0.0d0/)
        S(4,:) = (/vsq*gm1,-Vave(2)*gm1,-Vave(3)*gm1,gm1/)
        Cainv(1,:) = (/1.0d0,1.0d0/(2.0d0*a2),0.0d0,1.0d0/(2.0d0*a2)/)
        Cainv(2,:) = (/0.0d0,1.0d0/(2.0d0*Vave(1)*a),0.0d0,-1.0d0/(2.0d0*Vave(1)*a)/)
        Cainv(3,:) = (/0.0d0,0.0d0,1.0d0,0.0d0/)
        Cainv(4,:) = (/0.0d0,0.5d0,0.0d0,0.5d0/)
        Ca(1,:) = (/1.0d0,0.0d0,0.0d0,-1.0d0/a2/)
        Ca(2,:) = (/0.0d0,Vave(1)*a,0.0d0,1.0d0/)
        Ca(3,:) = (/0.0d0,0.0d0,1.0d0,0.0d0/)
        Ca(4,:) = (/0.0d0,-Vave(1)*a,0.0d0,1.0d0/)
        !Rotation into the face frame
        up = cosa*Vave(2) + sina*Vave(3)
        Rainv(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
        Rainv(2,:) = (/0.0d0,cosa,-sina,0.0d0/)
        Rainv(3,:) = (/0.0d0,sina,cosa,0.0d0/)
        Rainv(4,:) = (/0.0d0,0.0d0,0.0d0,1.0d0/)
        Ra(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
        Ra(2,:) = (/0.0d0,cosa,sina,0.0d0/)
        Ra(3,:) = (/0.0d0,-sina,cosa,0.0d0/)
        Ra(4,:) = (/0.0d0,0.0d0,0.0d0,1.0d0/)
        !Eigenvalues (times face length)
        Lambda(1,:,1) = (/da*up,0.0d0,0.0d0,0.0d0/)
        Lambda(2,:,1) = (/0.0d0,da*(up+a),0.0d0,0.0d0/)
        Lambda(3,:,1) = (/0.0d0,0.0d0,da*up,0.0d0/)
        Lambda(4,:,1) = (/0.0d0,0.0d0,0.0d0,da*(up-a)/)
        !forward (k=1) and backward (k=2) characteristics
        do i = 1,4
          Lambda(i,i,k) = (Lambda(i,i,1) + sgn*sqrt(Lambda(i,i,1)**2 + eps**2))/2.0d0
        enddo
        !A+ and A-
        dims = (/4,4,4,4/)
        call mat_mult(Sinv,Rainv,m1,dims)
        call mat_mult(m1,Cainv,m2,dims)
        call mat_mult(m2,Lambda(:,:,k),m3,dims)
        call mat_mult(m3,Ca,m4,dims)
        call mat_mult(m4,Ra,m5,dims)
        call mat_mult(m5,S,Af(:,:,k),dims)
      enddo

      if(ipass == 1) then
        !Face flux (always from the MSW split)
        F1 = U(1,:)
        F2 = U(2,:)
        F = matmul(Af(:,:,1),F1) + matmul(Af(:,:,2),F2)
      endif
    enddo
   end subroutine

end module
