        module solver
        use utilities
        implicit none


        contains
          subroutine solver_exp(dt,U1,U2)
          implicit none
          real*8, intent(inout) :: dt, U1(:,:,:), U2(:,:,:)
          integer*8 :: i,j,k,gid,fi,bcin
          real*8 :: U(2,4),np(2,2),V(4)
          real*8 :: f1(4),f2(4),f3(4),f4(4)
          real*8, allocatable :: flux(:,:),A(:,:,:,:),L(:,:,:),vf(:,:)
        
          !Looping over all faces
          allocate(flux(fnum,4),A(fnum,4,4,2),L(fnum,4,4),vf(fnum,4)) 
          call calc_fluxes(dt,U1,U2,A,flux)
          if(ivisc == 1) then
            call calc_visc(dt,U1,U2,L,vf) 
          else
            vf(:,:) = 0.0d0
          endif
          !Looping over all cells
          do i = 2,nx
          do j = 2,ny
             f1 = -flux(ncf(i,j,1),:)
             f2 =  flux(ncf(i,j,2),:)
             f3 = -(flux(ncf(i,j,3),:)-vf(ncf(i,j,3),:))
             f4 =  (flux(ncf(i,j,4),:)-vf(ncf(i,j,4),:))
             !Solving for conserved variables with fluxes
             U2(i,j,:) = U1(i,j,:) - (dt/dxy(i,j))*((f1+f2) + (f3+f4))
          enddo
          enddo
                    
          end subroutine
        
          subroutine solver_impl(dt,U1,U2)
          implicit none
          real*8, intent(inout) :: dt, U1(:,:,:), U2(:,:,:)
          integer*8 :: i,j,k,gid,fi,dims(4),jj,ii,bcin
          real*8 :: U(2,4),np(2,2),V(4),M(4,4)
          real*8 :: A1(4,4),A2(4,4),A3(4,4),A4(4,4),A5(4,4),A6(4,4)
          real*8 :: f1(4),f2(4),f3(4),f4(4),RHS(4),dUk(2,4)
          real*8 :: D(4),t1(4),t2(4),t3(4),t4(4),dUkt(4),t5(4),t6(4)
          real*8, allocatable :: A(:,:,:,:),L(:,:,:)
          real*8, allocatable :: flux(:,:),Rij(:,:,:),vf(:,:)
          real*8, allocatable :: Minv(:,:,:,:),dUk1(:,:,:),dUk2(:,:,:)

          !Looping over all faces
          allocate(flux(fnum,4),vf(fnum,4),L(fnum,4,4))
          allocate(Rij(nx+1,ny+1,4),A(fnum,4,4,2),Minv(nx+1,ny+1,4,4))
          allocate(dUk1(nx+1,ny+1,4),dUk2(nx+1,ny+1,4))
          call calc_fluxes(dt,U1,U2,A,flux)
          if(ivisc==1) then
             call calc_visc(dt,U1,U2,L,vf) 
          else
            vf(:,:) = 0.0d0
            L(:,:,:) = 0.0d0
          endif
 
          !Looping over all cells to solve for exp RHS
          do i = 2,nx
          do j = 2,ny
             f1 = -flux(ncf(i,j,1),:)
             f2 =  flux(ncf(i,j,2),:)
             f3 = -(flux(ncf(i,j,3),:)-vf(ncf(i,j,3),:))
             f4 =  (flux(ncf(i,j,4),:)-vf(ncf(i,j,4),:))
             !Solving for conserved variables with fluxes
             Rij(i,j,:) = (dt/dxy(i,j))*((f1+f2) + (f3+f4))

             !Calculating LHS matrix (M)
             A1 = A(ncf(i,j,1),:,:,2)
             A2 = A(ncf(i,j,2),:,:,1)
             A3 = A(ncf(i,j,3),:,:,2)
             A4 = A(ncf(i,j,4),:,:,1)
             if(ivisc==1) then
                !Viscous diagonal: bottom + top face Jacobians. At wall and
                !mirror faces the ghost state depends on this cell, so that
                !dependence is folded into M (see visc_diag)
                call visc_diag(L(ncf(i,j,3),:,:),bclist(ncf(i,j,3)),A5)
                call visc_diag(L(ncf(i,j,4),:,:),bclist(ncf(i,j,4)),A6)
                A5 = A5 + A6
             else
               A5(:,:) = 0.0d0
             endif
             M(:,:) = Id + (dt/dxy(i,j))*(-A1+A2-A3+A4+A5)
             call inverse_4x4(M(:,:),Minv(i,j,:,:))
          enddo
          enddo

          dUk1(:,:,:) = 0.0d0
          dUk2(:,:,:) = 0.0d0
          !Relation steps
          do k=1,kmax+1
            do i = 2,nx
            do j = 2,ny
             !Jacobians
             A1 = A(ncf(i,j,1),:,:,1)
             A2 = A(ncf(i,j,2),:,:,2)
             A3 = A(ncf(i,j,3),:,:,1)
             A4 = A(ncf(i,j,4),:,:,2)
             A5 = L(ncf(i,j,3),:,:)   !bottom face viscous Jacobian
             A6 = L(ncf(i,j,4),:,:)   !top face viscous Jacobian
            
             !F-1/2
             dims = (/4,4,4,1/)
             if(bclist(ncf(i,j,1)).NE.0 .AND. i==2) then
               fi = ncf(i,j,1)
               np(1,:) = (/cxy(i,j,1),cxy(i-1,j,1)/)
               np(2,:) = (/cxy(i,j,2),cxy(i-1,j,2)/)
               U(1,:) = U1(i-1,j,:)
               U(2,:) = U1(i,j,:)
               dUk(1,:) = dUk1(i-1,j,:)
               dUk(2,:) = dUk1(i,j,:)
               if(bclist(fi) == 3) then
                 t1(:) = 0.0d0   !wall: Jacobian is already in M
               else
                 call impl_bcs(dt,np,U,dUk,flux(fi,:),bclist(fi),1,fi)
                 call mat_mult(A1,dUk(1,:),t1,dims)
               endif
             else
               call mat_mult(A1,dUk1(i-1,j,:),t1,dims)
             endif

             !F+1/2
             if(bclist(ncf(i,j,2)).NE.0 .AND.i==nx) then
               fi = ncf(i,j,2)
               np(1,:) = (/cxy(i,j,1),cxy(i+1,j,1)/)
               np(2,:) = (/cxy(i,j,2),cxy(i+1,j,2)/)
               U(1,:) = U1(i,j,:)
               U(2,:) = U1(i+1,j,:)
               dUk(1,:) = dUk1(i,j,:)
               dUk(2,:) = dUk1(i+1,j,:)
               if(bclist(fi) == 3) then
                 t2(:) = 0.0d0   !wall: Jacobian is already in M
               else
                 call impl_bcs(dt,np,U,dUk,flux(fi,:),bclist(fi),2,fi)
                 call mat_mult(A2,dUk(2,:),t2,dims)
               endif
             else
               call mat_mult(A2,dUk1(i+1,j,:),t2,dims)
             endif

             !G-1/2
             if(bclist(ncf(i,j,3)).NE.0 .AND. j==2) then
               fi = ncf(i,j,3)
               np(1,:) = (/cxy(i,j,1),cxy(i,j-1,1)/)
               np(2,:) = (/cxy(i,j,2),cxy(i,j-1,2)/)
               U(1,:) = U1(i,j-1,:)
               U(2,:) = U1(i,j,:)
               dUk(1,:) = dUk1(i,j-1,:)
               dUk(2,:) = dUk1(i,j,:)
               if(bclist(fi) == 3 .OR. bclist(fi) == 5) then
                !Wall (inviscid part of the adiabatic wall is a slip wall):
                !its pressure-flux Jacobian is already in M
                t3(:) = 0.0d0
               else
                call impl_bcs(dt,np,U,dUk,flux(fi,:),bclist(fi),1,fi)
                call mat_mult(A3,dUk(1,:),t3,dims)
               endif
             else
               call mat_mult(A3,dUk1(i,j-1,:),t3,dims)
             endif

             !G+1/2
             if(bclist(ncf(i,j,4)).NE.0 .AND. j==ny) then
               fi = ncf(i,j,4)
               np(1,:) = (/cxy(i,j,1),cxy(i,j+1,1)/)
               np(2,:) = (/cxy(i,j,2),cxy(i,j+1,2)/)
               U(1,:) = U1(i,j,:)
               U(2,:) = U1(i,j+1,:)
               dUk(1,:) = dUk1(i,j,:)
               dUk(2,:) = dUk1(i,j+1,:)
               if(bclist(fi) == 3) then
                 t4(:) = 0.0d0   !wall: Jacobian is already in M
               else
                 call impl_bcs(dt,np,U,dUk,flux(fi,:),bclist(fi),2,fi)
                 call mat_mult(A4,dUk(2,:),t4,dims)
               endif
             else
                call mat_mult(A4,dUk1(i,j+1,:),t4,dims)
             endif

             !Visc G-1/2
             if(ivisc == 1) then
              if(bclist(ncf(i,j,3)).NE.0 .AND. j==2) then
                fi = ncf(i,j,3)
                np(1,:) = (/cxy(i,j,1),cxy(i,j-1,1)/)
                np(2,:) = (/cxy(i,j,2),cxy(i,j-1,2)/)
                U(1,:) = U1(i,j-1,:)
                U(2,:) = U1(i,j,:)
                dUk(1,:) = dUk1(i,j-1,:)
                dUk(2,:) = dUk1(i,j,:)
                if(bclist(fi) >= 3) then
                  !wall / mirror: ghost depends only on this cell, which
                  !is already accounted for in M (visc_diag)
                  t5(:) = 0.0d0
                else
                  call impl_bcs(dt,np,U,dUk,flux(fi,:),bclist(fi),1,fi)
                  call mat_mult(A5,dUk(1,:),t5,dims)
                endif
              else
                call mat_mult(A5,dUk1(i,j-1,:),t5,dims)
              endif

              !Visc G+1/2
              if(bclist(ncf(i,j,4)).NE.0 .AND. j==ny) then
                fi = ncf(i,j,4)
                np(1,:) = (/cxy(i,j,1),cxy(i,j+1,1)/)
                np(2,:) = (/cxy(i,j,2),cxy(i,j+1,2)/)
                U(1,:) = U1(i,j,:)
                U(2,:) = U1(i,j+1,:)
                dUk(1,:) = dUk1(i,j,:)
                dUk(2,:) = dUk1(i,j+1,:)
                if(bclist(fi) >= 3) then
                  t6(:) = 0.0d0   !wall / mirror: handled in M
                else
                  call impl_bcs(dt,np,U,dUk,flux(fi,:),bclist(fi),2,fi)
                  call mat_mult(A6,dUk(2,:),t6,dims)
                endif
              else
                 call mat_mult(A6,dUk1(i,j+1,:),t6,dims)
              endif
             endif


             !K-1 step inward facing fluxes
             if(k==1) then
               t1(:) = 0.0d0
               t2(:) = 0.0d0
               t3(:) = 0.0d0
               t4(:) = 0.0d0
               t5(:) = 0.0d0
               t6(:) = 0.0d0
             endif
             !Inviscid off-diagonals: +dF(minus faces) - dF(plus faces)
             !Viscous off-diagonals: +L_b*dU(j-1) + L_t*dU(j+1)
             D = (dt/dxy(i,j))*(t1-t2+t3-t4+t5+t6)

             !compute dUk
             RHS = -Rij(i,j,:) + D
             call mat_mult(Minv(i,j,:,:),RHS,dUk2(i,j,:),dims)

             if(k==kmax+1) then
               U2(i,j,:) = U1(i,j,:) + dUk2(i,j,:)
             endif

            enddo
            enddo
            dUk1 = dUk2
          enddo
          deallocate(Rij,dUk1,dUk2,A,Minv)

          end subroutine

          subroutine calc_fluxes(dt,U1,U2,A,flux)
          implicit none
          real*8, intent(inout) :: dt, U1(:,:,:), U2(:,:,:)
          integer*8 :: i,j,k,gid,fi,bcin
          real*8 :: U(2,4),np(2,2),V(4)
          real*8 :: f1(4),f2(4),f3(4),f4(4)
          real*8, intent(inout) :: flux(:,:),A(:,:,:,:)
          !Flux Jacobians: A(:,:,:,1) wrt U(1) side, A(:,:,:,2) wrt U(2) side.
          !Wall faces only get the pressure-flux Jacobian (wall_jac)
          A(:,:,:,:) = 0.0d0
          fi = 1
          do j = 1,ny-1
          do i = 1,nx
             !F face between cells (i,j+1) and (i+1,j+1)
             !np(1,:) = x of both cell centers, np(2,:) = y
             np(1,:) = (/cxy(i,j+1,1),cxy(i+1,j+1,1)/)
             np(2,:) = (/cxy(i,j+1,2),cxy(i+1,j+1,2)/)
             U(1,:) = U1(nfc(fi,1,1),nfc(fi,1,2),:)
             U(2,:) = U1(nfc(fi,2,1),nfc(fi,2,2),:)
             if(bclist(fi).NE.0 .AND. i==nx) then
               call exp_bcs(dt,np,U,flux(fi,:),bclist(fi),2,fi)
               if(bclist(fi) .NE. 3) then
                call calc_MSWexp(flux(fi,:),U,fi,A(fi,:,:,:))
               else
                call wall_jac(U(1,:),fi,A(fi,:,:,1))
               endif
               !Persist the updated outflow ghost state
               if(bclist(fi) == 2) U2(nfc(fi,2,1),nfc(fi,2,2),:) = U(2,:)
             elseif(bclist(fi).NE.0 .AND. i==1) then
               call exp_bcs(dt,np,U,flux(fi,:),bclist(fi),1,fi)
               if(bclist(fi) .NE. 3) then
                call calc_MSWexp(flux(fi,:),U,fi,A(fi,:,:,:))
               else
                call wall_jac(U(2,:),fi,A(fi,:,:,2))
               endif
               if(bclist(fi) == 2) U2(nfc(fi,1,1),nfc(fi,1,2),:) = U(1,:)
             else
                call calc_MSWexp(flux(fi,:),U,fi,A(fi,:,:,:))
             endif

             fi = fi + 1
          enddo
          enddo

          do i = 1,nx-1
          do j = 1,ny
             !G face between cells (i+1,j) and (i+1,j+1)
             np(1,:) = (/cxy(i+1,j,1),cxy(i+1,j+1,1)/)
             np(2,:) = (/cxy(i+1,j,2),cxy(i+1,j+1,2)/)
             U(1,:) = U1(nfc(fi,1,1),nfc(fi,1,2),:)
             U(2,:) = U1(nfc(fi,2,1),nfc(fi,2,2),:)
             if(bclist(fi).NE.0 .AND. j==ny) then
              call exp_bcs(dt,np,U,flux(fi,:),bclist(fi),2,fi)
              if(bclist(fi) .NE. 3) then
               call calc_MSWexp(flux(fi,:),U,fi,A(fi,:,:,:))
              else
               call wall_jac(U(1,:),fi,A(fi,:,:,1))
              endif
              if(bclist(fi) == 2) U2(nfc(fi,2,1),nfc(fi,2,2),:) = U(2,:)
             elseif(bclist(fi).NE.0 .AND. j==1) then
              if(bclist(fi) == 5) then
               bcin = 3
               call exp_bcs(dt,np,U,flux(fi,:),bcin,1,fi)
               call wall_jac(U(2,:),fi,A(fi,:,:,2))
              else
               call exp_bcs(dt,np,U,flux(fi,:),bclist(fi),1,fi)
               if(bclist(fi) .NE. 3) then
                call calc_MSWexp(flux(fi,:),U,fi,A(fi,:,:,:))
               else
                call wall_jac(U(2,:),fi,A(fi,:,:,2))
               endif
               if(bclist(fi) == 2) U2(nfc(fi,1,1),nfc(fi,1,2),:) = U(1,:)
              endif
             else
               call calc_MSWexp(flux(fi,:),U,fi,A(fi,:,:,:))
             endif
             fi = fi+1
          enddo
          enddo
          end subroutine

          subroutine wall_jac(Ui,dir,W)
          !Jacobian of the slip-wall flux (0, nx*p*S, ny*p*S, 0) with
          !respect to the conserved variables of the interior cell
          implicit none
          real*8, intent(in) :: Ui(:)
          integer*8, intent(in) :: dir
          real*8, intent(inout) :: W(:,:)
          real*8 :: Ut(4),V(4),dpdU(4)
           Ut = Ui
           call dfromp(Ut,V,-2)
           dpdU = gamm*(/0.5d0*(V(2)*V(2)+V(3)*V(3)),-V(2),-V(3),1.0d0/)
           W(:,:) = 0.0d0
           W(2,:) = ndf(dir,1)*ndf(dir,3)*dpdU
           W(3,:) = ndf(dir,2)*ndf(dir,3)*dpdU
          end subroutine

          subroutine visc_diag(Lf,bc,Ld)
          !Diagonal (own-cell) part of a face's viscous Jacobian.
          !Interior/inflow/outflow: Lf. Adiabatic wall (ghost = diag(1,-1,-1,1)
          !times the interior state): Lf*(I-Rw). Mirror (ghost = interior)
          !and slip wall: no viscous flux, so zero.
          implicit none
          real*8, intent(in) :: Lf(:,:)
          integer*8, intent(in) :: bc
          real*8, intent(inout) :: Ld(:,:)
           Ld(:,:) = 0.0d0
           if(bc == 5) then
             Ld(:,2) = 2.0d0*Lf(:,2)
             Ld(:,3) = 2.0d0*Lf(:,3)
           elseif(bc == 3 .OR. bc == 4) then
             Ld(:,:) = 0.0d0
           else
             Ld(:,:) = Lf(:,:)
           endif
          end subroutine

          subroutine calc_visc(dt,U1,U2,A,flux)
          implicit none
          real*8, intent(inout) :: dt, U1(:,:,:), U2(:,:,:)
          integer*8 :: i,j,k,gid,fi
          real*8 :: U(2,4),np(2,2),V(4)
          real*8 :: f1(4),f2(4),f3(4),f4(4)
          real*8, intent(inout) :: flux(:,:),A(:,:,:)

          !Only G (j-direction) faces carry a thin-layer viscous flux;
          !F faces and slip walls stay zero
          flux(:,:) = 0.0d0
          A(:,:,:) = 0.0d0
          fi = nx*(ny-1)+1
          do i = 1,nx-1
          do j = 1,ny
             !G face between cells (i+1,j) and (i+1,j+1)
             np(1,:) = (/cxy(i+1,j,1),cxy(i+1,j+1,1)/)
             np(2,:) = (/cxy(i+1,j,2),cxy(i+1,j+1,2)/)
             U(1,:) = U1(nfc(fi,1,1),nfc(fi,1,2),:)
             U(2,:) = U1(nfc(fi,2,1),nfc(fi,2,2),:)
             if(bclist(fi) == 3) then
              !Slip wall / symmetry plane: no shear stress, no heat flux
              flux(fi,:) = 0.0d0
             elseif(bclist(fi).NE.0 .AND. j==ny) then
              call exp_bcs(dt,np,U,flux(fi,:),bclist(fi),2,fi)
              call calc_vf(flux(fi,:),U,fi,A(fi,:,:),df(fi,:))
             elseif(bclist(fi).NE.0 .AND. j==1) then
              call exp_bcs(dt,np,U,flux(fi,:),bclist(fi),1,fi)
              call calc_vf(flux(fi,:),U,fi,A(fi,:,:),df(fi,:))
             else
               call calc_vf(flux(fi,:),U,fi,A(fi,:,:),df(fi,:))
             endif
             !Note: calc_vf already returns the total face flux (the
             !metric terms in Mnn include the face length), so no
             !extra multiplication by ndf(fi,3) here
             fi = fi+1
          enddo
          enddo
          end subroutine
        
          subroutine init_flow(U)
          implicit none
          real*8, intent(inout) :: U(:,:,:)
          real*8 :: a,ufs,V(4)
          integer*8 :: i,j,k
           a = sqrt(gam*fs(2)/fs(3))
           ufs = fs(1)*a
           do i = 1,nx+1
           do j = 1,ny+1
             if(switch == 1) then !initialize to freestream
              V(:) = (/fs(3),ufs,0.0d0,fs(2)/)
             elseif(switch == 2) then !initialize with normal shock
               if(j .LE. (ny-1)/2) then
                V(:) = (/2.0d0*fs(3),0.0d0,0.0d0,2.0d0*fs(2)/)
                !V(:) = (/fs(3),0.0d0,0.0d0,fs(2)/)
               else
                !V(:) = (/2.0d0*fs(3),0.0d0,0.0d0,2.0d0*fs(2)/)
                V(:) = (/fs(3),0.0d0,0.0d0,fs(2)/)
               endif
             elseif(switch == 3) then
               V(:) = (/fs(3),ufs,0.0d0,fs(2)/)
               if(i == (nx-1)/2) then
                V(:) = (/fs(3)+0.01*fs(3),ufs+0.01*ufs,0.0d0,fs(2)+0.01*fs(2)/)
               endif

             endif
             call dfromp(U(i,j,:),V(:),2)
           enddo
           enddo

          end subroutine

          subroutine impl_bcs(dt,np,U,dUk,F,bc,gid,dir)
          implicit none
          integer*8, intent(in) :: bc,dir
          integer*4, intent(in) :: gid
          real*8, intent(in) :: dt
          real*8, intent(inout) ::U(:,:),F(:),np(:,:),dUk(:,:)
          real*8 :: xf,yf,nxt,nyt,V(4),V1(4),V2(4),fdir
          real*8 :: r1,r2,r3,r4,ch(3),a,dxave,qs,beta,pb
          real*8 :: C1(4,4),lam(4,4),S(4,4),vsq,gm1
          real*8 :: q1(4,4),q2(4,4),q3(4,4),q4(4,4),q5(4,4)
          real*8 :: dV(4),RHSn(4),RHS(4),D2(4)
          real*8 :: D1(4,4),D1inv(4,4)
          real*8 :: dUk1(4),dUk2(4)
          real*8 :: C1n(4,4),Rp(4,4),un1
          integer*8 :: inid,dims1(4),dims2(4)
            dims1=(/4,4,4,4/)
            dims2=(/4,4,4,1/)
            if(gid == 1) then
               inid = 2
               fdir = -1.0d0
            else
               inid = 1
               fdir = 1.0d0
            endif
            !Supersonic inflow
            if(bc == 1) then
             U(gid,:) = U(gid,:)

            !Supersonic outflow
            elseif(bc == 2) then
             call dfromp(U(inid,:),V1,-2)
             call dfromp(U(gid,:),V2,-2)
             !Distance between the two cell centers
             dxave = sqrt((np(1,2)-np(1,1))**2 + (np(2,2)-np(2,1))**2)
             a = sqrt(gam*V1(4)/V1(1))
             !Outward unit normal of the boundary face
             nxt = fdir*ndf(dir,1)
             nyt = fdir*ndf(dir,2)
             un1 = V1(2)*nxt + V1(3)*nyt
             !Rotation of primitive vars into (normal, tangential) frame
             Rp(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
             Rp(2,:) = (/0.0d0,nxt,nyt,0.0d0/)
             Rp(3,:) = (/0.0d0,-nyt,nxt,0.0d0/)
             Rp(4,:) = (/0.0d0,0.0d0,0.0d0,1.0d0/)
             C1n(1,:) = (/1.0d0,0.0d0,0.0d0,-1.0d0/(a*a)/)
             C1n(2,:) = (/0.0d0,V1(1)*a,0.0d0,1.0d0/)
             C1n(3,:) = (/0.0d0,0.0d0,1.0d0,0.0d0/)
             C1n(4,:) = (/0.0d0,-V1(1)*a,0.0d0,1.0d0/)
             call mat_mult(C1n,Rp,C1,dims1)

             !Only outgoing characteristics (lambda > 0) are taken
             !from the interior; incoming ones are held fixed
             lam(:,:) = 0.0d0
             lam(1,1) = max(un1,0.0d0)
             lam(2,2) = max(un1+a,0.0d0)
             lam(3,3) = max(un1,0.0d0)
             lam(4,4) = max(un1-a,0.0d0)

             gm1 = gam-1.0d0
             vsq = 0.5d0*(V1(2)**2 + V1(3)**2) 
             s(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
             s(2,:) = (/-V1(2)/V1(1),1.0d0/V1(1),0.0d0,0.0d0/)
             s(3,:) = (/-V1(3)/V1(1),0.0d0,1.0d0/V1(1),0.0d0/)
             s(4,:) = (/vsq*gm1,-V1(2)*gm1,-V1(3)*gm1,gm1/)

              
             q1 = (dt/dxave)*lam(:,:)
             !calculating RHS
             call mat_mult(q1,C1,q2,dims1) 
             dV = V2(:) - V1(:)
             call mat_mult(q2,dV,RHSn,dims2)

             !contribution from the interior cell's dUk
             call mat_mult(q2,S,q3,dims1)
             call mat_mult(q3,dUk(inid,:),D2,dims2)

             !solving for the ghost cell's dUk
             q4 = Id + q1
             call mat_mult(q4,C1,q5,dims1)
             call mat_mult(q5,S,D1,dims1)
             call inverse_4x4(D1,D1inv)

             RHS = D2-RHSn
             call mat_mult(D1inv,RHS,dUk(gid,:),dims2)

            !Inviscid wall or symmetry surface
            elseif(bc == 3) then
             call dfromp(U(inid,:),V,-2)
             gm1 = gam-1.0d0
             vsq = 0.5d0*(V(2)**2 + V(3)**2)
             s(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
             s(2,:) = (/-V(2)/V(1),1.0d0/V(1),0.0d0,0.0d0/)
             s(3,:) = (/-V(3)/V(1),0.0d0,1.0d0/V(1),0.0d0/)
             s(4,:) = (/vsq*gm1,-V(2)*gm1,-V(3)*gm1,gm1/)
             !Pressure perturbation of the interior cell
             call mat_mult(s,dUk(inid,:),D2,dims2)
             
             F(:) = (/0.0d0,ndf(dir,1)*D2(4)*ndf(dir,3),ndf(dir,2)*D2(4)*ndf(dir,3),0.0d0/)
             U(gid,:) = (/0.0d0,0.0d0,0.0d0,0.0d0/) 

            !Mirror BC (for shock tube and testing)
            elseif(bc==4) then
              U(gid,:) = U(inid,:)
              dUk(gid,:) = dUk(inid,:)

            !Adiabatic no-slip wall: ghost increment is the reflection
            !of the interior increment
            elseif(bc==5) then
             dUk(gid,:) = (/dUk(inid,1),-dUk(inid,2),-dUk(inid,3),dUk(inid,4)/)
            endif
          end subroutine


          subroutine exp_bcs(dt,np,U,F,bc,gid,dir)
          implicit none
          integer*8, intent(in) :: bc,dir
          integer*4, intent(in) :: gid
          real*8, intent(in) :: dt
          real*8, intent(inout) ::U(:,:),F(:),np(:,:)
          real*8 :: xf,yf,nxt,nyt,V(4),V1(4),V2(4),Vw(4)
          real*8 :: Tin,r1,r2,r3,r4,ch(3),a,dxave,qs,beta,pb
          real*8 :: fdir,un1,ut1,un2,ut2,dpc,dun,un,ut
          integer*8 :: inid
            if(gid == 1) then
               inid = 2
               fdir = -1.0d0
            else
               inid = 1
               fdir = 1.0d0
            endif
            !Supersonic inflow
            if(bc == 1) then
             U(gid,:) = U(gid,:)

            !Supersonic outflow
            elseif(bc == 2) then
             call dfromp(U(inid,:),V1,-2)
             call dfromp(U(gid,:),V2,-2)
             !Distance between the two cell centers
             dxave = sqrt((np(1,2)-np(1,1))**2 + (np(2,2)-np(2,1))**2)
             a = sqrt(gam*V1(4)/V1(1))

             !Outward unit normal; velocities in (normal, tangential) frame
             nxt = fdir*ndf(dir,1)
             nyt = fdir*ndf(dir,2)
             un1 =  V1(2)*nxt + V1(3)*nyt
             ut1 = -V1(2)*nyt + V1(3)*nxt
             un2 =  V2(2)*nxt + V2(3)*nyt
             ut2 = -V2(2)*nyt + V2(3)*nxt

             !Characteristic CFL numbers. Only outgoing waves (>0) are
             !updated from the interior, and c/(1+c) (implicit in the
             !ghost cell) keeps the update bounded at large CFL
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
             call dfromp(U(gid,:),V,2)
            
            !Inviscid wall or symmetry surface
            elseif(bc == 3) then
             call dfromp(U(inid,:),V,-2)
             F(:) = (/0.0d0,ndf(dir,1)*V(4)*ndf(dir,3),ndf(dir,2)*V(4)*ndf(dir,3),0.0d0/)
             U(gid,:) = (/0.0d0,0.0d0,0.0d0,0.0d0/) 

            !Mirror BC (for shock tube and testing)
            elseif(bc==4) then
              U(gid,:) = U(inid,:)
           
            !Adiabatic wall 
            elseif(bc==5) then
             call dfromp(U(inid,:),V,-2)
             Tin = V(4)/(Rair*V(1))
             !Vw = (/V(1),-V(2),-V(3),V(4)/)
             !call dfromp(U(gid,:),Vw,2)
             U(gid,:) = (/V(1),-U(inid,2),-U(inid,3),U(inid,4)/)
            endif
          end subroutine

          subroutine blend_vars(VL,VR,Vave,Vp,Vm)
          implicit none
          real*8, intent(in) :: VL(:),VR(:)
          real*8, intent(inout) :: Vave(:),Vp(:),Vm(:)
          real*8 :: pg,den,AL,AR,wt
          !Finding min pressure b/w i,i+1
          if(VL(4) .LE. VR(4)) then
            den = VL(4)
          else
            den = VR(4)
          endif

          !calculating pressure gradient across face
          pg = (VR(4)-VL(4))/den
          
          !Calculating wt factor
          wt = 1.0d0/(1.0d0+pg**2)
     
          !Blending variables
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

          subroutine calc_MSWexp(F,U,dir,Af)
          implicit none
          integer*8, intent(in) :: dir
          real*8, intent(inout) :: U(:,:),F(:),Af(:,:,:)
          integer*8 :: i,j,k,dims(4)
          real*8 :: VL(4),VR(4),Vave(4),Uave(4),Vp(5),Vm(5)
          real*8 :: a,a2,vsq,gm1,eps,da,cosa,sina,sgn,up
          real*8 :: S(4,4),Ca(4,4),Ra(4,4)
          real*8 :: Sinv(4,4),Cainv(4,4),Rainv(4,4)
          real*8 :: Lambda(4,4,2)!,Af(2,4,4)
          real*8 :: m1(4,4),m2(4,4),m3(4,4),m4(4,4),m5(4,4)
          real*8 :: F1(4),F2(4),nv(2)
          integer*8 :: ipass,npass
           call dfromp(U(1,:),VL,-2)
           call dfromp(U(2,:),VR,-2)
           Uave(:) = (U(1,:) + U(2,:))/2.0d0
           call dfromp(Uave,Vave,-2)
           gm1 = gam-1.0d0
           call blend_vars(VL,VR,Vave,Vp,Vm)

           !Pass 1: MSW split (pressure-weighted blending) -> face flux.
           !Pass 2 (implicit only): Jacobians for the implicit operator from
           !an unblended split (A+ at the left state, A- at the right state,
           !i.e. Steger-Warming-like). The blended MSW Jacobians are not
           !consistent with the flux at shocks, which makes the implicit
           !scheme unstable at high CFL; this split is more dissipative.
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
             elseif(k==2) then
              if(iblend == 1) then
               Vave(2) = Vm(2)
               Vave(3) = Vm(3)
               a = Vm(5)
              else
               a = sqrt(gam*Vave(4)/Vave(1))
              endif

              sgn = -1.0d0
             endif

             !Kinetic energy from the (possibly blended) velocities so
             !that S and Sinv stay exact inverses
             vsq = 0.5d0*(Vave(2)**2 + Vave(3)**2)
             a2 = a*a
             eps = 0.001d0*a
             !Defining Sinv
             Sinv(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
             sinv(2,:) = (/Vave(2),Vave(1),0.0d0,0.0d0/)
             sinv(3,:) = (/Vave(3),0.0d0,Vave(1),0.0d0/)
             sinv(4,:) = (/vsq,Vave(1)*Vave(2),Vave(1)*Vave(3),1.0d0/gm1/)
             !Defining S
             s(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
             s(2,:) = (/-Vave(2)/Vave(1),1.0d0/Vave(1),0.0d0,0.0d0/)
             s(3,:) = (/-Vave(3)/Vave(1),0.0d0,1.0d0/Vave(1),0.0d0/)
             s(4,:) = (/vsq*gm1,-Vave(2)*gm1,-Vave(3)*gm1,gm1/)
             !Defining Cainv
             Cainv(1,:) = (/1.0d0,1.0d0/(2.0d0*a2),0.0d0,1.0d0/(2.0d0*a2)/)
             Cainv(2,:) = (/0.0d0,1.0d0/(2.0d0*Vave(1)*a),0.0d0,-1.0d0/(2.0d0*Vave(1)*a)/)
             Cainv(3,:) = (/0.0d0,0.0d0,1.0d0,0.0d0/)
             Cainv(4,:) = (/0.0d0,0.5d0,0.0d0,0.5d0/)
             !Defining Ca
             Ca(1,:) = (/1.0d0,0.0d0,0.0d0,-1.0d0/a2/)
             Ca(2,:) = (/0.0d0,Vave(1)*a,0.0d0,1.0d0/)
             Ca(3,:) = (/0.0d0,0.0d0,1.0d0,0.0d0/)
             Ca(4,:) = (/0.0d0,-Vave(1)*a,0.0d0,1.0d0/)


             !Defining Rotation Matrix
             cosa = ndf(dir,1)
             sina = ndf(dir,2)
             up = cosa*Vave(2) + sina*Vave(3)
             da = ndf(dir,3)
             !Defining Rainv
             Rainv(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
             Rainv(2,:) = (/0.0d0,cosa,-sina,0.0d0/)
             Rainv(3,:) = (/0.0d0,sina,cosa,0.0d0/)
             Rainv(4,:) = (/0.0d0,0.0d0,0.0d0,1.0d0/)
             !Defining Ra
             Ra(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
             Ra(2,:) = (/0.0d0,cosa,sina,0.0d0/)
             Ra(3,:) = (/0.0d0,-sina,cosa,0.0d0/)
             Ra(4,:) = (/0.0d0,0.0d0,0.0d0,1.0d0/)


             !Defining Lambda
             Lambda(1,:,1) = (/da*up,0.0d0,0.0d0,0.0d0/)
             Lambda(2,:,1) = (/0.0d0,da*(up+a),0.0d0,0.0d0/)
             Lambda(3,:,1) = (/0.0d0,0.0d0,da*up,0.0d0/)
             Lambda(4,:,1) = (/0.0d0,0.0d0,0.0d0,da*(up-a)/)

             !Finding forward and backward characteristics
             do i = 1,4
               Lambda(i,i,k) = (Lambda(i,i,1) + sgn*sqrt(Lambda(i,i,1)**2 + eps**2))/2.0d0
             enddo
             !Multiplying through the matrices to get A+ and A-
             dims = (/4,4,4,4/)
             call mat_mult(Sinv,Rainv,m1,dims)
             call mat_mult(m1,Cainv,m2,dims)
             call mat_mult(m2,Lambda(:,:,k),m3,dims)
             call mat_mult(m3,Ca,m4,dims)
             call mat_mult(m4,Ra,m5,dims)
             call mat_mult(m5,S,Af(:,:,k),dims)

          enddo

          if(ipass == 1) then
           !Finding Flux at the face (always from the MSW split)
           dims = (/4,4,4,1/)
           call mat_mult(Af(:,:,1),U(1,:),F1,dims)
           call mat_mult(Af(:,:,2),U(2,:),F2,dims)
           F = F1+F2
          endif
          enddo
          end subroutine
         
          subroutine calc_vf(F,U,dir,Af,df)
          implicit none
          integer*8, intent(in) :: dir
          real*8, intent(inout) :: U(:,:),F(:),Af(:,:),df(:)
          integer*8 :: i,j,k,dims(4) 
          real*8 :: Mxx(4,4),Mxy(4,4),Myx(4,4),Myy(4,4),S(4,4)
          real*8 :: Mee(4,4), Men(4,4), Mne(4,4), Mnn(4,4),d(4)
          real*8 :: VL(4),VR(4),Vd(4),gm1,vsq
            call dfromp(U(1,:),VL,-2)
            call dfromp(U(2,:),VR,-2)
            

            Mxx(1,:) = (/0.0d0,0.0d0,0.0d0,0.0d0/)
            Mxx(2,:) = (/0.0d0,lamu+2*mu,0.0d0,0.0d0/)
            Mxx(3,:) = (/0.0d0,0.0d0,mu,0.0d0/)
            Mxx(4,:) = (/0.0d0,VL(2)*(lamu+2*mu),VL(3)*mu,Kq/)

            Mxy(1,:) = (/0.0d0,0.0d0,0.0d0,0.0d0/)
            Mxy(2,:) = (/0.0d0,0.0d0,lamu,0.0d0/)
            Mxy(3,:) = (/0.0d0,mu,0.0d0,0.0d0/)
            Mxy(4,:) = (/0.0d0,VL(3)*mu,VL(2)*lamu,0.0d0/)

            Myx(1,:) = (/0.0d0,0.0d0,0.0d0,0.0d0/)
            Myx(2,:) = (/0.0d0,0.0d0,mu,00d0/)
            Myx(3,:) = (/0.0d0,lamu,0.0d0,0.0d0/)
            Myx(4,:) = (/0.0d0,VL(3)*lamu,VL(2)*mu,0.0d0/)
            
            Myy(1,:) = (/0.0d0,0.0d0,0.0d0,0.0d0/)
            Myy(2,:) = (/0.0d0,mu,0.0d0,0.0d0/)
            Myy(3,:) = (/0.0d0,0.0d0,lamu+2*mu,0.0d0/)
            Myy(4,:) = (/0.0d0,VL(2)*mu,VL(3)*(lamu+2*mu),Kq/)
            
            d = (/df(4)*df(4),df(2)*df(4),df(4)*df(2),df(2)*df(2)/)
            Mee =  (1.0d0/df(5))*(d(1)*Mxx - d(2)*Mxy - d(3)*Myx + d(4)*Myy)

            d = (/df(3)*df(4),df(1)*df(4),df(3)*df(2),df(1)*df(2)/)
            Men = -(1.0d0/df(5))*(d(1)*Mxx - d(2)*Mxy - d(3)*Myx + d(4)*Myy)

            d = (/df(4)*df(3),df(2)*df(3),df(4)*df(1),df(2)*df(1)/)
            Mne = -(1.0d0/df(5))*(d(1)*Mxx - d(2)*Mxy - d(3)*Myx + d(4)*Myy)

            d = (/df(3)*df(3),df(1)*df(3),df(3)*df(1),df(1)*df(1)/)
            Mnn =  (1.0d0/df(5))*(d(1)*Mxx - d(2)*Mxy - d(3)*Myx + d(4)*Myy)
            
            !Differences across the face; heat flux needs dT, not dp
            Vd = VR-VL
            Vd(4) = VR(4)/(Rair*VR(1)) - VL(4)/(Rair*VL(1))
            dims=(/4,4,4,1/)
            !Calculating TLNS viscous flux
            call mat_mult(Mnn,Vd,F,dims)

            !Calculating TLNS viscous jacobian
            gm1 = gam-1.0d0
            vsq = 0.5d0*(VL(2)*VL(2) + VL(3)*VL(3))
            s(1,:) = (/1.0d0,0.0d0,0.0d0,0.0d0/)
            s(2,:) = (/-VL(2)/VL(1),1.0d0/VL(1),0.0d0,0.0d0/)
            s(3,:) = (/-VL(3)/VL(1),0.0d0,1.0d0/VL(1),0.0d0/)
            s(4,:) = (/vsq*gm1,-VL(2)*gm1,-VL(3)*gm1,gm1/)
            !4th variable is T, not p: dT = (dp - (p/rho)*drho)/(rho*R)
            s(4,:) = (s(4,:) - (VL(4)/VL(1))*s(1,:))/(VL(1)*Rair)
 
            dims=(/4,4,4,4/)
            call mat_mult(Mnn,S,Af,dims)
          end subroutine
        end module
