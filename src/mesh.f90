module mesh
!-----------------------------------------------------------------------
! Unstructured 2D mesh: SU2 reader, face connectivity, geometry, checks
! and VTK output.
!
! Cells are triangles or quads (stored counter-clockwise). Faces are
! numbered interior first (1..nintface), then boundary (nintface+1..nface).
! Every face has a left and a right cell: for a boundary face the right
! "cell" is a ghost cell numbered ncell+k (k = 1..nbface), so cell arrays
! sized (ncell+nbface) hold one ghost state per boundary face.
! Face normals are unit vectors pointing from the left cell to the right.
!-----------------------------------------------------------------------
implicit none
private

!Sizes
integer, public :: nnode = 0     !nodes
integer, public :: ncell = 0     !real cells
integer, public :: nface = 0     !all faces
integer, public :: nintface = 0  !interior faces (1..nintface)
integer, public :: nbface = 0    !boundary faces (nintface+1..nface) = ghost cells
integer, public :: nmark = 0     !boundary markers

!Nodes
real*8, allocatable, public :: xy(:,:)          !(2,nnode)

!Cells
integer, allocatable, public :: nvert(:)        !(ncell) 3 or 4
integer, allocatable, public :: cnode(:,:)      !(4,ncell) node ids, CCW
integer, allocatable, public :: ncface(:)       !(ncell) number of faces
integer, allocatable, public :: cface(:,:)      !(4,ncell) face ids
real*8, allocatable, public :: cvol(:)          !(ncell) area
real*8, allocatable, public :: cctr(:,:)        !(2,ncell+nbface) centres, incl. ghosts

!Faces
integer, allocatable, public :: fnode(:,:)      !(2,nface) end nodes
integer, allocatable, public :: fcell(:,:)      !(2,nface) left, right (right > ncell: ghost)
integer, allocatable, public :: fmark(:)        !(nface) marker id, 0 for interior
real*8, allocatable, public :: fnrm(:,:)        !(2,nface) unit normal, left -> right
real*8, allocatable, public :: flen(:)          !(nface) length
real*8, allocatable, public :: fmid(:,:)        !(2,nface) midpoint

!Markers
character(len=64), allocatable, public :: markname(:)   !(nmark)

!Cell centre definition: the area centroid (needed for gradients and
!second-order reconstruction). Setting this to .true. uses the vertex
!average instead, as the structured v1 code did.
logical, public :: use_vertex_average = .false.

!Marker edges read from the file (used while building faces)
integer :: nmedge = 0
integer, allocatable :: medge(:,:)       !(2,nmedge) nodes
integer, allocatable :: medge_mark(:)    !(nmedge) marker id

public :: mesh_setup, read_su2, build_faces, compute_geometry, check_mesh
public :: write_mesh_vtk, write_boundary_vtk, marker_id

contains

   subroutine mesh_setup(fname)
   !Read an SU2 mesh and build everything the solver needs
   implicit none
   character(len=*), intent(in) :: fname
    call read_su2(fname)
    call build_faces()
    call compute_geometry()
    call check_mesh()
   end subroutine

   !--------------------------------------------------------------------
   ! SU2 reader
   !--------------------------------------------------------------------
   subroutine read_su2(fname)
   implicit none
   character(len=*), intent(in) :: fname
   integer, parameter :: un = 31
   character(len=512) :: line,key,val
   integer :: ios,ndime,n,k,itype,eq,m,nedge_m,iv(4),ie
   logical :: eof,got_elem,got_poin
   real*8 :: x,y
   integer, allocatable :: itmp(:,:),jtmp(:)

    open(unit=un,file=trim(fname),status='old',action='read',iostat=ios)
    if(ios /= 0) call mesh_error('cannot open mesh file '//trim(fname))
    ndime = 0
    got_elem = .false.
    got_poin = .false.
    nmark = 0
    nmedge = 0
    m = 0

    do
      call next_line(un,line,eof)
      if(eof) exit
      eq = index(line,'=')
      if(eq == 0) call mesh_error('unexpected line in mesh file: '//trim(line))
      key = adjustl(line(1:eq-1))
      val = adjustl(line(eq+1:))

      select case(trim(key))
      case('NDIME')
        read(val,*) ndime
        if(ndime /= 2) call mesh_error('only 2D meshes (NDIME= 2) are supported')

      case('NELEM')
        read(val,*) ncell
        allocate(nvert(ncell),cnode(4,ncell))
        cnode = 0
        do k = 1,ncell
          call next_line(un,line,eof)
          if(eof) call mesh_error('end of file inside NELEM block')
          read(line,*) itype
          select case(itype)
          case(5)       !triangle
            read(line,*) itype,iv(1:3)
            nvert(k) = 3
          case(9)       !quadrilateral
            read(line,*) itype,iv(1:4)
            nvert(k) = 4
          case default
            call mesh_error('unsupported element type in NELEM (only 5 = triangle, 9 = quad)')
          end select
          cnode(1:nvert(k),k) = iv(1:nvert(k)) + 1     !SU2 is 0-based
        enddo
        got_elem = .true.

      case('NPOIN')
        read(val,*) nnode
        allocate(xy(2,nnode))
        do k = 1,nnode
          call next_line(un,line,eof)
          if(eof) call mesh_error('end of file inside NPOIN block')
          read(line,*) x,y
          xy(1,k) = x
          xy(2,k) = y
        enddo
        got_poin = .true.

      case('NMARK')
        read(val,*) nmark
        allocate(markname(nmark))
        markname = ''
        m = 0

      case('MARKER_TAG')
        m = m+1
        if(m > nmark) call mesh_error('more MARKER_TAG entries than NMARK')
        markname(m) = trim(val)

      case('MARKER_ELEMS')
        if(m == 0) call mesh_error('MARKER_ELEMS before MARKER_TAG')
        read(val,*) nedge_m
        !grow the marker-edge list
        allocate(itmp(2,nmedge+nedge_m),jtmp(nmedge+nedge_m))
        if(nmedge > 0) then
          itmp(:,1:nmedge) = medge(:,1:nmedge)
          jtmp(1:nmedge) = medge_mark(1:nmedge)
        endif
        call move_alloc(itmp,medge)
        call move_alloc(jtmp,medge_mark)
        do k = 1,nedge_m
          call next_line(un,line,eof)
          if(eof) call mesh_error('end of file inside MARKER_ELEMS block')
          read(line,*) itype
          if(itype /= 3) call mesh_error('boundary elements must be lines (type 3)')
          read(line,*) itype,ie,n
          medge(1,nmedge+k) = ie+1
          medge(2,nmedge+k) = n+1
          medge_mark(nmedge+k) = m
        enddo
        nmedge = nmedge + nedge_m

      case default
        !ignore keywords not needed in 2D (e.g. NZONE, IZONE)
      end select
    enddo
    close(un)

    if(ndime /= 2) call mesh_error('NDIME= 2 not found')
    if(.NOT. got_elem) call mesh_error('NELEM block not found')
    if(.NOT. got_poin) call mesh_error('NPOIN block not found')
    if(m /= nmark) call mesh_error('number of MARKER_TAG entries does not match NMARK')
    if(any(cnode(1:3,:) < 1) .OR. any(cnode > nnode)) &
       call mesh_error('element refers to a node outside 0..NPOIN-1')
    if(nmedge > 0) then
      if(any(medge < 1) .OR. any(medge > nnode)) &
         call mesh_error('boundary edge refers to a node outside 0..NPOIN-1')
    endif

    call orient_cells()

    print*,'Mesh: ',trim(fname)
    print*,'   nodes =',nnode,'  cells =',ncell,'  (tri',count(nvert==3),', quad',count(nvert==4),')'
    print*,'   markers =',nmark,'  boundary edges =',nmedge
   end subroutine

   subroutine next_line(un,line,eof)
   !Next non-blank line with '%' comments removed
   implicit none
   integer, intent(in) :: un
   character(len=*), intent(out) :: line
   logical, intent(out) :: eof
   integer :: ios,ic
    eof = .false.
    do
      read(un,'(A)',iostat=ios) line
      if(ios /= 0) then
        eof = .true.
        line = ''
        return
      endif
      ic = index(line,'%')
      if(ic > 0) line(ic:) = ''
      line = adjustl(line)
      if(len_trim(line) > 0) return
    enddo
   end subroutine

   subroutine orient_cells()
   !Make every cell counter-clockwise (positive signed area)
   implicit none
   integer :: c,nv,nflip
    nflip = 0
    do c = 1,ncell
      nv = nvert(c)
      if(signed_area(c) < 0.0d0) then
        cnode(1:nv,c) = cnode(nv:1:-1,c)
        nflip = nflip+1
      endif
    enddo
    if(nflip > 0) print*,'   reoriented',nflip,' clockwise cells'
   end subroutine

   real*8 function signed_area(c)
   !Shoelace formula
   implicit none
   integer, intent(in) :: c
   integer :: k,a,b,nv
    nv = nvert(c)
    signed_area = 0.0d0
    do k = 1,nv
      a = cnode(k,c)
      b = cnode(mod(k,nv)+1,c)
      signed_area = signed_area + xy(1,a)*xy(2,b) - xy(1,b)*xy(2,a)
    enddo
    signed_area = 0.5d0*signed_area
   end function

   !--------------------------------------------------------------------
   ! Face connectivity
   !--------------------------------------------------------------------
   subroutine build_faces()
   !Faces are found by matching cell edges on the key (min node, max node).
   !Edges are bucketed by their min node (CSR), so no sorting is needed.
   implicit none
   integer :: nhe,c,k,nv,a,b,lo,hi,s,t,ib,ii,f,nmatch
   integer, allocatable :: hstart(:),hfill(:),hhi(:),hcell(:),hk(:),hpair(:)
   integer, allocatable :: mstart(:),mfill(:),mhi(:),mid(:)
   logical, allocatable :: mused(:)

    !--- half-edges (one per cell edge), bucketed by min node
    nhe = sum(nvert)
    allocate(hstart(nnode+1),hfill(nnode),hhi(nhe),hcell(nhe),hk(nhe),hpair(nhe))
    hstart = 0
    do c = 1,ncell
      nv = nvert(c)
      do k = 1,nv
        a = cnode(k,c)
        b = cnode(mod(k,nv)+1,c)
        if(a == b) call mesh_error('degenerate cell edge (repeated node)')
        lo = min(a,b)
        hstart(lo+1) = hstart(lo+1) + 1
      enddo
    enddo
    hstart(1) = 1
    do a = 1,nnode
      hstart(a+1) = hstart(a+1) + hstart(a)
    enddo
    hfill = 0
    do c = 1,ncell
      nv = nvert(c)
      do k = 1,nv
        a = cnode(k,c)
        b = cnode(mod(k,nv)+1,c)
        lo = min(a,b)
        s = hstart(lo) + hfill(lo)
        hfill(lo) = hfill(lo) + 1
        hhi(s) = max(a,b)
        hcell(s) = c
        hk(s) = k
      enddo
    enddo

    !--- pair up matching half-edges
    hpair = 0
    nintface = 0
    nbface = 0
    do lo = 1,nnode
      do s = hstart(lo),hstart(lo+1)-1
        if(hpair(s) /= 0) cycle
        nmatch = 0
        do t = s+1,hstart(lo+1)-1
          if(hhi(t) == hhi(s)) then
            nmatch = nmatch+1
            if(nmatch > 1) call mesh_error('edge shared by more than two cells')
            if(hpair(t) /= 0) call mesh_error('edge shared by more than two cells')
            hpair(s) = t
            hpair(t) = s
          endif
        enddo
        if(nmatch == 1) then
          nintface = nintface+1
        else
          hpair(s) = -1          !boundary half-edge
          nbface = nbface+1
        endif
      enddo
    enddo
    nface = nintface + nbface

    !--- marker edges, bucketed the same way
    allocate(mstart(nnode+1),mfill(nnode),mhi(max(nmedge,1)),mid(max(nmedge,1)),mused(max(nmedge,1)))
    mstart = 0
    do ii = 1,nmedge
      lo = min(medge(1,ii),medge(2,ii))
      mstart(lo+1) = mstart(lo+1) + 1
    enddo
    mstart(1) = 1
    do a = 1,nnode
      mstart(a+1) = mstart(a+1) + mstart(a)
    enddo
    mfill = 0
    do ii = 1,nmedge
      lo = min(medge(1,ii),medge(2,ii))
      s = mstart(lo) + mfill(lo)
      mfill(lo) = mfill(lo) + 1
      mhi(s) = max(medge(1,ii),medge(2,ii))
      mid(s) = medge_mark(ii)
    enddo
    mused = .false.

    !--- create the faces: interior first, then boundary
    allocate(fnode(2,nface),fcell(2,nface),fmark(nface))
    fmark = 0
    f = 0
    ib = nintface
    do lo = 1,nnode
      do s = hstart(lo),hstart(lo+1)-1
        c = hcell(s)
        k = hk(s)
        nv = nvert(c)
        if(hpair(s) > 0) then
          if(hpair(s) < s) cycle            !pair already created
          !Left cell = this half-edge's cell; nodes in its CCW order,
          !so the normal (dy,-dx) points out of the left cell
          f = f+1
          fnode(1,f) = cnode(k,c)
          fnode(2,f) = cnode(mod(k,nv)+1,c)
          fcell(1,f) = c
          fcell(2,f) = hcell(hpair(s))
        else
          ib = ib+1
          fnode(1,ib) = cnode(k,c)
          fnode(2,ib) = cnode(mod(k,nv)+1,c)
          fcell(1,ib) = c
          fcell(2,ib) = ncell + (ib-nintface)     !ghost cell
          !find this edge's marker
          hi = hhi(s)
          do t = mstart(lo),mstart(lo+1)-1
            if(mhi(t) == hi) then
              fmark(ib) = mid(t)
              mused(t) = .true.
              exit
            endif
          enddo
          if(fmark(ib) == 0) then
            print*,'Boundary edge with no marker: nodes',fnode(1,ib)-1,fnode(2,ib)-1,' (0-based)'
            call mesh_error('every boundary edge must belong to a marker')
          endif
        endif
      enddo
    enddo
    if(nmedge > 0) then
      if(.NOT. all(mused(1:nmedge))) then
        do t = 1,nmedge
          if(.NOT. mused(t)) then
            print*,'Marker edge that is not on the boundary: marker ',trim(markname(mid(t)))
            exit
          endif
        enddo
        call mesh_error('marker edges must lie on the domain boundary')
      endif
    endif

    !--- faces of each cell
    allocate(ncface(ncell),cface(4,ncell))
    ncface = 0
    cface = 0
    do f = 1,nface
      c = fcell(1,f)
      ncface(c) = ncface(c)+1
      cface(ncface(c),c) = f
      c = fcell(2,f)
      if(c <= ncell) then
        ncface(c) = ncface(c)+1
        cface(ncface(c),c) = f
      endif
    enddo

    deallocate(hstart,hfill,hhi,hcell,hk,hpair,mstart,mfill,mhi,mid,mused)
    print*,'   faces =',nface,'  (interior',nintface,', boundary',nbface,')'
   end subroutine

   !--------------------------------------------------------------------
   ! Geometry
   !--------------------------------------------------------------------
   subroutine compute_geometry()
   implicit none
   integer :: c,k,nv,a,b,f,g
   real*8 :: cr,ax,ay,dx,dy,dn

    if(allocated(cvol)) deallocate(cvol,cctr,fnrm,flen,fmid)
    allocate(cvol(ncell),cctr(2,ncell+nbface))
    allocate(fnrm(2,nface),flen(nface),fmid(2,nface))

    !Cells: area (shoelace) and centre
    do c = 1,ncell
      nv = nvert(c)
      cvol(c) = signed_area(c)
      if(use_vertex_average) then
        cctr(1,c) = sum(xy(1,cnode(1:nv,c)))/dble(nv)
        cctr(2,c) = sum(xy(2,cnode(1:nv,c)))/dble(nv)
      else
        !area centroid of a polygon
        ax = 0.0d0
        ay = 0.0d0
        do k = 1,nv
          a = cnode(k,c)
          b = cnode(mod(k,nv)+1,c)
          cr = xy(1,a)*xy(2,b) - xy(1,b)*xy(2,a)
          ax = ax + (xy(1,a)+xy(1,b))*cr
          ay = ay + (xy(2,a)+xy(2,b))*cr
        enddo
        cctr(1,c) = ax/(6.0d0*cvol(c))
        cctr(2,c) = ay/(6.0d0*cvol(c))
      endif
    enddo

    !Faces: midpoint, length and unit normal (left -> right)
    do f = 1,nface
      a = fnode(1,f)
      b = fnode(2,f)
      dx = xy(1,b) - xy(1,a)
      dy = xy(2,b) - xy(2,a)
      flen(f) = sqrt(dx*dx + dy*dy)
      fnrm(1,f) =  dy/flen(f)
      fnrm(2,f) = -dx/flen(f)
      fmid(1,f) = 0.5d0*(xy(1,a)+xy(1,b))
      fmid(2,f) = 0.5d0*(xy(2,a)+xy(2,b))
    enddo

    !Ghost-cell centres: left centre mirrored across the face line, so the
    !vector from a cell to its ghost is normal to the boundary
    do f = nintface+1,nface
      g = fcell(2,f)
      dn = dot_product(cctr(:,fcell(1,f)) - fmid(:,f), fnrm(:,f))
      cctr(:,g) = cctr(:,fcell(1,f)) - 2.0d0*dn*fnrm(:,f)
    enddo
   end subroutine

   !--------------------------------------------------------------------
   ! Checks
   !--------------------------------------------------------------------
   subroutine check_mesh()
   implicit none
   integer :: c,f,k,m,euler,worst
   real*8 :: s(2),err,errmax,perim,total
   integer, allocatable :: nfm(:)

    !1) positive areas
    if(any(cvol <= 0.0d0)) then
      print*,'Cell with non-positive area:',minloc(cvol)
      call mesh_error('mesh check failed: non-positive cell area')
    endif

    !2) closure: outward normals * lengths sum to zero for every cell
    errmax = 0.0d0
    worst = 0
    do c = 1,ncell
      s = 0.0d0
      perim = 0.0d0
      do k = 1,ncface(c)
        f = cface(k,c)
        if(fcell(1,f) == c) then
          s = s + fnrm(:,f)*flen(f)
        else
          s = s - fnrm(:,f)*flen(f)
        endif
        perim = perim + flen(f)
      enddo
      err = sqrt(s(1)**2 + s(2)**2)/perim
      if(err > errmax) then
        errmax = err
        worst = c
      endif
    enddo
    if(errmax > 1.0d-12) then
      print*,'Worst closure error',errmax,' in cell',worst
      call mesh_error('mesh check failed: cell faces do not close')
    endif

    !3) each cell has as many faces as vertices
    if(any(ncface /= nvert)) call mesh_error('mesh check failed: cell face count /= vertex count')

    !4) every boundary face has a marker, every marker has faces
    allocate(nfm(nmark))
    nfm = 0
    do f = nintface+1,nface
      if(fmark(f) < 1) call mesh_error('mesh check failed: boundary face without marker')
      nfm(fmark(f)) = nfm(fmark(f)) + 1
    enddo

    total = sum(cvol)
    euler = nnode - nface + ncell
    print*,'Mesh check passed'
    print*,'   total area       =',total
    print*,'   min / max area   =',minval(cvol),maxval(cvol)
    print*,'   max closure error=',errmax
    print*,'   nodes - faces + cells =',euler,'  (1 = no holes; holes =',1-euler,')'
    print*,'   boundary markers:'
    do m = 1,nmark
      print '(6X,I3,2X,A20,I8,A)',m,markname(m),nfm(m),' faces'
      if(nfm(m) == 0) print*,'      warning: marker ',trim(markname(m)),' has no faces'
    enddo
    deallocate(nfm)
   end subroutine

   integer function marker_id(name)
   !Marker number for a marker name (0 if not found)
   implicit none
   character(len=*), intent(in) :: name
   integer :: m
    marker_id = 0
    do m = 1,nmark
      if(trim(markname(m)) == trim(name)) then
        marker_id = m
        return
      endif
    enddo
   end function

   !--------------------------------------------------------------------
   ! VTK output (legacy ASCII, opens in ParaView)
   !--------------------------------------------------------------------
   subroutine write_mesh_vtk(fname)
   !Cells with their id and area
   implicit none
   character(len=*), intent(in) :: fname
   integer, parameter :: un = 32
   integer :: c,k
    open(unit=un,file=trim(fname),status='replace',action='write')
    write(un,'(A)') '# vtk DataFile Version 3.0'
    write(un,'(A)') 'SCUFFS mesh'
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
    write(un,'(A)') 'SCALARS cell_id int 1'
    write(un,'(A)') 'LOOKUP_TABLE default'
    do c = 1,ncell
      write(un,'(I0)') c
    enddo
    write(un,'(A)') 'SCALARS area double 1'
    write(un,'(A)') 'LOOKUP_TABLE default'
    do c = 1,ncell
      write(un,'(ES24.15)') cvol(c)
    enddo
    close(un)
    print*,'Wrote ',trim(fname)
   end subroutine

   subroutine write_boundary_vtk(fname)
   !Boundary faces as line cells, coloured by marker id
   implicit none
   character(len=*), intent(in) :: fname
   integer, parameter :: un = 33
   integer :: f,k
    open(unit=un,file=trim(fname),status='replace',action='write')
    write(un,'(A)') '# vtk DataFile Version 3.0'
    write(un,'(A)') 'SCUFFS boundary faces (marker ids are listed in the mesh check output)'
    write(un,'(A)') 'ASCII'
    write(un,'(A)') 'DATASET UNSTRUCTURED_GRID'
    write(un,'(A,I0,A)') 'POINTS ',nnode,' double'
    do k = 1,nnode
      write(un,'(3ES24.15)') xy(1,k),xy(2,k),0.0d0
    enddo
    write(un,'(A,I0,1X,I0)') 'CELLS ',nbface,3*nbface
    do f = nintface+1,nface
      write(un,'(A,I0,1X,I0)') '2 ',fnode(1,f)-1,fnode(2,f)-1
    enddo
    write(un,'(A,I0)') 'CELL_TYPES ',nbface
    do f = nintface+1,nface
      write(un,'(I0)') 3
    enddo
    write(un,'(A,I0)') 'CELL_DATA ',nbface
    write(un,'(A)') 'SCALARS marker int 1'
    write(un,'(A)') 'LOOKUP_TABLE default'
    do f = nintface+1,nface
      write(un,'(I0)') fmark(f)
    enddo
    write(un,'(A)') 'VECTORS normal double'
    do f = nintface+1,nface
      write(un,'(3ES24.15)') fnrm(1,f),fnrm(2,f),0.0d0
    enddo
    close(un)
    print*,'Wrote ',trim(fname)
   end subroutine

   subroutine mesh_error(msg)
   implicit none
   character(len=*), intent(in) :: msg
    print*,'Mesh error: ',trim(msg)
    stop 1
   end subroutine

end module mesh
