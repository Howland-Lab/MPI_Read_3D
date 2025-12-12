!===============================================================
! MPI_diag.F90
!===============================================================

module MPIR3D
  use mpi
  use decomp_2d_io
  use decomp_2d, only: decomp_info, decomp_2d_init, xstart, &
                       xend, ystart, yend, zstart, zend, &
                       decomp_2d_finalize, decomp_info_init, &
                       mytype
  implicit none

  integer, parameter :: rk = mytype
  integer :: myrank = -1
  integer :: nprocs = -1

  type :: FieldReader2Decomp
     private
     ! Global sizes and process grid
     integer :: nx = 0, ny = 0, nz = 0
     integer :: p_row = 0, p_col = 0
     logical :: is_init = .false.

     ! 2DECOMP&FFT descriptor (x-pencil)
     type(DECOMP_INFO) :: d

     ! Local extents for x-pencil
     integer :: xs=0, xe=0, ys=0, ye=0, zs=0, ze=0
     integer :: nxloc=0, nyloc=0, nzloc=0

   contains
     procedure :: init          => frd_init
     procedure :: read_field    => frd_read_field
     procedure :: global_shape  => frd_global_shape
     procedure :: local_shape   => frd_local_shape
     procedure :: indices       => frd_index
     procedure, private :: choose_proc_grid_ => frd_choose_proc_grid
     final       :: frd_finalize
  end type FieldReader2Decomp

contains

  ! Initialize the FieldReader2Decomp
  subroutine frd_init(this, nx, ny, nz)
    class(FieldReader2Decomp), intent(inout) :: this
    integer,               intent(in)    :: nx, ny, nz
    integer :: ierr
    character(256) :: buf
    integer :: loc(6)
    integer, allocatable :: recvbuf(:)
    integer :: r, off

    this%nx = nx; this%ny = ny; this%nz = nz
    if (myrank == 0) allocate(recvbuf(6*nprocs))
    
    ! Choose process grid and init 2DECOMP
    call this%choose_proc_grid_(nprocs, this%nx, this%ny, this%nz, this%p_row, this%p_col)
    
    if (myrank == 0) then
      write(buf,'("FRD_INIT: global size = ",I0,", ",I0,", ",I0)') this%nx, this%ny, this%nz
      call message(trim(buf))

      write(buf,'("FRD_INIT: proc grid = ",I0," x ",I0)') this%p_row, this%p_col
      call message(trim(buf))
    end if

    call decomp_2d_init(this%nx, this%ny, this%nz, this%p_row, this%p_col)
    call decomp_info_init(this%nx, this%ny, this%nz, this%d)
    
    ! Cache local extents
    this%xs = xstart(1); this%xe = xend(1)
    this%ys = xstart(2); this%ye = xend(2)
    this%zs = xstart(3); this%ze = xend(3)
    this%nxloc = this%xe - this%xs + 1
    this%nyloc = this%ye - this%ys + 1
    this%nzloc = this%ze - this%zs + 1    

    ! Store local extents in temp
    loc = [this%xs, this%xe, this%ys, this%ye, this%zs, this%ze]
    call MPI_Gather(loc, 6, MPI_INTEGER, recvbuf, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Print local extents
    if (myrank == 0) then      
      do r = 0, nprocs-1
        off = 6*r
        write(buf,'("FRD_INIT: Rank=",I0," X: ",I0,", ",I0," // Y: ",I0,", ",I0," // Z: ",I0,", ",I0)') &
            r, recvbuf(off+1), recvbuf(off+2), recvbuf(off+3), recvbuf(off+4), recvbuf(off+5), recvbuf(off+6)
        call message(trim(buf))
      end do
      deallocate(recvbuf)
    end if
    
    this%is_init = .true.

    if (myrank == 0) then
      call message('FRD_INIT: 2DECOMP reader initialized.')
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end subroutine frd_init

  ! Read a 3D field from file
  function frd_read_field(this, path) result(field)
    class(FieldReader2Decomp), intent(in) :: this
    character(*),               intent(in) :: path
    real(rk) :: field(this%nxloc, this%nyloc, this%nzloc)
    integer :: ierr

    if (.not. this%is_init) then
      if (myrank == 0) call message('ERROR(FRD): call init() before read_field().')
      call MPI_Abort(MPI_COMM_WORLD, 101, ierr)
    end if

    if (myrank == 0) call message('Reading from: '//path(:len_trim(path)))
    call decomp_2d_read_one(1, field, path, this%d)
  end function frd_read_field

  ! Get global shape
  subroutine frd_global_shape(this, nx, ny, nz)
    class(FieldReader2Decomp), intent(in)  :: this
    integer,                    intent(out) :: nx, ny, nz
    nx = this%nx; ny = this%ny; nz = this%nz
  end subroutine frd_global_shape

  ! Get local shape
  subroutine frd_local_shape(this, nxloc, nyloc, nzloc)
    class(FieldReader2Decomp), intent(in)  :: this
    integer,                    intent(out) :: nxloc, nyloc, nzloc
    nxloc = this%nxloc; nyloc = this%nyloc; nzloc = this%nzloc
  end subroutine frd_local_shape

  ! Get local indices
  subroutine frd_index(this, xs, xe, ys, ye, zs, ze)
    class(FieldReader2Decomp), intent(in)  :: this
    integer,                    intent(out) :: xs, xe, ys, ye, zs, ze
    xs=this%xs; xe=this%xe
    ys=this%ys; ye=this%ye
    zs=this%zs; ze=this%ze
  end subroutine frd_index

  ! Finalizer
  subroutine frd_finalize(this)
    type(FieldReader2Decomp), intent(inout) :: this
    if (this%is_init) then
      call decomp_2d_finalize()
      this%is_init = .false.
    end if
  end subroutine frd_finalize

  ! Choose process grid based on global sizes
  subroutine frd_choose_proc_grid(this, nprocs, nx, ny, nz, p_row, p_col)
    class(FieldReader2Decomp), intent(in)  :: this
    integer,                    intent(in)  :: nprocs, nx, ny, nz
    integer,                    intent(out) :: p_row, p_col
    integer :: i, br, bc
    real(rk) :: best, sc
    best = huge(1.0_rk); br = 1; bc = nprocs
    do i = 1, nprocs
      if (mod(nprocs, i) == 0) then
        sc = abs(real(nx, rk)/real(ny, rk) - real(nprocs/i, rk)/real(i, rk))
        if (sc < best) then
          best = sc; br = i; bc = nprocs/i
        end if
      end if
    end do
    p_row = br; p_col = bc
  end subroutine frd_choose_proc_grid

  ! Compute horizontal average speed and direction
  subroutine do_horizontal_average_ws_wd(nz, nxloc, nyloc, nzloc, zs, u, v, ws, wd)
    integer, intent(in)  :: nz, nxloc, nyloc, nzloc, zs
    real(rk), intent(in) :: u(nxloc, nyloc, nzloc), v(nxloc, nyloc, nzloc)
    real(rk), intent(out) :: ws(nz), wd(nz)
    integer :: k, kg, i, j
    integer :: ierr
    real(rk) :: sum_spd(nz), sum_c(nz), sum_s(nz), sum_cnt(nz), buff(nz)
    real(rk) :: uu, vv, mag, cbar, sbar, ang

    sum_spd= 0.0_rk; sum_c= 0.0_rk; sum_s= 0.0_rk; sum_cnt= 0.0_rk
    buff = 0.0_rk

    do k = 1, nzloc
      kg = zs + k - 1
      do j = 1, nyloc
        do i = 1, nxloc
          uu = u(i,j,k); vv = v(i,j,k)
          mag = sqrt(uu*uu + vv*vv)                  
          if (mag > 0.0_rk) then    
            sum_spd(kg) = sum_spd(kg) + mag        
            sum_c(kg) = sum_c(kg) + uu/mag
            sum_s(kg) = sum_s(kg) + vv/mag
            sum_cnt(kg) = sum_cnt(kg) + 1.0_rk
          end if
        end do
      end do
    end do

    call MPI_Allreduce(sum_spd, buff, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    sum_spd = buff
    call MPI_Allreduce(sum_c, buff, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    sum_c = buff
    call MPI_Allreduce(sum_s, buff, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    sum_s = buff
    call MPI_Allreduce(sum_cnt, buff, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    sum_cnt = buff

    do kg = 1, nz      
      if (sum_cnt(kg) > 0.5_rk) then
        ws(kg) = sum_spd(kg) / sum_cnt(kg)
        cbar = sum_c(kg) / sum_cnt(kg)
        sbar = sum_s(kg) / sum_cnt(kg)
        ang  = atan2(sbar, cbar) * (180.0_rk/acos(-1.0_rk))
        wd(kg) = ang
      else
        ws(kg) = 0.0_rk
        wd(kg) = 0.0_rk
      end if
    end do
  end subroutine do_horizontal_average_ws_wd

  ! for a generic field, perform horizontal average
  subroutine do_horizontal_average(nz, nxloc, nyloc, nzloc, zs, f, profile)
    integer, intent(in)  :: nz, nxloc, nyloc, nzloc, zs
    real(rk), intent(in) :: f(nxloc, nyloc, nzloc)
    real(rk), intent(out) :: profile(nz)
    integer :: k, kg, i, j
    integer :: ierr
    real(rk) :: sum_f(nz), sum_cnt(nz), buff(nz)

    sum_f= 0.0_rk; sum_cnt= 0.0_rk; buff=0.0_rk
    do k = 1, nzloc
      kg = zs + k - 1
      do j = 1, nyloc
        do i = 1, nxloc                    
          sum_f(kg) = sum_f(kg) + f(i,j,k)        
          sum_cnt(kg) = sum_cnt(kg) + 1.0_rk
        end do
      end do
    end do

    call MPI_Allreduce(sum_f, buff, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    sum_f = buff
    call MPI_Allreduce(sum_cnt, buff, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    sum_cnt = buff

    do kg = 1, nz      
      if (sum_cnt(kg) > 0.5_rk) then
        profile(kg) = sum_f(kg) / sum_cnt(kg)
      else
        profile(kg) = 0.0_rk
      end if
    end do
  end subroutine do_horizontal_average
  
  subroutine ha_driver(reader, Lz, runid, path, outdir, field)
    class(FieldReader2Decomp), intent(inout) :: reader
    integer, intent(in) :: runid
    character(*), intent(in) :: path, outdir, field
    real(rk), intent(in) :: Lz
    real(rk), allocatable :: f1(:,:,:), z(:)
    character(len=:), allocatable :: keys(:), sorted_keys(:)
    character(len=2) :: rc
    character(len=10):: f_
    integer :: k, nx, ny, nz, xs, xe, ys, ye, zs, ze, nxloc, nyloc, nzloc

    ! Fetch local and global sizes & start indices
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%global_shape(nx, ny, nz)
    call reader%indices(xs, xe, ys, ye, zs, ze)

    ! allocate the data holder    
    allocate(f1(nxloc, nyloc, nzloc))
    
    ! Convert runid to character
    write(rc, '(I2.2)') runid

    ! Convert field name to proper string
    f_ = field_to_name(trim(field))

    ! Vertical coordinate z
    z = linspace(0.0_rk, Lz, nz)   

    ! List files matching pattern
    call list_matching_keys(trim(path), 'Run'//trim(rc)//'_'//trim(f_)//'_t*.out', keys)
    
    ! Sort keys in time order (numeric)
    sorted_keys = sort_keys_numeric(keys)

    do k = 1, size(sorted_keys)

      ! If it is a single field, only f1 is read
      ! for WS and WD profiles, another field f2 (vVel) is also read
      f1 = reader%read_field(trim(path)//'/'//'Run'//trim(rc)//'_'//trim(f_)//'_t'//trim(sorted_keys(k))//'.out')
      if (field == 'S') then
        block
          real(rk) :: f2(nxloc, nyloc, nzloc), ws(nz), wd(nz)

          f2 = reader%read_field(trim(path)//'/'//'Run'//trim(rc)//'_vVel_t'//trim(sorted_keys(k))//'.out')  
          call do_horizontal_average_ws_wd(nz, nxloc, nyloc, nzloc, zs, f1, f2, ws, wd)
          call csvprofile(nz,trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//'_HA_WS.csv',z,ws)
          call csvprofile(nz,trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//'_HA_WD.csv',z,wd)
        end block
      else
        block
          real(rk) :: profile(nz)
          call do_horizontal_average(nz, nxloc, nyloc, nzloc, zs, f1, profile)
          call csvprofile(nz, trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//&
                '_HA_'//trim(f_)//'.csv',z,profile)
        end block
      end if
    end do
  end subroutine ha_driver

  subroutine slice_driver(reader, Lx, Ly, Lz, runid, path, outdir, field, axis, nslice, slice_coord)
    class(FieldReader2Decomp), intent(inout) :: reader
    integer,          intent(in)  :: runid, nslice
    character(*),     intent(in)  :: path, outdir, field
    real(rk),         intent(in)  :: Lx, Ly, Lz, slice_coord(nslice)
    character(*),     intent(in)  :: axis   ! e.g. 'z'

    ! Local 3D field
    real(rk), allocatable, target :: f1(:,:,:)

    ! Grid and indices
    integer :: nx, ny, nz
    integer :: nxloc, nyloc, nzloc
    integer :: xs, xe, ys, ye, zs, ze

    ! Global coordinates
    real(rk), allocatable :: x1(:), x2(:)

    ! Slice indices and interpolation
    integer :: k, k0, k1
    real(rk) :: alpha
    character(len=2) :: rc, term
    character(len=10) :: f_
    character(len=1) :: ax, x1name, x2name, budget, eax
    character(len=256) :: fname, msg
    character(len=2048) :: pattern

    ! Time keys
    character(len=:), allocatable :: keys(:), sorted_keys(:), stamps(:), sorted_stamps(:)

    ! Slice arrays (global shape on every rank)
    real(rk), allocatable :: local_slice0(:,:), local_slice1(:,:)
    real(rk), allocatable :: global_slice0(:,:), global_slice1(:,:), slice_interp(:,:)
    real(rk), pointer :: slice_ptr(:,:)

    integer :: i, j, isl, ierr, nx1, nx2, nax
    integer :: x1s, x1e, x2s, x2e, axs, axe, mode
    real(rk) :: L1, L2, Lax, slice_

    ! Handling reference to wind speed and wind direction
    if(trim(field) == 'S')then
      if(myrank == 0) call message('ERROR(Slice): Export u and v separately and calculate WS & WD offline.')
      call MPI_Abort(MPI_COMM_WORLD, 100, ierr)
    end if

    ! Convert runid to character
    write(rc, '(I2.2)') runid
    ax = to_lower(axis(1:1))

    ! Shapes and local indices
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%global_shape(nx, ny, nz)
    call reader%indices(xs, xe, ys, ye, zs, ze)

    ! Allocate local 3D field
    allocate(f1(nxloc, nyloc, nzloc))

    ! Get file list and sort by time
    mode = field_mode(trim(field))
    if (mode == 0)then
      f_ = field_to_name(trim(field))
      call list_matching_keys(trim(path), 'Run'//trim(rc)//'_'//trim(f_)//'_t*.out', keys)
      sorted_keys = sort_keys_numeric(keys)
    else
      f_ = trim(field)
      call define_budget(trim(field), budget, term) 
      pattern = 'Run'//trim(rc)
      if (mode == 2) pattern = trim(pattern)//'_deficit'
      pattern = trim(pattern)//'_budget'//budget//'_term'//term//'_t*_n~.s3D'
      if (myrank == 0)then
        call message('Pattern is '//trim(pattern))
      end if
      call list_matching_keys_budget(trim(path), trim(pattern), keys, stamps)
      call sort_keys_and_stamps_numeric(keys, stamps, sorted_keys, sorted_stamps)
    end if    

    ! Figure which coordinate names to use
    select case (ax)
    case ('x')
      x1name = 'y'; x2name = 'z'
      nx1 = ny; nx2 = nz; nax = nx
      L1 = Ly; L2 = Lz; Lax = Lx
      x1s = ys; x1e = ye; x2s = zs; x2e = ze
      axs = xs; axe = xe
      eax = 'i'
    case ('y')
      x1name = 'x'; x2name = 'z'
      nx1 = nx; nx2 = nz; nax = ny
      L1 = Lx; L2 = Lz; Lax = Ly
      x1s = xs; x1e = xe; x2s = zs; x2e = ze
      axs = ys; axe = ye
      eax = 'j'
    case ('z')
      x1name = 'x'; x2name = 'y' 
      nx1 = nx; nx2 = ny; nax = nz
      L1 = Lx; L2 = Ly; Lax = Lz
      x1s = xs; x1e = xe; x2s = ys; x2e = ye
      axs = zs; axe = ze
      eax = 'k'
    end select

    ! Global coordinate arrays (same on all ranks)
    x1 = linspace(0.0_rk, L1, nx1)
    x2 = linspace(0.0_rk, L2, nx2)

    ! Global slice arrays (same shape on all ranks)
    allocate(local_slice0(nx1, nx2), local_slice1(nx1, nx2))
    allocate(global_slice0(nx1, nx2), global_slice1(nx1, nx2))
    allocate(slice_interp(nx1, nx2))

    do isl = 1, nslice
      slice_ = slice_coord(isl)
      if (myrank == 0) then
        write(msg,'(A, I0, A, I0, A, ES12.4)') 'Slice ',isl,'/',nslice,', '//ax//' = ',slice_
        call message(trim(msg))
      end if

      ! Wipe clean slice arrays
      local_slice0 = 0.0_rk; local_slice1 = 0.0_rk
      global_slice0 = 0.0_rk; global_slice1 = 0.0_rk
      slice_interp = 0.0_rk

      ! Compute bracketing indices for z
      call find_bracket_uniform(nax, Lax, slice_, k0, k1, alpha)
      if (myrank == 0) then
        write(msg,'(A, I0, A, I0, A, ES12.4)')'k0 = ',k0,', k1 = ',k1,', alpha = ',alpha
        call message(trim(msg))
      end if

      ! Loop over time snapshots
      do k = 1, size(sorted_keys)

        ! 1) Read 3D field for this time step on each rank
        if(mode == 0)then
          fname = 'Run'//trim(rc)//'_'//trim(f_)//'_t'//trim(sorted_keys(k))//'.out'
        else if(mode == 1)then          
          fname = 'Run'//trim(rc)//'_budget'//budget//'_term'//term//'_t'//&
              trim(sorted_keys(k))//'_n'//trim(sorted_stamps(k))//'.s3D'
        else if(mode == 2)then          
          fname = 'Run'//trim(rc)//'_deficit_budget'//budget//'_term'//term//'_t'//&
              trim(sorted_keys(k))//'_n'//trim(sorted_stamps(k))//'.s3D'
        end if
        f1 = reader%read_field(trim(path)//'/'//trim(fname))

        ! 2) Build local contributions to the two bracketing planes
        local_slice0 = 0.0_rk
        local_slice1 = 0.0_rk

        block
          logical :: has0, has1
          integer :: k0_loc, k1_loc

          has0 = (k0 >= axs .and. k0 <= axe)
          has1 = (k1 >= axs .and. k1 <= axe)

          if (has0) then
            k0_loc = k0 - axs + 1
            select case(ax)
            case('x')
              slice_ptr => f1(k0_loc,:,:)
            case('y')
              slice_ptr => f1(:,k0_loc,:)
            case('z')
              slice_ptr => f1(:,:,k0_loc)
            end select
            local_slice0(x1s:x1e, x2s:x2e) = slice_ptr(:,:)
          end if

          if (has1) then
            k1_loc = k1 - axs + 1
            select case(ax)
            case('x')
              slice_ptr => f1(k1_loc,:,:)
            case('y')
              slice_ptr => f1(:,k1_loc,:)
            case('z')
              slice_ptr => f1(:,:,k1_loc)
            end select
            local_slice1(x1s:x1e, x2s:x2e) = slice_ptr(:,:)
          end if
        end block

        ! 3) Sum contributions from all ranks to get full global slices
        call MPI_Allreduce(local_slice0, global_slice0, nx1*nx2, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(local_slice1, global_slice1, nx1*nx2, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr)

        ! 4) On root: do linear interpolation and write CSV
        if (myrank == 0) then
          do j = 1, nx2
            do i = 1, nx1
              slice_interp(i,j) = (1.0_rk - alpha) * global_slice0(i,j) + alpha * global_slice1(i,j)
            end do
          end do

          ! Output file name
          fname = trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//&
                    '_SL_'//trim(f_)//'_'//ax
          if(slice_ <= -1) fname = trim(fname)//'_'//eax ! A direct index is given
          fname = trim(fname)//'='//trim(real2string(slice_))//'.csv'
          
          call writeslice(nx1, nx2, trim(fname), slice_interp)
        end if

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

      end do  ! time loop
    end do

    ! Cleanup
    if (allocated(f1))           deallocate(f1)
    if (allocated(x1))           deallocate(x1)
    if (allocated(x2))           deallocate(x2)
    if (allocated(local_slice0)) deallocate(local_slice0)
    if (allocated(local_slice1)) deallocate(local_slice1)
    if (allocated(global_slice0)) deallocate(global_slice0)
    if (allocated(global_slice1)) deallocate(global_slice1)
    if (allocated(slice_interp))  deallocate(slice_interp)
  end subroutine slice_driver

  pure function real2string(z) result(tag)
    ! Convert slice coordinate to a compact string.
    ! Modes:
    !   z >= 0.0  : treated as a physical location, e.g.
    !               0.250  -> "0p25"
    !               1.000  -> "1"
    !   z  < 0.0  : treated as an index flag, e.g.
    !               -2.0   -> "2"   (level 2)

    real(rk), intent(in) :: z
    character(len=32)    :: tag

    character(len=64) :: tmp
    integer :: i, n, idx

    ! --- Index mode for negative values ---
    if (z < 0.0_rk) then
      idx = nint(-z)                  ! -2.0 -> 2, -2.3 -> 2 (round to nearest)
      write(tag, '(I0)') idx
      if (len_trim(tag) < len(tag)) tag(len_trim(tag)+1:) = ' '
      return
    end if

    ! --- Original behaviour for non-negative values ---

    ! 1) Write with fixed decimals (adjust precision as you like)
    write(tmp, '(F15.8)') z      ! e.g. "      0.25000000"

    ! 2) Left-adjust and trim spaces on the right
    tmp = adjustl(tmp)
    n   = len_trim(tmp)          ! now e.g. "0.25000000"

    ! 3) Replace '.' with 'p'
    do i = 1, n
      if (tmp(i:i) == '.') tmp(i:i) = 'p'
    end do
    ! e.g. "0p25000000"

    ! 4) Remove trailing zeros
    do while (n > 0 .and. tmp(n:n) == '0')
      n = n - 1
    end do
    ! e.g. "0p25"

    ! 5) If we ended up with only integer part ("1p"), drop 'p'
    if (n > 0 .and. tmp(n:n) == 'p') then
      n = n - 1
    end if

    ! 6) Edge case: if everything vanished (z ~ 0), return "0"
    if (n <= 0) then
      tag = '0'
    else
      tag = tmp(1:n)
      if (n < len(tag)) tag(n+1:) = ' '
    end if

  end function real2string

  ! Utility function to convert character to lower case
  pure function to_lower(str) result(out)
    character(*), intent(in) :: str
    character(len(str))      :: out
    integer :: i, c

    do i = 1, len(str)
      c = iachar(str(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) then
        out(i:i) = achar(c + 32)   ! ASCII: 'A'..'Z' → 'a'..'z'
      else
        out(i:i) = str(i:i)
      end if
    end do
  end function to_lower

  ! Utility function to generate linearly spaced array
  pure function linspace(a, b, n) result(x)
    real(kind=8), intent(in) :: a, b
    integer,      intent(in) :: n
    real(kind=8), allocatable :: x(:)
    integer :: i
    if(allocated(x)) deallocate(x)
    if (n <= 0) then
      allocate(x(0))
    else if (n == 1) then
      allocate(x(1)); x = a
    else
      allocate(x(n))
      x = [( a + (b-a)*real(i-1,8)/real(n-1,8), i=1,n )]
    end if
  end function linspace

  ! Utility function to get proper field name
  function field_to_name(field) result(name)
    character(len=1), intent(in) :: field
    character(len=10)            :: name
    select case (field)
    case ('u'); name = 'uVel'
    case ('v'); name = 'vVel'
    case ('w'); name = 'wVel'
    case ('T'); name = 'potT'
    case ('p'); name = 'prss'
    case ('S'); name = 'uVel'
    case default; name = 'uVel'
    end select
  end function field_to_name

  function field_mode(field) result(mode)
    character(len=*), intent(in) :: field
    integer            :: mode
    select case (trim(field))
    case ('u', 'v', 'w', 'T', 'p'); mode = 0
    case ('delta_u', 'delta_v', 'delta_w', &
          'dup_dup','dvp_dvp','dwp_dwp',   &
          'dup_bup','dvp_bvp','dwp_bwp', &
          'delta_p'); mode = 2
    case default; mode = 1
    end select
  end function field_mode

  ! Utility function to get proper budget name
  subroutine define_budget(field, b, t)
    character(*), intent(in) :: field
    character(1), intent(out) :: b
    character(2), intent(out) :: t
    if(trim(field) == 'ubar')then
      b = '0'; t = '01'
    elseif(trim(field) == 'vbar')then
      b = '0'; t = '02'
    elseif(trim(field) == 'wbar')then
      b = '0'; t = '03'
    elseif(trim(field) == 'delta_u')then
      b = '0'; t = '01'
    elseif(trim(field) == 'delta_v')then
      b = '0'; t = '02'
    elseif(trim(field) == 'delta_w')then
      b = '0'; t = '03'
    elseif(trim(field) == 'dup_dup')then
      b = '0'; t = '05'
    elseif(trim(field) == 'dvp_dvp')then
      b = '0'; t = '08'
    elseif(trim(field) == 'dwp_dwp')then
      b = '0'; t = '10'
    elseif(trim(field) == 'dup_bup')then
      b = '0'; t = '11'
    elseif(trim(field) == 'dvp_bvp')then
      b = '0'; t = '16'
    elseif(trim(field) == 'dwp_bwp')then
      b = '0'; t = '19'
    elseif(trim(field) == 'R11')then
      b = '0'; t = '04'
    elseif(trim(field) == 'R22')then
      b = '0'; t = '07'
    elseif(trim(field) == 'R33')then
      b = '0'; t = '09'
    elseif(trim(field) == 'delta_p')then
      b = '0'; t = '04'
    end if
  end subroutine define_budget

  ! Utility function to convert integers or reals to strings
  pure function to_string(i) result(str)
    integer, intent(in) :: i
    character(len=32) :: str
    write(str, '(I0)') i
  end function to_string

  ! Utility subroutine to print messages
  subroutine message(msg)
    character(*), intent(in) :: msg
    write(*,*) trim(msg)
  end subroutine message

  ! Utility subroutine to write z-profile to CSV file
  subroutine csvprofile(nz, filename, z, profile)
    integer, intent(in) :: nz
    character(*), intent(in) :: filename
    real(rk),      intent(in) :: z(nz)
    real(rk),      intent(in) :: profile(nz)
    integer :: uo, n, k, ierr

    n = size(z)
    if (size(profile) /= n) then
      if(myrank == 0) call message('ERROR(CSV): z and profile size mismatch.')
      call MPI_Abort(MPI_COMM_WORLD, 100, ierr)
    end if

    open(newunit=uo, file=trim(filename), status='replace', action='write')
    write(uo,'(A)') 'z,profile'
    do k = 1, n
      write(uo,'(ES23.15, ",", ES23.15)') z(k), profile(k)
    end do
    close(uo)
  end subroutine csvprofile

  !---------------------------------------------------------------------------
  ! Write a 2D slice f(x1,x2) as CSV:
  !   - nx1 rows (i = 1..nx1)
  !   - nx2 columns (j = 1..nx2)
  !   - comma-separated, one row per line
  !---------------------------------------------------------------------------
  subroutine writeslice(nx1, nx2, filename, slice)
    integer,      intent(in) :: nx1, nx2
    character(*), intent(in) :: filename
    real(rk),     intent(in) :: slice(nx1, nx2)

    integer :: i, j, uo, ierr

    open(newunit=uo, file=trim(filename), status='replace', action='write', iostat=ierr)
    if (ierr /= 0) then
      if (myrank == 0) call message('ERROR(writeslice): cannot open file '//trim(filename))
      return
    end if

    do i = 1, nx1
      ! First column
      write(uo, '(ES23.15)', advance='no') slice(i, 1)

      ! Remaining columns with leading commas
      do j = 2, nx2
        write(uo, '(",",ES23.15)', advance='no') slice(i, j)
      end do

      ! End of line
      write(uo, *)   ! advance to next line
    end do

    close(uo)
  end subroutine writeslice

  !---------------------------------------------------------------------------
  ! Given a uniform grid in [0, L] with n points (1-based),
  ! find k0, k1 and alpha such that
  !   f(z0) ≈ (1 - alpha) * f(k0) + alpha * f(k1)
  !
  ! If z0 is outside [0,L], clamp it and set k0 = k1, alpha = 0.
  !---------------------------------------------------------------------------
  subroutine find_bracket_uniform(n, L, z0, k0, k1, alpha)
    integer,  intent(in)  :: n
    real(rk), intent(in)  :: L, z0
    integer,  intent(out) :: k0, k1
    real(rk), intent(out) :: alpha
    real(rk) :: dz, zz, s

    if (n <= 1) then
      k0    = 1
      k1    = 1
      alpha = 0.0_rk
      return
    end if

    if (z0 <= -1)then
      ! A direct index is provided
      k0 = int(abs(z0))
      k1 = k0
      alpha = 0.0_rk
      return
    end if

    dz = L / real(n - 1, rk)

    ! Clamp z0 to [0, L]
    zz = max(0.0_rk, min(L, z0))

    ! Left boundary
    if (zz <= 0.0_rk) then
      k0    = 1
      k1    = 1
      alpha = 0.0_rk
      return
    end if

    ! Right boundary
    if (zz >= L) then
      k0    = n
      k1    = n
      alpha = 0.0_rk
      return
    end if

    ! Now 0 < zz < L, so we’re between grid points
    ! s is between 0 and (n-1)
    s  = zz / dz                ! distance in "index" units from point 1 (0-based)
    k0 = int(floor(s)) + 1      ! 1 <= k0 <= n-1
    k1 = k0 + 1                 ! 2 <= k1 <= n

    ! local coordinate between k0 and k1
    alpha = (zz - dz * real(k0 - 1, rk)) / dz
  end subroutine find_bracket_uniform

  ! Escape single quotes so we can safely single-quote strings in the shell command
  pure function escape_single_quotes(s) result(t)
    character(*), intent(in) :: s
    character(len=:), allocatable :: t
    integer :: i, n, extra, pos
    n = len_trim(s)
    extra = 0
    do i = 1, n
      if (s(i:i) == "'") extra = extra + 3  ! "'" -> '\'' (3 extra chars)
    end do
    t = repeat(' ', n + extra)
    pos = 1
    do i = 1, n
      if (s(i:i) == "'") then
        t(pos:pos) = "'"; pos = pos + 1
        t(pos:pos) = "\"; pos = pos + 1
        t(pos:pos) = "'"; pos = pos + 1
        t(pos:pos) = "'"; pos = pos + 1
      else
        t(pos:pos) = s(i:i); pos = pos + 1
      end if
    end do
    if (pos <= len(t)) t = t(:pos-1)
  end function escape_single_quotes

  ! Split pattern with one '*' into prefix and suffix
  subroutine split_one_star(pattern, prefix, suffix, ok)
    character(*), intent(in)  :: pattern
    character(len=:), allocatable, intent(out) :: prefix, suffix
    logical, intent(out) :: ok
    integer :: p, q, n
    n = len_trim(pattern)
    p = index(pattern(:n), '*')
    if (p == 0) then
      ok = .false.; prefix = ''; suffix = ''; return
    end if
    q = index(pattern(p+1:n), '*')
    if (q /= 0) then
      ok = .false.; prefix = ''; suffix = ''; return
    end if
    prefix = pattern(:p-1)
    suffix = pattern(p+1:n)
    ok = .true.
  end subroutine split_one_star

  ! String utility functions
  logical pure function starts_with(s, pre) result(ok)
    character(*), intent(in) :: s, pre
    integer :: lp
    lp = len_trim(pre)
    if (lp == 0) then
      ok = .true.
    else
      ok = (len_trim(s) >= lp) .and. (s(1:lp) == pre(1:lp))
    end if
  end function starts_with

  logical pure function ends_with(s, suf) result(ok)
    character(*), intent(in) :: s, suf
    integer :: ls, ts
    ls = len_trim(suf); ts = len_trim(s)
    if (ls == 0) then
      ok = .true.
    else
      ok = (ts >= ls) .and. (s(ts-ls+1:ts) == suf(1:ls))
    end if
  end function ends_with

  ! Check if VAL is in LIST
  logical pure function in_list(list, n, val) result(found)
    character(len=*), intent(in) :: list(:)
    integer,          intent(in) :: n
    character(len=*), intent(in) :: val
    integer :: i
    found = .false.
    do i = 1, n
      if (list(i) == val) then
        found = .true.; return
      end if
    end do
  end function in_list

  subroutine list_matching_keys(dir, pattern, keys)
    ! Return unique substrings that fill the '*' in PATTERN.
    ! e.g., pattern "Run01_uVel_t*.out" -> keys like ["0001","0002",...]
    character(*), intent(in) :: dir
    character(*), intent(in) :: pattern
    character(len=:), allocatable, intent(out) :: keys(:)

    character(len=:), allocatable :: pre, suf, d_esc, p_esc, tmpfile, cmd
    character(len=4096) :: line
    integer :: istat, u, nlines, maxlen, klen, ts, lp, ls
    logical :: ok, ex

    ! Default empty result
    allocate(keys(0), mold='     ')

    call split_one_star(pattern, pre, suf, ok)
    if (.not. ok) then
      ! either no '*' or more than one '*'
      return
    end if

    d_esc = escape_single_quotes(trim(dir))
    p_esc = escape_single_quotes(trim(pattern))
    tmpfile = '/tmp/fortran_glob_'//to_string(getpid())//'_keys.txt'

    cmd = "find '"//d_esc//"' -maxdepth 1 -type f -name '"//p_esc//"' -printf '%f\n' > '"//tmpfile//"' 2>/dev/null"
    call execute_command_line(cmd, exitstat=istat)
    if (istat /= 0) return

    inquire(file=tmpfile, exist=ex); if (.not. ex) return

    ! Collect unique keys in a temporary list with conservative max count
    ! First count matches to size a temp array
    nlines = 0
    open(newunit=u, file=tmpfile, status='old', action='read', iostat=istat)
    if (istat /= 0) return
    do
      read(u,'(A)', iostat=istat) line
      if (istat /= 0) exit
      nlines = nlines + 1
    end do
    close(u)
    if (nlines == 0) then
      call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
      return
    end if

    ! Temp store (over-allocated), we’ll dedupe and then resize
    if (allocated(keys)) deallocate(keys)
    allocate(character(len=1024) :: keys(nlines))
    klen = 0; maxlen = 0

    open(newunit=u, file=tmpfile, status='old', action='read', iostat=istat)
    if (istat /= 0) then
      deallocate(keys); allocate(keys(0), mold='     ')
      call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
      return
    end if

    lp = len_trim(pre); ls = len_trim(suf)
    do
      read(u,'(A)', iostat=istat) line
      if (istat /= 0) exit
      ts = len_trim(line)
      if (.not. starts_with(line(:ts), pre) ) cycle
      if (.not. ends_with(  line(:ts), suf) ) cycle
      ! Extract middle between pre and suf
      if (ts < lp + ls) cycle
      ! Middle can be empty -> allow it
      line = line(lp+1 : ts-ls)
      if (.not. in_list(keys, klen, trim(line))) then
        klen = klen + 1
        keys(klen) = trim(line)
        maxlen = max(maxlen, len_trim(line))
      end if
    end do
    close(u)

    call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)

    ! Resize to exactly klen and right element length
    if (klen == 0) then
      deallocate(keys); allocate(keys(0), mold='     ')
    else
      block
        character(len=:), allocatable :: tmp(:)
        integer :: j
        allocate(character(len=maxlen) :: tmp(klen))
        do j = 1, klen
          tmp(j) = adjustl(keys(j)(:maxlen))
        end do
        call move_alloc(tmp, keys)
      end block
    end if
  end subroutine list_matching_keys

  subroutine list_matching_keys_budget(dir, pattern, keys, stamps)
    ! Variant of list_matching_keys to handle files like:
    !   Run06_budget0_term13_t*_n~.s3D
    !
    ! where:
    !   *  -> time stamp (returned in KEYS)
    !   ~  -> 6-digit stamp (returned in STAMPS)
    !
    ! Example filenames:
    !   Run06_budget0_term13_t000900_n123456.s3D
    !   Run06_budget0_term13_t001050_n654321.s3D
    !
    ! Result:
    !   keys   = ["000900","001050",...]
    !   stamps = ["123456","654321",...]
    !
    character(*), intent(in) :: dir
    character(*), intent(in) :: pattern
    character(len=:), allocatable, intent(out) :: keys(:)
    character(len=:), allocatable, intent(out) :: stamps(:)

    character(len=:), allocatable :: pre, suf
    character(len=:), allocatable :: d_esc, p_glob, tmpfile, cmd
    character(len=4096) :: line
    integer :: istat, u, nlines, maxlen_k, maxlen_s, klen
    integer :: ts, lp, pos_n, extpos
    logical :: ok, ex

    ! Default empty result
    allocate(keys(0),   mold='     ')
    allocate(stamps(0), mold='     ')

    ! Split pattern around the single '*' to get prefix PRE (up to 't')
    call split_one_star(pattern, pre, suf, ok)
    if (.not. ok) then
      ! either no '*' or more than one '*'
      return
    end if

    ! Escape directory name
    d_esc  = escape_single_quotes(trim(dir))

    ! Build a glob pattern for 'find':
    !   original:  Run06_budget0_term13_t*_n~.s3D
    !   glob:      Run06_budget0_term13_t*_n*.s3D
    !
    ! i.e. replace '~' with '*' so we ignore the 6-digit stamp in the shell.
    block
      integer :: i, L
      character(len=:), allocatable :: tmp
      L = len_trim(pattern)
      allocate(character(len=L) :: tmp)
      tmp = pattern
      do i = 1, L
        if (tmp(i:i) == '~') tmp(i:i) = '*'
      end do
      p_glob = escape_single_quotes(trim(tmp))
    end block

    tmpfile = '/tmp/fortran_glob_'//to_string(getpid())//'_keys.txt'

    cmd = "find '"//d_esc//"' -maxdepth 1 -type f -name '"//p_glob// &
          "' -printf '%f\n' > '"//tmpfile//"' 2>/dev/null"
    call execute_command_line(cmd, exitstat=istat)
    if (istat /= 0) return

    inquire(file=tmpfile, exist=ex); if (.not. ex) return

    ! Count matches first
    nlines = 0
    open(newunit=u, file=tmpfile, status='old', action='read', iostat=istat)
    if (istat /= 0) return
    do
      read(u,'(A)', iostat=istat) line
      if (istat /= 0) exit
      nlines = nlines + 1
    end do
    close(u)

    if (nlines == 0) then
      call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
      return
    end if

    ! Temp store (over-allocated), we'll dedupe then shrink
    if (allocated(keys))   deallocate(keys)
    if (allocated(stamps)) deallocate(stamps)
    allocate(character(len=1024) :: keys(nlines))
    allocate(character(len=1024) :: stamps(nlines))
    klen      = 0
    maxlen_k  = 0
    maxlen_s  = 0

    open(newunit=u, file=tmpfile, status='old', action='read', iostat=istat)
    if (istat /= 0) then
      deallocate(keys);   allocate(keys(0),   mold='     ')
      deallocate(stamps); allocate(stamps(0), mold='     ')
      call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
      return
    end if

    lp = len_trim(pre)

    do
      read(u,'(A)', iostat=istat) line
      if (istat /= 0) exit
      ts = len_trim(line)
      if (ts <= 0) cycle

      ! Must start with PRE (e.g. "Run06_budget0_term13_t")
      if (.not. starts_with(line(:ts), pre)) cycle

      ! Find the "_n" that comes after the timestamp
      pos_n = index(line(:ts), '_n')
      if (pos_n <= 0) cycle   ! no "_n" -> not our file

      ! Check extension ".s3D"
      if (ts < 4) cycle
      extpos = ts - 3          ! position of '.' in ".s3D"
      if (line(extpos:ts) /= '.s3D') cycle

      ! Extract timestamp between PRE and "_n"
      if (pos_n <= lp+1) cycle   ! nothing between prefix and "_n"
      ! time stamp (*)
      block
        character(len=1024) :: tstamp, sstamp
        integer :: lt, ls

        tstamp = line(lp+1 : pos_n-1)

        ! Extract the 6-digit stamp (~) between "n" and ".s3D"
        ! line: "..._n123456.s3D"
        ! pos_n: index of "_"
        ! 'n' is pos_n+1, stamp starts at pos_n+2, ends at extpos-1
        if (extpos <= pos_n+2) cycle
        sstamp = line(pos_n+2 : extpos-1)

        ! Deduplicate based on time stamp; if same time stamp appears twice
        ! we'll ignore duplicates (assuming 1-to-1 as you said).
        if (.not. in_list(keys, klen, trim(tstamp))) then
          klen = klen + 1
          keys(klen)   = trim(tstamp)
          stamps(klen) = trim(sstamp)
          lt = len_trim(tstamp)
          ls = len_trim(sstamp)
          maxlen_k = max(maxlen_k, lt)
          maxlen_s = max(maxlen_s, ls)
        end if
      end block
    end do

    close(u)
    call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)

    ! Resize KEYS and STAMPS to exactly klen and appropriate lengths
    if (klen == 0) then
      deallocate(keys);   allocate(keys(0),   mold='     ')
      deallocate(stamps); allocate(stamps(0), mold='     ')
    else
      block
        character(len=:), allocatable :: tmpk(:), tmps(:)
        integer :: j

        allocate(character(len=maxlen_k) :: tmpk(klen))
        allocate(character(len=maxlen_s) :: tmps(klen))

        do j = 1, klen
          tmpk(j) = adjustl(keys(j)(:maxlen_k))
          tmps(j) = adjustl(stamps(j)(:maxlen_s))
        end do

        call move_alloc(tmpk, keys)
        call move_alloc(tmps, stamps)
      end block
    end if

  end subroutine list_matching_keys_budget

  subroutine list_matching_filenames(dir, pattern, names)
    ! Return array of filenames (no path) in DIR matching PATTERN (e.g., Run01_uVel_t*.out)
    character(*), intent(in) :: dir
    character(*), intent(in) :: pattern
    character(len=:), allocatable, intent(out) :: names(:)

    character(len=:), allocatable :: d_esc, p_esc, tmpfile, cmd
    integer :: istat, u, nlines, maxlen, i
    character(len=4096) :: line

    ! Default: empty result
    allocate(names(0), mold='     ')

    ! Make a reasonably unique temp file in /tmp (pid-based)
    call random_seed()  ! not needed, but harmless
    tmpfile = '/tmp/fortran_glob_' // trim(adjustl(to_string(getpid()))) // '.txt'

    ! Escape single quotes for safe single-quoting
    d_esc = escape_single_quotes(trim(dir))
    p_esc = escape_single_quotes(trim(pattern))

    ! Use find to list matching files at depth 1; print only the basename (%f)
    cmd = "find '"//d_esc//"' -maxdepth 1 -type f -name '"//p_esc//"' -printf '%f\n' > '"//tmpfile//"' 2>/dev/null"

    call execute_command_line(cmd, exitstat=istat)
    if (istat /= 0) then
      ! Command failed; leave names(:) empty
      return
    end if

    ! First pass: count lines and max length
    nlines = 0; maxlen = 0
    open(newunit=u, file=tmpfile, status='old', action='read', iostat=istat)
    if (istat /= 0) then
      return
    end if
    do
      read(u, '(A)', iostat=istat) line
      if (istat /= 0) exit
      nlines = nlines + 1
      maxlen = max(maxlen, len_trim(line))
    end do
    close(u)

    if (nlines == 0) then
      call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
      return
    end if

    ! Allocate result with a unified length = max filename length
    if (allocated(names)) deallocate(names)
    allocate(character(len=maxlen) :: names(nlines))

    ! Second pass: fill names(:)
    open(newunit=u, file=tmpfile, status='old', action='read', iostat=istat)
    if (istat /= 0) then
      deallocate(names)
      allocate(names(0), mold='     ')
      call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
      return
    end if
    i = 0
    do
      read(u, '(A)', iostat=istat) line
      if (istat /= 0) exit
      i = i + 1
      names(i) = trim(line)
    end do
    close(u)

    ! Clean up temp file
    call execute_command_line("rm -f '"//tmpfile//"'", exitstat=istat)
  end subroutine list_matching_filenames

  ! Small helpers using ISO_C_BINDING to getpid()
  function getpid() result(pid)
    use iso_c_binding, only: c_int
    implicit none
    integer :: pid
    interface
      function c_getpid() bind(C, name="getpid") result(c_pid)
        import :: c_int
        integer(c_int) :: c_pid
      end function c_getpid
    end interface
    pid = int(c_getpid(), kind(pid))
  end function getpid

  function sort_keys_numeric(keys) result(sorted)
    !! Sort character keys by their integer value (ascending).
    !! - keys(:) may have leading zeros; they’re parsed with '(I0)'.
    !! - Non-numeric entries are placed at the end (treated as +huge()).
    character(len=*), intent(in) :: keys(:)
    character(len=:), allocatable :: sorted(:)
    character(len=:), allocatable :: s

    integer :: n, i, j, ios, val
    integer, allocatable :: vals(:), idx(:)
    integer :: maxlen

    n = size(keys)
    maxlen = 0
    do i = 1, n
      maxlen = max(maxlen, len_trim(keys(i)))
    end do

    ! Allocate outputs
    if (n == 0) then
      allocate(character(len=1) :: sorted(0))
      return
    end if
    allocate(vals(n), idx(n))
    allocate(character(len=maxlen) :: sorted(n))

    ! Parse integers; non-numeric => push to the end
    do i = 1, n
      s = trim(keys(i))
      read(s, *, iostat=ios) val
      if (ios == 0) then
        vals(i) = val
      else
        vals(i) = huge(1)    ! put non-numeric after all numeric
      end if
      idx(i) = i
    end do

    ! Indirect in-place sort of idx by vals (simple O(n^2) – fine for modest n)
    do i = 1, n-1
      do j = i+1, n
        if (vals(idx(j)) < vals(idx(i))) then
          call swap(idx(i), idx(j))
        end if
      end do
    end do

    ! Reorder keys according to idx
    do i = 1, n
      sorted(i) = adjustl(keys(idx(i))( : maxlen))
    end do
  end function sort_keys_numeric

  subroutine sort_keys_and_stamps_numeric(keys, stamps, sorted_keys, sorted_stamps)
    !! Sort KEYS (time stamps) by their integer value (ascending),
    !! and apply the same ordering to STAMPS.
    !!
    !! Input:
    !!   keys(:)   - character time stamps, e.g. "000900", "001050"
    !!   stamps(:) - corresponding "~" stamps, e.g. "123456", "654321"
    !!
    !! Output (allocatable):
    !!   sorted_keys(:), sorted_stamps(:) - reordered copies
    !!
    character(len=*), intent(in)  :: keys(:)
    character(len=*), intent(in)  :: stamps(:)
    character(len=:), allocatable, intent(out) :: sorted_keys(:)
    character(len=:), allocatable, intent(out) :: sorted_stamps(:)

    integer :: n, i, j, ios, val
    integer, allocatable :: vals(:), idx(:)
    integer :: maxlen_k, maxlen_s
    character(len=:), allocatable :: s

    ! Basic checks
    n = size(keys)
    if (n == 0 .or. size(stamps) /= n) then
      allocate(character(len=1) :: sorted_keys(0))
      allocate(character(len=1) :: sorted_stamps(0))
      return
    end if

    allocate(vals(n), idx(n))

    ! Parse integers from KEYS; non-numeric => sent to the end
    do i = 1, n
      s = trim(keys(i))
      read(s, *, iostat=ios) val
      if (ios == 0) then
        vals(i) = val
      else
        vals(i) = huge(1)    ! put non-numeric keys after numeric ones
      end if
      idx(i) = i
    end do

    ! Simple O(n^2) indirect sort of idx by vals
    do i = 1, n-1
      do j = i+1, n
        if (vals(idx(j)) < vals(idx(i))) then
          call swap(idx(i), idx(j))   ! your existing swap(int,int)
        end if
      end do
    end do

    ! Decide output lengths
    maxlen_k = 0
    maxlen_s = 0
    do i = 1, n
      maxlen_k = max(maxlen_k, len_trim(keys(i)))
      maxlen_s = max(maxlen_s, len_trim(stamps(i)))
    end do
    if (maxlen_k <= 0) maxlen_k = 1
    if (maxlen_s <= 0) maxlen_s = 1

    ! Allocate outputs with trimmed lengths
    allocate(character(len=maxlen_k) :: sorted_keys(n))
    allocate(character(len=maxlen_s) :: sorted_stamps(n))

    ! Fill outputs according to permutation idx
    do i = 1, n
      sorted_keys(i)   = adjustl(keys(idx(i))(1:maxlen_k))
      sorted_stamps(i) = adjustl(stamps(idx(i))(1:maxlen_s))
    end do

    deallocate(vals, idx)

  end subroutine sort_keys_and_stamps_numeric

  ! Simple integer swap
  pure subroutine swap(a, b)
    integer, intent(inout) :: a, b
    integer :: t
    t = a; a = b; b = t
  end subroutine swap
end module MPIR3D

!===============================================================
!===============================================================
program MPIR3D_
  use MPIR3D
  implicit none

  type(FieldReader2Decomp) :: reader
  integer :: ierr
  character(len=256) :: path, outdir
  character(len=10) :: field
  integer:: nx=1, ny=1, nz=1, runid=1, taskid=0, num_slice=1
  real(rk):: Lx=1.0_rk, Ly=1.0_rk, Lz=1.0_rk
  integer :: nlen, ioUnit=28
  character(:), allocatable :: inputfile
  character(len=1) :: slice_axis = 'z'
  real(rk) :: slice_coord(1000)
  real(rk), allocatable :: slice_coord_(:)
  namelist /SETUP/ nx, ny, nz, Lx, Ly, Lz, path, outdir, runid, taskid, field, &
                   slice_axis, num_slice, slice_coord
      
  ! Initiate MPI
  ! -------------------------------------------------------------------------!
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count() < 1) then
    if(myrank == 0) call message('Usage: MPIR3D <inputfile>')
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  else
    call get_command_argument(1, length=nlen)
  end if
  allocate(character(nlen) :: inputfile)
  call get_command_argument(1, value=inputfile)
  ! -------------------------------------------------------------------------!

  if(myrank == 0) call message('Inputfile is: '//trim(inputfile)) ! Remove later

  ! Read input namelist
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
  read(unit=ioUnit, NML=SETUP)
  close(ioUnit)
  if (taskid == 1)then
    if(allocated(slice_coord_))deallocate(slice_coord_)
    allocate(slice_coord_(num_slice))
    slice_coord_ = slice_coord(1:num_slice)
  end if

  ! Initiate reader
  call reader%init(nx, ny, nz)

  if (taskid == 0)then
    ! Horizontal average
    call ha_driver(reader, Lz, runid, path, outdir, field)
  else if (taskid == 1) then
    ! Slice
    call slice_driver(reader, Lx, Ly, Lz, runid, path, outdir, field, slice_axis, &
      num_slice, slice_coord_)
  end if
  
  if(myrank == 0) call message('Wrapping up ...')
  call MPI_Finalize(ierr)
end program MPIR3D_
