!===============================================================
! MPI_diag.F90
!===============================================================

module MPIR3D
  use mpi
  use netcdf
  use decomp_2d_io
  use decomp_2d, only: decomp_info, decomp_2d_init, xstart, &
                       xend, ystart, yend, zstart, zend, &
                       decomp_2d_finalize, decomp_info_init, &
                       mytype, transpose_x_to_y, transpose_y_to_z
  implicit none

  integer, parameter :: rk = mytype
  integer :: myrank = -1
  integer :: nprocs = -1
  integer :: instfieldsrc = 0
  integer :: timeavgsrc = 1
  integer :: mdgtsrc = 2
  integer :: mbdgtsrc = 3
  integer :: nxloc_,nyloc_,nzloc_

  interface to_string
    module procedure to_string_int
    module procedure to_string_real
  end interface to_string

  type :: FieldReader2Decomp
     private
     ! Global sizes and process grid
     integer :: nx = 0, ny = 0, nz = 0
     integer :: p_row = 0, p_col = 0
     logical :: is_init = .false.

     ! 2DECOMP&FFT descriptor (x-pencil; cell based)
     type(DECOMP_INFO) :: gpC
     
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

  type :: slice_packet_t
     integer :: nslices = 0
     character(len=1) :: axis = ' '
     real(rk), allocatable :: coords(:)
     integer :: integrate = 0
  end type slice_packet_t

  type :: rz_params
    real(rk) :: tm
    real(rk) :: a
    real(rk) :: b
    real(rk) :: l
    real(rk) :: d
    real(rk) :: sse
    integer  :: status
  end type rz_params

contains

  pure real(rk) function f_basis(x)
    real(rk), intent(in) :: x
    f_basis = 0.5_rk * (tanh(x) + 1.0_rk)
  end function f_basis

  pure real(rk) function log_2cosh_stable(x)
    real(rk), intent(in) :: x

    ! log(2 cosh x) = log(exp(x) + exp(-x))
    ! Stable implementation avoiding overflow.
    if (x >= 0.0_rk) then
        log_2cosh_stable = x + log(1.0_rk + exp(-2.0_rk*x))
    else
        log_2cosh_stable = -x + log(1.0_rk + exp(2.0_rk*x))
    end if
  end function log_2cosh_stable

  pure real(rk) function g_basis(x)
    real(rk), intent(in) :: x
    g_basis = 0.5_rk * (log_2cosh_stable(x) + x)
  end function g_basis

  subroutine fit_rz_profile(z, t, n, l0, d0, params, d_min, ridge, max_iter, tol)
    integer,  intent(in)  :: n
    real(rk), intent(in)  :: z(n), t(n)
    real(rk), intent(in)  :: l0, d0
    type(rz_params), intent(out) :: params

    real(rk), intent(in), optional :: d_min
    real(rk), intent(in), optional :: ridge
    integer,  intent(in), optional :: max_iter
    real(rk), intent(in), optional :: tol

    real(rk) :: dmin_loc, ridge_loc, tol_loc
    integer  :: max_iter_loc

    real(rk) :: x(2,3), fval(3)
    real(rk) :: centroid(2), xr(2), xe(2), xc(2)
    real(rk) :: fr, fe, fc
    real(rk) :: alpha, gamma, rho, sigma
    real(rk) :: scale_l, scale_d, z_range, simplex_size
    integer :: iter
    integer :: i_best, i_worst, i_mid
    integer :: stat

    dmin_loc = 1.0e-6_rk
    ridge_loc = 0.0_rk
    max_iter_loc = 500
    tol_loc = 1.0e-8_rk

    if (present(d_min))   dmin_loc = d_min
    if (present(ridge))   ridge_loc = ridge
    if (present(max_iter)) max_iter_loc = max_iter
    if (present(tol))     tol_loc = tol

    alpha = 1.0_rk   ! reflection
    gamma = 2.0_rk   ! expansion
    rho   = 0.5_rk   ! contraction
    sigma = 0.5_rk   ! shrink

    ! Initial simplex in the two nonlinear variables: x = [l, d].
    z_range = maxval(z) - minval(z)
    if (z_range <= 0.0_rk) then
        params%status = -2
        return
    end if
    scale_l = max(10.0_rk * dmin_loc, 0.05_rk * z_range)
    scale_d = max(10.0_rk * dmin_loc, 0.10_rk * z_range)
    
    x(:,1) = [l0, max(d0, dmin_loc)]
    x(:,2) = [l0 + scale_l, max(d0, dmin_loc)]
    x(:,3) = [l0, max(d0 + scale_d, dmin_loc)]

    do i_best = 1, 3
        fval(i_best) = objective(z, t, n, x(1,i_best), x(2,i_best), dmin_loc, ridge_loc, stat)
    end do

    do iter = 1, max_iter_loc

        call order_simplex(fval, i_best, i_mid, i_worst)

        simplex_size = maxval(abs(x - spread(x(:,i_best), 2, 3)))

        if (maxval(abs(fval - fval(i_best))) < tol_loc * (1.0_rk + abs(fval(i_best)))) exit
        if (simplex_size < tol_loc * max(1.0_rk, z_range)) exit

        centroid = 0.5_rk * (x(:,i_best) + x(:,i_mid))

        ! Reflection
        xr = centroid + alpha * (centroid - x(:,i_worst))
        xr(2) = max(xr(2), dmin_loc)
        fr = objective(z, t, n, xr(1), xr(2), dmin_loc, ridge_loc, stat)

        if (fr < fval(i_best)) then
            ! Expansion
            xe = centroid + gamma * (xr - centroid)
            xe(2) = max(xe(2), dmin_loc)
            fe = objective(z, t, n, xe(1), xe(2), dmin_loc, ridge_loc, stat)

            if (fe < fr) then
                x(:,i_worst) = xe
                fval(i_worst) = fe
            else
                x(:,i_worst) = xr
                fval(i_worst) = fr
            end if

        else if (fr < fval(i_mid)) then
            x(:,i_worst) = xr
            fval(i_worst) = fr

        else
            ! Contraction
            if (fr < fval(i_worst)) then
                xc = centroid + rho * (xr - centroid)
            else
                xc = centroid + rho * (x(:,i_worst) - centroid)
            end if

            xc(2) = max(xc(2), dmin_loc)
            fc = objective(z, t, n, xc(1), xc(2), dmin_loc, ridge_loc, stat)

            if (fc < fval(i_worst)) then
                x(:,i_worst) = xc
                fval(i_worst) = fc
            else
                ! Shrink toward best point
                x(:,i_mid)   = x(:,i_best) + sigma * (x(:,i_mid)   - x(:,i_best))
                x(:,i_worst) = x(:,i_best) + sigma * (x(:,i_worst) - x(:,i_best))

                x(2,i_mid)   = max(x(2,i_mid), dmin_loc)
                x(2,i_worst) = max(x(2,i_worst), dmin_loc)

                fval(i_mid) = objective(z, t, n, x(1,i_mid), x(2,i_mid), dmin_loc, ridge_loc, stat)
                fval(i_worst) = objective(z, t, n, x(1,i_worst), x(2,i_worst), dmin_loc, ridge_loc, stat)
            end if
        end if
    end do

    call order_simplex(fval, i_best, i_mid, i_worst)

    params%l = x(1,i_best)
    params%d = x(2,i_best)
    params%sse = fval(i_best)

    call solve_tmab(z, t, n, params%l, params%d, ridge_loc, params%tm, params%a, params%b, stat)

    params%status = stat
    if (iter > max_iter_loc) params%status = 1
  end subroutine fit_rz_profile

  real(rk) function objective(z, t, n, l, d, d_min, ridge, status)
    integer,  intent(in)  :: n
    real(rk), intent(in)  :: z(n), t(n)
    real(rk), intent(in)  :: l, d, d_min, ridge
    integer,  intent(out) :: status

    real(rk) :: tm, a, b
    real(rk) :: eta, t_fit, res
    integer :: i

    if (d <= d_min) then
        objective = huge(1.0_rk)
        status = -1
        return
    end if

    call solve_tmab(z, t, n, l, d, ridge, tm, a, b, status)

    if (status /= 0) then
        objective = huge(1.0_rk)
        return
    end if

    objective = 0.0_rk

    do i = 1, n
        eta = (z(i) - l) / d
        t_fit = tm + a * f_basis(eta) + b * g_basis(eta)
        res = t(i) - t_fit
        objective = objective + res * res
    end do
  end function objective

  subroutine solve_tmab(z, t, n, l, d, ridge, tm, a, b, status)
    integer,  intent(in)  :: n
    real(rk), intent(in)  :: z(n), t(n)
    real(rk), intent(in)  :: l, d, ridge
    real(rk), intent(out) :: tm, a, b
    integer,  intent(out) :: status

    real(rk) :: M(3,3), rhs(3), sol(3)
    real(rk) :: eta, fv, gv
    real(rk) :: S1, Sf, Sg, Sff, Sgg, Sfg
    real(rk) :: St, Sft, Sgt
    integer :: i

    S1  = real(n, rk)
    Sf  = 0.0_rk
    Sg  = 0.0_rk
    Sff = 0.0_rk
    Sgg = 0.0_rk
    Sfg = 0.0_rk
    St  = 0.0_rk
    Sft = 0.0_rk
    Sgt = 0.0_rk

    do i = 1, n
        eta = (z(i) - l) / d
        fv = f_basis(eta)
        gv = g_basis(eta)

        Sf  = Sf  + fv
        Sg  = Sg  + gv
        Sff = Sff + fv * fv
        Sgg = Sgg + gv * gv
        Sfg = Sfg + fv * gv

        St  = St  + t(i)
        Sft = Sft + fv * t(i)
        Sgt = Sgt + gv * t(i)
    end do

    M(1,:) = [S1, Sf,  Sg]
    M(2,:) = [Sf, Sff, Sfg]
    M(3,:) = [Sg, Sfg, Sgg]

    if (ridge > 0.0_rk) then
        M(1,1) = M(1,1) + ridge
        M(2,2) = M(2,2) + ridge
        M(3,3) = M(3,3) + ridge
    end if

    rhs = [St, Sft, Sgt]

    call solve_3x3(M, rhs, sol, status)

    if (status == 0) then
        tm = sol(1)
        a  = sol(2)
        b  = sol(3)
    else
        tm = 0.0_rk
        a  = 0.0_rk
        b  = 0.0_rk
    end if
  end subroutine solve_tmab

  subroutine solve_3x3(Ain, bin, x, status)
    real(rk), intent(in)  :: Ain(3,3), bin(3)
    real(rk), intent(out) :: x(3)
    integer,  intent(out) :: status

    real(rk) :: A(3,3), b(3)
    real(rk) :: factor, tmp, pivot_abs
    integer :: i, j, k, p

    A = Ain
    b = bin
    status = 0

    do k = 1, 2
        p = k
        pivot_abs = abs(A(k,k))

        do i = k+1, 3
            if (abs(A(i,k)) > pivot_abs) then
                p = i
                pivot_abs = abs(A(i,k))
            end if
        end do

        if (pivot_abs < 1.0e-14_rk) then
            status = -1
            x = 0.0_rk
            return
        end if

        if (p /= k) then
            do j = k, 3
                tmp = A(k,j)
                A(k,j) = A(p,j)
                A(p,j) = tmp
            end do

            tmp = b(k)
            b(k) = b(p)
            b(p) = tmp
        end if

        do i = k+1, 3
            factor = A(i,k) / A(k,k)
            A(i,k) = 0.0_rk

            do j = k+1, 3
                A(i,j) = A(i,j) - factor * A(k,j)
            end do

            b(i) = b(i) - factor * b(k)
        end do
    end do

    if (abs(A(3,3)) < 1.0e-14_rk) then
        status = -1
        x = 0.0_rk
        return
    end if

    x(3) = b(3) / A(3,3)
    x(2) = (b(2) - A(2,3) * x(3)) / A(2,2)
    x(1) = (b(1) - A(1,2) * x(2) - A(1,3) * x(3)) / A(1,1)
  end subroutine solve_3x3

  subroutine order_simplex(fval, i_best, i_mid, i_worst)
    real(rk), intent(in)  :: fval(3)
    integer,  intent(out) :: i_best, i_mid, i_worst

    integer :: idx(3), i, j, tmp

    idx = [1, 2, 3]

    do i = 1, 2
        do j = i+1, 3
            if (fval(idx(j)) < fval(idx(i))) then
                tmp = idx(i)
                idx(i) = idx(j)
                idx(j) = tmp
            end if
        end do
    end do

    i_best  = idx(1)
    i_mid   = idx(2)
    i_worst = idx(3)
  end subroutine order_simplex

  pure real(rk) function eval_rz_profile(z, params)
    real(rk), intent(in) :: z
    type(rz_params), intent(in) :: params

    real(rk) :: eta

    eta = (z - params%l) / params%d
    eval_rz_profile = params%tm + params%a * f_basis(eta) + params%b * g_basis(eta)
  end function eval_rz_profile

  subroutine read_slice_map(filename, packets, ierr)
    implicit none

    character(len=*), intent(in) :: filename
    type(slice_packet_t), allocatable, intent(out) :: packets(:)
    integer, intent(out), optional :: ierr

    integer :: unit
    integer :: ios
    integer :: npackets
    integer :: ip
    character(len=2048) :: line

    if (present(ierr)) ierr = 0

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       if (present(ierr)) then
          ierr = ios
          return
       else
          error stop "Could not open slice-map file."
       end if
    end if

    ! First line: number of slice packets
    call read_nonempty_line(unit, line, ios)
    if (ios /= 0) call fail(unit, ios, ierr, "Could not read number of packets.")

    read(line, *, iostat=ios) npackets
    if (ios /= 0) call fail(unit, ios, ierr, "Invalid number of packets.")

    allocate(packets(npackets))

    do ip = 1, npackets

       ! First line in packet: number of slices
       call read_nonempty_line(unit, line, ios)
       if (ios /= 0) call fail(unit, ios, ierr, "Could not read number of slices.")

       read(line, *, iostat=ios) packets(ip)%nslices
       if (ios /= 0) call fail(unit, ios, ierr, "Invalid number of slices.")

       allocate(packets(ip)%coords(packets(ip)%nslices))

       ! Second line in packet: axis, e.g. "x", "y", "z"
       call read_nonempty_line(unit, line, ios)
       if (ios /= 0) call fail(unit, ios, ierr, "Could not read slice axis.")

       packets(ip)%axis = parse_axis(line)

       if (.not. any(packets(ip)%axis == ['x', 'y', 'z'])) then
          call fail(unit, -1, ierr, "Invalid slice axis.")
       end if

       ! Third line in packet: coordinate list
       call read_nonempty_line(unit, line, ios)
       if (ios /= 0) call fail(unit, ios, ierr, "Could not read coordinate list.")

       call parse_real_list(line, packets(ip)%coords, ios)
       if (ios /= 0) call fail(unit, ios, ierr, "Invalid coordinate list.")

       ! Fourth line in packet: integration flag
       call read_nonempty_line(unit, line, ios)
       if (ios /= 0) call fail(unit, ios, ierr, "Could not read integration flag.")

       read(line, *, iostat=ios) packets(ip)%integrate
       if (ios /= 0) call fail(unit, ios, ierr, "Invalid integration flag.")

       if (.not. any(packets(ip)%integrate == [0, 1])) then
          call fail(unit, -2, ierr, "Integration flag must be 0 or 1.")
       end if

    end do

    close(unit)

  end subroutine read_slice_map

  subroutine read_nonempty_line(unit, line, ios)
    implicit none

    integer, intent(in) :: unit
    character(len=*), intent(out) :: line
    integer, intent(out) :: ios

    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) return
       if (len_trim(line) > 0) return
    end do

  end subroutine read_nonempty_line

  function parse_axis(line) result(axis)
    implicit none

    character(len=*), intent(in) :: line
    character(len=1) :: axis
    character(len=:), allocatable :: tmp
    integer :: i

    tmp = adjustl(trim(line))

    ! Remove quotes if present.
    do i = 1, len_trim(tmp)
       if (tmp(i:i) /= '"' .and. tmp(i:i) /= "'") then
          axis = tmp(i:i)
          return
       end if
    end do

    axis = ' '

  end function parse_axis

  subroutine parse_real_list(line, values, ios)
    implicit none

    character(len=*), intent(in) :: line
    real(rk), intent(out) :: values(:)
    integer, intent(out) :: ios

    character(len=:), allocatable :: tmp
    integer :: i

    tmp = trim(line)

    ! Convert comma-separated list to whitespace-separated list.
    do i = 1, len(tmp)
       if (tmp(i:i) == ',') tmp(i:i) = ' '
    end do

    read(tmp, *, iostat=ios) values

  end subroutine parse_real_list


  subroutine fail(unit, code, ierr, message)
    implicit none

    integer, intent(in) :: unit
    integer, intent(in) :: code
    integer, intent(out), optional :: ierr
    character(len=*), intent(in) :: message

    close(unit)

    if (present(ierr)) then
       ierr = code
       return
    else
       print *, trim(message)
       error stop
    end if

  end subroutine fail

  logical function within_range(istart, iend, tidx)
      implicit none
      character(*), intent(in) :: tidx
      integer, intent(in) :: istart, iend
      integer :: itime
      integer :: ios

      read(tidx, '(I6)', iostat=ios) itime
      if (ios /= 0) then
         within_range = .false.
         return
      end if
      within_range = (itime >= istart .and. itime <= iend)
   end function within_range

  ! Initialize the FieldReader2Decomp
  subroutine frd_init(this, nx, ny, nz)
    class(FieldReader2Decomp), intent(inout) :: this
    integer,               intent(in)    :: nx, ny, nz
    integer :: ierr
    character(256) :: buf
    integer :: loc(6)
    integer, allocatable :: recvbuf(:)
    !integer :: r, off

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
    call decomp_info_init(this%nx, this%ny, this%nz, this%gpC)
    
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
    ! if (myrank == 0) then      
    !   do r = 0, nprocs-1
    !     off = 6*r
    !     write(buf,'("FRD_INIT: Rank=",I0," X: ",I0,", ",I0," // Y: ",I0,", ",I0," // Z: ",I0,", ",I0)') &
    !         r, recvbuf(off+1), recvbuf(off+2), recvbuf(off+3), recvbuf(off+4), recvbuf(off+5), recvbuf(off+6)
    !     call message(trim(buf))
    !   end do
    !   deallocate(recvbuf)
    ! end if
    
    this%is_init = .true.

    if (myrank == 0) then
      call message('FRD_INIT: 2DECOMP reader initialized.')
      call message(' ')
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
    call decomp_2d_read_one(1, field, path, this%gpC)
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
  
  subroutine ha_driver(reader, Lz, runid, path, outdir, field, budget_source, start_idx, end_idx, filename)
    class(FieldReader2Decomp), intent(inout) :: reader
    integer, intent(in) :: runid
    integer, intent(in) :: start_idx, end_idx
    character(*), intent(in), optional :: filename
    character(*), intent(in) :: path, outdir, field
    integer,          intent(in)  :: budget_source
    real(rk), intent(in) :: Lz    
    real(rk), allocatable :: f1(:,:,:), x(:), y(:), z(:)
    character(len=:), allocatable :: sorted_keys(:), sorted_stamps(:)
    character(len=2) :: rc
    character(len=256):: f_
    integer :: k, nx, ny, nz, xs, xe, ys, ye, zs, ze, nxloc, nyloc, nzloc
    character(len=1024) :: filename_, outname
    logical :: break=.false.

    if(present(filename)) then
      filename_ = trim(filename)
    else
      filename_ = 'null'
    end if

    ! Fetch local and global sizes & start indices
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%global_shape(nx, ny, nz)
    call reader%indices(xs, xe, ys, ye, zs, ze)

    ! allocate the data holder    
    allocate(f1(nxloc, nyloc, nzloc))
    
    ! Convert runid to character
    write(rc, '(I2.2)') runid

    ! Vertical coordinate z
    allocate(x(nx), y(ny), z(nz))
    call create_grid(1._rk, 1._rk, Lz, nx, ny, nz, x, y, z)

    ! Get file list and sort by time
    call get_keys_stamps(trim(path), trim(rc), budget_source, trim(field), f_, sorted_keys, sorted_stamps)
    
    do k = 1, size(sorted_keys)

      if(.not. within_range(start_idx, end_idx, trim(sorted_keys(k)))) cycle

      if(break) exit  ! If we read from file, no need to loop over time snapshots

      if(trim(filename_) == 'null')then
        f1 = eval_field(trim(field), reader, budget_source, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
      else
        f1 = reader%read_field(trim(path)//'/'//trim(filename_))
        break = .true.
      end if
      
      if ((field == 'S') .and. (trim(filename_) == 'null')) then
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
          if(trim(filename_) == 'null')then
            outname = trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//&
                '_HA_'//trim(f_)//'.csv'
          else
            outname = trim(outdir)//'/'//trim(filename_)//'_HA_profile.csv'
          end if
          call csvprofile(nz, trim(outname) ,z,profile)
        end block
      end if
    end do
  end subroutine ha_driver

  function strip_extension(filename) result(basename)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=:), allocatable :: basename
    integer :: i, last_dot, n

    n = len_trim(filename)
    last_dot = 0

    ! Find last dot
    do i = n, 1, -1
      if (filename(i:i) == '.') then
        last_dot = i
        exit
      end if
    end do

    ! Case 1: no dot found OR dot is first character (hidden file)
    if (last_dot == 0 .or. last_dot == 1) then
      basename = filename(1:n)
    else
      basename = filename(1:last_dot-1)
    end if
  end function strip_extension

  subroutine max_time_change(reader, runid, path, field, budget_source, start_idx, end_idx)
    implicit none
    class(FieldReader2Decomp), intent(inout) :: reader
    integer,          intent(in)  :: runid
    integer,          intent(in)  :: start_idx, end_idx
    character(*),     intent(in)  :: path, field
    integer,          intent(in)  :: budget_source

    real(rk), allocatable :: f1(:,:,:), f2(:,:,:)

    integer :: nx, ny, nz
    integer :: nxloc, nyloc, nzloc
    integer :: xs, xe, ys, ye, zs, ze
    integer :: num_stamps
    integer :: k, ierr

    character(len=:), allocatable :: sorted_keys(:), sorted_stamps(:)
    character(len=256) :: f_
    character(len=16)  :: rc

    logical :: first

    first = .true.

    write(rc, '(I2.2)') runid

    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%global_shape(nx, ny, nz)
    call reader%indices(xs, xe, ys, ye, zs, ze)

    allocate(f1(nxloc, nyloc, nzloc))
    allocate(f2(nxloc, nyloc, nzloc))

    call get_keys_stamps(trim(path), trim(rc), budget_source, trim(field), f_, sorted_keys, sorted_stamps)
    num_stamps = size(sorted_keys)

    timeloop: do k = 1, num_stamps
      if (.not. within_range(start_idx, end_idx, trim(sorted_keys(k)))) cycle

      if (.not. first) f2 = f1

      f1 = eval_field(trim(field), reader, budget_source, trim(path), trim(rc), trim(rc), &
                      trim(sorted_keys(k)), trim(sorted_stamps(k)))

      if (first) then
        first = .false.
        cycle
      end if

      block
        real(rk) :: local_max, global_max
        integer  :: i, j, m

        local_max = 0.0_rk
        do m = 1, nzloc
          do j = 1, nyloc
            do i = 1, nxloc
              local_max = max(local_max, abs(f1(i,j,m) - f2(i,j,m)))
            end do
          end do
        end do

        call MPI_Allreduce(local_max, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

        if (myrank == 0) then
            write(*,'("Time: ",A," Max change: ",ES12.5)') trim(sorted_stamps(k)), global_max
        end if
      end block
    end do timeloop

    if (allocated(f1)) deallocate(f1)
    if (allocated(f2)) deallocate(f2)
  end subroutine max_time_change

  subroutine read_field_list(filename, field_list, ierr)
    implicit none

    character(len=*), intent(in) :: filename
    character(len=100), allocatable, intent(out) :: field_list(:)
    integer, intent(out), optional :: ierr

    integer :: unit
    integer :: ios
    integer :: nfields
    integer :: i
    character(len=2048) :: line

    if (present(ierr)) ierr = 0

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      if (present(ierr)) then
          ierr = ios
          return
      else
          error stop "Could not open field-list file."
      end if
    end if

    ! First pass: count nonempty lines.
    nfields = 0

    do
      call read_nonempty_line(unit, line, ios)
      if (ios /= 0) exit

      nfields = nfields + 1
    end do

    if (ios > 0) then
      close(unit)
      if (present(ierr)) then
          ierr = ios
          return
      else
          error stop "Error while reading field-list file."
      end if
    end if

    allocate(field_list(nfields))

    ! Second pass: read field names.
    rewind(unit)

    i = 0

    do
      call read_nonempty_line(unit, line, ios)
      if (ios /= 0) exit

      i = i + 1
      field_list(i) = adjustl(trim(line))
    end do

    close(unit)

  end subroutine read_field_list

  logical function has_s3d_extension(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: n

    n = len_trim(filename)

    has_s3d_extension = .false.

    if (n >= 4) then
        has_s3d_extension = filename(n-3:n) == ".s3D"
    end if
  end function has_s3d_extension

  subroutine slice_driver(reader, Lx, Ly, Lz, runid, path, outdir, fieldfile, slice_map, budget_source, start_idx, end_idx, filename)
    class(FieldReader2Decomp), intent(inout) :: reader
    integer,          intent(in)  :: runid, start_idx, end_idx
    character(*),     intent(in)  :: path, outdir, fieldfile, slice_map, filename
    real(rk),         intent(in)  :: Lx, Ly, Lz
    integer,          intent(in)  :: budget_source

    ! Local 3D field
    real(rk), allocatable, target :: f1(:,:,:)

    ! Grid and indices
    integer :: nx, ny, nz
    integer :: nxloc, nyloc, nzloc
    integer :: xs, xe, ys, ye, zs, ze
    
    ! Slice indices and interpolation
    integer :: k, k0, k1
    real(rk) :: alpha
    character(len=2) :: rc
    character(len=256) :: f_, field
    character(len=1) :: ax, x1name, x2name, eax
    character(len=256) :: fname, msg

    ! Time keys
    character(len=:), allocatable :: sorted_keys(:), sorted_stamps(:)
    integer :: num_stamps

    ! Slice arrays (global shape on every rank)
    real(rk), allocatable :: local_slice0(:,:), local_slice1(:,:)
    real(rk), allocatable :: global_slice0(:,:), global_slice1(:,:), slice_interp(:,:)
    real(rk), pointer :: x1(:), x2(:)
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), pointer :: slice_ptr(:,:)
    real(rk) :: delta

    integer :: i, j, isl, ierr, nx1, nx2, nax, jloc
    integer :: x1s, x1e, x2s, x2e, axs, axe
    real(rk) :: L1, L2, Lax, slice_
    logical :: filemode = .false.
    character(len=2) :: term
    character(len=1) :: budget
    logical :: calc=.false., integrate
    character(len=100), allocatable :: field_list(:)
    integer :: nfields=0, ifield, islice
    type(slice_packet_t), allocatable :: slice_packets(:)
    integer :: nslice

    filemode = .not.(trim(filename) == 'null')
    if(myrank == 0) then
      if(filemode)then
        call message('file mode is true')
      else
        call message('file mode is false')
      end if
    end if

    ! Read list of fields to be processed if not in file mode
    if(filemode)then
      nfields = 1
    else
      call read_field_list(trim(fieldfile), field_list, ierr)
      if (ierr /= 0) then
        print *, "Error reading field list. ierr = ", ierr
        call MPI_ABORT(MPI_COMM_WORLD, 1123, ierr)
      end if
      nfields=size(field_list)
    end if

    ! Read slice map to get axis and slice coordinates
    call read_slice_map(trim(slice_map), slice_packets, ierr)
    if (ierr /= 0) then
      print *, "Error reading slice map. ierr = ", ierr
      call MPI_ABORT(MPI_COMM_WORLD, 1124, ierr)
    end if

    ! Convert runid to character
    write(rc, '(I2.2)') runid

    ! Shapes and local indices
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%global_shape(nx, ny, nz)
    call reader%indices(xs, xe, ys, ye, zs, ze)

    ! Allocate local 3D field
    allocate(f1(nxloc, nyloc, nzloc))

    ! Get file list and sort by time
    if(.not. filemode)then
      call get_keys_stamps(trim(path), trim(rc), budget_source, trim(field_list(1)), f_, sorted_keys, sorted_stamps)
      num_stamps = size(sorted_keys)
    else
      num_stamps = 1
    end if

    ! Mesh
    allocate(x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)

      ! Loop over files first to avoid redundant multiple reads of the same file
      timeloop: do k = 1, num_stamps
        fields_loop : do ifield = 1, nfields

        ! Either read a field or directly from a filename
        if(filemode)then
          f1 = reader%read_field(trim(path)//'/'//trim(filename))
          field = "filebased"
        else
          field = trim(field_list(ifield))
          if(.not. within_range(start_idx, end_idx, trim(sorted_keys(k)))) cycle
          
          if(myrank == 0) then
            write(msg,'(A,I0,A,I0,A,A,A,A)') 'Reading field ', ifield, '/', nfields, &
              ': ', trim(field), ' at time ', trim(sorted_stamps(k))
            call message(trim(msg))
          end if
          if(has_s3d_extension(trim(field)))then
            ! It is a file name. Better to have a unified approach through an external map file only.
            f1 = reader%read_field(trim(path)//'/'//trim(field))
          else
            f1 = eval_field(trim(field), reader, budget_source, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          end if
        end if

        slicepackets_loop: do islice = 1, size(slice_packets)
          ax = slice_packets(islice)%axis
          integrate = (slice_packets(islice)%integrate == 1)
          nslice = size(slice_packets(islice)%coords)

          if (myrank == 0) then
            write(msg,'(A,I0,A,I0,A,A)') 'Processing slice packet ', islice, '/', &
              size(slice_packets), ' along axis ', ax
            call message(trim(msg))
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
            x1 => y
            x2 => z
          case ('y')
            x1name = 'x'; x2name = 'z'
            nx1 = nx; nx2 = nz; nax = ny
            L1 = Lx; L2 = Lz; Lax = Ly
            x1s = xs; x1e = xe; x2s = zs; x2e = ze
            axs = ys; axe = ye
            eax = 'j'
            x1 => x
            x2 => z
          case ('z')
            x1name = 'x'; x2name = 'y' 
            nx1 = nx; nx2 = ny; nax = nz
            L1 = Lx; L2 = Ly; Lax = Lz
            x1s = xs; x1e = xe; x2s = ys; x2e = ye
            axs = zs; axe = ze
            eax = 'k'
            x1 => x
            x2 => y
          end select
          delta = Lax/real(nax,rk)

          ! Global slice arrays (same shape on all ranks)
          if(allocated(local_slice0)) deallocate(local_slice0)
          allocate(local_slice0(nx1, nx2))
          if(allocated(global_slice0)) deallocate(global_slice0)
          allocate(global_slice0(nx1, nx2))
          if(.not. integrate)then
            if(allocated(local_slice1)) deallocate(local_slice1)
            allocate(local_slice1(nx1, nx2))
            if(allocated(global_slice1)) deallocate(global_slice1)
            allocate(global_slice1(nx1, nx2))
            if(allocated(slice_interp)) deallocate(slice_interp)
            allocate(slice_interp(nx1, nx2))
          end if

          mode: if(integrate)then
            ! Wipe clean slice arrays
            local_slice0 = 0.0_rk
            global_slice0 = 0.0_rk

            ! Loop over the axis we integrate along
            ! j is global index
            do j=axs, axe
              jloc = j - axs + 1
              select case(ax)
              case('x')
                slice_ptr => f1(jloc,:,:)
              case('y')
                slice_ptr => f1(:,jloc,:)
              case('z')
                slice_ptr => f1(:,:,jloc)
              end select
              local_slice0(x1s:x1e, x2s:x2e) = local_slice0(x1s:x1e, x2s:x2e) + slice_ptr(:,:)*delta
            end do

            ! 3) Sum contributions from all ranks to get full global slices
            call MPI_Reduce(local_slice0, global_slice0, nx1*nx2, MPI_DOUBLE_PRECISION, &
              MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            if (myrank == 0) then
              ! Output file name
              if(.not. filemode) then
                if(has_s3d_extension(trim(field))) then
                  fname = trim(outdir)//'/'//trim(strip_extension(trim(field)))
                else
                  call define_budget(trim(field), budget, term, budget_source, calc)
                  call create_file_name(budget_source, rc, trim(field), budget, term, trim(sorted_keys(k)), trim(sorted_stamps(k)), fname)
                  fname = trim(outdir)//'/'//trim(strip_extension(trim(fname)))
                end if
              else
                fname = trim(outdir)//'/'//trim(strip_extension(trim(filename)))    
              end if
              fname = trim(fname)//'_integ_'//trim(ax)//'.nc'
              
              call message("Exporting to: "//trim(fname))
              call export_slice_to_netcdf(trim(fname), trim(field), global_slice0, x1, x2, x1name, x2name)
            end if
          
          else

            sliceloop: do isl = 1, nslice
              slice_ = slice_packets(islice)%coords(isl)
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
              call MPI_Reduce(local_slice0, global_slice0, nx1*nx2, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierr)
              call MPI_Reduce(local_slice1, global_slice1, nx1*nx2, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierr)

              ! 4) On root: do linear interpolation
              if (myrank == 0) then
                do j = 1, nx2
                  do i = 1, nx1
                    slice_interp(i,j) = (1.0_rk - alpha) * global_slice0(i,j) + alpha * global_slice1(i,j)
                  end do
                end do

                ! Output file name
                if(.not. filemode) then
                  if(has_s3d_extension(trim(field))) then
                    fname = trim(outdir)//'/'//trim(strip_extension(trim(field)))//'_SL_'//ax 
                  else
                    call define_budget(trim(field), budget, term, budget_source, calc)
                    call create_file_name(budget_source, rc, trim(field), budget, term, trim(sorted_keys(k)), trim(sorted_stamps(k)), fname)
                    fname = trim(outdir)//'/'//trim(strip_extension(trim(fname)))//&
                            '_SL_'//ax
                  end if
                else
                  fname = trim(outdir)//'/'//trim(strip_extension(filename))//'_SL_'//ax            
                end if
                if(slice_ <= -1) fname = trim(fname)//'_'//eax ! A direct index is given
                fname = trim(fname)//'='//trim(real2string(slice_))//'.nc'      
                
                call message("Exporting to: "//trim(fname))
                call export_slice_to_netcdf(trim(fname), trim(field), slice_interp, x1, x2, x1name, x2name)

              end if
            end do sliceloop
            
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            if(myrank == 0) call message(' ')
          end if mode

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        if(myrank == 0) call message(' ')

        end do slicepackets_loop
      end do fields_loop
    end do timeloop

    ! Cleanup
    if (allocated(f1))           deallocate(f1)
    if (allocated(local_slice0)) deallocate(local_slice0)
    if (allocated(local_slice1)) deallocate(local_slice1)
    if (allocated(global_slice0)) deallocate(global_slice0)
    if (allocated(global_slice1)) deallocate(global_slice1)
    if (allocated(slice_interp))  deallocate(slice_interp)
  end subroutine slice_driver

  subroutine one_d_profile(reader, Lx, Ly, Lz, runid, budget_source, path, outdir, &
        field, start_idx, end_idx, filename, axis, x1, x2, y1, y2, z1, z2)
    implicit none
    class(FieldReader2Decomp), intent(inout) :: reader
    integer,          intent(in)  :: runid
    integer,          intent(in)  :: start_idx, end_idx, budget_source
    real(rk),         intent(in)  :: Lx, Ly, Lz, x1, x2, y1, y2, z1, z2
    character(*),     intent(in)  :: path, outdir, filename, field
    character(1),     intent(in)  :: axis  
    character(1) :: ax
    character(len=2) :: rc
    integer :: nx, ny, nz
    integer :: nxloc, nyloc, nzloc
    integer :: xs, xe, ys, ye, zs, ze  
    character(len=:), allocatable :: sorted_keys(:), sorted_stamps(:)
    real(rk), allocatable, target :: f(:,:,:)
    real(rk), allocatable :: mask(:,:,:), local_sum(:), global_sum(:)
    real(rk), dimension(:,:), pointer :: f_slice
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), pointer :: xi(:)
    character(len=512) :: f_
    logical :: filemode = .false.
    integer :: num_stamps, ni, i, j, k, kstart, kend, ierr, istamp
    character(len=512) :: outname
    real(rk) :: dx, dy, dz, x_, y_, z_
    integer :: kloc, kglob, nloc
    
    filemode = .not.(trim(filename) == 'null')

    ! Convert runid to character
    write(rc, '(I2.2)') runid
    ax = to_lower(axis(1:1))
    
    ! Shapes and local indices
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%global_shape(nx, ny, nz)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    
    ! Allocate local 3D field
    allocate(f(nxloc, nyloc, nzloc))
    allocate(mask(nxloc, nyloc, nzloc))

    ! Get file list and sort by time
    if(.not. filemode)then
      call get_keys_stamps(trim(path), trim(rc), budget_source, trim(field), f_, sorted_keys, sorted_stamps)
      num_stamps = size(sorted_keys)
    else
      num_stamps = 1
    end if

    allocate(x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)

    select case(ax)
    case('x')
      xi => x
      ni=nx
      nloc = nxloc
      kstart = xs; kend=xe
    case('y')
      xi => y
      ni=ny
      nloc = nyloc
      kstart = ys; kend=ye
    case('z')
      xi => z
      ni=nz
      nloc = nzloc
      kstart = zs; kend=ze
    case default
      call MPI_ABORT(MPI_COMM_WORLD, 322, ierr)
    end select
    allocate(local_sum(ni))
    allocate(global_sum(ni))

    dx = Lx/real(nx,rk)
    dy = Ly/real(ny,rk)
    dz = Lz/real(nz,rk)

    do k=1,nzloc
      do j=1,nyloc
        do i=1,nxloc
            x_ = real( xs + i - 1, rk) * dx - dx
            y_ = real( ys + j - 1, rk ) * dy - dy
            z_ = real( zs + k - 1, rk ) * dz + dz/2.0_rk - dz

            if(((x1 <= x_ ).and.(x_ <= x2)) .and. &
               ((y1 <= y_ ).and.(y_ <= y2)) .and. &
               ((z1 <= z_ ).and.(z_ <= z2)) )then
              mask(i,j,k) = 1._rk
            else
              mask(i,j,k) = 0._rk
            end if
        end do
      end do
    end do

    ! Loop over time snapshots
    do istamp = 1, num_stamps

      local_sum = 0.0_rk
      global_sum = 0.0_rk

      if(filemode)then
        f = reader%read_field(trim(path)//'/'//trim(filename))
      else
        if(.not. within_range(start_idx, end_idx, trim(sorted_keys(istamp)))) cycle
        f = eval_field(trim(field), reader, budget_source, trim(path), trim(rc), trim(rc), trim(sorted_keys(istamp)), trim(sorted_stamps(istamp)))
      end if

      f = f * mask
      
      do kloc = 1, nloc
        kglob = kstart + kloc - 1

        select case(ax)
        case('x')
          f_slice => f(kloc,:,:)
        case('y')
          f_slice => f(:,kloc,:)
        case('z')
          f_slice => f(:,:,kloc)
        case default
          f_slice => null()
        end select

        local_sum(kglob) = SUM(f_slice)
      end do

      call MPI_Allreduce(local_sum, global_sum, ni, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      
      if(myrank == 0)then
        if(filemode)then
          outname = trim(outdir)//'/'//trim(strip_extension(filename))//'_'//ax//'_profile.csv'
        else
          outname = trim(outdir)//'/Run'//trim(rc)//'_t'//trim(sorted_keys(istamp))//'_'//trim(ax)//'_profile_'//trim(f_)//'.csv'
        end if
        call csvprofile(ni, trim(outname) ,xi, global_sum)
      end if
    end do
  end subroutine

  subroutine create_grid(Lx, Ly, Lz, nx, ny, nz, xaxis, yaxis, zaxis)
    implicit none
    real(rk), intent(in) :: Lx, Ly, Lz
    integer, intent(in) :: nx, ny, nz
    real(rk), dimension(*), intent(out) :: xaxis, yaxis, zaxis
    real(rk) :: dx, dy, dz
    integer :: i
    
    dx = Lx/real(nx,rk)
    dy = Ly/real(ny,rk)
    dz = Lz/real(nz,rk)

    do i=1,nx
      xaxis(i) = real(i-1, rk) * dx
    end do
    do i=1,ny
      yaxis(i) = real(i-1, rk) * dy
    end do
    do i=1,nz
      zaxis(i) = real(i - 0.5, rk) * dz
    end do
  end subroutine

  subroutine compute_abl(reader, Lx, Ly, Lz, runid, baserunid, budget_source, abl_type, path, outdir, start_idx, end_idx, lengthscale)
    implicit none
    class(FieldReader2Decomp), intent(inout) :: reader
    integer,          intent(in)  :: runid, baserunid
    integer,          intent(in)  :: start_idx, end_idx, budget_source, abl_type
    real(rk),         intent(in)  :: Lx, Ly, Lz, lengthscale
    character(*),     intent(in)  :: path, outdir    
    character(len=2) :: rc, rcbase
    integer :: nx, ny, nz
    integer :: nxloc, nyloc, nzloc
    integer :: xs, xe, ys, ye, zs, ze  
    real(rk), allocatable :: uw(:,:,:), vw(:,:,:), buffer(:,:,:), zcross(:,:), zcross_global(:,:)
    real(rk), allocatable :: x(:), y(:), z(:)
    character(len=256) :: f_, fname, Lxc, Lyc, Lzc
    character(len=:), allocatable :: sorted_keys(:), sorted_stamps(:)
    character(1) :: method
    integer :: k, ierr
    
    ! Convert runid to character
    write(rc, '(I2.2)') runid
    write(rcbase, '(I2.2)') baserunid

    ! Shapes and local indices
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%global_shape(nx, ny, nz)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    allocate(uw(nxloc, nyloc, nzloc))
    allocate(vw(nxloc, nyloc, nzloc))
    allocate(buffer(nxloc, nyloc, nzloc))
    allocate(zcross(nx, ny))
    allocate(zcross_global(nx, ny))

    allocate(x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)
    
    if (myrank == 0) then
      write(Lxc, '(F10.3)') Lx
      call message('Domain length is '//trim(Lxc))
      write(Lyc, '(F10.3)') Ly
      call message('Domain width is '//trim(Lyc))
      write(Lzc, '(F10.3)') Lz
      call message('Domain height is '//trim(Lzc))
    end if

    if(budget_source == 1)then
      call get_keys_stamps(trim(path), trim(rc), 1, 'R13', f_, sorted_keys, sorted_stamps)
    else
      call get_keys_stamps(trim(path), trim(rc), 3, 'dup_dwp', f_, sorted_keys, sorted_stamps)
    end if   
    
    ! Loop over time snapshots
    do k = 1, size(sorted_keys)

      if(.not. within_range(start_idx, end_idx, trim(sorted_keys(k)))) cycle

      uw = 0.0_rk
      vw = 0.0_rk
      buffer = 0.0_rk
      zcross = 0.0_rk
      zcross_global = 0.0_rk

      if(abl_type == 0)then
        if(budget_source == 1)then
          ! u'w'
          buffer = eval_field('R13', reader, 1, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          uw = uw + buffer

          buffer = eval_field('tau13', reader, 1, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          uw = uw + buffer

          ! v'w'
          buffer = eval_field('R23', reader, 1, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          vw = vw + buffer

          buffer = eval_field('tau23', reader, 1, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          vw = vw + buffer

        else
          ! u'w'
          buffer = eval_field('R13', reader, 1, trim(path), trim(rcbase), trim(rcbase), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          uw = uw + buffer

          buffer = eval_field('dup_dwp', reader, 3, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          uw = uw + buffer

          buffer = eval_field('dup_bwp', reader, 3, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          uw = uw + buffer

          buffer = eval_field('dwp_bup', reader, 3, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          uw = uw + buffer

          buffer = eval_field('tau13', reader, 1, trim(path), trim(rcbase), trim(rcbase), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          uw = uw + buffer

          buffer = eval_field('delta_tau13', reader, 3, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          uw = uw + buffer

          ! v'w'
          buffer = eval_field('R23', reader, 1, trim(path), trim(rcbase), trim(rcbase), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          vw = vw + buffer

          buffer = eval_field('dvp_dwp', reader, 3, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          vw = vw + buffer

          buffer = eval_field('dvp_bwp', reader, 3, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          vw = vw + buffer

          buffer = eval_field('dwp_bvp', reader, 3, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          vw = vw + buffer

          buffer = eval_field('tau23', reader, 1, trim(path), trim(rcbase), trim(rcbase), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          vw = vw + buffer

          buffer = eval_field('delta_tau23', reader, 3, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))
          vw = vw + buffer
        end if

        ! Shear magnitude
        buffer = sqrt(uw**2 + vw**2)

        ! ABL height
        call find_threshold_crossing_z(buffer, z, zcross(xs:xe, ys:ye), 0.05_rk, 0.0_rk)

        ! MPI exchange
        call MPI_Reduce(zcross, zcross_global, nx*ny, MPI_DOUBLE_PRECISION, &
                  MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if(myrank == 0)then
          write(method, '(I1.1)')abl_type
          fname = trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//'_SL_BLH_M'//method//'.nc' 
          call message('Exporting to: '//trim(fname))
          call export_slice_to_netcdf(trim(fname), 'BLH', zcross_global, x, y, 'x', 'y')
        end if
      
      else if (abl_type == 1)then

        ! Read temperature vertical gradient
        ! if(budget_source == 1)then
        !   buffer = reader%read_field(trim(path)//'/ddz_Run'//trim(rc)//'_budget0_term26_t'//trim(sorted_keys(k))//'_n'//trim(sorted_stamps(k))//'.s3D')
        ! else
        !   ! Base
        !   uw = reader%read_field(trim(path)//'/ddz_Run'//trim(rcbase)//'_budget0_term26_t'//trim(sorted_keys(k))//'_n'//trim(sorted_stamps(k))//'.s3D')
        !   buffer = buffer + uw

        !   uw = reader%read_field(trim(path)//'/ddz_Run'//trim(rc)//'_comp_deficit_budget0_term05_t'//trim(sorted_keys(k))//'_n'//trim(sorted_stamps(k))//'.s3D')
        !   buffer = buffer + uw
        ! end if 

        ! ! Now buffer holds d(theta)/dz
        ! call find_capping_inversion(buffer, z, zcross(xs:xe, ys:ye), 0.15_rk)   

        ! Read temperature
        if(budget_source == 1)then
          buffer = reader%read_field(trim(path)//'/Run'//trim(rc)//'_budget0_term26_t'//trim(sorted_keys(k))//'_n'//trim(sorted_stamps(k))//'.s3D')
        else
          ! Base
          uw = reader%read_field(trim(path)//'/Run'//trim(rcbase)//'_budget0_term26_t'//trim(sorted_keys(k))//'_n'//trim(sorted_stamps(k))//'.s3D')
          buffer = uw

          uw = reader%read_field(trim(path)//'/Run'//trim(rc)//'_comp_deficit_budget0_term05_t'//trim(sorted_keys(k))//'_n'//trim(sorted_stamps(k))//'.s3D')
          buffer = buffer + uw ! Full
        end if 

        block
          real(rk), allocatable :: ytmp(:,:,:), ztmp(:,:,:)
          real(rk), allocatable :: h0map(:,:), h2map(:,:), GMAP(:,:)
          integer :: ic, jc, ig, jg
          type(rz_params) :: params
          real(rk) :: l0, d0, h0, h2
          real(rk) :: xi
          real(rk), allocatable :: tcol(:), z_in_m(:)
          
          allocate(ytmp(reader%gpC%ysz(1), reader%gpC%ysz(2), reader%gpC%ysz(3)))
          allocate(ztmp(reader%gpC%zsz(1), reader%gpC%zsz(2), reader%gpC%zsz(3)))
          allocate(h0map(nx, ny))
          allocate(h2map(nx, ny))
          allocate(GMAP(nx, ny))
          allocate(tcol(nz))
          allocate(z_in_m(nz))

          z_in_m = z * lengthscale

          call transpose_x_to_y(buffer, ytmp, reader%gpC)
          call transpose_y_to_z(ytmp, ztmp, reader%gpC)

          if (reader%gpC%zsz(3) /= nz) then
              if (myrank == 0) call message('ERROR: z-pencil third dimension is not full nz.')
              call MPI_Abort(MPI_COMM_WORLD, 222, ierr)
          end if

          l0 = 700.0_rk
          d0 = 200.0_rk
          xi = 1.3_rk

          h0map = 0.0_rk
          h2map = 0.0_rk

          do ic = 1, reader%gpC%zsz(1)
            do jc = 1, reader%gpC%zsz(2)
              ig = reader%gpC%zst(1) + ic - 1
              jg = reader%gpC%zst(2) + jc - 1
              tcol = ztmp(ic,jc,:)

              call fit_rz_profile( &
                  z_in_m, tcol, nz, l0, d0, params, &
                  d_min   = 1.0e-6_rk, &
                  ridge   = 1.0e-12_rk, &
                  max_iter = 1000, &
                  tol     = 1.0e-10_rk &
              )
              
              if (params%status == 0) then
                  h0 = params%l - xi * params%d
                  h2 = params%l + xi * params%d
              else
                  h0 = 0.0_rk
                  h2 = 0.0_rk
              end if

              h0map(ig,jg) = h0
              h2map(ig,jg) = h2
            end do
          end do

          ! Write h0
          GMAP = 0.0_rk
          ! MPI exchange
          call MPI_Reduce(h0map, GMAP, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          if(myrank == 0)then
            fname = trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//'_INVH0.nc' 
            call message('Exporting to: '//trim(fname))
            call export_slice_to_netcdf(trim(fname), 'INVH0', GMAP, x, y, 'x', 'y')
          end if

          ! Write h2
          GMAP = 0.0_rk
          ! MPI exchange
          call MPI_Reduce(h2map, GMAP, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          if(myrank == 0)then
            fname = trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//'_INVH2.nc' 
            call message('Exporting to: '//trim(fname))
            call export_slice_to_netcdf(trim(fname), 'INVH2', GMAP, x, y, 'x', 'y')
          end if
          deallocate(ytmp, ztmp, h0map, h2map, GMAP, z_in_m, tcol)
        end block    
      end if 
    end do

    if (allocated(uw)) deallocate(uw)
    if (allocated(vw)) deallocate(vw)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(zcross)) deallocate(zcross)
    if (allocated(zcross_global)) deallocate(zcross_global)
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(z)) deallocate(z)
    if (allocated(sorted_keys)) deallocate(sorted_keys)
    if (allocated(sorted_stamps)) deallocate(sorted_stamps)
  end subroutine

  subroutine find_capping_inversion(field, z, zcross, alpha)
    ! This suborutine should be used with a single process
    implicit none
    real(rk), intent(in)  :: field(:,:,:)      ! 3D field: (nx, ny, nz)
    real(rk), intent(in)  :: z(:)              ! z locations, size nz
    real(rk), intent(in), optional :: alpha  
    real(rk), intent(out) :: zcross(size(field,1), size(field,2))
    real(rk) :: alpha_, max_dthetadz
    integer :: nx, ny, nz
    integer :: i, j, k, izmax
    real(rk) :: zstop
    real(rk) :: Gbg, thr
    integer :: count

    nx = size(field,1)
    ny = size(field,2)
    nz = size(field,3)

    ! Basic sanity check
    if (size(z) /= nz) then
      stop "Error in find_capping_inversion: size(z) must equal size(field,3)"
    end if

    ! Defaults
    alpha_  = 0.15
    if (present(alpha)) alpha_ = alpha

    zcross = 0.0_rk
    do j = 1, ny
      do i = 1, nx
        ! find maximum
        max_dthetadz = field(i,j,2)
        izmax = 2
        do k=3, nz
          if(field(i,j,k) > max_dthetadz)then
            max_dthetadz=field(i,j,k)
            izmax = k
          end if
        end do
        
        ! Near-surface mean
        count = 1
        Gbg = field(i,j,2)
        zstop = z(izmax) * 0.3_rk
        do k=3,nz
          if(z(k) > zstop) exit
          count = count + 1
          Gbg = Gbg + field(i,j,k)
        end do
        Gbg = Gbg / real(count, rk)
        thr = Gbg + alpha_ * (max_dthetadz - Gbg)

        ! Find first crossing
        do k=3, izmax
          if(field(i,j,k) > thr)then
            zcross(i,j) = z(k)
            exit
          end if
        end do  
        
      end do
    end do
  end subroutine

  subroutine find_threshold_crossing_z(field, z, zcross, threshold, missing_value)
    implicit none

    !-----------------------------------------------------------------
    ! Inputs
    !-----------------------------------------------------------------
    real(rk), intent(in)  :: field(:,:,:)      ! 3D field: (nx, ny, nz)
    real(rk), intent(in)  :: z(:)              ! z locations, size nz

    ! Optional inputs
    real(rk), intent(in), optional :: threshold      ! default = 0.05
    real(rk), intent(in), optional :: missing_value  ! default = -huge(1.0)

    !-----------------------------------------------------------------
    ! Output
    !-----------------------------------------------------------------
    real(rk), intent(out) :: zcross(size(field,1), size(field,2))

    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------
    integer :: nx, ny, nz
    integer :: i, j, k
    real(rk)    :: thr, miss
    real(rk)    :: fref, fthr
    real(rk)    :: f1, f2
    real(rk)    :: z1, z2
    real(rk)    :: alpha
    logical :: found

    nx = size(field,1)
    ny = size(field,2)
    nz = size(field,3)

    ! Defaults
    thr  = 0.05
    if (present(threshold)) thr = threshold

    miss = 0.0_rk
    if (present(missing_value)) miss = missing_value

    ! Basic sanity check
    if (size(z) /= nz) then
      stop "Error in find_threshold_crossing_z: size(z) must equal size(field,3)"
    end if

    ! Initialize output
    zcross = miss

    do j = 1, ny
      do i = 1, nx
        fref = field(i,j,1)
        fthr = thr * fref
        found = .false.

        ! Handle case where the first point is already below threshold
        if (field(i,j,1) <= fthr) then
            zcross(i,j) = z(1)
            cycle
        end if

        ! Search for first crossing
        do k = 2, nz
          if (field(i,j,k) <= fthr) then
            f1 = field(i,j,k-1)
            f2 = field(i,j,k)
            z1 = z(k-1)
            z2 = z(k)

            ! Linear interpolation:
            ! f(zcross) = fthr
            if (abs(f2 - f1) > tiny(1.0)) then
                alpha = (fthr - f1) / (f2 - f1)
                zcross(i,j) = z1 + alpha * (z2 - z1)
            else
                ! Degenerate case: identical consecutive values
                zcross(i,j) = z1
            end if

            found = .true.
            exit
          end if
        end do

        if (.not. found) then
            zcross(i,j) = miss
        end if
      end do
    end do
  end subroutine find_threshold_crossing_z

  subroutine export_slice_to_netcdf(fname, varname, slice, x1, x2, x1_name, x2_name)
      !! Export a 2D slice on a uniform grid to a NetCDF file.
      !!
      !! Inputs
      !!   fname   : output NetCDF filename
      !!   varname : variable name for the 2D field in the NetCDF file
      !!   slice   : 2D real array, size (nx, ny)
      !!   x1      : coordinate array for dim-1, size nx
      !!   x2      : coordinate array for dim-2, size ny
      !!
      !! Notes
      !! - Writes dimensions (x1,x2) and variable (varname) with units-free metadata.
      !! - Uses netcdf-fortran (module netcdf).
      !!
      implicit none

      character(len=*), intent(in) :: fname
      character(len=*), intent(in) :: varname
      real(rk), intent(in) :: slice(:,:)   ! (nx, ny)
      real(rk), intent(in) :: x1(:)
      real(rk), intent(in) :: x2(:)
      character(len=*), intent(in), optional :: x1_name, x2_name

      integer :: ncid
      integer :: dimid_x1, dimid_x2
      integer :: varid_x1, varid_x2, varid_f
      integer :: nx, ny
      integer :: ierr
      integer :: dimids_f(2)
      character(len=64) :: cx1, cx2

      !-----------------------------------------
      ! Defaults for coordinate names
      !-----------------------------------------
      if (present(x1_name)) then
          cx1 = trim(x1_name)
      else
          cx1 = "x1"
      end if

      if (present(x2_name)) then
          cx2 = trim(x2_name)
      else
          cx2 = "x2"
      end if

      nx = size(slice, 1)
      ny = size(slice, 2)

      if (size(x1) /= nx) error stop "x1 size mismatch"
      if (size(x2) /= ny) error stop "x2 size mismatch"

      !-----------------------------------------
      ! Create file
      !-----------------------------------------
      ierr = nf90_create(trim(fname), ior(NF90_CLOBBER, NF90_NETCDF4), ncid)
      call nc_check(ierr, "nf90_create")

      !-----------------------------------------
      ! Define dimensions
      !-----------------------------------------
      ierr = nf90_def_dim(ncid, trim(cx1), nx, dimid_x1)
      call nc_check(ierr, "def_dim x1")

      ierr = nf90_def_dim(ncid, trim(cx2), ny, dimid_x2)
      call nc_check(ierr, "def_dim x2")

      !-----------------------------------------
      ! Define coordinate variables
      !-----------------------------------------
      ierr = nf90_def_var(ncid, trim(cx1), NF90_DOUBLE, dimid_x1, varid_x1)
      call nc_check(ierr, "def_var x1")

      ierr = nf90_def_var(ncid, trim(cx2), NF90_DOUBLE, dimid_x2, varid_x2)
      call nc_check(ierr, "def_var x2")

      !-----------------------------------------
      ! Define field variable
      !-----------------------------------------
      dimids_f = [dimid_x1, dimid_x2]
      ierr = nf90_def_var(ncid, trim(varname), NF90_DOUBLE, dimids_f, varid_f)
      call nc_check(ierr, "def_var field")

      ierr = nf90_put_att(ncid, varid_f, "long_name", trim(varname)//" slice")
      call nc_check(ierr, "put_att field")

      ierr = nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8")
      call nc_check(ierr, "put_att global")

      ierr = nf90_enddef(ncid)
      call nc_check(ierr, "enddef")

      !-----------------------------------------
      ! Write data
      !-----------------------------------------
      ierr = nf90_put_var(ncid, varid_x1, x1)
      call nc_check(ierr, "put_var x1")

      ierr = nf90_put_var(ncid, varid_x2, x2)
      call nc_check(ierr, "put_var x2")

      ierr = nf90_put_var(ncid, varid_f, slice)
      call nc_check(ierr, "put_var slice")

      ierr = nf90_close(ncid)
      call nc_check(ierr, "close")

    contains

      subroutine nc_check(status, where)
          integer, intent(in) :: status
          character(len=*), intent(in) :: where
          if (status /= nf90_noerr) then
            write(*,*) "NetCDF error in ", trim(where), ": ", trim(nf90_strerror(status))
            error stop
          end if
      end subroutine nc_check

    end subroutine export_slice_to_netcdf

  subroutine miscellaneous_driver(reader, runid, field, path, start_idx, end_idx)
    implicit none
    class(FieldReader2Decomp), intent(inout) :: reader
    integer, intent(in) :: runid, start_idx, end_idx
    character(*), intent(in) :: field, path

    call verify_budgets(reader, runid, trim(field), trim(path), start_idx, end_idx)

  end subroutine miscellaneous_driver

  subroutine verify_budgets(reader, runid, field, path, start_idx, end_idx)
    implicit none
    class(FieldReader2Decomp), intent(inout) :: reader
    integer, intent(in) :: runid, start_idx, end_idx
    character(*), intent(in) :: field, path
    character(len=256) :: field_, f_
    character(len=2) :: rc, rcbase
    integer :: k, ierr
    real(rk) :: error, global_error
    real(rk), allocatable :: errors(:)

    ! Grid and indices
    integer :: nx, ny, nz
    integer :: nxloc, nyloc, nzloc
    integer :: xs, xe, ys, ye, zs, ze

    ! Time keys
    character(len=:), allocatable :: sorted_keys(:), sorted_stamps(:)

    ! Local 3D field
    real(rk), allocatable, target :: f1(:,:,:), fprim(:,:,:), fpre(:,:,:), ferr(:,:,:)

    ! Convert runid to character
    write(rc, '(I2.2)') runid
    write(rcbase, '(I2.2)') (runid-1)

    ! Shapes and local indices
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%global_shape(nx, ny, nz)
    call reader%indices(xs, xe, ys, ye, zs, ze)

    ! Allocate local 3D field
    allocate(f1(nxloc, nyloc, nzloc))
    allocate(fprim(nxloc, nyloc, nzloc))
    allocate(fpre(nxloc, nyloc, nzloc))
    allocate(ferr(nxloc, nyloc, nzloc))

    call get_keys_stamps(trim(path), trim(rc), 3, trim(field), f_, sorted_keys, sorted_stamps)
    call deficit_field_to_time_field(trim(field), field_)
    if(myrank == 0) call message('Verifying '// trim(field)//' vs difference in '//trim(field_))
    
    allocate(errors(size(sorted_keys)))
    if(myrank == 0)call message(' ')

    do k = 1, size(sorted_keys)
      if(.not. within_range(start_idx, end_idx, trim(sorted_keys(k)))) cycle

      ! Deficit field
      f1 = eval_field(trim(field), reader, mbdgtsrc, trim(path), trim(rc), trim(rcbase), trim(sorted_keys(k)), trim(sorted_stamps(k)))

      ! Primary field
      fprim = eval_field(trim(field_), reader, timeavgsrc, trim(path), trim(rc), trim(rc), trim(sorted_keys(k)), trim(sorted_stamps(k)))

      ! Precursor field
      fpre = eval_field(trim(field_), reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(sorted_keys(k)), trim(sorted_stamps(k)))

      ! Error
      ferr = f1 - (fprim - fpre)
      ferr = abs(ferr)
      error = maxval(ferr)
      if(myrank == 0) print*, error

      ! Share with other ranks
      call MPI_Allreduce(error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      errors(k) = global_error
      
      if(myrank == 0) call message(' ')
    end do

    ! Export the max errors for each time step
    if (myrank == 0)then
      do k = 1, size(sorted_keys)
        call message('Time stamp: '//trim(sorted_keys(k))//', max abs error: '//trim(to_string(errors(k))))
      end do
      call message(' ')
    end if
  end subroutine verify_budgets

  function reverse_base_budget_sign(field) result(sgn)
    implicit none
    character(*), intent(in) :: field
    integer :: sgn
    select case(trim(field))
    case('adv_x', 'adv_y', 'adv_z','ddx_p','ddy_p','ddz_p', &
         'usgs','vsgs','wsgs','pcov','sgsTr', &
         'tkeTr')
      sgn = -1
    case default
      sgn = 1
    end select
  end function reverse_base_budget_sign

  subroutine deficit_field_to_time_field(field, field_)
    character(*), intent(in) :: field
    character(len=256), intent(out) :: field_

    select case (trim(field))
    case ('delta_u')
      field_ = 'ubar'
    case ('delta_v')
      field_ = 'vbar'
    case ('delta_w')
      field_ = 'wbar'
    case ('delta_T')
      field_ = 'tbar'
    case ('delta_p')
      field_ = 'pbar'
    case ('delta_tau11')
      field_ = 'tau11'
    case('delta_tau12')
      field_ = 'tau12'
    case('delta_tau13')
      field_ = 'tau13'
    case('delta_tau22')
      field_ = 'tau22'
    case('delta_tau23')
      field_ = 'tau23'
    case('delta_tau33')
      field_ = 'tau33'
    case('delta_ucor')
      field_ = 'ucor'
    case('delta_vcor')
      field_ = 'vcor'
    case('delta_bouyancy')
      field_ = 'bouyancy'
    case('delta_turbx')
      field_ = 'turbx'
    case('delta_turby')
      field_ = 'turby'
    case('delta_R11')
      field_ = 'R11'
    case('delta_R12')
      field_ = 'R12'
    case('delta_R13')
      field_ = 'R13'
    case('delta_R22')
      field_ = 'R22'
    case('delta_R23')
      field_ = 'R23'
    case('delta_R33')
      field_ = 'R33'
    case('delta_adv_x')
      field_ = 'adv_x'
    case('delta_adv_y')
      field_ = 'adv_y'
    case('delta_adv_z')
      field_ = 'adv_z'
    case('ddx_delta_p')
      field_ = 'ddx_p'
    case('ddy_delta_p')
      field_ = 'ddy_p'
    case('ddz_delta_p')
      field_ = 'ddz_p'
    case('delta_usgs')
      field_ = 'usgs'
    case('delta_vsgs')
      field_ = 'vsgs'
    case('delta_wsgs')
      field_ = 'wsgs'
    case('delta_pcov')
      field_ = 'pcov'
    case('delta_sgsTr')
      field_ = 'sgsTr'
    case('delta_sgsDiss')
      field_ = 'sgsDiss'
    case('delta_tcov')
      field_ = 'tcov'
    case('delta_tkeTr')
      field_ = 'tkeTr'
    end select

  end subroutine deficit_field_to_time_field

  subroutine get_keys_stamps(path, rc, budget_source, field, f_, sorted_keys, sorted_stamps)
    implicit none
    character(*), intent(in) :: path, rc
    integer, intent(in) :: budget_source
    character(*), intent(in) :: field
    character(len=:), allocatable, intent(inout) :: sorted_keys(:), sorted_stamps(:)
    character(len=256), intent(out):: f_
    character(len=:), allocatable :: keys(:), stamps(:)
    character(len=2048) :: pattern
    character(len=1) :: budget
    character(len=2) :: term    
    logical :: calc=.false.

    if (budget_source == instfieldsrc)then
      f_ = field_to_name(trim(field))
      pattern = 'Run'//trim(rc)//'_'//trim(f_)//'_t*.out'

      if (myrank == 0)call message('Pattern is '//trim(pattern))

      call list_matching_keys(trim(path), trim(pattern), keys)
      sorted_keys = sort_keys_numeric(keys)
    else
      f_ = trim(field)
      pattern = 'Run'//trim(rc)
      call define_budget(trim(field), budget, term, budget_source, calc) 
      if(budget_source == 2)then        
        pattern = trim(pattern)//'_deficit'
      elseif(budget_source == 3) then       
        pattern = trim(pattern)//'_comp_deficit'
      elseif(budget_source == 4) then       
        pattern = trim(pattern)//'_mdeficit'
      end if
      pattern = trim(pattern)//'_budget'//budget//'_term'//term//'_t*_n~.s3D'
      if (myrank == 0)then
        call message('Pattern is '//trim(pattern))
        if(calc) call message('The budget is evaluated using multiple files')
      end if
      call list_matching_keys_budget(trim(path), trim(pattern), keys, stamps)
      call sort_keys_and_stamps_numeric(keys, stamps, sorted_keys, sorted_stamps)
    end if
    if (myrank == 0) call message(' ')
  end subroutine

  function eval_field(field, reader, budget_source, path, rc, rcbase, key, stamp) result(f1)
    implicit none
    real(rk) :: f1(nxloc_,nyloc_,nzloc_)
    class(FieldReader2Decomp), intent(inout) :: reader
    character(*), intent(in) :: path, rc, rcbase, field, key, stamp
    integer, intent(in) :: budget_source
    character(len=256) :: filename, f_
    character(len=2) :: term
    character(len=1) :: budget
    logical :: calc=.false.
    integer :: sgn=1

    call define_budget(trim(field), budget, term, budget_source, calc)
    f_ = field_to_name(trim(field))

    if(calc)then
        ! Compute budget from multiple files
        call compute_budget_from_multiple_files(reader, trim(path), &
          trim(rc), trim(rcbase), trim(field), trim(key), trim(stamp), f1)
    else
        sgn = reverse_base_budget_sign(trim(field))
        if(myrank == 0 .and. sgn == -1) call message('Reversed sign of time-averaged budget '//trim(field))
    
        call create_file_name(budget_source, trim(rc), trim(f_), trim(budget), &
                    trim(term), trim(key), trim(stamp), filename)
        f1 = reader%read_field(trim(path)//'/'//trim(filename))
        f1 = f1 * sgn
    end if
  end function

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
    real(rk), intent(in) :: a, b
    integer,      intent(in) :: n
    real(rk), allocatable :: x(:)
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
    character(*), intent(in) :: field
    character(len=256)            :: name
    select case (trim(field))
    case ('u'); name = 'uVel'
    case ('v'); name = 'vVel'
    case ('w'); name = 'wVel'
    case ('T'); name = 'potT'
    case ('p'); name = 'prss'
    case ('S'); name = 'uVel'
    case default; name = trim(field)
    end select
  end function field_to_name

  ! Utility function to form a file name
  subroutine create_file_name(budget_source, rc, field, budget, term, key, stamp, fname)
    implicit none
    integer, intent(in) :: budget_source
    character(len=*), intent(in) :: rc, field, budget, term, key, stamp
    character(len=*), intent(out) :: fname

    if (budget_source == 0) then
        fname = 'Run'//trim(rc)//'_'//trim(field)//'_t'//trim(key)//'.out'
    else
      ! Budgets
      fname = 'Run' // trim(rc)

      select case (budget_source)
      case (1)
          fname = trim(fname) // '_budget'
      case (2)
          fname = trim(fname) // '_deficit_budget'
      case (3)
          fname = trim(fname) // '_comp_deficit_budget'
      case (4)
        fname = trim(fname) // '_mdeficit_budget'
      end select
      
      fname = trim(fname)//trim(budget)//'_term'//trim(term)// &
            '_t'//trim(key)//'_n'//trim(stamp)//'.s3D'
    end if   
  end subroutine create_file_name

  ! Utility function to compute a field from multiple files
  subroutine compute_budget_from_multiple_files(reader, path, rc, rcbase, field, key, stamp, f1)
  implicit none
  class(FieldReader2Decomp), intent(inout) :: reader
  character(*), intent(in) :: path, rc, rcbase, field, key, stamp
  real(rk), intent(out) :: f1(:,:,:)
  
  select case (trim(field))
  case('delta_R11')
    f1 = eval_field('dup_dup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
     2.0*eval_field('dup_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
    
  case('delta_R12')
    f1 = eval_field('dup_dvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('dup_bvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('dvp_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('delta_R13')
    f1 = eval_field('dup_dwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('dup_bwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('dwp_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('delta_R22')
    f1 = eval_field('dvp_dvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
     2.0*eval_field('dvp_bvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('delta_R23')
    f1 = eval_field('dvp_dwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('dvp_bwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('dwp_bvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('delta_R33')
    f1 = eval_field('dwp_dwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
     2.0*eval_field('dwp_bwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('dj_d1p_bjp')
    f1 = eval_field('ddx_dup_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddy_dup_bvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddz_dup_bwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('dj_d1p_djp')
    f1 = eval_field('ddx_dup_dup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddy_dup_dvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddz_dup_dwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('dj_b1p_djp')
    f1 = eval_field('ddx_dup_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddy_dvp_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddz_dwp_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
  
  case('delta_adv_x')
    block 
      real(rk), allocatable :: bf1(:,:,:), bf2(:,:,:), bf3(:,:,:)
      allocate(bf1(size(f1,1), size(f1,2), size(f1,3)))
      allocate(bf2(size(f1,1), size(f1,2), size(f1,3)))
      allocate(bf3(size(f1,1), size(f1,2), size(f1,3)))
            
      bf1 = eval_field('ddx_delta_u', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
      bf2 = eval_field('ddy_delta_u', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
      bf3 = eval_field('ddz_delta_u', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

      ! f1 = eval_field('delta_u',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) * bf1 + &
      !      eval_field('delta_v',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) * bf2 + &
      !      eval_field('delta_w',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) * bf3 

      f1 = eval_field('delta_u',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf1+eval_field('ddx_base_u', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
           eval_field('delta_v',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf2+eval_field('ddy_base_u', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
           eval_field('delta_w',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf3+eval_field('ddz_base_u', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
        bf1*eval_field('ubar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp))  + &
        bf2*eval_field('vbar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp))  + &
        bf3*eval_field('wbar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp)) &
         + eval_field('ddx_dup_dup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) &
         + eval_field('ddy_dup_dvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) &
         + eval_field('ddz_dup_dwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) &
         + eval_field('ddx_dup_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) &
         + eval_field('ddy_dup_bvp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) &
         + eval_field('ddz_dup_bwp', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) &
         + eval_field('ddy_dvp_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) &
         + eval_field('ddz_dwp_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) &
         + eval_field('ddx_dup_bup', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
      deallocate(bf1, bf2, bf3)
  end block
  case('delta_adv_y')
    block 
      real(rk), allocatable :: bf1(:,:,:), bf2(:,:,:), bf3(:,:,:)
      allocate(bf1(size(f1,1), size(f1,2), size(f1,3)))
      allocate(bf2(size(f1,1), size(f1,2), size(f1,3)))
      allocate(bf3(size(f1,1), size(f1,2), size(f1,3)))
            
      bf1 = eval_field('ddx_delta_v', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
      bf2 = eval_field('ddy_delta_v', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
      bf3 = eval_field('ddz_delta_v', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

      f1 = eval_field('delta_u',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf1+eval_field('ddx_base_v', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
           eval_field('delta_v',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf2+eval_field('ddy_base_v', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
           eval_field('delta_w',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf3+eval_field('ddz_base_v', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
        bf1*eval_field('ubar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp))  + &
        bf2*eval_field('vbar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp))  + &
        bf3*eval_field('wbar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp)) 
    deallocate(bf1, bf2, bf3)
  end block
  
  case('delta_adv_z')
    block 
      real(rk), allocatable :: bf1(:,:,:), bf2(:,:,:), bf3(:,:,:)
      allocate(bf1(size(f1,1), size(f1,2), size(f1,3)))
      allocate(bf2(size(f1,1), size(f1,2), size(f1,3)))
      allocate(bf3(size(f1,1), size(f1,2), size(f1,3)))
            
      bf1 = eval_field('ddx_delta_w', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
      bf2 = eval_field('ddy_delta_w', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
      bf3 = eval_field('ddz_delta_w', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

      f1 = eval_field('delta_u',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf1+eval_field('ddx_base_w', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
           eval_field('delta_v',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf2+eval_field('ddy_base_w', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
           eval_field('delta_w',    reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))   * &
      (bf3+eval_field('ddz_base_w', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)))  + &
        bf1*eval_field('ubar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp))  + &
        bf2*eval_field('vbar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp))  + &
        bf3*eval_field('wbar', reader, timeavgsrc, trim(path), trim(rcbase), trim(rcbase), trim(key), trim(stamp)) 
    deallocate(bf1, bf2, bf3)
  end block

  case('delta_div')
    f1 = eval_field('ddx_delta_u', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddy_delta_v', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddz_delta_w', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('base_div')  
    f1 = eval_field('ddx_base_u', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddy_base_v', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('ddz_base_w', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
  
  ! case('ddj_delta_tau1j')
  !   f1 = eval_field('ddx_delta_tau11', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
  !        eval_field('ddy_delta_tau12', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
  !        eval_field('ddz_delta_tau13', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('delta_pcov')
    f1 = eval_field('pcov_dd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('pcov_bd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('pcov_db', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))
  
  case('delta_sgsTr')
    f1 = eval_field('sgsTr_bd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('sgsTr_db', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('sgsTr_dd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('delta_sgsDiss')
    f1 = eval_field('sgsDiss_db', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('sgsDiss_bd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('sgsDiss_dd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('delta_tcov')
    f1 = eval_field('tcov_dd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('tcov_db', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('tcov_bd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

  case('delta_tkeTr')
    f1 = eval_field('tkeTr_dbb', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
    2.d0*eval_field('tkeTr_bbd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
    2.d0*eval_field('tkeTr_dbd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('tkeTr_bdd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp)) +&
         eval_field('tkeTr_ddd', reader, mbdgtsrc, trim(path), trim(rc), trim(rc), trim(key), trim(stamp))

    f1 = f1 / 2.d0

  end select
  end subroutine compute_budget_from_multiple_files

  ! Utility function to get proper budget name
  subroutine define_budget(field, b, t, budget_source, calc)
    character(*), intent(in) :: field
    integer, intent(in) :: budget_source
    logical, intent(inout) :: calc
    character(1), intent(out) :: b
    character(2), intent(out) :: t

    calc=.false.
    if(budget_source == 1)then
      ! Time averaged budgets
      select case(trim(field))
      case('ubar')
        b = '0'; t = '01'
      case('vbar')
        b = '0'; t = '02'
      case('wbar')
        b = '0'; t = '03'      
      case('R11')
        b = '0'; t = '04'
      case('R12')
        b = '0'; t = '05'
      case('R13')
        b = '0'; t = '06'
      case('R22')
        b = '0'; t = '07'
      case('R23')
        b = '0'; t = '08'
      case('R33')
        b = '0'; t = '09'
      case('pbar')
        b = '0'; t = '10'
      case('tau11')
        b = '0'; t = '11'
      case('tau12')
        b = '0'; t = '12'
      case('tau13')
        b = '0'; t = '13'
      case('tau22')
        b = '0'; t = '14'
      case('tau23')
        b = '0'; t = '15'
      case('tau33')
        b = '0'; t = '16'
      case('up_p')
        b = '0'; t = '17'
      case('vp_p')
        b = '0'; t = '18'
      case('wp_p')
        b = '0'; t = '19'
      case('tbar')
        b = '0'; t = '26'
      case('up_tp')
        b = '0'; t = '27'
      case('vp_tp')
        b = '0'; t = '28'
      case('wp_tp')
        b = '0'; t = '29'
      case('tp_tp')
        b = '0'; t = '30'
      case('bouyancy')
        b = '0'; t = '31'
      case('adv_x')
        b = '1'; t = '01'
      case('ddx_p')
        b = '1'; t = '02'
      case('ddy_p')
        b = '1'; t = '06'
      case('ddz_p')
        b = '1'; t = '09'
      case('usgs')
        b = '1'; t = '03'
      case('turbx')
        b = '1'; t = '04'
      case('adv_y')
        b = '1'; t = '05'
      case('vsgs')
        b = '1'; t = '07'
      case('adv_z')
        b = '1'; t = '08'
      case('wsgs')
        b = '1'; t = '10'
      case('ucor')
        b = '1'; t = '11'
      case('vcor')
        b = '1'; t = '13'
      case('turby')
        b = '1'; t = '15'  
      case('tkeTr')
        b = '3'; t = '03'
      case('pcov')  
        b = '3'; t = '04'
      case('sgsTr')
        b = '3'; t = '05'
      case('sgsDiss')
        b = '3'; t = '06'
      case('tcov')  
        b = '3'; t = '08'
      end select

    elseif(budget_source == 2)then
      ! Deficit budgets
      select case(trim(field))
      case('delta_u')
        b = '0'; t = '01'
      case('delta_v')
        b = '0'; t = '02'
      case('delta_w')
        b = '0'; t = '03'
      case('delta_p')
        b = '0'; t = '04'
      case('dup_dup')
        b = '0'; t = '05'
      case('dup_dvp')
        b = '0'; t = '06'
      case('dup_dwp')
        b = '0'; t = '07' 
      case('dvp_dvp')
        b = '0'; t = '08'
      case('dvp_dwp')
        b = '0'; t = '09'
      case('dwp_dwp')
        b = '0'; t = '10'
      case('dup_bup')
        b = '0'; t = '11'
      case('dup_bvp')
        b = '0'; t = '12'
      case('dvp_bup')
        b = '0'; t = '13' 
      case('dup_bwp')
        b = '0'; t = '14'
      case('dwp_bup')
        b = '0'; t = '15'
      case('dvp_bvp')
        b = '0'; t = '16'
      case('dvp_bwp')
        b = '0'; t = '17'
      case('dwp_bvp')
        b = '0'; t = '18'
      case('dwp_bwp')
        b = '0'; t = '19'
      case('delta_tau11')
        b = '0'; t = '20'
      case('delta_tau12')
        b = '0'; t = '21'
      case('delta_tau13')
        b = '0'; t = '22'
      case('delta_tau22')
        b = '0'; t = '23'
      case('delta_tau23')
        b = '0'; t = '24'
      case('delta_tau33')
        b = '0'; t = '25'
      case('delta_T')
        b = '0'; t = '26'
      case('dup_dTp')
        b = '0'; t = '27'
      case('dvp_dTp')
        b = '0'; t = '28'
      case('dwp_dTp')
        b = '0'; t = '29'
      case('b_dj_dup')
        b = '1'; t = '05'
      case('d_dj_dup')
        b = '1'; t = '06'
      case('d_dj_bup')
        b = '1'; t = '07'
      case('b_dj_dvp')
        b = '1'; t = '15'
      case('d_dj_dvp')
        b = '1'; t = '16'
      case('d_dj_bvp')
        b = '1'; t = '17'
      case('b_dj_dwp')
        b = '1'; t = '25'
      case('d_dj_dwp')
        b = '1'; t = '26'
      case('d_dj_bwp')
        b = '1'; t = '27'
      case('B2_14')
        b = '2'; t = '14'
      case('B2_18')
        b = '2'; t = '18'
      case('sgsTr_dd')    
        b = '3'; t = '05'
      case('sgsDiss_dd')    
        b = '3'; t = '06'
      case('tkeTr_ddd')
        b = '3'; t = '12'
      case('tkeTr_bdd')
        b = '3'; t = '13'
      case('tkeTr_ddb')
        b = '3'; t = '14'
      case('B3_16')
        b = '3'; t = '16'
      case('tkeTr_dbd')
        b = '3'; t = '22'
      end select

    elseif(budget_source == 3)then
      ! Compact deficit budgets
      select case(trim(field))
      case('delta_u')
        b = '0'; t = '01'
      case('delta_v')
        b = '0'; t = '02'
      case('delta_w')
        b = '0'; t = '03'
      case('delta_p')
        b = '0'; t = '04'
      case('delta_T')
        b = '0'; t = '05'
      case('delta_tau11')
        b = '0'; t = '06'
      case('delta_tau12')
        b = '0'; t = '07'
      case('delta_tau13')
        b = '0'; t = '08'
      case('delta_tau22')
        b = '0'; t = '09'
      case('delta_tau23')
        b = '0'; t = '10'
      case('delta_tau33')
        b = '0'; t = '11'
      case('delta_usgs')
        b = '0'; t = '12'
      case('delta_vsgs')
        b = '0'; t = '13'
      case('delta_wsgs')
        b = '0'; t = '14'
      case('delta_ucor')
        b = '0'; t = '15'
      case('delta_vcor')
        b = '0'; t = '16'
      case('delta_bouyancy')
        b = '0'; t = '17'
      case('ddx_delta_p')
        b = '0'; t = '18'
      case('ddy_delta_p')
        b = '0'; t = '19'
      case('ddz_delta_p')
        b = '0'; t = '20'
      case('delta_turbx')
        b = '0'; t = '21'
      case('delta_turby')
        b = '0'; t = '22'

      case('dup_dup')
        b = '1'; t = '01'
      case('dup_dvp')
        b = '1'; t = '02'
      case('dup_dwp')
        b = '1'; t = '03' 
      case('dvp_dvp')
        b = '1'; t = '04'
      case('dvp_dwp')
        b = '1'; t = '05'
      case('dwp_dwp')
        b = '1'; t = '06'
      case('dup_bup')
        b = '1'; t = '07'
      case('dup_bvp')
        b = '1'; t = '08'
      case('dvp_bup')
        b = '1'; t = '09' 
      case('dup_bwp')
        b = '1'; t = '10'
      case('dwp_bup')
        b = '1'; t = '11'
      case('dvp_bvp')
        b = '1'; t = '12'
      case('dvp_bwp')
        b = '1'; t = '13'
      case('dwp_bvp')
        b = '1'; t = '14'
      case('dwp_bwp')
        b = '1'; t = '15'

      case('d_dj_dup')
        b = '2'; t = '01'
      case('d_dj_dvp')
        b = '2'; t = '02'
      case('d_dj_dwp')
        b = '2'; t = '03' 
      case('d_dj_bup')
        b = '2'; t = '04'
      case('d_dj_bvp')
        b = '2'; t = '05'
      case('d_dj_bwp')
        b = '2'; t = '06'
      case('b_dj_dup')
        b = '2'; t = '07'
      case('b_dj_dvp')
        b = '2'; t = '08'
      case('b_dj_dwp')
        b = '2'; t = '09'
      case('b_dj_bup')
        b = '2'; t = '10'
      case('b_dj_bvp')
        b = '2'; t = '11'
      case('b_dj_bwp')
        b = '2'; t = '12'
      case('bup_dup_ddx')
        b = '2'; t = '13'
      case('bup_dup_ddy')
        b = '2'; t = '14'
      case('bup_dup_ddz')
        b = '2'; t = '15' 

      case('pcov_dd')   
        b = '3'; t = '01'
      case('pcov_bd')    
        b = '3'; t = '02'
      case('pcov_db')    
        b = '3'; t = '03'
      case('sgsTr_bd')    
        b = '3'; t = '04'
      case('sgsTr_db')    
        b = '3'; t = '05'
      case('sgsTr_dd')    
        b = '3'; t = '06'
      case('sgsDiss_db')    
        b = '3'; t = '07'
      case('sgsDiss_bd')    
        b = '3'; t = '08'
      case('sgsDiss_dd')    
        b = '3'; t = '09'
      case('tcov_dd')
        b = '3'; t = '10'
      case('tcov_db')
        b = '3'; t = '11'
      case('tcov_bd')
        b = '3'; t = '12'
      case('tkeTr_dbb')
        b = '3'; t = '13'
      case('tkeTr_bbd')
        b = '3'; t = '14'
      case('tkeTr_bdb')
        b = '3'; t = '15'
      case('tkeTr_dbd')
        b = '3'; t = '16'
      case('tkeTr_ddb')
        b = '3'; t = '17'
      case('tkeTr_bdd')
        b = '3'; t = '18'
      case('tkeTr_ddd')
        b = '3'; t = '19'
      case('turbcov_dd')
        b = '3'; t = '20'
      case('turbcov_db')
        b = '3'; t = '21'

      case('Adv_ddd')
        b = '4'; t = '01'
      case('Adv_ddb')
        b = '4'; t = '02'
      case('Adv_dbb')
        b = '4'; t = '03'
      case('Adv_bdd')
        b = '4'; t = '04'
      case('Adv_bdb')
        b = '4'; t = '05'
      case('prod_ddd')
        b = '4'; t = '06'
      case('prod_dbd')
        b = '4'; t = '07'
      case('prod_bdd')
        b = '4'; t = '08'
      case('prod_bbd')
        b = '4'; t = '09'
      case('prod_ddb')
        b = '4'; t = '10'
      case('prod_dbb')
        b = '4'; t = '11'
      case('prod_bdb')
        b = '4'; t = '12'

      case('adv_dudu')
        b = '5'; t = '01'
      case('adv_dvdu')
        b = '5'; t = '02'
      case('adv_dwdu')
        b = '5'; t = '03'
      case('adv_dubu')
        b = '5'; t = '04'
      case('adv_dvbu')
        b = '5'; t = '05'
      case('adv_dwbu')
        b = '5'; t = '06'
      case('adv_budu')
        b = '5'; t = '07'
      case('adv_bvdu')
        b = '5'; t = '08'
      case('adv_bwdu')
        b = '5'; t = '09'
      case('dx_dup_dup')
        b = '5'; t = '16'
      case('dy_dup_dvp')
        b = '5'; t = '17'
      case('dz_dup_dwp')
        b = '5'; t = '18'
      case('dx_dup_bup')
        b = '5'; t = '19'
      case('dy_dup_bvp')
        b = '5'; t = '20'
      case('dz_dup_bwp')
        b = '5'; t = '21'
      case('dx_bup_dup')
        b = '5'; t = '22'
      case('dy_bup_dvp')
        b = '5'; t = '23'
      case('dz_bup_dwp')
        b = '5'; t = '24' 

      case('adv_dudv')
        b = '6'; t = '01'
      case('adv_dvdv')
        b = '6'; t = '02'
      case('adv_dwdv')
        b = '6'; t = '03'
      case('adv_dubv')
        b = '6'; t = '04'
      case('adv_dvbv')
        b = '6'; t = '05'
      case('adv_dwbv')
        b = '6'; t = '06'
      case('adv_budv')
        b = '6'; t = '07'
      case('adv_bvdv')
        b = '6'; t = '08'
      case('adv_bwdv')
        b = '6'; t = '09'
      case('dx_dvp_dup')
        b = '6'; t = '16'
      case('dy_dvp_dvp')
        b = '6'; t = '17'
      case('dz_dvp_dwp')
        b = '6'; t = '18'
      case('dx_dvp_bup')
        b = '6'; t = '19'
      case('dy_dvp_bvp')
        b = '6'; t = '20'
      case('dz_dvp_bwp')
        b = '6'; t = '21'
      case('dx_bvp_dup')
        b = '6'; t = '22'
      case('dy_bvp_dvp')
        b = '6'; t = '23'
      case('dz_bvp_dwp')
        b = '6'; t = '24'  

      case('adv_dudw')
        b = '7'; t = '01'
      case('adv_dvdw')
        b = '7'; t = '02'
      case('adv_dwdw')
        b = '7'; t = '03'
      case('adv_dubw')
        b = '7'; t = '04'
      case('adv_dvbw')
        b = '7'; t = '05'
      case('adv_dwbw')
        b = '7'; t = '06'
      case('adv_budw')
        b = '7'; t = '07'
      case('adv_bvdw')
        b = '7'; t = '08'
      case('adv_bwdw')
        b = '7'; t = '09'
      case('dx_dwp_dup')
        b = '7'; t = '16'
      case('dy_dwp_dvp')
        b = '7'; t = '17'
      case('dz_dwp_dwp')
        b = '7'; t = '18'
      case('dx_dwp_bup')
        b = '7'; t = '19'
      case('dy_dwp_bvp')
        b = '7'; t = '20'
      case('dz_dwp_bwp')
        b = '7'; t = '21'
      case('dx_bwp_dup')
        b = '7'; t = '22'
      case('dy_bwp_dvp')
        b = '7'; t = '23'
      case('dz_bwp_dwp')
        b = '7'; t = '24'   
  
      case default
        b = '0'; t = '01'; calc = .true.     
      end select

    elseif(budget_source == 4)then
      ! Mini deficit budgets
      select case(trim(field))
      case('delta_u')
        b = '0'; t = '01'
      case('delta_v')
        b = '0'; t = '02'
      case('delta_w')
        b = '0'; t = '03'
      case('delta_p')
        b = '0'; t = '04'
      case('delta_T')
        b = '0'; t = '05'
      case('delta_tau11')
        b = '0'; t = '06'
      case('delta_tau12')
        b = '0'; t = '07'
      case('delta_tau13')
        b = '0'; t = '08'
      case('delta_tau22')
        b = '0'; t = '09'
      case('delta_tau23')
        b = '0'; t = '10'
      case('delta_tau33')
        b = '0'; t = '11'
      case('delta_usgs')
        b = '0'; t = '12'
      case('delta_vsgs')
        b = '0'; t = '13'
      case('delta_wsgs')
        b = '0'; t = '14'
      case('delta_ucor')
        b = '0'; t = '15'
      case('delta_vcor')
        b = '0'; t = '16'
      case('delta_bouyancy')
        b = '0'; t = '17'
      case('ddx_delta_p')
        b = '0'; t = '18'
      case('ddy_delta_p')
        b = '0'; t = '19'
      case('ddz_delta_p')
        b = '0'; t = '20'
      case('delta_turbx')
        b = '0'; t = '21'
      case('delta_turby')
        b = '0'; t = '22'

      case('dup_dup')
        b = '1'; t = '01'
      case('dup_dvp')
        b = '1'; t = '02'
      case('dup_dwp')
        b = '1'; t = '03' 
      case('dvp_dvp')
        b = '1'; t = '04'
      case('dvp_dwp')
        b = '1'; t = '05'
      case('dwp_dwp')
        b = '1'; t = '06'
      case('dup_bup')
        b = '1'; t = '07'
      case('dup_bvp')
        b = '1'; t = '08'
      case('dvp_bup')
        b = '1'; t = '09' 
      case('dup_bwp')
        b = '1'; t = '10'
      case('dwp_bup')
        b = '1'; t = '11'
      case('dvp_bvp')
        b = '1'; t = '12'
      case('dvp_bwp')
        b = '1'; t = '13'
      case('dwp_bvp')
        b = '1'; t = '14'
      case('dwp_bwp')
        b = '1'; t = '15'

      case('ddx_delta_u')
        b = '2'; t = '01'
      case('ddy_delta_u')
        b = '2'; t = '02' 
      case('ddz_delta_u')
        b = '2'; t = '03'
      case('ddx_delta_v')
        b = '2'; t = '04'
      case('ddy_delta_v')
        b = '2'; t = '05' 
      case('ddz_delta_v')
        b = '2'; t = '06' 
      case('ddx_delta_w')
        b = '2'; t = '07'
      case('ddy_delta_w')
        b = '2'; t = '08' 
      case('ddz_delta_w')
        b = '2'; t = '09'
      case('ddx_base_u')
        b = '2'; t = '10'
      case('ddy_base_u')
        b = '2'; t = '11' 
      case('ddz_base_u')
        b = '2'; t = '12'
      case('ddx_base_v')
        b = '2'; t = '13'
      case('ddy_base_v')
        b = '2'; t = '14' 
      case('ddz_base_v')
        b = '2'; t = '15' 
      case('ddx_base_w')
        b = '2'; t = '16'
      case('ddy_base_w')
        b = '2'; t = '17' 
      case('ddz_base_w')
        b = '2'; t = '18'

      case('ddx_dup_dup')
        b = '3'; t = '01'
      case('ddy_dup_dvp')
        b = '3'; t = '02'
      case('ddz_dup_dwp')
        b = '3'; t = '03'      
      case('ddx_dup_bup')
        b = '3'; t = '04'
      case('ddy_dup_bvp')
        b = '3'; t = '05'
      case('ddz_dup_bwp')
        b = '3'; t = '06'
      case('ddy_dvp_bup')
        b = '3'; t = '07'
      case('ddz_dwp_bup')
        b = '3'; t = '08' 

      case('pcov_dd')   
        b = '4'; t = '01'
      case('pcov_bd')    
        b = '4'; t = '02'
      case('pcov_db')    
        b = '4'; t = '03'
      case('sgsTr_bd')    
        b = '4'; t = '04'
      case('sgsTr_db')    
        b = '4'; t = '05'
      case('sgsTr_dd')    
        b = '4'; t = '06'
      case('sgsDiss_db')    
        b = '4'; t = '07'
      case('sgsDiss_bd')    
        b = '4'; t = '08'
      case('sgsDiss_dd')    
        b = '4'; t = '09'
      case('tcov_dd')
        b = '4'; t = '10'
      case('tcov_db')
        b = '4'; t = '11'
      case('tcov_bd')
        b = '4'; t = '12'
      case('tkeTr_dbb')
        b = '4'; t = '13'
      case('tkeTr_bbd')
        b = '4'; t = '14'
      case('tkeTr_dbd')
        b = '4'; t = '15'
      case('tkeTr_bdd')
        b = '4'; t = '16'
      case('tkeTr_ddd')
        b = '4'; t = '17'
      case('tkeAdv_ddd')
        b = '4'; t = '18'
      case('tkeAdv_ddb')
        b = '4'; t = '19'
      case('tkeAdv_dbb')
        b = '4'; t = '20'
      case('tkeAdv_bdd')
        b = '4'; t = '21'
      case('tkeAdv_bdb')
        b = '4'; t = '22'
      case default
        b = '0'; t = '01'; calc = .true.     
      end select
    end if
  end subroutine define_budget

  ! Utility function to convert integers to strings
  pure function to_string_int(i) result(str)
    integer, intent(in) :: i
    character(len=32) :: str
    write(str, '(I0)') i
  end function to_string_int

  pure function to_string_real(x) result(str)
    real(rk), intent(in) :: x
    character(len=32) :: str
    ! ES format avoids overflow and is portable
    write(str, '(ES16.8)') x
    str = adjustl(str)
  end function to_string_real

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
    ! cmd = "find '"//d_esc//"' -maxdepth 1 \( -type f -o -type l \) -name '"//p_esc//"' -printf '%f\n' > '"//tmpfile//"' 2>/dev/null"
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
  character(len=256) :: field, slice_map
  integer :: nx=1, ny=1, nz=1, runid=1, taskid=0, abl_type=0
  real(rk) :: Lx=1.0_rk, Ly=1.0_rk, Lz=1.0_rk
  real(rk) :: x1=-1._rk, x2=-1._rk, y1=-1._rk, y2=-1._rk, z1=-1._rk, z2=-1._rk
  real(rk) :: lengthscale=1.0_rk
  integer :: nlen, ioUnit=28
  character(:), allocatable :: inputfile
  integer :: budget_source, nxloc, nyloc, nzloc
  character(len=256) :: filename = 'null'
  integer :: start_idx=0, end_idx=huge(1)
  namelist /SETUP/ nx, ny, nz, Lx, Ly, Lz, path, outdir, runid, taskid, field, &
                   slice_map, budget_source, filename, &
                   start_idx, end_idx, abl_type, lengthscale
  namelist /BOX/ x1, x2, y1, y2, z1, z2
      
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
  read(unit=ioUnit, NML=BOX)
  close(ioUnit)

  ! Initiate reader
  call reader%init(nx, ny, nz)
  call reader%local_shape(nxloc, nyloc, nzloc)
  nxloc_=nxloc
  nyloc_=nyloc
  nzloc_=nzloc

  if (taskid == 0)then
    ! Horizontal average
    call ha_driver(reader, Lz, runid, trim(path), trim(outdir), trim(field), budget_source, &
      start_idx, end_idx, filename=trim(filename))
  else if (taskid == 1) then
    ! Slice
    call slice_driver(reader, Lx, Ly, Lz, runid, trim(path), trim(outdir), trim(field), &
      trim(slice_map), budget_source, start_idx, end_idx, trim(filename))
  else if (taskid == -1) then
    ! Miscellaneous tasks
    call miscellaneous_driver(reader, runid, trim(field), trim(path), start_idx, end_idx)
  else if (taskid == 2)then
    call compute_abl(reader, Lx, Ly, Lz, runid, runid-1, budget_source, abl_type, trim(path), trim(outdir), start_idx, end_idx, lengthscale)
  else if (taskid == 3)then
    ! call one_d_profile(reader, Lx, Ly, Lz, runid, budget_source, trim(path), &
    !      trim(outdir), trim(field), start_idx, end_idx, trim(filename), slice_axis, &
    !      x1, x2, y1, y2, z1, z2)
  else if (taskid == 4)then
    call max_time_change(reader, runid, trim(path), trim(field), &
      budget_source, start_idx, end_idx)
  end if
  
  if(myrank == 0) call message('Wrapping up ...')
  call MPI_Finalize(ierr)
end program MPIR3D_
