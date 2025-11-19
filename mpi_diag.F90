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
    real(rk) :: sum_spd(nz), sum_c(nz), sum_s(nz), sum_cnt(nz)
    real(rk) :: uu, vv, mag, cbar, sbar, ang

    sum_spd= 0.0_rk; sum_c= 0.0_rk; sum_s= 0.0_rk; sum_cnt= 0.0_rk

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

    call MPI_Allreduce(MPI_IN_PLACE, sum_spd, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, sum_c,   nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, sum_s,   nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, sum_cnt, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

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
    real(rk) :: sum_f(nz), sum_cnt(nz)

    sum_f= 0.0_rk; sum_cnt(nz)= 0.0_rk
    do k = 1, nzloc
      kg = zs + k - 1
      do j = 1, nyloc
        do i = 1, nxloc                    
          sum_f(kg) = sum_f(kg) + f(i,j,k)        
          sum_cnt(kg) = sum_cnt(kg) + 1.0_rk
        end do
      end do
    end do

    call MPI_Allreduce(MPI_IN_PLACE, sum_f, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, sum_cnt, nz, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

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
          call csvprofile(nz,trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//'HA_WS.csv',z,ws)
          call csvprofile(nz,trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//'HA_WD.csv',z,wd)
        end block
      else
        block
          real(rk) :: profile(nz)
          call do_horizontal_average(nz, nxloc, nyloc, nzloc, zs, f1, profile)
          call csvprofile(nz, trim(outdir)//'/'//'Run'//trim(rc)//'_t'//trim(sorted_keys(k))//&
                'HA_'//trim(f_)//'.csv',z,profile)
        end block
      end if
    end do
  end subroutine ha_driver

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
    case ('T'); name = 'potT'
    case ('p'); name = 'prss'
    case ('S'); name = 'uVel'
    case default; name = 'uVel'
    end select
  end function field_to_name

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
      call message('ERROR(CSV): z and profile size mismatch.')
      call MPI_Abort(MPI_COMM_WORLD, 100, ierr)
    end if

    open(newunit=uo, file=trim(filename), status='replace', action='write')
    write(uo,'(A)') 'z,profile'
    do k = 1, n
      write(uo,'(ES23.15, ",", ES23.15)') z(k), profile(k)
    end do
    close(uo)
  end subroutine csvprofile

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
  integer:: nx=1, ny=1, nz=1, runid=1, taskid=0
  real(rk):: Lx=1.0_rk, Ly=1.0_rk, Lz=1.0_rk
  integer :: nlen, ioUnit=28
  character(:), allocatable :: inputfile
  namelist /SETUP/ nx, ny, nz, Lx, Ly, Lz, path, outdir, runid, taskid, field 
      
  ! Initiate MPI
  ! -------------------------------------------------------------------------!
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count() < 1) then
    call message('Usage: MPIR3D <inputfile>')
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

  ! Initiate reader
  call reader%init(nx, ny, nz)

  if (taskid == 0)then
    ! Horizontal average
    call ha_driver(reader, Lz, runid, path, outdir, field)
  end if
  
  if(myrank == 0) call message('Wrapping up ...')
  call MPI_Finalize(ierr)
end program MPIR3D_
