!===============================================================
! mpi_diag_lean.F90
!
! A clean-input diagnostic program for PadeOps/2DECOMP fields.
!
! This file intentionally does not preserve the legacy keyword-driven
! interface in mpi_diag.F90.  All diagnostics read direct filenames and
! can assemble linear combinations of fields with expressions like:
!
!   {
!     (1.0)(file_a.s3D)
!   + (-2.0)(file_b.s3D)
!   }
!
! Relative filenames are resolved under SETUP%path.  Absolute filenames
! are used as-is.
!===============================================================

module MPIR3D_Lean
  use mpi
  use netcdf
  use decomp_2d_io
  use decomp_2d, only: decomp_info, decomp_2d_init, xstart, &
                       xend, ystart, yend, zstart, zend, &
                       decomp_2d_finalize, decomp_info_init, &
                       mytype, transpose_x_to_y, transpose_y_to_z
  implicit none

  integer, parameter :: rk = mytype
  integer, parameter :: mpi_rk = merge(MPI_DOUBLE_PRECISION, MPI_REAL, rk == kind(1.0d0))
  integer, parameter :: nc_rk = merge(NF90_DOUBLE, NF90_FLOAT, rk == kind(1.0d0))
  integer, parameter :: str_len = 2048
  integer, parameter :: name_len = 128
  integer :: myrank = -1
  integer :: nprocs = -1
  integer :: glob_counter = 0

  type :: FieldReader2Decomp
     private
     integer :: nx = 0, ny = 0, nz = 0
     integer :: p_row = 0, p_col = 0
     logical :: is_init = .false.
     type(DECOMP_INFO) :: gpC
     integer :: xs=0, xe=0, ys=0, ye=0, zs=0, ze=0
     integer :: nxloc=0, nyloc=0, nzloc=0
   contains
     procedure :: init         => frd_init
     procedure :: read_field   => frd_read_field
     procedure :: global_shape => frd_global_shape
     procedure :: local_shape  => frd_local_shape
     procedure :: indices      => frd_index
     procedure, nopass, private :: choose_proc_grid_ => frd_choose_proc_grid
     final :: frd_finalize
  end type FieldReader2Decomp

  type :: field_term_t
    real(rk) :: coeff = 1.0_rk
    character(len=str_len) :: filename = ''
  end type field_term_t

  type :: file_list_t
    integer :: n = 0
    character(len=str_len), allocatable :: names(:)
  end type file_list_t

  type :: field_case_set_t
    integer :: ncases = 0
    type(file_list_t), allocatable :: lists(:)
  end type field_case_set_t

  type :: field_expr_t
    character(len=name_len) :: name = ''
    character(len=32) :: derived = ''
    integer :: nterms = 0
    type(field_term_t), allocatable :: terms(:)
    type(field_term_t), allocatable :: u_terms(:)
    type(field_term_t), allocatable :: v_terms(:)
    integer :: nu_terms = 0
    integer :: nv_terms = 0
  end type field_expr_t

  type :: slice_spec_t
    character(len=1) :: axis = 'z'
    integer :: integrate = 0
    integer :: ncoords = 0
    real(rk), allocatable :: coords(:)
  end type slice_spec_t

  type :: rms_spec_t
    type(field_expr_t) :: field
    character(len=1) :: axis = 'z'
    real(rk) :: bounds_min(3) = [-huge(1.0_rk), -huge(1.0_rk), -huge(1.0_rk)]
    real(rk) :: bounds_max(3) = [ huge(1.0_rk),  huge(1.0_rk),  huge(1.0_rk)]
    logical :: has_min(3) = [.false., .false., .false.]
    logical :: has_max(3) = [.false., .false., .false.]
    logical :: subtract_reference = .false.
    logical :: has_ref_bounds = .false.
    real(rk) :: ref_min = -huge(1.0_rk)
    real(rk) :: ref_max =  huge(1.0_rk)
  end type rms_spec_t

  type :: abl_spec_t
    character(len=32) :: method = ''
    type(field_expr_t) :: theta
    type(field_expr_t) :: uw
    type(field_expr_t) :: vw
    real(rk) :: threshold = 0.05_rk
    real(rk) :: lengthscale = 1.0_rk
    real(rk) :: inversion_l0 = 700.0_rk
    real(rk) :: inversion_d0 = 200.0_rk
    real(rk) :: inversion_xi = 1.3_rk
  end type abl_spec_t

  type :: march_mode_t
    character(len=name_len) :: name = ''
    integer :: nremove = 0
    character(len=name_len), allocatable :: remove(:)
  end type march_mode_t

  type :: march_spec_t
    character(len=1) :: axis = 'x'
    character(len=name_len) :: reference = ''
    character(len=name_len) :: normalizer = ''
    real(rk) :: start = 0.0_rk
    real(rk) :: finish = 0.0_rk
    character(len=32) :: scheme = 'euler'
    logical :: export_rms = .false.
    logical :: export_linear_average = .false.
    integer :: nslices = 0
    real(rk), allocatable :: slices(:)
    integer :: nmodes = 0
    type(march_mode_t), allocatable :: modes(:)
  end type march_spec_t

  type :: rz_params
    real(rk) :: tm = 0.0_rk
    real(rk) :: a = 0.0_rk
    real(rk) :: b = 0.0_rk
    real(rk) :: l = 0.0_rk
    real(rk) :: d = 0.0_rk
    real(rk) :: sse = huge(1.0_rk)
    integer :: status = -1
  end type rz_params

  type :: diag_job_t
    character(len=32) :: driver = ''
    integer :: nfields = 0
    type(field_expr_t), allocatable :: fields(:)
    integer :: nslices = 0
    type(slice_spec_t), allocatable :: slices(:)
    integer :: nrms = 0
    type(rms_spec_t), allocatable :: rms(:)
    integer :: nprofiles = 0
    type(rms_spec_t), allocatable :: profiles(:)
    type(abl_spec_t) :: abl
    type(march_spec_t) :: march
  end type diag_job_t

contains

  subroutine message(msg)
    character(*), intent(in) :: msg
    if (myrank == 0) write(*,*) trim(msg)
  end subroutine message

  pure function lower(s) result(out)
    character(*), intent(in) :: s
    character(len(s)) :: out
    integer :: i, c
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) then
        out(i:i) = achar(c + 32)
      else
        out(i:i) = s(i:i)
      end if
    end do
  end function lower

  pure function strip_comments(line) result(out)
    character(*), intent(in) :: line
    character(len=len(line)) :: out
    integer :: p1, p2, p, first
    out = line
    p1 = index(out, '!')
    first = 1
    do while (first <= len_trim(out))
      if (out(first:first) /= ' ' .and. out(first:first) /= char(9)) exit
      first = first + 1
    end do
    p2 = 0
    if (first <= len_trim(out)) then
      if (out(first:first) == '#') p2 = first
    end if
    p = 0
    if (p1 > 0) p = p1
    if (p2 > 0) then
      if (p == 0) then
        p = p2
      else
        p = min(p, p2)
      end if
    end if
    if (p > 0) out(p:) = ' '
  end function strip_comments

  ! Expects a lower-cased, trimmed parser line.  The character after the key
  ! must be a delimiter, so "name_override" does not match "name".
  pure logical function key_is(line, key)
    character(*), intent(in) :: line, key
    character(len=len(line)) :: tmp
    integer :: n, lt
    tmp = adjustl(trim(line))
    n = len_trim(key)
    lt = len_trim(tmp)
    key_is = .false.
    if (lt < n) return
    if (tmp(1:n) /= key(1:n)) return
    if (lt == n) then
      key_is = .true.
    else
      key_is = tmp(n+1:n+1) == ' ' .or. tmp(n+1:n+1) == char(9) .or. &
               tmp(n+1:n+1) == '=' .or. tmp(n+1:n+1) == '{'
    end if
  end function key_is

  pure logical function is_abs_path(filename)
    character(*), intent(in) :: filename
    is_abs_path = len_trim(filename) > 0 .and. filename(1:1) == '/'
  end function is_abs_path

  function resolve_input_path(root, filename) result(fullpath)
    character(*), intent(in) :: root, filename
    character(len=:), allocatable :: fullpath
    if (is_abs_path(trim(filename))) then
      fullpath = trim(filename)
    else
      fullpath = trim(root)//'/'//trim(filename)
    end if
  end function resolve_input_path

  pure logical function has_wildcard(text)
    character(*), intent(in) :: text
    has_wildcard = index(text, '*') > 0 .or. index(text, '?') > 0 .or. index(text, '[') > 0
  end function has_wildcard

  function dirname_only(filename) result(dir)
    character(*), intent(in) :: filename
    character(len=:), allocatable :: dir
    integer :: i, last, n
    n = len_trim(filename)
    last = 0
    do i = n, 1, -1
      if (filename(i:i) == '/') then
        last = i
        exit
      end if
    end do
    if (last <= 0) then
      dir = '.'
    else if (last == 1) then
      dir = '/'
    else
      dir = filename(:last-1)
    end if
  end function dirname_only

  function shell_quote(text) result(quoted)
    character(*), intent(in) :: text
    character(len=:), allocatable :: quoted
    if (index(text, "'") > 0) then
      quoted = "''"
    else
      quoted = "'"//trim(text)//"'"
    end if
  end function shell_quote

  function basename_only(filename) result(base)
    character(*), intent(in) :: filename
    character(len=:), allocatable :: base
    integer :: i, last, n
    n = len_trim(filename)
    last = 0
    do i = n, 1, -1
      if (filename(i:i) == '/') then
        last = i
        exit
      end if
    end do
    if (last == 0) then
      base = filename(:n)
    else if (last < n) then
      base = filename(last+1:n)
    else
      base = ''
    end if
  end function basename_only

  function strip_extension(filename) result(stem)
    character(*), intent(in) :: filename
    character(len=:), allocatable :: stem
    integer :: i, last_dot, n
    n = len_trim(filename)
    last_dot = 0
    do i = n, 1, -1
      if (filename(i:i) == '.') then
        last_dot = i
        exit
      end if
    end do
    if (last_dot <= 1) then
      stem = filename(:n)
    else
      stem = filename(:last_dot-1)
    end if
  end function strip_extension

  pure logical function is_digit_char(c)
    character(len=1), intent(in) :: c
    is_digit_char = (c >= '0' .and. c <= '9')
  end function is_digit_char

  function legacy_run_time_stem(filename) result(stem)
    character(*), intent(in) :: filename
    character(len=:), allocatable :: stem
    character(len=:), allocatable :: leaf, run_part, time_part
    integer :: n, irun, run_end, itime, time_end, i

    leaf = strip_extension(trim(basename_only(trim(filename))))
    n = len_trim(leaf)
    irun = index(leaf, 'Run')
    if (irun == 0) then
      stem = leaf
      return
    end if

    run_end = n
    do i = irun, n
      if (leaf(i:i) == '_') then
        run_end = i - 1
        exit
      end if
    end do
    run_part = leaf(irun:run_end)

    itime = 0
    do i = 1, n - 2
      if (leaf(i:i+1) == '_t' .and. is_digit_char(leaf(i+2:i+2))) then
        itime = i + 2
      end if
    end do
    if (itime == 0) then
      stem = leaf
      return
    end if

    time_end = n
    do i = itime, n
      if (leaf(i:i) == '_') then
        time_end = i - 1
        exit
      end if
    end do
    time_part = leaf(itime:time_end)
    stem = trim(run_part)//'_t'//trim(time_part)
  end function legacy_run_time_stem

  function expr_output_stem(expr) result(stem)
    type(field_expr_t), intent(in) :: expr
    character(len=:), allocatable :: stem
    character(len=:), allocatable :: leaf
    integer :: i
    if (len_trim(expr%name) > 0) then
      stem = trim(expr%name)
    else if (expr%nterms == 1) then
      leaf = basename_only(trim(expr%terms(1)%filename))
      stem = strip_extension(trim(leaf))
    else
      stem = 'composed'
      do i = 1, expr%nterms
        leaf = basename_only(trim(expr%terms(i)%filename))
        stem = trim(stem)//'_'//trim(strip_extension(trim(leaf)))
      end do
    end if
  end function expr_output_stem

  subroutine append_field(fields, nfields, field)
    type(field_expr_t), allocatable, intent(inout) :: fields(:)
    integer, intent(inout) :: nfields
    type(field_expr_t), intent(in) :: field
    type(field_expr_t), allocatable :: tmp(:)
    integer :: i
    if (.not. allocated(fields)) allocate(fields(0))
    allocate(tmp(nfields + 1))
    do i = 1, nfields
      tmp(i) = fields(i)
    end do
    tmp(nfields + 1) = field
    call move_alloc(tmp, fields)
    nfields = nfields + 1
  end subroutine append_field

  subroutine append_slice(slices, nslices, slice)
    type(slice_spec_t), allocatable, intent(inout) :: slices(:)
    integer, intent(inout) :: nslices
    type(slice_spec_t), intent(in) :: slice
    type(slice_spec_t), allocatable :: tmp(:)
    integer :: i
    if (.not. allocated(slices)) allocate(slices(0))
    allocate(tmp(nslices + 1))
    do i = 1, nslices
      tmp(i) = slices(i)
    end do
    tmp(nslices + 1) = slice
    call move_alloc(tmp, slices)
    nslices = nslices + 1
  end subroutine append_slice

  subroutine append_rms(rms, nrms, item)
    type(rms_spec_t), allocatable, intent(inout) :: rms(:)
    integer, intent(inout) :: nrms
    type(rms_spec_t), intent(in) :: item
    type(rms_spec_t), allocatable :: tmp(:)
    integer :: i
    if (.not. allocated(rms)) allocate(rms(0))
    allocate(tmp(nrms + 1))
    do i = 1, nrms
      tmp(i) = rms(i)
    end do
    tmp(nrms + 1) = item
    call move_alloc(tmp, rms)
    nrms = nrms + 1
  end subroutine append_rms

  subroutine append_profile(profiles, nprofiles, item)
    type(rms_spec_t), allocatable, intent(inout) :: profiles(:)
    integer, intent(inout) :: nprofiles
    type(rms_spec_t), intent(in) :: item
    type(rms_spec_t), allocatable :: tmp(:)
    integer :: i
    if (.not. allocated(profiles)) allocate(profiles(0))
    allocate(tmp(nprofiles + 1))
    do i = 1, nprofiles
      tmp(i) = profiles(i)
    end do
    tmp(nprofiles + 1) = item
    call move_alloc(tmp, profiles)
    nprofiles = nprofiles + 1
  end subroutine append_profile

  subroutine append_march_mode(modes, nmodes, mode)
    type(march_mode_t), allocatable, intent(inout) :: modes(:)
    integer, intent(inout) :: nmodes
    type(march_mode_t), intent(in) :: mode
    type(march_mode_t), allocatable :: tmp(:)
    integer :: i
    if (.not. allocated(modes)) allocate(modes(0))
    allocate(tmp(nmodes + 1))
    do i = 1, nmodes
      tmp(i) = modes(i)
    end do
    tmp(nmodes + 1) = mode
    call move_alloc(tmp, modes)
    nmodes = nmodes + 1
  end subroutine append_march_mode

  subroutine append_term_list(terms, nterms, coeff, filename)
    type(field_term_t), allocatable, intent(inout) :: terms(:)
    integer, intent(inout) :: nterms
    real(rk), intent(in) :: coeff
    character(*), intent(in) :: filename
    type(field_term_t), allocatable :: tmp(:)
    integer :: i
    if (.not. allocated(terms)) allocate(terms(0))
    allocate(tmp(nterms + 1))
    do i = 1, nterms
      tmp(i) = terms(i)
    end do
    tmp(nterms + 1)%coeff = coeff
    tmp(nterms + 1)%filename = trim(filename)
    call move_alloc(tmp, terms)
    nterms = nterms + 1
  end subroutine append_term_list

  subroutine parse_field_expr(text, name, expr, ierr)
    character(*), intent(in) :: text, name
    type(field_expr_t), intent(out) :: expr
    integer, intent(out) :: ierr

    ierr = 0
    expr%name = trim(name)
    expr%nterms = 0
    if (allocated(expr%terms)) deallocate(expr%terms)
    call parse_term_list(text, expr%terms, expr%nterms, ierr)
  end subroutine parse_field_expr

  subroutine parse_term_list(text, terms, nterms, ierr)
    character(*), intent(in) :: text
    type(field_term_t), allocatable, intent(inout) :: terms(:)
    integer, intent(inout) :: nterms
    integer, intent(out) :: ierr
    integer :: pos, n, p1, p2
    character(len=str_len) :: coeff_part, file_part
    real(rk) :: coeff

    ierr = 0
    nterms = 0
    if (allocated(terms)) deallocate(terms)
    allocate(terms(0))

    n = len_trim(text)
    pos = 1
    do
      call skip_blanks_and_plus(text, pos, n)
      if (pos > n) exit

      if (text(pos:pos) /= '(') then
        ierr = 10
        return
      end if
      p1 = index(text(pos+1:n), ')')
      if (p1 <= 0) then
        ierr = 11
        return
      end if
      p1 = pos + p1
      coeff_part = adjustl(text(pos+1:p1-1))
      read(coeff_part, *, iostat=ierr) coeff
      if (ierr /= 0) return
      pos = p1 + 1

      call skip_spaces(text, pos, n)
      if (pos > n .or. text(pos:pos) /= '(') then
        ierr = 12
        return
      end if
      p2 = index(text(pos+1:n), ')')
      if (p2 <= 0) then
        ierr = 13
        return
      end if
      p2 = pos + p2
      file_part = adjustl(text(pos+1:p2-1))
      if (len_trim(file_part) == 0) then
        ierr = 14
        return
      end if
      call append_term_list(terms, nterms, coeff, trim(file_part))
      pos = p2 + 1
    end do

    if (nterms <= 0) ierr = 15
  end subroutine parse_term_list

  subroutine expand_file_pattern(pattern, files, nfiles, ierr)
    character(*), intent(in) :: pattern
    character(len=str_len), allocatable, intent(out) :: files(:)
    integer, intent(out) :: nfiles, ierr
    character(len=:), allocatable :: dir, leaf, qdir, qleaf, qtmp, cmd
    character(len=str_len), allocatable :: filebuf(:)
    character(len=str_len) :: line, tmpfile, tmpdir
    integer :: unit, ios, exitstat, cmdstat, ifile, mpi_ierr, clock_count
    logical :: opened

    ierr = 0
    nfiles = 0
    opened = .false.
    if (allocated(files)) deallocate(files)

    if (index(pattern, "'") > 0) then
      ierr = 610
    end if

    if (ierr == 0 .and. myrank == 0) then
      dir = dirname_only(trim(pattern))
      leaf = basename_only(trim(pattern))
      qdir = shell_quote(trim(dir))
      qleaf = shell_quote(trim(leaf))

      ! Only rank 0 creates temp files; glob_counter is meaningful only on rank 0.
      glob_counter = glob_counter + 1
      call get_environment_variable('TMPDIR', tmpdir)
      if (len_trim(tmpdir) == 0) tmpdir = '/tmp'
      call system_clock(count=clock_count)
      write(tmpfile,'(A,A,I0,A,I0,A,I0,A)') trim(tmpdir), '/mpir3d_glob_', myrank, '_', glob_counter, '_', clock_count, '.lst'
      if (index(tmpfile, "'") > 0) then
        ierr = 610
      else
        qtmp = shell_quote(trim(tmpfile))
        cmd = 'find '//trim(qdir)//' -maxdepth 1 -type f -name '//trim(qleaf)//' | sort > '//trim(qtmp)
        call execute_command_line(trim(cmd), exitstat=exitstat, cmdstat=cmdstat)
        if (cmdstat /= 0 .or. exitstat /= 0) ierr = 611
        if (ierr == 611 .and. cmdstat == 0) then
          open(newunit=unit, file=trim(tmpfile), status='old', action='readwrite', iostat=ios)
          if (ios == 0) close(unit, status='delete')
        end if
      end if

      if (ierr == 0) then
        open(newunit=unit, file=trim(tmpfile), status='old', action='read', iostat=ios)
        if (ios /= 0) then
          ierr = 612
        else
          opened = .true.
        end if
      end if

      if (ierr == 0) then
        do
          read(unit,'(A)',iostat=ios) line
          if (ios /= 0) exit
          if (len_trim(line) > 0) nfiles = nfiles + 1
        end do
        if (nfiles <= 0) ierr = 613
      end if

      if (ierr == 0) then
        rewind(unit)
        allocate(files(nfiles))
        do ifile = 1, nfiles
          read(unit,'(A)',iostat=ios) line
          if (ios /= 0) then
            ierr = 614
            exit
          end if
          files(ifile) = trim(line)
        end do
      end if

      if (opened .and. (ierr == 0 .or. ierr == 613 .or. ierr == 614)) then
        close(unit, status='delete')
      else if (opened) then
        ! Defensive fallback for future error paths after a successful open.
        close(unit)
      end if
    end if

    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    if (ierr /= 0) return

    call MPI_Bcast(nfiles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    allocate(filebuf(nfiles))
    if (myrank == 0) filebuf = files
    call MPI_Bcast(filebuf, nfiles * str_len, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi_ierr)
    if (myrank /= 0) then
      allocate(files(nfiles))
      files = filebuf
    end if
  end subroutine expand_file_pattern

  subroutine expand_terms(root, terms, nterms, lists, ncases, ierr)
    character(*), intent(in) :: root
    type(field_term_t), intent(in) :: terms(:)
    integer, intent(in) :: nterms
    type(file_list_t), allocatable, intent(out) :: lists(:)
    integer, intent(out) :: ncases, ierr
    character(len=:), allocatable :: infile
    integer :: iterm

    ierr = 0
    ncases = 1
    if (allocated(lists)) deallocate(lists)
    if (nterms <= 0) then
      ierr = 621
      return
    end if
    allocate(lists(nterms))

    do iterm = 1, nterms
      infile = resolve_input_path(trim(root), trim(terms(iterm)%filename))
      if (has_wildcard(trim(infile))) then
        call expand_file_pattern(trim(infile), lists(iterm)%names, lists(iterm)%n, ierr)
        if (ierr /= 0) return
      else
        lists(iterm)%n = 1
        allocate(lists(iterm)%names(1))
        lists(iterm)%names(1) = trim(infile)
      end if
      ncases = max(ncases, lists(iterm)%n)
    end do

    do iterm = 1, nterms
      if (lists(iterm)%n /= 1 .and. lists(iterm)%n /= ncases) then
        ierr = 620
        return
      end if
    end do
  end subroutine expand_terms

  subroutine abort_expand_error(ierr_in)
    integer, intent(in) :: ierr_in
    integer :: mpi_ierr
    select case(ierr_in)
    case(610)
      call message('ERROR: wildcard paths containing single quotes are not supported.')
    case(613)
      call message('ERROR: wildcard pattern matched no files.')
    case(620)
      call message('ERROR: wildcard terms in one expression must have either one match or the same number of matches.')
    case(621)
      call message('ERROR: field expression has no input terms.')
    case default
      call message('ERROR: failed to expand wildcard input pattern.')
    end select
    call MPI_Abort(MPI_COMM_WORLD, 800 + ierr_in, mpi_ierr)
  end subroutine abort_expand_error

  subroutine assemble_terms_case(reader, terms, nterms, lists, icase, f, scratch)
    class(FieldReader2Decomp), intent(inout) :: reader
    type(field_term_t), intent(in) :: terms(:)
    integer, intent(in) :: nterms, icase
    type(file_list_t), intent(in) :: lists(:)
    real(rk), intent(out) :: f(:,:,:)
    real(rk), intent(inout) :: scratch(:,:,:)
    integer :: iterm, idx

    f = 0.0_rk
    do iterm = 1, nterms
      idx = merge(1, icase, lists(iterm)%n == 1)
      call message('Reading '//trim(lists(iterm)%names(idx)))
      call reader%read_field(trim(lists(iterm)%names(idx)), scratch)
      f = f + terms(iterm)%coeff * scratch
    end do
  end subroutine assemble_terms_case

  function expr_output_stem_case(expr, lists, nlists, icase) result(stem)
    type(field_expr_t), intent(in) :: expr
    type(file_list_t), intent(in) :: lists(:)
    integer, intent(in) :: nlists, icase
    character(len=:), allocatable :: stem
    character(len=:), allocatable :: leaf
    integer :: i, idx

    do i = 1, nlists
      if (lists(i)%n > 1) then
        ! For multi-wildcard compositions, the first time-varying term names the output case.
        idx = min(icase, lists(i)%n)
        leaf = strip_extension(trim(basename_only(trim(lists(i)%names(idx)))))
        if (len_trim(expr%name) > 0) then
          stem = trim(expr%name)//'_'//trim(leaf)
        else
          stem = trim(leaf)
        end if
        return
      end if
    end do
    stem = expr_output_stem(expr)
  end function expr_output_stem_case

  subroutine skip_spaces(s, pos, n)
    character(*), intent(in) :: s
    integer, intent(inout) :: pos
    integer, intent(in) :: n
    do while (pos <= n)
      if (s(pos:pos) /= ' ' .and. s(pos:pos) /= char(9)) exit
      pos = pos + 1
    end do
  end subroutine skip_spaces

  subroutine skip_blanks_and_plus(s, pos, n)
    character(*), intent(in) :: s
    integer, intent(inout) :: pos
    integer, intent(in) :: n
    do while (pos <= n)
      if (s(pos:pos) /= ' ' .and. s(pos:pos) /= char(9) .and. &
          s(pos:pos) /= '+') exit
      pos = pos + 1
    end do
  end subroutine skip_blanks_and_plus

  function value_after_equals(line) result(value)
    character(*), intent(in) :: line
    character(len=str_len) :: value
    integer :: p
    value = ''
    p = index(line, '=')
    if (p > 0) value = adjustl(trim(line(p+1:)))
  end function value_after_equals

  subroutine parse_real_csv(line, values, nvalues, ierr)
    character(*), intent(in) :: line
    real(rk), allocatable, intent(out) :: values(:)
    integer, intent(out) :: nvalues, ierr
    character(len=str_len) :: tmp, token
    real(rk), allocatable :: vals(:), grown(:)
    integer :: i, start, ntmp
    real(rk) :: val

    ierr = 0
    nvalues = 0
    allocate(vals(0))
    tmp = trim(line)
    start = 1
    do i = 1, len_trim(tmp) + 1
      if (i > len_trim(tmp) .or. tmp(i:i) == ',') then
        token = adjustl(tmp(start:i-1))
        if (len_trim(token) > 0) then
          read(token, *, iostat=ierr) val
          if (ierr /= 0) return
          ntmp = nvalues + 1
          allocate(grown(ntmp))
          if (nvalues > 0) grown(1:nvalues) = vals
          grown(ntmp) = val
          call move_alloc(grown, vals)
          nvalues = ntmp
        end if
        start = i + 1
      end if
    end do
    call move_alloc(vals, values)
  end subroutine parse_real_csv

  subroutine parse_name_csv(line, names, nnames, ierr)
    character(*), intent(in) :: line
    character(len=name_len), allocatable, intent(out) :: names(:)
    integer, intent(out) :: nnames, ierr
    character(len=str_len) :: tmp, token
    character(len=name_len), allocatable :: vals(:), grown(:)
    integer :: i, start, ntmp

    ierr = 0
    nnames = 0
    allocate(vals(0))
    tmp = trim(line)
    if (len_trim(tmp) == 0) then
      call move_alloc(vals, names)
      return
    end if
    start = 1
    do i = 1, len_trim(tmp) + 1
      if (i > len_trim(tmp) .or. tmp(i:i) == ',') then
        token = adjustl(trim(tmp(start:i-1)))
        if (len_trim(token) > 0) then
          ntmp = nnames + 1
          allocate(grown(ntmp))
          if (nnames > 0) grown(1:nnames) = vals
          grown(ntmp) = trim(token)
          call move_alloc(grown, vals)
          nnames = ntmp
        end if
        start = i + 1
      end if
    end do
    call move_alloc(vals, names)
  end subroutine parse_name_csv

  subroutine read_diag_map(filename, job, ierr)
    character(*), intent(in) :: filename
    type(diag_job_t), intent(out) :: job
    integer, intent(out) :: ierr
    integer :: unit, ios
    character(len=str_len) :: line, clean, low

    ierr = 0
    job%nfields = 0; job%nslices = 0; job%nrms = 0; job%nprofiles = 0
    job%march%nmodes = 0; job%march%nslices = 0
    if (allocated(job%fields)) deallocate(job%fields)
    if (allocated(job%slices)) deallocate(job%slices)
    if (allocated(job%rms)) deallocate(job%rms)
    if (allocated(job%profiles)) deallocate(job%profiles)
    if (allocated(job%march%modes)) deallocate(job%march%modes)
    if (allocated(job%march%slices)) deallocate(job%march%slices)
    allocate(job%fields(0), job%slices(0), job%rms(0), job%profiles(0))
    allocate(job%march%modes(0))

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      ierr = ios
      return
    end if

    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      low = lower(clean)

      if (key_is(low, 'driver')) then
        job%driver = trim(lower(value_after_equals(clean)))
      else if (key_is(low, 'march')) then
        call read_march_block(unit, job, ierr)
      else if (key_is(low, 'fields')) then
        call read_fields_block(unit, job, ierr)
      else if (key_is(low, 'modes')) then
        call read_modes_block(unit, job, ierr)
      else if (key_is(low, 'slices')) then
        call read_slice_block(unit, job, ierr)
      else if (key_is(low, 'rms')) then
        call read_rms_block(unit, job, ierr)
      else if (key_is(low, 'profile')) then
        call read_profile_block(unit, job, ierr)
      else if (key_is(low, 'abl')) then
        call read_abl_block(unit, job, ierr)
      end if
      if (ierr /= 0) exit
    end do
    close(unit)

    if (ios > 0 .and. ierr == 0) ierr = ios
  end subroutine read_diag_map

  subroutine read_march_block(unit, job, ierr)
    integer, intent(in) :: unit
    type(diag_job_t), intent(inout) :: job
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean, low, val
    integer :: ios, nvals
    real(rk), allocatable :: vals(:)

    ierr = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      low = lower(clean)
      if (clean(1:1) == '}') exit

      if (key_is(low, 'axis')) then
        val = adjustl(trim(value_after_equals(clean)))
        val = lower(val)
        job%march%axis = val(1:1)
      else if (key_is(low, 'reference')) then
        job%march%reference = trim(value_after_equals(clean))
      else if (key_is(low, 'normalizer')) then
        job%march%normalizer = trim(value_after_equals(clean))
      else if (key_is(low, 'start')) then
        val = value_after_equals(clean)
        read(val, *, iostat=ierr) job%march%start
        if (ierr /= 0) return
      else if (key_is(low, 'end') .or. key_is(low, 'finish')) then
        val = value_after_equals(clean)
        read(val, *, iostat=ierr) job%march%finish
        if (ierr /= 0) return
      else if (key_is(low, 'scheme')) then
        job%march%scheme = trim(lower(value_after_equals(clean)))
      else if (key_is(low, 'rms')) then
        val = value_after_equals(clean)
        call parse_logical_value(trim(val), job%march%export_rms, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'linear_average') .or. key_is(low, 'average')) then
        val = value_after_equals(clean)
        call parse_logical_value(trim(val), job%march%export_linear_average, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'slices')) then
        val = value_after_equals(clean)
        call parse_real_csv(trim(val), vals, nvals, ierr)
        if (ierr /= 0) return
        job%march%nslices = nvals
        if (allocated(job%march%slices)) deallocate(job%march%slices)
        call move_alloc(vals, job%march%slices)
      end if
    end do

    if (.not. any(job%march%axis == ['x','y','z'])) ierr = 702
    if (trim(job%march%scheme) /= 'euler' .and. trim(job%march%scheme) /= 'trapezoid') ierr = 703
  end subroutine read_march_block

  subroutine read_modes_block(unit, job, ierr)
    integer, intent(in) :: unit
    type(diag_job_t), intent(inout) :: job
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean, low, val, remove_text
    integer :: ios
    type(march_mode_t) :: mode

    ierr = 0
    mode = march_mode_t()
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      low = lower(clean)
      if (clean(1:1) == '}') exit

      if (key_is(low, 'name')) then
        if (len_trim(mode%name) > 0) then
          call append_march_mode(job%march%modes, job%march%nmodes, mode)
          mode = march_mode_t()
        end if
        mode%name = trim(value_after_equals(clean))
      else if (key_is(low, 'remove') .or. key_is(low, 'fields')) then
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, remove_text, ierr)
        else
          val = value_after_equals(clean)
          remove_text = trim(val)
        end if
        if (ierr /= 0) return
        call parse_name_csv(trim(remove_text), mode%remove, mode%nremove, ierr)
        if (ierr /= 0) return
      end if
    end do

    if (len_trim(mode%name) > 0) then
      call append_march_mode(job%march%modes, job%march%nmodes, mode)
    end if
    if (job%march%nmodes <= 0) ierr = 704
  end subroutine read_modes_block

  subroutine read_fields_block(unit, job, ierr)
    integer, intent(in) :: unit
    type(diag_job_t), intent(inout) :: job
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean, low, name, expr_text, derived
    integer :: ios
    type(field_expr_t) :: expr

    ierr = 0
    name = ''
    derived = ''
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      low = lower(clean)
      if (clean(1:1) == '}') exit

      if (key_is(low, 'name')) then
        name = trim(value_after_equals(clean))
        derived = ''
        expr = field_expr_t()
      else if (key_is(low, 'derived')) then
        derived = trim(lower(value_after_equals(clean)))
      else if (key_is(low, 'expr')) then
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, expr_text, ierr)
        else
          call read_braced_expr(unit, expr_text, ierr)
        end if
        if (ierr /= 0) return
        call parse_field_expr(trim(expr_text), trim(name), expr, ierr)
        if (ierr /= 0) return
        call append_field(job%fields, job%nfields, expr)
        name = ''
        derived = ''
      else if (key_is(low, 'u')) then
        if (len_trim(derived) == 0) then
          call message('ERROR: u block requires a preceding derived = ... key.')
          ierr = 501
          return
        end if
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, expr_text, ierr)
        else
          call read_braced_expr(unit, expr_text, ierr)
        end if
        if (ierr /= 0) return
        expr%name = trim(name)
        expr%derived = trim(derived)
        call parse_term_list(trim(expr_text), expr%u_terms, expr%nu_terms, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'v')) then
        if (len_trim(derived) == 0) then
          call message('ERROR: v block requires a preceding derived = ... key.')
          ierr = 502
          return
        end if
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, expr_text, ierr)
        else
          call read_braced_expr(unit, expr_text, ierr)
        end if
        if (ierr /= 0) return
        expr%name = trim(name)
        expr%derived = trim(derived)
        call parse_term_list(trim(expr_text), expr%v_terms, expr%nv_terms, ierr)
        if (ierr /= 0) return
        if (trim(expr%derived) == 'ws_wd') then
          if (expr%nu_terms <= 0) then
            call message('ERROR: v block for derived = ws_wd must follow a u block.')
            ierr = 503
            return
          end if
          if (expr%nv_terms <= 0) then
            ierr = 503
            return
          end if
          if (len_trim(expr%name) == 0) then
            call message('WARNING: derived = ws_wd has no name; output stem will be "composed".')
          end if
          call append_field(job%fields, job%nfields, expr)
          name = ''
          derived = ''
          expr = field_expr_t()
        else
          ierr = 504
          return
        end if
      end if
    end do
  end subroutine read_fields_block

  subroutine read_slice_block(unit, job, ierr)
    integer, intent(in) :: unit
    type(diag_job_t), intent(inout) :: job
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean, low, val
    integer :: ios, nvals
    real(rk), allocatable :: vals(:)
    type(slice_spec_t) :: slice

    ierr = 0
    slice%axis = 'z'
    slice%integrate = 0
    slice%ncoords = 0
    if (allocated(slice%coords)) deallocate(slice%coords)

    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      low = lower(clean)
      if (clean(1:1) == '}') exit

      if (key_is(low, 'axis')) then
        val = adjustl(trim(value_after_equals(clean)))
        val = lower(val)
        slice%axis = val(1:1)
      else if (key_is(low, 'coords')) then
        val = value_after_equals(clean)
        call parse_real_csv(trim(val), vals, nvals, ierr)
        if (ierr /= 0) return
        slice%ncoords = nvals
        call move_alloc(vals, slice%coords)
      else if (key_is(low, 'integrate')) then
        val = value_after_equals(clean)
        read(val, *, iostat=ierr) slice%integrate
        if (ierr /= 0) return
      end if
    end do

    if (.not. any(slice%axis == ['x','y','z'])) then
      ierr = 200
      return
    end if
    if (slice%integrate == 0 .and. slice%ncoords == 0) then
      ierr = 201
      return
    end if
    if (slice%integrate > 0 .and. slice%ncoords > 0) then
      call message('WARNING: coords ignored for integrated slice block.')
    end if
    call append_slice(job%slices, job%nslices, slice)
  end subroutine read_slice_block

  subroutine read_rms_block(unit, job, ierr)
    integer, intent(in) :: unit
    type(diag_job_t), intent(inout) :: job
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean, low, val, expr_text, name
    integer :: ios
    type(rms_spec_t) :: item

    ierr = 0
    name = ''
    item%axis = 'z'

    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      low = lower(clean)
      if (clean(1:1) == '}') exit

      if (key_is(low, 'name')) then
        name = trim(value_after_equals(clean))
      else if (key_is(low, 'expr')) then
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, expr_text, ierr)
        else
          call read_braced_expr(unit, expr_text, ierr)
        end if
        if (ierr /= 0) return
        call parse_field_expr(trim(expr_text), trim(name), item%field, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'axis')) then
        val = adjustl(trim(value_after_equals(clean)))
        val = lower(val)
        item%axis = val(1:1)
      else if (key_is(low, 'bounds')) then
        val = value_after_equals(clean)
        call parse_bounds_csv(trim(val), item, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'subtract_reference')) then
        val = value_after_equals(clean)
        call parse_logical_value(trim(val), item%subtract_reference, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'ref_bounds')) then
        val = value_after_equals(clean)
        call parse_ref_bounds_csv(trim(val), item, ierr)
        if (ierr /= 0) return
      end if
    end do

    if (item%field%nterms <= 0) then
      ierr = 301
      return
    end if
    if (.not. any(item%axis == ['x','y','z'])) then
      ierr = 302
      return
    end if
    if (item%subtract_reference .and. .not. item%has_ref_bounds) then
      ierr = 303
      return
    end if
    call append_rms(job%rms, job%nrms, item)
  end subroutine read_rms_block

  subroutine read_profile_block(unit, job, ierr)
    integer, intent(in) :: unit
    type(diag_job_t), intent(inout) :: job
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean, low, val, expr_text, name
    integer :: ios
    type(rms_spec_t) :: item

    ierr = 0
    name = ''
    item%axis = 'z'

    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      low = lower(clean)
      if (clean(1:1) == '}') exit

      if (key_is(low, 'name')) then
        name = trim(value_after_equals(clean))
      else if (key_is(low, 'expr')) then
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, expr_text, ierr)
        else
          call read_braced_expr(unit, expr_text, ierr)
        end if
        if (ierr /= 0) return
        call parse_field_expr(trim(expr_text), trim(name), item%field, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'axis')) then
        val = adjustl(trim(value_after_equals(clean)))
        val = lower(val)
        item%axis = val(1:1)
      else if (key_is(low, 'bounds')) then
        val = value_after_equals(clean)
        call parse_bounds_csv(trim(val), item, ierr)
        if (ierr /= 0) return
      end if
    end do

    if (item%field%nterms <= 0) then
      ierr = 311
      return
    end if
    if (.not. any(item%axis == ['x','y','z'])) then
      ierr = 312
      return
    end if
    call append_profile(job%profiles, job%nprofiles, item)
  end subroutine read_profile_block

  subroutine read_abl_block(unit, job, ierr)
    integer, intent(in) :: unit
    type(diag_job_t), intent(inout) :: job
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean, low, val, expr_text
    integer :: ios

    ierr = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      low = lower(clean)
      if (clean(1:1) == '}') exit

      if (key_is(low, 'method')) then
        job%abl%method = trim(lower(value_after_equals(clean)))
      else if (key_is(low, 'threshold')) then
        val = value_after_equals(clean)
        read(val, *, iostat=ierr) job%abl%threshold
        if (ierr /= 0) return
      else if (key_is(low, 'lengthscale')) then
        val = value_after_equals(clean)
        read(val, *, iostat=ierr) job%abl%lengthscale
        if (ierr /= 0) return
      else if (key_is(low, 'inversion_l0') .or. key_is(low, 'l0')) then
        val = value_after_equals(clean)
        read(val, *, iostat=ierr) job%abl%inversion_l0
        if (ierr /= 0) return
      else if (key_is(low, 'inversion_d0') .or. key_is(low, 'd0')) then
        val = value_after_equals(clean)
        read(val, *, iostat=ierr) job%abl%inversion_d0
        if (ierr /= 0) return
      else if (key_is(low, 'inversion_xi') .or. key_is(low, 'xi')) then
        val = value_after_equals(clean)
        read(val, *, iostat=ierr) job%abl%inversion_xi
        if (ierr /= 0) return
      else if (key_is(low, 'theta')) then
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, expr_text, ierr)
        else
          call read_braced_expr(unit, expr_text, ierr)
        end if
        if (ierr /= 0) return
        call parse_field_expr(trim(expr_text), 'theta', job%abl%theta, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'uw')) then
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, expr_text, ierr)
        else
          call read_braced_expr(unit, expr_text, ierr)
        end if
        if (ierr /= 0) return
        call parse_field_expr(trim(expr_text), 'uw', job%abl%uw, ierr)
        if (ierr /= 0) return
      else if (key_is(low, 'vw')) then
        if (index(clean, '{') > 0) then
          call read_expr_body(unit, expr_text, ierr)
        else
          call read_braced_expr(unit, expr_text, ierr)
        end if
        if (ierr /= 0) return
        call parse_field_expr(trim(expr_text), 'vw', job%abl%vw, ierr)
        if (ierr /= 0) return
      end if
    end do
  end subroutine read_abl_block

  subroutine read_braced_expr(unit, expr_text, ierr)
    integer, intent(in) :: unit
    character(len=str_len), intent(out) :: expr_text
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean
    integer :: ios

    ierr = 0
    expr_text = ''
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      if (clean(1:1) == '{') exit
    end do

    call read_expr_body(unit, expr_text, ierr)
  end subroutine read_braced_expr

  subroutine read_expr_body(unit, expr_text, ierr)
    integer, intent(in) :: unit
    character(len=str_len), intent(out) :: expr_text
    integer, intent(out) :: ierr
    character(len=str_len) :: line, clean
    integer :: ios

    ierr = 0
    expr_text = ''
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) then
        ierr = ios
        return
      end if
      clean = adjustl(trim(strip_comments(line)))
      if (len_trim(clean) == 0) cycle
      if (clean(1:1) == '}') exit
      expr_text = trim(expr_text)//' '//trim(clean)
    end do
  end subroutine read_expr_body

  subroutine parse_bounds_csv(line, item, ierr)
    character(*), intent(in) :: line
    type(rms_spec_t), intent(inout) :: item
    integer, intent(out) :: ierr
    character(len=str_len) :: tmp, token
    real(rk) :: value
    integer :: i, start, ntok, iax

    ierr = 0
    item%bounds_min = [-huge(1.0_rk), -huge(1.0_rk), -huge(1.0_rk)]
    item%bounds_max = [ huge(1.0_rk),  huge(1.0_rk),  huge(1.0_rk)]
    item%has_min = [.false., .false., .false.]
    item%has_max = [.false., .false., .false.]
    tmp = trim(line)
    start = 1
    ntok = 0
    do i = 1, len_trim(tmp) + 1
      if (i > len_trim(tmp) .or. tmp(i:i) == ',') then
        token = adjustl(trim(tmp(start:i-1)))
        if (len_trim(token) > 0) then
          ntok = ntok + 1
          if (ntok > 6) then
            ierr = 300
            return
          end if
          if (trim(token) /= '*') then
            read(token, *, iostat=ierr) value
            if (ierr /= 0) return
            iax = (ntok + 1) / 2
            if (mod(ntok, 2) == 1) then
              item%bounds_min(iax) = value
              item%has_min(iax) = .true.
            else
              item%bounds_max(iax) = value
              item%has_max(iax) = .true.
            end if
          end if
        end if
        start = i + 1
      end if
    end do
    if (ntok /= 6) then
      ierr = 300
      return
    end if
    do iax = 1, 3
      if (item%has_min(iax) .and. item%has_max(iax)) then
        if (item%bounds_min(iax) > item%bounds_max(iax)) ierr = 400 + iax
      end if
    end do
  end subroutine parse_bounds_csv

  subroutine parse_ref_bounds_csv(line, item, ierr)
    character(*), intent(in) :: line
    type(rms_spec_t), intent(inout) :: item
    integer, intent(out) :: ierr
    character(len=str_len) :: tmp, token
    real(rk) :: value
    integer :: i, start, ntok

    ierr = 0
    item%ref_min = -huge(1.0_rk)
    item%ref_max =  huge(1.0_rk)
    item%has_ref_bounds = .true.
    tmp = trim(line)
    start = 1
    ntok = 0
    do i = 1, len_trim(tmp) + 1
      if (i > len_trim(tmp) .or. tmp(i:i) == ',') then
        token = adjustl(trim(tmp(start:i-1)))
        if (len_trim(token) > 0) then
          ntok = ntok + 1
          if (ntok > 2) then
            ierr = 304
            return
          end if
          if (trim(token) /= '*') then
            read(token, *, iostat=ierr) value
            if (ierr /= 0) return
            if (ntok == 1) then
              item%ref_min = value
            else
              item%ref_max = value
            end if
          end if
        end if
        start = i + 1
      end if
    end do
    if (ntok /= 2) then
      ierr = 304
      return
    end if
    if (item%ref_min > item%ref_max) ierr = 305
  end subroutine parse_ref_bounds_csv

  subroutine parse_logical_value(text, value, ierr)
    character(*), intent(in) :: text
    logical, intent(out) :: value
    integer, intent(out) :: ierr
    character(len=str_len) :: val

    ierr = 0
    val = trim(lower(adjustl(text)))
    select case(trim(val))
    case('true', 't', '.true.', '1', 'yes', 'y', 'on')
      value = .true.
    case('false', 'f', '.false.', '0', 'no', 'n', 'off')
      value = .false.
    case default
      ierr = 306
    end select
  end subroutine parse_logical_value

  pure logical function coord_in_bounds(coord, item, iax)
    real(rk), intent(in) :: coord
    type(rms_spec_t), intent(in) :: item
    integer, intent(in) :: iax
    coord_in_bounds = coord >= item%bounds_min(iax) .and. coord <= item%bounds_max(iax)
  end function coord_in_bounds

  pure logical function coord_in_ref_bounds(coord, item)
    real(rk), intent(in) :: coord
    type(rms_spec_t), intent(in) :: item
    coord_in_ref_bounds = coord >= item%ref_min .and. coord <= item%ref_max
  end function coord_in_ref_bounds

  subroutine create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)
    real(rk), intent(in) :: Lx, Ly, Lz
    integer, intent(in) :: nx, ny, nz
    real(rk), intent(out) :: x(nx), y(ny), z(nz)
    integer :: i
    do i = 1, nx
      x(i) = real(i-1, rk) * Lx / real(nx, rk)
    end do
    do i = 1, ny
      y(i) = real(i-1, rk) * Ly / real(ny, rk)
    end do
    do i = 1, nz
      z(i) = real(i, rk) * Lz / real(nz, rk) - 0.5_rk * Lz / real(nz, rk)
    end do
  end subroutine create_grid

  subroutine do_horizontal_average(nx, ny, nz, nxloc, nyloc, nzloc, zs, f, profile)
    integer, intent(in) :: nx, ny, nz, nxloc, nyloc, nzloc, zs
    real(rk), intent(in) :: f(nxloc,nyloc,nzloc)
    real(rk), intent(out) :: profile(nz)
    real(rk), allocatable :: local(:)
    integer :: i, j, k, kg, ierr
    allocate(local(nz))
    local = 0.0_rk
    do k = 1, nzloc
      kg = zs + k - 1
      do j = 1, nyloc
        do i = 1, nxloc
          local(kg) = local(kg) + f(i,j,k)
        end do
      end do
    end do
    call MPI_Allreduce(local, profile, nz, mpi_rk, MPI_SUM, MPI_COMM_WORLD, ierr)
    profile = profile / real(nx * ny, rk)
  end subroutine do_horizontal_average

  subroutine run_horizontal_average(reader, Lx, Ly, Lz, path, outdir, job)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: path, outdir
    type(diag_job_t), intent(in) :: job
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze, i, k
    integer :: icase, ncases, ncases_u, ncases_v, ierr
    real(rk), allocatable :: f(:,:,:), f2(:,:,:), scratch(:,:,:), x(:), y(:), z(:), profile(:), profile_u(:), profile_v(:)
    type(file_list_t), allocatable :: lists(:), u_lists(:), v_lists(:)
    character(len=:), allocatable :: outname, stem
    real(rk) :: pi

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    allocate(f(nxloc,nyloc,nzloc), f2(nxloc,nyloc,nzloc), scratch(nxloc,nyloc,nzloc), &
      x(nx), y(ny), z(nz), profile(nz), profile_u(nz), profile_v(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)
    pi = acos(-1.0_rk)

    do i = 1, job%nfields
      if (trim(job%fields(i)%derived) == 'ws_wd') then
        call expand_terms(path, job%fields(i)%u_terms, job%fields(i)%nu_terms, u_lists, ncases_u, ierr)
        if (ierr /= 0) call abort_expand_error(ierr)
        call expand_terms(path, job%fields(i)%v_terms, job%fields(i)%nv_terms, v_lists, ncases_v, ierr)
        if (ierr /= 0) call abort_expand_error(ierr)
        ncases = max(ncases_u, ncases_v)
        if ((ncases_u /= 1 .and. ncases_u /= ncases) .or. (ncases_v /= 1 .and. ncases_v /= ncases)) then
          call message('ERROR: ws_wd wildcard u/v series must have the same number of matches.')
          call MPI_Abort(MPI_COMM_WORLD, 821, ierr)
        end if

        do icase = 1, ncases
          call assemble_terms_case(reader, job%fields(i)%u_terms, job%fields(i)%nu_terms, u_lists, icase, f, scratch)
          call assemble_terms_case(reader, job%fields(i)%v_terms, job%fields(i)%nv_terms, v_lists, icase, f2, scratch)
          stem = expr_output_stem_case(job%fields(i), u_lists, job%fields(i)%nu_terms, icase)

          ! WS profile is <sqrt(u*u + v*v)>, not sqrt(<u>**2 + <v>**2).
          scratch = sqrt(f*f + f2*f2)
          call do_horizontal_average(nx, ny, nz, nxloc, nyloc, nzloc, zs, scratch, profile)
          if (myrank == 0) then
            outname = trim(outdir)//'/'//trim(stem)//'_WS_HA_z.csv'
            call csv_profile(nz, trim(outname), 'z', z, profile)
          end if

          call do_horizontal_average(nx, ny, nz, nxloc, nyloc, nzloc, zs, f, profile_u)
          call do_horizontal_average(nx, ny, nz, nxloc, nyloc, nzloc, zs, f2, profile_v)
          ! Mathematical direction: counter-clockwise degrees from +x, not meteorological direction.
          do k = 1, nz
            profile(k) = atan2(profile_v(k), profile_u(k)) * 180.0_rk / pi
          end do
          if (myrank == 0) then
            outname = trim(outdir)//'/'//trim(stem)//'_WD_HA_z.csv'
            call csv_profile(nz, trim(outname), 'z', z, profile)
          end if
        end do
      else
        call expand_terms(path, job%fields(i)%terms, job%fields(i)%nterms, lists, ncases, ierr)
        if (ierr /= 0) call abort_expand_error(ierr)
        do icase = 1, ncases
          call assemble_terms_case(reader, job%fields(i)%terms, job%fields(i)%nterms, lists, icase, f, scratch)
          call do_horizontal_average(nx, ny, nz, nxloc, nyloc, nzloc, zs, f, profile)
          if (myrank == 0) then
            stem = expr_output_stem_case(job%fields(i), lists, job%fields(i)%nterms, icase)
            outname = trim(outdir)//'/'//trim(stem)//'_HA_z.csv'
            call csv_profile(nz, trim(outname), 'z', z, profile)
          end if
        end do
      end if
    end do
  end subroutine run_horizontal_average

  subroutine run_rms(reader, Lx, Ly, Lz, path, outdir, job)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: path, outdir
    type(diag_job_t), intent(in) :: job
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze
    integer :: ip, i, j, k, ig, jg, kg, ni, profile_axis, ierr, icase, ncases
    integer :: ref_size, ref_idx
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), pointer :: coord(:)
    real(rk), allocatable :: f(:,:,:), scratch(:,:,:), l2norm(:)
    real(rk), allocatable :: local_ref_sum(:), global_ref_sum(:), local_ref_count(:), global_ref_count(:), ref(:)
    type(file_list_t), allocatable :: lists(:)
    character(len=:), allocatable :: outname, stem, rms_label

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    allocate(f(nxloc,nyloc,nzloc), scratch(nxloc,nyloc,nzloc), x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)

    do ip = 1, job%nrms
      call expand_terms(path, job%rms(ip)%field%terms, job%rms(ip)%field%nterms, lists, ncases, ierr)
      if (ierr /= 0) call abort_expand_error(ierr)
      select case(job%rms(ip)%axis)
      case('x')
        ni = nx; coord => x; profile_axis = 1
      case('y')
        ni = ny; coord => y; profile_axis = 2
      case default
        ni = nz; coord => z; profile_axis = 3
      end select
      allocate(l2norm(ni))
      do icase = 1, ncases
        call assemble_terms_case(reader, job%rms(ip)%field%terms, job%rms(ip)%field%nterms, lists, icase, f, scratch)
        if (job%rms(ip)%subtract_reference) then
          select case(job%rms(ip)%axis)
          case('x')
            ref_size = ny * nz
          case('y')
            ref_size = nx * nz
          case default
            ref_size = nx * ny
          end select
          allocate(local_ref_sum(ref_size), global_ref_sum(ref_size), &
            local_ref_count(ref_size), global_ref_count(ref_size), ref(ref_size))
          local_ref_sum = 0.0_rk
          local_ref_count = 0.0_rk

          select case(job%rms(ip)%axis)
          case('x')
            do i = 1, nxloc
              ig = xs + i - 1
              if (.not. coord_in_ref_bounds(x(ig), job%rms(ip))) cycle
              do j = 1, nyloc
                jg = ys + j - 1
                if (.not. coord_in_bounds(y(jg), job%rms(ip), 2)) cycle
                do k = 1, nzloc
                  kg = zs + k - 1
                  if (.not. coord_in_bounds(z(kg), job%rms(ip), 3)) cycle
                  ref_idx = (kg - 1) * ny + jg
                  local_ref_sum(ref_idx) = local_ref_sum(ref_idx) + f(i,j,k)
                  local_ref_count(ref_idx) = local_ref_count(ref_idx) + 1.0_rk
                end do
              end do
            end do
          case('y')
            do j = 1, nyloc
              jg = ys + j - 1
              if (.not. coord_in_ref_bounds(y(jg), job%rms(ip))) cycle
              do i = 1, nxloc
                ig = xs + i - 1
                if (.not. coord_in_bounds(x(ig), job%rms(ip), 1)) cycle
                do k = 1, nzloc
                  kg = zs + k - 1
                  if (.not. coord_in_bounds(z(kg), job%rms(ip), 3)) cycle
                  ref_idx = (kg - 1) * nx + ig
                  local_ref_sum(ref_idx) = local_ref_sum(ref_idx) + f(i,j,k)
                  local_ref_count(ref_idx) = local_ref_count(ref_idx) + 1.0_rk
                end do
              end do
            end do
          case default
            do k = 1, nzloc
              kg = zs + k - 1
              if (.not. coord_in_ref_bounds(z(kg), job%rms(ip))) cycle
              do i = 1, nxloc
                ig = xs + i - 1
                if (.not. coord_in_bounds(x(ig), job%rms(ip), 1)) cycle
                do j = 1, nyloc
                  jg = ys + j - 1
                  if (.not. coord_in_bounds(y(jg), job%rms(ip), 2)) cycle
                  ref_idx = (jg - 1) * nx + ig
                  local_ref_sum(ref_idx) = local_ref_sum(ref_idx) + f(i,j,k)
                  local_ref_count(ref_idx) = local_ref_count(ref_idx) + 1.0_rk
                end do
              end do
            end do
          end select

          call MPI_Allreduce(local_ref_sum, global_ref_sum, ref_size, mpi_rk, MPI_SUM, MPI_COMM_WORLD, ierr)
          call MPI_Allreduce(local_ref_count, global_ref_count, ref_size, mpi_rk, MPI_SUM, MPI_COMM_WORLD, ierr)
          if (maxval(global_ref_count) <= 0.0_rk) then
            call message('ERROR: subtract_reference ref_bounds did not select any cells.')
            call MPI_Abort(MPI_COMM_WORLD, 307, ierr)
          end if
          ref = 0.0_rk
          do ref_idx = 1, ref_size
            if (global_ref_count(ref_idx) > 0.0_rk) then
              ref(ref_idx) = global_ref_sum(ref_idx) / global_ref_count(ref_idx)
            end if
          end do
          call subtract_crossplane_reference(reader, f, job%rms(ip), ref)
        end if

        ! This driver exports the cross-plane L2 norm, sqrt(integral f**2 dA).
        ! Area normalization to convert this profile to RMS is done offline.
        call l2_profile_bounded(reader, x, y, z, Lx, Ly, Lz, f, job%rms(ip), l2norm)
        if (myrank == 0) then
          stem = expr_output_stem_case(job%rms(ip)%field, lists, job%rms(ip)%field%nterms, icase)
          if (job%rms(ip)%subtract_reference) then
            rms_label = 'delta_rms'
          else
            rms_label = 'rms'
          end if
          outname = trim(outdir)//'/'//trim(stem)//'_'//trim(rms_label)//'_'//job%rms(ip)%axis//'.csv'
          call csv_profile_bounded(ni, trim(outname), job%rms(ip)%axis, coord, l2norm, &
            job%rms(ip), profile_axis, 'l2_norm')
        end if
        if (job%rms(ip)%subtract_reference) then
          deallocate(local_ref_sum, global_ref_sum, local_ref_count, global_ref_count, ref)
        end if
      end do
      deallocate(l2norm)
    end do
  end subroutine run_rms

  subroutine run_linear_profiles(reader, Lx, Ly, Lz, path, outdir, job)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: path, outdir
    type(diag_job_t), intent(in) :: job
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze
    integer :: ip, ni, profile_axis, ierr, icase, ncases
    real(rk) :: dx, dy, dz, area_weight
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), pointer :: coord(:)
    real(rk), allocatable :: f(:,:,:), scratch(:,:,:)
    real(rk), allocatable :: avg(:)
    type(file_list_t), allocatable :: lists(:)
    character(len=:), allocatable :: outname, stem

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    allocate(f(nxloc,nyloc,nzloc), scratch(nxloc,nyloc,nzloc), x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)
    dx = Lx/real(nx,rk); dy = Ly/real(ny,rk); dz = Lz/real(nz,rk)

    do ip = 1, job%nprofiles
      call expand_terms(path, job%profiles(ip)%field%terms, job%profiles(ip)%field%nterms, lists, ncases, ierr)
      if (ierr /= 0) call abort_expand_error(ierr)
      select case(job%profiles(ip)%axis)
      case('x')
        ni = nx; coord => x; area_weight = dy*dz; profile_axis = 1
      case('y')
        ni = ny; coord => y; area_weight = dx*dz; profile_axis = 2
      case default
        ni = nz; coord => z; area_weight = dx*dy; profile_axis = 3
      end select
      do icase = 1, ncases
        call assemble_terms_case(reader, job%profiles(ip)%field%terms, job%profiles(ip)%field%nterms, lists, icase, f, scratch)
        allocate(avg(ni))
        call compute_bounded_linear_average(reader, x, y, z, Lx, Ly, Lz, f, job%profiles(ip), avg)

        if (myrank == 0) then
          stem = expr_output_stem_case(job%profiles(ip)%field, lists, job%profiles(ip)%field%nterms, icase)
          outname = trim(outdir)//'/'//trim(stem)//'_avg_'//job%profiles(ip)%axis//'.csv'
          call csv_profile_bounded(ni, trim(outname), job%profiles(ip)%axis, coord, avg, &
            job%profiles(ip), profile_axis, 'average')
        end if
        deallocate(avg)
      end do
    end do
  end subroutine run_linear_profiles

  subroutine run_slices(reader, Lx, Ly, Lz, path, outdir, job)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: path, outdir
    type(diag_job_t), intent(in) :: job
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, ifield, ispec, ierr, icase, ncases
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), allocatable :: f(:,:,:), scratch(:,:,:)
    type(file_list_t), allocatable :: lists(:)
    character(len=:), allocatable :: stem

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    allocate(f(nxloc,nyloc,nzloc), scratch(nxloc,nyloc,nzloc), x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)

    do ifield = 1, job%nfields
      if (len_trim(job%fields(ifield)%derived) > 0) then
        call message('ERROR: derived fields are not supported by the slice driver.')
        call MPI_Abort(MPI_COMM_WORLD, 811, ierr)
      end if
      call expand_terms(path, job%fields(ifield)%terms, job%fields(ifield)%nterms, lists, ncases, ierr)
      if (ierr /= 0) call abort_expand_error(ierr)
      do icase = 1, ncases
        call assemble_terms_case(reader, job%fields(ifield)%terms, job%fields(ifield)%nterms, lists, icase, f, scratch)
        stem = expr_output_stem_case(job%fields(ifield), lists, job%fields(ifield)%nterms, icase)
        do ispec = 1, job%nslices
          call process_slice_spec(reader, f, x, y, z, Lx, Ly, Lz, outdir, trim(stem), job%slices(ispec))
        end do
      end do
    end do
  end subroutine run_slices

  subroutine process_slice_spec(reader, f, x, y, z, Lx, Ly, Lz, outdir, stem, spec)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in), target :: f(:,:,:)
    real(rk), intent(in), target :: x(:), y(:), z(:)
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: outdir, stem
    type(slice_spec_t), intent(in) :: spec
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze
    integer :: nx1, nx2, nax, x1s, x1e, x2s, x2e, axs, axe
    integer :: icoord, j, jloc, k0, k1, k0loc, k1loc, ierr
    real(rk) :: delta, alpha, coord_value
    character(len=1) :: x1name, x2name, eax
    real(rk), pointer :: x1(:), x2(:), axis_coord(:)
    real(rk), pointer :: plane(:,:)
    real(rk), allocatable :: local0(:,:), local1(:,:), global0(:,:), global1(:,:), interp(:,:)
    character(len=:), allocatable :: fname

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)

    select case(spec%axis)
    case('x')
      x1name='y'; x2name='z'; nx1=ny; nx2=nz; nax=nx
      x1s=ys; x1e=ye; x2s=zs; x2e=ze; axs=xs; axe=xe; eax='i'
      x1 => y; x2 => z; axis_coord => x
    case('y')
      x1name='x'; x2name='z'; nx1=nx; nx2=nz; nax=ny
      x1s=xs; x1e=xe; x2s=zs; x2e=ze; axs=ys; axe=ye; eax='j'
      x1 => x; x2 => z; axis_coord => y
    case default
      x1name='x'; x2name='y'; nx1=nx; nx2=ny; nax=nz
      x1s=xs; x1e=xe; x2s=ys; x2e=ye; axs=zs; axe=ze; eax='k'
      x1 => x; x2 => y; axis_coord => z
    end select

    allocate(local0(nx1,nx2), global0(nx1,nx2))

    if (spec%integrate > 0) then
      local0 = 0.0_rk
      select case(spec%axis)
      case('x'); delta = Lx / real(nax, rk)
      case('y'); delta = Ly / real(nax, rk)
      case default; delta = Lz / real(nax, rk)
      end select
      do j = axs, axe
        jloc = j - axs + 1
        select case(spec%axis)
        case('x'); plane => f(jloc,:,:)
        case('y'); plane => f(:,jloc,:)
        case default; plane => f(:,:,jloc)
        end select
        local0(x1s:x1e,x2s:x2e) = local0(x1s:x1e,x2s:x2e) + plane**spec%integrate * delta
      end do
      call MPI_Reduce(local0, global0, nx1*nx2, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (myrank == 0) then
        global0 = global0**(1.0_rk / real(spec%integrate, rk))
        fname = trim(outdir)//'/'//trim(stem)//'_integ_'//spec%axis//'_'//trim(int_to_string(spec%integrate))//'.nc'
        call export_slice_to_netcdf(trim(fname), trim(stem), global0, x1, x2, x1name, x2name)
      end if
      return
    end if

    allocate(local1(nx1,nx2), global1(nx1,nx2), interp(nx1,nx2))
    do icoord = 1, spec%ncoords
      coord_value = spec%coords(icoord)
      local0 = 0.0_rk; local1 = 0.0_rk
      global0 = 0.0_rk; global1 = 0.0_rk; interp = 0.0_rk
      call find_bracket_uniform(axis_coord, coord_value, k0, k1, alpha)

      if (k0 >= axs .and. k0 <= axe) then
        k0loc = k0 - axs + 1
        select case(spec%axis)
        case('x'); plane => f(k0loc,:,:)
        case('y'); plane => f(:,k0loc,:)
        case default; plane => f(:,:,k0loc)
        end select
        local0(x1s:x1e,x2s:x2e) = plane
      end if
      if (k1 >= axs .and. k1 <= axe) then
        k1loc = k1 - axs + 1
        select case(spec%axis)
        case('x'); plane => f(k1loc,:,:)
        case('y'); plane => f(:,k1loc,:)
        case default; plane => f(:,:,k1loc)
        end select
        local1(x1s:x1e,x2s:x2e) = plane
      end if

      call MPI_Reduce(local0, global0, nx1*nx2, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(local1, global1, nx1*nx2, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      if (myrank == 0) then
        interp = (1.0_rk - alpha) * global0 + alpha * global1
        fname = trim(outdir)//'/'//trim(stem)//'_SL_'//spec%axis
        if (coord_value <= -1.0_rk) fname = trim(fname)//'_'//eax
        fname = trim(fname)//'='//trim(real_to_tag(coord_value))//'.nc'
        call export_slice_to_netcdf(trim(fname), trim(stem), interp, x1, x2, x1name, x2name)
      end if
    end do
  end subroutine process_slice_spec

  subroutine run_abl(reader, Lx, Ly, Lz, path, outdir, job)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: path, outdir
    type(diag_job_t), intent(in) :: job
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, icase, ncases, ncases_uw, ncases_vw, ierr
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), allocatable :: uw(:,:,:), vw(:,:,:), buffer(:,:,:), scratch(:,:,:), yscratch(:,:,:), zscratch(:,:,:)
    type(file_list_t), allocatable :: lists(:), uw_lists(:), vw_lists(:)
    character(len=:), allocatable :: stem

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    allocate(x(nx), y(ny), z(nz), buffer(nxloc,nyloc,nzloc), scratch(nxloc,nyloc,nzloc), &
      yscratch(reader%gpC%ysz(1), reader%gpC%ysz(2), reader%gpC%ysz(3)), &
      zscratch(reader%gpC%zsz(1), reader%gpC%zsz(2), reader%gpC%zsz(3)))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)

    select case(trim(lower(job%abl%method)))
    case('stress')
      allocate(uw(nxloc,nyloc,nzloc), vw(nxloc,nyloc,nzloc))
      call expand_terms(path, job%abl%uw%terms, job%abl%uw%nterms, uw_lists, ncases_uw, ierr)
      if (ierr /= 0) call abort_expand_error(ierr)
      call expand_terms(path, job%abl%vw%terms, job%abl%vw%nterms, vw_lists, ncases_vw, ierr)
      if (ierr /= 0) call abort_expand_error(ierr)
      ncases = max(ncases_uw, ncases_vw)
      if ((ncases_uw /= 1 .and. ncases_uw /= ncases) .or. (ncases_vw /= 1 .and. ncases_vw /= ncases)) then
        call message('ERROR: ABL stress wildcard uw/vw series must have the same number of matches.')
        call MPI_Abort(MPI_COMM_WORLD, 831, ierr)
      end if
      do icase = 1, ncases
        call assemble_terms_case(reader, job%abl%uw%terms, job%abl%uw%nterms, uw_lists, icase, uw, scratch)
        call assemble_terms_case(reader, job%abl%vw%terms, job%abl%vw%nterms, vw_lists, icase, vw, scratch)
        uw = sqrt(uw**2 + vw**2)
        if (ncases > 1) then
          stem = 'ABL_stress_'//trim(expr_output_stem_case(job%abl%uw, uw_lists, job%abl%uw%nterms, icase))
        else
          stem = 'ABL_stress'
        end if
        call export_abl_threshold(reader, uw, x, y, z, outdir, trim(stem), job%abl%threshold, yscratch, zscratch)
      end do
    case('inversion')
      call expand_terms(path, job%abl%theta%terms, job%abl%theta%nterms, lists, ncases, ierr)
      if (ierr /= 0) call abort_expand_error(ierr)
      do icase = 1, ncases
        call assemble_terms_case(reader, job%abl%theta%terms, job%abl%theta%nterms, lists, icase, buffer, scratch)
        stem = legacy_run_time_stem(lists(1)%names(merge(1, icase, lists(1)%n == 1)))
        call export_abl_inversion(reader, buffer, x, y, z, outdir, trim(stem), job%abl%lengthscale, &
          job%abl%inversion_l0, job%abl%inversion_d0, job%abl%inversion_xi, yscratch, zscratch)
      end do
    case default
      call message('ERROR: ABL method must be stress or inversion.')
    end select
  end subroutine run_abl

  subroutine export_abl_threshold(reader, field, x, y, z, outdir, stem, threshold, ytmp, ztmp)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: field(:,:,:), x(:), y(:), z(:), threshold
    character(*), intent(in) :: outdir, stem
    real(rk), intent(inout) :: ytmp(:,:,:), ztmp(:,:,:)
    real(rk), allocatable :: local(:,:), global(:,:)
    integer :: nx, ny, nz, ierr
    character(len=:), allocatable :: fname

    call reader%global_shape(nx, ny, nz)
    if (nz /= size(z)) call message('WARNING: ABL threshold z-grid length differs from global nz.')
    allocate(local(nx,ny), global(nx,ny))
    local = 0.0_rk; global = 0.0_rk
    call require_shape3('ABL threshold x-pencil field', shape(field), reader%gpC%xsz)
    call require_shape3('ABL threshold y-pencil scratch', shape(ytmp), reader%gpC%ysz)
    call require_shape3('ABL threshold z-pencil scratch', shape(ztmp), reader%gpC%zsz)
    call transpose_x_to_y(field, ytmp, reader%gpC)
    call transpose_y_to_z(ytmp, ztmp, reader%gpC)
    call find_threshold_crossing_z(ztmp, z, &
      local(reader%gpC%zst(1):reader%gpC%zst(1)+reader%gpC%zsz(1)-1, &
            reader%gpC%zst(2):reader%gpC%zst(2)+reader%gpC%zsz(2)-1), threshold, -huge(1.0_rk))
    call MPI_Reduce(local, global, nx*ny, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) then
      fname = trim(outdir)//'/'//trim(stem)//'.nc'
      call export_slice_to_netcdf(trim(fname), trim(stem), global, x, y, 'x', 'y')
    end if
  end subroutine export_abl_threshold

  subroutine export_abl_inversion(reader, theta, x, y, z, outdir, stem, lengthscale, l0, d0, xi, ytmp, ztmp)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: theta(:,:,:), x(:), y(:), z(:), lengthscale
    real(rk), intent(in) :: l0, d0, xi
    character(*), intent(in) :: outdir, stem
    integer :: nx, ny, nz, ic, jc, ig, jg, ierr
    real(rk), intent(inout) :: ytmp(:,:,:), ztmp(:,:,:)
    real(rk), allocatable :: h0(:,:), h2(:,:), global(:,:)
    real(rk), allocatable :: tcol(:), z_m(:)
    type(rz_params) :: params

    call reader%global_shape(nx, ny, nz)
    allocate(h0(nx,ny), h2(nx,ny), global(nx,ny), tcol(nz), z_m(nz))
    h0 = 0.0_rk; h2 = 0.0_rk; z_m = z * lengthscale
    call require_shape3('ABL inversion x-pencil theta', shape(theta), reader%gpC%xsz)
    call require_shape3('ABL inversion y-pencil scratch', shape(ytmp), reader%gpC%ysz)
    call require_shape3('ABL inversion z-pencil scratch', shape(ztmp), reader%gpC%zsz)
    call transpose_x_to_y(theta, ytmp, reader%gpC)
    call transpose_y_to_z(ytmp, ztmp, reader%gpC)

    do ic = 1, reader%gpC%zsz(1)
      do jc = 1, reader%gpC%zsz(2)
        ig = reader%gpC%zst(1) + ic - 1
        jg = reader%gpC%zst(2) + jc - 1
        tcol = ztmp(ic,jc,:)
        call fit_rz_profile(z_m, tcol, nz, l0, d0, params, &
          d_min    = 1.0e-6_rk,  &
          ridge    = 1.0e-12_rk, &
          max_iter = 1000,        &
          tol      = 1.0e-10_rk)
        if (params%status == 0) then
          h0(ig,jg) = params%l - xi * params%d
          h2(ig,jg) = params%l + xi * params%d
        end if
      end do
    end do

    global = 0.0_rk
    call MPI_Reduce(h0, global, nx*ny, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) call export_slice_to_netcdf(trim(outdir)//'/'//trim(stem)//'_INVH0.nc', 'INVH0', global, x, y, 'x', 'y')
    global = 0.0_rk
    call MPI_Reduce(h2, global, nx*ny, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) call export_slice_to_netcdf(trim(outdir)//'/'//trim(stem)//'_INVH2.nc', 'INVH2', global, x, y, 'x', 'y')
  end subroutine export_abl_inversion

  pure real(rk) function f_basis(x)
    real(rk), intent(in) :: x
    f_basis = 0.5_rk * (tanh(x) + 1.0_rk)
  end function f_basis

  pure real(rk) function log_2cosh_stable(x)
    real(rk), intent(in) :: x
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

  subroutine run_march(reader, Lx, Ly, Lz, path, outdir, job)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: path, outdir
    type(diag_job_t), intent(in) :: job
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys_dum, ye_dum, zs_dum, ze_dum
    integer :: norm_field_idx, ref_field_idx, ierr, imode, field_idx, icase, ncases
    integer :: istart, iend, islice, idx
    integer, allocatable :: slice_idx(:)
    real(rk) :: dx
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), allocatable :: normalizer(:,:,:), ref(:,:,:), rhs(:,:,:)
    real(rk), allocatable :: recon(:,:,:), scratch(:,:,:)
    real(rk), allocatable :: profile(:), avg(:)
    type(field_case_set_t), allocatable :: case_sets(:)
    type(rms_spec_t) :: range_spec
    character(len=str_len) :: fname, ref_name
    character(len=name_len) :: mode_name
    character(len=32) :: case_suffix

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys_dum, ye_dum, zs_dum, ze_dum)
    if (xs /= 1 .or. xe /= nx) then
      call message('ERROR: march driver requires full x-lines on each rank.')
      call MPI_Abort(MPI_COMM_WORLD, 720, ierr)
    end if

    allocate(x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)
    dx = Lx / real(nx, rk)

    if (job%march%axis /= 'x') then
      call message('ERROR: march driver currently supports axis = x only.')
      call MPI_Abort(MPI_COMM_WORLD, 721, ierr)
    end if

    istart = nearest_grid_index(x, job%march%start)
    iend = nearest_grid_index(x, job%march%finish)
    range_spec = rms_spec_t()
    range_spec%axis = 'x'
    range_spec%bounds_min(1) = min(x(istart), x(iend))
    range_spec%bounds_max(1) = max(x(istart), x(iend))
    range_spec%has_min(1) = .true.
    range_spec%has_max(1) = .true.
    if (myrank == 0) then
      call message('MARCH: requested start='//trim(real_to_tag(job%march%start))// &
        ', using x('//trim(int_to_string(istart))//')='//trim(real_to_tag(x(istart))))
      call message('MARCH: requested end='//trim(real_to_tag(job%march%finish))// &
        ', using x('//trim(int_to_string(iend))//')='//trim(real_to_tag(x(iend))))
      call message('MARCH: scheme='//trim(job%march%scheme))
    end if
    if (job%march%nslices > 0) then
      allocate(slice_idx(job%march%nslices))
      do islice = 1, job%march%nslices
        slice_idx(islice) = nearest_grid_index(x, job%march%slices(islice))
        if (myrank == 0) then
          call message('MARCH: requested slice='//trim(real_to_tag(job%march%slices(islice)))// &
            ', using x('//trim(int_to_string(slice_idx(islice)))//')='// &
            trim(real_to_tag(x(slice_idx(islice)))))
        end if
      end do
    end if

    ref_name = trim(job%march%reference)
    if (len_trim(ref_name) == 0) then
      call message('ERROR: march driver requires reference = <field_name>.')
      call MPI_Abort(MPI_COMM_WORLD, 727, ierr)
    end if
    ref_field_idx = find_field(job, trim(ref_name))
    if (ref_field_idx <= 0) then
      call message('ERROR: march reference field not found: '//trim(ref_name))
      call MPI_Abort(MPI_COMM_WORLD, 722, ierr)
    end if
    norm_field_idx = 0
    if (len_trim(job%march%normalizer) > 0) then
      norm_field_idx = find_field(job, trim(job%march%normalizer))
      if (norm_field_idx <= 0) then
        call message('ERROR: march normalizer field not found: '//trim(job%march%normalizer))
        call MPI_Abort(MPI_COMM_WORLD, 726, ierr)
      end if
    end if
    if (job%march%nmodes <= 0) then
      call message('ERROR: march driver requires at least one mode.')
      call MPI_Abort(MPI_COMM_WORLD, 723, ierr)
    end if

    call prepare_field_case_sets(path, job, case_sets, ncases)
    call validate_march_modes(job, ref_field_idx, norm_field_idx)

    allocate(ref(nxloc,nyloc,nzloc), rhs(nxloc,nyloc,nzloc), recon(nxloc,nyloc,nzloc), &
      scratch(nxloc,nyloc,nzloc), profile(nx))
    if (job%march%export_linear_average) allocate(avg(nx))
    if (norm_field_idx > 0) allocate(normalizer(nxloc,nyloc,nzloc))

    do icase = 1, ncases
      case_suffix = march_case_suffix(icase, ncases)
      call assemble_field_from_case_set(reader, job%fields(ref_field_idx), case_sets(ref_field_idx), icase, ref, scratch)
      if (norm_field_idx > 0) then
        call assemble_field_from_case_set(reader, job%fields(norm_field_idx), case_sets(norm_field_idx), icase, normalizer, scratch)
      end if

      if (job%march%nslices > 0) then
        do islice = 1, job%march%nslices
          idx = slice_idx(islice)
          fname = trim(outdir)//'/'//trim(ref_name)//trim(case_suffix)// &
            '_SL_x='//trim(real_to_tag(x(idx)))//'.nc'
          call export_x_plane_to_netcdf(reader, trim(fname), trim(ref_name), ref, idx, y, z)
        end do
      end if

      do imode = 1, job%march%nmodes
        mode_name = trim(job%march%modes(imode)%name)
        rhs = 0.0_rk
        do field_idx = 1, job%nfields
          if (field_idx == ref_field_idx .or. field_idx == norm_field_idx) cycle
          if (mode_removes_field(job%march%modes(imode), trim(job%fields(field_idx)%name))) cycle
          call assemble_field_from_case_set(reader, job%fields(field_idx), case_sets(field_idx), icase, recon, scratch)
          rhs = rhs + recon
        end do

        if (norm_field_idx > 0) then
          rhs = rhs / normalizer
        end if
        call integrate_march_x(ref, rhs, recon, dx, istart, iend, trim(job%march%scheme))
        ! Reuse scratch as the 3D error field after file assembly is complete.
        scratch = recon - ref

        if (job%march%export_rms) then
          call l2_profile_bounded(reader, x, y, z, Lx, Ly, Lz, recon, range_spec, profile)
          fname = trim(outdir)//'/'//trim(mode_name)//trim(case_suffix)//'_'//trim(ref_name)//'_rms_x.csv'
          if (myrank == 0) then
            call csv_profile_bounded(nx, trim(fname), 'x', x, profile, range_spec, 1, 'l2_norm')
          end if
          call l2_profile_bounded(reader, x, y, z, Lx, Ly, Lz, scratch, range_spec, profile)
          fname = trim(outdir)//'/'//trim(mode_name)//trim(case_suffix)//'_error_rms_x.csv'
          if (myrank == 0) then
            call csv_profile_bounded(nx, trim(fname), 'x', x, profile, range_spec, 1, 'l2_norm')
          end if
        end if
        if (job%march%export_linear_average) then
          call compute_bounded_linear_average(reader, x, y, z, Lx, Ly, Lz, recon, range_spec, avg)
          fname = trim(outdir)//'/'//trim(mode_name)//trim(case_suffix)//'_'//trim(ref_name)//'_avg_x.csv'
          if (myrank == 0) then
            call csv_profile_bounded(nx, trim(fname), 'x', x, avg, range_spec, 1, 'average')
          end if
          call compute_bounded_linear_average(reader, x, y, z, Lx, Ly, Lz, scratch, range_spec, avg)
          fname = trim(outdir)//'/'//trim(mode_name)//trim(case_suffix)//'_error_avg_x.csv'
          if (myrank == 0) then
            call csv_profile_bounded(nx, trim(fname), 'x', x, avg, range_spec, 1, 'average')
          end if
        end if

        if (job%march%nslices > 0) then
          do islice = 1, job%march%nslices
            idx = slice_idx(islice)
            fname = trim(outdir)//'/'//trim(mode_name)//trim(case_suffix)//'_recon_SL_x='// &
              trim(real_to_tag(x(idx)))//'.nc'
            call export_x_plane_to_netcdf(reader, trim(fname), trim(mode_name)//'_'//trim(ref_name), recon, idx, y, z)
            fname = trim(outdir)//'/'//trim(mode_name)//trim(case_suffix)//'_error_SL_x='// &
              trim(real_to_tag(x(idx)))//'.nc'
            call export_x_plane_to_netcdf(reader, trim(fname), trim(mode_name)//'_error', scratch, idx, y, z)
          end do
        end if
      end do

      if (job%march%export_rms) then
        call l2_profile_bounded(reader, x, y, z, Lx, Ly, Lz, ref, range_spec, profile)
        fname = trim(outdir)//'/'//trim(ref_name)//trim(case_suffix)//'_rms_x.csv'
        if (myrank == 0) then
          call csv_profile_bounded(nx, trim(fname), 'x', x, profile, range_spec, 1, 'l2_norm')
        end if
      end if
      if (job%march%export_linear_average) then
        call compute_bounded_linear_average(reader, x, y, z, Lx, Ly, Lz, ref, range_spec, avg)
        fname = trim(outdir)//'/'//trim(ref_name)//trim(case_suffix)//'_avg_x.csv'
        if (myrank == 0) then
          call csv_profile_bounded(nx, trim(fname), 'x', x, avg, range_spec, 1, 'average')
        end if
      end if
    end do
    if (allocated(avg)) deallocate(avg)
  end subroutine run_march

  function march_case_suffix(icase, ncases) result(suffix)
    integer, intent(in) :: icase, ncases
    character(len=:), allocatable :: suffix
    if (ncases > 1) then
      suffix = '_case'//trim(int_to_string(icase))
    else
      suffix = ''
    end if
  end function march_case_suffix

  integer function find_field(job, name) result(idx)
    type(diag_job_t), intent(in) :: job
    character(*), intent(in) :: name
    integer :: i
    idx = 0
    do i = 1, job%nfields
      if (trim(lower(job%fields(i)%name)) == trim(lower(name))) then
        idx = i
        return
      end if
    end do
  end function find_field

  logical function mode_removes_field(mode, field_name) result(removes)
    type(march_mode_t), intent(in) :: mode
    character(*), intent(in) :: field_name
    integer :: i
    removes = .false.
    do i = 1, mode%nremove
      if (trim(lower(mode%remove(i))) == trim(lower(field_name))) then
        removes = .true.
        return
      end if
    end do
  end function mode_removes_field

  subroutine validate_march_modes(job, ref_field_idx, norm_field_idx)
    type(diag_job_t), intent(in) :: job
    integer, intent(in) :: ref_field_idx, norm_field_idx
    integer :: imode, iremove, remove_idx
    do imode = 1, job%march%nmodes
      do iremove = 1, job%march%modes(imode)%nremove
        remove_idx = find_field(job, trim(job%march%modes(imode)%remove(iremove)))
        if (remove_idx <= 0) then
          call message('WARNING: march mode '//trim(job%march%modes(imode)%name)// &
            ' removes unknown field '//trim(job%march%modes(imode)%remove(iremove))//'; ignoring.')
        else if (remove_idx == ref_field_idx .or. remove_idx == norm_field_idx) then
          call message('WARNING: march mode '//trim(job%march%modes(imode)%name)// &
            ' removes reference/normalizer field '//trim(job%march%modes(imode)%remove(iremove))// &
            '; this has no effect.')
        end if
      end do
    end do
  end subroutine validate_march_modes

  subroutine prepare_field_case_sets(path, job, case_sets, ncases_total)
    character(*), intent(in) :: path
    type(diag_job_t), intent(in) :: job
    type(field_case_set_t), allocatable, intent(out) :: case_sets(:)
    integer, intent(out) :: ncases_total
    integer :: ifield, ierr
    ncases_total = 1
    if (allocated(case_sets)) deallocate(case_sets)
    allocate(case_sets(job%nfields))
    do ifield = 1, job%nfields
      if (job%fields(ifield)%nterms <= 0) then
        call message('ERROR: march driver does not support derived fields (field: '// &
          trim(job%fields(ifield)%name)//').')
        call MPI_Abort(MPI_COMM_WORLD, 725, ierr)
      end if
      call expand_terms(path, job%fields(ifield)%terms, job%fields(ifield)%nterms, &
        case_sets(ifield)%lists, case_sets(ifield)%ncases, ierr)
      if (ierr /= 0) call abort_expand_error(ierr)
      ncases_total = max(ncases_total, case_sets(ifield)%ncases)
    end do
    do ifield = 1, job%nfields
      if (case_sets(ifield)%ncases /= 1 .and. case_sets(ifield)%ncases /= ncases_total) then
        call message('ERROR: march wildcard fields must have one case or the same number of cases.')
        call MPI_Abort(MPI_COMM_WORLD, 724, ierr)
      end if
    end do
  end subroutine prepare_field_case_sets

  subroutine assemble_field_from_case_set(reader, field, case_set, icase, f, scratch)
    class(FieldReader2Decomp), intent(inout) :: reader
    type(field_expr_t), intent(in) :: field
    type(field_case_set_t), intent(in) :: case_set
    integer, intent(in) :: icase
    real(rk), intent(out) :: f(:,:,:)
    real(rk), intent(inout) :: scratch(:,:,:)
    call assemble_terms_case(reader, field%terms, field%nterms, case_set%lists, icase, f, scratch)
  end subroutine assemble_field_from_case_set

  integer function nearest_grid_index(coord, value) result(idx)
    real(rk), intent(in) :: coord(:), value
    integer :: i
    real(rk) :: best, dist
    idx = 1
    best = abs(coord(1) - value)
    do i = 2, size(coord)
      dist = abs(coord(i) - value)
      if (dist < best) then
        best = dist
        idx = i
      end if
    end do
  end function nearest_grid_index

  subroutine integrate_march_x(ref, dfdx, recon, dx, istart, iend, scheme)
    real(rk), intent(in) :: ref(:,:,:), dfdx(:,:,:), dx
    real(rk), intent(out) :: recon(:,:,:)
    integer, intent(in) :: istart, iend
    character(*), intent(in) :: scheme
    integer :: i, idir
    recon = ref
    idir = merge(1, -1, iend >= istart)
    if (istart == iend) return
    do i = istart, iend - idir, idir
      if (trim(scheme) == 'trapezoid') then
        recon(i+idir,:,:) = recon(i,:,:) + real(idir, rk) * dx * &
          0.5_rk * (dfdx(i,:,:) + dfdx(i+idir,:,:))
      else
        recon(i+idir,:,:) = recon(i,:,:) + real(idir, rk) * dx * dfdx(i,:,:)
      end if
    end do
  end subroutine integrate_march_x

  subroutine subtract_crossplane_reference(reader, f, item, ref)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(inout) :: f(:,:,:)
    type(rms_spec_t), intent(in) :: item
    real(rk), intent(in) :: ref(:)
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze
    integer :: i, j, k, ig, jg, kg, ref_idx
    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)

    select case(item%axis)
    case('x')
      do i = 1, nxloc
        do j = 1, nyloc
          jg = ys + j - 1
          do k = 1, nzloc
            kg = zs + k - 1
            ref_idx = (kg - 1) * ny + jg
            f(i,j,k) = f(i,j,k) - ref(ref_idx)
          end do
        end do
      end do
    case('y')
      do j = 1, nyloc
        do i = 1, nxloc
          ig = xs + i - 1
          do k = 1, nzloc
            kg = zs + k - 1
            ref_idx = (kg - 1) * nx + ig
            f(i,j,k) = f(i,j,k) - ref(ref_idx)
          end do
        end do
      end do
    case default
      do k = 1, nzloc
        do i = 1, nxloc
          ig = xs + i - 1
          do j = 1, nyloc
            jg = ys + j - 1
            ref_idx = (jg - 1) * nx + ig
            f(i,j,k) = f(i,j,k) - ref(ref_idx)
          end do
        end do
      end do
    end select
  end subroutine subtract_crossplane_reference

  subroutine l2_profile_bounded(reader, x, y, z, Lx, Ly, Lz, f, item, profile)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: x(:), y(:), z(:), Lx, Ly, Lz, f(:,:,:)
    type(rms_spec_t), intent(in) :: item
    real(rk), intent(out) :: profile(:)
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze
    integer :: i, j, k, ig, jg, kg, ni, ierr
    real(rk) :: dx, dy, dz, area_weight
    real(rk), allocatable :: local_sum(:), global_sum(:)

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    dx = Lx / real(nx, rk)
    dy = Ly / real(ny, rk)
    dz = Lz / real(nz, rk)

    select case(item%axis)
    case('x')
      ni = nx
      area_weight = dy * dz
    case('y')
      ni = ny
      area_weight = dx * dz
    case default
      ni = nz
      area_weight = dx * dy
    end select

    allocate(local_sum(ni), global_sum(ni))
    local_sum = 0.0_rk

    select case(item%axis)
    case('x')
      do i = 1, nxloc
        ig = xs + i - 1
        if (.not. coord_in_bounds(x(ig), item, 1)) cycle
        do j = 1, nyloc
          jg = ys + j - 1
          if (.not. coord_in_bounds(y(jg), item, 2)) cycle
          do k = 1, nzloc
            kg = zs + k - 1
            if (.not. coord_in_bounds(z(kg), item, 3)) cycle
            local_sum(ig) = local_sum(ig) + f(i,j,k)**2 * area_weight
          end do
        end do
      end do
    case('y')
      do j = 1, nyloc
        jg = ys + j - 1
        if (.not. coord_in_bounds(y(jg), item, 2)) cycle
        do i = 1, nxloc
          ig = xs + i - 1
          if (.not. coord_in_bounds(x(ig), item, 1)) cycle
          do k = 1, nzloc
            kg = zs + k - 1
            if (.not. coord_in_bounds(z(kg), item, 3)) cycle
            local_sum(jg) = local_sum(jg) + f(i,j,k)**2 * area_weight
          end do
        end do
      end do
    case default
      do k = 1, nzloc
        kg = zs + k - 1
        if (.not. coord_in_bounds(z(kg), item, 3)) cycle
        do i = 1, nxloc
          ig = xs + i - 1
          if (.not. coord_in_bounds(x(ig), item, 1)) cycle
          do j = 1, nyloc
            jg = ys + j - 1
            if (.not. coord_in_bounds(y(jg), item, 2)) cycle
            local_sum(kg) = local_sum(kg) + f(i,j,k)**2 * area_weight
          end do
        end do
      end do
    end select

    call MPI_Allreduce(local_sum, global_sum, ni, mpi_rk, MPI_SUM, MPI_COMM_WORLD, ierr)
    profile = sqrt(global_sum)
    deallocate(local_sum, global_sum)
  end subroutine l2_profile_bounded

  subroutine compute_bounded_linear_average(reader, x, y, z, Lx, Ly, Lz, f, item, profile)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: x(:), y(:), z(:), Lx, Ly, Lz, f(:,:,:)
    type(rms_spec_t), intent(in) :: item
    real(rk), intent(out) :: profile(:)
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze
    integer :: i, j, k, ig, jg, kg, ni, ierr
    real(rk) :: dx, dy, dz, area_weight
    real(rk), allocatable :: local_sum(:), global_sum(:), local_area(:), global_area(:)

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    dx = Lx / real(nx, rk)
    dy = Ly / real(ny, rk)
    dz = Lz / real(nz, rk)

    select case(item%axis)
    case('x')
      ni = nx
      area_weight = dy * dz
    case('y')
      ni = ny
      area_weight = dx * dz
    case default
      ni = nz
      area_weight = dx * dy
    end select

    allocate(local_sum(ni), global_sum(ni), local_area(ni), global_area(ni))
    local_sum = 0.0_rk
    local_area = 0.0_rk

    select case(item%axis)
    case('x')
      do i = 1, nxloc
        ig = xs + i - 1
        if (.not. coord_in_bounds(x(ig), item, 1)) cycle
        do j = 1, nyloc
          jg = ys + j - 1
          if (.not. coord_in_bounds(y(jg), item, 2)) cycle
          do k = 1, nzloc
            kg = zs + k - 1
            if (.not. coord_in_bounds(z(kg), item, 3)) cycle
            local_sum(ig) = local_sum(ig) + f(i,j,k) * area_weight
            local_area(ig) = local_area(ig) + area_weight
          end do
        end do
      end do
    case('y')
      do j = 1, nyloc
        jg = ys + j - 1
        if (.not. coord_in_bounds(y(jg), item, 2)) cycle
        do i = 1, nxloc
          ig = xs + i - 1
          if (.not. coord_in_bounds(x(ig), item, 1)) cycle
          do k = 1, nzloc
            kg = zs + k - 1
            if (.not. coord_in_bounds(z(kg), item, 3)) cycle
            local_sum(jg) = local_sum(jg) + f(i,j,k) * area_weight
            local_area(jg) = local_area(jg) + area_weight
          end do
        end do
      end do
    case default
      do k = 1, nzloc
        kg = zs + k - 1
        if (.not. coord_in_bounds(z(kg), item, 3)) cycle
        do i = 1, nxloc
          ig = xs + i - 1
          if (.not. coord_in_bounds(x(ig), item, 1)) cycle
          do j = 1, nyloc
            jg = ys + j - 1
            if (.not. coord_in_bounds(y(jg), item, 2)) cycle
            local_sum(kg) = local_sum(kg) + f(i,j,k) * area_weight
            local_area(kg) = local_area(kg) + area_weight
          end do
        end do
      end do
    end select

    call MPI_Allreduce(local_sum, global_sum, ni, mpi_rk, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_area, global_area, ni, mpi_rk, MPI_SUM, MPI_COMM_WORLD, ierr)
    profile = 0.0_rk
    do i = 1, ni
      if (global_area(i) > 0.0_rk) profile(i) = global_sum(i) / global_area(i)
    end do
    deallocate(local_sum, global_sum, local_area, global_area)
  end subroutine compute_bounded_linear_average

  subroutine export_x_plane_to_netcdf(reader, fname, varname, f, idx, y, z)
    class(FieldReader2Decomp), intent(inout) :: reader
    character(*), intent(in) :: fname, varname
    real(rk), intent(in) :: f(:,:,:), y(:), z(:)
    integer, intent(in) :: idx
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze, ierr
    real(rk), allocatable :: local(:,:), global(:,:)
    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    allocate(local(ny,nz))
    if (myrank == 0) then
      allocate(global(ny,nz))
    else
      allocate(global(1,1))
    end if
    local = 0.0_rk
    if (idx >= xs .and. idx <= xe) then
      local(ys:ye,zs:ze) = f(idx-xs+1,:,:)
    end if
    call MPI_Reduce(local, global, ny*nz, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) call export_slice_to_netcdf(trim(fname), trim(varname), global, y, z, 'y', 'z')
    deallocate(local, global)
  end subroutine export_x_plane_to_netcdf

  subroutine fit_rz_profile(z, t, n, l0, d0, params, d_min, ridge, max_iter, tol)
    integer, intent(in) :: n
    real(rk), intent(in) :: z(n), t(n), l0, d0
    type(rz_params), intent(out) :: params
    real(rk), intent(in), optional :: d_min, ridge, tol
    integer, intent(in), optional :: max_iter

    real(rk) :: dmin_loc, ridge_loc, tol_loc
    integer :: max_iter_loc
    real(rk) :: x(2,3), fval(3)
    real(rk) :: centroid(2), xr(2), xe(2), xc(2)
    real(rk) :: fr, fe, fc
    real(rk) :: alpha, gamma, rho, sigma
    real(rk) :: scale_l, scale_d, z_range, simplex_size
    integer :: iter, i_best, i_worst, i_mid, stat

    dmin_loc = 1.0e-6_rk
    ridge_loc = 0.0_rk
    max_iter_loc = 500
    tol_loc = 1.0e-8_rk
    if (present(d_min)) dmin_loc = d_min
    if (present(ridge)) ridge_loc = ridge
    if (present(max_iter)) max_iter_loc = max_iter
    if (present(tol)) tol_loc = tol

    alpha = 1.0_rk
    gamma = 2.0_rk
    rho = 0.5_rk
    sigma = 0.5_rk

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
      xr = centroid + alpha * (centroid - x(:,i_worst))
      xr(2) = max(xr(2), dmin_loc)
      fr = objective(z, t, n, xr(1), xr(2), dmin_loc, ridge_loc, stat)

      if (fr < fval(i_best)) then
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
          x(:,i_mid) = x(:,i_best) + sigma * (x(:,i_mid) - x(:,i_best))
          x(:,i_worst) = x(:,i_best) + sigma * (x(:,i_worst) - x(:,i_best))
          x(2,i_mid) = max(x(2,i_mid), dmin_loc)
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
    integer, intent(in) :: n
    real(rk), intent(in) :: z(n), t(n), l, d, d_min, ridge
    integer, intent(out) :: status
    real(rk) :: tm, a, b, eta, t_fit, res
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
    integer, intent(in) :: n
    real(rk), intent(in) :: z(n), t(n), l, d, ridge
    real(rk), intent(out) :: tm, a, b
    integer, intent(out) :: status
    real(rk) :: M(3,3), rhs(3), sol(3)
    real(rk) :: eta, fv, gv
    real(rk) :: S1, Sf, Sg, Sff, Sgg, Sfg, St, Sft, Sgt
    integer :: i

    S1 = real(n, rk)
    Sf = 0.0_rk; Sg = 0.0_rk
    Sff = 0.0_rk; Sgg = 0.0_rk; Sfg = 0.0_rk
    St = 0.0_rk; Sft = 0.0_rk; Sgt = 0.0_rk
    do i = 1, n
      eta = (z(i) - l) / d
      fv = f_basis(eta)
      gv = g_basis(eta)
      Sf = Sf + fv
      Sg = Sg + gv
      Sff = Sff + fv * fv
      Sgg = Sgg + gv * gv
      Sfg = Sfg + fv * gv
      St = St + t(i)
      Sft = Sft + fv * t(i)
      Sgt = Sgt + gv * t(i)
    end do

    M(1,:) = [S1, Sf, Sg]
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
      a = sol(2)
      b = sol(3)
    else
      tm = 0.0_rk
      a = 0.0_rk
      b = 0.0_rk
    end if
  end subroutine solve_tmab

  subroutine solve_3x3(Ain, bin, x, status)
    real(rk), intent(in)  :: Ain(3,3), bin(3)
    real(rk), intent(out) :: x(3)
    integer, intent(out)  :: status
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
    real(rk), intent(in) :: fval(3)
    integer, intent(out) :: i_best, i_mid, i_worst
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
    i_best = idx(1)
    i_mid = idx(2)
    i_worst = idx(3)
  end subroutine order_simplex

  subroutine find_threshold_crossing_z(field, z, zcross, threshold, missing_value)
    real(rk), intent(in) :: field(:,:,:), z(:), threshold, missing_value
    real(rk), intent(out) :: zcross(size(field,1), size(field,2))
    integer :: i, j, k, nx, ny, nz
    real(rk) :: fref, fthr, alpha
    nx = size(field,1); ny = size(field,2); nz = size(field,3)
    zcross = missing_value
    do j = 1, ny
      do i = 1, nx
        fref = field(i,j,1)
        fthr = threshold * fref
        if (field(i,j,1) <= fthr) then
          zcross(i,j) = z(1)
          cycle
        end if
        do k = 2, nz
          if (field(i,j,k) <= fthr) then
            if (abs(field(i,j,k)-field(i,j,k-1)) > tiny(1.0_rk)) then
              alpha = (fthr - field(i,j,k-1)) / (field(i,j,k)-field(i,j,k-1))
              zcross(i,j) = z(k-1) + alpha * (z(k)-z(k-1))
            else
              zcross(i,j) = z(k-1)
            end if
            exit
          end if
        end do
      end do
    end do
  end subroutine find_threshold_crossing_z

  subroutine csv_profile(n, filename, axis, coord, profile)
    integer, intent(in) :: n
    character(*), intent(in) :: filename, axis
    real(rk), intent(in) :: coord(n), profile(n)
    integer :: u, i
    call message('Writing file: '//trim(filename))
    open(newunit=u, file=trim(filename), status='replace', action='write')
    write(u,'(A)') trim(axis)//',profile'
    do i = 1, n
      write(u,'(ES23.15,",",ES23.15)') coord(i), profile(i)
    end do
    close(u)
  end subroutine csv_profile

  subroutine csv_profile_bounded(n, filename, axis, coord, profile, item, iax, header)
    integer, intent(in) :: n, iax
    character(*), intent(in) :: filename, axis, header
    real(rk), intent(in) :: coord(n), profile(n)
    type(rms_spec_t), intent(in) :: item
    integer :: u, i
    call message('Writing file: '//trim(filename))
    open(newunit=u, file=trim(filename), status='replace', action='write')
    write(u,'(A)') trim(axis)//','//trim(header)
    do i = 1, n
      if (.not. coord_in_bounds(coord(i), item, iax)) cycle
      write(u,'(ES23.15,",",ES23.15)') coord(i), profile(i)
    end do
    close(u)
  end subroutine csv_profile_bounded

  subroutine export_slice_to_netcdf(fname, varname, slice, x1, x2, x1_name, x2_name)
    character(*), intent(in) :: fname, varname, x1_name, x2_name
    real(rk), intent(in) :: slice(:,:), x1(:), x2(:)
    integer :: ncid, dimid_x1, dimid_x2, varid_x1, varid_x2, varid_f, ierr
    integer :: dimids_f(2)
    call message('Writing file: '//trim(fname))
    ierr = nf90_create(trim(fname), ior(NF90_CLOBBER, NF90_NETCDF4), ncid)
    call nc_check(ierr, 'nf90_create')
    ierr = nf90_def_dim(ncid, trim(x1_name), size(x1), dimid_x1); call nc_check(ierr, 'def_dim x1')
    ierr = nf90_def_dim(ncid, trim(x2_name), size(x2), dimid_x2); call nc_check(ierr, 'def_dim x2')
    ierr = nf90_def_var(ncid, trim(x1_name), nc_rk, dimid_x1, varid_x1); call nc_check(ierr, 'def_var x1')
    ierr = nf90_def_var(ncid, trim(x2_name), nc_rk, dimid_x2, varid_x2); call nc_check(ierr, 'def_var x2')
    dimids_f = [dimid_x1, dimid_x2]
    ierr = nf90_def_var(ncid, trim(varname), nc_rk, dimids_f, varid_f); call nc_check(ierr, 'def_var field')
    ierr = nf90_enddef(ncid); call nc_check(ierr, 'enddef')
    ierr = nf90_put_var(ncid, varid_x1, x1); call nc_check(ierr, 'put x1')
    ierr = nf90_put_var(ncid, varid_x2, x2); call nc_check(ierr, 'put x2')
    ierr = nf90_put_var(ncid, varid_f, slice); call nc_check(ierr, 'put field')
    ierr = nf90_close(ncid); call nc_check(ierr, 'close')
  contains
    subroutine nc_check(status, where)
      integer, intent(in) :: status
      character(*), intent(in) :: where
      integer :: mpi_ierr
      if (status /= nf90_noerr) then
        if (myrank == 0) then
          write(*,*) 'NetCDF error in ', trim(where), ': ', trim(nf90_strerror(status))
        end if
        call MPI_Abort(MPI_COMM_WORLD, 900, mpi_ierr)
      end if
    end subroutine nc_check
  end subroutine export_slice_to_netcdf

  subroutine find_bracket_uniform(coord, z0, k0, k1, alpha)
    real(rk), intent(in) :: coord(:), z0
    integer, intent(out) :: k0, k1
    real(rk), intent(out) :: alpha
    integer :: n
    real(rk) :: denom, s

    n = size(coord)
    if (n <= 1) then
      k0 = 1
      k1 = 1
      alpha = 0.0_rk
      return
    end if

    ! Negative slice coordinates are retained as an index shortcut:
    ! -1 selects index 1, -2 selects index 2, etc.  Use physical
    ! coordinates for domains where negative positions are meaningful.
    if (z0 <= -1.0_rk) then
      k0 = max(1, min(n, int(abs(z0))))
      k1 = k0; alpha = 0.0_rk
      return
    end if

    if (z0 <= coord(1)) then
      k0 = 1; k1 = 1; alpha = 0.0_rk
      return
    end if
    if (z0 >= coord(n)) then
      k0 = n; k1 = n; alpha = 0.0_rk
      return
    end if

    ! O(1) bracket lookup assumes coord(:) is uniformly spaced.
    denom = coord(2) - coord(1)
    if (abs(denom) <= tiny(1.0_rk)) then
      k0 = 1; k1 = 1; alpha = 0.0_rk
      return
    end if

    s = (z0 - coord(1)) / denom
    k0 = max(1, min(n - 1, int(floor(s)) + 1))
    k1 = k0 + 1
    denom = coord(k1) - coord(k0)
    if (abs(denom) > tiny(1.0_rk)) then
      alpha = (z0 - coord(k0)) / denom
    else
      alpha = 0.0_rk
    end if
  end subroutine find_bracket_uniform

  function real_to_tag(z) result(tag)
    real(rk), intent(in) :: z
    character(len=64) :: tmp
    character(len=:), allocatable :: tag
    integer :: n, i
    write(tmp,'(F20.8)') z
    tmp = adjustl(tmp)
    n = len_trim(tmp)
    do i = 1, n
      if (tmp(i:i) == '.') tmp(i:i) = 'p'
    end do
    do while (n > 1 .and. tmp(n:n) == '0')
      n = n - 1
    end do
    if (tmp(n:n) == 'p') n = n - 1
    tag = tmp(:n)
  end function real_to_tag

  function int_to_string(i) result(s)
    integer, intent(in) :: i
    character(len=:), allocatable :: s
    character(len=32) :: tmp
    write(tmp,'(I0)') i
    s = trim(tmp)
  end function int_to_string

  subroutine require_shape3(label, actual, expected)
    character(*), intent(in) :: label
    integer, intent(in) :: actual(3), expected(3)
    integer :: ierr
    if (all(actual == expected)) return
    call message('ERROR: '//trim(label)//' shape is '//int_to_string(actual(1))//' x '// &
      int_to_string(actual(2))//' x '//int_to_string(actual(3))//'; expected '// &
      int_to_string(expected(1))//' x '//int_to_string(expected(2))//' x '//int_to_string(expected(3)))
    call MPI_Abort(MPI_COMM_WORLD, 910, ierr)
  end subroutine require_shape3

  subroutine frd_init(this, nx, ny, nz)
    class(FieldReader2Decomp), intent(inout) :: this
    integer, intent(in) :: nx, ny, nz
    if (this%is_init) return
    this%nx = nx; this%ny = ny; this%nz = nz
    call this%choose_proc_grid_(nprocs, nx, ny, this%p_row, this%p_col)
    if (myrank == 0) then
      call message('FRD_INIT: proc grid = '//int_to_string(this%p_row)//' x '//int_to_string(this%p_col))
    end if
    call decomp_2d_init(nx, ny, nz, this%p_row, this%p_col)
    call decomp_info_init(nx, ny, nz, this%gpC)
    this%xs=this%gpC%xst(1); this%xe=this%gpC%xen(1)
    this%ys=this%gpC%xst(2); this%ye=this%gpC%xen(2)
    this%zs=this%gpC%xst(3); this%ze=this%gpC%xen(3)
    this%nxloc=this%gpC%xsz(1)
    this%nyloc=this%gpC%xsz(2)
    this%nzloc=this%gpC%xsz(3)
    this%is_init = .true.
  end subroutine frd_init

  subroutine frd_read_field(this, path, field)
    class(FieldReader2Decomp), intent(inout) :: this
    character(*), intent(in) :: path
    real(rk), intent(out) :: field(this%nxloc,this%nyloc,this%nzloc)
    integer :: ierr
    character(len=str_len) :: read_path
    logical :: exists
    read_path = trim(path)
    if (.not. this%is_init) then
      call message('ERROR: FieldReader2Decomp is not initialized before read_field().')
      call MPI_Abort(MPI_COMM_WORLD, 101, ierr)
    end if
    inquire(file=trim(read_path), exist=exists)
    if (.not. exists) then
      call message('ERROR: input field file does not exist: '//trim(read_path))
      call MPI_Abort(MPI_COMM_WORLD, 102, ierr)
    end if
    call decomp_2d_read_one(1, field, read_path, this%gpC)
  end subroutine frd_read_field

  subroutine frd_global_shape(this, nx, ny, nz)
    class(FieldReader2Decomp), intent(in) :: this
    integer, intent(out) :: nx, ny, nz
    nx=this%nx; ny=this%ny; nz=this%nz
  end subroutine frd_global_shape

  subroutine frd_local_shape(this, nxloc, nyloc, nzloc)
    class(FieldReader2Decomp), intent(in) :: this
    integer, intent(out) :: nxloc, nyloc, nzloc
    nxloc=this%nxloc; nyloc=this%nyloc; nzloc=this%nzloc
  end subroutine frd_local_shape

  subroutine frd_index(this, xs, xe, ys, ye, zs, ze)
    class(FieldReader2Decomp), intent(in) :: this
    integer, intent(out) :: xs, xe, ys, ye, zs, ze
    xs=this%xs; xe=this%xe; ys=this%ys; ye=this%ye; zs=this%zs; ze=this%ze
  end subroutine frd_index

  subroutine frd_finalize(this)
    type(FieldReader2Decomp), intent(inout) :: this
    if (this%is_init) then
      call decomp_2d_finalize
      this%is_init = .false.
    end if
  end subroutine frd_finalize

  subroutine frd_choose_proc_grid(nproc, nx, ny, p_row, p_col)
    integer, intent(in) :: nproc, nx, ny
    integer, intent(out) :: p_row, p_col
    integer :: r, best_row, best_col
    real(rk) :: best_score, score

    best_score = huge(1.0_rk)
    best_row = 1
    best_col = nproc
    do r = 1, nproc
      if (mod(nproc, r) /= 0) cycle
      score = abs(real(nx, rk) / real(max(1, ny), rk) - real(nproc / r, rk) / real(r, rk))
      if (score < best_score) then
        best_score = score
        best_row = r
        best_col = nproc / r
      end if
    end do
    p_row = best_row
    p_col = best_col
  end subroutine frd_choose_proc_grid

end module MPIR3D_Lean

program MPIR3D_Lean_Main
  use MPIR3D_Lean
  implicit none

  type(FieldReader2Decomp) :: reader
  type(diag_job_t) :: job
  integer :: ierr, nlen, ioUnit
  character(:), allocatable :: inputfile
  character(len=256) :: path='.', outdir='.', diag_map=''
  character(len=32) :: driver=''
  integer :: nx=1, ny=1, nz=1
  real(rk) :: Lx=1.0_rk, Ly=1.0_rk, Lz=1.0_rk
  namelist /SETUP/ nx, ny, nz, Lx, Ly, Lz, path, outdir, driver, diag_map

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  if (command_argument_count() < 1) then
    call message('Usage: MPIR3D_lean <inputfile>')
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  call get_command_argument(1, length=nlen)
  allocate(character(nlen) :: inputfile)
  call get_command_argument(1, value=inputfile)

  open(newunit=ioUnit, file=trim(inputfile), status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    call message('ERROR: could not open input file.')
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
  end if
  read(ioUnit, nml=SETUP, iostat=ierr)
  close(ioUnit)
  if (ierr /= 0) then
    call message('ERROR: could not read SETUP namelist.')
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
  end if

  call read_diag_map(trim(diag_map), job, ierr)
  if (ierr /= 0) then
    call message('ERROR: could not parse diagnostic map; ierr='//trim(int_to_string(ierr)))
    call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
  end if

  if (len_trim(driver) > 0) job%driver = trim(lower(driver))
  call reader%init(nx, ny, nz)

  select case(trim(lower(job%driver)))
  case('ha', 'horizontal_average')
    call run_horizontal_average(reader, Lx, Ly, Lz, trim(path), trim(outdir), job)
  case('slice', 'slices')
    call run_slices(reader, Lx, Ly, Lz, trim(path), trim(outdir), job)
  case('rms')
    call run_rms(reader, Lx, Ly, Lz, trim(path), trim(outdir), job)
  case('profile', 'profiles', 'linear_profile')
    call run_linear_profiles(reader, Lx, Ly, Lz, trim(path), trim(outdir), job)
  case('abl')
    call run_abl(reader, Lx, Ly, Lz, trim(path), trim(outdir), job)
  case('march')
    call run_march(reader, Lx, Ly, Lz, trim(path), trim(outdir), job)
  case default
    call message('ERROR: unknown driver '//trim(job%driver))
    call MPI_Abort(MPI_COMM_WORLD, 5, ierr)
  end select

  call message('Wrapping up ...')
  call MPI_Finalize(ierr)
end program MPIR3D_Lean_Main
