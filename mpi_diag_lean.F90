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
    type(abl_spec_t) :: abl
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

  subroutine read_diag_map(filename, job, ierr)
    character(*), intent(in) :: filename
    type(diag_job_t), intent(out) :: job
    integer, intent(out) :: ierr
    integer :: unit, ios
    character(len=str_len) :: line, clean, low

    ierr = 0
    job%nfields = 0; job%nslices = 0; job%nrms = 0
    if (allocated(job%fields)) deallocate(job%fields)
    if (allocated(job%slices)) deallocate(job%slices)
    if (allocated(job%rms)) deallocate(job%rms)
    allocate(job%fields(0), job%slices(0), job%rms(0))

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
      else if (key_is(low, 'fields')) then
        call read_fields_block(unit, job, ierr)
      else if (key_is(low, 'slices')) then
        call read_slice_block(unit, job, ierr)
      else if (key_is(low, 'rms')) then
        call read_rms_block(unit, job, ierr)
      else if (key_is(low, 'abl')) then
        call read_abl_block(unit, job, ierr)
      end if
      if (ierr /= 0) exit
    end do
    close(unit)

    if (ios > 0 .and. ierr == 0) ierr = ios
  end subroutine read_diag_map

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
    call append_rms(job%rms, job%nrms, item)
  end subroutine read_rms_block

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

  pure logical function coord_in_bounds(coord, item, iax)
    real(rk), intent(in) :: coord
    type(rms_spec_t), intent(in) :: item
    integer, intent(in) :: iax
    coord_in_bounds = coord >= item%bounds_min(iax) .and. coord <= item%bounds_max(iax)
  end function coord_in_bounds

  subroutine assemble_field(reader, root, expr, f)
    class(FieldReader2Decomp), intent(inout) :: reader
    character(*), intent(in) :: root
    type(field_expr_t), intent(in) :: expr
    real(rk), intent(out) :: f(:,:,:)
    real(rk), allocatable :: tmp(:,:,:)
    character(len=:), allocatable :: infile
    integer :: iterm, ierr

    if (expr%nterms <= 0) then
      call message('ERROR: field expression has no direct terms. Derived fields are only supported by selected drivers.')
      call MPI_Abort(MPI_COMM_WORLD, 810, ierr)
    end if
    f = 0.0_rk
    do iterm = 1, expr%nterms
      infile = resolve_input_path(trim(root), trim(expr%terms(iterm)%filename))
      call message('Reading '//trim(infile))
      tmp = reader%read_field(trim(infile))
      f = f + expr%terms(iterm)%coeff * tmp
    end do
  end subroutine assemble_field

  subroutine assemble_terms(reader, root, terms, nterms, f)
    class(FieldReader2Decomp), intent(inout) :: reader
    character(*), intent(in) :: root
    type(field_term_t), intent(in) :: terms(:)
    integer, intent(in) :: nterms
    real(rk), intent(out) :: f(:,:,:)
    real(rk), allocatable :: tmp(:,:,:)
    character(len=:), allocatable :: infile
    integer :: iterm

    f = 0.0_rk
    do iterm = 1, nterms
      infile = resolve_input_path(trim(root), trim(terms(iterm)%filename))
      call message('Reading '//trim(infile))
      tmp = reader%read_field(trim(infile))
      f = f + terms(iterm)%coeff * tmp
    end do
  end subroutine assemble_terms

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
    real(rk), allocatable :: f(:,:,:), f2(:,:,:), f3(:,:,:), x(:), y(:), z(:), profile(:), profile_u(:), profile_v(:)
    character(len=:), allocatable :: outname, stem
    real(rk) :: pi

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    allocate(f(nxloc,nyloc,nzloc), f2(nxloc,nyloc,nzloc), f3(nxloc,nyloc,nzloc), &
      x(nx), y(ny), z(nz), profile(nz), profile_u(nz), profile_v(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)
    pi = acos(-1.0_rk)

    do i = 1, job%nfields
      if (trim(job%fields(i)%derived) == 'ws_wd') then
        call assemble_terms(reader, path, job%fields(i)%u_terms, job%fields(i)%nu_terms, f)
        call assemble_terms(reader, path, job%fields(i)%v_terms, job%fields(i)%nv_terms, f2)
        stem = expr_output_stem(job%fields(i))

        ! WS profile is <sqrt(u*u + v*v)>, not sqrt(<u>**2 + <v>**2).
        f3 = sqrt(f*f + f2*f2)
        call do_horizontal_average(nx, ny, nz, nxloc, nyloc, nzloc, zs, f3, profile)
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
      else
        call assemble_field(reader, path, job%fields(i), f)
        call do_horizontal_average(nx, ny, nz, nxloc, nyloc, nzloc, zs, f, profile)
        if (myrank == 0) then
          stem = expr_output_stem(job%fields(i))
          outname = trim(outdir)//'/'//trim(stem)//'_HA_z.csv'
          call csv_profile(nz, trim(outname), 'z', z, profile)
        end if
      end if
    end do
  end subroutine run_horizontal_average

  subroutine run_rms(reader, Lx, Ly, Lz, path, outdir, job)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: path, outdir
    type(diag_job_t), intent(in) :: job
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, xs, xe, ys, ye, zs, ze
    integer :: ip, i, j, k, ig, jg, kg, ni, profile_axis, ierr
    real(rk) :: dx, dy, dz, area_weight
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), pointer :: coord(:)
    real(rk), allocatable :: f(:,:,:), local_sum(:), global_sum(:), l2norm(:)
    character(len=:), allocatable :: outname, stem

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    call reader%indices(xs, xe, ys, ye, zs, ze)
    allocate(f(nxloc,nyloc,nzloc), x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)
    dx = Lx/real(nx,rk); dy = Ly/real(ny,rk); dz = Lz/real(nz,rk)

    do ip = 1, job%nrms
      call assemble_field(reader, path, job%rms(ip)%field, f)
      select case(job%rms(ip)%axis)
      case('x')
        ni = nx; coord => x; area_weight = dy*dz; profile_axis = 1
      case('y')
        ni = ny; coord => y; area_weight = dx*dz; profile_axis = 2
      case default
        ni = nz; coord => z; area_weight = dx*dy; profile_axis = 3
      end select
      allocate(local_sum(ni), global_sum(ni), l2norm(ni))
      local_sum = 0.0_rk

      select case(job%rms(ip)%axis)
      case('x')
        do i = 1, nxloc
          ig = xs + i - 1
          if (.not. coord_in_bounds(x(ig), job%rms(ip), 1)) cycle
          do j = 1, nyloc
            jg = ys + j - 1
            if (.not. coord_in_bounds(y(jg), job%rms(ip), 2)) cycle
            do k = 1, nzloc
              kg = zs + k - 1
              if (.not. coord_in_bounds(z(kg), job%rms(ip), 3)) cycle
              local_sum(ig) = local_sum(ig) + f(i,j,k)**2 * area_weight
            end do
          end do
        end do
      case('y')
        do j = 1, nyloc
          jg = ys + j - 1
          if (.not. coord_in_bounds(y(jg), job%rms(ip), 2)) cycle
          do i = 1, nxloc
            ig = xs + i - 1
            if (.not. coord_in_bounds(x(ig), job%rms(ip), 1)) cycle
            do k = 1, nzloc
              kg = zs + k - 1
              if (.not. coord_in_bounds(z(kg), job%rms(ip), 3)) cycle
              local_sum(jg) = local_sum(jg) + f(i,j,k)**2 * area_weight
            end do
          end do
        end do
      case default
        do k = 1, nzloc
          kg = zs + k - 1
          if (.not. coord_in_bounds(z(kg), job%rms(ip), 3)) cycle
          do i = 1, nxloc
            ig = xs + i - 1
            if (.not. coord_in_bounds(x(ig), job%rms(ip), 1)) cycle
            do j = 1, nyloc
              jg = ys + j - 1
              if (.not. coord_in_bounds(y(jg), job%rms(ip), 2)) cycle
              local_sum(kg) = local_sum(kg) + f(i,j,k)**2 * area_weight
            end do
          end do
        end do
      end select

      call MPI_Allreduce(local_sum, global_sum, ni, mpi_rk, MPI_SUM, MPI_COMM_WORLD, ierr)
      ! This driver exports the cross-plane L2 norm, sqrt(integral f**2 dA).
      ! Area normalization to convert this profile to RMS is done offline.
      l2norm = sqrt(global_sum)
      if (myrank == 0) then
        stem = expr_output_stem(job%rms(ip)%field)
        outname = trim(outdir)//'/'//trim(stem)//'_rms_'//job%rms(ip)%axis//'.csv'
        call csv_profile_bounded(ni, trim(outname), job%rms(ip)%axis, coord, l2norm, job%rms(ip), profile_axis)
      end if
      deallocate(local_sum, global_sum, l2norm)
    end do
  end subroutine run_rms

  subroutine run_slices(reader, Lx, Ly, Lz, path, outdir, job)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: Lx, Ly, Lz
    character(*), intent(in) :: path, outdir
    type(diag_job_t), intent(in) :: job
    integer :: nx, ny, nz, nxloc, nyloc, nzloc, ifield, ispec, ierr
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), allocatable :: f(:,:,:)
    character(len=:), allocatable :: stem

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    allocate(f(nxloc,nyloc,nzloc), x(nx), y(ny), z(nz))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)

    do ifield = 1, job%nfields
      if (len_trim(job%fields(ifield)%derived) > 0) then
        call message('ERROR: derived fields are not supported by the slice driver.')
        call MPI_Abort(MPI_COMM_WORLD, 811, ierr)
      end if
      call assemble_field(reader, path, job%fields(ifield), f)
      stem = expr_output_stem(job%fields(ifield))
      do ispec = 1, job%nslices
        call process_slice_spec(reader, f, x, y, z, Lx, Ly, Lz, outdir, trim(stem), job%slices(ispec))
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
    integer :: nx, ny, nz, nxloc, nyloc, nzloc
    real(rk), allocatable, target :: x(:), y(:), z(:)
    real(rk), allocatable :: uw(:,:,:), vw(:,:,:), buffer(:,:,:)

    call reader%global_shape(nx, ny, nz)
    call reader%local_shape(nxloc, nyloc, nzloc)
    allocate(x(nx), y(ny), z(nz), buffer(nxloc,nyloc,nzloc))
    call create_grid(Lx, Ly, Lz, nx, ny, nz, x, y, z)

    select case(trim(lower(job%abl%method)))
    case('stress')
      allocate(uw(nxloc,nyloc,nzloc), vw(nxloc,nyloc,nzloc))
      call assemble_field(reader, path, job%abl%uw, uw)
      call assemble_field(reader, path, job%abl%vw, vw)
      buffer = sqrt(uw**2 + vw**2)
      call export_abl_threshold(reader, buffer, x, y, z, outdir, 'ABL_stress', job%abl%threshold)
    case('inversion')
      call assemble_field(reader, path, job%abl%theta, buffer)
      call export_abl_inversion(reader, buffer, x, y, z, outdir, job%abl%lengthscale, &
        job%abl%inversion_l0, job%abl%inversion_d0, job%abl%inversion_xi)
    case default
      call message('ERROR: ABL method must be stress or inversion.')
    end select
  end subroutine run_abl

  subroutine export_abl_threshold(reader, field, x, y, z, outdir, stem, threshold)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: field(:,:,:), x(:), y(:), z(:), threshold
    character(*), intent(in) :: outdir, stem
    real(rk), allocatable :: ytmp(:,:,:), ztmp(:,:,:), local(:,:), global(:,:)
    integer :: nx, ny, nz, ierr
    character(len=:), allocatable :: fname

    call reader%global_shape(nx, ny, nz)
    if (nz /= size(z)) call message('WARNING: ABL threshold z-grid length differs from global nz.')
    allocate(ytmp(reader%gpC%ysz(1), reader%gpC%ysz(2), reader%gpC%ysz(3)))
    allocate(ztmp(reader%gpC%zsz(1), reader%gpC%zsz(2), reader%gpC%zsz(3)))
    allocate(local(nx,ny), global(nx,ny))
    local = 0.0_rk; global = 0.0_rk
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

  subroutine export_abl_inversion(reader, theta, x, y, z, outdir, lengthscale, l0, d0, xi)
    class(FieldReader2Decomp), intent(inout) :: reader
    real(rk), intent(in) :: theta(:,:,:), x(:), y(:), z(:), lengthscale
    real(rk), intent(in) :: l0, d0, xi
    character(*), intent(in) :: outdir
    integer :: nx, ny, nz, ic, jc, ig, jg, ierr
    real(rk), allocatable :: ytmp(:,:,:), ztmp(:,:,:), h0(:,:), h2(:,:), global(:,:)
    real(rk), allocatable :: tcol(:), z_m(:)
    type(rz_params) :: params

    call reader%global_shape(nx, ny, nz)
    allocate(ytmp(reader%gpC%ysz(1), reader%gpC%ysz(2), reader%gpC%ysz(3)))
    allocate(ztmp(reader%gpC%zsz(1), reader%gpC%zsz(2), reader%gpC%zsz(3)))
    allocate(h0(nx,ny), h2(nx,ny), global(nx,ny), tcol(nz), z_m(nz))
    h0 = 0.0_rk; h2 = 0.0_rk; z_m = z * lengthscale
    call transpose_x_to_y(theta, ytmp, reader%gpC)
    call transpose_y_to_z(ytmp, ztmp, reader%gpC)

    do ic = 1, reader%gpC%zsz(1)
      do jc = 1, reader%gpC%zsz(2)
        ig = reader%gpC%zst(1) + ic - 1
        jg = reader%gpC%zst(2) + jc - 1
        tcol = ztmp(ic,jc,:)
        call fit_rz_profile(z_m, tcol, nz, l0, d0, params)
        if ((params%status == 0 .or. params%status == 2) .and. params%sse < huge(1.0_rk)) then
          h0(ig,jg) = params%l - xi * params%d
          h2(ig,jg) = params%l + xi * params%d
        end if
      end do
    end do

    global = 0.0_rk
    call MPI_Reduce(h0, global, nx*ny, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) call export_slice_to_netcdf(trim(outdir)//'/ABL_inversion_h0.nc', 'ABL_h0', global, x, y, 'x', 'y')
    global = 0.0_rk
    call MPI_Reduce(h2, global, nx*ny, mpi_rk, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (myrank == 0) call export_slice_to_netcdf(trim(outdir)//'/ABL_inversion_h2.nc', 'ABL_h2', global, x, y, 'x', 'y')
  end subroutine export_abl_inversion

  subroutine fit_rz_profile(z, t, n, l0, d0, params)
    integer, intent(in) :: n
    real(rk), intent(in) :: z(n), t(n), l0, d0
    type(rz_params), intent(out) :: params
    integer :: iter, i_best
    real(rk) :: l, d, best, trial, step_l, step_d

    l = l0
    d = max(d0, 1.0e-6_rk)
    step_l = max(0.05_rk * (maxval(z)-minval(z)), 1.0_rk)
    step_d = max(0.10_rk * (maxval(z)-minval(z)), 1.0_rk)
    best = rz_objective(z, t, n, l, d)
    do iter = 1, 80
      i_best = 0
      trial = rz_objective(z, t, n, l + step_l, d)
      if (trial < best) then; best = trial; l = l + step_l; i_best = 1; end if
      trial = rz_objective(z, t, n, l - step_l, d)
      if (trial < best) then; best = trial; l = l - step_l; i_best = 1; end if
      trial = rz_objective(z, t, n, l, max(d + step_d, 1.0e-6_rk))
      if (trial < best) then; best = trial; d = max(d + step_d, 1.0e-6_rk); i_best = 1; end if
      trial = rz_objective(z, t, n, l, max(d - step_d, 1.0e-6_rk))
      if (trial < best) then; best = trial; d = max(d - step_d, 1.0e-6_rk); i_best = 1; end if
      if (i_best == 0) then
        step_l = 0.5_rk * step_l
        step_d = 0.5_rk * step_d
      end if
      if (max(step_l, step_d) < 1.0e-6_rk) exit
    end do
    best = rz_objective(z, t, n, l, d, params)
    params%l = l
    params%d = d
    params%sse = best
    if (max(step_l, step_d) < 1.0e-6_rk) then
      params%status = 0
    else
      params%status = 2
    end if
  end subroutine fit_rz_profile

  real(rk) function rz_objective(z, t, n, l, d, params)
    integer, intent(in) :: n
    real(rk), intent(in) :: z(n), t(n), l, d
    type(rz_params), intent(inout), optional :: params
    real(rk) :: A(3,3), bvec(3), sol(3), x, f, g, pred
    integer :: i, status
    A = 0.0_rk; bvec = 0.0_rk
    do i = 1, n
      x = (z(i) - l) / d
      f = 0.5_rk * (tanh(x) + 1.0_rk)
      g = 0.5_rk * (log(2.0_rk*cosh(max(min(x,50.0_rk),-50.0_rk))) + x)
      A(1,1)=A(1,1)+1.0_rk; A(1,2)=A(1,2)+f; A(1,3)=A(1,3)+g
      A(2,2)=A(2,2)+f*f; A(2,3)=A(2,3)+f*g; A(3,3)=A(3,3)+g*g
      bvec(1)=bvec(1)+t(i); bvec(2)=bvec(2)+t(i)*f; bvec(3)=bvec(3)+t(i)*g
    end do
    A(2,1)=A(1,2); A(3,1)=A(1,3); A(3,2)=A(2,3)
    call solve_3x3(A, bvec, sol, status)
    if (status /= 0) then
      rz_objective = huge(1.0_rk)
      return
    end if
    rz_objective = 0.0_rk
    do i = 1, n
      x = (z(i) - l) / d
      f = 0.5_rk * (tanh(x) + 1.0_rk)
      g = 0.5_rk * (log(2.0_rk*cosh(max(min(x,50.0_rk),-50.0_rk))) + x)
      pred = sol(1) + sol(2)*f + sol(3)*g
      rz_objective = rz_objective + (pred - t(i))**2
    end do
    if (present(params)) then
      params%tm = sol(1); params%a = sol(2); params%b = sol(3)
    end if
  end function rz_objective

  subroutine solve_3x3(Ain, bin, x, status)
    real(rk), intent(in) :: Ain(3,3), bin(3)
    real(rk), intent(out) :: x(3)
    integer, intent(out) :: status
    real(rk) :: A(3,4), factor, pivot, pivot_tol
    integer :: i, j, k, p
    A(:,1:3) = Ain; A(:,4) = bin
    pivot_tol = 1.0e-10_rk * max(1.0_rk, maxval(abs(Ain)))
    status = 0
    do k = 1, 3
      p = k
      do i = k+1, 3
        if (abs(A(i,k)) > abs(A(p,k))) p = i
      end do
      if (abs(A(p,k)) <= pivot_tol) then
        status = 1
        x = 0.0_rk
        return
      end if
      if (p /= k) A([k,p],:) = A([p,k],:)
      pivot = A(k,k)
      A(k,k:4) = A(k,k:4) / pivot
      do i = 1, 3
        if (i == k) cycle
        factor = A(i,k)
        do j = k, 4
          A(i,j) = A(i,j) - factor * A(k,j)
        end do
      end do
    end do
    x = A(:,4)
  end subroutine solve_3x3

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
    open(newunit=u, file=trim(filename), status='replace', action='write')
    write(u,'(A)') trim(axis)//',profile'
    do i = 1, n
      write(u,'(ES23.15,",",ES23.15)') coord(i), profile(i)
    end do
    close(u)
  end subroutine csv_profile

  subroutine csv_profile_bounded(n, filename, axis, coord, profile, item, iax)
    integer, intent(in) :: n, iax
    character(*), intent(in) :: filename, axis
    real(rk), intent(in) :: coord(n), profile(n)
    type(rms_spec_t), intent(in) :: item
    integer :: u, i
    open(newunit=u, file=trim(filename), status='replace', action='write')
    write(u,'(A)') trim(axis)//',l2_norm'
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

  subroutine frd_init(this, nx, ny, nz)
    class(FieldReader2Decomp), intent(inout) :: this
    integer, intent(in) :: nx, ny, nz
    if (this%is_init) return
    this%nx = nx; this%ny = ny; this%nz = nz
    call this%choose_proc_grid_(nprocs, this%p_row, this%p_col)
    call decomp_2d_init(nx, ny, nz, this%p_row, this%p_col)
    call decomp_info_init(nx, ny, nz, this%gpC)
    this%xs=xstart(1); this%xe=xend(1)
    this%ys=ystart(1); this%ye=yend(1)
    this%zs=zstart(1); this%ze=zend(1)
    this%nxloc=this%xe-this%xs+1
    this%nyloc=this%ye-this%ys+1
    this%nzloc=this%ze-this%zs+1
    this%is_init = .true.
  end subroutine frd_init

  function frd_read_field(this, path) result(field)
    class(FieldReader2Decomp), intent(inout) :: this
    character(*), intent(in) :: path
    real(rk), allocatable :: field(:,:,:)
    integer :: ierr
    if (.not. this%is_init) then
      call message('ERROR: FieldReader2Decomp is not initialized before read_field().')
      call MPI_Abort(MPI_COMM_WORLD, 101, ierr)
    end if
    allocate(field(this%nxloc,this%nyloc,this%nzloc))
    call decomp_2d_read_one(1, field, trim(path), this%gpC)
  end function frd_read_field

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

  subroutine frd_choose_proc_grid(nproc, p_row, p_col)
    integer, intent(in) :: nproc
    integer, intent(out) :: p_row, p_col
    integer :: r
    p_row = 1; p_col = nproc
    do r = nint(sqrt(real(nproc, rk))), 1, -1
      if (mod(nproc, r) == 0) then
        p_row = r
        p_col = nproc / r
        exit
      end if
    end do
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
  case('abl')
    call run_abl(reader, Lx, Ly, Lz, trim(path), trim(outdir), job)
  case default
    call message('ERROR: unknown driver '//trim(job%driver))
    call MPI_Abort(MPI_COMM_WORLD, 5, ierr)
  end select

  call message('Wrapping up ...')
  call MPI_Finalize(ierr)
end program MPIR3D_Lean_Main
