program mandelbrot
  use mpi
  use write_image
  use read_opt

  implicit none

  abstract interface
     integer function fn(c, max_iter)
       complex, intent(in) :: c
       integer, intent(in) :: max_iter
     end function fn
  end interface

  integer, allocatable, dimension(:,:) :: image

  integer :: gridsizex, gridsizey
  integer :: iter
  integer :: ierr
  integer :: rank
  integer :: comm
  complex :: cmin
  complex :: cmax

  procedure(fn), pointer :: fp

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call read_options(gridsizex, gridsizey, iter, cmin, cmax)

  if ( rank .eq. 0 ) then
     call initialise_image(image, gridsizex, gridsizey)
  end if

  ! Set fp to point to the correct function (point_in_mandelbrot_set
  ! or point_in_julia_set)

  call compute_set(fp, image, cmin, cmax, gridsizex, gridsizey, iter)

  if ( rank .eq. 0 ) then
     call write_ppm("output.ppm", image, iter)
     deallocate(image)
  end if

  call MPI_Finalize(ierr)

contains

  ! Initialise the data for the image array.  Data is stored in
  ! "scanline order", i.e. x dimension varies fastest.  You get an
  ! array with shape (grid_size_x, grid_size_y) from this function
  subroutine initialise_image(image, grid_size_x, grid_size_y)
    integer, intent(out), allocatable, dimension(:,:) :: image
    integer, intent(in) :: grid_size_x, grid_size_y

    allocate(image(grid_size_x, grid_size_y))
    image = -1
  end subroutine initialise_image

  subroutine calc_slice_bounds(slice, nslice, grid_size_y, &
       start, end)
    integer, intent(in) :: slice, nslice, grid_size_y
    integer, intent(out) :: start, end

    start = 1 + ((slice - 1) * grid_size_y / nslice)
    end = 1 + (slice * grid_size_y / nslice)
  end subroutine calc_slice_bounds

  subroutine send_image_slice(image_slice, start, end, dest, comm)
    integer, intent(in), dimension(:,:) :: image_slice
    integer, intent(in) :: dest
    integer, intent(in) :: comm
    integer, intent(in) :: start, end

    integer :: ierr

    call MPI_Send(start, 1, MPI_INTEGER, dest, 0, comm, ierr)
    call MPI_Send(end, 1, MPI_INTEGER, dest, 0, comm, ierr)
    call MPI_Send(image_slice, size(image_slice, 1) * (end - start), &
         MPI_INTEGER, dest, 0, comm, ierr)
  end subroutine send_image_slice

  subroutine recv_image_slice(image, source, comm)
    integer, intent(inout), dimension(:,:) :: image
    integer, intent(in) :: source
    integer, intent(in) :: comm
    
    integer :: start, end, ierr
    integer, dimension(MPI_STATUS_SIZE) :: status

    call MPI_Recv(start, 1, MPI_INTEGER, source, 0, comm, status, ierr)
    call MPI_Recv(end, 1, MPI_INTEGER, source, 0, comm, status, ierr)
    call MPI_Recv(image(1, start), size(image, 1) * (end - start), &
         MPI_INTEGER, source, 0, comm, status, ierr)

  end subroutine recv_image_slice

  subroutine copy_slice_to_image(image_slice, image, grid_size_y, &
       slice, nslice)
    integer, dimension(:,:), intent(in) :: image_slice
    integer, dimension(:,:), intent(inout) :: image
    integer, intent(in) :: slice, nslice, grid_size_y

    integer :: comm, rank, nproc, ierr
    integer :: start, end, i

    comm = MPI_COMM_WORLD
    call MPI_Comm_rank(comm, rank, ierr)

    if ( rank .eq. 0 ) then
       ! Copy local slices into global image
       ! Receive remote slices into global image
    else
       ! Send slices to rank 0
    end if
  end subroutine copy_slice_to_image

  pure function point_in_mandelbrot_set(c, max_iter) result(ret)
    complex, intent(in) :: c
    integer, intent(in) :: max_iter
    integer :: ret
    ret = 0
  end function point_in_mandelbrot_set

  pure function point_in_julia_set(c, max_iter) result(ret)
    complex, intent(in) :: c
    integer, intent(in) :: max_iter
    integer :: ret
    ret = 0
  end function point_in_julia_set

  subroutine compute_slice(in_set_fn, image_slice, slice, nslice, &
       min, max, grid_size_x, grid_size_y, max_iter)
    procedure(fn) :: in_set_fn
    integer, allocatable, dimension(:,:), intent(out) :: image_slice
    integer, intent(in) :: slice, nslice
    complex, intent(in) :: min, max
    integer, intent(in) :: grid_size_x, grid_size_y
    integer, intent(in) :: max_iter

    ! Compute the set given by IN_SET_FN in a SLICE of the whole domain.
    ! The total number of slices is given by NSLICE.
    !
    ! The utility function calc_slice_bounds is provided to work out
    ! the bounds of a slice.  It divides the image in the y direction.
    !
    ! This subroutine should return a newly initialised image slice.
  end subroutine compute_slice

  subroutine compute_set(in_set_fn, image, min, max, grid_size_x, &
       grid_size_y, max_iter)
    integer, dimension(:,:), intent(inout) :: image
    procedure(fn) :: in_set_fn
    complex, intent(in) :: min, max
    integer, intent(in) :: grid_size_x, grid_size_y
    integer, intent(in) :: max_iter
    integer :: rank, nproc, comm, ierr
    integer :: slice
    integer :: nslice
    integer, allocatable, dimension(:,:) :: image_slice

    comm = MPI_COMM_WORLD
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nproc, ierr)

    nslice = nproc
    slice = rank + 1
    call compute_slice(in_set_fn, image_slice, slice, nslice, &
         min, max, grid_size_x, grid_size_y, max_iter)

    call copy_slice_to_image(image_slice, image, &
         grid_size_y, slice, nslice)
    if ( allocated(image_slice) ) deallocate(image_slice)
  end subroutine compute_set

end program mandelbrot
