program mandelbrot
  use mpi
  use write_image
  use read_opt
  implicit none

  integer, allocatable, dimension(:,:) :: image

  integer :: gridsizex, gridsizey
  integer :: iter
  integer :: ierr
  integer :: rank
  integer :: comm
  complex :: cmin
  complex :: cmax

  call MPI_Init(ierr)
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call read_options(gridsizex, gridsizey, iter, cmin, cmax)

  if ( rank .eq. 0 ) then
     call initialise_image(image, gridsizex, gridsizey)
  end if

  call compute_mandelbrot_set(image, cmin, cmax, gridsizex, gridsizey, iter)

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

    ! We send these values to the master process for writing because
    ! we don't assume that the mapping from rank to slice is
    ! one-to-one.  For example, we might divide the domain into twice
    ! as many slices as there are processes
    call calc_slice_bounds(slice, nslice, grid_size_y, start, end)

    if ( rank .eq. 0 ) then
       ! Copy local slice into global image
       image(:, start:end-1) = image_slice
       call MPI_Comm_size(comm, nproc, ierr)
       do i = 1, nproc-1
          call recv_image_slice(image, i, comm)
       end do
    else
       call send_image_slice(image_slice, start, end, 0, comm)
    end if
  end subroutine copy_slice_to_image

  pure function point_in_mandelbrot_set(c, max_iter) result(ret)
    complex, intent(in) :: c
    integer, intent(in) :: max_iter
    integer :: ret
    integer :: i
    complex :: z

    z = c

    do i = 1, max_iter
       if ( real(z * conjg(z)) .gt. 4.0 ) then
          ret = i-1
          return
       else
          z = z**2 + c
       end if
    end do
    ret = max_iter
  end function point_in_mandelbrot_set

  subroutine compute_mandelbrot_slice(image_slice, slice, nslice, &
       min, max, grid_size_x, grid_size_y, max_iter)
    integer, allocatable, dimension(:,:), intent(out) :: image_slice
    integer, intent(in) :: slice, nslice
    complex, intent(in) :: min, max
    integer, intent(in) :: grid_size_x, grid_size_y
    integer, intent(in) :: max_iter

    integer :: i, j
    integer :: start, end
    complex :: c
    real :: rpart
    real :: ipart

    call calc_slice_bounds(slice, nslice, grid_size_y, start, end)

    call initialise_image(image_slice, grid_size_x, end - start)
    rpart = real(max - min) / grid_size_x
    ipart = aimag(max - min) / grid_size_y

    do j = 1, end - start
       do i = 1, grid_size_x
          ! Calculate coordinate value of current pixel
          c = min + cmplx(rpart * real(i-1), ipart * real(j + start - 1))
          image_slice(i, j) = point_in_mandelbrot_set(c, max_iter)
       end do
    end do
  end subroutine compute_mandelbrot_slice

  subroutine compute_mandelbrot_set(image, min, max, grid_size_x, &
       grid_size_y, max_iter)
    integer, dimension(:,:), intent(inout) :: image
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
    call compute_mandelbrot_slice(image_slice, slice, nslice, &
         min, max, grid_size_x, grid_size_y, max_iter)

    call copy_slice_to_image(image_slice, image, grid_size_y, slice, nslice)
    if ( allocated(image_slice) ) deallocate(image_slice)
  end subroutine compute_mandelbrot_set

end program mandelbrot
