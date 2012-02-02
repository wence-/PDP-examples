module fnptr
  implicit none
  abstract interface
     integer function fn(c, max_iter)
       complex, intent(in) :: c
       integer, intent(in) :: max_iter
     end function fn
  end interface

contains

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

  pure function point_in_julia_set(c, max_iter) result(ret)
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
          z = z**2 + cmplx(0.285, 0.01)
       end if
    end do
    ret = max_iter
  end function point_in_julia_set

end module fnptr

module routines
  use mpi
  use fnptr
  implicit none
  integer, parameter :: EMPTY_TASK = -1

contains

  subroutine build_task(lines, task, grid_size_y)
    integer, dimension(2), intent(out) :: lines
    integer, intent(inout) :: task
    integer, intent(in) :: grid_size_y

    integer :: lines_in_task = 10

    if ( task .le. grid_size_y ) then
       lines(1) = task
       lines(2) = min(task + lines_in_task, grid_size_y)
    else
       lines = EMPTY_TASK
    end if

    task = task + lines_in_task
  end subroutine build_task

  ! Initialise the data for the image array.  Data is stored in
  ! "scanline order", i.e. x dimension varies fastest.  You get an
  ! array with shape (grid_size_x, grid_size_y) from this function
  subroutine initialise_image(image, grid_size_x, grid_size_y)
    integer, intent(out), allocatable, dimension(:,:) :: image
    integer, intent(in) :: grid_size_x, grid_size_y

    allocate(image(grid_size_x, grid_size_y))
    image = -1
  end subroutine initialise_image

  subroutine send_image_lines(image_lines, lines, comm)
    integer, intent(in), dimension(:,:) :: image_lines
    integer, intent(in), dimension(2) :: lines
    integer, intent(in) :: comm

    integer :: ierr

    call MPI_Send(lines, 2, MPI_INTEGER, 0, 0, comm, ierr)
    call MPI_Send(image_lines, size(image_lines, 1) * (lines(2) - lines(1)), &
         MPI_INTEGER, 0, 0, comm, ierr)
  end subroutine send_image_lines

  subroutine recv_image_lines(image, worker, comm)
    integer, intent(inout), dimension(:,:) :: image
    integer, intent(out) :: worker
    integer, intent(in) :: comm
    integer :: ierr, lines(2)
    integer, dimension(MPI_STATUS_SIZE) :: status

    call MPI_Recv(lines, 2, MPI_INTEGER, MPI_ANY_SOURCE, 0, &
         comm, status, ierr)
    call MPI_Recv(image(1, lines(1)), &
         size(image, 1) * (lines(2) - lines(1)), &
         MPI_INTEGER, status(MPI_SOURCE), 0, comm, status, ierr)
    worker = status(MPI_SOURCE)
  end subroutine recv_image_lines

  subroutine compute_lines(in_set_fn, data, lines, &
       min, max, grid_size_x, grid_size_y, max_iter)
    procedure(fn), pointer :: in_set_fn
    integer, allocatable, dimension(:,:), intent(out) :: data
    integer, intent(in), dimension(2) :: lines
    complex, intent(in) :: min, max
    integer, intent(in) :: grid_size_x, grid_size_y
    integer, intent(in) :: max_iter

    integer :: i, j
    integer :: start, end
    complex :: c
    real :: rpart
    real :: ipart

    start = lines(1)
    end = lines(2)
    call initialise_image(data, grid_size_x, end - start)
    rpart = real(max - min) / grid_size_x
    ipart = aimag(max - min) / grid_size_y

    do j = 1, end - start
       do i = 1, grid_size_x
          ! Calculate coordinate value of current pixel
          c = min + cmplx(rpart * real(i-1), ipart * real(j + start - 1))
          data(i, j) = in_set_fn(c, max_iter)
       end do
    end do
  end subroutine compute_lines

  subroutine master_loop(image, min, max, grid_size_x, grid_size_y, &
       max_iter, comm)
    integer, dimension(:,:), intent(inout) :: image
    complex, intent(in) :: min, max
    integer, intent(in) :: grid_size_x, grid_size_y, max_iter
    integer, intent(in) :: comm
    integer :: ierr, nproc
    integer :: i, outstanding_tasks, current_task, worker
    integer, dimension(2) :: lines
    complex, dimension(2) :: ctmp
    integer, dimension(3) :: itmp

    call MPI_Comm_size(comm, nproc, ierr)
    ctmp = (/ min, max /)
    call MPI_Bcast(ctmp, 2, MPI_COMPLEX, 0, comm, ierr)
    itmp = (/ grid_size_x, grid_size_y, max_iter /)
    call MPI_Bcast(itmp, 3, MPI_INTEGER, 0, comm, ierr)

    outstanding_tasks = 0
    current_task = 1

    do i = 1, nproc - 1
       call build_task(lines, current_task, grid_size_y)
       call MPI_Send(lines, 2, MPI_INTEGER, i, 0, comm, ierr)
       if ( lines(1) .ne. EMPTY_TASK ) then
          outstanding_tasks = outstanding_tasks + 1
       end if
    end do

    do
       if ( outstanding_tasks .eq. 0 ) exit
       call recv_image_lines(image, worker, comm)
       outstanding_tasks = outstanding_tasks - 1
       call build_task(lines, current_task, grid_size_y)
       call MPI_Send(lines, 2, MPI_INTEGER, worker, 0, comm, ierr)
       if ( lines(1) .ne. EMPTY_TASK ) then
          outstanding_tasks = outstanding_tasks + 1
       end if
    end do
  end subroutine master_loop

  subroutine worker_loop(in_set_fn, comm)
    procedure(fn), pointer :: in_set_fn
    integer, intent(in) :: comm
    complex :: min, max, ctmp(2)
    integer :: grid_size_x, grid_size_y, max_iter, itmp(3)
    integer :: ierr
    integer, dimension(2) :: lines
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(:,:), allocatable :: data

    call MPI_Bcast(ctmp, 2, MPI_COMPLEX, 0, comm, ierr)
    min = ctmp(1)
    max = ctmp(2)
    call MPI_Bcast(itmp, 3, MPI_INTEGER, 0, comm, ierr)
    grid_size_x = itmp(1)
    grid_size_y = itmp(2)
    max_iter = itmp(3)

    do
       call MPI_Recv(lines, 2, MPI_INTEGER, 0, 0, comm, status, ierr)
       if ( lines(1) .eq. EMPTY_TASK ) then
          return
       else
          call compute_lines(in_set_fn, data, lines, min, max, &
               grid_size_x, grid_size_y, max_iter)
          call send_image_lines(data, lines, comm)
          deallocate(data)
       end if
    end do
  end subroutine worker_loop

  subroutine compute_set(in_set_fn, image, min, max, grid_size_x, &
       grid_size_y, max_iter, comm)
    integer, dimension(:,:), allocatable, intent(out) :: image
    procedure(fn), pointer :: in_set_fn
    complex, intent(in) :: min, max
    integer, intent(in) :: grid_size_x, grid_size_y
    integer, intent(in) :: max_iter
    integer, intent(in) :: comm
    integer :: rank, ierr

    call MPI_Comm_rank(comm, rank, ierr)

    if ( rank .eq. 0 ) then
       call initialise_image(image, grid_size_x, grid_size_y)
       call master_loop(image, min, max, grid_size_x, grid_size_y, &
            max_iter, comm)
    else
       call worker_loop(in_set_fn, comm)
    end if
  end subroutine compute_set

end module routines

program fractal_set
  use mpi
  use write_image
  use read_opt
  use fnptr
  use routines
  implicit none

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

  fp => point_in_julia_set
  call compute_set(fp, image, cmin, cmax, gridsizex, gridsizey, iter, comm)

  if ( rank .eq. 0 ) then
     call write_ppm("output.ppm", image, iter)
     deallocate(image)
  end if

  call MPI_Finalize(ierr)

end program fractal_set
