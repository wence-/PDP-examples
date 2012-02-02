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
  end function point_in_mandelbrot_set

  pure function point_in_julia_set(c, max_iter) result(ret)
    complex, intent(in) :: c
    integer, intent(in) :: max_iter
    integer :: ret
  end function point_in_julia_set

end module fnptr

module routines
  use mpi
  use fnptr
  implicit none

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

  subroutine master_loop(image, min, max, grid_size_x, grid_size_y, &
       max_iter, comm)
    integer, dimension(:,:), intent(inout) :: image
    complex, intent(in) :: min, max
    integer, intent(in) :: grid_size_x, grid_size_y, max_iter
    integer, intent(in) :: comm
  end subroutine master_loop

  subroutine worker_loop(in_set_fn, comm)
    procedure(fn), pointer :: in_set_fn
    integer, intent(in) :: comm
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
