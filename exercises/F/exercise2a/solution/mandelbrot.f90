program mandelbrot

  use write_image
  use read_opt
  implicit none

  integer, allocatable, dimension(:,:) :: image

  integer :: gridsizex, gridsizey
  integer :: iter
  complex :: cmin
  complex :: cmax

  call read_options(gridsizex, gridsizey, iter, cmin, cmax)

  call initialise_image(image, gridsizex, gridsizey)

  call compute_mandelbrot_set(image, cmin, cmax, gridsizex, gridsizey, iter)

  call write_ppm("output.ppm", image, iter)

  deallocate(image)

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

  subroutine copy_slice_to_image(image_slice, image, &
       slice, nslice)
    integer, dimension(:,:), intent(in) :: image_slice
    integer, dimension(:,:), intent(inout) :: image
    integer, intent(in) :: slice, nslice

    integer :: start, end

    call calc_slice_bounds(slice, nslice, size(image, 2), start, end)

    image(:, start:end-1) = image_slice
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
          ret = i
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

    integer :: slice
    integer :: nslice
    integer, allocatable, dimension(:,:) :: image_slice

    ! Arbitrary number of slices
    nslice = 3

    do slice = 1, nslice
       call compute_mandelbrot_slice(image_slice, slice, nslice, &
            min, max, grid_size_x, grid_size_y, max_iter)
       call copy_slice_to_image(image_slice, image, slice, nslice)
       if ( allocated(image_slice) ) deallocate(image_slice)
    end do
  end subroutine compute_mandelbrot_set

end program mandelbrot
