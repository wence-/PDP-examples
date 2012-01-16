program mandelbrot

  use write_image
  use read_opt
  implicit none

  integer, allocatable, dimension(:,:) :: image

  integer :: gridsizex, gridsizey
  integer :: max_iter
  complex :: cmin
  complex :: cmax

  call read_options(gridsizex, gridsizey, max_iter, cmin, cmax)

  call initialise_image(image, gridsizex, gridsizey)

  call compute_mandelbrot_set(image, cmin, cmax, gridsizex, gridsizey, max_iter)

  call write_ppm("output.ppm", image, max_iter)

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

  subroutine compute_mandelbrot_set(image, min, max, grid_size_x, &
       grid_size_y, max_iter)
    integer, dimension(:,:), intent(inout) :: image
    complex, intent(in) :: min, max
    integer, intent(in) :: grid_size_x, grid_size_y
    integer, intent(in) :: max_iter

    integer :: i
    integer :: j
    complex :: c
    real :: rpart
    real :: ipart

    rpart = real(max - min) / grid_size_x
    ipart = aimag(max - min) / grid_size_y
    do j = 1, grid_size_y
       do i = 1, grid_size_x
          ! Calculate coordinate value of current pixel
          c = min + cmplx(rpart * real(i-1), ipart * real(j-1))
          image(i, j) = point_in_mandelbrot_set(c, max_iter)
       end do
    end do
  end subroutine compute_mandelbrot_set
end program mandelbrot
