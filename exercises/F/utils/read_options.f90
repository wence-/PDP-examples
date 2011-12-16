module read_opt
  implicit none
  complex, parameter :: MIN = (-2.0, -1.2)
  complex, parameter :: MAX = (1, 1.2)

  integer, parameter :: GRIDSIZE_X = 768
  integer, parameter :: GRIDSIZE_Y = int(GRIDSIZE_X * &
       aimag(MAX - MIN) / real(MAX - MIN))

  integer, parameter :: ITERATIONS = 5000
contains

  subroutine read_options(gridsizex, gridsizey, iter, cmin, cmax)
    integer, intent(out) :: gridsizex
    integer, intent(out) :: gridsizey
    integer, intent(out) :: iter
    complex, intent(out) :: cmin
    complex, intent(out) :: cmax

    integer :: i
    integer :: nargs
    real :: tmp
    character(len=128) :: arg
    character(len=128) :: opt
    nargs = command_argument_count()
    i = 1

    iter = ITERATIONS
    gridsizex = GRIDSIZE_X
    cmin = MIN
    cmax = MAX

    do
       call get_command_argument(i, value=arg)
       call get_command_argument(i+1, value=opt)

       select case(trim(arg(2:)))
       case ('S')
          read(opt, *)gridsizex
          if ( gridsizex < 4 ) gridsizex = GRIDSIZE_X
       case ('i')
          read(opt, *)iter
       case ('x')
          read(opt, *)tmp
          cmin = cmplx(tmp, aimag(cmin))
       case ('X')
          read(opt, *)tmp
          cmax = cmplx(tmp, aimag(cmax))
       case ('y')
          read(opt, *)tmp
          cmin = cmplx(real(cmin), tmp)
       case ('Y')
          read(opt, *)tmp
          cmax = cmplx(real(cmax), tmp)
       case default
          write(*, *)'Usage:'
          write(*, *)'mandelbrot [-SixXyY]'
          write(*, *)'   -S NPIXEL    Set number of pixels in X dimension of image'
          write(*, *)'                Y dimension is scaled to ensure pixels are square'
          write(*, *)'   -i ITS       Set max number of iterations for a point to be inside'
          write(*, *)'   -x XMIN      Set xmin coordinate'
          write(*, *)'   -X XMAX      Set xmax coordinate'
          write(*, *)'   -y YMIN      Set ymin coordinate'
          write(*, *)'   -Y YMAX      Set ymax coordinate'
          write(*, *)'   -h           Show this help'
          stop
       end select
       i = i + 2
       if ( i .gt. nargs ) exit
    end do

    if ( real(cmin) .ge. real(cmax) ) cmin = cmplx(real(MIN), aimag(cmin))
    if ( aimag(cmin) .ge. aimag(cmax) ) cmin = cmplx(real(cmin), aimag(MIN))
    gridsizey = int(gridsizex * aimag(cmax - cmin) / real(cmax - cmin))

  end subroutine read_options
end module read_opt
