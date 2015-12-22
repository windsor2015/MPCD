#ifdef __INTEL_COMPILER
include 'mkl_vsl.f90'
#endif

module shape_all
#ifdef __INTEL_COMPILER
    use ifport
    USE MKL_VSL_TYPE
    USE MKL_VSL
#endif
    ! 下标约定: p - polymer, s - solution, b - boundary / phantom

    !结构
    integer, parameter :: n_p=40, n_cell_x=10, n_cell_y=20, n_cell_z=50

    integer n_b, n_s
        !n是单体数目，n0是单体、壁、孔数目之和

    real(8), parameter :: kB = 1.38064852d-23
    real(8), parameter :: pi = 3.141592653589793238462643383279

    real(8), parameter :: Ek_fac = 1.5d0, ratio_y=0.3, ratio_z=0.4

    real(8), parameter :: time_step_p=1d-4, time_step_s=5d-3, mass_p=1, mass_s=0.2, T_set=1, v_gradient=0.2

    ! polymer 位置 速度 力 上一次力
    real(8), dimension(3,n_p) :: x_p, v_p, f_p, f0_p, x0_p
    ! solution
    real(8), allocatable, dimension(:,:) :: x_s, v_s, f_s, x0_s, x_s0
    ! boundaaries, 1~nb-up, nb+1~2nb-down
    real(8), allocatable, dimension(:,:) :: x_b, v_b, f_b

    real(8) box_size_unit,  box_size(3)    !!, half_box_size(3)   !!half_box_size_unit,

    real(8) aver_v(3), distant

    real(8), parameter :: sigma=1, epson=1

    real(8) :: density_s, gama

    real(8) U, U_LJ, U_FENE, U_BEND, U_WALL

    ! pointer - 每个格子中的粒子编号
    ! count - 每个格子中的粒子计数，1 - p，2 - s
    ! momentum - 每个格子中的总动量
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z, n_p) :: pointer_cell_p
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z, 100) :: pointer_cell_s, pointer_cell_b
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z) :: count_cell_p, count_cell_s, count_cell_b
    real(8), dimension(3,0:n_cell_x, 0:n_cell_y, 0:n_cell_z) :: momentum_cell

    integer :: desk_interval_step
    integer :: equili_step
    integer :: equili_interval_step
    integer :: total_step
    integer :: output_interval_step

    integer :: string_form = 0

    integer :: thermostat_method, thermostat_interval
    real(8) thermostat_B_parameter, thermostat_A_parameter


    real(8) :: BEND_b=10d0

    character(50) :: equili_filename='dump.equili.lammpstrj'
    character(50) :: production_filename='dump.production.lammpstrj'

    namelist /basic/ gama, density_s, &
        desk_interval_step, equili_step, equili_interval_step, total_step, output_interval_step, &
        thermostat_method, thermostat_interval, thermostat_B_parameter,thermostat_A_parameter, &
        string_form, BEND_b, equili_filename, production_filename

#ifdef __INTEL_COMPILER
    TYPE (VSL_STREAM_STATE) :: vsl_stream
#endif
    integer :: debug=0
    real(8) :: time0=0

    real(8) sum_v(3,0:100,0:400),sum_grid_v(3,0:100,0:100,0:400)
    integer n(0:100,0:400),n_grid(0:100,0:100,0:400)


    contains

end module


