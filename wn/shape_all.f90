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
    integer, parameter :: n_cell_x=10, n_cell_y=20, n_cell_z=40, n_p = 40

    integer :: string_form = 0,  n_b, n_s

    real(8) :: density_s, gama, ratio_y, ratio_z, BEND_b

    real(8), parameter :: kB = 1.38064852d-23
    real(8), parameter :: pi = 3.141592653589793238462643383279

    real(8), parameter :: Ek_fac = 1.5d0

    real(8), parameter :: time_step_p=1d-4, time_step_s=5d-3, mass_p=1d0, mass_s=0.2d0, T_set=1d0, v_gradient=0.2d0

    ! polymer 位置 速度 力 上一次力
    real(8), dimension(3,n_p) :: x_p, v_p, f_p, f0_p, x0_p
    ! solution
    real(8), allocatable, dimension(:,:) :: x_s, v_s, x0_s, x_s0
    ! boundaaries, 1~nb-up, nb+1~2nb-down
    real(8), allocatable, dimension(:,:) :: x_b, v_b

    real(8) box_size_unit,  box_size(3), aver_v(3), distant

    real(8), parameter :: sigma=1, epson=1

    real(8) U, U_LJ, U_FENE, U_BEND, U_WALL
    ! pointer - 每个格子中的粒子编号
    ! count - 每个格子中的粒子计数，1 - p，2 - s
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z, n_p) :: pointer_cell_p
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z, 100) :: pointer_cell_s, pointer_cell_b
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z) :: count_cell_p, count_cell_s, count_cell_b
    real(8), dimension(3,0:n_cell_x, 0:n_cell_y, 0:n_cell_z) :: momentum_cell

    integer :: desk_interval_step
    integer :: equili_step
    integer :: equili_interval_step
    integer :: total_step
    integer :: produ_interval_step

    integer :: thermostat_method, thermostat_interval
    real(8) thermostat_B_parameter, thermostat_A_parameter

    real(8), dimension(2)::radius=(/4, 4/)

    character(50) :: equili_filename='dump.equili.lammpstrj'
    character(50) :: production_filename='dump.production.lammpstrj'

    namelist /basic/ gama, density_s, &
        desk_interval_step, equili_step, equili_interval_step, total_step, produ_interval_step, &
        thermostat_method, thermostat_interval, thermostat_B_parameter,thermostat_A_parameter, &
        string_form, BEND_b, equili_filename, production_filename,radius,ratio_z,ratio_y

#ifdef __INTEL_COMPILER
    TYPE (VSL_STREAM_STATE) :: vsl_stream
#endif
    integer :: debug=0
    real(8) :: time0=0

    contains

end module


