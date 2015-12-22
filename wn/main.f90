program Poisellie_field
    use parameters

    implicit none
    integer :: cur_step,output_file,energy_file,production_file,velocity_file,coord_velo_file
    integer i,j,k,h_p
    real(8) :: EK_scaled,T_scaled,r,t,t0

    output_file=912
    energy_file=913
    production_file=914
    velocity_file=915
    coord_velo_file=916
    !gama=0.001

    call report()

    call readin()

    call get_time(t0)
    time0=t0
    write(*,*) 'begin at'
    call output_date()
    write(*,*)

    box_size = [n_cell_x, n_cell_y, n_cell_z]
    box_size_unit=1.0

    !!!读链的大小 改成1个文件
    open(output_file,file=equili_filename)
    open(energy_file,file='energy.out')
    open(production_file,file=production_filename)
    open(velocity_file,file='velocity_radius')
    open(coord_velo_file,file='coordinate_velocity')

    call init()
    call thermostat_init()
    call output(output_file,0,equili_interval_step)

    !!!polymer的初速度
    call random_number(v_p)
    !!!solution的初速度
    call random_number(v_s)
    v_p=v_p-0.5
    v_s=v_s-0.5
    ! write(*,*)v_p(:,1)
    call thermostat_I()
    !call cal_collision_velocity(0)
    call Ek_T(EK_scaled,T_scaled)

    write(*,*) ''
    write(*,*) 'Initial scaled kinetics Ek_scaled: ',Ek_scaled
    write(*,*) 'Initial setted temperature T_set: ',T_set
    write(*,*) 'Initial scaled temperature T_scaled: ',T_scaled

    ! 没有外场时，polymer和solution达到平衡
    call update_force(0)
    f_s=0
    t=0
    write(*,*) ''
    write(*,*)'Equilibrium begin:'
    write(*,'(A7,6A12)') 'step', 'BEND','FENE','LJ','total','T_scaled','time'
    write(*,*) '-------------------------------------------------------------------------------'
    write(*,'(I7,6F12.3)') 0, U_BEND, U_FENE, U_LJ, U,T_scaled,t

    call clear_stat()
    do cur_step=1,equili_step
        ! write(*,*) U
        call one_step(cur_step,equili_interval_step, output_file)
        call stat_velocity(cur_step)
    enddo
    call output_velocity(1,velocity_file,coord_velo_file, equili_step)

    write(*,*)
    write(*,*)'Production begin:'

    !!! compute a(t-dt)

    call update_force(0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,'(A7,6A12)') 'step', 'BEND','FENE','LJ', 'total','T_scaled','time'

    write(*,*) '-------------------------------------------------------------------------------'

    call clear_stat()
    do cur_step=1,total_step
        do i=1,n_s
            if (x_s(i,3)>=-25 .and. x_s(i,3)<-23) then
                v_s(3,i) = v_s(3,i) + gama !- gama*(x_s(1,:)**2+x_s(2,:)**2)/radius**2
            end if
        end do

        call one_step(cur_step, output_interval_step,production_file)
        call stat_velocity(cur_step)
    enddo
    call output_velocity(2,velocity_file,coord_velo_file,total_step)


    close(output_file)
    !close(energy_file)
    close(production_file)
    close(velocity_file)
    close(coord_velo_file)
    write(*,*)
    write(*,*) 'end at'
    call output_date()
    call get_time(t)
    write(*,*) 'elapsed time is', t-t0, 's'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program Poisellie_field
