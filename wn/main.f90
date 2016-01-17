program Poisellie_field
    use parameters

    implicit none
    integer :: cur_step,equi_file,energy_file,produ_file,velocity_file,coord_velo_file
    integer i,j,k,h_p
    real(8) :: EK_scaled,T_scaled,r,t,t0,tc0,tc1


    equi_file=912
    !energy_file=913
    produ_file=914
    velocity_file=915
    coord_velo_file=916

    call report()
    call readin()

    call get_time(t0)
    call cpu_time(tc0)
    time0=t0
    write(*,*) 'begin at'
    call output_date()
    write(*,*)

    box_size = [n_cell_x, n_cell_y, n_cell_z]
    box_size_unit=1.0

    !!!读链的大小 改成1个文件
    open(equi_file,file=equili_filename)
    !open(energy_file,file='energy.out')
    open(produ_file,file=production_filename)
    open(velocity_file,file='velocity_radius')
    open(coord_velo_file,file='coordinate_velocity')

    call init()
    call thermostat_init()
    call output(equi_file,0,equili_interval_step)

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
    call update_force(0,0)

    t=0
    write(*,*) ''
    write(*,*)'Equilibrium begin:'
    write(*,'(A7,6A12)') 'step', 'BEND','FENE','LJ','total','T_scaled','time'
    write(*,*) '-------------------------------------------------------------------------------'
    write(*,'(I7,6F12.3)') 0, U_BEND, U_FENE, U_LJ, U,T_scaled,t

    call clear_stat()
    do cur_step=1,equili_step
        ! write(*,*) U
        call one_step(cur_step,equili_interval_step, equi_file,0)
        call stat_velocity(cur_step,equili_interval_step)
    enddo
    call output_velocity(1,velocity_file,coord_velo_file, equili_step,equili_interval_step)

    write(*,*)
    write(*,*)'Production begin:'

    !!! compute a(t-dt)

    call update_force(0,1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,'(A7,6A12)') 'step', 'BEND','FENE','LJ', 'total','T_scaled','time'

    write(*,*) '-------------------------------------------------------------------------------'

    call clear_stat()
    knot_flag=0
    do cur_step=1,total_step

        if (field_interval==0 .or. mod(cur_step, field_interval*2)<=field_interval) then
            v_s(3,:) = v_s(3,:) + gama !- gama*(x_s(1,:)**2+x_s(2,:)**2)/radius**2
        end if

        call one_step(cur_step, produ_interval_step,produ_file,1)
        call stat_velocity(cur_step,produ_interval_step)

        if (knot_flag) then
            write(*,*) 'the knot has untied'
            exit
        end if

    enddo
    call output_velocity(2,velocity_file,coord_velo_file,total_step,produ_interval_step)

    close(equi_file)
    !close(energy_file)
    close(produ_file)
    close(velocity_file)
    close(coord_velo_file)
    write(*,*)
    write(*,*) 'end at'
    call output_date()
    call get_time(t)
    call cpu_time(tc1)
    write(*,'(A,$)') ' wall time is'
    call output_time_format(t-t0)
    write(*,'(A,$)') ' cpu  time is'
    call output_time_format(tc1-tc0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program Poisellie_field
