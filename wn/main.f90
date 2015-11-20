!DIR$ ATTRIBUTES FORCEINLINE :: FENE, LJ, BEND, inbox, scale_v
module parameters
    !use ifport
    implicit none
    ! 下标约定: p - polymer, s - solution, b - boundary

    !结构
    integer, parameter :: n_p=40, n_cell_x=30, n_cell_y=30, n_cell_z=20

    integer n_b, n_s

    !n是单体数目，n0是单体、壁、孔数目之和

    real(8), parameter :: kB = 1.38064852d-23
    real(8), parameter :: pi = 3.141592653589793238462643383279

    real(8), parameter :: time_step_p=1d-4, time_step_s=time_step_p, mass_p=1, mass_s=1, T_set=1, v_gradient=0.2

    ! polymer 位置 速度 力 上一次力
    real(8), dimension(3,n_p) :: x_p, v_p, f_p, f0_p
    ! solution
    real(8), allocatable, dimension(:,:) :: x_s, v_s, f_s, x0_s
    ! boundaaries, 1~nb-up, nb+1~2nb-down
    real(8), allocatable, dimension(:,:) :: x_b, v_b, f_b

    real(8) box_size_unit, half_box_size_unit, box_size(3), half_box_size(3)

    real(8) aver_v(3), distant

    real(8), parameter::sigma=1, epson=1, radius=4

contains

    function distance(x1,x2)
        implicit none
        real(8) distance,x1(3),x2(3)

        x1(3)=x1(3)-nint(x1(3)/n_cell_z)*n_cell_z
        x2(3)=x2(3)-nint(x2(3)/n_cell_z)*n_cell_z
        distance=norm2(x1-x2)
        return
    endfunction

    subroutine init()
        implicit none
        integer i,j,k,h,count_number
        real(8) distant,dx(2,4), temp(3,10000)
        logical success
        integer,parameter :: seed = 86456

        !!!!!!!!!以下是cylinder channel初值!!!!!!!!!!!!!!
        n_b=0
        k=nint(2.0*pi*radius)
        do j=0,2*n_cell_z
            do i=0,2*k-1
                n_b=n_b+1
                if(mod(j,2)==0)then
                    temp(1,n_b)=radius*cos(i*pi/k)
                    temp(2,n_b)=radius*sin(i*pi/k)
                else
                    temp(1,n_b)=radius*cos(i*pi/k+pi/(2d0*k))
                    temp(2,n_b)=radius*sin(i*pi/k+pi/(2d0*k))
                endif
                temp(3,n_b)=(j-n_cell_z)*0.5*box_size_unit
                ! write(output_file,'(2I6,3F13.4)') n_b,1,x_b(:,n_b)
            enddo
        enddo
        allocate(x_b(3,n_b),v_b(3,n_b),f_b(3,n_b))
        x_b=temp(:,1:n_b)
        write(*,*)"Cylinder particle number: ", n_b

        !!!!!!!!!以下是polymer chain初值!!!!!!!!!!!!!!
        x_p(1,1)=0d0
        x_p(2,1)=0d0
        x_p(3,1)=-5d0

        dx(1,1)=0
        dx(2,1)=1d0

        dx(1,2)=0
        dx(2,2)=-1d0

        dx(1,3)=1d0
        dx(2,3)=0

        dx(1,4)=-1d0
        dx(2,4)=0

        call srand(seed)
        do i=2,n_p
            success=.false.
            do count_number=1,11116

                h=1+int(4*rand())
                ! 当步数是5的倍数直接z轴加一
                if(mod(i,5)==0)then
                    x_p(1,i)=x_p(1,i-1)+0
                    x_p(2,i)=x_p(2,i-1)+0
                    x_p(3,i)=x_p(3,i-1)+1.0
                else
                    x_p(1,i)=x_p(1,i-1)+dx(1,h)
                    x_p(2,i)=x_p(2,i-1)+dx(2,h)
                    x_p(3,i)=x_p(3,i-1)
                endif

                ! 和壁保持距离
                if(sqrt(x_p(1,i)**2+x_p(2,i)**2)>radius-1.0) then
                    cycle
                endif

                ! 和除一级近邻外之前全部要保持一定距离
                success = .true.
                do k=1,i-2
                    if(distance(x_p(:,i),x_p(:,k))<1.1)then
                        success = .false.
                        exit
                    endif
                enddo

                if (success) then
                    exit
                endif
            enddo

            if (.not. success) then
                x_p(1,i)=x_p(1,i-1)+0
                x_p(2,i)=x_p(2,i-1)+0
                x_p(3,i)=x_p(3,i-1)+1.0
            endif
        enddo

        call periodic_p()

        write(*,*)"Polymer monomer number: ", n_p

        !!!!!!!!!!!!!!!!以下是solution粒子的初值!!!!!!!!!!!!!!

        n_s=0
        do i=0,n_cell_x
            do j=0,n_cell_y
                do k=0,n_cell_z
                    n_s=n_s+1
                    temp(1,n_s)=(i-n_cell_x/2d0)*box_size_unit
                    temp(2,n_s)=(j-n_cell_y/2d0)*box_size_unit
                    temp(3,n_s)=(k-n_cell_z/2d0)*box_size_unit
                    distant=sqrt(temp(1,n_s)**2+temp(2,n_s)**2)
                    if(distant>radius-1.0 .or. temp(3,n_s)>n_cell_z/2.0 .or. temp(3,n_s)<-n_cell_z/2.0)then
                        n_s=n_s-1
                        cycle
                    endif

                    !write(output_file,'(2I6,3F13.4)') n_b+n_p+n_s,3,x_s(:,n_s)
                enddo
            enddo
        enddo
        allocate(x_s(3,n_s),v_s(3,n_s), f_s(3,n_s), x0_s(3,n_s))
        x_s=temp(:,1:n_s)
        write(*,*)"Solvent particle number: ", n_s
        write(*,*)"Total particle number: ", n_b+n_p+n_s

    endsubroutine

    SUBROUTINE connect_z()
        INTEGER f,i,j
        REAL(8) distant_z,distant_max
        PARAMETER(distant_max=7)

        DO i=2,n_p
            distant_z=x_p(3,i)-x_p(3,i-1)
            IF(abs(distant_z)>distant_max)THEN
                IF(distant_z>0)THEN
                    f=-1
                ELSE
                    f=1
                ENDIF
                DO j=i,n_p
                    x_p(3,j)=x_p(3,j)+n_cell_z*f
                ENDDO
            ENDIF
        ENDDO
        RETURN
    END SUBROUTINE

    subroutine periodic_p()
        implicit none
        x_p(3,:)=x_p(3,:)-nint(x_p(3,:)/n_cell_z)*n_cell_z
    endsubroutine

    subroutine periodic_s()
        implicit none
        x_s(3,:)=x_s(3,:)-nint(x_s(3,:)/n_cell_z)*n_cell_z
    endsubroutine

    subroutine FENE(f,U,rx)
        implicit none
        real(8) f(3), U, rx(3)
        real(8), parameter :: FENE_rc=1.5*sigma, FENE_k=30
        real(8) temp,r

        r=norm2(rx)
        if(r/=0 .and. r<=FENE_rc) then
            temp=(FENE_k*FENE_rc**2)/(FENE_rc**2-r**2)
            f=f-rx*temp
            U=U-(FENE_k*FENE_rc**2)*log(1-(r/FENE_rc)**2)/2
        endif
    end subroutine

    subroutine LJ(f,U,rx)
        implicit none
        real(8) f(3), U, rx(3)
        real(8), parameter :: LJ_rc=sigma*2d0**(1d0/6d0)
        real(8) temp,r

        r=norm2(rx)
        if(r/=0 .and. r<=LJ_rc)then
            temp=48d0*epson*((sigma/r)**12)/(r**2)
            f=f+rx*temp
            U=U+4*epson*(sigma/r)**12
        endif
    end subroutine

    subroutine BEND(f,U,rx1,rx2)
        implicit none
        real(8) f(3), U, rx1(3), rx2(3)
        real(8), parameter :: BEND_b = 0

        f=f+BEND_b*(rx1-rx2)
        U=U-BEND_b*dot_product(rx1,rx2)
    end subroutine

    subroutine update_force(mode)
        !    use parameters, only : x_p, f_p, n_p,f0_p
        implicit none
        integer mode, i, j
        real(8) U, temp(3)
        call connect_z()
        f_p=0
        U=0
        temp=0
        ! 链两端
        call FENE(f_p(:,1),U,x_p(:,1)-x_p(:,2))
        call FENE(f_p(:,n_p),U,x_p(:,n_p)-x_p(:,n_p-1))
!write(*,*) 'force0',f_p(:,1)
        !$omp parallel do private(j,temp) reduction(+:U,f_p)
        do i=1,n_p
            if (i>1 .and. i<n_p) then
                ! UFENE(r) force
                call FENE(f_p(:,i),U,x_p(:,i)-x_p(:,j+1))
                call FENE(f_p(:,i),U,x_p(:,i)-x_p(:,j-1))
                ! bend energy
                call BEND(f_p(:,i),U,x_p(:,i+1)-x_p(:,i),x_p(:,i)-x_p(:,i-1))
            endif
            ! ULJ(r) force

            do j=1,n_p
                if (j/=i) then
                    call LJ(f_p(:,i),U,x_p(:,i)-x_p(:,j))
                endif
            enddo

!            temp(2)=x_p(2,i)
!            call LJ(f_p(:,i),U,temp)
!            temp(2)=n_cell_y-x_p(2,i)
!            call LJ(f_p(:,i),U,temp)

        enddo
        !$omp end parallel do
        !write(*,*) 'force1',f_p(:,1)
        if (mode==0) then
            f0_p=f_p
        endif

        call periodic_p()

    end subroutine

    logical function inbox(x,ix,iy,iz)
        implicit none
        real(8) x(3), temp(3)
        integer ix,iy,iz

        temp=x/box_size_unit-[ix,iy,iz]*1d0
        inbox = temp(1)>=0 .and. temp(1)<1 &
            .and. temp(2)>=0 .and. temp(2)<1 &
            .and. temp(3)>=0 .and. temp(3)<1
    end function

    subroutine cal_collision_velocity()
        !    use parameters
        implicit none
        integer ix,iy,iz,k,count_p,count_s,i
        real(8) momentum(3), matrix(3,3), l(3), fai, theta
        real(8), parameter :: alpha = 130*pi/180, s=sin(alpha), c=cos(alpha)
        real(8) v_aver_p(3), v_aver_s(3), v_aver(3), temp(3)
        logical mask_p(n_p), mask_s(n_s)
        ! calculate velocity of all particles in each cell
        k=0
        do ix=0,n_cell_x-1
            do iy=0,n_cell_y-1
                do iz=0,n_cell_x-1

                    k=k+1

                    fai=2.0*pi*rand(0)
                    theta=2.0*rand(0)-1

                    l(1)=cos(fai)*SQRT(1-theta**2)
                    l(2)=sin(fai)*SQRT(1-theta**2)
                    l(3)=theta

                    matrix(1,1) = l(1)*l(1)*(1-c) + c
                    matrix(1,2) = l(1)*l(2)*(1-c) - s*l(3)
                    matrix(1,3) = l(1)*l(3)*(1-c) + s*l(2)

                    matrix(2,1) = l(2)*l(1)*(1-c) + s*l(3)
                    matrix(2,2) = l(2)*l(2)*(1-c) + c
                    matrix(2,3) = l(2)*l(3)*(1-c) - s*l(1)

                    matrix(3,1) = l(3)*l(1)*(1-c) - s*l(2)
                    matrix(3,2) = l(3)*l(2)*(1-c) + s*l(1)
                    matrix(3,3) = l(3)*l(3)*(1-c) + c

                    momentum=0

                    count_p=0
                    v_aver_p=0d0
                    mask_p=.false.
                    do i=1,n_p
                        if(inbox(x_p(:,i),ix,iy,iz))then
                            momentum = momentum + v_p(:,i)*mass_p
                            count_p=count_p+1
                            mask_p(i)=.true.
                            v_aver_p=v_aver_p+v_p(:,i)
                        endif
                    enddo
!write(*,*)ix,iy,iz,count_p
                    count_s=0
                    v_aver_s=0d0
                    mask_s=.false.
                    do i=1,n_s
                        if(inbox(x_s(:,i),ix,iy,iz))then
                            momentum = momentum + v_s(:,i)*mass_s
                            count_s=count_s+1
                            mask_s(i)=.true.
                            v_aver_s=v_aver_s+v_s(:,i)
                        endif
                    enddo
!write(*,*)ix,iy,iz,count_s
                    if((count_s==0).and.(count_p==0))then
                        v_aver=0
                    else
                        v_aver=momentum/(count_p*mass_p+count_s*mass_s)
                    endif

                    do i=1,n_p
                        if(mask_p(i)) then
                            temp=v_p(:,i)-v_aver
                            v_p(1,i)=v_aver(1) + matrix(1,1)*temp(1) + matrix(1,2)*temp(2) + matrix(1,3)*temp(3)
                            v_p(2,i)=v_aver(2) + matrix(2,1)*temp(1) + matrix(2,2)*temp(2) + matrix(2,3)*temp(3)
                            v_p(3,i)=v_aver(3) + matrix(3,1)*temp(1) + matrix(3,2)*temp(2) + matrix(3,3)*temp(3)
                        endif
                    enddo

                    do i=1,n_s
                        if(mask_s(i)) then
                            temp=v_s(:,i)-v_aver
                            v_s(1,i)=v_aver(1) + matrix(1,1)*temp(1) + matrix(1,2)*temp(2) + matrix(1,3)*temp(3)
                            v_s(2,i)=v_aver(2) + matrix(2,1)*temp(1) + matrix(2,2)*temp(2) + matrix(2,3)*temp(3)
                            v_s(3,i)=v_aver(3) + matrix(3,1)*temp(1) + matrix(3,2)*temp(2) + matrix(3,3)*temp(3)
                        endif
                    enddo

                enddo
            enddo
        enddo
    end subroutine

    !!!!以下是Isokinetics thermostat，可以试一下Berendsen thermostat
    subroutine scale_v(Ek,T,T_out)
        implicit none
        integer i
        real(8), parameter :: Ek_fac = 1.5d0
        real(8) v, Ek,T, scalar, Ek1, T_out, T1


        do i=1,3
            v=sum(v_p(i,:)*mass_p+v_s(i,:)*mass_s)/(n_p*mass_p+n_s*mass_s)
            v_p(i,:) = v_p(i,:)-v
            v_s(i,:) = v_s(i,:)-v
        enddo

        Ek1=0.5*(mass_p*sum(v_p**2)+mass_s*sum(v_s**2))
        T1=Ek1/(Ek_fac*(n_p+n_s))
        scalar=sqrt(T/T1)
        v_p=v_p*scalar
        v_s=v_s*scalar
        Ek=0.5*(mass_p*sum(v_p**2)+mass_s*sum(v_s**2))
        T_out=Ek/(Ek_fac*(n_p+n_s))

    end subroutine

    subroutine output(output_file,cur_step,step)
        implicit none
        integer cur_step,step,output_file,k

        if(mod(cur_step,step)==0)then
        write(output_file,*)'ITEM:TIMESTEP'
        write(output_file,'(I9)')cur_step
        write(output_file,*)'ITEM:NUMBER OF ATOMS'
        write(output_file,'(I6)')n_b+n_p+n_s
        write(output_file,*)'ITEM:BOX BOUNDS'
        write(output_file,'(2F7.1)')-radius,radius
        write(output_file,'(2F7.1)')-radius,radius
        write(output_file,'(2F7.1)')-n_cell_z/2.0,n_cell_z/2.0
        write(output_file,*)'ITEM:ATOMS id type x y z'
        do k=1,n_b
            write(output_file,'(2I6,3F13.4)') k,1,x_b(:,k)
        enddo
        do k=1,n_p
            write(output_file,'(2I6,3F13.4)') n_b+k,2,x_p(:,k)
        enddo
        do k=1,n_s
            write(output_file,'(2I6,3F13.4)') n_b+n_p+k,3,x_s(:,k)
        enddo
        endif

    end subroutine

    subroutine output_U(energy_file,cur_step,step)
        implicit none
        integer cur_step, step, energy_file
        real(8) U
        energy_file=13
        if(mod(cur_step,step)==0)then
        write(energy_file,*) cur_step,U
        endif
    end subroutine

end module


program Poissonfield
    use parameters
    implicit none
    real(8) :: Ek, EK_scaled, T_scaled, density
    integer :: equili_step,equili_interval_step,total_step,output_interval_step,cur_step_per_rot,total_step_per_rot, &
        cur_step,total_rot_step,output_file

   ! real(8), allocatable :: randnum_group(:,:)

    integer i, j
    real(8) randx, randz, ppx, ppz

    output_file=12

    equili_step=500000
    equili_interval_step=10
    total_step=3000000
    output_interval_step=10000
    total_rot_step=500000
    total_step_per_rot=200

    box_size = [n_cell_x, n_cell_y, n_cell_z]
    half_box_size(3) = n_cell_z/2d0
    density=3.0
    box_size_unit=(1/density)**(1d0/3)
    half_box_size_unit=box_size_unit/2

    !!!读链的大小 改成1个文件
    open(output_file,file='dump.cylinder.lammpstrj')
    call init()
    call output(output_file,0,equili_interval_step)

!    allocate(randnum_group(3,n_s))
!    call random_number(randnum_group)
!    x_s(:,:) = x_s(:,:) + (randnum_group(:,:)-0.5)*box_size_unit
!    deallocate(randnum_group)

!    do i=1,n_s
!        x_s(3,i) = x_s(3,i) - box_size(3)*nint((x_s(3,i)-half_box_size(3))/box_size(3))
!    enddo

     !!!polymer的初速度
    call random_number(v_p)
    !!!solution的初速度
    call random_number(v_s)
    v_p=v_p-0.5
    v_s=v_s-0.5
    call scale_v(Ek_scaled, T_set, T_scaled)
        write(*,*) 'Polymer: '
        write(*,*) 'Initial scaled kinetics Ek_scaled: ',Ek_scaled
        write(*,*) 'Initial setted temperature T_set: ',T_set
        write(*,*) 'Initial scaled temperature T_scaled: ',T_scaled

        write(*,*) 'Solution: '
        write(*,*) 'Initial scaled kinetics Ek_scaled=',Ek_scaled
        write(*,*) 'Initial setted temperature T_set=',T_set
        write(*,*) 'Initial scaled temperature T_scaled=',T_scaled

    ! 没有外场时，polymer和solution达到平衡
        call update_force(0)

    do cur_step=1,equili_step
                if(mod(cur_step,10000)==0)then
            !write(*,*)v_s(:,1)
        endif
        x_p = x_p + v_p*time_step_p + 0.5*f0_p*time_step_p**2
        x_s = x_s + v_s*time_step_s
        call update_force(1)
        !write(*,*) v_p(:,2)
        v_p = v_p + 0.5*(f0_p+f_p)*time_step_p
        !write(*,*) v_p(:,2),f0_p(:,2),f_p(:,2)
        call scale_v(Ek,T_set,T_scaled)
        call cal_collision_velocity()
        call scale_v(Ek,T_set,T_scaled)
        f0_p=f_p
        call output(output_file,cur_step,equili_interval_step)
        if(mod(cur_step,10)==0)then
            !write(*,*)v_s(:,1)
        endif
    enddo

    call cal_collision_velocity()

    write(*,*)'循环开始'
    !!! compute a(t-dt)
    call update_force(0)
    !    do i=1,256
    !        x_up1(i)=x_bu(i)
    !        z_up1(i)=z_bu(i)
    !        x_down1(i)=x_bd(i)
    !        z_down1(i)=z_bd(i)
    !    enddo

    x0_s=x_s

    do cur_step=1,total_step

        do cur_step_per_rot=1,total_step_per_rot

            !!!!!! compute a(t+dt)
            x_p = x_p + v_p*time_step_p + 0.5*f0_p*time_step_p**2
            x_p(1,:) = x_p(1,:)-0.1*(x_p(2,:)**2-15.0*x_p(2,:))*v_gradient*time_step_p
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call update_force(1)

            v_p = v_p + 0.5*(f0_p+f_p)*time_step_p
            call scale_v(Ek,T_set,T_scaled)
            f0_p=f_p
            !   call output(cur_step_pri)

        enddo        !!!!st3

        randx=rand(0)
        randz=rand(0)
        x_p(1,:)=x_p(1,:)+(randx-0.5)*box_size_unit
        x_p(3,:)=x_p(3,:)+(randz-0.5)*box_size_unit

        ppx=0.05
        ppz=0.05

        do j=1,n_p

            if(box_size(1)*nint((x_p(1,j)-half_box_size(1))/box_size(1))/=ppx &
                    .or.box_size(3)*nint((x_p(3,j)-half_box_size(3))/box_size(3))/=ppz)then
                !                do i=1,256
                !                    x_bu(i)=x_up1(i)
                !                    z_bu(i)=z_up1(i)
                !                    x_bd(i)=x_down1(i)
                !                    z_bd(i)=z_down1(i)
                !                enddo

                x_s=x0_s

                ppx=box_size(1)*nint((x_p(1,j)-half_box_size(1))/box_size(1))
                ppz=box_size(3)*nint((x_p(3,j)-half_box_size(3))/box_size(3))

                !                do i=1,256
                !                    x_bu(i)=x_bu(i)+ppx
                !                    z_bu(i)=z_bu(i)+ppz
                !                    x_bd(i)=x_bd(i)+ppx
                !                    z_bd(i)=z_bd(i)+ppz
                !                enddo

                x_s(1,:)=x_s(1,:)+ppx
                x_s(3,:)=x_s(3,:)+ppz

                x_s = x_s + v_s*time_step_s
                x_s(1,:) = x_s(1,:) -0.1*(x_s(2,:)**2-15.0*x_s(2,:))*v_gradient*time_step_s
                ! 弹回, 看起来不对
                do i=1,n_s
                    !y12(i)=y_s(i)
                    !if(y_s(i)>=box_size_y.or.y_s(i)<=0)then
                    !y_s(i)=y12(i)
                    !endif
                enddo

                x_s(1,:)=x_s(1,:)+(randx-0.5)*box_size_unit
                x_s(3,:)=x_s(3,:)+(randz-0.5)*box_size_unit

                x_s(1,:)=x_s(1,:)-box_size(1)*nint((x_s(1,:)-(ppx+half_box_size(1)))/box_size(1))
                x_s(3,:)=x_s(3,:)-box_size(3)*nint((x_s(3,:)-(ppz+half_box_size(3)))/box_size(3))

                !if(cur_step>=st0.and.cur_step<st00.and.MOD(cur_step,st22)==0)then
                !    write(100,*)cur_step,j
                !endif

                call cal_collision_velocity()

                call scale_v(Ek,T_set,T_scaled)

            endif

        enddo



    enddo        !!!!!step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    close(70)
    close(80)
    close(81)
    close(821)
    close(82)
    close(83)
    close(84)
    close(88)
    close(90)
    close(91)
    close(100)
    close(110)
    close(111)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end program Poissonfield
