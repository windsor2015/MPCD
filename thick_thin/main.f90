#ifdef __INTEL_COMPILER
include 'mkl_vsl.f90'
#endif

module parameters
#ifdef __INTEL_COMPILER
    use ifport
    USE MKL_VSL_TYPE
    USE MKL_VSL
#endif
    implicit none
    ! 下标约定: p - polymer, s - solution, b - boundary / phantom

    !结构
    integer, parameter :: n_p=0, n_cell_x=12, n_cell_y=12, n_cell_z=40

    integer n_b(3), n_s(2)

    !n是单体数目，n0是单体、壁、孔数目之和

    real(8), parameter :: kB = 1.38064852d-23
    real(8), parameter :: pi = 3.141592653589793238462643383279

    real(8), parameter :: Ek_fac = 1.5d0

    real(8), parameter :: time_step_p=1d-4, time_step_s=5d-3, mass_p=1, mass_s=0.2, T_set=1, v_gradient=0.2

    ! polymer 位置 速度 力 上一次力
    real(8), dimension(3,n_p) :: x_p, v_p, f_p, f0_p, x0_p
    ! solution
    real(8), allocatable, dimension(:,:) :: x_s, v_s, f_s, x0_s, x_s0
    ! boundaaries, 1~nb-up, nb+1~2nb-down
    real(8), allocatable, dimension(:,:) :: x_b, v_b, f_b

    real(8) box_size_unit, half_box_size_unit, box_size(3), half_box_size(3)

    real(8) aver_v(3), distant

    real(8), parameter :: sigma=1, epson=1

    real(8) :: density_s=5, gama=0.01

    real(8), dimension(2)::radius=(/4, 2/)

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

    integer :: thermostat_method, thermostat_interval
    real(8) thermostat_B_parameter, thermostat_A_parameter

    namelist /basic/ radius, gama, density_s, &
        desk_interval_step, equili_step, equili_interval_step, total_step, output_interval_step, &
        thermostat_method, thermostat_interval, thermostat_B_parameter,thermostat_A_parameter

#ifdef __INTEL_COMPILER
    TYPE (VSL_STREAM_STATE) :: vsl_stream
#endif

contains

    function distance(x1,x2)
        implicit none
        real(8) distance,x1(3),x2(3)

        x1(3)=x1(3)-nint(x1(3)/n_cell_z)*n_cell_z
        x2(3)=x2(3)-nint(x2(3)/n_cell_z)*n_cell_z
        distance=norm2(x1-x2)
        return
    endfunction

    subroutine readin()
        implicit none
        integer input_file
        input_file=9000
        open(input_file, file='input_file')
        read(input_file, nml=basic)
        close(input_file)
    endsubroutine

    subroutine init()
        implicit none
        integer i,j,k,count_number
        real(8) dx(2,4), theta, r
        logical success
        integer,parameter :: seed = 11111

        !!!!!!!!!以下是polymer chain初值!!!!!!!!!!!!!!
        x_p(1,1)=0d0
        x_p(2,1)=0d0
        x_p(3,1)=-15d0

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

                j=1+int(4*rand())
                ! 当步数是5的倍数直接z轴加一
                if(mod(i,5)==0)then
                    x_p(1,i)=x_p(1,i-1)+0
                    x_p(2,i)=x_p(2,i-1)+0
                    x_p(3,i)=x_p(3,i-1)+1.0
                else
                    x_p(1,i)=x_p(1,i-1)+dx(1,j)
                    x_p(2,i)=x_p(2,i-1)+dx(2,j)
                    x_p(3,i)=x_p(3,i-1)
                endif

                ! 和壁保持距离
                if((sqrt(x_p(1,i)**2+x_p(2,i)**2)>radius(1).and.x_p(3,i)<-5d0).or.(sqrt(x_p(1,i)**2+x_p(2,i)**2)>radius(2).and.x_p(3,i)>-5d0.and.x_p(3,i)<5d0)&
                                                 .or.(sqrt(x_p(1,i)**2+x_p(2,i)**2)>radius(1).and.x_p(3,i)>5d0)) then
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

        n_s(1)=nint(pi*radius_thick**2*3d0*n_cell_z/8d0*density_s)
        n_s(2)=n_s(1)+nint(pi*radius_thin**2*n_cell_z/4d0*density_s)
        n_s(3)=n_s(1)+n_s(2)
        allocate(x_s(3,n_s(3)),v_s(3,n_s(3)),f_s(3,n_s(3)),x0_s(3,n_s(3)),x_s0(3,n_s(3)))
        i=0
            do while(.true.)
                i=i+1
                if(i>n_s(1))exit
                theta=pi*2*rand()
                r=sqrt(rand())*radius(1)
                x_s(1,i)=r*cos(theta)
                x_s(2,i)=r*sin(theta)
                x_s(3,i)=(rand()-1d0)*n_cell_z/2d0
                if(x_s(3,i)>-5d0)then
                    i=i-1
                    cycle
                end if
            enddo
            i=i-1

            do while(.true.)
                i=i+1
                if(i>n_s(2))exit
                theta=pi*2*rand()
                r=sqrt(rand())*radius(2)
                x_s(1,i)=r*cos(theta)
                x_s(2,i)=r*sin(theta)
                x_s(3,i)=(rand()-5d-1)*n_cell_z/4d0
            enddo
            i=i-1

            do while(.true.)
                i=i+1
                if(i>n_s(3))exit
                theta=pi*2*rand()
                r=sqrt(rand())*radius(1)
                x_s(1,i)=r*cos(theta)
                x_s(2,i)=r*sin(theta)
                x_s(3,i)=rand()*n_cell_z/2d0
                if(x_s(3,i)<5d0)then
                    i=i-1
                    cycle
                end if
            enddo
            i=i-1
        !!!!!!!!!!phantom particles
        n_b(1)=nint(pi*((radius(1)+sqrt(5d-1))**2-radius(1)**2)*n_cell_z/2d0*density_s)
        n_b(2)=n_b(1)+nint(pi*((radius(2)+sqrt(5d-1))**2-radius(2)**2)*n_cell_z/2d0*density_s)
        n_b(3)=n_b(2)+nint(pi*((radius(1)+sqrt(5d-1))**2-(radius(2)+sqrt(5d-1))**2)*sqrt(5d-1)*density_s)

        allocate(x_b(3,n_b(5)),v_b(3,n_b(5)))
        i=0
            do while(.true.)
                i=i+1
                if(i>n_b(1))exit
                theta=pi*2*rand()
                r=sqrt(rand())*(radius(1)+sqrt(5d-1))
                if(r<=radius(1))then
                    i=i-1
                else
                    x_b(1,i)=r*cos(theta)
                    x_b(2,i)=r*sin(theta)
                    x_b(3,i)=(rand()-1d0)*n_cell_z/2d0
                if(x_b(3,i)>-5d0)then
                    i=i-1
                    cycle
                end if
                end if
            enddo
            i=i-1

        do while(.true.)
            i=i+1
            if(i>n_b(2))exit
            theta=pi*2*rand()
            r=sqrt(rand())*radius(1)
            if(r<=radius(2))then
                i=i-1
            else
                x_b(1,i)=r*cos(theta)
                x_b(2,i)=r*sin(theta)
                x_b(3,i)=rand()*(-5d0+sqrt(5d-1))
            end if
        enddo
        i=i-1

            do while(.true.)
                i=i+1
                if(i>n_b(3))exit
                theta=pi*2*rand()
                r=sqrt(rand())*(radius(2)+sqrt(5d-1))
                if(r<=radius(2))then
                    i=i-1
                else
                    x_b(1,i)=r*cos(theta)
                    x_b(2,i)=r*sin(theta)
                    x_b(3,i)=(rand()-5d-1)*n_cell_z/4d0
                if(x_b(3,i)>5d0-sqrt(5d-1).or.x_b(3,i)<-5d0+sqrt(5d-1))then
                    i=i-1
                    cycle
                end if
                end if
            enddo
            i=i-1
        write(*,*)i-1
        write(*,*)"Solvent particle number: ", n_s(2)
        write(*,*)"Total particle number: ", n_p+n_s(2)
        write(*,*)"Phantom particle number: ", n_b(3)

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
        x_p(3,:)=x_p(3,:)-n_cell_z*nint((x_p(3,:))/n_cell_z)
    endsubroutine

    subroutine periodic_s()
        implicit none
        x_s(3,:)=x_s(3,:)-n_cell_z*nint((x_s(3,:))/n_cell_z)
    endsubroutine

    subroutine FENE(f,U,rx)
        implicit none
        real(8) f(3), U, rx(3)
        real(8), parameter :: FENE_rc=1.5*sigma, FENE_k=20
        real(8) temp,r

        rx(3)=rx(3)-n_cell_z*nint(rx(3)/n_cell_z)
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

        rx(3)=rx(3)-n_cell_z*nint(rx(3)/n_cell_z)
        r=norm2(rx)
        if(r/=0 .and. r<=LJ_rc)then
            temp=24d0*epson*(2d0*(sigma/r)**12-(sigma/r)**6)/(r**2)
            f=f+rx*temp
            U=U+4*epson*((sigma/r)**12-(sigma/r)**6)
        endif
    end subroutine

    subroutine BEND(f,U,rx1,rx2)
        implicit none
        real(8) f(3), U, rx1(3), rx2(3), c, r1, r2
        real(8), parameter :: BEND_b = 10

        rx1(3)=rx1(3)-n_cell_z*nint(rx1(3)/n_cell_z)
        rx2(3)=rx2(3)-n_cell_z*nint(rx2(3)/n_cell_z)

        r1=norm2(rx1)
        r2=norm2(rx2)
        c=dot_product(rx1,rx2)/r1/r2
        f=f-BEND_b*((rx1+rx2)/(r1*r2)-c*rx1/r1/r1-c*rx2/r2/r2)
        U=U+BEND_b*(1+c)

        !c=dot_product(rx1,rx2)
        !f=f+BEND_b*(rx1-rx2)
        !U=U-BEND_b*c
    end subroutine

    subroutine update_force(mode)
        !    use parameters, only : x_p, f_p, n_p,f0_p
        implicit none
        integer mode, i, j
        real(8) temp(3)

        f_p=0
        U=0
        U_BEND=0
        U_FENE=0
        U_LJ=0
        U_WALL=0
        temp=0
        !call connect_z()
        ! 链两端
        call FENE(f_p(:,1),U_FENE,x_p(:,1)-x_p(:,2))
        call FENE(f_p(:,n_p),U_FENE,x_p(:,n_p)-x_p(:,n_p-1))
        !if (u>10000) write(*,*) 'force0',U

        !$omp parallel do private(j) reduction(+:U_FENE,U_LJ,U_BEND,f_p)
        do i=1,n_p
            if (i>1 .and. i<n_p) then
                ! UFENE(r) force
                call FENE(f_p(:,i),U_FENE,x_p(:,i)-x_p(:,i+1))
                call FENE(f_p(:,i),U_FENE,x_p(:,i)-x_p(:,i-1))

                ! bend energy
                call BEND(f_p(:,i),U_BEND,-x_p(:,i+1)+x_p(:,i),x_p(:,i)-x_p(:,i-1))
            endif

            do j=1,n_p
                if (j/=i) then
                    call LJ(f_p(:,i),U_LJ,x_p(:,i)-x_p(:,j))
                endif
            enddo

        enddo
        !$omp end parallel do
        !if (u>10000) write(*,*) 'force1',U
        !call periodic_p()

        !        !$omp parallel do private(j) reduction(+:U_LJ,f_p,U_WALL)
        !        !ULJ(r) force
        !        do i=1,n_p
        !            do j=1,n_p
        !                if (j/=i) then
        !                    call LJ(f_p(:,i),U_LJ,x_p(:,i)-x_p(:,j))
        !                endif
        !            enddo
        !
        !!            do j=1,n_b
        !!                call LJ(f_p(:,i),U_WALL,x_p(:,i)-x_b(:,j))
        !!            enddo
        !        enddo
        !        !$omp end parallel do
        U=U_FENE+U_BEND+U_LJ+U_WALL
        !if (u>10000) write(*,*) 'force2',U
        if (mode==0) then
            f0_p=f_p
        endif

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

    subroutine get_cell_xyz(r,x,y,z,s)
        implicit none
        real(8) r(3)
        integer x,y,z
        logical s

        x=floor(r(1)+n_cell_x/2d0)
        y=floor(r(2)+n_cell_y/2d0)
        z=floor(r(3)-n_cell_z*nint(r(3)/n_cell_z)+n_cell_z/2d0)

        s= x>=0 .and. x<=n_cell_x &
            .and. y>=0 .and. y<=n_cell_y &
            .and. z>=0 .and. z<=n_cell_z

        !if (.not.s) write(*,*) r

    end subroutine

    subroutine cal_collision_velocity(cur_step)
        !    use parameters
        implicit none
        integer ix,iy,iz,k,count_p,count_s,count_b,i,pointer_p(n_p),pointer_s(1000),cur_step
        real(8)  matrix(3,3), l(3), fai, theta
        real(8), parameter :: alpha = 130*pi/180, s=sin(alpha), c=cos(alpha)
        real(8)  v_aver(3), temp(3), randx, randy, randz, Ek, scalar, randr(3), sigma1
        logical  check

        pointer_cell_p=0
        pointer_cell_s=0
        count_cell_p=0
        count_cell_s=0
        count_cell_b=0
        momentum_cell=0d0
        scalar=1d0
        ! calculate velocity of all particles in each cell
        !randx=(rand()-0.5)*box_size_unit
        !randy=(rand()-0.5)*box_size_unit
        !randz=(rand()-0.5)*box_size_unit

        call random_number(randr)
        randr=(randr-0.5)*box_size_unit

        !        x0_p(1,:)=x_p(1,:)+randx
        !        x0_s(1,:)=x_s(1,:)+randx
        !        x0_p(2,:)=x_p(2,:)+randy
        !        x0_s(2,:)=x_s(2,:)+randy
        !        x0_p(3,:)=x_p(3,:)+randz
        !        x0_s(3,:)=x_s(3,:)+randz
        !        x0_p(3,:)=x0_p(3,:)-n_cell_z*nint((x0_p(3,:))/n_cell_z)
        !        x0_s(3,:)=x0_s(3,:)-n_cell_z*nint((x0_s(3,:))/n_cell_z)

        !$omp parallel do private(ix,iy,iz,check)
        do i=1,n_p
            call get_cell_xyz(x_p(:,i)+randr,ix,iy,iz,check)
            if (.not.check) cycle
            count_cell_p(ix,iy,iz)=count_cell_p(ix,iy,iz)+1
            pointer_cell_p(ix,iy,iz,count_cell_p(ix,iy,iz))=i
            momentum_cell(:,ix,iy,iz)=momentum_cell(:,ix,iy,iz)+mass_p*v_p(:,i)
        enddo
        !$omp end parallel do

        !$omp parallel do private(ix,iy,iz,check)
        do i=1,n_s
            call get_cell_xyz(x_s(:,i)+randr,ix,iy,iz,check)
            if (check) then
                count_cell_s(ix,iy,iz)=count_cell_s(ix,iy,iz)+1
                pointer_cell_s(ix,iy,iz,count_cell_s(ix,iy,iz))=i
                momentum_cell(:,ix,iy,iz)=momentum_cell(:,ix,iy,iz)+mass_s*v_s(:,i)
            endif
        enddo
        !$omp end parallel do
        !write(*,*) maxval(count_cell_s)

        ! phantom particle, the velocity
        sigma1=sqrt(T_set)
        i=vdrnggaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, vsl_stream, n_b*3, v_b,0d0,sigma1)
        !$omp parallel do private(ix,iy,iz,check,temp)
        do i=1,n_b
            call get_cell_xyz(x_b(:,i)+randr,ix,iy,iz,check)
            if (check) then
                count_cell_b(ix,iy,iz)=count_cell_b(ix,iy,iz)+1
                pointer_cell_b(ix,iy,iz,count_cell_b(ix,iy,iz))=i
                !temp(1)=rand_gaussian(sigma1)
                !temp(2)=rand_gaussian(sigma1)
                !temp(3)=rand_gaussian(sigma1)
                !v_b(:,i)=temp/mass_s
                ! gaussian should be used
                momentum_cell(:,ix,iy,iz)=momentum_cell(:,ix,iy,iz)+mass_s*v_b(:,i)
            endif
        enddo
        !$omp end parallel do

        !$omp parallel do private(iy,iz,i,matrix,l,fai,theta,v_aver,temp,k,scalar,count_p,count_s,count_b,Ek)
        do ix=0,n_cell_x
            do iy=0,n_cell_y
                do iz=0,n_cell_z

                    if (count_cell_p(ix,iy,iz)+count_cell_s(ix,iy,iz)<=1) cycle

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

                    v_aver=momentum_cell(:,ix,iy,iz) &
                        / (mass_p*count_cell_p(ix,iy,iz)+mass_s*count_cell_s(ix,iy,iz)+mass_s*count_cell_b(ix,iy,iz))
                    count_p=count_cell_p(ix,iy,iz)
                    count_s=count_cell_s(ix,iy,iz)
                    count_b=count_cell_b(ix,iy,iz)
                    ! need two loops, 1 calculate Ek with delta v, 2 rotate and scale
                    Ek=0
                    do i=1,count_p
                        k=pointer_cell_p(ix,iy,iz,i)
                        !pointer_p(i)=pointer_cell_p(ix,iy,iz,i)
                        v_p(:,k)=v_p(:,k)-v_aver
                        Ek=Ek+mass_p*sum(v_p(:,k)**2)
                    enddo
                    do i=1,count_s
                        k=pointer_cell_s(ix,iy,iz,i)
                        !pointer_s(i)=pointer_cell_s(ix,iy,iz,i)
                        v_s(:,k)=v_s(:,k)-v_aver
                        Ek=Ek+mass_s*sum(v_s(:,k)**2)
                    enddo
                    do i=1,count_b
                        k=pointer_cell_b(ix,iy,iz,i)
                        !pointer_s(i)=pointer_cell_s(ix,iy,iz,i)
                        v_b(:,k)=v_b(:,k)-v_aver
                        Ek=Ek+mass_s*sum(v_b(:,k)**2)
                    enddo
                    Ek=Ek/2d0

                    scalar=1d0
                    if (thermostat_method>=10) then
                        scalar=thermostat_cal_scalar(cur_step, count_cell_p(ix,iy,iz)+count_cell_s(ix,iy,iz)+count_cell_b(ix,iy,iz), Ek)
                    end if
                    !write(*,*)scalar

                    do i=1,count_p
                        k=pointer_cell_p(ix,iy,iz,i)
                        temp=v_p(:,k)*scalar
                        v_p(1,k)= v_aver(1) + matrix(1,1)*temp(1) + matrix(1,2)*temp(2) + matrix(1,3)*temp(3)
                        v_p(2,k)= v_aver(2) + matrix(2,1)*temp(1) + matrix(2,2)*temp(2) + matrix(2,3)*temp(3)
                        v_p(3,k)= v_aver(3) + matrix(3,1)*temp(1) + matrix(3,2)*temp(2) + matrix(3,3)*temp(3)
                    enddo
                    do i=1,count_s
                        k=pointer_cell_s(ix,iy,iz,i)
                        temp=v_s(:,k)*scalar
                        v_s(1,k)= v_aver(1) + matrix(1,1)*temp(1) + matrix(1,2)*temp(2) + matrix(1,3)*temp(3)
                        v_s(2,k)= v_aver(2) + matrix(2,1)*temp(1) + matrix(2,2)*temp(2) + matrix(2,3)*temp(3)
                        v_s(3,k)= v_aver(3) + matrix(3,1)*temp(1) + matrix(3,2)*temp(2) + matrix(3,3)*temp(3)
                    enddo

                enddo
            enddo
        enddo
        !$omp end parallel do

    end subroutine

    real(8) function thermostat_cal_scalar(cur_step, n_c, Ek) result(scalar)
        implicit none
        integer cur_step, n_c
        real(8) Ek, T1
        scalar=1d0
        if (mod(cur_step, thermostat_interval)/=0) return

        select case (thermostat_method)
            case(10)
                scalar=sqrt(1.5d0*(n_c-1)*T_set/Ek)
            case(11)
                call thermostat_MCS(Ek, n_c, scalar)
            case(12)
                call thermostat_MBS(Ek, n_c, scalar)
        end select
    end function

    subroutine thermostat_init()
        implicit none
        integer s
#ifdef __INTEL_COMPILER
        s=vslnewstream(vsl_stream,VSL_BRNG_MT19937,77777)
#endif
        if (thermostat_method/=0) thermostat_interval=1

    end subroutine

    subroutine Ek_T(Ek,T_out)
        implicit none
        real(8) Ek, T_out

        Ek=0.5*(mass_p*sum(v_p**2)+mass_s*sum(v_s**2))
        T_out=Ek/(Ek_fac*(n_p+n_s-1))
    end subroutine

    subroutine thermostat(cur_step)
        implicit none
        integer cur_step

        if (mod(cur_step, thermostat_interval)/=0) return

        select case (thermostat_method)
            case(0)
                call thermostat_I()
            case(1)
                call thermostat_B()
            case(2)
                call thermostat_A()
            case(3)
                !call thermostat_MBS(count_p,pointer_p,count_s,pointer_s)
            case(4)
                !call thermostat_MCS(count_p,pointer_p,count_s,pointer_s)
            case(5)
                !call thermostat_L(Ek,T,T_out)
        end select

    end subroutine

    !!!!以下是Isokinetics thermostat
    subroutine thermostat_I()
        implicit none
        real(8) Ek1, T1, scalar
        Ek1=0.5*(mass_p*sum(v_p**2)+mass_s*sum(v_s**2))
        T1=Ek1/(Ek_fac*(n_p+n_s))
        scalar=sqrt(T_set/T1)
        v_p=v_p*scalar
        v_s=v_s*scalar
    end subroutine

    !Berendsen thermostat
    subroutine thermostat_B()
        implicit none
        real(8) v, Ek,T, Ek1, T_out, T1,lambda,tau,scalar
        scalar=thermostat_B_parameter
        Ek1=0.5*(mass_p*sum(v_p**2)+mass_s*sum(v_s**2))
        T1=Ek1/(Ek_fac*(n_p+n_s))
        lambda = sqrt(1+scalar*(T_set/T1-1.0))
        v_p=v_p*lambda
        v_s=v_s*lambda
    end subroutine

    function rand_gaussian(sigma) result(r)
        implicit none
        real(8) sigma, r, x, y, r2

        do while(.true.)
            x=-1+2*rand()
            y=-1+2*rand()
            r2=x**2+y**2
            if (r2<=1d0 .and. r2/=0d0) exit
        enddo
        r=sigma*y*sqrt(-2d0*log(r2)/r2)

    end function


    subroutine thermostat_A()
        implicit none
        real(8) nu_plus_dt,sigma1,r(6)
        integer i,k,s
        nu_plus_dt=thermostat_A_parameter
        sigma1=sqrt(m_s*T_set)    !!？？？？？？？？？？？？？？？？？？？可以吗？
        #ifdef __INTEL_COMPILER
        s=vdrnggaussian(vsl_rng_method_gaussian_boxmuller,vsl_stream,6,r,0d0,sigma1)
        #endif

        do i=1,n_p
            if (rand()<nu_plus_dt) then
                v_p(1,i)=r(1)/mass_p
                v_p(2,i)=r(2)/mass_p
                v_p(3,i)=r(3)/mass_p
            endif
        enddo
        do i=1,n_s
            if (rand()<nu_plus_dt) then
                v_s(1,i)=r(4)/mass_s
                v_s(2,i)=r(5)/mass_s
                v_s(3,i)=r(6)/mass_s
            endif
        enddo

        !        do i=1,n_p
        !            if (rand()<nu_plus_dt) then
        !                v_p(1,i)=rand_gaussian(sigma1)/mass_p
        !                v_p(2,i)=rand_gaussian(sigma1)/mass_p
        !                v_p(3,i)=rand_gaussian(sigma1)/mass_p
        !            endif
        !        enddo
        !        do i=1,n_s
        !            if (rand()<nu_plus_dt) then
        !                v_s(1,i)=rand_gaussian(sigma1)/mass_s
        !                v_s(2,i)=rand_gaussian(sigma1)/mass_s
        !                v_s(3,i)=rand_gaussian(sigma1)/mass_s
        !            endif
        !        enddo
        !    call cal_Ek_T(Ek,T_out)
    end subroutine

    !    subroutine thermostat_L(Ek,T,T_out)
    !        implicit none
    !        real(8) Ek,T,T_out
    !        real(8) noise, gfric, gama
    !        integer i
    !
    !        gfric=1d0-thermostat_parameter1/2d0;
    !        noise=sqrt(6.0*thermostat_parameter1*T/time_step_s**2)
    !
    !        do i=1,n_p
    !            f0_p(1,i)=f0_p(1,i)+2*noise*(rand()-0.5)
    !            f0_p(2,i)=f0_p(2,i)+2*noise*(rand()-0.5)
    !            f0_p(3,i)=f0_p(3,i)+2*noise*(rand()-0.5)
    !        enddo
    !        do i=1,n_s
    !            f_s(1,i)=2*noise*(rand()-0.5)
    !            f_s(2,i)=2*noise*(rand()-0.5)
    !            f_s(3,i)=2*noise*(rand()-0.5)
    !        enddo
    !        !call cal_Ek_T(Ek,T_out)
    !    end subroutine

    subroutine thermostat_MCS(Ek, n_c, scalar)
        implicit none
        real(8) Ek, scalar
        integer n_c

        real(8) :: zeta = 0.03
        real(8) epsilon0, A

        epsilon0=rand()*zeta
        if(rand()>0.5)then
            scalar=1+epsilon0
        else
            scalar=1/(1d0+epsilon0)
        endif
        A=scalar**(3*(n_c-1))*exp(-(scalar*2-1)*Ek/T_set)
        if(rand()>A)then
            scalar=1d0
        endif
    end subroutine

    subroutine thermostat_MBS(Ek, n_c, scalar)
        implicit none
        real(8) Ek, scalar
        integer n_c, s
        real(8) p, x, y, r(1) !, gamma(100000)
#ifdef __INTEL_COMPILER
        s=vdrnggamma(VSL_RNG_METHOD_GAMMA_GNORM,vsl_stream,1,r,3d0/2*(n_c-1),0d0,T_set)
#else
        r(1)=1d0
#endif
        scalar=sqrt(r(1)/Ek)
        !write(*,*)scalar
        !    return

        !        do while(.true.)
        !            x=rand()*10
        !            p=1/(x*gamma(3d0/2*(n_c-1)))*(x/T_set)**(3d0/2*(n_c-1))*exp(-x/T_set)
        !            write(*,*)p
        !            y=rand()
        !            if(p>y)then
        !                scalar=sqrt(x/Ek)
        !                exit
        !            endif
        !        enddo
    end subroutine

    subroutine output(output_file,cur_step,step)
        implicit none
        integer cur_step,step,output_file,k

        if(mod(cur_step,step)==0)then
            write(output_file,*)'ITEM:TIMESTEP'
            write(output_file,'(I9)')cur_step
            write(output_file,*)'ITEM:NUMBER OF ATOMS'
            write(output_file,'(I6)')n_p+n_s
            write(output_file,*)'ITEM:BOX BOUNDS'
            write(output_file,'(2F7.1)')-radius,radius
            write(output_file,'(2F7.1)')-radius,radius
            write(output_file,'(2F7.1)')-n_cell_z/2.0,n_cell_z/2.0
            write(output_file,*)'ITEM:ATOMS id type x y z'
            !            do k=1,n_b
            !                write(output_file,'(2I6,3F13.4)') k,1,x_b(:,k)
            !            enddo
            do k=1,n_p
                write(output_file,'(2I6,3F13.4)') n_b+k,1,x_p(:,k)
            enddo
            do k=1,n_s
                write(output_file,'(2I6,3F13.4)') n_b+n_p+k,2,x_s(:,k)
            enddo
        endif

    end subroutine

    subroutine output_U(energy_file,cur_step,step)
        implicit none
        integer cur_step, step, energy_file
        if(mod(cur_step,step)==0)then
            write(energy_file,*) cur_step,U
        endif
    end subroutine

    subroutine cross_border(x,cross_cylinder,cross_plane)
        implicit none
        real(8), dimension(3):: x(3)
        real(8) r
        logical cross_cylinder,cross_plane
        r=sqrt(x(1)**2+x(2)**2)
        cross_cylinder=(x(3)<0.and.r>radius_thick).or.(x(3)>0.and.r>radius_thin)
        cross_plane=r>radius_thin.and.r<radius_thick.and.x(3)>0
    endsubroutine

    subroutine cross(x,x_0,xc)
        implicit none
        integer i,n
        real(8),dimension(3):: x, x_0, xc
        real(8),dimension(2):: ss, cc, xm
        real(8) sc, t, t1, a, b, c, d, det
        logical check_cylinder,check_plane
        !!越界与界面的交点

        call cross_border(x,check_cylinder,check_plane)
        if(check_cylinder)then
            if(sqrt(x(1)**2+x(2)**2)>radius_thick)then
                r=radius_thick
            elseif(sqrt(x(1)**2+x(2)**2)>radius_thin.and.sqrt(x(1)**2+x(2)**2)<radius_thick)then
                r=radius_thin
            endif
            xm=x(1:2)-x_0(1:2)
            ! solve equation sum((t*xm-x0)**2)=radius**2
            ! 对于a*t^2+b*t+c=0，c必定小于0，因此解必有一正一负，仅取正值
            a=xm(1)**2+xm(2)**2
            b=2*x_0(1)*xm(1)+2*x_0(2)*xm(2)
            c=x_0(1)**2+x_0(2)**2-r**2
            delta=b**2-4*a*c
            if (delta<0) cycle
            t=(-b+sqrt(delta))/2/a
            if (t<0 .or. t>1) cycle
            ! 旋转，最好参看配图
            xc(1:2)=x_0(1:2)+(x(1:2)-x_0(1:2))*t
            xc(3)=(xc(1)*x(3)-x_0(1)*x(3)+x(1)*x_0(3)-xc(1)*x_0(3))/(x(1)-x_0(1))
        elseif(check_plane)then
            xc(3)=0d0
        end if
    end subroutine

    subroutine bounce_back(x,x_0,v,n)
        implicit none
        integer i,n
        real(8), dimension(2):: x0,v1,xc
        real(8), dimension(3,n):: x, x_0, v
        real(8) cc,ss,t,det,a,b,c,d

        !$omp parallel do private(x,v1,v,cc,ss,t,det)
        do i=1,n
            ! 越界则回弹
            call cross(x(:,i),x_0(:,i),xc)
            if(xc(3)/=0)then
                det=xc(1)**2+xc(2)**2
                if (det==0) cycle
                cc=x(1,i)*xc(1)
                ss=x(2,i)*xc(2)
                t=1-(cc+ss)/det
                x(1:2,i)=2d0*(x(1:2,i)+xc(1:2)*t)-x(1:2,i)
                ! x(3,i)=x1(3)
                v1=(x(1:2,i)-xc(1:2))*norm2(v(1:2,i))/norm2(x(1:2,i)-xc(1:2))
                v(1:2,i)=v1
            else
                x(3,i)=2*xc(3)-x(3,i)
            endif
        enddo
        !$omp end parallel do
    end subroutine

    subroutine one_step(cur_step,output_file)
        implicit none

        integer cur_step, output_file, i
        real(8) :: EK_scaled,T_scaled
        ! solvent
        x0_s=x_s
        x_s = x_s + v_s*time_step_s
        !v_s = v_s + 0.5*f_s*time_step_s
        call bounce_back(x_s,x0_s,v_s,n_s)
        call periodic_s()
        ! polymer chain
        do i=1,int(time_step_s/time_step_p)

            x0_p=x_p
            x_p = x_p + v_p*time_step_p + 0.5*f0_p*time_step_p**2

            call bounce_back(x_p,x0_p,v_p,n_p)
            call periodic_p()
            call update_force(1)
            v_p = v_p + 0.5*(f0_p+f_p)*time_step_p
            f0_p=f_p
        enddo

        call cal_collision_velocity(cur_step)
        if (thermostat_method<10) call thermostat(cur_step)

        if(mod(cur_step,desk_interval_step)==0)then
            call Ek_T(EK_scaled, T_scaled)
            write(*,'(I7,5F12.3)') cur_step, U_BEND, U_FENE, U_LJ, U,T_scaled
        endif

        call output(output_file,cur_step,equili_interval_step)
        !call output_U(energy_file,cur_step,equili_interval_step)
    end subroutine


end module parameters


program Poisellie_field
    use parameters
    implicit none
    integer :: cur_step,output_file,energy_file,production_file,velocity_file
    integer i,j, h_p,n(0:20)
    real(8) :: EK_scaled,T_scaled,sum_v(0:20),r

    output_file=12
    energy_file=13
    production_file=14
    velocity_file=15

    !gama=0.001

    call readin()

    box_size = [n_cell_x, n_cell_y, n_cell_z]
    half_box_size(3) = n_cell_z/2d0
    box_size_unit=1.0
    half_box_size_unit=box_size_unit/2

    !!!读链的大小 改成1个文件
    open(output_file,file='dump.cylinder.lammpstrj')
    open(energy_file,file='energy.out')
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
    write(*,*) ''
    write(*,*)'Equilibrium begin:'
    write(*,'(A7,5A12)') 'step', 'BEND','FENE','LJ','total','T_scaled'
    write(*,*) '--------------------------------------------------------------------'
    write(*,'(I7,5F12.3)') 0, U_BEND, U_FENE, U_LJ, U,T_scaled
    open(velocity_file,file='velocity_radius')
    sum_v=0
    n=0
    do cur_step=1,equili_step
        ! write(*,*) U
        call one_step(cur_step, output_file)
        !        call date_and_time(TIME=time0)
        !        write(*,*) time0, f_p(:,20)
        if(mod(cur_step,output_interval_step)==0)then
            do i=1,n_s
                if(x_s(3,i)<2.0 .and. x_s(3,i)>-2.0)then
                    r=sqrt(x_s(1,i)**2+x_s(2,i)**2)
                    j=floor(r*5)
                    sum_v(j)=sum_v(j)+v_s(3,i)
                    n(j)=n(j)+1
                end if
            enddo

        end if
    enddo
    do j=0,20
        sum_v(j)=sum_v(j)/(total_step/output_interval_step)
        write(velocity_file,*)j,sum_v(j)/n(j)
    enddo


    write(*,*)
    write(*,*)'Production begin'
    open(production_file,file='dump.production.lammpstrj')
    !!! compute a(t-dt)
    call update_force(0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,'(A7,5A12)') 'step', 'BEND','FENE','LJ', 'total','T_scaled'

    write(*,*) '--------------------------------------------------------------------'

    sum_v=0
    n=0
    do cur_step=1,total_step
        v_s(3,:) = v_s(3,:) + gama !- gama*(x_s(1,:)**2+x_s(2,:)**2)/radius**2
        call one_step(cur_step, production_file)
        if(mod(cur_step,output_interval_step)==0)then
            do i=1,n_s
                if(x_s(3,i)<2.0 .and. x_s(3,i)>-2.0)then
                    r=sqrt(x_s(1,i)**2+x_s(2,i)**2)
                    j=floor(r*5)
                    sum_v(j)=sum_v(j)+v_s(3,i)
                    n(j)=n(j)+1
                end if
            enddo

        end if

    enddo
    do j=0,20
        sum_v(j)=sum_v(j)/(total_step/output_interval_step)
        write(velocity_file,*)j,sum_v(j)/n(j)
    enddo
    close(output_file)
    close(energy_file)
    close(production_file)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program Poisellie_field
