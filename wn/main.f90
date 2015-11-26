
module parameters
#ifdef __INTEL_COMPILER
    use ifport
#endif
    implicit none
    ! 下标约定: p - polymer, s - solution, b - boundary

    !结构
    integer, parameter :: n_p=40, n_cell_x=12, n_cell_y=12, n_cell_z=40

    integer n_b, n_s

    !n是单体数目，n0是单体、壁、孔数目之和

    real(8), parameter :: kB = 1.38064852d-23
    real(8), parameter :: pi = 3.141592653589793238462643383279

    real(8), parameter :: time_step_p=1d-4, time_step_s=5d-3, mass_p=1, mass_s=1, T_set=1, v_gradient=0.2

    ! polymer 位置 速度 力 上一次力
    real(8), dimension(3,n_p) :: x_p, v_p, f_p, f0_p, x0_p
    ! solution
    real(8), allocatable, dimension(:,:) :: x_s, v_s, f_s, x0_s, x_s0
    ! boundaaries, 1~nb-up, nb+1~2nb-down
    real(8), allocatable, dimension(:,:) :: x_b, v_b, f_b

    real(8) box_size_unit, half_box_size_unit, box_size(3), half_box_size(3)

    real(8) aver_v(3), distant

    real(8), parameter :: sigma=1, epson=1

    real(8) :: density_s=5, radius=4, gama=0.001

    real(8) U, U_LJ, U_FENE, U_BEND, U_WALL

    ! pointer - 每个格子中的粒子编号
    ! count - 每个格子中的粒子计数，1 - p，2 - s
    ! momentum - 每个格子中的总动量
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z, n_p) :: pointer_cell_p
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z, 100) :: pointer_cell_s
    integer, dimension(0:n_cell_x, 0:n_cell_y, 0:n_cell_z) :: count_cell_p, count_cell_s
    real(8), dimension(3,0:n_cell_x, 0:n_cell_y, 0:n_cell_z) :: momentum_cell

    integer :: desk_interval_step=10000
    integer :: equili_step=50000
    integer :: equili_interval_step=1000
    integer :: total_step=3000000
    integer :: output_interval_step=1000

    namelist /basic/ radius, gama, density_s, &
    desk_interval_step, equili_step, equili_interval_step, total_step, output_interval_step

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
        integer i,j,k,h,count_number
        real(8) distant,dx(2,4), theta, r
        logical success
        integer,parameter :: seed = 86456

        !!!!!!!!!以下是cylinder channel初值!!!!!!!!!!!!!!
!        n_b=0
!        k=nint(2.0*pi*radius)
!        do j=0,2*n_cell_z
!            do i=0,2*k-1
!                n_b=n_b+1
!                if(mod(j,2)==0)then
!                    temp(1,n_b)=radius*cos(i*pi/k)
!                    temp(2,n_b)=radius*sin(i*pi/k)
!                else
!                    temp(1,n_b)=radius*cos(i*pi/k+pi/(2d0*k))
!                    temp(2,n_b)=radius*sin(i*pi/k+pi/(2d0*k))
!                endif
!                temp(3,n_b)=(j-n_cell_z)*0.5
!                ! write(output_file,'(2I6,3F13.4)') n_b,1,x_b(:,n_b)
!            enddo
!        enddo
!        allocate(x_b(3,n_b),v_b(3,n_b),f_b(3,n_b))
!        x_b=temp(:,1:n_b)
!        write(*,*)"Cylinder particle number: ", n_b

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

        n_s=nint(pi*radius**2*n_cell_z*density_s)

        allocate(x_s(3,n_s),v_s(3,n_s),f_s(3,n_s),x0_s(3,n_s),x_s0(3,n_s))
        do i=1,n_s
            theta=pi*2*rand()
            r=sqrt(rand())*radius
            x_s(1,i)=r*cos(theta)
            x_s(2,i)=r*sin(theta)
            x_s(3,i)=(rand()-0.5)*n_cell_z
        enddo

        write(*,*)"Solvent particle number: ", n_s
        write(*,*)"Total particle number: ", n_p+n_s

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
        real(8), parameter :: FENE_rc=1.5*sigma, FENE_k=100
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
        real(8), parameter :: BEND_b = 100

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
        z=floor(r(3)+n_cell_z/2d0)

        s= x>=0 .and. x<=n_cell_x &
            .and. y>=0 .and. y<=n_cell_y &
            .and. z>=0 .and. z<=n_cell_z

        !if (.not.s) write(*,*) r

    end subroutine

    subroutine cal_collision_velocity()
        !    use parameters
        implicit none
        integer ix,iy,iz,k,count_p,count_s,i
        real(8) momentum(3), matrix(3,3), l(3), fai, theta
        real(8), parameter :: alpha = 130*pi/180, s=sin(alpha), c=cos(alpha)
        real(8) v_aver_p(3), v_aver_s(3), v_aver(3), temp(3), randx, randy, randz, R(2,2), s1, c1
        logical mask_p(n_p), mask_s(n_s), check

        pointer_cell_p=0
        pointer_cell_s=0
        count_cell_p=0
        count_cell_s=0
        momentum_cell=0d0
        ! calculate velocity of all particles in each cell
        randx=(rand()-0.5)*box_size_unit
        randy=(rand()-0.5)*box_size_unit
        randz=(rand()-0.5)*box_size_unit

        !randx=rand()*pi*2
!        c1=cos(randx)
!        s1=sin(randx)

!        R(1,1)=c1
!        R(1,2)=-s1
!        R(2,1)=s1
!        R(2,2)=c1

!        x0_p(1:2,:)=matmul(R,x_p(1:2,:))
!        x0_s(1:2,:)=matmul(R,x_s(1:2,:))
!        x0_p(1:2,:)=x_p(1:2,:)
!        x0_s(1:2,:)=x_s(1:2,:)
        x0_p(1,:)=x_p(1,:)+randx
        x0_s(1,:)=x_s(1,:)+randx
        x0_p(2,:)=x_p(2,:)+randy
        x0_s(2,:)=x_s(2,:)+randy
        x0_p(3,:)=x_p(3,:)+randz
        x0_s(3,:)=x_s(3,:)+randz
        x0_p(3,:)=x0_p(3,:)-n_cell_z*nint((x0_p(3,:))/n_cell_z)
        x0_s(3,:)=x0_s(3,:)-n_cell_z*nint((x0_s(3,:))/n_cell_z)
        !$omp parallel do private(ix,iy,iz,check)
        do i=1,n_p
            call get_cell_xyz(x0_p(:,i),ix,iy,iz,check)
            if (.not.check) cycle
            count_cell_p(ix,iy,iz)=count_cell_p(ix,iy,iz)+1
            pointer_cell_p(ix,iy,iz,count_cell_p(ix,iy,iz))=i
            momentum_cell(:,ix,iy,iz)=momentum_cell(:,ix,iy,iz)+mass_p*v_p(:,i)
        enddo
        !$omp end parallel do
        !$omp parallel do private(ix,iy,iz,check)
        do i=1,n_s
            call get_cell_xyz(x0_s(:,i),ix,iy,iz,check)
            if (check) then
                count_cell_s(ix,iy,iz)=count_cell_s(ix,iy,iz)+1
                pointer_cell_s(ix,iy,iz,count_cell_s(ix,iy,iz))=i
                momentum_cell(:,ix,iy,iz)=momentum_cell(:,ix,iy,iz)+mass_s*v_s(:,i)
            endif
        enddo
        !$omp end parallel do
        !$omp parallel do private(iy,iz,i,matrix,l,fai,theta,v_aver,temp,k)
        do ix=0,n_cell_x
            do iy=0,n_cell_y
                do iz=0,n_cell_z

                    if (count_cell_p(ix,iy,iz)+count_cell_s(ix,iy,iz)==0) cycle

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

                    v_aver=momentum_cell(:,ix,iy,iz)/(mass_p*count_cell_p(ix,iy,iz)+mass_s*count_cell_s(ix,iy,iz))
                    do i=1,count_cell_p(ix,iy,iz)
                        k=pointer_cell_p(ix,iy,iz,i)
                        temp=v_p(:,k)-v_aver
                        v_p(1,k)=v_aver(1) + matrix(1,1)*temp(1) + matrix(1,2)*temp(2) + matrix(1,3)*temp(3)
                        v_p(2,k)=v_aver(2) + matrix(2,1)*temp(1) + matrix(2,2)*temp(2) + matrix(2,3)*temp(3)
                        v_p(3,k)=v_aver(3) + matrix(3,1)*temp(1) + matrix(3,2)*temp(2) + matrix(3,3)*temp(3)
                    enddo
                    do i=1,count_cell_s(ix,iy,iz)
                        k=pointer_cell_s(ix,iy,iz,i)
                        temp=v_s(:,k)-v_aver
                        v_s(1,k)=v_aver(1) + matrix(1,1)*temp(1) + matrix(1,2)*temp(2) + matrix(1,3)*temp(3)
                        v_s(2,k)=v_aver(2) + matrix(2,1)*temp(1) + matrix(2,2)*temp(2) + matrix(2,3)*temp(3)
                        v_s(3,k)=v_aver(3) + matrix(3,1)*temp(1) + matrix(3,2)*temp(2) + matrix(3,3)*temp(3)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do

    end subroutine

    !!!!以下是Isokinetics thermostat，可以试一下Berendsen thermostat
    subroutine scale_v(Ek,T,T_out)
        implicit none
        integer i
        real(8), parameter :: Ek_fac = 1.5d0
        real(8) v, Ek,T, scalar, Ek1, T_out, T1


        !do i=1,3
        !v=sum(v_p(i,:)*mass_p+v_s(i,:)*mass_s)/(n_p*mass_p+n_s*mass_s)
        !v_p(i,:) = v_p(i,:)-v
        !v_s(i,:) = v_s(i,:)-v
        !enddo
        !write(*,*)v

        Ek1=0.5*(mass_p*sum(v_p**2)+mass_s*sum(v_s**2))
        T1=Ek1/(Ek_fac*(n_p+n_s))
        scalar=sqrt(T/T1)
        v_p=v_p*scalar
        v_s=v_s*scalar
        Ek=0.5*(mass_p*sum(v_p**2)+mass_s*sum(v_s**2))
        T_out=Ek/(Ek_fac*(n_p+n_s))

    end subroutine

    subroutine scale_v_p(Ek,T,T_out)
        implicit none
        integer i
        real(8), parameter :: Ek_fac = 1.5d0
        real(8) v, Ek,T, scalar, Ek1, T_out, T1


        !do i=1,3
        !v=sum(v_p(i,:)*mass_p+v_s(i,:)*mass_s)/(n_p*mass_p+n_s*mass_s)
        !v_p(i,:) = v_p(i,:)-v
        !v_s(i,:) = v_s(i,:)-v
        !enddo
        !write(*,*)v

        Ek1=0.5*mass_p*sum(v_p**2)
        T1=Ek1/(Ek_fac*n_p)
        scalar=sqrt(T/T1)
        v_p=v_p*scalar
        Ek=0.5*mass_p*sum(v_p**2)
        T_out=Ek/(Ek_fac*n_p)

    end subroutine
    !
    !        subroutine scale_v_s(Ek,T,T_out)
    !        implicit none
    !        integer i
    !        real(8), parameter :: Ek_fac = 1.5d0
    !        real(8) v, Ek,T, scalar, Ek1, T_out, T1
    !
    !
    !        !do i=1,3
    !            !v=sum(v_p(i,:)*mass_p+v_s(i,:)*mass_s)/(n_p*mass_p+n_s*mass_s)
    !            !v_p(i,:) = v_p(i,:)-v
    !            !v_s(i,:) = v_s(i,:)-v
    !        !enddo
    !        !write(*,*)v
    !
    !        Ek1=0.5*mass_s*sum(v_s**2)
    !        T1=Ek1/(Ek_fac*n_s)
    !        scalar=sqrt(T/T1)
    !        v_s=v_s*scalar
    !        Ek=0.5*mass_s*sum(v_s**2)
    !        T_out=Ek/(Ek_fac*n_s)
    !
    !    end subroutine

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

    subroutine bounce_back(x,v,n,time_step, r)
        implicit none
        integer i,n
        real(8), dimension(2):: x0,x1,v0,v1,xc,xm
        real(8) a,b,c,t,delta,norm_rest,s,det,norm_cs,x(3,n),v(3,n),time_step,r
        !$omp parallel do private(x0,x1,v0,v1,xc,xm,a,b,c,t,delta,norm_rest,det,norm_cs)
        do i=1,n
            ! 越界则回弹
            if (x(1,i)**2+x(2,i)**2>=r**2) then
                x1=x(1:2,i)
                v0=v(1:2,i)

                if (n==n_p) then

                    x0 = x1 - v0*time_step-0.5*f0_p(1:2,i)*time_step**2
                    else
                    x0 = x1 - v0*time_step
                endif
                xm=x1-x0
                ! solve equation sum((t*xm-x0)**2)=radius**2
                ! 对于a*t^2+b*t+c=0，c必定小于0，因此解必有一正一负，仅取正值
                a=xm(1)**2+xm(2)**2
                b=2*x0(1)*xm(1)+2*x0(2)*xm(2)
                c=x0(1)**2+x0(2)**2-r**2
                delta=b**2-4*a*c
                if (delta<0) cycle
                t=(-b+sqrt(delta))/2/a
                !write(*,*) t
                if (t<0 .or. t>1) cycle
                ! 旋转，最好参看配图
                xc=x0+(x1-x0)*t
                det=xm(1)**2+xm(2)**2
                if (det==0) cycle
                c=(xc(1)*xm(1)+xc(2)*xm(2))/det
                s=(xc(2)*xm(1)-xc(1)*xm(2))/det

                norm_cs=sqrt(c**2+s**2)
                c=c/norm_cs
                s=s/norm_cs

                norm_rest=norm2(x1-xc)
                x1(1)=-(c*xc(1)-s*xc(2))
                x1(2)=-(s*xc(1)+c*xc(2))

                v1=x1*norm2(v0)/r
                x1=xc+x1*norm_rest/r

                x(1:2,i)=x1
                v(1:2,i)=v1

            endif
        enddo
        !$omp end parallel do
    end subroutine

end module


program Poissonfield
    use parameters
    implicit none
    real(8) :: Ek, EK_scaled,EK_scaled_p,T_scaled,T_scaled_p,density
    integer :: cur_step,output_file,energy_file,production_file

    character(10) :: time0

    integer i, j, h_p
    real(8) randx, randz

    output_file=12
    energy_file=13
    production_file=14
    gama=0.001
    h_p=50

call readin()

    box_size = [n_cell_x, n_cell_y, n_cell_z]
    half_box_size(3) = n_cell_z/2d0
    box_size_unit=1.0
    half_box_size_unit=box_size_unit/2

    !!!读链的大小 改成1个文件
    open(output_file,file='dump.cylinder.lammpstrj')
    open(energy_file,file='energy.out')
    call init()
    call output(output_file,0,equili_interval_step)

    !!!polymer的初速度
    call random_number(v_p)
    !!!solution的初速度
    call random_number(v_s)
    v_p=v_p-0.5
    v_s=v_s-0.5
    call scale_v(EK_scaled,T_set,T_scaled)
    !    call scale_v_p(Ek_scaled_p, T_set, T_scaled_p)
    !    call scale_v_s(Ek_scaled_s, T_set, T_scaled_s)
    write(*,*) ''
    write(*,*) 'Initial scaled kinetics Ek_scaled: ',Ek_scaled
    write(*,*) 'Initial setted temperature T_set: ',T_set
    write(*,*) 'Initial scaled temperature T_scaled: ',T_scaled
    !    write(*,*) 'Solution: '
    !    write(*,*) 'Initial scaled kinetics Ek_scaled: ',Ek_scaled
    !    write(*,*) 'Initial setted temperature T_set: ',T_set
    !    write(*,*) 'Initial scaled temperature T_scaled: ',T_scaled
    ! 没有外场时，polymer和solution达到平衡
    call update_force(0)
    write(*,*) ''
    write(*,*)'Equilibrium begin:'
    write(*,'(A7,4A12)') 'step', 'BEND','FENE','LJ','total'
    write(*,*) '--------------------------------------------------------------------'
    write(*,'(I7,4F12.3)') 0, U_BEND, U_FENE, U_LJ, U
    do cur_step=1,equili_step
        ! write(*,*) U
        ! solvent
        x_s = x_s + v_s*time_step_s
        call bounce_back(x_s,v_s,n_s,time_step_s, radius)
        call periodic_s()
        do i=1,h_p
        ! polymer chain
        x_p = x_p + v_p*time_step_p + 0.5*f0_p*time_step_p**2
        call bounce_back(x_p,v_p,n_p,time_step_p, radius)
        call periodic_p()
        call update_force(1)
        v_p = v_p + 0.5*(f0_p+f_p)*time_step_p
        !call scale_v(EK_scaled,T_set,T_scaled)
        f0_p=f_p
        enddo
        !call scale_v(EK_scaled,T_set,T_scaled)
        call cal_collision_velocity()
        if (mod(cur_step,10)==0)call scale_v(EK_scaled,T_set,T_scaled)
        ! call scale_v_p(EK_scaled_p,T_set,T_scaled_p)
        ! call scale_v_s(EK_scaled_s,T_set,T_scaled_s)

        if(mod(cur_step,desk_interval_step)==0)then
            write(*,'(I7,4F12.3)') cur_step, U_BEND, U_FENE, U_LJ, U
        endif

        call output(output_file,cur_step,equili_interval_step)
        call output_U(energy_file,cur_step,equili_interval_step)
        !        call date_and_time(TIME=time0)
        !        write(*,*) time0, f_p(:,20)
    enddo

    write(*,*)'Production begin'
    open(production_file,file='dump.production.lammpstrj')
    !!! compute a(t-dt)
    call update_force(0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,'(A7,4A12)') 'step', 'BEND','FENE','LJ', 'total'
    write(*,*) '--------------------------------------------------------------------'

    do cur_step=1,total_step
        v_s(3,:) = v_s(3,:) + gama !- gama*(x_s(1,:)**2+x_s(2,:)**2)/radius**2
        x_s0(1:2,:)=x_s(1:2,:)
        x_s = x_s + v_s*time_step_s
        call bounce_back(x_s,v_s,n_s,time_step_s,radius)
        call periodic_s()
        !v_s(3,:) = v_s(3,:) - gama + gama*(x_s0(1,:)**2+x_s0(2,:)**2)/radius**2
        !x_s0(1:2,:)=x_s(1:2,:)

        do i=1,h_p
        ! polymer chain
        x_p = x_p + v_p*time_step_p + 0.5*f0_p*time_step_p**2
        call bounce_back(x_p,v_p,n_p,time_step_p,radius)
        call periodic_p()
        call update_force(1)
        v_p = v_p + 0.5*(f0_p+f_p)*time_step_p
        !call scale_v_p(EK_scaled_p,T_set,T_scaled_p)
        f0_p=f_p
        enddo
        !call scale_v(EK_scaled,T_set,T_scaled)
        !v_s(3,:) = v_s(3,:) + gama - gama*(x_s(1,:)**2+x_s(2,:)**2)/radius**2

        call cal_collision_velocity()
        if (mod(cur_step,10)==0)call scale_v(EK_scaled,T_set,T_scaled)
        !v_s(3,:) = v_s(3,:) - gama + gama*(x_s(1,:)**2+x_s(2,:)**2)/radius**2

        ! call scale_v_p(EK_scaled_p,T_set,T_scaled_p)
        ! call scale_v_s(EK_scaled_s,T_set,T_scaled_s

        if(mod(cur_step,desk_interval_step)==0)then
            write(*,'(I7,5F12.3)') cur_step, U_BEND, U_FENE, U_LJ, U
        endif

        call output(production_file,cur_step,output_interval_step)

    enddo

    close(output_file)
    close(energy_file)
    close(production_file)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program Poissonfield
