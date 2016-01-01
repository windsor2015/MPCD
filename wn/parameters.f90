
module parameters

#if defined (_T_TUBE)
    use shape_t_tube
#elif defined (_T_TUBE1)
    use shape_t_tube1
#elif defined (_FUNNEL)
    use shape_funnel
#else
    use shape_cylinder
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
        real(8) dx(2,4), theta, r, d, t
        logical success
        integer,parameter :: seed = 11111

        select case (string_form)
            case (0)
                !!!!!!!!!以下是polymer chain初值!!!!!!!!!!!!!!
                x_p(1,1)=0d0
#if defined (_T_TUBE)
                x_p(2,1)=n_cell_y*(1d0-ratio_y)/2d0+n_cell_y*ratio_y/4d0
                x_p(3,1)=-n_cell_z/2d0+1d0
#elif defined (_T_TUBE1)
                x_p(2,1)=n_cell_y*(1d0-ratio_y)/2d0+n_cell_y*ratio_y/4d0
                x_p(3,1)=-n_cell_z/2d0+1d0
#elif defined (_FUNNEL)
                x_p(2,1)=0
                x_p(3,1)=-n_cell_z/2d0+2d0
#else
                x_p(2,1)=0d0
                x_p(3,1)=-n_cell_z/2d0+2d0
#endif
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
#if defined (_T_TUBE)
                        if(((x_p(2,i)<0.or.x_p(2,i)>n_cell_y/2d0).and.abs(x_p(3,i))<=n_cell_z*ratio_z/2d0) &
                                .or.((x_p(2,i)<n_cell_y*(1d0-ratio_y)/2d0.or.x_p(2,i)>n_cell_y/2d0).and. &
                                abs(x_p(3,i))>n_cell_z*ratio_z/2d0)) then
                            cycle
                        endif
#elif defined (_T_TUBE1)
                        if(abs(x_p(1,i))>n_cell_x/2d0 .or. &
                                ((x_p(2,i)<0.or.x_p(2,i)>n_cell_y/2d0).and.abs(x_p(3,i))<=n_cell_z*ratio_z/2d0) &
                                .or.((x_p(2,i)<n_cell_y*(1d0-ratio_y)/2d0.or.x_p(2,i)>n_cell_y/2d0).and. &
                                abs(x_p(3,i))>n_cell_z*ratio_z/2d0)) then
                            cycle
                        endif
#else
                        if((sqrt(x_p(1,i)**2+x_p(2,i)**2)>radius(1).and.abs(x_p(3,i))>=n_cell_z*ratio_z/2d0) &
                                .or.(sqrt(x_p(1,i)**2+x_p(2,i)**2)>radius(2).and.abs(x_p(3,i))<n_cell_z*ratio_z/2d0)) then
                            cycle
                        endif
#endif
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
            case (1)
                do i=1,n_p
                    t=-pi+2*pi*i/n_p
                    x_p(1,i)=sin(t)+2*sin(2*t)
#if defined (_T_TUBE)
                    x_p(2,i)=-sin(3*t)+(1d0-ratio_y)*n_cell_y/2d0+ratio_y*n_cell_y/4d0
                    x_p(3,i)=cos(t)-2*cos(2*t)-n_cell_z*ratio_z
#elif defined (_T_TUBE1)
                    x_p(2,i)=-sin(3*t)+(1d0-ratio_y)*n_cell_y/2d0+ratio_y*n_cell_y/4d0
                    x_p(3,i)=cos(t)-2*cos(2*t)-n_cell_z*ratio_z
#elif defined (_FUNNEL)
                    x_p(2,i)=cos(t)-2*cos(2*t)
                    x_p(3,i)=-sin(3*t)-n_cell_z/2+5d0
#else
                    x_p(2,i)=cos(t)-2*cos(2*t)
                    x_p(3,i)=-sin(3*t)-n_cell_z*ratio_z
#endif
                end do
                call periodic_p()
        end select



        write(*,*)"Polymer monomer number: ", n_p

        !!!!!!!!!!!!!!!!以下是solution粒子的初值!!!!!!!!!!!!!!

        n_s=get_pipe_volume()*density_s

        allocate(x_s(3,n_s),v_s(3,n_s),x0_s(3,n_s),x_s0(3,n_s))
        i=0
        do while(.true.)
            i=i+1
            if(i>n_s)exit
            x_s(1,i)=(rand()-5d-1)*n_cell_x
            x_s(2,i)=(rand()-5d-1)*n_cell_y
            x_s(3,i)=(rand()-5d-1)*n_cell_z
            if(.not. in_pipe(x_s(:,i)))then
                i=i-1
            end if
        enddo

        !!!!!!!!!!phantom particles

        n_b=get_pipe_phantom_volume()*density_s

        allocate(x_b(3,n_b),v_b(3,n_b))
        i=0
        do while(.true.)
            i=i+1
            if(i>n_b)exit
            d=sqrt(5d-1)
            x_b(1,i)=(rand()-5d-1)*(n_cell_x+2)
            x_b(2,i)=(rand()-5d-1)*(n_cell_y+2)
            x_b(3,i)=(rand()-5d-1)*n_cell_z
            if(.not. in_pipe_phantom(x_b(:,i)))then
                i=i-1
            end if
        enddo

        min_y=minval(int(x_s(2,:)))
        max_y=maxval(int(x_s(2,:)))

        write(*,*)"Solvent particle number: ", n_s
        write(*,*)"Total particle number: ", n_p+n_s
        write(*,*)"Phantom particle number: ", n_b

    endsubroutine

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
        real(8), parameter :: FENE_rc=1.5*sigma, FENE_k=30
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
        !rx(1)=rx(1)-n_cell_x*nint(rx(1)/n_cell_x)
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
        integer i
        real(8) f(3,3), U, rx1(3), rx2(3), c, s, ss, r1, r2, small, d(3,3)
        small=0.001
        !real(8), parameter :: BEND_b = 10
        rx1(3)=rx1(3)-n_cell_z*nint(rx1(3)/n_cell_z)
        rx2(3)=rx2(3)-n_cell_z*nint(rx2(3)/n_cell_z)

        r1=norm2(rx1)
        r2=norm2(rx2)
        c=-dot_product(rx1,rx2)/(r1*r2)

        if((dble(1)-c**2)<=0)then
            s=0
        else
            s=sqrt(dble(1)-c**2)
        endif
        if(c>dble(1))then
            c=dble(1)
        elseif(c<-dble(1))then
            c=-dble(1)
        endif
        if(s<small)then
            s=small
        endif

        s=1/s
        ss=BEND_b*(-1/s)

        d(:,1)=(-c*rx1*r2/r1-rx2)*s/(r2*r1)

        d(:,2)=(c*(rx1*r2/r1-rx2*r1/r2)-rx1+rx2)*s/(r1*r2)

        d(:,3)=(c*rx2*r1/r2+rx1)*s/(r1*r2)

        f=f-ss*d

        !        c=dot_product(rx1,rx2)/r1/r2
        !        f=f-BEND_b*((rx1+rx2)/(r1*r2)-c*rx1/r1/r1-c*rx2/r2/r2)
        U=U+BEND_b*(1+c)

    end subroutine

    subroutine exter(f,flag)
        implicit none
        real(8) f
        integer flag
        if(flag==0)then
            f=f+0
        elseif(flag==1)then
            f=f+0
        end if

    end subroutine

    subroutine update_force(mode,flag)
        implicit none
        integer mode, i, j,flag
        real(8) temp(3)

        f_p=0
        U=0
        U_BEND=0
        U_FENE=0
        U_LJ=0
        temp=0
        ! 链两端
        call FENE(f_p(:,1),U_FENE,x_p(:,1)-x_p(:,2))
        call FENE(f_p(:,n_p),U_FENE,x_p(:,n_p)-x_p(:,n_p-1))

        !$omp parallel do private(j) reduction(+:U_FENE,U_LJ,U_BEND,f_p)
        do i=1,n_p
            if (i>1 .and. i<n_p) then
                ! UFENE(r) force
                call FENE(f_p(:,i),U_FENE,x_p(:,i)-x_p(:,i+1))
                call FENE(f_p(:,i),U_FENE,x_p(:,i)-x_p(:,i-1))

                ! bend energy
                call BEND(f_p(:,i-1:i+1),U_BEND,x_p(:,i)-x_p(:,i-1),x_p(:,i+1)-x_p(:,i))
            endif

            do j=1,n_p
                if (j/=i) then
                    call LJ(f_p(:,i),U_LJ,x_p(:,i)-x_p(:,j))
                endif
            enddo

        enddo
        !$omp end parallel do
        call exter(f_p(3,1),flag)

        U=U_FENE+U_BEND+U_LJ

        if (mode==0) then
            f0_p=f_p
        endif

    end subroutine

    logical function inbox(x,ix,iy,iz)
        implicit none
        real(8) x(3), temp(3)
        integer ix,iy,iz

        temp=x/box_size_unit-[ix,iy,iz]*1d0
        inbox = temp(1)>=0 .and. temp(1)<1d0 &
            .and. temp(2)>=0 .and. temp(2)<1d0 &
            .and. temp(3)>=0 .and. temp(3)<1d0
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
    end subroutine

    subroutine cal_collision_velocity(cur_step)
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

        call random_number(randr)
        randr=(randr-0.5)*box_size_unit

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

        ! phantom particle, the velocity
        sigma1=sqrt(T_set)
        i=vdrnggaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, vsl_stream, n_b*3, v_b, 0d0, sigma1)
        !$omp parallel do private(ix,iy,iz,check,temp)
        do i=1,n_b
            call get_cell_xyz(x_b(:,i)+randr,ix,iy,iz,check)
            if (check) then
                count_cell_b(ix,iy,iz)=count_cell_b(ix,iy,iz)+1
                pointer_cell_b(ix,iy,iz,count_cell_b(ix,iy,iz))=i
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
                        v_p(:,k)=v_p(:,k)-v_aver
                        Ek=Ek+mass_p*sum(v_p(:,k)**2)
                    enddo
                    do i=1,count_s
                        k=pointer_cell_s(ix,iy,iz,i)
                        v_s(:,k)=v_s(:,k)-v_aver
                        Ek=Ek+mass_s*sum(v_s(:,k)**2)
                    enddo
                    do i=1,count_b
                        k=pointer_cell_b(ix,iy,iz,i)
                        v_b(:,k)=v_b(:,k)-v_aver
                        Ek=Ek+mass_s*sum(v_b(:,k)**2)
                    enddo
                    Ek=Ek/2d0

                    scalar=1d0
                    if (thermostat_method>=10) then
                        scalar=thermostat_cal_scalar(cur_step, count_cell_p(ix,iy,iz)+count_cell_s(ix,iy,iz)+count_cell_b(ix,iy,iz), Ek)
                    end if

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
        sigma1=sqrt(mass_s*T_set)    !!？？？？？？？？？？？？？？？？？？？可以吗？
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
        !        enddo5d0-x_0(3,i)
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

        !if(mod(cur_step,step)==0)then
        write(output_file,*)'ITEM:TIMESTEP'
        write(output_file,'(I9)')cur_step
        write(output_file,*)'ITEM:NUMBER OF ATOMS'
        write(output_file,'(I6)')n_p+n_s!+n_b
        write(output_file,*)'ITEM:BOX BOUNDS'
        write(output_file,'(2F7.1)')-n_cell_x/2d0-1,n_cell_x/2d0+1
        write(output_file,'(2F7.1)')-n_cell_y/2d0-1,n_cell_y/2d0+1
        write(output_file,'(2F7.1)')-n_cell_z/2d0,n_cell_z/2d0
        write(output_file,*)'ITEM:ATOMS id type x y z'

        do k=1,n_p
            write(output_file,'(I6,I3,3F9.4)') k,1,x_p(:,k)
        enddo
        do k=1,n_s
            write(output_file,'(I6,I3,3F9.4)') n_p+k,2,x_s(:,k)
        enddo
        do k=1,n_b
            !write(output_file,'(I6,I3,3F9.4)') n_p+n_s+k,3,x_b(:,k)
        enddo
        !endif

    end subroutine

    subroutine output_U(energy_file,cur_step,step)
        implicit none
        integer cur_step, step, energy_file
        if(mod(cur_step,step)==0)then
            write(energy_file,*) cur_step,U
        endif
    end subroutine


    subroutine one_step(cur_step,interval_step,output_file,flag)
        implicit none

        integer cur_step, output_file, i, interval_step,flag
        real(8) :: EK_scaled,T_scaled,t,min_z,min_z0,count_z

        ! solvent
        x0_s=x_s
        x_s = x_s + v_s*time_step_s

        call bounce_back(x_s,x0_s,v_s,n_s)
        call periodic_s()

        ! polymer chain
        do i=1,int(time_step_s/time_step_p)

            x0_p=x_p
            x_p = x_p + v_p*time_step_p + 0.5*f0_p*time_step_p**2

            if (count(x_p(3,:)>0)>0 .and. count(x0_p(3,:)>0)>0) then
                min_z=minval(x_p(3,:),x_p(3,:)>0)
                min_z0=minval(x0_p(3,:),x0_p(3,:)>0)
                if (min_z>=1 .and. min_z0<1) then
                    cross_flag=1
                endif
            endif

            call bounce_back(x_p,x0_p,v_p,n_p)

            call periodic_p()
            call update_force(1,flag)
            v_p = v_p + 0.5*(f0_p+f_p)*time_step_p
            f0_p=f_p
        enddo

        call cal_collision_velocity(cur_step)
        if (thermostat_method<10) call thermostat(cur_step)

        if(mod(cur_step,desk_interval_step)==0)then
            call Ek_T(EK_scaled, T_scaled)
            call get_time(t)
            write(*,'(I7,6F12.3)') cur_step,U_BEND,U_FENE,U_LJ,U,T_scaled,t-time0
            if (output_sketch/=0) then
                call output_sketch_sub()
            end if
            time0=t
        endif
        if (mod(cur_step,interval_step)==0)then
            call output(output_file,cur_step,interval_step)
            if (cross_flag==1 .and. min_z>n_cell_z/3d0) cross_flag=2
        endif
        !call output_U(energy_file,cur_step,interval_step)
    end subroutine


    subroutine output_sketch_sub()
        implicit none
        integer x_p_2d(n_p)
        integer i,j,k,l,x_p_int(2,n_p)

        x_p_int(1,:)=int(-x_p(2,:))
        !min_y=minval(x_p_int(1,:))
        !max_y=maxval(x_p_int(1,:))
        x_p_int(1,:)=x_p_int(1,:)-min_y+1
        x_p_int(2,:)=int((x_p(3,:)+n_cell_z/2d0)*80d0/n_cell_z+1)
        x_p_2d=x_p_int(1,:)*10000+x_p_int(2,:)

        call quick_sort(x_p_2d,n_p,1,n_p)
        !do i=1,n_p-1
        !    do j=i+1,n_p
        !        if (x_p_2d(i)>x_p_2d(j)) call swap(x_p_2d(i),x_p_2d(j))
        !    end do
        !end do

        x_p_int(1,:)=x_p_2d/10000
        x_p_int(2,:)=mod(x_p_2d,10000)

        k=1
        do i=1,max_y-min_y+1
            do j=1,80
                if (x_p_int(1,k)==i .and. x_p_int(2,k)==j) then
                    write(*,'(A,$)') '*'
                    !k=k+1
                    do l=k,n_p
                        if (x_p_int(2,l)/=x_p_int(2,k)) then
                            k=l
                            exit
                        end if
                    end do

                    if (k>n_p) exit
                else
                    write(*,'(A,$)') ' '
                endif
            enddo
            write(*,*)
        enddo

        !write(*,*) x_p_int
        !stop

    end subroutine

    ! 交换
    subroutine swap(a,b)
        integer a,b,c
        c=a
        a=b
        b=c
    end subroutine

    ! 快速排序
    recursive subroutine quick_sort(v,n,s,e)
        integer v(n),key
        integer n,s,e,l,r,m

        l=s
        r=e
        m=(s+e)/2
        if (l>=r) return
        key=v(m)
        do while (l<r)
            do while(v(l)<key)
                l=l+1
            enddo
            do while(v(r)>key)
                r=r-1
            enddo
            if (l<r) then
                call swap(v(l),v(r))
                if (v(r)==key) l=l+1
                if (v(l)==key) r=r-1
            endif
        enddo

        call quick_sort(v,n,s,l-1)
        call quick_sort(v,n,r+1,e)
    end subroutine


    subroutine get_time(t)
        implicit none
        real(8) t
        integer(8) t1, clock_rate, clock_max
        call system_clock(t1,clock_rate,clock_max)
        t=1d0*t1/clock_rate
        !time0=t1
    end subroutine

    subroutine output_date()
        implicit none
        integer d(8)

        call date_and_time(values=d)
        write (*,'(I5,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I3.3)') &
            d(1),'-',d(2),'-',d(3),'/',d(5),':',d(6),':',d(7),'.',d(8)
    end subroutine

end module parameters
