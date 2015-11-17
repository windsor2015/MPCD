!DIR$ ATTRIBUTES FORCEINLINE :: FENE, LJ, BEND, inbox, scale_v
module parameters
    !use ifport
    implicit none
    ! 下标约定: p - polymer, s - solution, b - boundary

    !结构
    integer, parameter :: n_p=40, n_cell_x=5, n_cell_y=5, n_cell_z=15, n_s=n_cell_x*n_cell_y*n_cell_z

    integer n_b

    !n是单体数目，n0是单体、壁、孔数目之和

    real(8), parameter :: kB = 1.38064852d-23
    real(8), parameter :: pi = 3.141592653589793238462643383279

    real(8), parameter :: time_step_p=0.0001, time_step_s=0.0001, mass_p=1, mass_s=0.1, T_set=1, v_gradient=0.2

    ! polymer 位置 速度 力 上一次力
    real(8), dimension(3,n_p) :: x_p, v_p, f_p, f0_p
    ! solution
    real(8), dimension(3,n_s) :: x_s, v_s, f_s, x0_s
    ! boundaaries, 1~nb-up, nb+1~2nb-down
    real(8), dimension(3,2*n_b) :: x_b, v_b, f_b

    real(8) box_size_unit, half_box_size_unit, box_size(3), half_box_size(3)

    real(8) randnum_group(3,n_s)

    real(8) aver_v(3)

    real(8), parameter::sigma=1, epson=1

contains

    subroutine test_rand
        integer,parameter :: seed = 86456

        call srand(seed)
        print *, rand(), rand(), rand(), rand()
        print *, rand(seed), rand(), rand(), rand()
    end subroutine test_rand

    subroutine init()
        implicit none
        integer i,j,k,count_number
        real(8) radius,lz,distant

        open(12,file='dump.cylinder.lammpstrj')
        write(12,*)'ITEM:TIMESTEP'
        write(12,'(I9)')0
        write(12,*)'ITEM:NUMBER OF ATOMS'
        write(12,'(I6)')40300

        write(12,*)'ITEM:BOX BOUNDS'
        write(12,'(2F7.1)')-radius,radius
        write(12,'(2F7.1)')-radius,radius
        write(12,'(2F7.1)')-n_cell_z/2,n_cell_z/2.0
        write(12,*)'ITEM:ATOMS id type x y z'

        !!!!!!!!!以下是cylinder channel初值!!!!!!!!!!!!!!
        n_b=0
        k=nint(2.0*pi*radius)
        do j=0,2*n_cell_z-1
            do i=0,2*k-1
                n_b=n_b+1
                if(mod(j,2)==0)then
                    x_b(1,n_b)=radius*cos(dble(i)*pi/dble(k))
                    x_b(2,n_b)=radius*sin(dble(i)*pi/dble(k))
                else
                    x_b(1,n_b)=radius*cos(dble(i)*pi/dble(k)+pi/dble(2*k))
                    x_b(2,n_b)=radius*sin(dble(i)*pi/dble(k)+pi/dble(2*k))
                endif
                x_b(3,n_b)=dble(j-n_cell_z)*0.5
                write(12,'(2I6,3F13.4)') k,1,x_b(:,n_b)
            enddo
        enddo

        write(*,*)"cylinder单体数目："，n_b

        !!!!!!!!!以下是polymer chain初值!!!!!!!!!!!!!!
        x_p(1,1)=0.0
        x_p(2,1)=0.0
        x_p(3,1)=-5.0

        dx(1,1)=0
        dx(2,1)=1.0

        dx(1,2)=0
        dx(2,2)=-1.0

        dx(1,3)=1.0
        dx(2,3)=0

        dx(1,4)=-1.0
        dx(2,4)=0

        do i=2,n_p

            print*,i
            count_number=0
            call random_seed()
            call random_number(j)
            200  h=1+int(4*rand(j))

            if(count_number>11116)then

                x_p(1,i)=x_p(1,i)+0
                x_p(2,i)=x_p(2,i)+0
                x_p(3,i)=x_p(3,i)+1.0

                print*,(x_p(1,i)**2+x_p(2,i)**2)**(0.5)
                cycle
            endif

            if(mod(i,5)==0)then
                x_p(1,i)=x_p(1,i-1)+0
                x_p(2,i)=x_p(2,i-1)+0
                x_p(3,i)=x_p(3,i-1)+1.0

            else
                x_p(1,i)=x_p(1,i-1)+dx(1,h)
                x_p(2,i)=x_p(2,i-1)+dx(2,h)
                x_p(3,i)=x_p(3,i-1)

            endif


            if(x_p(3,i)>n_cell_z/2.0.or.x_p(3,i)<-n_cell_z/2.0)then
                x_p(3,i)=x_p(3,i)-nint(x_p(3,i)/n_cell_z)*n_cell_z
            endif

            r=sqrt(x_p(1,i)**2+x_p(2,i)**2)
            if(r>radius-1.0)then
                count_number=count_number+1
                goto 200
            endif

            if(i>1)then
                rz=x_p(3,i)-x_p(3,i-1)
                rz=rz-n_cell_z*nint(rz/n_cell_z)
                distant=((x_p(1,i)-x_p(1,(i-1)))**2+(x_p(2,i)-x_p(2,(i-1)))**2+rz**2)**(0.5)
                if(distant<0.9.or.distant>1.1)then
                    count_number=count_number+1
                    goto 200
                endif
            endif

            if(i>2)then
                do k=i-2,1,-1
                    rz=x_p(3,i)-x_p(3,k)
                    rz=rz-n_cell_z*nint(rz/n_cell_z)
                    distant=((x_p(1,i)-x_p(1,k))**2+(x_p(2,i)-x_p(2,k))**2+rz**2)**(0.5)
                    if(distant<1.1)then
                        exit
                    endif
                enddo

                if(k>0)then
                    count_number=count_number+1
                    goto 200
                endif
            endif

            distant=(x_p(1,i)**2+x_p(2,i)**2)**(0.5)

            if(distant>radius-1.0)then
                count_number=count_number+1
                goto 200
            endif

            write(output_file,'(2I6,3F13.4)') i,2,x_p(:,i)

        enddo

        write(*,*)"polymer单体数目："，n_p

        !!!!!!!!!!!!!!!!以下是solution粒子的初值!!!!!!!!!!!!!!

        l=0
        do i=0,int(n_cell_x)
            do j=0,int(n_cell_y)
                do k=0,int(n_cell_z)
                    l=l+1
                    x_s(1,l)=dble(i-int(n_cell_x)/2)*sid
                    x_s(2,l)=dble(j-int(n_cell_y)/2)*sid
                    x_s(3,l)=dble(k-int(n_cell_z)/2)*sid
                    write(output_file,'(2I6,3F13.4)') i,3,x_s(:,l)
                enddo
            enddo
        enddo
        write(*,*)"solvent单体数目："，n_s
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
        use parameters, only : x_p, f_p, n_p,f0_p
        implicit none
        integer mode, i, j
        real(8) U, temp(3)

        f_p=0
        U=0
        temp=0
        ! 链两端
        call FENE(f_p(:,1),U,x_p(:,1)-x_p(:,2))
        call FENE(f_p(:,n_p),U,x_p(:,n_p)-x_p(:,n_p-1))
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

            temp(2)=x_p(2,i)
            call LJ(f_p(:,i),U,temp)
            temp(2)=n_cell_y-x_p(2,i)
            call LJ(f_p(:,i),U,temp)
        enddo
        !$omp end parallel do

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

    subroutine cal_collision_velocity()
        use parameters
        implicit none
        integer ix,iy,iz,k,count_p,count_s,i
        real(8) momentum(3), matrix(3,3), l(3), fai, theta
        real(8), parameter :: alpha = 130*pi/180, s=sin(alpha), c=cos(alpha)
        real(8) v_aver_p(3), v_aver_s(3), v_aver(3), temp(3)
        logical mask_p(n_p), mask_p(n_s)
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
    subroutine scale_v(v,n,mass,Ek,T,T_out)
        implicit none
        integer n, i
        real(8), parameter :: Ek_fac = 1.5d0
        real(8) v(3,n), mass, Ek,T, scalar, Ek1, T_out, T1

        do i=1,3
            v(i,:) = v(i,:)-sum(v(i,:))/n
        enddo

        Ek1=0.5*mass*sum(v**2)
        T1=Ek1/(Ek_fac*n)
        scalar=sqrt(T/T1)
        v=v*scalar
        Ek=0.5*mass*sum(v**2)
        T_out=Ek1/(Ek_fac*n)
        write(*,*) '初始标度后动能Ek=',Ek
        write(*,*) '初始标度后温度T_out=',T_out
    end subroutine

    subroutine output(step)
        implicit none
        integer cur_step,output_interval_step,equilibrium_interval_step,output_file,energy_file
        !  energy_file=11
        !  output_file=12
        if(MOD(cur_step,output_interval_step)==0)then

            !   open(output_file,file='dump.CNT.lammpstrj')
            write(output_file,*)'ITEM:TIMESTEP'
            write(output_file,'(I9)')cur_step
            write(output_file,*)'ITEM:NUMBER OF ATOMS'
            write(output_file,'(I6)')n_b+n_p+n_s

            write(output_file,*)'ITEM:BOX BOUNDS'
            write(output_file,'(2F7.1)')-radius,radius
            write(output_file,'(2F7.1)')-radius,radius
            write(output_file,'(2F7.1)')-lz/2.0,lz/2.0
            write(output_file,*)'ITEM:ATOMS id type x y z'
            write(output_file,'(2I6,3F13.4)') k,1,x_b
            write(output_file,'(2I6,3F13.4)') k,2,x_p
            write(output_file,'(2I6,3F13.4)') k,3,x_s
            !write(13,'(3I6,3F13.4)') k,1,1,x1(k),y1(k),z1(k)
        endif

        if(MOD(cur_step,equilibrium_interval_step)==0)then
            write(energy_file,*) cur_step,U
        endif
    end subroutine

end module


program Poissonfield
    use parameters
    implicit none
    real(8) :: Ek, EK_scaled, T_Ek, T_scaled, density
    integer :: equili_step,equili_interval_step,total_step,output_interval_step,step,cur_step_per_rot, total_step_per_rot, &
        cur_step_pri, total_step_pri, cur_step, total_step

    integer i, ix, iy, iz, k, j
    real(8) randx, randz, ppx, ppz

    character(len=20) filename

    call random_seed()
    equili_step=500000
    equili_interval_step=10000
    total_step=3000000
    output_interval_step=10000
    total_rot_step=500000
    total_step_per_rot=200

    box_size = [n_cell_x, n_cell_y, n_cell_z]
    half_box_size = box_size/2
    density=0.85
    box_size_unit=(n_s/density)**(1d0/3)/n_cell_x
    half_box_size_unit=box_size_unit/2.0

    !!!d读链的大小 改成1个文件

    call test_rand

    call init()

    call random_number(randnum_group)
    x_s = x_s + (randnum_group-0.5)*box_size_unit
    do i=1,n_s
        x_s(:,k) = x_s(:,k) - box_size*nint((x_s(:,k)-half_box_size)/box_size)
    enddo


    open(70,file='qu0.out')
    open(80,file='qu.out')
    open(81,file='pcc.out')
    open(821,file='gxy.out')
    open(82,file='gzr.out')
    open(83,file='masscenter.out')
    open(84,file='yta.out')
    open(88,file='numy.out')
    open(90,file='r22.out')
    open(91,file='sf.out')
    open(100,file='vxyz1.out')
    open(110,file='pt.out')
    open(111,file='cc.out')
    open(200,file='vxyz0.out')
    open(210,file='xyz.out')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!assign initial velocity
    !!!链的初速度

    call random_number(v_p)
    v_p=v_p-0.5

    call scale_v(v_p, n_p, mass_p, Ek_scaled, T_set, T_scaled)

    write(*,*) '初始标度后动能kin2=',Ek_scaled
    write(*,*) '初始标度后温度tmp=',T_set

    !!!溶剂的初速度
    call random_number(v_s)
    v_s=v_s-0.5

    call scale_v(v_s, n_s, mass_s, Ek_scaled, T_set, T_scaled)

    write(*,*) '初始标度后动能kin2=',Ek_scaled
    write(*,*) '初始标度后温度tmp=',T_set

    call cal_average_momentum()

    ! 预处理
    !!! compute a(t-dt)
    call update_force(0)

    do cur_step_pri=1,total_step_pri
        x_p = x_p + v_p*time_step_p + 0.5*f0_p*time_step_p**2
        !!!!!! compute a(t+dt)
        call update_force(1)
        v_p = v_p + 0.5*(f0_p+f_p)*time_step_p
        call scale_v(v_p,n_p,mass_p,Ek,T_set,T_scaled)
        f0_p=f_p
        !call output(cur_step_pri)
    enddo

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
            call scale_v(v_p,n_p,mass_p,Ek,T_set,T_scaled)
            f0_p=f_p
            call output(cur_step_pri)

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

                call cal_average_momentum()

                call scale_v(v_p,n_p,mass_p,Ek,T_set,T_scaled)
                call scale_v(v_s,n_s,mass_s,Ek,T_set,T_scaled)

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


end program translocation
