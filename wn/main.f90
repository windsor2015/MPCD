module parameters
    implicit none
    ! 下标约定: p - particle, s - solution

    !结构
    integer, parameter :: n_p=40, n_cell_x=15, n_cell_y=15, n_cell_z=15, n_unit_p=10, n_s=n_unit_p*n_cell_x*n_cell_y*n_cell_z

    integer, parameter :: n_b = (n_cell_x+1)*(n_cell_y+1)*(n_cell_z+1)

    !n是单体数目，n0是单体、壁、孔数目之和

    real(8), parameter :: kB = 1.38064852d-23
    real(8), parameter :: pi = 3.141592653589793238462643383279

    real(8), parameter :: time_step_p=0.0001, time_step_s=0.0001, mass_p=1, mass_s=0.1, T_set=1, v_probability=0.2


    ! polymer 位置 速度 力 上一次力
    real(8), dimension(3,n_p) :: x_p, v_p, f_p, f0_p
    ! solution
    real(8), dimension(3,n_s) :: x_s, v_s, f_s, x0_s
    ! boundaaries, 1~nb-up, nb+1~2nb-down
    real(8), dimension(3,2*n_b) :: x_b, v_b, f_b

    real(8) box_size_unit, half_box_size_unit, box_size(3), half_box_size(3)


    real(8) randnum_group(3,n_s)

    real(8) aver_v(3)


    !    real(8) x_bu(n_wall),y_bu(n_wall),z_bu(n_wall),x_bd(n_wall),y_bd(n_wall),z_bd(n_wall)  ! 上下边界坐标, b-boundary, u-up, d-down
    !    real(8) y_p(n_p),z_p(n_p),vr_p(n_p),vy_p(n_p),vz_p(n_p),r0_s(n_s),y11(n_s),z11(n_s)
    !    real(8) dx1(n_s),dy1(n_s),dz1(n_s),x12(n_s),y12(n_s),z12(n_s),y13(n_s)
    !    real(8) dx2(n_s),dy2(n_s),dz2(n_s),x0(n_s),y0(n_s),z0(n_s),x10(n_s),y10(n_s),z10(n_s)
    !    real(8) y_s(n_s),z_s(n_s),vr_s(n_s),vy_s(n_s),vz_s(n_s)
    !    real(8) vxx(4500),vyy(4500),vzz(4500),d1(n_p),dd0,dd1,dd
    !    real(8) vx10(n_s),vy10(n_s),vz10(n_s),g(n_s),gee,vu(n_s),vvu(n_s),Ek,U,ppx,ppy,ppz,x30(n_p),y30(n_p)
    !    real(8) forcex1_p(n_p),forcey1_p(n_p),forcez1_p(n_p),forcex2_p(n_p),forcey2_p(n_p),forcez2_p(n_p),kx3(n_p),ky3(n_p),kz3(n_p),f
    !    real(8) rannx(n_p),ranny(n_p),rannz(n_p),ran1,ran2,rco,lx,ly,lz
    !    real(8) vp,scalx,scaly,scalar_T,scalz,vp2,ranx,rany,ranz,rx,ry,rz,rr,rr1,rl
    !    real(8) sum_x,sum_y,sum_z,mx1,my1,mz1,Ek_real,Ek_scaled,kin11,kin22,Ek0,Ek1,Ek2
    !    real(8) rx1,ry1,rz1,rx2,ry2,rz2,rx3,ry3,rz3,r1,r2,r3,alpha,cos_alpha,sin_alpha
    !    real(8) v_probability,dlty,R22,xcm,ycm,zcm,S2mm(n_p),gg(n_p),ge
    !    real(8) x_up1(n_wall),y_up1(n_wall),z_up1(n_wall),x_down1(n_wall),y_down1(n_wall),z_down1(n_wall)
    !    real(8) rx00,ry00,rz00,rxyz,ddx1,ddy1,ddz1,ddx2,ddy2,ddz2,ddx3,ddy3,ddz3
    !    real(8) dt,sf,Ek_fac,dx0(n_s),dy0(n_s),dz0(n_s),theta,fai
    !    real(8) rrx1,rry1,rrx2,rry2,rrz2,rr2,rrx3,rry3,rrz3,rr3,rr4
    !    real(8) t0,t1,t2,cc,ss,sss,c1,c2,c3,c12,c22,c32,pc,pcc,pyx,yta
    !    real(8) s(3,3),d(3),vv(n_p),vvv(n_p),temp1,T_Ek,r,ll
    !    real(8) LJ_rb,FENE_max_distance,ld1,LJ_rc    !,   small
    !    integer count_s,count_p,j,sub,jj,p,n0,l,t,q1,cur_step_pri,cur_step,qq,total_step_pri,total_step,st01
    !    integer st0,st00,st1,st11,step_per_run,st22,a,n2,rand1_num,gama,num(n_cell_y)
end module



program translocation
    use parameters
    implicit none
    real(8) :: Ek, EK_scaled, T_Ek, T_scaled
    integer :: cur_step_per_rot, total_step_per_rot, cur_step_pri, total_step_pri, cur_step, total_step

    integer i, ix, iy, iz, k, j
    real(8) randx, randz, ppx, ppz

    character(len=20) filename

    call random_seed()
    total_step_pri=3800000
    !    st1=3800000
    !    st11=10000
    total_step=500000
    !    st0=199000
    !    st01=1000
    !    st00=200000
    !    step_per_run=1000
    !    st22=10000
    total_step_per_rot=200
    !    alpha=130*pi/180
    !
    !    rand1_num=-4500
    !
    !
    !    temp1=1

    box_size = [n_cell_x, n_cell_y, n_cell_z]
    half_box_size = box_size/2

    box_size_unit=box_size(1)/n_cell_x
    half_box_size_unit=box_size_unit/2.0


    !cos_alpha=cos(alpha)
    !sin_alpha=sin(alpha)
    !small=dble(0.001)
    !!!d读链的大小 改成1个文件

    open(11,file='1.txt')
    do i=1,n_p
        read(11,*) x_p(:,i)
    end do
    close(11)

    !!!!!!!!!!!!!!!!!建立初始平板边界
    ! 这个好像后面并没用上
    !    do ix=1,n_cell_x+1
    !        do iz=1,n_cell_z+1
    !            sub=ix+(iz-1)*(n_cell_x+1)
    !            x_bu(sub)=(i-1)*box_size_unit
    !            z_bu(sub)=(j-1)*box_size_unit
    !            x_bd(sub)=(i-1)*box_size_unit
    !            z_bd(sub)=(j-1)*box_size_unit
    !        enddo
    !    enddo
    !    y_bu(:)=n_cell_y
    !    y_bd(:)=0

    !OPEN(2,FILE='rongji.xyz')
    k=0
    do ix=0,n_cell_x-1
        do iy=0,n_cell_y-1
            do iz=0,n_cell_z-1
                do i=1,n_unit_p
                    k=k+1
                    x_s(:,k) = [iz,iy,ix] * box_size_unit + i*box_size_unit/(n_unit_p+1)
                enddo
            enddo
        enddo
    enddo

    call random_number(randnum_group)
    x_s = x_s + (randnum_group-0.5)*box_size_unit
    do i=1,n_s
        x_s(:,k) = x_s(:,k) - box_size*nint((x_s(:,k)-half_box_size)/box_size)
    enddo

    !OPEN(2,file='2.xyz')
    !write(2,*) mn+512+n
    !write(2,*)'cohar'
    !			   do i=1,n
    !            	write(2,'(a4,f17.7,f13.7,f13.7)')'p',x(i),y(i),z(i)
    !               enddo
    !			   do i=1,256
    !            	write(2,'(a4,f17.7,f13.7,f13.7)')'o',x_up(i),y_up(i),z_up(i)
    !               enddo
    !               do i=1,256
    !            	write(2,'(a4,f17.7,f13.7,f13.7)')'o',x_down(i),y_down(i),z_down(i)
    !               enddo
    !       do i=1,mn
    !        write(2,'(a4,f17.7,f13.7,f13.7)') 'f',x1(i),y1(i),z1(i)
    !	   end do
    !CLOSE(2)

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
            x_p(1,:) = x_p(1,:)-0.1*(x_p(2,:)**2-15.0*x_p(2,:))*v_probability*time_step_p
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
                x_s(1,:) = x_s(1,:) -0.1*(x_s(2,:)**2-15.0*x_s(2,:))*v_probability*time_step_s
                ! 弹回, 看起来不对
                do i=1,n_s
                    !y12(i)=y_s(i)
                    !if(y_s(i)>=box_size_y.or.y_s(i)<=0)then
                    !y_s(i)=y12(i)
                    !endif
                enddo

                x_s(1,:)=x_s(1,:)+(randx-0.5)*box_size_unit
                x_s(3,:)=x_s(3,:)+(randz-0.5)*box_size_unit

                x_s(1,:)=x_s(1,i)-box_size(1)*nint((x_s(1,:)-(ppx+half_box_size(1)))/box_size(1))
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

contains

    subroutine tran(q,string)
        implicit none
        integer q,i,j,temp,k
        character(len=20) string
        do i=1,20
            string(i:i)=' '
        enddo
        temp=q
        i=1
        do while(temp/10/=0)
            i=i+1
            temp=temp/10
        end do
        temp=q
        do j=i,1,-1
            k=mod(temp,10)
            temp=temp/10
            string(j:j)=achar(k+48)
        end do
        string((i+1):(i+4))=".xyz"
        string=trim(string)
        return
    end subroutine tran

    subroutine FENE(f,U,rx)
        implicit none
        real(8) f(3), U, rx(3)
        real(8), parameter :: FENE_rc=1.5, FENE_k=30
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
        real(8), parameter :: LJ_sigma = 1, LJ_rc=LJ_sigma*2.0**(1.0/6.0)
        real(8), parameter :: LJ_epsilon = 30
        real(8) temp,r

        r=norm2(rx)
        if(r/=0 .and. r<=LJ_rc)then
            temp=48.0*LJ_epsilon*((LJ_sigma/r)**6**2)/(r**2)
            f=f+rx*temp
            U=U+4*LJ_epsilon*(LJ_sigma/r)**12
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
            if (i<n_p) then
                do j=i+1,n_p
                    call LJ(f_p(:,i),U,x_p(:,i)-x_p(:,j))
                    call LJ(f_p(:,j),U,x_p(:,j)-x_p(:,i))
                enddo
            endif
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
        real(8) x(3)
        integer ix,iy,iz

        inbox=.true.
    end function

    subroutine cal_average_momentum()
        use parameters
        implicit none
        integer ix,iy,iz,k,count_p,count_s,i,j
        real(8) momentum(3), matrix(3,3), aver_momentum(3),l(3),suml,fai,theta
        real(8), parameter :: alpha = 130*pi/180, s=sin(alpha), c=cos(alpha)
        ! calculate average momentum of each cell
        k=0
        do ix=0,n_cell_x-1
            do iy=0,n_cell_y-1
                do iz=0,n_cell_x-1

                    k=k+1

                    fai=2.0*pi*rand(0)
                    theta=2.0*rand(0)-1
                    !c=cos(alpha)
                    !s=sin(alpha)
                    l(1)=cos(fai)*SQRT(1-theta**2)
                    l(2)=sin(fai)*SQRT(1-theta**2)
                    l(3)=theta
                    suml=sum(l)
                    matrix=0

                    do i=1,3
                        do j=1,3
                            matrix(i,j)=l(i)*l(j)
                            if (i/=j) then
                                if (i<j) then
                                    matrix(i,j)=matrix(i,j)*(1-c)+(-1)**(i+j)*(suml-l(i)-l(j))*s
                                else
                                    matrix(i,j)=matrix(i,j)*(1-c)-(-1)**(i+j)*(suml-l(i)-l(j))*s
                                endif
                            else
                                matrix(i,j)=matrix(i,j)+(1-l(i)**2)*c
                            endif

                        enddo
                    enddo

                    momentum=0
                    count_p=0
                    do i=1,n_p
                        if(inbox(x_p(:,i),ix,iy,iz))then
                            momentum = momentum + v_p(:,i)*mass_p
                            count_p=count_p+1
                            v_p(:,i)=matmul(matrix,v_p(:,i))
                        endif
                    enddo
                    count_s=0
                    do i=1,n_s
                        if(inbox(x_s(:,i),ix,iy,iz))then
                            momentum = momentum + v_s(:,i)*mass_s
                            count_s=count_s+1
                            v_s(:,i)=matmul(matrix,v_s(:,i))
                        endif
                    enddo

                    if((count_s==0).and.(count_p==0))then
                        aver_momentum=0
                    else
                        aver_momentum=momentum/(count_s*mass_s+count_p*mass_p)
                    endif
                enddo
            enddo
        enddo
    end subroutine

    subroutine scale_v(v,n,mass,Ek,T,T_out)
        implicit none
        integer n, i
        real(8), parameter :: Ek_fac = 1.5d0
        real(8) v(3,n),mass,Ek,T,scalar,Ek1,T_out,T1

        do i=1,3
            v(i,:) = v(i,:)-sum(v(i,:))/n
        enddo
        Ek1=0.5*mass*sum(v**2)
        T1=Ek1/(Ek_fac*n)
        scalar=sqrt(T/T1)
        v=v*scalar
        Ek=0.5*mass*sum(v**2)
        T_out=Ek1/(Ek_fac*n)

    end subroutine

    subroutine output(step)
        implicit none
        integer step
        !        if(MOD(cur_step_pri,st1)==0)then
        !
        !            call tran(cur_step_pri+5,filename)
        !            open(31,file=filename)
        !            write(31,*) n_s+512+n_p
        !            write(31,*)'cohar'
        !
        !            do i=1,n_p
        !                write(31,'(a4,f17.7,f13.7,f13.7)')'p',r_p(i),y_p(i),z_p(i)
        !            enddo
        !
        !            do i=1,512
        !                write(31,'(a4,f17.7,f13.7,f13.7)')'o',x_bu(i),y_bu(i),z_bu(i)
        !            enddo
        !
        !            do i=1,512
        !                write(31,'(a4,f17.7,f13.7,f13.7)')'o',x_bd(i),y_bd(i),z_bd(i)
        !            enddo
        !            do i=1,n_s
        !                write(31,'(a4,f17.7,f13.7,f13.7)') 'f',r_s(i),y_s(i),z_s(i)
        !            end do
        !            close(31)
        !
        !        endif
        !
        !        if(MOD(cur_step_pri,st11)==0)then
        !            write(70,*) cur_step_pri,U1
        !        endif
    end subroutine

end program translocation
