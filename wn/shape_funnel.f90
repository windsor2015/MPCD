module shape_funnel
    use shape_all
    real(8) sum_v(0:100,0:400),sum_grid_v(2,0:100,0:400)
    integer n(0:100,0:400),n_grid(0:100,0:400)
contains

    subroutine report()
        write (*,*) 'funnel'
    end subroutine

    logical function in_pipe(x)
        implicit none
        real(8) x(3),r,theta

        r=sqrt(x(1)**2+x(2)**2)
        if (abs(x(3))<=n_cell_z*ratio_z/2d0) then
            in_pipe=r<=radius(2)
        elseif(abs(x(3))<n_cell_z*ratio_y/2d0.and.abs(x(3))>n_cell_z*ratio_z/2d0) then
            theta=(radius(1)-radius(2))/(n_cell_z*ratio_y/2d0-n_cell_z*ratio_z/2d0)
            in_pipe=r<=(abs(x(3))-n_cell_z*ratio_z/2d0)*theta + radius(2)
        else
            in_pipe=r<=radius(1)
        end if
    end function

    logical function in_pipe_phantom(x)
        implicit none
        real(8) x(3),r,theta
        real d
        d=sqrt(5d-1)
        r=sqrt(x(1)**2+x(2)**2)
        if (abs(x(3))<n_cell_z*ratio_z/2d0) then
            in_pipe_phantom=r<=radius(2)+d
        elseif(abs(x(3))<n_cell_z*ratio_y/2d0.and.abs(x(3))>n_cell_z*ratio_z/2d0) then
            theta=(radius(1)-radius(2))/(n_cell_z*ratio_y/2d0-n_cell_z*ratio_z/2d0)
            in_pipe_phantom=r<=(abs(x(3))-n_cell_z*ratio_z/2d0)*theta + radius(2) + d
        else
            in_pipe_phantom=r<=radius(1)+d
        end if
        in_pipe_phantom=in_pipe_phantom .and. (.not. in_pipe(x))
    end function

    real(8) function get_pipe_volume()
        implicit none
        get_pipe_volume=pi*radius(1)**2*(1d0-ratio_y)*n_cell_z+pi*radius(2)**2*ratio_z*n_cell_z+ &
            pi*n_cell_z*(ratio_y-ratio_z)*(radius(1)**2+radius(2)**2+radius(1)*radius(2))/6d0
    end function

    real(8) function get_pipe_phantom_volume()
        implicit none
        real d
        d=sqrt(5d-1)
        if(ratio_z>0d0)then
            get_pipe_phantom_volume=pi*(radius(1)+d)**2*(1d0-ratio_z)*n_cell_z+pi*(radius(2)+d)**2*ratio_z*n_cell_z+ &
                pi*n_cell_z*(ratio_y-ratio_z)*((radius(1)+d)**2+(radius(2)+d)**2+(radius(1)+d)*(radius(2)+d))/6d0-get_pipe_volume()
        elseif(ratio_z==0d0)then
            get_pipe_phantom_volume=pi*(radius(1)+d)**2*(1d0-ratio_z)*n_cell_z+pi*(radius(2)+d)**2*ratio_z*n_cell_z- &
                get_pipe_volume()
        end if

    end function

    integer function get_region(x)
        implicit none
        real(8) x(3),r,z,theta
        r=sqrt(x(1)**2+x(2)**2)
        z=x(3)
        theta=(radius(1)-radius(2))/(n_cell_z*ratio_y/2d0-n_cell_z*ratio_z/2d0)
        if (r>radius(1).or.(r>radius(2).and.abs(z)<n_cell_z*ratio_z/2d0).or. &
                (r>((abs(x(3))-n_cell_z*ratio_z/2d0)*theta + radius(2)).and.abs(z)<n_cell_z*ratio_y/2d0.and.abs(z)>n_cell_z*ratio_z/2d0))then
            get_region=0
            return
        endif
        if (ratio_z>0d0.and.z<=-n_cell_z*ratio_z/2d0)then
            get_region=2
            return
        endif
        if (ratio_z>0d0.and.z>=n_cell_z*ratio_z/2d0)then
            get_region=3
            return
        endif
        get_region=1
    end function

    subroutine do_ba(x,x0,v)
        implicit none
        real(8), dimension(3)::x,x0,xc,normal,v
        real(8) d,t
        ! 获取交点和法线
        call cross(x,x0,xc,normal)

        d=-dot_product(xc,normal)
        t=(dot_product(normal,x)+d)/norm2(normal)

        x=x-2*t*normal
        x0=xc+(x-xc)*1d-8
        v=norm2(v)/norm2(x-xc)*(x-xc)

        !return
        !这里没看太懂
        !        if(xc(3)/=0)then
        !            det=xc(1)**2+xc(2)**2
        !            if (det==0) cycle
        !            cc=x(1)*xc(1)
        !            ss=x(2)*xc(2)
        !            t=1-(cc+ss)/det
        !            x(1:2)=2d0*(x(1:2)+xc(1:2)*t)-x(1:2)
        !            ! x(3,i)=x1(3)
        !            !v1=(x(1:2,i)-xc(1:2))*norm2(v(1:2,i))/norm2(x(1:2,i)-xc(1:2))
        !            !v(1:2,i)=v1
        !        else
        !            x(3)=2*xc(3)-x(3)
        !        endif

    end subroutine

    subroutine cross(x,x0,xc,normal)
        implicit none
        real(8),dimension(3):: x, x0, xc, normal, xm
        real(8) t(4), a, b, c, r, delta, xc_list(3,4), normal_list(3,4), mint,z,z0,a0,b0,c0,cc,theta,t1,t2
        logical check(4)
        integer i
        !!越界与界面的交点
        ! 第1类，大圆交点
        normal_list=0
        t=-1
        xm=x-x0
        a=xm(1)**2+xm(2)**2
        b=2*x0(1)*xm(1)+2*x0(2)*xm(2)
        c=x0(1)**2+x0(2)**2-radius(1)**2

        delta=b**2-4*a*c
        if (delta>=0) t(1)=(-b+sqrt(delta))/(2*a)
        if(ratio_z==0d0)then
            xc=x0+xm*t(1)
            normal(1:2)=-xc(1:2)
        endif

        if(ratio_z>0d0)then

            theta=(radius(1)-radius(2))/n_cell_z/(ratio_y-ratio_z)*2d0
            z=radius(1)/theta!/(radius(1)-radius(2))
            z0=n_cell_z*ratio_y/2d0-z

            ! 第2类，右侧锥面
            a0=xm(1)**2+xm(2)**2-xm(3)**2*theta**2
            b0=2d0*(xm(1)*x0(1)+xm(2)*x0(2)-xm(3)*(x0(3)-z0)*theta**2)
            c0=x0(1)**2+x0(2)**2-(x0(3)-z0)**2*theta**2
            delta=b0**2-4*a0*c0
            ! 应取较小的正解
            if (delta>=0) then
                t1=(-b0-sqrt(delta))/(2*a0)
                t2=(-b0+sqrt(delta))/(2*a0)
                t(2)=min(t1,t2)
                if (t(2)<=0) t(2)=max(t1,t2)
            endif
            ! 第3类，左侧锥面
            b0=2d0*(xm(1)*x0(1)+xm(2)*x0(2)-xm(3)*(x0(3)+z0)*theta**2)
            c0=x0(1)**2+x0(2)**2-(x0(3)+z0)**2*theta**2
            delta=b0**2-4*a0*c0
            ! 应取较小的正解
            if (delta>=0) then
                t1=(-b0-sqrt(delta))/(2*a0)
                t2=(-b0+sqrt(delta))/(2*a0)
                t(3)=min(t1,t2)
                if (t(3)<=0) t(3)=max(t1,t2)
            endif

            ! 第4类，小圆交点
            c=x0(1)**2+x0(2)**2-radius(2)**2
            delta=b**2-4*a*c
            if (delta>=0) t(4)=(-b+sqrt(delta))/(2*a)

            do i=1,4
                xc_list(:,i)=x0+xm*t(i)
            enddo

            ! 3种法线
            normal_list(1:2,1)=xc_list(1:2,1)
            normal_list(1:2,4)=xc_list(1:2,4)
            normal_list(1:2,2:3)=xc_list(1:2,2:3)
            normal_list(3,2)=xc_list(3,2)-(xc_list(1,2)**2+xc_list(2,2)**2+(xc_list(3,2)-z0)**2)/(xc_list(3,2)-z0)-z0
            normal_list(3,3)=xc_list(3,3)-(xc_list(1,3)**2+xc_list(2,3)**2+(xc_list(3,3)-z0)**2)/(xc_list(3,3)-z0)-z0

            ! 3类边界的范围
            check(1)=abs(xc_list(3,1))>=n_cell_z*ratio_y/2d0

            check(2)=xc_list(3,2)<n_cell_z*ratio_y/2d0 .and. &
                xc_list(3,2)>n_cell_z*ratio_z/2d0

            check(3)=xc_list(3,3)>-n_cell_z*ratio_y/2d0 .and. &
                xc_list(3,3)<-n_cell_z*ratio_z/2d0

            check(4)=abs(xc_list(3,4))<=n_cell_z*ratio_z/2d0

            mint=1
            do i=1,4
                if (t(i)>=0 .and. t(i)<1 .and. check(i) .and. t(i)<mint) then
                    xc=xc_list(:,i)
                    normal=normal_list(:,i)
                    mint=t(i)
                endif
            enddo
            if(debug==1)then
                write(*,*) 'x0,x,xc,normal'
                write(*,*) x0,x,xc,normal
                !write(*,*) xc_list
                write(*,*) 't',t
                write(*,*)
            end if

        endif
        !return
        !    call cross_border(x,check_cylinder,check_plane)
        !    if(check_cylinder)then
        !        if(sqrt(x(1)**2+x(2)**2)>radius(1))then
        !            r=radius(1)
        !        elseif(sqrt(x(1)**2+x(2)**2)>radius(2).and.sqrt(x(1)**2+x(2)**2)<radius(2))then
        !            r=radius(2)
        !        endif
        !        xm=x(1:2)-x_0(1:2)
        !        ! solve equation sum((t*xm-x0)**2)=radius**2
        !        ! 对于a*t^2+b*t+c=0，c必定小于0，因此解必有一正一负，仅取正值
        !        a=xm(1)**2+xm(2)**2
        !        b=2*x_0(1)*xm(1)+2*x_0(2)*xm(2)
        !        c=x_0(1)**2+x_0(2)**2-r**2
        !        delta=b**2-4*a*c
        !        !if (delta<0) cycle
        !        t=(-b+sqrt(delta))/2/a
        !        !if (t<0 .or. t>1) cycle
        !        ! 旋转，最好参看配图
        !        xc(1:2)=x_0(1:2)+(x(1:2)-x_0(1:2))*t
        !        xc(3)=(xc(1)*x(3)-x_0(1)*x(3)+x(1)*x_0(3)-xc(1)*x_0(3))/(x(1)-x_0(1))
        !    elseif(check_plane)then
        !        xc(3)=0d0
        !    end if
    end subroutine

    subroutine bounce_back(x,x0,v,n)
        implicit none
        integer i,n,c
        real(8), dimension(3,n):: x, x0, v

        !$omp parallel do private(c)
        do i=1,n
            ! 越界则回弹
            if (get_region(x(:,i))/=get_region(x0(:,i))) then
                c=0
                debug=0
                !if (i==3843) debug=1
                do while(.not. in_pipe(x(:,i)))
                    call do_ba(x(:,i),x0(:,i),v(:,i))
                    c=c+1
!                    if (debug==1) then
!                        write(*,*) x(:,i)
!                    endif
!                    if (c>2) then
!                        write(*,*)i
!                        debug=1
!                        if (c>=10) stop
!                    end if
                enddo
            endif
        enddo
        !$omp end parallel do
    end subroutine

    subroutine clear_stat()
        sum_v=0
        n=0
        sum_grid_v=0
        n_grid=0
    end subroutine

    subroutine stat_velocity(cur_step,interval_step)
        implicit none
        integer cur_step,interval_step
        integer i,j,k
        real(8) r
        if(mod(cur_step,1)==0)then
            do i=1,n_s
                if(abs(x_s(2,i))<0.25)then
                    j=floor(x_s(1,i)*2)+floor(radius(1)*2)
                    k=floor(x_s(3,i)*2)+n_cell_z
                    sum_grid_v(1,j,k)=sum_grid_v(1,j,k)+v_s(1,i)
                    sum_grid_v(2,j,k)=sum_grid_v(2,j,k)+v_s(3,i)
                    n_grid(j,k)=n_grid(j,k)+1
                end if
            enddo
        endif

        if(mod(cur_step,interval_step)==0)then
            do k=0,n_cell_z-1
                do i=1,n_s
                    if(x_s(3,i)<(k-n_cell_z/2+1)*1d0 .and. x_s(3,i)>=(k-n_cell_z/2)*1d0)then
                        r=sqrt(x_s(1,i)**2+x_s(2,i)**2)
                        j=floor(r*5)
                        sum_v(j,k)=sum_v(j,k)+v_s(3,i)
                        n(j,k)=n(j,k)+1
                    end if
                enddo
            enddo
        end if
    end subroutine

    subroutine output_velocity(num,velocity_file,coord_velo_file,step,interval_step)
        implicit none
        integer num, velocity_file,coord_velo_file,step,interval_step
        integer k,j
        do k=0,n_cell_z-1
            do j=0,floor(radius(1)*5)-1
                sum_v(j,k)=sum_v(j,k)/(step/interval_step)
                write(velocity_file,'(3I6,ES18.4)')num,k,j,sum_v(j,k)/n(j,k)
            enddo
        enddo

        do j=0,floor(radius(1)*4)-1
            do k=0,2*n_cell_z-1
                sum_grid_v(:,j,k)=sum_grid_v(:,j,k)/(step/1)
                write(coord_velo_file,'(3I6,2ES18.4)')num,j,k,sum_grid_v(:,j,k)/n_grid(j,k)
            enddo
        enddo
    end subroutine

end module

