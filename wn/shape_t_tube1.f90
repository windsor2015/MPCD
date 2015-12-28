module shape_t_tube1
    use shape_all
    integer n(0:100,0:400),n_grid(0:100,0:100,0:400)
    real(8) sum_v(3,0:100,0:400),sum_grid_v(3,0:100,0:100,0:400)
contains

    subroutine report()
        write (*,*) 't-tube-1'
    end subroutine

    logical function in_pipe(x)
        implicit none
        real(8) x(3),z
        z=x(3)-n_cell_z*nint((x(3))/n_cell_z)
        in_pipe=.false.
        if (abs(z)<=n_cell_z*ratio_z/2d0) then
            in_pipe=x(2)<=n_cell_y/2d0 .and. x(2)>=0 .and. abs(x(1))<=n_cell_x/2d0
        elseif ((z<=n_cell_z/2d0 .and. z>n_cell_z*ratio_z/2d0).or. &
                (z<-n_cell_z*ratio_z/2d0 .and. z>=-n_cell_z/2d0)) then
            in_pipe=x(2)<=n_cell_y/2d0 .and. x(2)>=(1d0-ratio_y)*n_cell_y/2d0 .and. abs(x(1))<=n_cell_x/2d0
        end if
    end function

    logical function in_pipe_phantom(x)
        implicit none
        real(8) x(3),r
        real d
        d=sqrt(5d-1)

        if (abs(x(3))<=n_cell_z*ratio_z/2d0+d) then
            in_pipe_phantom=x(2)<=n_cell_y/2d0+d.and.x(2)>=-d

        elseif ((x(3)<=n_cell_z/2d0.and.x(3)>n_cell_z*ratio_z/2d0+d).or. &
                (x(3)<-n_cell_z*ratio_z/2d0-d.and.x(3)>=-n_cell_z/2d0)) then

            in_pipe_phantom=x(2)<=n_cell_y/2d0+d.and.x(2)>=(1d0-ratio_y)*n_cell_y/2d0-d

        end if

        in_pipe_phantom=in_pipe_phantom .and. (.not. in_pipe(x))
    end function

    real(8) function get_pipe_volume()
        implicit none
        get_pipe_volume=n_cell_x*n_cell_y/2d0*n_cell_z*ratio_z+ &
            n_cell_x*ratio_y*n_cell_y/2d0*n_cell_z*(1d0-ratio_z)
    end function

    real(8) function get_pipe_phantom_volume()
        implicit none
        real d
        d=sqrt(5d-1)
        get_pipe_phantom_volume=(n_cell_x+2d0*d)*(n_cell_y/2d0+2d0*d)*(n_cell_z*ratio_z+2d0*d)+ &
            (n_cell_x+2d0*d)*(ratio_y*n_cell_y/2d0+2d0*d)*(n_cell_z*(1-ratio_z)-2d0*d)-get_pipe_volume()
    end function

    integer function get_region(x)
        implicit none
        real(8) x(3),x_,y,z
        x_=x(1)
        y=x(2)
        z=x(3)-n_cell_z*nint((x(3))/n_cell_z)
        if (x_>n_cell_x/2d0.or.x_<-n_cell_x/2d0.or.y>n_cell_y/2d0.or.(y<0.and.abs(z)<=n_cell_z*ratio_z/2d0).or. &
                (y<n_cell_y*(1d0-ratio_y)/2d0.and.(z<-n_cell_z*ratio_z/2d0.or.z>n_cell_z*ratio_z/2d0)))then
            get_region=0
            return
        endif
        if (z<-n_cell_z*ratio_z/2d0)then
            get_region=1
            return
        endif
        if (z>n_cell_z*ratio_z/2d0)then
            get_region=3
            return
        endif
        get_region=2
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
        real(8) t(7), a, b, c, r, delta, xc_list(3,7), normal_list(3,7), mint, xc_list2(3,7)
        logical check(7)
        integer i
        !!越界与界面的交点
        t=-1
        xm=x-x0
        ! 第1类，x2=n_cell_y/2d0
        t(1)=(n_cell_y/2d0-x0(2))/xm(2)

        ! 第2类，x2=(1d0-ratio_y)*n_cell_y/2d0
        t(2)=((1d0-ratio_y)*n_cell_y/2d0-x0(2))/xm(2)

        ! 第3类，x2=0
        t(3)=(0-x0(2))/xm(2)

        ! 第4类，x3=-n_cell_z*ratio_z/2d0
        t(4)=(-n_cell_z*ratio_z/2d0-x0(3))/xm(3)

        ! 第5类，x3=n_cell_z*ratio_z/2d0
        t(5)=(n_cell_z*ratio_z/2d0-x0(3))/xm(3)

        ! 第6类，x1=-n_cell_x/2d0
        t(6)=(-n_cell_x/2d0-x0(1))/xm(1)

        ! 第7类，x1=n_cell_x/2d0
        t(7)=(n_cell_x/2d0-x0(1))/xm(1)

        do i=1,7
            xc_list(:,i)=x0+xm*t(i)
        enddo
        xc_list2=xc_list

        xc_list2(3,:)=xc_list2(3,:)-n_cell_z*nint((xc_list2(3,:))/n_cell_z)

        ! 7种法线
        normal_list=0
        normal_list(2,1)=1
        normal_list(2,2)=-1
        normal_list(2,3)=-1
        normal_list(3,4)=-1
        normal_list(3,5)=1
        normal_list(1,6)=-1
        normal_list(1,7)=1

        ! 7类边界的范围
        check(1)=abs(xc_list2(3,1))<=n_cell_z/2d0 .and. abs(xc_list2(1,1))<=n_cell_x/2d0

        check(2)=((xc_list2(3,2)>=-n_cell_z/2d0 .and. xc_list2(3,2)<-n_cell_z*ratio_z/2d0) .or. &
            (xc_list2(3,2)<=n_cell_z/2d0 .and. xc_list2(3,2)>n_cell_z*ratio_z/2d0)) .and. abs(xc_list2(1,2))<=n_cell_x/2d0

        check(3)=abs(xc_list2(3,3))<=n_cell_z*ratio_z/2d0 .and. abs(xc_list2(1,3))<=n_cell_x/2d0

        check(4:5)=(xc_list2(2,4:5)<=n_cell_y*(1d0-ratio_y)/2d0 .and. xc_list2(2,4:5)>0) .and. abs(xc_list2(1,4:5))<=n_cell_x/2d0

        check(6:7)=(((xc_list2(3,6:7)>=-n_cell_z/2d0 .and. xc_list2(3,6:7)<-n_cell_z*ratio_z/2d0) .or. &
            (xc_list2(3,6:7)<=n_cell_z/2d0 .and. xc_list2(3,6:7)>n_cell_z*ratio_z/2d0)) .and. &
            (xc_list2(2,6:7)>=n_cell_y*(1d0-ratio_y)/2d0 .and. xc_list2(2,6:7)<=n_cell_y/2d0)) .or. &
            (abs(xc_list2(3,6:7))<=n_cell_z*ratio_z/2d0 .and.(xc_list2(2,6:7)<=n_cell_y/2d0 .and. xc_list2(2,6:7)>=0))


        mint=1
        do i=1,7
            if (t(i)>=0 .and. t(i)<1 .and. check(i) .and. t(i)<mint) then
                xc=xc_list(:,i)
                normal=normal_list(:,i)
                mint=t(i)
            endif
        enddo

!        if (debug==1)then
!            write(*,*) 't',t
!            write(*,*) 'min_t',mint,check
!        end if

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
        integer i,n,c, cur_step

        real(8), dimension(3,n):: x, x0, v

        !$omp parallel do private(c)
        do i=1,n
            ! 越界则回弹
            if (get_region(x(:,i))/=get_region(x0(:,i))) then
                c=0
                !                debug=0
                do while(get_region(x(:,i))==0)
                    !                    if (i==1979) debug=1
                    !                    if (debug==1)then
                    !                        write(*,*)x0(:,i)
                    !                        write(*,*)x(:,i)
                    !                        debug=1
                    !                    endif
                    call do_ba(x(:,i),x0(:,i),v(:,i))
                    c=c+1
                    !                    if (c>2 .or. debug==1) then
                    !                        write(*,*) 'c',c,i,cur_step
                    !                        write(*,*) 'x',x(:,i)
                    !                        if (c>10) stop
                    !                        debug=1
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
        integer i,j,k,l

        if(mod(cur_step,1)==0)then
            do i=1,n_s
                !if (.not. in_pipe(x_s(:,i))) cycle
                !if(x_s(2,i)>-0.25.and.x_s(2,i)<0.25)then
                l=floor(x_s(1,i)*2+n_cell_x)
                j=floor(x_s(2,i)*2+n_cell_y)
                k=floor(x_s(3,i)*2+n_cell_z)
                sum_grid_v(1,k,j,l)=sum_grid_v(1,k,j,l)+v_s(1,i)
                sum_grid_v(2,k,j,l)=sum_grid_v(2,k,j,l)+v_s(2,i)
                sum_grid_v(3,k,j,l)=sum_grid_v(3,k,j,l)+v_s(3,i)
                n_grid(k,j,l)=n_grid(k,j,l)+1
                !if (j==33.and.k==74) write(*,*) x_s(:,i),v_s(:,i), i
                !end if
                !write(*,*)j,k
            enddo
        endif

    end subroutine

    subroutine output_velocity(num,velocity_file,coord_velo_file,step,interval_step)
        implicit none
        integer num, velocity_file,coord_velo_file,step,interval_step
        integer k,j,l

        do l=0,2*n_cell_x-1
            do j=n_cell_y,2*n_cell_y-1
                do k=0,2*n_cell_z-1
                    sum_grid_v(:,k,j,l)=sum_grid_v(:,k,j,l)/(step/1)
                    if (n_grid(k,j,l)==1) then
                        n_grid(k,j,l)=0
                        sum_grid_v(:,k,j,l)=0
                    end if
                    write(coord_velo_file,'(4I6,3ES18.4)')num,l,j,k,sum_grid_v(:,k,j,l)/n_grid(k,j,l)
                    !write(*,*)sum_grid_v(:,j,k),n_grid(j,k)
                enddo
            enddo
        enddo
    end subroutine

end module
