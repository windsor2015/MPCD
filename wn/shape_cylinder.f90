module shape_cylinder
    use shape_all
    real(8), dimension(2)::radius=(/4, 2/)

contains

    subroutine report()
        write (*,*) 'cylinder'
    end subroutine

    logical function in_pipe(x)
        implicit none
        real(8) x(3),r

        r=sqrt(x(1)**2+x(2)**2)
        if (abs(x(3))<n_cell_z*ratio_z/2d0) then
            in_pipe=r<=radius(2)
        else
            in_pipe=r<=radius(1)
        end if
    end function

    logical function in_pipe_phantom(x)
        implicit none
        real(8) x(3),r
        real d
        d=sqrt(5d-1)
        r=sqrt(x(1)**2+x(2)**2)
        if (abs(x(3))<n_cell_z*ratio_z/2d0-d) then
            in_pipe_phantom=r<=radius(2)+d
        else
            in_pipe_phantom=r<=radius(1)+d
        end if
        in_pipe_phantom=in_pipe_phantom .and. (.not. in_pipe(x))
    end function

    real(8) function get_pipe_volume()
        implicit none
        get_pipe_volume=pi*radius(1)**2*(1d0-ratio_z)*n_cell_z+pi*radius(2)**2*ratio_z*n_cell_z
    end function

    real(8) function get_pipe_phantom_volume()
        implicit none
        real d
        d=sqrt(5d-1)
        if(ratio_z>0)then
            get_pipe_phantom_volume=pi*(radius(1)+d)**2*(1d0-ratio_z)*n_cell_z+pi*(radius(2)+d)**2*ratio_z*n_cell_z+ &
                (2*pi*((radius(1)+d)**2-(radius(2)+d)**2)*d)-get_pipe_volume()
        elseif(ratio_z==0)then
            get_pipe_phantom_volume=pi*(radius(1)+d)**2*(1d0-ratio_z)*n_cell_z+pi*(radius(2)+d)**2*ratio_z*n_cell_z- &
                get_pipe_volume()
        end if

    end function

    integer function get_region(x)
        implicit none
        real(8) x(3),r,z
        r=sqrt(x(1)**2+x(2)**2)
        z=x(3)
        if (r>radius(1).or.(r>radius(2).and.r<radius(1).and.abs(z)<n_cell_z*ratio_z/2d0))then
            get_region=0
            return
        endif
        if (ratio_z>0.and.z<=-n_cell_z*ratio_z/2d0)then
            get_region=2
            return
        endif
        if (ratio_z>0.and.z>=n_cell_z*ratio_z/2d0)then
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
        real(8) t(4), a, b, c, r, delta, xc_list(3,4), normal_list(3,4), mint
        logical check(4)
        integer i
        !!越界与界面的交点
        ! 第1类，大圆交点

        t=-1
        xm=x-x0
        a=xm(1)**2+xm(2)**2
        b=2*x0(1)*xm(1)+2*x0(2)*xm(2)
        c=x0(1)**2+x0(2)**2-radius(1)**2

        delta=b**2-4*a*c
        t(1)=(-b+sqrt(delta))/(2*a)

        ! 第2类，x3=-5
        t(2)=(-5-x0(3))/xm(3)

        ! 第3类，x3=5
        t(3)=(5-x0(3))/xm(3)

        ! 第4类，小圆交点
        c=x0(1)**2+x0(2)**2-radius(2)**2

        delta=b**2-4*a*c
        if (delta>=0) t(4)=(-b+sqrt(delta))/(2*a)

        do i=1,4
            xc_list(:,i)=x0+xm*t(i)
        enddo

        ! 4种法线
        normal_list=0
        normal_list(1:2,1)=-xc_list(1:2,1)
        normal_list(1:2,4)=-xc_list(1:2,4)
        normal_list(3,2)=1
        normal_list(3,3)=-1

        ! 4类边界的范围
        check(1)=abs(xc_list(3,1))>=n_cell_z*ratio_z/2d0

        r=sqrt(xc_list(1,2)**2+xc_list(2,2)**2)
        check(2)=r<radius(1) .and. r>radius(2)

        r=sqrt(xc_list(1,3)**2+xc_list(2,3)**2)
        check(3)=r<radius(1) .and. r>radius(2)

        check(4)=abs(xc_list(3,4))<n_cell_z*ratio_z/2d0

        if(ratio_z==0)then
            normal=normal_list(:,1)
        elseif(ratio_z>0)then
            mint=1
            do i=1,4
                if (t(i)>=0 .and. t(i)<1 .and. check(i) .and. t(i)<mint) then
                    xc=xc_list(:,i)
                    normal=normal_list(:,i)
                    mint=t(i)
                endif
            enddo
        end if


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
                do while(.not. in_pipe(x(:,i)))
                    call do_ba(x(:,i),x0(:,i),v(:,i))
                    c=c+1
                    !                    if (c>2) then
                    !                        write(*,*)c,i
                    !                        write(*,*)x0(:,i)
                    !                        write(*,*)x(:,i)
                    !                        write(*,*)norm2(x0(1:2,i)),norm2(x(1:2,i))
                    !                        debug=1
                    !                    end if
                enddo
            endif
        enddo
        !$omp end parallel do
    end subroutine

end module
