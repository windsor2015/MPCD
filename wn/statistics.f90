module statistics
    use shape_all , only: ratio_z, n_cell_z, n_p,desk_interval_step
    integer :: n_p_const=100
    real(8) x(4,100), z0(3000000)

    integer :: trans_end_t=-1, trans_begin_t=-1
    integer :: unknot_t=-1
    interface
        subroutine homfly(content,n,string) bind(C, name='homfly')
            import
            integer, value :: n
            integer, dimension(n) :: content
            character(1) string
        end subroutine homfly
    endinterface
contains

    subroutine stat_main(cur_step,x_p,n_p0,n_p,string,Rg,c_axis,std_deviation)

        integer n_p,n_p0,cur_step,i,trans_t
        character(80) string
        real(8) x_p(3,n_p0),Rg,c_axis,std_deviation

        n_p_const=n_p0

        x(1:3,1:n_p0)=x_p(1:3,1:n_p0)
        do i=1,n_p0
            x(4,i)=i
        end do

        n_p=n_p0
        call stat_gyration(x_p,n_p,Rg,c_axis,std_deviation)

        !call translocation_t(x_p,n_p,cur_step)
        !if(trans_begin_t/=0.and.trans_end_t/=0)then
        !    trans_t=trans_end_t-trans_begin_t
        !end if
        !call cal_t(x,n_p,t,hom)
        !测试时只看这一步,用31结测试一下
        call cal_t(x,n_p,string)
        call removeparticle(x,n_p,cur_step)

    end subroutine

    subroutine stat_gyration(x_p,n_p,Rg,c_axis,std_deviation)   !result(Rg,c_axis,std_deviation)
        real(8) Rg,xc(3),x_p(3,n_p),c_axis,coord(n_p),std_deviation
        integer i,n_p
        xc=0d0
        Rg=0d0

        xc(1)=sum(x_p(1,:))
        xc(2)=sum(x_p(2,:))
        xc(3)=sum(x_p(3,:))
        xc=xc/n_p

        Rg=sum((x_p(1,:)-xc(1))**2)+sum((x_p(2,:)-xc(2))**2)+sum((x_p(3,:)-xc(3))**2)
        Rg=sqrt(Rg/n_p)

        coord(:)=sqrt(x_p(1,:)**2+x_p(2,:)**2)
        c_axis=sum(coord(:))
        c_axis=c_axis/n_p
        std_deviation=sum((coord(:)-c_axis)**2)
        std_deviation=sqrt(std_deviation/n_p)
    end subroutine   !function

    subroutine translocation_t(x_p,n_p,cur_step)
        real(8) x_p(3,n_p)
        integer cur_step,n_p

        if(count(x_p(3,:)<-n_cell_z*ratio_z/2d0)==n_p-1 .and. trans_begin_t==-1)then
            trans_begin_t=cur_step
        endif
        if(count(x_p(3,:)>n_cell_z*ratio_z/2d0)==n_p-1 .and. trans_end_t==-1)then
            trans_end_t=cur_step
        end if
    end subroutine

    integer function return_translocation_t()
        return_translocation_t=trans_end_t-trans_begin_t
    end function

    !!!!!!!!!穿孔散结类型!!!!!!!!!!
    integer function trans_unknot_type() !result(trans_unknot_type)
    integer unknot_before
            unknot_before=unknot_t-desk_interval_step
        if (trans_begin_t>0 .and. trans_end_t>0) then
            !            if(unknot_t<trans_begin_t)then
            !                trans_unknot_type=1
            !            elseif(unknot_t<trans_end_t .and. unknot_t>trans_begin_t)then
            !                trans_unknot_type=2
            !            elseif(unknot_t>trans_end_t)then
            !                trans_unknot_type=3
            !            elseif(unknot_t==-1)then
            !                trans_unknot_type=4
            !            end if
            if(z0(unknot_before)<-n_cell_z*ratio_z/2d0)then
                trans_unknot_type=1
            elseif(z0(unknot_before)<=n_cell_z*ratio_z/2d0 .and. z0(unknot_before)>=-n_cell_z*ratio_z/2d0)then
                trans_unknot_type=2
            elseif(z0(unknot_before)>n_cell_z*ratio_z/2d0)then
                trans_unknot_type=3
            elseif(unknot_t==-1)then
                trans_unknot_type=4
            end if
        elseif(trans_begin_t>0 .and. trans_end_t==-1)then
            trans_unknot_type=5
        elseif(trans_begin_t==-1 .and. trans_end_t==-1)then
            trans_unknot_type=6
        end if
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!传字串给homfly,貌似没用到，传了整数数组代替
    subroutine append_string(str0,str)
        character(30) str
        character(1000) str0
        integer l1,l2
        l1=len(trim(str0))
        l2=len(trim(str))
        !write(*,*)l1,l2
        str0(l1+1:l1+l2)=str(1:l2)
    end subroutine

    function inside_triangle(x1,x2,x3,t1,t2)
        implicit none
        real(8), dimension(3) :: x1,x2,x3,t1,t2,normal,tm,xc
        real(8) t,s1,s2,s3,s
        logical inside_triangle

        normal=cross_product(x1-x2,x3-x2)
        tm=t2-t1
        t=dot_product(x2-t1,normal)/dot_product(normal,tm)
        xc=t1+tm*t

        if (t<=0 .or. t>=1) then
            inside_triangle=.false.
            return
        end if

        s=area(x1,x2,x3)
        s1=area(xc,x2,x3)
        s2=area(xc,x1,x3)
        s3=area(xc,x2,x1)

        inside_triangle = abs(s1+s2+s3-s)<1e-8
    end function

    function area(x1,x2,x3)
        implicit none
        real(8), dimension(3) :: x1,x2,x3
        real(8) :: area
        area = norm2(cross_product(x1-x2,x3-x2))/2
    end function

    function cross_product(v1, v2)
        implicit none
        real(8), dimension(3) :: v1, v2
        real(8), dimension(3) :: cross_product

        cross_product(1) = v1(2)*v2(3) - v1(3)*v2(2)
        cross_product(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cross_product(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end function

    subroutine connect(x,n)
        integer f,i,j,n
        real(8) rz,x(4,n)
        real(8),parameter::r0=7
        !n_cell_z=40
        do i=2,n
            rz=x(3,i)-x(3,i-1)
            if(abs(rz)>r0)then
                if(rz>0)then
                    f=-1
                else
                    f=1
                endif
                do j=i,n
                    x(3,j)=x(3,j)+n_cell_z*f
                enddo
            endif
        enddo
        return
    end subroutine


    subroutine output(output_file,cur_step)
        implicit none
        integer cur_step,output_file,k

        write(output_file,*)'ITEM:TIMESTEP'
        write(output_file,'(I9)')cur_step
        write(output_file,*)'ITEM:NUMBER OF ATOMS'
        write(output_file,'(I6)')n_p
        write(output_file,*)'ITEM:BOX BOUNDS'
        write(output_file,'(2F7.1)')-6d0,6d0
        write(output_file,'(2F7.1)')-6d0,6d0
        write(output_file,'(2F7.1)')-20d0,20d0
        write(output_file,*)'ITEM:ATOMS id type x y z'

        do k=1,n_p
            write(output_file,'(I6,I3,3F9.4)') k,1,x(1:3,k)
        enddo

    end subroutine

    !input prmeter never chned
    subroutine removeparticle1(x0,n_p0,cur_step,mode,n_reult,z)
        integer n_p0, n_p, i, j, t, i1, i2, i3, minj, maxj,output_file,q,mode
        integer i_begin,i_end,i_step,select_i,cur_step,n_reult
        real(8) x(4,n_p0), x0(4,n_p0)
        logical removed, t_f
        real(8) maxcos, curcos
        real(8) r1(3),r2(3),z

        removed=.true.

        x=x0
        n_p=n_p0
        !call output(output_file,q)

        do while (removed)
            removed=.false.
            ! 检查粒子
            !call output(output_file,q)
            select_i=-1
            select case(mode)
                case(1,2)
                    if(mode==1)then
                        i_begin=2
                        i_end=n_p-1
                        i_step=1
                    else
                        i_begin=n_p-1
                        i_end=2
                        i_step=-1
                    endif
                    !write(*,*)i_ein,i_end,i_tep
                    do i=i_begin,i_end,i_step
                        t=0
                        ! 线段
                        do j=1,n_p-1
                            if ( j+1<i-1 .or. j>i+1 ) then
                                t_f = inside_triangle(x(1:3,i),x(1:3,i+1),x(1:3,i-1),x(1:3,j),x(1:3,j+1))
                                if (t_f) then
                                    t=t+1
                                    !write(*,*) i,j,n_p,t_f
                                    exit
                                endif
                            endif
                        enddo
                        if (t==0)then
                            select_i=i
                            exit
                        end if
                        !如果t为0, 可以被移除, 如果t不是0则不能移除
                    enddo
                case(3)
                    removed=.false.
                    maxcos=-10
                    select_i=-1
                    do i=2,n_p-1
                        r1=x(1:3,i)-x(1:3,i-1)
                        r2=x(1:3,i)-x(1:3,i+1)
                        curcos=dot_product(r1,r2)/norm2(r1)/norm2(r2)
                        !write(*,*) curcos
                        t_f=.false.
                        do j=1,n_p-1
                            if ( j+1<i-1 .or. j>i+1 ) then
                                t_f = inside_triangle(x(1:3,i),x(1:3,i+1),x(1:3,i-1),x(1:3,j),x(1:3,j+1))
                            endif
                            if (t_f) exit
                        end do
                        if ((.not. t_f) .and. curcos>maxcos) then !
                            select_i=i
                            maxcos=curcos
                        end if
                    end do
                case(4)
                    if (n_p==n_p0) then
                        select_i=2
                    else
                        select_i=select_i+1
                    endif
                    if (.not.removed)select_i=select_i+1
                    if (select_i>=n_p) then
                        select_i=2
                    endif
                    t_f=.false.
                    i=select_i
                    do j=1,n_p-1
                        if ( j+1<i-1 .or. j>i+1 ) then
                            t_f = inside_triangle(x(1:3,i),x(1:3,i+1),x(1:3,i-1),x(1:3,j),x(1:3,j+1))
                        endif
                        if (t_f) then
                            select_i=-1
                            exit
                        end if
                    end do
            end select

            if(select_i>0)then
                !write(*,*)select_i,n_p
                x(1:4,select_i:n_p-1) = x(1:4,select_i+1:n_p)
                n_p=n_p-1
                removed=.true.
            endif
            ! write(*,*)n_p
            if (n_p<=2) exit
        enddo

        call meat_on_bone(x,n_p,x0,n_p0,n_reult,z)
        !write(*,*) n_p,n_p0,n_reult
        !call output(output_file,q)

    end subroutine


    subroutine meat_on_bone(x_b,n_b,x_m,n_m,measure,coord)
        integer n_b,n_m,measure
        real(8) x_b(4,n_b),x_m(4,n_m),coord
        integer i1,i2,i3,minj,maxj,j

        if (n_b>2) then
            minj=n_m
            maxj=0
            ! write(*,*)n_p,n_p0
            do i1=2,n_b-3; do i2=i1+1,n_b-2; do i3=i2+1,n_b-1
                do j=1,n_m-1
                    if ( j+1<x_b(4,i1) .or. j>x_b(4,i3) ) then
                        if (inside_triangle(x_b(:,i1),x_b(:,i2),x_b(:,i3),x_m(:,j),x_m(:,j+1))) then
                            minj=min(j,minj,int(x_b(4,i1))-1)
                            maxj=max(j+1,maxj,int(x_b(4,i3))+1)
                            !write(*,*) i1,i2,i3,j
                        end if
                    endif
                end do
            enddo; enddo; enddo
            !write(*,*) minj, maxj
            !x(:,1:maxj-minj+1)=x0(:,minj:maxj)
            measure=maxj-minj+1
        else
            measure=0
        end if
        if(measure==0)then
            coord=0
        else
            coord=sum(x_m(3,minj:maxj))/measure
        endif
        !write(*,*)poition,ize
    endsubroutine

    subroutine removeparticle(x,n_p,q)
        integer n_p, n_p0,i,q,n_reult
        real(8) x(4,n_p),x0(4,n_p),zk(3),z
        integer n_pk(3)
        integer:: method=0

        call removeparticle1(x,n_p,q,1,n_pk(1),zk(1))
        call removeparticle1(x,n_p,q,2,n_pk(2),zk(2))
        call removeparticle1(x,n_p,q,3,n_pk(3),zk(3))
        !write(*,*)n_pk,zk
        n_p=minval(n_pk, n_pk>=0)
        z0(q)=zk(minloc(n_pk,1,n_pk>=0))

    end subroutine

    subroutine cal_t(x,n_p,string)
        real(8) x(4,n_p),t(6,n_p*2)
        integer n_p
        real(8) xp(3,n_p),x0(2),x1(2),y0(2),y1(2),v1(2),v2(2), s
        real(8) a11,a12,a21,a22,b1,b2,t1,t2,deta
        integer i,j,k,c,l
        integer content(1000),n
        character(80) string

        string=' '
        string='0'

        xp(1,:)=x(2,:)
        xp(2,:)=x(3,:)
        xp(3,:)=x(1,:)
        t=-1
        k=0
        do i=1,n_p-2
            do j=i+2,n_p-1
                x0=xp(1:2,i)
                x1=xp(1:2,i+1)
                y0=xp(1:2,j)
                y1=xp(1:2,j+1)
                b1=-x0(1)+y0(1)
                b2=-x0(2)+y0(2)
                a11=x1(1)-x0(1)
                a21=x1(2)-x0(2)
                a12=-y1(1)+y0(1)
                a22=-y1(2)+y0(2)
                deta=a11*a22-a12*a21
                t1=(b1*a22-a12*b2)/deta
                t2=(b2*a11-a21*b1)/deta
                if (t1>=0 .and. t1<1 .and. t2>=0 .and. t2<1) then
                    t(1,k+1)=i+t1
                    t(3,k+1)=xp(3,i)+t1*(xp(3,i+1)-xp(3,i))
                    t(2,k+1)=j+t2
                    t(1,k+2)=j+t2
                    t(3,k+2)=xp(3,j)+t2*(xp(3,j+1)-xp(3,j))
                    t(2,k+2)=i+t1
                    v1=x1-x0
                    v2=y1-y0
                    s=sign(1d0,v1(1)*v2(2)-v1(2)*v2(1))
                    if (t(3,k+1)<t(3,k+2)) then
                        s=-s
                    end if
                    t(6,k+1)=s
                    t(6,k+2)=s
                    !write(*,*) i,j,xp(:,i),xp(:,j)
                    k=k+2
                end if
            end do
        end do

        call quick_sort(t,k,1,k)

        do i=1,k
            !t(6,i)=i
        end do

        do i=1,k-1
            do j=i+1,k
                if (t(1,i)==t(2,j)) then
                    t(4,i)=j
                    t(4,j)=i
                    t(5,i)=sign(1d0,t(3,i)-t(3,j))
                    t(5,j)=-t(5,i)
                    exit
                end if
            end do
        end do

        !write(*,'(6F8.3)') t(:,1:k)

        !write(*,*) 1
        content(1)=1
        !write(*,*) k
        content(2)=k

        if (k==0) return

        c=1
        n=3
        do i=1,k
            l=-1
            do j=1,i-1
                if (int(t(4,j))==i) then
                    l=t(3,j)
                end if
            end do
            !找到了新的点
            if (l==-1) then
                !write(*,'(2I3,$)') c-1, int(t(5,i))
                content(n)=c-1
                content(n+1)=int(t(5,i))
                t(3,i)=c
                c=c+1
            else
                !write(*,'(2I3,$)') l-1, int(t(5,i))
                content(n)=l-1
                content(n+1)=int(t(5,i))
            end if
            n=n+2
        end do
        !write(*,*)
        !左右旋
        c=1
        do i=1,k
            l=-1
            do j=1,i-1
                if (int(t(4,j))==i) then
                    l=i
                end if
            end do
            !找到了新的点
            if (l==-1) then
                !write(*,*) c-1, int(t(6,i))
                content(n)=c-1
                content(n+1)=int(t(6,i))
                n=n+2
                c=c+1
            end if
        end do
        n=n-1
        !write(*,*) content(1:n)

        call homfly(content, n, string)
        do i=1,80
            if (iachar(string(i:i))<32) then
                string(i:i)=' '
            end if
        end do
        !write(*,*) string
        !stop
    end subroutine

    ! 交换
    subroutine swap_r8(a,b)
        real(8), dimension(6) :: a,b,c
        c=a
        a=b
        b=c
    end subroutine

    ! 快速排序
    recursive subroutine quick_sort(v,n,s,e)
        real(8) v(6,n),key
        integer n,s,e,l,r,m

        l=s
        r=e
        m=(s+e)/2
        if (l>=r) return
        key=v(1,m)
        do while (l<r)
            do while(v(1,l)<key)
                l=l+1
            enddo
            do while(v(1,r)>key)
                r=r-1
            enddo
            if (l<r) then
                call swap_r8(v(:,l),v(:,r))
                ! 与key相等则多走一位, 注意是交换后的值
                if (v(1,r)==key) l=l+1
                if (v(1,l)==key) r=r-1
            endif
        enddo
        call quick_sort(v,n,s,l-1)
        call quick_sort(v,n,r+1,e)
    end subroutine

    function is_knot(x_,n_)
        implicit none
        integer n_, n
        real(8) x_(3,n_), x(3,n_)
        logical removed, is_knot, t_f
        integer i,j,t

        removed=.true.
        is_knot=.true.
        x=x_
        n=n_
        call connect(x,n)
        do while (removed)
            removed=.false.
            ! 检查粒子
            !call output(output_file,q)
            do i=1,n-2,1
                t=0
                ! 线段
                do j=1,n-1
                    if ( j+1<i .or. j>i+2 ) then
                        t_f = inside_triangle(x(1:3,i),x(1:3,i+1),x(1:3,i+2),x(1:3,j),x(1:3,j+1))
                        if (t_f) then
                            t=t+1
                            !write(*,*) i,j,t_f
                        endif
                    endif
                enddo
                if (t==0)exit
                !如果t为0, 可以被移除, 如果t不是0则不能移除
            enddo
            if(t==0)then
                x(:,i+1:n-1) = x(:,i+2:n)
                n=n-1
                removed=.true.
            endif
            ! write(*,*)n_p
            if (n<=2) then
                is_knot=.false.
                exit
            endif
        enddo
    end function

end module

