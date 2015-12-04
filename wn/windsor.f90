PROGRAM translocation
    PARAMETER(n=40,ncell=15,ncelly=15,nn=10,mn=nn*ncell*ncelly*ncell) !!!n是单体数目，n0是单体、壁、孔数目之和
    REAL(8) x(n),y(n),z(n),vx(n),vy(n),vz(n),x11(mn),y11(mn),z11(mn),dx1(mn),dy1(mn),dz1(mn),x12(mn),y12(mn),z12(mn),y13(mn),dx2(mn),dy2(mn),dz2(mn),x0(mn),y0(mn),z0(mn),x10(mn),y10(mn),z10(mn)
    REAL(8) x1(mn),y1(mn),z1(mn),vx1(mn),vy1(mn),vz1(mn),vxx(4500),vyy(4500),vzz(4500),d1(n),dd0,dd1,dd,vx10(mn),vy10(mn),vz10(mn),g(mn),gee,vu(mn),vvu(mn),Ek,U,ppx,ppy,ppz,x30(n),y30(n)
    REAL(8) kx1(n),ky1(n),kz1(n),kx2(n),ky2(n),kz2(n),kx3(n),ky3(n),kz3(n),f,rannx(n),ranny(n),rannz(n),ran1,ran2,rco,lx,ly,lz
    REAL(8) vp,scalx,scaly,scl,scalz,vp2,ranx,rany,ranz,rx,ry,rz,rr,rr1,rl,mx,my,mz,mx1,my1,mz1,kin1,kin2,kin11,kin22,Ek0,Ek1,Ek2
    REAL(8) rx1,ry1,rz1,rx2,ry2,rz2,rx3,ry3,rz3,r1,r2,r3,alpha,c,s1,uu,dlty,R22,xcm,ycm,zcm,S2mm(n),gg(n),ge,x_up1(256),y_up1(256),z_up1(256),x_down1(256),y_down1(256),z_down1(256)
    REAL(8) rx00,ry00,rz00,rxyz,ddx1,ddy1,ddz1,ddx2,ddy2,ddz2,ddx3,ddy3,ddz3,tol,tol1,dt,sf,fac,dx0(mn),dy0(mn),dz0(mn),theta,fai
    REAL(8) rrx1,rry1,rrx2,rry2,rrz2,rr2,rrx3,rry3,rrz3,rr3,rr4,t0,t1,t2,cc,ss,sss,c1,c2,c3,c12,c22,c32,pc,pcc,pyx,yta,sf0(90000000)
    REAL(8) s(3,3),d(3),vv(n),vvv(n),kB,temp,temp1,tmp,m,ms,r,ll,d0,r00,ld1,ld2,small,x_up(256),y_up(256),z_up(256),x_down(256),y_down(256),z_down(256)
    INTEGER nc,ncc,i,j,ii,jj,p,n0,l,t,q1,q0,q,qq,step0,step,st01,st0,st00,st1,st11,st2,st22,st3,a,n2,int,idum,gama,num(ncelly)            !!! kB为玻尔兹曼常数

    REAL(8) s0,e0,k,kk,pi,sum                               !!! s0是LJ势的深度, k是弹簧劲度系数,柔性连对应kk为0，半柔性连对应kk为10
    character(len=20) filename

    step0=3800000
    st1=3800000
    st11=10000
    step=500000
    st0=199000
    st01=1000
    st00=200000
    st2=1000
    st22=10000
    st3=200

    idum=-4500

    pi=dble(3.14159265)
    kB=1.380662D-23
    m=dble(1)
    ms=dble(0.1)
    temp=dble(1)
    temp1=dble(1)
    fac=dble(1.5)

    d0=dble(1)
    r00=dble(1.5)
    s0=dble(1)
    k=dble(30)
    kk=dble(0)
    ld2=d0*2.0**(1.0/6.0)         !!! LJ势排斥力的作用范围2.0**(1.0/6.0)
    tol=dble(0.0001)
    tol1=dble(0.0001)
    boxx=dble(ncell)
    boxy=dble(ncelly)
    boxz=dble(ncell)
    sid=boxx/dble(ncell)
    hboxx=boxx/dble(2.0)
    hboxz=boxz/dble(2.0)
    hboxy=boxy/dble(2.0)
    hsid=sid/dble(2.0)

    uu=dble(0.2)

    alpha=dble(130*pi/180)

    c=cos(alpha)
    s1=sin(alpha)
    small=dble(0.001)
    !!!!!! assign initial coordinate


    OPEN(11,FILE='1.txt')
    do i=1,n
        READ(11,*) x(i)
    end do
    CLOSE(11)

    OPEN(12,FILE='2.txt')
    do i=1,n
        READ(12,*) y(i)
    end do
    CLOSE(12)

    OPEN(13,FILE='3.txt')
    do i=1,n
        READ(13,*) z(i)
    end do
    CLOSE(13)

    !!!!!!!!!!!!!!!!!建立初始平板
    do i=1,ncell+1
        do j=1,ncell+1
            ii=i+(j-1)*(ncell+1)
            x_up(ii)=(i-1)*sid
            z_up(ii)=(j-1)*sid
            x_down(ii)=(i-1)*sid
            z_down(ii)=(j-1)*sid
        enddo
    enddo
    y_up(:)=dble(ncelly)

    y_down(:)=0

    !OPEN(2,FILE='rongji.xyz')
    i=0
    do l=0,ncell-1
        do j=0,ncelly-1
            do ii=0,ncell-1
                do lg=1,nn
                    i=i+1
                    x1(i)=ii*sid+lg*sid/dble(nn+1)
                    y1(i)=j*sid+lg*sid/dble(nn+1)
                    z1(i)=l*sid+lg*sid/dble(nn+1)
                enddo
            enddo
        enddo
    enddo

    DO i=1,mn
        x1(i)=x1(i)+(rand(idum)-0.5)*sid
        y1(i)=y1(i)+(rand(idum)-0.5)*sid
        z1(i)=z1(i)+(rand(idum)-0.5)*sid
    ENDDO

    DO i=1,mn
        x1(i)=x1(i)-boxx*nint((x1(i)-hboxx)/boxx)
        y1(i)=y1(i)-boxy*nint((y1(i)-hboxy)/boxy)
        z1(i)=z1(i)-boxz*nint((z1(i)-hboxz)/boxz)
    ENDDO

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

    OPEN(70,FILE='qU0.out')
    OPEN(80,FILE='qU.out')
    OPEN(81,FILE='pcc.out')
    OPEN(821,FILE='gxy.out')
    OPEN(82,FILE='gzr.out')
    OPEN(83,FILE='masscenter.out')
    OPEN(84,FILE='yta.out')
    OPEN(88,FILE='numy.out')
    OPEN(90,FILE='R22.out')
    OPEN(91,FILE='sf.out')
    OPEN(100,FILE='vxyz1.out')
    OPEN(110,FILE='pt.out')
    OPEN(111,FILE='cc.out')
    OPEN(200,FILE='vxyz0.out')
    OPEN(210,FILE='xyz.out')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!assign initial velocity
    !!!链的初速度
    mx=0.0
    my=0.0
    mz=0.0
    call random_number(vx)
    call random_number(vy)
    call random_number(vz)
    do i=1,n
        vx(i)=vx(i)-0.5
        vy(i)=vy(i)-0.5
        vz(i)=vz(i)-0.5
        mx=mx+vx(i)
        my=my+vy(i)
        mz=mz+vz(i)
    end do
    mx=mx/dble(n)
    my=my/dble(n)
    mz=mz/dble(n)


    kin1=0.0
    do i=1,n
        vx(i)=vx(i)-mx
        vy(i)=vy(i)-my
        vz(i)=vz(i)-mz
        kin1=kin1+dble(0.5)*m*(vx(i)**2+vy(i)**2+vz(i)**2)
    end do

    tmp=kin1/(fac*dble(n))
    scl=sqrt(temp/tmp)

    kin2=0.0
    do i=1,n
        vx(i)=scl*vx(i)
        vy(i)=scl*vy(i)
        vz(i)=scl*vz(i)
        kin2=kin2+dble(0.5)*m*(vx(i)**2+vy(i)**2+vz(i)**2)
    end do

    tmp=kin2/(fac*dble(n))
    write(*,*) '初始标度后动能kin2=',kin2
    write(*,*) '初始标度后温度tmp=',tmp

    !!!溶剂的初速度
    mx=0.0
    my=0.0
    mz=0.0
    call random_number(vx1)
    call random_number(vy1)
    call random_number(vz1)
    do i=1,mn
        vx1(i)=vx1(i)-0.5
        vy1(i)=vy1(i)-0.5
        vz1(i)=vz1(i)-0.5
        mx=mx+vx1(i)
        my=my+vy1(i)
        mz=mz+vz1(i)
    end do
    mx=mx/dble(mn)
    my=my/dble(mn)
    mz=mz/dble(mn)


    kin1=0.0
    do i=1,mn
        vx1(i)=vx1(i)-mx
        vy1(i)=vy1(i)-my
        vz1(i)=vz1(i)-mz
        kin1=kin1+dble(0.5)*ms*(vx1(i)**2+vy1(i)**2+vz1(i)**2)
    end do

    tmp=kin1/(fac*dble(mn)*ms)
    scl=sqrt(temp1/tmp)

    kin2=0.0
    do i=1,mn
        vx1(i)=scl*vx1(i)
        vy1(i)=scl*vy1(i)
        vz1(i)=scl*vz1(i)
        kin2=kin2+dble(0.5)*ms*(vx1(i)**2+vy1(i)**2+vz1(i)**2)
    end do

    tmp=kin2/(fac*dble(mn)*ms)
    write(*,*) '初始标度后动能kin2=',kin2
    write(*,*) '初始标度后温度tmp=',tmp

    DO l=0,ncell-1
        DO j=0,ncelly-1
            DO ii=0,ncell-1

                ki=ii+1+j*ncell+l*ncelly*ncell

                x10(ki)=ii*sid+sid/dble(2)
                y10(ki)=j*sid+sid/dble(2)
                z10(ki)=l*sid+sid/dble(2)

                vx0=0.0
                vy0=0.0
                vz0=0.0
                ncc=0
                DO i=1,n
                    IF(x(i)>=ii*sid.AND.x(i)<(ii+1)*sid.AND.y(i)>=j*sid.AND.y(i)<(j+1)*sid.AND.z(i)>=l*sid.AND.z(i)<(l+1)*sid)THEN
                        vx0=vx0+vx(i)*m
                        vy0=vy0+vy(i)*m
                        vz0=vz0+vz(i)*m
                        ncc=ncc+1
                        gg(ncc)=i
                    ENDIF
                ENDDO
                nc=0
                DO i=1,mn
                    IF(x1(i)>=ii*sid.AND.x1(i)<(ii+1)*sid.AND.y1(i)>=j*sid.AND.y1(i)<(j+1)*sid.AND.z1(i)>=l*sid.AND.z1(i)<(l+1)*sid)THEN
                        vx0=vx0+vx1(i)*ms
                        vy0=vy0+vy1(i)*ms
                        vz0=vz0+vz1(i)*ms
                        nc=nc+1
                        g(nc)=i
                    ENDIF
                ENDDO

                ki=ii+1+j*ncell+l*ncell*ncell

                if((nc==0).AND.(ncc==0))then
                    vxx(ki)=0.0
                    vyy(ki)=0.0
                    vzz(ki)=0.0
                else
                    vxx(ki)=vx0/dble(nc*ms+ncc*m)
                    vyy(ki)=vy0/dble(nc*ms+ncc*m)
                    vzz(ki)=vz0/dble(nc*ms+ncc*m)
                endif

                write(200,*)vxx(ki),vyy(ki),vzz(ki)
                write(210,*)x10(ki),y10(ki),z10(ki)
            ENDDO
        ENDDO
    ENDDO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! compute a(t-dt)

    DO i=1,n
        kx2(i)=dble(0)
        ky2(i)=dble(0)
        kz2(i)=dble(0)
    ENDDO

    !!!!!! UFENE(r) force

    DO i=2,n-1
        DO j=i+1,i-1,-2
            rx=x(i)-x(j)
            ry=y(i)-y(j)
            rz=z(i)-z(j)

            r=SQRT(rx**2+ry**2+rz**2)

            IF(r/=0.AND.r<=r00) THEN
                rr1=(k*r00**2)/(r00**2-r**2)
                kx2(i)=kx2(i)-rx*rr1
                ky2(i)=ky2(i)-ry*rr1
                kz2(i)=kz2(i)-rz*rr1
            ENDIF
        ENDDO
    ENDDO

    rx=x(1)-x(2)
    ry=y(1)-y(2)
    rz=z(1)-z(2)

    r=SQRT(rx**2+ry**2+rz**2)

    IF(r/=0.AND.r<=r00) THEN
        rr1=(k*r00**2)/(r00**2-r**2)
        kx2(1)=kx2(1)-rx*rr1
        ky2(1)=ky2(1)-ry*rr1
        kz2(1)=kz2(1)-rz*rr1
    ENDIF

    rx=x(n)-x(n-1)
    ry=y(n)-y(n-1)
    rz=z(n)-z(n-1)

    r=SQRT(rx**2+ry**2+rz**2)

    IF(r/=0.AND.r<=r00) THEN
        rr1=(k*r00**2)/(r00**2-r**2)
        kx2(n)=kx2(n)-rx*rr1
        ky2(n)=ky2(n)-ry*rr1
        kz2(n)=kz2(n)-rz*rr1
    ENDIF
    !!!!!! end UFENE(r) force


    !!!!!! ULJ(r) force

    DO i=1,n
        DO j=1,n
            IF(i/=j)THEN
                rx=x(i)-x(j)
                ry=y(i)-y(j)
                rz=z(i)-z(j)
                r=SQRT(rx**2+ry**2+rz**2)

                IF(r<=ld2)THEN
                    rr=(d0/r)**6
                    rr1=dble(48)*s0*(rr**2)/(r**2)
                    kx2(i)=kx2(i)+rx*rr1
                    ky2(i)=ky2(i)+ry*rr1
                    kz2(i)=kz2(i)+rz*rr1
                ENDIF
            ENDIF
        ENDDO
    ENDDO

    !!!!!! end ULJ(r) force


    !!!!!!!bend energy!!!!!!!!!

    DO i=2,n-1
        rx1=x(i+1)-x(i)
        ry1=y(i+1)-y(i)
        rz1=z(i+1)-z(i)

        rx2=x(i)-x(i-1)
        ry2=y(i)-y(i-1)
        rz2=z(i)-z(i-1)

        kx2(i)=kx2(i)+kk*(rx1-rx2)
        ky2(i)=ky2(i)+kk*(ry1-ry2)
        kz2(i)=kz2(i)+kk*(rz1-rz2)
    ENDDO

    !!!!!!!!!!!!!!!!!!!end bend energy

    DO i=1,n
        ry=y(i)
        if(ry<=ld2)then
            rr=(d0/ry)**6
            rr1=dble(48)*s0*(rr**2)/ry
            ky2(i)=ky2(i)+rr1
        endif
        ry=dble(ncelly)-y(i)
        if(ry<=ld2)then
            rr=(d0/ry)**6
            rr1=dble(48)*s0*(rr**2)/ry
            ky2(i)=ky2(i)-rr1
        endif
    ENDDO


    DO q0=1,step0

        DO i=1,n
            x(i)=x(i)+vx(i)*tol+dble(0.5)*kx2(i)*tol**2
            y(i)=y(i)+vy(i)*tol+dble(0.5)*ky2(i)*tol**2
            z(i)=z(i)+vz(i)*tol+dble(0.5)*kz2(i)*tol**2
        ENDDO

        !!!!!! compute a(t+dt)

        DO i=1,n
            kx1(i)=dble(0)
            ky1(i)=dble(0)
            kz1(i)=dble(0)
        ENDDO

        !!!!!! UFENE(r) force


        U1=dble(0)

        DO i=2,n-1
            DO j=i+1,i-1,-2
                rx=x(i)-x(j)
                ry=y(i)-y(j)
                rz=z(i)-z(j)
                r=SQRT(rx**2+ry**2+rz**2)
                IF(r/=0.AND.r<=r00) THEN
                    rr1=(k*r00**2)/(r00**2-r**2)
                    kx1(i)=kx1(i)-rx*rr1
                    ky1(i)=ky1(i)-ry*rr1
                    kz1(i)=kz1(i)-rz*rr1
                    U1=U1-(k*r00**2)*log(1-(r/r00)**2)/2
                ENDIF
            ENDDO
        ENDDO

        rx=x(1)-x(2)
        ry=y(1)-y(2)
        rz=z(1)-z(2)
        r=SQRT(rx**2+ry**2+rz**2)
        IF(r/=0.AND.r<=r00) THEN
            rr1=(k*r00**2)/(r00**2-r**2)
            kx1(1)=kx1(1)-rx*rr1
            ky1(1)=ky1(1)-ry*rr1
            kz1(1)=kz1(1)-rz*rr1
            U1=U1-(k*r00**2)*log(1-(r/r00)**2)/2
        ENDIF

        rx=x(n)-x(n-1)
        ry=y(n)-y(n-1)
        rz=z(n)-z(n-1)
        r=SQRT(rx**2+ry**2+rz**2)
        IF(r/=0.AND.r<=r00) THEN
            rr1=(k*r00**2)/(r00**2-r**2)
            kx1(n)=kx1(n)-rx*rr1
            ky1(n)=ky1(n)-ry*rr1
            kz1(n)=kz1(n)-rz*rr1
            U1=U1-(k*r00**2)*log(1-(r/r00)**2)/2
        ENDIF
        !!!!!! end UFENE(r) force


        !!!!!! ULJ(r) force

        DO i=1,n
            DO j=1,n
                IF(i/=j)THEN
                    rx=x(i)-x(j)
                    ry=y(i)-y(j)
                    rz=z(i)-z(j)
                    r=SQRT(rx**2+ry**2+rz**2)
                    IF(r<=ld2)THEN
                        rr=(d0/r)**6
                        rr1=dble(48)*s0*(rr**2)/(r**2)
                        kx1(i)=kx1(i)+rx*rr1
                        ky1(i)=ky1(i)+ry*rr1
                        kz1(i)=kz1(i)+rz*rr1
                        U1=U1+4*s0*((rr)**2-rr)
                    ENDIF
                ENDIF
            ENDDO
        ENDDO

        !!!!!! end ULJ(r) force


        !!!!!!!bend energy!!!!!!!!!

        DO i=2,n-1
            rx1=x(i+1)-x(i)
            ry1=y(i+1)-y(i)
            rz1=z(i+1)-z(i)

            rx2=x(i)-x(i-1)
            ry2=y(i)-y(i-1)
            rz2=z(i)-z(i-1)

            kx1(i)=kx1(i)+kk*(rx1-rx2)
            ky1(i)=ky1(i)+kk*(ry1-ry2)
            kz1(i)=kz1(i)+kk*(rz1-rz2)
            U1=U1-kk*(rx1*rx2+ry1*ry2+rz1*rz2)
        ENDDO

        !!!!!!!!!!!!!!!!!!!end bend energy

        DO i=1,n
            ry=y(i)
            if(ry<=ld2)then
                rr=(d0/ry)**6
                rr1=dble(48)*s0*(rr**2)/ry
                ky1(i)=ky1(i)+rr1
            endif
            ry=dble(ncelly)-y(i)
            if(ry<=ld2)then
                rr=(d0/ry)**6
                rr1=dble(48)*s0*(rr**2)/ry
                ky1(i)=ky1(i)-rr1
            endif
        ENDDO


        DO i=1,n
            vx(i)=vx(i)+dble(0.5)*(kx1(i)+kx2(i))*tol
            vy(i)=vy(i)+dble(0.5)*(ky1(i)+ky2(i))*tol
            vz(i)=vz(i)+dble(0.5)*(kz1(i)+kz2(i))*tol
        ENDDO

        Ek1=0.0
        DO i=1,n
            Ek1=Ek1+dble(0.5)*m*(vx(i)**2+vy(i)**2+vz(i)**2)
        ENDDO

        t1=Ek1/(fac*dble(n))
        scl=SQRT(temp/t1)

        DO i=1,n
            vx(i)=scl*vx(i)
            vy(i)=scl*vy(i)
            vz(i)=scl*vz(i)
        ENDDO

        DO i=1,n
            kx2(i)=kx1(i)
            ky2(i)=ky1(i)
            kz2(i)=kz1(i)
        ENDDO

        IF(MOD(q0,st1)==0)THEN

            call tran(q0+5,filename)
            open(31,file=filename)
            write(31,*) mn+512+n
            write(31,*)'cohar'

            do i=1,n

                write(31,'(a4,f17.7,f13.7,f13.7)')'p',x(i),y(i),z(i)

            enddo

            do i=1,512

                write(31,'(a4,f17.7,f13.7,f13.7)')'o',x_up(i),y_up(i),z_up(i)

            enddo

            do i=1,512

                write(31,'(a4,f17.7,f13.7,f13.7)')'o',x_down(i),y_down(i),z_down(i)

            enddo
            do i=1,mn
                write(31,'(a4,f17.7,f13.7,f13.7)') 'f',x1(i),y1(i),z1(i)
            end do
            close(31)

        ENDIF

        IF(MOD(q0,st11)==0)THEN

            write(70,*)q0,U1

        ENDIF

    ENDDO        !!!!!step0


    !!! compute a(t-dt)

    DO i=1,n
        kx2(i)=dble(0)
        ky2(i)=dble(0)
        kz2(i)=dble(0)
    ENDDO

    !!!!!! UFENE(r) force

    DO i=2,n-1
        DO j=i+1,i-1,-2
            rx=x(i)-x(j)
            ry=y(i)-y(j)
            rz=z(i)-z(j)
            r=SQRT(rx**2+ry**2+rz**2)
            IF(r/=0.AND.r<=r00) THEN
                rr1=(k*r00**2)/(r00**2-r**2)
                kx2(i)=kx2(i)-rx*rr1
                ky2(i)=ky2(i)-ry*rr1
                kz2(i)=kz2(i)-rz*rr1
            ENDIF
        ENDDO
    ENDDO

    rx=x(1)-x(2)
    ry=y(1)-y(2)
    rz=z(1)-z(2)
    r=SQRT(rx**2+ry**2+rz**2)
    IF(r/=0.AND.r<=r00) THEN
        rr1=(k*r00**2)/(r00**2-r**2)
        kx2(1)=kx2(1)-rx*rr1
        ky2(1)=ky2(1)-ry*rr1
        kz2(1)=kz2(1)-rz*rr1
    ENDIF

    rx=x(n)-x(n-1)
    ry=y(n)-y(n-1)
    rz=z(n)-z(n-1)
    r=SQRT(rx**2+ry**2+rz**2)
    IF(r/=0.AND.r<=r00) THEN
        rr1=(k*r00**2)/(r00**2-r**2)
        kx2(n)=kx2(n)-rx*rr1
        ky2(n)=ky2(n)-ry*rr1
        kz2(n)=kz2(n)-rz*rr1
    ENDIF
    !!!!!! end UFENE(r) force


    !!!!!! ULJ(r) force

    DO i=1,n
        DO j=1,n
            IF(i/=j)THEN
                rx=x(i)-x(j)
                ry=y(i)-y(j)
                rz=z(i)-z(j)
                r=SQRT(rx**2+ry**2+rz**2)
                IF(r<=ld2)THEN
                    rr=(d0/r)**6
                    rr1=dble(48.0)*s0*(rr**2)/(r**2)
                    kx2(i)=kx2(i)+rx*rr1
                    ky2(i)=ky2(i)+ry*rr1
                    kz2(i)=kz2(i)+rz*rr1
                ENDIF
            ENDIF
        ENDDO
    ENDDO

    !!!!!! end ULJ(r) force


    !!!!!!!bend energy!!!!!!!!!

    DO i=2,n-1
        rx1=x(i+1)-x(i)
        ry1=y(i+1)-y(i)
        rz1=z(i+1)-z(i)

        rx2=x(i)-x(i-1)
        ry2=y(i)-y(i-1)
        rz2=z(i)-z(i-1)

        kx2(i)=kx2(i)+kk*(rx1-rx2)
        ky2(i)=ky2(i)+kk*(ry1-ry2)
        kz2(i)=kz2(i)+kk*(rz1-rz2)
    ENDDO

    !!!!!!!!!!!!!!!!!!!end bend energy

    DO i=1,n
        ry=y(i)
        if(ry<=ld2)then
            rr=(d0/ry)**6
            rr1=dble(48)*s0*(rr**2)/ry
            ky2(i)=ky2(i)+rr1
        endif
        ry=dble(ncelly)-y(i)
        if(ry<=ld2)then
            rr=(d0/ry)**6
            rr1=dble(48)*s0*(rr**2)/ry
            ky2(i)=ky2(i)-rr1
        endif
    ENDDO

    write(*,*)'循环开始'

    DO i=1,256
        x_up1(i)=x_up(i)
        z_up1(i)=z_up(i)
        x_down1(i)=x_down(i)
        z_down1(i)=z_down(i)
    ENDDO

    DO i=1,mn
        x11(i)=x1(i)
        y11(i)=y1(i)
        z11(i)=z1(i)
    ENDDO

    q1=1

    DO q=1,step

        DO st=1,st3

            !!!!!! compute a(t+dt)

            DO i=1,n
                x(i)=x(i)+vx(i)*tol+dble(0.5)*kx2(i)*tol**2-0.1*(y(i)**2-15.0*y(i))*uu*tol
                y(i)=y(i)+vy(i)*tol+dble(0.5)*ky2(i)*tol**2
                z(i)=z(i)+vz(i)*tol+dble(0.5)*kz2(i)*tol**2
            ENDDO

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            DO i=1,n
                kx1(i)=dble(0)
                ky1(i)=dble(0)
                kz1(i)=dble(0)
            ENDDO

            !!!!!! UFENE(r) force

            U=dble(0)

            DO i=2,n-1
                DO j=i+1,i-1,-2
                    rx=x(i)-x(j)
                    ry=y(i)-y(j)
                    rz=z(i)-z(j)
                    r=SQRT(rx**2+ry**2+rz**2)
                    IF(r/=0.AND.r<=r00) THEN
                        rr1=(k*r00**2)/(r00**2-r**2)
                        kx1(i)=kx1(i)-rx*rr1
                        ky1(i)=ky1(i)-ry*rr1
                        kz1(i)=kz1(i)-rz*rr1
                        U=U-(k*r00**2)*log(1-(r/r00)**2)/2
                    ENDIF
                ENDDO
            ENDDO

            rx=x(1)-x(2)
            ry=y(1)-y(2)
            rz=z(1)-z(2)
            r=SQRT(rx**2+ry**2+rz**2)
            IF(r/=0.AND.r<=r00) THEN
                rr1=(k*r00**2)/(r00**2-r**2)
                kx1(1)=kx1(1)-rx*rr1
                ky1(1)=ky1(1)-ry*rr1
                kz1(1)=kz1(1)-rz*rr1
                U=U-(k*r00**2)*log(1-(r/r00)**2)/2
            ENDIF

            rx=x(n)-x(n-1)
            ry=y(n)-y(n-1)
            rz=z(n)-z(n-1)
            r=SQRT(rx**2+ry**2+rz**2)
            IF(r/=0.AND.r<=r00) THEN
                rr1=(k*r00**2)/(r00**2-r**2)
                kx1(n)=kx1(n)-rx*rr1
                ky1(n)=ky1(n)-ry*rr1
                kz1(n)=kz1(n)-rz*rr1
                U=U-(k*r00**2)*log(1-(r/r00)**2)/2
            ENDIF

            !!!!!! end UFENE(r) force


            !!!!!! ULJ(r) force

            DO j=1,n
                DO i=1,n
                    IF(i/=j)THEN
                        rx=x(i)-x(j)
                        ry=y(i)-y(j)
                        rz=z(i)-z(j)
                        r=SQRT(rx**2+ry**2+rz**2)
                        IF(r<=ld2)THEN
                            rr=(d0/r)**6
                            rr1=dble(48.0)*s0*(rr**2)/(r**2)
                            kx1(i)=kx1(i)+rx*rr1
                            ky1(i)=ky1(i)+ry*rr1
                            kz1(i)=kz1(i)+rz*rr1
                            U=U+4*s0*(d0/r)**12
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO

            !!!!!! end ULJ(r) force


            !!!!!!!bend energy!!!!!!!!!

            DO i=2,n-1
                rx1=x(i+1)-x(i)
                ry1=y(i+1)-y(i)
                rz1=z(i+1)-z(i)

                rx2=x(i)-x(i-1)
                ry2=y(i)-y(i-1)
                rz2=z(i)-z(i-1)

                kx1(i)=kx1(i)+kk*(rx1-rx2)
                ky1(i)=ky1(i)+kk*(ry1-ry2)
                kz1(i)=kz1(i)+kk*(rz1-rz2)
                U=U-kk*(rx1*rx2+ry1*ry2+rz1*rz2)
            ENDDO

            !!!!!!!!!!!!!!!!!!!end bend energy

            DO i=1,n
                ry=y(i)
                if(ry<=ld2)then
                    rr=(d0/ry)**6
                    rr1=dble(48)*s0*(rr**2)/ry
                    ky1(i)=ky1(i)+rr1
                endif
                ry=dble(ncelly)-y(i)
                if(ry<=ld2)then
                    rr=(d0/ry)**6
                    rr1=dble(48)*s0*(rr**2)/ry
                    ky1(i)=ky1(i)-rr1
                endif
            ENDDO

            DO i=1,n
                vx(i)=vx(i)+dble(0.5)*(kx1(i)+kx2(i))*tol
                vy(i)=vy(i)+dble(0.5)*(ky1(i)+ky2(i))*tol
                vz(i)=vz(i)+dble(0.5)*(kz1(i)+kz2(i))*tol
            ENDDO

            Ek1=0.0
            DO i=1,n
                Ek1=Ek1+dble(0.5)*m*(vx(i)**2+vy(i)**2+vz(i)**2)
            ENDDO
            t1=Ek1/(fac*dble(n)*m)
            scl=sqrt(temp/t1)
            do i=1,n
                vx(i)=scl*vx(i)
                vy(i)=scl*vy(i)
                vz(i)=scl*vz(i)
            end do

            DO i=1,n
                kx2(i)=kx1(i)
                ky2(i)=ky1(i)
                kz2(i)=kz1(i)
            ENDDO

        ENDDO        !!!!st3

        ranx=rand(idum)
        ranz=rand(idum)

        do i=1,n
            x(i)=x(i)+(ranx-0.5)*sid
            z(i)=z(i)+(ranz-0.5)*sid
        end do

        ppx=0.05
        ppz=0.05
        DO p=1,n

            IF(boxx*nint((x(p)-hboxx)/boxx)/=ppx.OR.boxz*nint((z(p)-hboxz)/boxz)/=ppz)THEN
                DO i=1,256
                    x_up(i)=x_up1(i)
                    z_up(i)=z_up1(i)
                    x_down(i)=x_down1(i)
                    z_down(i)=z_down1(i)
                ENDDO
                DO i=1,mn
                    x1(i)=x11(i)
                    y1(i)=y11(i)
                    z1(i)=z11(i)
                ENDDO

                ppx=boxx*nint((x(p)-hboxx)/boxx)
                ppz=boxz*nint((z(p)-hboxz)/boxz)

                DO i=1,256
                    x_up(i)=x_up(i)+ppx
                    z_up(i)=z_up(i)+ppz
                    x_down(i)=x_down(i)+ppx
                    z_down(i)=z_down(i)+ppz
                ENDDO

                DO i=1,mn
                    x1(i)=x1(i)+ppx
                    z1(i)=z1(i)+ppz
                ENDDO

                DO i=1,mn
                    y12(i)=y1(i)
                    x1(i)=x1(i)+vx1(i)*tol1-0.1*(y1(i)**2-15.0*y1(i))*uu*tol1
                    y1(i)=y1(i)+vy1(i)*tol1
                    z1(i)=z1(i)+vz1(i)*tol1
                    if(y1(i)>=boxy.OR.y1(i)<=0)then
                        y1(i)=y12(i)
                    endif
                ENDDO

                do i=1,mn
                    x1(i)=x1(i)+(ranx-0.5)*sid
                    z1(i)=z1(i)+(ranz-0.5)*sid
                end do

                do i=1,mn
                    x1(i)=x1(i)-boxx*nint((x1(i)-(ppx+hboxx))/boxx)
                    z1(i)=z1(i)-boxz*nint((z1(i)-(ppz+hboxz))/boxz)
                end do

                IF(q>=st0.AND.q<st00.AND.MOD(q,st22)==0)then
                    write(100,*)q,p
                ENDIF

                DO l=0,ncell-1
                    DO j=0,ncelly-1
                        DO ii=0,ncell-1

                            ki=ii+1+j*ncell+l*ncelly*ncell

                            !x10(ki)=ii*sid+sid/dble(2)
                            y10(ki)=j*sid+sid/dble(2)
                            !z10(ki)=l*sid+sid/dble(2)

                            vx0=0.0
                            vy0=0.0
                            vz0=0.0
                            ncc=0
                            DO i=1,n
                                IF(x(i)>=(ii*sid+ppx).AND.x(i)<((ii+1)*sid+ppx).AND.y(i)>=j*sid.AND.y(i)<(j+1)*sid.AND.z(i)>=(l*sid+ppz).AND.z(i)<((l+1)*sid+ppz))THEN
                                    vx0=vx0+vx(i)*m
                                    vy0=vy0+vy(i)*m
                                    vz0=vz0+vz(i)*m
                                    ncc=ncc+1
                                    gg(ncc)=i
                                ENDIF
                            ENDDO
                            nc=0
                            DO i=1,mn
                                IF(x1(i)>=(ii*sid+ppx).AND.x1(i)<((ii+1)*sid+ppx).AND.y1(i)>=j*sid.AND.y1(i)<(j+1)*sid.AND.z1(i)>=(l*sid+ppz).AND.z1(i)<((l+1)*sid+ppz))THEN
                                    vx0=vx0+vx1(i)*ms
                                    vy0=vy0+vy1(i)*ms
                                    vz0=vz0+vz1(i)*ms
                                    nc=nc+1
                                    g(nc)=i
                                ENDIF
                            ENDDO

                            if((nc==0).AND.(ncc==0))then
                                vxx(ki)=0.0
                                vyy(ki)=0.0
                                vzz(ki)=0.0
                            else
                                vxx(ki)=vx0/dble(nc*ms+ncc*m)
                                vyy(ki)=vy0/dble(nc*ms+ncc*m)
                                vzz(ki)=vz0/dble(nc*ms+ncc*m)
                            endif

                            IF(q>=st0.AND.q<st00.AND.MOD(q,st22)==0)then

                                write(100,*)vxx(ki)-0.1*(y10(ki)**2-15.0*y10(ki))*uu,vyy(ki),vzz(ki)

                            ENDIF

                            fai=2.0*pi*rand(idum)
                            theta=2.0*rand(idum)-1

                            lx=cos(fai)*SQRT(1-theta**2)
                            ly=sin(fai)*SQRT(1-theta**2)
                            lz=theta

                            DO t=1,ncc
                                i=gg(t)
                                vx(i)=vxx(ki)+(lx**2+(1-lx**2)*c)*(vx(i)-vxx(ki))+(lx*ly*(1-c)-lz*s1)*(vy(i)-vyy(ki))+(lx*lz*(1-c)+ly*s1)*(vz(i)-vzz(ki))
                                vy(i)=vyy(ki)+(lx*ly*(1-c)+lz*s1)*(vx(i)-vxx(ki))+(ly**2+(1-ly**2)*c)*(vy(i)-vyy(ki))+(ly*lz*(1-c)-lx*s1)*(vz(i)-vzz(ki))
                                vz(i)=vzz(ki)+(lx*lz*(1-c)-ly*s1)*(vx(i)-vxx(ki))+(ly*lz*(1-c)+lx*s1)*(vy(i)-vyy(ki))+(lz**2+(1-lz**2)*c)*(vz(i)-vzz(ki))
                            ENDDO
                            DO t=1,nc
                                i=g(t)
                                vx1(i)=vxx(ki)+(lx**2+(1-lx**2)*c)*(vx1(i)-vxx(ki))+(lx*ly*(1-c)-lz*s1)*(vy1(i)-vyy(ki))+(lx*lz*(1-c)+ly*s1)*(vz1(i)-vzz(ki))
                                vy1(i)=vyy(ki)+(lx*ly*(1-c)+lz*s1)*(vx1(i)-vxx(ki))+(ly**2+(1-ly**2)*c)*(vy1(i)-vyy(ki))+(ly*lz*(1-c)-lx*s1)*(vz1(i)-vzz(ki))
                                vz1(i)=vzz(ki)+(lx*lz*(1-c)-ly*s1)*(vx1(i)-vxx(ki))+(ly*lz*(1-c)+lx*s1)*(vy1(i)-vyy(ki))+(lz**2+(1-lz**2)*c)*(vz1(i)-vzz(ki))
                            ENDDO

                        ENDDO
                    ENDDO
                ENDDO

                Ek1=0.0
                DO i=1,n
                    Ek1=Ek1+dble(0.5)*m*(vx(i)**2+vy(i)**2+vz(i)**2)
                ENDDO
                t1=Ek1/(fac*dble(n)*m)
                scl=sqrt(temp/t1)
                do i=1,n
                    vx(i)=scl*vx(i)
                    vy(i)=scl*vy(i)
                    vz(i)=scl*vz(i)
                end do
                Ek2=0.0
                DO i=1,mn
                    Ek2=Ek2+dble(0.5)*ms*(vx1(i)**2+vy1(i)**2+vz1(i)**2)
                ENDDO
                t2=Ek2/(fac*dble(mn)*ms)
                scl=sqrt(temp/t2)
                do i=1,mn
                    vx1(i)=scl*vx1(i)
                    vy1(i)=scl*vy1(i)
                    vz1(i)=scl*vz1(i)
                end do
            ENDIF

        ENDDO

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!以下为所求!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF(q>=st0.AND.q<st00.AND.MOD(q,st01)==0)then

            call tran(q,filename)
            OPEN(203,file=filename)
            write(203,*) mn+512+n
            write(203,*)'cohar'

            do i=1,n

                write(203,'(a4,f17.7,f13.7,f13.7)')'p',x(i),y(i),z(i)

            enddo

            do i=1,256

                write(203,'(a4,f17.7,f13.7,f13.7)')'o',x_up(i),y_up(i),z_up(i)

            enddo

            do i=1,256

                write(203,'(a4,f17.7,f13.7,f13.7)')'o',x_down(i),y_down(i),z_down(i)

            enddo
            do i=1,mn
                write(203,'(a4,f17.7,f13.7,f13.7)') 'f',x1(i),y1(i),z1(i)
            end do
            CLOSE(203)


        ENDIF


        IF(MOD(q,st2)==0)THEN

            write(111,*)q
            DO i=2,n-1
                rx1=x(i)-x(i-1)
                ry1=y(i)-y(i-1)
                rz1=z(i)-z(i-1)
                rx2=x(i+1)-x(i)
                ry2=y(i+1)-y(i)
                rz2=z(i+1)-z(i)
                r1=SQRT(rx1**2+ry1**2+rz1**2)
                r2=SQRT(rx2**2+ry2**2+rz2**2)
                cc=abs(rx1*rx2+ry1*ry2+rz1*rz2)/(r1*r2)
                write(111,*)i,acos(cc)
            ENDDO

            write(80,*)q,U

            pcc1=0.0
            pcc2=0.0
            pcc3=0.0
            DO i=1,n-1
                rx=x(i)-x(i+1)
                ry=y(i)-y(i+1)
                rz=z(i)-z(i+1)

                r=SQRT(rx**2+ry**2+rz**2)
                c1=abs(rx)/r
                c2=abs(ry)/r
                c3=abs(rz)/r

                c12=c1**2
                c22=c2**2
                c32=c3**2
                pc1=3.0*c12/2.0-dble(0.5)
                pc2=3.0*c22/2.0-dble(0.5)
                pc3=3.0*c32/2.0-dble(0.5)
                pcc1=pcc1+pc1
                pcc2=pcc2+pc2
                pcc3=pcc3+pc3
            ENDDO

            pcc1=pcc1/dble(n-1)
            pcc2=pcc2/dble(n-1)
            pcc3=pcc3/dble(n-1)
            write(81,*)pcc1,pcc2,pcc3

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            DO j=1,ncelly
                num(j)=0
                DO i=1,n
                    if(y(i)>(j-1)*sid.AND.y(i)<=j*sid)then
                        num(j)=num(j)+1
                    endif
                ENDDO
                write(88,*)j,num(j)
            ENDDO

            xc=0.0
            yc=0.0
            zc=0.0
            DO i=1,n
                xc=xc+x(i)
                yc=yc+y(i)
                zc=zc+z(i)
            ENDDO
            xc=xc/dble(n)
            yc=yc/dble(n)
            zc=zc/dble(n)
            gxx=0.0
            gyy=0.0
            gzz=0.0
            rg=0.0
            DO i=1,n
                gxx=gxx+(x(i)-xc)**2
                gyy=gyy+(y(i)-yc)**2
                gzz=gzz+(z(i)-zc)**2
                rg=rg+(x(i)-xc)**2+(y(i)-yc)**2+(z(i)-zc)**2
            ENDDO
            gxx=gxx/dble(n)
            gyy=gyy/dble(n)
            gzz=gzz/dble(n)
            rg=rg/dble(n)
            write(821,*)q,gxx,gyy
            write(82,*)gzz,rg
            write(83,*)xc,yc,zc

            R22=(x(n)-x(1))**2+(y(n)-y(1))**2+(z(n)-z(1))**2
            write(90,*)q,R22

            ! compute shape factor
            s(:,:)=0.0
            DO i=1,n
                s(1,1)=s(1,1)+(x(i)-xc)**2
                s(1,2)=s(1,2)+(x(i)-xc)*(y(i)-yc)
                s(1,3)=s(1,3)+(x(i)-xc)*(z(i)-zc)
                s(2,1)=s(2,1)+(y(i)-yc)*(x(i)-xc)
                s(2,2)=s(2,2)+(y(i)-yc)**2
                s(2,3)=s(2,3)+(y(i)-yc)*(z(i)-zc)
                s(3,1)=s(3,1)+(z(i)-zc)*(x(i)-xc)
                s(3,2)=s(3,2)+(z(i)-zc)*(y(i)-yc)
                s(3,3)=s(3,3)+(z(i)-zc)**2
            ENDDO
            s(:,:)=s(:,:)/float(n)
            call jacobi(s,d)
            sf=1-3*(d(1)*d(2)+d(2)*d(3)+d(1)*d(3))/((d(1)+d(2)+d(3))**2)
            write(91,*)q,sf

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            pyx=dble(0.0)
            DO i=1,n
                pyx=pyx+m*abs(vy(i)*vx(i))
            ENDDO
            DO i=1,mn
                pyx=pyx+ms*abs(vy1(i)*vx1(i))
            ENDDO
            pyx=pyx/dble(ncell*ncelly*ncell)
            yta=pyx/dble(uu)
            write(84,*)q,yta    !!!粘滞系数84

        ENDIF


        ! compute shape factor
        s(:,:)=0.0
        DO i=1,n
            s(1,1)=s(1,1)+(x(i)-xc)**2
            s(1,2)=s(1,2)+(x(i)-xc)*(y(i)-yc)
            s(1,3)=s(1,3)+(x(i)-xc)*(z(i)-zc)
            s(2,1)=s(2,1)+(y(i)-yc)*(x(i)-xc)
            s(2,2)=s(2,2)+(y(i)-yc)**2
            s(2,3)=s(2,3)+(y(i)-yc)*(z(i)-zc)
            s(3,1)=s(3,1)+(z(i)-zc)*(x(i)-xc)
            s(3,2)=s(3,2)+(z(i)-zc)*(y(i)-yc)
            s(3,3)=s(3,3)+(z(i)-zc)**2
        ENDDO
        s(:,:)=s(:,:)/float(n)
        call jacobi(s,d)
        sf0(q)=1-3*(d(1)*d(2)+d(2)*d(3)+d(1)*d(3))/((d(1)+d(2)+d(3))**2)

        IF(q>1)THEN
            DO qq=q1,q-1
                if(((sf0(q-1)-sf0(qq))<-dble(0.5).AND.sf0(q)>sf0(q-1)).OR.((sf0(q-1)-sf0(qq))>dble(0.5).AND.sf0(q)<sf0(q-1)))then
                    write(110,*)q-1,q-1-q1
                    q1=q-1
                    exit
                endif
            ENDDO
        ENDIF

    ENDDO        !!!!!step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CLOSE(70)
    CLOSE(80)
    CLOSE(81)
    CLOSE(821)
    CLOSE(82)
    CLOSE(83)
    CLOSE(84)
    CLOSE(88)
    CLOSE(90)
    CLOSE(91)
    CLOSE(100)
    CLOSE(110)
    CLOSE(111)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
    FUNCTION rand(idum)
        INTEGER idum,ia,im,iq,ir,ntab,ndiv
        REAL rand,am,eps,rnmx
        PARAMETER (ia=16807,im=2147483647,am=1./im, iq=127773,ir=2836,  &
            ntab=32,ndiv=1+(im-1)/ntab, eps=1.2e-7,rnmx=1.-eps)
        INTEGER j,k,iv(ntab),iy
        SAVE iv,iy
        DATA iv /ntab*0/, iy /0/
        if (idum<=0.or.iy==0) then
            idum=max(-idum,1)
            do j=ntab+8,1,-1
                k=idum/iq
                idum=ia*(idum-k*iq)-ir*k
                if (idum<0) idum=idum+im
                if (j<=ntab) iv(j)=idum
            end do
            iy=iv(1)
        endif
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum<0) idum=idum+im
        j=1+iy/ndiv
        iy=iv(j)
        iv(j)=idum
        rand=min(am*iy,rnmx)
    END FUNCTION rand

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!形状因子子程序
    SUBROUTINE jacobi(s,d)
        real(8)::s(3,3),d(3)
        real(8)::bb(500),zz(500)
        integer nrot,ip,i,j,iq
        real(8) sm,tresh,gg,hh,tt,theta,cc,si,tau
        do ip=1,3
            bb(ip)=s(ip,ip)
            d(ip)=bb(ip)
            zz(ip)=0.0
        enddo
        nrot=0
        do i=1,50
            sm=0.0d0
            do ip=1,2
                do iq=ip+1,3
                    sm=sm+abs(s(ip,iq))
                enddo
            enddo
            if(sm.eq.0.0)return
            if(i.lt.4)then
                tresh=0.2*sm/9.0
            else
                tresh=0.0
            endif
            do ip=1,2
                do iq=ip+1,3
                    gg=100.0*abs(s(ip,iq))
                    if((i.gt.4).and.((abs(d(ip))+gg).eq.abs(d(ip))).and.((abs(d(iq))+gg).eq.abs(d(iq))))then
                        s(ip,iq)=0.0
                    else if(abs(s(ip,iq)).gt.tresh)then
                        hh=d(iq)-d(ip)
                        if((abs(hh)+gg).eq.abs(hh))then
                            tt=s(ip,iq)/hh
                        else
                            theta=0.5*hh/s(ip,iq)
                            tt=1.0/(abs(theta)+sqrt(1.0+theta**2))
                            if(theta.lt.0.0)tt=-tt
                        endif
                        cc=1.0/sqrt(1+tt**2)
                        si=tt*cc
                        tau=si/(1.0+cc)
                        hh=tt*s(ip,iq)
                        zz(ip)=zz(ip)-hh
                        zz(iq)=zz(iq)+hh
                        d(ip)=d(ip)-hh
                        d(iq)=d(iq)+hh
                        s(ip,iq)=0.0
                        do j=1,ip-1
                            gg=s(j,ip)
                            hh=s(j,iq)
                            s(j,ip)=gg-si*(hh+gg*tau)
                            s(j,iq)=hh+si*(gg-hh*tau)
                        enddo
                        do j=ip+1,iq-1
                            gg=s(ip,j)
                            hh=s(j,iq)
                            s(ip,j)=gg-si*(hh+gg*tau)
                            s(j,iq)=hh+si*(gg-hh*tau)
                        enddo
                        do j=iq+1,3
                            gg=s(ip,j)
                            hh=s(iq,j)
                            s(ip,j)=gg-si*(hh+gg*tau)
                            s(iq,j)=hh+si*(gg-hh*tau)
                        enddo
                        nrot=nrot+1
                    endif
                enddo
            enddo
            do ip=1,3
                bb(ip)=bb(ip)+zz(ip)
                d(ip)=bb(ip)
                zz(ip)=0.0
            enddo
        enddo
    END SUBROUTINE jacobi

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE xyz0(x,y,z,n)
        INTEGER n,f,i,j
        REAL(8) x(n),y(n),z(n),rx,ry,rz,box,r0
        PARAMETER(box=10,r0=7)

        DO i=2,n
            rx=x(i)-x(i-1)
            IF(abs(rx)>r0)THEN
                IF(rx>0)THEN
                    f=-1
                ELSE
                    f=1
                ENDIF
                DO j=i,n
                    x(j)=x(j)+box*f
                ENDDO
            ENDIF

            ry=y(i)-y(i-1)
            IF(abs(ry)>r0)THEN
                IF(ry>0)THEN
                    f=-1
                ELSE
                    f=1
                ENDIF
                DO j=i,n
                    y(j)=y(j)+box*f
                ENDDO
            ENDIF

            rz=z(i)-z(i-1)
            IF(abs(rz)>r0)THEN
                IF(rz>0)THEN
                    f=-1
                ELSE
                    f=1
                ENDIF
                DO j=i,n
                    z(j)=z(j)+box*f
                ENDDO
            ENDIF
        ENDDO
        RETURN
    END SUBROUTINE

    SUBROUTINE tran(q,string)
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
    End SUBROUTINE

end



