module typedef
    implicit none
    integer,parameter::filep=10
    integer,parameter::filef=11
    integer,parameter::filefv=12
    integer,parameter::filee=13
    integer,parameter::filev=14
    integer,parameter::np=80
    integer,parameter::avgn=5
    integer,parameter::a0=8
    integer,parameter::b0=5
    integer,parameter::a=31
    integer,parameter::b=2
    integer,parameter::c=8
    integer,parameter::nnc=(a0*b0*c)
    integer,parameter::nc=(2*a0*b0*c+a*b*c)
    integer,parameter::nf=(nc*avgn)
    integer,parameter::step=100000000
    integer,parameter::hn=5
    integer,parameter::nf0=19167
    real,parameter::len=1.0D-9 !(l=1)
    real,parameter::fr0=(0.2*len)
    real,parameter::x0=3.0
    real,parameter::x1=17.0
    real,parameter::y0=-b0
    real,parameter::z0=0.0
    real,parameter::hmd=1.0D-14 !(hp=2.152D-5)
    real,parameter::kb=1.380662D-23
    real,parameter::m=29.8D-27  !(m=1)
    real,parameter::t=1.0D-2
    real,parameter::e0=1.380662D-25  !(e0=kt=1)
    real,parameter::rp0=(0.5*len)
    real,parameter::r10=(1.5*rp0)
    real,parameter::eb0=1.380662D-22
    real,parameter::rb0=(0.5*len)
    real,parameter::k=4.14D-6 !(k=30)
    real,parameter::r00=(0.2*len)
    real,parameter::pi=3.14159
    real,parameter::rotang=(2*pi/3)
    real,parameter::mp=(m*avgn)
    real,parameter::hmpc=(hmd*hn)
    real,external::rand_num
    type r
        real rx
        real ry
        real rz
        integer ni
    end type
    type v
        real vx
        real vy
        real vz
    end type
    type f
        real fx
        real fy
        real fz
    end type
    type cell
        real x
        real y
        real z
        real avgvx
        real avgvy
        real avgvz
        real nm
        real ang
        real ranp
    end type
    type frv
        real frx,fry
        real vm
        real fvx,fvy
    end type
end module

module initializer
    use typedef
    implicit none
    type(r)::rp(np)
    type(v)::vp(np)
    type(r)::rf(nf)
    type(v)::vf(nf)
    type(v)::vff(nf)
    type(v)::vf0(nf)
    type(cell)::ce(nc)
    type(frv)::ff(nf0)
    integer::nz(np)
    integer flag,idum
contains
    subroutine initializerfp()
        implicit none
        integer i,j,k,count,temp,p,f
        real randn,rmin,dr,dr0,rr
        rmin=0.1
        dr=0.5
        dr0=1.0
        count=1
        do i=1,a0
            do j=1,b0
                do k=1,c
                    ce(count)%x=(x0+(i-1))*len
                    ce(count)%y=(y0+j-1)*len
                    ce(count)%z=(z0+k-1)*len
                    ce(count)%nm=0
                    ce(count)%avgvx=0.0
                    ce(count)%avgvy=0.0
                    ce(count)%avgvz=0.0
                    ce(count)%ang=0.0
                    ce(count)%ranp=0.0
                    count=count+1
                end do
            end do
        end do
        do i=1,a0
            do j=1,b0
                do k=1,c
                    ce(count)%x=(x1+(i-1))*len
                    ce(count)%y=(y0+j-1)*len
                    ce(count)%z=(z0+k-1)*len
                    ce(count)%nm=0
                    ce(count)%avgvx=0.0
                    ce(count)%avgvy=0.0
                    ce(count)%avgvz=0.0
                    ce(count)%ang=0.0
                    ce(count)%ranp=0.0
                    count=count+1
                end do
            end do
        end do
        do i=1,a
            do j=1,b
                do k=1,c
                    ce(count)%x=(0.0+i-1)*len
                    ce(count)%y=(0.0+j-1)*len
                    ce(count)%z=(0.0+k-1)*len
                    ce(count)%nm=0
                    ce(count)%avgvx=0.0
                    ce(count)%avgvy=0.0
                    ce(count)%avgvz=0.0
                    ce(count)%ang=0.0
                    ce(count)%ranp=0.0
                    count=count+1
                end do
            end do
        end do
        temp=(nnc)*(avgn)
        do i=1,temp
            randn=rand_num(idum)
            rf(i)%rx=(x0+randn*a0)*len
            randn=rand_num(idum)
            rf(i)%ry=(y0+randn*b0)*len
            randn=rand_num(idum)
            rf(i)%rz=(z0+randn*c)*len
            rf(i)%ni=0
        end do
        count=(nnc)*(avgn)
        temp=(nnc)*(avgn)
        do i=1,temp
            randn=rand_num(idum)
            rf(i+count)%rx=(x1+randn*a0)*len
            randn=rand_num(idum)
            rf(i+count)%ry=(y0+randn*b0)*len
            randn=rand_num(idum)
            rf(i+count)%rz=(z0+randn*c)*len
            rf(i+count)%ni=0
        end do
        count=2*(nnc)*(avgn)
        temp=((nc)-2*(nnc))*(avgn)
        do i=1,temp
            randn=rand_num(idum)
            rf(i+count)%rx=(0.0+randn*a)*len
            randn=rand_num(idum)
            rf(i+count)%ry=(0.0+randn*b)*len
            randn=rand_num(idum)
            rf(i+count)%rz=(0.0+randn*c)*len
            rf(i+count)%ni=0
        end do
        do i=1,nf
            vf(i)%vx=0.0
            vf(i)%vy=0.0
            vf(i)%vz=0.0
        end do
        rp(1)%rx=x0+4.0
        rp(1)%ry=y0+2.0
        rp(1)%rz=z0+4.0
        do i=2,np
            f=1
            do while(f==1)
                rp(i)%rx=rp(i-1)%rx
                rp(i)%ry=rp(i-1)%ry
                rp(i)%rz=rp(i-1)%rz
                randn=rand_num(idum)
                p=int(randn*6.0)
                select case(p)
                    case(0)
                        rp(i)%rx=rp(i)%rx-dr
                    case(1)
                        rp(i)%rx=rp(i)%rx+dr
                    case(2)
                        rp(i)%ry=rp(i)%ry-dr
                    case(3)
                        rp(i)%ry=rp(i)%ry+dr
                    case(4)
                        rp(i)%rz=rp(i)%rz-dr
                    case(5)
                        rp(i)%rz=rp(i)%rz+dr
                    case default
                        rp(i)%rx=rp(i)%rx
                        rp(i)%ry=rp(i)%ry
                        rp(i)%rz=rp(i)%rz
                end select
                if((rp(i)%rx<(x0+dr0)).or.(rp(i)%rx>(x0+a0-dr0)).or.(rp(i)%ry<(y0+dr0)).or.&
                        (rp(i)%ry>(y0+b0-dr0)).or.(rp(i)%rz<(z0+dr0)).or.(rp(i)%rz>(z0+c-dr0))) then
                    cycle
                end if
                do j=i-1,1,-1
                    rr=sqrt((rp(i)%rx-rp(j)%rx)**2+(rp(i)%ry-rp(j)%ry)**2+(rp(i)%rz-rp(j)%rz)**2)
                    if(rr<rmin) exit
                end do
                if(j>0) cycle
                f=0
            end do
        end do
        do i=1,np
            rp(i)%rx=rp(i)%rx*len
            rp(i)%ry=rp(i)%ry*len
            rp(i)%rz=rp(i)%rz*len
            rp(i)%ni=0
            nz(i)=0
            vp(i)%vx=0.0
            vp(i)%vy=0.0
            vp(i)%vz=0.0
        end do
        return
    end subroutine

    subroutine initialvelocity(vi,n,mm)
        implicit none
        integer i,n
        type(v)::vi(n)
        real sumvx,sumvy,sumvz,tempv,tempv2,randn,tempx,tempy,tempz,mm
        tempv2=kb*t/mm
        tempv=sqrt(tempv2)
        do i=1,n
            randn=rand_num(idum)
            if(randn>=0.50) then
                vi(i)%vx=2*tempv*sqrt(randn-0.50)
            else
                vi(i)%vx=-2*tempv*sqrt(0.50-randn)
            end if
            randn=rand_num(idum)
            if(randn>=0.50) then
                vi(i)%vy=2*tempv*sqrt(randn-0.50)
            else
                vi(i)%vy=-2*tempv*sqrt(0.50-randn)
            end if
            randn=rand_num(idum)
            if(randn>=0.50) then
                vi(i)%vz=2*tempv*sqrt(randn-0.50)
            else
                vi(i)%vz=-2*tempv*sqrt(0.50-randn)
            end if
        end do
        sumvx=0.0
        sumvy=0.0
        sumvz=0.0
        do i=1,n
            sumvx=sumvx+vi(i)%vx
            sumvy=sumvy+vi(i)%vy
            sumvz=sumvz+vi(i)%vz
        end do
        do i=1,n
            vi(i)%vx=vi(i)%vx-sumvx/(n*1.0)
            vi(i)%vy=vi(i)%vy-sumvy/(n*1.0)
            vi(i)%vz=vi(i)%vz-sumvz/(n*1.0)
        end do
        sumvx=0.0
        sumvy=0.0
        sumvz=0.0
        do i=1,n
            sumvx=sumvx+(vi(i)%vx)**2
            sumvy=sumvy+(vi(i)%vy)**2
            sumvz=sumvz+(vi(i)%vz)**2
        end do
        tempx=sqrt(tempv2*n*1.0/sumvx)
        tempy=sqrt(tempv2*n*1.0/sumvy)
        tempz=sqrt(tempv2*n*1.0/sumvz)
        do i=1,n
            vi(i)%vx=(vi(i)%vx)*tempx
            vi(i)%vy=(vi(i)%vy)*tempy
            vi(i)%vz=(vi(i)%vz)*tempz
        end do
        return
    end subroutine

    subroutine fluidvelocity(tempr,tempv)
        implicit none
        type(r)::tempr
        type(v)::tempv
        integer j,num,nn
        real(8)::rr
        tempv%vx=0.0
        tempv%vy=0.0
        tempv%vz=0.0
        num=0
        do j=1,nf0
            rr=sqrt((tempr%rx-ff(j)%frx)**2+(tempr%ry-ff(j)%fry)**2)
            if(rr<=fr0) then
                num=num+1
                tempv%vx=tempv%vx+ff(j)%fvx
                tempv%vy=tempv%vy+ff(j)%fvy
            end if
        end do
        if(num==0) then
            tempv%vx=0
            tempv%vy=0
        else
            tempv%vx=tempv%vx/num
            tempv%vy=tempv%vy/num
        end if
        return
    end subroutine

    subroutine fvelocity(tempr1,tempv1)
        implicit none
        type(r)::tempr1
        type(v)::tempv1
        integer i
        real randn
        tempv1%vx=0.0
        tempv1%vy=0.0
        tempv1%vz=0.0
        call fluidvelocity(tempr1,tempv1)
        randn=rand_num(idum)
        i=int(randn*nf)
        tempv1%vx=tempv1%vx+vf0(i)%vx
        tempv1%vy=tempv1%vy+vf0(i)%vy
        tempv1%vz=tempv1%vz+vf0(i)%vz
        return
    end subroutine
end module

module boundary
    use initializer
    implicit none
    real::oldry(nf)
    type(f)::fb(np)
    real(8) Eb
contains
    subroutine boundaryconditionf()
        implicit none
        integer i
        do i=1,nf
            if((rf(i)%rz)<0.0) then
                rf(i)%rz=(c*1.0)*len+(rf(i)%rz)
                vf(i)%vz=vf(i)%vz
            else if((rf(i)%rz)>(c*1.0)*len) then
                rf(i)%rz=(rf(i)%rz)-c*1.0*len
                vf(i)%vz=vf(i)%vz
            else
                rf(i)%rz=rf(i)%rz
                vf(i)%vz=vf(i)%vz
            end if
            if(oldry(i)>=0.0) then
                if((rf(i)%rx)<0.0) then
                    rf(i)%rx=a*len+rf(i)%rx
                else if((rf(i)%rx)>(a*1.0)*len) then
                    rf(i)%rx=(rf(i)%rx)-a*1.0*len
                else
                    rf(i)%rx=rf(i)%rx
                end if
                if((rf(i)%ry)>(b*1.0*len)) then
                    rf(i)%ry=b*len-((rf(i)%ry)-b*len)
                    vf(i)%vy=-vf(i)%vy
                else if(((rf(i)%ry)<0.0).and.(.not.((((rf(i)%rx)>=x0*len).and.((rf(i)%rx)<=(x0+a0)*len))&
                        .or.(((rf(i)%rx)>=x1*len).and.((rf(i)%rx)<=(x1+a0)*len))))) then
                    rf(i)%ry=-rf(i)%ry
                    vf(i)%vy=-vf(i)%vy
                else
                    rf(i)%ry=rf(i)%ry
                    vf(i)%vy=vf(i)%vy
                end if
            else
                if((rf(i)%rx)<x0*len) then
                    rf(i)%rx=x0*len+(x0*len-(rf(i)%rx))
                    vf(i)%vx=-vf(i)%vx
                else if(((rf(i)%rx)>(x0+a0)*len).and.((rf(i)%rx)<len*(x1+x0+a0)/2)) then
                    rf(i)%rx=(x0+a0)*len-((rf(i)%rx)-(x0+a0)*len)
                    vf(i)%vx=-vf(i)%vx
                else if(((rf(i)%rx)<x1*len).and.((rf(i)%rx)>=len*(x1+x0+a0)/2)) then
                    rf(i)%rx=x1*len+(x1*len-(rf(i)%rx))
                    vf(i)%vx=-vf(i)%vx
                else if((rf(i)%rx)>(x1+a0)*len) then
                    rf(i)%rx=(x1+a0)*len-((rf(i)%rx)-(x1+a0)*len)
                    vf(i)%vx=-vf(i)%vx
                else
                    rf(i)%rx=rf(i)%rx
                    vf(i)%vx=vf(i)%vx
                end if
                if((rf(i)%ry)<y0*len) then
                    rf(i)%ry=y0*len+(y0*len-(rf(i)%ry))
                    vf(i)%vy=-vf(i)%vy
                else
                    rf(i)%ry=rf(i)%ry
                    vf(i)%vy=vf(i)%vy
                end if
            end if
        end do
        return
    end subroutine

    subroutine boundaryconditionff()
        implicit none
        integer i
        do i=1,nf
            if((rf(i)%rz)<0.0) then
                rf(i)%rz=(c*1.0)*len+(rf(i)%rz)
                vf(i)%vz=vf(i)%vz
            else if((rf(i)%rz)>(c*1.0)*len) then
                rf(i)%rz=(rf(i)%rz)-c*1.0*len
                vf(i)%vz=vf(i)%vz
            else
                rf(i)%rz=rf(i)%rz
                vf(i)%vz=vf(i)%vz
            end if
            if(oldry(i)>=0.0) then
                if((rf(i)%ry)>(b*1.0*len)) then
                    rf(i)%ry=b*len-((rf(i)%ry)-b*len)
                    vf(i)%vy=-vf(i)%vy
                else if(((rf(i)%ry)<0.0).and.(.not.((((rf(i)%rx)>=x0*len).and.((rf(i)%rx)<=(x0+a0)*len))&
                        .or.(((rf(i)%rx)>=x1*len).and.((rf(i)%rx)<=(x1+a0)*len))))) then
                    rf(i)%ry=-rf(i)%ry
                    vf(i)%vy=-vf(i)%vy
                else
                    rf(i)%ry=rf(i)%ry
                    vf(i)%vy=vf(i)%vy
                end if
                if((rf(i)%rx)<0.0) then
                    rf(i)%rx=-rf(i)%rx
                    call fvelocity(rf(i),vf(i))
                else if((rf(i)%rx)>(a*1.0)*len) then
                    rf(i)%rx=(rf(i)%rx)-a*1.0*len
                    call fvelocity(rf(i),vf(i))
                else
                    rf(i)%rx=rf(i)%rx
                    vf(i)%vx=vf(i)%vx
                end if
            else
                if((rf(i)%rx)<x0*len) then
                    rf(i)%rx=x0*len+(x0*len-(rf(i)%rx))
                    vf(i)%vx=-vf(i)%vx
                else if(((rf(i)%rx)>(x0+a0)*len).and.((rf(i)%rx)<len*(x1+x0+a0)/2)) then
                    rf(i)%rx=(x0+a0)*len-((rf(i)%rx)-(x0+a0)*len)
                    vf(i)%vx=-vf(i)%vx
                else if(((rf(i)%rx)<x1*len).and.((rf(i)%rx)>=len*(x1+x0+a0)/2)) then
                    rf(i)%rx=x1*len+(x1*len-(rf(i)%rx))
                    vf(i)%vx=-vf(i)%vx
                else if((rf(i)%rx)>(x1+a0)*len) then
                    rf(i)%rx=(x1+a0)*len-((rf(i)%rx)-(x1+a0)*len)
                    vf(i)%vx=-vf(i)%vx
                else
                    rf(i)%rx=rf(i)%rx
                    vf(i)%vx=vf(i)%vx
                end if
                if((rf(i)%ry)<y0*len) then
                    rf(i)%ry=y0*len+(y0*len-(rf(i)%ry))
                    vf(i)%vy=-vf(i)%vy
                else
                    rf(i)%ry=rf(i)%ry
                    vf(i)%vy=vf(i)%vy
                end if
            end if
        end do
        return
    end subroutine

    subroutine boundaryconditionp()
        implicit none
        integer i
        real co,si,rr
        do i=1,np
            fb(i)%fx=0.0
            fb(i)%fy=0.0
            fb(i)%fz=0.0
        end do
        Eb=0.0
        do i=1,np
            if(rp(i)%ry>=0) then
                if((rp(i)%ry<2**(1.0/6)*rb0).and.(.not.((((rp(i)%rx)>=x0*len).and.&
                        ((rp(i)%rx)<=(x0+a0)*len)).or.(((rp(i)%rx)>=x1*len).and.((rp(i)%rx)<=(x1+a0)*len))))) then
                    fb(i)%fy=(eb0/rp(i)%ry)*(12*(rb0/rp(i)%ry)**12-6*(rb0/rp(i)%ry)**6)
                    Eb=Eb+(eb0)*((rb0/rp(i)%ry)**12-(rb0/rp(i)%ry)**6)
                else if((rp(i)%ry<2**(1.0/6)*rb0).and.((((rp(i)%rx)>=x0*len).and.&
                        ((rp(i)%rx)<=(x0+a0)*len)).or.(((rp(i)%rx)>=x1*len).and.((rp(i)%rx)<=(x1+a0)*len)))) then
                    if((rp(i)%rx-x0*len)**2+rp(i)%ry**2<=(2**(1.0/6)*rb0)**2) then
                        rr=sqrt((rp(i)%rx-x0*len)**2+rp(i)%ry**2)
                        co=(rp(i)%rx-x0*len)/rr
                        si=rp(i)%ry/rr
                        fb(i)%fx=co*(eb0/rr)*(12*(rb0/rr)**12-6*(rb0/rr)**6)
                        fb(i)%fy=si*(eb0/rr)*(12*(rb0/rr)**12-6*(rb0/rr)**6)
                        Eb=Eb+(eb0)*((rb0/rr)**12-(rb0/rr)**6)
                    else if((rp(i)%rx-(x0+a0)*len)**2+rp(i)%ry**2<=(2**(1.0/6)*rb0)**2) then
                        rr=sqrt(((x0+a0)*len-rp(i)%rx)**2+rp(i)%ry**2)
                        co=((x0+a0)*len-rp(i)%rx)/rr
                        si=rp(i)%ry/rr
                        fb(i)%fx=-co*(eb0/rr)*(12*(rb0/rr)**12-6*(rb0/rr)**6)
                        fb(i)%fy=si*(eb0/rr)*(12*(rb0/rr)**12-6*(rb0/rr)**6)
                        Eb=Eb+(eb0)*((rb0/rr)**12-(rb0/rr)**6)
                    else if((rp(i)%rx-x1*len)**2+rp(i)%ry**2<=(2**(1.0/6)*rb0)**2) then
                        rr=sqrt((rp(i)%rx-x1*len)**2+rp(i)%ry**2)
                        co=(rp(i)%rx-x1*len)/rr
                        si=rp(i)%ry/rr
                        fb(i)%fx=co*(eb0/rr)*(12*(rb0/rr)**12-6*(rb0/rr)**6)
                        fb(i)%fy=si*(eb0/rr)*(12*(rb0/rr)**12-6*(rb0/rr)**6)
                        Eb=Eb+(eb0)*((rb0/rr)**12-(rb0/rr)**6)
                    else if((rp(i)%rx-(x1+a0)*len)**2+rp(i)%ry**2<=(2**(1.0/6)*rb0)**2) then
                        rr=sqrt(((x1+a0)*len-rp(i)%rx)**2+rp(i)%ry**2)
                        co=((x1+a0)*len-rp(i)%rx)/rr
                        si=rp(i)%ry/rr
                        fb(i)%fx=-co*(eb0/rr)*(12*(rb0/rr)**12-6*(rb0/rr)**6)
                        fb(i)%fy=si*(eb0/rr)*(12*(rb0/rr)**12-6*(rb0/rr)**6)
                        Eb=Eb+(eb0)*((rb0/rr)**12-(rb0/rr)**6)
                    else
                        fb(i)%fx=0.0
                        fb(i)%fy=0.0
                    end if
                else if(rp(i)%ry>(b*len-2**(1.0/6)*rb0)) then
                    fb(i)%fy=-(eb0/(b*len-rp(i)%ry))*(12*(rb0/(b*len-rp(i)%ry))**12-6*(rb0/(b*len-rp(i)%ry))**6)
                    Eb=Eb+(eb0)*((rb0/(b*len-rp(i)%ry))**12-(rb0/(b*len-rp(i)%ry))**6)
                else
                    fb(i)%fy=0.0
                end if
            else
                if(rp(i)%rx<(x0*len+2**(1.0/6)*rb0)) then
                    fb(i)%fx=(eb0/(rp(i)%rx-x0*len))*(12*(rb0/(rp(i)%rx-x0*len))**12-6*(rb0/(rp(i)%rx-x0*len))**6)
                    Eb=Eb+(eb0)*((rb0/(rp(i)%rx-x0*len))**12-(rb0/(rp(i)%rx-x0*len))**6)
                else if(rp(i)%rx>((x0+a0)*len-2**(1.0/6)*rb0)) then
                    fb(i)%fx=-(eb0/((x0+a0)*len-rp(i)%rx))*(12*(rb0/((x0+a0)*len-rp(i)%rx))**12-6*(rb0/((x0+a0)*len-rp(i)%rx))**6)
                    Eb=Eb+(eb0)*((rb0/((x0+a0)*len-rp(i)%rx))**12-(rb0/((x0+a0)*len-rp(i)%rx))**6)
                else if(rp(i)%rx<(x1*len+2**(1.0/6)*rb0)) then
                    fb(i)%fx=(eb0/(rp(i)%rx-x1*len))*(12*(rb0/(rp(i)%rx-x1*len))**12-6*(rb0/(rp(i)%rx-x1*len))**6)
                    Eb=Eb+(eb0)*((rb0/(rp(i)%rx-x1*len))**12-(rb0/(rp(i)%rx-x1*len))**6)
                else if(rp(i)%rx>((x1+a0)*len-2**(1.0/6)*rb0)) then
                    fb(i)%fx=-(eb0/((x1+a0)*len-rp(i)%rx))*(12*(rb0/((x1+a0)*len-rp(i)%rx))**12-6*(rb0/((x1+a0)*len-rp(i)%rx))**6)
                    Eb=Eb+(eb0)*((rb0/((x1+a0)*len-rp(i)%rx))**12-(rb0/((x1+a0)*len-rp(i)%rx))**6)
                else
                    fb(i)%fx=0.0
                end if
                if(rp(i)%ry<(y0*len+2**(1.0/6)*rb0)) then
                    fb(i)%fy=(eb0/(rp(i)%ry-y0*len))*(12*(rb0/(rp(i)%ry-y0*len))**12-6*(rb0/(rp(i)%ry-y0*len))**6)
                    Eb=Eb+(eb0)*((rb0/(rp(i)%ry-y0*len))**12-(rb0/(rp(i)%ry-y0*len))**6)
                else
                    fb(i)%fy=0.0
                end if
            end if
        end do
        return
    end subroutine

    subroutine boundaryconditionprz()
        implicit none
        integer i
        do i=1,np
            if((rp(i)%rz)<0.0) then
                rp(i)%rz=c*1.0*len+(rp(i)%rz)
                nz(i)=nz(i)-1
            else if((rp(i)%rz)>(c*1.0*len)) then
                rp(i)%rz=(rp(i)%rz)-c*1.0*len
                nz(i)=nz(i)+1
            else
                rp(i)%rz=rp(i)%rz
                nz(i)=nz(i)
            end if
        end do
        return
    end subroutine

    subroutine mpccondition(rs,n,ty)
        implicit none
        integer n
        type(r)::rs(n)
        real::ty(n)
        integer i
        do i=1,n
            if(rs(i)%rz<0) then
                rs(i)%rz=c*len+rs(i)%rz
            else if(rs(i)%rz>(c*len)) then
                rs(i)%rz=rs(i)%rz-c*len
            else
                rs(i)%rz=rs(i)%rz
            end if
            if(ty(i)>=0) then
                if(rs(i)%rx<0) then
                    rs(i)%rx=a*len+rs(i)%rx
                else if(rs(i)%rx>(a*len)) then
                    rs(i)%rx=rs(i)%rx-a*len
                else
                    rs(i)%rx=rs(i)%rx
                end if
                if(rs(i)%ry<0) then
                    rs(i)%ry=b*len+rs(i)%ry
                else if(rs(i)%ry>(b*len)) then
                    rs(i)%ry=rs(i)%ry-b*len
                else
                    rs(i)%ry=rs(i)%ry
                end if
            else
                if(rs(i)%rx<(x0*len)) then
                    rs(i)%rx=rs(i)%rx+a0*len
                else if((rs(i)%rx>((x0+a0)*len)).and.(rs(i)%rx<((x1+x0+a0)*len/2))) then
                    rs(i)%rx=rs(i)%rx-a0*len
                else if((rs(i)%rx<(x1*len)).and.(rs(i)%rx>((x1+x0+a0)*len/2))) then
                    rs(i)%rx=rs(i)%rx+a0*len
                else if(rs(i)%rx>((x1+a0)*len)) then
                    rs(i)%rx=rs(i)%rx-a0*len
                else
                    rs(i)%rx=rs(i)%rx
                end if
                if(rs(i)%ry<(y0*len)) then
                    rs(i)%ry=rs(i)%ry+b0*len
                else if(rs(i)%ry>((y0+b0)*len)) then
                    rs(i)%ry=rs(i)%ry-b0*len
                else
                    rs(i)%ry=rs(i)%ry
                end if
            end if
        end do
        return
    end subroutine
end module

module mpc
    use initializer
    use boundary
    implicit none
contains
    subroutine mpcmethod()
        implicit none
        integer i,j
        real randn,lx,ly,lz,co,si,rx0,ry0,rz0
        real::tempv(3,1),rm(3,3),tempvm(3,1),tempm(3,1),tempv1(3,1)
        real::tempy(nf),tempy1(np)
        do i=1,nc
            ce(i)%nm=0.0
            ce(i)%avgvx=0.0
            ce(i)%avgvy=0.0
            ce(i)%avgvz=0.0
            ce(i)%ang=0.0
            ce(i)%ranp=0.0
        end do
        randn=rand_num(idum)
        rx0=-0.5+randn
        randn=rand_num(idum)
        ry0=-0.5+randn
        randn=rand_num(idum)
        rz0=-0.5+randn
        do i=1,nf
            tempy(i)=rf(i)%ry
        end do
        do i=1,nf
            rf(i)%rx=rf(i)%rx+rx0*len
            rf(i)%ry=rf(i)%ry+ry0*len
            rf(i)%rz=rf(i)%rz+rz0*len
        end do
        call mpccondition(rf,nf,tempy)
        !do i=1,np
        !tempy1(i)=rp(i)%ry
        !end do
        !do i=1,np
        !rp(i)%rx=rp(i)%rx+rx0*len
        !rp(i)%ry=rp(i)%ry+ry0*len
        !rp(i)%rz=rp(i)%rz+rz0*len
        !end do
        !call mpccondition(rp,np,tempy1)
        do i=1,nc
            do j=1,nf
                if((rf(j)%rx>=ce(i)%x).and.(rf(j)%rx<(ce(i)%x+1.0*len)).and.(rf(j)%ry>=ce(i)%y).and.&
                        (rf(j)%ry<(ce(i)%y+1.0*len)).and.(rf(j)%rz>=ce(i)%z).and.(rf(j)%rz<(ce(i)%z+1.0*len))) then
                    rf(j)%ni=i
                    ce(i)%nm=ce(i)%nm+m
                    ce(i)%avgvx=ce(i)%avgvx+m*vf(j)%vx
                    ce(i)%avgvy=ce(i)%avgvy+m*vf(j)%vy
                    ce(i)%avgvz=ce(i)%avgvz+m*vf(j)%vz
                end if
            end do
            do j=1,np
                if((rp(j)%rx>=ce(i)%x).and.(rp(j)%rx<(ce(i)%x+1.0*len)).and.(rp(j)%ry>=ce(i)%y).and.&
                        (rp(j)%ry<(ce(i)%y+1.0*len)).and.(rp(j)%rz>=ce(i)%z).and.(rp(j)%rz<(ce(i)%z+1.0*len))) then
                    rp(j)%ni=i
                    ce(i)%nm=ce(i)%nm+mp
                    ce(i)%avgvx=ce(i)%avgvx+mp*vp(j)%vx
                    ce(i)%avgvy=ce(i)%avgvy+mp*vp(j)%vy
                    ce(i)%avgvz=ce(i)%avgvz+mp*vp(j)%vz
                end if
            end do
            if(ce(i)%nm==0) then
                ce(i)%avgvx=0
                ce(i)%avgvy=0
                ce(i)%avgvz=0
            else
                ce(i)%avgvx=ce(i)%avgvx/ce(i)%nm
                ce(i)%avgvy=ce(i)%avgvy/ce(i)%nm
                ce(i)%avgvz=ce(i)%avgvz/ce(i)%nm
            end if
        end do
        do i=1,nc
            randn=rand_num(idum)
            ce(i)%ang=2*pi*randn
        end do
        do i=1,nc
            randn=rand_num(idum)
            ce(i)%ranp=-1.0+2*randn
        end do
        do i=1,nf
            tempv(1,1)=vf(i)%vx
            tempv(2,1)=vf(i)%vy
            tempv(3,1)=vf(i)%vz
            tempvm(1,1)=ce(rf(i)%ni)%avgvx
            tempvm(2,1)=ce(rf(i)%ni)%avgvy
            tempvm(3,1)=ce(rf(i)%ni)%avgvz
            lx=cos(ce(rf(i)%ni)%ang)*sqrt(1-(ce(rf(i)%ni)%ranp)**2)
            ly=sin(ce(rf(i)%ni)%ang)*sqrt(1-(ce(rf(i)%ni)%ranp)**2)
            lz=ce(rf(i)%ni)%ranp
            co=cos(rotang)
            si=sin(rotang)
            rm(1,1)=lx**2+(1-lx**2)*co
            rm(2,1)=lx*ly*(1-co)+lz*si
            rm(3,1)=lx*lz*(1-co)-ly*si
            rm(1,2)=lx*ly*(1-co)-lz*si
            rm(2,2)=ly**2+(1-ly**2)*co
            rm(3,2)=ly*lz*(1-co)+lx*si
            rm(1,3)=lx*lz*(1-co)+ly*si
            rm(2,3)=ly*lz*(1-co)-lx*si
            rm(3,3)=lz**2+(1-lz**2)*co
            tempv1=tempv-tempvm
            tempm=tempvm+matmul(rm,tempv1)
            vf(i)%vx=tempm(1,1)
            vf(i)%vy=tempm(2,1)
            vf(i)%vz=tempm(3,1)
        end do
        do i=1,np
            tempv(1,1)=vp(i)%vx
            tempv(2,1)=vp(i)%vy
            tempv(3,1)=vp(i)%vz
            tempvm(1,1)=ce(rp(i)%ni)%avgvx
            tempvm(2,1)=ce(rp(i)%ni)%avgvy
            tempvm(3,1)=ce(rp(i)%ni)%avgvz
            lx=cos(ce(rp(i)%ni)%ang)*sqrt(1-(ce(rp(i)%ni)%ranp)**2)
            ly=sin(ce(rp(i)%ni)%ang)*sqrt(1-(ce(rp(i)%ni)%ranp)**2)
            lz=ce(rp(i)%ni)%ranp
            co=cos(rotang)
            si=sin(rotang)
            rm(1,1)=lx**2+(1-lx**2)*co
            rm(2,1)=lx*ly*(1-co)+lz*si
            rm(3,1)=lx*lz*(1-co)-ly*si
            rm(1,2)=lx*ly*(1-co)-lz*si
            rm(2,2)=ly**2+(1-ly**2)*co
            rm(3,2)=ly*lz*(1-co)+lx*si
            rm(1,3)=lx*lz*(1-co)+ly*si
            rm(2,3)=ly*lz*(1-co)-lx*si
            rm(3,3)=lz**2+(1-lz**2)*co
            tempv1=tempv-tempvm
            tempm=tempvm+matmul(rm,tempv1)
            vp(i)%vx=tempm(1,1)
            vp(i)%vy=tempm(2,1)
            vp(i)%vz=tempm(3,1)
        end do
        call ntv(vp,np,mp)
        return
    end subroutine
end module

module force
    use boundary
    implicit none
    type(f)::ftot(np)
    real E
contains
    subroutine totalforce()
        implicit none
        integer i,j
        real rxx,ryy,rzz,rr,rr0,rr1,temp
        real::x(np),y(np),z(np)
        do i=1,np
            x(i)=rp(i)%rx
            y(i)=rp(i)%ry
            z(i)=rp(i)%rz+c*nz(i)*len
        end do
        do i=1,np
            ftot(i)%fx=0.0
            ftot(i)%fy=0.0
            ftot(i)%fz=0.0
        end do
        E=0.0
        !!!!!! UFENE(r) force

        DO i=2,np-1
            DO j=i+1,i-1,-2
                rxx=x(i)-x(j)
                ryy=y(i)-y(j)
                rzz=z(i)-z(j)

                rr=SQRT(rxx**2+ryy**2+rzz**2)

                IF(rr/=0.AND.rr<=r10) THEN
                    rr1=(k*r00**2)/(r10**2-rr**2)
                    ftot(i)%fx=ftot(i)%fx-rxx*rr1
                    ftot(i)%fy=ftot(i)%fy-ryy*rr1
                    ftot(i)%fz=ftot(i)%fz-rzz*rr1
                ENDIF
            ENDDO
        ENDDO

        rxx=x(1)-x(2)
        ryy=y(1)-y(2)
        rzz=z(1)-z(2)

        rr=SQRT(rxx**2+ryy**2+rzz**2)

        IF(rr/=0.AND.rr<=r10) THEN
            rr1=(k*r10**2)/(r10**2-rr**2)
            E=E-0.5*k*r10**2*log(1-(rr/r10)**2)
            ftot(1)%fx=ftot(1)%fx-rxx*rr1
            ftot(1)%fy=ftot(1)%fy-ryy*rr1
            ftot(1)%fz=ftot(1)%fz-rzz*rr1
        ENDIF

        rxx=x(np)-x(np-1)
        ryy=y(np)-y(np-1)
        rzz=z(np)-z(np-1)

        rr=SQRT(rxx**2+ryy**2+rzz**2)

        IF(rr/=0.AND.rr<=r10) THEN
            rr1=(k*r10**2)/(r10**2-rr**2)
            ftot(np)%fx=ftot(np)%fx-rxx*rr1
            ftot(np)%fy=ftot(np)%fy-ryy*rr1
            ftot(np)%fz=ftot(np)%fz-rzz*rr1
        ENDIF
        !!!!!! end UFENE(r) force

        !!!!!! ULJ(r) force

        DO i=1,np
            DO j=1,np
                IF(i/=j)THEN
                    rxx=x(i)-x(j)
                    ryy=y(i)-y(j)
                    rzz=z(i)-z(j)
                    rr=SQRT(rxx**2+ryy**2+rzz**2)

                    IF(rr<=2**(1.0/6)*rp0)THEN
                        rr0=(rp0/rr)**6
                        rr1=dble(48)*e0*(rr0**2)/(rr**2)
                        E=E+4.0*e0*rr0**12
                        ftot(i)%fx=ftot(i)%fx+rxx*rr1
                        ftot(i)%fy=ftot(i)%fy+ryy*rr1
                        ftot(i)%fz=ftot(i)%fz+rzz*rr1
                    ENDIF
                ENDIF
            ENDDO
        ENDDO

        !!!!!! end ULJ(r) force

        call boundaryconditionp()
        do i=1,np
            ftot(i)%fx=ftot(i)%fx+fb(i)%fx
            ftot(i)%fy=ftot(i)%fy+fb(i)%fy
            ftot(i)%fz=ftot(i)%fx+fb(i)%fz
        end do
        E=E+Eb
        return
    end subroutine
end module

program md
    use initializer
    use mpc
    use force
    implicit none
    integer i,j,s,s1,nn
    real(8) oldE
    type(f)::f1(np),f2(np),f3(np)
    idum=-2000
    open(filep,file="polymer.txt")
    open(filef,file="fluid.txt")
    open(filee,file="engery.txt")
    open(filev,file="fv.txt")
    open(filefv,file="velocity.txt",status="old")
    call initializerfp()
    call initialvelocity(vf,nf,m)
    call initialvelocity(vf0,nf,m)
    call initialvelocity(vp,np,mp)
    flag=0
    do i=1,nf0
        read(filefv,*)nn,ff(i)%frx,ff(i)%fry,ff(i)%vm,ff(i)%fvx,ff(i)%fvy
    end do
    close(filefv)
    do i=1,nf0
        ff(i)%frx=(ff(i)%frx-7.0)*len
        ff(i)%fry=(ff(i)%fry)*len
    end do
    do i=1,np
        write(filep,*) rp(i)%rx/len,rp(i)%ry/len,rp(i)%rz/len,nz(i)
    end do
    do i=1,nf
        !write(filef,*) rf(i)%rx,rf(i)%ry,rf(i)%rz
        write(filev,*) vf(i)
    end do
    call xyzout(0)
    do i=1,np
        f1(i)%fx=0.0
        f1(i)%fy=0.0
        f1(i)%fz=0.0
    end do
    call totalforce()
    do i=1,np
        f2(i)%fx=ftot(i)%fx
        f2(i)%fy=ftot(i)%fy
        f2(i)%fz=ftot(i)%fz
    end do
    oldE=E
    do s=1,step
        do s1=1,hn
            do i=1,np
                rp(i)%rx=rp(i)%rx+vp(i)%vx*hmd+((4*f2(i)%fx-f1(i)%fx)/mp)*(hmd**2/6)
                rp(i)%ry=rp(i)%ry+vp(i)%vy*hmd+((4*f2(i)%fy-f1(i)%fy)/mp)*(hmd**2/6)
                rp(i)%rz=rp(i)%rz+vp(i)%vz*hmd+((4*f2(i)%fz-f1(i)%fz)/mp)*(hmd**2/6)
            end do
            call boundaryconditionprz()
            do i=1,nf
                oldry(i)=rf(i)%ry
            end do
            do i=1,nf
                rf(i)%rx=rf(i)%rx+vf(i)%vx*hmd
                rf(i)%ry=rf(i)%ry+vf(i)%vy*hmd
                rf(i)%rz=rf(i)%rz+vf(i)%vz*hmd
            end do
            if(flag==0) then
                call boundaryconditionf()
            else
                call boundaryconditionff()
            end if
            call totalforce()
            do i=1,np
                f3(i)%fx=ftot(i)%fx
                f3(i)%fy=ftot(i)%fy
                f3(i)%fz=ftot(i)%fz
            end do
            do i=1,np
                vp(i)%vx=vp(i)%vx+((2*f3(i)%fx+5*f2(i)%fx-f1(i)%fx)/mp)*(hmd/6)
                vp(i)%vy=vp(i)%vy+((2*f3(i)%fy+5*f2(i)%fy-f1(i)%fy)/mp)*(hmd/6)
                vp(i)%vz=vp(i)%vz+((2*f3(i)%fz+5*f2(i)%fz-f1(i)%fz)/mp)*(hmd/6)
            end do
            call ntv(vp,np,mp)
            do i=1,np
                f1(i)%fx=f2(i)%fx
                f1(i)%fy=f2(i)%fy
                f1(i)%fz=f2(i)%fz
            end do
            do i=1,np
                f2(i)%fx=f3(i)%fx
                f2(i)%fy=f3(i)%fy
                f2(i)%fz=f3(i)%fz
            end do
            do i=1,np
                E=E+0.5*mp*(vp(i)%vx**2+vp(i)%vy**2+vp(i)%vz**2)
            end do
            if(flag==0) then
                if(s>=100000) then
                    flag=1
                    do i=1,nf
                        call fluidvelocity(rf(i),vff(i))
                        vf(i)%vx=vf(i)%vx+vff(i)%vx
                        vf(i)%vy=vf(i)%vy+vff(i)%vy
                        vf(i)%vz=vf(i)%vz+vff(i)%vz
                        write(filev,*) i,vff(i)
                    end do
                end if
            end if
            oldE=E
        end do
        call mpcmethod()
        !if(mod(s,500000)==0) then
        !write(*,"('step='I9)") s
        !end if
        if(mod(s,1000)==0) then
            write(filee,*) oldE/e0
        end if
        if((mod(s,5000)==0)) then
            write(filep,"('step='I9)") s
            do i=1,np
                write(filep,*)  rp(i)%rx/len,rp(i)%ry/len,rp(i)%rz/len,nz(i)
            end do
        end if
        if((mod(s,10000)==0)) then
            !write(filef,"('step='I12)") s
            !do i=1,nf
            write(filef,*) rf(1000)%rx,rf(1000)%ry,rf(1000)%rz,vf(1000)
            !write(filev,*) vf(i)
            !end do
        end if
        if(mod(s,100000)==0) then
            call xyzout(s)
        end if
    end do
    close(filep)
    close(filef)
    close(filee)
    close(filev)
    write(*,*) "program over!"
    stop
end

subroutine ntv(vs,n,mm)
    use initializer
    implicit none
    integer i,n
    real sumvx,sumvy,sumvz,scalx,scaly,scalz,temp,mm
    type(v)::vs(n)
    temp=kb*t/mm
    sumvx=0.0
    sumvy=0.0
    sumvz=0.0
    do i=1,n
        sumvx=sumvx+vs(i)%vx**2
        sumvy=sumvy+vs(i)%vy**2
        sumvz=sumvz+vs(i)%vz**2
    end do
    scalx=sqrt(1.0+2.5D-3*(temp*n/sumvx-1.0))
    scaly=sqrt(1.0+2.5D-3*(temp*n/sumvy-1.0))
    scalz=sqrt(1.0+2.5D-3*(temp*n/sumvz-1.0))
    do i=1,n
        vs(i)%vx=scalx*vs(i)%vx
        vs(i)%vy=scalx*vs(i)%vy
        vs(i)%vz=scalx*vs(i)%vz
    end do
    return
end subroutine

subroutine xyzout(n)
    use initializer
    implicit none
    integer num,n,fileid,i
    character t1,t2
    character(len=20) filename
    fileid=20
    t1='P'
    t2='F'
    num=nf+np
    call integertostring(n,filename)
    open(fileid,file=filename)
    write(fileid,"(I5)") num
    write(fileid,"(a6)")'Atoms'
    do i=1,np
        write(fileid,'(a4,F12.5,F12.5,F12.5)')t1,rp(i)%rx/len,rp(i)%ry/len,rp(i)%rz/len
    end do
    do i=1,nf
        write(fileid,'(a4,F12.5,F12.5,F12.5)')t2,rf(i)%rx/len,rf(i)%ry/len,rf(i)%rz/len
    end do
    close(fileid)
    return
end subroutine

subroutine integertostring(n,s)
    implicit none
    integer n,i,j,k,temp
    character(len=20) s
    temp=n
    i=1
    do while(temp/10/=0)
        i=i+1
        temp=temp/10
    end do
    temp=n
    do k=i,1,-1
        j=mod(temp,10)
        temp=temp/10
        s(k:k)=achar(j+48)
    end do
    s((i+1):(i+4))=".xyz"
    s=trim(s)
    return
end subroutine

FUNCTION rand_num(idum)
    INTEGER idum,ia,im,iq,ir,ntab,ndiv
    REAL rand_num,am,eps,rnmx
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
    rand_num=min(am*iy,rnmx)
END FUNCTION rand_num
