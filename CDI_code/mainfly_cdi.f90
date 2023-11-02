

    program main

    integer:: ntc,nolayer,ns    ! number of time channel / layers / parameter per point / points
    integer:: i,j, k  ! index
    integer:: num_sig

    Real *8 t_st1,t_ed1,xr1, hr1,rt1, rr1
    integer::nturn1, nturn11

    real *8 para1,para2,para3,para4

    real *8 pi,htt
    real *8 delta
    real *8 rho_st,rho_ed,rholog
    real *8,allocatable:: point1set(:), point2set(:), point3set(:), point4set(:),ht(:)

    real *8,allocatable::hz1obs(:),rho(:)      ! dBz/dt HCP
    real *8,allocatable::hz2obs(:)      ! dBz/dt VCX
    real *8,allocatable::time(:)

    real *8,allocatable:: rho_true(:), hh_true(:), mu_true(:)
    real*8 :: rho_n = 0.98, rho_tao = 1.0, rho_c = 0.25
    ! real*8 :: rho_n = 0.01, rho_tao = 0.01, rho_c = 0.3
    Data pi/3.1415926D0/
    Common /complex_rho/rho_n, rho_tao, rho_c 

    open(1, file='parameter_in.txt', status='old')

    read(1,*)para1,para2,para3,para4,t_st1,t_ed1,xr1, hr1,rt1, rr1

    ntc = floor(para1)
    ns = floor(para2)
    nturn1 = floor(para3)
    nturn11 = floor(para4)


    !ntc = 56          ! 时间道
    nolayer = 1
    num_sig = 256*2

    ns = 1

    allocate(rho(num_sig))
    allocate(hz2obs(ntc))
    allocate(hz1obs(ntc),time(ntc))


    allocate(rho_true(nolayer), hh_true(nolayer), mu_true(nolayer))

    allocate(point1set(ns*2), point2set(ns*2), point3set(ns*2), point4set(ns*2),ht(ns))


    !******************************************************************
    Open (6, File='log1.dat', Status='unknown')
    Open (7, File='rho1.dat', Status='unknown')
    Open (10, File='point1set.txt', Status='old')
    Open (11, File='point2set.txt', Status='old')
    Open (12, File='point3set.txt', Status='old')
    Open (13, File='point4set.txt', Status='old')


    Open (15, File='flag.dat', Status='unknown')



    do i=1,ns*2
        read(10,*)point1set(i)
        read(11,*)point2set(i)
        read(12,*)point3set(i)
        read(13,*)point4set(i)
    end do


    rho_st = 1.d-05
    rho_ed = 1.d05

    delta = (dlog10(rho_ed)-dlog10(rho_st))/num_sig

    do i=1,num_sig
        rholog=dlog10(rho_st)+i*delta
        rho(i)=10**rholog
    end do

    Write(7,100)rho


    do k=1,ns
        htt = hr1

        do i=1,num_sig  ! 测点 从0m出开始
            do j=1,nolayer
                rho_true(j) = rho(i)
                mu_true(j) = pi*4d-07
                hh_true(j) = 10
            end do

            call forwardprocess(rho_true,hh_true,mu_true,hz1obs,nolayer,ntc,time,k,&
            & point1set,point2set,point3set,point4set,htt,ns,t_st1,t_ed1,xr1,hr1,rt1,rr1,nturn1,nturn11) ! observe value

            Write(6,10)hz1obs
        end do
    end do


    Write(15,*)1.0

10  Format(512E14.6)
100 Format(1E14.6)

    end program


    Subroutine forwardprocess(rho, hh, mu, hz1, nlayer,nt,t,ns_id,p1,p2,p3,p4,ht,ns1,t_st,t_ed,xr,hr,rt,rr,nturn,nturn1)
    integer::nlayer, nt, npls, ns_id
    Real *8 t_st,t_ed,delta_t
    Real *8 t(nt), hz1(nt),tlog(nt)

    Real *8 s2, ss2
    Real *8 r, rplus, tt1, tt2, tt3, tt4, pi
    Real *8 rho(nlayer), hh(nlayer), frq(67), mu(nlayer)

    Real *8 xr, yr, ht, hr, zplus, zminus

    real *8 slope1, slope2, slope3, imax, imin

    real *8 p1(ns1*2), p2(ns1*2), p3(ns1*2), p4(ns1*2)

    real *8  rt, rr
    integer::nturn, nturn1

    real *8 R_matrix(3,3), moment(3,1)

    Complex *16 func(9, 67), hs_mat(3,3),fun_calibrate(3,1), fun_calibrate_z(67)


    Data pi/3.1415926D0/
    Data ngau/10000/

    Common /para/r

    Common /funn/frq, nfrq, func, fun_calibrate_z
    Data nfrq, (frq(i), i=1, 67)/67,&
        & 0.10000000d-02,0.14677993d-02,0.21544347d-02,0.31622777d-02,&
        & 0.46415888d-02,0.68129207d-02,0.10000000d-01,0.14677993d-01,&
        & 0.21544347d-01,0.31622777d-01,0.46415888d-01,0.68129207d-01,&
        & 0.10000000d+00,0.14677993d+00,0.21544347d+00,0.31622777d+00,&
        & 0.46415888d+00,0.68129207d+00,0.10000000d+01,0.14677993d+01,&
        & 0.21544347d+01,0.31622777d+01,0.46415888d+01,0.68129207d+01,&
        & 0.10000000d+02,0.14677993d+02,0.21544347d+02,0.31622777d+02,&
        & 0.46415888d+02,0.68129207d+02,0.10000000d+03,0.14677993d+03,&
        & 0.21544347d+03,0.31622777d+03,0.46415888d+03,0.68129207d+03,&
        & 0.10000000d+04,0.14677993d+04,0.21544347d+04,0.31622777d+04,&
        & 0.46415888d+04,0.68129207d+04,0.10000000d+05,0.14677993d+05,&
        & 0.21544347d+05,0.31622777d+05,0.46415888d+05,0.68129207d+05,&
        & 0.10000000d+06,0.14677993d+06,0.21544347d+06,0.31622777d+06,&
        & 0.46415888d+06,0.68129207d+06,0.10000000d+07,0.14677993d+07,&
        & 0.21544347d+07,0.31622777d+07,0.46415888d+07,0.68129207d+07,&
        & 0.10000000d+08,0.14677993d+08,0.21544347d+08,0.31622777d+08,&
        & 0.46415888d+08,0.68129207d+08,0.1000000d+09/
    !**************************************************************
    Call filter
    !**************************************************************

    !hr = 0.01
    !xr = 0.58
    yr = 0.

    zplus = ht - hr
    zminus = ht + hr
    r = dsqrt(xr*xr+yr*yr)
    rplus = dsqrt(r*r+zplus*zplus)

    !**************************************************************
    ! 发射线圈
    !rt = 0.5          ! 线圈的半径
    !nturn = 3.         ! 线圈的匝数

    !rr = 0.25       ! 半径
    !nturn1 = 20     ! 接收线圈的匝数


    npls = 1
    moment(1,1) = 0
    moment(2,1) = 0
    moment(3,1) = 1
    delta_t = (dlog10(t_ed)-dlog10(t_st))/(nt-1)


    do i=1,nt
        tlog(i)=dlog10(t_st)+(i-1)*delta_t
        t(i)=10**tlog(i)
    end do


    ic = 3

    If (ic==2) Then   ! 单极性方波
        tt1 = 2e-3
        tt2 = t_ed-tt1
    End If

    If (ic==3) Then   ! 单极性方波

        tt1 = p2(1+(ns_id-1)*2)-p1(1+(ns_id-1)*2)
        tt2 = p3(1+(ns_id-1)*2)-p2(1+(ns_id-1)*2)
        tt3 = p4(1+(ns_id-1)*2)-p3(1+(ns_id-1)*2)
        tt4 = t_ed-tt1-tt2-tt3

        slope1 = (p2(2+(ns_id-1)*2)-p1(2+(ns_id-1)*2))/tt1
        slope2 = (p2(2+(ns_id-1)*2)-p3(2+(ns_id-1)*2))/tt2
        slope3 = (p3(2+(ns_id-1)*2)-p4(2+(ns_id-1)*2))/tt3

    End If

    If (ic==4) Then   ! 单极性方波
        tt1 = 0.5e-03
        tt2 = 1e-03
        tt3 = t_ed-2.0d0*tt1-tt2
    End If

    Do i = 1, nfrq
        Call forward(rho, hh, frq(i), hs_mat, zplus, zminus, nlayer) ! mu0*H = B
        fun_calibrate = matmul(hs_mat,moment)   ! 3x1
        fun_calibrate_z(i) = fun_calibrate(3,1)
    End Do

    !********************* impulse and step wave ************************
    If (ic==0 .Or. ic==1) Then
        ik = 0
        Do i = 1,nt

            Call frt(rho, hh, t(i), hz1(i), 2, zplus, zminus, ic, ik, nlayer)

            ik = ik + 1

            ! Transformation of cylindrical coordinate system
            hz1(i) = (4*rt*rt*nturn)*hz1(i)*(pi*rr**2*nturn1)

        End Do

    !******************* for rectangle wave(multi-pulse)***********************
    Else If (ic==3) Then
        Do i = 1, nt

            ss2 = 0.D0

            ik = 0

            Do ip = 1, npls
                If (t(i)>(ip-1)*(tt1+tt2+tt3+tt4) .And. t(i)<ip*(tt1+tt2+tt3+tt4)) kpls = ip
            End Do

            Do ip = 1, kpls - 1

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2+tt3+tt4), s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + s2*slope1

                ik = ik + 1

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2+tt3+tt4)-tt1, s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + (-1)*s2*(slope1+slope2)

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2+tt3+tt4)-tt1-tt2, s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + (-1)*s2*(slope3-slope2)

                Call frt(rho, hh, t(i)-(ip-1)*(tt1+tt2+tt3+tt4)-tt1-tt2-tt3, s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + s2*slope3
            End Do

            !****************** contribution from resting pulses**************************
            If (t(i)>(kpls-1)*(tt1+tt2+tt3+tt4))Then

                Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2+tt3+tt4), s2, 2, zplus, zminus, 1, ik, nlayer)

                ss2 = ss2 + s2*slope1

                ik = ik + 1

                If (t(i)>(kpls-1)*(tt1+tt2+tt3+tt4)+tt1)Then

                    Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2+tt3+tt4)-tt1, s2, 2, zplus, zminus, 1, ik, nlayer)

                    ss2 = ss2 + (-1)*s2*(slope1+slope2)


                    If (t(i)>(kpls-1)*(tt1+tt2+tt3+tt4)+tt1+tt2)Then

                        Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2+tt3+tt4)-tt1-tt2, s2, 2, zplus, zminus, 1, ik, nlayer)

                        ss2 = ss2 + (-1)*s2*(slope3-slope2)

                        If (t(i)>(kpls-1)*(tt1+tt2+tt3+tt4)+tt1+tt2+tt3)Then

                            Call frt(rho, hh, t(i)-(kpls-1)*(tt1+tt2+tt3+tt4)-tt1-tt2-tt3, s2, 2, zplus, zminus, 1, ik, nlayer)

                            ss2 = ss2 + s2*slope3

                        End If
                    End If
                End If
            End If

            ! dB/dt for trapezoid
            hz1(i) = 4*rt*rt*nturn*ss2*(pi*rr**2*nturn1)
        End Do

    End If


    End Subroutine forwardprocess

    !        gauleg(0.d0,t(i),x,w,ngau)
    Subroutine gauleg(x1, x2, x, w, n)
    Implicit Real *8(A-H, O-Z)
    Real *8 x1, x2, xm, xl, x(n), w(n)
    Parameter (eps=3.D-14)
    m = (n+1)/2
    xm = 0.5D0*(x2+x1)
    xl = 0.5D0*(x2-x1)
    Do i = 1, m
        z = cos(3.141592654D0*(i-0.25D0)/(n+0.5D0))
1       Continue
        p1 = 1.D0
        p2 = 0.D0
        Do j = 1, n
            p3 = p2
            p2 = p1
            p1 = ((2.D0*j-1.D0)*z*p2-(j-1.D0)*p3)/j
        End Do
        pp = n*(z*p1-p2)/(z*z-1.D0)
        z1 = z
        z = z1 - p1/pp
        If (abs(z-z1)>eps) Goto 1
        x(i) = xm - xl*z
        x(n+1-i) = xm + xl*z
        w(i) = 2.D0*xl/((1.D0-z*z)*pp*pp)
        w(n+1-i) = w(i)
    End Do
    Return
    End Subroutine gauleg

    !        spl(nfrq, nc, frq, funr0, f, funr1)
    Subroutine spl(nx, n2, x, fx, x2, fx2)
    Real *8 x(nx), fx(nx), c(3, nx), x2(n2), fx2(n2), xint, xl1, xl2
    Call splin1(nx, fx, c)
    xl1 = dlog10(x(1))
    xl2 = dlog10(x(nx))
    Do ix = 1, n2
        xint = dlog10(x2(ix))
        Call splin2(nx, xint, xl1, xl2, c, fx2(ix))
    End Do
    End Subroutine spl

    Subroutine splin1(n, y, c)
    Real *8 y(n), c(3, n), p
    n1 = n - 1
    Do i = 2, n1
        c(1, i) = y(i+1) - 2.*y(i) + y(i-1)
    End Do
    c(2, 1) = 0.D0
    c(3, 1) = 0.D0
    Do i = 2, n1
        p = 4. + c(2, i-1)
        c(2, i) = -1.D0/p
        c(3, i) = (c(1,i)-c(3,i-1))/p
    End Do
    c(1, n) = 0.D0
    Do ii = 2, n1
        i = n + 1 - ii
        c(1, i) = c(2, i)*c(1, i+1) + c(3, i)
    End Do
    c(1, 1) = 0.
    Do i = 1, n1
        c(2, i) = y(i+1) - y(i) - c(1, i+1) + c(1, i)
        c(3, i) = y(i) - c(1, i)
    End Do
    c(3, n) = y(n)
    Return
    End Subroutine splin1

    Subroutine splin2(n, xint, x1, x2, c, yint)
    Real *8 c(3, n), xint, x1, x2, yint, h, u, p, q
    h = (x2-x1)/dble(float(n-1))
    If (xint<x1) Goto 10
    If (xint>=x2) Goto 20
    u = (xint-x1)/h
    i = 1 + int(u)
    p = u - i + 1
    q = 1.D0 - p
    yint = c(1, i)*q**3 + c(1, i+1)*p**3 + c(2, i)*p + c(3, i)
    Return
10  p = (xint-x1)/h
    yint = c(2, 1)*p + c(3, 1)
    Return
20  p = (xint-x2)/h
    yint = c(2, n-1)*p + c(3, n)
    Return
    End Subroutine splin2

Subroutine frt(rho1, hh1, t, ft, item, zplus, zminus, ic, ik, nlayer)
  Complex *16 fun, iomega, func(9, 67)
  Real *8 t, ft, zplus, zminus, pi, q, rho1(nlayer), hh1(nlayer)
  Complex *16 fun_calibrate_z(67)
    Real *8 frq(67), funr0(67), funi0(67)
  Real *8 f(160), omega(160), funr1(160), funi1(160), h(200)
  Common /funn/frq, nfrq, func, fun_calibrate_z

  Data pi, q/3.141592654D0, 1.258925412D0/
  Data ncnull, nc, ndec, (h(i), i=1, 160)/80,160,10,&
     & 2.59511139938829d-13,3.66568771323555d-13,5.17792876616242d-13,&
     & 7.31400730405791d-13,1.03313281156235d-12,1.45933600088387d-12,&
     & 2.06137146234699d-12,2.91175733962418d-12,4.11297804457870d-12,&
     & 5.80971771117984d-12,8.20647323099742d-12,1.15919058389365d-11,&
     & 1.63740746547780d-11,2.31288803930431d-11,3.26705938902288d-11,&
     & 4.61481520721098d-11,6.51864545047052d-11,9.20775899532545d-11,&
     & 1.30064200980219d-10,1.83718747396255d-10,2.59512512377884d-10,&
     & 3.66566596154242d-10,5.17796324027279d-10,7.31395266627501d-10,&
     & 1.03314147106736d-09,1.45932227649333d-09,2.06139321404013d-09,&
     & 2.91172286551380d-09,4.11303268236158d-09,5.80963111612975d-09,&
     & 8.20661047490285d-09,1.15916883220051d-08,1.63744193958818d-08,&
     & 2.31283340152144d-08,3.26714598407299d-08,4.61467796330556d-08,&
     & 6.84744728867720d-08,5.46574677490374d-08,1.13319898777493d-07,&
     & 2.16529974157527d-07,2.88629942214140d-07,3.42872728051125d-07,&
     & 4.79119488706262d-07,7.42089418889752d-07,1.07736520535271d-06,&
     & 1.46383231306575d-06,2.01727682134668d-06,2.89058197617431d-06,&
     & 4.15237808867022d-06,5.84448989361742d-06,8.18029430348419d-06,&
     & 1.15420854481494d-05,1.63897017145322d-05,2.31769096113890d-05,&
     & 3.26872676331330d-05,4.60786866701851d-05,6.51827321351636d-05,&
     & 9.20862589540037d-05,1.30169142615951d-04,1.83587481111627d-04,&
     & 2.59595544393723d-04,3.66324383719323d-04,5.18210697462501d-04,&
     & 7.30729969562531d-04,1.03385239132389d-03,1.45738764044730d-03,&
     & 2.06298256402732d-03,2.90606401578959d-03,4.11467957883740d-03,&
     & 5.79034253321120d-03,8.20005721235220d-03,1.15193892333104d-02,&
     & 1.63039398900789d-02,2.28256810984487d-02,3.22248555163692d-02,&
     & 4.47865101670011d-02,6.27330674874545d-02,8.57058672847471d-02,&
     & 1.17418179407605d-01,1.53632645832305d-01,1.97718111895102d-01,&
     & 2.28849924263247d-01,2.40310905012422d-01,1.65409071929404d-01,&
     & 2.84709685167114d-03,-2.88015846269687d-01,-3.69097391853225d-01,&
     & -2.50109865922601d-02,5.71811109500426d-01,-3.92261390212769d-01,&
     & 7.63282774297327d-02,5.16233692927851d-02,-6.48015160576432d-02,&
     & 4.89045522502552d-02,-3.26934307794750d-02,2.10542570949745d-02,&
     & -1.33862848934736d-02,8.47098801479259d-03,-5.35134515919751d-03,&
     & 3.37814023806349d-03,-2.13157364002470d-03,1.34506352474558d-03,&
     & -8.48929743771803d-04,5.35521822356713d-04,-3.37744799986382d-04,&
     & 2.13268792633204d-04,-1.34629969723156d-04,8.47737416679279d-05,&
     & -5.34940635827096d-05,3.3904416298191d-05,-2.13315638358794d-05,&
     & 1.33440911625019d-05,-8.51629073825634d-06,5.44362672273211d-06,&
     & -3.32112278417896d-06,2.07147190852386d-06,-1.42009412555511d-06,&
     & 8.78247754998004d-07,-4.5566280473703d-07,3.38598103040009d-07,&
     & -2.87407830772251d-07,1.07866150545699d-07,-2.4724024185358d-08,&
     & 5.35535110396030d-08,-3.3789981131378d-08,2.13200367531820d-08,&
     & -1.34520337740075d-08,8.48765950790546d-09,-5.35535110396018d-09,&
     & 3.37899811131383d-09,-2.13200367531819d-09,1.34520337740075d-09,&
     & -8.48765950790576d-10,5.35535110396015d-10,-3.37899811131382d-10,&
     & 2.13200367531811d-10,-1.34520337740079d-10,8.48765950790572d-11,&
     & -5.35535110396034d-11,3.37899811131381d-11,-2.13200367531818d-11,&
     & 1.34520337740074d-11,-8.48765950790571d-12,5.35535110396031d-12,&
     & -3.37899811131379d-12,2.13200367531817d-12,-1.34520337740073d-12,&
     & 8.48765950790567d-13,-5.35535110396029d-13,3.37899811131377d-13,&
     & -2.13200367531816d-13,1.34520337740078d-13,-8.48765950790596d-14,&
     & 5.35535110396007d-14,-3.37899811131377d-14,2.13200367531816d-14,&
     & -1.34520337740083d-14,8.4876550790558d-15,-5.35535110396025d-15,&
     & 3.37899811131389d-15/

  Do nn = 1, nc
    n = -nc + ncnull + nn
    omega(nn) = q**(-(n-1))/t
    f(nn) = omega(nn)/(2.D0*pi)
  End Do

  Do i = 1, nfrq
      funr0(i) = dreal(fun_calibrate_z(i))
      funi0(i) = dimag(fun_calibrate_z(i))
  End Do
  Call spl(nfrq, nc, frq, funr0, f, funr1)  !interpolation
  Call spl(nfrq, nc, frq, funi0, f, funi1)
  ft = 0.D0
  Do nn = 1, nc
    If (ic==0) & ! impulse
      Then
      iomega = (1.D0, 0.D0)
    Else If (ic==1) Then
      iomega = 1.D0/((0.,-1.D0)*omega(nn)) ! -1/iw
    End If
    fun = dcmplx(funr1(nn), funi1(nn))*iomega
    ita = max0(1, nn-nc+1)
    ite = min0(1, nn)
    Do it = ita, ite
      itn = nc - nn + it
      ft = ft + dimag(fun)*dsqrt(omega(nn))*h(itn) ! primary field st
    End Do
  End Do
  ft = -ft*dsqrt(2.D0/pi/t)
  Return
End Subroutine frt

Subroutine forward(rho2, hh2, f, fun, zplus, zminus, nlayer)
  Complex *16 t3, t5, t6, hf, fun(3,3)
  Real *8 pi, f, zplus, zminus
  Real *8 r, rho2(nlayer), hh2(nlayer)
  Common /para/r
  pi = 3.1415926D0

    fun(1,1) = 0
    fun(1,2) = 0
    fun(1,3) = 0
    fun(2,1) = 0
    fun(2,2) = 0
    fun(2,3) = 0
    fun(3,1) = 0
    fun(3,2) = 0
    fun(3,3) = -t5(rho2, hh2, f, zplus, zminus, nlayer)*4.d-7*pi

  Return
End Subroutine forward

    Complex *16 Function t5(rho5, hh5, f, zplus, zminus, nlayer)
  Complex *16 s, s1, b
    Real *8 r, rho5(nlayer), hh5(nlayer)
    Real *8 h0, h1, f, zplus, zminus, u, fac, expc

  Common /para/r
  Common /hankel/nc, ncnull, h0(100), h1(100)
  fac = 0.1D0*dlog(10.D0)
  s = (0.D0, 0.D0)
  Do nn = 1, nc
    nu = nn
    mn = nc - nn + 1
    nnn = ncnull - nc + nu
    u = expc(-(nnn-1)*fac)/r
    ! z+h, z-h, 其中z是接收线圈的z坐标，z=-ht，-h是发射线圈的坐标, h=hr
    ! 前文 zminus = ht + hr = (h + h), 此处使用-zminus取反
    ! zplus = ht - hr = -z - h = - (z+h)，所以要使用 -zplus 对应论文中的 z+h
    ! print*, "r=",r
    s1 = u * bessj0(u*r) * &
    ( (b(rho5, hh5,f,u, nlayer)-u)/(b(rho5, hh5,f,u, nlayer)+u)*expc(-u*zminus) )
    ! ( expc(-u*abs(zplus)) + (b(rho5, hh5,f,u, nlayer)-u)/(b(rho5, hh5,f,u, nlayer)+u)*expc(-u*zminus) )
    s = s + s1*h1(mn)
  End Do
  t5 = s/r
  Return
End Function t5

    FUNCTION bessj0(x)
        REAL bessj0
        real*8 x
        REAL ax,xx,z
        DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,&
        r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,y
        SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,&
        s1,s2,s3,s4,s5,s6
        DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,&
        -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5&
        /-.1562499995d-1,.1430488765d-3,-.6911147651d-5,&
        .7621095161d-6,-.934945152d-7/
        DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,&
        651619640.7d0,-11214424.18d0,77392.33017d0,&
        -184.9052456d0/,s1,s2,s3,s4,s5,s6/57568490411.d0,&
        1029532985.d0,9494680.718d0,59272.64853d0,&
        267.8532712d0,1.d0/
        if(abs(x)<8.) then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/&
        (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
        else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*&
        (p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*&
        (q3+y*(q4+y*q5)))))
        endif
    END FUNCTION bessj0

! output: BE
! input: f - frequency
! input: u - k
! input: mu - 磁导率
! nlayer -  no. of layers
Complex *16 Function b(rho4, hh4, f, u, nlayer)
  Complex *16 alpha, s1, s2
  Real *8 f, u, pi,gam
  Real *8 r, rho4(nlayer), hh4(nlayer)
  real *8 mu0
    real *8 :: rho_n, rho_tao, rho_c
    Complex *16 cplxRho
  Common /para/r
    Common /complex_rho/rho_n, rho_tao, rho_c 

  pi = 3.1415926D0
  mu0 = pi*4d-07

  ! print*, rhoNs
  b = cdsqrt(u*u+(0.D0,1.D0)*mu0*2*pi*f/ &
  cplxRho(rho4(nlayer), rho_n, rho_tao, rho_c, 2*pi*f))  ! alpha nlayer omega = 2*pi*f
  ! b = cdsqrt(u*u+(0.D0,1.D0)*mu0*2*pi*f/ rho4(nlayer))

  If (nlayer==1) Return
  Do i = nlayer-1, 1, -1

      alpha = cdsqrt(u*u+(0.D0,1.D0)*mu0*2*pi*f/ &
      cplxRho(rho4(i), rho_n, rho_tao, rho_c, 2*pi*f))

      s1 = (0.D0, 0.D0)
      If (dreal(2.D0*alpha*hh4(i))<400.D0) s1 = cdexp(-2.D0*alpha*hh4(i))
      s2 = (1.D0-s1)/(1.D0+s1)  ! tanh
      b = alpha*(b+alpha*s2)/(alpha+b*s2)
  End Do
End Function b

complex *16 Function cplxRho(rho_0, rho_n, rho_tao, rho_c, omega)
! implicit none
real *8 :: rho_0, rho_n, rho_tao, rho_c, omega

! print*, rho_n
cplxRho = dcmplx(rho_0, 0) * &
(1 - rho_n * (1 - 1 / (1 + dcmplx(0, omega*rho_tao) ** rho_c)))
! print*, 'cplxRho', cplxRho
return
End Function cplxRho

Real *8 Function expc(x)
  Real *8 x, x1
  x1 = x
  If (dabs(x1)>650.D0) x1 = dsign(650.D0, x1)
  expc = dexp(x1)
  Return
End Function expc

Subroutine filter
  Common /hankel/nc, ncnull, h0(100), h1(100)
  Real *8 h0, h1

  Data(h0(i),i=1,48)/&
  & 2.89878288d-07,3.64935144d-07,4.59426126d-07,5.78383226d-07,&
  & 7.28141338d-07,9.16675639d-07,1.15402625d-06,1.45283298d-06,&
  & 1.82900834d-06,2.30258511d-06,2.89878286d-06,3.64935148d-06,&
  & 4.59426119d-06,5.78383236d-06,7.28141322d-06,9.16675664d-06,&
  & 1.15402621d-05,1.45283305d-05,1.82900824d-05,2.30258527d-05,&
  & 2.89878259d-05,3.64935186d-05,4.59426051d-05,5.78383329d-05,&
  & 7.28141144d-05,9.16675882d-05,1.15402573d-04,1.45283354d-04,&
  & 1.82900694d-04,2.30258630d-04,2.89877891d-04,3.64935362d-04,&
  & 4.59424960d-04,5.78383437d-04,7.28137738d-04,9.16674828d-04,&
  & 1.15401453d-03,1.45282561d-03,1.82896826d-03,2.30254535d-03,&
  & 2.89863979d-03,3.64916703d-03,4.59373308d-03,5.78303238d-03,&
  & 7.27941497d-03,9.16340705d-03,1.15325691d-02,1.45145832d-02/

      Data(h0(i),i=49,100)/&
     & 1.82601199d-02,2.29701042d-02,2.88702619d-02,3.62691810d-02,&
     & 4.54794031d-02,5.69408192d-02,7.09873072d-02,8.80995426d-02,&
     & 1.08223889d-01,1.31250483d-01,1.55055715d-01,1.76371506d-01,&
     & 1.85627738d-01,1.69778044d-01,1.03405245d-01,-3.02583233d-02,&
     & -2.27574393d-01,-3.62173217d-01,-2.05500446d-01,3.37394873d-01,&
     & 3.17689897d-01,-5.13762160d-01,3.09130264d-01,-1.26757592d-01,&
     & 4.61967890d-02,-1.80968674d-02,8.35426050d-03,-4.47368304d-03,&
     & 2.61974783d-03,-1.60171357d-03,9.97717882d-04,-6.26275815d-04,&
     & 3.94338818d-04,-2.48606354d-04,1.56808604d-04,-9.89266288d-05,&
     & 6.24152398d-05,-3.93805393d-05,2.48472358d-05,-1.56774945d-05,&
     & 9.89181741d-06,-6.24131160d-06,3.93800058d-06,-2.48471018d-06,&
     & 1.56774609d-06,-9.89180896d-07,6.24130948d-07,-3.93800005d-07,&
     & 2.48471005d-07,-1.56774605d-07,9.89180888d-08,-6.24130946d-08/

      Data(h1(i),i=1,48)/&
     & 1.84909557d-13,2.85321327d-13,4.64471808d-13,7.16694771d-13,&
     & 1.16670043d-12,1.80025587d-12,2.93061898d-12,4.52203829d-12,&
     & 7.36138206d-12,1.13588466d-11,1.84909557d-11,2.85321327d-11,&
     & 4.64471808d-11,7.16694771d-11,1.16670043d-10,1.80025587d-10,&
     & 2.93061898d-10,4.52203829d-10,7.36138206d-10,1.13588466d-09,&
     & 1.84909557d-09,2.85321326d-09,4.64471806d-09,7.16694765d-09,&
     & 1.16670042d-08,1.80025583d-08,2.93061889d-08,4.52203807d-08,&
     & 7.36138149d-08,1.13588452d-07,1.84909521d-07,2.85321237d-07,&
     & 4.64471580d-07,7.16694198d-07,1.16669899d-06,1.80025226d-06,&
     & 2.93060990d-06,4.52201549d-06,7.36132477d-06,1.13587027d-05,&
     & 1.84905942d-05,2.85312247d-05,4.64449000d-05,7.16637480d-05,&
     & 1.16655653d-04,1.79989440d-04,2.92971106d-04,4.51975783d-04/

      Data(h1(i),i=49,100)/&
     & 7.35565435d-04,1.13444615d-03,1.84548306d-03,2.84414257d-03,&
     & 4.62194743d-03,7.10980590d-03,1.15236911d-02,1.76434485d-02,&
     & 2.84076233d-02,4.29770596d-02,6.80332569d-02,9.97845929d-02,&
     & 1.51070544d-01,2.03540581d-01,2.71235377d-01,2.76073871d-01,&
     & 2.16691977d-01,-7.83723737d-02,-3.40675627d-01,-3.60693673d-01,&
     & 5.13024526d-01,-5.94724729d-02,-1.95117123d-01,1.99235600d-01,&
     & -1.38521553d-01,8.79320859d-02,-5.50697146d-02,3.45637848d-02,&
     & -2.17527180d-02,1.37100291d-02,-8.64656417d-03,5.45462758d-03,&
     & -3.44138864d-03,2.17130686d-03,-1.36998628d-03,8.64398952d-04,&
     & -5.45397874d-04,3.44122545d-04,-2.17126585d-04,1.36997597d-04,&
     & -8.64396364d-05,5.45397224d-05,-3.44122382d-05,2.17126544d-05,&
     & -1.36997587d-05,8.64396338d-06,-5.45397218d-06,3.44122380d-06,&
     & -2.17126543d-06,1.36997587d-06,-8.64396337d-07,5.45397218d-07/

  nc = 100
  ncnull = 60
  Return
End Subroutine filter
