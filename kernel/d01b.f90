
subroutine d01bcf(itype,aa,bb,cc,dd,npnts,weight,abscis,ifail)
! mark 8 release. nag copyright 1979.
! mark 9c revised. ier-370 (jun 1982).
! mark 11.5(f77) revised. (sept 1985.)
! mark 13 revised. use of mark 12 x02 functions (apr 1988).
! mark 14a revised. ier-677 (dec 1989).
! mark 14b revised. ier-840 (mar 1990).
! subroutine for the determination of gaussian quadrature rules
! **************************************************************

! input parameters

! itype  integer which specifies the rule type chosen
! weight w(x)          interval     restrictions
! 0            1                    a,b          b.gt.a
! 1    (b-x)**c*(x-a)**d            a,b     b.gt.a,c,d.gt.-1
! 2   abs(x-0.5*(a+b))**c           a,b     c.gt.-1,b.gt.a
! 3  abs(x-a)**c*exp(-b*x)          a,inf   c.gt.-1,b.gt.0
! 3  abs(x-a)**c*exp(-b*x)       -inf,a     c.gt.-1,b.lt.0
! 4 abs(x-a)**c*exp(-b*(x-a)**2) -inf,inf   c.gt.-1,b.gt.0
! 5  abs(x-a)**c/abs(x+b)**d      a,inf a.gt.-b,c.gt.-1,d.gt.c+1
! 5  abs(x-a)**c/abs(x+b)**d     -inf,a a.lt.-b,c.gt.-1,d.gt.c+1
! abs(itype) must be less than 6. if itype is given less than
! zero then the adjusted weights are calculated. if npnts is
! odd and itype equals -2 or -4 and c is not zero, there may be
! problems.

! aa     real parameter used to specify rule type.  see itype.

! bb     real parameter used to specify rule type.  see itype.

! cc     real parameter used to specify rule type.  see itype.

! dd     real parameter used to specify rule type.  see itype.

! ifail  nag failure parameter.  see nag documentation.

! npnts  integer that determines dimension of weight and abscis

! output parameters

! weight  real array of dimension npnts which contains
! rule weights

! abscis  real array of dimension npnts which contains
! rule abscissae

! ifail integer nag failure parameter
! ifail=0 for normal exit
! ifail=1 for failure in nag routine f02avf
! ifail=2 for parameter npnts or itype out of range
! ifail=3 for parameter aa or bb or cc or dd out of
! allowed range
! ifail=4 for overflow in calculation of weights
! ifail=5 for underflow in calculation of weights
! ifail=6 for itype=-2 or -4, npnts odd, c not zero

! **************************************************************
! .. parameters ..
    character(6) ::       srname
    parameter         (srname='d01bcf')
! .. scalar arguments ..
    double precision ::  aa, bb, cc, dd
    integer ::           ifail, itype, npnts
! .. array arguments ..
    double precision ::  abscis(npnts), weight(npnts)
! .. local scalars ..
    double precision ::  a, abspnc, b, bn, c, cn, cno, d, facn, fn, four, &
    gamma, gammab, gammb, half, one, pna, pnb, pnc, &
    ponorm, psqrd, realmx, small, sqrtcn, store, &
    twnapb, two, wtsum, y, zero
    integer ::           ierror, isub, j, mitype, n, nbug, nfac, nhalf
! .. local arrays ..
    character(1) ::       p01rec(1)
! .. external functions ..
    double precision ::  s14aaf, x02ajf, x02alf
    integer ::           p01abf
    external          s14aaf, x02ajf, x02alf, p01abf
! .. external subroutines ..
    external          f02avf
! .. intrinsic functions ..
    intrinsic         abs, mod, log, exp, dble, sqrt, int
! .. data statements ..
    data              zero, one, two, four/0.0d0, 1.0d0, 2.0d0, 4.0d0/
    data              half/0.5d0/
! .. executable statements ..

! initialisation and parameter checking

    small = x02ajf()
    if (npnts <= 0) go to 780
    do 20 j = 1, npnts
        abscis(j) = zero
        weight(j) = zero
    20 ENDDO
    mitype = abs(itype) + 1
    if (mitype > 6) go to 780
    a = aa
    b = bb
    c = cc
    d = dd
    go to (40,60,100,120,140,160) mitype
    40 c = zero
    d = zero
    60 if (c <= -one .OR. d <= -one) go to 800
    if (b <= a) go to 800
    ponorm = (half*(b-a))**(c+d+one)
    if (itype < 0) ponorm = ponorm/(half*(b-a))**(c+d)
    80 ierror = 1
    gamma = s14aaf(c+one,ierror)
    if (ierror > 0) go to 800
    ierror = 1
    gammb = s14aaf(d+one,ierror)
    if (ierror > 0) go to 800
    ierror = 1
    gammab = s14aaf(c+d+two,ierror)
    if (ierror > 0) go to 800
    ponorm = ponorm*two**(c+d+one)*gamma*gammb/gammab
    abscis(1) = (d-c)/(c+d+two)
    go to 180
    100 if (c <= -one .OR. b <= a) go to 800
    ponorm = two*(half*(b-a))**(c+one)/(c+one)
    if (itype < 0) ponorm = ponorm/(half*(b-a))**c
    go to 180
    120 if (c <= -one .OR. b == zero) go to 800
    ierror = 1
    ponorm = s14aaf(c+one,ierror)*exp(-b*a)/abs(b)**(c+one)
    if (itype < 0) ponorm = ponorm/exp(-b*a)*abs(b)**c
    if (ierror > 0) go to 800
    abscis(1) = c + one
    go to 180
    140 if (c <= -one .OR. b <= zero) go to 800
    ierror = 1
    ponorm = s14aaf((c+one)/two,ierror)/b**((c+one)/two)
    if (itype < 0) ponorm = ponorm*b**(c/two)
    if (ierror > 0) go to 800
    go to 180
    160 if (a+b == zero) go to 800
    if (c <= -one .OR. d <= c+one) go to 800
    d = d - c - two
    ponorm = one/(two**(c+d+one))/(abs(a+b)**(d+one))
    if (itype < 0) ponorm = ponorm*(two**(c+d+two))*(abs(a+b) &
    **(d+two))
    go to 80

! compute diagonal and off-diagonal of symmetric tri-diagonal
! matrix which has abscissae as eigenvalues

    180 if (npnts == 1) go to 320
    do 300 n = 2, npnts
        fn = n - 1
        go to (200,200,220,240,260,200) mitype
        200 twnapb = fn + fn + c + d
        abscis(n) = (d+c)*(d-c)/(twnapb*(twnapb+two))
        cn = four*(fn+c)*(fn+d)*fn/(twnapb**2*(twnapb+one))
        if (n > 2) cn = cn*((c+d+fn)/(twnapb-one))
        go to 280
        220 abscis(n) = zero
        cn = (fn+c*mod(fn,two))**2/((fn+fn+c)**2-one)
        go to 280
        240 abscis(n) = c + fn + fn + one
        cn = fn*(c+fn)
        go to 280
        260 abscis(n) = zero
        cn = (fn+c*mod(fn,two))/two
        280 weight(n) = sqrt(cn)
    300 ENDDO

! use nag routine to find eigenvalues which are abscissae

    320 ierror = 1
    call f02avf(npnts,x02ajf(),abscis,weight,ierror)
    if (ierror > 0) go to 760

! loop to determine weights
! evaluate each orthonormal polynomial of degree
! less than npnts at abscis(j) and sum squares of
! results to determine weight(j)
    ierror = 0
    realmx = x02alf()
    do 700 j = 1, npnts
    
    ! initialise inner loop and scale weight(j) and abscis(j)
    ! divide exponential terms into factors that don't underflow
    
        weight(j) = zero
        y = abscis(j)
        pna = zero
        cno = zero
        nfac = 0
        pnb = one/sqrt(ponorm)
        go to (340,340,360,400,420,440) mitype
        340 abscis(j) = y*(half*(b-a)) + (half*(a+b))
        if (itype > 0) go to 460
        pnb = pnb*(one-y)**(c*half)*(one+y)**(d*half)
        go to 460
        360 abscis(j) = y*(half*(b-a)) + (half*(a+b))
        if (itype > 0 .OR. c == zero) go to 460
        if (y == zero .AND. c > zero) go to 660
        if (c > zero) go to 380
        if (ponorm >= one) go to 380
        if (abs(y) <= (one/(realmx*ponorm))**(-one/c)) go to 680
        380 pnb = pnb*abs(y)**(c*half)
        go to 460
        400 abscis(j) = y/b + a
        if (itype > 0) go to 460
        pnb = pnb*y**(c*half)
        nfac = int(y/log(half*realmx)) + 1
        facn = exp(-half*y/dble(nfac))
        go to 460
        420 abscis(j) = y/sqrt(b) + a
        if (itype > 0) go to 460
        nfac = int(y*y/log(half*realmx)) + 1
        facn = exp(-half*y*y/dble(nfac))
        if (c == zero) go to 460
        if (y == zero .AND. c > zero) go to 660
        if (y == zero .AND. c < zero) go to 680
        pnb = pnb*abs(y)**(c*half)
        go to 460
        440 abscis(j) = two*(a+b)/(y+one) - b
        if (itype > 0) go to 460
        pnb = pnb*(one-y)**(c*half)*(one+y)**(half*(d+two))
        460 wtsum = pnb*pnb
        if (npnts == 1) go to 640
    
    ! loop to evaluate orthonormal polynomials using three
    ! term recurrence relation.
    
        do 620 n = 2, npnts
            fn = n - 1
            go to (480,480,500,520,540,480) mitype
            480 twnapb = fn + fn + c + d
            bn = (d-c)/twnapb
            if (n > 2) bn = bn*(c+d)/(twnapb-two)
            cn = four*fn*(c+fn)*(d+fn)/(twnapb**2*(twnapb+one))
            if (n > 2) cn = cn*((c+d+fn)/(twnapb-one))
            go to 560
            500 bn = zero
            cn = (fn+c*mod(fn,two))**2/((fn+fn+c)**2-one)
            go to 560
            520 bn = c + fn + fn - one
            cn = fn*(fn+c)
            go to 560
            540 bn = zero
            cn = (fn+c*mod(fn,two))/two
            560 sqrtcn = sqrt(cn)
            pnc = ((y-bn)*pnb-cno*pna)/sqrtcn
            cno = sqrtcn
            abspnc = abs(pnc)
            if (abspnc <= one) go to 580
            if (abspnc <= realmx/abspnc) go to 580
            if (itype > 0) go to 680
            if (nfac <= 0) go to 680
            pnb = pnb*facn
            pnc = pnc*facn
            wtsum = wtsum*facn*facn
            nfac = nfac - 1
            580 psqrd = pnc*pnc
            if (wtsum <= realmx-psqrd) go to 600
            if (itype > 0) go to 680
            if (nfac <= 0) go to 680
            pnb = pnb*facn
            pnc = pnc*facn
            wtsum = wtsum*facn*facn
            psqrd = psqrd*facn*facn
            nfac = nfac - 1
            600 wtsum = wtsum + psqrd
            pna = pnb
            pnb = pnc
        620 ENDDO
    
    ! end loop for polynomial evaluation
    
    ! richard brankin - nag, oxford - 26th july 1989
    ! replaced the following line ....
    
    ! 640    if (nfac.gt.0) wtsum = wtsum*facn**(2*nfac)
    
    ! so as not to get needless underflow to zero when powering up facn
    ! for 0.0 < facn << 1.0. the error was brought to light in a vax
    ! double precision implementation when a user tried to compute modified
    ! laguerre weights (itype = -3) for more than 25 abscissae (n > 25).
    ! as a result, before the assignment in the above line
    ! wtsum = o(1.0e+38), facn = o(1.0e-10), nfac = 2
    ! wtsum was assigned a value of 0.0 since o(1.0e-10)**4 underflows
    ! although wtsum should have been assigned o(1.0e+2). this correction
    ! also applies for other values of itype.
    
        640 if (nfac > 0) then
            do 650 nbug = 1, 2*nfac
                wtsum = wtsum*facn
            650 ENDDO
        end if
    
    ! end of correction
    
        if (wtsum == zero) go to 660
        weight(j) = one/wtsum
        go to 700
        660 ierror = 4
        weight(j) = realmx
        go to 700
        680 ierror = 5
    700 ENDDO

! end loop for weights

! reverse rational or laguerre points

    if ((mitype /= 6 .OR. a+b < zero) &
     .AND. (mitype /= 4 .OR. b > zero)) go to 740
    nhalf = npnts/2
    if (nhalf <= 1) go to 740
    do 720 j = 1, nhalf
        isub = npnts + 1 - j
        store = abscis(j)
        abscis(j) = abscis(isub)
        abscis(isub) = store
        store = weight(j)
        weight(j) = weight(isub)
        weight(isub) = store
    720 ENDDO

! assignment of ifail parameter

    740 if ((itype == -2 .OR. itype == -4) .AND. mod(npnts,2) &
     == 1 .AND. c /= zero) ierror = 6
    go to 820
    760 ierror = 1
    go to 820
    780 ierror = 2
    go to 820
    800 ierror = 3
    820 ifail = p01abf(ifail,ierror,srname,0,p01rec)
    return
 end subroutine d01bcf
 subroutine f02avf(n,acheps,d,e,ifail)
! mark 2 release. nag copyright 1972
! mark 3 revised.
! mark 4 revised.
! mark 4.5 revised
! mark 9 revised. ier-326 (sep 1981).
! mark 11.5(f77) revised. (sept 1985.)

! tql1
! this subroutine finds the eigenvalues of a tridiagonal
! matrix,
! t, given with its diagonal elements in the array d(n) and
! its subdiagonal elements in the last n - 1 stores of the
! array e(n), using ql transformations. the eigenvalues are
! overwritten on the diagonal elements in the array d in
! ascending order. the subroutine will fail if all
! eigenvalues take more than 30*n iterations.
! 1st april 1972

! .. parameters ..
    character(6) ::       srname
    parameter         (srname='f02avf')
! .. scalar arguments ..
    double precision ::  acheps
    integer ::           ifail, n
! .. array arguments ..
    double precision ::  d(n), e(n)
! .. local scalars ..
    double precision ::  b, c, f, g, h, p, r, s
    integer ::           i, i1, ii, isave, j, l, m, m1
! .. local arrays ..
    character(1) ::       p01rec(1)
! .. external functions ..
    integer ::           p01abf
    external          p01abf
! .. intrinsic functions ..
    intrinsic         abs, sqrt
! .. executable statements ..
    isave = ifail
    if (n == 1) go to 40
    do 20 i = 2, n
        e(i-1) = e(i)
    20 ENDDO
    40 e(n) = 0.0d0
    b = 0.0d0
    f = 0.0d0
    j = 30*n
    do 340 l = 1, n
        h = acheps*(abs(d(l))+abs(e(l)))
        if (b < h) b = h
    ! look for small sub diagonal element
        do 60 m = l, n
            if (abs(e(m)) <= b) go to 80
        60 ENDDO
        80 if (m == l) go to 260
        100 if (j <= 0) go to 360
        j = j - 1
    ! form shift
        g = d(l)
        h = d(l+1) - g
        if (abs(h) >= abs(e(l))) go to 120
        p = h*0.5d0/e(l)
        r = sqrt(p*p+1.0d0)
        h = p + r
        if (p < 0.0d0) h = p - r
        d(l) = e(l)/h
        go to 140
        120 p = 2.0d0*e(l)/h
        r = sqrt(p*p+1.0d0)
        d(l) = e(l)*p/(1.0d0+r)
        140 h = g - d(l)
        i1 = l + 1
        if (i1 > n) go to 180
        do 160 i = i1, n
            d(i) = d(i) - h
        160 ENDDO
        180 f = f + h
    ! ql transformation
        p = d(m)
        c = 1.0d0
        s = 0.0d0
        m1 = m - 1
        do 240 ii = l, m1
            i = m1 - ii + l
            g = c*e(i)
            h = c*p
            if (abs(p) < abs(e(i))) go to 200
            c = e(i)/p
            r = sqrt(c*c+1.0d0)
            e(i+1) = s*p*r
            s = c/r
            c = 1.0d0/r
            go to 220
            200 c = p/e(i)
            r = sqrt(c*c+1.0d0)
            e(i+1) = s*e(i)*r
            s = 1.0d0/r
            c = c/r
            220 p = c*d(i) - s*g
            d(i+1) = h + s*(c*g+s*d(i))
        240 ENDDO
        e(l) = s*p
        d(l) = c*p
        if (abs(e(l)) > b) go to 100
        260 p = d(l) + f
    ! order eigenvalue
        if (l == 1) go to 300
        do 280 ii = 2, l
            i = l - ii + 2
            if (p >= d(i-1)) go to 320
            d(i) = d(i-1)
        280 ENDDO
        300 i = 1
        320 d(i) = p
    340 ENDDO
    ifail = 0
    return
    360 ifail = p01abf(isave,1,srname,0,p01rec)
    return
 	end subroutine f02avf
 	
 	double precision function s14aaf(x,ifail)
! mark 7 release. nag copyright 1978.
! mark 7c revised ier-184 (may 1979)
! mark 11.5(f77) revised. (sept 1985.)
! gamma function

! **************************************************************

! to extract the correct code for a particular machine-range,
! activate the statements contained in comments beginning  cdd ,
! where  dd  is the approximate number of significant decimal
! digits represented by the machine
! delete the illegal dummy statements of the form
! * expansion (nnnn) *

! also insert appropriate data statements to define constants
! which depend on the range of numbers represented by the
! machine, rather than the precision (suitable statements for
! some machines are contained in comments beginning crd where
! d is a digit which simply distinguishes a group of machines).
! delete the illegal dummy data statements with values written
! *value*

! **************************************************************

! .. parameters ..
    character(6) ::                      srname
    parameter                        (srname='s14aaf')
! .. scalar arguments ..
    double precision ::                 x
    integer ::                          ifail
! .. local scalars ..
    double precision ::                 g, gbig, t, xbig, xminv, xsmall, y
    integer ::                          i, m
! .. local arrays ..
    character(1) ::                      p01rec(1)
! .. external functions ..
    integer ::                          p01abf
    external                         p01abf
! .. intrinsic functions ..
    intrinsic                        abs, sign, dble
! .. data statements ..
! 8   data xsmall/1.0d-8/
! 9   data xsmall/3.0d-9/
! 2   data xsmall/1.0d-12/
! 5   data xsmall/3.0d-15/
    data xsmall/1.0d-17/
! 9    data xsmall/1.7d-18/

    data xbig,gbig,xminv/ 1.70d+2,4.3d+304,2.23d-308 /
! xbig = largest x such that  gamma(x) .lt. maxreal
! and  1.0/gamma(x+1.0) .gt. minreal
! (rounded down to an integer)
! gbig = gamma(xbig)
! xminv = max(1.0/maxreal,minreal)  (rounded up)
! for ieee single precision
! 0   data xbig,gbig,xminv /33.0e0,2.6e+35,1.2e-38/
! for ibm 360/370 and similar machines
! 1   data xbig,gbig,xminv /57.0d0,7.1d+74,1.4d-76/
! for dec-10, honeywell, univac 1100 (s.p.)
! 2   data xbig,gbig,xminv /34.0d0,8.7d+36,5.9d-39/
! for icl 1900
! 3   data xbig,gbig,xminv /58.0d0,4.0d+76,1.8d-77/
! for cdc 7600/cyber
! 4   data xbig,gbig,xminv /164.0d0,2.0d+291,3.2d-294/
! for univac 1100 (d.p.)
! 5   data xbig,gbig,xminv /171.0d0,7.3d+306,1.2d-308/
! for ieee double precision
! 7   data xbig,gbig,xminv /170.0d0,4.3d+304,2.3d-308/
! .. executable statements ..

! error 1 and 2 test
    t = abs(x)
    if (t > xbig) go to 160
! small range test
    if (t <= xsmall) go to 140
! main range reduction
    m = x
    if (x < 0.0d0) go to 80
    t = x - dble(m)
    m = m - 1
    g = 1.0d0
    if (m) 20, 120, 40
    20 g = g/x
    go to 120
    40 do 60 i = 1, m
        g = (x-dble(i))*g
    60 ENDDO
    go to 120
    80 t = x - dble(m-1)
! error 4 test
    if (t == 1.0d0) go to 220
    m = 1 - m
    g = x
    do 100 i = 1, m
        g = (dble(i)+x)*g
    100 ENDDO
    g = 1.0d0/g
    120 t = 2.0d0*t - 1.0d0

! * expansion (0026) *

! expansion (0026) evaluated as y(t)  --precision 08e.09
! 8   y = ((((((((((((+1.88278283d-6)*t-5.48272091d-6)
! 8  *    *t+1.03144033d-5)*t-3.13088821d-5)*t+1.01593694d-4)
! 8  *    *t-2.98340924d-4)*t+9.15547391d-4)*t-2.42216251d-3)
! 8  *    *t+9.04037536d-3)*t-1.34119055d-2)*t+1.03703361d-1)
! 8  *    *t+1.61692007d-2)*t + 8.86226925d-1

! expansion (0026) evaluated as y(t)  --precision 09e.10
! 9   y = (((((((((((((-6.463247484d-7)*t+1.882782826d-6)
! 9  *    *t-3.382165478d-6)*t+1.031440334d-5)*t-3.393457634d-5)
! 9  *    *t+1.015936944d-4)*t-2.967655076d-4)*t+9.155473906d-4)
! 9  *    *t-2.422622002d-3)*t+9.040375355d-3)*t-1.341184808d-2)
! 9  *    *t+1.037033609d-1)*t+1.616919866d-2)*t + 8.862269255d-1

! expansion (0026) evaluated as y(t)  --precision 12e.13
! 2   y = (((((((((((((((-7.613347676160d-8)*t+2.218377726362d-7)
! 2  *    *t-3.608242105549d-7)*t+1.106350622249d-6)
! 2  *    *t-3.810416284805d-6)*t+1.138199762073d-5)
! 2  *    *t-3.360744031186d-5)*t+1.008657892262d-4)
! 2  *    *t-2.968993359366d-4)*t+9.158021574033d-4)
! 2  *    *t-2.422593898516d-3)*t+9.040332894085d-3)
! 2  *    *t-1.341185067782d-2)*t+1.037033635205d-1)
! 2  *    *t+1.616919872669d-2)*t + 8.862269254520d-1

! expansion (0026) evaluated as y(t)  --precision 15e.16
! 5   y = (((((((((((((((-1.243191705600000d-10
! 5  *    *t+3.622882508800000d-10)*t-4.030909644800000d-10)
! 5  *    *t+1.265236705280000d-9)*t-5.419466096640000d-9)
! 5  *    *t+1.613133578240000d-8)*t-4.620920340480000d-8)
! 5  *    *t+1.387603440435200d-7)*t-4.179652784537600d-7)
! 5  *    *t+1.253148247777280d-6)*t-3.754930502328320d-6)
! 5  *    *t+1.125234962812416d-5)*t-3.363759801664768d-5)
! 5  *    *t+1.009281733953869d-4)*t-2.968901194293069d-4)
! 5  *    *t+9.157859942174304d-4)*t-2.422595384546340d-3
! 5   y = ((((y*t+9.040334940477911d-3)*t-1.341185057058971d-2)
! 5  *    *t+1.037033634220705d-1)*t+1.616919872444243d-2)*t +
! 5  *     8.862269254527580d-1

! expansion (0026) evaluated as y(t)  --precision 17e.18
    y = (((((((((((((((-1.46381209600000000d-11 &
    *t+4.26560716800000000d-11)*t-4.01499750400000000d-11) &
    *t+1.27679856640000000d-10)*t-6.13513953280000000d-10) &
    *t+1.82243164160000000d-9)*t-5.11961333760000000d-9) &
    *t+1.53835215257600000d-8)*t-4.64774927155200000d-8) &
    *t+1.39383522590720000d-7)*t-4.17808776355840000d-7) &
    *t+1.25281466396672000d-6)*t-3.75499034136576000d-6) &
    *t+1.12524642975590400d-5)*t-3.36375833240268800d-5) &
    *t+1.00928148823365120d-4)*t-2.96890121633200000d-4
    y = ((((((y*t+9.15785997288933120d-4)*t-2.42259538436268176d-3) &
    *t+9.04033494028101968d-3)*t-1.34118505705967765d-2) &
    *t+1.03703363422075456d-1)*t+1.61691987244425092d-2)*t + &
    8.86226925452758013d-1

! expansion (0026) evaluated as y(t)  --precision 19e.20
! 9   y = (((((((((((((((+6.7108864000000000000d-13
! 9  *    *t-1.6777216000000000000d-12)*t+6.7108864000000000000d-13)
! 9  *    *t-4.1523609600000000000d-12)*t+2.4998051840000000000d-11)
! 9  *    *t-6.8985815040000000000d-11)*t+1.8595971072000000000d-10)
! 9  *    *t-5.6763875328000000000d-10)*t+1.7255563264000000000d-9)
! 9  *    *t-5.1663077376000000000d-9)*t+1.5481318277120000000d-8)
! 9  *    *t-4.6445740523520000000d-8)*t+1.3931958370304000000d-7)
! 9  *    *t-4.1782339907584000000d-7)*t+1.2528422549504000000d-6)
! 9  *    *t-3.7549858152857600000d-6)*t+1.1252456510305280000d-5
! 9   y = (((((((((y*t-3.3637584239226880000d-5)
! 9  *    *t+1.0092815021080832000d-4)
! 9  *    *t-2.9689012151880000000d-4)*t+9.1578599714350784000d-4)
! 9  *    *t-2.4225953843706897600d-3)*t+9.0403349402888779200d-3)
! 9  *    *t-1.3411850570596516480d-2)*t+1.0370336342207529018d-1)
! 9  *    *t+1.6169198724442506740d-2)*t + 8.8622692545275801366d-1

    s14aaf = y*g
    ifail = 0
    go to 240

! error 3 test
    140 if (t < xminv) go to 200
    s14aaf = 1.0d0/x
    ifail = 0
    go to 240

! error exits
    160 if (x < 0.0d0) go to 180
    ifail = p01abf(ifail,1,srname,0,p01rec)
    s14aaf = gbig
    go to 240

    180 ifail = p01abf(ifail,2,srname,0,p01rec)
    s14aaf = 0.0d0
    go to 240

    200 ifail = p01abf(ifail,3,srname,0,p01rec)
    t = x
    if (x == 0.0d0) t = 1.0d0
    s14aaf = sign(1.0d0/xminv,t)
    go to 240

    220 ifail = p01abf(ifail,4,srname,0,p01rec)
    s14aaf = gbig

    240 return
 	end function s14aaf

    double precision function x02ajf()
! mark 12 release. nag copyright 1986.

! returns  (1/2)*b**(1-p)  if rounds is .true.
! returns  b**(1-p)  otherwise

    double precision :: x02con
    data x02con /1.11022302462516d-16 /
! .. executable statements ..
    x02ajf = x02con
    return
    end function x02ajf
    
    double precision function x02alf()
! mark 12 release. nag copyright 1986.

! returns  (1 - b**(-p)) * b**emax  (the largest positive model
! number)

    double precision :: x02con
    data x02con /1.79769313486231d+308 /
! .. executable statements ..
    x02alf = x02con
    return
    end function x02alf
    integer function p01abf(ifail,ierror,srname,nrec,rec)
! mark 11.5(f77) release. nag copyright 1986.
! mark 13 revised. ier-621 (apr 1988).
! mark 13b revised. ier-668 (aug 1988).

! p01abf is the error-handling routine for the nag library.

! p01abf either returns the value of ierror through the routine
! name (soft failure), or terminates execution of the program
! (hard failure). diagnostic messages may be output.

! if ierror = 0 (successful exit from the calling routine),
! the value 0 is returned through the routine name, and no
! message is output

! if ierror is non-zero (abnormal exit from the calling routine),
! the action taken depends on the value of ifail.

! ifail =  1: soft failure, silent exit (i.e. no messages are
! output)
! ifail = -1: soft failure, noisy exit (i.e. messages are output)
! ifail =-13: soft failure, noisy exit but standard messages from
! p01abf are suppressed
! ifail =  0: hard failure, noisy exit

! for compatibility with certain routines included before mark 12
! p01abf also allows an alternative specification of ifail in which
! it is regarded as a decimal integer with least significant digits
! cba. then

! a = 0: hard failure  a = 1: soft failure
! b = 0: silent exit   b = 1: noisy exit

! except that hard failure now always implies a noisy exit.

! s.hammarling, m.p.hooper and j.j.du croz, nag central office.

! .. scalar arguments ..
    integer ::                 ierror, ifail, nrec
    character*(*)           srname
! .. array arguments ..
    character*(*)           rec(*)
! .. local scalars ..
    integer ::                 i, nerr
    character(72) ::            mess
! .. external subroutines ..
    external                p01abz, x04aaf, x04baf
! .. intrinsic functions ..
    intrinsic               abs, mod
! .. executable statements ..
    if (ierror /= 0) then
    ! abnormal exit from calling routine
        if (ifail == -1 .OR. ifail == 0 .OR. ifail == -13 .OR. &
        (ifail > 0 .AND. mod(ifail/10,10) /= 0)) then
        ! noisy exit
            call x04aaf(0,nerr)
            do 20 i = 1, nrec
                call x04baf(nerr,rec(i))
            20 ENDDO
            if (ifail /= -13) then
                write (mess,fmt=99999) srname, ierror
                call x04baf(nerr,mess)
                if (abs(mod(ifail,10)) /= 1) then
                ! hard failure
                    call x04baf(nerr, &
                    ' ** nag hard failure - execution terminated' &
                    )
                    call p01abz
                else
                ! soft failure
                    call x04baf(nerr, &
                    ' ** nag soft failure - control returned')
                end if
            end if
        end if
    end if
    p01abf = ierror
    return

    99999 format (' ** abnormal exit from nag library routine ',a,': ifail', &
    ' =',i6)
    end function p01abf
    subroutine p01abz
! mark 11.5(f77) release. nag copyright 1986.

! terminates execution when a hard failure occurs.

! ******************** implementation note ********************
! the following stop statement may be replaced by a call to an
! implementation-dependent routine to display a message and/or
! to abort the program.
! *************************************************************
! .. executable statements ..
    stop
    end subroutine p01abz
    subroutine x04aaf(i,nerr)
! mark 7 release. nag copyright 1978
! mark 7c revised ier-190 (may 1979)
! mark 11.5(f77) revised. (sept 1985.)
! mark 14 revised. ier-829 (dec 1989).
! if i = 0, sets nerr to current error message unit number
! (stored in nerr1).
! if i = 1, changes current error message unit number to
! value specified by nerr.

! .. scalar arguments ..
    integer ::           i, nerr
! .. local scalars ..
    integer ::           nerr1
! .. save statement ..
    save              nerr1
! .. data statements ..
    data              nerr1/0/
! .. executable statements ..
    if (i == 0) nerr = nerr1
    if (i == 1) nerr1 = nerr
    return
    end subroutine x04aaf
    subroutine x04baf(nout,rec)
! mark 11.5(f77) release. nag copyright 1986.

! x04baf writes the contents of rec to the unit defined by nout.

! trailing blanks are not output, except that if rec is entirely
! blank, a single blank character is output.
! if nout.lt.0, i.e. if nout is not a valid fortran unit identifier,
! then no output occurs.

! .. scalar arguments ..
    integer ::           nout
    character*(*)     rec
! .. local scalars ..
    integer ::           i
! .. intrinsic functions ..
    intrinsic         len
! .. executable statements ..
    if (nout >= 0) then
    ! remove trailing blanks
        do 20 i = len(rec), 2, -1
            if (rec(i:i) /= ' ') go to 40
        20 ENDDO
    ! write record to external file
        40 write (nout,fmt=99999) rec(1:i)
    end if
    return

    99999 format (a)
    end subroutine x04baf
