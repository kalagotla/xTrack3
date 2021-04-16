      SUBROUTINE V3_STREAM(IXP, IYP, ID, CELL, IC, IXY, DTRP, XY, XYZ,
     &       V, CEL1, CCEL1, CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &       NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK)
C
C     CALCULATE THE STREAM LINE AT POSITION IXP IYP IN 2D WINDOW
C
C	Copyright 1990 - 2005, Massachusetts Institute of Technology.
C
      INCLUDE 'Visual3.inc'
C
      REAL    XYZ(3,*),    V(3,*),     DTRP(*),     XY(2,*)
      INTEGER CELL(3,*),   IXY(3,*),   IC(*)
      INTEGER CEL1(4,*),   CEL2(5,*),  CEL3(6,*),   CEL4(8,*)
      INTEGER CCEL1(4,*),  CCEL2(5,*), CCEL3(5,*),  CCEL4(6,*)
      INTEGER NPOLYT(8,*), POLYT(*),   CPOLYT(2,*), BLOCKS(6,*)
      INTEGER CBLOCK(*), N
      REAL    rhop, dp, nu, sigma, delta1
      INTEGER dpInt
      CHARACTER*32 name
      CHARACTER*200 location
C
      REAL    X0(3)
C
      Common /Orkwis/ rhop, dp
      CALL CPU_TIME(start)
      DO N=1,189
        rhop = 819
        dp = 0.281e-6
        KSTREAM = 0
        CALL V3_2DW23D(IXP, IYP, KC, X0, CELL, IC, IXY, DTRP, XY, XYZ)
        IF(KC .EQ. 0) RETURN
C
        X0(1) = 13.3*1e-3
        X0(2) = 0.409091*1e-3
        X0(3) = (0.35*N*1e-3)/100.0

        dpInt = dp*1e9
        write(location,*) '../../output/umData/'
        write(name,'(I0,A,I0,A,I0)') dpInt, '/Z', dpInt, '.', N
        open(unit=77,file=TRIM(ADJUSTL(location))//name,
     &                  form='formatted',status='unknown')
        Print*, 'Particle number=', N
        Print*, 'Spawned at', X0
      	CALL V3_CalcSL(X0, KC, ID, XYZ, V, CEL1, CCEL1, CEL2, CCEL2,
     &               CEL3, CCEL3, CEL4, CCEL4, NPOLYT, POLYT, CPOLYT,
     &               BLOCKS, CBLOCK)
        close(77)
      ENDDO
      CALL CPU_TIME(finish)
      Print*, 'ELAPSED CPU TIME = ', finish-start
 77   RETURN
      END


      SUBROUTINE V3_CalcSL(XOLD, KC, ID, XYZ, V, CEL1, CCEL1, CEL2,
     &                     CCEL2, CEL3, CCEL3, CEL4, CCEL4, NPOLYT, 
     &                     POLYT, CPOLYT, BLOCKS, CBLOCK)
C
C     CALCULATE THE STREAM LINE AT POSITION XOLD (IN CELL KC)
C
C	Copyright 1990 - 2005, Massachusetts Institute of Technology.
C
      PARAMETER (DELMAX = 30.0*0.017453292)
      PARAMETER (MCOUNT = 100, mnod=1e7)
      INCLUDE 'XFtn.inc'
      INCLUDE 'Visual3.inc'
C
      REAL    XYZ(3,*)
      INTEGER CEL1(4,*),   CEL2(5,*),  CEL3(6,*),   CEL4(8,*)
      INTEGER CCEL1(4,*),  CCEL2(5,*), CCEL3(5,*),  CCEL4(6,*)
      INTEGER NPOLYT(8,*), POLYT(*),   CPOLYT(2,*), BLOCKS(6,*)
      INTEGER CBLOCK(*)
      REAL XOLD(3), VpOLD(3), Vp0(3), V(3,*)
C
      REAL  W(8), XSAV(3)
      REAL  UCURL(3), UCURL0(3), UTDIV, UTDIV0
      REAL  Vpa(3), Av(3), Au(3)
      REAL  rho(8), Vf(3,8), ork(3), wis(3), UOLD(3), rhor, Kmc
      REAL Cd, Re, knd, mfp, df, nu, vps, ufs, vrs, T(8), mu, mu0
      INTEGER Di, Lip, Dilip, KCOLD
      LOGICAL MAXOUT, DOUBLE, PASS, PASSI, HALVED
      REAL X(3), DX1(3), DX2(3), DX3(3), DX4(3)
      REAL U0(3), U(3),Vp(3),A(3),A0(3)
      REAL DVp1(3),DVp2(3),DVp3(3),DVp4(3)
	  
C
      DATA EPLUS/0.025/, EMINUS/1e-6/
      DATA MHIT/24/
      Common /Orkwis/ rhop, dp
      common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
C
      mu0 = 1.716e-5
      Kmc = 0.75/(rhop*dp)
      xsav(1) = xold(1)
      xsav(2) = xold(2)
      xsav(3) = xold(3)
      ork(1) = 0
      ork(2) = 0
      ork(3) = 0
      tosl = 0.0
      aold = 0.0
      asav = aold
      dold = 0.0
      dsav = dold
      KC0  = KC
      Dilip = 0
C
      CALL V3_UX(KC, XOLD, U0, XYZ, V, W, CEL1, CCEL1,
     &           CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &           NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &           %val(piblank), UCURL0, UTDIV0, PASSI)
      VpOLD(1)= U0(1)
      VpOLD(2)= U0(2)
      VpOLD(3)= U0(3)
      Vp0(1)= U0(1)
      Vp0(2)= U0(2)
      Vp0(3)= U0(3)
      UOLD(1) = U0(1)
      UOLD(2) = U0(2)
      UOLD(3) = U0(3)
      do j = 1, 3
        do i = 1, 3
          ugrad(i,j) = 0.0
        enddo
      enddo
      TIME  = 0.0
      ATIM  = 0.0
      DTlim = V3_DTLIM(KC, XYZ, V, CEL1, CEL2, CEL3, CEL4,
     &                 NPOLYT, POLYT, BLOCKS, X2, U2)
      IF(U2 .EQ. 0.0) THEN
       WRITE(*,*) 'StreamLine Error - Zero Velocity in Starting Cell!'
       RETURN
      ENDIF
      KCOUNT  = 0
      DT = 0.5*X2/U2
      IF(ID .NE. 1) DT = -DT
      MAXOUT = .FALSE.
      DOUBLE = .TRUE.
      PASS   = .FALSE.
      NHIT   = 0
      LASTKC = KC
      Print*, 'Particle run started'
C
C---- Main loop
C
 2    CALL V3_UX(KC, XOLD, U0, XYZ, V, W, CEL1, CCEL1,
     &           CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &           NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &           %val(piblank), UCURL0, UTDIV0, PASSI)
      IF(KC.EQ.0) GOTO 4
      IF(KC.LT.0) THEN
        XOLD(1) = XSAV(1)
        XOLD(2) = XSAV(2)
        XOLD(3) = XSAV(3)
        AOLD    = ASAV
        DOLD    = DSAV
        TIME    = ATIM
        DT      = 0.5*DT
        KC      = KC0
        DOUBLE  = .FALSE.
        NHIT    = NHIT + 1
        IF(NHIT .GT. MHIT) THEN
            Print*, 'NHIT>MHIT; Terminated'
            GOTO 1
        ENDIF
        GOTO 3
      ENDIF
      IF(KSTREAM .GT. 2) THEN
       DO I = 0, 2
         IF(XOLD(1) .NE. XYZSTREAM(1,KSTREAM-I)   .OR.
     &      XOLD(2) .NE. XYZSTREAM(2,KSTREAM-I)   .OR.
     &      XOLD(3) .NE. XYZSTREAM(3,KSTREAM-I)) GOTO 5
       ENDDO
       TIME = TSTREAM(KSTREAM)
       Print*, 'Terminated from KSTREAM'
       GOTO 1
      ENDIF
 5    KSTREAM = KSTREAM + 1
      CSTREAM(KSTREAM)     = KC
      TSTREAM(KSTREAM)     = TIME
      XYZSTREAM(1,KSTREAM) = XOLD(1)
      XYZSTREAM(2,KSTREAM) = XOLD(2)
      XYZSTREAM(3,KSTREAM) = XOLD(3)
      ASTREAM(KSTREAM)     = AOLD
      DSTREAM(KSTREAM)     = DOLD
      DO I = 1, 8
        FSTREAM(I,KSTREAM) = W(I)
      ENDDO
      DOUBLE = .TRUE.
      IF(KC .NE. LASTKC) THEN
        LASTKC = KC
        NHIT   = 0
      ENDIF
      IF(KSTREAM .EQ. MSTREAM) THEN
       IF(MAXOUT) THEN
       ENDIF
       MAXOUT = .TRUE.
       KS = 2
       DO J = 4, KSTREAM, 2
         KS = KS + 1
         CSTREAM(KS)     = CSTREAM(J)
         TSTREAM(KS)     = TSTREAM(J)
         XYZSTREAM(1,KS) = XYZSTREAM(1,J) 
         XYZSTREAM(2,KS) = XYZSTREAM(2,J) 
         XYZSTREAM(3,KS) = XYZSTREAM(3,J) 
         ASTREAM(KS)     = ASTREAM(J)
         DSTREAM(KS)     = DSTREAM(J)
         DO I = 1, 8
           FSTREAM(I,KS) = FSTREAM(I,J)
         ENDDO
       ENDDO
       KSTREAM = KS
      ENDIF
      KC0 = KC
      DTlim = V3_DTLIM(KC, XYZ, V, CEL1, CEL2, CEL3, CEL4,
     &                 NPOLYT, POLYT, BLOCKS, X2M, U2M)
      IF(U2M .EQ. 0.0) GO TO 1
      DTMAX = 3.0*X2M/U2M
      DTPMAX = X2M/(VpOLD(1)**2+VpOLD(2)**2+VpOLD(3)**2)
      IF(ABS(DT) .GT. DTPMAX) THEN
       DT = SIGN(DTPMAX, DT)
       DOUBLE = .FALSE.
      ENDIF
      IF(ABS(DT) .GT. DTMAX) THEN
       DT = SIGN(DTMAX, DT)
       DOUBLE = .FALSE.
      ENDIF
      IF(ABS(DT) .GT. DTlim) then
        DT = SIGN(DTlim, DT)
        DOUBLE = .FALSE.
      ENDIF
      DT = 1e-9
C
C     CALCULATE A0, INPUT VOLD, K
3     IF(ID .NE. 1) DT = -DT
      Dilip=Dilip+1
      Do N=1,8
        KN = CEL4(N,KC)
        rho(N) = q(1,KN)
        T(N) = q(5,KN)
      ENDDO
      Km = Kmc*(SUM(rho)/SIZE(rho))
      rhor = (SUM(rho)/SIZE(rho))/rhop
      KCOLD = KC
      mu = mu0*(((SUM(T)/SIZE(T))/273.15)**1.5)
     &            *(384.15/(((SUM(T)/SIZE(T))+111)))
      nu = mu/(SUM(rho)/SIZE(rho))
      vps = SQRT(SUM((VpOLD)**2))
      ufs = SQRT(SUM(U0**2))
      vrs = SQRT(SUM((VpOLD-U0)**2))
      Re = (vrs)*dp/nu
      Cd = C_DRAG(Re, rho, mu, T)
      DO N = 1, 3
        Di = U0(N)*1000
        Lip = UOLD(N)*1000
        A0(N)  = Km*Cd*(VpOLD(N)-U0(N))*vrs
        DVp1(N) = DT*A0(N)
        Vp(N)  = VpOLD(N) - 0.5*DVp1(N)
        DX1(N) = DT*VpOLD(N)
        X(N)   = XOLD(N) + 0.5*DX1(N)
      ENDDO
C
C     CHECK FOR ESCAPE
C
      IF(MOD(KCOUNT, MCOUNT) .EQ. 0) THEN
        IF(XFtnGetEvent(display, key_window, ButtonRelease, IE, JE,
     &                  ISTATE) .EQ. 0) THEN
          IF(IE.GE.10.AND.IE.LE.26.AND.JE.GE.67.AND.JE.LE.83) THEN
            KSTREAM = -1
            RETURN
          ENDIF
        ENDIF
      ENDIF
      KCOUNT = KCOUNT + 1
C
C     INTEGRATE RIBBON ANGLE AND TUBE SIZE
C       note: rotation rate is 1/2 the curl
C
      UMAG = SQRT(U0(1)**2 + U0(2)**2 + U0(3)**2)
      IF(UMAG .EQ. 0.0) THEN
        DA1 = 0.0
      ELSE
        DA1 = 0.5*DT*(UCURL0(1)*U0(1)+UCURL0(2)*U0(2)+
     &                    UCURL0(3)*U0(3))/UMAG
      ENDIF
      DD1 = DT*UTDIV0
C
      CALL V3_UX(KC, X, U, XYZ, V, W, CEL1, CCEL1,
     &           CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &           NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &           %val(piblank), UCURL, UTDIV, PASSI)
      IF(REENTER .AND. NHIT.EQ.0 .AND.KC.LT.0.AND.KC.NE.-99999999) THEN
        CALL V3_GetBCell(KC, X, XYZ, %val(pscel),
     &                       %val(psurf), %val(ppsurf),
     &                       %val(pblock), %val(piblank))
        IF(KC .GT. 0)
     &    CALL V3_UX(KC, X, U, XYZ, V, W, CEL1, CCEL1,
     &           CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &           NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &           %val(piblank), UCURL, UTDIV, PASSI)
        PASS = KC .GT. 0
      ENDIF
      PASS = PASS .OR. PASSI
      IF(KC.LE.0) GOTO 44
C
C     CALCULATE A(N)
      Do N=1,8
        KN = CEL4(N,KC)
        rho(N) = q(1,KN)
        T(N) = q(5,KN)
      ENDDO
      Km = Kmc*(SUM(rho)/SIZE(rho))
      rhor = (SUM(rho)/SIZE(rho))/rhop
      mu = mu0*(((SUM(T)/SIZE(T))/273.15)**1.5)
     &            *(384.15/(((SUM(T)/SIZE(T))+111)))
      nu = mu/(SUM(rho)/SIZE(rho))
      vps = SQRT(SUM((Vp)**2))
      ufs = SQRT(SUM(U**2))
      vrs = SQRT(SUM((Vp-U)**2))
      Re = (vrs)*dp/nu
      Cd = C_DRAG(Re, rho, mu, T)
      DO N = 1, 3
        A(N)   = Km*Cd*(Vp(N)-U(N))*vrs
        DX2(N) = DT*Vp(N)
        DVp2(N) = DT*A(N)
        Vp(N)   = VpOLD(N) - 0.5*DVp2(N)
        Di = (Vp(N)*1000)
        Lip = (U(N)*1000)
        X(N)   = XOLD(N) + 0.5*DX2(N)
      ENDDO
C
      UMAG = SQRT(U(1)**2 + U(2)**2 + U(3)**2)
      IF(UMAG .EQ. 0.0) THEN
        DA2 = 0.0
      ELSE
        DA2 = 0.5*DT*(UCURL(1)*U(1)+UCURL(2)*U(2)+UCURL(3)*U(3))/UMAG
      ENDIF
      DD2 = DT*UTDIV
C
      CALL V3_UX(KC, X, U, XYZ, V, W, CEL1, CCEL1,
     &           CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &           NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &           %val(piblank), UCURL, UTDIV, PASSI)
      IF(REENTER .AND. NHIT.EQ.0 .AND.KC.LT.0.AND.KC.NE.-99999999) THEN
        CALL V3_GetBCell(KC, X, XYZ, %val(pscel),
     &                       %val(psurf), %val(ppsurf),
     &                       %val(pblock), %val(piblank))
        IF(KC .GT. 0)
     &    CALL V3_UX(KC, X, U, XYZ, V, W, CEL1, CCEL1,
     &           CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &           NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &           %val(piblank), UCURL, UTDIV, PASSI)
        PASS = KC .GT. 0
      ENDIF
      PASS = PASS .OR. PASSI
      IF(KC.LE.0) GOTO 44
C
C     CALCULATE A
      Do N=1,8
        KN = CEL4(N,KC)
        rho(N) = q(1,KN)
        T(N) = q(5,KN)
      ENDDO
      Km = Kmc*(SUM(rho)/SIZE(rho))
      rhor = (SUM(rho)/SIZE(rho))/rhop
      mu = mu0*(((SUM(T)/SIZE(T))/273.15)**1.5)
     &            *(384.15/(((SUM(T)/SIZE(T))+111)))
      nu = mu/(SUM(rho)/SIZE(rho))
      vps = SQRT(SUM((Vp)**2))
      ufs = SQRT(SUM(U**2))
      vrs = SQRT(SUM((Vp-U)**2))
      Re = (vrs)*dp/nu
      Cd = C_DRAG(Re, rho, mu, T)
      DO N = 1, 3
        A(N)   = Km*Cd*(Vp(N)-U(N))*vrs
        DX3(N) = DT*Vp(N)
        DVp3(N) = DT*A(N)
        Vp(N)   = VpOLD(N) - DVp3(N)
        Di = (Vp(N)*1000)
        Lip = (U(N)*1000)
        X(N)   = XOLD(N) + DX3(N)
      ENDDO
C
      UMAG = SQRT(U(1)**2 + U(2)**2 + U(3)**2)
      IF(UMAG .EQ. 0.0) THEN
        DA3 = 0.0
      ELSE
        DA3 = 0.5*DT*(UCURL(1)*U(1)+UCURL(2)*U(2)+UCURL(3)*U(3))/UMAG
      ENDIF
      DD3 = DT*UTDIV
C
      CALL V3_UX(KC, X, U, XYZ, V, W, CEL1, CCEL1,
     &           CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &           NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &           %val(piblank), UCURL, UTDIV, PASSI)
      IF(REENTER .AND. NHIT.EQ.0 .AND.KC.LT.0.AND.KC.NE.-99999999) THEN
        CALL V3_GetBCell(KC, X, XYZ, %val(pscel),
     &                       %val(psurf), %val(ppsurf),
     &                       %val(pblock), %val(piblank))
        IF(KC .GT. 0)
     &    CALL V3_UX(KC, X, U, XYZ, V, W, CEL1, CCEL1,
     &           CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &           NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &           %val(piblank), UCURL, UTDIV, PASSI)
        PASS = KC .GT. 0
      ENDIF
      PASS = PASS .OR. PASSI
      IF(KC.LE.0) GOTO 44
C
C     CALCULATE A
      Do N=1,8
        KN = CEL4(N,KC)
        rho(N) = q(1,KN)
        T(N) = q(5,KN)
      ENDDO
      Km = Kmc*(SUM(rho)/SIZE(rho))
      rhor = (SUM(rho)/SIZE(rho))/rhop
      mu = mu0*(((SUM(T)/SIZE(T))/273.15)**1.5)
     &            *(384.15/(((SUM(T)/SIZE(T))+111)))
      nu = mu/(SUM(rho)/SIZE(rho))
      vps = SQRT(SUM((Vp)**2))
      ufs = SQRT(SUM(U**2))
      vrs = SQRT(SUM((Vp-U)**2))
      Re = (vrs)*dp/nu
      Cd = C_DRAG(Re, rho, mu, T)
      DO N = 1, 3
        Di = U(N)*1000
        Lip = UOLD(N)*1000
        A(N)   = Km*Cd*(Vp(N)-U(N))*vrs
        DX4(N) = DT*Vp(N)
        DVp4(N) = DT*A(N)
        Vp(N)   = VpOLD(N)
     &    - 0.166666*( DVp1(N) + 2.0*DVp2(N) + 2.0*DVp3(N) + DVp4(N))
        Di = (Vp(N)*1000)
        Lip = (U(N)*1000)
        X(N)   = XOLD(N)
     &    + 0.166666*( DX1(N) + 2.0*DX2(N) + 2.0*DX3(N) + DX4(N))
      ENDDO
C
      UMAG = SQRT(U(1)**2 + U(2)**2 + U(3)**2)
      IF(UMAG .EQ. 0.0) THEN
        DA4 = 0.0
      ELSE
        DA4 = 0.5*DT*(UCURL(1)*U(1)+UCURL(2)*U(2)+UCURL(3)*U(3))/UMAG
      ENDIF
      DD4 = DT*UTDIV
      ACUR = AOLD + 0.166666*( DA1 + 2.0*DA2 + 2.0*DA3 + DA4)
      DCUR = DOLD + 0.166666*( DD1 + 2.0*DD2 + 2.0*DD3 + DD4)
C
C
C     Don't look at errors if passed a domain boundary
      IF(PASS) THEN
        DTOLD = DT
        XSAV(1) = XOLD(1)
        XSAV(2) = XOLD(2)
        XSAV(3) = XOLD(3)
        XOLD(1) = X(1)
        XOLD(2) = X(2)
        XOLD(3) = X(3)
        ork(1) = XOLD(1)-XSAV(1)
        ork(2) = XOLD(2)-XSAV(2)
        ork(3) = XOLD(3)-XSAV(3)
        Au(1)  = (U0(1) - UOLD(1))/DTOLD
        Au(2)  = (U0(2) - UOLD(2))/DTOLD
        Au(3)  = (U0(3) - UOLD(3))/DTOLD
        UOLD(1) = U0(1)
        UOLD(2) = U0(2)
        UOLD(3) = U0(3)
        Av(1)  = (VpOLD(1) - Vp0(1))/DTOLD
        Av(2)  = (VpOLD(2) - Vp0(2))/DTOLD
        Av(3)  = (VpOLD(3) - Vp0(3))/DTOLD
        Vp0(1) = VpOLD(1)
        Vp0(2) = VpOLD(2)
        Vp0(3) = VpOLD(3)
        write(77,*) XSAV(1), XSAV(2), XSAV(3),
     &   VpOLD(1), VpOLD(2), VpOLD(3), UOLD(1),
     &   UOLD(2), UOLD(3), A0(1), A0(2), A0(3), 
     &   Av(1), Av(2), Av(3), Au(1), Au(2), Au(3), DTOLD, KCOLD
        VpOLD(1) = Vp(1)
        VpOLD(2) = Vp(2)
        VpOLD(3) = Vp(3)
        ASAV = AOLD
        AOLD = ACUR
        DSAV = DOLD
        DOLD = DCUR
        ATIM = TIME
        TIME = TIME + DT
        PASS = .FALSE.
        GO TO 2
      ENDIF
C
C---- Check timesteps
C
C---- Used for making streamribbons(curl) better
      EPSA = ABS(ACUR-AOLD)
      IF(EPSA .GT. DELMAX) THEN
       DT = 0.5*DT
       KC = KC0
       DOUBLE = .FALSE.
       GOTO 3
      ENDIF
C----
C
C---- Error check for Adaptive RK4 method
      DX1S = DX1(1)*DX1(1) + DX1(2)*DX1(2) + DX1(3)*DX1(3)
      DX2S = DX2(1)*DX2(1) + DX2(2)*DX2(2) + DX2(3)*DX2(3)
      DX3S = DX3(1)*DX3(1) + DX3(2)*DX3(2) + DX3(3)*DX3(3)
      DX4S = DX4(1)*DX4(1) + DX4(2)*DX4(2) + DX4(3)*DX4(3)
      EPS2 = 1.0 - (DX1(1)*DX2(1) + DX1(2)*DX2(2) + DX1(3)*DX2(3))**2
     &           / DX1S / DX2S
      EPS3 = 1.0 - (DX1(1)*DX3(1) + DX1(2)*DX3(2) + DX1(3)*DX3(3))**2
     &           / DX1S / DX3S
      EPS4 = 1.0 - (DX1(1)*DX4(1) + DX1(2)*DX4(2) + DX1(3)*DX4(3))**2
     &           / DX1S / DX4S
      IF(EPS2.GT.EPLUS .OR. EPS3.GT.EPLUS .OR. EPS4.GT.EPLUS) THEN
        DT = 0.5*DT
        KC = KC0
        DOUBLE = .FALSE.
C       GOTO 3
      ENDIF
C
      DTOLD = DT
      XSAV(1) = XOLD(1)
      XSAV(2) = XOLD(2)
      XSAV(3) = XOLD(3)
      XOLD(1) = X(1)
      XOLD(2) = X(2)
      XOLD(3) = X(3)
      ork(1) = XOLD(1)-XSAV(1)
      ork(2) = XOLD(2)-XSAV(2)
      ork(3) = XOLD(3)-XSAV(3)
      Au(1)  = (U0(1) - UOLD(1))/DTOLD
      Au(2)  = (U0(2) - UOLD(2))/DTOLD
      Au(3)  = (U0(3) - UOLD(3))/DTOLD
      UOLD(1) = U0(1)
      UOLD(2) = U0(2)
      UOLD(3) = U0(3)
      Av(1)  = (VpOLD(1) - Vp0(1))/DTOLD
      Av(2)  = (VpOLD(2) - Vp0(2))/DTOLD
      Av(3)  = (VpOLD(3) - Vp0(3))/DTOLD
      Vp0(1) = VpOLD(1)
      Vp0(2) = VpOLD(2)
      Vp0(3) = VpOLD(3)
      write(77,*) XSAV(1), XSAV(2), XSAV(3),
     &   VpOLD(1), VpOLD(2), VpOLD(3), UOLD(1),
     &   UOLD(2), UOLD(3), A0(1), A0(2), A0(3), 
     &   Av(1), Av(2), Av(3), Au(1), Au(2), Au(3), DTOLD, KCOLD
      VpOLD(1) = Vp(1)
      VpOLD(2) = Vp(2)
      VpOLD(3) = Vp(3)
      ASAV = AOLD
      AOLD = ACUR
      DSAV = DOLD
      DOLD = DCUR
      ATIM = TIME
      TIME = TIME + DT
      IF(MAX(EPS2,EPS3,EPS4).LT.EMINUS .AND. EPSA*2.0 .LT. DELMAX
     &                                 .AND. DOUBLE) THEN
        DT = 2.0*DT
      ENDIF
      GO TO 2
C
C---- Hit Boundary - 1/2 TimeStep
C
 4    IF(KC .EQ. 0) THEN
       WRITE(*,*) 'StreamLine Integration Error - Lost Position!'
       GOTO 1
      ENDIF
 44   DT = 0.5*DT
      KC = KC0
      DOUBLE = .FALSE.
      PASS   = .FALSE.
      NHIT   = NHIT + 1
      IF(NHIT .LE. MHIT) GO TO 3
C
C---- terminate streamline
C
 1    IF(ID.NE.1 .AND. TIME.LT.0.0 .AND. KSTREAM.GT.1) THEN
       tosl = -time
       DO KS = 1, KSTREAM/2
         KT = KSTREAM-KS+1
         XYZS            = XYZSTREAM(1,KT)
         XYZSTREAM(1,KT) = XYZSTREAM(1,KS)
         XYZSTREAM(1,KS) = XYZS
         XYZS            = XYZSTREAM(2,KT)
         XYZSTREAM(2,KT) = XYZSTREAM(2,KS)
         XYZSTREAM(2,KS) = XYZS
         XYZS            = XYZSTREAM(3,KT)
         XYZSTREAM(3,KT) = XYZSTREAM(3,KS)
         XYZSTREAM(3,KS) = XYZS
         ICS         = CSTREAM(KT)
         CSTREAM(KT) = CSTREAM(KS)
         CSTREAM(KS) = ICS
         TS          = TSTREAM(KT) - TIME
         TSTREAM(KT) = TSTREAM(KS) - TIME
         TSTREAM(KS) = TS
         TS          = ASTREAM(KT)
         ASTREAM(KT) = ASTREAM(KS)
         ASTREAM(KS) = TS
         TS          = DSTREAM(KT)
         DSTREAM(KT) = DSTREAM(KS)
         DSTREAM(KS) = TS
         DO I = 1, 8
           FS            = FSTREAM(I,KT)
           FSTREAM(I,KT) = FSTREAM(I,KS)
           FSTREAM(I,KS) = FS
         ENDDO
       ENDDO
       KS = KSTREAM/2 + 1
       IF(MOD(KSTREAM,2).EQ.1) TSTREAM(KS) = TSTREAM(KS) - TIME
      ENDIF
C
      IF(ID.EQ.0 .AND. TIME.LT.0.0 .AND. KSTREAM.NE.MSTREAM .AND.
     &                                   KSTREAM.NE.0) THEN
       MAXOUT  = .FALSE.
       DT      = 0.5*X2/U2
       TIME    = -TIME
       KC      = CSTREAM(KSTREAM)
       XOLD(1) = XYZSTREAM(1,KSTREAM)
       XOLD(2) = XYZSTREAM(2,KSTREAM)
       XOLD(3) = XYZSTREAM(3,KSTREAM)
       AOLD    = ASTREAM(KSTREAM)
       DOLD    = DSTREAM(KSTREAM)
       CALL V3_UX(KC, XOLD, U0, XYZ, V, W, CEL1, CCEL1,
     &            CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &            NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &            %val(piblank), UCURL0, UTDIV0, PASSI)
       KC0    = KC
       LASTKC = KC
       XSAV(1) = XOLD(1)
       XSAV(2) = XOLD(2)
       XSAV(3) = XOLD(3)
       ASAV = AOLD
       DSAV = DOLD
       ATIM = TIME
       GOTO 3
      ENDIF
C
C     Remove Duplicate Points
      IF(KSTREAM .EQ. 1) THEN
        KSTREAM = 0
      ELSE
        ii = 1
        do j = 2, kstream
          xmag = sqrt( (xyzstream(1,ii)-xyzstream(1,j))**2 +
     &                 (xyzstream(2,ii)-xyzstream(2,j))**2 +
     &                 (xyzstream(3,ii)-xyzstream(3,j))**2 )
          if(xmag .ne. 0.0) then
            ii = ii + 1
            xyzstream(1,ii) = xyzstream(1,j)
            xyzstream(2,ii) = xyzstream(2,j)
            xyzstream(3,ii) = xyzstream(3,j)
            cstream(ii) = cstream(j)
            tstream(ii) = tstream(j)
            astream(ii) = astream(j)
            dstream(ii) = dstream(j)
            do k = 1,8
              fstream(k,ii) = fstream(k,j)
            enddo
          endif
        enddo
        if(ii .eq. 1) then
          kstream = 0
        else
          kstream = ii
        endif
      ENDIF
C
      RETURN
      END


      FUNCTION V3_DTLIM(KC, XYZ, V, CEL1, CEL2, CEL3, CEL4,
     &                  NPOLYT, POLYT, BLOCKS, X2, U2)
C
C     RETURN CELL LIMIT ON DELTA TIME
C
C     Copyright 1990 - 2005, Massachusetts Institute of Technology.
C
      INCLUDE 'Visual3.inc'
C
      REAL    XYZ(3,*),  V(3,*), XA(3), XN(3,8), UA(3), VN(3,8)
      INTEGER CEL1(4,*), CEL2(5,*), CEL3(6,*), CEL4(8,*)
      INTEGER NPOLYT(8,*), POLYT(*), BLOCKS(6,*)
C
      REAL    XT(3), UT(3), UGABS(3,3)
      INTEGER CEL(8)
      LOGICAL FLAG
C
C     There are two coefficients which are needed:
C
C     (1) eta = maximum fraction of cell which should be crossed in
C               one iteration
C
C     (2) ekeps = maximum eigenvalue*dt which is needed for accuracy
C
C     These two coefficients can be found in Eq. (11)
C
C     eta = 1.0 (for lack of anything better). 
C     And ekeps=0.6 which gives an amplification error of under 0.001
C     for the RK4 scheme (see figure 4).
C
      DATA Ceta/1.0/, Cekeps/0.6/
C
      XA(1) = 0.0
      XA(2) = 0.0
      XA(3) = 0.0
      UA(1) = 0.0
      UA(2) = 0.0
      UA(3) = 0.0
C
C     TETRAHEDRA
      IF(KC .LE. KCEL31) THEN
       NCN = 4
       KCR = KC
       DO NN = 1, NCN
         KN = CEL1(NN,KCR)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
         VN(1,NN) = V(1,KN)
         VN(2,NN) = V(2,KN)
         VN(3,NN) = V(3,KN)
         XA(1) = XA(1) + XN(1,NN)
         XA(2) = XA(2) + XN(2,NN)
         XA(3) = XA(3) + XN(3,NN)
         UA(1) = UA(1) + VN(1,NN)
         UA(2) = UA(2) + VN(2,NN)
         UA(3) = UA(3) + VN(3,NN)
       ENDDO
C     PYRAMIDS
      ELSEIF(KC .LE. KCEL31+KCEL32) THEN
       NCN = 5
       KCR = KC - KCEL31
       DO NN = 1, NCN
         KN = CEL2(NN,KCR)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
         VN(1,NN) = V(1,KN)
         VN(2,NN) = V(2,KN)
         VN(3,NN) = V(3,KN)
         XA(1) = XA(1) + XN(1,NN)
         XA(2) = XA(2) + XN(2,NN)
         XA(3) = XA(3) + XN(3,NN)
         UA(1) = UA(1) + VN(1,NN)
         UA(2) = UA(2) + VN(2,NN)
         UA(3) = UA(3) + VN(3,NN)
       ENDDO
C     PRISMS
      ELSEIF(KC .LE. KCEL31+KCEL32+KCEL33) THEN
       NCN = 6
       KCR = KC - KCEL31-KCEL32
       DO NN = 1, NCN
         KN = CEL3(NN,KCR)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
         VN(1,NN) = V(1,KN)
         VN(2,NN) = V(2,KN)
         VN(3,NN) = V(3,KN)
         XA(1) = XA(1) + XN(1,NN)
         XA(2) = XA(2) + XN(2,NN)
         XA(3) = XA(3) + XN(3,NN)
         UA(1) = UA(1) + VN(1,NN)
         UA(2) = UA(2) + VN(2,NN)
         UA(3) = UA(3) + VN(3,NN)
       ENDDO
C     HEXAHEDRA
      ELSEIF(KC .LE. KCEL31+KCEL32+KCEL33+KCEL34) THEN
       NCN = 8
       KCR = KC - KCEL31-KCEL32-KCEL33
       DO NN = 1, NCN
         KN = CEL4(NN,KCR)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
         VN(1,NN) = V(1,KN)
         VN(2,NN) = V(2,KN)
         VN(3,NN) = V(3,KN)
         XA(1) = XA(1) + XN(1,NN)
         XA(2) = XA(2) + XN(2,NN)
         XA(3) = XA(3) + XN(3,NN)
         UA(1) = UA(1) + VN(1,NN)
         UA(2) = UA(2) + VN(2,NN)
         UA(3) = UA(3) + VN(3,NN)
       ENDDO
C     POLY-TETRAS
      ELSEIF(KC .LE. KCELLD) THEN
       NCN = 4
       KCR = KC - KCEL31-KCEL32-KCEL33-KCEL34
       IHEAD = 0
       CEL(4) = POLYT(KCR)
       IF(CEL(4) .LT. 0) THEN
         KP = -CEL(4)
         IF(NPOLYT(1,KP) .EQ. KCR) THEN
           CEL(4) = NPOLYT(6,KP)
         ELSE
           CEL(4) = NPOLYT(5,KP)
           IHEAD = 5
         ENDIF
       ENDIF
       DO L = 3, 1, -1
         KCR = KCR - 1
         IF(IHEAD .EQ. 0) THEN
           CEL(L) = POLYT(KCR)
           IF(CEL(L) .LT. 0) THEN
             KP = -CEL(L)
             CEL(L) = NPOLYT(5,KP)
             IHEAD = 5
           ENDIF
         ELSE
           IHEAD = IHEAD - 1
           CEL(L) = NPOLYT(IHEAD,KP)
         ENDIF
       ENDDO
       DO NN = 1, NCN
         KN = CEL(NN)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
         VN(1,NN) = V(1,KN)
         VN(2,NN) = V(2,KN)
         VN(3,NN) = V(3,KN)
         XA(1) = XA(1) + XN(1,NN)
         XA(2) = XA(2) + XN(2,NN)
         XA(3) = XA(3) + XN(3,NN)
         UA(1) = UA(1) + VN(1,NN)
         UA(2) = UA(2) + VN(2,NN)
         UA(3) = UA(3) + VN(3,NN)
       ENDDO
C     STRUCTURED BLOCKS
      ELSE
       NCN = 8
       KCR = KC
       include 'getblock.inc'
       N1  = BLOCKS(1,KB)
       N2  = BLOCKS(2,KB)
       KCR = KC - KCS
       I   = MOD( KCR-1,         N1-1) + 1
       J   = MOD((KCR-1)/(N1-1), N2-1) + 1
       K   =     (KCR-1)/((N1-1)*(N2-1)) + 1
       KN  = KNS + I + (J-1)*N1 + (K-1)*N1*N2
       CEL(1) = KN
       CEL(2) = KN + 1
       CEL(3) = KN + 1 + N1
       CEL(4) = KN     + N1
       CEL(5) = KN          + N1*N2
       CEL(6) = KN + 1      + N1*N2
       CEL(7) = KN + 1 + N1 + N1*N2
       CEL(8) = KN     + N1 + N1*N2
       DO NN = 1, NCN
         KN = CEL(NN)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
         VN(1,NN) = V(1,KN)
         VN(2,NN) = V(2,KN)
         VN(3,NN) = V(3,KN)
         XA(1) = XA(1) + XN(1,NN)
         XA(2) = XA(2) + XN(2,NN)
         XA(3) = XA(3) + XN(3,NN)
         UA(1) = UA(1) + VN(1,NN)
         UA(2) = UA(2) + VN(2,NN)
         UA(3) = UA(3) + VN(3,NN)
       ENDDO
      ENDIF
C
      XA(1) = XA(1)/FLOAT(NCN)
      XA(2) = XA(2)/FLOAT(NCN)
      XA(3) = XA(3)/FLOAT(NCN)
      UA(1) = UA(1)/FLOAT(NCN)
      UA(2) = UA(2)/FLOAT(NCN)
      UA(3) = UA(3)/FLOAT(NCN)
C
      X2 = 0.0
      U2 = 0.0
      UX = 0.0
      DTVEL = 0.0
      DO NN = 1, NCN
        X2 = MAX(X2, (XN(1,NN)-XA(1))**2 +
     &               (XN(2,NN)-XA(2))**2 +
     &               (XN(3,NN)-XA(3))**2 )
        U2 = MAX(U2, (VN(1,NN)-UA(1))**2 +
     &               (VN(2,NN)-UA(2))**2 +
     &               (VN(3,NN)-UA(3))**2 )
        UX = MAX(UX, ABS(UA(1)*(XN(1,NN)-XA(1)) +
     &                   UA(2)*(XN(2,NN)-XA(2)) +
     &                   UA(3)*(XN(3,NN)-XA(3))))
C
        XT(1) = XN(1,NN)-XA(1)
        XT(2) = XN(2,NN)-XA(2)
        XT(3) = XN(3,NN)-XA(3)
	XT2   = XT(1)**2 + XT(2)**2 + XT(3)**2
        UT(1) = 0.5*(VN(1,NN)+UA(1))
        UT(2) = 0.5*(VN(2,NN)+UA(2))
        UT(3) = 0.5*(VN(3,NN)+UA(3))
        DTVEL = MAX(DTVEL, 
     &              ABS(UT(1)*XT(1) + UT(2)*XT(2) + UT(3)*XT(3))/XT2)
      ENDDO
C
      X2 = UX + SQRT(U2*X2)
      U2 = U2 + UA(1)*UA(1) + UA(2)*UA(2) + UA(3)*UA(3)
      DTVEL = 2.0/DTVEL
C
C     Now, calculate the timestep based on Ugrad
C
      DO II = 1,3
        DO JJ = 1,3
	  UGABS(II,JJ) = ABS(UGRAD(II,JJ))
        ENDDO
      ENDDO
      UGMAX = 0.0
      DUMAX = 0.0
      DO II = 1,3
	UGMAX = MAX(UGMAX, UGABS(1,II) +
     &                     UGABS(2,II) +
     &                     UGABS(3,II) )
	DUMAX = MAX(DUMAX, UGABS(II,1) +
     &                     UGABS(II,2) +
     &                     UGABS(II,3) )
      ENDDO
      EIGMAX = MAX(MIN(UGMAX,DUMAX),0.00001)
C
      V3_DTLIM = MIN(Ceta*DTVEL, Cekeps/EIGMAX)
      RETURN
      END


      SUBROUTINE V3_UX(KC, X0, U0, XYZ, V, W, CEL1, CCEL1,
     &                 CEL2, CCEL2, CEL3, CCEL3, CEL4, CCEL4,
     &                 NPOLYT, POLYT, CPOLYT, BLOCKS, CBLOCK,
     &                 IBLANK, UCURL0, UTDIV, PASSI)
C
C     RETURNS THE VECTOR QUANTITY AT THE REQUESTED POSITION
C       USES NEWTON-RAPHSON ITERATION TO GET TO POSITION
C
C	Copyright 1990 - 2005, Massachusetts Institute of Technology.
C
      INCLUDE 'Visual3.inc'
C
      INTEGER  CEL1(4,*),  CEL2(5,*),  CEL3(6,*),  CEL4(8,*)
      INTEGER CCEL1(4,*), CCEL2(5,*), CCEL3(5,*), CCEL4(6,*)
      INTEGER NPOLYT(8,*), POLYT(*), CPOLYT(2,*)
      INTEGER BLOCKS(6,*), CBLOCK(*), IBLANK(*)
      REAL    XYZ(3,*), V(3,*), W(8), UCURL0(3),  UTDIV
      LOGICAL PASSI
C
      DIMENSION XN(3,8), X0(3), U0(3), X(3), DXI(3)
      DIMENSION AA(3), BB(3), CC(3), DD(3), EE(3), FF(3), GG(3), HH(3)
      REAL      JAC(3,3), XIOLD(3), XI(3), CLOLD(6), CL(6), VJAC(3,3)
      REAL      II(3)
      INTEGER   KNE(8), FACE(4,6), KSAV(4)
      LOGICAL   BOUND, FLAG
C
      DATA FACE/1,2,3,4, 2,3,7,6, 3,4,8,7, 4,8,5,1, 5,6,7,8, 6,5,1,2/
C	  Common /VUX/ XN(3,8)
C
      NITER = 0
      PASSI = .FALSE.
 3    NITER = NITER + 1
      IF(NITER .GT. 64) THEN
       KC = 0
       RETURN
      ENDIF
C
C---- Load in nodal values
C
      KTYPE = 0
      IF(KC .LE. KCEL31) THEN
C
C----- Tetrahedra
       NCN = 4
       NCF = 4
       KCR = KC
       DO NN = 1, NCN
         KN = CEL1(NN,KCR)
         KNE(NN) = KN
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = XN(N,1)
         BB(N) = XN(N,2) - XN(N,1)
         CC(N) = XN(N,4) - XN(N,1)
         DD(N) = XN(N,3) - XN(N,1)
         EE(N) = 0.0
         FF(N) = 0.0
         GG(N) = 0.0
         HH(N) = 0.0
         II(N) = 0.0
       ENDDO
       XI(1) = 0.0
       XI(2) = 0.0
       XI(3) = 0.0
      ELSEIF(KC .LE. KCEL31+KCEL32) THEN
C
C----- Pyramids
       NCN = 5
       NCF = 5
       KCR = KC - KCEL31
       DO NN = 1, NCN
         KN = CEL2(NN,KCR)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = 0.25*(+ XN(N,1) + XN(N,2) + XN(N,3) + XN(N,4) )
         BB(N) = 0.25*(- XN(N,1) + XN(N,2) + XN(N,3) - XN(N,4) )
         CC(N) = 0.25*(- XN(N,1) - XN(N,2) + XN(N,3) + XN(N,4) )
         DD(N) = 2.0*(XN(N,5) - AA(N))
         EE(N) = 0.25*(+ XN(N,1) - XN(N,2) + XN(N,3) - XN(N,4) )
         FF(N) = -CC(N)
         GG(N) = -BB(N)
         HH(N) = 0.0
         II(N) = AA(N) - XN(N,5)
       ENDDO
       XI(1) = 0.0
       XI(2) = 0.0
       XI(3) = 0.25
      ELSEIF(KC .LE. KCEL31+KCEL32+KCEL33) THEN
C
C----- Prisms
       NCN = 6
       NCF = 5
       KCR = KC - KCEL31-KCEL32
       DO NN = 1, NCN
         KN = CEL3(NN,KCR)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = XN(N,1)
         BB(N) = XN(N,4) - XN(N,1)
         CC(N) = XN(N,6) - XN(N,1)
         DD(N) = XN(N,2) - XN(N,1)
         EE(N) = 0.0
         FF(N) = XN(N,5) - XN(N,6) - DD(N)
         GG(N) = XN(N,3) - XN(N,4) - DD(N)
         HH(N) = 0.0
         II(N) = 0.0
       ENDDO
       XI(1) = 0.3333333
       XI(2) = 0.3333333
       XI(3) = 0.5
      ELSEIF(KC .LE. KCEL31+KCEL32+KCEL33+KCEL34) THEN
C
C----- Hexahedra
       NCN = 8
       NCF = 6
       KCR = KC - KCEL31-KCEL32-KCEL33
       DO NN = 1, NCN
         KN = CEL4(NN,KCR)
         KNE(NN) = KN
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = 0.125*(+ XN(N,1) + XN(N,2) + XN(N,3) + XN(N,4)
     &                  + XN(N,5) + XN(N,6) + XN(N,7) + XN(N,8) )
         BB(N) = 0.125*(- XN(N,1) + XN(N,2) + XN(N,3) - XN(N,4)
     &                  - XN(N,5) + XN(N,6) + XN(N,7) - XN(N,8) )
         CC(N) = 0.125*(- XN(N,1) - XN(N,2) + XN(N,3) + XN(N,4)
     &                  - XN(N,5) - XN(N,6) + XN(N,7) + XN(N,8) )
         DD(N) = 0.125*(- XN(N,1) - XN(N,2) - XN(N,3) - XN(N,4)
     &                  + XN(N,5) + XN(N,6) + XN(N,7) + XN(N,8) )
         EE(N) = 0.125*(+ XN(N,1) - XN(N,2) + XN(N,3) - XN(N,4)
     &                  + XN(N,5) - XN(N,6) + XN(N,7) - XN(N,8) )
         FF(N) = 0.125*(+ XN(N,1) + XN(N,2) - XN(N,3) - XN(N,4)
     &                  - XN(N,5) - XN(N,6) + XN(N,7) + XN(N,8) )
         GG(N) = 0.125*(+ XN(N,1) - XN(N,2) - XN(N,3) + XN(N,4)
     &                  - XN(N,5) + XN(N,6) + XN(N,7) - XN(N,8) )
         HH(N) = 0.125*(- XN(N,1) + XN(N,2) - XN(N,3) + XN(N,4)
     &                  + XN(N,5) - XN(N,6) + XN(N,7) - XN(N,8) )
         II(N) = 0.0
       ENDDO
       XI(1) = 0.0
       XI(2) = 0.0
       XI(3) = 0.0
      ELSEIF(KC .LE. KCELLD) THEN
C
C----- Poly-Tetrahedral Strips
       KTYPE = 5
       NCN = 4
       NCF = 4
       KCR = KC - KCEL31-KCEL32-KCEL33-KCEL34
       KCX = KCR
       IHEAD = 0
       KHEAD = 0
       KNE(4) = POLYT(KCR)
       IF(KNE(4) .LT. 0) THEN
         KP = -KNE(4)
         IF(NPOLYT(1,KP) .EQ. KCR) THEN
           KNE(4) = NPOLYT(6,KP)
           KHEAD  = 2
         ELSE
           KNE(4) = NPOLYT(5,KP)
           IHEAD = 5
           KHEAD = 1
         ENDIF
       ENDIF
       DO L = 3, 1, -1
         KCX = KCX - 1
         IF(IHEAD .EQ. 0) THEN
           KNE(L) = POLYT(KCX)
           IF(KNE(L) .LT. 0) THEN
             KP = -KNE(L)
             KNE(L) = NPOLYT(5,KP)
             IHEAD = 5
           ENDIF
         ELSE
           IHEAD = IHEAD - 1
           KNE(L) = NPOLYT(IHEAD,KP)
         ENDIF
       ENDDO
       DO NN = 1, NCN
         KN = KNE(NN)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = XN(N,1)
         BB(N) = XN(N,2) - XN(N,1)
         CC(N) = XN(N,4) - XN(N,1)
         DD(N) = XN(N,3) - XN(N,1)
         EE(N) = 0.0
         FF(N) = 0.0
         GG(N) = 0.0
         HH(N) = 0.0
         II(N) = 0.0
       ENDDO
       XI(1) = 0.0
       XI(2) = 0.0
       XI(3) = 0.0
      ELSE
C
C----- Structured Blocks
       KTYPE = 6
       NCN = 8
       NCF = 6
       KCR = KC
       include 'getblock.inc'
       N1  = BLOCKS(1,KB)
       N2  = BLOCKS(2,KB)
       N3  = BLOCKS(3,KB)
       KCR = KC - KCS
       I   = MOD( KCR-1,         N1-1) + 1
       J   = MOD((KCR-1)/(N1-1), N2-1) + 1
       K   =     (KCR-1)/((N1-1)*(N2-1)) + 1
       KN  = KNS + I + (J-1)*N1 + (K-1)*N1*N2
       KNE(1) = KN
       KNE(2) = KN + 1
       KNE(3) = KN + 1 + N1
       KNE(4) = KN     + N1
       KNE(5) = KN          + N1*N2
       KNE(6) = KN + 1      + N1*N2
       KNE(7) = KN + 1 + N1 + N1*N2
       KNE(8) = KN     + N1 + N1*N2
       DO NN = 1, NCN
         KN = KNE(NN)
         XN(1,NN) = XYZ(1,KN)
         XN(2,NN) = XYZ(2,KN)
         XN(3,NN) = XYZ(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = 0.125*(+ XN(N,1) + XN(N,2) + XN(N,3) + XN(N,4)
     &                  + XN(N,5) + XN(N,6) + XN(N,7) + XN(N,8) )
         BB(N) = 0.125*(- XN(N,1) + XN(N,2) + XN(N,3) - XN(N,4)
     &                  - XN(N,5) + XN(N,6) + XN(N,7) - XN(N,8) )
         CC(N) = 0.125*(- XN(N,1) - XN(N,2) + XN(N,3) + XN(N,4)
     &                  - XN(N,5) - XN(N,6) + XN(N,7) + XN(N,8) )
         DD(N) = 0.125*(- XN(N,1) - XN(N,2) - XN(N,3) - XN(N,4)
     &                  + XN(N,5) + XN(N,6) + XN(N,7) + XN(N,8) )
         EE(N) = 0.125*(+ XN(N,1) - XN(N,2) + XN(N,3) - XN(N,4)
     &                  + XN(N,5) - XN(N,6) + XN(N,7) - XN(N,8) )
         FF(N) = 0.125*(+ XN(N,1) + XN(N,2) - XN(N,3) - XN(N,4)
     &                  - XN(N,5) - XN(N,6) + XN(N,7) + XN(N,8) )
         GG(N) = 0.125*(+ XN(N,1) - XN(N,2) - XN(N,3) + XN(N,4)
     &                  - XN(N,5) + XN(N,6) + XN(N,7) - XN(N,8) )
         HH(N) = 0.125*(- XN(N,1) + XN(N,2) - XN(N,3) + XN(N,4)
     &                  + XN(N,5) - XN(N,6) + XN(N,7) - XN(N,8) )
         II(N) = 0.0
       ENDDO
       XI(1) = 0.0
       XI(2) = 0.0
       XI(3) = 0.0
      ENDIF
C
C---- MAIN LOOP
C
      ITER  = 0
      BOUND = .FALSE.
C
C---- evaluate x, and dx/dxi
C
 1    ITER = ITER + 1
      IF(ITER .GT. 10) GOTO 999
C
      DO N = 1, 3
        X(N) = AA(N)
     &       + BB(N)*XI(1)       + CC(N)*XI(2)       + DD(N)*XI(3) 
     &       + EE(N)*XI(1)*XI(2) + FF(N)*XI(2)*XI(3) + GG(N)*XI(3)*XI(1)
     &       + HH(N)*XI(1)*XI(2)*XI(3) + II(N)*XI(3)*XI(3)
     &       - X0(N)
        JAC(N,1) = BB(N) + EE(N)*XI(2) + GG(N)*XI(3) + HH(N)*XI(2)*XI(3)
        JAC(N,2) = CC(N) + EE(N)*XI(1) + FF(N)*XI(3) + HH(N)*XI(1)*XI(3)
        JAC(N,3) = DD(N) + FF(N)*XI(2) + GG(N)*XI(1) + HH(N)*XI(1)*XI(2)
     &                   + 2.0*II(N)*XI(3)
      ENDDO
C
C---- newton raphson step
C
C---- -1.0/|JAC|
      DTI = -1.0 / ( JAC(1,1)*(JAC(2,2)*JAC(3,3)-JAC(2,3)*JAC(3,2))
     &             + JAC(2,1)*(JAC(3,2)*JAC(1,3)-JAC(3,3)*JAC(1,2))
     &             + JAC(3,1)*(JAC(1,2)*JAC(2,3)-JAC(1,3)*JAC(2,2)) )
C
      DXI(1) = DTI* (   X(1)  *(JAC(2,2)*JAC(3,3)-JAC(2,3)*JAC(3,2))
     &              +   X(2)  *(JAC(3,2)*JAC(1,3)-JAC(3,3)*JAC(1,2))
     &              +   X(3)  *(JAC(1,2)*JAC(2,3)-JAC(1,3)*JAC(2,2)) )
C
      DXI(2) = DTI* ( JAC(1,1)*(  X(2)  *JAC(3,3)-JAC(2,3)*  X(3)  )
     &              + JAC(2,1)*(  X(3)  *JAC(1,3)-JAC(3,3)*  X(1)  )
     &              + JAC(3,1)*(  X(1)  *JAC(2,3)-JAC(1,3)*  X(2)  ) )
C
      DXI(3) = DTI*( JAC(1,1)*(JAC(2,2)*  X(3)  -  X(2)  *JAC(3,2))
     &             + JAC(2,1)*(JAC(3,2)*  X(1)  -  X(3)  *JAC(1,2))
     &             + JAC(3,1)*(JAC(1,2)*  X(2)  -  X(1)  *JAC(2,2)) )
C
      XIOLD(1) = XI(1)
      XIOLD(2) = XI(2)
      XIOLD(3) = XI(3)
C
      XI(1) = XI(1) + DXI(1)
      XI(2) = XI(2) + DXI(2)
      XI(3) = XI(3) + DXI(3)
C
C---- clip
C
      IF(NCN .EQ. 4) THEN
       CLOLD(1) = XIOLD(2)
       CLOLD(2) = 1.0 - XIOLD(1) - XIOLD(2) - XIOLD(3)
       CLOLD(3) = XIOLD(1)
       CLOLD(4) = XIOLD(3)
       CL(1)    = XI(2)
       CL(2)    = 1.0 - XI(1) - XI(2) - XI(3)
       CL(3)    = XI(1)
       CL(4)    = XI(3)
C
      ELSE IF(NCN .EQ. 5) THEN
       CLOLD(1) = XIOLD(3)
       CLOLD(2) = 1.0 - XIOLD(1) - XIOLD(3)
       CLOLD(3) = 1.0 - XIOLD(2) - XIOLD(3)
       CLOLD(4) = 1.0 + XIOLD(1) - XIOLD(3)
       CLOLD(5) = 1.0 + XIOLD(2) - XIOLD(3)
       CL(1)    = XI(3)
       CL(2)    = 1.0 - XI(1) - XI(3)
       CL(3)    = 1.0 - XI(2) - XI(3)
       CL(4)    = 1.0 + XI(1) - XI(3)
       CL(5)    = 1.0 + XI(2) - XI(3)
C
      ELSE IF(NCN .EQ. 6) THEN
       CLOLD(1) = XIOLD(2)
       CLOLD(2) = XIOLD(1)
       CLOLD(3) = 1.0 - XIOLD(1) - XIOLD(2)
       CLOLD(4) = XIOLD(3)
       CLOLD(5) = 1.0 - XIOLD(3)
       CL(1)    = XI(2)
       CL(2)    = XI(1)
       CL(3)    = 1.0 - XI(1) - XI(2)
       CL(4)    = XI(3)
       CL(5)    = 1.0 - XI(3)
C
      ELSE
       CLOLD(1) = 1.0 + XIOLD(3)
       CLOLD(2) = 1.0 - XIOLD(1)
       CLOLD(3) = 1.0 - XIOLD(2)
       CLOLD(4) = 1.0 + XIOLD(1)
       CLOLD(5) = 1.0 - XIOLD(3)
       CLOLD(6) = 1.0 + XIOLD(2)
       CL(1)    = 1.0 + XI(3)
       CL(2)    = 1.0 - XI(1)
       CL(3)    = 1.0 - XI(2)
       CL(4)    = 1.0 + XI(1)
       CL(5)    = 1.0 - XI(3)
       CL(6)    = 1.0 + XI(2)
      ENDIF
C
      IFLAG = 0
      RFLAG = 1.0
C
      DO NF = 1, NCF
        IF(CL(NF) .LT. -0.0001) THEN
         RATIO = CLOLD(NF) / (CLOLD(NF)-CL(NF))
         IF(RATIO .LT. RFLAG) THEN
          RFLAG = RATIO
          IFLAG = NF
         ENDIF
        ENDIF
      ENDDO
C
      IF(IFLAG .EQ. 0) THEN
       IF(NCN .EQ. 4) THEN
        GOTO 2
       ELSE IF( DXI(1)**2 + DXI(2)**2 + DXI(3)**2 .GT. 0.001) THEN
        GOTO 1
       ELSE
        GOTO 2
       ENDIF
C
      ELSE
       IF(BOUND) THEN
        if(sup_con) then
          n = kc
          call v3connect(n, iflag, kc)
          if(kc .le. 0) then
            if(kc .eq. 0) kc = -99999999
            return
          endif
          goto 3
        endif
        IF(NCN .EQ. 4 .AND. KTYPE .EQ. 0) THEN
         KC = CCEL1(IFLAG,KCR)
        ELSEIF(NCN .EQ. 5) THEN
         KC = CCEL2(IFLAG,KCR)
        ELSEIF(NCN .EQ. 6) THEN
         KC = CCEL3(IFLAG,KCR)
        ELSEIF(NCN .EQ. 8 .AND. KTYPE .EQ. 0) THEN
         KC = CCEL4(IFLAG,KCR)
        ELSEIF(KTYPE .EQ. 5) THEN
         IF(IFLAG .EQ. 1) THEN
           IF(KHEAD .NE. 1) THEN
             KC = KC-1
           ELSE
             KC = NPOLYT(7,KP)
           ENDIF
         ELSEIF(IFLAG .EQ. 2) THEN
           IF(KHEAD .NE. 2) THEN
             KC = KC+1
           ELSE
             KC = NPOLYT(8,KP)
           ENDIF
         ELSE
           KC = CPOLYT(IFLAG-2,KCR)
         ENDIF
        ELSE
         IF(IFLAG .EQ. 1) THEN
           IF(K .EQ. 1) THEN
             IND = BLOCKS(6,KB) + I + (J -1)*(N1-1)
             KC  = CBLOCK(IND)
             passi = .true.
           ELSE
             KC = KC - (N1-1)*(N2-1)
           ENDIF
         ELSEIF(IFLAG .EQ. 2) THEN
           IF(I .EQ. N1-1) THEN
             IND = BLOCKS(6,KB) + J + (K -1)*(N2-1) +   (N1-1)*(N2-1)
             KC  = CBLOCK(IND)
             passi = .true.
           ELSE
             KC = KC + 1
           ENDIF
         ELSEIF(IFLAG .EQ. 3) THEN
           IF(J .EQ. N2-1) THEN
             IND = BLOCKS(6,KB) + I + (K -1)*(N1-1) +   (N1-1)*(N2-1) +
     &                                (N2-1)*(N3-1)
             KC  = CBLOCK(IND)
             passi = .true.
           ELSE
             KC = KC + (N1-1)
           ENDIF
         ELSEIF(IFLAG .EQ. 4) THEN
           IF(I .EQ. 1) THEN
             IND = BLOCKS(6,KB) + J + (K -1)*(N2-1) +   (N1-1)*(N2-1) +
     &                                (N2-1)*(N3-1) +   (N1-1)*(N3-1)
             KC  = CBLOCK(IND)
             passi = .true.
           ELSE
             KC = KC - 1
           ENDIF
         ELSEIF(IFLAG .EQ. 5) THEN
           IF(K .EQ. N3-1) THEN
             IND = BLOCKS(6,KB) + I + (J -1)*(N1-1) +   (N1-1)*(N2-1) +
     &                              2*(N2-1)*(N3-1) +   (N1-1)*(N3-1)
             KC  = CBLOCK(IND)
             passi = .true.
           ELSE
             KC = KC + (N1-1)*(N2-1)
           ENDIF
         ELSE
           IF(J .EQ. 1) THEN
             IND = BLOCKS(6,KB) + I + (K -1)*(N1-1) + 2*(N1-1)*(N2-1) +
     &                              2*(N2-1)*(N3-1) +   (N1-1)*(N3-1)
             KC  = CBLOCK(IND)
             passi = .true.
           ELSE
             KC = KC - (N1-1)
           ENDIF
         ENDIF
        ENDIF
	IF(PIBLANK .NE. 0 .AND. KC .GT. KCELLD) THEN
	  kcr = kc
	  kc  = 0
	  kbs = kb
	  include 'getblock.inc'
	  n1 = blocks(1,kb)
	  n2 = blocks(2,kb)
	  kb = kbs
	  kz = kcr - kcs
	  i  = mod( kz-1,        n1-1) + 1
	  j  = mod((kz-1)/(n1-1),n2-1) + 1
	  k  =     (kz-1)/((n1-1)*(n2-1)) + 1
	  kn = kns + i + (j-1)*n1 + (k-1)*n1*n2 - knode3d
	  if(iblank(kn                 ) .eq. 0 .or.
     &       iblank(kn + 1             ) .eq. 0 .or.
     &       iblank(kn + 1 + n1        ) .eq. 0 .or.
     &       iblank(kn     + n1        ) .eq. 0 .or.
     &       iblank(kn          + n1*n2) .eq. 0 .or.
     &       iblank(kn + 1      + n1*n2) .eq. 0 .or.
     &       iblank(kn + 1 + n1 + n1*n2) .eq. 0 .or.
     &       iblank(kn     + n1 + n1*n2) .eq. 0) then
	    if(ktype .ne. 6) return
	    ksav(1) = kne(face(1,iflag))
	    ksav(2) = kne(face(2,iflag))
	    ksav(3) = kne(face(3,iflag))
	    ksav(4) = kne(face(4,iflag))
	    flag  = passi
	    passi = .true.
	    do 5 i = 1, 4
	      kbt = -iblank(ksav(i)-knode3d)
	      if(kbt .le. 0) go to 5
	      do j = 1, i-1
	        if(kbt .eq. -iblank(ksav(j)-knode3d)) go to 5
	      enddo
	      call V3_XIB(kc, ksav, kb, kbt, X0, XYZ, BLOCKS, IBLANK)
	      if(kc .ne. 0) go to 3
 5	    continue
	    passi = flag
	    return
	  endif
	  kc = kcr
	ENDIF
	IF(KC.LT.0 .AND. KTYPE.EQ.6 .AND. PIBLANK.NE.0) THEN
	  kcr = kc
	  ksav(1) = kne(face(1,iflag))
	  ksav(2) = kne(face(2,iflag))
	  ksav(3) = kne(face(3,iflag))
	  ksav(4) = kne(face(4,iflag))
	  do 6 i = 1, 4
	    kbt = -iblank(ksav(i)-knode3d)
	    if(kbt .le. 0) go to 6
	    do j = 1, i-1
	      if(kbt .eq. -iblank(ksav(j)-knode3d)) go to 6
	    enddo
	    call V3_XIB(kc, ksav, kb, kbt, X0, XYZ, BLOCKS, IBLANK)
	    if(kc .eq. 0) go to 6
	    if(iabs(iopt) .lt. 2) then
	      do j = 1, 4
	        if(iblank(ksav(j)-knode3d) .ne. -kbt) go to 3
	      enddo
	      cblock(ind) = kc
	    endif
	    go to 3
 6	  continue
	  kc = kcr
	ENDIF
	IF(KC.LE.0) THEN
	  IF(KC .EQ. 0) KC = -99999999
	  RETURN
	ENDIF
	GOTO 3
       ELSE
        BOUND = .TRUE.
        XI(1) = XIOLD(1) + RFLAG*DXI(1)
        XI(2) = XIOLD(2) + RFLAG*DXI(2)
        XI(3) = XIOLD(3) + RFLAG*DXI(3)
        GOTO 1
       ENDIF
      ENDIF
C
 999  WRITE(*,*) 'StreamLine Integration Error - convergence failure'
C
C---- Load in weights and nodal velocity values
C
 2    IF(NCN .EQ. 4) THEN
C
C----- Tetrahedra or Poly-tetrahedral strips
       W(1) = 1.0 - XI(1) - XI(2) - XI(3)
       W(2) = XI(1)
       W(3) = XI(3)
       W(4) = XI(2)
       W(5) = 0.0
       W(6) = 0.0
       W(7) = 0.0
       W(8) = 0.0
       DO NN = 1, NCN
         KN       = KNE(NN)
         XN(1,NN) = V(1,KN)
         XN(2,NN) = V(2,KN)
         XN(3,NN) = V(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = XN(N,1)
         BB(N) = XN(N,2) - XN(N,1)
         CC(N) = XN(N,4) - XN(N,1)
         DD(N) = XN(N,3) - XN(N,1)
         EE(N) = 0.0
         FF(N) = 0.0
         GG(N) = 0.0
         HH(N) = 0.0
         II(N) = 0.0
       ENDDO
C
      ELSEIF(NCN .EQ. 5) THEN
C
C----- Pyramids
       W(1) = 0.25*(1.0 - XI(1) - XI(3))*(1.0 - XI(2) - XI(3))
       W(2) = 0.25*(1.0 + XI(1) - XI(3))*(1.0 - XI(2) - XI(3))
       W(3) = 0.25*(1.0 + XI(1) - XI(3))*(1.0 + XI(2) - XI(3))
       W(4) = 0.25*(1.0 - XI(1) - XI(3))*(1.0 + XI(2) - XI(3))
       W(5) = 1.0 - (1.0 - XI(3))*(1.0 - XI(3))
       W(6) = 0.0
       W(7) = 0.0
       W(8) = 0.0
       KNE(1) = CEL2(1,KCR)
       KNE(2) = CEL2(2,KCR)
       KNE(3) = CEL2(3,KCR)
       KNE(4) = CEL2(4,KCR)
       KNE(5) = CEL2(5,KCR)
       DO NN = 1, NCN
         KN       = KNE(NN)
         XN(1,NN) = V(1,KN)
         XN(2,NN) = V(2,KN)
         XN(3,NN) = V(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = 0.25*(+ XN(N,1) + XN(N,2) + XN(N,3) + XN(N,4) )
         BB(N) = 0.25*(- XN(N,1) + XN(N,2) + XN(N,3) - XN(N,4) )
         CC(N) = 0.25*(- XN(N,1) - XN(N,2) + XN(N,3) + XN(N,4) )
         DD(N) = 2.0*(XN(N,5) - AA(N))
         EE(N) = 0.25*(+ XN(N,1) - XN(N,2) + XN(N,3) - XN(N,4) )
         FF(N) = -CC(N)
         GG(N) = -BB(N)
         HH(N) = 0.0
         II(N) = AA(N) - XN(N,5)
       ENDDO
C
      ELSEIF(NCN .EQ. 6) THEN
C
C----- Prisms
       W(1) = (1.0 - XI(3))*(1.0 - XI(1) - XI(2))
       W(2) =        XI(3) *(1.0 - XI(1) - XI(2))
       W(3) = XI(1)*XI(3)
       W(4) = XI(1)*(1.0 - XI(3))
       W(5) = XI(2)*XI(3)
       W(6) = XI(2)*(1.0 - XI(3))
       W(7) = 0.0
       W(8) = 0.0
       KNE(1) = CEL3(1,KCR)
       KNE(2) = CEL3(2,KCR)
       KNE(3) = CEL3(3,KCR)
       KNE(4) = CEL3(4,KCR)
       KNE(5) = CEL3(5,KCR)
       KNE(6) = CEL3(6,KCR)
       DO NN = 1, NCN
         KN       = KNE(NN)
         XN(1,NN) = V(1,KN)
         XN(2,NN) = V(2,KN)
         XN(3,NN) = V(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = XN(N,1)
         BB(N) = XN(N,4) - XN(N,1)
         CC(N) = XN(N,6) - XN(N,1)
         DD(N) = XN(N,2) - XN(N,1)
         EE(N) = 0.0
         FF(N) = XN(N,5) - XN(N,6) - DD(N)
         GG(N) = XN(N,3) - XN(N,4) - DD(N)
         HH(N) = 0.0
         II(N) = 0.0
       ENDDO
C
      ELSE
C
C----- Hexahedra or blocks
       W(1) = 0.125*(1.0 - XI(1))*(1.0 - XI(2))*(1.0 - XI(3))
       W(2) = 0.125*(1.0 + XI(1))*(1.0 - XI(2))*(1.0 - XI(3))
       W(3) = 0.125*(1.0 + XI(1))*(1.0 + XI(2))*(1.0 - XI(3))
       W(4) = 0.125*(1.0 - XI(1))*(1.0 + XI(2))*(1.0 - XI(3))
       W(5) = 0.125*(1.0 - XI(1))*(1.0 - XI(2))*(1.0 + XI(3))
       W(6) = 0.125*(1.0 + XI(1))*(1.0 - XI(2))*(1.0 + XI(3))
       W(7) = 0.125*(1.0 + XI(1))*(1.0 + XI(2))*(1.0 + XI(3))
       W(8) = 0.125*(1.0 - XI(1))*(1.0 + XI(2))*(1.0 + XI(3))
       DO NN = 1, NCN
         KN       = KNE(NN)
         XN(1,NN) = V(1,KN)
         XN(2,NN) = V(2,KN)
         XN(3,NN) = V(3,KN)
       ENDDO
       DO N = 1, 3
         AA(N) = 0.125*(+ XN(N,1) + XN(N,2) + XN(N,3) + XN(N,4)
     &                  + XN(N,5) + XN(N,6) + XN(N,7) + XN(N,8) )
         BB(N) = 0.125*(- XN(N,1) + XN(N,2) + XN(N,3) - XN(N,4)
     &                  - XN(N,5) + XN(N,6) + XN(N,7) - XN(N,8) )
         CC(N) = 0.125*(- XN(N,1) - XN(N,2) + XN(N,3) + XN(N,4)
     &                  - XN(N,5) - XN(N,6) + XN(N,7) + XN(N,8) )
         DD(N) = 0.125*(- XN(N,1) - XN(N,2) - XN(N,3) - XN(N,4)
     &                  + XN(N,5) + XN(N,6) + XN(N,7) + XN(N,8) )
         EE(N) = 0.125*(+ XN(N,1) - XN(N,2) + XN(N,3) - XN(N,4)
     &                  + XN(N,5) - XN(N,6) + XN(N,7) - XN(N,8) )
         FF(N) = 0.125*(+ XN(N,1) + XN(N,2) - XN(N,3) - XN(N,4)
     &                  - XN(N,5) - XN(N,6) + XN(N,7) + XN(N,8) )
         GG(N) = 0.125*(+ XN(N,1) - XN(N,2) - XN(N,3) + XN(N,4)
     &                  - XN(N,5) + XN(N,6) + XN(N,7) - XN(N,8) )
         HH(N) = 0.125*(- XN(N,1) + XN(N,2) - XN(N,3) + XN(N,4)
     &                  + XN(N,5) - XN(N,6) + XN(N,7) - XN(N,8) )
         II(N) = 0.0
       ENDDO
C
      ENDIF
C
C	*****diags
C     U0(1) = 0.0
C     U0(2) = 0.0
C     U0(3) = 0.0
C     DO NN = 1, NCN
C       U0(1) = U0(1) + W(NN)*XN(1,NN)
C       U0(2) = U0(2) + W(NN)*XN(2,NN)
C       U0(3) = U0(3) + W(NN)*XN(3,NN)
C     ENDDO
C	dis = sqrt((x0(1)-u0(1))**2+(x0(2)-u0(2))**2+(x0(3)-u0(3))**2)
C       if(dis .gt. 1.e-5) write(*,*) 'not at target',x0,u0,ncn
C	*********
C
      U0(1) = 0.0
      U0(2) = 0.0
      U0(3) = 0.0
      DO NN = 1, NCN
        KN    = KNE(NN)
        U0(1) = U0(1) + W(NN)*V(1,KN)
        U0(2) = U0(2) + W(NN)*V(2,KN)
        U0(3) = U0(3) + W(NN)*V(3,KN)
      ENDDO
C
C
      DO N = 1, 3
        VJAC(N,1)=BB(N) + EE(N)*XI(2) + GG(N)*XI(3) + HH(N)*XI(2)*XI(3)
        VJAC(N,2)=CC(N) + EE(N)*XI(1) + FF(N)*XI(3) + HH(N)*XI(1)*XI(3)
        VJAC(N,3)=DD(N) + FF(N)*XI(2) + GG(N)*XI(1) + HH(N)*XI(1)*XI(2)
     &                  + 2.0*II(N)*XI(3)
      ENDDO
C
      DTI = 1.0 / ( JAC(1,1)*(JAC(2,2)*JAC(3,3)-JAC(2,3)*JAC(3,2))
     &            + JAC(2,1)*(JAC(3,2)*JAC(1,3)-JAC(3,3)*JAC(1,2))
     &            + JAC(3,1)*(JAC(1,2)*JAC(2,3)-JAC(1,3)*JAC(2,2)) )
C
      UGRAD(1,1) = DTI*(
     &        VJAC(1,1)*(JAC(2,2)*JAC(3,3)-JAC(3,2)*JAC(2,3)) +
     &        VJAC(1,2)*(JAC(2,3)*JAC(3,1)-JAC(3,3)*JAC(2,1)) +
     &        VJAC(1,3)*(JAC(2,1)*JAC(3,2)-JAC(3,1)*JAC(2,2)) )
      UGRAD(2,1) = DTI*(
     &        VJAC(2,1)*(JAC(2,2)*JAC(3,3)-JAC(3,2)*JAC(2,3)) +
     &        VJAC(2,2)*(JAC(2,3)*JAC(3,1)-JAC(3,3)*JAC(2,1)) +
     &        VJAC(2,3)*(JAC(2,1)*JAC(3,2)-JAC(3,1)*JAC(2,2)) )
      UGRAD(3,1) = DTI*(
     &        VJAC(3,1)*(JAC(2,2)*JAC(3,3)-JAC(3,2)*JAC(2,3)) +
     &        VJAC(3,2)*(JAC(2,3)*JAC(3,1)-JAC(3,3)*JAC(2,1)) +
     &        VJAC(3,3)*(JAC(2,1)*JAC(3,2)-JAC(3,1)*JAC(2,2)) )
C
      UGRAD(1,2) = DTI*(
     &        VJAC(1,1)*(JAC(3,2)*JAC(1,3)-JAC(3,3)*JAC(1,2)) +
     &        VJAC(1,2)*(JAC(3,3)*JAC(1,1)-JAC(3,1)*JAC(1,3)) +
     &        VJAC(1,3)*(JAC(3,1)*JAC(1,2)-JAC(3,2)*JAC(1,1)) )
      UGRAD(2,2) = DTI*(
     &        VJAC(2,1)*(JAC(3,2)*JAC(1,3)-JAC(3,3)*JAC(1,2)) +
     &        VJAC(2,2)*(JAC(3,3)*JAC(1,1)-JAC(3,1)*JAC(1,3)) +
     &        VJAC(2,3)*(JAC(3,1)*JAC(1,2)-JAC(3,2)*JAC(1,1)) )
      UGRAD(3,2) = DTI*(
     &        VJAC(3,1)*(JAC(3,2)*JAC(1,3)-JAC(3,3)*JAC(1,2)) +
     &        VJAC(3,2)*(JAC(3,3)*JAC(1,1)-JAC(3,1)*JAC(1,3)) +
     &        VJAC(3,3)*(JAC(3,1)*JAC(1,2)-JAC(3,2)*JAC(1,1)) )
C
      UGRAD(1,3) = DTI*(
     &        VJAC(1,1)*(JAC(1,2)*JAC(2,3)-JAC(2,2)*JAC(1,3)) +
     &        VJAC(1,2)*(JAC(2,1)*JAC(1,3)-JAC(1,1)*JAC(2,3)) +
     &        VJAC(1,3)*(JAC(1,1)*JAC(2,2)-JAC(2,1)*JAC(1,2)) )
      UGRAD(2,3) = DTI*(
     &        VJAC(2,1)*(JAC(1,2)*JAC(2,3)-JAC(2,2)*JAC(1,3)) +
     &        VJAC(2,2)*(JAC(2,1)*JAC(1,3)-JAC(1,1)*JAC(2,3)) +
     &        VJAC(2,3)*(JAC(1,1)*JAC(2,2)-JAC(2,1)*JAC(1,2)) )
      UGRAD(3,3) = DTI*(
     &        VJAC(3,1)*(JAC(1,2)*JAC(2,3)-JAC(2,2)*JAC(1,3)) +
     &        VJAC(3,2)*(JAC(2,1)*JAC(1,3)-JAC(1,1)*JAC(2,3)) +
     &        VJAC(3,3)*(JAC(1,1)*JAC(2,2)-JAC(2,1)*JAC(1,2)) )
C
C     Calculate ucurl
C
      UCURL0(1) = UGRAD(3,2) - UGRAD(2,3)
      UCURL0(2) = UGRAD(1,3) - UGRAD(3,1)
      UCURL0(3) = UGRAD(2,1) - UGRAD(1,2)
C
C     Calculate transverse div
C
      TMAG  = SQRT(U0(1)*U0(1) + U0(2)*U0(2) + U0(3)*U0(3))
      IF(TMAG .EQ. 0.0) TMAG = 1.0
      TX    = U0(1)/TMAG
      TY    = U0(2)/TMAG
      TZ    = U0(3)/TMAG
      WPX   = UGRAD(1,1)*TX + UGRAD(2,1)*TY + UGRAD(3,1)*TZ
      WPY   = UGRAD(1,2)*TX + UGRAD(2,2)*TY + UGRAD(3,2)*TZ
      WPZ   = UGRAD(1,3)*TX + UGRAD(2,3)*TY + UGRAD(3,3)*TZ
      UTDIV = 0.5*(UGRAD(1,1)    + UGRAD(2,2)    + UGRAD(3,3)    -
     &                   (WPX*TX +        WPY*TY +        WPZ*TZ) )
C
      RETURN
      END


	subroutine V3_SStream(ixp, iyp, id, xyz, u2d, scel, sccel,
     &                        ixy, cell, xy, ic, iblank)
c
c	calculates surface streamline starting form position ixp iyp
c
c	Copyright 1990 - 2005, Massachusetts Institute of Technology.
c
	include 'Visual3.inc'
c
c	NOTE: ixy and u2d share storage!
	real    xyz(3,*),  u2d(3,*),   xy(2,*)
	integer scel(4,*), sccel(5,*), cell(3,*), ic(*), ixy(3,*)
	integer iblank(*)
c
	data rad/0.017453292/
c
c------ Initialize position and cell
c
	ca   = cos(angw*rad)
	sa   = sin(angw*rad)
	xc   = 2.0*s2x*float(ixp)/float(npixx-1) - s2x
	yc   = 2.0*s2y*float(npixy-iyp-1)/float(npixy-1) - s2y
	xp   = xc*halfw
	yp   = yc*abs(halfw)
	xcur = xp*ca + yp*sa + xpc
	ycur =-xp*sa + yp*ca + ypc
	call V3_Intgr8SSl(xcur, ycur, id, xyz, u2d, scel, sccel,
     &                    ixy, cell, xy, ic, iblank)
	return
	end


	subroutine V3_CalcSSL(w, kc, id, xyz, u2d, scel, sccel,
     &                        ixy, cell, xy, ic, iblank)
c
c	calculates surface streamline starting form position in weights
c
c	Copyright 1990 - 2005, Massachusetts Institute of Technology.
c
	include 'Visual3.inc'
c
c	NOTE: ixy and u2d share storage!
	real    w(*), xyz(3,*), u2d(3,*), xy(2,*)
	integer scel(4,*), sccel(5,*), cell(3,*), ic(*), ixy(3,*)
	integer iblank(*)
c
c------ Initialize position and cell
c
	ks = 0
	do i = 1, kcell
	  if(ks .eq. 0) then
	    if (ic(i) .eq. -kc) ks = i
	  endif
	enddo
	if (ks .eq. 0) return
c
	kn1 = cell(1,ks)
	kn2 = cell(2,ks)
	kn3 = cell(3,ks)
	kn4 = kn1
	if(ks .ne. kcell) then
	  if(ic(ks) .eq. ic(ks+1)) kn4 = cell(2,ks+1)
	endif
	xcur = xy(1,kn1)*w(1) + xy(1,kn2)*w(2) + 
     &         xy(1,kn3)*w(3) + xy(1,kn4)*w(4)
	ycur = xy(2,kn1)*w(1) + xy(2,kn2)*w(2) + 
     &         xy(2,kn3)*w(3) + xy(2,kn4)*w(4)
	call V3_Intgr8SSl(xcur, ycur, id, xyz, u2d, scel, sccel,
     &                    ixy, cell, xy, ic, iblank)
	return
	end


	subroutine V3_Intgr8SSl(xcur, ycur, id, xyz, u2d, scel, sccel,
     &                          ixy, cell, xy, ic, iblank)
c
c	calculates surface streamline starting form position xcur ycur
c
c	Copyright 1990 - 2005, Massachusetts Institute of Technology.
c
	parameter (mcount = 20)
	include 'XFtn.inc'
	include 'Visual3.inc'
c
c	NOTE: ixy and u2d share storage!
	real    xyz(3,*),  u2d(3,*),   xy(2,*)
	integer scel(4,*), sccel(5,*), cell(3,*), ic(*), ixy(3,*)
	integer iblank(*)
c
	dimension frac(3), f2stream(2,1)
	dimension xy1(2), xy2(2), xy3(2), u1(2), u2(2), u3(2)
	equivalence (f2stream, fstream)
c
c------ Initialize position and cell
c
	kstream = 0
	kcount  = 0
	tosl    = 0.0
	call V3_GetTri(xcur, ycur, 1, cell, xy, kcc, frac)
	if(kcc .eq. 0) return
	              sense =  1.0
	if(id .ne. 1) sense = -1.0
c
 2	xc =  xcur
	yc =  ycur
	kc =  kcc
	kx = -ic(kc)
	if(id .ne. 0 .or. sense .ne. 1.0) then
	  kstream       = 1
	  f2stream(1,1) = frac(1)
	  f2stream(2,1) = frac(2)
	  cstream(1)    = kc
	  tstream(1)    = 0.0
	endif
c
c	setup local data for the tri
 1	                                          itype =  0
	if(kc .ne.     1 .and. kx .eq. -ic(kc-1)) itype = -1
	if(kc .ne. kcell .and. kx .eq. -ic(kc+1)) itype =  1
c
	ki1    = cell(1,kc)
	xy1(1) =        xy(1,ki1)
	xy1(2) =        xy(2,ki1)
	u1(1)  = sense*u2d(2,ki1)
	u1(2)  = sense*u2d(3,ki1)
	umag   = u1(1) + u1(2)
	if(umag/umag .ne. 1.0) go to 3
c
	ki2    = cell(2,kc)
	xy2(1) =        xy(1,ki2)
	xy2(2) =        xy(2,ki2)
	u2(1)  = sense*u2d(2,ki2)
	u2(2)  = sense*u2d(3,ki2)
	umag   = u2(1) + u2(2)
	if(umag/umag .ne. 1.0) go to 3
c
	ki3    = cell(3,kc)
	xy3(1) =        xy(1,ki3)
	xy3(2) =        xy(2,ki3)
	u3(1)  = sense*u2d(2,ki3)
	u3(2)  = sense*u2d(3,ki3)
	umag   = u3(1) + u3(2)
	if(umag/umag .ne. 1.0) go to 3
c
c	integrate surface streamline for this tri
	xyswork(1,1) = xc
	xyswork(2,1) = yc
	astream(1)   = tstream(kstream)
	np = min(mwork, mstream-kstream)/2
	if(np .lt. 4) go to 3
	call V3_TriStream(xy1, xy2, xy3, u1, u2, u3,
     &                    xyswork, astream, np, nf)
	if(np .eq. 0) go to 3
	np  = np + 1
	dx1 = xy1(1) - xy3(1)
	dy1 = xy1(2) - xy3(2)
	dx2 = xy2(1) - xy3(1)
	dy2 = xy2(2) - xy3(2)
	a3  = dx1*dy2-dy1*dx2
	if(a3 .eq. 0) go to 3
	do i = 2, np
	  if(xyswork(1,i-1) .eq. xyswork(1,i) .and.
     &       xyswork(2,i-1) .eq. xyswork(2,i)) go to 3
	  dxx = xyswork(1,i) - xy3(1)
	  dyy = xyswork(2,i) - xy3(2)
	  w1  = (dxx*dy2-dyy*dx2) / a3
	  w2  =-(dxx*dy1-dyy*dx1) / a3
	  if(w1    .lt. -0.001 .or. w2 .lt. -0.001 .or.
     &       w1+w2 .gt.  1.001) go to 3
	  kstream             = kstream + 1
	  cstream(kstream)    = kc
	  f2stream(1,kstream) = w1
	  f2stream(2,kstream) = w2
	  tstream(kstream)    = astream(i)
	  if(kstream .ge. mstream) go to 3
	enddo
	if(nf .eq. 0) go to 3
	if(kstream .gt. 2) then
	  if(cstream(kstream) .eq. cstream(kstream-2)) then
	    if(cstream(kstream) .ne. cstream(kstream-1)) go to 3
	  endif
	endif
c
c	move to next facet
	xc = xyswork(1,np)
	yc = xyswork(2,np)
	if(iabs(itype) .eq. 1 .and. nf .eq. 2) then
	  kc    = kc + itype
	  itype = -itype
	else
	  if(itype .eq. 0) then
	    iface = nf + 1
	    if(iface .eq. 4) iface = 1
	    ifac1 = iface+1
	    if(ifac1 .eq. 4) ifac1 = 1
	  else
	    if(itype .eq. 1) then
	      iface = 2
              if(nf .eq. 3) iface = 1
	    else
	      iface = 4
	      if(nf .eq. 1) iface = 3
	    endif
	    ifac1 = iface+1
	    if(ifac1 .eq. 5) ifac1 = 1
	  endif
	  kn1 = scel(iface,kx)
	  kn2 = scel(ifac1,kx)
	  kx = sccel(iface,kx)
	  if(kx .le. 0) go to 3
	  if(piblank .ne. 0) then
	    kn = scel(1,kx)
	    if(kn .gt. knode3d) then
	      if(iblank(kn-knode3d) .eq. 0) go to 3
	    endif
	    kn = scel(2,kx)
	    if(kn .gt. knode3d) then
	      if(iblank(kn-knode3d) .eq. 0) go to 3
	    endif
	    kn = scel(3,kx)
	    if(kn .gt. knode3d) then
	      if(iblank(kn-knode3d) .eq. 0) go to 3
	    endif
	    kn = scel(4,kx)
	    if(kn .gt. knode3d) then
	      if(iblank(kn-knode3d) .eq. 0) go to 3
	    endif
	  endif
	  kc = sccel(5,kx)
c	  check for quad face and pick appropriate tri
	  itype = 0
	  if(kc .ne. kcell .and. kx .eq. -ic(kc+1)) then
	    itype = 1
	    if(kn1 .eq. scel(4,kx) .or. kn2 .eq. scel(4,kx)) then
	      itype = -1
	      kc    = kc + 1
	    endif
	  endif
	endif
c
c	check for escape
c
	if(mod(kcount, mcount) .eq. 0) then
	  IF(XFtnGetEvent(display, key_window, ButtonRelease,
     &                    ie, je, istate) .eq. 0) then
	    if(ie.ge.10 .and. ie.le.26 .and.
     &         je.ge.67 .and. je.le.83) then
	      kstream = -1
	      return
	    endif
	  endif
	endif
	kcount = kcount + 1
	go to 1
c
 3	if(sense .eq. -1.0) then
	  tosl = tstream(kstream)
c$dir no_recurrence
	  do i = 1, kstream/2
	    j = kstream - i + 1
	    k          = cstream(i)
	    cstream(i) = cstream(j)
	    cstream(j) = k
	    sav           = f2stream(1,i)
	    f2stream(1,i) = f2stream(1,j)
	    f2stream(1,j) = sav
	    sav           = f2stream(2,i)
	    f2stream(2,i) = f2stream(2,j)
	    f2stream(2,j) = sav
	    sav        = tosl - tstream(i)
	    tstream(i) = tosl - tstream(j)
	    tstream(j) = sav
	  enddo
	  ks = kstream/2 + 1
	  if(mod(kstream,2).eq.1) tstream(ks) = tosl - tstream(ks)
	  if(id .eq. 0) then
	    sense = 1.0
	    if(kstream .ge. mstream) then
	      ks = 2
	      do j = 4, kstream, 2
	        ks = ks + 1
	        cstream(ks)    = cstream(j)
	        f2stream(1,ks) = f2stream(1,j)
	        f2stream(2,ks) = f2stream(2,j)
	        tstream(ks)    = tstream(j)
	      enddo
	      kstream = ks
	    endif
	    go to 2
	  endif
	endif
c
	if(kstream .eq. 1) kstream = 0
	return
	end


      SUBROUTINE V3_TriStream(X1, X2, X3, U1, U2, U3, XP, TS, NP, NF)
C
C     X1,X2,X3  coordinates of corners of triangle
C     U1,U2,U3  velocities at corners of triangle
C     XP        array of streamline points (first is input; others output)
C     NP        input: max number of new points to be added
C               output: number of new points added
C     NF        number of exiting face (for identification of new triangle)
C		1 -> exits 2-3; 2 -> exits 1-3; 3 -> exits 1-2
C
C	Copyright 1990 - 2005, Massachusetts Institute of Technology.
C
      DIMENSION X1(2), X2(2), X3(2), U1(2), U2(2), U3(2), XP(2,*)
      DIMENSION TS(*)
C
      DIMENSION A(2,2), AI(2,2), B(2,2), C(2,2), C2(2,2),
     &          V(2), VP(2), VP2(2), Z(2), ZP(2)
C
      NPMAX = NP
      NP = 0
      TSTART = TS(1)
C
C---- Setup A,B
C
      A(1,1) = X2(1) - X1(1)
      A(1,2) = X3(1) - X1(1)
      A(2,1) = X2(2) - X1(2)
      A(2,2) = X3(2) - X1(2)
C
      B(1,1) = U2(1) - U1(1)
      B(1,2) = U3(1) - U1(1)
      B(2,1) = U2(2) - U1(2)
      B(2,2) = U3(2) - U1(2)
C
C---- Calc A^{-1}
C
      DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      IF(DET .EQ. 0.0) RETURN
      DETINV = 1.0 / DET
      AI(1,1) =  A(2,2)*DETINV
      AI(1,2) = -A(1,2)*DETINV
      AI(2,1) = -A(2,1)*DETINV
      AI(2,2) =  A(1,1)*DETINV
C
C---- Calc C, V
C
      C(1,1) = AI(1,1)*B(1,1) + AI(1,2)*B(2,1)
      C(1,2) = AI(1,1)*B(1,2) + AI(1,2)*B(2,2)
      C(2,1) = AI(2,1)*B(1,1) + AI(2,2)*B(2,1)
      C(2,2) = AI(2,1)*B(1,2) + AI(2,2)*B(2,2)
C
      V(1) = AI(1,1)*U1(1) + AI(1,2)*U1(2) 
      V(2) = AI(2,1)*U1(1) + AI(2,2)*U1(2) 
C
C---- Calc magnitude of C (square root of largest eigenvalue of C^T C)
C
      C2(1,1) = C(1,1)*C(1,1) + C(2,1)*C(2,1)
      C2(1,2) = C(1,1)*C(1,2) + C(2,1)*C(2,2)
      C2(2,2) = C(1,2)*C(1,2) + C(2,2)*C(2,2)
      CMAG = 0.5*( C2(1,1)+C2(2,2) + 
     &             SQRT((C2(1,1)-C2(2,2))**2 + 4.0*C2(1,2)**2) )
C
C---- Calc maximum acceptable time step 
C     (0.2 value is important and subject to future change!)
C
      TMAX  = 0.2 / (SQRT(CMAG) + 1.0E-20)
C
C---- Calc starting value of Z
C
      Z(1) = AI(1,1)*(XP(1,1)-X1(1)) + AI(1,2)*(XP(2,1)-X1(2))
      Z(2) = AI(2,1)*(XP(1,1)-X1(1)) + AI(2,2)*(XP(2,1)-X1(2))
C
C---- Begin main loop
C
 1    T  = TMAX
      NF = 0
C
C---- Calc VP, VP2
C
      VP(1) = V(1) + C(1,1)*Z(1) + C(1,2)*Z(2)
      VP(2) = V(2) + C(2,1)*Z(1) + C(2,2)*Z(2)
C
      VP2(1) = C(1,1)*VP(1) + C(1,2)*VP(2)
      VP2(2) = C(2,1)*VP(1) + C(2,2)*VP(2)
C
C---- Check for crossing boundary 1
C
      AA = -0.5*(VP2(1)+VP2(2))
      BB = -(VP(1)+VP(2))
      CC = -(Z(1)+Z(2)-1.0)
      CALL V3_Quad(AA,BB,CC,T,NF,1)
C
C---- Check for crossing boundary 2
C
      AA = 0.5*VP2(1)
      BB = VP(1)
      CC = Z(1)
      CALL V3_Quad(AA,BB,CC,T,NF,2)
C
C---- Check for crossing boundary 3
C
      AA = 0.5*VP2(2)
      BB = VP(2)
      CC = Z(2)
      CALL V3_Quad(AA,BB,CC,T,NF,3)
C
C---- Check for streamline failure
C
      IF(T.LT.0.) THEN
       WRITE(*,*) 'Streamline failure: point is outside triangle'
       NF = 0
       RETURN
      ENDIF
C
C---- Add new streamline points
C
      NMAX = NINT(3.0*T/TMAX) + 1
      DT   = T / FLOAT(NMAX)
      TSAV = T
C
      DO 2 N = 1, NMAX
        T  = FLOAT(N)*DT
        T2 = 0.5*T*T
        ZP(1) = Z(1) + T*VP(1) + T2*VP2(1)
        ZP(2) = Z(2) + T*VP(2) + T2*VP2(2)
        NP = NP+1
        XP(1,NP+1) = X1(1) + A(1,1)*ZP(1) + A(1,2)*ZP(2)
        XP(2,NP+1) = X1(2) + A(2,1)*ZP(1) + A(2,2)*ZP(2)
        TS(NP+1) = TSTART + T
 2    CONTINUE
      TSTART = TSTART + TSAV
C
C---- Loop back if necessary and possible
C
      IF(NF.EQ.0 .AND. NP.LT.NPMAX) THEN
       Z(1) = ZP(1)
       Z(2) = ZP(2)
       GOTO 1
      ENDIF
C
      RETURN
      END



      SUBROUTINE V3_Quad(A, B, C, TOLD, NFOLD, NF)
C
C	Copyright 1990 - 2005, Massachusetts Institute of Technology.
C
C---- Find positive crossing time T for C>0
C
      IF(C.GT.0.) THEN
C
       IF(A.GE.0.) THEN
        IF(B.GE.0.) THEN
         RETURN
        ELSE
         D = B*B - 4.0*A*C
         IF(D.LT.0.) RETURN
         T = (2.0*C) / (-B+SQRT(D))
        ENDIF
       ELSE
        IF(B.GE.0.) THEN
         D = B*B - 4.0*A*C
         T = (-B-SQRT(D))/(2.0*A)
        ELSE
         D = B*B - 4.0*A*C
         T = (2.0*C) / (-B+SQRT(D))
        ENDIF
       ENDIF
C
C---- Find positive crossing time T for C<0, and checks for errors
C
      ELSE
C
       IF(B.LE.0.) THEN
        T = -999.0
        RETURN
C        STOP 'ERROR: C<0 and B<0'
       ELSE
        IF(A.GE.0) THEN
         RETURN
        ELSE
         D = B*B - 4.0*A*C
         IF(D.LE.0) THEN
          T = -999.0
          RETURN
C         STOP 'ERROR: C<0 and D<0'
         ELSE
          T = (-B-SQRT(D))/(2.0*A)
         ENDIF
        ENDIF
       ENDIF
C
      ENDIF
C
C---- Check if T is less than old T
C
      IF(T.LT.TOLD) THEN
       TOLD  = T
       NFOLD = NF
      ENDIF
C
      RETURN
      END


      FUNCTION C_DRAG(Rer, rho, mu, T)
        Real Rer, lf, Knd, rho(8), T(8), mu
        Real rho_av, T_av
        PARAMETER (kB=1.38e-23, R=8.314, df=4e-10)
        Common /Orkwis/ rhop, dp

        rho_av = SUM(rho)/SIZE(rho)
        T_av = SUM(T)/SIZE(T)
        Knd = 1.26*(SQRT(mu*1.4/(rho_av*dp*SQRT(1.4*R*T_av))))
c	print*, Knd, rho_av, T_av, mu
        IF(ABS(Rer) .LE. 1e-7) THEN
            C_DRAG = 0
        ELSE
            C_DRAG = 24/Rer
        ENDIF
      RETURN
      END
