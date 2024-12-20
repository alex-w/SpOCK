C***********************************************************************
C
C
C     SUBROUTINE GEOMAG (GEOMAGNETIC FIELD COMPUTATION)
C
C
C***********************************************************************
C
C     GEOMAG IS A NATIONAL GEOSPATIAL INTELLIGENCE AGENCY (NGA) STANDARD
C     PRODUCT.  IT IS COVERED UNDER NGA MILITARY SPECIFICATION:
C     MIL-W-89500 (1993).
C
C***********************************************************************
C     Contact Information
C
C     Software and Model Support
C     	National Geophysical Data Center
C     	NOAA EGC/2
C     	325 Broadway
C     	Boulder, CO 80303 USA
C     	Attn: Susan McLean or Stefan Maus
C     	Phone:  (303) 497-6478 or -6522
C     	Email:  Susan.McLean@noaa.gov or Stefan.Maus@noaa.gov
C		Web: http://www.ngdc.noaa.gov/seg/WMM/
C
C     Sponsoring Government Agency
C	   National Geospatial-Intelligence Agency
C    	   PRG / CSAT, M.S. L-41
C    	   3838 Vogel Road
C    	   Arnold, MO 63010
C    	   Attn: Craig Rollins
C    	   Phone:  (314) 263-4186
C    	   Email:  Craig.M.Rollins@Nga.Mil
C
C      Original Program By:
C        Dr. John Quinn
C        FLEET PRODUCTS DIVISION, CODE N342
C        NAVAL OCEANOGRAPHIC OFFICE (NAVOCEANO)
C        STENNIS SPACE CENTER (SSC), MS 39522-5001
C
C***********************************************************************
C
C     PURPOSE:  THIS ROUTINE COMPUTES THE DECLINATION (DEC),
C               INCLINATION (DIP), TOTAL INTENSITY (TI) AND
C               GRID VARIATION (GV - POLAR REGIONS ONLY, REFERENCED
C               TO GRID NORTH OF A STEREOGRAPHIC PROJECTION) OF THE
C               EARTH'S MAGNETIC FIELD IN GEODETIC COORDINATES
C               FROM THE COEFFICIENTS OF THE CURRENT OFFICIAL
C               DEPARTMENT OF DEFENSE (DOD) SPHERICAL HARMONIC WORLD
C               MAGNETIC MODEL (WMM.COF).  THE WMM SERIES OF MODELS IS
C               UPDATED EVERY 5 YEARS ON JANUARY 1ST OF THOSE YEARS
C               WHICH ARE DIVISIBLE BY 5 (I.E. 2000, 2005, 2010 ETC.)
C               BY NOAA'S NATIONAL GEOPHYSICAL DATA CENTER IN
C               COOPERATION WITH THE BRITISH GEOLOGICAL SURVEY (BGS).
C               THE MODEL IS BASED ON GEOMAGNETIC FIELD MEASUREMENTS
C               FROM SATELLITE AND GROUND OBSERVATORIES.
C
C***********************************************************************
C
C     MODEL:  THE WMM SERIES GEOMAGNETIC MODELS ARE COMPOSED
C             OF TWO PARTS:  THE MAIN FIELD MODEL, WHICH IS
C             VALID AT THE BASE EPOCH OF THE CURRENT MODEL AND
C             A SECULAR VARIATION MODEL, WHICH ACCOUNTS FOR SLOW
C             TEMPORAL VARIATIONS IN THE MAIN GEOMAGNETIC FIELD
C             FROM THE BASE EPOCH TO A MAXIMUM OF 5 YEARS BEYOND
C             THE BASE EPOCH.  FOR EXAMPLE, THE BASE EPOCH OF
C             THE WMM-2005 MODEL IS 2005.0.  THIS MODEL IS THEREFORE
C             CONSIDERED VALID BETWEEN 2005.0 AND 2010.0. THE
C             COMPUTED MAGNETIC PARAMETERS ARE REFERENCED TO THE
C             WGS-84 ELLIPSOID.
C
C***********************************************************************
C
C     ACCURACY:  IN OCEAN AREAS AT THE EARTH'S SURFACE OVER THE
C                ENTIRE 5 YEAR LIFE OF THE DEGREE AND ORDER 12
C                SPHERICAL HARMONIC MODEL WMM-2005, THE ESTIMATED
C                MAXIMUM RMS ERRORS FOR THE VARIOUS MAGNETIC COMPONENTS
C                ARE:
C
C                DEC  -   0.5 Degrees
C                DIP  -   0.5 Degrees
C                TI   - 280.0 nanoTeslas (nT)
C                GV   -   0.5 Degrees
C
C                OTHER MAGNETIC COMPONENTS THAT CAN BE DERIVED FROM
C                THESE FOUR BY SIMPLE TRIGONOMETRIC RELATIONS WILL
C                HAVE THE FOLLOWING APPROXIMATE ERRORS OVER OCEAN AREAS:
C
C                X    - 140 nT (North)
C                Y    - 140 nT (East)
C                Z    - 200 nT (Vertical) Positive is down
C                H    - 200 nT (Horizontal)
C
C                OVER LAND THE MAXIMUM RMS ERRORS ARE EXPECTED TO BE
C                HIGHER, ALTHOUGH THE RMS ERRORS FOR DEC, DIP, AND GV
C                ARE STILL ESTIMATED TO BE LESS THAN 1.0 DEGREE, FOR
C                THE ENTIRE 5-YEAR LIFE OF THE MODEL AT THE EARTH's
C                SURFACE.  THE OTHER COMPONENT ERRORS OVER LAND ARE
C                MORE DIFFICULT TO ESTIMATE AND SO ARE NOT GIVEN.
C
C                THE ACCURACY AT ANY GIVEN TIME FOR ALL OF THESE
C                GEOMAGNETIC PARAMETERS DEPENDS ON THE GEOMAGNETIC
C                LATITUDE.  THE ERRORS ARE LEAST FROM THE EQUATOR TO
C                MID-LATITUDES AND GREATEST NEAR THE MAGNETIC POLES.
C
C                IT IS VERY IMPORTANT TO NOTE THAT A DEGREE AND
C                ORDER 12 MODEL, SUCH AS WMM-2005, DESCRIBES ONLY
C                THE LONG WAVELENGTH SPATIAL MAGNETIC FLUCTUATIONS
C                DUE TO EARTH'S CORE.  NOT INCLUDED IN THE WMM SERIES
C                MODELS ARE INTERMEDIATE AND SHORT WAVELENGTH
C                SPATIAL FLUCTUATIONS OF THE GEOMAGNETIC FIELD
C                WHICH ORIGINATE IN THE EARTH'S MANTLE AND CRUST.
C                CONSEQUENTLY, ISOLATED ANGULAR ERRORS AT VARIOUS
C                POSITIONS ON THE SURFACE (PRIMARILY OVER LAND, IN
C                CONTINENTAL MARGINS AND OVER OCEANIC SEAMOUNTS,
C                RIDGES AND TRENCHES) OF SEVERAL DEGREES MAY BE
C                EXPECTED. ALSO NOT INCLUDED IN THE MODEL ARE
C                NONSECULAR TEMPORAL FLUCTUATIONS OF THE GEOMAGNETIC
C                FIELD OF MAGNETOSPHERIC AND IONOSPHERIC ORIGIN.
C                DURING MAGNETIC STORMS, TEMPORAL FLUCTUATIONS CAN
C                CAUSE SUBSTANTIAL DEVIATIONS OF THE GEOMAGNETIC
C                FIELD FROM MODEL VALUES.  IN ARCTIC AND ANTARCTIC
C                REGIONS, AS WELL AS IN EQUATORIAL REGIONS, DEVIATIONS
C                FROM MODEL VALUES ARE BOTH FREQUENT AND PERSISTENT.
C
C                IF THE REQUIRED DECLINATION ACCURACY IS MORE
C                STRINGENT THAN THE WMM SERIES OF MODELS PROVIDE, THEN
C                THE USER IS ADVISED TO REQUEST SPECIAL (REGIONAL OR
C                LOCAL) SURVEYS BE PERFORMED AND MODELS PREPARED.
C                REQUESTS OF THIS NATURE SHOULD BE MADE TO NIMA
C                AT THE ADDRESS ABOVE.
C
C***********************************************************************
C
C     USAGE:  THIS ROUTINE IS BROKEN UP INTO TWO PARTS:
C
C             A) AN INITIALIZATION MODULE, WHICH IS CALLED ONLY
C                ONCE AT THE BEGINNING OF THE MAIN (CALLING)
C                PROGRAM
C             B) A PROCESSING MODULE, WHICH COMPUTES THE MAGNETIC
C                FIELD PARAMETERS FOR EACH SPECIFIED GEODETIC
C                POSITION (ALTITUDE, LATITUDE, LONGITUDE) AND TIME
C
C             INITIALIZATION IS MADE VIA A SINGLE CALL TO THE MAIN
C             ENTRY POINT (GEOMAG), WHILE SUBSEQUENT PROCESSING
C             CALLS ARE MADE THROUGH THE SECOND ENTRY POINT (GEOMG1).
C             ONE CALL TO THE PROCESSING MODULE IS REQUIRED FOR EACH
C             POSITION AND TIME.
C
C             THE VARIABLE MAXDEG IN THE INITIALIZATION CALL IS THE
C             MAXIMUM DEGREE TO WHICH THE SPHERICAL HARMONIC MODEL
C             IS TO BE COMPUTED.  IT MUST BE SPECIFIED BY THE USER
C             IN THE CALLING ROUTINE.  NORMALLY IT IS 12 BUT IT MAY
C             BE SET LESS THAN 12 TO INCREASE COMPUTATIONAL SPEED AT
C             THE EXPENSE OF REDUCED ACCURACY.
C
C             THE PC VERSION OF THIS SUBROUTINE MUST BE COMPILED
C             WITH A FORTRAN 77 COMPATIBLE COMPILER SUCH AS THE
C             MICROSOFT OPTIMIZING FORTRAN COMPILER VERSION 4.1
C             OR LATER.
C
C**********************************************************************
C
C     REFERENCES:
C
C       JOHN M. QUINN, DAVID J. KERRIDGE AND DAVID R. BARRACLOUGH,
C            WORLD MAGNETIC CHARTS FOR 1985 - SPHERICAL HARMONIC
C            MODELS OF THE GEOMAGNETIC FIELD AND ITS SECULAR
C            VARIATION, GEOPHYS. J. R. ASTR. SOC. (1986) 87,
C            PP 1143-1157
C
C       DEFENSE MAPPING AGENCY TECHNICAL REPORT, TR 8350.2:
C            DEPARTMENT OF DEFENSE WORLD GEODETIC SYSTEM 1984,
C            SEPT. 30 (1987)
C
C       JOHN M. QUINN, RACHEL J. COLEMAN, MICHAEL R. PECK, AND
C            STEPHEN E. LAUBER; THE JOINT US/UK 1990 EPOCH
C            WORLD MAGNETIC MODEL, TECHNICAL REPORT NO. 304,
C            NAVAL OCEANOGRAPHIC OFFICE (1991)
C
C       JOHN M. QUINN, RACHEL J. COLEMAN, DONALD L. SHIEL, AND
C            JOHN M. NIGRO; THE JOINT US/UK 1995 EPOCH WORLD
C            MAGNETIC MODEL, TECHNICAL REPORT NO. 314, NAVAL
C            OCEANOGRAPHIC OFFICE (1995)
C
C            SUSAN AMCMILLAN, DAVID R. BARRACLOUGH, JOHN M. QUINN, AND
C            RACHEL J. COLEMAN;  THE 1995 REVISION OF THE JOINT US/UK
C            GEOMAGNETIC FIELD MODELS - I. SECULAR VARIATION, JOURNAL OF
C            GEOMAGNETISM AND GEOELECTRICITY, VOL. 49, PP. 229-243
C            (1997)
C
C            JOHN M. QUINN, RACHEL J. COELMAN, SUSAM MACMILLAN, AND
C            DAVID R. BARRACLOUGH;  THE 1995 REVISION OF THE JOINT
C            US/UK GEOMAGNETIC FIELD MODELS: II. MAIN FIELD,JOURNAL OF
C            GEOMAGNETISM AND GEOELECTRICITY, VOL. 49, PP. 245 - 261
C            (1997)
C
C***********************************************************************
C
C     PARAMETER DESCRIPTIONS:
C
C       A      - SEMIMAJOR AXIS OF WGS-84 ELLIPSOID (KM)
C       B      - SEMIMINOR AXIS OF WGS-84 ELLIPSOID (KM)
C       RE     - MEAN RADIUS OF IAU-66 ELLIPSOID (KM)
C       SNORM  - SCHMIDT NORMALIZATION FACTORS
C       C      - GAUSS COEFFICIENTS OF MAIN GEOMAGNETIC MODEL (NT)
C       CD     - GAUSS COEFFICIENTS OF SECULAR GEOMAGNETIC MODEL (NT/YR)
C       TC     - TIME ADJUSTED GEOMAGNETIC GAUSS COEFFICIENTS (NT)
C       OTIME  - TIME ON PREVIOUS CALL TO GEOMAG (YRS)
C       OALT   - GEODETIC ALTITUDE ON PREVIOUS CALL TO GEOMAG (YRS)
C       OLAT   - GEODETIC LATITUDE ON PREVIOUS CALL TO GEOMAG (DEG.)
C       TIME   - COMPUTATION TIME (YRS)                        (INPUT)
C                (EG. 1 JULY 1995 = 1995.500)
C       ALT    - GEODETIC ALTITUDE (KM)                        (INPUT)
C       GLAT   - GEODETIC LATITUDE (DEG.)                      (INPUT)
C       GLON   - GEODETIC LONGITUDE (DEG.)                     (INPUT)
C       EPOCH  - BASE TIME OF GEOMAGNETIC MODEL (YRS)
C       DTR    - DEGREE TO RADIAN CONVERSION
C       SP(M)  - SINE OF (M*SPHERICAL COORD. LONGITUDE)
C       CP(M)  - COSINE OF (M*SPHERICAL COORD. LONGITUDE)
C       ST     - SINE OF (SPHERICAL COORD. LATITUDE)
C       CT     - COSINE OF (SPHERICAL COORD. LATITUDE)
C       R      - SPHERICAL COORDINATE RADIAL POSITION (KM)
C       CA     - COSINE OF SPHERICAL TO GEODETIC VECTOR ROTATION ANGLE
C       SA     - SINE OF SPHERICAL TO GEODETIC VECTOR ROTATION ANGLE
C       BR     - RADIAL COMPONENT OF GEOMAGNETIC FIELD (NT)
C       BT     - THETA COMPONENT OF GEOMAGNETIC FIELD (NT)
C       BP     - PHI COMPONENT OF GEOMAGNETIC FIELD (NT)
C       P(N,M) - ASSOCIATED LEGENDRE POLYNOMIALS (UNNORMALIZED)
C       PP(N)  - ASSOCIATED LEGENDRE POLYNOMIALS FOR M=1 (UNNORMALIZED)
C       DP(N,M)- THETA DERIVATIVE OF P(N,M) (UNNORMALIZED)
C       BX     - NORTH GEOMAGNETIC COMPONENT (NT)
C       BY     - EAST GEOMAGNETIC COMPONENT (NT)
C       BZ     - VERTICALLY DOWN GEOMAGNETIC COMPONENT (NT)
C       BH     - HORIZONTAL GEOMAGNETIC COMPONENT (NT)
C       DEC    - GEOMAGNETIC DECLINATION (DEG.)                (OUTPUT)
C                  EAST=POSITIVE ANGLES
C                  WEST=NEGATIVE ANGLES
C       DIP    - GEOMAGNETIC INCLINATION (DEG.)                (OUTPUT)
C                  DOWN=POSITIVE ANGLES
C                    UP=NEGATIVE ANGLES
C       TI     - GEOMAGNETIC TOTAL INTENSITY (NT)              (OUTPUT)
C       GV     - GEOMAGNETIC GRID VARIATION (DEG.)             (OUTPUT)
C                REFERENCED TO GRID NORTH
C                GRID NORTH REFERENCED TO 0 MERIDIAN
C                OF A POLAR STEREOGRAPHIC PROJECTION
C                (ARCTIC/ANTARCTIC ONLY)
C       MAXDEG - MAXIMUM DEGREE OF SPHERICAL HARMONIC MODEL    (INPUT)
C       MOXORD - MAXIMUM ORDER OF SPHERICAL HARMONIC MODEL
C
C***********************************************************************
C
C     NOTE:  THIS VERSION OF GEOMAG USES A WMM SERIES GEOMAGNETIC
C            FIELS MODEL REFERENCED TO THE WGS-84 GRAVITY MODEL
C            ELLIPSOID
C
C***********************************************************************
C
C
C                        INITIALIZATION MODULE
C
C
C***********************************************************************
C
C
      SUBROUTINE GEOMAG(MAXDEG,EPOCH,WMMCOF)
C
      IMPLICIT NONE
C
      INTEGER MAXDEG
      REAL EPOCH
      CHARACTER*200 WMMCOF

C
      INTEGER N, M, MAXORD, J
      REAL PI, DTR, A, B, RE, A2, B2, C2, A4, B4, C4
      REAL GNM, HNM, DGNM, DHNM, FLNMJ
      REAL OTIME, OALT, OLAT
C
      REAL C(0:12,0:12), CD(0:12,0:12)
      REAL P(0:12,0:12), DP(0:12,0:12), SNORM(0:12,0:12)
      REAL SP(0:12), CP(0:12), FN(0:12), FM(0:12), PP(0:12)
      REAL K(0:12,0:12)
C
      COMMON /INIT/ MAXORD, PI, DTR, A, B, RE, A2, B2, C2, A4, 
     1     B4, C4, OTIME, OALT, OLAT,
     2     C, CD, P, DP, SP, CP, FN, FM, PP, K
C
      CHARACTER*20 MODEL
C      EQUIVALENCE (SNORM,P)
      !OPEN (UNIT=9,FILE='data/magneticfield/WMM2020/WMM.COF')
      !WRITE(*,'(A,A)') 'File = ', WMMCOF
      OPEN (UNIT=9,FILE=WMMCOF)
C
C
C
C
C
C
C        INITIALIZE CONSTANTS
C
C
C
      READ(9,1) EPOCH,MODEL
    1 FORMAT(4X,F6.1,5X,A)
      IF (MAXDEG .GT. 12) MAXDEG=12
      MAXORD=MAXDEG
      PI=3.14159265359
      DTR=PI/180.0
      SP(0)=0.
      CP(0)=1.
      P(0,0)=1.
      PP(0)=1.
      DP(0,0)=0.
      A=6378.137
      B=6356.7523142
      RE=6371.2
      A2=A**2
      B2=B**2
      C2=A2-B2
      A4=A2**2
      B4=B2**2
      C4=A4-B4
C
C      READ WORLD MAGNETIC MODEL SPHERICAL HARMONIC COEFFICIENTS
C
      C(0,0)=0.0
      CD(0,0)=0.0
    3 CONTINUE
      READ(9,4) N,M,GNM,HNM,DGNM,DHNM
    4 FORMAT(2(1X,I2),2(1X,F9.1),2(1X,F10.1))
      IF (N .EQ. 99) GO TO 5
      IF (M .LE. N) THEN
      C(N,M)=GNM
      CD(N,M)=DGNM
      IF (M .NE. 0) THEN
      C(M-1,N)=HNM
      CD(M-1,N)=DHNM
      ENDIF
      ENDIF
      GO TO 3
    5 CONTINUE
      CLOSE (UNIT=9)
C
C        CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED
C
      SNORM(0,0)=1.
      DO 20 N=1,MAXORD
      SNORM(N,0)=SNORM(N-1,0)*FLOAT(2*N-1)/FLOAT(N)
      J=2
      DO 10 M=0,N
      K(N,M)=FLOAT((N-1)**2-M**2)/FLOAT((2*N-1)*(2*N-3))
      IF (M .GT. 0) THEN
      FLNMJ=FLOAT((N-M+1)*J)/FLOAT(N+M)
      SNORM(N,M)=SNORM(N,M-1)*SQRT(FLNMJ)
      J=1
      C(M-1,N)=SNORM(N,M)*C(M-1,N)
      CD(M-1,N)=SNORM(N,M)*CD(M-1,N)
      ENDIF
      C(N,M)=SNORM(N,M)*C(N,M)
      CD(N,M)=SNORM(N,M)*CD(N,M)
   10 CONTINUE
      FN(N)=FLOAT(N+1)
      FM(N)=FLOAT(N)
   20 CONTINUE
      K(1,1)=0.
C
C
      OTIME=-1000.
      OALT=-1000.
      OLAT=-1000.
C
C
      RETURN
      END
C
C
C***********************************************************************
C
C
C                        PROCESSING MODULE
C
C
C***********************************************************************
C
C
C      SUBROUTINE GEOMG1(ALT,GLAT,GLON,TIME,DEC,DIP,TI,GV,EPOCH)
      SUBROUTINE GEOMG1(ALT,GLAT,GLON,TIME,BX,BY,BZ,GV,EPOCH)
C
      IMPLICIT NONE
C
      INTEGER N, M, MAXORD
      REAL ALT, GLAT, GLON, TIME, DEC, DIP, TI, GV, EPOCH
      REAL DT, RLON, RLAT, SRLON, SRLAT, CRLON, CRLAT
      REAL SRLON2, SRLAT2, CRLON2, CRLAT2
      REAL Q, Q1, Q2, CT, ST, R2, R, D, CA, SA, AOR, AR
      REAL BR, BT, BP, BPP, PAR, TEMP1, TEMP2, PARP
      REAL BX, BY, BZ, BH
C
      REAL PI, DTR, A, B, RE, A2, B2, C2, A4, B4, C4
      REAL OTIME, OALT, OLAT
c
      REAL C(0:12,0:12), CD(0:12,0:12), TC(0:12,0:12)
      REAL P(0:12,0:12), DP(0:12,0:12), SNORM(0:12,0:12)
      REAL SP(0:12), CP(0:12), FN(0:12), FM(0:12), PP(0:12)
      REAL K(0:12,0:12)

      COMMON /INIT/ MAXORD, PI, DTR, A, B, RE, A2, B2, C2, A4, 
     1     B4, C4, OTIME, OALT, OLAT,
     2     C, CD, P, DP, SP, CP, FN, FM, PP, K

C
      DT = TIME - EPOCH
      RLON=GLON*DTR
      RLAT=GLAT*DTR
      SRLON=SIN(RLON)
      SRLAT=SIN(RLAT)
      CRLON=COS(RLON)
      CRLAT=COS(RLAT)
      SRLON2=SRLON**2
      SRLAT2=SRLAT**2
      CRLON2=CRLON**2
      CRLAT2=CRLAT**2
      SP(1)=SRLON
      CP(1)=CRLON
      
      OTIME=-1000.
      OALT=-1000.
      OLAT=-1000.
C
C        CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS.
C
      Q=SQRT(A2-C2*SRLAT2)
      Q1=ALT*Q
      Q2=((Q1+A2)/(Q1+B2))**2
      CT=SRLAT/SQRT(Q2*CRLAT2+SRLAT2)
      ST=SQRT(1.0-CT**2)
      R2=ALT**2+2.0*Q1+(A4-C4*SRLAT2)/Q**2
      R=SQRT(R2)
      D=SQRT(A2*CRLAT2+B2*SRLAT2)
      CA=(ALT+D)/R
      SA=C2*CRLAT*SRLAT/(R*D)
C
C
      DO 40 M=2,MAXORD
      SP(M)=SP(1)*CP(M-1)+CP(1)*SP(M-1)
      CP(M)=CP(1)*CP(M-1)-SP(1)*SP(M-1)
   40 CONTINUE
C
C
      AOR=RE/R
      AR=AOR**2
C
C
      BR=0.
      BT=0.
      BP=0.
      BPP=0.
      
      !WRITE(*,'(A,F8.3)') 'OTIME = ', OTIME
      !WRITE(*,'(A,F8.3)') 'ALT = ', ALT
      !WRITE(*,'(A,F8.3)') 'GLAT = ', GLAT
      !WRITE(*,'(A,F8.3)') 'GLON = ', GLON
      !WRITE(*,'(A,F8.3)') 'TIME = ', TIME
      !WRITE(*,'(A,F8.3)') 'EPOCH = ', EPOCH
      !WRITE(*,'(A,F8.3)') 'CD(10,3) = ', CD(10,3)
C
C
      DO 70 N=1,MAXORD
         AR=AR*AOR
         DO 60 M=0,N
C     
C     COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
C     AND DERIVATIVES VIA RECURSION RELATIONS
C     
            IF (ALT .NE. OALT .OR. GLAT .NE. OLAT) THEN
               IF (N .EQ. M) THEN
                  P(N,M)=ST*P(N-1,M-1)
                  DP(N,M)=ST*DP(N-1,M-1)+CT*P(N-1,M-1)
                  GO TO 50
               ENDIF
               IF (N .EQ. 1 .AND. M .EQ. 0) THEN
                  P(N,M)=CT*P(N-1,M)
                  DP(N,M)=CT*DP(N-1,M)-ST*P(N-1,M)
                  GO TO 50
               ENDIF
               IF (N .GT. 1 .AND. N .NE. M) THEN
                  IF (M .GT. N-2) P(N-2,M)=0.0
                  IF (M .GT. N-2) DP(N-2,M)=0.0
                  P(N,M)=CT*P(N-1,M)-K(N,M)*P(N-2,M)
                  DP(N,M)=CT*DP(N-1,M)-ST*P(N-1,M)-K(N,M)*DP(N-2,M)
               ENDIF
            ENDIF
 50         CONTINUE
C     
C        TIME ADJUST THE GAUSS COEFFICIENTS
C
            IF (TIME .NE. OTIME) THEN
               TC(N,M)=C(N,M)+DT*CD(N,M)
               IF (M .NE. 0) THEN
                  TC(M-1,N)=C(M-1,N)+DT*CD(M-1,N)
               ENDIF
            ENDIF
C     
C        ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
C
            PAR=AR*P(N,M)
            IF (M .EQ. 0) THEN
               TEMP1=TC(N,M)*CP(M)
               TEMP2=TC(N,M)*SP(M)
            ELSE
               TEMP1=TC(N,M)*CP(M)+TC(M-1,N)*SP(M)
               TEMP2=TC(N,M)*SP(M)-TC(M-1,N)*CP(M)
            ENDIF
            BT=BT-AR*TEMP1*DP(N,M)
            BP=BP+FM(M)*TEMP2*PAR
            BR=BR+FN(N)*TEMP1*PAR
C     
C     SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
C     
            IF (ST .EQ. 0.0 .AND. M .EQ. 1) THEN
               IF (N .EQ. 1) THEN
                  PP(N)=PP(N-1)
               ELSE
                  PP(N)=CT*PP(N-1)-K(N,M)*PP(N-2)
               ENDIF
               PARP=AR*PP(N)
               BPP=BPP+FM(M)*TEMP2*PARP
            ENDIF
C     
C     
 60      CONTINUE
 70   CONTINUE
C     
C
      IF (ST .EQ. 0.0) THEN
      BP=BPP
      ELSE
      BP=BP/ST
      ENDIF
C
C        ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
C        GEODETIC COORDINATES
C
      BX=-BT*CA-BR*SA
      BY=BP
      BZ=BT*SA-BR*CA
C
C        COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
C        TOTAL INTENSITY (TI)
C
      BH=SQRT(BX**2+BY**2)
      TI=SQRT(BH**2+BZ**2)
      DEC=ATAN2(BY,BX)/DTR
      DIP=ATAN2(BZ,BH)/DTR
C
C        COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
C        GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
C        (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)
C
C        OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
C
      GV=-999.0
      IF (ABS(GLAT) .GE. 55.) THEN
      IF (GLAT .GT. 0. .AND. GLON .GE. 0.) GV=DEC-GLON
      IF (GLAT .GT. 0. .AND. GLON .LT. 0.) GV=DEC+ABS(GLON)
      IF (GLAT .LT. 0. .AND. GLON .GE. 0.) GV=DEC+GLON
      IF (GLAT .LT. 0. .AND. GLON .LT. 0.) GV=DEC-ABS(GLON)
      IF (GV .GT. +180.) GV=GV-360.
      IF (GV .LT. -180.) GV=GV+360.
      ENDIF
C
C
      OTIME=TIME
      OALT=ALT
      OLAT=GLAT
C
C
      RETURN
      END























