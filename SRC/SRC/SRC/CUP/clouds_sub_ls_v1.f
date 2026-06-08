C***************************************************************************
C*****                   SUBROUTINE CLOUDS_SUB_LS                      *****
C*****                                                                 *****
C*****                          VERSION 1.0                            *****
C*****                                                                 *****
C*****                        16 November, 2001                        *****
C*****                                                                 *****
C***************************************************************************
C
C Purpose:
C --------
C
C   This parameterization predicts the cloudiness (cloud amount, 
C   cloud water content) associated with cumulus convection. 
C
C   It is based upon the idea that the convection scheme predicts 
C   the local concentration of condensed water (the in-cloud water 
C   content) produced at the subgrid-scale, and that a statistical 
C   cloud scheme predicts how this condensed water is spatially 
C   distributed within the domain.
C
C   The statistical cloud scheme uses a generalized log-normal 
C   Probability Distribution Function (GNO PDF) of the total water
C   whose variance and skewness coefficient are diagnosed from
C 	the amount of condensed water produced at the subgrid-scale
C	by cumulus convection and at the large-scale by supersaturation,
C   from the degree of saturation of the environment, and from the
C 	lower bound of the total water distribution that is taken 
C   equal to 0.
C
C Remarks:
C --------
C
C   * this parameterization is intended to be used in association
C   with the Emanuel convection scheme CONVECT (Emanuel 1991, 
C   Emanuel and Zivkovic-Rothman 1999).
C
C   * in the present version of the scheme, the only source of 
C   subgrid-scale condensation that is considered is cumulus 
C   convection; clouds associated with other processes (e. g. 
C   boundary-layer turbulence) are not predicted by this 
C   parameterization.
C
C   * if no subgrid-scale condensation occurs within the domain
C   (i.e, for the moment, in non-convective regions), the scheme 
C   becomes equivalent to an "all-or-nothing" large-scale 
C   saturation scheme. 
C
C	* the subgrid-scale condensation being predicted by the 
C   convection scheme, this is the convection scheme that takes 
C   care of the temperature/water tendencies and precipitation 
C   associated with subgrid-scale clouds. The T/q tendencies and 
C   precipitation computed in this parameterization are related
C   to the occurence of large-scale supersaturation only.
C   On the other hand, the cloud fraction and water content in
C   output of this routine are for subgrid-scale+large-scale clouds.
C
C Authors:
C --------
C
C  Sandrine Bony (LMD/CNRS, bony@lmd.jussieu.fr) 
C  & Kerry Emanuel (MIT, emanuel@texmex.mit.edu)
C
C Reference:
C ----------
C
C  Bony, S and K A Emanuel, 2001: A parameterization of the cloudiness
C  associated with cumulus convection; Evaluation using TOGA COARE data.
C  J. Atmos. Sci., 58, No 21, 3158-3183.
C
C============================================================================

        SUBROUTINE CLOUDS_SUB_LS(ND,R,RS,T,P,PH,LV0,DT,QSUBGRID
     :                          ,CLDF,CLDQ,PRADJ,FTADJ,FRADJ,IFC)

        IMPLICIT NONE

C--------------------------------------------------------------------------------
C
C Inputs:
C
C  ND----------: Number of vertical levels
C  R--------ND-: Grid-box average of the total water mixing ratio [kg/kg]
C  RS-------ND-: Mean saturation humidity mixing ratio within the gridbox [kg/kg]
C  T--------ND-: Grid-box average temperature [K]
C  P-------ND+1: Pressure at mid-levels [mb]
C  PH------ND+1: Pressure at interface levels [mb]
C  LV0-----ND  : Specific heat of condensation/deposition (J/kg)
C  DT----------: Timestep [seconds]
C  QSUBGRID-ND-: in-cloud mixing ratio of cloud condensate [kg/kg] from CONVECT
C
C Outputs:
C
C  CLDF-----ND-: cloud fraction [0-1]
C  CLDQ-----ND-: in-cloud mixing ratio of condensed water [kg/kg]
C  IFC------ND-: flag: =1 if numerical convergence in clouds_gno, 2 otherwise
C  PRADJ----ND-: precipitation associated with the LS super-saturation [mm/day]
C  FTADJ----ND-: temperature tendency associated with the LS adjustment [K/s]
C  FRADJ----ND-: total water tendency associated with the LS adjustment [kg/kg/s]
C
C--------------------------------------------------------------------------------

      integer ND,I
      integer IFC(ND)
      real DT, PRADJ
      real TCA,ELACRIT,ALV,CPN,TNEW,RNEW,EP,ELCRITI
      real R(ND),RS(ND),T(ND),P(ND),PH(ND+1)
     :    ,QSUBGRID(ND),CLDF(ND),CLDQ(ND),LV0(ND)
     :    ,FTADJ(ND),FRADJ(ND),RNEWLS(ND),TNEWLS(ND),QLSP(ND)

C Thermodynamical constants:
C
        REAL CPD,CPV,CL,RD,RV,G,ROWL,EPS,EPSI,CPVMCL
        PARAMETER (CPD=1004.64, CPV=1870.0, CL=4190.0, RD=287.04)
        PARAMETER (RV=461.5, G=9.79764, ROWL=1000.0 )
        PARAMETER (EPS=RD/RV, EPSI=1./EPS, CPVMCL=CL-CPV)

C Microphysical parameters:
C (here, we use the same values as in the convection scheme)
C
        REAL TLCRIT,ELCRIT,EPMAX
        real om

        TLCRIT= -55.0
        ELCRIT = 0.0011 
        ELCRITI = 0.0001  ! for ice
        EPMAX = 0.999

C
C ***  Initializations   ***
C
        PRADJ=0.0
        DO I = 1, ND
         FRADJ(I)  = 0.0
         FTADJ(I)  = 0.0
         TNEWLS(I) = T(I)
         RNEWLS(I) = R(I)
         QLSP(I)   = 0.0
         IFC(I)    = 0
        ENDDO

C
C ***   Call the subgrid-cloud parameterization   ***
C

         CALL CLOUDS_GNO(ND,R,RS,QSUBGRID,CLDF,CLDQ,IFC)

C
C ***      When the atmosphere is super-saturated       ***
C ***   at the large-scale, we precipitate a fraction   ***
C ***   of the water content exceeding saturation       ***
C

        DO 9999 I = ND, 1, -1

         IF(R(I).GT.RS(I))THEN

C  Compute a LS precipitation efficiency as done in CONVECT:

         EP=0.0
         TCA=T(I)-273.15
         IF(TCA.GE.0.0)THEN
          ELACRIT=ELCRIT
         ELSE
       !   om = max(0.,min(1., -TCA/20.))
       !   ELACRIT=ELCRITI*om + ELCRIT*(1.-om)
          ELACRIT=ELCRIT*(1.0-TCA/TLCRIT)
          ELACRIT=MAX(ELACRIT,0.0)**1.05
         END IF
         EP= EPMAX * (1.0-ELACRIT/MAX(R(I)-RS(I),1.0E-8))
         EP=MAX(EP,0.0)
         EP=MIN(EP,EPMAX)

         QLSP(I) = EP*(R(I)-RS(I)) 
         CLDQ(I) = CLDQ(I) - QLSP(I)

C  Adjust temperature and humidity profiles:

          TCA=T(I)-273.15
          ALV=LV0(I)-CPVMCL*TCA
          CPN=CPD*(1.-R(I))+CPV*R(I)

          TNEW=(ALV*(EP*R(I)+RS(I)*(ALV/(RV*T(I))-EP))
     :       +CPN*T(I))/
     1     (CPN+ALV*ALV*RS(I)/(RV*T(I)*T(I)))
          RNEW=RS(I)*(1.+(TNEW-T(I))*ALV/(RV*T(I)*T(I)))
     :       + (1.-EP)*(R(I)-RS(I))

          TNEWLS(I) = TNEW
          RNEWLS(I) = RNEW

          FRADJ(I)=FRADJ(I)+(RNEW-R(I))/DT
          FTADJ(I)=FTADJ(I)-ALV*(RNEW-R(I))/DT/CPN

          PRADJ=PRADJ-100.0*(PH(I)-PH(I+1))*(RNEW-R(I))/DT
     :            *1000.0*3600.0*24.0/(ROWL*G)

         ENDIF 

9999    CONTINUE

         RETURN
         END

C
C================================================================================
C
      SUBROUTINE CLOUDS_GNO(ND,R,RS,QSUB,CLDF,CLDQ,IFC)
      IMPLICIT NONE
C     
C--------------------------------------------------------------------------------
C
C Inputs:
C
C  ND----------: Number of vertical levels
C  R--------ND-: Domain-averaged mixing ratio of total water 
C  RS-------ND-: Mean saturation humidity mixing ratio within the gridbox
C  QSUB-----ND-: Mixing ratio of condensed water within clouds associated
C                with SUBGRID-SCALE condensation processes (here, it is
C                predicted by the convection scheme)
C Outputs:
C
C  CLDF-----ND-: cloud fractional area [0-1]
C  CLDQ-----ND-: in-cloud mixing ratio of condensed water [kg/kg]
C  IFC------ND-: flag = 1 when numerical convergence, 2 otherwise.
C
C--------------------------------------------------------------------------------

      INTEGER ND
      INTEGER IFC(ND)
      REAL R(ND),  RS(ND), QSUB(ND), CLDF(ND), CLDQ(ND)

c -- parameters controlling the iteration:
c --    nmax    : maximum nb of iterations (hopefully never reached)
c --    epsilon : accuracy of the numerical resolution 
c --    vmax    : v-value above which we use an asymptotic expression for ERF(v)

      INTEGER nmax, niter
      PARAMETER ( nmax = 10) 
      REAL epsilon, vmax0, vmax
      PARAMETER ( epsilon = 0.02, vmax0 = 2.0 ) 

      REAL min_mu, min_Q
      PARAMETER ( min_mu =  1.e-12, min_Q=1.e-12 )
     
      INTEGER K, n, m
      REAL*8 mu, qsat, delta, beta 
      REAL*8 xx, aux, coeff, block, dist, fprime, det
      REAL*8 pi, u, v, erfu, erfv, xx1, xx2
      LOGICAL lconv
      real*8 derf


      pi = ACOS(-1.)


      DO 500 K = 1, ND

      mu = R(K)
      mu = MAX(mu,min_mu)
      qsat = RS(K) 
      qsat = MAX(qsat,min_mu)
      delta = log(mu/qsat)

C
C ***          There is no subgrid-scale condensation;        ***
C ***   the scheme becomes equivalent to an "all-or-nothing"  *** 
C ***             large-scale condensation scheme.            ***
C

      IF ( QSUB(K) .lt. min_Q ) THEN

        CLDQ(K) = MAX( 0.0, mu-qsat )
        CLDF(K) = CLDQ(K) / MAX( CLDQ(K), min_mu )
        IFC(K) = 1

      ELSE 

C
C ***     Some condensation is produced at the subgrid-scale       ***
C ***                                                              ***
C ***       PDF = generalized log-normal distribution (GNO)        ***
C ***   (k<0 because a lower bound is considered for the PDF)      ***
C ***                                                              ***
C ***  -> Determine x (the parameter k of the GNO PDF) such        ***
C ***  that the contribution of subgrid-scale processes to         ***
C ***  the in-cloud water content is equal to QSUB(K)              ***
C ***  (equations (13), (14), (15) + Appendix B of the paper)      ***
C ***                                                              ***
C ***    Here, an iterative method is used for this purpose        ***
C ***    (other numerical methods might be more efficient)         ***
C ***                                                              ***
C ***          NB: the "error function" is called ERF              ***
C ***                 (DERF in double precision)                   ***
C

        lconv  = .FALSE. 
        niter = 0
        vmax = vmax0

        beta = QSUB(K)/mu + EXP( -MIN(0.0,delta) )

        IF ( .NOT. lconv ) then 
        
c --  roots of equation v > vmax:

        det = delta + vmax**2.
        if (det.LE.0.0) vmax = vmax0 + 1.0
        det = delta + vmax**2.

        if (det.LE.0.) then
          xx = -0.0001
        else 
         xx1 = -SQRT(2.)*vmax*(1.0-SQRT(1.0+delta/(vmax**2.)))
         xx2 = -SQRT(2.)*vmax*(1.0+SQRT(1.0+delta/(vmax**2.)))
         xx = 1.01 * xx1
         if ( xx1 .GE. 0.0 ) xx = 0.5*xx2
        endif
        if (delta.LT.0.) xx = -0.5*SQRT(log(2.)) 

        DO n = 1, nmax 

          u = delta/(xx*sqrt(2.)) + xx/(2.*sqrt(2.))
          v = delta/(xx*sqrt(2.)) - xx/(2.*sqrt(2.))

          IF ( v .GT. vmax ) THEN 

            IF (     ABS(u)  .GT. vmax 
     :          .AND.  delta .LT. 0. ) THEN

c -- use asymptotic expression of erf for u and v large:
c ( -> analytic solution for xx )

             aux = 2.0*delta*(1.-beta*EXP(delta))
     :                       /(1.+beta*EXP(delta))
             xx = -SQRT(aux)
             block = EXP(-v*v) / v / SQRT(pi)
             dist = 0.0
             fprime = 1.0

            ELSE

c -- erfv -> 1.0, use an asymptotic expression of erfv for v large:

             erfu = DERF(u)
             aux = SQRT(pi) * (1.0-erfu) * EXP(v*v)
             coeff = 1.0 - 1./2./(v**2.) + 3./4./(v**4.)
             block = coeff * EXP(-v*v) / v / SQRT(pi)
             dist = v * aux / coeff - beta
             fprime = 2.0 / xx * (v**2.)
     :           * ( coeff*EXP(-delta) - u * aux )
     :           / coeff / coeff
            
            ENDIF ! ABS(u)

          ELSE

c -- general case:

           erfu = DERF(u)
           erfv = DERF(v)
           block = 1.0-erfv
           dist = (1.0 - erfu) / (1.0 - erfv) - beta
           fprime = 2. /SQRT(pi) /xx /(1.0-erfv)**2.
     :           * (   (1.0-erfv)*v*EXP(-u*u)
     :               - (1.0-erfu)*u*EXP(-v*v) )

          ENDIF ! x

c -- test numerical convergence:

        if ( ABS(dist/beta) .LT. epsilon ) then 

        lconv = .TRUE. 
        CLDF(K) = 0.5 * block
        CLDQ(K) = QSUB(K) + MAX(mu-qsat,0.0)
        IFC(K)  = 1

        GOTO 100

        else
           xx = xx - dist/fprime
        endif

        ENDDO ! n

 100    continue

        if (.NOT. lconv) then 

C use a all-or-nothing scheme in that (rare) case:
C (may be improved later on)

          CLDQ(K) = MAX( 0.0, mu-qsat )
          CLDF(K) = CLDQ(K) / MAX( CLDQ(K), min_mu )
          IFC(K)  = 2

        endif

        ENDIF ! lconv

      ENDIF ! qsub

500   CONTINUE  ! K

       RETURN
       END



