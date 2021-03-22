      subroutine nad_eos_dt(var,xxMass,AA,ZZ,term_var)
      include 'implno.dek'
      include 'vector_eos.dek'

! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer          nrow
      integer          ionmax
      parameter        (ionmax=3)
      parameter        (nrow=1)
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),abar,zbar
      double precision dens,temp
      double precision var(3)
      double precision xxMass(ionmax),AA(ionmax),ZZ(ionmax)
      double precision term_var(6)

! set the mass fractions, z's and a's of the composition
! hydrogen, heliu, and carbon
      !xmass(1) = 0.75d0 ; aion(1)  = 1.0d0  ; zion(1)  = 1.0d0
      !xmass(2) = 0.23d0 ; aion(2)  = 4.0d0  ; zion(2)  = 2.0d0
      !xmass(3) = 0.02d0 ; aion(3)  = 12.0d0 ; zion(3)  = 6.0d0
      xmass(1) = xxMass(1); aion(1) = AA(1); zion(1) = ZZ(1)
      xmass(2) = xxMass(2); aion(2) = AA(2); zion(2) = ZZ(2)
      xmass(3) = xxMass(3); aion(3) = AA(3); zion(3) = ZZ(3)

! average atomic weight and charge
      abar   = 1.0d0/sum(xmass(1:ionmax)/aion(1:ionmax))
      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))

      temp_row(1) = var(2); den_row(1)  = var(1); abar_row(1) = abar ; zbar_row(1) = zbar
      jlo_eos = 1 ; jhi_eos = 1

! call the eos
      call nados(term_var)

      end   




      subroutine nad_eos_dp(var,xxMass,AA,ZZ,term_var)
      include 'implno.dek'
      include 'vector_eos.dek'

! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer          nrow
      integer          ionmax
      parameter        (ionmax=3)
      parameter        (nrow=1)
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),abar,zbar
      double precision var(3)
      double precision xxMass(ionmax),AA(ionmax),ZZ(ionmax)
      double precision term_var(6)
      double precision D(nrow),P(nrow),E(nrow),tguess(nrow)
      double precision A(nrow),Z(nrow)

! set the mass fractions, z's and a's of the composition
! hydrogen, heliu, and carbon
      !xmass(1) = 0.75d0 ; aion(1)  = 1.0d0  ; zion(1)  = 1.0d0
      !xmass(2) = 0.23d0 ; aion(2)  = 4.0d0  ; zion(2)  = 2.0d0
      !xmass(3) = 0.02d0 ; aion(3)  = 12.0d0 ; zion(3)  = 6.0d0
      xmass(1) = xxMass(1); aion(1) = AA(1); zion(1) = ZZ(1)
      xmass(2) = xxMass(2); aion(2) = AA(2); zion(2) = ZZ(2)
      xmass(3) = xxMass(3); aion(3) = AA(3); zion(3) = ZZ(3)

! average atomic weight and charge
      abar   = 1.0d0/sum(xmass(1:ionmax)/aion(1:ionmax))
      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))

      D(1) = var(1)
      P(1) = var(2)
      A(1) = abar
      Z(1) = zbar
      tguess(1) = 1.0d7
      call call_helmeos_DP(nrow,D,P,A,Z,tguess,term_var)

      end   




      subroutine nad_eos_de(var,xxMass,AA,ZZ,term_var)
      include 'implno.dek'
      include 'vector_eos.dek'

! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer          nrow
      integer          ionmax
      parameter        (ionmax=3)
      parameter        (nrow=1)
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),abar,zbar
      double precision var(3)
      double precision xxMass(ionmax),AA(ionmax),ZZ(ionmax)
      double precision term_var(6)
      double precision D(nrow),P(nrow),E(nrow),tguess(nrow)
      double precision A(nrow),Z(nrow)

! set the mass fractions, z's and a's of the composition
! hydrogen, heliu, and carbon
      !xmass(1) = 0.75d0 ; aion(1)  = 1.0d0  ; zion(1)  = 1.0d0
      !xmass(2) = 0.23d0 ; aion(2)  = 4.0d0  ; zion(2)  = 2.0d0
      !xmass(3) = 0.02d0 ; aion(3)  = 12.0d0 ; zion(3)  = 6.0d0
      xmass(1) = xxMass(1); aion(1) = AA(1); zion(1) = ZZ(1)
      xmass(2) = xxMass(2); aion(2) = AA(2); zion(2) = ZZ(2)
      xmass(3) = xxMass(3); aion(3) = AA(3); zion(3) = ZZ(3)

! average atomic weight and charge
      abar   = 1.0d0/sum(xmass(1:ionmax)/aion(1:ionmax))
      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))

      D(1) = var(1)
      E(1) = var(3)
      A(1) = abar
      Z(1) = zbar
      tguess(1) = 1.0d7
      call call_helmeos_DE(nrow,D,E,A,Z,tguess,term_var)

      end   




      subroutine teos(temp,dens,xxMass,AA,ZZ,pre)
      include 'implno.dek'
      include 'vector_eos.dek'

! tests the eos routine
! 
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer          nrow
      integer          ionmax
      parameter        (ionmax=3)
      parameter        (nrow=1)
      double precision term_var(6)
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),abar,zbar
      double precision xxMass(ionmax),AA(ionmax),ZZ(ionmax),temp,dens,pre
      double precision density(nrow),pres(nrow),ener(nrow),tguess(nrow),A(nrow),Z(nrow)


! set the mass fractions, z's and a's of the composition
! hydrogen, heliu, and carbon
      !xmass(1) = 0.75d0 ; aion(1)  = 1.0d0  ; zion(1)  = 1.0d0
      !xmass(2) = 0.23d0 ; aion(2)  = 4.0d0  ; zion(2)  = 2.0d0
      !xmass(3) = 0.02d0 ; aion(3)  = 12.0d0 ; zion(3)  = 6.0d0
      xmass(1) = xxMass(1); aion(1) = AA(1); zion(1) = ZZ(1)
      xmass(2) = xxMass(2); aion(2) = AA(2); zion(2) = ZZ(2)
      xmass(3) = xxMass(3); aion(3) = AA(3); zion(3) = ZZ(3)

! average atomic weight and charge
      abar   = 1.0d0/sum(xmass(1:ionmax)/aion(1:ionmax))
      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))

! set the input vector. pipeline is only 1 element long
      temp_row(1) = temp ; den_row(1)  = dens ; abar_row(1) = abar ; zbar_row(1) = zbar
      jlo_eos = 1 ; jhi_eos = 1

! call the eos
      call nados(term_var)

! write out the results
      call pretty_eos_out('nados:  ')

      density(1) = 1.0d7
      tguess(1) = 1.0d7
      pres(1)   = 2.68840907d24
      A(1) = abar
      Z(1) = zbar
      call call_helmeos_DP(nrow, density, pres, A, Z, tguess, term_var)

      density(1) = 1.0d7
      tguess(1) = 1.0d7
      ener(1)   = 5.07259918d24
      A(1) = abar
      Z(1) = zbar
      call call_helmeos_DE(nrow, density, ener, A, Z, tguess, term_var)

      end   




! here is the tabular helmholtz free energy eos:
!
! routine read_helm_table reads an electron helm free energy table
! routine helmeos computes the pressure, energy and entropy via tables


      subroutine read_helm_table
      include 'implno.dek'
      include 'helm_table_storage.dek'

! this routine reads the helmholtz eos file, and
! must be called once before the helmeos routine is invoked.

! declare local variables
      integer          i,j
      double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                       dd,dd2,ddi,dd2i,dd3i


! open the file (use softlinks to input the desired table)

       open(unit=19,file='helm_table.dat',status='old')


! for standard table limits
       tlo   = 3.0d0
       thi   = 13.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -12.0d0
       dhi   = 15.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp

! read the helmholtz free energy and its derivatives
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                  fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo
!       write(6,*) 'read main table'


! read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo
!       write(6,*) 'read dpdd table'

! read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo
!       write(6,*) 'read eta table'

! read the number density table
       do j=1,jmax
        do i=1,imax
         read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo
!       write(6,*) 'read xne table'

! close the file
      close(unit=19)


! construct the temperature and density deltas and their inverses
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt3i        = dt2i*dti
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
        dt3i_sav(j) = dt3i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd3i        = dd2i*ddi
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
        dd3i_sav(i) = dd3i
       enddo



!      write(6,*)
!      write(6,*) 'finished reading eos table'
!      write(6,04) 'imax=',imax,' jmax=',jmax
!04    format(1x,4(a,i4))
!      write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
!      write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
!03    format(1x,4(a,1pe11.3))
!      write(6,*)

      return
      end






      subroutine read_helm_iontable
      include 'implno.dek'
      include 'helm_table_storage.dek'

! this routine reads the helmholtz eos file, and
! must be called once before the helmeos routine is invoked.

! declare local variables
      integer          i,j
      double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                       dd,dd2,ddi,dd2i,dd3i


! open the file (use softlinks to input the desired table)

       open(unit=19,file='helm_iontable.dat',status='old')


! for the standard table
       tion_lo   = 3.0d0
       tion_hi   = 13.0d0
       tion_stp  = (thi - tlo)/float(jmax-1)
       tion_stpi = 1.0d0/tstp
       dion_lo   = -12.0d0
       dion_hi   = 15.0d0
       dion_stp  = (dhi - dlo)/float(imax-1)
       dion_stpi = 1.0d0/dstp

! read the helmholtz free energy and its derivatives
       do j=1,jmax
        tsav = tion_lo + (j-1)*tion_stp
        tion(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dion_lo + (i-1)*dion_stp
         dion(i) = 10.0d0**(dsav)
         read(19,*) fion(i,j),fiond(i,j),fiont(i,j),fiondd(i,j), &
                    fiontt(i,j),fiondt(i,j),fionddt(i,j),fiondtt(i,j), &
                    fionddtt(i,j)
        enddo
       enddo


! read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(19,*) dpiondf(i,j),dpiondfd(i,j), &
                    dpiondft(i,j),dpiondfdt(i,j)
        enddo
       enddo

! read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(19,*) efion(i,j),efiond(i,j),efiont(i,j),efiondt(i,j)
        enddo
       enddo

! read the number density table
       do j=1,jmax
        do i=1,imax
         read(19,*) xfion(i,j),xfiond(i,j),xfiont(i,j),xfiondt(i,j)
        enddo
       enddo

! close the file
      close(unit=19)


! construct the temperature and density deltas and their inverses
       do j=1,jmax-1
        dth             = t(j+1) - t(j)
        dt2             = dth * dth
        dti             = 1.0d0/dth
        dt2i            = 1.0d0/dt2
        dt3i            = dt2i*dti
        dt_sav_ion(j)   = dth
        dt2_sav_ion(j)  = dt2
        dti_sav_ion(j)  = dti
        dt2i_sav_ion(j) = dt2i
        dt3i_sav_ion(j) = dt3i
       end do
       do i=1,imax-1
        dd              = d(i+1) - d(i)
        dd2             = dd * dd
        ddi             = 1.0d0/dd
        dd2i            = 1.0d0/dd2
        dd3i            = dd2i*ddi
        dd_sav_ion(i)   = dd
        dd2_sav_ion(i)  = dd2
        ddi_sav_ion(i)  = ddi
        dd2i_sav_ion(i) = dd2i
        dd3i_sav_ion(i) = dd3i
       enddo


!      write(6,*)
!      write(6,*) 'finished reading eos ion table'
!      write(6,04) 'imax=',imax,' jmax=',jmax
!04    format(1x,4(a,i4))
!      write(6,03) 'temp(1)     =',tion(1),' temp(jmax)     =',tion(jmax)
!      write(6,03) 'ytot*den(1) =',dion(1),' ytot*den(imax) =',dion(imax)
!03    format(1x,4(a,1pe11.3))
!      write(6,*)

      return
      end




! this file contains nados's equation of state
!
! for density, temperature and composition as primary thermo variables:
! routine nados is a units interface routine
! routine epeos returns the various thermodynamic quantities
!
!
!
      subroutine nados(term_var)


!          ***  version 1.1 santa cruz, august 2, 1992  ***
!***********************************************************************
!  *** equation of state for completely ionized matter
!  *** electron & positron component --- fermi-dirac statistics using
!                       various asymptotic expansions where possible
!  *** ion component --- a perfect gas approximation
!  *** black-body radiation
!***********************************************************************
!                            references
!   1. nadyozhin d.k. 1974, "naucnye informatsii", nos. 32, 33
!                            (ucsc science library: qb 1 a4)
!***********************************************************************
      implicit double precision (a-h,o-z)
      save
      include 'vector_eos.dek'


      double precision term_var(6)
      double precision chit,chid,cv,cp,gam1,gam2,gam3,nabad,sound,zz,z, &
                       prad,pion,dpiondd,dpiondt,erad,eion,srad,sion
      double precision clight,con2,me,mecc,kerg,avo
      parameter        (clight = 2.997925d10, &
                        con2   = clight * clight, &
                        me     = 9.1093897d-28, &
                        mecc   = me * con2, &
                        kerg   = 1.380658d-16, &
                        avo    = 6.0221367d23)


!
!  *** the arguments

!      double precision t,den,psi,pl
      common/arg/t,den,psi
      equivalence(den,pl)
!      integer     lst,kentr,kpar,jurs,jkk
      common /jarg/ lst,kentr,kpar,jurs,jkk

!***********************************************************************
!
!  ***   t --- temperature in 10**9 k
!  *** den --- density in 10**7 g/ccm
!  *** psi --- parameter of degeneracy. works as an argument only when
!              one enters entry fd12f to get fermi-dirac functions,
!              otherwise it is calculated as a function of t and den.
!  *** lst, kentr, kpar --- the keys aimed to make the calculations
!               faster where possible (default: lst, kentr, kpar = 0)
!  *** lst=0 --- epeos calculates only thermodynamic functions p,e,s
!      lst=1 --- the first temperature-derivatives are calculated
!                in addition to the thermodynamic functions
!      lst=2 --- extends calculations to get density-derivatives as well
!  *** kentr=0 --- turns the calculation of entropy off and suspends
!            the calculation of psi in a perfect gas asymptotics (nz=1).
!      kentr=1 --- turns the calculation of entropy on.
!  *** kpar=0 --- when in relativistic asymptotics (nz=4), turns off the
!        calculation of total number of pair (hpr),(kpar=1 --- turns on)
!      kpar is inactive for other asymptotics.
!
!  *** jkk --- the current number of mesh point, is inactive in this
!        version of epeos. however, it appears in epeos error messages
!        and, if defined, may be useful for locating of errors in
!        a program calling epeos.
!***********************************************************************
      common/nz/nz
!***********************************************************************
!  *** nz --- specifies the operational mode (default: nz=0)
!      nz=0 --- calculations with the overall search for
!               the appropriate working region
!      for 0<nz<6 epeos works within one of five following modes
!      independent of the values of temperature and density specified
!      nz=1 --- perfect gas approximation with the first order
!               corrections for degeneracy and pairs
!      nz=2 --- expansion over half-integer fermi--dirac functions
!      nz=3 --- chandrasekhar's expansion for degenerate gas
!      nz=4 --- relativistic asymptotics
!      nz=5 --- quadratures taken with the gauss method
!***********************************************************************
      common/az/as,zs,scn
!***********************************************************************
!  ***  as --- mean mass number, zs --- mean nuclear charge
!          emue=as/zs --- mean electron molecular weight: the total
!          number of "atomic" electrons in a unit volume, nea,
!          amounts to (density)/(mu*emue), mu is atomic mass unit.
!          for a mixture:   as=1/sum{y_i},  zs=as*sum{z_i*y_i}, where
!          y_i=x_i/a_i,  and  x_i  being a mass fraction of 'i' species.
!  *** scn --- additive entropy constant for the ion component
!          for a gas of nuclei with mass number a, and spin i
!          scn=ln[(2i+1)*a**2.5], for iron-56, scn=2.5*ln(56)=10.063379.
!          for a mixture:   scn=as*sum{y_i*ln[(2i_i+1)*(a_i)**2.5}.
!  *** jurs --- if jurs=0 then common-block  /az/ is to be preliminary
!            filled with all necessary information. (default values
!            of as,zs,scn are specified in block data for pure iron-56).
! if jurs=1, epeos applies to subroutine chemic for as,zs,scn
!***********************************************************************
!
!                         diagnostics
!***********************************************************************
! epeos opens file 'epeos.mes', writes the epeos version in it, analyses
! the values of arguments in common-blocks /arg/, /jarg/, /nz/, /az/ and
! then writes additional information in 'epeos.mes' if something wrong
! with the arguments: in particular, epeos stops when  t = or < 0.
! in case  den < 0 , epeos sends a warning in 'epeos.mes' and continues
! to calculate with den=0 (black body radiation only).
!***********************************************************************
!
!  *** the results of calculations
      common/result/p,e,s,sk,pt,et,st
      common/str/ppl,epl,spl,cp,gam,da,dpe,dse,dsp,beta
      common/resel/pe,ee,se,sek,hpr
      common/nzr/nzr
!***********************************************************************
!  *** p --- total pressure in units 10**24 erg/ccm
!  *** e --- total specific energy in units 10**17 erg/g
!  *** s --- total entropy in units  10**8 erg/(g k)
!  *** sk --- total dimensionless entropy per nucleon: sk=s*mu/kb
!             mu -- atomic mass unit, kb -- boltzmann's constant
!  *** pt,et,st --- temperature derivatives of p,e,s at constant density
!  *** ppl,spl --- density derivatives of p and s at constant temperature
!  *** epl---- density derivatives of e multiplied by density den
!  *** pe,ee,se,sek ----  pressure, specific energy and entropy
!                         of the electron-positron gas component
!  *** hpr --- the total number of the electron-positron pairs
!              per "atomic" electron (hpr=mu*emue*np/den), where
!              'np' is the number of pairs per unit volume.
!  *** cp --- specific heat at constant pressure
!  *** gam --- adiabatic index = {d log(p)/d log(den)} at s=const
!  *** da --- logarithmic adiabatic temperature gradient
!  *** dpe = (den/p)(epl+t*pt/den)-1 -- thermodynamic identity: dpe=0
!  *** dse = t*st/et-1 -- thermodynamic identity: dse=0
!  *** dsp = -spl*(den/pt)-1 -- thermodynamic identity: dsp=0
!  *** beta --- ratio of black-body radiation pressure to total pressure
!  *** nzr --- identificator of working region on t-den plane when nz=0
!***********************************************************************
      common/fdf/f12,f32,f52,f72,f12s,f32s,f52s,f72s
!***********************************************************************
!  *** f12,f32,f52,f72 --- half-integer fermi-dirac functions
!  *** f12s,f32s,f52s,f72s --- the first derivatives of f-d functions
!  *** psi (in common/arg/t,den,psi) is the argument of f-d functions
!      there exists special entry to get these functions separately --
!      use command call fd12f after specifying psi in common/arg/
!***********************************************************************
      dimension cu(62),ck(3),a(17),c(8),b(22),b1(4),b2(4),b3(4),b4(4), &
      b5(4),c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),d1(4),d2(4),d3(4),d4(4), &
      d5(4),d6(4),d(4),a1(4),a2(4),a3(4),a4(4),df1(4),df2(4)
      dimension uio(5),ui1(5),ui2(5),cio(5),ci1(5),ci2(5),aio(5),ai1(5)
      dimension ai2(5),xxi(5),aai(5),cci(5),bbi(5),wk1(5),wk2(5),wk3(5)
      dimension uac(95),wk4(5),wk5(5),wk6(5)
      dimension cpp(5),abc(85),ado(5),ad1(5),ad2(5),bdo(5),bd1(5),fgs(8)
      dimension bd2(5),cdo(5),cd1(5),cd2(5),gdo(5),gd1(5),gd2(5)
      dimension ggsi(5),zzi(5),vvi(5),hhi(5),ggi(5)
      dimension asp(3),bsp(3),csp(3),gsp(3),aspa(3),bspa(3),cspa(3),gspa &
      (3),abcg(24),wk7(5)
      equivalence(psi,pc)
      equivalence (uio(1),uac(1)),(ui1(1),uac(6)),(ui2(1),uac(11)), &
      (cio(1),uac(16)),(ci1(1),uac(21)),(ci2(1),uac(26)),(aio(1),uac(31)),&
      (ai1(1),uac(36)),(ai2(1),uac(41)),(xxi(1),uac(46)),(aai(1),uac(51)),&
      (cci(1),uac(56)),(bbi(1),uac(61)),(wk1(1),uac(66)),(wk2(1),uac(71)),&
      (wk3(1),uac(76)),(wk4(1),uac(81)),(wk5(1),uac(86)),(wk6(1),uac(91))
      equivalence(abc(1),ado(1)),(abc(6),ad1(1)),(abc(11),ad2(1)), &
      (abc(16),bdo(1)),(abc(21),bd1(1)),(abc(26),bd2(1)),(abc(31),cdo(5)),&
      (abc(36),cd1(1)),(abc(41),cd2(1)),(abc(46),gdo(1)),(abc(51),gd1(1)),&
      (abc(56),gd2(1)),(abc(61),ggsi(1)),(abc(66),zzi(1)), &
      (abc(71),vvi(1)),(abc(76),hhi(1)),(abc(81),ggi(1)),(abcg(1),asp(1)),&
      (abcg(4),bsp(1)),(abcg(7),csp(1)),(abcg(10),gsp(1)),(abcg(13),aspa(1)),&
      (abcg(16),bspa(1)),(abcg(19),cspa(1)),(abcg(22),gspa(1))


! for the coulomb corrections
      double precision lami,inv_lami,lamidd,ktinv, &
                       plasg,plasgdd,plasgdt, &
                       aa1,bb1,cc1,dd1,ee1,aa2,bb2,cc2, &
                       ecoul,decouldd,decouldt, &
                       pcoul,dpcouldd,dpcouldt, &
                       scoul,dscouldd,dscouldt, &
                       third,qe,esqu,forth,pie
      parameter        (third = 1.0d0/3.0d0, &
                        qe    = 4.8032068d-10, &
                        esqu  = qe * qe, &
                        forth = 4.0d0/3.0d0, &
                        pie   = 3.1415926535897932384d0)

! see yakovlev & shalybkov 1989
      data aa1,bb1,cc1,dd1,ee1,aa2,bb2,cc2 &
           /-0.898004d0, 0.96786d0, 0.220703d0, -0.86097d0, &
             2.5269d0  , 0.29561d0, 1.9885d0,    0.288675d0/



      data eit/1.d-6/,dst/1.d-3/,grcf/1.3d0/,grif/0.7d0/,tpar/0.3d0/, &
      tt1/0.8d0/,tt2/0.96d0/,tpar1/0.32d0/,ro1/0.12d0/,ro2/3.4884680d-2/ &
      ep2/0.7d0/,ps1/7.2427389d-5/,ro1l/-2.1202635d0/,ro2l/-3.3557075d0/ &
      xep2/0.3d0/,rotp/3.1220424d-10/,pss1/0.1662467d0/, &
       pss2/0.1881164d0/
      data cu/5.93013d0,0.831434d0,1.5d0,1.d0,2.8125d0,5.5957031d0, &
      3.75d0,5.15625d0,1.25d0,0.9375d0,1.8652344d0,5.25d0,1.265625d1, &
      0.5625d0,2.d0,0.5d0,3.d0,2.25d0,0.75d0,2.984375d0,1.3125d1, &
      1.6875d1,3.5d0,2.5d0,1.6787109d1,1.40625d0,2.5d0,0.25d0, &
      5.6568542d0,5.75d0,1.453125d1,4.63125d1,5.8125d1,1.35d1, &
      1.013429d-3,1.33333333333d0,4.d0,2.8613184d-2,0.94051952d-1, &
      0.71428571d0,0.39650038026d-1,0.21875d0,0.66666666667d0,0.6d0, &
      1.05d0,0.18007375d0,0.292178d0,0.3601475d0,0.33333333333d0, &
      5.d0,8.d0,2.66666666667d0,24.d0,15.d0,6.d0,9.d0,1.2d1, &
      1.02677135171d+1,4.9305117064d0,4.80195683116d-1, &
      8.09755744169d-2,7.99196885645d-1/
      data a1/6.78196814d-1,6.78183752d-1,6.94779549d-1,7.60042563d-1/
      data a2/5.37738527d-1,5.37346659d-1,4.87559266d-1,4.09243650d-1/
      data a3/1.73981666d-1,1.70062988d-1,2.19850381d-1,0.251176627d0/
      data a4/2.38231503d-2,1.07608902d-2,-5.83490747d-3,-.0100117403d0/
      data df1/3.47963332d-1,3.40125976d-1,4.39700762d-1,0.502353254d0/
      data df2/7.14694509d-2,3.22826706d-2,-1.75047241d-2, &
               -3.00352209d-2/
      data  b1/1.15292705d0,1.15292656d0,1.14670313d0,1.08551906d0/
      data  b2/1.01729522d0,1.01727563d0,1.04216932d0,1.14006384d0/
      data  b3/4.03303895d-1,4.03009994d-1,3.65669449d-1,3.06932737d-1/
      data b4/8.69908330d-2,8.50314940d-2,1.09925190d-1,1.25588313d-1/
      data b5/8.93368136d-3,4.03533382d-3,-2.18809030d-3,-3.75440261d-3/
      data c1/3.08286343d0,3.08286341d0,3.08597512d0,3.16245521d0/
      data c2/2.88231761d0,2.88231639d0,2.86675783d0,2.71379764d0/
      data c3/1.27161903d0,1.27159453d0,1.30271165d0,1.42507981d0/
      data c4/3.36086579d-1,3.35841662d-1,3.04724541d-1,2.55777281d-1/
      data c5/5.43692706d-2,5.31446837d-2,6.87032441d-2,7.84926959d-2/
      data c6/4.46684068d-3,2.01766691d-3,-1.09404515d-3,-1.87720131d-3/
      data d1/1.11841602d1,1.11841602d1,1.11823451d1,1.10708116d1/
      data d2/1.07900220d1,1.07900219d1,1.08009129d1,1.10685932d1/
      data d3/5.04405582d0,5.04405368d0,5.01682620d0,4.74914588d0/
      data d4/1.48355553d0,1.48352696d0,1.51983026d0,1.66259311d0/
      data d5/2.94075757d-1,2.93861454d-1,2.66633974d-1,2.23805121d-1/
      data d6/3.80584894d-2,3.72012786d-2,4.80922708d-2,5.49448872d-2/
      data d/2.60565706d-3,1.17697237d-3,-6.38193005d-4,-1.09503410d-3/
      data a/-3.53553400d-1,1.92450000d-1,-1.25000000d-1,8.94427000d-2 &
            ,-1.76777000d-1,6.41500000d-2,-3.12500000d-2,1.78885000d-2 &
            ,-8.83883000d-2,2.13833000d-2,-7.81250000d-3,3.57771000d-3 &
            ,1.16317000d1,-4.41942000d-2,7.12778000d-3,-1.95313000d-3 &
            ,7.15541750d-4/
      data b/6.666667d-1,8.22467d-1,7.102746d-1,6.467679d0 &
            ,4.00000000d-1,2.46740000d0,-7.10275000d-1,-2.77186000d0 &
            ,2.85714000d-1,4.11234000d0,3.55137000d0,2.77186000d0 &
            ,2.22222200d-1,5.75727000d0,2.48596000d1,-6.4676800d0 &
            ,7.79429075d-3,-4.94746507d-2,1.94857269d-2,3.56600189d-2 &
            ,-1.73161277d-1,3.41000220d-2/
      data c/-7.07106800d-1,5.77350000d-1,-5.00000000d-1,4.47213500d-1 &
            ,1.00000015d0,-4.11233500d-1,-1.77568650d0,-2.91045555d1/
      data ck/8.86227d-1,1.32934d0,3.32335d0/
      data uio/0.43139881d0,1.7597537d0,4.1044654d0,7.7467038d0, &
               1.3457678d1/
      data ui1/0.81763176d0,2.4723339d0,5.1160061d0,9.0441465d0, &
               1.5049882d1/
      data ui2/1.2558461d0,3.2070406d0,6.1239082d0,1.0316126d1, &
               1.6597079d1/
      data cio/0.37045057d0,0.41258437d0,9.777982d-2,5.3734153d-3, &
               3.8746281d-5/
      data ci1/0.39603109d0,0.69468795d0,0.2232276d0,1.5262934d-2, &
               1.3081939d-4/
      data ci2/0.76934619d0,1.7891437d0,0.70754974d0,5.6755672d-2, &
               5.557148d-4/
      data aio/0.64959979d0,0.17208724d0,0.016498837d0,4.321647d-4, &
               1.4302261d-6/
      data ai1/0.44147594d0,8.4387677d-2,5.9999383d-3,1.180802d-4, &
               2.9101763d-7/
      data ai2/0.28483475d0,4.0476222d-2,2.1898807d-3,3.3095078d-5, &
               6.194128d-8/
      data xxi/7.265351d-2,0.2694608d0,0.5331220d0,0.7868801d0, &
               0.9569313d0/
      data aai/3.818735d-2,0.1256732d0,0.1986308d0,0.1976334d0, &
               0.1065420d0/
      data cci/0.26356032d0,1.4134031d0,3.5964258d0,7.0858100d0, &
               1.2640801d1/
      data bbi/0.29505869d0,0.32064856d0,7.391557d-2,3.6087389d-3, &
               2.3369894d-5/
      data pc1/0.5d0/,pc2/0.7d0/
      data cpp/5.d0,1.d1,1.5d1,2.5d1,2.5d0/
      data fgs/0.571428571d0,0.33333333d0,0.2272727d0,0.168269d0, &
       0.142857143d0,5.5555555d-2,2.840909d-2,1.68269d-2/
      data (abc(k),k=1,76)/7.72885519d1, &
        1.42792768d2, 4.30552910d1, 2.43440537d0, &
        1.75674547d-2, 9.99400362d1, 2.73430265d2, 1.00130386d2, &
        6.91871969d0, 5.93132645d-2, 2.30460043d2, 7.56122303d2, &
        3.19543255d2, 2.57313963d1, 2.51960145d-1,-2.35425685d1, &
       -4.38697759d1,-1.32985534d1,-7.52438243d-1,-5.42994019d-3, &
       -3.05287674d1,-8.42357074d1,-3.09412790d1,-2.13850223d0, &
       -1.83331906d-2,-7.06062732d1,-2.33317365d2,-9.87584116d1, &
       -7.95332909d0,-7.78785901d-2, 1.42401164d-1, 4.12778210d-1, &
        1.52786427d-1, 8.84665279d-3, 6.38815164d-5, 2.18702630d-1, &
        8.82651141d-1, 3.60865258d-1, 2.51545288d-2, 2.15684504d-4, &
        5.87304073d-1, 2.59226969d0, 1.15817403d0, 9.35640728d-2, &
        9.16218624d-4, 2.94914091d-1, 5.29893251d-1, 1.56942521d-1, &
        8.85295620d-3, 6.38816670d-5, 3.77889882d-1, 1.00545595d0, &
        3.64435019d-1, 2.51594259d-2, 2.15684607d-4, 8.63109766d-1, &
        2.76526224d0, 1.16235562d0, 9.35691781d-2, 9.16218718d-4, &
        7.22774549d-1, 6.91700407d-1, 6.36940508d-1, 5.69038300d-1, &
        5.14846835d-1, 9.63560320d-1, 2.11340310d0, 4.29642580d0, &
        7.78581000d0, 1.33408010d1, 5.08574570d-2, 1.88622560d-1, &
        3.73185400d-1, 5.50816070d-1, 6.69851910d-1, 2.89632880d-1/
        data (abc(k),k=77,85)/ &
        4.66144392d-1, 1.53210873d-1, 1.00694874d-2, 8.53586810d-5, &
        1.46896384d-2, 4.60107809d-2, 6.75861819d-2, 6.21820743d-2, &
        3.16690579d-2/
          data nitm/60/
          data pi2/9.8696044011d0/
          data t5/1.3d1/,t4/1.5d1/,ro3/3.d1/,ro4/3.9d1/
          data nfil/1/

! ***
!         if(nfil.ne.1) go to 4000
!         open(unit=101,file='epeos.mes')
!         write(6,5001)
!         nfil=0
!4000  if((nz.ge.0).and.(nz.le.5)) go to 4004
!       write(6,5002) nz
!       print 4005,nz
!       print 4006
!       stop


4005  format('  illegal value of nz:  nz=',i5)
4006  format(' allowed values are nz=0,1,2,3,4,5  *stop in epeos*')
5001    format(10x,'epeos  ***  version 1.1 august 2, 1992  ***  epeos')
5002    format(20x,'epeos  ***  error in   nz   ***  epeos'/ &
       1x,'illegal value of nz:  nz=',i5, &
       1x,'! allowed values are nz=0,1,2,3,4,5 *stop*')
5012    format(20x,'epeos  ***  error in   t    ***  epeos'/ &
       1x,1p,'temperature t must be positive: t=',d13.6, &
       1x,'jkk=',i4,' *stop*')
5013    format('  temperature t must be positive: t=', &
       1p,d13.6,' jkk=',i4,' *stop in epeos*')
5022    format(20x,'epeos  ***     warning!     ***  epeos'/ &
       1x,1p,'negative or zero density: den=',d13.6,' jkk=',i4/ &
       1x,'calculations are going on with den=0.')
5023    format('  negative or zero density den=', &
        d13.6,' jkk=',i4,' *warning in epeos* ')
5032    format(20x,'epeos  *** error in as or zs ***  epeos'/ &
       1x,1p,'as and zs must be positive: as=',d10.3,' zs=',d10.3, &
       1x,' jkk=',i4,' *stop*')
5033    format(' as and zs must be positive: as=',1p,d10.3, &
       1x,'zs=',d10.3,' jkk=',i4,'*stop in epeos*')
!
4004  continue


! use the given abar & zbar, get derivatives, entropy, include pairs

       jurs  = 0
       lst   = 2
       kentr = 1
       kpar  = 1
       eosfail = .false.


! start the pipeline loop
!      write(6,*) 'starting pipeline'

      do j = jlo_eos,jhi_eos
       t   = temp_row(j) * 1.0d-9
       den = den_row(j) * 1.0d-7
       as  = abar_row(j)
       zs  = zbar_row(j)
       scn = 2.5d0 * log(abar_row(j))

! flag set to find the working region
       nz  = 0

! check the input
      if(t.gt.0.d0) go to 4010
       write(6,5012) t,jkk
       stop
4010  continue
      if(pl.gt.0.d0) go to 102
      write(6,5022) pl,jkk
      eosfail = .true.
      stop

!      pt=3.025884d-2*t**3/cu(17)
!      p=0.25d0*t*pt
!      ppl=0.d0
!      gam=cu(36)
!      da=0.25d0
!      go to 9


  102 if(jurs.eq.1) stop 'tried a call to chemic'
! call chemic
! *** subroutine chemic calculates as,zs,scn
      if((as.gt.0.d0).and.(zs.gt.0.d0)) go to 4030
      write(6,5032) as,zs,jkk
      eosfail = .true.
      stop
!
! start normal calculations
4030  emue=as/zs
      rg=cu(2)/emue
      ki=1
      hpr=0.d0
   90 alf=cu(1)/t
      al1=cu(4)/alf
      plm=pl/emue
      sqe=0.d0
      ei=t*rg
      pi=ei*pl
      eg=cu(3)*ei
  590  continue
!
! *** search for required working region
      if(nz.ne.0) go to(1,2,3,4,5),nz
      if(ki.ne.1) go to  123
      if(plm.le.ro3) go to 310
      if(plm.ge.ro4) go to 4
      if(t.le.t5) go to 550
      if(t.ge.t4) go to 4
! *** searching around the triangle
          x=(ro4-ro3)/(t4-t5)
          y=ro3+(t4-t)*x
          if(y.gt.plm) go to 800
          go to 4
  310  if(t.le.t5) go to 123
          if(t.ge.t4) go to 4
! *** interpolation over t in region 45
          t1=t5
          t2=t4
          nz2=4
          nz1=5
          go to 128
  550  continue
! *** interpolation over density for density < ro4
      nzp1=0
      if(t.lt.tt2) nzp1=3
  583  kin=ki
       nz=nzp1
       ki=6
       go to 590

  577  pn1=pe
       en1=ee
       sn1=se
       psn1=psi
       hprn=hpr
       pnp1=ppl
       enp1=epl
       snp1=spl
       pnt=pt
       ent=et
       snt=st
       nzr1=nzr
       nz=4
       ki=7
       go to 590
  578  z1=ro4-ro3
       x2=(plm-ro3)/z1
       x=x2**2*(3.d0-2.d0*x2)
       x1=1.d0-x
       z1=z1*emue
       x3=6.d0*x2*(1.d0-x2)/z1
       p=pe
       e=ee
       s=se
       pe=pe*x+pn1*x1
       ee=ee*x+en1*x1
       if(kentr.eq.0) go to 591
       se=se*x+sn1*x1
  591 if(lst.eq.0) go to 592
      pt=pt*x+pnt*x1
      et=et*x+ent*x1
      ppl=ppl*x+pnp1*x1+(p-pn1)*x3
      epl=epl*x+enp1*x1+(e-en1)*x3
      if(kentr.eq.0) go to 592
      st=st*x+snt*x1
      spl=spl*x+snp1*x1+(s-sn1)*x3
  592  psi=psi*x+psn1*x1
      hpr=hpr*x+hprn*x1
      nzr=10*nzr1+nzr
      ki=kin
      go to 134
  800 continue
! ***** the triangle
          nz2=4
         nz1=0
          t1=t5
          t2=t4-(plm-ro3)/x
          go to 128
  123 if(ki.ne.4) go to 136
      if(nz2.eq.5) go to 111
      nzp1=5
          go to 583
  136 if(plm.lt.ps1) go to 110
      if(t.lt.tt1) go to 111
      kzz=2
         go to 121
  115 y=x*grcf
      if(plm.gt.y) go to 3
      kzz=3
         z1=t
      if(t.lt.tt2) go to 112
      if(plm.lt.x) go to 5
  113 kkz=0
  116 go to 121
  114 dl=(x-plm)/xz
      x1=1.d0
      if(dl.lt.0.d0) x1=-1.d0
      x2=dl*x1
      if(x2.gt.0.3d0) dl=0.3d0*x1
      t=t*(1.d0-dl)
      if(x2.gt.eit)go to 116
      if(kkz.eq.1) go to 118
      t2=t
      t=z1
  138 z2=plm
         plm=plm/grcf
         kkz=1
      go to 116
  118 t1=t
         nz1=3
         plm=z2
         t=z1
         nz2=5
      go to 128
  112 if(plm.gt.pss1) go to 117
      t1=tt1
         t2=tt2
        nz1=0
        nz2=5
      go to 128
  117 if(plm.gt.pss2) go to 113
      t2=tt2
      go to 138
  110 if(t.lt.tpar) go to 111
      if(t.gt.tt2) go to 5
      sqe=exp(-alf)
      if(t.gt.tt1) go to 119
      if(plm.gt.ro1*sqe) go to 1
      if(t.gt.tpar1) go to 119
      if(plm.gt.rotp) go to 120
      t1=tpar
      t2=tpar1
  122 nz1=1
      nz2=5
      sqe=0.d0
      go to 128
  119 if(plm.lt.ro2*sqe) go to 5
  120 t1=log(plm)
      t2=cu(1)/(ro2l-t1)
      t1=cu(1)/(ro1l-t1)
      go to 122
  111 sq=t*sqrt(t)
      y=2.095d-2*sq
      x=y*grif
      if(plm-x)44,44,51
!
! perfect gas with corrections for degeneracy and pairs (nz=1)
    1 sq=t*sqrt(t)
   44 nzr=1
      qa=6.9712909d0*plm/sq
      x=qa*al1
      x2=qa-x*(cu(5)-cu(6)*al1)
      pe=cu(4)+x2
      epl=qa-x*(cu(10)+cu(11)*al1)
      x8=cu(9)*al1
      x3=x8*(cu(4)+(cu(14)*al1-cu(4))*al1)
      ee=cu(4)+epl+x3
      if(lst.eq.0) go to 6
      pt=cu(4)-cu(16)*x2-x*(cu(5)-cu(15)*cu(6)*al1)
      ppl=cu(4)+cu(15)*x2
      et=cu(4)-cu(16)*epl+x8*(cu(15)+al1*(cu(18)*al1-cu(17))-qa* &
         (cu(19)+cu(20)*al1))
    6 if(t.lt.tpar) go to 8
      if(sqe.eq.0.d0) sqe=exp(-alf)
      x4=0.268192d0*(al1*sqe/plm)**2
      x5=x4*al1
      hpr=x5*(cu(4)+al1*(cu(7)+cu(8)*al1))
      pe=pe+hpr
      x6=(x4+x5*(cu(12)+cu(13)*al1))/cu(3)
      ee=ee+x6
      if(lst.eq.0) go to 8
      ppl=ppl-hpr
      pt=pt+hpr*cu(15)*(cu(15)+alf)+x5*(cu(7)+cu(21)*al1)*al1
      x8=x6*cu(15)
      et=et+x8*(cu(3)+alf)+x5*(cu(23)+cu(22)*al1)
      epl=epl-x8
    8 if((kentr+ki).eq.1) go to 56
      x7=cu(16)*(qa+x*(cu(5)-cu(25)*al1))
      psi=log(cu(29)*qa)
      se=cu(24)-psi+x7+al1*(cu(7)+al1*(al1*cu(26)-cu(5)))
      psi=psi+2.d0*x2+0.5d0*hpr-1.875d0*al1* &
       (1.d0+al1*(al1*0.1875d0-0.5d0))
      if(kentr.eq.0) go to 56
      if(lst.eq.0) go to 53
      st=cu(3)*(cu(4)+al1*(cu(27)+al1*(cu(5)*al1-cu(7))))
      st=st-cu(28)*qa*(cu(17)+al1*(cu(5)+cu(25)*al1))
      spl=x7-cu(4)
   53 if(t.lt.tpar) go to 101
      x8=x4+x5*(cu(30)+cu(31)*al1)
      se=se+x8
      if(lst.eq.0) go to 10
      spl=spl-cu(15)*x8
      st=st+x4*(cu(15)*alf+cu(34)+al1*(cu(32)+cu(33)*al1))
  101 if(lst.eq.0) go to 10
      st=st*rg/t
      spl=spl*rg
   10 se=se*rg
   56 pe=pi*pe
      ee=eg*ee
      if(lst.eq.0) go to 50
      pt=rg*pl*pt
      et=eg*et/t
      ppl=ei*ppl
      epl=eg*epl
   50 go to 135
!
! *** addition the ion and black-body radiation components to eos
   57 continue
! *********************************************************************
      x=pi/zs
      x1=eg/zs
      v=7.56471d-3*t**4
!      v = 7.56590624067000303d-3 * t**4
      x3=v/cu(17)
      p=x+pe+x3

      pion = x
      dpiondd = x/pl
      dpiondt = x/t


      prad = x3


      beta=x3/p
      x3=v/pl
      e=x1+x3+ee

      eion = x1
      erad = x3

      if(kentr.eq.0) go to 7
      x6=cu(2)/as
      x4=pl*(cu(35)/(t*sqrt(t)))
      x4=log(x4)
      x4=x6*(cu(24)-x4+scn)
      x5=cu(36)*x3/t
      s=x5+x4+se

      sion = x4/cu(2)
      srad = x5/cu(2)

      sek=se/cu(2)
      sk=s/cu(2)
      if(lst.eq.0) go to 9

! store the electron-positron entropy derivatives before adding rad+ion
      dsepdd = spl
      dsepdt = st
      st=st+(cu(17)*x5+x1/t)/t
      spl=spl-x6-x5

! dsraddt and dsraddd
!      write(6,*) cu(17)*x5/t*1.0d-1,-x5*1.0d1/den

! dsiondt and dsiondd
!      write(6,*) x1/t/t * 1.0d-1,-x6*1.0d1/den


      go to 45
    7 if(lst.eq.0) go to 9
   45 continue

! store the electron-positron pressure derivatives before adding rad+ion
      dpepdd = ppl
      dpepdt = pt
      pt=pt+(x+cu(36)*v)/t
      ppl=ppl+x/pl

! store the electron-positron energy derivatives before adding rad+ion
      deepdd = epl
      deepdt = et
      et=et+(x1+cu(37)*x3)/t
      epl=epl-x3

! *********************************************************************
      x4=pt/(pl*et)
      x5=pl/p
      gam=x5*(ppl+t*x4*pt/pl)
      da=x4/gam
      cp=gam*et*(p/ppl)
      dpe=x5*(epl+t*pt/pl)-1.d0
      if(kentr.ne.0) then
      dse=t*st/et-1.d0
      dsp=-spl*(pl/pt)-1.d0
      endif


! normal exit from epeos
! set the return arguments in cgs units
   9  continue

!  *** dpe = (den/p)(epl+t*pt/den)-1 -- thermodynamic identity: dpe=0
!  *** dse = t*st/et-1 -- thermodynamic identity: dse=0
!  *** dsp = -spl*(den/pt)-1 -- thermodynamic identity: dsp=0
!       dse = t*st/et - 1.0d0
!       dpe = (epl*den + t*pt)/p - 1.0d0
!       dsp = -spl*(den/pt) -1.0d0
!       write(6,696) dse,dpe,dsp

      p       = p       * 1.0d24
      prad    = prad    * 1.0d24
      pion    = pion    * 1.0d24
      pe      = pe      * 1.0d24

      pt      = pt      * 1.0d15
      dpiondt = dpiondt * 1.0d15
      dpepdt  = dpepdt  * 1.0d15

      ppl     = ppl     * 1.0d17
      dpiondd = dpiondd * 1.0d17
      dpepdd  = dpepdd  * 1.0d17

      e       = e       * 1.0d17
      erad    = erad    * 1.0d17
      eion    = eion    * 1.0d17
      ee      = ee      * 1.0d17

      et      = et      * 1.0d8
      deepdt  = deepdt  * 1.0d8

      epl     = epl     * 1.0d10/den
      deepdd  = deepdd  * 1.0d10/den

      s       = s       * 1.0d8
      sion    = sion    * 1.0d8
      srad    = srad    * 1.0d8
      se      = se      * 1.0d8

      st      = st      * 1.0d-1
      dsepdt  = dsepdt  * 1.0d-1

      spl     = spl     * 1.0d1/den
      dsepdd  = dsepdd  * 1.0d1/den


! coulomb section:
! initialize

       pcoul    = 0.0d0
       dpcouldd = 0.0d0
       dpcouldt = 0.0d0
       ecoul    = 0.0d0
       decouldd = 0.0d0
       decouldt = 0.0d0
       scoul    = 0.0d0
       dscouldd = 0.0d0
       dscouldt = 0.0d0


! uniform background corrections only
! see yakovlev & shalybkov 1989
! lami is the average ion seperation
! plasg is the plasma coupling parameter
       z        = forth * pie
       xni      = avo/as * den_row(j)
       dxnidd   = avo/as
       sss      = z * xni
       dsdd     = z * dxnidd

       lami     = 1.0d0/sss**third
       inv_lami = 1.0d0/lami
       z        = -third * lami/sss
       lamidd   = z * dsdd
       ktinv    = 1.0d0/(kerg*temp_row(j))

       plasg    = zs*zs*esqu*ktinv*inv_lami
       z        = -plasg * inv_lami
       plasgdd  = z * lamidd
       plasgdt  = -plasg*ktinv * kerg


! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg .ge. 1.0) then
         x        = plasg**(0.25d0)
         y        = avo/as * kerg
         ecoul    = y * temp_row(j) * (aa1*plasg + bb1*x + cc1/x + dd1)
         pcoul    = third * den_row(j) * ecoul
         scoul    = -y * (3.0d0*bb1*x - 5.0d0*cc1/x &
                    + dd1 * (log(plasg) - 1.0d0) - ee1)

         y        = avo/as*kerg*temp_row(j) &
                    * (aa1 + 0.25d0/plasg*(bb1*x - cc1/x))
         decouldd = y * plasgdd
         decouldt = y * plasgdt + ecoul/temp_row(j)

         y        = third * den_row(j)
         dpcouldd = third * ecoul + y*decouldd
         dpcouldt = y * decouldt

         y        = -avo*kerg/(as*plasg) &
                      *(0.75d0*bb1*x +1.25d0*cc1/x +dd1)
         dscouldd = y * plasgdd
         dscouldt = y * plasgdt


! yakovlev & shalybkov 1989 equations 102, 103, 104
        else if (plasg .lt. 1.0) then
         x        = plasg*sqrt(plasg)
         y        = plasg**bb2
         z        = cc2 * x - third * aa2 * y
         pcoul    = -pion * z
         ecoul    = 3.0d0 * pcoul/den_row(j)
         scoul    = -avo/as*kerg*(cc2*x -aa2*(bb2-1.0d0)/bb2*y)

         sss      = 1.5d0*cc2*x/plasg - third*aa2*bb2*y/plasg
         dpcouldd = -dpiondd*z - pion * sss * plasgdd
         dpcouldt = -dpiondt*z - pion * sss * plasgdt

         sss      = 3.0d0/den_row(j)
         decouldd = sss * dpcouldd - ecoul/den_row(j)
         decouldt = sss * dpcouldt

         sss      = -avo*kerg/(as*plasg) &
                      * (1.5d0*cc2*x -aa2*(bb2-1.0d0)*y)
         dscouldd = sss * plasgdd
         dscouldt = sss * plasgdt
        end if


! bomb proof
       x   = p + pcoul
       if (x .le. 0.0) then
        write(6,*)
        write(6,*)
        write(6,*) 'coulomb corrections are causing a negative pressure'
        write(6,*) 'setting all coulomb corrections to zero'

        pcoul    = 0.0d0
        dpcouldd = 0.0d0
        dpcouldt = 0.0d0
        ecoul    = 0.0d0
        decouldd = 0.0d0
        decouldt = 0.0d0
        scoul    = 0.0d0
        dscouldd = 0.0d0
        dscouldt = 0.0d0

        write(6,*)
        write(6,*)
       end if



! and add in the coulomb corrections
      p   = p + pcoul
      pt  = pt + dpcouldt
      ppl = ppl + dpcouldd

      e   = e + ecoul
      et  = et + decouldt
      epl = epl + decouldd

      sk  = sk + scoul/(kerg*avo)
      s   = s  + scoul
      st  = st + dscouldt
      spl = spl + dscouldd



! the temperature and density exponents (c&g 9.81 9.82)
! the specific heat at constant volume (c&g 9.92)
! the third adiabatic exponent (c&g 9.93)
! the first adiabatic exponent (c&g 9.97)
! the second adiabatic exponent (c&g 9.105)
! the specific heat at constant pressure (c&g 9.98)
! and relativistic formula for the sound speed (c&g 14.29)
      zz    = p/den_row(j)
      chit  = temp_row(j)/p * pt
      chid  = ppl/zz
      cv    = et
      x     = zz * chit/(temp_row(j) * cv)
      gam3  = x + 1.0d0
      gam1  = chit*x + chid
      nabad = x/gam1
      gam2  = 1.0d0/(1.0d0 - nabad)
      cp    = cv * gam1/chid
      z     = 1.0d0 + (e + con2)/zz
      sound = clight * sqrt(gam1/z)

! these are the 3 thermodynamic consistency checks
! each is zero if the consistency is perfect
        dse = temp_row(j)*st/et - 1.0d0
        dpe = (epl*den_row(j)*den_row(j) + temp_row(j)*pt)/p - 1.0d0
        dsp = -spl*(den_row(j)**2/pt) - 1.0d0

! store this row and go to the end of the do loop
        ptot_row(j)   = p
        dpt_row(j)    = pt
        dpd_row(j)    = ppl
        dpa_row(j)    = 0.0d0
        dpz_row(j)    = 0.0d0

        etot_row(j)   = e
        det_row(j)    = et
        ded_row(j)    = epl
        dea_row(j)    = 0.0d0
        dez_row(j)    = 0.0d0

        stot_row(j)   = s
        dst_row(j)    = st
        dsd_row(j)    = spl
        dsa_row(j)    = 0.0d0
        dsz_row(j)    = 0.0d0

        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad

        pion_row(j)   = pion
        eion_row(j)   = eion
        sion_row(j)   = sion
        xni_row(j)    = avo/as * den_row(j)

        pele_row(j)   = pe
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = 0.0d0
        dpepz_row(j)  = 0.0d0

        eele_row(j)   = ee
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = 0.0d0
        deepz_row(j)  = 0.0d0

        sele_row(j)   = se
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt
        dsepd_row(j)  = dsepdd
        dsepa_row(j)  = 0.0d0
        dsepz_row(j)  = 0.0d0

        xnem_row(j)   = zs * xni_row(j)
        xne_row(j)    = xnem_row(j) + 0.5d0 * hpr * xnem_row(j)
        dxnet_row(j)  = 0.0d0
        dxned_row(j)  = avo/emue * (1.0d0 + hpr)
        dxnea_row(j)  = 0.0d0
        dxnez_row(j)  = 0.0d0
        xnp_row(j)    = 0.5d0 * hpr * xnem_row(j)

        etaele_row(j) = psi
        detat_row(j)  = 0.0d0
        detad_row(j)  = 0.0d0
        detaa_row(j)  = 0.0d0
        detaz_row(j)  = 0.0d0
        etapos_row(j) = 0.0d0
        beta          = kerg * temp_row(j)/mecc
        if (beta .gt. 0.02) etapos_row(j) = -psi - 2.0d0/beta

        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul
        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound

      goto 6969

!
! *** further search for working region
   51 if(plm.lt.y) go to 52
      kzz=1
  121 z=t/ep2
      if(z.gt.0.57142857d0) go to 64
      x=xep2*sq
      if(kzz.eq.3) xz=cu(3)*x
      go to 125
   64 if(z.gt.1.07142857d0) go to 54
      xz=7.306392d-2*z
      x=xz+3.414385d-2
      go to 125
   54 if(z.gt.1.d1) go to 58
      x=((4.77856d-2*z-1.41839d-2)*z+7.2001d-2)*z-7.20897d-3
      if(kzz.eq.3) xz=((1.433568d-1*z-2.83678d-2)*z+7.2001d-2)*z
      go to 125
   58 x=4.708d-2*z**3
      if(kzz.eq.3) xz=cu(17)*x
  125 go to (55,115,114),kzz
   55 if(plm-x) 59,59,61
!
! *** expansion over half-integer fermi--dirac functions (nz=2)
    2 sq=t*sqrt(t)
   59 nzr=2
      kenf=1
      zp=plm/(cu(38)*sq)
      x=psi
      nit=0
      kf=1
      x1=cu(9)*al1
      al2=al1**2
      al3=0.21875d0*al2
      go to 26
   11 z=zp-f12-x1*f32-al3*f52
      v=f12s+x1*f32s+al3*f52s
      dl=z/v
      z1=x
      if(z1.lt.0.d0)z1=-x
      v=1.d0
      if(dl.lt.0.d0) v=-1.d0
      z=v*dl
      if(z1.lt.0.1d0) go to 41
      z=z/z1
      if(z.gt.0.3d0) dl=0.3d0*v*z1
   41 x=x+dl
      nit=nit+1
      if(z.lt.eit) go to 73
      if(nit.lt.nitm) go to 26
   66  write(6,5072) nzr
       write(6,5073) dl,x,eit,nitm,jkk
!        print 72,nzr
!        print 4072,dl,x,eit,nitm,jkk
       eosfail = .true.
       stop
   72  format('  iterations in region nzr=',i1, &
       1x,'do not converge!')
4072   format('  dx=',1p,d10.3,' x=',d10.3,' eit=',d10.3, &
       1x,'nitm=',i3,' jkk=',i4,' *stop in epeos*')
5072   format(10x,'epeos  ***  runaway of iterations in region nzr=', &
       i1,' *** stop*')
5073   format(10x,1p,'dx=',d10.3,' x=',d10.3,' eit=',d10.3,' nitm=',i3, &
       1x,'jkk=',i4)
!
   73 kf=2
      go to 26
  300 go to(31,9),kenf
   31 psi=x
      v=sq*al1
      x=cu(19)*al1
      z=cu(39)*v
      z1=cu(3)*z/pl
      z3=x1*f32
      z4=x*f52
      al4=9.375d-2*al2
      y=al4*f72
      pe=z*(f32+z4+y)
      z5=x1*f52
      y1=al3*f72
      ee=z1*(f32+z5+y1)
      if(kentr.eq.0) go to 34
      z2=cu(41)*sq/pl
      dl=cu(40)*psi
      z10=cu(44)*psi
      z11=cu(45)*al1
      z6=z11*(f52-dl*f32)
      al5=0.16875d0*al2
      y3=0.77777777778d0*psi
      y2=al5*(f72-y3*f52)
      se=z2*(f32-z10*f12+z6+y2)
   34 if(lst.eq.0) go to 35
      pal=f12s+x1*f32s+al3*f52s
      pap=zp/pal
      pal=(cu(3)*zp+z3+cu(15)*al3*f52)/pal
      z8=f32s+x1*f52s+al3*f72s
      z7=z5+cu(15)*y1-pal*z8
      et=(cu(27)*ee+z1*z7)/t
      epl=z1*z8*pap-ee
      z9=f32s+x*f52s+al4*f72s
      pt=(cu(27)*pe+z*(z4+cu(15)*y-pal*z9))/t
      ppl=z*z9*pap/pl
      if(kentr.eq.0) go to 35
      z9=f32s-z10*f12s-cu(44)*f12+z11*(f52s-dl*f32s-cu(40)*f32)
      z9=z9+al5*(f72s-0.77777777778d0*f52-y3*f52s)
      st=(cu(3)*se+z2*(z6+cu(15)*y2-pal*z9))/t
      spl=z2*pap*z9-se
   35 go to 135
!
! *** procedure of calculation of half-integer fermi--dirac functions
!       entry fd12f
!      kenf=2
!      kf=2
!       x=psi
   26 if(x.ge.-1.d0) go to 21
      z=exp(x)
      v=ck(1)*z
      f12=v*(1.d0+z*(a(1)+z*(a(2)+z*(a(3)+z*a(4)))))
      f12s=v*(1.d0+z*(c(1)+z*(c(2)+z*(c(3)+z*c(4)))))
      f32=ck(2)*z*(1.d0+z*(a(5)+z*(a(6)+z*(a(7)+z*a(8)))))
      f52=ck(3)*z*(1.d0+z*(a(9)+z*(a(10)+z*(a(11)+z*a(12)))))
      go to(30,12),kf
   12 f72=a(13)*z*(1.d0+z*(a(14)+z*(a(15)+z*(a(16)+z*a(17)))))
      go to 30
   21 if(x.ge.-.1d0) go to 22
      n=1
      go to 14
   22 if(x.ge.1.d0) go to 23
      n=2
      go to 14
   23 if(x.ge.2.5d0) go to 24
      n=3
      go to 14
   24 if(x.ge.4.5d0) go to 25
      n=4
   14 f12=a1(n)+x*(a2(n)+x*(a3(n)+x*a4(n)))
      f12s=a2(n)+x*(df1(n)+x*df2(n))
      f32=b1(n)+x*(b2(n)+x*(b3(n)+x*(b4(n)+x*b5(n))))
      f52=c1(n)+x*(c2(n)+x*(c3(n)+x*(c4(n)+x*(c5(n)+x*c6(n)))))
      go to(30,13),kf
   13 f72=d1(n)+x*(d2(n)+x*(d3(n)+x*(d4(n)+x*(d5(n)+x*(d6(n)+x*d(n))))))
   30 f32s=1.5d0*f12
      f52s=2.5d0*f32
      if(kf.eq.1) go to 11
      f72s=3.5d0*f52
      go to 300
   25 z=sqrt(x)
      z1=z*x
      z2=1.d0/x**2
      f12=z1*(b(1)+(b(2)+(b(3)+b(4)*z2)*z2)*z2)
      f12s=z*(c(5)+(c(6)+(c(7)+c(8)*z2)*z2)*z2)
      z=z1*x
      f32=z*(b(5)+(b(6)+(b(7)+b(8)*z2)*z2)*z2)
      f52=x*z*(b(9)+(b(10)+(b(11)+b(12)*z2)*z2)*z2)
      f32=f32+b(17)
      f52=f52+b(18)+b(19)*x
      go to(30,17),kf
   17 f72=z*(b(13)+(b(14)+(b(15)+b(16)*z2)*z2)*z2)/z2
      f72=f72+b(20)+x*(b(21)+b(22)*x)
      go to 30
!
! *** search for working region is continued
   61 y=x*grcf
      if(plm.lt.y) go to 62
!
! *** chandrasekhar's expansion for degenerate gas (nz=3)
    3 nzr=3
      x1=plm/cu(47)
      if(psi.lt.3.d0) psi=3.d0
      nit=0
      x=psi*al1
      x=sqrt(x*(x+cu(15)))
      al2=al1**2
      z2=1.644934d0*al2
      z4=1.894066d0*al2**2
   74 x2=x**2
      x3=x2*x
      x5=cu(17)*(z4/x3)/x2
      x6=z2/x
      x8=cu(15)*z2*x
      dl=(x3*cu(49)+x8+x6+x5-x1)/(x3+x8-x6-cu(50)*x5)
      x6=1.d0
      if(dl.lt.0.d0) x6=-1.d0
      z1=x6*dl
      if(z1.gt.0.9d0) dl=0.9d0*x6
      x=x*(cu(4)-dl)
      nit=nit+1
      if(z1.lt.eit) go to 71
      if(nit.lt.nitm) go to 74
      go to 66
   71 x2=x**2
      x4=cu(4)+x2
      z=sqrt(x4)
      y1=cu(4)+z
         z1=x2/y1
      psi=alf*z1
      z3=x*z
      x3=x*x2
      z5=cu(15)*x2
      x7=z5+cu(4)
      if(x.gt.0.1d0) go to 174
      x5=x2*x3
      f0=x5*(1.6d0-x2*(fgs(1)-x2*(fgs(2)-x2*(fgs(3)-x2*fgs(4)))))*cu(49)
      g0=x5*(0.8d0-x2*(fgs(5)-x2*(fgs(6)-x2*(fgs(7)-x2*fgs(8)))))
      go to 175
  174 x6=log(x+z)
      f0=z3*(z5-cu(17))*cu(49)+x6
      g0=x7*z3-x6-x3*cu(52)
  175 x5=z4/x3
      y5=z5-cu(4)
      f2=z2*cu(51)*z3
      f4=x5*cu(51)*y5*z
      pe=cu(46)*(f0+f2+f4)
      y2=cu(51)/y1
      y4=z*y1
      g2=z2*y2*x*(x7+y4)
      g4=x5*cu(17)*y2*(cu(4)+y5*y4)
      ee=cu(46)*(g0+g2+g4)/pl
      if(kentr.eq.0) go to 75
      y6=cu(48)/(t*pl)
      se=y6*(f2+cu(15)*f4)
   75 if(lst.eq.0) go to 76
      z6=cu(54)*x5/x3
      z7=x7-cu(15)
      pap=x2+z2*z7/x2-z6
      pal=cu(15)*(z2*x7+cu(55)*x5/x)/(x*pap)
      pap=x1/pap
      z9=cu(51)*x/z
      z10=x1*z9
      z11=cu(46)*pap/pl
      ppl=z11*z10
      y3=cu(46)/t
      pt=y3*(cu(15)*(f2+cu(15)*f4)-z10*pal)
      g0=z1*x2*cu(51)
      v=cu(15)*(g2+cu(15)*g4)
      g2=cu(51)*z2*((cu(4)+cu(17)*x7*y1)/y4-cu(15))
      g4=cu(53)*x5*(cu(37)-z)/(x*y4)
      g0=g0+g2+g4
      epl=z11*g0-ee
      et=y3*(v-pal*g0)/pl
      if(kentr.eq.0) go to 76
      g4=(z2*x7+cu(55)*x5/x)*z9/x
      spl=y6*pap*g4-se
      st=(se+y6*(cu(37)*f4-pal*g4))/t
   76 go to 135
!
! *** quadratures taken with the gauss method (nz=5)
    5 nzr=5
      kkk=1
         kk1=1
      al3=alf**3
         z11=plm*al3
         x2=cu(15)*alf
      z10=z11/cu(47)
      kw=1
         ku=10
      go to 151
  169 z=(g1m-z10)/g1mp
      z3=1.d0
      if(z.lt.0.d0) z3=-1.d0
      z2=z*z3
      z1=pc
      if(pc.lt.0.d0) z1=-pc
      if(z1.lt.1.d0) go to 181
      z2=z2/z1
      if(z2.gt.0.3d0) z=0.3d0*z1*z3
  181 pc=pc-z
      if(z2.gt.eit) go to 151
      ku=15
         kw=2
      go to 151
  170 z=1.44059d0/(al3*alf)
      z1=z/pl
        z2=z/cu(17)
      pe=z2*gp
         ee=z1*ge
         hpr=cu(15)*g1/z10
      if(kentr.eq.0) go to 182
      z3=pc+alf
         y3=g3+g31+cu(49)*gp-z3*g1m
      se=z1*y3/t
  182 if(lst.eq.0) go to 183
      pap=g1m/g1mp
         x=g1a1-g1a
      pal=(cu(17)*g1m-alf*(x+cu(15)*g1p))/g1mp
      x1=g2p1-g2p
         y2=g2a+g2a1-cu(15)*g2p
      ppl=ei*cu(49)*x1/g1mp
      pt=pe*(cu(37)-(alf*y2+x1*pal)/gp)/t
      y1=g4p1-g4p-x2*g1p
      epl=ee*(pap*y1/ge-cu(4))
      et=ee*(cu(37)-(alf*(g4a1+g4a)+x2*(g1-g4p+alf*g1a-x2*g1p)+y1*pal)/ &
      ge)/t
      if(kentr.eq.0) go to 183
      y1=g3p1-g3p+cu(49)*x1-g1m-z3*g1mp
      spl=se*(pap*y1/y3-cu(4))
      st=g3a1+g3a+y2*cu(49)-g1m-z3*(g1a1-g1a+cu(15)*g1p)-cu(15)*g3p
      st=se*(cu(17)-(pal*y1+st*alf)/y3)/t
  183 go to 135
  151 kpg=0
  152 wo=0.d0
         w1=0.d0
         w2=0.d0
      wop=0.d0
         w1p=0.d0
        w2p=0.d0
      woa=0.d0
         w1a=0.d0
        w2a=0.d0
      if(pc.gt.pc2) go to 155
      if(kkk.eq.0) go to 158
      do 157 k=1,15
  157 uac(k+65)=sqrt(uac(k)+x2)
      kkk=0
  158 if(pc.gt.pc1) go to 156
      if(pc.lt.-4.4d1) go to 163
      x=exp(pc)
      do 161 k=1,ku
  161 uac(k+80)=uac(k+15)/(uac(k+30)*x+1.d0)
      do 162 k=1,5
      z=wk1(k)*wk4(k)
      wo=wo+z
      wop=wop+z/(1.d0+aio(k)*x)
      z=wk2(k)*wk5(k)
      w1=w1+z
      w1p=w1p+z/(1.d0+ai1(k)*x)
      if(kw.eq.1) go to 162
      z=wk6(k)*wk3(k)
      w2=w2+z
      if(lst.eq.0) go to 162
      woa=woa+wk4(k)/wk1(k)
      w1a=w1a+wk5(k)/wk2(k)
      w2p=w2p+z/(1.d0+ai2(k)*x)
      w2a=w2a+wk6(k)/wk3(k)
  162 continue
      wop=wop*x
         w1p=w1p*x
      wo=wo*x
        w1=w1*x
      if(kw.eq.1) go to 163
      w2=w2*x
      if(lst.eq.0) go to 163
      w1a=w1a*x
         woa=woa*x
      w2a=w2a*x
         w2p=w2p*x
  163 g1=w1+alf*wo
      g1p=w1p+alf*wop
      if(kw.eq.1) go to 164
      g2=w2+x2*w1
         g3=w2+alf*(w1+g1)
        g4=w2+alf*w1
      if(lst.eq.0) go to 164
      g1a=w1a+wo+alf*woa
      g2p=w2p+x2*w1p
         g3p=w2p+alf*(w1p+g1p)
         g4p=w2p+alf*w1p
      g2a=w2a+2.d0*w1+x2*w1a
         g3a=w2a+w1+g1+alf*(w1a+g1a)
      g4a=w2a+w1+alf*w1a
  164 if(kpg.eq.1) go to 166
      g11=g1
         g1p1=g1p
      if(kw.eq.1) go to 168
      g21=g2
         g31=g3
        g41=g4
      if(lst.eq.0) go to 168
      g1a1=g1a
         g2p1=g2p
         g3p1=g3p
         g4p1=g4p
      g2a1=g2a
         g3a1=g3a
         g4a1=g4a
  168 pc=-pc-x2
         kpg=1
      if(pc.gt.-4.4d1) go to 152
         g1=0.d0
         g2=0.d0
         g3=0.d0
         g4=0.d0
         g1p=0.d0
         g2p=0.d0
         g3p=0.d0
         g4p=0.d0
         g1a=0.d0
         g2a=0.d0
         g3a=0.d0
         g4a=0.d0
  166 pc=-pc-x2
      g1m=g11-g1
      g1mp=g1p1+g1p
      gp=g2+g21
      ge=g4+g41+x2*g1
  167 go to(169,170),kw
  155 do 171 k=1,5
      z4=xxi(k)-1.d0
      z1=exp(pc*z4)
      y1=pc*xxi(k)
      z2=1.d0+z1
      z3=x2+y1
      y2=pc*aai(k)*sqrt(pc*z3)/z2
      y4=cci(k)+pc
      z=y4+x2
      y6=bbi(k)*sqrt(y4*z)
      wo=wo+y2+y6
      if((lst.eq.0).and.(kw.ne.1)) go to 172
      z5=1./pc
         y3=0.5d0*xxi(k)/z3-z4*z1/z2+1.5d0*z5
      z6=1.d0/y4
         y5=0.5d0*(1.d0/z+z6)
      wop=wop+y2*y3+y6*y5
      z3=y2/z3
         z=y6/z
  172 y2=y2*y1
         y6=y6*y4
      w1=w1+y2+y6
      y3=y3+z5
        y5=y5+z6
      w1p=w1p+y2*y3+y6*y5
      if(kw.eq.1) go to 171
      if(lst.eq.0) go to 173
      woa=woa+z3+z
      z3=z3*y1
         z=z*y4
      w1a=w1a+z3+z
  173 y2=y2*y1
         y6=y6*y4
      w2=w2+y2+y6
      if(lst.eq.0) go to 171
      w2p=w2p+y2*(y3+z5)+y6*(y5+z6)
      w2a=w2a+z3*y1+z*y4
  171 continue
      go to 163
  156 if(kk1.eq.0) go to 191
      do 197 k=1,5
      wk6(k)=hhi(k)
         wk7(k)=ggi(k)
      wk4(k)=sqrt(vvi(k)+x2)
  197 wk5(k)=sqrt(zzi(k)+x2)
      do 195 k=1,24
  195 abcg(k)=0.d0
      k1=0
      do 198 i=1,3
      x=i
         x1=x-0.5d0
      do 190 k=1,5
      k2=k1+k
        k3=k2+65
      z=(x+ggsi(k))/pc2
         y=x1/zzi(k)
      z1=wk7(k)*(z-cpp(2))*cpp(4)
         z2=wk6(k)*(y-cpp(2))*cpp(4)
      z4=xxi(k)*wk7(k)/wk4(k)
        z5=wk6(k)/wk5(k)
         z6=cpp(5)*(z4+z5)
      asp(i)=asp(i)+abc(k2)*uac(k3)+z1*wk4(k)+z2*wk5(k)+z6*cpp(1)
      z8=cpp(5)*(z4/(vvi(k)+x2)+z5/(zzi(k)+x2))
      aspa(i)=aspa(i)+abc(k2)/uac(k3)+z1/wk4(k)+z2/wk5(k)-z8*cpp(1)
      k4=k2+15
      z1=wk7(k)*(cpp(3)-z)*cpp(1)
         z2=wk6(k)*(cpp(3)-y)*cpp(1)
      bsp(i)=bsp(i)+abc(k4)*uac(k3)+z1*wk4(k)+z2*wk5(k)-z6
      bspa(i)=bspa(i)+abc(k4)/uac(k3)+z1/wk4(k)+z2/wk5(k)+z8
      k4=k4+15
      csp(i)=csp(i)+abc(k4)*uac(k3)
      cspa(i)=cspa(i)+abc(k4)/uac(k3)
      k4=k4+15
      gsp(i)=gsp(i)+abc(k4)*uac(k3)
      gspa(i)=gspa(i)+abc(k4)/uac(k3)
      wk6(k)=wk6(k)*zzi(k)
  190 wk7(k)=wk7(k)*vvi(k)
  198 k1=k1+5
      kk1=0
  191 z=pc-pc1
         z1=2.d0*z
         z2=1.5d0*z
      wo=gsp(1)+z*(csp(1)+z*(bsp(1)+z*asp(1)))
      w1=gsp(2)+z*(csp(2)+z*(bsp(2)+z*asp(2)))
      wop=csp(1)+z1*(bsp(1)+z2*asp(1))
      w1p=csp(2)+z1*(bsp(2)+z2*asp(2))
      if(kw.eq.1) go to 163
      w2=gsp(3)+z*(csp(3)+z*(bsp(3)+z*asp(3)))
      if(lst.eq.0) go to 163
      w2p=csp(3)+z1*(bsp(3)+z2*asp(3))
      woa=gspa(1)+z*(cspa(1)+z*(bspa(1)+z*aspa(1)))
      w1a=gspa(2)+z*(cspa(2)+z*(bspa(2)+z*aspa(2)))
      w2a=gspa(3)+z*(cspa(3)+z*(bspa(3)+z*aspa(3)))
      go to 163
!
! *** relativistic asymptotics (nz=4)
    4 nzr=4
  520  ro=pl*cu(58)
          hi=1.d0/emue
          r1=ro*0.5d0*hi
          pit=pi2*al1
          pt2=pit*al1
          pa=pt2-1.5d0
      hu=psi*al1+1.d0
      do 525 it=1,4
          hu1=hu
          hu2=hu**2
          hu =2.d0*(hu2*hu+r1)/(3.d0*hu2+pa)
      if(abs(hu1/hu-1.d0).le.eit) go to 527
  525  continue
          r=r1**2+(pa*.3333333333d0)**3
      x=(r1+sqrt(r))**(1.d0/3.d0)
          hu=x-pa/(3.d0*x)
  527  continue
          hu2=hu**2
          psi=-7.77d2
          if(al1.gt.1.d-8) psi=(hu-1.d0)/al1
          pe=.25d0*(hu2**2+2.d0*pa*hu2+.46666666667d0*pt2**2-pt2)
          ee=(3.d0*pe+0.5d0*(3.d0*hu2+pt2))/ro-hi
          pe=pe+1.1550858d0
          ee=ee-1.1550858d0/ro
          if(lst.eq.0) go to 555
          r=1.d0/(3.d0*hu2+pa)
          r1=hu2+pa
          pt=pit*(hu2-0.5d0+0.46666666667d0*pt2)-2.d0*pit*hu2*r1*r
          r2=hi*hu*r
          ppl=r2*r1
          et=pit*(hu2*(1.d0-4.d0*pt2*r)+1.4d0*pt2-0.5d0)/ro
          epl=(3.d0*(ppl+r2)-ee-hi)
  555  ee=ee*cu(59)
       pe=pe*cu(60)
      if(kpar.ne.1) go to 558
      hpr=0.d0
      if(t.lt.5.96d0) go to 558
      eta=alf*hu
      if(eta.gt.6.d1) go to 558
      hu1=exp(-eta)
      al2=al1**2
      hpr=hu1*(1.2d1*al2-3.d0+hu1*((0.444444444d0*al2-1.d0) &
       *hu1-1.5d0*(al2-1.d0)))/(eta*r1)
  558 continue
          if(lst.eq.0)go to 556
          pt=pt*cu(61)
          ppl=ppl*cu(59)
          epl=epl*cu(59)
          et=et*cu(2)
 556   if(kentr.eq.0) go to 557
          y=cu(62)
      se=y*al1*(hu2+.466666666667d0*pt2-.5d0)/pl
          if(lst.eq.0) go to 557
      spl=-se+2.d0*pi2*cu(2)*al1*r2
          st=se/t+2.d0*y*pt2*(0.46666666667d0-2.d0*hu2*r)/(pl*cu(1))
  557  go to 135
!
! *** interpolation between perfect gas and expansion over
! *** half-integer f-d functions
   52 pl1=x*emue
         pl2=y*emue
        nzp1=1
         nzp2=2
!
! *** interpolation over density
   83 kin=ki
      lst1=lst
         lst=1
         kk=0
      dni=pl
      if(lst1.eq.0) go to 81
      kk=1
      tni=t
         t=t*(cu(4)+dst)
   81 pl=pl1
         nz=nzp1
        ki=2
      go to 90
   77 pn1=pe
         en1=ee
        sn1=se
         psn1=psi
         hprn1=hpr
      pnp1=ppl
         enp1=epl/pl
         snp1=spl/pl
      pl=pl2
         nz=nzp2
        ki=3
         nzr1=nzr
      go to 90
   78 if(kk.eq.2) go to 92
      wv=pl2-pl1
         wv1=cu(4)/wv
      wv4=cu(15)*wv1
         wv3=dni-pl1
         wv5=wv1*wv3
   92 x=epl/pl
      x1=pe-pn1
         x2=x1*wv4
         x3=x1*wv1
      z1=pnp1+ppl-x2
         z2=x3-z1-pnp1
         z1=wv1*z1
      pe=pn1+wv3*(pnp1+wv5*(z2+wv3*z1))
      x1=ee-en1
         x2=x1*wv4
         x3=x1*wv1
      y1=enp1+x-x2
         y2=x3-y1-enp1
         y1=wv1*y1
      ee=en1+wv3*(enp1+wv5*(y2+wv3*y1))
      if(kentr.eq.0) go to 91
      x1=se-sn1
         x2=x1*wv4
         x3=x1*wv1
      v1=snp1+spl/pl-x2
         v2=x3-v1-snp1
         v1=wv1*v1
      se=sn1+wv3*(snp1+wv5*(v2+wv3*v1))
   91 if(kk.eq.1) go to 82
      if(kk.eq.2) go to 124
  127 x1=(dni-pl1)*wv1
      psi=(psi-psn1)*x1+psn1
      hpr=(hpr-hprn1)*x1+hprn1
      nzr=10*nzr1+nzr
       pl=dni
        plm=pl/emue
      ki=kin
         lst=lst1
  134 nz=0
         pi=ei*pl
  135  go to(57,77,78,79,80,577,578), ki
   82 pn2=pe
         en2=ee
        sn2=se
      kk=2
         t=tni
      go to 81
  124 x1=t*dst
      pt=(pn2-pe)/x1
         et=(en2-ee)/x1
      if(kentr.eq.0) go to 126
      st=(sn2-se)/x1
      spl=snp1+wv5*(cu(15)*v2+cu(17)*wv3*v1)
      spl=dni*spl
  126 ppl=pnp1+wv5*(cu(15)*z2+cu(17)*wv3*z1)
      epl=enp1+wv5*(cu(15)*y2+cu(17)*wv3*y1)
      epl=dni*epl
      go to 127
!
! *** interpolation between degenerate gas and expansion over
! *** half-integer fermi-dirac functions
   62 pl1=x*emue
         pl2=y*emue
        nzp1=2
         nzp2=3
      go to 83
!
! *** interpolation over temperature
  128 kit=ki
      lst2=lst
         lst=1
         kkt=0
      dnt=t
      if(lst2.eq.0) go to 129
      kkt=1
         plni=pl
        pl=pl*(cu(4)+dst)
  129 t=t1
         nz=nz1
         ki=4
      go to 90
   79 pnt=pe
         ent=ee
         snt=se
         psnt=psi
         hprnt=hpr
         pntt=pt
         entt=et
         sntt=st
         t=t2
         nz=nz2
         ki=5
         nzrt=nzr
      go to 90
   80 if(kkt.eq.2) go to 93
      vw=t2-t1
         vw1=cu(4)/vw
      vw4=cu(15)*vw1
         vw3=dnt-t1
        vw5=vw1*vw3
   93 x1=pe-pnt
         x2=x1*vw4
         x3=x1*vw1
      z1=pntt+pt-x2
         z2=x3-z1-pntt
         z1=vw1*z1
      pe=pnt+vw3*(pntt+vw5*(z2+vw3*z1))
      x1=ee-ent
         x2=x1*vw4
         x3=x1*vw1
      y1=entt+et-x2
         y2=x3-y1-entt
         y1=vw1*y1
      ee=ent+vw3*(entt+vw5*(y2+vw3*y1))
      if(kentr.eq.0) go to 94
      x1=se-snt
         x2=x1*vw4
         x3=x1*vw1
      v1=sntt+st-x2
         v2=x3-v1-sntt
         v1=vw1*v1
      se=snt+vw3*(sntt+vw5*(v2+vw3*v1))
   94 if(kkt.eq.1) go to 130
      if(kkt.eq.2) go to 131
  133 x1=(dnt-t1)*vw1
      psi=(psi-psnt)*x1+psnt
      hpr=(hpr-hprnt)*x1+hprnt
      nzr=10*nzrt+nzr
      t=dnt
         alf=cu(1)/t
         al1=cu(4)/alf
         ei=t*rg
      eg=cu(3)*ei
         ki=kit
         lst=lst2
      go to 134
  130 pnt2=pe
        ent2=ee
         snt2=se
      kkt=2
         pl=plni
      go to 129
  131 x1=cu(4)/dst
      ppl=x1*(pnt2-pe)/pl
         epl=x1*(ent2-ee)
      if(kentr.eq.0) go to 132
      spl=(snt2-se)*x1
      st=sntt+vw5*(cu(15)*v2+cu(17)*vw3*v1)
  132 pt=pntt+vw5*(cu(15)*z2+cu(17)*vw3*z1)
      et=entt+vw5*(cu(15)*y2+cu(17)*vw3*y1)
      go to 133


! finish the pipeline loop
6969  continue
      enddo

      term_var(1) = den_row(1)
      term_var(2) = ptot_row(1)
      term_var(3) = etot_row(1)
      term_var(4) = stot_row(1)
      term_var(5) = temp_row(1)
      term_var(6) = cs_row(1)
      end







      subroutine pretty_eos_out(whose)
      include 'implno.dek'
      include 'vector_eos.dek'

! writes a pretty output for the eos tester


! declare the pass
      character*(*) whose


! local variables
      integer          j
      double precision ye,xcess,avo,kerg,xka
      parameter        (avo     = 6.0221417930d23, &
                        kerg    = 1.380650424d-16, &
                        xka = kerg*avo)


! popular formats
!01    format(1x,t2,a,t11,a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a)
!02    format(1x,t2,a,1p7e16.8)
!03    format(1x,t2,a7,1pe12.4,t22,a7,1pe12.4, &
!               t42,a7,1pe12.4,t62,a7,1pe12.4)
!04    format(1x,t2,a,t11,'total',t24,'ion',t34,'e- + e+', &
!             t58,'radiation',t70,'coulomb')
!05    format(1x,t2,a,1p3e12.4,t56,1p2e12.4)
!06    format(1x,t2,a,a,1pe12.4, &
!                t30,a,a,1pe12.4, &
!                t58,a,a,1pe12.4)



! loop over the pipeline
      do j=jlo_eos,jhi_eos


! the input
!      write(6,*)  ' '
!      write(6,03) 'temp  =',temp_row(j),'den   =',den_row(j), &
!                  'abar  =',abar_row(j),'zbar  =',zbar_row(j)

      ye = zbar_row(1)/abar_row(1)
      xcess = 1.0d0 - 2.0d0*ye
!      write(6,03) 'ye    =',ye,'xcess =',xcess
!      write(6,*) ' '


! and the output

!       write(6,01)  whose,'value','d/dd','d/dt','d/da','d/dz'

!       write(6,02) 'p tot=',ptot_row(j), &
!                    dpd_row(j),dpt_row(j),dpa_row(j),dpz_row(j)
!       write(6,02) 'p gas=',pgas_row(j), &
!                 dpgasd_row(j),dpgast_row(j),dpgasa_row(j),dpgasz_row(j)
!       write(6,02) 'p rad=',prad_row(j), &
!                 dpradd_row(j),dpradt_row(j),dprada_row(j),dpradz_row(j)
!       write(6,02) 'p ion=',pion_row(j), &
!                dpiond_row(j),dpiont_row(j),dpiona_row(j),dpionz_row(j)
!       write(6,02) 'p  e-=',pele_row(j), &
!                dpepd_row(j),dpept_row(j),dpepa_row(j),dpepz_row(j)
!       write(6,02) 'p  e+=',ppos_row(j)
!       write(6,02) 'p cou=',pcou_row(j), &
!                dpcoud_row(j),dpcout_row(j),dpcoua_row(j),dpcouz_row(j)


!       write(6,*)  ' '
!       write(6,02) 'e tot=',etot_row(j), &
!                    ded_row(j),det_row(j),dea_row(j),dez_row(j)
!       write(6,02) 'e gas=',egas_row(j), &
!                 degasd_row(j),degast_row(j),degasa_row(j),degasz_row(j)
!       write(6,02) 'e rad=',erad_row(j), &
!                deradd_row(j),deradt_row(j),derada_row(j),deradz_row(j)
!       write(6,02) 'e ion=',eion_row(j), &
!                deiond_row(j),deiont_row(j),deiona_row(j),deionz_row(j)
!       write(6,02) 'e  e-=',eele_row(j), &
!                deepd_row(j),deept_row(j),deepa_row(j),deepz_row(j)
!       write(6,02) 'e  e+=',epos_row(j)
!       write(6,02) 'e cou=',ecou_row(j), &
!                decoud_row(j),decout_row(j),decoua_row(j),decouz_row(j)

!       write(6,*)  ' '
!       write(6,02) 's tot=',stot_row(j), &
!                    dsd_row(j),dst_row(j),dsa_row(j),dsz_row(j)
!       write(6,02) 's/xka=',stot_row(j)/xka, &
!             dsd_row(j)/xka,dst_row(j)/xka,dsa_row(j)/xka,dsz_row(j)/xka
!       write(6,02) 's gas=',sgas_row(j), &
!                 dsgasd_row(j),dsgast_row(j),dsgasa_row(j),dsgasz_row(j)
!       write(6,02) 's rad=',srad_row(j), &
!                dsradd_row(j),dsradt_row(j),dsrada_row(j),dsradz_row(j)
!       write(6,02) 's ion=',sion_row(j), &
!                dsiond_row(j),dsiont_row(j),dsiona_row(j),dsionz_row(j)
!       write(6,02) 's  e-=',sele_row(j), &
!                dsepd_row(j),dsept_row(j),dsepa_row(j),dsepz_row(j)
!       write(6,02) 's  e+=',spos_row(j)
!       write(6,02) 's cou=',scou_row(j), &
!                dscoud_row(j),dscout_row(j),dscoua_row(j),dscouz_row(j)


! specific heats, and ratio of electostatic to thermal energy
! the 3 gammas and the sound speed for both the gas and the total
!       write(6,*)  ' '
!       write(6,02) 'cv  =',cv_row(j)/(kerg*avo)*abar_row(1), &
!                    dcvdd_row(j),dcvdt_row(j), &
!                    dcvda_row(j),dcvdz_row(j)
!       write(6,02) 'cp  =',cp_row(j), &
!                    dcpdd_row(j),dcpdt_row(j), &
!                    dcpda_row(j),dcpdz_row(j)
!       write(6,02) 'gam1=',gam1_row(j), &
!                    dgam1dd_row(j),dgam1dt_row(j), &
!                    dgam1da_row(j),dgam1dz_row(j)
!       write(6,02) 'gam2=',gam2_row(j), &
!                    dgam2dd_row(j),dgam2dt_row(j), &
!                    dgam2da_row(j),dgam2dz_row(j)
!       write(6,02) 'gam3=',gam3_row(j), &
!                    dgam3dd_row(j),dgam3dt_row(j), &
!                    dgam3da_row(j),dgam3dz_row(j)
!       write(6,02) 'cs  =',cs_row(j), &
!                    dcsdd_row(j),dcsdt_row(j), &
!                    dcsda_row(j),dcsdz_row(j)

!       write(6,*)  ' '
!       write(6,02) 'cvgas=',cv_gas_row(j)/(kerg*avo)*abar_row(1), &
!                    dcv_gasdd_row(j),dcv_gasdt_row(j), &
!                    dcv_gasda_row(j),dcv_gasdz_row(j)
!       write(6,02) 'cpgas=',cp_gas_row(j), &
!                    dcp_gasdd_row(j),dcp_gasdt_row(j), &
!                    dcp_gasda_row(j),dcp_gasdz_row(j)
!       write(6,02) 'g1gas=',gam1_gas_row(j), &
!                    dgam1_gasdd_row(j),dgam1_gasdt_row(j), &
!                    dgam1_gasda_row(j),dgam1_gasdz_row(j)
!       write(6,02) 'g2gas=',gam2_gas_row(j), &
!                    dgam2_gasdd_row(j),dgam2_gasdt_row(j), &
!                    dgam2_gasda_row(j),dgam2_gasdz_row(j)
!       write(6,02) 'g3gas=',gam3_gas_row(j), &
!                    dgam3_gasdd_row(j),dgam3_gasdt_row(j), &
!                    dgam3_gasda_row(j),dgam3_gasdz_row(j)
!       write(6,02) 'csgas=',cs_gas_row(j), &
!                    dcs_gasdd_row(j),dcs_gasdt_row(j), &
!                    dcs_gasda_row(j),dcs_gasdz_row(j)


! the thermodynamic consistency relations, these should all be
! at the floating point limit of zero
!       write(6,*) ' '
!       write(6,03) 'maxw1 =',dse_row(j),'maxw2 =',dpe_row(j), &
!                   'maxw3 =',dsp_row(j)

! number density of ions and its derivatives
!       write(6,03) 'xni   =',xni_row(j),  'xnim  =',xnim_row(j)
!       write(6,03) 'dxnidd=',dxned_row(j),'dxnidt=',dxnet_row(j), &
!                   'dxnida=',dxnea_row(j),'dxnidz=',dxnez_row(j)

! ion chemical potential and its derivatives
!       write(6,03) 'etaion=',etaion_row(j)
!       write(6,03) 'detaid=',detaid_row(j),'detait=',detait_row(j), &
!                   'detaia=',detaia_row(j),'detaiz=',detaiz_row(j)


! number density of electrons+positrons and its derivatives
!       write(6,03) 'xnele =',xne_row(j),'xnpos =',xnp_row(j), &
!                   'xnem  =',xnem_row(j)
!       write(6,03) 'dxnedd=',dxned_row(j),'dxnedt=',dxnet_row(j), &
!                   'dxneda=',dxnea_row(j),'dxnedz=',dxnez_row(j)


! electron chemical potential, positron chemical potential and its derivatives
!       write(6,03) 'etaele=',etaele_row(j),'etapos=',etapos_row(j)
!       write(6,03) 'detadd=',detad_row(j),'detadt=',detat_row(j), &
!                   'detada=',detaa_row(j),'detadz=',detaz_row(j)

!       write(6,03) 'zeff  =',zeff_row(j), &
!                   'ionzd =',zeff_row(j)/zbar_row(j), &
!                   'plasg =',plasg_row(j)

! end of pipeline loop
      enddo

      return
      end




      subroutine call_helmeos_DP(nrow, den, pres, abar, zbar, tguess, term_var)
      include 'implno.dek'
      include 'vector_eos.dek'

! tests the eos routine
!
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer, intent(in) :: nrow
      double precision :: term_var(6)
      double precision, intent(in), dimension(nrow) :: den, pres, abar, zbar, tguess
      !f2py INTEGER, INTENT(hide) :: nrow
      !f2py DOUBLE PRECISION, DIMENSION(nrow), INTENT(in) :: den, pres, abar, zbar, tguess
      double precision, dimension(nrow) :: rerr_P, rerr_T
      double precision, dimension(nrow) :: delta_P, delta_T
      double precision, dimension(nrow) :: T_lower, T_upper
      double precision, dimension(nrow) :: Pgoal
      logical, dimension(nrow) :: NR_converged

      double precision, parameter :: temp_floor = 1e4
      double precision, parameter :: rtol = 1e-6

      integer :: i, iter
      integer, parameter :: max_iter  = 100

! don't try and overfill the array
!      if (ninput.gt.nrowmax) then stop

! read the data table
!      call read_helm_table

! set the input vector. pipeline is only 1 element long

      jlo_eos = 1 ; jhi_eos = nrow
      do i = 1, nrow
         Pgoal(i) = pres(i)
         temp_row(i) = tguess(i)
         den_row(i)  = den(i)
         abar_row(i) = abar(i)
         zbar_row(i) = zbar(i)
      end do

      T_lower = temp_floor
      T_upper = 1e12
      NR_converged = .FALSE.

      ! do the NR iteration

      do iter = 1, max_iter

         call nados(term_var)

         do i = 1, nrow

            ! if this point is converged, go to the next one
            if (NR_converged(i)) cycle

            ! energy difference
            delta_P(i) = Pgoal(i) - ptot_row(i)

            ! keep things safe with bisect-limits
            if (delta_P(i).gt.0) then
               t_lower(i) = temp_row(i)
            else
               t_upper(i) = temp_row(i)
            end if

            ! update temperature
            delta_T(i) = delta_P(i) / dpt_row(i)
            temp_row(i) = temp_row(i) + delta_T(i)

            ! if this took us out of bounds, don't let it happen
            ! choose a new point inside the interval [t_lower, t_upper]
            ! the point is in the middle of the interval (logarthmically)

            if ((temp_row(i).gt.t_upper(i)).OR.(temp_row(i).lt.t_lower(i))) then
               temp_row(i) = sqrt(t_lower(i) * t_upper(i))
            end if

            ! calculate relative errors
            rerr_P(i) = delta_P(i) / Pgoal(i)
            rerr_T(i) = delta_T(i) / temp_row(i)

            ! if we're at tolerances, end this
            if ((abs(rerr_P(i)).LE.rtol).AND.(abs(rerr_T(i)).LE.rtol)) then
               NR_converged(i) = .TRUE.
            endif

            !allow points at the temperature floor to "converge"
            if (t_upper(i).le.temp_floor * (1d0 + rtol)) then
               NR_converged(i) = .TRUE.
               temp_row(i) = temp_floor
            end if

         end do

         if (ALL(NR_converged)) exit

      end do

      ! once more, with feeling
      NR_converged = .FALSE.

      call nados(term_var)

      end subroutine call_helmeos_DP

      subroutine call_helmeos_DE(nrow, den, ener, abar, zbar, tguess, term_var)
      include 'implno.dek'
      include 'vector_eos.dek'

! tests the eos routine
!
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer, intent(in) :: nrow
      double precision, intent(in), dimension(nrow) :: den, ener, abar, zbar, tguess
      double precision :: term_var(6)
      !f2py INTEGER, INTENT(hide) :: nrow
      !f2py DOUBLE PRECISION, DIMENSION(nrow), INTENT(in) :: den, ener, abar, zbar, tguess
      double precision, dimension(nrow) :: rerr_e, rerr_T
      double precision, dimension(nrow) :: delta_e, delta_T
      double precision, dimension(nrow) :: T_lower, T_upper
      double precision, dimension(nrow) :: egoal
      logical, dimension(nrow) :: NR_converged

      double precision, parameter :: temp_floor = 1e4
      double precision, parameter :: rtol = 1e-4

      integer :: i, iter
      integer, parameter :: max_iter  = 100

! don't try and overfill the array
!      if (ninput.gt.nrowmax) then stop

! read the data table
!      call read_helm_table

! set the input vector. pipeline is only 1 element long

      jlo_eos = 1 ; jhi_eos = nrow
      do i = 1, nrow
         egoal(i) = ener(i) / den(i) ! eos works on specific internal energy
         temp_row(i) = tguess(i)
         den_row(i)  = den(i)
         abar_row(i) = abar(i)
         zbar_row(i) = zbar(i)
      end do

      T_lower = temp_floor
      T_upper = 1e12
      NR_converged = .FALSE.

      ! do the NR iteration

      do iter = 1, max_iter

         call nados(term_var)

         do i = 1, nrow

            ! if this point is converged, go to the next one
            if (NR_converged(i)) cycle

            ! energy difference
            delta_E(i) = egoal(i) - etot_row(i)

            ! keep things safe with bisect-limits
            if (delta_E(i).gt.0) then
               t_lower(i) = temp_row(i)
            else
               t_upper(i) = temp_row(i)
            end if

            ! update temperature
            delta_T(i) = delta_E(i) / det_row(i)
            temp_row(i) = temp_row(i) + delta_T(i)

            ! if this took us out of bounds, don't let it happen
            ! choose a new point inside the interval [t_lower, t_upper]
            ! the point is in the middle of the interval (logarthmically)

            if ((temp_row(i).gt.t_upper(i)).OR.(temp_row(i).lt.t_lower(i))) then
               temp_row(i) = sqrt(t_lower(i) * t_upper(i))
            end if

            ! calculate relative errors
            rerr_e(i) = delta_e(i) / egoal(i)
            rerr_T(i) = delta_T(i) / temp_row(i)

            ! if we're at tolerances, end this
            if ((abs(rerr_e(i)).LE.rtol).AND.(abs(rerr_T(i)).LE.rtol)) then
               NR_converged(i) = .TRUE.
            endif

            !allow points at the temperature floor to "converge"
            if (t_upper(i).le.temp_floor * (1d0 + rtol)) then
               NR_converged(i) = .TRUE.
               temp_row(i) = temp_floor
            end if

         end do

         if (ALL(NR_converged)) exit

      end do

      ! once more, with feeling
      NR_converged = .FALSE.

      call nados(term_var)

    end subroutine call_helmeos_DE
