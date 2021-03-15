c	-----------------------------------------------------------------------------
c	Lemaitre 3D Hypoelastic-Plasticity 
c     with combined isotropic hardening of type:
c		POWER 
c		POWER_PIRES
c		POWER_LUDWIK 
c		SATURATION
c		DOUBLE_SATURATION
c		SWIFT
c		DISCRETE 
c     and Kinematic (nonlinear Armstrong Frederick type) Hardening
c	and Hill48 type initial anisotropy
c	-----------------------------------------------------------------------------
c
c	code development by C. Soyarslan
c	property of Department of Applied Mechanics in Forming Technologies
c	(formed in April 2011)
c	IUL TU Dortmund
c
c	-----------------------------------------------------------------------------
c
c	EXISTING IMPROVEMENTS
c	------+--------+-------------+----------------+-----------------------+-----+
c	Lemaitre with quasi unilateral effects
c	Shear factor introduced to the damage evolution
c	------+--------+-------------+----------------+-----------------------+-----+
c
c	------+--------+-------------+----------------+-----------------------+-----+
c     MODIFICATION LIST
c	------+--------+-------------+----------------+-----------------------+-----+
c
c	Done by Soyarslan, 21.AUG.2012:
c	 1) Overall theoretical check is done
c	 2) Code is optimized, compute Y is put in NR scheme, this way computed only
c		for the yielded elements
c	 3) Outputs are corrected, now triaxiality and max_shear_norm works
c	 4) A crucial correction for 3D code in use of vsprinc and a copy paste error
c		in computation of seffective is corrected, do not use older versions of 
c	    3D code
c	Done by Soyarslan, 6.AUG.2012:
c	 1) Isotropic hardening models are enriched to involve
c			POWER, POWER_PIRES, POWER_LUDWIK 
c			SATURATION, DOUBLE_SATURATION, SWIFT
c		    DISCRETE
c	Done by Soyarslan, 3.AUG.2012:
c	 1) Araf problem (see TODO.txt) is solved
c		if (yld_fnctn_NPO_KPO.lt.zero)then is replaced by
c		if (yld_fnctn_NPO_KPO.lt.TOLF)then is replaced by
c
c	-----------------------------------------------------------------------------
c
c	References:
c
c	See following internal reports for implementation details:
c	 1) FOSTA report internal 1
c	 2) FOSTA report internal 2
c
c	-----------------------------------------------------------------------------
c

c	TO DO LIST
c	------+--------+-------------+----------------+-----------------------+-----+
c      1) Clean code for the following warnings if still exist
c         a) D_0 initiation, sacma bir bicimde var...sdv1den aliyor hem de
c	    b) i_flow_index var...kullaniliyor mu?
c         c) compute_yld_function'da stress_eq_mises 2 checkinde 1 refer ediliyor ... 
c            triaxiality 2 icin giderken
c         d) compute dmg fonksiyonunda karsilastirma var ama EPSILON yok!!! Basa ekle!
c         e) Ayni yerde, no output flags
c         f) ABAQUSun principal bulucusu calisiyor gibi
c	 2) ...
c	------+--------+-------------+----------------+-----------------------+-----+
c
c	WARNING 1:	Should be used with double precision code,
c				otherwise will not work properly!
c	WARNING 2:	Due to iterative code segments with convergence check,
c	            do not use with parallel programming with multiple CPUS
c				If you have to do, supply required checks.	
c     WARNING 3:	code of POWER_LUDWIK isotropic hardening till now can not work 
c                 perfectly. need further modification
c	-----------------------------------------------------------------------------
c	UTILIZES NOTATION SIMO HUGHES PAGES 168,169,170
c
c     E STRAIN TENSOR       E11, E22, E33, 2*E12, 2*E23, 2*E13
c     A STRAIN TENSOR       A11, A22, A33, 2*A12, 2*A23, 2*A13
c
c     T STRESS TENSOR       S11, S22, S33,   S12,   S23,   S13
c     X STRESS TENSOR       X11, X22, X33,   X12,   X23,   X13
c
c     CONSTITUTIVE TENSOR   C1111, C1122, C1133, C1112, C1123, C1113
c                           C2211, C2222, C2233, C2212, C2223, C2213
c                           C3311, C3322, C3333, C3312, C3323, C3313
c                           C1211, C1222, C1233, C1212, C1223, C1213
c                           C2311, C2322, C2333, C2312, C2323, C2313
c                           C1311, C1322, C1333, C1312, C1323, C1313
c	-----------------------------------------------------------------------------

      subroutine vumat (
c	Read only -
     1     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
c     Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew )
c
      include 'vaba_param.inc'
c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname

      parameter (
     1     i_umt_nblock = 1,
     2     i_umt_npt    = 2,
     3     i_umt_layer  = 3,
     4     i_umt_kspt   = 4,
     5     i_umt_noel   = 5 )

      call  vumatXtrArg ( jblock(i_umt_nblock),
     1     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
     7     stressNew, stateNew, enerInternNew, enerInelasNew,
     8     jblock(i_umt_noel), jblock(i_umt_npt),
     9     jblock(i_umt_layer), jblock(i_umt_kspt))

      return
      end

c	*****************************************************************************
c	*****************************************************************************

	subroutine vumatXtrArg (
c	read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, timeinc, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
c	write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew,
c	read only extra arguments -
     6  nElement, nMatPoint, nLayer, nSecPoint)
c
      include 'vaba_param.inc'
c	
c	all arrays dimensioned by (*) are not used in this algorithm
	dimension props(nprops), density(nblock), 
	1  strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), defgradOld(nblock,9),
     4  stressOld(nblock,ndir+nshr),
     5  stateOld(nblock,nstatev), enerInternOld(nblock),
     6  enerInelasOld(nblock),
     7  stretchNew(nblock,ndir+nshr), defgradNew(nblock,9), 
     8  stressNew(nblock,ndir+nshr) 
     
      dimension enerInelasNew(nblock),stateNew(nblock,nstatev),
	1  enerInternNew(nblock)

	dimension nElement(nblock),nMatPoint(nblock),nLayer(nblock),
	1  nSecPoint(nblock)

c	new set for 3D applications
c	-----------------------------------------------------------------------------
      dimension delta_E_v6(6)
      
c	stress like (conjugate) internal variables (Ep, A)
c	-----------------------------------------------------------------------------
	dimension delta_stress_eff_tri_v6(6)

	dimension stress_hom_N_v6(6),stress_hom_NPO_v6(6)
	dimension delta_stress_hom_tri_v6(6),stress_hom_tri_v6(6)

	dimension X_N_v6(6),X_NPO_v6(6) ! X_tri_v6=X_N_v6, no need to use additional memory

	dimension C_e_eff_m66(6,6)
	
	character*80 cmname

c	additional tensors added for Simo-Hughes approach	
c	-----------------------------------------------------------------------------
	COMMON /constants1/  zero,one,two,three,four
	COMMON /constants2/  one_third,half,two_thirds,three_halves
	COMMON /constants3/  sqrt_two_thirds, sqrt_three_halves
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	COMMON /tolerances/  ALF,TOLX,TOLF,TOLMIN,EPSILON,MAXITS
	COMMON /HILL/ x_F,x_G,x_H,x_M,x_N,x_L
	COMMON /LEMAITRE/ Beta, Big_S, Small_s, Y_0, D_cri, D_del
	COMMON /KINHARD/ C_X, Q_X ! Armstrong Frederic type nonlinear kinematic hardening
	COMMON /ISOHARD/ C_q, Q_q 
	COMMON /VISCOSITY/ ETA

c	IF SOME VALUES ARE DIFFERENT FOR SOME ELEMENTS PUT IN NBLOCK CYCLE
c	-----------------------------------------------------------------------------
	call initialize_constants()
	call initialize_computational_parameters()
	call initialize_analysis_flags()
	call select_analysis_flags()
	call compute_elastic_constants(props)
	call compute_inelastic_constants(props) ! props is not "actually" used
	
c	-----------------------------------------------------------------------------
c     *INITIAL CONDITIONS, TYPE=SOLUTION
c      Set-1,1.0,0.,0.,0.,0.,0.,0., ! First one: Element deletion, Second one: Initial Damage, ! 7 vars per this line
c      0.,0.,0.,0.,0.,0.,0., 0.,    ! 8 vars per this line
c	-----------------------------------------------------------------------------
c	same privilaged directions for all elements
c	-----------------------------------------------------------------------------
c     definition of orthogonal priviliged directions
c	CHECK OK! added to dim list
c	-----------------------------------------------------------------------------
	
c	make computations for every material point
c	-----------------------------------------------------------------------------
      do 1000 nblck = 1,nblock
	
c	computation of elastic compliance tensor
c	-----------------------------------------------------------------------------
	call compute_eff_elastic_stiffness_matrix_6_by_6(C_e_eff_m66)
	
	Nelem=nElement(nblck)

c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "STEP TIME..............:",stepTime
	write(6,*) "TOTAL TIME.............:",totalTime
	write(6,*) "TIME INCREMENT.........:",timeinc
	write(6,*) ""
	write(6,*) "ELEMENT................:",Nelem
	write(6,*) ""
	endif
c	-----------------------------------------------------------------------------

c     NOTE THAT VUMAT CONVENTION IS DIFFERENT AND IT PASSES CONTINUUM COMPONENTS
c     YOU SHOULD CORRECT THIS OPTION WHILE CONVERTING TO UMAT
c	-----------------------------------------------------------------------------
      delta_E_v6(1)=    strainInc(nblck,1)
      delta_E_v6(2)=    strainInc(nblck,2)
      delta_E_v6(3)=    strainInc(nblck,3)
      delta_E_v6(4)=two*strainInc(nblck,4) !the shear 12 strain component of VUMAT is continuum component!
      delta_E_v6(5)=two*strainInc(nblck,5)
      delta_E_v6(6)=two*strainInc(nblck,6)
       
      stress_hom_N_v6(1)=stressOld(nblck,1)
      stress_hom_N_v6(2)=stressOld(nblck,2)      
      stress_hom_N_v6(3)=stressOld(nblck,3)
      stress_hom_N_v6(4)=stressOld(nblck,4)      
      stress_hom_N_v6(5)=stressOld(nblck,5)      
      stress_hom_N_v6(6)=stressOld(nblck,6) 
      
c	take damage from initial step  
c	-----------------------------------------------------------------------------
c     STATE VAR 1 IS FOR ELEMENT DELETION
c	-----------------------------------------------------------------------------
      D_N     =stateOld(nblck,2)           
	alpha_N =stateOld(nblck,3)
	      
c	state variables of the problem, 
c	-----------------------------------------------------------------------------
	X_N_v6(1)  =stateOld(nblck,4)
	X_N_v6(2)  =stateOld(nblck,5)
	X_N_v6(3)  =stateOld(nblck,6)
	X_N_v6(4)  =stateOld(nblck,7)
	X_N_v6(5)  =stateOld(nblck,8)
	X_N_v6(6)  =stateOld(nblck,9)
	
c	state variables of the problem, 
c	-----------------------------------------------------------------------------
      i_flow_index_tracker    =stateOld(nblck,10)
		
c	compute trial stress increment
c	-----------------------------------------------------------------------------
	call compute_stress_from_strain(delta_E_v6,C_e_eff_m66,
	1                                delta_stress_eff_tri_v6) 

      x_coef=(one-D_N)
      
      delta_stress_hom_tri_v6(1)=x_coef*delta_stress_eff_tri_v6(1)
      delta_stress_hom_tri_v6(2)=x_coef*delta_stress_eff_tri_v6(2)
      delta_stress_hom_tri_v6(3)=x_coef*delta_stress_eff_tri_v6(3)
      delta_stress_hom_tri_v6(4)=x_coef*delta_stress_eff_tri_v6(4)
      delta_stress_hom_tri_v6(5)=x_coef*delta_stress_eff_tri_v6(5)
      delta_stress_hom_tri_v6(6)=x_coef*delta_stress_eff_tri_v6(6)
      
      stress_hom_tri_v6(1)=stress_hom_N_v6(1)+delta_stress_hom_tri_v6(1)
	stress_hom_tri_v6(2)=stress_hom_N_v6(2)+delta_stress_hom_tri_v6(2)
      stress_hom_tri_v6(3)=stress_hom_N_v6(3)+delta_stress_hom_tri_v6(3)
	stress_hom_tri_v6(4)=stress_hom_N_v6(4)+delta_stress_hom_tri_v6(4)
      stress_hom_tri_v6(5)=stress_hom_N_v6(5)+delta_stress_hom_tri_v6(5)
	stress_hom_tri_v6(6)=stress_hom_N_v6(6)+delta_stress_hom_tri_v6(6)
      
c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "BEFORE THE CORE:"
	write(6,*) ""
	write(6,*) "stress_hom_tri_v6(1)...:",stress_hom_tri_v6(1)
	write(6,*) "stress_hom_tri_v6(2)...:",stress_hom_tri_v6(2)
	write(6,*) "stress_hom_tri_v6(3)...:",stress_hom_tri_v6(3)
	write(6,*) "stress_hom_tri_v6(4)...:",stress_hom_tri_v6(4)
	write(6,*) "stress_hom_tri_v6(5)...:",stress_hom_tri_v6(5)
	write(6,*) "stress_hom_tri_v6(6)...:",stress_hom_tri_v6(6)
	write(6,*) ""
	endif
c	-----------------------------------------------------------------------------

c	call the return mapping algorithm
c	-----------------------------------------------------------------------------
	call return_map(
c     INPUTS	 
c	-----------------------------------------------------------------------------
     1                stepTime,i_flow_index_tracker,
	2                stress_hom_tri_v6,
	3                alpha_N,D_N,
     5                X_N_v6,
	4                C_e_eff_m66,
c     OUTPUTS	 
c	-----------------------------------------------------------------------------
     8	 		    stress_hom_NPO_v6,
     9                x_gamma_NPO,
     1                alpha_NPO,D_NPO,
     5                X_NPO_v6,
	6                delta_alpha,
	1                delta_D,
     3	  			NUMITERATIONS,x_NORM_R,i_delete,
     5                Y_NPO,triaxiality,shear_max_norm)
c	-----------------------------------------------------------------------------
c	call the return mapping algorithm over
c	-----------------------------------------------------------------------------
      
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "OUT OF THE CORE:"
	write(6,*) ""
	write(6,*) "alpha_NPO..............:",alpha_NPO
	write(6,*) "D_NPO..................:",D_NPO
	write(6,*) ""
	write(6,*) "stress_hom_NPO_v6(1)...:",stress_hom_NPO_v6(1)
	write(6,*) "stress_hom_NPO_v6(2)...:",stress_hom_NPO_v6(2)
	write(6,*) "stress_hom_NPO_v6(3)...:",stress_hom_NPO_v6(3)
	write(6,*) "stress_hom_NPO_v6(4)...:",stress_hom_NPO_v6(4)
	write(6,*) "stress_hom_NPO_v6(5)...:",stress_hom_NPO_v6(5)
	write(6,*) "stress_hom_NPO_v6(6)...:",stress_hom_NPO_v6(6)
	write(6,*) ""
	endif

c	update stresses
c	-----------------------------------------------------------------------------
	stressNew(nblck,1)=stress_hom_NPO_v6(1)
	stressNew(nblck,2)=stress_hom_NPO_v6(2)
	stressNew(nblck,3)=stress_hom_NPO_v6(3)
	stressNew(nblck,4)=stress_hom_NPO_v6(4)
	stressNew(nblck,5)=stress_hom_NPO_v6(5)
	stressNew(nblck,6)=stress_hom_NPO_v6(6)

c	update other variables
c	-----------------------------------------------------------------------------
	stateNew(nblck,2) =D_NPO
	stateNew(nblck,3) =alpha_NPO
	stateNew(nblck,4) =X_NPO_v6(1)
	stateNew(nblck,5) =X_NPO_v6(2)
	stateNew(nblck,6) =X_NPO_v6(3)
	stateNew(nblck,7) =X_NPO_v6(4)
	stateNew(nblck,8) =X_NPO_v6(5)
	stateNew(nblck,9) =X_NPO_v6(6)
	stateNew(nblck,10)=i_flow_index_tracker
	stateNew(nblck,11)=delta_alpha
	stateNew(nblck,12)=delta_D
	stateNew(nblck,13)=NUMITERATIONS
	stateNew(nblck,14)=x_NORM_R
	stateNew(nblck,15)=Y_NPO 
	stateNew(nblck,16)=triaxiality 
	stateNew(nblck,17)=shear_max_norm
	stateNew(nblck,18)=0.0d0
	stateNew(nblck,19)=0.0d0
	stateNew(nblck,20)=0.0d0

	stateNew(nblck,1)=1.0d0  ! DELETION VARIABLE
	if(stateNew(nblck,2).gt.D_del) stateNew(nblck,1)=0.0d0

1000  continue

      return
      end

c	*****************************************************************************
c	*****************************************************************************

	subroutine compute_stress_from_strain(e,C_e,s) 
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
	include 'vaba_param.inc'

	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME

	dimension s(6),e(6),C_e(6,6)

c	compute trial stress increment
c	-----------------------------------------------------------------------------
      s(1)=C_e(1,1)*e(1)+C_e(1,2)*e(2)+C_e(1,3)*e(3)
      s(2)=C_e(2,1)*e(1)+C_e(2,2)*e(2)+C_e(2,3)*e(3)
      s(3)=C_e(3,1)*e(1)+C_e(3,2)*e(2)+C_e(3,3)*e(3)
      s(4)=C_e(4,4)*e(4)
      s(5)=C_e(5,5)*e(5)
      s(6)=C_e(6,6)*e(6)

	return
	end

c	*****************************************************************************
c	*****************************************************************************
	
	subroutine compute_K_alpha_power
	1	                             (alpha,x_K_alpha,x_K_alpha_prime)
c	-----------------------------------------------------------------------------
c	last update : 06 AUGUST 2012, 14:58, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
	include 'vaba_param.inc'

	COMMON /constants1/ zero,one,two,three,four
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /plasticity_power     /  x_C,x_pow,x_add
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	
c	K' (alpha)=x_C*alpha^x_pow
c	K''(alpha)=x_pow*x_C*alpha^(x_pow-1)
c	-----------------------------------------------------------------------------
c	PREVENT ZERO DIVISION, which may occur for x_pow<1 at first step
c	-----------------------------------------------------------------------------
	x_K_alpha      =x_C*alpha**x_pow		
	x_K_alpha_prime=x_pow*x_C*alpha**(x_pow-one)

	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "Inside POWER:"
	write(6,*) ""
	write(6,*) "alpha is.............:",alpha
	write(6,*) "x_K_alpha is.........:",x_K_alpha
	write(6,*) "x_K_alpha_prime is...:",x_K_alpha_prime
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif

      return
	end

c	*****************************************************************************
c	*****************************************************************************
	
	subroutine compute_K_alpha_power_ludwik
	1	                             (alpha,x_K_alpha,x_K_alpha_prime)

c	-----------------------------------------------------------------------------
c	last update : 06 AUGUST 2012, 14:58, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
	include 'vaba_param.inc'

	COMMON /constants1/ zero,one,two,three,four
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /plasticity_power     /  x_C,x_pow,x_add
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	
c	K' (alpha)=x_C*alpha^x_pow-sigma0
c	K''(alpha)=x_pow*x_C*alpha^(x_pow-1)
c	-----------------------------------------------------------------------------
c	PREVENT ZERO DIVISION, which may occur for x_pow<1 at first step
c	-----------------------------------------------------------------------------
	x_K_alpha      =x_C*alpha**x_pow-sigma_y_0_11	
	x_K_alpha_prime=x_pow*x_C*alpha**(x_pow-one)

	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "Inside POWER LUDWIK:"
	write(6,*) ""
	write(6,*) "alpha is.............:",alpha
	write(6,*) "x_K_alpha is.........:",x_K_alpha
	write(6,*) "x_K_alpha_prime is...:",x_K_alpha_prime
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif

      return
	end

c	*****************************************************************************
c	*****************************************************************************
	
	subroutine compute_K_alpha_power_pires
	1	                             (alpha,x_K_alpha,x_K_alpha_prime)

c	-----------------------------------------------------------------------------
c	last update : 06 AUGUST 2012, 14:58, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
	include 'vaba_param.inc'

	COMMON /constants1/ zero,one,two,three,four
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /plasticity_power     /  x_C,x_pow,x_add
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	
	x_K_alpha      =x_C*(alpha+x_add)**x_pow		
	x_K_alpha_prime=x_pow*x_C*(alpha+x_add)**(x_pow-one)

	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "Inside POWER PIRES:"
	write(6,*) ""
	write(6,*) "alpha is.............:",alpha
	write(6,*) "x_K_alpha is.........:",x_K_alpha
	write(6,*) "x_K_alpha_prime is...:",x_K_alpha_prime
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif

      return
	end

c	*****************************************************************************
c	*****************************************************************************
	
	subroutine compute_K_alpha_saturation
     1                                 (alpha,x_K_alpha,x_K_alpha_prime) 
c	-----------------------------------------------------------------------------
c	last update : 06 AUGUST 2012, 14:58, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
	include 'vaba_param.inc'

	COMMON /constants1/ zero,one,two,three,four
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /plasticity_saturation/  x_K,x_K_inf,x_K_0,x_delta
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	
	x_K_alpha=x_K*alpha+(x_K_inf-x_K_0)*(one-exp(-x_delta*alpha))
	x_K_alpha_prime=x_K+x_delta*(x_K_inf-x_K_0)*exp(-x_delta*alpha)

	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "Inside SATURATION:"
	write(6,*) ""
	write(6,*) "alpha is.............:",alpha
	write(6,*) "x_K_alpha is.........:",x_K_alpha
	write(6,*) "x_K_alpha_prime is...:",x_K_alpha_prime
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif

      return
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine compute_K_alpha_double_saturation
     1                                 (alpha,x_K_alpha,x_K_alpha_prime) 
c	-----------------------------------------------------------------------------
c	last update : 06 AUGUST 2012, 14:58, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
	include 'vaba_param.inc'
c
	COMMON /constants1/ zero,one,two,three,four
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /plasticity_saturation/  x_K,x_K_inf,x_K_0,x_delta
	COMMON /plasticity_saturation_double/  x_K_inf1,x_K_inf2,
	1                                       x_delta1,x_delta2
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
c	

	call precision_exp(-x_delta1*alpha,exp_minus_x_delta_alpha1)
      call precision_exp(-x_delta2*alpha,exp_minus_x_delta_alpha2)
c     
	x_K_alpha=x_K*alpha+(x_K_inf1-x_K_0)*(one-exp_minus_x_delta_alpha1)
     1	               +(x_K_inf2-x_K_0)*(one-exp_minus_x_delta_alpha2)
	x_K_alpha_prime=x_K
	1               +x_delta1*(x_K_inf1-x_K_0)*exp_minus_x_delta_alpha1
     2               +x_delta2*(x_K_inf2-x_K_0)*exp_minus_x_delta_alpha2

	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "Inside DOUBLE SATURATION:"
	write(6,*) ""
	write(6,*) "alpha is.............:",alpha
	write(6,*) "x_K_alpha is.........:",x_K_alpha
	write(6,*) "x_K_alpha_prime is...:",x_K_alpha_prime
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif
 
      return
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine compute_K_alpha_swift
     1                                 (alpha,x_K_alpha,x_K_alpha_prime) 
c	-----------------------------------------------------------------------------
c	last update : 06 AUGUST 2012, 14:58, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
	include 'vaba_param.inc'
c
	COMMON /constants1/ zero,one,two,three,four
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /plasticity_swift/ x_K_a,x_K_b,x_K_n
	COMMON /Lueders/ sigma_y_1_11,x_alpha_offset,x_K_offset
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
c	
	
      if(alpha.lt.(x_alpha_offset))then
        x_K_alpha=alpha*x_K_offset
        x_K_alpha_prime=x_K_offset
      return
      endif
      
      x_K_alpha=x_K_a*((x_K_b+alpha)**x_K_n)-sigma_y_0_11        
     	x_K_alpha_prime=x_K_a*x_K_n*((x_K_b+alpha)**(x_K_n-one))

	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "Inside SWIFT:"
	write(6,*) ""
	write(6,*) "alpha is.............:",alpha
	write(6,*) "x_alpha_offset is....:",x_alpha_offset
	write(6,*) "x_K_alpha is.........:",x_K_alpha
	write(6,*) "x_K_alpha_prime is...:",x_K_alpha_prime
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif
 
      return
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine compute_K_alpha_discrete(index_tracker,alpha,
     2                                    x_K_alpha,x_K_alpha_prime) 
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
	include 'vaba_param.inc'

	COMMON /constants1/ zero,one,two,three,four
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	COMMON /FLOW_CURVE/flow_curve(66),i_dim_flow_curve
	COMMON /tolerances/  ALF,TOLX,TOLF,TOLMIN,EPSILON,MAXITS	

      if(alpha.lt.EPSILON)then	
        x_K_alpha   	=zero
        x_K_alpha_prime	=one
        x_K_alpha_prime	=zero
        return
      endif

      i_flow_curve_ended=1

      xnom              =zero
      xdenom            =zero
      dstress           =zero
      stress_inc        =zero

      x_K_alpha   	    =zero
      x_K_alpha_prime	=one

	index_tracker=1

      do i=index_tracker,i_dim_flow_curve
        if(flow_curve(2*i).gt.alpha) then
		    i_flow_curve_ended=0 ! means flow curve has not come to an end
		    goto 10 
	  endif
      enddo
      
10    index_tracker=i-1

      if(i_flow_curve_ended.ne.1)then
      
           xdstrain  =alpha-flow_curve(2*i-2)
           dstrain	 =flow_curve(2*i)-flow_curve(2*i-2)
           dstress	 =flow_curve(2*i-1)-flow_curve(2*i-3)
           slope	 =dstress/dstrain
          
           stress_inc	   =xdstrain*slope
           x_K_alpha	   =flow_curve(2*i-3)+stress_inc-sigma_y_0_11
           x_K_alpha_prime =slope
      
      else
      
           x_K_alpha=flow_curve(2*i-3)-sigma_y_0_11
           x_K_alpha_prime=zero
          
      endif
      
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "Inside DISCRETE:"
	write(6,*) ""
	write(6,*) "index_tracker     ",index_tracker
	write(6,*) "i_dim_flow_curve  ",i_dim_flow_curve
	write(6,*) "alpha             ",alpha
	write(6,*) ""
	write(6,*) "xnom              ",xnom
	write(6,*) "xdenom            ",xdenom
	write(6,*) "dstress           ",dstress
	write(6,*) ""
	write(6,*) "x_K_alpha         ",x_K_alpha
	write(6,*) "x_K_alpha_prime   ",x_K_alpha_prime
	write(6,*) ""
	write(6,*) "*****************************************************"	
	write(6,*) "*****************************************************"
	endif
	

	return
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine swap_these(a,b)
c	-----------------------------------------------------------------------------
c	last update   : 08 AUGUST 2011, 16:34, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 08 AUGUST 2011
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'
	                  
      dummy=a
      a=b
      b=dummy

	return
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine order_my_eigenvalues(s_eff_1,s_eff_2,s_eff_3)
c	-----------------------------------------------------------------------------
c	last update   : 08 AUGUST 2011, 16:34, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 08 AUGUST 2011
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'
	                  
      if(s_eff_2.gt.s_eff_1) call swap_these(s_eff_1,s_eff_2)
      if(s_eff_3.gt.s_eff_2) call swap_these(s_eff_2,s_eff_3)
      if(s_eff_2.gt.s_eff_1) call swap_these(s_eff_1,s_eff_2)

	return
	end

c	*****************************************************************************
c	*****************************************************************************
          
      subroutine compute_damage_fnctn(stress_hom_v6,D,dDdot_dgamma,
	1                                Y,triaxiality,shear_max_norm)
c	-----------------------------------------------------------------------------
c	last update   : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

      dimension stress_hom_v6(6)
      dimension stress_eff_v6(6)
      
      dimension eigVal(1,3) 
      
c	IF MY EIGENVALUE CODE IS USED
c	-----------------------------------------------------------------------------
c     dimension eigValue(3)     ! vector containing the eigenvalues of b
c     dimension princ(3,3)      ! matrix containing the three principal row vectors, not used here           
c	-----------------------------------------------------------------------------
      
      COMMON /constants1/  zero,one,two,three,four     
	COMMON /constants2/  one_third,half,two_thirds,three_halves  
	COMMON /LEMAITRE/ Beta, Big_S, Small_s, Y_0, D_cri, D_del	     
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /CRACK_CLOSURE/ h_crack_closure		
	COMMON /tolerances/  ALF,TOLX,TOLF,TOLMIN,EPSILON,MAXITS	
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME	
	COMMON /MAX_SHEAR_POWER/ shear_max_power
      
c     compute effective components of stress tensor      
c	-----------------------------------------------------------------------------   
      x_coef=one/(one-D)
      
      stress_eff_v6(1)=x_coef*stress_hom_v6(1)
      stress_eff_v6(2)=x_coef*stress_hom_v6(2)
      stress_eff_v6(3)=x_coef*stress_hom_v6(3)
      stress_eff_v6(4)=x_coef*stress_hom_v6(4)
      stress_eff_v6(5)=x_coef*stress_hom_v6(5)
      stress_eff_v6(6)=x_coef*stress_hom_v6(6)
          
      p_eff=(stress_eff_v6(1)
     1      +stress_eff_v6(2)
     2      +stress_eff_v6(3))/three      
      
      s_dev_eff_v6_1=stress_eff_v6(1)-p_eff
      s_dev_eff_v6_2=stress_eff_v6(2)-p_eff
      s_dev_eff_v6_3=stress_eff_v6(3)-p_eff
      s_dev_eff_v6_4=stress_eff_v6(4)
      s_dev_eff_v6_5=stress_eff_v6(5)
      s_dev_eff_v6_6=stress_eff_v6(6)
      
c     -----------------------------------------------------------------------------
      seq_eff_Mises=s_dev_eff_v6_1*s_dev_eff_v6_1
     1             +s_dev_eff_v6_2*s_dev_eff_v6_2
     2             +s_dev_eff_v6_3*s_dev_eff_v6_3
     3         +two*s_dev_eff_v6_4*s_dev_eff_v6_4
     4         +two*s_dev_eff_v6_5*s_dev_eff_v6_5
     5         +two*s_dev_eff_v6_6*s_dev_eff_v6_6
      seq_eff_Mises=sqrt(three_halves*seq_eff_Mises)

      if(abs(seq_eff_Mises).lt.EPSILON)seq_eff_Mises=EPSILON

      triaxiality=p_eff/seq_eff_Mises ! for output

c    -----------------------------------------------------------------------------
c    FIND EIGENVALUES OF THE STRESS TENSOR
c    -----------------------------------------------------------------------------
c     call vsprinc(nblock,stress_eff_v6,eigVal,ndir,nshr)
      call vsprinc(1,stress_eff_v6, eigVal, 3, 3)
c    -----------------------------------------------------------------------------

c    MY CODE TO FIND EIGENVALUES OF THE STRESS TENSOR, JUST IN CASE!!!
c    -----------------------------------------------------------------------------
c    subroutine jacobi(stress_eff_v6,eigValue,princ)
c    -----------------------------------------------------------------------------
c     Evaluates the principal values and principal directions given any symmetric
c	second order b tensor using the Jacobi iteration. Adapted from numerical recipies
c
c     stress_eff_v6(6)-->  any second order tensor in vector form
c     eigValue(3)     -->  vector containing the eigenvalues of b
c     princ(3,3)      -->  matrix containing the three principal row vectors, not used here
c    -----------------------------------------------------------------------------

c	although no need to order, we order them for possible future use
c    -----------------------------------------------------------------------------

c     if((eig1.gt.eig2.gt.eig3))then
      eig1=eigVal(1,1)
      eig2=eigVal(1,2)
      eig3=eigVal(1,3)

c	IF MY EIGENVALUE CODE IS USED
c	-----------------------------------------------------------------------------
c     eig1=eigValue(1)
c     eig2=eigValue(2)
c     eig3=eigValue(3)      
c	-----------------------------------------------------------------------------
        
      s_eff_1=eig1
      s_eff_2=eig2
      s_eff_3=eig3
                  
c    -----------------------------------------------------------------------------
c	ordering is needed to decide on the maximum shear stress
c    -----------------------------------------------------------------------------
	call order_my_eigenvalues(s_eff_1,s_eff_2,s_eff_3)           

c     -----------------------------------------------------------------------------
	s_eff_1_PLUS = s_eff_1
	s_eff_2_PLUS = s_eff_2
	s_eff_3_PLUS = s_eff_3

	s_eff_1_MINUS=-s_eff_1
	s_eff_2_MINUS=-s_eff_2
	s_eff_3_MINUS=-s_eff_3

      if(s_eff_1_PLUS.lt.zero) s_eff_1_PLUS  =zero
      if(s_eff_2_PLUS.lt.zero) s_eff_2_PLUS  =zero
      if(s_eff_3_PLUS.lt.zero) s_eff_3_PLUS  =zero

      if(s_eff_1_MINUS.lt.zero)s_eff_1_MINUS =zero
      if(s_eff_2_MINUS.lt.zero)s_eff_2_MINUS =zero
      if(s_eff_3_MINUS.lt.zero)s_eff_3_MINUS =zero
      
      p_eff_PLUS = p_eff
      p_eff_MINUS=-p_eff

      if(p_eff_PLUS .lt.zero) p_eff_PLUS =zero
      if(p_eff_MINUS.lt.zero) p_eff_MINUS=zero     

	Y_PLUS=(one+poisson)/(two*E)*(s_eff_1_PLUS*s_eff_1_PLUS
	1                             +s_eff_2_PLUS*s_eff_2_PLUS
	2                             +s_eff_3_PLUS*s_eff_3_PLUS)
	3	  -(9.0d0*poisson)/(two*E)*p_eff_PLUS*p_eff_PLUS

	Y_MINUS=(one+poisson)/(two*E)*(s_eff_1_MINUS*s_eff_1_MINUS
	1                              +s_eff_2_MINUS*s_eff_2_MINUS
	2                              +s_eff_3_MINUS*s_eff_3_MINUS)
	3	  -(9.0d0*poisson)/(two*E)*p_eff_MINUS*p_eff_MINUS

c	-----------------------------------------------------------------------------     
      Y=Y_PLUS+Y_MINUS*h_crack_closure   ! for quasi-unilateral damage growth
c	-----------------------------------------------------------------------------     
	      
      Y=Y-Y_0

      dDdot_dgamma=zero
     
      if(Y.gt.zero)then
        dDdot_dgamma=(Y/Big_S)**Small_s
        dDdot_dgamma=dDdot_dgamma/(one-D)**Beta
      endif 
      
c	-----------------------------------------------------------------------------     
c     apply Yahnshan Lou type shear correction      
c	-----------------------------------------------------------------------------     
c     compute shear max normalized with Mises stress        
c     -----------------------------------------------------------------------------
      shear_max_norm=(s_eff_1-s_eff_3)/seq_eff_Mises    
      dDdot_dgamma=shear_max_norm**shear_max_power*dDdot_dgamma
     
      return
      end
      
c	*****************************************************************************
c	*****************************************************************************

	subroutine return_map(
c     INPUTS	 
c	-----------------------------------------------------------------------------
     1                stepTime,i_flow_index_tracker,
	2                stress_hom_tri_v6,
	3                alpha_N,D_N,
     5                X_N_v6,
	4                C_e_eff_m66,
c     OUTPUTS	 
c	-----------------------------------------------------------------------------
     8	 		      stress_hom_NPO_v6,
     9                x_gamma_NPO,
     1                alpha_NPO,D_NPO,
     5                X_NPO_v6,
	6                delta_alpha,
	1                delta_D,
     3	  			  NUMITERATIONS,x_NORM_R,i_delete,
     5                Y_NPO,triaxiality,shear_max_norm)
c	-----------------------------------------------------------------------------
c	last update   : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

c	"TRIAL" HOMOGENIZED STRESS AT STEP N+1
c	-----------------------------------------------------------------------------
      dimension stress_hom_tri_v6(6)
      dimension X_N_v6(6)

c	VARIABLES AT STEP N+1 (OUTPUTS of this SUBROUTINE)
c	-----------------------------------------------------------------------------
      dimension stress_hom_NPO_v6(6)
      dimension X_NPO_v6(6)
      
c	for CUTTING PLANE iterations with K
c	-----------------------------------------------------------------------------
	dimension stress_hom_NPO_K_v6(6), stress_hom_NPO_KPO_v6(6)
	dimension X_NPO_K_v6(6), X_NPO_KPO_v6(6)
	
	dimension r_NPO_K_v6(6)   !d_yld_fnctn_d_stress_hom_NPO_K_v6(6)
	dimension r_NPO_KPO_v6(6) !d_yld_fnctn_d_stress_hom_NPO_KPO_v6(6)
	dimension s_NPO_K_v6(6)   !d_plastic_potential_d_X_NPO_K_v6(6)
	dimension s_NPO_KPO_v6(6) !d_plastic_potential_d_X_NPO_KPO_v6(6)
	
	dimension C_e_eff_m66(6,6)
	dimension C_e_hom_m66(6,6)
	
      COMMON /constants1/  zero,one,two,three,four
	COMMON /constants2/  one_third,half,two_thirds,three_halves
	COMMON /constants3/  sqrt_two_thirds, sqrt_three_halves
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	COMMON /tolerances/  ALF,TOLX,TOLF,TOLMIN,EPSILON,MAXITS
	COMMON /HILL/ x_F,x_G,x_H,x_M,x_N,x_L
	COMMON /LEMAITRE/ Beta, Big_S, Small_s, Y_0, D_cri, D_del
	COMMON /KINHARD/ C_X, Q_X ! Armstrong Frederic type nonlinear kinematic hardening
	COMMON /ISOHARD/ C_q, Q_q 
	COMMON /VISCOSITY/ ETA	
	
c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_YES) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) " VISITING  return map"
	write(6,*) "*****************************************************"	
	endif
c	-----------------------------------------------------------------------------

c	just in case
c	-----------------------------------------------------------------------------
      stress_eq_NPO_K   =zero
      stress_eq_NPO_KPO =zero

	Y_NPO_K  =zero
	Y_NPO_KPO=zero 

c	initialize the values
c	-----------------------------------------------------------------------------
	alpha_NPO_KPO     =alpha_N
	D_NPO_KPO         =D_N
      
      x_gamma_NPO_KPO   =zero     ! as an initial value

	stress_hom_NPO_KPO_v6(1)=stress_hom_tri_v6(1)
	stress_hom_NPO_KPO_v6(2)=stress_hom_tri_v6(2)
	stress_hom_NPO_KPO_v6(3)=stress_hom_tri_v6(3)
	stress_hom_NPO_KPO_v6(4)=stress_hom_tri_v6(4)
	stress_hom_NPO_KPO_v6(5)=stress_hom_tri_v6(5)
	stress_hom_NPO_KPO_v6(6)=stress_hom_tri_v6(6)
c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "KONTROL b"
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif
c	-----------------------------------------------------------------------------	  	
      X_NPO_KPO_v6(1)=X_N_v6(1)
      X_NPO_KPO_v6(2)=X_N_v6(2)
      X_NPO_KPO_v6(3)=X_N_v6(3)
      X_NPO_KPO_v6(4)=X_N_v6(4)
      X_NPO_KPO_v6(5)=X_N_v6(5)
      X_NPO_KPO_v6(6)=X_N_v6(6)
c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "KONTROL c"
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif
c	-----------------------------------------------------------------------------	  
      chapa_NPO_K   = one
      chapa_NPO_KPO = one
      x_kappa_NPO_K   = one
      x_kappa_NPO_KPO = one
      
      dDdot_dgamma_NPO_K   =zero
      dDdot_dgamma_NPO_KPO =zero
      
      x_K_alpha_prime_NPO_K=zero
      x_K_alpha_prime_NPO_KPO=zero
      
	x_NORM_R		=zero
	NUMITERATIONS	=0
c	-----------------------------------------------------------------------------	      
      triaxiality_1=-1000.0d0 ! Elastic step
      triaxiality_2=-1000.0d0 ! Elastic step
c	-----------------------------------------------------------------------------	      
      triaxiality		=-1000.0d0 ! Elastic step
	shear_max_norm  =-1000.0d0 ! Elastic step
c	-----------------------------------------------------------------------------	      
c	 A material defined in user subroutine VUMAT must be defined as purely elastic 
c	(using the initial elastic modulus) at the beginning of the analysis 
c	(stepTime=0). This is an informative message. It does not necessarily 
c	indicate that user subroutine VUMAT is incorrectly defined.
c	-----------------------------------------------------------------------------
	if (stepTime.eq.0) goto 1000 ! elastic initial step ABAQUS SPECIFIC CODE SEGMENT!!!
c	-----------------------------------------------------------------------------
      call compute_yld_fnctn(stress_hom_NPO_KPO_v6,X_NPO_KPO_v6,
     1                       i_flow_index_tracker,
	2                       alpha_NPO_KPO,D_NPO_KPO, !main problem variables
	3                       r_NPO_KPO_v6,
     4                       s_NPO_KPO_v6,
     5                       chapa_NPO_KPO, ! dyld_dq
     6                       x_kappa_NPO_KPO, ! dyld_dD
     7                       yld_fnctn_NPO_KPO,
     8                       stress_eq_NPO_KPO
     9                       x_K_alpha_prime_NPO_KPO)
c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "YIELD FUNCTION:"
	write(6,*) ""
	write(6,*) "yld_fnctn is...:",yld_fnctn_NPO_KPO
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif
c	-----------------------------------------------------------------------------
           
	if (yld_fnctn_NPO_KPO.lt.TOLF) goto 1000 ! NO PLASTIC FLOW, ELASTIC STEP!!!

	x_NORM_R=yld_fnctn_NPO_KPO
	
c	-----------------------------------------------------------------------------
	i_damage_compute_once=0
c	-----------------------------------------------------------------------------
c	START NEWTON ITERATIONS for LOCAL INTEGRATION (PLASTIC-DAMAGE CORRECTION)!!!
c	-----------------------------------------------------------------------------
	do 5000 while(x_NORM_R.gt.TOLF.and.NUMITERATIONS.lt.MAXITS)!=>START OF NR LOOP

	NUMITERATIONS=NUMITERATIONS+1

	if(i_damage_compute_once.eq.0)then
		call compute_damage_fnctn(stress_hom_NPO_KPO_v6,
     1                         D_NPO_KPO,dDdot_dgamma_NPO_KPO,! d(damage rate)/d(plastic multiplier)
     2                         Y_NPO_KPO,triaxiality,
     3                         shear_max_norm) 
		i_damage_compute_once=1 ! change value for not to be revisited
	endif

c	Transfer iteration varibles (KPO to K)
c	-----------------------------------------------------------------------------
	D_NPO_K            =D_NPO_KPO
	dDdot_dgamma_NPO_K =dDdot_dgamma_NPO_KPO 

	alpha_NPO_K    =alpha_NPO_KPO
	x_gamma_NPO_K	 =x_gamma_NPO_KPO
		
	stress_hom_NPO_K_v6(1)=stress_hom_NPO_KPO_v6(1)
	stress_hom_NPO_K_v6(2)=stress_hom_NPO_KPO_v6(2)
	stress_hom_NPO_K_v6(3)=stress_hom_NPO_KPO_v6(3)
	stress_hom_NPO_K_v6(4)=stress_hom_NPO_KPO_v6(4)
	stress_hom_NPO_K_v6(5)=stress_hom_NPO_KPO_v6(5)
	stress_hom_NPO_K_v6(6)=stress_hom_NPO_KPO_v6(6)
	
      X_NPO_K_v6(1)=X_NPO_KPO_v6(1)
      X_NPO_K_v6(2)=X_NPO_KPO_v6(2)
      X_NPO_K_v6(3)=X_NPO_KPO_v6(3)
      X_NPO_K_v6(4)=X_NPO_KPO_v6(4)
      X_NPO_K_v6(5)=X_NPO_KPO_v6(5)
      X_NPO_K_v6(6)=X_NPO_KPO_v6(6)	

	r_NPO_K_v6(1)=r_NPO_KPO_v6(1)
	r_NPO_K_v6(2)=r_NPO_KPO_v6(2)
	r_NPO_K_v6(3)=r_NPO_KPO_v6(3)
	r_NPO_K_v6(4)=r_NPO_KPO_v6(4)
	r_NPO_K_v6(5)=r_NPO_KPO_v6(5)
	r_NPO_K_v6(6)=r_NPO_KPO_v6(6)

	s_NPO_K_v6(1)=s_NPO_KPO_v6(1)
	s_NPO_K_v6(2)=s_NPO_KPO_v6(2)
	s_NPO_K_v6(3)=s_NPO_KPO_v6(3)
	s_NPO_K_v6(4)=s_NPO_KPO_v6(4)
	s_NPO_K_v6(5)=s_NPO_KPO_v6(5)
	s_NPO_K_v6(6)=s_NPO_KPO_v6(6)
	
	x_K_alpha_prime_NPO_K=x_K_alpha_prime_NPO_KPO
	
	Y_NPO_K=Y_NPO_KPO

      chapa_NPO_K       =chapa_NPO_KPO
      x_kappa_NPO_K     =x_kappa_NPO_KPO
	
	yld_fnctn_NPO_K   =yld_fnctn_NPO_KPO
	stress_eq_NPO_K   =stress_eq_NPO_KPO

c	-----------------------------------------------------------------------------
 	if (i_WANT_OUTPUT.eq.i_YES) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "INSIDE THE CHAMBER -> phase I ---------------------->"
	write(6,*) ""
	write(6,*) "alpha_NPO_K is................:", alpha_NPO_K
	write(6,*) "D_NPO_K is....................:", D_NPO_K
	write(6,*) ""
	write(6,*) "x_gamma_NPO_K is..............:", x_gamma_NPO_K
	write(6,*) "chapa_NPO_K is................:", chapa_NPO_K
	write(6,*) "x_kappa_NPO_K is..............:", x_kappa_NPO_K
	write(6,*) "yld_fnctn_NPO_K is............:", yld_fnctn_NPO_K	
	write(6,*) "stress_eq_NPO_K is............:", stress_eq_NPO_K						
	write(6,*) ""	                            
	write(6,*) "stress_hom_NPO_K_v6(1) is.....:", stress_hom_NPO_K_v6(1)
	write(6,*) "stress_hom_NPO_K_v6(2) is.....:", stress_hom_NPO_K_v6(2)
	write(6,*) "stress_hom_NPO_K_v6(3) is.....:", stress_hom_NPO_K_v6(3)
	write(6,*) "stress_hom_NPO_K_v6(4) is.....:", stress_hom_NPO_K_v6(4)
	write(6,*) "stress_hom_NPO_K_v6(5) is.....:", stress_hom_NPO_K_v6(5)
	write(6,*) "stress_hom_NPO_K_v6(6) is.....:", stress_hom_NPO_K_v6(6)				
	write(6,*) ""	
	write(6,*) "r_NPO_K_v6(1) is..............:", r_NPO_K_v6(1)
	write(6,*) "r_NPO_K_v6(2) is..............:", r_NPO_K_v6(2)
	write(6,*) "r_NPO_K_v6(3) is..............:", r_NPO_K_v6(3)
	write(6,*) "r_NPO_K_v6(4) is..............:", r_NPO_K_v6(4)
	write(6,*) "r_NPO_K_v6(5) is..............:", r_NPO_K_v6(5)
	write(6,*) "r_NPO_K_v6(6) is..............:", r_NPO_K_v6(6)				
	write(6,*) ""	
	write(6,*) "s_NPO_K_v6(1) is..............:", s_NPO_K_v6(1)
	write(6,*) "s_NPO_K_v6(2) is..............:", s_NPO_K_v6(2)
	write(6,*) "s_NPO_K_v6(3) is..............:", s_NPO_K_v6(3)
	write(6,*) "s_NPO_K_v6(4) is..............:", s_NPO_K_v6(4)
	write(6,*) "s_NPO_K_v6(5) is..............:", s_NPO_K_v6(5)
	write(6,*) "s_NPO_K_v6(6) is..............:", s_NPO_K_v6(6)				
	write(6,*) ""	
	write(6,*) "*****************************************************"
	endif	
c	-----------------------------------------------------------------------------

      call compute_hom_elastic_stiffness_matrix_6_by_6
     1                                 (C_e_eff_m66,C_e_hom_m66,D_NPO_K)
     	
      x_nom=yld_fnctn_NPO_K
      
      x_denom=zero      
	do i=1,6
      	do j=1,6
            x_denom=x_denom+r_NPO_K_v6(i)*C_e_hom_m66(i,j)*r_NPO_K_v6(j)
        enddo
      enddo

c	-----------------------------------------------------------------------------
 	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "INSIDE THE CHAMBER -> phase II----------------------->"
	write(6,*) ""
	write(6,*) "x_nom is   ",x_nom
	write(6,*) "x_denom is ",x_denom
	write(6,*) ""	
	write(6,*) "*****************************************************"
	endif
	
c	-----------------------------------------------------------------------------
c     SECOND COMPONENT AT THE DENOMINATOR, FROM LINEARIZATION OF THE STRESS
c     ELIMINATED BY DAMAGE RELATED FORM     
c	-----------------------------------------------------------------------------
c     sum=zero
c     do i=1,6
c         sum=sum+r_NPO_K_v6(i)*stress_hom_NPO_K_v6(i)  
c     enddo
c     
c     sum=sum*dDdot_dgamma_NPO_K/(one-D_NPO_K)
c     
c     x_denom=x_denom+sum
c      
c	-----------------------------------------------------------------------------
c 	if (i_WANT_OUTPUT.eq.i_SOME) then
c	write(6,*) "*****************************************************"
c	write(6,*) ""
c 	write(6,*) "INSIDE THE CHAMBER -> phase III----------------------->"
c	write(6,*) ""
c	write(6,*) "x_denom is ",x_denom
c	write(6,*) ""	
c	write(6,*) "*****************************************************"
c	endif
c	-----------------------------------------------------------------------------
      
      sum=zero
      do i=1,6
         sum=sum+r_NPO_K_v6(i)*s_NPO_K_v6(i) ! vector multiplication energetic (stress x strain) Voigt  
      enddo
 
      sum=sum*two_thirds*(one-D_NPO_K)*Q_X
 
      x_denom=x_denom-sum
 
c	-----------------------------------------------------------------------------
 	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "INSIDE THE CHAMBER -> phase IV----------------------->"
	write(6,*) ""
	write(6,*) "x_denom is ",x_denom
	write(6,*) ""	
	write(6,*) "*****************************************************"
	endif
c	-----------------------------------------------------------------------------
  
      x_denom=x_denom-chapa_NPO_K*x_K_alpha_prime_NPO_K
c     x_denom=x_denom-chapa_NPO_K*Q_q*(one-C_q*alpha_NPO_K)

c	-----------------------------------------------------------------------------
 	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "INSIDE THE CHAMBER -> phase V----------------------->"
	write(6,*) ""
	write(6,*) "x_denom is ",x_denom
	write(6,*) ""	
	write(6,*) "*****************************************************"
	endif
c	-----------------------------------------------------------------------------
c     FIFTH AT THE DENOMINATOR, FROM LINEARIZATION OF THE STRESS
c     ELIMINATES SECOND COMPONENT
c     EQUIVALENCE OF EFFECTIVE AND HOMOGENIZED STRESS SPACE DEFINITION OF 
c     THE YIELD FUNCTION
c     check definition of x_kappa_NPO_K (d_yld_function_d_D)
c	-----------------------------------------------------------------------------
c     
c     x_denom=x_denom-x_kappa_NPO_K(i)*dDdot_dgamma_NPO_K
c      
c	-----------------------------------------------------------------------------
c 	if (i_WANT_OUTPUT.eq.i_SOME) then
c	write(6,*) "*****************************************************"
c	write(6,*) ""
c 	write(6,*) "INSIDE THE CHAMBER -> phase III----------------------->"
c	write(6,*) ""
c	write(6,*) "x_denom is ",x_denom
c	write(6,*) ""	
c	write(6,*) "*****************************************************"
c	endif
c	-----------------------------------------------------------------------------      

 	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "INSIDE THE CHAMBER -> phase V----------------------->"
	write(6,*) ""
	write(6,*) "x_denom is ",x_denom
	write(6,*) ""	
	write(6,*) "*****************************************************"
	endif
c	-----------------------------------------------------------------------------
      
      delta_x_gamma_NPO_KPO=x_nom/x_denom
      x_gamma_NPO_KPO=x_gamma_NPO_K+delta_x_gamma_NPO_KPO

c	-----------------------------------------------------------------------------
 	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "INSIDE THE CHAMBER -> phase VII--------------------->"
	write(6,*) ""
	write(6,*) "delta_x_gamma_NPO_KPO is ",delta_x_gamma_NPO_KPO
	write(6,*) ""	
	write(6,*) "*****************************************************"
	endif
c	-----------------------------------------------------------------------------
	
	stress_hom_NPO_KPO_v6(1)=stress_hom_NPO_K_v6(1)
	stress_hom_NPO_KPO_v6(2)=stress_hom_NPO_K_v6(2)
	stress_hom_NPO_KPO_v6(3)=stress_hom_NPO_K_v6(3)
	stress_hom_NPO_KPO_v6(4)=stress_hom_NPO_K_v6(4)
	stress_hom_NPO_KPO_v6(5)=stress_hom_NPO_K_v6(5)
	stress_hom_NPO_KPO_v6(6)=stress_hom_NPO_K_v6(6)
	
      do i=1,6
        do j=1,6
            stress_hom_NPO_KPO_v6(i)=stress_hom_NPO_KPO_v6(i)
     1                                -delta_x_gamma_NPO_KPO*
     2                         C_e_hom_m66(i,j)*r_NPO_K_v6(j)
        enddo
      enddo

	stress_hom_NPO_KPO_v6(1)=stress_hom_NPO_KPO_v6(1)
	1 -delta_x_gamma_NPO_KPO*dDdot_dgamma_NPO_K
	2 *stress_hom_NPO_K_v6(1)/(one-D_NPO_K)
	stress_hom_NPO_KPO_v6(2)=stress_hom_NPO_KPO_v6(2)
	1 -delta_x_gamma_NPO_KPO*dDdot_dgamma_NPO_K
	2 *stress_hom_NPO_K_v6(2)/(one-D_NPO_K)
	stress_hom_NPO_KPO_v6(3)=stress_hom_NPO_KPO_v6(3)
	1 -delta_x_gamma_NPO_KPO*dDdot_dgamma_NPO_K
	2 *stress_hom_NPO_K_v6(3)/(one-D_NPO_K)
	stress_hom_NPO_KPO_v6(4)=stress_hom_NPO_KPO_v6(4)
	1 -delta_x_gamma_NPO_KPO*dDdot_dgamma_NPO_K
	2 *stress_hom_NPO_K_v6(4)/(one-D_NPO_K)
	stress_hom_NPO_KPO_v6(5)=stress_hom_NPO_KPO_v6(5)
	1 -delta_x_gamma_NPO_KPO*dDdot_dgamma_NPO_K
	2 *stress_hom_NPO_K_v6(5)/(one-D_NPO_K)
	stress_hom_NPO_KPO_v6(6)=stress_hom_NPO_KPO_v6(6)
	1 -delta_x_gamma_NPO_KPO*dDdot_dgamma_NPO_K
	2 *stress_hom_NPO_K_v6(6)/(one-D_NPO_K)

c	-----------------------------------------------------------------------------
 	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "INSIDE THE CHAMBER -> phase IV---------------------->"
	write(6,*) ""
	write(6,*) "stress_hom_NPO_KPO_v6(1) is ",stress_hom_NPO_KPO_v6(1)
	write(6,*) "stress_hom_NPO_KPO_v6(2) is ",stress_hom_NPO_KPO_v6(2)
	write(6,*) "stress_hom_NPO_KPO_v6(3) is ",stress_hom_NPO_KPO_v6(3)
	write(6,*) "stress_hom_NPO_KPO_v6(4) is ",stress_hom_NPO_KPO_v6(4)
	write(6,*) "stress_hom_NPO_KPO_v6(5) is ",stress_hom_NPO_KPO_v6(5)
	write(6,*) "stress_hom_NPO_KPO_v6(6) is ",stress_hom_NPO_KPO_v6(6)				
	write(6,*) ""	
	write(6,*) "*****************************************************"
	endif
c	-----------------------------------------------------------------------------

      X_NPO_KPO_v6(1)=X_NPO_K_v6(1)
     1               -two_thirds*Q_X*delta_x_gamma_NPO_KPO*s_NPO_K_v6(1)
      X_NPO_KPO_v6(2)=X_NPO_K_v6(2)
     1               -two_thirds*Q_X*delta_x_gamma_NPO_KPO*s_NPO_K_v6(2)
      X_NPO_KPO_v6(3)=X_NPO_K_v6(3)
     1               -two_thirds*Q_X*delta_x_gamma_NPO_KPO*s_NPO_K_v6(3)
      X_NPO_KPO_v6(4)=X_NPO_K_v6(4)
     1               -two_thirds*Q_X*delta_x_gamma_NPO_KPO*s_NPO_K_v6(4) ! continuum vs voigt form
      X_NPO_KPO_v6(5)=X_NPO_K_v6(5)
     1               -two_thirds*Q_X*delta_x_gamma_NPO_KPO*s_NPO_K_v6(5) ! continuum vs voigt form
      X_NPO_KPO_v6(6)=X_NPO_K_v6(6)
     1               -two_thirds*Q_X*delta_x_gamma_NPO_KPO*s_NPO_K_v6(6) ! continuum vs voigt form
      

      alpha_NPO_KPO=alpha_NPO_K+delta_x_gamma_NPO_KPO
     
    
	D_NPO_KPO    =D_NPO_K+delta_x_gamma_NPO_KPO*dDdot_dgamma_NPO_K
	
c	-----------------------------------------------------------------------------
 	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "INSIDE THE CHAMBER -> phase VII--------------------->"
	write(6,*) ""
	write(6,*) "alpha_NPO_KPO is.................",alpha_NPO_KPO
	write(6,*) "D_NPO_KPO is.....................",D_NPO_KPO
	write(6,*) ""	
	write(6,*) "*****************************************************"
	endif
c	-----------------------------------------------------------------------------
	      
      call compute_yld_fnctn(stress_hom_NPO_KPO_v6,X_NPO_KPO_v6,
     1                       i_flow_index_tracker,
	2                       alpha_NPO_KPO,D_NPO_KPO, !main problem variables
	3                       r_NPO_KPO_v6,
     4                       s_NPO_KPO_v6,
     5                       chapa_NPO_KPO, ! dyld_dq
     6                       x_kappa_NPO_KPO, ! dyld_dD
     7                       yld_fnctn_NPO_KPO,
     8                       stress_eq_NPO_KPO
     9                       x_K_alpha_prime_NPO_KPO)

      call compute_damage_fnctn(stress_hom_NPO_KPO_v6,
     1                         D_NPO_KPO,dDdot_dgamma_NPO_KPO,! d(damage rate)/d(plastic multiplier)
     2                         Y_NPO_KPO,triaxiality,
     3                         shear_max_norm)   
     
	call precision_abs(yld_fnctn_NPO_KPO,x_NORM_R)
	
c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
 	write(6,*) "END OF ITERATION..........:",NUMITERATIONS
	write(6,*) ""
	write(6,*) "yld_fnctn is..............:",yld_fnctn_NPO_KPO
	write(6,*) ""	
	write(6,*) "stress_hom_NPO_KPO_v6(1) is ",stress_hom_NPO_KPO_v6(1)
	write(6,*) "stress_hom_NPO_KPO_v6(2) is ",stress_hom_NPO_KPO_v6(2)
	write(6,*) "stress_hom_NPO_KPO_v6(3) is ",stress_hom_NPO_KPO_v6(3)
	write(6,*) "stress_hom_NPO_KPO_v6(4) is ",stress_hom_NPO_KPO_v6(4)
	write(6,*) "stress_hom_NPO_KPO_v6(5) is ",stress_hom_NPO_KPO_v6(5)
	write(6,*) "stress_hom_NPO_KPO_v6(6) is ",stress_hom_NPO_KPO_v6(6)
	write(6,*) ""
	write(6,*) "alpha_NPO_KPO            is ",alpha_NPO_KPO
	write(6,*) "D_NPO_KPO                is ",D_NPO_KPO
	write(6,*) ""
	write(6,*) "*****************************************************"
	endif
c	-----------------------------------------------------------------------------
	

5000	continue !=============================================> END OF NR LOOP
	
1000	alpha_NPO		          =alpha_NPO_KPO
	x_gamma_NPO	          =x_gamma_NPO_KPO     ! as an initial value
      D_NPO                   =D_NPO_KPO
      
c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "KONTROL 3"
	write(6,*) ""
	write(6,*) "**********************************************************"
	write(6,*) "**********************************************************"
	write(6,*) ""
	endif
c	-----------------------------------------------------------------------------	     

    	stress_hom_NPO_v6(1)=stress_hom_NPO_KPO_v6(1)
	stress_hom_NPO_v6(2)=stress_hom_NPO_KPO_v6(2)
	stress_hom_NPO_v6(3)=stress_hom_NPO_KPO_v6(3)
	stress_hom_NPO_v6(4)=stress_hom_NPO_KPO_v6(4)
	stress_hom_NPO_v6(5)=stress_hom_NPO_KPO_v6(5)
	stress_hom_NPO_v6(6)=stress_hom_NPO_KPO_v6(6)
	
	X_NPO_v6(1)=X_NPO_KPO_v6(1)
	X_NPO_v6(2)=X_NPO_KPO_v6(2)
	X_NPO_v6(3)=X_NPO_KPO_v6(3)
	X_NPO_v6(4)=X_NPO_KPO_v6(4)
	X_NPO_v6(5)=X_NPO_KPO_v6(5)
	X_NPO_v6(6)=X_NPO_KPO_v6(6)

      delta_alpha=alpha_NPO-alpha_N
      delta_D=D_NPO-D_N
          	
	if (i_WANT_OUTPUT.eq.i_SOME) then
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "END OF ITERATIONS:"
	write(6,*) ""
	write(6,*) "stress_hom_NPO_v6(1) is ",stress_hom_NPO_v6(1)
	write(6,*) "stress_hom_NPO_v6(2) is ",stress_hom_NPO_v6(2)
	write(6,*) "stress_hom_NPO_v6(3) is ",stress_hom_NPO_v6(3)
	write(6,*) "stress_hom_NPO_v6(4) is ",stress_hom_NPO_v6(4)
	write(6,*) "stress_hom_NPO_v6(5) is ",stress_hom_NPO_v6(5)
	write(6,*) "stress_hom_NPO_v6(6) is ",stress_hom_NPO_v6(6)
	write(6,*) ""
	write(6,*) "X_NPO_v6(1) is          ",X_NPO_v6(1)
	write(6,*) "X_NPO_v6(2) is          ",X_NPO_v6(2)
	write(6,*) "X_NPO_v6(3) is          ",X_NPO_v6(3)
	write(6,*) "X_NPO_v6(4) is          ",X_NPO_v6(4)
	write(6,*) "X_NPO_v6(5) is          ",X_NPO_v6(5)
	write(6,*) "X_NPO_v6(6) is          ",X_NPO_v6(6)
	write(6,*) ""
	write(6,*) "alpha_NPO        is     ",alpha_NPO
	write(6,*) "D_NPO            is     ",D_NPO
	write(6,*) ""
	write(6,*) "*****************************************************"
	endif

	return 
	end
            
c	*****************************************************************************
c	*****************************************************************************

	subroutine precision_abs(x_input,x_output)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

c	computational clarity
	if(j_sys_Dimension.eq.1)then
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=abs(x_input)
	else
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=abs(x_input)
	endif
	
	return 
	end
	       
c	*****************************************************************************
c	*****************************************************************************

	subroutine precision_sinh(x_input,x_output)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

c	computational clarity
	if(j_sys_Dimension.eq.1)then
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=sinh(x_input)
	else
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=sinh(x_input)
	endif
	
	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine precision_cosh(x_input,x_output)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

c	computational clarity
	if(j_sys_Dimension.eq.1)then
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=cosh(x_input)
	else
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=cosh(x_input)
	endif
	
	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine precision_sqrt(x_input,x_output)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

c	computational clarity
	if(j_sys_Dimension.eq.1)then
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=sqrt(x_input)
	else
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=sqrt(x_input)
	endif
	
	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine precision_exp(x_input,x_output)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

c	computational clarity
	if(j_sys_Dimension.eq.1)then
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=exp(x_input)
	else
c		write(6,*) "PROBLEM AREA j_sys_Dimension", j_sys_Dimension
		x_output=exp(x_input)
	endif
	
	return 
	end
	            
c	*****************************************************************************
c	*****************************************************************************
     
      subroutine compute_yld_fnctn(stress_hom_v6,X_v6,
     1                             i_flow_index_tracker,
	2                             alpha,D, !main problem variables
	3                             r_v6,
     4                             s_v6,
     5                             chapa,     ! dyld_dq
     6                             x_kappa,   ! dyld_dD
     7                             yld_fnctn,
     8                             stress_eq
     9                             x_K_alpha_prime)
c	-----------------------------------------------------------------------------
c	last update   : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'
	
      dimension stress_hom_v6(6)
      dimension X_v6(6)
      dimension r_v6(6)
      dimension s_v6(6)

c     Internally used 
c	-----------------------------------------------------------------------------
      dimension stress_eff_minus_X_v6(6)
	dimension H_Hill_m66(6,6)	

	COMMON /ISO_HARD_FLAGS   /  i_ISO_HARD_MODEL,i_POWER,i_SATURATION
	COMMON /ISO_HARD_FLAGS2  /  i_POWER_PIRES,i_POWER_LUDWIK
	COMMON /ISO_HARD_FLAGS3  /  i_DOUBLE_SATURATION,i_SWIFT,i_DISCRETE	

	COMMON /constants1/  zero,one,two,three,four
	COMMON /constants2/  one_third,half,two_thirds,three_halves	
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /tolerances/  ALF,TOLX,TOLF,TOLMIN,EPSILON,MAXITS	
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	COMMON /HILL/ x_F,x_G,x_H,x_M,x_N,x_L
	COMMON /LEMAITRE/ Beta, Big_S, Small_s, Y_0, D_cri, D_del
	COMMON /KINHARD/ C_X, Q_X ! Armstrong Frederic type nonlinear kinematic h.
	COMMON /ISOHARD/ C_q, Q_q 
	COMMON /VISCOSITY/ ETA		
      
c     -----------------------------------------------------------------------------  
	if (i_WANT_OUTPUT.eq.i_YES) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) " VISITING  compute_yld_fnctn"
	write(6,*) "*****************************************************"	
	endif	
c     -----------------------------------------------------------------------------  

c     -----------------------------------------------------------------------------  
c     FOLLOWING DOIG AND FRENCH SCHOOL
c     -----------------------------------------------------------------------------  
c      x_K_alpha      =Q_q*alpha ! q
c      x_K_alpha_prime=Q_q       ! d_q/d_alpha

c     PLASTIC ISOTROPIC HARDENING
c     -----------------------------------------------------------------------------  
	if(i_ISO_HARD_MODEL.eq.i_SATURATION)then
	    call compute_K_alpha_saturation
	1								 (alpha,x_K_alpha,x_K_alpha_prime)
c     -----------------------------------------------------------------------------  
	else if(i_ISO_HARD_MODEL.eq.i_DOUBLE_SATURATION)then
		call compute_K_alpha_double_saturation
	1	                             (alpha,x_K_alpha,x_K_alpha_prime)
c     -----------------------------------------------------------------------------  
	else if(i_ISO_HARD_MODEL.eq.i_POWER)then
		call compute_K_alpha_power
	1								 (alpha,x_K_alpha,x_K_alpha_prime)
c     -----------------------------------------------------------------------------  
	else if(i_ISO_HARD_MODEL.eq.i_POWER_PIRES)then
		call compute_K_alpha_power_pires
	1								 (alpha,x_K_alpha,x_K_alpha_prime)
c     -----------------------------------------------------------------------------  
	else if(i_ISO_HARD_MODEL.eq.i_POWER_LUDWIK)then
		call compute_K_alpha_power_ludwik
	1	                             (alpha,x_K_alpha,x_K_alpha_prime)
c     -----------------------------------------------------------------------------  
	else if(i_ISO_HARD_MODEL.eq.i_SWIFT)then
		call compute_K_alpha_swift
	1	                             (alpha,x_K_alpha,x_K_alpha_prime)
c     -----------------------------------------------------------------------------  
	else if(i_ISO_HARD_MODEL.eq.i_DISCRETE)then
		call compute_K_alpha_discrete(i_flow_index_tracker,alpha,
     1                             x_K_alpha,x_K_alpha_prime) 
c     -----------------------------------------------------------------------------  
	endif
c     -----------------------------------------------------------------------------  

c	compute yield function
c	-----------------------------------------------------------------------------
      x_coef=one/(one-D)
      stress_eff_minus_X_v6(1)=x_coef*stress_hom_v6(1)-X_v6(1)
      stress_eff_minus_X_v6(2)=x_coef*stress_hom_v6(2)-X_v6(2)
      stress_eff_minus_X_v6(3)=x_coef*stress_hom_v6(3)-X_v6(3)
      stress_eff_minus_X_v6(4)=x_coef*stress_hom_v6(4)-X_v6(4)
      stress_eff_minus_X_v6(5)=x_coef*stress_hom_v6(5)-X_v6(5)
      stress_eff_minus_X_v6(6)=x_coef*stress_hom_v6(6)-X_v6(6)

	call compute_Hill48_structural_matrix_6_by_6(H_Hill_m66)
	      
      sum=zero      
	do i=1,6
      	do j=1,6
            sum=sum+stress_eff_minus_X_v6(i)
     1                *H_Hill_m66(i,j)*stress_eff_minus_X_v6(j)
        enddo
      enddo      

	call precision_sqrt(sum,stress_eq)
	if(stress_eq.lt.EPSILON)stress_eq=EPSILON ! to avoid zero division		
	
	yld_stress=sigma_y_0_11+x_K_alpha
		
	yld_fnctn=stress_eq-yld_stress
	
	sum=zero      
      
c	initialize r and s vectors
c	-----------------------------------------------------------------------------
      do i=1,6
      	   r_v6(i)=zero
      	   s_v6(i)=zero
      enddo   

      x_coef=one/(one-D)/stress_eq
	do i=1,6
      	do j=1,6
           r_v6(i)=r_v6(i)+H_Hill_m66(i,j)*stress_eff_minus_X_v6(j)
      	enddo
      	r_v6(i)=x_coef*r_v6(i)
      enddo   
	
	x_coef=three_halves*C_X/Q_X
	
      s_v6(1)=     -(one-D)*r_v6(1)+x_coef*X_v6(1)
      s_v6(2)=     -(one-D)*r_v6(2)+x_coef*X_v6(2)
      s_v6(3)=     -(one-D)*r_v6(3)+x_coef*X_v6(3)
      s_v6(4)=-half*(one-D)*r_v6(4)+x_coef*X_v6(4)
      s_v6(5)=-half*(one-D)*r_v6(5)+x_coef*X_v6(5)
      s_v6(6)=-half*(one-D)*r_v6(6)+x_coef*X_v6(6)
      
      chapa=-one
     
c	needs actual computation, not done yet since eliminated in further stages
c	-----------------------------------------------------------------------------
      x_kappa=one
     		     	
	return
	end

c	*****************************************************************************
c	*****************************************************************************
c
c	AUX MATHEMATICAL ROUTINES AND INITIALIZATIONS
c
c	*****************************************************************************
c	*****************************************************************************

	subroutine compute_eff_elastic_stiffness_matrix_6_by_6(C_e)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

      dimension C_e(6,6)
	
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /constants1/  zero,one,two,three,four
	COMMON /constants2/  one_third,half,two_thirds,three_halves	
	
	x_diag   =E*(one-poisson)/(one+poisson)/(one-two*poisson)
	x_nondiag=E/(one+poisson)/(one-two*poisson)*poisson
	x_shear  =half*E/(one+poisson)

	C_e(1,1)=x_diag
	C_e(1,2)=x_nondiag
	C_e(1,3)=x_nondiag
	C_e(1,4)=zero
	C_e(1,5)=zero
	C_e(1,6)=zero

	C_e(2,1)=x_nondiag
	C_e(2,2)=x_diag
	C_e(2,3)=x_nondiag
	C_e(2,4)=zero
	C_e(2,5)=zero
	C_e(2,6)=zero

	C_e(3,1)=x_nondiag
	C_e(3,2)=x_nondiag
     	C_e(3,3)=x_diag
	C_e(3,4)=zero
	C_e(3,5)=zero
	C_e(3,6)=zero
      
	C_e(4,1)=zero
	C_e(4,2)=zero
     	C_e(4,3)=zero
	C_e(4,4)=x_shear
	C_e(4,5)=zero
	C_e(4,6)=zero

	C_e(5,1)=zero
	C_e(5,2)=zero
     	C_e(5,3)=zero
	C_e(5,4)=zero
	C_e(5,5)=x_shear
	C_e(5,6)=zero
      
	C_e(6,1)=zero
	C_e(6,2)=zero
     	C_e(6,3)=zero
	C_e(6,4)=zero
	C_e(6,5)=zero
	C_e(6,6)=x_shear

	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine compute_eff_elastic_compliance_matrix_6_by_6(C_e_inv)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------

      include 'vaba_param.inc'

      dimension C_e_inv(6,6)
	
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /constants1/  zero,one,two,three,four
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME

	diag=one/E
	offdiag=-poisson/E
	diagshear=one/x_mu

	C_e_inv(1,1)=diag
	C_e_inv(1,2)=offdiag
	C_e_inv(1,3)=offdiag
	C_e_inv(1,4)=zero
	C_e_inv(1,5)=zero
	C_e_inv(1,6)=zero

	C_e_inv(2,1)=offdiag
	C_e_inv(2,2)=diag
	C_e_inv(2,3)=offdiag
	C_e_inv(2,4)=zero
	C_e_inv(2,5)=zero
	C_e_inv(2,6)=zero

	C_e_inv(3,1)=offdiag
	C_e_inv(3,2)=offdiag
	C_e_inv(3,3)=diag
	C_e_inv(3,4)=zero
	C_e_inv(3,5)=zero
	C_e_inv(3,6)=zero

	C_e_inv(4,1)=zero
	C_e_inv(4,2)=zero
	C_e_inv(4,3)=zero
	C_e_inv(4,4)=diagshear
	C_e_inv(4,5)=zero
	C_e_inv(4,6)=zero

	C_e_inv(5,1)=zero
	C_e_inv(5,2)=zero
	C_e_inv(5,3)=zero
	C_e_inv(5,4)=zero
	C_e_inv(5,5)=diagshear
	C_e_inv(5,6)=zero

	C_e_inv(6,1)=zero
	C_e_inv(6,2)=zero
	C_e_inv(6,3)=zero
	C_e_inv(6,4)=zero
	C_e_inv(6,5)=zero
	C_e_inv(6,6)=diagshear

	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine compute_hom_elastic_stiffness_matrix_6_by_6(C_eff,C_hom,D)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

      dimension C_eff(6,6),C_hom(6,6)

	COMMON /constants1/  zero,one,two,three,four
	
	x_coef=one-D
	
	C_hom(1,1)=x_coef*C_eff(1,1)
	C_hom(1,2)=x_coef*C_eff(1,2)
	C_hom(1,3)=x_coef*C_eff(1,3)
	C_hom(1,4)=x_coef*C_eff(1,4)
	C_hom(1,5)=x_coef*C_eff(1,5)
	C_hom(1,6)=x_coef*C_eff(1,6)
               
	C_hom(2,1)=x_coef*C_eff(2,1)
	C_hom(2,2)=x_coef*C_eff(2,2)
	C_hom(2,3)=x_coef*C_eff(2,3)
	C_hom(2,4)=x_coef*C_eff(2,4)
	C_hom(2,5)=x_coef*C_eff(2,5)
	C_hom(2,6)=x_coef*C_eff(2,6)
               
	C_hom(3,1)=x_coef*C_eff(3,1)
	C_hom(3,2)=x_coef*C_eff(3,2)
     	C_hom(3,3)=x_coef*C_eff(3,3)
	C_hom(3,4)=x_coef*C_eff(3,4)
	C_hom(3,5)=x_coef*C_eff(3,5)
	C_hom(3,6)=x_coef*C_eff(3,6)
               
	C_hom(4,1)=x_coef*C_eff(4,1)
	C_hom(4,2)=x_coef*C_eff(4,2)
     	C_hom(4,3)=x_coef*C_eff(4,3)
	C_hom(4,4)=x_coef*C_eff(4,4)
	C_hom(4,5)=x_coef*C_eff(4,5)
	C_hom(4,6)=x_coef*C_eff(4,6)
                
	C_hom(5,1)=x_coef*C_eff(5,1)
	C_hom(5,2)=x_coef*C_eff(5,2)
     	C_hom(5,3)=x_coef*C_eff(5,3)
	C_hom(5,4)=x_coef*C_eff(5,4)
	C_hom(5,5)=x_coef*C_eff(5,5)
	C_hom(5,6)=x_coef*C_eff(5,6)
                
	C_hom(6,1)=x_coef*C_eff(6,1)
	C_hom(6,2)=x_coef*C_eff(6,2)
     	C_hom(6,3)=x_coef*C_eff(6,3)
	C_hom(6,4)=x_coef*C_eff(6,4)
	C_hom(6,5)=x_coef*C_eff(6,5)
	C_hom(6,6)=x_coef*C_eff(6,6)

	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine compute_Hill48_structural_matrix_6_by_6(H_Hill)
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 09APRIL2009
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

      dimension H_Hill(6,6)
	
	COMMON /constants1/  zero,one,two,three,four
	COMMON /constants2/  one_third,half,two_thirds,three_halves	
	COMMON /HILL/ x_F,x_G,x_H,x_M,x_N,x_L	
	
	H_Hill(1,1)= x_G+x_H
	H_Hill(1,2)=-x_H
	H_Hill(1,3)=-x_G
	H_Hill(1,4)=zero
	H_Hill(1,5)=zero
	H_Hill(1,6)=zero

	H_Hill(2,1)=-x_H
	H_Hill(2,2)= x_H+x_F
	H_Hill(2,3)=-x_F
	H_Hill(2,4)=zero
	H_Hill(2,5)=zero
	H_Hill(2,6)=zero

	H_Hill(3,1)=-x_G
	H_Hill(3,2)=-x_F
     	H_Hill(3,3)= x_F+x_G
	H_Hill(3,4)=zero
	H_Hill(3,5)=zero
	H_Hill(3,6)=zero
      
	H_Hill(4,1)=zero
	H_Hill(4,2)=zero
     	H_Hill(4,3)=zero
	H_Hill(4,4)=two*x_N
	H_Hill(4,5)=zero
	H_Hill(4,6)=zero

	H_Hill(5,1)=zero
	H_Hill(5,2)=zero
     	H_Hill(5,3)=zero
	H_Hill(5,4)=zero
	H_Hill(5,5)=two*x_M
	H_Hill(5,6)=zero
      
	H_Hill(6,1)=zero
	H_Hill(6,2)=zero
     	H_Hill(6,3)=zero
	H_Hill(6,4)=zero
	H_Hill(6,5)=zero
	H_Hill(6,6)=two*x_L
	
	return 
	end
		
c	*****************************************************************************
c	*****************************************************************************

      subroutine jacobi(b,stret,princ)
c-----------------------------------------------------------------------
c
c     Evaluates the principal values and principal directions given any symmetric
c	second order b tensor using the Jacobi iteration. Adapted from numerical recipies
c
c     b(6)         -->  any second order tensor in vector form, VUMAT order 11,22,33,12,23,13
c     stret(3)     -->  vector containing the eigenvalues of b
c     princ(3,3)   -->  matrix containing the three principal row vectors
c
c-----------------------------------------------------------------------
      include 'vaba_param.inc'

      dimension b(6),btens(3,3),stret(3),princ(3,3)

	COMMON /tolerances/  ALF,TOLX,TOLF,TOLMIN,EPSILON,MAXITS
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME

	btens(1,1)=b(1)
	btens(2,2)=b(2)
	btens(3,3)=b(3)
	btens(1,2)=b(4)
	btens(2,1)=b(4)
	btens(3,2)=b(5)
	btens(2,3)=b(5)
	btens(1,3)=b(6)
	btens(3,1)=b(6)
c
c     Initialise princ to the identity
c     
      do i=1,3
         do j=1,3
            princ(i,j)=0.0d0
         enddo
         princ(i,i)=1.0d0
         stret(i)=btens(i,i)
      enddo
c
c     Starts sweeping. 
c
      do is=1,50
         sum=0.0d0
         do ip=1,2
            do iq=ip+1,3
               sum=sum+abs(btens(ip,iq))
            enddo
         enddo
c
c     If the sum of off-diagonal terms is zero evaluates the 
c     stretches and returns
c
         if(sum.lt.EPSILON) return
c
c     Performs the sweep in three rotations. One per off diagonal term     
c
         do ip=1,2
            do iq=ip+1,3
               od=100.*abs(btens(ip,iq))
               if((od+abs(stret(ip)).ne.abs(stret(ip))).and.
     &            (od+abs(stret(iq)).ne.abs(stret(iq)))) then
                  hd=stret(iq)-stret(ip)
c
c    Evaluates the rotation angle 
c
                  if(abs(hd)+od.eq.abs(hd)) then
                    t=btens(ip,iq)/hd
                  else
                    theta=0.5d0*hd/btens(ip,iq)
                     t=1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
                     if(theta.lt.0.) t=-t
                  endif
c
c     Re-evaluates the diagonal terms
c
                  c=1.0d0/sqrt(1.0d0+t**2)
                  s=t*c
                  tau=s/(1.0d0+c)
                  h=t*btens(ip,iq)
                  stret(ip)=stret(ip)-h
                  stret(iq)=stret(iq)+h
c
c     Re-evaluates the remaining off-diagonal terms         
c
                  ir=6-ip-iq
                  g=btens(min(ir,ip),max(ir,ip))
                  h=btens(min(ir,iq),max(ir,iq))
                  btens(min(ir,ip),max(ir,ip))=g-s*(h+g*tau)
                  btens(min(ir,iq),max(ir,iq))=h+s*(g-h*tau)
c
c     Rotates the eigenvectors
c
                  do ir=1,3
                     g=princ(ip,ir)
                     h=princ(iq,ir)
                     princ(ip,ir)=g-s*(h+g*tau)
                     princ(iq,ir)=h+s*(g-h*tau)
                  enddo
               endif
               btens(ip,iq)=0.0d0
            enddo
         enddo
      enddo
c
c     If convergence is not achieved stops
c
	return
      end
      		
c	*****************************************************************************
c	*****************************************************************************

	subroutine initialize_constants()
c	-----------------------------------------------------------------------------
c	last update : 25 APRIL 2009, 01:13, Celal Soyarslan
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

	COMMON /constants1/ zero,one,two,three,four
	COMMON /constants2/ one_third,half,two_thirds,three_halves
	COMMON /constants3/ sqrt_two_thirds, sqrt_three_halves
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME

c	computational clarity

      zero				= 0.0D0
	one				= 1.0D0
	two				= 2.0D0
	three				= 3.0D0
	four				= 4.0D0
	half				= 0.5D0
	three_halves		= 1.5D0
	one_third			= one/three
	two_thirds		= two/three
	call precision_sqrt(two/three,sqrt_two_thirds)
	call precision_sqrt(three/two,sqrt_three_halves)

	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine initialize_computational_parameters()
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

	COMMON /tolerances/  ALF,TOLX,TOLF,TOLMIN,EPSILON,MAXITS
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME

c	tolerances	
c	-----------------------------------------------------------------------------
	TOLX	    =1.e-10
	TOLF	    =1.e-10
	ALF		=1.e-10
	TOLMIN    =1.e-10
	EPSILON   =1.e-10  !!! controls the yield band near zero !!!
	MAXITS    =1000

	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine initialize_analysis_flags()
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

	COMMON /ISO_HARD_FLAGS   /  i_ISO_HARD_MODEL,i_POWER,i_SATURATION
	COMMON /ISO_HARD_FLAGS2  /  i_POWER_PIRES,i_POWER_LUDWIK
	COMMON /ISO_HARD_FLAGS3  /  i_DOUBLE_SATURATION,i_SWIFT,i_DISCRETE
	
	COMMON /OUTPUT_FLAGS     /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
	COMMON /SOLN_SCHEME_FLAGS/  i_SOLN_SCHEME,i_SIMUL,i_STAG

c	NONLINEAR SOLUTION SCHEME ENUMERATIONS
c	-----------------------------------------------------------------------------
	i_SIMUL			=1
	i_STAG       		=2

c	ISOTROPIC HARDENING MODEL ENUMERATIONS
c	-----------------------------------------------------------------------------
	i_POWER			    =1
	i_SATURATION		    =2
	i_POWER_PIRES		    =3
	i_POWER_LUDWIK        =4
	i_DOUBLE_SATURATION   =5
	i_SWIFT				=6
	i_DISCRETE			=7

c	ADDITIONAL ENUMERATIONS 
c	-----------------------------------------------------------------------------
	i_YES				=1
	i_NO				=2
	i_SOME			=3

	return
	end

c	*****************************************************************************
c	*****************************************************************************	

	subroutine compute_elastic_constants(props)
c	-----------------------------------------------------------------------------
c	Input:
c	props(3)			is the properties array
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

	dimension props(3)

	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
	COMMON /constants1/ zero,one,two,three,four
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME

c	TEST PARAMETERS
c	-----------------------------------------------------------------------------
	sigma_y_0_11=sigma_y_0_11_par
	poisson=0.30d0
	E=210000.0d0
	x_mu=E/two/(one+poisson)
	x_hh=E/three/(one-two*poisson)
	two_mu=two*x_mu
	
	return
	end

c	*****************************************************************************
c	*****************************************************************************

		subroutine compute_inelastic_constants(props)
c	-----------------------------------------------------------------------------
c	Input:
c	props(3)			is the properties array
c	-----------------------------------------------------------------------------
c	Required Output:
c	     				NONE
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS, 06JUNE07
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /elasticity/  E,x_hh,x_mu,two_mu,poisson,sigma_y_0_11
c	------+--------+-------------+----------------+-----------------------+-----+
c	PLASTIC ISOTROPIC HARDENING
c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /plasticity_saturation/  x_K,x_K_inf,x_K_0,x_delta
	COMMON /plasticity_saturation_double/  x_K_inf1,x_K_inf2,
	1                                       x_delta1,x_delta2
	COMMON /plasticity_swift/ x_K_a,x_K_b,x_K_n
	COMMON /Lueders/ sigma_y_1_11,x_alpha_offset,x_K_offset
	COMMON /plasticity_power     /  x_C,x_pow,x_add
	COMMON /FLOW_CURVE/flow_curve(66),i_dim_flow_curve	
c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /KINHARD/ C_X, Q_X ! Armstrong Frederic type nonlinear kinematic hardening
c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /constants1/ zero,one,two,three,four
	COMMON /constants2/ one_third,half,two_thirds,three_halves
c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /HILL/ x_F,x_G,x_H,x_M,x_N,x_L
c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /LEMAITRE/ Beta, Big_S, Small_s, Y_0, D_cri, D_del
	COMMON /CRACK_CLOSURE/ h_crack_closure
	COMMON /MAX_SHEAR_POWER/ shear_max_power	
c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /VISCOSITY/ ETA
	COMMON /OUTPUT_FLAGS    /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
c	------+--------+-------------+----------------+-----------------------+-----+

	dimension props(8)

	i_dim_flow_curve=33

c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+
c	PLASTIC ANISOTROPY PARAMETERS
c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+

c     Lankford's Coefficients: R0, R45, R90: FOR ISOTROPY ALL ARE IDENTITY
c	-----------------------------------------------------------------------------
      R0 =one
      R45=one
      R90=one  

c	COMMON /HILL/ x_F,x_G,x_H,x_M,x_N,x_L
c	-----------------------------------------------------------------------------
c     Conversion from Lankford to Hill follows:
c	-----------------------------------------------------------------------------
      x_F=R0/(R90*(one+R0))
      x_G=one/(one+R0)
      x_H=R0/(one+R0)
      x_M=half*(R0+R90)*(one+two*R45)/(R90*(one+R0))
      x_N=x_M
      x_L=x_M
      
c	-----------------------------------------------------------------------------
	if (i_WANT_OUTPUT.eq.i_YES) then
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	write(6,*) "Inside compute inelastic:"
	write(6,*) ""
	write(6,*) "x_F is...:",x_F
	write(6,*) "x_G is...:",x_G
	write(6,*) "x_H is...:",x_H
	write(6,*) "x_N is...:",x_N
	write(6,*) ""
	write(6,*) "*****************************************************"
	write(6,*) "*****************************************************"
	write(6,*) ""
	endif
c	-----------------------------------------------------------------------------

c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+
c	PLASTIC ISOTROPIC HARDENING PARAMETERS
c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+

c	ISOTROPIC HARDENING Discrete Data Type
c	-----------------------------------------------------------------------------
c     INPUT YOUR FLOW CURVE, ABAQUS NOTATION, {STRESS,STRAIN}
c	-----------------------------------------------------------------------------
      flow_curve = (/709.160100d0,	0.000000d0,
     1	 814.160092d0,	0.000500d0,
     2	 870.053986d0,	0.001000d0,
     3	 931.377426d0,	0.002000d0,
     4	 969.699121d0,	0.003000d0,
     5	 997.883244d0,	0.004000d0,
     6	1020.236739d0,	0.005000d0,
     7	1091.420281d0,	0.010000d0,
     8	1133.044454d0,	0.015000d0,
     9	1161.706040d0,	0.020000d0,
     1	1182.997728d0,	0.025000d0,
     2	1199.536990d0,	0.030000d0,
     3	1212.774065d0,	0.035000d0,
     4	1223.601180d0,	0.040000d0,
     5	1232.607924d0,	0.045000d0,
     6	1240.204575d0,	0.050000d0,
     7	1246.687965d0,	0.055000d0,
     8	1252.279413d0,	0.060000d0,
     9	1257.147908d0,	0.065000d0,
     1	1261.424969d0,	0.070000d0,
     2	1265.214553d0,	0.075000d0,
     3	1268.599880d0,	0.080000d0,
     4	1271.648275d0,	0.085000d0,
     5	1274.414681d0,	0.090000d0,
     6	1276.944263d0,	0.095000d0,
     7	1279.274377d0,	0.100000d0,
     8	1296.488467d0,	0.150000d0,
     9	1309.594665d0,	0.200000d0,
     1	1321.854402d0,	0.250000d0,
     2	1333.936616d0,	0.300000d0,
     3	1358.018439d0,	0.400000d0,
     4	1382.088725d0,	0.500000d0,
     5	3428.017500d0,	9.000000d0/)	    

c	ISOTROPIC HARDENING Double Saturation Type
c	-----------------------------------------------------------------------------
 	x_K		  =0.0D0                                  
	x_K_inf1  =sigma_y_0_11  ! x_K_inf=sigma_y_0_11+Q_q  (MPa) INPRO 
	x_K_inf2  =sigma_y_0_11
	x_K_0	  =sigma_y_0_11  ! DO NOT CHANGE!!!
	x_delta1  =6.9d0         ! C_q (-)                  INPRO notation
      x_delta2  =12.0d0
      
c	ISOTROPIC HARDENING POWER TYPE
c	-----------------------------------------------------------------------------
c	x_K_prime=x_C*alpha^x_pow
c	-----------------------------------------------------------------------------
	x_C			=1042.0D0
	x_pow		=0.14D0

c	ISOTROPIC HARDENING POWER PIRES TYPE
c	-----------------------------------------------------------------------------
c	x_K_prime=x_C*alpha^x_pow
c	-----------------------------------------------------------------------------
	x_C			=1042.0D0
	x_pow		=0.14D0
	x_add		=1.0e-4

c	ISOTROPIC HARDENING POWER LUDWIK TYPE
c	-----------------------------------------------------------------------------
c	x_K_prime=x_C*alpha^x_pow
c	-----------------------------------------------------------------------------
	x_C			=1042.0D0
	x_pow		=0.14D0

c	ISOTROPIC HARDENING Saturation Type
c	-----------------------------------------------------------------------------
	x_K		=0.00000000001D0
	x_K_inf	=sigma_y_0_11+500.0d0  ! for only linear hardening
	x_K_0	=sigma_y_0_11  ! DO NOT CHANGE !!!
	x_delta	=9.0D0       ! for only linear hardening

c	ISOTROPIC HARDENING Swift Type
c	-----------------------------------------------------------------------------
      x_K_a=x_K_a_par
      x_K_b=x_K_b_par
      x_K_n=x_K_n_par

c	ISOTROPIC HARDENING Requiring Offset for Lueders
c	-----------------------------------------------------------------------------
      x_alpha_offset=0.0072D0
      sigma_y_1_11=x_K_a*((x_K_b+x_alpha_offset)**x_K_n)
      x_K_offset=(sigma_y_1_11-sigma_y_0_11)/x_alpha_offset 

c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+
c	PLASTIC KINEMATIC HARDENING PARAMETERS
c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+

c	COMMON /KINHARD/ C_X, Q_X 
c     Armstrong Frederic type nonlinear kinematic hardening
c	-----------------------------------------------------------------------------
c     C_X =props(1)    !(-)
c     Q_X =props(2)    !(MPa)
      
      C_X= 0.0d0       !(-)
      Q_X= 0.000000001d0     !(MPa) 

c     Armstrong Frederic type nonlinear kinematic hardening
c	-----------------------------------------------------------------------------
      C_X =0.0d0           !(-)
      Q_X =100000.0d0    !(MPa)

c     Armstrong Frederic type nonlinear kinematic hardening
c	-----------------------------------------------------------------------------
      C_X =30.26d0           !(-)
      Q_X =8408.3d0    !(MPa)
      
      C_X= 0.0d0       !(-)
      Q_X= 0.000000001d0     !(MPa) 

c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+
c	DAMAGE PARAMETERS
c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+

c	COMMON /LEMAITRE/ Beta, Big_S, Small_s, Y_0, D_cri, D_del
c	-----------------------------------------------------------------------------
      Beta    =Beta_par         !(-)
      Big_S   =Big_S_par         !(MPa)
      Small_s =Small_s_par           !(-)
      Y_0     =0.0D0    !(MPa)
      D_cri   =0.99D0   !(-)
      D_del   =D_cri   !(-)

c     CRACK CLOSURE PARAMETER (for 1 identical to classical Lemaitre damage model)
c	-----------------------------------------------------------------------------   
      h_crack_closure=0.15d0

c	COMMON /MAX_SHEAR_POWER/ shear_max_power=0 for classical Lemaitre
c	                         shear_max_power=1 OK
c	                         shear_max_power=limits should be understood
c	-----------------------------------------------------------------------------
      shear_max_power=0.0d0

c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+
c	RATE DEPENDENCE PARAMETERS - DOES NOT WORK CURRENTLY!!!
c	------+--------+-------------+----------------+-----------------------+-----+
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c	------+--------+-------------+----------------+-----------------------+-----+

c	COMMON /VISCOSITY/ ETA
c	-----------------------------------------------------------------------------
      ETA=zero         ! when taken zero NO rate effects are due!
      	                    
	return 
	end

c	*****************************************************************************
c	*****************************************************************************

	subroutine select_analysis_flags()
c	-----------------------------------------------------------------------------
c	Just activate the flag you wish to select
c	-----------------------------------------------------------------------------
c	TO DO LIST	: NOTHING AT THIS LEVEL, CHECKED FOR CORRECTNESS
c	-----------------------------------------------------------------------------
      include 'vaba_param.inc'

c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /SOLN_SCHEME_FLAGS/  i_SOLN_SCHEME,i_SIMUL,i_STAG
c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /ISO_HARD_FLAGS   /  i_ISO_HARD_MODEL,i_POWER,i_SATURATION
	COMMON /ISO_HARD_FLAGS2  /  i_POWER_PIRES,i_POWER_LUDWIK
	COMMON /ISO_HARD_FLAGS3  /  i_DOUBLE_SATURATION,i_SWIFT,i_DISCRETE
c	------+--------+-------------+----------------+-----------------------+-----+
	COMMON /OUTPUT_FLAGS     /  i_WANT_OUTPUT,i_YES,i_NO,i_SOME
c	------+--------+-------------+----------------+-----------------------+-----+

c	------+--------+-------------+----------------+-----------------------+-----+
c	SOLUTION FLAGS - DOES NOT WORK CURRENTLY!!!
c	------+--------+-------------+----------------+-----------------------+-----+
	i_SOLN_SCHEME   =i_SIMUL
c	i_SOLN_SCHEME   =i_STAG

c	------+--------+-------------+----------------+-----------------------+-----+
c	HARDENING FLAGS - WORKS!!!
c	------+--------+-------------+----------------+-----------------------+-----+
c	i_ISO_HARD_MODEL=i_POWER 
c	i_ISO_HARD_MODEL=i_POWER_PIRES
c	i_ISO_HARD_MODEL=i_POWER_LUDWIK 
c	i_ISO_HARD_MODEL=i_SATURATION
c	i_ISO_HARD_MODEL=i_DOUBLE_SATURATION
	i_ISO_HARD_MODEL=i_SWIFT
c	i_ISO_HARD_MODEL=i_DISCRETE
	
c	------+--------+-------------+----------------+-----------------------+-----+
c	OUTPUT FLAGS - DOES NOT WORK CURRENTLY!!!
c	------+--------+-------------+----------------+-----------------------+-----+
c	i_WANT_OUTPUT=i_YES
	i_WANT_OUTPUT=i_NO
c     i_WANT_OUTPUT=i_SOME

	return
	end

c	*****************************************************************************
c	*****************************************************************************

