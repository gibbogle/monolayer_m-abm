!----------------------------------------------------------------------------------
! Note: The value of spcrad was first determined by writing out the value computed in rkc.
! Later it was just determined by trial, then made into a run parameter.
!----------------------------------------------------------------------------------
double precision function spcrad(neqn,t,y)
!DEC$ ATTRIBUTES DLLEXPORT :: spcrad
use global
integer :: neqn
double precision :: t, y(neqn)
spcrad = spcrad_value
end function

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
module ode_diffuse

use chemokine
use metabolism
use rkc_90

implicit none

integer :: ivdbug

!real(REAL_KIND) :: work_rkc(8+5*2*MAX_CHEMO)
real(REAL_KIND) :: work_rkc(8+5*3*(N1D+1))
logical :: chemo_active(2*MAX_CHEMO)    ! flags necessity to solve for the constituent
real(REAL_KIND) :: CO2_rkc				! O2 concentration for f_rkc_drug
integer :: idrug_rkc					! drug number for f_rkc_drug

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine CheckDrugConcs
integer :: ndrugs_present, drug_present(3*MAX_DRUGTYPES), drug_number(3*MAX_DRUGTYPES)
integer :: idrug, iparent, im, kcell, ichemo, i
type(cell_type), pointer :: cp

ndrugs_present = 0
drug_present = 0
do idrug = 1,ndrugs_used
	iparent = DRUG_A + 3*(idrug-1)
	if (chemo(iparent)%present) then		! simulation with this drug has started
	    do im = 0,2
	        ichemo = iparent + im
	        ndrugs_present = ndrugs_present + 1
	        drug_present(ndrugs_present) = ichemo
	        drug_number(ndrugs_present) = idrug
	    enddo
	endif
enddo

do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
    cp => cell_list(kcell)
	do i = 1,ndrugs_present
	    ichemo = drug_present(i)
	    idrug = drug_number(i)
	    if (cp%Cin(ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
!	    if (cp%Cex(ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
	enddo
enddo
do i = 1,ndrugs_present
    ichemo = drug_present(i)
    idrug = drug_number(i)
    if (Caverage(MAX_CHEMO + ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
enddo
end subroutine

!----------------------------------------------------------------------------------
! For constituent ichemo, the extracellular concentration is:
! Cex = chemo(ichemo)%conc
! In the case of oxygen this is determined from: 
!   depth, Kdiff, chemo(OXYGEN)%flux, chemo(OXYGEN)%bdry_conc
! where:
!   depth = depth of medium in the well
!   %flux = total O2 uptake rate by cells
!   %bdry_conc = specified O2 concentration at the medium-air boundary
! For other constituents %conc is the average concentration in the medium,
! i.e. the medium is considered to be fully mixed.  In this case:
!   dC/dt = -flux/V
! where:
!   V = total volume of medium
!   C = medium concentration %Cin
!
! neqn = 2*ncvars = 2*number of constituents present
! ic > ncvars implies a medium concentration
! chemo_active(ic) = false means we do not solve for it (only medium variables)
! NOT USED NOW - we use f_rkc_OGL
!----------------------------------------------------------------------------------
subroutine f_rkc(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: ic, ichemo, idrug, im, ict, Ng, ncvars
real(REAL_KIND) :: dCsum, dCdiff, dCreact, vol_cm3, val, Cin(MAX_CHEMO), Cmedium(MAX_CHEMO), Cex
real(REAL_KIND) :: G_rate,PP_rate,P_rate
real(REAL_KIND) :: decay_rate, C, membrane_kin, membrane_kout, membrane_flux, area_factor, n_O2(0:2)
logical :: metabolised(MAX_CELLTYPES,0:2)
real(REAL_KIND) :: metab, cell_flux, dMdt, KmetC, vcell_actual, Kd(0:2), dC, C0
type(drug_type), pointer :: dp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: average_volume = 1.2
logical :: use_average_volume = .true.

write(*,*)
write(*,'(a,f9.4)') 'f_rkc: t: ',t
write(*,'(a)') '-------------------'
if (use_average_volume) then
    vol_cm3 = Vcell_cm3*average_volume	  ! not accounting for cell volume change
    area_factor = (average_volume)**(2./3.)
endif
!Vcell_actual = Vcell_cm3*cell_list(kcell)%volume
!vol_cm3 = Vcell_actual	            ! accounting for cell volume change
!Cin = cell_list(kcell)%Cin
!ict = cell_list(kcell)%celltype
ncvars = neqn/2
do ic = 1,ncvars
	ichemo = chemomap(ic)
    Cin(ichemo) = y(ic)
    Cmedium(ichemo) = y(ncvars+ic)
enddo
!write(*,'(a,8e12.3)') 'Cin: ',Cin(1:ncvars)
!write(*,'(a,8e12.3)') 'Cmedium: ',Cmedium(1:ncvars)
ict = icase
if (use_metabolism) then
	mp => metabolic(ict)
	call get_metab_rates(ict,mp,Cin)
endif
!write(*,*) 'icase, neqn: ',icase,neqn

do ic = 1,neqn
    if (ic <= ncvars) then
    	ichemo = chemomap(ic)
    else
    	ichemo = chemomap(ic-ncvars)
    endif
    if (ichemo == GLUCOSE) then
	    Ng = chemo(GLUCOSE)%Hill_N
    elseif (ichemo == LACTATE) then
	    Ng = chemo(LACTATE)%Hill_N
    endif
    if (ichemo > TRACER) then
        idrug = (ichemo - TRACER - 1)/3 + 1
        im = ichemo - TRACER - 1 - 3*(idrug-1)		! 0 = drug, 1 = metab1, 2 = metab2
        dp => drug(idrug)
        metabolised(:,:) = (dp%Kmet0(:,:) > 0)
        if (idrug > 0) then
            n_O2(:) = dp%n_O2(ict,:)
        endif
    endif

    decay_rate = chemo(ichemo)%decay_rate
    membrane_kin = chemo(ichemo)%membrane_diff_in
    membrane_kout = chemo(ichemo)%membrane_diff_out

    Cex = Cmedium(ichemo)
	C = Cin(ichemo)     ! = y(ic)
	membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
	dydt(ic) = 0
	if (ic <= ncvars .and. chemo_active(ic)) then      ! cell variable
	    dCreact = 0
	    if (use_metabolism) then
!			call get_metab_rates(ict,mp,Cin)
			if (ichemo == OXYGEN) then
				dCreact = (-mp%O_rate + membrane_flux)/vol_cm3
			elseif (ichemo == GLUCOSE) then
				dCreact = (-mp%G_rate + membrane_flux)/vol_cm3
			elseif (ichemo == LACTATE) then
				dCreact = (mp%L_rate + membrane_flux)/vol_cm3
			endif
	    else
			if (ichemo == OXYGEN) then
				metab = O2_metab(C)
				dCreact = (-metab*chemo(ichemo)%max_cell_rate + membrane_flux)/vol_cm3	! convert mass rate (mumol/s) to concentration rate (mM/s)
	!		    write(*,'(a,6e12.3)') 'O2: ',C,metab,chemo(ichemo)%max_cell_rate,membrane_flux,vol_cm3,dCreact
			elseif (ichemo == GLUCOSE) then
				metab = C**Ng/(chemo(ichemo)%MM_C0**Ng + C**Ng)
				cell_flux = metab*chemo(ichemo)%max_cell_rate
				dCreact = (-cell_flux + membrane_flux)/vol_cm3	! convert mass rate (mumol/s) to concentration rate (mM/s)
	!		    write(*,'(a,6e11.3)') 'glucose: ',C,metab,chemo(ichemo)%max_cell_rate,membrane_flux,vol_cm3,dCreact
			endif
		endif
		if (ichemo > TRACER) then
			if (im == 0) then
				if (metabolised(ict,0) .and. C > 0) then
					KmetC = dp%Kmet0(ict,0)*C
					if (dp%Vmax(ict,0) > 0) then
						KmetC = KmetC + dp%Vmax(ict,0)*C/(dp%Km(ict,0) + C)
					endif
					dCreact = -(1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2(0)/(dp%KO2(ict,0)**n_O2(0) + Cin(OXYGEN)**n_O2(0)))*KmetC
				endif
				dCreact = dCreact + membrane_flux/vol_cm3
			elseif (im == 1) then	! ichemo-1 is the PARENT drug
				if (metabolised(ict,0) .and. Cin(ichemo-1) > 0) then
					dCreact = (1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2(0)/(dp%KO2(ict,0)**n_O2(0) + Cin(OXYGEN)**n_O2(0)))*dp%Kmet0(ict,0)*Cin(ichemo-1)
				endif
				if (metabolised(ict,1) .and. C > 0) then
					dCreact = dCreact - (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)**n_O2(1)/(dp%KO2(ict,1)**n_O2(1) + Cin(OXYGEN)**n_O2(1)))*dp%Kmet0(ict,1)*C
				endif
				dCreact = dCreact + membrane_flux/vol_cm3
			elseif (im == 2) then	! ichemo-1 is the METAB1
				if (metabolised(ict,1) .and. Cin(ichemo-1) > 0) then
					dCreact = (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)**n_O2(1)/(dp%KO2(ict,1)**n_O2(1) + Cin(OXYGEN)**n_O2(1)))*dp%Kmet0(ict,1)*Cin(ichemo-1)
				endif
				if (metabolised(ict,2) .and. C > 0) then
					dCreact = dCreact - (1 - dp%C2(ict,2) + dp%C2(ict,2)*dp%KO2(ict,2)**n_O2(2)/(dp%KO2(ict,2)**n_O2(2) + Cin(OXYGEN)**n_O2(2)))*dp%Kmet0(ict,2)*C
				endif
				dCreact = dCreact + membrane_flux/vol_cm3
			endif
		endif
        dydt(ic) = dCreact - C*decay_rate
        if (ichemo <= LACTATE) then
			write(*,'(a,i2,e12.3)') 'cell dC/dt: ',ichemo,dydt(ic)
		endif
    else    ! medium variable 
        if (chemo_active(ic)) then
            dydt(ic) = -Ncells*membrane_flux/total_volume - Cex*decay_rate
        else
            dydt(ic) = 0
        endif
        if (ichemo <= LACTATE) then
			write(*,'(a,i2,e12.3)') 'medium dC/dt: ',ichemo,dydt(ic)
		endif
    endif
	if (isnan(dydt(ic))) then
		write(nflog,*) 'f_rkc: dydt isnan: ',ic,ichemo,dydt(ic)
		write(*,*) 'f_rkc: dydt isnan: ',ic,ichemo,dydt(ic)
		stop
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine f_rkc_drug(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: k, kk, i, ichemo, idrug, iparent, im, ict, Nlivecells
real(REAL_KIND) :: dCsum, dCdiff, dCreact, vol_cm3, Cex
real(REAL_KIND) :: decay_rate, C, membrane_kin, membrane_kout, membrane_flux, area_factor, n_O2(0:2)
logical :: metabolised(MAX_CELLTYPES,0:2)
real(REAL_KIND) :: metab, cell_flux, dMdt, KmetC, vcell_actual, dC, CO2, A, d, dX, dV, Kd, KdAVX
type(drug_type), pointer :: dp
real(REAL_KIND) :: average_volume = 1.2
logical :: use_average_volume = .true.
logical :: is_metab1

ict = icase
CO2 = CO2_rkc
idrug = idrug_rkc
A = well_area
d = total_volume/A
dX = d/N1D
dV = A*dX
if (use_average_volume) then
    vol_cm3 = Vcell_cm3*average_volume	  ! not accounting for cell volume change
    area_factor = (average_volume)**(2./3.)
endif
Nlivecells = Ncells - (Ncells_dying(1) + Ncells_dying(2))
iparent = DRUG_A + 3*(idrug-1)
dp => drug(idrug)
metabolised(:,:) = (dp%Kmet0(:,:) > 0)
n_O2(:) = dp%n_O2(ict,:)

k = 0
do im = 0,2
	! First process IC reactions
	k = k+1
	C = y(k)
	Cex = y(k+1)
	ichemo = iparent + im
	Kd = chemo(ichemo)%medium_diff_coef
    decay_rate = chemo(ichemo)%decay_rate
    membrane_kin = chemo(ichemo)%membrane_diff_in
    membrane_kout = chemo(ichemo)%membrane_diff_out
	membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
	dCreact = 0
    if (im == 0) then
        if (metabolised(ict,0) .and. C > 0) then
		    KmetC = dp%Kmet0(ict,0)*C
		    if (dp%Vmax(ict,0) > 0) then
			    KmetC = KmetC + dp%Vmax(ict,0)*C/(dp%Km(ict,0) + C)
		    endif
		    dCreact = -(1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2(0)/(dp%KO2(ict,0)**n_O2(0) + CO2**n_O2(0)))*KmetC
	    endif
!		write(nflog,'(a,3e12.3)') 'dCreact, flux, decay: ',dCreact,membrane_flux/vol_cm3,-C*decay_rate
	    dCreact = dCreact + membrane_flux/vol_cm3
    elseif (im == 1) then	! kk=1 is the PARENT drug
		kk = 1
	    if (metabolised(ict,0) .and. y(kk) > 0) then
		    dCreact = (1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2(0)/(dp%KO2(ict,0)**n_O2(0) + CO2**n_O2(0)))*dp%Kmet0(ict,0)*y(kk)
	    endif
	    if (metabolised(ict,1) .and. C > 0) then
		    dCreact = dCreact - (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)**n_O2(1)/(dp%KO2(ict,1)**n_O2(1) + CO2**n_O2(1)))*dp%Kmet0(ict,1)*C
	    endif
	    dCreact = dCreact + membrane_flux/vol_cm3
    elseif (im == 2) then	! kk=N1D+2 is the METAB1
		kk = N1D+2
	    if (metabolised(ict,1) .and. y(kk) > 0) then
		    dCreact = (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)**n_O2(1)/(dp%KO2(ict,1)**n_O2(1) + CO2**n_O2(1)))*dp%Kmet0(ict,1)*y(kk)
	    endif
	    if (metabolised(ict,2) .and. C > 0) then
		    dCreact = dCreact - (1 - dp%C2(ict,2) + dp%C2(ict,2)*dp%KO2(ict,2)**n_O2(2)/(dp%KO2(ict,2)**n_O2(2) + CO2**n_O2(2)))*dp%Kmet0(ict,2)*C
	    endif
	    dCreact = dCreact + membrane_flux/vol_cm3
    endif
	dydt(k) = dCreact - C*decay_rate
!	write(nflog,'(a,i4,e12.3)') 'dydt: ',im,dydt(k)
	if (isnan(dydt(k))) then
		write(nflog,*) 'f_rkc_drug: dydt isnan: ',im,dydt(k)
		write(*,*) 'f_rkc_drug: dydt isnan: ',im,dydt(k)
		stop
	endif
	
	! Next process grid cell next to the cell layer - note that membrane _flux has already been computed
	k = k+1
	C = y(k)
	dydt(k) = (-Nlivecells*membrane_flux - Kd*A*(C - y(k+1))/dX)/dV - C*decay_rate
	
	! Next compute diffusion and decay on the FD grid
	KdAVX = Kd*A/(dV*dX)
	do i = 2,N1D
		k = k+1
		C = y(k)
		if (i < N1D) then
			dydt(k) = KdAVX*(y(k+1) - 2*C + y(k-1)) - C*decay_rate
		else
			dydt(k) = KdAVX*(-C + y(k-1)) - C*decay_rate
		endif
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine f_rkc_OGL(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: k, kk, i, ichemo, ict, Nlivecells
real(REAL_KIND) :: dCsum, dCdiff, dCreact, vol_cm3, Cex, Cin(3)
real(REAL_KIND) :: C, membrane_kin, membrane_kout, membrane_flux, area_factor, Cbnd
real(REAL_KIND) :: A, d, dX, dV, Kd, KdAVX
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: average_volume = 1.2
logical :: use_average_volume = .true.

ict = icase
A = well_area
d = total_volume/A
dX = d/N1D
dV = A*dX
if (use_average_volume) then
    vol_cm3 = Vcell_cm3*average_volume	  ! not accounting for cell volume change
    area_factor = (average_volume)**(2./3.)
endif
Nlivecells = Ncells - (Ncells_dying(1) + Ncells_dying(2))
Cin(1) = y(1)
Cin(2) = y(N1D+2)
Cin(3) = y(2*N1D+3)
mp => metabolic(ict)
call get_metab_rates(ict,mp,Cin)

k = 0
do ichemo = 1,3
	! First process IC reactions
	k = k+1
	C = y(k)
	Cex = y(k+1)
	Kd = chemo(ichemo)%medium_diff_coef
    membrane_kin = chemo(ichemo)%membrane_diff_in
    membrane_kout = chemo(ichemo)%membrane_diff_out
	membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
    if (ichemo == OXYGEN) then
		dCreact = (-mp%O_rate + membrane_flux)/vol_cm3
    elseif (ichemo == GLUCOSE) then	! 
		dCreact = (-mp%G_rate + membrane_flux)/vol_cm3
    elseif (ichemo == LACTATE) then
		dCreact = (mp%L_rate + membrane_flux)/vol_cm3
		dCreact = dCreact*f_MM(C,chemo(LACTATE)%MM_C0, int(chemo(LACTATE)%Hill_N))
    endif
	dydt(k) = dCreact
!	write(nflog,'(a,i4,e12.3)') 'dydt: ',im,dydt(k)
	if (isnan(dydt(k))) then
		write(nflog,*) 'f_rkc_drug: dydt isnan: ',ichemo,dydt(k)
		write(*,*) 'f_rkc_drug: dydt isnan: ',ichemo,dydt(k)
		stop
	endif
	
	! Next process grid cell next to the cell layer - note that membrane _flux has already been computed
	k = k+1
	C = y(k)
	dydt(k) = (-Nlivecells*membrane_flux - Kd*A*(C - y(k+1))/dX)/dV
	
	! Next compute diffusion and decay on the FD grid
	! Need special treatment for oxygen at air boundary
	KdAVX = Kd*A/(dV*dX)
	do i = 2,N1D
		k = k+1
		C = y(k)
		if (i < N1D) then
			dydt(k) = KdAVX*(y(k+1) - 2*C + y(k-1))
		else
			if (ichemo == OXYGEN) then
				Cbnd = chemo(OXYGEN)%bdry_conc
				dydt(k) = KdAVX*(Cbnd - 2*C + y(k-1))
			else
				dydt(k) = KdAVX*(-C + y(k-1))
			endif
		endif
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine Solver1(it,tstart,dt,nc,ok)
integer :: it, nc
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, ic, k, ict, ncvars, neqn, kcell
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: C(2*MAX_CHEMO)
real(REAL_KIND) :: timer1, timer2
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)
logical :: solve_O2 = .true.
logical :: use_drugsolver = .true.

k = 0
do ic = 1,nchemo
	ichemo = chemomap(ic)
	k = k + 1
    chemo_active(k) = .not.chemo(ichemo)%constant
    if (ichemo == OXYGEN) then
        if (.not.solve_O2) then
            ! Suppress solving for cell oxygen
!	        Caverage(OXYGEN) = getCin(OXYGEN,Caverage(MAX_CHEMO+OXYGEN))
            chemo_active(k) = .false.
        endif
    endif
	if (use_drugsolver .and. ichemo >= DRUG_A) then
        chemo_active(k) = .false.
	endif	
	C(k) = Caverage(ichemo)                ! average cell concentration
enddo
ncvars = k
! Note: ncvars = nchemo
do ic = 1,nchemo
	ichemo = chemomap(ic)
	k = k + 1
	C(k) = Caverage(MAX_CHEMO + ichemo)      ! average medium concentration
    chemo_active(k) = .not.chemo(ichemo)%constant
    if (ichemo == OXYGEN) then
        ! Suppress solving for medium oxygen
        chemo_active(k) = .false.
    endif
enddo
neqn = k
! Note: neqn = 2*ncvars

!write(*,*) 'solver: nchemo,neqn: ',nchemo,neqn
!write(*,'(10f7.3)') C(1:neqn)
!write(*,'(a,3f8.5)') 'solver: metab1: ',Caverage(MAX_CHEMO+4:MAX_CHEMO+6)
!if (chemo(DRUG_A)%present) then
!	write(*,'(a,3e12.3)') 'medium drug conc: ',Caverage(MAX_CHEMO +DRUG_A:MAX_CHEMO +DRUG_A+2)
!endif

ict = 1 ! for now just a single cell type

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-2
atol = rtol

idid = 0
t = tstart
tend = t + dt
call rkc(comm_rkc(1),neqn,f_rkc,C,t,tend,rtol,atol,info,work_rkc,idid,ict)
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif

! This determines average cell concentrations, assumed the same for all cells
! Now put the concentrations into the cells

k = 0
do ic = 1,nchemo
    ichemo = chemomap(ic)
    k = k + 1
    Caverage(ichemo) = C(k)
    do kcell = 1,nlist
        if (cell_list(kcell)%state == DEAD) cycle
        cell_list(kcell)%Cin(ichemo) = Caverage(ichemo)
    enddo
enddo
k = ncvars
do ic = 1,nchemo
    ichemo = chemomap(ic)
    k = k + 1
    Caverage(MAX_CHEMO + ichemo) = C(k)
enddo

!if (chemo(DRUG_A)%present) then
!	write(*,'(a,3e12.3)') 'did solver: drug conc: ',Caverage(DRUG_A:DRUG_A+2)
!endif
!write(*,'(a,3f8.5)') 'did solver: metab1: ',Caverage(MAX_CHEMO+4:MAX_CHEMO+6)
! Note: medium oxygen is unchanged

end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine Solver2(it,tstart,dt,nc,ok)
integer :: it, nc
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, ic, k, ict, ncvars, neqn, kcell
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: C(2*MAX_CHEMO)
real(REAL_KIND) :: timer1, timer2
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)
logical :: solve_O2 = .true.
logical :: use_drugsolver = .true.

k = 0
do ic = 1,nchemo
	ichemo = chemomap(ic)
	k = k + 1
    chemo_active(k) = .not.chemo(ichemo)%constant
    if (ichemo == OXYGEN) then
        if (.not.solve_O2) then
            ! Suppress solving for cell oxygen
!	        Caverage(OXYGEN) = getCin(OXYGEN,Caverage(MAX_CHEMO+OXYGEN))
            chemo_active(k) = .false.
        endif
    endif
	if (use_drugsolver .and. ichemo >= DRUG_A) then
        chemo_active(k) = .false.
	endif	
	C(k) = Caverage(ichemo)                ! average cell concentration
enddo
ncvars = k
! Note: ncvars = nchemo
do ic = 1,nchemo
	ichemo = chemomap(ic)
	k = k + 1
	C(k) = Caverage(MAX_CHEMO + ichemo)      ! average medium concentration
    chemo_active(k) = .not.chemo(ichemo)%constant
    if (ichemo == OXYGEN) then
        ! Suppress solving for medium oxygen
        chemo_active(k) = .false.
    endif
	if (use_drugsolver .and. ichemo >= DRUG_A) then
        chemo_active(k) = .false.
	endif	
enddo
neqn = k
! Note: neqn = 2*ncvars

!write(*,*) 'solver: nchemo,neqn: ',nchemo,neqn
!write(*,'(10f7.3)') C(1:neqn)
!write(*,'(a,3f8.5)') 'solver: metab1: ',Caverage(MAX_CHEMO+4:MAX_CHEMO+6)

ict = 1 ! for now just a single cell type

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-2
atol = rtol

idid = 0
t = tstart
tend = t + dt
call rkc(comm_rkc(1),neqn,f_rkc,C,t,tend,rtol,atol,info,work_rkc,idid,ict)
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif

! This determines average cell concentrations, assumed the same for all cells
! Now put the concentrations into the cells

k = 0
do ic = 1,nchemo
    ichemo = chemomap(ic)
    k = k + 1
    if (.not.chemo_active(k)) cycle
    Caverage(ichemo) = C(k)
    do kcell = 1,nlist
        if (cell_list(kcell)%state == DEAD) cycle
        cell_list(kcell)%Cin(ichemo) = Caverage(ichemo)
    enddo
enddo
k = ncvars
do ic = 1,nchemo
    ichemo = chemomap(ic)
    k = k + 1
    if (.not.chemo_active(k)) cycle
    Caverage(MAX_CHEMO + ichemo) = C(k)
enddo

if (.not.use_drugsolver) return
if (chemo(DRUG_A)%present) then
	call DrugSolver(DRUG_A,tstart,dt,1,ok)
endif
if (chemo(DRUG_B)%present) then
	call DrugSolver(DRUG_B,tstart,dt,2,ok)
endif
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine Solver(it,tstart,dt,nc,ok)
integer :: it, nc
real(REAL_KIND) :: tstart, dt
logical :: ok
logical :: use_drugsolver = .true.

call OGLSolver(tstart,dt,ok)
if (.not.use_drugsolver) return
if (chemo(DRUG_A)%present) then
	call DrugSolver(DRUG_A,tstart,dt,1,ok)
endif
if (chemo(DRUG_B)%present) then
	call DrugSolver(DRUG_B,tstart,dt,2,ok)
endif
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine DrugSolver(iparent,tstart,dt,idrug,ok)
integer :: iparent, idrug
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, k, ict, neqn, i, kcell, im
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: C(3*N1D+3), Csum
real(REAL_KIND) :: timer1, timer2
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)

!write(nflog,*) 'DrugSolver: ',istep
ict = 1 ! for now just a single cell type
idrug_rkc = idrug

k = 0
do im = 0,2
	ichemo = iparent + im
	if (.not.chemo(ichemo)%present) cycle
	k = k+1
	C(k) = Caverage(ichemo)		! IC 
	do i = 1,N1D
		k = k+1
!		C(k) = Cdrug(im,i)		! EC
		C(k) = chemo(ichemo)%Cmedium(i)
	enddo
enddo

neqn = k

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-5
atol = rtol

idid = 0
t = tstart
tend = t + dt
call rkc(comm_rkc(1),neqn,f_rkc_drug,C,t,tend,rtol,atol,info,work_rkc,idid,ict)
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif
!write(nflog,'(a,3e12.3)') 'IC: ',C(1),C(N1D+2),C(2*N1D+3)

! This determines average cell concentrations, assumed the same for all cells
! Now put the concentrations into the cells 

do im = 0,2
    ichemo = iparent + im
	if (.not.chemo(ichemo)%present) cycle
    k = im*(N1D+1) + 1
    Caverage(ichemo) = C(k)
	Csum = 0
    do i = 1,N1D
		Csum = Csum + C(k+i)
	enddo
	Cmediumave(ichemo) = Csum/N1D
    do kcell = 1,nlist
        if (cell_list(kcell)%state == DEAD) cycle
        cell_list(kcell)%Cin(ichemo) = Caverage(ichemo)
    enddo
    Caverage(MAX_CHEMO + ichemo) = C(k+1)	! not really average, this is medium at the cell layer, i.e. EC
!	write(nflog,'(a,i3,5e12.3)') 'Cdrug: im: ',im,Cdrug(im,1:5)
enddo
!write(*,'(a,3e12.3)') 'Cell drug conc: ',(Caverage(DRUG_A+k),k=0,2)

end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine OGLSolver(tstart,dt,ok)
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, k, ict, neqn, i, kcell
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: C(3*N1D+3), Cin(3), Csum, Itotal, I2Divide
real(REAL_KIND) :: timer1, timer2
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)
type(metabolism_type), pointer :: mp
type(cell_type), pointer :: cp
real(REAL_KIND) :: Cic,Cex,average_volume,area_factor,membrane_kin,membrane_kout,membrane_flux

!write(nflog,*) 'DrugSolver: ',istep
ict = selected_celltype ! for now just a single cell type

k = 0
do ichemo = 1,3
!	if (.not.chemo(ichemo)%present) cycle
	k = k+1
	C(k) = Caverage(ichemo)		! IC 
	do i = 1,N1D
		k = k+1
		C(k) = C_OGL(ichemo,i)	! EC
	enddo
enddo

neqn = k

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-5
atol = rtol

idid = 0
t = tstart
tend = t + dt
call rkc(comm_rkc(1),neqn,f_rkc_OGL,C,t,tend,rtol,atol,info,work_rkc,idid,ict)
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif
!write(nflog,'(a,3e12.3)') 'IC: ',C(1),C(N1D+2),C(2*N1D+3)
!write(*,'(a,3e12.3)') 'IC: ',C(1),C(N1D+2),C(2*N1D+3)

! This determines average cell concentrations, assumed the same for all cells
! Now put the concentrations into the cells 

do ichemo = 1,3
	if (.not.chemo(ichemo)%present) cycle
    k = (ichemo-1)*(N1D+1) + 1
    Caverage(ichemo) = C(k)
    Csum = 0
    do i = 1,N1D
		C_OGL(ichemo,i) = C(k+i)
		Csum = Csum + C(k+i)
		chemo(ichemo)%Cmedium(i) = C(k+i)
	enddo
	Cmediumave(ichemo) = Csum/N1D
    do kcell = 1,nlist
        if (cell_list(kcell)%state == DEAD) cycle
        cell_list(kcell)%Cin(ichemo) = Caverage(ichemo)
    enddo
    Caverage(MAX_CHEMO + ichemo) = C(k+1)	! not really average, this is medium at the cell layer, i.e. EC
!	write(nflog,'(a,i3,5e12.3)') 'Cdrug: im: ',im,Cdrug(im,1:5)
enddo
!write(*,'(a,3e12.3)') 'C O2: ',C(1),C(2),C(N1D+1)
Cin(1) = C(1)
Cin(2) = C(N1D+2)
Cin(3) = C(2*N1D+3)
mp => metabolic(ict)
call get_metab_rates(ict,mp,Cin)
do kcell = 1,nlist
	cp => cell_list(kcell)
    if (cp%state == DEAD .or. cp%state == DYING) cycle
! First back up cell metabolism parameters that we need to preserve
	Itotal = cp%metab%Itotal
	I2Divide = cp%metab%I2Divide
    cp%metab = metabolic(ict)
    cp%metab%Itotal = Itotal
    cp%metab%I2Divide = I2Divide
enddo

! Check Lactate flux balance
ichemo = LACTATE
Cic = Caverage(LACTATE)
Cex = Caverage(MAX_CHEMO + LACTATE)
average_volume = 1.2
area_factor = (average_volume)**(2./3.)
membrane_kin = chemo(ichemo)%membrane_diff_in
membrane_kout = chemo(ichemo)%membrane_diff_out
membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*Cic)
!write(*,'(a,3e12.3)') 'Lactate flux: ',mp%L_rate,membrane_flux,2*(1-N_GI(1))*mp%G_rate-mp%P_rate
! Checks OK

end subroutine

!----------------------------------------------------------------------------------
! Note: This computes a rate of change of concentration! mM/s
! Currently only for O2!!! 
! There are two options: use_Cex_Cin = true/false
!
! use_Cex_Cin = true
! ------------------
! The idea is that the speed of the intracellular reactions, compared with other
! processes, is so fast that effectively the intracellular concentration is always
! in equilibrium with the extracellular value.  This means that the rate of consumption
! in the cell matches the rate of transport across the cell membrane: both these rates 
! depend on Cin, therefore we can solve for Cin given Cex then deduce uptake rate
!
! use_Cex_Cin = false
! -------------------
! In this case we just use Cin = Cex to calculate the consumption rate - no
! dependence on chemo(OXYGEN)%membrane_diff
!----------------------------------------------------------------------------------
real(REAL_KIND) function UptakeRate(ichemo,Cex)
integer :: ichemo
real(REAL_KIND) :: Cex
real(REAL_KIND) :: vol, K1, Cin, flux
integer :: n, i

if (ichemo == OXYGEN) then
!    vol = Vsite_cm3
!    vol = Vsite_cm3 - Vextra_cm3	! this was used in the RKC solution
    vol = Vextra_cm3	! the current extracellular volume should be used I think !!!!!!!!!!!!!!!
	if (use_Cex_Cin) then
		Cin = getCin(ichemo,Cex)
!		flux = chemo(ichemo)%membrane_diff*(Cex - Cin)
		flux = (chemo(ichemo)%membrane_diff_in*Cex - chemo(ichemo)%membrane_diff_out*Cin)
	else	! 
		flux = O2_metab(Cex)*chemo(ichemo)%max_cell_rate
	endif
	if (dbug) write(nfout,'(a,2e12.4)') 'Cex, flux: ',Cex,flux
	UptakeRate = flux/vol	! concentration rate (mM/s)
else
	write(logmsg,*) 'ERROR: UptakeRate: currently only for OXYGEN'
	call logger(logmsg)
	stop
endif
end function

!----------------------------------------------------------------------------------
! Computes intracellular O2 concentration as a function of the extracellular level C,
! assuming equilibrium, i.e. rate of consumption = rate of membrane transport.
! Note that the cell's O2 uptake rate is taken to be independent of any other factors,
! e.g. independent of cell size.
! NOTE: Currently only for OXYGEN and GLUCOSE - OK because membrane_diff_in = membrane_diff_out
! Note: needs to be amended to account for HIF-1
!----------------------------------------------------------------------------------
!real(REAL_KIND) function getCinO2(C)
real(REAL_KIND) function getCin(ichemo,C)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: K1, K2, K2K1, C0, a, b, cc, D, r(3), Cin
integer :: i, n

if (ichemo >= DRUG_A) then
	write(logmsg,*) 'ERROR: getCin: currently only for OXYGEN, GLUCOSE, LACTATE'
	call logger(logmsg)
	stop
endif
!ichemo = OXYGEN
!K1 = chemo(OXYGEN)%membrane_diff*(Vsite_cm3 - Vextra_cm3)
K1 = chemo(ichemo)%membrane_diff_in
K2 = chemo(ichemo)%max_cell_rate
K2K1 = K2/K1
C0 = chemo(ichemo)%MM_C0
if (chemo(ichemo)%Hill_N == 2) then
	a = K2K1 - C
	b = C0*C0
	cc = -b*C
	call cubic_roots(a,b,cc,r,n)
	if (n == 1) then
		Cin = r(1)
	else
		n = 0
		do i = 1,3
			if (r(i) > 0) then
				n = n+1
				Cin = r(i)
			endif
		enddo
		if (n > 1) then
			write(nflog,*) 'getCin: two roots > 0: ',r
			stop
		endif
	endif
elseif (chemo(ichemo)%Hill_N == 1) then
	b = K2K1 + C0 - C
	cc = -C0*C
	D = sqrt(b*b - 4*cc)
	Cin = (D - b)/2
endif
getCin = Cin
end function


!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine UpdateCbnd_1D
integer :: kpar = 0
integer :: i, ic, ichemo
real(REAL_KIND) :: tnow, alpha_Cbnd = 0.3
real(REAL_KIND) :: t_buffer = 3600	! one hour delay before applying smoothing to Cbnd
integer :: ndrugs_present, drug_present(3*MAX_DRUGTYPES), drug_number(3*MAX_DRUGTYPES)
integer :: idrug, iparent, im
logical :: present

tnow = istep*DELTA_T
ndrugs_present = 0
drug_present = 0
do idrug = 1,ndrugs_used
	iparent = DRUG_A + 3*(idrug-1)
	if (chemo(iparent)%present) then		! simulation with this drug has started
	    do im = 0,2
	        ichemo = iparent + im
	        ndrugs_present = ndrugs_present + 1
	        drug_present(ndrugs_present) = ichemo
	        drug_number(ndrugs_present) = idrug
	    enddo
	endif
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! 1D FD solution
! uptake_rate is the total rate summed over all cells in mumol/s
! flux/vol_cm3	to convert mass rate (mumol/s) to concentration rate (mM/s)
! NOT USED NOW
!-----------------------------------------------------------------------------------------
subroutine SolveMediumGlucose(dt)
real(REAL_KIND) :: dt
real(REAL_KIND) :: A, d, dX, Kd, dV, area_factor, membrane_kin, membrane_kout
real(REAL_KIND) :: C, Cex, membrane_flux, uptake_rate, F(N1D+1)
integer :: ichemo, i, k
integer :: ndt = 20
real(REAL_KIND) :: average_volume = 1.2
real(REAL_KIND), dimension(:), pointer :: Cglucose

!write(*,*) 'SolveMediumGlucose: ',dt
ichemo = GLUCOSE
!Cglucose(:) = chemo(ichemo)%Cmedium(:)
Cglucose => chemo(ichemo)%Cmedium
Kd = chemo(ichemo)%medium_diff_coef
membrane_kin = chemo(ichemo)%membrane_diff_in
membrane_kout = chemo(ichemo)%membrane_diff_out
A = well_area
d = total_volume/A
dX = d/N1D
dV = A*dX
area_factor = (average_volume)**(2./3.)
Cex = Caverage(MAX_CHEMO + ichemo)
do k = 1,ndt
!	Cex = Cglucose(1)
!	C = Caverage(ichemo)
	C = getCin(ichemo,Cex)
	membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
	uptake_rate = Ncells*membrane_flux
	F(1) = -uptake_rate
	do i = 2,N1D
		F(i) = Kd*A*(Cglucose(i-1) - Cglucose(i))/dX
	enddo
	F(N1D+1) = 0
	do i = 1,N1D
		Cglucose(i) = Cglucose(i) + (F(i) - F(i+1))*(dt/ndt)/dV
	enddo
	Cex = Cglucose(1)
enddo
!write(nflog,'(6e12.3)') F(1),C,(Cglucose(i),i=1,4)
!write(nflog,'(10e12.3)') (Cglucose(i),i=1,N1D)
Caverage(MAX_CHEMO + ichemo) = Cex
!chemo(ichemo)%Cmedium(:) = Cglucose(:)
end subroutine

end module

