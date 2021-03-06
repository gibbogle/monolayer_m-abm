! Cancer cell state development

module cellstate
use global
use chemokine
use metabolism
use cycle_mod
implicit none

integer :: kcell_dividing = 0
logical :: first

contains

!-----------------------------------------------------------------------------------------
! Need to initialize site and cell concentrations when a cell divides and when there is
! cell death.
! dt seconds
! radiation dose has been separated out
!-----------------------------------------------------------------------------------------
subroutine GrowCells(dt,t_simulation,ok)
real(REAL_KIND) :: dt, t_simulation
logical :: ok
logical :: changed
integer :: kcell, idrug, ichemo
type(cell_type),pointer :: cp

!call logger('GrowCells: ')
!tnow = istep*DELTA_T
tnow = t_simulation		! now = time at the start of the timestep
!write(*,*) 'tnow, t_simulation: ',tnow,t_simulation
ok = .true.
call new_grower(dt,changed,ok)
if (.not.ok) return
!if (dose > 0) then
!	call Irradiation(dose, ok)
!	if (.not.ok) return
!endif
call CellDeath(dt,ok)
if (.not.ok) return
cp => cell_list(kcell_test)
!write(nflog,'(a,2i6,2f8.3,e12.3)') 'kcell_test: ',kcell_test,cp%phase,tnow/3600,cp%mitosis,cp%dVdt
end subroutine

!-----------------------------------------------------------------------------------------
! The O2 concentration to use with cell kcell is either the intracellular concentration,
! or if use_extracellular_O2, the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getO2conc(cp, C_O2)
type(cell_type), pointer :: cp
real(REAL_KIND) :: Cex, C_O2

if (use_extracellular_O2 .and. istep > 1) then		! fix 30/04/2015
    C_O2 = chemo(OXYGEN)%conc
else
    C_O2 = cp%Cin(OXYGEN)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The glucose concentration to use with cell kcell is either the intracellular concentration,
! or if use_extracellular_O2 (!), the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getGlucoseconc(cp, C_glucose)
type(cell_type), pointer :: cp
real(REAL_KIND) :: C_glucose

C_glucose = cp%Cin(GLUCOSE)
end subroutine

!-----------------------------------------------------------------------------------------
! Irradiate cells with dose (and duration tmin to be added)
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, site(3), iv, ityp, idrug, im, ichemo, kpar=0
real(REAL_KIND) :: C_O2, SER, p_death, p_recovery, R, kill_prob, tmin
real(REAL_KIND) :: SER_OER(2)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

ok = .true.
call logger('Irradiation')
!tnow = istep*DELTA_T	! seconds
if (use_volume_method) then
    do kcell = 1,nlist
        if (colony_simulation) then
            cp => ccell_list(kcell)
        else
            cp => cell_list(kcell)
        endif
	    if (cp%state == DEAD) cycle
	    if (cp%radiation_tag) cycle	! we do not tag twice (yet)
	    call getO2conc(cp,C_O2)
	    ! Compute sensitisation SER
	    if (Ndrugs_used > 0) then
	        SER = getSER(cp,C_O2)
	    else
    	    SER = 1.0
    	endif
	    ityp = cp%celltype
	    call get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
	    kill_prob = 1 - p_recovery
	    R = par_uni(kpar)
	    if (R < kill_prob) then
		    cp%radiation_tag = .true.
		    Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
		    cp%p_rad_death = p_death
		    if (LQ(ityp)%growth_delay_N > 0 .and. cp%Iphase) then
			    cp%growth_delay = .true.
			    cp%dt_delay = LQ(ityp)%growth_delay_factor*dose
			    cp%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
		    else
			    cp%growth_delay = .false.
		    endif
	    elseif (use_radiation_growth_delay_all .and. LQ(ityp)%growth_delay_N > 0) then
		    cp%growth_delay = .true.
		    cp%dt_delay = LQ(ityp)%growth_delay_factor*dose
		    cp%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
	    else
		    cp%growth_delay = .false.
	    endif
    enddo
else
    ccp => cc_parameters
    tmin = 1.0      ! for now...
    do kcell = 1,nlist
        if (colony_simulation) then
            cp => ccell_list(kcell)
        else
            cp => cell_list(kcell)
        endif
	    if (cp%state == DEAD) cycle
	    call getO2conc(cp,C_O2)
	    ! Compute sensitisation SER
	    if (Ndrugs_used > 0) then
	        SER = getSER(cp,C_O2)
	    else
    	    SER = 1.0
    	endif
        SER_OER(1) = SER*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)      ! OER_alpha
        SER_OER(2) = SER*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)      ! OER_beta
        call radiation_damage(cp, ccp, dose, SER_OER, tmin)
    enddo
endif
call logger('did Irradiation')
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function getSER(cp,C_O2) result(SER)
type(cell_type), pointer :: cp
real(REAL_KIND) ::C_O2, SER
real(REAL_KIND) :: Cs							! concentration of radiosensitising drug
real(REAL_KIND) :: SER_max0, SER_Km, SER_KO2	! SER parameters of the drug
real(REAL_KIND) :: SERmax						! max sensitisation at the drug concentration
integer :: ityp, idrug, iparent, ichemo, im

ityp = cp%celltype
SER = 1.0
do idrug = 1,Ndrugs_used
!    ichemo = 4 + 3*(idrug-1)
    iparent = DRUG_A + 3*(idrug-1)
    if (.not.chemo(ichemo)%present) cycle
    do im = 0,2
!	    ichemo = 4 + 3*(idrug-1) + im
	    ichemo = iparent + im
	    if (drug(idrug)%sensitises(ityp,im)) then
		    Cs = cp%Cin(ichemo)	! concentration of drug/metabolite in the cell
		    SER_max0 = drug(idrug)%SER_max(ityp,im)
		    SER_Km = drug(idrug)%SER_Km(ityp,im)
		    SER_KO2 = drug(idrug)%SER_KO2(ityp,im)
		    SERmax = (Cs*SER_max0 + SER_Km)/(Cs + SER_Km)
		    SER = SER*(C_O2 + SER_KO2*SERmax)/(C_O2 + SER_KO2)
	    endif
    enddo
enddo
end function

!-----------------------------------------------------------------------------------------
! A cell that receives a dose of radiation either recovers completely before reaching 
! mitosis or retains damage that has a probability of causing cell death during mitosis.
! A damaged cell that does not die at this point passes the damage on to the progeny cells.
! The probability of complete recovery = p_recovery = p_r
! The probability of death for a damaged cell at mitosis = p_death = p_d
! To ensure that the short-term death probability is consistent with the previous
! LQ formulation, we require p_d(1-p_r) = kill_prob as previously calculated.
! If p_d is determined (currently it is fixed), then 1-p_r = kill_prob/p_d,
! therefore p_r = 1 - kill_prob/p_d
!-----------------------------------------------------------------------------------------
subroutine get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, SER, p_recovery, p_death
real(REAL_KIND) :: OER_alpha_d, OER_beta_d, expon, kill_prob_orig

OER_alpha_d = dose*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
OER_beta_d = dose*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)

OER_alpha_d = OER_alpha_d*SER
OER_beta_d = OER_beta_d*SER

expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_beta_d**2
p_recovery = exp(-expon)	! = SF
p_death = LQ(ityp)%death_prob
end subroutine

!-----------------------------------------------------------------------------------------
! This is the probability of death at time of division of cell that received a dose of 
! radiation and did not recover.
! In general one would expect this to depend of damage, i.e. on dose and C_O2, but
! for now it is a constant for a cell type
!-----------------------------------------------------------------------------------------
subroutine get_pdeath(ityp,dose,C_O2,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, p_death

p_death = LQ(ityp)%death_prob
end subroutine


!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia or aglucosia, or they can be tagged 
! for death at division time if the drug is effective.
! Note: if simulating colony, no tagging, no death from anoxia, aglucosia
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3), i, ichemo, idrug, im, ityp, killmodel, kpar=0 
real(REAL_KIND) :: C_O2, C_glucose, Cdrug, n_O2, kmet, Kd, dMdt, kill_prob, dkill_prob, death_prob,survival_prob
logical :: anoxia_death, aglucosia_death
real(REAL_KIND) :: delayed_death_prob(MAX_CELLTYPES)
type(drug_type), pointer :: dp
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

!call logger('CellDeath')
ok = .true.
if (colony_simulation) then
    return
endif

!tnow = istep*DELTA_T	! seconds
anoxia_death = chemo(OXYGEN)%controls_death
aglucosia_death = chemo(GLUCOSE)%controls_death
ccp => cc_parameters
delayed_death_prob = ccp%apoptosis_rate*DELTA_T/3600
nlist0 = nlist
first = .true.
do kcell = 1,nlist
    cp => cell_list(kcell)
	ityp = cp%celltype
	if (cp%state == DEAD) cycle
	if (cp%state == DYING) then
	    if (par_uni(kpar) < delayed_death_prob(ityp)) then	
			call CellDies(kcell)
		endif
		cycle
	endif	
	call getO2conc(cp,C_O2)
	if (use_metabolism) then
!		if (cp%metab%A_rate*cp%V < cp%ATP_rate_factor*ATPs(ityp)*Vcell_cm3) then
		if (cp%metab%A_rate < ATPs) then
!			if (kcell == 1) write(nflog,'(a,2e12.3)') 'insufficient A_rate: ',cp%metab%A_rate,ATPs
!			write(*,'(a,i8,2e12.3)') 'insufficient A_rate: ',kcell,cp%metab%A_rate,ATPs
!			call CellDies(kcell)
			cp%state = DYING
			cp%dVdt = 0
			Ncells_dying(ityp) = Ncells_dying(ityp) + 1
			cycle
		endif			
	else
		if (cp%anoxia_tag) then
			if (tnow >= cp%t_anoxia_die) then
				call CellDies(kcell)
				Nanoxia_dead(ityp) = Nanoxia_dead(ityp) + 1
				cycle
			endif
		else
			if (anoxia_death .and. C_O2 < anoxia_threshold) then
				cp%t_anoxia = cp%t_anoxia + dt
				if (cp%t_anoxia > t_anoxia_limit) then
					cp%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
					if (.not.cp%anoxia_tag) then	! don't retag
						Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
					endif
					cp%anoxia_tag = .true.						! tagged to die later
				endif
			else
				cp%t_anoxia = 0
			endif
		endif
		
		call getGlucoseconc(cp,C_glucose)
		if (cp%aglucosia_tag) then
			if (tnow >= cp%t_aglucosia_die) then
				call CellDies(kcell)
				Naglucosia_dead(ityp) = Naglucosia_dead(ityp) + 1
				cycle
			endif
		else
			if (aglucosia_death .and. C_glucose < aglucosia_threshold) then
				cp%t_aglucosia = cp%t_aglucosia + dt
				if (cp%t_aglucosia > t_aglucosia_limit) then
					cp%t_aglucosia_die = tnow + aglucosia_death_delay	! time that the cell will die
					if (.not.cp%aglucosia_tag) then	! don't retag
						Naglucosia_tag(ityp) = Naglucosia_tag(ityp) + 1
					endif
					cp%aglucosia_tag = .true.						! tagged to die later
				endif
			else
				cp%t_aglucosia = 0
			endif
		endif
	endif
	
	do idrug = 1,ndrugs_used	
		ichemo = DRUG_A + 3*(idrug-1)
		if (.not.chemo(ichemo)%present) cycle
		if (cp%drug_tag(idrug)) cycle	! don't tag more than once
		dp => drug(idrug)
		kill_prob = 0
		death_prob = 0
		survival_prob = 1
		do im = 0,2
			if (.not.dp%kills(ityp,im)) cycle
			killmodel = dp%kill_model(ityp,im)		! could use %drugclass to separate kill modes
			Cdrug = cp%Cin(ichemo + im)
			Kd = dp%Kd(ityp,im)
			n_O2 = dp%n_O2(ityp,im)
			kmet = (1 - dp%C2(ityp,im) + dp%C2(ityp,im)*dp%KO2(ityp,im)**n_O2/(dp%KO2(ityp,im)**n_O2 + C_O2**n_O2))*dp%Kmet0(ityp,im)
			dMdt = kmet*Cdrug
			if (first) write(nflog,'(a,4e12.3)') 'Kmet0,C_O2,kmet,Cdrug: ',dp%Kmet0(ityp,im),C_O2,kmet,Cdrug
			call getDrugKillProb(killmodel,Kd,dMdt,Cdrug,dt,dkill_prob)
!			kill_prob = kill_prob + dkill_prob
			survival_prob = survival_prob*(1 - dkill_prob)
			death_prob = max(death_prob,dp%death_prob(ityp,im))
			first = .false.
		enddo
		kill_prob = 1 - survival_prob
	    if (par_uni(kpar) < kill_prob) then	
			cp%p_drug_death(idrug) = death_prob
			cp%drug_tag(idrug) = .true.
            Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) + 1
		endif
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getDrugKillProb(kill_model,Kd,dMdt,Cdrug,dt,dkill_prob)
integer :: kill_model
real(REAL_KIND) :: Kd, dMdt, Cdrug, dt, dkill_prob
real(REAL_KIND) :: SF, dtstep, kill_prob, c
integer :: Nsteps, istep

if (kill_model == 1) then
	c = Kd*dMdt
elseif (kill_model == 2) then
	c = Kd*dMdt*Cdrug
elseif (kill_model == 3) then
	c = Kd*dMdt**2
elseif (kill_model == 4) then
	c = Kd*Cdrug
elseif (kill_model == 5) then
	c = Kd*Cdrug**2
endif
SF = exp(-c*dt)
if (first) write(nflog,'(a,5e12.3)') 'Cdrug,dMdt,c,dt,SF: ',Cdrug,dMdt,c,dt,SF
dkill_prob = 1 - SF
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDies
integer :: kcell, i, kpar=0

do i = 1,10
	kcell = random_int(1,nlist,kpar)
	call CellDies(kcell)
enddo
stop
end subroutine

!-----------------------------------------------------------------------------------------
! A cell that dies must be subtracted from any counts it is on.
!-----------------------------------------------------------------------------------------
subroutine CellDies(kcell)
integer :: kcell
integer :: site(3), ityp, idrug
type(cell_type), pointer :: cp

cp => cell_list(kcell)
ityp = cp%celltype
if (cp%state == DYING) then
	Ncells_dying(ityp) = Ncells_dying(ityp) - 1
endif
Ncells = Ncells - 1
Ncells_type(ityp) = Ncells_type(ityp) - 1
if (cp%anoxia_tag) then
	Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
endif
if (cp%aglucosia_tag) then
	Naglucosia_tag(ityp) = Naglucosia_tag(ityp) - 1
endif
do idrug = 1,ndrugs_used
	if (cp%drug_tag(idrug)) then
		Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
	endif
enddo
if (cp%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
endif
cp%state = DEAD

ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(logmsg,'(a,i6,i6)') 'CellDies: ngaps > max_ngaps: ',ngaps,max_ngaps
    call logger(logmsg)
    stop
endif
gaplist(ngaps) = kcell
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AddToMedium(cp,site)
integer :: kcell, site(3)
integer :: ichemo
real(REAL_KIND) :: V, Cex(MAX_CHEMO), Cin(MAX_CHEMO)
type(cell_type),pointer :: cp
return

!Cex = occupancy(site(1),site(2),site(3))%C
Cin = cp%Cin
V = cp%V
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	chemo(ichemo)%medium_M = chemo(ichemo)%medium_M + V*Cin(ichemo) + (Vsite_cm3 - V)*Cex(ichemo)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine RemoveFromMedium
integer :: ichemo

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (chemo(ichemo)%constant) cycle
	chemo(ichemo)%medium_M = chemo(ichemo)%medium_M - Vsite_cm3*chemo(ichemo)%medium_Cbnd
enddo
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDivision(ok)
logical :: ok
integer :: kcell, kpar=0

kcell = random_int(1,nlist,kpar)
call divider(kcell,ok)
end subroutine

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  Division occurs when cell volume 
! exceeds the divide volume. 
! As the cell grows we need to adjust both Cin and Cex to maintain mass conservation.
! GROWTH DELAY
! When a cell has received a dose of radiation (or possibly drug - not yet considered)
! the cycle time is increased by an amount that depends on the dose.  The delay may be
! transmitted to progeny cells.
! 
! NOTE: now the medium concentrations are not affected by cell growth
!-----------------------------------------------------------------------------------------
subroutine new_grower(dt, changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: k, kcell, nlist0, ityp, idrug, prev_phase, kpar=0
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: rr(3), c(3), rad, d_desired, R, rrsum
integer, parameter :: MAX_DIVIDE_LIST = 100000
integer :: ndivide, divide_list(MAX_DIVIDE_LIST)
logical :: drugkilled
logical :: mitosis_entry, in_mitosis, divide, tagged

ok = .true.
changed = .false.
ccp => cc_parameters
nlist0 = nlist
ndivide = 0
!tnow = istep*DELTA_T !+ t_fmover
!if (colony_simulation) write(*,*) 'grower: ',nlist0,use_volume_method,tnow

do kcell = 1,nlist0
	kcell_now = kcell
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	divide = .false.
	mitosis_entry = .false.
	in_mitosis = .false.
	tagged = cp%anoxia_tag .or. cp%aglucosia_tag .or. (cp%state == DYING)
	if (tagged) then
		cp%dVdt = 0
	endif
	if (use_volume_method) then
!        if (colony_simulation) then
!            write(*,'(a,i6,L2,2e12.3)') 'kcell: ',kcell,cp%Iphase,cp%V,cp%divide_volume
!        endif
	    if (cp%Iphase) then
		    call growcell(cp,dt)
		    if (cp%V > cp%divide_volume) then	! time to enter mitosis
!    	        mitosis_entry = .true.
!				cp%Iphase = .false.
	            in_mitosis = .true.
!				cp%mitosis = 0
!				cp%t_start_mitosis = tnow
				
				cp%Iphase = .false.
				cp%mitosis = 0
				cp%t_start_mitosis = tnow
				ncells_mphase = ncells_mphase + 1
				
! The following are applicable when cell-cell forces are relevant
!				cp%nspheres = 2
!				call get_random_vector3(rr)	! set initial axis direction
!				cp%d = 0.1*small_d
!				c = cp%centre(:,1)
!				cp%centre(:,1) = c + (cp%d/2)*rr
!				cp%centre(:,2) = c - (cp%d/2)*rr
!				cp%d_divide = 2.0**(2./3)*cp%radius(1)
	        endif
	    else
	        in_mitosis = .true.
	    endif
	else
	    prev_phase = cp%phase
!	    if (cp%dVdt == 0) then
!			write(nflog,*) 'dVdt=0: kcell, phase: ',kcell,cp%phase 
!		endif
		if (cp%dVdt > 0) then
	        call timestep(cp, ccp, dt)
	    endif
        if (.not.cp%radiation_tag .and.(cp%NL2(1) > 0 .or. cp%NL2(2) > 0)) then	! irrepairable damage
			! For now, tag for death at mitosis
			cp%radiation_tag = .true.	! new_grower
		    Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
!			write(*,*) 'Tagged with NL2: ',istep,cp%NL2,Nradiation_tag(ityp)
		endif
        if (cp%phase >= M_phase) then
            if (prev_phase == Checkpoint2) then		! this is mitosis entry
!				write(*,*) 'mitosis_entry: kcell,istep: ',kcell,istep
                if (.not.cp%radiation_tag .and. cp%NL1 > 0) then		! lesions still exist, no time for repair
					cp%radiation_tag = .true.
				    Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
!				    write(*,*) 'Tagged with NL1: ',istep,cp%NL1,Nradiation_tag(ityp)
				endif
				
				cp%Iphase = .false.
				cp%mitosis = 0
				cp%t_start_mitosis = tnow
				ncells_mphase = ncells_mphase + 1
!				cp%nspheres = 2
!				call get_random_vector3(rr)	! set initial axis direction
!				rrsum = rrsum + rr
!				cp%d = 0.1*small_d
!				c = cp%centre(:,1)
!				cp%centre(:,1) = c + (cp%d/2)*rr
!				cp%centre(:,2) = c - (cp%d/2)*rr
!				cp%d_divide = 2.0**(2./3)*cp%radius(1)
            endif
            in_mitosis = .true.
        endif
		if (cp%phase < Checkpoint2 .and. cp%phase /= Checkpoint1) then
		    call growcell(cp,dt)
		endif	
	endif
!	if (mitosis_entry) then
	if (in_mitosis) then
		drugkilled = .false.
		do idrug = 1,ndrugs_used
			if (cp%drug_tag(idrug)) then
				call CellDies(kcell)
				changed = .true.
				Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
				drugkilled = .true.
				exit
			endif
		enddo
		if (drugkilled) cycle

		cp%mitosis = (tnow - cp%t_start_mitosis)/mitosis_duration
!		d_desired = max(cp%mitosis*cp%d_divide,small_d)
!		rr = cp%centre(:,2) - cp%centre(:,1)
!		cp%d = sqrt(dot_product(rr,rr))
!		c = (cp%centre(:,1) + cp%centre(:,2))/2
!		cp%site = c/DELTA_X + 1
!		call cubic_solver(d_desired,cp%V,rad)
!		cp%radius = rad
!		if (cp%d > 2*rad) then	! completion of mitosis - note that this overrides cp%phase
!			write(logmsg,'(a,i6,2e12.3,f7.3)') 'divides: d > 2*rad: ',kcell,cp%d,rad,cp%mitosis
!			call logger(logmsg)
!			divide = .true.
!		endif
			
		if (use_volume_method) then
			if (cp%growth_delay) then
				if (cp%G2_M) then
					if (tnow > cp%t_growth_delay_end) then
						cp%G2_M = .false.
					else
						cycle
					endif
				else
					cp%t_growth_delay_end = tnow + cp%dt_delay
					cp%G2_M = .true.
					cycle
				endif
			endif
			! try moving death prob test to here
			if (cp%radiation_tag) then
				R = par_uni(kpar)
				if (R < cp%p_rad_death) then
					call CellDies(kcell)
					changed = .true.
					Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
					cycle
				endif
			endif		
		else
		    ! Check for cell death by radiation lesions
!		    if (cp%NL1 > 0 .or. cp%NL2(2) > 0) then
			if (cp%radiation_tag) then
				call CellDies(kcell)
				changed = .true.
				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
				cycle
			endif		        
		endif
		
!		write(*,*) 'istep, kcell, mitosis: ',istep,kcell,cp%mitosis
        if (cp%mitosis >= 1) then
			divide = .true.
		endif
	endif
	if (divide) then
		ndivide = ndivide + 1
		if (ndivide > MAX_DIVIDE_LIST) then
		    write(logmsg,*) 'Error: growcells: MAX_DIVIDE_LIST exceeded: ',MAX_DIVIDE_LIST
		    call logger(logmsg)
		    ok = .false.
		    return
		endif
		divide_list(ndivide) = kcell
	endif
enddo
do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	call divider(kcell, ok)
	if (.not.ok) return
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  Division occurs when cell volume 
! exceeds the divide volume. 
! As the cell grows we need to adjust both Cin and Cex to maintain mass conservation.
! GROWTH DELAY
! When a cell has received a dose of radiation (or possibly drug - not yet considered)
! the cycle time is increased by an amount that depends on the dose.  The delay may be
! transmitted to progeny cells.
! 
! NOTE: now the medium concentrations are not affected by cell growth
!-----------------------------------------------------------------------------------------
subroutine grower(dt, changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: k, kcell, nlist0, ityp, idrug, prev_phase, kpar=0
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: R
integer, parameter :: MAX_DIVIDE_LIST = 10000
integer :: ndivide, divide_list(MAX_DIVIDE_LIST)
logical :: drugkilled
logical :: mitosis_entry, in_mitosis, divide

ok = .true.
changed = .false.
ccp => cc_parameters
nlist0 = nlist
ndivide = 0
!tnow = istep*DELTA_T !+ t_fmover
!if (colony_simulation) write(*,*) 'grower: ',nlist0,use_volume_method,tnow

do kcell = 1,nlist0
	kcell_now = kcell
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	ityp = cp%celltype
	divide = .false.
	mitosis_entry = .false.
	in_mitosis = .false.
	if (use_volume_method) then
!        if (colony_simulation) then
!            write(*,'(a,i6,L2,2e12.3)') 'kcell: ',kcell,cp%Iphase,cp%V,cp%divide_volume
!        endif
	    if (cp%Iphase) then
		    call growcell(cp,dt)
		    if (cp%V > cp%divide_volume) then	! time to enter mitosis
 !   	        mitosis_entry = .true.
				cp%Iphase = .false.
				in_mitosis = .true.
				cp%mitosis = 0
				cp%t_start_mitosis = tnow
	        endif
	    else
	        in_mitosis = .true.
	    endif
	else
	    prev_phase = cp%phase
!	    if (cp%dVdt == 0) then
!			write(nflog,*) 'dVdt=0: kcell, phase: ',kcell,cp%phase
!		endif
        call timestep(cp, ccp, dt)
!		if (istep > 1560 .and. kcell == 6) write(*,*) 'kcell: ',kcell,cp%phase
        if (cp%phase >= M_phase) then
            if (prev_phase == Checkpoint2) then
!                mitosis_entry = .true.
!				write(*,*) 'mitosis_entry: ',kcell,cp%drug_tag(1)
				cp%Iphase = .false.
				cp%mitosis = 0
				cp%t_start_mitosis = tnow
	            in_mitosis = .true.
!            else
!                in_mitosis = .true.
            endif
            in_mitosis = .true.
        endif
		if (cp%phase < Checkpoint2 .and. cp%phase /= Checkpoint1) then
		    call growcell(cp,dt)
		endif	
	endif
!	if (mitosis_entry) then
	if (in_mitosis) then
		drugkilled = .false.
!		write(*,*) 'mitosis_entry: ',kcell,cp%drug_tag(1)
		do idrug = 1,ndrugs_used
			if (cp%drug_tag(idrug)) then
				call CellDies(kcell)
				changed = .true.
				Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
				drugkilled = .true.
				exit
			endif
		enddo
		if (drugkilled) cycle
			
		if (use_volume_method) then
			if (cp%growth_delay) then
				if (cp%G2_M) then
					if (tnow > cp%t_growth_delay_end) then
						cp%G2_M = .false.
					else
						cycle
					endif
				else
					cp%t_growth_delay_end = tnow + cp%dt_delay
					cp%G2_M = .true.
					cycle
				endif
			endif
			! try moving death prob test to here
			if (cp%radiation_tag) then
				R = par_uni(kpar)
				if (R < cp%p_rad_death) then
					call CellDies(kcell)
					changed = .true.
					Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
					cycle
				endif
			endif		
		else
		    ! Check for cell death by radiation lesions
		    ! For simplicity: 
		    !   cell death occurs only at mitosis entry
		    !   remaining L1 lesions and L2c misrepair (non-reciprocal translocation) are treated the same way
		    !   L2a and L2b are treated as non-fatal
		    if (cp%NL1 > 0 .or. cp%NL2(2) > 0) then
				call CellDies(kcell)
				changed = .true.
				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
				cycle
			endif		        
		endif
		
!		cp%Iphase = .false.
!		cp%mitosis = 0
!		cp%t_start_mitosis = tnow
!	elseif (in_mitosis) then
		cp%mitosis = (tnow - cp%t_start_mitosis)/mitosis_duration
        if (cp%mitosis >= 1) then
			divide = .true.
		endif
	endif
	if (divide) then
		ndivide = ndivide + 1
		if (ndivide > MAX_DIVIDE_LIST) then
		    write(logmsg,*) 'Error: growcells: MAX_DIVIDE_LIST exceeded: ',MAX_DIVIDE_LIST
		    call logger(logmsg)
		    ok = .false.
		    return
		endif
		divide_list(ndivide) = kcell
	endif
enddo
do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	call divider(kcell, ok)
	if (.not.ok) return
enddo
end subroutine


!-----------------------------------------------------------------------------------------
! Need to check phase-dependence of growth
!-----------------------------------------------------------------------------------------
subroutine growcell(cp, dt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt
real(REAL_KIND) :: Cin_0(NCONST), Cex_0(NCONST)		! at some point NCONST -> MAX_CHEMO
real(REAL_KIND) :: dVdt,  Vin_0, dV,  metab_O2, metab_glucose, metab
integer :: ityp
logical :: oxygen_growth, glucose_growth, tagged

if (use_cell_cycle .and. .not.(cp%phase == G1_phase .or. cp%phase == S_phase .or. cp%phase == G2_phase)) then
	write(*,*) 'no growth - phase'
	return
endif
if (use_metabolism .and. cp%metab%A_rate < ATPg) then
	cp%dVdt = 0
!	write(*,*) 'A_rate < ATPg: ', kcell_now
	return
endif
ityp = cp%celltype
tagged = cp%anoxia_tag .or. cp%aglucosia_tag .or. (cp%state == DYING)
if (tagged) then
	cp%dVdt = 0
	write(*,*) 'tagged: ',kcell_now
	return
endif
if (colony_simulation) then
	if (use_metabolism) then
!		cp%metab%Itotal = cp%metab%Itotal + dt*cp%metab%I_rate
!		dVdt = max_growthrate(ityp)*cp%metab%I_rate/cp%metab%I_rate_max	! ***** Convert %I_rate to %dVdt *****
!		dVdt = cp%growth_rate_factor*dVdt
		metab = cp%metab%I_rate/cp%metab%I_rate_max
		dVdt = get_dVdt(cp,metab)
	else
		metab = 1
		dVdt = get_dVdt(cp,metab)
	endif
else
	if (use_metabolism) then
!		cp%metab%Itotal = cp%metab%Itotal + dt*cp%metab%I_rate
		! need to set cp%dVdt from cp%metab%I_rate
!		dVdt = max_growthrate(ityp)*cp%metab%I_rate/cp%metab%I_rate_max	! ***** Convert %I_rate to %dVdt *****
!		dVdt = cp%growth_rate_factor*cp%dVdt
		metab = cp%metab%I_rate/cp%metab%I_rate_max
		dVdt = get_dVdt(cp,metab)
	else
		oxygen_growth = chemo(OXYGEN)%controls_growth
		glucose_growth = chemo(GLUCOSE)%controls_growth
		Cin_0 = cp%Cin
		metab_O2 = O2_metab(Cin_0(OXYGEN))	! Note that currently growth depends only on O2
		metab_glucose = glucose_metab(Cin_0(GLUCOSE))   
		if (oxygen_growth .and. glucose_growth) then
			metab = metab_O2*metab_glucose
		elseif (oxygen_growth) then
			metab = metab_O2
		elseif (glucose_growth) then
			metab = metab_glucose
		endif
		dVdt = get_dVdt(cp,metab)
		if (suppress_growth) then	! for checking solvers
			dVdt = 0
		endif
	endif
endif
cp%dVdt = dVdt
Vin_0 = cp%V
dV = dVdt*dt
if (use_cell_cycle .and. .not.(cp%phase == G1_phase .or. cp%phase == S_phase .or. cp%phase == G2_phase)) then
    dV = 0
endif
cp%V = Vin_0 + dV
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetInitialGrowthRate
integer :: kcell, ityp
real(REAL_KIND) :: C_O2, C_glucose, metab, metab_O2, metab_glucose, dVdt
logical :: oxygen_growth, glucose_growth
type(cell_type), pointer :: cp

oxygen_growth = chemo(OXYGEN)%controls_growth
glucose_growth = chemo(GLUCOSE)%controls_growth
do kcell = 1,nlist
    if (colony_simulation) then
        cp => ccell_list(kcell)
    else
        cp => cell_list(kcell)
    endif
	if (cp%state == DEAD) cycle
	metab = 1
	dVdt = get_dVdt(cp,metab)
!	C_O2 = chemo(OXYGEN)%bdry_conc
!	C_glucose = chemo(GLUCOSE)%bdry_conc
!	metab_O2 = 1
!	metab_glucose = 1
!	if (oxygen_growth .and. glucose_growth) then
!	    metab_O2 = O2_metab(C_O2)
!		metab_glucose = glucose_metab(C_glucose)
!		metab = metab_O2*metab_glucose
!	elseif (oxygen_growth) then
!	    metab_O2 = O2_metab(C_O2)
!		metab = metab_O2
!	elseif (glucose_growth) then
!		metab_glucose = glucose_metab(C_glucose)
!		metab = metab_glucose
!	endif
	if (suppress_growth) then	! for checking solvers
		dVdt = 0
	endif
	cp%dVdt = dVdt
!    if (cp%dVdt == 0) then
!        write(*,*) 'setinitialgrowthrate: dVdt: = 0: kcell: ',kcell,metab_O2,metab_glucose
!        stop
!    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Need to account for time spent in mitosis.  Because there is no growth during mitosis,
! the effective divide_time must be reduced by mitosis_duration.
! Note that TERMINAL_MITOSIS is the only option.
!-----------------------------------------------------------------------------------------
function get_dVdt(cp, metab) result(dVdt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: metab, dVdt
integer :: ityp
real(REAL_KIND) :: r_mean, c_rate

ityp = cp%celltype
if (use_cell_cycle) then
    if (use_constant_growthrate) then
        dVdt = max_growthrate(ityp)
    else
        dVdt = metab*max_growthrate(ityp)
    endif
else
    if (use_V_dependence) then
	    if (use_constant_divide_volume) then
		    dVdt = metab*log(2.0)*cp%V/(cp%divide_time - mitosis_duration)
	    else
		    c_rate = log(2.0)/(divide_time_mean(ityp) - mitosis_duration)
		    dVdt = c_rate*cp%V*metab
	    endif
    else
	    if (use_constant_divide_volume) then
		    dVdt = 0.5*metab*Vdivide0/(cp%divide_time  - mitosis_duration)
	    else
		    ityp = cp%celltype
		    r_mean = 0.5*Vdivide0/(divide_time_mean(ityp) - mitosis_duration)
		    dVdt = r_mean*metab
	    endif
    endif
endif
!if (cp%dVdt == 0) then
!    write(nflog,*) 'get_dVdt: = 0'
!    write(*,*) 'get_dVdt: dVdt = 0'
!    stop
!endif
end function


!-----------------------------------------------------------------------------------------
! New version
!-----------------------------------------------------------------------------------------
subroutine divider(kcell1, ok)
integer :: kcell1
logical :: ok
integer :: kcell2, ityp, nbrs0, im, imax, ipdd
real(REAL_KIND) :: r(3), c(3), cfse0, cfse2, V0, Tdiv, gfactor
type(cell_type), pointer :: cp1, cp2
type(cycle_parameters_type), pointer :: ccp

!write(*,*) 'divider:'
!write(logmsg,*) 'divider: ',kcell1 
!call logger(logmsg)
ok = .true.
ndivided = ndivided + 1
ccp => cc_parameters
if (colony_simulation) then
    cp1 => ccell_list(kcell1)
else
	cp1 => cell_list(kcell1)
endif
if (ngaps > 0) then
    kcell2 = gaplist(ngaps)
    ngaps = ngaps - 1
else
	nlist = nlist + 1
	if (nlist > MAX_NLIST) then
		ok = .false.
		call logger('Error: Maximum number of cells MAX_NLIST has been exceeded.  Increase and rebuild.')
		return
	endif
	kcell2 = nlist
endif
ncells = ncells + 1
ityp = cp1%celltype
ncells_type(ityp) = ncells_type(ityp) + 1
ncells_mphase = ncells_mphase - 1
if (colony_simulation) then
    cp2 => ccell_list(kcell2)
else
	cp2 => cell_list(kcell2)
endif

cp1%state = ALIVE
cp1%generation = cp1%generation + 1
V0 = cp1%V/2
cp1%V = V0
cp1%birthtime = tnow
cp1%divide_volume = get_divide_volume(ityp,V0,Tdiv, gfactor)
cp1%divide_time = Tdiv
cp1%fg = gfactor
!if (use_metabolism) then	! Fraction of I needed to divide = fraction of volume needed to divide NOT USED
!	cp1%metab%I2Divide = get_I2Divide(cp1)
!	cp1%metab%Itotal = 0
!endif
cp1%mitosis = 0
cfse0 = cp1%CFSE
cp1%CFSE = generate_CFSE(cfse0/2)
cfse2 = cfse0 - cp1%CFSE
! Halve drug levels - OK for normal drugs?
do ipdd = 0,1
	if (chemo(DRUG_A + 3*ipdd)%used) then
		if (drug(ipdd+1)%phase_dependent) then
			imax = 1
		else
			imax = 2
		endif
		do im = 0,imax
			cp1%Cin(DRUG_A + 3*ipdd + im) = cp1%Cin(DRUG_A + 3*ipdd + im)/2
		enddo
	endif
enddo

cp1%drug_tag = .false.
cp1%anoxia_tag = .false.
cp1%t_anoxia = 0
cp1%aglucosia_tag = .false.
cp1%t_aglucosia = 0

if (cp1%growth_delay) then
	cp1%N_delayed_cycles_left = cp1%N_delayed_cycles_left - 1
	cp1%growth_delay = (cp1%N_delayed_cycles_left > 0)
endif
cp1%G2_M = .false.
!if (use_metabolism) then
!    cp1%G1_time = tnow + (cp1%metab%I_rate_max/cp1%metab%I_rate)*cp1%fg*ccp%T_G1(ityp)
!endif
cp1%G1_time = tnow + (max_growthrate(ityp)/cp1%dVdt)*cp1%fg*ccp%T_G1(ityp)    ! time spend in G1 varies inversely with dV/dt
cp1%Iphase = .true.
cp1%phase = G1_phase

!if (max_growthrate(ityp)/cp1%dVdt > 2) then
!    write(*,*) 'dVdt: ',kcell1,max_growthrate(ityp)/cp1%dVdt
!endif

ndoublings = ndoublings + 1
doubling_time_sum = doubling_time_sum + tnow - cp1%t_divide_last
if (kcell1 < 100) then
	write(nflog,'(a,2i6,3e12.3)') 'kcell1, ndoublings, doubling_time_sum, tnow, cp1%t_divide_last: ',kcell1,ndoublings,doubling_time_sum, tnow, cp1%t_divide_last
endif
cp1%t_divide_last = tnow 

! Second cell
cp2 = cp1

! These are the variations from cp1
cp2%divide_volume = get_divide_volume(ityp,V0,Tdiv, gfactor)
cp2%divide_time = Tdiv
cp2%fg = gfactor
if (use_metabolism) then	! Fraction of I needed to divide = fraction of volume needed to divide
    cp2%G1_time = tnow + (cp2%metab%I_rate_max/cp2%metab%I_rate)*cp2%fg*ccp%T_G1(ityp)
	cp2%growth_rate_factor = get_growth_rate_factor()
	cp2%ATP_rate_factor = get_ATP_rate_factor()
else
	cp2%G1_time = tnow + (max_growthrate(ityp)/cp2%dVdt)*cp2%fg*ccp%T_G1(ityp)    ! time spend in G1 varies inversely with dV/dt
endif
cp2%CFSE = cfse2
if (cp2%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
endif
!if (colony_simulation) write(*,'(a,i6,2e12.3)') 'new cell: ',kcell2,cp2%V,cp2%divide_volume
end subroutine

#if 0
!----------------------------------------------------------------------------------
! Makes a slight modification to the Michaelis-Menten function to create a
! "soft landing" as C -> 0
! f1(C) = (C-d)/(C0 + C-d)  is a shifted version of the MM curve
! f0(C) = kC^2             is a function with derivative -> 0 as C -> 0
! At C=e, we want f0(e) = f1(e), and f0'(e) = f1'(e)
! =>
! ke^2 = (e-d)/(C0 + e-d)
! 2ke = C0/(Co + e-d)^2
! => e/2 = (e-d)(C0 + e-d)/C0
! Set x = e-d
! C0(x+d) = 2x(C0+x)
! => x = e-d = (sqrt(Co^2 + 8dC0) - C0)/4
! k = (e-d)/(e^2(C0 + e-d))
! We fix d (as small as possible, by trial and error) then deduce e, k.
! Note: This has really been superceded by the option of a Hill function with N = 2.
!----------------------------------------------------------------------------------
subroutine AdjustMM
real(REAL_KIND) :: deltaC, C0, C1

C0 = chemo(OXYGEN)%MM_C0
if (MM_THRESHOLD > 0) then
	deltaC = MM_THRESHOLD
	C1 = deltaC + (sqrt(C0*C0 + 8*C0*deltaC) - C0)/4
	ODEdiff%k_soft = (C1-deltaC)/(C1*C1*(C0+C1-deltaC))
	ODEdiff%C1_soft = C1
	ODEdiff%deltaC_soft = deltaC
else
	ODEdiff%k_soft = 0
	ODEdiff%C1_soft = 0
	ODEdiff%deltaC_soft = 0
endif
!write(logmsg,'(a,4e12.4)') 'AdjustMM: C0, deltaC, C1, k: ',C0, ODEdiff%deltaC_soft, ODEdiff%C1_soft, ODEdiff%k_soft
!call logger(logmsg)
end subroutine
#endif

!----------------------------------------------------------------------------------
! This is used to adjust a cell's growth rate to introduce some random variability
! into divide times.
!----------------------------------------------------------------------------------
function get_growth_rate_factor() result(r)
real(REAL_KIND) :: r
real(REAL_KIND) :: dr = 0.2
integer :: kpar = 0

r = 1 + (par_uni(kpar) - 0.5)*dr
end function

!----------------------------------------------------------------------------------
! This is used to adjust a cell's ATP rate to introduce some random variability
! into transitions to no-growth and death.
!----------------------------------------------------------------------------------
function get_ATP_rate_factor() result(r)
real(REAL_KIND) :: r
real(REAL_KIND) :: dr = 0.2		! was 0.2
integer :: kpar = 0

r = 1 + (par_uni(kpar) - 0.5)*dr
end function


end module
