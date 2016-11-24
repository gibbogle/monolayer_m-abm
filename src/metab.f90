! To test simple cell metabolism model
! Concentration of ATP varies by < 10%  https://en.wikipedia.org/wiki/Glycolysis#Intermediates_for_other_pathways

! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM

! Question: How do the results of this model translate into cell rate of volume growth?
!----------------------------------------------------------------------------------------------------------------
module metabolism
use real_kind_mod
use global
use chemokine
implicit none

! From spheroid-abm, the max rates of consumption of oxygen and glucose are:
!   oxygen:  6.25e-17 moles/cell/s
!   glucose: 3.80e-17 moles/cell/s
! We work with mumol/cell/sec, and convert these values by scaling by 1.0e6, to give
!   oxygen:  6.25e-11 mumol/cell/s
!   glucose: 3.80e-11 mumol/cell/s

real(REAL_KIND) :: Hill_Km_O2
real(REAL_KIND) :: Hill_N_O2
real(REAL_KIND) :: Hill_Km_G     ! Hill Km for dependence of glycolysis rate on glucose
real(REAL_KIND) :: Hill_N_G      ! Hill N for dependence of glycolysis rate on glucose
real(REAL_KIND) :: Hill_Km_P     ! Hill Km for dependence of pyruvate oxidation rate on pyruvate
real(REAL_KIND) :: Hill_N_P      ! Hill N for dependence of pyruvate oxidation rate on pyruvate
real(REAL_KIND) :: K_H1(MAX_CELLTYPES)     ! HIF-1 k1
real(REAL_KIND) :: K_H2(MAX_CELLTYPES)     ! HIF-1 k2
real(REAL_KIND) :: K_Ha     ! HIF-1 ka
real(REAL_KIND) :: K_Hb(MAX_CELLTYPES)     ! HIF-1 kb
real(REAL_KIND) :: K_PDK(MAX_CELLTYPES)    ! K_PDK
real(REAL_KIND) :: PDKmin(MAX_CELLTYPES)    ! PDKmin
real(REAL_KIND) :: P_conc_min
!real(REAL_KIND) :: Vcell     ! cell volume
real(REAL_KIND) :: K_PL(MAX_CELLTYPES)     ! P -> L
real(REAL_KIND) :: K_LP(MAX_CELLTYPES)     ! L -> P

real(REAL_KIND) :: f_G_norm, f_P_norm, r_P_norm, r_G_norm, r_A_norm, r_I_norm

real(REAL_KIND) :: A_rate_base(MAX_CELLTYPES)	! total rate of production of ATP under full nutrition

!real(REAL_KIND) :: C_O, C_G, C_P, C_L, H, M_I

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_glycosis_rate(ityp, H, G) result(rate)
integer :: ityp
real(REAL_KIND) :: H, G, rate

!write(*,*) 'Hill_Km_G,Hill_N_G: ',Hill_Km_G,Hill_N_G
rate = (1 + K_Hb(ityp)*H)*G**Hill_N_G/(G**Hill_N_G + Hill_Km_G**Hill_N_G)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_Poxidation_rate(ityp, PDK, P) result(rate)
integer :: ityp
real(REAL_KIND) :: PDK, P, rate

!write(*,*) 'Hill_Km_G,Hill_N_G: ',Hill_Km_G,Hill_N_G
rate = PDK*P**Hill_N_P/(P**Hill_N_P + Hill_Km_P**Hill_N_P)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine SetupMetabolism
real(REAL_KIND) :: f_PO, f_PA
integer :: ityp

Hill_Km_O2 = chemo(OXYGEN)%MM_C0
Hill_N_O2 = chemo(OXYGEN)%Hill_N
Hill_Km_G = chemo(GLUCOSE)%MM_C0
Hill_N_G = chemo(GLUCOSE)%Hill_N
Hill_Km_P = Hill_Km_G/10	! just because as yet we have no idea 
Hill_N_P = Hill_N_G
N_PO = 3
K_Ha = chemo(GLUCOSE)%max_cell_rate	! mumol/cell/s = 3.80e-17*1.0e6

do ityp = 1,2
	f_G_norm = N_GI(ityp)
	f_P_norm = N_PI(ityp)
	f_PO = N_PO(ityp)
	f_PA = N_PA(ityp)
	r_P_norm = f_MM(0.18d0,Hill_Km_O2,int(Hill_N_O2))*chemo(OXYGEN)%max_cell_rate/f_PO
	r_G_norm = K_Ha*get_glycosis_rate(ityp,0.0d0,5.5d0)
	r_A_norm = 2*(1-f_G_norm)*r_G_norm + f_PA*(1-f_P_norm)*r_P_norm
	r_I_norm = f_G_norm*r_G_norm + f_P_norm*r_P_norm
	ATPg(ityp) = f_ATPg(ityp)*r_A_norm
	ATPs(ityp) = f_ATPs(ityp)*r_A_norm
	metabolic(ityp)%I_rate_max = f_G_norm*r_G_norm + f_P_norm*r_P_norm
	write(*,'(a,4e12.3)') 'f_G_norm,r_G_norm,f_P_norm,r_P_norm: ',f_G_norm,r_G_norm,f_P_norm,r_P_norm
enddo
metabolic%HIF1 = 0
metabolic%PDK1 = 1
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_I2Divide(cp) result(I2div)
type(cell_type), pointer :: cp
real(REAL_KIND) :: I2div
integer :: ityp

ityp = cp%celltype
I2div = cp%divide_time*metabolic(ityp)%I_rate_max
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_HIF1steadystate(ityp,C_O) result(H)
integer :: ityp
real(REAL_KIND) :: C_O, H

H = exp(-K_H1(ityp)*C_O)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine oldSetHIF1(ityp, C_O, H, dt_big)
integer :: ityp
real(REAL_KIND) :: C_O, H, dt_big
integer :: it, ntlim = 10, nt
real(REAL_KIND) :: tlim, dt, dHdt

if (C_O > 0.01) then
    H = get_HIF1steadystate(ityp,C_O)
else
    dHdt = K_H2(ityp)*(1 - H*exp(K_H1(ityp)*C_O))
    if (dHdt > 0) then
		tlim = (1 - H)/dHdt
	else
		tlim = H/dHdt
	endif
	dt = min(dt_big, tlim/ntlim)
	nt = max(ntlim,int(dt_big/dt))
	dt = dt_big/nt
	do it = 1,nt
        dHdt = K_H2(ityp)*(1 - H*exp(K_H1(ityp)*C_O))
        H = H + dHdt*dt
        H = max(H,0.0)
        H = min(H,1.0)
	enddo
endif
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine analyticSetHIF1(ityp, C_O, H, dt)
integer :: ityp
real(REAL_KIND) :: C_O, H, dt
real(REAL_KIND) :: a, b, c, H0

a = K_H2(ityp)
b = exp(K_H1(ityp)*C_O)
H0 = H
c = 1 - b*H0
H = (1 - c*exp(-a*b*dt))/b
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine analyticSetPDK1(ityp, H, P, dt)
integer :: ityp
real(REAL_KIND) :: P, H, dt
real(REAL_KIND) :: a, b, c, d, P0

a = K_PDK(ityp)
d = 1 - PDKmin(ityp)
b = (1 - d*H)
P0 = P
c = P0 - b
P = b + c*exp(-a*dt)
end subroutine

!--------------------------------------------------------------------------
! The fraction of pyruvate that is converted to acetyl-CoA depends on the
! rate of glycolysis.  The normalised rate has a maximum value of 1 under
! normoxic conditions. This corresponds to the minimum pyruvate oxidation fraction
! (equivalently, the maximum lactate fraction).  The minimum pyruvate oxidation
! fraction is the fraction of pyruvate that is directed to the Krebs cycle
! when both glucose and oxygen are in ample supply.  
! In fact this is the nominal minimum, and the fraction can be even less 
! if hypoxia elevates glycosis rate.
! The basic idea is that over a range of glycolysis rate, the total rate of
! production of ATP is a constant.
! The total is the sum of ATP produced by glycolysis and by pyruvate oxidation.
! If the glycolysis rate is reduced, the intermediate production rate can be maintained 
! by increasing the pyruvate oxidation fraction (reducing lactate production), to the limit
! of fraction = 1.  Further reductions in glycolysis will then reduce ATP production.
! REVISED
! Glycolysis:
! A fraction N_GI goes to make intermediates, I_rate = N_GI*G_rate
! the remainder makes pyruvate, PP_rate = 2*(1 - N_GI)*G_rate
! and A_rate = PP_rate
! Pyruvate oxidation:
! A fraction N_PI goes to make intermediates via the TCA cycle, I_rate = N_PI*PO_rate
! the remainder (1 - N_PI) is fully oxidised, producing N_PA ATP/mole, (N_PA = 18)
! and A_rate from pyruvate oxidation = N_PA*(1 - N_PI)*PO_rate
! Bill:
! The intermediates used for anabolism are downstream of acetyl-CoA (including acetyl-CoA itself).
! Yes, PDK1 reduces pyruvate utilisation (r) but this is upstream of acetyl-CoA. So the formalism 
! needs to have the intermediates coming from acetyl-CoA rather than pyruvate itself.
!--------------------------------------------------------------------------
!subroutine get_metab_rates(ityp, H, PDK, Cin, G_rate, PP_rate, PO_rate)
subroutine get_metab_rates(ityp, M, Cin)
integer :: ityp
type(metabolism_type), pointer :: M
real(REAL_KIND) :: Cin(:)
real(REAL_KIND) :: G_norm       ! normalised glycolysis rate
real(REAL_KIND) :: C_P			! steady-state pyruvate concentration
logical :: use_new_metab = .true.

if (Cin(GLUCOSE) == 0) then
	write(*,*) 'Glucose concentration = 0'
	stop
endif

call new_metab2(ityp,M,Cin(OXYGEN),Cin(GLUCOSE),Cin(LACTATE))
end subroutine

!--------------------------------------------------------------------------
! This is used in the monolayer model, where all cells have the same IC
! concentrations, and therefore the same metabolism.
! We may not need G_rate, PP_rate, PO_rate
! We update only by cell type here, because in the monolayer all cells of the
! same type have the same metabolic state.
!--------------------------------------------------------------------------
subroutine update_metabolism
type(metabolism_type), pointer :: M
integer :: kcell, ityp
real(REAL_KIND) :: Itotal, I2Divide
type(cell_type), pointer :: cp

do ityp = 1,Ncelltypes
	M => metabolic(ityp)
!	call get_metab_rates(ityp, M%HIF1, M%PDK1, Caverage(1:MAX_CHEMO), M%G_rate, M%PP_rate, M%PO_rate)
	call get_metab_rates(ityp, M, Caverage(1:MAX_CHEMO))
! Revised
!	M%L_rate = M%PP_rate - M%PO_rate						! rate of production of lactate
!	M%A_rate = M%PP_rate + N_PA(ityp)*(1 - N_PI(ityp))*M%PO_rate
!	M%I_rate = N_GI(ityp)*M%G_rate + N_PI(ityp)*M%PO_rate	! rate of production of intermediates
!	M%O_rate = N_PO(ityp)*(1 - N_PI(ityp))*M%PO_rate		! rate of consumption of oxygen
enddo

do kcell = 1,nlist
    if (colony_simulation) then
        cp => ccell_list(kcell)
    else
        cp => cell_list(kcell)
    endif
    if (cp%state == DEAD) cycle
    Itotal = cp%metab%Itotal	! The only fields we don't want to change are Itotal and I2Divide, which are cell-specific
    I2Divide = cp%metab%I2Divide
    cp%metab = metabolic(cp%celltype)
    cp%metab%Itotal = Itotal
    cp%metab%I2Divide = I2Divide
enddo

end subroutine

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
! Use the "soft landing" option for Hill_N = 1 if MM_threshold = 0
! TO BE REMOVED - NEEDED for metab only
!----------------------------------------------------------------------------------
!real(REAL_KIND) function O2_metab(C)
!implicit none
!integer :: ichemo
!real(REAL_KIND) :: C
!
!ichemo = OXYGEN
!if (C > 0) then
!	O2_metab = C/(chemo(ichemo)%MM_C0 + C)
!else
!	O2_metab = 0
!endif
!end function

!--------------------------------------------------------------------------
! Use:
!	r_G for dG/dt, the rate of glycolysis = rate of glucose consumption
!	f_G for N_GI, the fraction of dG/dt that goes to intermediates
!	r_P for dPO/dt, rate of oxidation of pyruvate
!	f_P for N_PI, the fraction of dPO/dt that goes to intermediates
!	r_A for dA/dt, rate of ATP production
!	r_I for dI/dt, rate of intermediates production
!
!	r_G_norm for dG/dt under normal conditions, no nutrient constraints, H = 1
!	f_G_norm for the f_G under normal conditions (upper bound of f_G)
!	r_P_norm for dPO/dt under normal conditions
!	f_P_norm for f_P under normal conditions (upper bound of f_P)
!	r_A_norm for r_A under normal conditions
!	r_I_norm for r_I under normal conditions
!	alpha = r_P as a fraction of dP/dt under normal conditions
!
!	P_conc for IC pyruvate concentration
!	L_conc for IC lactate concentration
!	P_conc_min for minimumvalue of P_conc
!	
! r_P = 2*(1 - f_G)*r_G + V*(K2*C_L - K1*C_P - dC_P/dt) = fPDK*r_P_max*MM(O2)*MM(C_P)
! with the constraint that C_P >= 0
! Steady-state approach may not be feasible, because if there is 
! plenty of glucose but O2 is very low, r_G will be high but r_P will tend
! towards 0.  This must lead to an increase in C_P.
!--------------------------------------------------------------------------
subroutine new_metab2(ityp, M, C_O2, C_G, C_L)
integer :: ityp
type(metabolism_type), pointer :: M
real(REAL_KIND) :: C_O2, C_G, C_L
real(REAL_KIND) :: r_G, fPDK
real(REAL_KIND) :: f_G, f_P, r_P, r_A, r_I, r_L, f_PO, f_PA
real(REAL_KIND) :: K1, K2, C_P
real(REAL_KIND) :: r_GP, r_GA, r_PA, r_P_fac, V, Km_O2, Km_P, a, b, c, e, MM_P, r_A_target
real(REAL_KIND) :: base_O_rate = 2.0e-11
integer :: N_O2, N_P, it
logical :: use_MM_P = .true.

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_P = 1
Km_P = Hill_Km_P		! not used if .not.use_MM_P
V = Vcell_cm3

M%G_rate = K_Ha*get_glycosis_rate(ityp,M%HIF1,C_G)
r_G = M%G_rate
if (r_G == 0) then
	write(*,'(a,3e12.3)') 'r_G=0: K_Ha,M%HIF1,C_G: ',K_Ha,M%HIF1,C_G
	stop
endif
fPDK = M%PDK1

f_G_norm = N_GI(ityp)
f_P_norm = N_PI(ityp)
f_PO = N_PO(ityp)
f_PA = N_PA(ityp)
K1 = K_PL(ityp)
K2 = K_LP(ityp)
f_G = f_G_norm
f_P = f_P_norm
r_GP = 2*(1-f_G)*r_G	! rate of production of pyruvate by glycolysis, depends on f_G
r_GA = r_GP				! rate of production of ATP by glycolysis

!write(*,'(a,5e12.3)') 'f_G,r_GA: ',f_G,r_GA
r_P_fac = fPDK*f_MM(C_O2,Km_O2,N_O2)*chemo(OXYGEN)%max_cell_rate/f_PO
!write(*,'(a,2e12.3)') 'C_O2,f_MM(C_O2,Km_O2,N_O2): ',C_O2,f_MM(C_O2,Km_O2,N_O2)
!write(*,'(a,5e12.3)') 'r_P_norm,f_MM(C_O2,Km_O2,N_O2),r_Ptemp: ',r_P_norm,f_MM(C_O2,Km_O2,N_O2),r_Ptemp
! (1-f_P)*(f_PA*r_Ptemp) = (r_A_norm - r_GA) = r_PA

MM_P = 1
r_A_target = r_A_norm

do it = 1,2
	r_PA = (r_A_target - r_GA)
	f_P = 1 - r_PA/(f_PA*r_P_fac*MM_P)
	if (f_P < 0) then
		f_P = 0
		r_A_target = r_GA + f_PA*r_P_fac*MM_P
	endif
	r_PA = (1-f_P)*f_PA*r_P_fac*MM_P
	r_GA = r_A_target - r_PA
	if (r_GA < 0) then
		write(*,'(a,i2,5e12.3)') 'r_GA < 0: ',it,f_P,r_G,r_A_target,r_PA,r_GA
		r_GA = 0
		r_A_target = r_PA
!		stop
	endif
	f_G = 1 - r_GA/(2*r_G)
	if (f_G < 0) then
		write(*,*) 'f_G < 0: ',f_G
		f_G = 0
	endif
	r_GP = 2*(1-f_G)*r_G
	! Now from f_G, solve for steady-state C_P and r_P
	! Note: this assumes that N_P = 1
	if (use_MM_P) then
		e = r_GP + V*K2*C_L
		a = V*K1
		b = r_P_fac - e + V*K1*Km_P
		c = -e*Km_P
		C_P = solve_C_P(a,b,c)
		!	write(*,'(a,4e12.3)') 'r_G,f_G,r_GP,C_L: ',r_G,f_G,r_GP,C_L
		!	write(*,'(a,5e12.3)') 'a,b,c,e,C_P: ',a,b,c,e,C_P
		MM_P = f_MM(C_P,Km_P,N_P)
		r_P = r_P_fac*MM_P	! can check this against r_GP + V*(K2*C_L - K1*C_P)
		!write(*,'(a,2e12.3)') 'r_P check: ',r_P,r_GP + V*(K2*C_L - K1*C_P)
		!if (abs(r_P-(r_GP + V*(K2*C_L - K1*C_P))) > 1.0e-6*r_P) then
		!	write(*,*) 'r_P check failed'
		!	stop
		!endif
		!write(*,'(a,3e12.3)') 'r_P: ',r_P_fac,MM_P,r_P
	else
		! Simply let r_P = r_P_fac = r_GP + V*(K2*C_L - K1*C_P)
		r_P = r_P_fac
		C_P = (K2*C_L + (r_GP - r_P)/V)/K1
		if (C_P < 0) then
			write(*,'(a,3e12.3)') 'C_P < 0: C_O2, C_P, C_L: ', C_O2,C_P,C_L
			stop
!			C_P = 0
!			r_P = r_GP + V*K2*C_L
		endif
	endif
	r_GA = r_GP
	r_PA = f_PA*(1-f_P)*r_P
	r_A = r_GA + r_PA
	r_I = f_G*r_G + f_P*r_P
	r_L = V*(K1*C_P - K2*C_L)
	if (f_P == 0) exit
enddo
!write(*,'(a,4e12.3)') 'C_L,C_P,Km_P,MM_P: ',C_L,C_P,Km_P,MM_P
!write(*,'(a,e12.3,2f8.3)') 'r_G,f_G,f_P: ',r_G,f_G,f_P

if (r_I < 0) stop
M%A_rate = r_A				! production
M%I_rate = r_I				! production
M%PO_rate = r_P
M%O_rate = f_PO*r_P			! consumption
M%L_rate = r_L				! production
M%GA_rate = r_GA			! production

!if (M%O_rate < 0.37e-10) write(*,'(a,5e12.3)') 'r_P_fac,MM_P,C_P: ',r_P_fac,MM_P,C_P,r_P,f_PO
! Add base rate correction
M%O_rate = M%O_rate + base_O_rate
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function solve_C_P(a,b,c) result(x)
real(REAL_KIND) :: a, b, c, x
real(REAL_KIND) :: d

d = b*b - 4*a*c
if (d < 0) then
	write(*,*) 'Error: solve_C_P: a,b,c,d: ',a,b,c,d
!	write(*,'(a,3e12.3)') 'a,b,c: ',a,b,c
!	write(*,'(a,e12.3)') '-b/2a: ',-b/(2*a)
	x = 0
	return
else
	d = sqrt(d)
endif
x = (-b + d)/(2*a)
if (x < 0) then
	write(*,*) 'solve_C_P: x < 0: ',x
	stop
endif
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function f_MM(C,Km,N) result(v)
real(REAL_KIND) :: C, Km, v
integer :: N

v = C**N/(Km**N + C**N)
end function

end module

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!program main
!use metabolism
!integer :: nt = 10
!real(REAL_KIND) :: Hss, dGdt, dOdt
!
!open(nfout,file='metab.out',status='replace')
!call setup
!C_O = 0.01
!C_G = 3.0
!C_L = 1.0
!C_P = C_L*K_LP/K_PL
!M_I = 0
!
!Hss = get_HIF1steadystate(C_O)
!write(*,'(a,2e12.3)') 'steady-state HIF-1: ',C_O,Hss
!H = 1.0*Hss
!H = min(H, 1.0)
!H = max(H, 0.0)
!dGdt = -0.00004    ! test rate of change of ambient glucose
!dOdt = -0.0000004    ! test rate of change of ambient O2
!write(nfout,'(a)') '   hour C_G C_O C_L C_P H G_rate PP_rate PO_rate L_rate A_rate I_rate O_rate'
!do istep = 1,60*24
!    C_G = C_G + dGdt*dt
!    C_G = max(C_G,0.0)
!    C_G = min(C_G,5.5)
!!    C_O = C_O + dOdt*dt
!!    C_O = max(C_O,0.0)
!!    C_O = min(C_O,0.10)
!    call timestep(nt)
!!    write(*,'(a,e12.3)') 'M_I: ', M_I
!enddo
!!write(*,*) 'Hss: ',Hss
!end program
