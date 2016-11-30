module optimise
use metabolism
implicit none

contains

!--------------------------------------------------------------------------
! Given:
! r_G = glycolysis rate
! C_L = IC lactate concentration
! r_Pm = max pyruvate oxidation rate accounting for O2 and fPDK, but not f_G, C_P
! solves for f_G, f_P and hence C_P and all the rates.
! Note: f_G always = its maximum value.
! This suggests a simpler solution method that treats f_G as fixed.
!
! NOTE:
! It's OK to change this file, but metab.f90 should not be altered except to
! comment out 'use chemokine'.  The code should remain identical to that
! used by monolayer_m.
! 
! For now assume that f_G is always at its maximum, = f_G_norm.
!--------------------------------------------------------------------------
subroutine optimiser(ityp, r_G, C_L, r_Pm, M, x, y, C_P)
integer :: ityp
type(metabolism_type), pointer :: M
real(REAL_KIND) :: r_G, r_Pm, C_L, x, y, C_P
real(REAL_KIND) :: f_G, f_P, r_P, r_A, r_I, r_L, f_PO, f_PA
real(REAL_KIND) :: K1, K2
real(REAL_KIND) :: r_GP, r_GA, r_PA, V, Km_O2, Km_P, a, b, c, e, MM_P, r_A_target
integer :: N_O2, N_P, it
logical :: use_MM_P = .true.

N_P = 1
Km_P = Hill_Km_P(ityp)		! not used if .not.use_MM_P
N_O2 = chemo(OXYGEN)%Hill_N
Km_O2 = chemo(OXYGEN)%MM_C0
V = Vcell_cm3

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

MM_P = 1
r_A_target = r_A_norm

do it = 1,2
	r_PA = (r_A_target - r_GA)
	f_P = 1 - r_PA/(f_PA*r_Pm*MM_P)
	if (f_P < 0) then
		f_P = 0
		r_A_target = r_GA + f_PA*r_Pm*MM_P
	endif
	r_PA = (1-f_P)*f_PA*r_Pm*MM_P
	r_GA = r_A_target - r_PA
	if (r_GA < 0) then
		write(*,'(a,i2,5e12.3)') 'r_GA < 0: ',it,f_P,r_G,r_A_target,r_PA,r_GA
		r_GA = 0
		r_A_target = r_PA
!		stop
	endif
	if (r_G == 0) then
		f_G = f_G_norm
	else
		f_G = 1 - r_GA/(2*r_G)
	endif
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
		b = r_Pm - e + V*K1*Km_P
		c = -e*Km_P
		C_P = solve_C_P(a,b,c)
		!	write(*,'(a,4e12.3)') 'r_G,f_G,r_GP,C_L: ',r_G,f_G,r_GP,C_L
		!	write(*,'(a,5e12.3)') 'a,b,c,e,C_P: ',a,b,c,e,C_P
		MM_P = f_MM(C_P,Km_P,N_P)
		r_P = r_Pm*MM_P	! can check this against r_GP + V*(K2*C_L - K1*C_P)
!		write(*,'(a,2e12.3)') 'r_P check: ',r_P,r_GP + V*(K2*C_L - K1*C_P)
!		if (abs(r_P-(r_GP + V*(K2*C_L - K1*C_P))) > 1.0e-6*r_P) then
!			write(*,*) 'r_P check failed'
!			stop
!		endif
!		write(*,'(a,3e12.3)') 'r_P: ',r_Pm,MM_P,r_P
	else
		! Simply let r_P = r_Pm = r_GP + V*(K2*C_L - K1*C_P)
		r_P = r_Pm
		C_P = (K2*C_L + (r_GP - r_P)/V)/K1
		if (C_P < 0) then
			write(*,'(a,3e12.3)') 'C_P < 0: C_P, C_L: ', C_P,C_L
			stop
!			C_P = 0 
!			r_P = r_GP + V*K2*C_L
		endif
	endif
!	write(*,'(a,i2,2f8.4,e12.3)') 'it,f_G,f_P,r_P: ',it,f_G,f_P,r_P
	r_GA = r_GP
	r_PA = f_PA*(1-f_P)*r_P
	r_A = r_GA + r_PA
	r_I = f_G*r_G + f_P*r_P
	r_L = V*(K1*C_P - K2*C_L)
	if (f_P == 0) then
		write(*,*) 'f_P = 0'
		exit
	endif
enddo
if (r_I < 0) stop
M%A_rate = r_A				! production
M%I_rate = r_I				! production
M%PO_rate = r_P
M%O_rate = f_PO*r_P			! consumption
M%L_rate = r_L				! production
x = f_G
y = f_P
!write(*,'(a,3e12.3)') 'A_rate, I_rate, max: ',M%A_rate,M%I_rate,M%I_rate_max
!write(*,'(a,4e12.3)') 'r_P,_norm, r_G,_norm: ',r_P,r_P_norm,r_G,r_G_norm
end subroutine

!--------------------------------------------------------------------------
! After creating lookup tables for f_G, f_P, C_P based on running optimiser:
! HIF1 is updated with analyticSetHIF1
! PDK1 is updated with analyticSetPDK1
! HIF1, C_G      -> r_G 
! PDK1, C_O2     -> r_Pm (analyticSetPDK1)
! r_Pm, C_L, r_G -> f_G, f_P, C_P (table lookup)
! f_G, f_P, C_P  -> r_P, r_A, r_I, r_L, r_O
!--------------------------------------------------------------------------
subroutine run_optimiser
real(REAL_KIND) :: HIF1, fPDK, f_PO, f_PA, V, K1, K2, Km_O2, C_O2, C_G, r_G_max, C_L_max
real(REAL_KIND) :: r_G, r_Pm, C_L, alfa
real(REAL_KIND) :: x, y, C_P
real(REAL_KIND) :: r_P, r_A, r_I, r_O, r_L
integer :: ityp, N_O2, N, i
type(metabolism_type), target :: met
type(metabolism_type), pointer :: M

M => met
V = Vcell_cm3
ityp = 1
met = metabolic(ityp)
HIF1 = 0
C_O2 = 0.18
C_G = 5.5
fPDK = 1
N_O2 = chemo(OXYGEN)%Hill_N
Km_O2 = chemo(OXYGEN)%MM_C0
K1 = K_PL(ityp)
K2 = K_LP(ityp)
f_PO = N_PO(ityp)
f_PA = N_PA(ityp)
r_G_max = K_Ha*get_glycosis_rate(ityp,HIF1,C_G)
r_Pm = fPDK*f_MM(C_O2,Km_O2,N_O2)*chemo(OXYGEN)%max_cell_rate/f_PO
C_L_max = 1.0
C_L = C_L_max
r_G = r_G_max
N = 21
do i = 1,N
	alfa = (i-1)*1.0/(N-1)
	C_L = alfa*C_L_max
	call optimiser(ityp, r_G, C_L, r_Pm, M, x, y, C_P)
	write(*,'(a,4f8.4,e12.3)') 'x, y, C_L, C_P: ',x,y,C_L,C_P
	!write(*,'(a,4e12.3)') 'A, I, O, L: ',M%A_rate, M%I_rate, M%O_rate, M%L_rate
	r_P = 2*(1-x)*r_G - V*(K1*C_P - K2*C_L)
	r_A = 2*(1-x)*r_G + f_PA*(1-y)*r_P
	r_I = x*r_G + y*r_P
	r_O = f_PO*r_P
	r_L = V*(K1*C_P - K2*C_L)
!	write(*,'(a,4e12.3)') 'r_A, r_I: ',r_A,r_I,r_O,r_L
enddo
end subroutine

end module