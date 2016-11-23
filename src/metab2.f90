! To test simple cell metabolism model
! Concentration of ATP varies by < 10%  https://en.wikipedia.org/wiki/Glycolysis#Intermediates_for_other_pathways

! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM

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

real(REAL_KIND) :: Hill_Km_G     ! Hill Km for dependence of glycolysis rate on glucose
real(REAL_KIND) :: Hill_N_G      ! Hill N for dependence of glycolysis rate on glucose
real(REAL_KIND) :: K_H1(MAX_CELLTYPES)     ! HIF-1 k1
real(REAL_KIND) :: K_H2(MAX_CELLTYPES)     ! HIF-1 k2
real(REAL_KIND) :: K_Ha     ! HIF-1 ka
real(REAL_KIND) :: K_Hb(MAX_CELLTYPES)     ! HIF-1 kb
real(REAL_KIND) :: K_PDK(MAX_CELLTYPES)    ! K_PDK

!real(REAL_KIND) :: Vcell     ! cell volume
real(REAL_KIND) :: K_PL(MAX_CELLTYPES)     ! P -> L
real(REAL_KIND) :: K_LP(MAX_CELLTYPES)     ! L -> P

!real(REAL_KIND) :: N_GA     ! number of ATP molecules generated per glucose molecule in glycosis
!real(REAL_KIND) :: N_GI     ! number of intermediate molecules generated per glucose molecule in glycosis
!real(REAL_KIND) :: F_PO_BASE	! base level of pyruvate oxidation (fraction of glycolysis rate)
!real(REAL_KIND) :: N_PA     ! number of ATP molecules generated per pyruvate molecule in pyruvate oxidation
!real(REAL_KIND) :: N_PI     ! number of intermediate molecules generated per pyruvate molecule in pyruvate oxidation
!real(REAL_KIND) :: N_PO     ! number of O2 molecules consumed per pyruvate molecule in pyruvate oxidation

real(REAL_KIND) :: A_rate_base(MAX_CELLTYPES)	! total rate of production of ATP under full nutrition

!real(REAL_KIND) :: C_O, C_G, C_P, C_L, H, M_I

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_glycosis_rate(ityp, H, G) result(rate)
integer :: ityp
real(REAL_KIND) :: H, G, rate

rate = (1 + K_Hb(ityp)*H)*G**Hill_N_G/(G**Hill_N_G + Hill_Km_G**Hill_N_G)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine SetupMetabolism
real(REAL_KIND) :: G_rate, PO_rate(MAX_CELLTYPES)
integer :: ityp

!Hill_Km_G = 220e-3		! 220 uM -> mM
!Hill_N_G = 1
Hill_Km_G = chemo(GLUCOSE)%MM_C0
Hill_N_G = chemo(GLUCOSE)%Hill_N
!N_GA = 2
!N_GI = 1
!N_PA = 18
!N_PI = 4
!N_PO = 3
!K_H1 = 200
!K_H2 = 0.001
!K_Ha = 3.80e-11     ! mumol/cell/s = 3.80e-17*1.0e6
K_Ha = chemo(GLUCOSE)%max_cell_rate
!K_Hb = 0.2

!K_PL = 1.0e-2		! /s to convert conc in mM to mM/s, then mM/s.Vcell_cm3 -> mumol/s
!K_LP = 0.2e-2

!F_PO_BASE = 0.1

A_rate_base = K_Ha*(N_GA + F_PO_BASE*(2 - N_GI)*N_PA)	! NOT normalised
ATPg = ATPg*A_rate_base
ATPs = ATPs*A_rate_base

metabolic%HIF1 = 0
metabolic%PDK1 = 1
G_rate = K_Ha
!A_rate_base = N_GA*G_rate + N_PA*PO_rate
PO_rate = (A_rate_base - N_GA*G_rate)/N_PA 
do ityp = 1,Ncelltypes
	metabolic(ityp)%I_rate_max = N_GI(ityp)*G_rate + N_PI(ityp)*PO_rate(ityp)
enddo
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
real(REAL_KIND) :: a, b, c, P0

a = K_PDK(ityp)
b = (1 - H)
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
!--------------------------------------------------------------------------
subroutine get_metab_rates(ityp, H, PDK, Cin, G_rate, PP_rate, PO_rate)
integer :: ityp
real(REAL_KIND) :: H, PDK, Cin(:)
real(REAL_KIND) :: G_rate       ! glycolysis rate
real(REAL_KIND) :: PP_rate      ! pyruvate production rate by glycolysis
real(REAL_KIND) :: PO_rate      ! pyruvate oxidation rate
real(REAL_KIND) :: G_norm       ! normalised glycolysis rate
real(REAL_KIND) :: pyruvate
!real(REAL_KIND) :: L_rate       ! lactate production rate
!real(REAL_KIND) :: A_rate       ! total rate of production of ATP
!real(REAL_KIND) :: I_rate       ! total rate of production of anabolic intermediates ( <= 3.31E-10 )
!real(REAL_KIND) :: O_rate       ! rate of consumption of O2 ( <= 7.35E-11)

!call get_HIF1(cp%Cin(OXYGEN),cp%HIF1)
G_norm = get_glycosis_rate(ityp,H,Cin(GLUCOSE))
G_rate = K_Ha*G_norm
PP_rate = (2 - N_GI(ityp))*G_rate
PO_rate = (A_rate_base(ityp) - G_rate*N_GA(ityp))/N_PA(ityp)
PO_rate = min(PO_rate,PP_rate + Vcell_cm3*K_LP(ityp)*Cin(LACTATE))
PO_rate = PDK*PO_rate
pyruvate = (PP_rate - PO_rate + Vcell_cm3*K_LP(ityp)*Cin(LACTATE))/Vcell_cm3
G_rate = G_rate							! rate of consumption of glucose
PP_rate = PP_rate
PO_rate = PO_rate
!L_rate = PP_rate - PO_rate				! rate of production of lactate
!A_rate = N_GA*G_rate + N_PA*PO_rate	! rate of production of ATP
!I_rate = N_GI*G_rate + N_PI*PO_rate	! rate of production of intermediates
!O_rate = N_PO*PO_rate					! rate of consumption of oxygen
end subroutine

!--------------------------------------------------------------------------
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
	call get_metab_rates(ityp, M%HIF1, M%PDK1, Caverage(1:MAX_CHEMO), M%G_rate, M%PP_rate, M%PO_rate)
	M%L_rate = M%PP_rate - M%PO_rate			! rate of production of lactate
	M%A_rate = N_GA(ityp)*M%G_rate + N_PA(ityp)*M%PO_rate	! rate of production of ATP
	M%I_rate = N_GI(ityp)*M%G_rate + N_PI(ityp)*M%PO_rate	! rate of production of intermediates
	M%O_rate = N_PO(ityp)*M%PO_rate					! rate of consumption of oxygen
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
