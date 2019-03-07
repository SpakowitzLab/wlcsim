#include "../defines.inc"

module energies
    use precision, only: dp, epsapprox

    implicit none
    public
    integer, parameter :: NUMBER_OF_ENERGY_TYPES = 16

    integer, parameter :: chi_ = 1
    integer, parameter :: mu_ = 2
    integer, parameter :: field_ = 3
    integer, parameter :: couple_ = 4
    integer, parameter :: kap_ = 5
    integer, parameter :: bend_ = 6
    integer, parameter :: stretch_ = 7
    integer, parameter :: shear_ = 8
    integer, parameter :: maierSaupe_ = 9
    integer, parameter :: external_ = 10
    integer, parameter :: twoBody_ = 11
    integer, parameter :: twist_ = 12
    integer, parameter :: bind_ = 13
    integer, parameter :: self_ = 14
    integer, parameter :: explicitBinding_ = 15
    integer, parameter :: confine_ = 16

    type MC_energy
        real(dp) E
        real(dp) dE
        real(dp) x
        real(dp) dx
        real(dp) cof
        character(len = 8) name_str
        logical(dp) parallel_temper
    end type

    type(MC_energy), dimension(NUMBER_OF_ENERGY_TYPES) :: energyOf

contains
    subroutine set_up_energyOf()
        implicit none
        energyOf(1)%name_str='chi     '
        energyOf(2)%name_str='mu      '
        energyOf(3)%name_str='field   '
        energyOf(4)%name_str='couple  '
        energyOf(5)%name_str='kap     '
        energyOf(6)%name_str='bend    '
        energyOf(7)%name_str='stretch '
        energyOf(8)%name_str='shear   '
        energyOf(9)%name_str='maierSp '
        energyOf(10)%name_str='external'
        energyOf(11)%name_str='twoBody '
        energyOf(12)%name_str='twist   '
        energyOf(13)%name_str='bind    '
        energyOf(14)%name_str='self    '
        energyOf(15)%name_str='expl.Bnd'
        energyOf(16)%name_str='confine '

        energyOf(1)%parallel_temper = WLC_P__PT_CHI
        energyOf(2)%parallel_temper = WLC_P__PT_MU
        energyOf(3)%parallel_temper = WLC_P__PT_H
        energyOf(4)%parallel_temper = WLC_P__PT_COUPLE
        energyOf(5)%parallel_temper = WLC_P__PT_KAP
        energyOf(6)%parallel_temper = .FALSE.
        energyOf(7)%parallel_temper = .FALSE.
        energyOf(8)%parallel_temper = .FALSE.
        energyOf(9)%parallel_temper = WLC_P__PT_MAIERSAUPE
        energyOf(10)%parallel_temper = .FALSE.
        energyOf(11)%parallel_temper = WLC_P__PT_A2B
        energyOf(12)%parallel_temper = .FALSE.
        energyOf(13)%parallel_temper = .FALSE.
        energyOf(14)%parallel_temper = .FALSE.
        energyOf(15)%parallel_temper = .FALSE.
        energyOf(16)%parallel_temper = .FALSE.

        energyOf(chi_)%cof      = WLC_P__CHI
        energyOf(mu_)%cof       = WLC_P__MU
        energyOf(field_)%cof       = WLC_P__HA
        energyOf(couple_)%cof = WLC_P__HP1_BIND
        energyOf(kap_)%cof      = WLC_P__KAP
        energyOf(maierSaupe_)%cof   = WLC_P__CHI_L2
        energyOf(external_)%cof       = WLC_P__AmplitudeExternalField
        energyOf(twoBody_)%cof      = WLC_P__Amplitude2beadPotential
        energyOf(confine_)%cof = 1.0_dp
        energyOf(explicitBinding_)%cof = WLC_P__EXPLICIT_BIND_ENERGY
        energyOf(self_)%cof = 1.0_dp
        energyOf(bind_)%cof = 1.0_dp
        energyOf(bend_)%cof = 1.0_dp
        energyOf(stretch_)%cof = 1.0_dp
        energyOf(shear_)%cof = 1.0_dp
        energyOf(twist_)%cof = 1.0_dp
    end subroutine
    subroutine set_all_energy_to_zero()
        implicit none
        integer ii
        do ii = 1,NUMBER_OF_ENERGY_TYPES
            energyOf(ii)%dE=0.0_dp
            energyOf(ii)%E=0.0_dp
            energyOf(ii)%x=0.0_dp
            energyOf(ii)%dx=0.0_dp
        enddo
    end subroutine
    subroutine set_all_dEnergy_to_zero()
        implicit none
        integer ii
        do ii = 1,NUMBER_OF_ENERGY_TYPES
            energyOf(ii)%dE=0.0_dp
            energyOf(ii)%dx=0.0_dp
        enddo
    end subroutine
    subroutine sum_all_dEnergies(total_energy)
        implicit none
        integer ii
        real(dp), intent(out) :: total_energy
        total_energy=0.0_dp
        do ii = 1,NUMBER_OF_ENERGY_TYPES
            total_energy = total_energy + energyOf(ii)%dE
        enddo
    end subroutine
    subroutine accept_all_energies()
        implicit none
        integer ii
        do ii = 1,NUMBER_OF_ENERGY_TYPES
            energyOf(ii)%E=energyOf(ii)%E + energyOf(ii)%dE
            energyOf(ii)%x=energyOf(ii)%x + energyOf(ii)%dx
        enddo
    end subroutine
    subroutine calc_all_dE_from_dx()
        implicit none
        integer ii
        do ii = 1,NUMBER_OF_ENERGY_TYPES
            energyOf(ii)%dE=energyOf(ii)%dx * energyOf(ii)%cof
        enddo
    end subroutine

end module
