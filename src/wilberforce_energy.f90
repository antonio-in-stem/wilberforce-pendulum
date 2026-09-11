!> \file wilberforce_energy.f90
!! \brief Energy diagnostics for the coupled mechanical system.
!!
!! Conservation of the total energy is the package's primary correctness
!! indicator for the numerical integration. The five contributions reproduce the
!! original `e_kz, e_pz, e_ka, e_pa, e_ac` expressions verbatim.
module wilberforce_energy
    use wilberforce_kinds, only: qp
    use wilberforce_types, only: params_t, state_t
    implicit none
    private

    public :: energy_t, energy_of

    !> Decomposition of the instantaneous mechanical energy.
    type :: energy_t
        real(qp) :: kin_long = 0.0_qp   !< ½ m ż²
        real(qp) :: pot_long = 0.0_qp   !< ½ k z²
        real(qp) :: kin_ang  = 0.0_qp   !< ½ I θ̇²
        real(qp) :: pot_ang  = 0.0_qp   !< ½ δ θ²
        real(qp) :: coupling = 0.0_qp   !< ½ ε z θ
    contains
        procedure :: longitudinal => energy_longitudinal
        procedure :: angular      => energy_angular
        procedure :: total        => energy_total
    end type energy_t

contains

    !> Energy contributions for a given state and parameter set.
    function energy_of(u, p) result(en)
        type(state_t),  intent(in) :: u
        type(params_t), intent(in) :: p
        type(energy_t) :: en
        real(qp) :: m, k, w, s, e
        m = real(p%m,       qp)
        k = real(p%k,       qp)
        w = real(p%inertia, qp)
        s = real(p%delta,   qp)
        e = real(p%eps,     qp)
        en%kin_long = 0.5_qp*m*u%v**2
        en%pot_long = 0.5_qp*k*u%z**2
        en%kin_ang  = 0.5_qp*w*u%omega**2
        en%pot_ang  = 0.5_qp*s*u%theta**2
        en%coupling = 0.5_qp*e*u%z*u%theta
    end function energy_of

    pure function energy_longitudinal(self) result(r)
        class(energy_t), intent(in) :: self
        real(qp) :: r
        r = self%kin_long + self%pot_long
    end function energy_longitudinal

    pure function energy_angular(self) result(r)
        class(energy_t), intent(in) :: self
        real(qp) :: r
        r = self%kin_ang + self%pot_ang
    end function energy_angular

    pure function energy_total(self) result(r)
        class(energy_t), intent(in) :: self
        real(qp) :: r
        r = self%kin_long + self%pot_long + self%kin_ang + self%pot_ang + self%coupling
    end function energy_total

end module wilberforce_energy
