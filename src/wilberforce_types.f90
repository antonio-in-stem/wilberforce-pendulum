!> \file wilberforce_types.f90
!! \brief Derived types describing the physical system and its dynamical state.
!!
!! The original program packed everything into two bare arrays, `constantes(7)`
!! and `dsolve(99)`, and relied on magic indices (constantes(1) is k,
!! constantes(4) is m, ...). This module preserves that data *content* but gives
!! every quantity a name, which removes an entire class of index bugs and makes
!! the rest of the package self-documenting.
module wilberforce_types
    use wilberforce_kinds, only: dp, qp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: params_t, state_t, valid_params

    !> Physical parameters of a Wilberforce pendulum plus the initial amplitudes.
    !! Stored in double precision, exactly as the original `constantes` array.
    type :: params_t
        real(dp) :: k     = 0.0_dp   !< longitudinal spring constant   (constantes(1))
        real(dp) :: delta = 0.0_dp   !< torsional constant             (constantes(2))
        real(dp) :: eps   = 0.0_dp   !< coupling constant              (constantes(3))
        real(dp) :: m     = 0.0_dp   !< bob mass                       (constantes(4))
        real(dp) :: inertia = 0.0_dp !< moment of inertia I            (constantes(5))
        real(dp) :: z0    = 0.0_dp   !< initial longitudinal position  (constantes(6))
        real(dp) :: theta0 = 0.0_dp  !< initial angular position       (constantes(7))
    contains
        procedure :: enforce_resonance
    end type params_t

    !> Instantaneous dynamical state u = (z, z', theta, theta').
    !! Integrated in quadruple precision (original `real*16`).
    type :: state_t
        real(qp) :: z     = 0.0_qp   !< longitudinal position   (dsolve(2), "y")
        real(qp) :: v     = 0.0_qp   !< longitudinal velocity   (dsolve(3), "v")
        real(qp) :: theta = 0.0_qp   !< angular position        (dsolve(6), "x")
        real(qp) :: omega = 0.0_qp   !< angular velocity        (dsolve(7), "b")
    end type state_t

contains

    !> Impose the resonance condition k = m*delta/I on this parameter set,
    !! reproducing `constantes(1)=constantes(4)*constantes(2)/constantes(5)`.
    subroutine enforce_resonance(self)
        class(params_t), intent(inout) :: self
        if (.not. ieee_is_finite(self%inertia) .or. self%inertia <= 0.0_dp) &
            error stop 'inertia must be finite and positive'
        self%k = self%m * self%delta / self%inertia
    end subroutine enforce_resonance

    !> Supported domain: finite, positive-definite, resonant linear oscillator.
    pure logical function valid_params(p) result(ok)
        type(params_t), intent(in) :: p
        real(qp) :: k, d, e, m, j, scale
        ok = .false.
        if (.not. all(ieee_is_finite([p%k,p%delta,p%eps,p%m,p%inertia,p%z0,p%theta0]))) return
        if (min(p%k,p%delta,p%m,p%inertia) <= 0.0_dp) return
        k=real(p%k,qp); d=real(p%delta,qp); e=real(p%eps,qp)
        m=real(p%m,qp); j=real(p%inertia,qp)
        if (e*e >= 4.0_qp*k*d) return
        scale=max(abs(k/m),abs(d/j))
        if (abs(k/m-d/j) > 64.0_qp*epsilon(1.0_dp)*scale) return
        ok = .true.
    end function valid_params

end module wilberforce_types
