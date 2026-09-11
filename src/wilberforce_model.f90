!> \file wilberforce_model.f90
!! \brief Right-hand side of the coupled equations of motion.
!!
!! These four functions are the direct, behaviour-preserving descendants of the
!! original auxiliary functions `d1y`, `d2y`, `d1x`, `d2x`. The arithmetic is
!! reproduced verbatim (same operations, same order, same 0.5 literal) so the
!! refactored solver is bit-for-bit identical to the original.
!!
!!   z̈ = −(k/m) z − (ε/2m) θ
!!   θ̈ = −(δ/I) θ − (ε/2I) z
module wilberforce_model
    use wilberforce_kinds, only: qp
    use wilberforce_types, only: params_t
    implicit none
    private

    public :: deriv_z, accel_z, deriv_theta, accel_theta

contains

    !> Kinematic relation ż = v  (original d1y).
    pure function deriv_z(v) result(dz)
        real(qp), intent(in) :: v
        real(qp) :: dz
        dz = v
    end function deriv_z

    !> Longitudinal acceleration z̈  (original d2y).
    pure function accel_z(z, theta, p) result(a)
        real(qp),       intent(in) :: z, theta
        type(params_t), intent(in) :: p
        real(qp) :: a
        real(qp) :: k, m, e
        k = real(p%k,   qp)
        m = real(p%m,   qp)
        e = real(p%eps, qp)
        a = -k/m * z - 0.5_qp*e/m * theta
    end function accel_z

    !> Kinematic relation θ̇ = b  (original d1x).
    pure function deriv_theta(b) result(dth)
        real(qp), intent(in) :: b
        real(qp) :: dth
        dth = b
    end function deriv_theta

    !> Angular acceleration θ̈  (original d2x).
    pure function accel_theta(z, theta, p) result(a)
        real(qp),       intent(in) :: z, theta
        type(params_t), intent(in) :: p
        real(qp) :: a
        real(qp) :: s, w, e
        s = real(p%delta,   qp)
        w = real(p%inertia, qp)
        e = real(p%eps,     qp)
        a = -s/w * theta - 0.5_qp*e/w * z
    end function accel_theta

end module wilberforce_model
