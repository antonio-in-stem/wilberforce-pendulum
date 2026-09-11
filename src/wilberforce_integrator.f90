!> \file wilberforce_integrator.f90
!! \brief Classical fourth-order Runge–Kutta integrator for the 4-D system
!!        u = (z, ż, θ, θ̇), plus the simulation driver.
!!
!! `rk4_step` retains the four stages extracted from the inner loop
!! of the original `RK` subroutine, with a corrected RK4 weight. The driver retains the
!! loop: advance the state, evaluate the analytic solution at the same instant,
!! accumulate the energy budget, and stream every step to disk.
module wilberforce_integrator
    use wilberforce_kinds,    only: qp, dp
    use wilberforce_types,    only: params_t, state_t
    use wilberforce_model,    only: deriv_z, accel_z, deriv_theta, accel_theta
    use wilberforce_analytic, only: analytic_t, analytic_init
    use wilberforce_energy,   only: energy_t, energy_of
    use wilberforce_io,       only: datastore_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: rk4_step, run_simulation, step_count

contains

    !> Full steps within the requested horizon (no extra step at an exact endpoint).
    !! Ratios within double-precision roundoff of an integer are snapped to it,
    !! because the supplied timestep literals originate in double precision.
    function step_count(t0, t_final, h) result(n)
        real(qp), intent(in) :: t0, t_final, h
        integer :: n
        real(qp) :: ratio, nearest
        if (.not. all(ieee_is_finite([t0,t_final,h]))) error stop 'times must be finite'
        if (h <= 0.0_qp .or. t_final < t0) error stop 'invalid time interval or timestep'
        ratio=(t_final-t0)/h
        if (ratio > real(huge(n)-1,qp)) error stop 'step count exceeds integer capacity'
        nearest=anint(ratio)
        if (abs(ratio-nearest) <= 16.0_qp*epsilon(1.0_dp)*max(1.0_qp,ratio)) ratio=nearest
        n=int(ratio)
    end function step_count

    !> Classical RK4. The weight is evaluated in the state precision.
    !! The RK4 weight is evaluated in the state precision to avoid a
    !! single-precision phase-error floor.
    pure function rk4_step(p, u, h) result(un)
        type(params_t), intent(in) :: p
        type(state_t),  intent(in) :: u
        real(qp),       intent(in) :: h
        type(state_t) :: un
        real(qp) :: z, v, x, b
        real(qp) :: k1y, k1v, k2y, k2v, k3y, k3v, k4y, k4v
        real(qp) :: j1x, j1b, j2x, j2b, j3x, j3b, j4x, j4b

        z = u%z;  v = u%v;  x = u%theta;  b = u%omega

        ! Stage 1
        k1y = deriv_z(v)
        k1v = accel_z(z, x, p)
        j1x = deriv_theta(b)
        j1b = accel_theta(z, x, p)

        ! Stage 2  (half step, slopes k1/j1)
        k2y = deriv_z(v + h*k1v/2.)
        k2v = accel_z(z + h*k1y/2., x + h*j1x/2., p)
        j2x = deriv_theta(b + h*j1b/2.)
        j2b = accel_theta(z + h*k1y/2., x + h*j1x/2., p)

        ! Stage 3  (half step, slopes k2/j2)
        k3y = deriv_z(v + h*k2v/2.)
        k3v = accel_z(z + h*k2y/2., x + h*j2x/2., p)
        j3x = deriv_theta(b + h*j2b/2.)
        j3b = accel_theta(z + h*k2y/2., x + h*j2x/2., p)

        ! Stage 4  (full step, slopes k3/j3)
        k4y = deriv_z(v + h*k3v)
        k4v = accel_z(z + h*k3y, x + h*j3x, p)
        j4x = deriv_theta(b + h*j3b)
        j4b = accel_theta(z + h*k3y, x + h*j3x, p)

        ! Weighted combination
        un%z     = z + (1.0_qp/6.0_qp)*h*(k1y + 2.*k2y + 2.*k3y + k4y)
        un%v     = v + (1.0_qp/6.0_qp)*h*(k1v + 2.*k2v + 2.*k3v + k4v)
        un%theta = x + (1.0_qp/6.0_qp)*h*(j1x + 2.*j2x + 2.*j3x + j4x)
        un%omega = b + (1.0_qp/6.0_qp)*h*(j1b + 2.*j2b + 2.*j3b + j4b)
    end function rk4_step

    !> Run the full simulation: integrate, compare against the analytic solution,
    !! track the energy budget, and write every step to `store`.
    subroutine run_simulation(p, u0, t0, h, n, store)
        type(params_t),    intent(in)    :: p
        type(state_t),     intent(in)    :: u0
        real(qp),          intent(in)    :: t0, h
        integer,           intent(in)    :: n
        type(datastore_t), intent(inout) :: store
        type(analytic_t) :: ana
        type(state_t)    :: u, ua
        type(energy_t)   :: en
        real(qp)         :: t
        integer          :: i

        if (.not. all(ieee_is_finite([t0,h]))) error stop 'times must be finite'
        if (h <= 0.0_qp .or. n < 0) error stop 'invalid integration request'
        ana = analytic_init(p, u0)
        if (h*ana%w1 > sqrt(8.0_qp)) error stop 'RK4 timestep outside oscillatory stability region' 
        u   = u0
        do i = 1, n
            t  = t0 + h*real(i, qp)
            u  = rk4_step(p, u, h)
            ua = ana%eval(t-t0)
            if (.not. all(ieee_is_finite([u%z,u%v,u%theta,u%omega]))) &
                error stop 'non-finite integrated state'
            en = energy_of(u, p)
            call store%write_record(t, u, ua, en)
        end do
    end subroutine run_simulation

end module wilberforce_integrator
