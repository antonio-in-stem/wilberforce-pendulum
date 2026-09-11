!> \file wilberforce_analytic.f90
!! \brief Closed-form normal-mode solution at resonance.
!!
!! Reproduces the analytic branch of the original `RK` subroutine: the two
!! normal-mode frequencies w1, w2, the natural frequency w_bar = sqrt(delta/I),
!! the mode amplitudes B, D fixed by the initial conditions, and the evaluation
!! of z(t), z'(t), theta(t), theta'(t). See README.md and paper/en.pdf.
module wilberforce_analytic
    use wilberforce_kinds, only: qp
    use wilberforce_types, only: params_t, state_t, valid_params
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: analytic_t, analytic_init

    !> Pre-computed constants of the analytic solution for one parameter set.
    type :: analytic_t
        real(qp) :: w1    = 0.0_qp  !< normal-mode frequency  sqrt(w0^2 + eps/(2 sqrt(mI)))
        real(qp) :: w2    = 0.0_qp  !< normal-mode frequency  sqrt(w0^2 - eps/(2 sqrt(mI)))
        real(qp) :: w_bar = 0.0_qp  !< natural frequency w0 = sqrt(delta/I)
        real(qp) :: B     = 0.0_qp  !< amplitude of mode 1
        real(qp) :: D     = 0.0_qp  !< amplitude of mode 2
        real(qp) :: w     = 0.0_qp  !< inertia I (needed to rebuild z)
        real(qp) :: e     = 0.0_qp  !< coupling eps (needed to rebuild z)
        real(qp) :: scale = 1.0_qp !< signed sqrt(I/m), finite at zero coupling
        real(qp) :: Bv = 0.0_qp, Dv = 0.0_qp !< initial modal velocities
    contains
        procedure :: eval => analytic_eval
    end type analytic_t

contains

    !> Build the analytic solution from a (resonant) parameter set and the
    !! initial amplitudes z0 = p%z0, theta0 = p%theta0.
    function analytic_init(p, u0) result(a)
        type(params_t), intent(in) :: p
        type(state_t), intent(in), optional :: u0
        type(analytic_t) :: a
        type(state_t) :: initial
        real(qp) :: s, e, m, w, sign_e

        if (.not. valid_params(p)) error stop 'analytic solution requires stable resonant parameters'
        s=real(p%delta,qp); e=real(p%eps,qp); m=real(p%m,qp); w=real(p%inertia,qp)
        initial%z=real(p%z0,qp); initial%theta=real(p%theta0,qp)
        if (present(u0)) initial=u0
        if (.not. all(ieee_is_finite([initial%z,initial%v,initial%theta,initial%omega]))) &
            error stop 'initial state must be finite'
        sign_e=1.0_qp
        if (e < 0.0_qp) sign_e=-1.0_qp
        a%w1=sqrt(s/w+abs(e)/(2.0_qp*sqrt(m*w)))
        a%w2=sqrt(s/w-abs(e)/(2.0_qp*sqrt(m*w)))
        a%w_bar=sqrt(s/w)
        a%scale=sign_e*sqrt(w/m)
        a%B=0.5_qp*(initial%theta+initial%z/a%scale)
        a%D=0.5_qp*(initial%theta-initial%z/a%scale)
        a%Bv=0.5_qp*(initial%omega+initial%v/a%scale)
        a%Dv=0.5_qp*(initial%omega-initial%v/a%scale)
        a%w=w; a%e=e
    end function analytic_init

    !> Evaluate the analytic state at time t.
    function analytic_eval(self, t) result(u)
        class(analytic_t), intent(in) :: self
        real(qp),          intent(in) :: t
        type(state_t) :: u
        real(qp) :: q1,q2,v1,v2
        q1=self%B*cos(self%w1*t)+self%Bv*sin(self%w1*t)/self%w1
        q2=self%D*cos(self%w2*t)+self%Dv*sin(self%w2*t)/self%w2
        v1=-self%w1*self%B*sin(self%w1*t)+self%Bv*cos(self%w1*t)
        v2=-self%w2*self%D*sin(self%w2*t)+self%Dv*cos(self%w2*t)
        u%theta=q1+q2; u%omega=v1+v2
        u%z=self%scale*(q1-q2); u%v=self%scale*(v1-v2)
    end function analytic_eval

end module wilberforce_analytic
