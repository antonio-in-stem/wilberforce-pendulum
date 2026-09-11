!> \file test_wilberforce.f90
!! \brief Regression and physics test suite for the Wilberforce package.
!!
!! Exercises the numerical core without touching the filesystem:
!!   1. resonance condition k = mδ/I is enforced for every preset;
!!   2. total energy is conserved along a full RK4 trajectory;
!!   3. the RK4 solution tracks the closed-form normal-mode solution;
!!   4. the periodicity condition θ0 = z0√(m/I) makes the second mode vanish;
!!   5. a numerical regression anchor pins the first integrated point of
!!      preset 14 to its published value.
!!
!! Exit status is non-zero if any check fails (suitable for CI / ctest).
program test_wilberforce
    use wilberforce_kinds,      only: qp, dp
    use wilberforce_types,      only: params_t, state_t, valid_params
    use wilberforce_presets,    only: sim_config_t, preset, initial_state, N_PRESETS
    use wilberforce_analytic,   only: analytic_t, analytic_init
    use wilberforce_integrator, only: rk4_step, step_count
    use wilberforce_energy,     only: energy_t, energy_of
    implicit none

    integer :: n_fail = 0

    write(*,'(A)') '== Wilberforce test suite =='
    call test_resonance()
    call test_energy_conservation()
    call test_analytic_vs_numeric()
    call test_periodicity_condition()
    call test_regression_anchor()
    call test_step_count()
    call test_zero_coupling()
    call test_initial_velocities()
    call test_convergence()

    write(*,'(A)') '----------------------------------------'
    if (n_fail == 0) then
        write(*,'(A)') 'ALL TESTS PASSED'
    else
        write(*,'(I0,A)') n_fail, ' TEST(S) FAILED'
        error stop 1
    end if

contains

    subroutine check(name, ok)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: ok
        if (ok) then
            write(*,'(A,A)') '  [PASS] ', name
        else
            write(*,'(A,A)') '  [FAIL] ', name
            n_fail = n_fail + 1
        end if
    end subroutine check

    !> (1) Every preset must satisfy the physical resonance condition
    !!     ω_z² = ω_θ², i.e. k/m = δ/I, to double-precision tolerance.
    subroutine test_resonance()
        type(sim_config_t) :: cfg
        integer :: id
        logical :: ok
        real(dp) :: wz2, wth2
        ok = .true.
        do id = 1, N_PRESETS
            cfg = preset(id)
            wz2  = cfg%params%k / cfg%params%m           ! ω_z²
            wth2 = cfg%params%delta / cfg%params%inertia ! ω_θ²
            if (abs(wz2 - wth2) > 1.0e-13_dp*abs(wth2)) ok = .false.
        end do
        call check('resonance ω_z²=ω_θ² (k/m=δ/I) for all 14 presets', ok)
    end subroutine test_resonance

    !> (2) Total energy is conserved along the trajectory of preset 14.
    subroutine test_energy_conservation()
        type(sim_config_t) :: cfg
        type(state_t)      :: u
        type(energy_t)     :: en
        real(qp) :: e0, emin, emax, drift
        integer  :: i, n
        cfg = preset(14)
        u   = initial_state(cfg%params)
        n   = step_count(cfg%t0, cfg%t_final, cfg%h)
        en  = energy_of(u, cfg%params)
        e0  = en%total(); emin = e0; emax = e0
        do i = 1, n
            u  = rk4_step(cfg%params, u, cfg%h)
            en = energy_of(u, cfg%params)
            emin = min(emin, en%total()); emax = max(emax, en%total())
        end do
        drift = (emax - emin)/abs(e0)
        write(*,'(A,ES10.3)') '         relative energy drift = ', real(drift)
        call check('energy conserved to < 1e-6 (preset 14, 4000 steps)', drift < 1.0e-6_qp)
    end subroutine test_energy_conservation

    !> (3) RK4 tracks the analytic normal-mode solution.
    subroutine test_analytic_vs_numeric()
        type(sim_config_t) :: cfg
        type(analytic_t)   :: ana
        type(state_t)      :: u, ua
        real(qp) :: t, dmax, amp
        integer  :: i, n
        cfg = preset(14)
        ana = analytic_init(cfg%params)
        u   = initial_state(cfg%params)
        n   = step_count(cfg%t0, cfg%t_final, cfg%h)
        dmax = 0.0_qp; amp = abs(real(cfg%params%theta0, qp))
        do i = 1, n
            t  = cfg%t0 + cfg%h*real(i, qp)
            u  = rk4_step(cfg%params, u, cfg%h)
            ua = ana%eval(t)
            dmax = max(dmax, abs(u%theta - ua%theta))
        end do
        write(*,'(A,ES10.3)') '         max|θ_num - θ_ana| = ', real(dmax)
        call check('numeric matches analytic θ(t) to < 1e-3', dmax < 1.0e-3_qp)
    end subroutine test_analytic_vs_numeric

    !> (4) Periodicity condition θ0 = z0√(m/I) kills the second normal mode.
    subroutine test_periodicity_condition()
        type(sim_config_t) :: cfg
        type(analytic_t)   :: ana
        cfg = preset(12)   ! θ0 set to z0*sqrt(m/I)
        ana = analytic_init(cfg%params)
        write(*,'(A,ES10.3)') '         |D| (2nd-mode amplitude) = ', real(abs(ana%D))
        call check('periodicity condition suppresses mode 2 (|D| < 1e-3)', abs(ana%D) < 1.0e-3_qp)
    end subroutine test_periodicity_condition

    !> (5) Regression anchor: first integrated point of preset 14.
    subroutine test_regression_anchor()
        type(sim_config_t) :: cfg
        type(state_t)      :: u
        real(qp) :: z_expected, az0, h, wz2, wt2
        real(qp) :: relerr
        cfg = preset(14)
        u   = initial_state(cfg%params)
        h=cfg%h
        az0=-real(cfg%params%eps,qp)*u%theta/(2.0_qp*real(cfg%params%m,qp))
        wz2=real(cfg%params%k,qp)/real(cfg%params%m,qp)
        wt2=real(cfg%params%delta,qp)/real(cfg%params%inertia,qp)
        ! Independent degree-four Taylor polynomial for this zero-velocity IC.
        z_expected=az0*h*h/2.0_qp-az0*(wz2+wt2)*h**4/24.0_qp
        u   = rk4_step(cfg%params, u, cfg%h)
        relerr = abs(u%z - z_expected)/abs(z_expected)
        write(*,'(A,ES10.3)') '         first-step z relerr = ', real(relerr)
        call check('RK4 first step matches independent matrix Taylor polynomial', relerr < 1.0e-28_qp)
    end subroutine test_regression_anchor

    subroutine test_step_count()
        call check('exact endpoint has no extra step', step_count(0.0_qp,1.0_qp,0.25_qp)==4)
        call check('empty interval has no steps', step_count(2.0_qp,2.0_qp,0.25_qp)==0)
        call check('nonintegral horizon uses full steps only', step_count(0.0_qp,1.0_qp,0.3_qp)==3)
    end subroutine

    subroutine test_zero_coupling()
        type(sim_config_t) :: c
        type(analytic_t) :: a
        type(state_t) :: u
        c=preset(7); c%params%eps=0.0_dp
        a=analytic_init(c%params); u=a%eval(1.0_qp)
        call check('zero coupling gives finite independent oscillator', &
            abs(u%z-real(c%params%z0,qp)*cos(a%w_bar))<1.0e-28_qp .and. abs(u%theta)<1.0e-28_qp)
        c%params%eps=-0.8_dp
        call check('negative stable coupling accepted',valid_params(c%params))
        c%params%eps=10.0_dp
        call check('unstable coupling rejected',.not.valid_params(c%params))
        c%params%inertia=0.0_dp
        call check('zero inertia rejected',.not.valid_params(c%params))
    end subroutine

    subroutine test_initial_velocities()
        type(sim_config_t) :: c
        type(analytic_t) :: a
        type(state_t) :: u,v
        c=preset(7); u=state_t(0.2_qp,0.3_qp,-0.4_qp,0.5_qp)
        a=analytic_init(c%params,u); v=a%eval(0.0_qp)
        call check('analytic reference uses supplied initial positions and velocities', &
            maxval(abs([u%z-v%z,u%v-v%v,u%theta-v%theta,u%omega-v%omega]))<1.0e-28_qp)
    end subroutine

    subroutine test_convergence()
        type(sim_config_t) :: c
        type(analytic_t) :: a
        type(state_t) :: u,v
        real(qp) :: h,err(3)
        integer :: j,i
        c=preset(14); a=analytic_init(c%params)
        do j=1,3
            h=0.04_qp/real(2**(j-1),qp); u=initial_state(c%params); err(j)=0.0_qp
            do i=1,nint(40.0_qp/h)
                u=rk4_step(c%params,u,h); v=a%eval(i*h)
                err(j)=max(err(j),abs(u%theta-v%theta))
            end do
        end do
        write(*,'(A,2F10.4)') '         error ratios h/h2: ',err(1)/err(2),err(2)/err(3)
        call check('fourth-order convergence, not single-weight phase plateau', &
            all(err(1:2)/err(2:3)>14.0_qp) .and. all(err(1:2)/err(2:3)<18.0_qp))
    end subroutine

end program test_wilberforce
