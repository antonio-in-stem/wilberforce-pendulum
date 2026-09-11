!> \file wilberforce_presets.f90
!! \brief The fourteen built-in parameter/initial-condition sets.
!!
!! These are the exact configurations from the `select case (p)` block of the
!! original program. To keep the refactored solver bit-for-bit identical:
!!
!!   * the physical constants are written as DEFAULT-REAL literals (no `_dp`
!!     suffix), reproducing the original single→double promotion that happens
!!     when e.g. `constantes(2)=2.2` is assigned to a double array;
!!   * the step size uses the `d0` (double) suffix, as in `dsolve(4)=0.01d0`;
!!   * initial angles keep the `atan(1.0)` single-precision form (π/4 ≈ ...).
!!
!! Change nothing here casually: these literals define the published numbers.
module wilberforce_presets
    use wilberforce_kinds, only: qp, dp
    use wilberforce_types, only: params_t, state_t
    implicit none
    private

    public :: sim_config_t, preset, initial_state, N_PRESETS

    integer, parameter :: N_PRESETS = 14

    !> A complete, ready-to-run simulation specification.
    type :: sim_config_t
        type(params_t)    :: params
        real(qp)          :: t0      = 0.0_qp
        real(qp)          :: h       = 0.01_qp
        real(qp)          :: t_final = 100.0_qp
        character(len=64) :: label   = ''
        logical           :: valid   = .true.
    end type sim_config_t

contains

    !> Initial dynamical state implied by a parameter set: released from rest
    !! (ż = θ̇ = 0) at (z0, θ0). Conversions reproduce the original qp values.
    pure function initial_state(p) result(u0)
        type(params_t), intent(in) :: p
        type(state_t) :: u0
        u0%z     = real(p%z0,     qp)
        u0%v     = 0.0_qp
        u0%theta = real(p%theta0, qp)
        u0%omega = 0.0_qp
    end function initial_state

    !> Return configuration number `id` (1..N_PRESETS); `valid=.false.` otherwise.
    function preset(id) result(cfg)
        integer, intent(in) :: id
        type(sim_config_t) :: cfg

        cfg%t0 = 0.0_qp
        cfg%h  = 0.01d0          ! matches dsolve(4)=0.01d0

        select case (id)
        case (1)
            cfg%params%delta = 2.2;  cfg%params%eps = 0.05
            cfg%params%m = 0.5;      cfg%params%inertia = 0.15
            cfg%params%z0 = 0.1;     cfg%params%theta0 = 4*atan(1.0)   ! π
            cfg%t_final = 100.0;     cfg%label = 'Set 1'
        case (2)
            cfg%params%delta = 7.0;  cfg%params%eps = 1.2
            cfg%params%m = 10.0;     cfg%params%inertia = 0.9
            cfg%params%z0 = 0.01;    cfg%params%theta0 = 2*atan(1.0)   ! π/2
            cfg%t_final = 100.0;     cfg%label = 'Set 2'
        case (3)
            cfg%params%delta = 5.4;  cfg%params%eps = 0.3
            cfg%params%m = 0.4;      cfg%params%inertia = 0.2
            cfg%params%z0 = 0.15;    cfg%params%theta0 = 3*atan(1.0)   ! 3π/4
            cfg%t_final = 10.0;      cfg%label = 'Set 3'
        case (4)
            cfg%params%delta = 0.4;  cfg%params%eps = 0.03
            cfg%params%m = 0.1;      cfg%params%inertia = 0.02
            cfg%params%z0 = 0.1;     cfg%params%theta0 = atan(1.0)     ! π/4
            cfg%t_final = 100.0;     cfg%label = 'Set 4'
        case (5)
            cfg%params%delta = 0.9;  cfg%params%eps = 0.6
            cfg%params%m = 1.0;      cfg%params%inertia = 0.1
            cfg%params%z0 = 0.2;     cfg%params%theta0 = 0.0
            cfg%t_final = 100.0;     cfg%label = 'Set 5'
        case (6)
            cfg%params%delta = 2.3;  cfg%params%eps = 1.6
            cfg%params%m = 0.1;      cfg%params%inertia = 0.8
            cfg%params%z0 = 0.08;    cfg%params%theta0 = 0.0
            cfg%t_final = 100.0;     cfg%label = 'Set 6 (eps=1.6, strong)'
        case (7)
            cfg%params%delta = 2.3;  cfg%params%eps = 0.8
            cfg%params%m = 0.1;      cfg%params%inertia = 0.8
            cfg%params%z0 = 0.08;    cfg%params%theta0 = 0.0
            cfg%t_final = 100.0;     cfg%label = 'Set 7 (eps=0.8)'
        case (8)
            cfg%params%delta = 2.3;  cfg%params%eps = 0.2
            cfg%params%m = 0.1;      cfg%params%inertia = 0.8
            cfg%params%z0 = 0.08;    cfg%params%theta0 = 0.0
            cfg%t_final = 100.0;     cfg%label = 'Set 8 (eps=0.2)'
        case (9)
            cfg%params%delta = 2.3;  cfg%params%eps = 0.05
            cfg%params%m = 0.1;      cfg%params%inertia = 0.8
            cfg%params%z0 = 0.08;    cfg%params%theta0 = 0.0
            cfg%t_final = 100.0;     cfg%label = 'Set 9 (eps=0.05)'
        case (10)
            cfg%params%delta = 2.3;  cfg%params%eps = 0.01
            cfg%params%m = 0.1;      cfg%params%inertia = 0.8
            cfg%params%z0 = 0.08;    cfg%params%theta0 = 0.0
            cfg%t_final = 100.0;     cfg%label = 'Set 10 (eps=0.01, weak)'
        case (11)
            cfg%params%delta = 1.6;  cfg%params%eps = 0.075
            cfg%params%m = 0.07;     cfg%params%inertia = 0.5
            cfg%params%z0 = 1.0;     cfg%params%theta0 = 0.0
            cfg%t_final = 100.0;     cfg%label = 'Set 11 (non-periodic)'
        case (12)
            cfg%params%delta = 1.6;  cfg%params%eps = 0.075
            cfg%params%m = 0.07;     cfg%params%inertia = 0.5
            cfg%params%z0 = 1.0;     cfg%params%theta0 = 0.374165738   ! z0*sqrt(m/I)
            cfg%t_final = 100.0;     cfg%label = 'Set 12 (periodicity condition)'
        case (13)
            cfg%params%delta = 1.6;  cfg%params%eps = 0.075
            cfg%params%m = 0.07;     cfg%params%inertia = 0.5
            cfg%params%z0 = 1.0;     cfg%params%theta0 = 0.374165738
            cfg%h = 0.001d0
            cfg%t_final = 5.0;       cfg%label = 'Set 13 (fine step, periodic)'
        case (14)
            cfg%params%delta = 7.44288e-4;  cfg%params%eps = 9.27e-3
            cfg%params%m = 0.4905;          cfg%params%inertia = 1.39e-4
            cfg%params%z0 = 0.0;            cfg%params%theta0 = 2*4*atan(1.0) ! 2π
            cfg%t_final = 40.0;             cfg%label = 'Set 14 (Berg-Marshall / Hill)'
        case default
            cfg%valid = .false.
            return
        end select

        ! Enforce the resonance condition k = m*delta/I, exactly as the original
        ! `constantes(1)=constantes(4)*constantes(2)/constantes(5)`.
        call cfg%params%enforce_resonance()
    end function preset

end module wilberforce_presets
