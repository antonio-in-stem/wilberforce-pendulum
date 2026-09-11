!> \file main.f90
!! \brief Driver program for the Wilberforce pendulum solver.
!!
!! Usage:
!!     wilberforce [PRESET] [--plot] [--outdir DIR]
!!
!!   PRESET     integer 1..14 selecting a built-in configuration (default: 14,
!!              the Berg–Marshall / Christian-Hill comparison set, matching the
!!              original program's hard-coded `p=14`).
!!   --plot     after writing the data, run Gnuplot head-less to render the PNGs.
!!   --outdir   parent directory for the `data/`, `gpl/`, `graficas/` folders
!!              (default: generated).
!!
!! This is the modular replacement for the original monolithic `program main`
!! plus `sets`. Column formats are preserved; corrected RK4 values differ from the historical version.
program wilberforce_main
    use wilberforce_kinds,      only: qp, PI
    use wilberforce_types,      only: state_t
    use wilberforce_presets,    only: sim_config_t, preset, initial_state, N_PRESETS
    use wilberforce_analytic,   only: analytic_t, analytic_init
    use wilberforce_integrator, only: run_simulation, step_count
    use wilberforce_io,         only: datastore_t, validate_path
    use wilberforce_gnuplot,    only: plot_all
    implicit none

    type(sim_config_t) :: cfg
    type(state_t)      :: u0
    type(analytic_t)   :: ana
    type(datastore_t)  :: store
    integer            :: id, n, nargs, i, ios
    logical            :: do_plot, seen_preset, seen_outdir
    character(len=:), allocatable :: arg, outdir, data_dir, gpl_dir, out_dir
    real(qp)           :: w0, tbeat

    ! ---- defaults -------------------------------------------------------------
    id      = 14
    do_plot = .false.
    seen_preset=.false.; seen_outdir=.false.
    outdir  = 'generated'

    ! ---- parse command line ---------------------------------------------------
    nargs = command_argument_count()
    i = 1
    do while (i <= nargs)
        call read_arg(i, arg)
        select case (trim(arg))
        case ('--plot')
            do_plot = .true.
        case ('--outdir')
            if (i + 1 > nargs) then
                write(*,'(A)') 'Error: --outdir requires a directory argument'
                call print_usage(); error stop 1
            end if
            i = i + 1
            call read_arg(i, outdir)
            if (seen_outdir) error stop 'duplicate --outdir'
            seen_outdir=.true.
            call validate_path(outdir)
        case ('-h', '--help')
            call print_usage(); stop
        case default
            if (len_trim(arg)==0) error stop 'empty argument'
            if (verify(arg,'0123456789')/=0 .or. seen_preset) error stop 'invalid or duplicate preset'
            seen_preset=.true.
            read(arg, *, iostat=ios) id
            if (ios /= 0) then
                write(*,'(A)') 'Unrecognised argument: '//trim(arg)
                call print_usage(); error stop 1
            end if
        end select
        i = i + 1
    end do

    call validate_path(outdir)
    data_dir = trim(outdir)//'/data'
    gpl_dir  = trim(outdir)//'/gpl'
    out_dir  = trim(outdir)//'/graficas'

    ! ---- banner ---------------------------------------------------------------
    write(*,'(A)') ''
    write(*,'(A)') '----------------------------------------------------------------------------'
    write(*,'(A)') ' Wilberforce Pendulum  ·  Física Computacional  ·  RK4 (quad precision)'
    write(*,'(A)') '----------------------------------------------------------------------------'
    write(*,'(A)') ''

    ! ---- load configuration ---------------------------------------------------
    cfg = preset(id)
    if (.not. cfg%valid) then
        write(*,'(A,I0,A,I0)') ' Preset ', id, ' is out of range 1..', N_PRESETS
        error stop 1
    end if
    u0 = initial_state(cfg%params)
    n  = step_count(cfg%t0, cfg%t_final, cfg%h)

    ! ---- report the physics ---------------------------------------------------
    ana = analytic_init(cfg%params)
    w0  = sqrt(real(cfg%params%delta,qp)/real(cfg%params%inertia,qp))
    write(*,'(A,I0,A,A)')      ' Preset            : ', id, '  ', trim(cfg%label)
    write(*,'(A,ES13.6)')      ' k (=mδ/I)         : ', cfg%params%k
    write(*,'(A,ES13.6,A,ES13.6,A,ES13.6)') ' δ, ε, m           : ', &
        cfg%params%delta, ',', cfg%params%eps, ',', cfg%params%m
    write(*,'(A,ES13.6)')      ' I                 : ', cfg%params%inertia
    write(*,'(A,F12.6)')       ' ω0  = √(δ/I)      : ', real(w0)
    write(*,'(A,F12.6)')       ' ω1  (mode 1)      : ', real(ana%w1)
    write(*,'(A,F12.6)')       ' ω2  (mode 2)      : ', real(ana%w2)
    if (ana%w1 /= ana%w2) then
        tbeat = 2.0_qp*PI/(ana%w1 - ana%w2)
        write(*,'(A,F12.4,A)') ' T_beat = 2π/(ω1-ω2): ', real(tbeat), ' s'
    end if
    write(*,'(A,F9.5,A,ES10.3,A,I0,A)') ' Integration       :  t∈[0,', real(cfg%t_final), &
        '], h=', real(cfg%h), ', ', n, ' steps'
    write(*,'(A)') ''

    ! ---- integrate ------------------------------------------------------------
    call store%open_all(trim(data_dir))
    call run_simulation(cfg%params, u0, cfg%t0, cfg%h, n, store)
    call store%close_all()
    write(*,'(A)') ' Data written to '//trim(data_dir)//'/'

    ! ---- plots ----------------------------------------------------------------
    call plot_all(cfg%params, do_plot, trim(data_dir), trim(gpl_dir), trim(out_dir))
    if (do_plot) then
        write(*,'(A)') ' Figures rendered to '//trim(out_dir)//'/'
    else
        write(*,'(A)') ' Gnuplot scripts written to '//trim(gpl_dir)//'/  (run with --plot to render)'
    end if
    write(*,'(A)') ''
    write(*,'(A)') ' Done (!)'
    write(*,'(A)') ''

contains

    subroutine read_arg(index, value)
        integer, intent(in) :: index
        character(len=:), allocatable, intent(out) :: value
        integer :: length, status
        call get_command_argument(index,length=length,status=status)
        if (status/=0) error stop 'could not read CLI argument length'
        allocate(character(len=length) :: value)
        call get_command_argument(index,value,status=status)
        if (status/=0) error stop 'could not read CLI argument'
    end subroutine read_arg

    subroutine print_usage()
        write(*,'(A)') 'Usage: wilberforce [PRESET] [--plot] [--outdir DIR]'
        write(*,'(A)') '  PRESET    1..14  (default 14)'
        write(*,'(A)') '  --plot    render PNGs with Gnuplot after writing data'
        write(*,'(A)') '  --outdir  output parent directory (default: generated)'
    end subroutine print_usage

end program wilberforce_main
