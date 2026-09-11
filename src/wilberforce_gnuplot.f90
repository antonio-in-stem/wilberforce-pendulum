!> \file wilberforce_gnuplot.f90
!! \brief Generation of the Gnuplot scripts, preserved from the original program.
!!
!! The original code emitted seven near-identical Gnuplot scripts, each 90+ lines,
!! by copy-paste. This module keeps the *output* identical (same titles, colours,
!! parameter side-panel, and plotted columns) but removes the duplication: the
!! shared header and parameter-label block are written once, and each figure only
!! specifies what is unique to it.
!!
!! Two behavioural improvements over the original, both opt-in and harmless:
!!   * scripts run head-less with `gnuplot <script>` (writing PNGs) instead of
!!     `start gnuplot -p <script>` (which pops seven windows);
!!   * the data/script/output directories are parameters, not hard-coded.
module wilberforce_gnuplot
    use wilberforce_types, only: params_t
    use wilberforce_io,    only: ensure_dir, validate_path
    implicit none
    private

    public :: generate_all_scripts, plot_all

    integer, parameter :: GPL = 20241   !< scratch unit for writing scripts

contains

    !> Format a double as the original `strd_mod` did: F0.4, with a leading '0'
    !! prepended when |x|<1 so that ".5000" prints as "0.5000".
    function strd(x) result(s)
        double precision, intent(in) :: x
        character(len=50) :: s
        write(s, '(F0.4)') x
        if (s(1:1)=='.') s='0'//trim(s)
        if (s(1:2)=='-.') s='-0'//trim(s(2:))
    end function strd

    !> Common preamble: terminal, output file, title, axis labels, key.
    subroutine write_header(unit, outfile, title, xlabel, ylabel)
        integer,          intent(in) :: unit
        character(len=*), intent(in) :: outfile, title, xlabel, ylabel
        write(unit,'(A)') "# Configuración del gráfico"
        write(unit,'(A)') "set terminal pngcairo size 2000,1200 enhanced font 'Verdana,22'"
        write(unit,'(A)') "set output '"//outfile//"'"
        write(unit,'(A)')
        write(unit,'(A)') "# Títulos y etiquetas"
        write(unit,'(A)') "set title '"//title//"'"
        write(unit,'(A)') "set xlabel '"//xlabel//"'"
        write(unit,'(A)') "set ylabel '"//ylabel//"'"
        write(unit,'(A)') "set key outside bottom right"
        write(unit,'(A)')
    end subroutine write_header

    !> The parameter side-panel (k, δ, ε, m, I, z0, θ0), identical to the original.
    subroutine write_param_labels(unit, p)
        integer,        intent(in) :: unit
        type(params_t), intent(in) :: p
        write(unit,'(A)') "set label 1 'Parámetros' at screen 0.88,0.9 left font 'Verdana,24' front"
        write(unit,'(A)') "set label 2 'k="//trim(adjustl(strd(p%k)))//"'&
            & at screen 0.88,0.85 left font 'Verdana,22' front"
        write(unit,'(A)') "set label 3 'δ="//trim(adjustl(strd(p%delta)))//"'&
            & at screen 0.88,0.8 left font 'Verdana,22' front"
        write(unit,'(A)') "set label 4 'ε="//trim(adjustl(strd(p%eps)))//"'&
            & at screen 0.88,0.75 left font 'Verdana,22' front"
        write(unit,'(A)') "set label 5 'm="//trim(adjustl(strd(p%m)))//"'&
            & at screen 0.88,0.7 left font 'Verdana,22' front"
        write(unit,'(A)') "set label 6 'I="//trim(adjustl(strd(p%inertia)))//"'&
            & at screen 0.88,0.65 left font 'Verdana,22' front"
        write(unit,'(A)') "set label 7 'z(0)="//trim(adjustl(strd(p%z0)))//"'&
            & at screen 0.88,0.6 left font 'Verdana,22' front"
        write(unit,'(A)') "set label 8 'θ(0)="//trim(adjustl(strd(p%theta0)))//"'&
            & at screen 0.88,0.55 left font 'Verdana,22' front"
        write(unit,'(A)') "show label"
        write(unit,'(A)')
        write(unit,'(A)') "# Adecuar los limites"
        write(unit,'(A)') "set grid lw 4 lc rgb '#808080'"
        write(unit,'(A)')
        write(unit,'(A)') "# Graficar datos"
    end subroutine write_param_labels

    !> Write one complete "position + velocity vs time" script.
    subroutine write_series_script(fname, outfile, title, ylabel, p, dat_pos, dat_vel, &
                                   col, c_pos, c_vel)
        character(len=*), intent(in) :: fname, outfile, title, ylabel
        type(params_t),   intent(in) :: p
        character(len=*), intent(in) :: dat_pos, dat_vel   !< data files
        integer,          intent(in) :: col                !< column (2=numeric, 3=analytic)
        character(len=*), intent(in) :: c_pos, c_vel       !< line colours
        character(len=1) :: cc
        write(cc,'(I1)') col
        open(GPL, file=fname)
        call write_header(GPL, outfile, title, 'Tiempo [s]', ylabel)
        call write_param_labels(GPL, p)
        write(GPL,'(A)') "plot 0 with lines lw 2 lt rgb '#FFFFFF00' notitle, \"
        write(GPL,'(A)') "     '"//dat_pos//"' u 1:"//cc//" w lp&
            & ps 0.5 pt 7 lc rgb'"//c_pos//"' title 'Posición', \"
        write(GPL,'(A)') "     '"//dat_vel//"' u 1:"//cc//" w lp&
            & ps 0.5 pt 7 lc rgb'"//c_vel//"' title 'Velocidad', \"
        write(GPL,'(A)') "     0 with lines lw 2 lt rgb 'black' notitle behind"
        close(GPL)
    end subroutine write_series_script

    !> Generate all seven Gnuplot scripts for a parameter set.
    subroutine generate_all_scripts(p, data_dir, gpl_dir, out_dir)
        type(params_t),   intent(in) :: p
        character(len=*), intent(in) :: data_dir, gpl_dir, out_dir
        character(len=:), allocatable :: d, g, o

        call validate_path(data_dir)
        call ensure_dir(gpl_dir)
        call ensure_dir(out_dir)
        d = data_dir; g = gpl_dir; o = out_dir

        ! --- longitudinal (numeric col 2, analytic col 3) ---
        call write_series_script(g//'/longitudinal.gpl', o//'/numerical_longitudinal.png', &
            '[Método Numérico] Longitudinal vs Tiempo', 'Posición[m] Velocidad[m/s]', p, &
            d//'/posicion_longitudinal.txt', d//'/velocidad_longitudinal.txt', 2, '#4778FF', '#A8BFFF')
        call write_series_script(g//'/analytic_longitudinal.gpl', o//'/analytic_longitudinal.png', &
            '[Método Analítico] Longitudinal vs Tiempo', 'Posición[m] Velocidad[m/s]', p, &
            d//'/posicion_longitudinal.txt', d//'/velocidad_longitudinal.txt', 3, '#4A4954', '#B0ADC7')

        ! --- angular ---
        call write_series_script(g//'/angular.gpl', o//'/numerical_angular.png', &
            '[Método Numérico] Angular vs Tiempo', 'Posición[rad] Velocidad[rad/s]', p, &
            d//'/posicion_angular.txt', d//'/velocidad_angular.txt', 2, '#7645FF', '#BDA8FF')
        call write_series_script(g//'/analytic_angular.gpl', o//'/analytic_angular.png', &
            '[Método Analítico] Ángular vs Tiempo', 'Posición[rad] Velocidad[rad/s]', p, &
            d//'/posicion_angular.txt', d//'/velocidad_angular.txt', 3, '#4A4954', '#B0ADC7')

        ! --- phase portraits ---
        call write_phase_script(g//'/fase.gpl', o//'/numerical_fase.png', &
            '[Método Numérico] Longitudinal vs Angular', p, d//'/fase.txt', 1, 2, '#C7329A')
        call write_phase_script(g//'/analytic_fase.gpl', o//'/analytic_fase.png', &
            '[Método Analítico] Longitudinal vs Angular', p, d//'/fase.txt', 3, 4, '#B0ADC7')

        ! --- energy budget ---
        call write_energy_script(g//'/energias.gpl', o//'/numerical_energias.png', p, d//'/energias.txt')
    end subroutine generate_all_scripts

    !> Phase-portrait script (angular vs longitudinal).
    subroutine write_phase_script(fname, outfile, title, p, datfile, cx, cy, colour)
        character(len=*), intent(in) :: fname, outfile, title, datfile, colour
        type(params_t),   intent(in) :: p
        integer,          intent(in) :: cx, cy
        character(len=1) :: sx, sy
        write(sx,'(I1)') cx;  write(sy,'(I1)') cy
        open(GPL, file=fname)
        call write_header(GPL, outfile, title, 'Posición Angular [rad]', 'Posición longitudinal [m]')
        call write_param_labels(GPL, p)
        write(GPL,'(A)') "plot 0 with lines lw 2 lt rgb '#FFFFFF00' notitle, \"
        write(GPL,'(A)') "     '"//datfile//"' u "//sx//":"//sy//" w lp&
            & ps 0.5 pt 7 lc rgb'"//colour//"' title 'Fase', \"
        write(GPL,'(A)') "     0 with lines lw 2 lt rgb 'black' notitle behind"
        close(GPL)
    end subroutine write_phase_script

    !> Energy-budget script (four contributions vs time).
    subroutine write_energy_script(fname, outfile, p, datfile)
        character(len=*), intent(in) :: fname, outfile, datfile
        type(params_t),   intent(in) :: p
        open(GPL, file=fname)
        call write_header(GPL, outfile, '[Método Numérico] Energia vs Tiempo', 'Tiempo [s]', 'Energía[J]')
        call write_param_labels(GPL, p)
        write(GPL,'(A)') "plot 0 with lines lw 2 lt rgb '#FFFFFF00' notitle, \"
        write(GPL,'(A)') "     '"//datfile//"' u 1:2 w lp ps 0.5 pt 7 lc rgb'#4778FF' title 'Longitudinal', \"
        write(GPL,'(A)') "     '"//datfile//"' u 1:3 w lp ps 0.5 pt 7 lc rgb'#7645FF' title 'Angular', \"
        write(GPL,'(A)') "     '"//datfile//"' u 1:4 w lp ps 0.5 pt 7 lc rgb'#4E59BA' title 'Acoplamiento', \"
        write(GPL,'(A)') "     '"//datfile//"' u 1:5 w lp ps 0.5 pt 7 lc rgb'#1CD9FF' title 'Energía Total', \"
        write(GPL,'(A)') "     0 with lines lw 2 lt rgb 'black' notitle behind"
        close(GPL)
    end subroutine write_energy_script

    !> Generate the scripts and (optionally) run Gnuplot on each to produce PNGs.
    subroutine plot_all(p, run, data_dir, gpl_dir, out_dir)
        type(params_t),   intent(in) :: p
        logical,          intent(in) :: run
        character(len=*), intent(in), optional :: data_dir, gpl_dir, out_dir
        character(len=:), allocatable :: d, g, o
        integer :: stat, cstat
        character(len=*), parameter :: names(7) = [character(len=24) :: &
            'longitudinal', 'analytic_longitudinal', 'angular', 'analytic_angular', &
            'fase', 'analytic_fase', 'energias']
        integer :: i

        d = 'generated/data';     if (present(data_dir)) d = data_dir
        g = 'generated/gpl';      if (present(gpl_dir))  g = gpl_dir
        o = 'generated/graficas'; if (present(out_dir))  o = out_dir

        call generate_all_scripts(p, d, g, o)
        if (run) then
            do i = 1, size(names)
                call execute_command_line('gnuplot "'//g//'/'//trim(names(i))//'.gpl"', &
                                          wait=.true., exitstat=stat, cmdstat=cstat)
                if (cstat/=0) error stop 'could not start Gnuplot'
                if (stat/=0) error stop 'Gnuplot failed'
            end do
        end if
    end subroutine plot_all

end module wilberforce_gnuplot
