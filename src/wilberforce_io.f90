!> \file wilberforce_io.f90
!! \brief Portable directory creation and time-series data output.
!!
!! Replaces the original Windows-only `folder.bat` trick (which wrote a `.bat`
!! file and shelled out to it) with a single portable routine that works under
!! both Windows `cmd` and POSIX shells. The column layout of every output file
!! is preserved exactly so downstream tooling (Gnuplot scripts, Python) is
!! unaffected.
module wilberforce_io
    use wilberforce_kinds, only: qp
    use wilberforce_types, only: state_t
    use wilberforce_energy, only: energy_t
    implicit none
    private

    public :: datastore_t, ensure_dir, validate_path

    !> Bundles the six output units and the directory they live in.
    type :: datastore_t
        integer :: u_pos_long = 0   !< posicion_longitudinal.txt : t, z_num, z_ana
        integer :: u_vel_long = 0   !< velocidad_longitudinal.txt: t, v_num, v_ana
        integer :: u_pos_ang  = 0   !< posicion_angular.txt      : t, th_num, th_ana
        integer :: u_vel_ang  = 0   !< velocidad_angular.txt     : t, om_num, om_ana
        integer :: u_phase    = 0   !< fase.txt      : th_num, z_num, th_ana, z_ana
        integer :: u_energy   = 0   !< energias.txt  : t, E_long, E_ang, E_coup, E_tot
        character(len=:), allocatable :: dir
    contains
        procedure :: open_all
        procedure :: write_record
        procedure :: close_all
    end type datastore_t

contains

    !> Paths enter both an OS shell and Gnuplot. Reject their metacharacters;
    !! use forward slashes on Windows. Spaces, Unicode and drive colons are allowed.
    subroutine validate_path(path)
        character(len=*), intent(in) :: path
        integer :: i, c
        if (len_trim(path)==0) error stop 'empty output path'
        if (path(1:1)=='-') error stop 'output path cannot start with a dash'
        do i=1,len_trim(path)
            c=iachar(path(i:i))
            if (c<32 .or. c==127 .or. c==34 .or. c==39 .or. c==36 .or. c==96 .or. &
                c==37 .or. c==33 .or. c==38 .or. c==124 .or. c==60 .or. c==62 .or. &
                c==94 .or. c==59 .or. c==92) error stop 'unsafe output path character'
        end do
    end subroutine validate_path

    subroutine ensure_dir(path)
        character(len=*), intent(in) :: path
        character(len=512) :: osname
        character(len=:), allocatable :: cmd, winpath
        integer :: length, stat, cstat, i
        call validate_path(path)
        osname=''
        call get_environment_variable('OS', osname, length, status=stat)
        if (stat==0 .and. index(osname, 'Windows')>0) then
            winpath=path
            do i=1,len(winpath)
                if (winpath(i:i)=='/') winpath(i:i)=achar(92)
            end do
            cmd='if not exist "'//winpath//'" mkdir "'//winpath//'"'
        else
            cmd='mkdir -p -- "'//path//'"'
        end if
        call execute_command_line(cmd,wait=.true.,exitstat=stat,cmdstat=cstat)
        if (cstat/=0) error stop 'could not start directory creation command'
        if (stat/=0) error stop 'could not create output directory'
    end subroutine ensure_dir

    !> Open all six data files inside `dir` (default "generated/data").
    subroutine open_all(self, dir)
        class(datastore_t), intent(inout) :: self
        character(len=*), intent(in), optional :: dir
        call self%close_all()
        if (present(dir)) then
            self%dir = dir
        else
            self%dir = 'generated/data'
        end if
        call ensure_dir(self%dir)
        open(newunit=self%u_pos_long, file=self%dir//'/posicion_longitudinal.txt', status='replace', action='write')
        open(newunit=self%u_vel_long, file=self%dir//'/velocidad_longitudinal.txt', status='replace', action='write')
        open(newunit=self%u_pos_ang, file=self%dir//'/posicion_angular.txt', status='replace', action='write')
        open(newunit=self%u_vel_ang, file=self%dir//'/velocidad_angular.txt', status='replace', action='write')
        open(newunit=self%u_phase, file=self%dir//'/fase.txt', status='replace', action='write')
        open(newunit=self%u_energy, file=self%dir//'/energias.txt', status='replace', action='write')
    end subroutine open_all

    !> Write one time step to every file. The column order matches the original
    !! program exactly (list-directed output).
    subroutine write_record(self, t, num, ana, en)
        class(datastore_t), intent(in) :: self
        real(qp),           intent(in) :: t
        type(state_t),      intent(in) :: num  !< numerical (RK4) state
        type(state_t),      intent(in) :: ana  !< analytic state
        type(energy_t),     intent(in) :: en

        if (any([self%u_pos_long,self%u_vel_long,self%u_pos_ang,self%u_vel_ang, &
                 self%u_phase,self%u_energy]==0)) error stop 'datastore is not open'
        write(self%u_pos_long, *) t, num%z,     ana%z
        write(self%u_vel_long, *) t, num%v,     ana%v
        write(self%u_pos_ang,  *) t, num%theta, ana%theta
        write(self%u_vel_ang,  *) t, num%omega, ana%omega
        write(self%u_phase,    *) num%theta, num%z, ana%theta, ana%z
        write(self%u_energy,   *) t, en%longitudinal(), en%angular(), &
                                  en%coupling, en%total()
    end subroutine write_record

    !> Close only owned connections; safe before open, after close, and on reopen.
    subroutine close_all(self)
        class(datastore_t), intent(inout) :: self
        if (self%u_pos_long/=0) close(self%u_pos_long)
        self%u_pos_long=0
        if (self%u_vel_long/=0) close(self%u_vel_long)
        self%u_vel_long=0
        if (self%u_pos_ang/=0) close(self%u_pos_ang)
        self%u_pos_ang=0
        if (self%u_vel_ang/=0) close(self%u_vel_ang)
        self%u_vel_ang=0
        if (self%u_phase/=0) close(self%u_phase)
        self%u_phase=0
        if (self%u_energy/=0) close(self%u_energy)
        self%u_energy=0
    end subroutine close_all

end module wilberforce_io
