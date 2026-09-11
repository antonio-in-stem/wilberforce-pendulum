!> \file wilberforce_kinds.f90
!! \brief Real-kind parameters used throughout the package.
!!
!! Centralising the working precisions in one module means the numerical
!! precision of the whole solver can be inspected — or changed — in a single
!! place, which is the standard practice in maintained scientific Fortran.
module wilberforce_kinds
    implicit none
    private

    !> Double precision (~15 significant digits). Used for the physical
    !! parameters, exactly as in the original program (`double precision`).
    integer, parameter, public :: dp = selected_real_kind(15, 307)

    !> Quadruple precision (~33 significant digits). The dynamical state is
    !! integrated in this kind, reproducing the original `real*16` declarations
    !! so that round-off stays far below the RK4 truncation error.
    integer, parameter, public :: qp = selected_real_kind(30)

    !> Pi in quadruple precision (matches 4*atan(1.0q0)).
    real(qp), parameter, public :: PI = 4.0_qp * atan(1.0_qp)

end module wilberforce_kinds
