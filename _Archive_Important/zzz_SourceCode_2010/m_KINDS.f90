! CVS info:
! $RCSfile: m_KINDS.F90,v $
! $Revision: 1.1.1.1 $
! $Date: 2008/03/20 08:07:58 $
!
module KINDS
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Selects the KIND for REAL data (compiler dependent)
!
! Revision record:
!-----------------------------------------------------------------------------------------------------------------------------------
! No.   Date          Programmer     Description of changes
!-----------------------------------------------------------------------------------------------------------------------------------
! 0     2002/07/15    E Z Muller     Original
!
  implicit none
  public
! for IEEE 754, Intel, Lahey
!  integer ,parameter :: sp = selected_real_kind(p = 6 ,r = 37)  ! IEEE 754 single-precision
!  integer ,parameter :: dp = selected_real_kind(p = 15,r = 307) ! IEEE 754 double-precision
  integer ,parameter :: sp = selected_real_kind(p = 6 ,r = 37)            ! single-precision KIND
  integer ,parameter :: dp = selected_real_kind(p = 13,r = 200)           ! double-precision KIND
  integer, parameter :: qp = MAX(dp, selected_real_kind(p = 26,r = 300))  ! quadruple-precision kind
! pseudonyms
  integer, parameter :: singel    = sp  ! single-precision
  integer, parameter :: doubel    = dp  ! double-precision
  integer, parameter :: quadrupel = qp  ! quadruple-precision
  integer, parameter :: reals     = sp  ! standard precision
  integer, parameter :: high_p    = dp  ! high precision
!
contains
  subroutine test_kinds
    real(KIND=sp) :: ps=0.0_sp
    real(KIND=dp) :: pd=0.0_dp
    real(KIND=qp) :: pq=0.0_qp
    write(*,'(" sp has precision",i3," and range",i5)') precision(ps), range(ps)
    write(*,'(" sp has machine epsilon",es14.6)') epsilon(ps)
    write(*,'(" dp has precision",i3," and range",i5)') precision(pd), range(pd)
    write(*,'(" dp has machine epsilon",es14.6)') epsilon(pd)
    write(*,'(" qp has precision",i3," and range",i5)') precision(pq), range(pq)
    write(*,'(" qp has machine epsilon",es14.6)') epsilon(pq)
  end subroutine test_kinds
end module KINDS
