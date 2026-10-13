Module FileHandle

! Module for central control of file opening / closing / error control and unit number control

IMPLICIT NONE

PRIVATE
PUBLIC                                          :: OpenFile, CloseFile


integer        :: initial_Unit = 200       ! starting unit number for file control
integer        :: fileCount = 0            ! number of files in unit
integer        :: internal_unit=199
integer        :: unitlimit=10000



!---------------------------------------------------------------------------
CONTAINS



!--------------------------------------------------------------------
Function OpenFile (FileName, FileStatus,FileAction)
! Subroutine for generic file opening and unit number allocation
!--------------------------------------------------------------------
character (*), INTENT (IN)     :: FileName              ! name of file to open
character (*), INTENT (IN)     :: FileStatus,FileAction ! Flags for file-open status/action
integer                        :: ierror,i                ! status flag for file opening
integer                        :: OpenFile,cur          ! File Unit to return
LOGICAL                        :: legal,used
character(len=80)              :: fln


   if (fileCount.eq.0) then
      OPEN (UNIT=internal_unit, FILE='filelog.txt', STATUS='REPLACE', ACTION='WRITE',IOSTAT=ierror)
   end if
   do i=initial_Unit-1,unitlimit
     inquire(UNIT=i,EXIST=legal,NAME=fln,OPENED=used)
     if (legal .and. .not. used) then
      cur=i
      exit
!     else
!       write(internal_unit,*) 'Unit ',i,' not available - used by ',fln

     end if
   end do
   !Attempt to open the target file
   OPEN (UNIT=cur, FILE=trim(adjustl(FileName)), STATUS=trim(adjustl(FileStatus)), ACTION=trim(adjustl(FileAction)),IOSTAT=ierror)
   if (ierror.eq.0) then

      OpenFile=cur
      fileCount=fileCount+1
      write(internal_unit,*) 'File ',trim(adjustl(FileName)),' opened on unit ',cur
   else
      write (*,*) 'Error opening file ',trim(adjustl(FileName)),' on unit ',initial_Unit+fileCount,' :',ierror
      OpenFile=-1
   end if



End FUNCTION OpenFile


!--------------------------------------------------------------
Subroutine CloseFile(UnitNumber)
! Generic Routine to close file on a given unit number
!--------------------------------------------------------------
integer,INTENT(IN)                 :: UnitNumber           ! Unit Number of file to close

   write(internal_unit,*) 'File closed on unit ',UnitNumber

close (unitnumber)

End Subroutine CloseFile


!--------------------------------------------------------------------------------------------
END Module FileHandle
