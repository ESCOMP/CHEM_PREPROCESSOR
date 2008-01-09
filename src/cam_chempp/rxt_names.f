
      subroutine make_rxt_name_mod
!--------------------------------------------------------------------------------
!	... Makes a module of parameter reaction names
!--------------------------------------------------------------------------------

      use rxt_mod, only : rxtcnt => rxntot, gascnt, phtcnt, rxt_alias, rxt_has_alias
      use io,      only : temp_path

      implicit none

!--------------------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------------------
      integer :: i
      character(len=80) :: buff
      character(len=5)  :: num
      logical :: lexist

!--------------------------------------------------------------------------------
!	... Check mod file existence; remove if found
!--------------------------------------------------------------------------------
      inquire( file = trim( temp_path ) // 'rxt_names.mod', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'rxt_names.mod' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'rxt_names.mod' )
             
      buff = ''
      write(30,'(a)') buff
      buff(7:) = 'module m_rxt_id'
      write(30,'(a)') buff
      buff = ''
      write(30,'(a)') buff
      buff(7:) = 'implicit none'
      write(30,'(a)') buff
      buff = ''
      write(30,'(a)') buff

      do i = 1,rxtcnt
         if( rxt_alias(i) /= ' ' ) then
            rxt_has_alias(i) = .true.
            write(buff(7:),'(''integer, parameter :: rid_'',a,1x,''='',1x,i4)') &
                             rxt_alias(i)(:len_trim(rxt_alias(i))), i
            write(30,'(a)') buff
         end if
      end do

      if( any( rxt_has_alias(:rxtcnt) ) ) then
         buff = ''
         write(30,'(a)') buff
      end if

      do i = 1,rxtcnt
         if( .not. rxt_has_alias(i) ) then
            write(num,'(i5)') i+10000
            if( i <= phtcnt ) then
               write(buff(7:),'(''integer, parameter :: rid_j'',a,1x,''='',1x,i4)') &
                   num(2:5), i
               write(rxt_alias(i)(:5),'(''j'',a)') num(2:5)
            else
               write(buff(7:),'(''integer, parameter :: rid_r'',a,1x,''='',1x,i4)') &
                        num(2:5), i
               write(rxt_alias(i)(:5),'(''r'',a)') num(2:5)
            end if
            write(30,'(a)') buff
         end if
      end do

      buff = ''
      write(30,'(a)') buff
      buff(7:) = 'end module m_rxt_id'
      write(30,'(a)') buff
      close( 30 )

      end subroutine make_rxt_name_mod
