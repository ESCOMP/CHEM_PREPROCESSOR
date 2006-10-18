      
      subroutine make_sim_dat( model, sparse )
!-------------------------------------------------------------------
!	... write the simulation data routine; only for CAM
!-------------------------------------------------------------------

      use io,      only : temp_path
      use sp_mods, only : sparsity
      use var_mod, only : clscnt, clsmap, permute, new_nq, new_solsym
      use var_mod, only : nq, newind, mass, temp_mass
      use rxt_mod, only : cls_rxt_cnt
      use rxt_mod, only : usrcnt, usrmap, frc_from_dataset
      use rxt_mod, only : rxt_has_alias, rxt_alias
!!      use var_mod, only : dvel_cnt, dvel_map
      use var_mod, only : nfs, fixsym
      use rxt_mod

      implicit none

!-------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------
      character(len=8), intent(in) :: model
      type(sparsity), intent(in)   :: sparse(2)

!-------------------------------------------------------------------
!	... local variables
!-------------------------------------------------------------------
      integer, parameter :: max_len= 128

      character(len=8) ::  wrk_chr(5)
      integer, allocatable   ::  ndx(:)
      integer  :: m1, l
      integer  ::  i, m, n, n1
      integer  ::  lpos
      integer  ::  lstrt
      character(len=max_len) :: line
      character(len=64)      :: frmt
      character(len=24)      :: number
      logical  ::  flush
      logical  ::  lexist

      inquire( file = trim( temp_path ) // 'mo_sim_dat.F', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'mo_sim_dat.F' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'mo_sim_dat.F' )

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'module mo_sim_dat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'private'
      write(30,100) trim(line)
      line(7:) = 'public :: set_sim_dat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'contains'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'subroutine set_sim_dat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'use chem_mods,   only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass'
      write(30,100) trim(line)
      if( clscnt(4) > 0 ) then
        line(7:) = 'use chem_mods,   only : diag_map'
        write(30,100) trim(line)
      endif
      line(7:) = 'use chem_mods,   only : het_lst, extfrc_lst, inv_lst'
      write(30,100) trim(line)
      line(7:) = 'use mo_tracname, only  : solsym'
      write(30,100) trim(line)
      line(7:) = 'use chem_mods,   only : frc_from_dataset'
      write(30,100) trim(line)
      line(7:) = 'use chem_mods,   only : rxt_alias_map, rxt_alias_lst, rxt_alias_cnt'
      write(30,100) trim(line)
!      line(7:) = 'use chem_mods,   only : drydep_lst, drydep_cnt'
!      write(30,100) trim(line)
      line(7:) = 'use abortutils,  only : endrun'
      write(30,100) trim(line)
      line(7:) = 'use shr_kind_mod, only : r8 => shr_kind_r8'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'implicit none '
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!      ... local variables'
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------'
      write(30,100) trim(line)
      line = '      integer :: ios'
      write(30,100) trim(line)
!-------------------------------------------------------------------
!	... set the simulation chemical mechanism data
!	    class species count
!-------------------------------------------------------------------
      line = '      clscnt(:) = (/'
      write(line(len_trim(line)+2:),*) clscnt(1),',',clscnt(2),',',clscnt(3),',',clscnt(4),',',clscnt(5)
      line(len_trim(line)+2:) = '/)'
      write(30,'(a)') trim(line)
      line = ' '
      write(30,100) trim(line)

!-------------------------------------------------------------------
!	... class reaction count
!-------------------------------------------------------------------
      do i = 1,5
         if( clscnt(i) > 0 ) then
            line = '      cls_rxt_cnt(:,'
            write(line(len_trim(line)+1:),'(i1,") = (/")') i
            m = len_trim(line) + 2
            write(line(m:),*) cls_rxt_cnt(1,i),',',cls_rxt_cnt(2,i),',',cls_rxt_cnt(3,i),',',cls_rxt_cnt(4,i),' /)'
            write(30,'(a)') trim(line)
         end if
      end do

!-------------------------------------------------------------------
!	... species symbols
!-------------------------------------------------------------------
      line = ' '
      write(30,100) trim(line)
      write(line,'("      solsym(:",i3,") = (/")') new_nq
      m = len_trim(line) + 2
      do n = 1,new_nq,5
         n1 = min( n+4,new_nq )
         if( n1 /= new_nq ) then
            write(line(m:),'(5("''",a8,"'',")," &")') new_solsym(n:n1)
         else
            if( n1 > n ) then
               write(frmt,'("(",i1)') n1 - n
               frmt(len_trim(frmt)+1:) = '("''",a8,"'',"),"''",a8,"'' /)")'
            else
               frmt = '("''",a8,"'' /)")'
            end if
            write(line(m:),trim(frmt)) new_solsym(n:n1)
         end if
         write(30,'(a)') trim(line)
         line = ' '
      end do

!-------------------------------------------------------------------
!	... species mass
!-------------------------------------------------------------------
      if( nq > 0 ) then
         line = ' '
         write(30,100) trim(line)
         temp_mass(:) = 0.
         do n = 1,nq
            if( newind(n) /= 0 ) then
               temp_mass(newind(n)) = mass(n)
            end if
         end do
         line = '      adv_mass(:'
         write(line(len_trim(line)+1:),'(i3,") = (/")') new_nq
         m     = len_trim(line) + 2
         lstrt = m
         do n = 1,new_nq
            number = ' '
            write(number,'(g15.9)') temp_mass(n)
            number = adjustl( number )
            lpos   = scan( number, '0123456789', back=.true. ) + 1
            if( n < new_nq ) then
               if( mod(n,5) /= 0 ) then
                  number(lpos:) = '_r8,'
                  flush = .false.
               else
                  number(lpos:) = '_r8, &'
                  flush = .true.
               end if
            else
               number(lpos:) = '_r8 /)'
               flush = .true.
            end if
            line(m:) = trim( number )
            if( .not. flush ) then
               m = len_trim(line) + 2
            else
               write(30,'(a)') trim(line)
               line = ' '
               m = lstrt
            end if
         end do
      end if
      ! fixed species masses
      if( nfs > 0 ) then
         line = ' '
         write(30,100) trim(line)
!!$         temp_mass(:) = 0.
!!$         do n = 1,nfs
!!$            if( newind(n) /= 0 ) then
!!$               temp_mass(newind(n)) = mass(nqs)
!!$            end if
!!$         end do
         temp_mass(newind(n)) = mass(nfs)

         line = '      fix_mass(:'
         write(line(len_trim(line)+1:),'(i3,") = (/")') nfs
         m     = len_trim(line) + 2
         lstrt = m

         do n = 1,nfs
            number = ' '
!!$            write(number,'(g15.9)') temp_mass(n+new_nq)
            write(number,'(g15.9)') mass(n+new_nq)
            number = adjustl( number )
            lpos   = scan( number, '0123456789', back=.true. ) + 1
            if( n < nfs ) then
               if( mod(n,5) /= 0 ) then
                  number(lpos:) = '_r8,'
                  flush = .false.
               else
                  number(lpos:) = '_r8, &'
                  flush = .true.
               end if
            else
               number(lpos:) = '_r8 /)'
               flush = .true.
            end if
            line(m:) = trim( number )
            if( .not. flush ) then
               m = len_trim(line) + 2
            else
               write(30,'(a)') trim(line)
               line = ' '
               m = lstrt
            end if
         end do
      end if

!-------------------------------------------------------------------
!	... class map
!-------------------------------------------------------------------
      line = ' '
      write(30,100) trim(line)
      do i = 1,5
         if( clscnt(i) > 0 ) then
            write(line,'("      clsmap(:",i3,",",i1,") = (/")') clscnt(i),i
            m = len_trim(line) + 2
            do n = 1,clscnt(i),10
               n1 = min( n+9,clscnt(i) )
               if( n1 /= clscnt(i) ) then
                  write(line(m:),'(10(i4,",")," &")') clsmap(n:n1,i,2)
               else
                  if( n1 > n ) then
                     write(frmt,'("(",i1)') n1 - n
                     frmt(len_trim(frmt)+1:) = '(i4,","),i4," /)")'
                  else
                     frmt = '(i4," /)")'
                  end if
                  write(line(m:),trim(frmt)) clsmap(n:n1,i,2)
               end if
               write(30,'(a)') trim(line)
               line = ' '
            end do
         end if
      end do

!-------------------------------------------------------------------
!	... class permutation map
!-------------------------------------------------------------------
      line = ' '
      write(30,100) trim(line)
      do i = 2,5
         if( clscnt(i) > 0 ) then
            write(line,'("      permute(:",i3,",",i1,") = (/")') clscnt(i),i
            m = len_trim(line) + 2
            do n = 1,clscnt(i),10
               n1 = min( n+9,clscnt(i) )
               if( n1 /= clscnt(i) ) then
                  write(line(m:),'(10(i4,",")," &")') permute(n:n1,i)
               else
                  if( n1 > n ) then
                     write(frmt,'("(",i1)') n1 - n
                     frmt(len_trim(frmt)+1:) = '(i4,","),i4," /)")'
                  else
                     frmt = '(i4," /)")'
                  end if
                  write(line(m:),trim(frmt)) permute(n:n1,i)
               end if
               write(30,'(a)') trim(line)
               line = ' '
            end do
         end if
      end do

!-------------------------------------------------------------------
!	... class diagonal indicies
!-------------------------------------------------------------------
      line = ' '
      write(30,100) trim(line)
      do i = 4,4
         if( clscnt(i) > 0 ) then
            write(line,'("      diag_map(:",i3,") = (/")') clscnt(i)
            m = len_trim(line) + 2
            do n = 1,clscnt(i),10
               n1 = min( n+9,clscnt(i) )
               if( n1 /= clscnt(i) ) then
                  write(line(m:),'(10(i4,",")," &")') sparse(i-3)%diag_map(n:n1)
               else
                  if( n1 > n ) then
                     write(frmt,'("(",i1)') n1 - n
                     frmt(len_trim(frmt)+1:) = '(i4,","),i4," /)")'
                  else
                     frmt = '(i4," /)")'
                  end if
                  write(line(m:),trim(frmt)) sparse(i-3)%diag_map(n:n1)
               end if
               write(30,'(a)') trim(line)
               line = ' '
            end do
         end if
      end do

!-----------------------------------------------------------------------
!        ... Write the wet removal species
!-----------------------------------------------------------------------
    if( hetcnt > 0 ) then
      line = ' '
      write(30,100) trim(line)
      write(line,'("      het_lst(:",i3,") = (/")') hetcnt
      m = len_trim(line) + 2
      do n = 1,hetcnt,5
	 wrk_chr(:) = ' '
         n1 = min( n+4,hetcnt )
	 do i = 1,n1-n+1 !!n,n1
	    wrk_chr(i) = new_solsym(hetmap(i+n-1,1))
         end do        
         if( n1 /= hetcnt ) then
            write(line(m:),'(5("''",a8,"'',")," &")') wrk_chr(1:n1-n+1)
         else
            if( n1 > n ) then
               write(frmt,'("(",i1)') n1 - n
               frmt(len_trim(frmt)+1:) = '("''",a8,"'',"),"''",a8,"'' /)")'
            else
               frmt = '("''",a8,"'' /)")'
            end if
            write(line(m:),trim(frmt)) wrk_chr(1:n1-n+1)
         end if
         write(30,'(a)') trim(line)
         line = ' '
      end do
    endif

!-----------------------------------------------------------------------
!        ... Write the ext frcing species
!-----------------------------------------------------------------------
   if( usrcnt > 0 ) then

      line = ' '
      write(30,100) trim(line)
      write(line,'("      extfrc_lst(:",i3,") = (/")') usrcnt
      m = len_trim(line) + 2
      do n = 1,usrcnt,5
	 wrk_chr(:) = ' '
         n1 = min( n+4,usrcnt )
	 do i = 1,n1-n+1 !!n,n1
	    wrk_chr(i) = new_solsym(usrmap(i+n-1))
         end do        
         if( n1 /= usrcnt ) then
            write(line(m:),'(5("''",a8,"'',")," &")') wrk_chr(1:n1-n+1)
         else
            if( n1 > n ) then
               write(frmt,'("(",i1)') n1 - n
               frmt(len_trim(frmt)+1:) = '("''",a8,"'',"),"''",a8,"'' /)")'
            else
               frmt = '("''",a8,"'' /)")'
            end if
            write(line(m:),trim(frmt)) wrk_chr(1:n1-n+1)
         end if
         write(30,'(a)') trim(line)
         line = ' '
      end do

         line = ' '
         write(30,100) trim(line)
         write(line,'("      frc_from_dataset(:",i3,") = (/")') usrcnt
         m1 = len_trim(line) + 2
         do n = 1,usrcnt,5
            n1 = min( n+4,usrcnt )
            m = m1
            do l = n,n1
               if( l /= usrcnt ) then
                  if( l /= n1 ) then
                     if( frc_from_dataset(l) ) then
                        write(line(m:),'(".true.,")')
                     else
                        write(line(m:),'(".false.,")')
                     end if
                  else
                     if( frc_from_dataset(l) ) then
                        write(line(m:),'(".true., &")')
                     else
                        write(line(m:),'(".false., &")')
                     end if
                  end if
               else
                  if( frc_from_dataset(l) ) then
                     write(line(m:),'(".true. /)")')
                  else
                     write(line(m:),'(".false. /)")')
                  end if
               end if
               m = len_trim(line) + 2
            end do
            write(30,'(a)') trim(line)
            line = ' '
         end do
    end if


!-------------------------------------------------------------------
!	... fixed species
!-------------------------------------------------------------------
      if( nfs > 0 ) then
         line = ' '
         write(30,100) trim(line)
         write(line,'("      inv_lst(:",i3,") = (/")') nfs
         m1 = len_trim(line) + 2
         do n = 1,nfs,5
            n1 = min( n+4,nfs )
            m = m1
            do l = n,n1
               if( l /= nfs ) then
                  if( l /= n1 ) then
                     write(line(m:),'("''",a8,"'',")') fixsym(l)
                  else
                     write(line(m:),'("''",a8,"'', &")') fixsym(l)
                  end if
               else
                  write(line(m:),'("''",a8,"'' /)")') fixsym(l)
               end if
               m = len_trim(line) + 2
            end do
            write(30,'(a)') trim(line)
            line = ' '
         end do
      end if


!-------------------------------------------------------------------
!	... reaction aliases
!-------------------------------------------------------------------
      i = count( rxt_has_alias(:rxntot) )
      if( i > 0 ) then
         allocate( ndx(i) )
         l = 0
         do m = 1,rxntot
            if( rxt_has_alias(m) ) then
               l = l + 1
               ndx(l) = m
            end if
         end do 
         line = ' '
         write(30,100) trim(line)
         write(line,'("      rxt_alias_cnt = ",i4)') l
         write(30,100) trim(line)
         line = ' '
         line(7:) = 'if( allocated( rxt_alias_lst ) ) then'
         write(30,100) trim(line)
         line(7:) = '   deallocate( rxt_alias_lst )'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line(7:) = 'allocate( rxt_alias_lst(rxt_alias_cnt),stat=ios )'
         write(30,100) trim(line)
         line(7:) = 'if( ios /= 0 ) then'
         write(30,100) trim(line)
         line = ' '
         line(10:) = 'write(*,*) ''set_sim_dat: failed to allocate rxt_alias_lst; error = '',ios'
         write(30,100) trim(line)
         line(10:) = 'call endrun'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line(7:) = 'if( allocated( rxt_alias_map ) ) then'
         write(30,100) trim(line)
         line(7:) = '   deallocate( rxt_alias_map )'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line(7:) = 'allocate( rxt_alias_map(rxt_alias_cnt),stat=ios )'
         write(30,100) trim(line)
         line(7:) = 'if( ios /= 0 ) then'
         write(30,100) trim(line)
         line = ' '
         line(10:) = 'write(*,*) ''set_sim_dat: failed to allocate rxt_alias_map; error = '',ios'
         write(30,100) trim(line)
         line(10:) = 'call endrun'
         write(30,100) trim(line)
         line(7:) = 'end if'
         write(30,100) trim(line)
         line = '      rxt_alias_lst(:rxt_alias_cnt) = (/ '
         m1 = len_trim(line) + 2
         do n = 1,i,4
            n1 = min( n+3,i )
            m = m1
            do l = n,n1
               if( l /= i ) then
                  if( l /= n1 ) then
                     write(line(m:),'("''",a16,"'',")') rxt_alias(ndx(l))
                  else
                     write(line(m:),'("''",a16,"'', &")') rxt_alias(ndx(l))
                  end if
               else
                  write(line(m:),'("''",a16,"'' /)")') rxt_alias(ndx(l))
               end if
               m = len_trim(line) + 2
            end do
            write(30,'(a)') trim(line)
            line = ' '
         end do

         line = '      rxt_alias_map(:rxt_alias_cnt) = (/'
         m = len_trim(line) + 2
         do n = 1,i,10
            n1 = min( n+9,i )
            if( n1 /= i ) then
               write(line(m:),'(10(i4,",")," &")') ndx(n:n1)
            else
               if( n1 > n ) then
                  write(frmt,'("(",i1)') n1 - n
                  frmt(len_trim(frmt)+1:) = '(i4,","),i4," /)")'
               else
                  frmt = '(i4," /)")'
               end if
               write(line(m:),trim(frmt)) ndx(n:n1)
            end if
            write(30,'(a)') trim(line)
            line = ' '
         end do
         deallocate( ndx )
      end if


!!$!-------------------------------------------------------------------
!!$!	... dry deposition
!!$!-------------------------------------------------------------------
!!$      if( dvel_cnt > 0 ) then
!!$         line = ' '
!!$         write(30,100) trim(line)
!!$         write(line,'("      drydep_cnt = ",i4)') dvel_cnt
!!$         write(30,100) trim(line)
!!$         line = ' '
!!$         line(7:) = 'if( allocated( drydep_lst ) ) then'
!!$         write(30,100) trim(line)
!!$         line(7:) = '   deallocate( drydep_lst )'
!!$         write(30,100) trim(line)
!!$         line(7:) = 'end if'
!!$         write(30,100) trim(line)
!!$         line(7:) = 'allocate( drydep_lst(drydep_cnt),stat=ios )'
!!$         write(30,100) trim(line)
!!$         line(7:) = 'if( ios /= 0 ) then'
!!$         write(30,100) trim(line)
!!$         line = ' '
!!$         line(10:) = 'write(*,*) ''set_sim_dat: failed to allocate drydep_lst; error = '',ios'
!!$         write(30,100) trim(line)
!!$         line(10:) = 'call endrun'
!!$         write(30,100) trim(line)
!!$         line(7:) = 'end if'
!!$         write(30,100) trim(line)
!!$         line = '      drydep_lst(:drydep_cnt) = (/ '
!!$         m1 = len_trim(line) + 2
!!$         do n = 1,dvel_cnt,5
!!$            n1 = min( n+4,dvel_cnt )
!!$            m = m1
!!$            do l = n,n1
!!$               if( l /= dvel_cnt ) then
!!$                  if( l /= n1 ) then
!!$                     write(line(m:),'("''",a8,"'',")') new_solsym(dvel_map(l))
!!$                  else
!!$                     write(line(m:),'("''",a8,"'', &")') new_solsym(dvel_map(l))
!!$                  end if
!!$               else
!!$                  write(line(m:),'("''",a8,"'' /)")') new_solsym(dvel_map(l))
!!$               end if
!!$               m = len_trim(line) + 2
!!$            end do
!!$            write(30,'(a)') trim(line)
!!$            line = ' '
!!$         end do
!!$      end if

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end subroutine set_sim_dat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end module mo_sim_dat'
      write(30,100) trim(line)

100   format(a)

      end subroutine make_sim_dat
