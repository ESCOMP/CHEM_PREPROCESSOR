
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,   only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass
      use chem_mods,   only : het_lst, extfrc_lst, inv_lst
      use mo_tracname, only  : solsym
      use chem_mods,   only : frc_from_dataset
      use chem_mods,   only : rxt_alias_map, rxt_alias_lst, rxt_alias_cnt
      use abortutils,  only : endrun
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios
      clscnt(:) = (/             1 ,            0 ,            0 ,            0 ,            0 /)

      cls_rxt_cnt(:,1) = (/             0 ,            1 ,            0 ,            0  /)

      solsym(:  1) = (/ 'CO      ' /)

      adv_mass(:  1) = (/ 28.0104008_r8 /)

      fix_mass(:  5) = (/ 0.00000000_r8, 28.0134792_r8, 31.9988003_r8, 18.0142002_r8, 17.0067997_r8 /)

      clsmap(:  1,1) = (/    1 /)



      extfrc_lst(:  1) = (/ 'CO      ' /)

      frc_from_dataset(:  1) = (/ .true. /)

      inv_lst(:  5) = (/ 'M       ', 'N2      ', 'O2      ', 'H2O     ', 'OH      ' /)

      rxt_alias_cnt =    1
      if( allocated( rxt_alias_lst ) ) then
         deallocate( rxt_alias_lst )
      end if
      allocate( rxt_alias_lst(rxt_alias_cnt),stat=ios )
      if( ios /= 0 ) then
         write(*,*) 'set_sim_dat: failed to allocate rxt_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_alias_map ) ) then
         deallocate( rxt_alias_map )
      end if
      allocate( rxt_alias_map(rxt_alias_cnt),stat=ios )
      if( ios /= 0 ) then
         write(*,*) 'set_sim_dat: failed to allocate rxt_alias_map; error = ',ios
         call endrun
      end if
      rxt_alias_lst(:rxt_alias_cnt) = (/ 'usr8            ' /)
      rxt_alias_map(:rxt_alias_cnt) = (/    1 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
