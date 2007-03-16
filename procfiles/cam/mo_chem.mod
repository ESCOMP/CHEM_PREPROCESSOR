
      module chem_mods
!--------------------------------------------------------------
!     	... Basic chemistry parameters and arrays
!--------------------------------------------------------------

      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

      save

      integer, parameter :: hetcnt    = HETCNT, &    ! number of heterogeneous processes
                            phtcnt    = PHTCNT, &    ! number of photolysis reactions
                            rxntot    = RXNCNT, &    ! number of total reactions
                            gascnt    = GASCNT, &    ! number of gas phase reactions
                            nabscol   = NCOL, &      ! number of absorbing column densities
                            gas_pcnst = PCNST, &     ! number of "gas phase" species
                            nfs       = NFS, &       ! number of "fixed" species
                            relcnt    = RELCNT, &    ! number of relationship species
                            grpcnt    = GRPCNT, &    ! number of group members
                            nzcnt     = IMP_NZCNT, & ! number of non-zero matrix entries
                            extcnt    = EXTCNT, &    ! number of species with external forcing
                            clscnt1   = CLSCNT1, &   ! number of species in explicit class
                            clscnt2   = CLSCNT2, &   ! number of species in hov class
                            clscnt3   = CLSCNT3, &   ! number of species in ebi class
                            clscnt4   = CLSCNT4, &   ! number of species in implicit class
                            clscnt5   = CLSCNT5, &   ! number of species in rodas class
                            indexm    = INDEXM, &    ! index of total atm density in invariant array
                            indexh2o  = INDEXH2O, &  ! index of water vapor density
                            clsze      = CLSZE, &       ! loop length for implicit chemistry
                            rxt_alias_cnt = ALIASCNT

      integer   :: clscnt(5)
      integer   :: cls_rxt_cnt(4,5)
      integer   :: clsmap(gas_pcnst,5)
      integer   :: permute(gas_pcnst,5)
# if CLSCNT4 != 0
      integer   :: diag_map(clscnt4)
# elif CLSCNT5 != 0
      integer   :: diag_map(clscnt5)
# endif
      real(r8)  :: adv_mass(gas_pcnst)
      real(r8)  :: fix_mass(max(1,nfs))
# if GRPCNT != 0
      real(r8)  :: nadv_mass(grpcnt)
# endif

      character(len=8)               :: het_lst(max(1,hetcnt))
      character(len=8)               :: inv_lst(max(1,nfs))
      character(len=8)               :: extfrc_lst(max(1,extcnt))
      logical                        :: frc_from_dataset(max(1,extcnt))
      integer,           allocatable :: rxt_alias_map(:)
      character(len=16), allocatable :: rxt_alias_lst(:)

      end module chem_mods
