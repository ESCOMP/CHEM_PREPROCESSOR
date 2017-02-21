
module mo_imp_sol

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : clscnt4, gas_pcnst, clsmap
  use cam_logfile,  only : iulog
  use ppgrid,       only : veclen

  implicit none
  private
  public :: imp_slv_inti, imp_sol
  save

  real(r8), parameter :: rel_err = 1.e-3_r8
  real(r8), parameter :: high_rel_err = 1.e-4_r8
  !-----------------------------------------------------------------------
  ! Newton-Raphson iteration limits
  !-----------------------------------------------------------------------
  integer, parameter  :: itermax   = 11
  integer, parameter  :: cut_limit = 5
  real(r8), parameter :: sol_min   = 1.e-20_r8
  real(r8), parameter :: small     = 1.e-40_r8

  real(r8) :: epsilon(clscnt4)
  logical  :: factor(itermax)

contains

  subroutine imp_slv_inti
    !-----------------------------------------------------------------------
    ! ... Initialize the implict solver
    !-----------------------------------------------------------------------
    use mo_chem_utls, only : get_spc_ndx

    implicit none

    !-----------------------------------------------------------------------
    ! ... Local variables
    !-----------------------------------------------------------------------
    integer  :: m, ox_ndx, o3a_ndx
    real(r8) :: eps(gas_pcnst)

    factor(:) = .true.
    eps(:) = rel_err

    ox_ndx = get_spc_ndx( 'OX' )
    if( ox_ndx < 1 ) then
       ox_ndx = get_spc_ndx( 'O3' )
    end if
    if( ox_ndx > 0 ) then
       eps(ox_ndx) = high_rel_err
    end if
    m = get_spc_ndx( 'NO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'N2O5' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'OH' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    o3a_ndx = get_spc_ndx( 'O3A' )
    if( o3a_ndx > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    do m = 1,clscnt4
       epsilon(m) = eps(clsmap(m,4))
    end do

  end subroutine imp_slv_inti

  subroutine imp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, &
                      ncol, nlev, lchnk, prod_out, loss_out  )
    !-----------------------------------------------------------------------
    !      	... imp_sol advances the volumetric mixing ratio
    !           forward one time step via the fully implicit euler scheme.
    !           this source is meant for vector architectures such as the
    !           nec sx6 and cray x1
    !-----------------------------------------------------------------------

    use chem_mods,     only : rxntot, extcnt, nzcnt, permute, cls_rxt_cnt
    use mo_tracname,   only : solsym
    use mo_lin_matrix, only : linmat
    use mo_nln_matrix, only : nlnmat
    use mo_lu_factor,  only : lu_fac
    use mo_lu_solve,   only : lu_slv
    use mo_prod_loss,  only : imp_prod_loss
    use mo_indprd,     only : indprd
    use time_manager,  only : get_nstep
    use perf_mod,      only : t_startf, t_stopf

    implicit none

    !-----------------------------------------------------------------------
    ! ... dummy args
    !-----------------------------------------------------------------------
    integer,  intent(in) :: ncol ! columns in chunck
    integer, intent(in)  :: nlev
    integer,  intent(in) :: lchnk ! chunk id
    real(r8), intent(in) :: delt ! time step (s)

    real(r8), intent(in) :: reaction_rates(ncol*nlev,max(1,rxntot)) ! rxt rates (1/cm^3/s)
    real(r8), intent(in) :: extfrc(ncol*nlev,max(1,extcnt))         ! external in-situ forcing (1/cm^3/s)
    real(r8), intent(in) :: het_rates(ncol*nlev,max(1,gas_pcnst))   ! washout rates (1/s)
    real(r8), intent(inout) :: base_sol(ncol*nlev,gas_pcnst)        ! species mixing ratios (vmr)
    real(r8), intent(out) :: prod_out(ncol*nlev,max(1,clscnt4))
    real(r8), intent(out) :: loss_out(ncol*nlev,max(1,clscnt4))
    !-----------------------------------------------------------------------
    ! ... local variables
    !-----------------------------------------------------------------------
    integer  ::  nr_iter
    integer  ::  ofl
    integer  ::  ofu
    integer  ::  avec_len
    integer  ::  bndx                                                      ! base index
    integer  ::  cndx                                                      ! class index
    integer  ::  pndx                                                      ! permuted class index
    integer  ::  i,m
    integer  ::  fail_cnt
    integer  ::  cut_cnt
    integer  ::  stp_con_cnt
    integer  ::  nstep
    real(r8) ::  interval_done
    real(r8) ::  dt
    real(r8) ::  dti
    real(r8) ::  max_delta(max(1,clscnt4))
    real(r8) :: ind_prd(ncol*nlev,max(1,clscnt4))
    logical  ::  convergence
    logical  ::  converged(max(1,clscnt4))
    integer  :: chnkpnts ! total spatial points in chunk; ncol*ncol
    logical :: diags_out(ncol*nlev,max(1,clscnt4))

    real(r8) :: sys_jac_blk(veclen,max(1,nzcnt))
    real(r8) :: lin_jac_blk(veclen,max(1,nzcnt))
    real(r8) :: solution_blk(veclen,max(1,clscnt4))
    real(r8) :: forcing_blk(veclen,max(1,clscnt4))
    real(r8) :: iter_invariant_blk(veclen,max(1,clscnt4))
    real(r8) :: prod_blk(veclen,max(1,clscnt4))
    real(r8) :: loss_blk(veclen,max(1,clscnt4))
    real(r8) :: ind_prd_blk(veclen,max(1,clscnt4))
    real(r8) :: sbase_sol_blk(veclen,gas_pcnst)
    real(r8) :: wrk_blk(veclen)
    logical :: spc_conv_blk(veclen,max(1,clscnt4))
    logical :: cls_conv_blk(veclen)
    real(r8) :: reaction_rates_blk(veclen,max(1,rxntot))
    real(r8) :: extfrc_blk(veclen,max(1,extcnt))
    real(r8) :: het_rates_blk(veclen,max(1,gas_pcnst))
    real(r8) :: base_sol_blk(veclen,gas_pcnst)

    chnkpnts = ncol*nlev
    prod_out = 0._r8
    loss_out = 0._r8
    diags_out = .false.

    !-----------------------------------------------------------------------      
    !        ... class independent forcing
    !-----------------------------------------------------------------------      
    if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
       call indprd( 4, ind_prd, clscnt4, base_sol, extfrc, &
            reaction_rates, chnkpnts )
    else
       do m = 1,clscnt4
          ind_prd(:,m) = 0._r8
       end do
    end if

    ofl = 1
    chnkpnts_loop : do 
       ofu = min( chnkpnts,ofl + veclen - 1 )
       avec_len = (ofu - ofl) + 1

       reaction_rates_blk(1:avec_len,:) = reaction_rates(ofl:ofu,:)
       extfrc_blk(1:avec_len,:) = extfrc(ofl:ofu,:)
       het_rates_blk(1:avec_len,:) = het_rates(ofl:ofu,:)
       ind_prd_blk(1:avec_len,:) = ind_prd(ofl:ofu,:)
       base_sol_blk(1:avec_len,:) = base_sol(ofl:ofu,:)

       do m = 1,gas_pcnst
         do i = 1, avec_len
           sbase_sol_blk(i,m) = base_sol_blk(i,m)
         end do
       end do
       !-----------------------------------------------------------------------      
       !        ... time step loop
       !-----------------------------------------------------------------------      
       dt            = delt
       cut_cnt       = 0
       fail_cnt      = 0
       stp_con_cnt   = 0
       interval_done = 0._r8
       time_step_loop : do
          dti = 1._r8 / dt
          !-----------------------------------------------------------------------      
          !        ... transfer from base to class array
          !-----------------------------------------------------------------------      
          do cndx = 1,clscnt4
             bndx = clsmap(cndx,4)
             pndx = permute(cndx,4)
             do i = 1, avec_len
                solution_blk(i,pndx) = base_sol_blk(i,bndx)
             end do
          end do
          !-----------------------------------------------------------------------      
          !        ... set the iteration invariant part of the function f(y)
          !-----------------------------------------------------------------------      
          if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
             do m = 1,clscnt4
                do i = 1, avec_len
                   iter_invariant_blk(i,m) = dti * solution_blk(i,m) + ind_prd_blk(i,m)
                end do
             end do
          else
             do m = 1,clscnt4
                do i = 1, avec_len
                    iter_invariant_blk(i,m) = dti * solution_blk(i,m)
                end do
             end do
          end if
          !-----------------------------------------------------------------------      
          !        ... the linear component
          !-----------------------------------------------------------------------      
          if( cls_rxt_cnt(2,4) > 0 ) then
             call t_startf( 'lin_mat' )
             call linmat( avec_len, lin_jac_blk, base_sol_blk, &
                  reaction_rates_blk, het_rates_blk )
             call t_stopf( 'lin_mat' )
          end if
          !=======================================================================
          !        the newton-raphson iteration for f(y) = 0
          !=======================================================================
          do i = 1, avec_len
             cls_conv_blk(i) = .false.
          end do
          iter_loop : do nr_iter = 1,itermax
             !-----------------------------------------------------------------------      
             !        ... the non-linear component
             !-----------------------------------------------------------------------      
             if( factor(nr_iter) ) then
                call t_startf( 'nln_mat' )
                call nlnmat( avec_len, sys_jac_blk, base_sol_blk, &
                     reaction_rates_blk, lin_jac_blk, dti )
                call t_stopf( 'nln_mat' )
                !-----------------------------------------------------------------------      
                !         ... factor the "system" matrix
                !-----------------------------------------------------------------------      
                call t_startf( 'lu_fac' )
                call lu_fac( avec_len, sys_jac_blk )
                call t_stopf( 'lu_fac' )
             end if
             !-----------------------------------------------------------------------      
             !   	... form f(y)
             !-----------------------------------------------------------------------      
             call t_startf( 'prod_loss' )
             call imp_prod_loss( avec_len, prod_blk, loss_blk, &
                  base_sol_blk, reaction_rates_blk, het_rates_blk )
             call t_stopf( 'prod_loss' )
             do m = 1,clscnt4
                do i = 1, avec_len
                   forcing_blk(i,m) = solution_blk(i,m)*dti &
                                    - (iter_invariant_blk(i,m) + prod_blk(i,m) - loss_blk(i,m))
                end do
             end do
             !-----------------------------------------------------------------------      
             !         ... solve for the mixing ratio at t(n+1)
             !-----------------------------------------------------------------------      
             call t_startf( 'lu_slv' )
             call lu_slv( avec_len, sys_jac_blk, forcing_blk )
             call t_stopf( 'lu_slv' )
             do m = 1,clscnt4
                do i = 1, avec_len
                   if( .not. cls_conv_blk(i) )then
                      solution_blk(i,m) = solution_blk(i,m) + forcing_blk(i,m)
                   else
                      forcing_blk(i,m)  = 0._r8
                   endif
                end do
             end do
             !-----------------------------------------------------------------------      
             !    	... convergence measures and test
             !-----------------------------------------------------------------------      
             conv_chk : if( nr_iter > 1 ) then
                !-----------------------------------------------------------------------      
                !    	... check for convergence
                !-----------------------------------------------------------------------      
                do cndx = 1,clscnt4
                   pndx = permute(cndx,4)
                   bndx = clsmap(cndx,4)
                   do i = 1, avec_len
                     if ( abs( solution_blk(i,pndx) ) > sol_min ) then
                        wrk_blk(i) = abs( forcing_blk(i,pndx)/solution_blk(i,pndx) )
                     else
                      wrk_blk(i) = 0._r8
                     endif
                   enddo

                   max_delta(cndx) = maxval( wrk_blk(1:avec_len) )
                   do i = 1, avec_len
                     solution_blk(i,pndx) = max( 0._r8,solution_blk(i,pndx) )
                     base_sol_blk(i,bndx) = solution_blk(i,pndx)
                     if ( abs( forcing_blk(i,pndx) ) > small ) then
                       spc_conv_blk(i,cndx) = abs(forcing_blk(i,pndx)) <= epsilon(cndx)*abs(solution_blk(i,pndx))
                     else
                       spc_conv_blk(i,cndx) = .true.
                     endif
                   enddo
                   where( spc_conv_blk(1:avec_len,cndx) .and. .not.diags_out(ofl:ofu,cndx) )
                      ! capture output production and loss diagnostics at converged ponits
                      prod_out(ofl:ofu,cndx) = prod_blk(1:avec_len,cndx) + ind_prd_blk(1:avec_len,cndx)
                      loss_out(ofl:ofu,cndx) = loss_blk(1:avec_len,cndx)
                      diags_out(ofl:ofu,cndx) = .true.
                   endwhere
                   converged(cndx) = all( spc_conv_blk(1:avec_len,cndx) )
                end do

                convergence = all( converged(:) )
                if( convergence ) then
                   base_sol(ofl:ofu,:) = base_sol_blk(1:avec_len,:)
                   ind_prd(ofl:ofu,:) = ind_prd_blk(1:avec_len,:)
                   exit iter_loop
                end if
                do m = 1, avec_len
                   if( .not. cls_conv_blk(m) ) then
                      cls_conv_blk(m) = all( spc_conv_blk(m,:) )
                   end if
                end do
             else conv_chk
!-----------------------------------------------------------------------
! ... limit iterate
!-----------------------------------------------------------------------
                do m = 1,clscnt4
                  do i = 1, avec_len
                    solution_blk(i,m) = max( 0._r8,solution_blk(i,m) )
                  end do
                end do
!-----------------------------------------------------------------------
! ... transfer latest solution back to base array
!-----------------------------------------------------------------------
                do cndx = 1,clscnt4
                   pndx = permute(cndx,4)
                   bndx = clsmap(cndx,4)
                   do i = 1, avec_len
                     base_sol_blk(i,bndx) = solution_blk(i,pndx)
                   end do
                end do
             end if conv_chk
          end do iter_loop

          !-----------------------------------------------------------------------      
          !    	... check for newton-raphson convergence
          !-----------------------------------------------------------------------      
          non_conv :  if( .not. convergence ) then
             !-----------------------------------------------------------------------      
             !   	... non-convergence
             !-----------------------------------------------------------------------      
             fail_cnt = fail_cnt + 1
             nstep    = get_nstep()
             write(iulog,'('' imp_sol: time step '',1p,g15.7,'' failed to converge @ (lchnk,ofl,ofu,nstep) = '',4i6)') &
                  dt,lchnk,ofl,ofu,nstep
             stp_con_cnt = 0
             step_reduction : if( cut_cnt < cut_limit ) then
                cut_cnt = cut_cnt + 1
                if( cut_cnt < cut_limit ) then
                   dt = .5_r8 * dt
                else
                   dt = .1_r8 * dt
                end if
                do m = 1,gas_pcnst
                   base_sol_blk(:,m) = sbase_sol_blk(:,m)
                end do
                cycle time_step_loop
             else step_reduction
                write(iulog,'('' imp_sol: step failed to converge @ (lchnk,ofl,ofu,nstep,dt,time) = '',4i6,1p,2g15.7)') &
                     lchnk,ofl,ofu,nstep,dt,interval_done+dt
                do m = 1,clscnt4
                   if( .not. converged(m) ) then
                      write(iulog,'(1x,a8,1x,1pe10.3)') solsym(clsmap(m,4)), max_delta(m)
                   end if
                end do
             end if step_reduction
          end if non_conv
          !-----------------------------------------------------------------------      
          !   	... check for interval done
          !-----------------------------------------------------------------------      
          interval_done = interval_done + dt
          time_step_done : if( abs( delt - interval_done ) <= .0001_r8 ) then
             if( fail_cnt > 0 ) then
                write(iulog,'(''imp_sol : @ (lchnk,ofl,ofu) = '',3i6,'' failed '',i2,'' times'')') &
                     lchnk,ofl,ofu,fail_cnt
             end if
             base_sol(ofl:ofu,:) = base_sol_blk(1:avec_len,:)
             ind_prd(ofl:ofu,:) = ind_prd_blk(1:avec_len,:)
             exit time_step_loop
          else time_step_done
             !-----------------------------------------------------------------------      
             !   	... transfer latest solution back to base array
             !-----------------------------------------------------------------------      
             if( convergence ) then
                stp_con_cnt = stp_con_cnt + 1
             end if
             do m = 1,gas_pcnst
                do i = 1, avec_len
                   sbase_sol_blk(i,m) = base_sol(i,m)
                end do
             end do
             if( stp_con_cnt >= 2 ) then
                dt = 2._r8*dt
                stp_con_cnt = 0
             end if
             dt = min( dt,delt-interval_done )
          end if time_step_done
       end do time_step_loop

       ofl = ofu + 1
       if( ofl > chnkpnts ) then
          base_sol(ofl:ofu,:) = base_sol_blk(1:avec_len,:)
          ind_prd(ofl:ofu,:) = ind_prd_blk(1:avec_len,:)
          exit chnkpnts_loop
       end if

       base_sol(ofl:ofu,:) = base_sol_blk(1:avec_len,:)
       ind_prd(ofl:ofu,:) = ind_prd_blk(1:avec_len,:)

    end do chnkpnts_loop

  end subroutine imp_sol

end module mo_imp_sol
