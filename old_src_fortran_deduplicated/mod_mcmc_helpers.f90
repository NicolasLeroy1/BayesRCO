module mod_mcmc_helpers
    use parz
    use RDistributions
    implicit none

contains

    subroutine sample_distribution_component(j, zz, zz_vare, rhs)
       integer, intent(in) :: j
       double precision, intent(in) :: zz, zz_vare, rhs
       integer :: kk, l
       double precision :: ssculm
       logical :: overflow

       s(1)=log_p(1,j)
       do kk=2,ndist
          logdetV=dlog(gp(kk)*zz_vare+1.0d0)
          uhat=rhs/(zz+vare_gp(kk))
          s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(kk,j)
       enddo
       
       stemp=0.0d0
       do kk=1,ndist
          skk=s(kk)
          sk=0.0d0
          overflow=.false.
          do l=1,ndist
             if(l==kk) cycle
             clike=s(l)-skk
             if(clike .lt. -700) then !undeflow
                cycle
             else if (clike .gt. 700) then 
                overflow=.true.
                exit
             endif
             sk=sk+dexp(clike)
          enddo
          if (overflow .eqv. .true.) then
             stemp(kk) = 0.0
          else
             stemp(kk)=1.0d0/(1.0d0+sk)
          endif
       enddo
       
       ssculm=0.0d0
       call random_number(r)
       indistflag=1
       do kk=1,ndist
          ssculm=ssculm+stemp(kk)
          if (r<ssculm) then
             indistflag=kk
             exit
          endif
       enddo
    end subroutine sample_distribution_component

    subroutine update_log_p(j)
       integer, intent(in) :: j
       integer :: i
       log_p(1,j)=dlog(p(1,j))
       do i=2,ndist
          log_p(i,j)=dlog(p(i,j))
       enddo
    end subroutine update_log_p

end module mod_mcmc_helpers
