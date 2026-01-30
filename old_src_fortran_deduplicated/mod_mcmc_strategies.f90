module mod_mcmc_strategies
    use parz
    use RDistributions
    use routinez
    use mod_mcmc_helpers
    implicit none

contains

    subroutine update_effects_mixture()
       integer :: k, j, kk, l, i, snploc, annotflag
       logical :: overflow
       double precision :: maxtemp, ssculm
    
       snptracker=0d0
       
       do k=1,nloci         
          snploc=permvec(k)
          ! 1. Allocation Step (Sample Category 'a')
          if (nannot(snploc)>1) THEN
             z => X(:,snploc)
             zz=xpx(snploc)
             zz_vare=zz/vare
             gk=g(snploc)
             atemp=0d0
             if (rep /=1) THEN
                atemp(a(snploc))=1
             endif
             dira=C(snploc,1:ncat)+atemp
             pia=rdirichlet2(ncat,dira)
             if (vsnptrack(snploc)>1) THEN
                ytemp=yadj+z*gk
             endif
             rhs=dot_product(ytemp,z)
             
             ! Compute log likelihoods for each category
             ss=0d0
             maxs=0d0

             do i=2,ndist
                uhat=rhs/(zz+vare_gp(i))
                maxtemp=0.5d0*uhat*rhs/vare
                if (maxtemp > maxs) maxs=maxtemp
             enddo
             
             do j=1,ncat
                if (C(snploc,j)==1) THEN
                   ss(j)=p(1,j)*exp(-maxs)
                   do kk=2,ndist 
                      detV=gp(kk)*zz_vare+1.0d0
                      uhat=rhs/(zz+vare_gp(kk))
                      ss(j)=ss(j)+p(kk,j)*detV**(-0.5d0)*exp(0.5d0*uhat*rhs/vare-maxs)
                   enddo
                   ss(j)=dlog(pia(j))+dlog(ss(j))
                endif
             enddo
             
             ! Softmax / Sampling
             sstemp=0d0
             do kk=1,ncat
                 skk=0.0d0
                 if (C(snploc,kk)==1) THEN
                    skk=ss(kk)
                    sk=0.0d0
                    overflow=.false.
                    do l=1,ncat
                       if (C(snploc,l)==1) THEN
                          if(l==kk) cycle
                          clike=ss(l)-skk
                          if(clike .lt. -700) then !undeflow
                             cycle
                          else if (clike .gt. 700) then 
                             overflow=.true.
                             exit
                          endif
                          sk=sk+dexp(clike)
                       endif
                    enddo
                    if (overflow .eqv. .true.) then
                       sstemp(kk) = 0.0
                    else
                       sstemp(kk)=1.0d0/(1.0d0+sk)
                    endif
                 endif
             enddo
             
             ssculm=0.0d0
             call random_number(r)
             annotflag=1
             do kk=1,ncat
                if (C(snploc,kk)==1) THEN
                   ssculm=ssculm+sstemp(kk)
                   if (r<ssculm) then
                      annotflag=kk
                      exit
                   endif
                endif
             enddo
             a(snploc)=annotflag
          endif
       enddo
    
       ! 2. Sample Effect 'g' for selected category 'a(snploc)'
       do k=1,nloci
          snploc=permvec(k)
          j=a(snploc)
          z => X(:,snploc)
          zz=xpx(snploc)
          zz_vare=zz/vare
          gk=g(snploc)
          if(vsnptrack(snploc) > 1) then
             yadj=yadj+z*gk
          endif
          rhs= dot_product(yadj,z)
          
          call sample_distribution_component(j, zz, zz_vare, rhs)
          
          snptracker(snploc,j)=indistflag
          vsnptrack(snploc)=indistflag
          
          if(indistflag==1) then
             gk=0.0d0
          else
             v1=zz+vare/gp(indistflag)
             gk=rand_normal(rhs/v1, dsqrt(vare/v1))
             yadj=yadj-z*gk  
             included=included+1
          endif
          g(snploc)=gk
          if(msize>0 .and. rep>mrep) then
             if(included>=msize) exit
          endif
       enddo
       
       ! Stats
       do j=1,ncat
          do i=1,ndist
             snpindist(i,j)=count(snptracker(:,j)==i)
             varindist(i,j)=sum(g*g, mask= snptracker(:,j)==i)
          enddo
       enddo
       included=nloci-sum(snpindist(1,:))
       
    end subroutine update_effects_mixture
    
    subroutine update_effects_additive_loci_major()
       integer :: k, kcat, j, i, snploc
       
       call permutate(permannot,ncat)
       
       do k=1,nloci
          snploc=permvec(k)
          gk=g(snploc)
          z => X(:,snploc)
          zz=xpx(snploc)
          zz_vare=zz/vare
          if(vsnptrack(snploc) > 1) then
             yadj=yadj+z*gk
          endif
          
          do kcat=1,ncat
             j=permannot(kcat)
             call update_log_p(j) 
             
             if (C(snploc,j) == 1) THEN
                 rhs= dot_product(yadj,z)
                 call sample_distribution_component(j, zz, zz_vare, rhs)
                 
                 snptracker(snploc,j)=indistflag
                 
                 if(indistflag==1) then
                    gk=0.0d0
                 else
                    v1=zz+vare/gp(indistflag)
                    gk=rand_normal(rhs/v1, dsqrt(vare/v1))
                    yadj=yadj-z*gk 
                 endif
                 gannot(snploc,j)=gk
                 if(msize>0 .and. rep>mrep) then
                   ! Logic for breaking early?
                 endif
             endif
          enddo 
       enddo
       
       ! Sum loci effects
       g=sum(gannot,dim=2)
       vsnptrack=maxval(snptracker,dim=2)
       
       do j=1,ncat
          do i=1,ndist
             snpindist(i,j)=count(snptracker(:,j)==i)
             varindist(i,j)=sum(gannot(:,j)*gannot(:,j), mask= snptracker(:,j)==i)
          enddo
       enddo
       
       ! Included count
       includedloci=0d0
       do i=1,nloci
          if (vsnptrack(i)>1) includedloci(i)=1
       enddo
       included=sum(includedloci)
    
    end subroutine update_effects_additive_loci_major
    
    subroutine update_effects_additive_category_major()
       integer :: j, k, i, snploc
       
       do j=1,ncat
          call update_log_p(j)
          do k=1,nloci
             snploc=permvec(k)
             if (C(snploc,j) == 1) THEN
                z => X(:,snploc)
                zz=xpx(snploc)
                zz_vare=zz/vare
                gk=gannot(snploc,j)
                
                if(snptracker(snploc,j) > 1) then
                   yadj=yadj+z*gk
                endif
                
                rhs= dot_product(yadj,z)
                
                call sample_distribution_component(j, zz, zz_vare, rhs)
                
                snptracker(snploc,j)=indistflag
                vsnptrack(snploc)=indistflag
                
                if(indistflag==1) then
                   gk=0.0d0
                else
                   v1=zz+vare/gp(indistflag)
                   gk=rand_normal(rhs/v1, dsqrt(vare/v1))
                   yadj=yadj-z*gk  
                   included=included+1
                endif
                gannot(snploc,j)=gk
                
             endif
          enddo
       enddo
    
       ! Sum loci effects
       g=sum(gannot,dim=2)
       
       do j=1,ncat
          do i=1,ndist
             snpindist(i,j)=count(snptracker(:,j)==i)
             varindist(i,j)=sum(gannot(:,j)*gannot(:,j), mask= snptracker(:,j)==i)
          enddo
       enddo
       
       ! Included count
       includedloci=0d0
       do j=1,ncat
          do i=1,nloci
             if (snptracker(i,j)>1) includedloci(i)=1
          enddo
       enddo
       included=sum(includedloci)
    
    end subroutine update_effects_additive_category_major

end module mod_mcmc_strategies
