module mod_mcmc_core
    use parz
    use RDistributions
    use routinez
    use mod_mcmc_helpers
    use mod_mcmc_strategies
    implicit none

contains

    subroutine run_mcmc_deduplicated()
       implicit none
       integer :: i, j
    
       call mcmc_init_common()
    
       if (nobayesCpi) THEN
          if(mixture) THEN
             call mcmc_init_mixture()
             call mcmc_loop(strategy=1) ! 1 = Mixture
          else ! Additive
             call mcmc_init_additive()
             call mcmc_loop(strategy=2) ! 2 = Additive Loci-Major
          endif
       else ! BayesCpi
          call mcmc_init_bayesCpi()
          call mcmc_loop(strategy=3) ! 3 = Additive Category-Major
       endif
    
       call mcmc_finalize()
    
    end subroutine run_mcmc_deduplicated
    
    !-------------------------------------------------------------------------------
    ! INITIALIZATION ROUTINES
    !-------------------------------------------------------------------------------
    subroutine mcmc_init_common()
       implicit none
       integer :: i, j
       character (len=10) :: ci, ca, cj
    
       nnind=dble(nt)
       if(snpout) open(unit=14,file=locfil,status='unknown',action='write')
       if(cat)    open(unit=100,file=catfil,status='unknown',action='write')
       if(beta)   open(unit=101,file=betafil,status='unknown',action='write')
       
       open(unit=25,file=hypfil,status='unknown',form='formatted')
       write(25,'(2(A10,1x),3(A12,1x),A7)',advance='no') 'Replicate','Nsnp','mu','Va','Ve',' '
       do j=1,ncat
          write(cj,'(I8)') j
          cj=adjustl(cj)
          do i=1,ndist
             write(ci,'(I8)') i
             ci=adjustl(ci)
             ca="Nk"//trim(ci)//"_"//trim(cj)
             write(25,'(A10,1x)',advance="no") ca
          enddo
       end do
       do j=1,ncat
          write(cj,'(I8)') j
          cj=adjustl(cj)
          do i=1,ndist
             write(ci,'(I8)') i
             ci=adjustl(ci)
             ca="Vk"//trim(ci)//"_"//trim(cj)
             write(25,'(A12)',advance="no") ca
          enddo
       end do
       write(25,*)
    
       ! Vara/Vare priors handling
       if(dfvara < -2) then
          VCE=.false.
          yhat=sum(why, mask=trains==0)/nnind
          vary= sum((why-yhat)*(why-yhat),mask=trains==0)/(nnind-1.0d0)
          vara=vara*vary
       else
          VCE=.true.
          vara_ap=vara
          vare_ap=vare
          if(dfvara == -2) vara_ap=0d0
          if(dfvare == -2) vare_ap=0d0
       endif
    
       ! Initialize stores
       pstore=0d0
       gstore=0d0
       varistore=0d0
       mu_vare_store=0
       snpstore=0d0
       indiststore=0d0
       varustore=0d0
       
       ! Initialize counters and matrices
       xpx=0d0
       do i=1,nloci
          xpx(i)=dot_product(X(:,i),X(:,i))
       enddo
       nannot=sum(C(:,1:ncat),dim=2)
       
       ! Starting values
       mu=1.0d0
       yadj=0.0d0
       yhat=sum(why, mask=trains==0)/nnind
       vary= sum((why-yhat)*(why-yhat),mask=trains==0)/(nnind-1.0d0)
       gp=gpin*vara
       scale=0.0d0
       
       do j=1,ncat
          p(1,j)=0.5d0
          p(2:ndist,j)=1.0d0/gpin(2:ndist)
          p(2:ndist,j)=0.5*p(2:ndist,j)/sum(p(2:ndist,j))
       enddo
       g=dsqrt(vara/(0.5*dble(nloci)))
       
       do i=1,nloci
          permvec(i)=i
       enddo
    
       call compute_residuals
    
    end subroutine mcmc_init_common
    
    subroutine mcmc_init_mixture()
       integer :: k, j
       print *, 'mixture'
       vsnptrack=2d0
       snptracker=0d0
       ! Fill the annotation vector for 1 annotation SNP
       do k=1,nloci
          if (nannot(k)==1) THEN
             do j=1,ncat
                if (C(k,j)==1) THEN
                   a(k)=j
                endif
             enddo
          endif
       enddo
    end subroutine mcmc_init_mixture
    
    subroutine mcmc_init_additive()
       integer :: kcat
       print *, 'additive'
       vsnptrack=2d0
       snptracker=0d0
       do kcat=1,ncat
          permannot(kcat)=kcat
       enddo
    end subroutine mcmc_init_additive
    
    subroutine mcmc_init_bayesCpi()
       integer :: j, k
       print *, 'bayesCpi'
       annotstore=C(:,1:ncat)
       gannot=0d0
       snptracker=0d0
       do j=1,ncat
          do k=1,nloci
             if (C(k,j)==1) then
                snptracker(k,j)=2
             else
                snptracker(k,j)=0
             endif
          enddo
       enddo
    end subroutine mcmc_init_bayesCpi
    
    
    !-------------------------------------------------------------------------------
    ! MAIN LOOP
    !-------------------------------------------------------------------------------
    subroutine mcmc_loop(strategy)
       implicit none
       integer, intent(in) :: strategy
       integer :: i
       
       do rep=1,numit
          included=0
          
          ! 1. Sample Residual Variance (if not VCE)
          if(.not. VCE) then
             vare=dot_product(yadj,yadj)/rand_chi_square(nnind+3.0d0)
          endif
          
          ! 2. Sample Mu
          yadj=yadj+mu
          mu=rand_normal(sum(yadj)/nnind, dsqrt(vare/nnind))
          yadj=yadj-mu
          
          ! 3. Update cached variance quantities
          do i=2,ndist
             log_gp(i)=dlog(gp(i))
             vare_gp(i)=vare/gp(i)
          enddo
          
          ! 4. Permutation
          if(permute) then
             call permutate(permvec,nloci)
          endif
          
          ! 5. Update Effects (Strategy specific)
          if (strategy == 1) then
             call update_effects_mixture()
          elseif (strategy == 2) then
             call update_effects_additive_loci_major()
          elseif (strategy == 3) then
             call update_effects_additive_category_major()
          endif
          
          ! 6. Compute stats for iteration
          call compute_iteration_stats(strategy)
          
          ! 7. Save Samples
          if(mod(rep,thin)==0 .and. rep>burnin) then
             call save_samples(strategy)
          endif
          
       enddo 
    end subroutine mcmc_loop
    
    
    !-------------------------------------------------------------------------------
    ! UPDATE HYPERPARAMS & STATS
    !-------------------------------------------------------------------------------
    subroutine compute_iteration_stats(strategy)
       integer, intent(in) :: strategy
       integer :: j
       
       if(VCE) then
          scale=(dble(included)*sum(g**2) + vara_ap*dfvara)/(dfvara+dble(included))
          vara=rand_scaled_inverse_chi_square(dble(included)+dfvara,scale)
          if (strategy==3) then
               gp(2)=vara/included ! Specific to BayesCpi block 3 ???
               ! Wait, Block 3 line 809: gp(2)=vara/included
               ! Block 1 & 2: gp=gpin*vara
               ! This is a real difference.
               ! Wait, why would gp(2) be set like that?
               ! Let's check original file.
               ! Line 809 of original: gp(2)=vara/included
               ! Line 581 of original: gp=gpin*vara
          else
               gp=gpin*vara
          endif
          
          vare=(dot_product(yadj,yadj)+vare_ap*dfvare)/ (nnind+dfvare)
          vare=rand_scaled_inverse_chi_square(nnind+dfvare,vare)
       endif
    
       do j=1,ncat
          dirx=dble(snpindist(:,j))+delta
          p(:,j)=rdirichlet(ndist,dirx)
          call update_log_p(j)
       enddo
    end subroutine compute_iteration_stats
    
    
    !-------------------------------------------------------------------------------
    ! SAVE SAMPLES
    !-------------------------------------------------------------------------------
    subroutine save_samples(strategy)
       integer, intent(in) :: strategy
       integer :: i, j, jj, aa
       
       counter=counter+1
       gstore=gstore+g
       pstore=pstore+p
       varistore=varistore+g**2
       mu_vare_store(1)=mu_vare_store(1)+mu
       mu_vare_store(2)=mu_vare_store(2)+included
       mu_vare_store(3)=mu_vare_store(3)+vara
       mu_vare_store(4)=mu_vare_store(4)+vare
       varstore=varstore+varindist
       snpstore=snpstore+snpindist
       if(counter>1) then
          varustore=varustore+(counter*g-gstore)**2/(counter*(counter-1))
       endif
       
       do i=1,nloci
          if (strategy == 1) then
             jj=vsnptrack(i)
             indiststore(i,jj)=indiststore(i,jj)+1
             aa=a(i)
             annotstore(i,aa)=annotstore(i,aa)+1
          else
             do j=1,ncat
                jj=snptracker(i,j)
                if (jj>0) THEN
                   indiststore(i,jj)=indiststore(i,jj)+1
                endif
             enddo
          endif
       enddo
       
       write(25,'(i10,1x,i10,1x,3(E15.7,1x),100(i10,1x))',advance='no')  rep, included , & 
            mu, vara, vare, snpindist
       write(25,'(100E15.7,1x)') varindist
       call flush(25)
       
       if(beta) call output_beta
    end subroutine save_samples
    
    subroutine mcmc_finalize()
       integer :: i
       !posterior means
       gstore=gstore/counter
       pstore=pstore/counter
       mu_vare_store=mu_vare_store/counter
       varstore=varstore/counter
       snpstore=snpstore/counter
       varustore=varustore/counter
       varistore=varistore/counter
       do i=1,nloci
          indiststore(i,:)=indiststore(i,:)/(counter*nannot(i))
          if (nobayesCpi .and. mixture) annotstore(i,:)=annotstore(i,:)/counter
       enddo
       call output_model
       mu=mu_vare_store(1)
       call compute_dgv
       call write_dgv
    end subroutine mcmc_finalize

end module mod_mcmc_core
