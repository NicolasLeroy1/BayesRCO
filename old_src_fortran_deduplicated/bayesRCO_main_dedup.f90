program bayesR_dedup
use parz
use cmd_parser
use routinez
use RDistributions
use mod_mcmc_core

implicit none
integer :: i, j, k, kk, jj,snploc, l, aa, kcat
character (len=8)  :: cdate
character (len=10) :: ctime, ci, ca, cj
logical :: overflow

call date_and_time(date=cdate,time=ctime)
call parse

call get_size
call load_phenos_plink
call allocate_data
call parse_priors
call load_categories

if(mcmc) then
   open(unit=22,file=logfil,status='unknown',form='formatted')
   write(22,901) 'Program BayesRCO Deduplicated'
   write(22,908) 'Run started at',cdate(1:4),cdate(5:6),cdate(7:8),ctime(1:2),ctime(3:4),ctime(5:6)
   write(22,902) 'Prefix for input files',trim(inprefix)
   write(22,902) 'Prefix for output files',trim(outprefix)
   write(22,903) 'Phenotype column',trait_pos
   write(22,903) 'No. of loci',nloci
   write(22,903) 'No. of individuals',nind
   write(22,903) 'No. of training individuals',nt
   write(22,906) 'Prior Vara', vara, dfvara
   write(22,906) 'Prior Vare', vare, dfvare
   write(22,903) 'Model size',msize
   write(22,903) 'No. of cycles',numit
   write(22,903) 'Burnin ',burnin
   write(22,903) 'Thinning rate',thin
   write(22,903) 'No. of mixtures',ndist
   write(22,905) 'Variance of dist ', gpin
   write(22,905) 'Dirichlet prior', delta
   write(22,903) 'Seed ', seed1
   write(22,909) 'SNP output ', snpout
   write(22,909) 'Cat output', cat
   write(22,909) 'Beta output', beta
   write(22,903) 'No. of SNP categories ', ncat
   call flush(22)
endif

call load_snp_binary
call xcenter 
call init_random_seed


if(mcmc) then
   call run_mcmc_deduplicated()
else   ! end mcmc
 
   open(unit=21,file=logfil,status='unknown',form='formatted')
   write(21,901) 'Program BayesR'
   write(21,908) 'Run started at',cdate(1:4),cdate(5:6),cdate(7:8),ctime(1:2),ctime(3:4),ctime(5:6)
   write(21,902) 'Prefix for input files',trim(inprefix)
   write(21,902) 'Prefix for output files',trim(outprefix)
   write(21,903) 'Phenotype column',trait_pos
   write(21,903) 'No. of loci',nloci
   write(21,903) 'No. of individuals',nind
   write(21,903) 'No. of individuals to predict',nt
   call load_param
   call compute_dgv
   call write_dgv
end if

call date_and_time(date=cdate,time=ctime)
write(21,908) 'Run ended at',cdate(1:4),cdate(5:6),cdate(7:8),ctime(1:2),ctime(3:4),ctime(5:6)
close(21)

901 format(a)
902 format(a,t30,': ',a)
903 format(a,t30,'= ',i8)
904 format(a,t30,'= ',f20.6)
905 format(a,t30,'= ',10f10.5)
906 format(a,t30,'= ',2f10.6)
907 format(a,t30,'= ',f10.2,a)
908 format(a20,1x,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a2)
909 format(a,t30,'= ',l)

end program bayesR_dedup
