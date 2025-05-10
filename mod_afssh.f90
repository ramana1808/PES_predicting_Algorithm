Module mod_afssh
!! Hammes-Schiffer, Tully, JCP 101, 4657 (1994)
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Potential
integer nquant
real*8 e_alpha_C(3),e_alpha_i(3),r0,l0,charge_H
real*8 a_solute,b_solute,c_solute,Diss_A,d_A,d_B,n_A,n_B
real*8 mu,dmu_drAB,k_coup,omg_s
real*8,allocatable :: mass(:),omg(:),charge(:)
real*8 temperature,gamma_B,beta
real*8 rAH,rAB
real*8 a_a,b_b,c_c,dA,x_min(2),ll,e_A_cov,e_A_ion,e_B_cov,e_B_ion
real*8 r_cutoff

!! Input/Output
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8 r_avg
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:)
complex*16,allocatable :: rho(:,:,:)

!! Classical
integer nclass,idistribution
real*8,allocatable :: x(:),v(:),acc(:)
real*8,allocatable :: x_old(:),v_old(:),acc_old(:),x_hop(:)
real*8 tim_hop
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_hop,flag_ortho
complex*16,allocatable :: delr(:,:,:),delp(:,:,:),delacc(:,:,:)
complex*16,allocatable :: delr_old(:,:,:),delp_old(:,:,:) 

!! Quantized vibration
integer ncl_site
integer nb_vib,n_dvr
real*8,allocatable :: ke_dvr(:,:),x_dvr(:)
real*8,allocatable ::si_sho(:,:,:),sho_overlap(:,:,:,:),q_exp(:,:,:)
real*8,allocatable ::qsq_exp(:,:,:)
real*8,allocatable ::en_sho(:,:),fc_init(:)
real*8,allocatable ::Hamil_diab_0(:,:)

!! Fitting
real*8 V0(2),cc(2),dd(2),q0(2),qr,qp,qt,Vc
real*8 parameters(10)

!! Quantum
integer state,nbasis,state_tentative,no_sites
integer state_old
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:),vdotd(:,:),V_k_old(:)
real*8,allocatable :: Hamil_site(:,:),Hamil_diab(:,:),delH_dels(:,:,:),delH_dels_ad(:,:,:)
real*8,allocatable :: pot(:,:),force(:,:,:),force_old(:,:,:),delF(:,:,:)
complex*16,allocatable :: ci(:),ci_old(:),sigma(:,:)
real*8,allocatable :: si_adiab_prev(:,:)
complex*16,allocatable :: mat(:,:),mat_adiab(:,:)
real*8,allocatable :: hop_prob(:),W_overlap(:,:),hop_prob_net(:)

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,total_time,curr_time,traj_num,tim_eq
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold,cnt_rate
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2
real*8 tim_T_jk
integer,allocatable:: seed(:)
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

! RAMANA

!! Classical MD Params 
real*8 :: md_sampling_t, md_total_t
real*8, allocatable :: x_r(:), x_p(:), x_ts(:)


!! Number of points in q_grid, test and train 
integer :: n_q_grid, ndims, test_points, train_points 

!! Model Parameters
real*8, allocatable :: q_grid(:), theta_c(:,:)
real*8, allocatable :: test_coords(:,:,:), train_coords(:,:,:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,j,size_seed,seed2(2)
  real*8 rnd,c_0,c_e,kt

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="AFSSH.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) total_time
  read(10,*) iwrite
  read(10,*) nstep_write
  read(10,*) nstep_avg
  read(10,*) idistribution
  read(10,*) flag_frust
  read(10,*) flag_ortho
  read(10,*) energy_cutoff
  read(10,*) nclass
  read(10,*) n_dvr
  read(10,*) nquant
  read(10,*) temperature
  read(10,*) icollapse
  read(10,*) ifriction
  read(10,*) seed2
  read(10,*) n_q_grid
  read(10,*) md_total_t
  read(10,*) md_sampling_t
  read(10,*) test_points
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  nbasis=n_dvr

  energy_cutoff=energy_cutoff*wave_to_J
  kt=kb*temperature
  beta=1.d0/(kb*temperature)
  nsteps=nint(total_time/dtc)+1

  ! RAMANA
  train_points = 2*(md_total_t/md_sampling_t+1) - test_points +(md_total_t/(2*md_sampling_t)+1)
  ndims = (nclass+1)*(nclass+2)/2

  !-----------------------------------------------------------------  
  i=nsteps/nstep_avg+1
  allocate(pop(nquant,i),pop_surf(nquant,i),pop_amp(nquant,i))
  allocate(rho(nbasis,nbasis,i))
  allocate(x(nclass),v(nclass),acc(nclass))
  allocate(x_old(nclass),v_old(nclass),acc_old(nclass),x_hop(nclass))
  allocate(mass(nclass),omg(nclass),charge(nclass))
  allocate(delr(nquant,nquant,nclass),delp(nquant,nquant,nclass),delacc(nquant,nquant,nclass))
  allocate(delr_old(nquant,nquant,nclass),delp_old(nquant,nquant,nclass))
  allocate(si_adiab(nbasis,nquant),ci(nquant),V_k(nquant),V_k_old(nquant),sigma(nquant,nquant))
  allocate(Hamil_site(no_sites,no_sites),Hamil_diab(nbasis,nbasis))
  allocate(delH_dels(nbasis,nbasis,nclass),delH_dels_ad(nquant,nquant,nclass))
  allocate(pot(nquant,nquant),force(nquant,nquant,nclass),force_old(nquant,nquant,nclass),delf(nquant,nquant,nclass))
  allocate(mat(nbasis,nbasis),mat_adiab(nquant,nquant))
  allocate(d_ij(nquant,nquant,nclass),vdotd(nquant,nquant),hop_prob(nquant),W_overlap(nquant,nquant))
  allocate(hop_prob_net(nquant))
  allocate(ci_old(nquant),si_adiab_prev(nbasis,nquant))
  allocate(si_sho(n_dvr,nb_vib,no_sites))
  allocate(ke_dvr(n_dvr,n_dvr),x_dvr(n_dvr))
  allocate(sho_overlap(nb_vib,nb_vib,no_sites,no_sites))
  allocate(q_exp(nb_vib,nb_vib,no_sites),qsq_exp(nb_vib,nb_vib,no_sites))
  allocate(en_sho(nb_vib,no_sites),fc_init(nb_vib))
  allocate(Hamil_diab_0(nbasis,nbasis))
  !RAMANA - model params 
  allocate(x_r(nclass),x_p(nclass),x_ts(nclass))
  allocate(q_grid(n_q_grid),theta_c(ndims, n_q_grid))
  allocate(test_coords(n_q_grid, test_points, nclass))
  allocate(train_coords(n_q_grid,train_points,nclass))

  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0
  
end subroutine setup

! --------------------------------------------------------
subroutine main 
  implicit none
  integer i,j,k,n
  real*8 t1,t2 

  call files(0)

  call cpu_time(t1)

  call setup_parameters
  call initialize_averages
  call evaluate_variables(0) 

  !call test_coeffs
  call rmse_model
  !stop
  !call draw_elec_pot
  call draw_pes
  !call check_acceleration

  do i=1,N_traj
    traj_num=i
    call init_cond
    call evolve(nsteps)
    call average_end
  enddo
  call write_average

  call cpu_time(t2);tim_tot=tim_tot+t2-t1
  call files(1)

end subroutine main
!---------------------------------------------------------- 

subroutine files(flag)
  implicit none
  integer,intent(in)::flag

  if(flag==0) then
    open(10,file="output")
    open(11,file="output_cl")
    open(12,file="output_qm")
    open(13,file="output_hop")
    open(14,file="output_overlap")
    open(15,file="output_dec")

    open(100,file="pop.out")
    open(101,file="cnts.out")
    open(16, file='check_acc')
  else
    write(10,*)
    write(10,*)"Total time=",tim_tot
    close(10);close(11);close(12);close(13);close(14);close(15)
    close(100);close(101)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine initialize_averages
  implicit none

  cnt_frust=0.d0
  cnt_collapse=0.d0
  cnt_init=0.d0
  cnt_term=0.d0
  pop=0.d0
  rho=0.d0
  pop_surf=0.d0
  pop_amp=0.d0

end subroutine initialize_averages
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i
  real*8 sig_x,sig_p,rnd,ak,su
  real*8 energy_0

  do i=1,nclass
    !ak=2/(hbar*omg(i))*dtanh(beta*hbar*omg(i)/2.d0) !! Wigner
    ak=beta    !! Classical
    sig_x=1.d0/dsqrt(ak*mass(i)*omg(i)**2)
    sig_p=dsqrt(mass(i)/ak)
    call gaussian_random_number(rnd)
    x(i)=rnd*sig_x
    call gaussian_random_number(rnd)
    v(i)=(1.d0/mass(i)*(rnd*sig_p))
  enddo

  x(1)=x(1)+2.7d-10
  x(2)=x(2)+0.05d-10

  state=1

  call evaluate_variables(0)
  call evaluate_variables(1)

  !! quantum state initialized on diabat 1
  !ci=si_adiab(state,:)
  ci(state)=1.d0

  call random_number(rnd)
  su=0.d0
  do i=1,nquant
    su=su+cdabs(ci(i))**2
    if(rnd<su) then
      state=i
      exit
    endif
  enddo

  delr=0.d0
  delp=0.d0

  ihop=1
  iaverage=1
  iterminate=0
  flag_terminate=0

  curr_time=0.d0
  call evaluate_variables(0)
  call evaluate_variables(1)
  call compute_mat_diab

  !! to compute the standard deviation of the energy of the trajectory
  en_avg=0.d0;ensq_avg=0.d0
  nst_av=0

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in) :: nsteps
  integer i,j,nstep_sm,iflag_coll,i_do_something
  real*8 t1,t2
  integer iterm

  !call cpu_time(t1)

  call write_output(1,1)
  iterm=0
  do i=1,nsteps
    call write_output(i,0)
    call average(i)
    call save_old_state
    call evolve_classical(dtc)
    call evolve_quantum_small_dtq
    if(ihop==1)call hop
    !if(i_do_something.ne.1.and.icollapse==1)call
    !collapse(dtc,iflag_coll)
    if(flag_terminate==1) call traj_terminate(iterm)
      if(iterm==1)exit

    curr_time=curr_time+dtc
  enddo
  call write_output(1,1)

  !call cpu_time(t2)
  !tim_evolve=tim_evolve+t2-t1

end subroutine evolve
!-----------------------------------------------------------------  

subroutine do_something(i_do_something)
  !! subroutine to maintain energy conservation
  implicit none
  integer,intent(in)::i_do_something
  real*8 acc_tent(nclass),dt,dtm(3)
  integer i,nstep

  if(i_do_something==1) then
    !! On first pass, check if hop happens; if yes, check if evolution
    !using the acceleration of hopped surface conservses energy
    !! Useful for very sharp crossings
    call evolve_quantum_small_dtq
    if(flag_hop==1) then
      state=state_tentative
      call evaluate_variables(0)
      v=v_old+0.5*(acc_old+acc)*dtc
      call evaluate_variables(1)
    endif
  else
    !! If the first pass did not work, reduce time-steps untill energy
    !is conserved
    dtm=1.d0 !! some large value
    dtm(1)=0.1d0/maxval(vdotd)
    dtm(2)=0.5*dtc*dsqrt(energy_cutoff/dabs(energy-energy_old))
    dtm(3)=dtc

    dt=minval(dtm)

    dt=dt/real(i_do_something)
    nstep=nint(dtc/dt)
    dt=dtc/real(nstep)
    call revert_state
    do i=1,nstep
      call evolve_classical(dt)
    enddo
  endif


end subroutine do_something
!-----------------------------------------------------------------  

subroutine average(i)
  implicit none
  integer,intent(in) :: i
  integer j,i1,j1,k,kp
  complex*16 ci_diab(nbasis),rho_ad(nquant,nquant)
  real*8 U(nbasis,nquant),U_exc(nquant,nquant)
  integer if_reactant
  real*8 t1,t2

  !call cpu_time(t1)

  if(iwrite==1) then
    en_avg=en_avg+energy
    ensq_avg=ensq_avg+energy*energy
    nst_av=nst_av+1
  endif

  if(iaverage==1.and.(mod(i,nstep_avg)==1.or.nstep_avg==1)) then
    if(nstep_avg==1) then
      j=i
    else
      j=i/nstep_avg+1
    endif

    !! Diabatic population
    !! J. Chem. Phys. 139, 211101 (2013)
    !U_exc(1,1)=-0.88449142962491789
    !U_exc(1,2)=-0.46655643915829631
    !U_exc(2,1)=-0.46655643915829631
    !U_exc(2,2)=0.88449142962491778
    !U=matmul(U_exc,si_adiab)
    !U=si_adiab
    !rho_ad=0.d0
    !rho_ad(state,state)=1.d0
    !do i1=1,nquant
    !  do j1=1,nquant
    !    if(i1.ne.j1) rho_ad(i1,j1)=ci(i1)*dconjg(ci(j1))
    !  enddo
    !enddo
    !rho(:,:,j)=rho(:,:,j)+matmul(U,matmul(rho_ad,transpose(U)))

    r_avg=sum(si_adiab(:,1)*si_adiab(:,1)*x_dvr)
    !if(state==1.and.r_avg<1.4d-10) pop(1,j)=pop(1,j)+1
    if(state==1.and.x(2)<0.14d-10) pop(1,j)=pop(1,j)+1


    !pop(:,j)=pop(:,j)+si_adiab(:,state)**2
    !pop_surf(:,j)=pop_surf(:,j)+si_adiab(:,state)**2
    !ci_diab=matmul(si_adiab,ci)
    !pop_amp(:,j)=pop_amp(:,j)+cdabs(ci_diab)**2
    !do j1=2,nquant
    !  do i1=1,j1-1
    !    pop(:,j)=pop(:,j)+2*real(ci(i1)*dconjg(ci(j1)))*si_adiab(:,i1)*si_adiab(:,j1)
    !  enddo
    !enddo
  endif

  !call cpu_time(t2)
  !tim_coll=tim_coll+t2-t1

end subroutine average
!-----------------------------------------------------------------  

subroutine average_end
  implicit none

end subroutine average_end
!-----------------------------------------------------------------  

subroutine save_old_state
  implicit none

  x_old=x
  v_old=v
  acc_old=acc
  ci_old=ci
  state_old=state
  !ci2_old=ci2
  si_adiab_prev=si_adiab
  V_k_old=V_k
  force_old=force
  energy_old=energy
  delr_old=delr
  delp_old=delp

end subroutine save_old_state
!-----------------------------------------------------------------  

subroutine revert_state
  implicit none

  x=x_old
  v=v_old
  state=state_old
  ci=ci_old
  delr=delr_old
  delp=delp_old
  force=force_old
  !ci2=ci2_old
  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine revert_state
!-----------------------------------------------------------------  

subroutine evolve_quantum_small_dtq
  implicit none
  integer i,nstep_qm
  real*8 dtq,dtq1,dtq2
  real*8 V_k_hold(nquant),dVk_dt(nquant)
  real*8 dforce_dt(nquant,nclass)
  complex*16 ci_prev(nquant),dci_dt(nquant)

  call compute_vdotd
  dVk_dt=(V_k-V_k_old)/dtc
  if(icollapse==1) then
    call compute_delH_dels_ad
  endif

  dtq1=0.02/maxval(vdotd)
  dtq2=0.02*hbar/maxval(V_k-sum(V_k)/real(nquant))
  dtq=dtq1
  if(dtq>dtq2)dtq=dtq2

  if(dtq>dtc)dtq=dtc
  nstep_qm=nint(dtc/dtq)
  dtq=dtc/real(nstep_qm)
  hop_prob=0.d0
  hop_prob_net=0.d0
  V_k_hold=V_k
  V_k=V_k_old
  call compute_mat_adiab

  flag_hop=0
  do i=1,nstep_qm
    call compute_hop_prob(dtq)
    if(flag_hop==0)call check_hop(i*dtq)
    call rk4(ci,dtq,dVk_dt)
    if(icollapse==1)call rk4_decoherence(dtq)
  enddo
  !if(icollapse==1)call vv_decoherence(dtc)

  !if(icollapse==1) then
  !  call verlet_decoherence(dtc,W_overlap,V_k_old,dvk_dt)
  !endif

  do i=1,nquant
    if(hop_prob_net(i)<0.d0)hop_prob_net=0.d0
    hop_prob_net(i)=1.d0-dexp(-hop_prob_net(i))
  enddo

end subroutine evolve_quantum_small_dtq
!-----------------------------------------------------------------  

subroutine compute_hop_prob(dtq)
  implicit none
  real*8,intent(in)::dtq
  integer i
  real*8 pr

  do i=1,nquant
    if(i.ne.state) then
      pr=-2*real(ci(i)*dconjg(ci(state)))*vdotd(i,state)
      pr=pr*dtq/cdabs(ci(state))**2
      if(pr<0.d0)pr=0.d0     !!!! CAUTION AMBER CHECK !!!!
      hop_prob(i)=pr
      hop_prob_net(i)=hop_prob_net(i)+pr
    endif
  enddo

end subroutine compute_hop_prob
!-----------------------------------------------------------------  

subroutine check_hop(tim)
  implicit none
  real*8,intent(in)::tim
  integer i
  real*8 rnd,pr

  call random_number(rnd)
  pr=0.d0
  flag_hop=0
  do i=1,nquant
    if(i.ne.state) then
      pr=pr+hop_prob(i)
      if(rnd<pr) then
        state_tentative=i
        flag_hop=1
        exit
      endif
    endif
  enddo

end subroutine check_hop
!-----------------------------------------------------------------  

subroutine rk4(ci,dtq,dVk_dt)
  implicit none
  complex*16,intent(inout)::ci(nquant)
  real*8,intent(in) :: dtq,dVk_dt(nquant)
  complex*16,dimension(1:nquant):: k1,k2,k3,k4

  k1=matmul(mat_adiab,ci)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k2=matmul(mat_adiab,ci+0.5*dtq*k1)
  k3=matmul(mat_adiab,ci+0.5*dtq*k2)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k4=matmul(mat_adiab,ci+dtq*k3)

  ci=ci+dtq/6.d0*(k1+2*k2+2*k3+k4)

end subroutine rk4
!-----------------------------------------------------------------  

subroutine vv_decoherence(dtc)
  implicit none
  integer i
  real*8,intent(in)::dtc
  real*8 delacc_old(nquant,nquant,nclass)

  delr=delr+delp/mass(1)*dtc+0.5*delacc*dtc**2/mass(1)
  delacc_old=delacc
  call compute_delacc
  delp=delp+0.5*(delacc+delacc_old)*dtc

  do i=1,nclass
    delr(:,:,i)=matmul(W_overlap,matmul(delr(:,:,i),W_overlap))
    delp(:,:,i)=matmul(W_overlap,matmul(delp(:,:,i),W_overlap))
  enddo

end subroutine vv_decoherence
!-----------------------------------------------------------------

subroutine compute_delacc
  implicit none
  integer i

  do i=1,nclass
    delacc(:,:,i)=0.5*anti_commute(delF(:,:,i),sigma,0)
  enddo

end subroutine compute_delacc
!-----------------------------------------------------------------  

subroutine rk4_decoherence(dtq)
  implicit none
  real*8,intent(in)::dtq
  complex*16,dimension(2,nquant,nquant,nclass):: kd1,kd2,kd3,kd4,vec

  vec(1,:,:,:)=delr
  vec(2,:,:,:)=delp

  call compute_T_jk(kd1,vec)
  call compute_T_jk(kd2,vec+0.5*dtq*kd1)
  call compute_T_jk(kd3,vec+0.5*dtq*kd2)
  call compute_T_jk(kd4,vec+dtq*kd3)

  vec=vec+dtq/6.d0*(kd1+2*kd2+2*kd3+kd4)
  delr=vec(1,:,:,:)
  delp=vec(2,:,:,:)

end subroutine rk4_decoherence
!-----------------------------------------------------------------  

subroutine compute_T_jk(T_jk,vec)
  implicit none
  complex*16,intent(in):: vec(2,nquant,nquant,nclass)
  complex*16,intent(out):: T_jk(2,nquant,nquant,nclass)
  complex*16 delr(nquant,nquant,nclass),delp(nquant,nquant,nclass)
  complex*16 Tr(nquant,nquant,nclass)
  complex*16 tmp1(nclass),tmp2(nclass)
  integer i
  real*8 t1,t2,t11,t12,t21,t22

  !call cpu_time(t1)

  delr=vec(1,:,:,:)
  delp=vec(2,:,:,:)

!  call cpu_time(t11)
  call compute_T_jk_R(Tr,delr,delp)
!  call cpu_time(t21);tim_check=tim_check+(t21-t11)

  T_jk(1,:,:,:)=Tr
!  call cpu_time(t12)
  call compute_T_jk_P(Tr,delr,delp)
!  call cpu_time(t22);tim_check2=tim_check2+(t22-t12)

  T_jk(2,:,:,:)=Tr

  tmp1=T_jk(1,state,state,:)
  tmp2=T_jk(2,state,state,:)

  do i=1,nquant
    T_jk(1,i,i,:)=T_jk(1,i,i,:)-tmp1
    T_jk(2,i,i,:)=T_jk(2,i,i,:)-tmp2
  enddo

  !call cpu_time(t2)
  !tim_T_jk=tim_T_jk+t2-t1

end subroutine compute_T_jk
!-----------------------------------------------------------------  

subroutine compute_T_jk_R(T_jk,delr,delp)
  !! Eq. 14 of JCP 137, 22A513
  implicit none
  complex*16,intent(in) :: delr(nquant,nquant,nclass),delp(nquant,nquant,nclass)
  complex*16,intent(out) :: T_jk(nquant,nquant,nclass)
  integer i1

  do i1=1,nclass
      T_jk(:,:,i1)=-iota/hbar*commute(pot,delr(:,:,i1),0)+delp(:,:,i1)/mass(i1)
      T_jk(:,:,i1)=T_jk(:,:,i1)-commute(vdotd,delr(:,:,i1),0)
  enddo

end subroutine compute_T_jk_R
!-----------------------------------------------------------------  

subroutine compute_T_jk_P(T_jk,delr,delp)
  !! Eq. 16 of JCP 137, 22A513
  implicit none
  complex*16,intent(in) :: delr(nquant,nquant,nclass),delp(nquant,nquant,nclass)
  complex*16,intent(out) :: T_jk(nquant,nquant,nclass)
  real*8 delF(nquant,nquant,nclass)
  integer i1,j1

  delF=force
  do i1=1,nquant
    delF(i1,i1,:)=delF(i1,i1,:)-force(state,state,:)
  enddo

  do i1=1,nclass
      T_jk(:,:,i1)=-iota/hbar*commute(pot,delp(:,:,i1),0)
      T_jk(:,:,i1)=T_jk(:,:,i1)+0.5*anti_commute(delF(:,:,i1),sigma,0)
      T_jk(:,:,i1)=T_jk(:,:,i1)-commute(vdotd,delp(:,:,i1),0)
  enddo

end subroutine compute_T_jk_P
!-----------------------------------------------------------------  

!subroutine verlet_decoherence(dt,W_mat,V_k0,dvk_dt)
! implicit none
! real*8,intent(in)::
! dt,W_mat(nquant,nquant),V_k0(nquant),dvk_dt(nquant)
! real*8 acc_dec(nquant,nclass),delf(nquant,nclass),temp(nclass)
! complex*16 temp_delr(nquant,nclass),temp_delp(nquant,nclass)
! !complex*16 ci_diab(nquant)
! integer i,j,k
!
! delF=force_old
! temp=delF(state,:)
! do i=1,nquant
!   delF(i,:)=delF(i,:)-temp
!   acc_dec(i,:)=delF(i,:)*cdabs(ci_old(i))**2/mass
! enddo
!
! do i=1,nquant
!   delr(i,:)=delr(i,:)+delp(i,:)/mass*dt+0.5*acc_dec(i,:)*dt**2
!   delp(i,:)=delp(i,:)+0.5*mass*acc_dec(i,:)*dt
! enddo
!
! !ci_diab=cdexp(iota*V_k0*dt/hbar)*cdexp(0.5*iota*dvk_dt*dt**2/hbar)*ci
! !ci_diab=matmul_lap(W_mat,ci_diab)
! delF=0.d0
! do j=1,nquant
!   do k=1,nquant
!     delF(j,:)=delF(j,:)+dabs(W_mat(j,k)**2)*(force(k,:)-force(state,:))
!   enddo
! enddo
! !temp=delF(state,:)
! do i=1,nquant
! !  delF(i,:)=delF(i,:)-temp
! !  !acc_dec(i,:)=delF(i,:)*cdabs(ci_diab(i))**2/mass
!   acc_dec(i,:)=delF(i,:)*cdabs(ci_old(i))**2/mass
! enddo
!
! do i=1,nquant
!   delp(i,:)=delp(i,:)+0.5*mass*acc_dec(i,:)*dt
! enddo
!
! temp_delr=0.d0;temp_delp=0.d0
! do j=1,nquant
!   do k=1,nquant
!     temp_delr(j,:)=temp_delr(j,:)+dabs(W_mat(k,j)**2)*delr(k,:)
!     temp_delp(j,:)=temp_delp(j,:)+dabs(W_mat(k,j)**2)*delp(k,:)
!   enddo
! enddo
! delr=temp_delr
! delp=temp_delp
!
! !do i=1,nclass
! !  delr(:,i)=delr(:,i)-delr(state,i)
! !  delp(:,i)=delp(:,i)-delp(state,i) !enddo
!
!end subroutine verlet_decoherence
!-----------------------------------------------------------------  

subroutine evolve_classical(dt)
  !! Velocity Verlet
  implicit none
  integer i
  real*8,intent(in) :: dt
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(nclass),delta_v(nclass),acc_sav(nclass)
  real*8 t1,t2

  !call cpu_time(t1)

  if(ifriction==0) then
    x=x+v*dt+0.5*acc*dt*dt
    v=v+0.5*acc*dt
    call evaluate_variables(0)
    v=v+0.5*acc*dt
    call evaluate_variables(1)
  endif

  if(ifriction==1) then
    gama_dt=gamma_B*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)

    call stochastic_force(delta_r,delta_v,dt)

    x=x+c1*dt*v+c2*dt*dt*acc+delta_r
    acc_sav=acc
    call evaluate_variables(0)
    v=c0*v+(c1-c2)*dt*acc_sav+c2*dt*acc+delta_v
    call evaluate_variables(1)
  endif

  !call cpu_time(t2);tim_ev_cl=tim_ev_cl+t2-t1

end subroutine evolve_classical
!-----------------------------------------------------------------  

subroutine traj_terminate(iterm)
  implicit none
  integer,intent(out) :: iterm

  iterm=0

end subroutine traj_terminate
!-----------------------------------------------------------------  

subroutine compute_mat_diab
  implicit none
  integer i,j
  real*8 t1,t2

  !call cpu_time(t1)

  mat=0.d0
  do i=1,nbasis
    do j=1,nbasis
      mat(i,j)=-iota/hbar*sum(si_adiab(i,:)*si_adiab(j,:)*V_k(1:nquant))
    enddo
  enddo

  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1

end subroutine compute_mat_diab
!-----------------------------------------------------------------  

subroutine compute_mat_adiab
  implicit none
  integer i,j
  real*8 t1,t2
  real*8 V_avg
  
  !call cpu_time(t1)

  mat_adiab=-vdotd
  V_avg=sum(V_k)/real(nquant)
  do i=1,nquant
    mat_adiab(i,i)=mat_adiab(i,i)-iota/hbar*(V_k(i)-V_avg)
  enddo
      
  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1
  
end subroutine compute_mat_adiab
!-----------------------------------------------------------------  

subroutine hop
  implicit none
  integer ifrust

  if(flag_hop==1) then
    call velocity_adjust(state_tentative,ifrust)
  endif

end subroutine hop
!-----------------------------------------------------------------  

subroutine velocity_adjust(state_tentative,ifrust)
  implicit none
  integer,intent(in)::state_tentative
  integer,intent(out)::ifrust
  real*8 gij,gama,aa,bb,cc,discr,dp(nclass),vd,f1,f2
  integer i,j,k,kp

  k=state;kp=state_tentative
  cc=V_k(state)-V_k(state_tentative)

  call compute_dij_2state(x,k,kp,dp)
  dp=dp/dsqrt(sum(dp*dp))

  aa=0.d0
  bb=0.d0
  do i=1,nclass

    aa=aa+0.5/mass(i)*(dp(i)*dp(i))
    bb=bb+(v(i)*dp(i))

  enddo

  discr=bb**2+4*aa*cc
  if(discr<0.d0) then
    ifrust=1
    cnt_frust=cnt_frust+1.d0
    if(flag_frust==0)then
      gama=0.d0
      call compute_delH_dels_ad
      f1=sum(force(k,k,:)*dp)
      f2=sum(force(kp,kp,:)*dp)
      vd=sum(v*dp)
      !! reverse velocity based on Truhlar's ideas
      if(f1*f2<0.d0.and.vd*f2<0.d0) then
      !if(f1*f2<0.d0) then
        gama=bb/aa
      endif
    endif
    if(flag_frust>0)gama=0.d0
  else
    ifrust=0
    if(bb>=0.d0) gama=(bb-dsqrt(discr))/(2*aa)
    if(bb<0.d0)  gama=(bb+dsqrt(discr))/(2*aa)
    state=state_tentative
    delr=0.d0
    delp=0.d0
  endif

  do i=1,nclass
    v(i)=v(i)-gama*dp(i)/mass(i)
  enddo

!write(20,*)curr_time*1.d15,dp/dsqrt(sum(dp*dp)),x(1),ifrust
!write(21,*)curr_time*1.d15,k,kp,gama

  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine velocity_adjust
!-----------------------------------------------------------------  

subroutine collapse(dt,iflag_coll)
  implicit none
  real*8,intent(in) :: dt
  integer,intent(out) :: iflag_coll
  real*8 rnd,gama_collapse,gama_reset
  complex*16 su1
  integer n,i,j

  i=state

  if(icollapse==1) then

    iflag_coll=0
    do n=1,nquant
      if(n.ne.state) then
        gama_reset=0.d0
        gama_reset=gama_reset+sum((force(n,n,:)-force(i,i,:))*dble(delr(n,n,:)))/(2*hbar)
        su1=0.d0
        su1=su1+sum(force(i,n,:)*delr(n,n,:))
        gama_collapse=gama_reset-2/hbar*cdabs(su1)
        gama_collapse=gama_collapse*dt
        gama_reset=-gama_reset*dt
        call random_number(rnd)

        if(rnd<gama_collapse) then
          iflag_coll=1
          cnt_collapse=cnt_collapse+1
          if(icollapse==1) then
            !do j=1,nquant
            !  if(j.ne.n) ci(j)=ci(j)/dsqrt(1-cdabs(ci(n)**2))
            !enddo
            !! Erratum: Landry, Subotnik JCP 137, 229901 (2012)
            ci(i)=ci(i)/cdabs(ci(i))*dsqrt(cdabs(ci(i))**2+cdabs(ci(n))**2)
            ci(n)=0.d0

          endif
        endif
        if(rnd<gama_collapse.or.rnd<gama_reset) then
          if(icollapse==1) then
            do j=1,nquant
              delr(j,n,:)=0.d0;delr(n,j,:)=0.d0
              delp(j,n,:)=0.d0;delp(n,j,:)=0.d0
            enddo
          endif
        endif
      endif
    enddo

  endif

end subroutine collapse
!-----------------------------------------------------------------  

subroutine write_output(n,nflag)
  !! nflag=0: Writes various variables as a function of time
  !! nflag=1: writes minimal useful information at the start and end of
  !trajectory
  implicit none
  integer,intent(in)::nflag,n
  integer i
  real*8 t1,t2
  real*8 phase

  !call cpu_time(t1)

  if(nflag==0) then
    if(iwrite==1) then
      if(mod(n,nstep_write)==1.or.nstep_write==1) then
        write(10,'(4es17.7,i5)')curr_time*1.d15,energy/wave_to_J,sum(cdabs(ci)**2),r_avg,state
        write(11,'(es15.5$)')curr_time*1.d15
        write(12,'(5f15.5)')curr_time*1.d15,cdabs(ci(1:2))**2
        write(13,'(5es15.5)')curr_time*1.d15,vdotd(1,2),dasin(W_overlap(1,2))/dtc,hop_prob_net(3-state),state*1.d0
        write(14,'(6f15.5)')curr_time*1.d15,W_overlap(1,1:2),W_overlap(2,1:2),determinant(W_overlap,nquant)
        write(15,'(6es15.5)')curr_time*1.d15,delr(1,1,1)*1.d10,delr(2,2,1)*1.d10
        do i=1,nclass
          write(11,'(2es15.5$)')x(i)*1.d10,v(i)
        enddo
        write(11,*)
      endif
    endif
  endif

  if(nflag==1) then
    if(iwrite==0)then
      write(10,'(5es15.5)')traj_num,energy/wave_to_J,sum(cdabs(ci)**2),temperature
      write(11,*) traj_num
      write(11,'(es15.5$)')curr_time*1.d15
      do i=1,nclass
        write(11,'(2es15.5$)')x(i)*1.d10,v(i)
      enddo
      write(11,*)
      write(11,*)
    endif
    if(iwrite==1) then
      write(10,*)"traj num=",traj_num
      write(10,*)"standard deviation=", dsqrt((ensq_avg-en_avg**2/dfloat(nst_av))/dfloat(nst_av))/wave_to_J
      write(10,*)"ci**2=",sum(cdabs(ci)**2)
      write(10,*);write(10,*)
      write(11,*);write(11,*)
      write(12,*);write(12,*)
      write(13,*);write(13,*)
      write(14,*);write(14,*)
      write(15,*);write(15,*)
    endif
  endif

  !call cpu_time(t2)
  !tim_wr_out=tim_wr_out+t2-t1

end subroutine write_output
!-----------------------------------------------------------------  

subroutine write_average
  !! Writes the final useful output
  implicit none
  integer i,j,i1,k
  real*8 nf,pop_el(3)

  nf=dfloat(n_traj)
  cnt_frust=cnt_frust/nf
  cnt_collapse=cnt_collapse/nf

  pop=pop/nf
  rho=rho/nf
  pop_surf=pop_surf/nf
  pop_amp=pop_amp/nf

  do i=1,nsteps/nstep_avg
    write(100,'(21f15.7)')(i-1)*nstep_avg*dtc*1.d15,pop(1,i)
  enddo

  write(101,*) cnt_frust,cnt_collapse

end subroutine write_average
!-----------------------------------------------------------------  

subroutine evaluate_variables(flag)
  implicit none
  integer,intent(in):: flag
  integer i,j

  if(flag==0) then
    !! position dependant variables only
    call tise
    do i=1,nquant
      do j=1,nquant
        sigma(i,j)=ci(i)*dconjg(ci(j))
      enddo
    enddo
  endif

  if(flag==1) then
    KE_en=0.d0
    do i=1,nclass
      KE_en=KE_en+0.5*mass(i)*v(i)*v(i)
    enddo

    energy=pot_en+KE_en
    !temperature=2*KE_en/(nclass*kb)

    !vdotd=0.d0
    !do i=1,nclass
    !  vdotd=vdotd+v(i)*d_ij(:,:,i)
    !enddo
    !call compute_vdotd
    
  endif

end subroutine evaluate_variables
!-----------------------------------------------------------------  

subroutine tise
  !! time independent schrodinger equation
  !! Output - pot_en,acc
  !! Output - V_k,d_ij
  implicit none
  integer i,j,k
  real*8 Hamil(nbasis,nbasis),ens(nbasis),vect(nbasis,nquant)
  real*8 pot_cl,acc_cl(nclass),acc_qm(nclass),dpotcl_dx(nclass)
  real*8 si_adiab_old(nquant,nbasis)
  real*8 t1,t2

  !call cpu_time(t1)

  call compute_potential(Hamil,delH_dels)
  Hamil_diab=Hamil
  call diag(Hamil,nbasis,ens,vect,nquant)

  do i=1,nquant
    si_adiab(:,i)=vect(:,i)
    if(sum(si_adiab(:,i)*si_adiab_prev(:,i))<0.d0)si_adiab(:,i)=-si_adiab(:,i)
  enddo

  do i=1,nclass
    delH_dels_ad(state,state,i)=sum(si_adiab(:,state)*matmul(delH_dels(:,:,i),si_adiab(:,state)))
  enddo
  !call cpu_time(t2);tim_diag=tim_diag+(t2-t1)

  !call cpu_time(t1)

  !call potential_classical(pot_cl,dpotcl_dx)
  pot_cl=0.d0
  dpotcl_dx=0.d0
  acc_qm=-1.d0/mass*delH_dels_ad(state,state,:)
  acc_cl=-1.d0/mass*dpotcl_dx

  pot_en=pot_cl+ens(state)
  V_k=pot_cl+ens(1:nquant)
  acc=acc_cl+acc_qm

  !call cpu_time(t2);tim_cl=tim_cl+(t2-t1)

end subroutine tise
!-----------------------------------------------------------------  

subroutine compute_delH_dels_ad
  implicit none
  integer i,k,kp,i1

  force=0.d0
  pot=0.d0
  do k=1,nquant
    do kp=k,nquant
      do i=1,nclass
        delH_dels_ad(k,kp,i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,kp)))
      enddo
      force(k,kp,:)=-delH_dels_ad(k,kp,:)
      force(kp,k,:)=-delH_dels_ad(k,kp,:)
      delH_dels_ad(kp,k,i)=delH_dels_ad(k,kp,i)
    enddo
    pot(k,k)=V_k(k)
  enddo

  delF=force
  do i1=1,nquant
    delF(i1,i1,:)=delF(i1,i1,:)-force(state,state,:)
  enddo

end subroutine compute_delH_dels_ad
!-----------------------------------------------------------------  

subroutine compute_dij
  implicit none
  integer i,k,kp

  do k=1,nquant-1
    do kp=k+1,nquant
      do i=1,nclass
        d_ij(k,kp,i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,kp)))
      enddo
      d_ij(k,kp,:)=d_ij(k,kp,:)/(V_k(kp)-V_k(k))
      d_ij(kp,k,:)=-d_ij(k,kp,:)
    enddo
  enddo

end subroutine compute_dij
!-----------------------------------------------------------------  

subroutine compute_dij_2state(x_hop,k,kp,dp)
  implicit none
  integer,intent(in):: k,kp
  real*8,intent(in):: x_hop(nclass)
  real*8,intent(out):: dp(nclass)
  real*8 x_sav(nclass)
  integer i

  x_sav=x
  x=x_hop
  call evaluate_variables(0)

  do i=1,nclass
    dp(i)=sum(si_adiab(:,k)*matmul(delH_dels(:,:,i),si_adiab(:,kp)))
  enddo
  dp=dp/(V_k(kp)-V_k(k))

  x=x_sav
  call evaluate_variables(0)

end subroutine compute_dij_2state
!-----------------------------------------------------------------  

subroutine compute_vdotd
  !! T matrix computation
  implicit none
  integer i,j,k
  real*8,dimension(nquant,nquant) :: W,ci_W,si_W
  real*8 A,B,C,D,E
  real*8 Wlj,Wlk

  !Method 
  !call compute_dij
  !vdotd=0.d0
  !do i=1,nclass
  !  vdotd=vdotd+v(i)*d_ij(:,:,i)
  !enddo

  !Method 2
  ! Meek, Levine, JPCL 5, 2351 (2014). Look at Supp info.
!  do j=1,nquant
!    do k=1,nquant
!      W(j,k)=sum(si_adiab_prev(:,j)*si_adiab(:,k))
!      ci_W(j,k)=dacos(W(j,k))
!      si_W(j,k)=dasin(W(j,k))
!    enddo
!  enddo
!
!  vdotd=0.d0
!  do k=1,nquant-1
!    do j=k+1,nquant
!      A=-sinx_x(ci_W(j,j)-si_W(j,k))
!      B=sinx_x(ci_W(j,j)+si_W(j,k))
!      C=sinx_x(ci_W(k,k)-si_W(k,j))
!      D=sinx_x(ci_W(k,k)+si_W(k,j))
!      Wlj=dsqrt(1.d0-W(j,j)**2-W(k,j)**2)
!      if(Wlj==0.d0.or.nquant==2) then
!        E=0.d0
!      else
!        Wlk=(-W(j,k)*W(j,j)-W(k,k)*W(k,j))/Wlj
!        E=2*dasin(Wlj)/(dasin(Wlj)**2-dasin(Wlk)**2)
!        E=E*(Wlj*Wlk*dasin(Wlj)+dasin(Wlk)*(dsqrt((1-Wlj**2)*(1-Wlk**2))-1.d0))
!      endif
!      vdotd(k,j)=0.5/dtc*(ci_W(j,j)*(A+B)+si_W(k,j)*(C+D)+E)
!      vdotd(j,k)=-vdotd(k,j)
!    enddo
!  enddo

  !Method 3
  do i=1,nquant
    do j=1,nquant
      W_overlap(i,j)=sum(si_adiab_prev(:,i)*si_adiab(:,j))
    enddo
  enddo

  if(flag_ortho==1)call orthoganalize(W_overlap,nquant)
  call logm(W_overlap,vdotd,nquant)
  vdotd=vdotd/dtc

end subroutine compute_vdotd
!-----------------------------------------------------------------  

subroutine orthoganalize(mat,n)
  integer,intent(in)::n
  real*8,intent(inout)::mat(n,n)
  real*8 S_mat(n,n)

  S_mat=matmul(transpose(mat),mat)
  call inverse_squareroot(S_mat,n)
  mat=matmul(mat,S_mat)

end subroutine orthoganalize
!-----------------------------------------------------------------  

subroutine setup_parameters
  implicit none
  integer :: i 

  r0=1.43d-10;l0=0.125d-10
  e_alpha_C(1)=-0.5d0
  e_alpha_C(2)=0.5d0
  e_alpha_C(3)=0.d0
  e_alpha_i(1)=-1.d0
  e_alpha_i(2)=0.5d0
  e_alpha_i(3)=0.5d0
  e_alpha_C=e_alpha_C*1.602d-19
  e_alpha_i=e_alpha_i*1.602d-19

    a_solute=11.2d10
  b_solute=7.1d13*kcal_to_J
  c_solute=0.776
  Diss_A=110.d0*kcal_to_J
  d_A=0.95d-10
  d_b=0.97d-10
  n_A=9.26d10
  n_B=11.42d10

  omg_s=300.d0*2*pi*clight
  k_coup=-12.d22 * wave_to_J/1.d-20
  gamma_B=omg_s

  omg=omg_s

  ! Mass - classial d.o.fs 
  mass(1) = 36.098*amu2kg
  do i = 2,nclass
    mass(i) = 50.49*amu2kg
  end do 

  ! Values of frequencies of classical vibrations
  omg(1) = 310.d0*2*pi*clight
  omg(2) = 300.d0*2*pi*clight

  if (nclass == 3) then 
    omg(3) = 200.d0*2*pi*clight
  else if (nclass > 3) then  
    call generate_pts_grid(200.d0*2*pi*clight, 700.d0*2*pi*clight, nclass-2, omg(3:))  
  end if 

  k_coup=-15.d22*wave_to_J*1.d20

  a_a=11.2*1.d10  
  !! Parameters needed for calculation of Potential VHAB
  b_b=7.1d13*kcal_to_J
  c_c=0.776
  d_A=0.95*1.d-10
  d_B=0.97*1.d-10
  DA=110*kcal_to_J              !! kcal/mol to J
  n_A=9.26*1.d10
  n_B=11.42*1.d10

  x_min(1)=2.68d-10             !! x(1) RAB
  x_min(2)=0.08d-10             !! x(2) Solvent co-ordinate
  !x_min=0.d0

  ll=0.125*1.d-10
  !! parameters needed to calculate charge for dipole
  r0=1.43*1.d-10
  e_A_cov=-0.5d0*1.602d-19!elementary_charge_to_coul
  e_A_ion=-1.d0*1.602d-19!elementary_charge_to_coul
  e_B_cov=0.d0
  e_B_ion=0.5d0*1.602d-19!elementary_charge_to_coul

  r_cutoff=1.4d-10
  gamma_B=omg(2)

  qr=1.d-10
  qp=1.6d-10

  ! RAMANA
  ! Point close to minima of reactant, TS and product - for initiating MD
 
  !Reactant well minima :
  x_r(1) = 2.69d-10
  x_r(2) = 0.09d-10
  do i = 3,nclass
    if (i <= nclass) then 
      x_r(i) = 0.0d-10
    end if 
  end do

  !Product well minima :
  x_p(1) = 2.74d-10
  x_p(2) = 0.37d-10
  do i = 3,nclass
    if (i <= nclass) then
      x_p(i) = 0.0d-10
    end if
  end do

  !Transition State minima :
  x_ts(1) = 2.70d-10
  x_ts(2) = 0.12d-10
  do i = 3,nclass
    if (i <= nclass) then
      x_ts(i) = 0.0d-10
    end if
  end do

  !RAMANA
  call setup_pred_model(q_grid, theta_c, test_coords, train_coords) !gives model_params

  call setup_quantized_vib

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine setup_quantized_vib
  implicit none
  integer i
  real*8 delq,mass_H

  mass_H=1.d0*amu2kg

  do i=1,n_dvr
    x_dvr(i)=0.8d-10+1.2d-10*(i-1)/real(n_dvr-1)
  enddo
  delq=x_dvr(2)-x_dvr(1)

  call compute_KE_matrix_dvr(KE_dvr,n_dvr,delq,mass_H)

end subroutine setup_quantized_vib
!-----------------------------------------------------------------  

subroutine compute_potential(H_diab,delV_dels)
  implicit none
  real*8,intent(out) :: H_diab(:,:),delV_dels(:,:,:)
  real*8 pot_AHB,dpot_AHB_dx(nclass)
  real*8 pot_sol,dpot_sol_dx(nclass)
  real*8 V,dv_dq,d2v_dq2,dv_dx(nclass),q
  integer i

  H_diab=KE_dvr
  delv_dels=0.d0

  rAB=x(1)

  do i=1,n_dvr
    rAH=x_dvr(i)
    q=x_dvr(i)

    !call full_pot(q,x,V,dv_dx)
    call predicted_pot(q,x,V,dv_dx)

    H_diab(i,i)=H_diab(i,i)+V
    delv_dels(i,i,:)=dv_dx
  enddo

end subroutine compute_potential
!-----------------------------------------------------------------  

subroutine full_pot(q,x,V,dv_dx)
  real*8,intent(in) :: q,x(nclass)
  real*8,intent(out) :: V,dV_dx(nclass)
  real*8 V1,dv1_dx(nclass),dV1_dq,d2V1_dq2
  real*8 V2,dv2_dx(nclass),dV2_dq,d2V2_dq2
  real*8 pot_cl,dpotcl_dx(nclass)
  real*8 pot_unc,dpotunc_dx(nclass)

  call pot_q(V1,dv1_dx,q,x,dV1_dq,d2V1_dq2)
  call pot_coup(V2,dv2_dx,q,x,dV2_dq,d2V2_dq2)             !! V_HAB
  call potential_classical(x, pot_cl,dpotcl_dx)

  ! BATH 
  !call potential_uncoupled(x, pot_unc,dpotunc_dx) ! harmonic bath 
  call pot_morse_uncoupled(x, pot_unc,dpotunc_dx) ! morse bath 
 
  V = V1 + V2 + pot_cl + pot_unc
  dV_dx = dV1_dx + dv2_dx + dpotcl_dx + dpotunc_dx

end subroutine full_pot
!-----------------------------------------------------------------  

subroutine scan_pot
  implicit none
  real*8 pot_AHB,dpot_AHB_dx(nclass)
  real*8 pot_sol,dpot_sol_dx(nclass)
  real*8 V,dv_dq,d2v_dq2,dv_dx(nclass)
  integer i
  real*8 V1,dv1_dx(nclass),q,dV1_dq,d2V1_dq2
  real*8 V2,dv2_dx(nclass),dV2_dq,d2V2_dq2
  real*8 pot_cl,dpotcl_dx(nclass)

  rAB=x(1)

  do i=1,100
    rAH=0.8d-10+1.2d-10*i/100.d0
    q=rAH
!    call compute_charges(rAH)
!    call comp_pot_AHB(pot_AHB,dpot_AHB_dx,rAH,x(1))
!    call comp_pot_sol(pot_sol,dpot_sol_dx)
  call pot_q(V1,dv1_dx,q,x,dV1_dq,d2V1_dq2)
  call pot_coup(V2,dv2_dx,q,x,dV2_dq,d2V2_dq2)             !! V_HAB
  call potential_classical(x, pot_cl,dpotcl_dx)

  write(40,*) q*1.d10,(V1+V2+pot_cl)/wave_to_J
  enddo

end subroutine scan_pot
!-----------------------------------------------------------------  

subroutine compute_charges(rAH)
  implicit none
  real*8,intent(in)::rAH
  real*8 f,r

  f(r)=0.5d0*(1+(r-r0)/dsqrt((r-r0)**2+l0**2))

  charge(1)=e_alpha_c(1)+f(rAH)*(e_alpha_i(1)-e_alpha_c(1))
  charge(2)=e_alpha_c(3)+f(rAH)*(e_alpha_i(3)-e_alpha_c(3))
  charge_H=e_alpha_c(2)+f(rAH)*(e_alpha_i(2)-e_alpha_c(2))

  mu=-charge(1)*rAH+charge(2)*(rAB-rAH)
  dmu_drAB = charge(2)

end subroutine compute_charges
!-----------------------------------------------------------------

subroutine comp_pot_AHB(pot,dpot_dx,rAH,rAB)
  implicit none
  real*8,intent(out):: pot,dpot_dx(nclass)
  real*8,intent(in):: rAH,rAB
  real*8 morse1,dM1_dr
  real*8 morse2,dM2_dr
  real*8 dV_drAB,dV_drAH
  real*8 tmp

  call morse(morse1,dM1_dr,rAH,Diss_A,n_A,d_A)
  call morse(morse2,dM2_dr,(RAB-rAH),c_solute*Diss_A,n_B,d_B)
  tmp=b_solute*dexp(-a_solute*rAB)

  dV_drAB=tmp*(-a_solute)+dM2_dr

  pot=morse1+morse2+tmp
  dpot_dx=0.d0
  dpot_dx(1)=dv_drAB
  dpot_dx(2)=0.d0

end subroutine comp_pot_AHB
!-----------------------------------------------------------------

subroutine morse(pot,dpot_dr,r,Diss_A,n_A,d_A)
  implicit none
  real*8,intent(out)::pot,dpot_dr
  real*8,intent(in)::r,diss_A,n_A,d_A
  real*8 tmp

  tmp=Diss_A*dexp(-n_A*(r-d_A)**2/(2*r))

  pot=Diss_A-tmp
  dpot_dr=tmp*n_A*0.5*(1-(d_A/r)**2)

end subroutine morse
!-----------------------------------------------------------------  

subroutine comp_pot_sol(pot,dpot_dx)
  implicit none
  real*8,intent(out) :: pot,dpot_dx(nclass)

  pot = 0.5*mass(2)*omg_s**2*x(2)**2 + k_coup*mu*x(2)
  dpot_dx(1)=k_coup*x(2)*dmu_drAB
  dpot_dx(2)=mass(2)*omg_s**2*x(2) + k_coup*mu

end subroutine comp_pot_sol
!-----------------------------------------------------------------

!subroutine potential_classical(pot_cl,acc_cl)
!  implicit none
!  real*8,intent(out) :: pot_cl,acc_cl(nclass)
!  integer n
!  real*8 q1,q3
!
!  pot_cl=0.d0!0.5*mass(1)*sum(omg**2*x**2)
!  acc_cl=0.d0!mass(1)*(omg**2*x)
!
!end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine check_acceleration
  !! A test subroutine that compares analytical accelerations with
  !numerical
  !accelerations
  implicit none
  integer i,nflag
  real*8 delx,en_old,acc_sav(nclass)
  real*8 q0,rnd

  qr=1.d-10
  qp=1.6d-10

  delx=1.d-17
  state=1

  do i=1,nclass
    call random_number(rnd)
    x(i)=(rnd*2-1.d0)*1.d-10
  enddo

  x = 0.d0  

  x(1)=2.74d-10
  x(2)=0.13d-10
  
  call evaluate_variables(0)
  en_old=pot_en;acc_sav=acc

  write(16,*) "delx=",delx
  write(16,*)

  do i=1,nclass
      x(i)=x(i)+delx
      call evaluate_variables(0)
      acc(i)=-(pot_en-en_old)/delx/mass(i)
      write(16,*)"Analytical acceleration =",acc_sav(i)
      write(16,*)"Numerical acceleration  =",acc(i)
      write(16,*)"Error =",(acc(i)-acc_sav(i))/acc(i)*100.d0
      write(16,*)
      x(i)=x(i)-delx
  enddo

  stop

end subroutine check_acceleration
!---------------------------------------------------------- 

subroutine stochastic_force(delr,delv,dt)
  !! stoachastic forces for langevin equation
  !! Not used for the Holstein model results 
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass),delv(nclass)!f(nclass)
  integer i
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma_B*dt

  do i=1,nclass
    sig_r=dt*dsqrt(kb*temperature/mass(i)*1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
    sig_v=dsqrt(kb*temperature/mass(i)*(1-dexp(-2*gdt)))
    sig_rv=(dt*kb*temperature/mass(i)* 1.d0/gdt*(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

    call gaussian_random_number(rnd1)
    call gaussian_random_number(rnd2)
    delr(i)=sig_r*rnd1
    delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------  

function commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  complex*16 tmp
  integer j,k

  if(iflag==0) commute=matmul(A,B)-matmul(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do j=1,nquant
      do k=1,nquant
        commute(j,k)=B(j,k)*(A(j,j)-A(k,k))
      enddo
    enddo
  endif

  if(iflag==2) then
    !! Assume A is tridiagonal, with a_ii=0, and a_ij=-a_ji (a is
    !assumed to be d_ij)
    do j=1,nquant
      do k=1,nquant
        tmp=0.d0
        if(j<nquant) tmp=tmp+A(j,j+1)*B(j+1,k)
        if(j>1) tmp=tmp-A(j-1,j)*B(j-1,k)
        if(k>1) tmp=tmp-A(k-1,k)*B(j,k-1)
        if(k<nquant) tmp=tmp+A(k,k+1)*B(j,k+1)
      enddo
    enddo
  endif

end function commute
!-----------------------------------------------------------------  

function anti_commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 anti_commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  integer i,j

  if(iflag==0) anti_commute=matmul(A,B)+matmul(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do i=1,nquant
      do j=1,nquant
       anti_commute(i,j)=B(i,j)*(A(i,i)+A(j,j))
      enddo
    enddo
  endif

end function anti_commute
!-----------------------------------------------------------------  

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and
  !eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork
  !are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n), eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 

subroutine logm(mat,log_mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(in):: mat(n,n)
  real*8,intent(out):: log_mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
  enddo

  log_mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine logm
!-----------------------------------------------------------------  

subroutine inverse_squareroot(mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(inout):: mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=1.d0/t(i,i)**0.5d0
  enddo

  mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine inverse_squareroot
!-----------------------------------------------------------------  

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and
  !eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork
  !are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n
  integer,intent(inout) :: nold
  complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
  real*8,intent(in) :: mat(n,n)
  complex*16,intent(out) :: T(n,n)
  complex*16,allocatable,intent(inout):: cwork(:)
  real*8 rwork(n)
  complex*16 mat_c(n,n)

  integer lwork
  logical:: select
  logical bwork(n)
  integer sdim,info,AllocateStatus

  T=mat

  info=0
  sdim=0

  if(nold.ne.n .or. .not.allocated(cwork)) then
  !if(nold.ne.n) then
    lwork=-1
    if(allocated(cwork))deallocate(cwork)
    allocate(cwork(n))
    call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
    lwork=int(cwork(1))
    deallocate(cwork)
    allocate(cwork(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(cwork)
  call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine schur
!---------------------------------------------------------- 

REAL FUNCTION determinant(matrix, n)
    !!http://web.hku.hk/~gdli/UsefulFiles/Example-Fortran-program.html
    IMPLICIT NONE
    REAL*8, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL*8 :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                determinant= 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    determinant= l
    DO i = 1, n
        determinant= determinant* matrix(i,i)
    END DO
    
END FUNCTION determinant
!-----------------------------------------------------------------  

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!---------------------------------------------------------- 

subroutine compute_KE_matrix_dvr(KE_matrix,ngrid,delq,mass)
  !! computes KE matrix in DVR basis
  !! Appendix A of JCP 96, 1982 (1991)
  implicit none
  integer,intent(in) :: ngrid       !! size of DVR grid
  real*8,intent(inout) :: KE_matrix(ngrid,ngrid)
  real*8,intent(in) :: delq         !! step size in position of DVR basis
  real*8,intent(in) :: mass         !! mass
  integer i,j
  real*8 pi,hbar

  pi=dacos(-1.d0);hbar=1.05457266D-34

  KE_matrix=hbar**2/(2*mass*delq**2)
  do i=1,ngrid
    do j=1,ngrid
      KE_matrix(i,j)=KE_matrix(i,j)*(-1.d0)**(i-j)
      if(i==j) then
        KE_matrix(i,j)=KE_matrix(i,j)*pi**2/3.d0
      else
        KE_matrix(i,j)=KE_matrix(i,j)*2.d0/real(i-j)**2
      endif
    enddo
  enddo
end subroutine compute_KE_matrix_dvr
!---------------------------------------------------------- 

subroutine draw_pes
  implicit none
  integer i,j,k

  x(1)=2.6d-10
  x(2)=0.0d-10
  !  call scan_pot
  k = 2
  x = 0.d0

  open(20,file='pes_vib.out')
  do j=1,30
    x(1)=2.6d-10+0.3d-10*j/30.d0
    do i=1,100
      x(k)=-0.1d-10+0.6d-10*i/100.d0
      call evaluate_variables(0)
      write(20,*) x(1)*1.d10, x(k)*1d10, V_K(1)/wave_to_J
    enddo
    write(20,*)
  enddo
  close(20)

  write(6,*) "Vibrational PES files generated"
  write(6,*) '------------------------------------------------'

  stop

end subroutine draw_pes
!-----------------------------------------------------------------  

subroutine draw_pes_1ES
  implicit none
  integer i,j

  x(1)=2.6d-10
  x(2)=0.0d-10
  !  call scan_pot

  open(31,file='pes_vib_ES.out')
  do j=1,30
    x(1)=2.6d-10+0.3d-10*j/30.d0
    do i=1,100
      x(2)=-0.1d-10+0.6d-10*i/100.d0
      call evaluate_variables(0)
      write(31,*) x(1)*1.d10, x(2)*1d10, V_K(2)/wave_to_J
    enddo
    write(31,*)
  enddo
  close(31)
  stop

end subroutine draw_pes_1ES
!--------------------------------------------------------------------

!! From Simran
!! Potential and its first and second derivative....!!

subroutine exact_potential(q,x,V,dV_dq,d2V_dq2)
  implicit none
  real*8,intent(in):: q,x(nclass)
  real*8,intent(out)::V,dV_dq,d2V_dq2
  real*8 V1, V2, pot_cl,dV1_dq, dV2_dq,d2V1_dq2,d2V2_dq2
  integer i,j,k
  real*8 dv1_dx(nclass), dv2_dx(nclass),dpotcl_dx(nclass)

  call pot_q(V1,dv1_dx,q,x,dV1_dq,d2V1_dq2)
  call pot_coup(V2,dv2_dx,q,x,dV2_dq,d2V2_dq2)             !! V_HAB
  call potential_classical(x, pot_cl,dpotcl_dx)
  V=V1+V2+pot_cl
  dV_dq=dV1_dq+dV2_dq
  d2V_dq2=d2V1_dq2+d2V2_dq2

end subroutine exact_potential
!-----------------------------------------------------------------------------------------

subroutine pot_q(V,dv_dx,q,x,dV_dq,d2V_dq2)             !! V_HAB
  implicit none
  real*8,intent(in)::q,x(nclass)                        !!q=r, x(1)=RAB
  real*8,intent(out)::V,dv_dx(nclass),dV_dq,d2V_dq2
  real*8 fac1, fac2, fac3
  real*8 f1,f2,f3,f4,f5,f6

  fac1=b_b*dexp(-a_a*x(1))
  fac2=DA*(1.d0-dexp((-n_A*(q-d_A)**2)/(2.d0*q)))
  fac3=c_c*DA*(1.d0-dexp((-n_B*(x(1)-q-d_B)**2)/(2.d0*(x(1)-q))))
  V=fac1+fac2+fac3

  f1=(n_A*(q-d_A)*(q+d_A))/(2.d0*q**2)
  f2=dexp((-n_A*(q-d_A)**2)/(2.d0*q))
  f3=DA*f1*f2
  f4=(n_B*(x(1)-q-d_B)*(x(1)-q+d_B))/(2.d0*(x(1)-q)**2)
  f5=dexp((-n_B*(x(1)-q-d_B)**2)/(2.d0*(x(1)-q)))
  f6=-c_c*DA*f4*f5

  dV_dq=f3+f6                                  !! First derivative
  dv_dx = 0.d0
  dv_dx(1)=-a_a*fac1+(c_c*DA-fac3)*n_B*(x(1)-q-d_B)*(x(1)-q+d_B)/(2.d0*(x(1)-q)**2)
  dv_dx(2)=0.d0
  dv_dx(3)=0.d0

  call second_derivative_potq(q,x,d2V_dq2)      !! Second derivative

end subroutine pot_q
!-----------------------------------------------------------------------------------------

subroutine pot_coup(V,dv_dx,q,x,dV_dq,d2V_dq2)
  implicit none
  real*8,intent(in)::q, x(nclass)
  real*8,intent(out)::V,dv_dx(nclass),dV_dq,d2V_dq2
  real*8 dipole, e_A, e_B,de_dq(2),d2e_dq2(2)
  call evaluate_charge(q,x,e_A,e_B,de_dq,d2e_dq2)
  dipole=-e_A*q + e_B*(x(1)-q)
  V=k_coup*dipole*x(2)!**3/1.d-20 !!_!! 
  dV_dq=k_coup*x(2)*(-e_A-(q*de_dq(1))-e_B+((x(1)-q)*de_dq(2)))    !!First derivative
  d2V_dq2=k_coup*x(2)*((-2.d0*de_dq(1))-(2.d0*de_dq(2))-(q*d2e_dq2(1))+((x(1)-q)*d2e_dq2(2))) !! Second derivative

  dv_dx = 0.d0
  dv_dx(1) = k_coup*e_B*x(2)!**3/1.d-22 
  dv_dx(2) = k_coup*dipole!*3.d0*x(2)**2/1.d-22
  dv_dx(3)=0.d0

end subroutine pot_coup
!-----------------------------------------------------------------  

subroutine potential_classical(x, pot_cl,dv_dx)
  implicit none
  real*8, intent(in) :: x(nclass)
  real*8,intent(out) :: pot_cl,dv_dx(nclass)

  real*8 :: pot_cl_x2, pot_cl_x1, a, De

  De = 10000.d0*wave_to_J ! Dissociation energy
  a = dsqrt(mass(2)*(omg(2)**2)*0.5d0/De) 

  pot_cl_x2 = 0.5d0*mass(2)*(omg(2)**2)*(x(2)**2) ! harmonic x2
  !pot_cl_x2 = De*(1.d0 - dexp(-a*x(2)))**2 ! morse x2
  pot_cl_x1 = 0.5d0*mass(1)*(omg(1)**2)*((x(1)-2.8d-10)**2)
  pot_cl = pot_cl_x2 + pot_cl_x1
  dv_dx = 0.d0
  dv_dx(1)=mass(1)*(omg(1)**2)*(x(1)-2.8d-10)
  dv_dx(2)=mass(2)*(omg(2)**2)*x(2) ! harmonic x2
  !dv_dx(2)= 2.d0*De*a*dexp(-a*x(2))*(1.d0-dexp(-a*x(2))) ! morse x2
  dv_dx(3)=0.d0

end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine potential_uncoupled(x, pot_unc,dv_dx)
  implicit none
  real*8, intent(in) :: x(nclass)
  real*8,intent(out) :: pot_unc,dv_dx(nclass)
  integer :: i 

  pot_unc =0.d0
  dv_dx = 0.d0
  dv_dx(1)=0.d0
  dv_dx(2)=0.d0
  if (nclass > 2) then
  do i = 3,nclass 
      pot_unc = pot_unc + 0.5d0*mass(i)*(omg(i)**2)*(x(i)**2)
      dv_dx(i)=mass(i)*(omg(i)**2)*x(i)
  end do
  end if  
end subroutine potential_uncoupled
!!------------------------------------------------------------------

subroutine pot_morse_uncoupled(x, pot_unc, dv_dx)
  implicit none 
  real*8, intent(in) :: x(nclass)
  real*8, intent(out) :: pot_unc, dv_dx(nclass)

  integer :: i 
  real*8 :: a(nclass), De

  De = Diss_A !50000.d0*wave_to_J ! dissociation energy
  pot_unc =0.d0
  dv_dx = 0.d0
  dv_dx(1)=0.d0
  dv_dx(2)=0.d0
  if (nclass > 2) then
  do i = 3,nclass
      a(i) = dsqrt(mass(i)*(omg(i)**2)*0.5d0/De) 
      pot_unc = pot_unc + De*(1.d0 - dexp(-a(i)*x(i)))**2
      dv_dx(i)= 2.d0*De*a(i)*dexp(-a(i)*x(i))*(1.d0-dexp(-a(i)*x(i)))
  end do
  end if  

end subroutine pot_morse_uncoupled
!!-----------------------------------------------------------------
subroutine second_derivative_potq(q,x,d2V_dq2)
  implicit none
  real*8,intent(in):: q, x(nclass)
  real*8,intent(out)::d2V_dq2
  real*8 fac1, fac2,fac3
  real*8 fac4,fac7,fac8,fac9,fac10,f1,f2,f3
  real*8 fac11,fac12,fac6

  fac1=(-n_A*(q-d_A)**2)/(2.d0*q)
  fac2=(n_A*d_A**2)/q**3
  f1=n_A*n_A*(q**2-d_A**2)**2
  f2=4.d0*q**4
  f3=f1/f2
  fac3=fac2-f3
  fac6=DA*(dexp(fac1)*fac3)

  fac4=(-n_B*(x(1)-q-d_B)**2)/(2.d0*(x(1)-q))
  fac8=(-n_B*d_B**2)/(x(1)-q)**3
  fac9=n_B*n_B*((x(1)-q)**2-d_B**2)**2
  fac10=fac9/(4.d0*(x(1)-q)**4)

  fac12=-c_c*DA*(fac8+fac10)*dexp(fac4)
  d2V_dq2=fac12+fac6

 end subroutine second_derivative_potq
!-----------------------------------------------------------------  

subroutine evaluate_charge(q,x,e_A,e_B,de_dq,d2e_dq2)
  implicit none
  real*8, intent(in):: q, x(nclass)
  real*8, intent(out):: e_A, e_B,de_dq(2),d2e_dq2(2)
  real*8 num, den, func_q, dfunc_dq,d2func_dq2
  real*8 de_A,de_B,d2eA_dq2,d2eB_dq2
  num=q-r0
  den=((q-r0)**2+ll**2)
  func_q=0.5d0*(1.d0+(num/dsqrt(den)))
  dfunc_dq=0.5d0*(ll**2/(den*dsqrt(den)))
  d2func_dq2=-1.5d0*((ll*ll*(q-r0))/(den**2*dsqrt(den)))
  e_A=(1-func_q)*e_A_cov + func_q*e_A_ion                    !!e_A_ion/e_B_ion, e_B_cov/e_B_ion are constants
  e_B=(1-func_q)*e_B_cov + func_q*e_B_ion
  de_A=dfunc_dq*(-e_A_cov+e_A_ion)
  de_B=dfunc_dq*(-e_B_cov+e_B_ion)
  d2eA_dq2=d2func_dq2*(-e_A_cov+e_A_ion)
  d2eB_dq2=d2func_dq2*(-e_B_cov+e_B_ion)
  de_dq(1)=de_A
  de_dq(2)=de_B
  d2e_dq2(1)=d2eA_dq2
  d2e_dq2(2)=d2eB_dq2
end subroutine evaluate_charge

!---------------------------------------------------------------------------------  
!! RAMANA
!! Subroutines for safeguards, metrics and plotting related to predicted potential 
!----------------------------------------------------------------------------------

subroutine test_coeffs
  implicit none

  real*8 :: flattened(ndims)
  real*8 :: coeff_mat(nclass,nclass), eigen_value(nclass)
  real*8 :: eigen_vect(nclass,nclass), pot_output
  integer :: k, i, j

  coeff_mat = 0.d0

  do i = 1,ndims
    call rbf_model_predictions(q_grid(15), q_grid, theta_c(i,:), flattened(i))
  end do

  k = 1
  do i = 1,nclass
    do j = i,nclass
      coeff_mat(i,j) = flattened(k)
      coeff_mat(j,i) = flattened(k)
      k = k+1
    end do
  end do

  call diag(coeff_mat, nclass, eigen_value, eigen_vect, nclass)

  pot_output = sum(eigen_value)

  if (pot_output <= 0) then
    write(*,*) "------------------------------------------------"
    write(*,*) "Possible Dissociative Potential, Increase points"
    write(*,*) "------------------------------------------------"
    !stop
  end if

end subroutine test_coeffs
!----------------------------------------------------------------------

subroutine rmse_model
  implicit none

  real*8 :: predicted_pots_test(n_q_grid, test_points), dummy(nclass)
  real*8 :: actual_pots_test(n_q_grid, test_points)
  real*8 :: predicted_pots_train(n_q_grid, train_points)
  real*8 :: actual_pots_train(n_q_grid, train_points)
  real*8 :: error_square_sum_test, rmse_test
  real*8 :: error_square_sum_train, rmse_train
  real*8 :: absolute_error_test, mae_test
  real*8 :: absolute_error_train, mae_train
  integer :: i,j

  error_square_sum_test = 0.d0 ! initialising value 
  error_square_sum_train = 0.d0 ! initialising value 
  absolute_error_test = 0.d0 ! initialising value 
  absolute_error_train = 0.d0 ! initialising value 

  

  do i = 1,n_q_grid
    do j = 1,test_points
      call predicted_pot(q_grid(i), test_coords(i,j,:), predicted_pots_test(i,j), dummy)
      call full_pot(q_grid(i), test_coords(i,j,:), actual_pots_test(i,j), dummy)
      error_square_sum_test = error_square_sum_test + (predicted_pots_test(i,j)-actual_pots_test(i,j))**2
      absolute_error_test = absolute_error_test + dabs(predicted_pots_test(i,j)-actual_pots_test(i,j))
      end do
  end do

  do i = 1,n_q_grid
    do j = 1,train_points
      call predicted_pot(q_grid(i), train_coords(i,j,:), predicted_pots_train(i,j), dummy)
      call full_pot(q_grid(i), train_coords(i,j,:), actual_pots_train(i,j),dummy)
      error_square_sum_train = error_square_sum_train + (predicted_pots_train(i,j)-actual_pots_train(i,j))**2
      absolute_error_train = absolute_error_train + dabs(predicted_pots_train(i,j)-actual_pots_train(i,j))
   end do
  end do

  rmse_test = dsqrt(error_square_sum_test/real(test_points*n_q_grid))
  rmse_train = dsqrt(error_square_sum_train/real(train_points*n_q_grid))

  mae_test = absolute_error_test/real(test_points*n_q_grid)
  mae_train = absolute_error_train/real(train_points*n_q_grid)

  write(6,*) '------------------------------------------------'
  write(6,*) 'Classical D.O.F.s =', nclass
  write(6,*) '------------------------------------------------'
  write(6,*) 'MD time run =', md_total_t
  write(6,*) '------------------------------------------------'
  write(6,*) 'Training Points =', train_points*n_q_grid
  write(6,*) 'Testing  Points =', test_points*n_q_grid
  write(6,*) '------------------------------------------------'
  write(6,*) 'Training RMSE =', rmse_train/wave_to_J, "cm^-1"
  write(6,*) 'Testing  RMSE =', rmse_test/wave_to_J, "cm^-1"
  write(6,*) '------------------------------------------------'
  !write(6,*) 'Training MAE  =', mae_train/wave_to_J, "cm^-1"
  !write(6,*) 'Testing  MAE  =', mae_test/wave_to_J, "cm^-1"
  !write(6,*) '------------------------------------------------' 

end subroutine rmse_model
!----------------------------------------------------------------------------

subroutine draw_elec_pot
  implicit none

  integer :: i, j
  real*8 :: q_grid(100), V_pred, V_orig, dV_dx(nclass), x(nclass)
  real*8 :: x_grid(100, nclass)

  call generate_pts_grid(0.7d-10, 2.d-10, 100, q_grid)
  open(77,file='pes_elec.out')
  do i=1,100
    call predicted_pot(q_grid(i), x_ts, V_pred, dV_dx)
    call full_pot(q_grid(i), x_ts, V_orig, dV_dx)
    write(77,*) q_grid(i)*1.d10, V_orig/wave_to_J, V_pred/wave_to_J
  enddo
  close(77)

  x_grid = 0.d0
  call generate_pts_grid(2.6d-10, 2.9d-10, 100, x_grid(:,1))
  open(87, file='pes_elec_x1.out')
  do i=1,100 
    call predicted_pot(2.d-10, x_grid(i,:), V_pred, dV_dx)
    call full_pot(2.d-10, x_grid(i,:), V_orig, dV_dx)
    write(87,*) x_grid(i,1)*1.d10, V_orig/wave_to_J, V_pred/wave_to_J
  end do
  close(87)

  j = 2
  x_grid = 0.d0
  x_grid(:,1) = 2.7d-10
  call generate_pts_grid(-0.5d-10, 1d-10, 100, x_grid(:,j)) ! 0.1 to 0.5d-10
  open(88, file='pes_elec_xj.out')
  do i=1,100
    call predicted_pot(2.d-10, x_grid(i,:), V_pred, dV_dx)
    call full_pot(2.d-10, x_grid(i,:), V_orig, dV_dx)
    write(88,*) x_grid(i,j)*1.d10, V_orig/wave_to_J, V_pred/wave_to_J
  end do
  close(88)

  write(6,*) "Electronic PES files generated"
  write(6,*) '------------------------------------------------'

  stop 

end subroutine draw_elec_pot

!-----------------------------------------------------------------------
!!! subroutines below inserted by RAMANA
!!! subroutines below are for defining the new_algo prediction potential 
!-------------------------------------------------------------------

subroutine predicted_pot(q, x, V, dV_dx)
  implicit none 

  real*8, intent(in) :: q, x(nclass)
  real*8, intent(out) :: V, dV_dx(nclass)

  integer :: i, j 
  real*8 :: x_scaled(nclass), basis_coeffs(ndims)

  x_scaled = x*1.d10 ! making the input dimensionless 

  ! Bath coupling energy - coeffs of bath energy 
  do i = 1,ndims
    call rbf_model_predictions(q, q_grid, theta_c(i,:), basis_coeffs(i))
  end do

  !do i = 1,ndims
  !  basis_coeffs(i) = theta_c(i,1)
  !end do  

  ! pot_energy and derivative   
  call ho_model_predictions(x_scaled, basis_coeffs, V, dV_dx) 
  
end subroutine predicted_pot
!---------------------------------------------------------------

subroutine setup_pred_model(q_grid, theta_c, test_coords, train_coords)
  implicit none

  real*8, intent(out) :: q_grid(n_q_grid)
  real*8, intent(out) :: theta_c(ndims, n_q_grid)
  real*8, intent(out) :: test_coords(n_q_grid, test_points, nclass)
  real*8, intent(out) :: train_coords(n_q_grid, train_points, nclass)

  real*8, allocatable :: coords_x(:,:,:), pots_x(:,:) 
  real*8 :: theta_B(n_q_grid, ndims), test_pots(n_q_grid,test_points)
  real*8 :: train_pots(n_q_grid,train_points)
  integer :: i, n, n_ts, st, en, m

  n = md_total_t/md_sampling_t
  n_ts = md_total_t/(md_sampling_t*2)
  m = 2*(n+1)+(n_ts+1)

  allocate(coords_x(n_q_grid,m,nclass), pots_x(n_q_grid,m))

  ! Calling setup parameters to get model coeffs 
  call model_setup_all_q(m, q_grid, coords_x, pots_x)

  ! Performing train_test split 
  call train_test_split(coords_x, pots_x, train_coords, test_coords, train_pots, test_pots)

  ! Getting coefficients of ho_bath  
  call ho_model_optimised_params(train_coords*1.d10, train_pots, theta_B) 
 
  do i = 1,ndims  
    call rbf_model_optimised_params(q_grid, theta_B(:,i), theta_c(i,:))
    !theta_c(i,:) = theta_B(:,i)
  end do  
 
  deallocate(coords_x, pots_x)

end subroutine setup_pred_model 
!-------------------------------------------------------------------

subroutine model_setup_all_q(m, q_grid, coords_x, pots_x)
  implicit none 

  integer, intent(in) :: m ! is the total number of points from MD  
  real*8, intent(out) :: q_grid(n_q_grid) 
  real*8, intent(out) :: coords_x(n_q_grid,m,nclass)
  real*8, intent(out) :: pots_x(n_q_grid,m)
 
  integer :: i,j

  call generate_pts_grid(0.7d-10, 2.d-10, n_q_grid, q_grid)
  
  !call generate_pts_grid(1.3275862068965519d-10,1.3275862068965519d-10,n_q_grid,q_grid)
 
  !q_grid(1) = 1.3275862068965519d-10 

  do i = 1,n_q_grid
    call model_setup_given_q(m, q_grid(i), coords_x(i,:,:),pots_x(i,:))  
  end do  

end subroutine model_setup_all_q
!----------------------------------------------------------------

subroutine model_setup_given_q(m, q, coords_needed_all, pots_needed_all)
  implicit none

  integer, intent(in) :: m ! is the total number of points from MD
  real*8, intent(in) :: q
  real*8, intent(out) :: coords_needed_all(m,nclass)
  real*8, intent(out) :: pots_needed_all(m)

  real*8, allocatable :: coords_needed_r(:,:), coords_needed_p(:,:)
  real*8, allocatable :: coords_needed_ts(:,:)
  real*8, allocatable :: pots_needed_r(:), pots_needed_p(:), pots_needed_ts(:) 
  integer :: i, j, n, n_ts
 
  n = md_total_t/md_sampling_t
  n_ts = md_total_t/(md_sampling_t*2)

  call md(q, x_r, md_sampling_t, md_total_t, coords_needed_all(1:n+1,:), pots_needed_all(1:n+1))
  call md(q, x_p, md_sampling_t, md_total_t, coords_needed_all(n+2:2*n+2,:), pots_needed_all(n+2:2*n+2))
  call md_ts(q, x_ts, md_sampling_t, md_total_t, coords_needed_all(2*n+3:m,:), pots_needed_all(2*n+3:m))

end subroutine model_setup_given_q
!-----------------------------------------------------------------

subroutine generate_pts_grid(start_pt, end_pt, num_pt, grid)
  implicit none

  real*8, intent(in) :: start_pt, end_pt
  integer, intent(in) :: num_pt
  real*8, intent(out) :: grid(num_pt)

  integer :: i
  real*8 :: step

  step = (end_pt - start_pt)/(num_pt-1)

  do i = 1,num_pt
    grid(i) = start_pt + (i-1)*step
  end do

end subroutine generate_pts_grid
!-------------------------------------------------------------------

subroutine rbf_model_predictions(q, q_train, theta_rbf, C_q_pred)
  implicit none

  real*8, intent(in) :: q, q_train(n_q_grid), theta_rbf(n_q_grid)
  real*8, intent(out) :: C_q_pred

  real*8 :: kernelised_q(n_q_grid), gamma_param, std_q
  integer :: i

  call standard_deviation(q_train, std_q)
  gamma_param = 1.d0/(2.d0*(std_q**2.d0))

  do i = 1,n_q_grid
    kernelised_q(i) = dexp(-1.d0*gamma_param*((q - q_train(i))**2))
  end do

  C_q_pred = dot_product(theta_rbf, kernelised_q)

end subroutine rbf_model_predictions
!------------------------------------------------------------------

subroutine rbf_model_optimised_params(q_train, c_train, theta_rbf)
  implicit none 

  real*8, intent(in) :: q_train(n_q_grid), c_train(n_q_grid)
  real*8, intent(out) :: theta_rbf(n_q_grid)

  real*8 :: kernel(n_q_grid,n_q_grid), std_q, gamma_param 
  real*8 :: k2(n_q_grid, n_q_grid), k2_inv(n_q_grid, n_q_grid) 
  real*8 :: identity(n_q_grid, n_q_grid), alpha_rbf, k2_inv_k(n_q_grid, n_q_grid)
  integer :: i, j

  call standard_deviation(q_train, std_q)
  gamma_param = 1.d0/(2.d0*(std_q**2))

  alpha_rbf = 1.d-10 !find a way to automate this choice 

  ! forming the kernel
  do i = 1,n_q_grid
    do j = 1,n_q_grid
      kernel(i,j) = dexp(-1.d0*gamma_param*(q_train(i)-q_train(j))**2) 
    end do 
  end do 

  ! Making an identity matrix
  identity = 0.d0
  do i = 1,n_q_grid 
    identity(i,i) = 1.d0
  end do  

  k2 = matmul(transpose(kernel), kernel) - alpha_rbf*identity
  call sq_matrix_inverse(k2, k2_inv, n_q_grid)
  k2_inv_k = matmul(k2_inv, transpose(kernel)) 
  theta_rbf = matmul(k2_inv_k, c_train)   

end subroutine rbf_model_optimised_params 
!--------------------------------------------------------------------

subroutine ho_model_predictions(x, theta_Vc, V_c_pred, dVc_dx)
  implicit none

  real*8, intent(in) :: x(nclass), theta_Vc(ndims)
  real*8, intent(out) :: V_c_pred, dVc_dx(nclass)

  real*8 :: kernelised_x(ndims)
  real*8 :: der_k(nclass,ndims)
  integer :: i 

  call ho_kernelise(x, kernelised_x)
  call ho_kernelise_derivative(x, der_k)
  
  V_c_pred = dot_product(theta_Vc, kernelised_x)

  do i=1,nclass
    dVc_dx(i) = dot_product(theta_Vc, der_k(i,:))
  end do 

  dVc_dx = dVc_dx/1d-10 ! converting acceleration from unitless to SI 

end subroutine ho_model_predictions
!-------------------------------------------------------------------

subroutine ho_model_optimised_params(train_coords, train_pots, theta_ho)
  implicit none 

  real*8, intent(in) :: train_coords(n_q_grid,train_points,nclass) 
  real*8, intent(in) :: train_pots(n_q_grid,train_points)
  real*8, intent(out) :: theta_ho(n_q_grid,ndims)

  real*8, allocatable :: kernel(:,:,:), k2(:,:), k2_inv(:,:), k2_inv_k(:,:) 
  real*8, allocatable :: identity(:,:)
  real*8 :: idid(ndims, ndims)
  real*8 :: alpha_ho, logk2(ndims, ndims), logk2_inv(ndims,ndims)
  integer :: i, j, p   

  allocate(kernel(n_q_grid,train_points,ndims), identity(ndims,ndims))
  allocate(k2(ndims,ndims), k2_inv(ndims, ndims), k2_inv_k(ndims, train_points)) 

  !alpha_ho = 1.d-50

  ! Making an identity matrix
  identity = 0.d0
  do i = 1,ndims
    identity(i,i) = 1.d0
  end do

  ! Finding optimised coefficients 
  do i = 1,n_q_grid
    do j = 1,train_points
      call ho_kernelise(train_coords(i,j,:), kernel(i,j,:))
    end do 
    k2 = matmul(transpose(kernel(i,:,:)), kernel(i,:,:)) !- alpha_ho*identity
    call sq_matrix_inverse(k2, k2_inv, ndims)
    k2_inv_k = matmul(k2_inv, transpose(kernel(i,:,:))) 
    theta_ho(i,:) = matmul(k2_inv_k, train_pots(i,:))
  end do 

  deallocate(kernel, k2, k2_inv, k2_inv_k)

end subroutine ho_model_optimised_params
!--------------------------------------------------------------------

subroutine ho_kernelise(x_input, ker_ut_flat)
  implicit none

  real*8, intent(in) :: x_input(nclass)
  real*8, intent(out) :: ker_ut_flat(ndims)

  real*8 :: ker_mat(nclass+1,nclass+1), x_appended(nclass+1,1)
  integer :: i, j, k

  x_appended(:nclass,1) = x_input
  x_appended(nclass+1,1) = 1.d0   

  ker_mat = matmul(x_appended,transpose(x_appended))

  k = 1
  do i = 1,nclass+1
    do j = 1,nclass+1
      if (i <= j) then 
        ker_ut_flat(k) = ker_mat(i,j)
        k = k+1
      end if
    end do
  end do

end subroutine ho_kernelise
!----------------------------------------------------------------

subroutine ho_kernelise_derivative(x_input, der_k)
  implicit none

  real*8, intent(in) :: x_input(nclass)
  real*8, intent(out) :: der_k(nclass, ndims)

  real*8 :: ker_mat(nclass+1,nclass+1), x_appended(nclass+1,1)
  integer :: i, j, k, l

  x_appended(:nclass,1) = x_input
  x_appended(nclass+1,1) = 1.d0
  ker_mat = matmul(x_appended, transpose(x_appended)) 

  do l = 1,nclass
    k = 1
    do i = 1,nclass+1
      do j = 1,nclass+1
        if (i<=j) then
          if(i==l .and. j==l) then
            der_k(l,k) = 2.d0*x_appended(l,1)
          else if(i==l .or. j==l) then
            if (x_appended(l,1) == 0.d0) then 
              der_k(l,k) = 0.d0
            else 
              der_k(l,k) = ker_mat(i,j)/x_appended(l,1)
            end if 
          else 
            der_k(l,k) = 0.d0
          end if
          k = k+1 
        end if
      end do
    end do 
  end do

end subroutine ho_kernelise_derivative 
!--------------------------------------------------------------------

subroutine md(q, x0, md_sampling_t, md_total_t, coords_needed, pots_needed)
  implicit none 

  real*8, intent(in) :: q, x0(nclass), md_total_t, md_sampling_t
  real*8, intent(out) :: coords_needed(int(md_total_t/md_sampling_t)+1,nclass)
  real*8, intent(out) :: pots_needed(int(md_total_t/md_sampling_t)+1)
   
  integer :: iterations, n, step_n, i, j 
  real*8 :: v_mean(nclass), v_sd(nclass), V, dV(nclass), x_sd(nclass)
  real*8, allocatable :: coords_all(:,:), pots_all(:), acc_all(:,:), vel_all(:,:)

  do i = 1,nclass
    v_mean(i) = 0.d0 
    v_sd(i) = dsqrt(kb*temperature/mass(i))
    x_sd(i) = dsqrt(kb*temperature/(mass(i)*(omg(i)**2)))
  end do  
 
  iterations = md_total_t/dtc
  n = md_total_t/md_sampling_t
  step_n = iterations/n

  allocate(coords_all(iterations+1, nclass), pots_all(iterations+1))
  allocate(acc_all(iterations+1, nclass), vel_all(iterations+1, nclass))
  
  call normal_distribution(v_mean, v_sd, nclass, vel_all(1,:))
  call normal_distribution(x0, x_sd, nclass, coords_all(1,:))

  call acceleration(q, x0, acc_all(1,:)) 

  do i = 1,iterations 
    call positions_friction(coords_all(i,:), vel_all(i,:), acc_all(i,:), dtc, coords_all(i+1,:))
    call acceleration(q, coords_all(i+1,:), acc_all(i+1,:))
    call velocity_friction(vel_all(i,:), acc_all(i,:), acc_all(i+1,:), dtc, vel_all(i+1,:))
  end do

  coords_needed(1,:) = x0
  call full_pot(q, x0, pots_needed(1), dV)
  do j = 1,n 
    coords_needed(j+1,:) = coords_all((j*step_n)+1,:)
    call full_pot(q, coords_needed(j+1,:), pots_needed(j+1), dV)
  end do   

  deallocate(coords_all, pots_all, acc_all, vel_all)

end subroutine md 
!-----------------------------------------------------------------

subroutine md_ts(q, x0, md_sampling_t, md_total_t, coords_needed, pots_needed)
  implicit none

  real*8, intent(in) :: q, x0(nclass), md_total_t, md_sampling_t
  real*8, intent(out) :: coords_needed(int(md_total_t/(md_sampling_t*2))+1,nclass)
  real*8, intent(out) :: pots_needed(int(md_total_t/(md_sampling_t*2))+1) 

  integer :: iterations, n, step_n, i, j
  real*8 :: v_mean(nclass), v_sd(nclass), x_sd(nclass), V, dV(nclass)
  real*8, allocatable :: coords_all(:,:), pots_all(:), acc_all(:,:), vel_all(:,:)

  do i = 1,nclass
    v_mean(i) = 0.d0
    v_sd(i) = dsqrt(kb*temperature/mass(i))
    x_sd(i) = dsqrt(kb*temperature/(mass(i)*(omg(i)**2)))
  end do 

  iterations = md_total_t/dtc
  n = md_total_t/(md_sampling_t*2)
  step_n = iterations/n

  allocate(coords_all(iterations+1, nclass), pots_all(iterations+1))
  allocate(acc_all(iterations+1, nclass), vel_all(iterations+1, nclass))

  call normal_distribution(v_mean, v_sd, nclass, vel_all(1,:))
  call normal_distribution(x0, x_sd, nclass, coords_all(1,:))

  vel_all(1,2) = 0.d0
  coords_all(1,2) = x0(2)

  call acceleration_ts(q, x0, acc_all(1,:))

  do i = 1,iterations
    call positions_friction(coords_all(i,:), vel_all(i,:), acc_all(i,:), dtc,coords_all(i+1,:))
    call acceleration_ts(q, coords_all(i+1,:), acc_all(i+1,:))
    call velocity_friction(vel_all(i,:), acc_all(i,:), acc_all(i+1,:), dtc, vel_all(i+1,:))
  end do

  coords_needed(1,:) = x0
  call full_pot(q, x0, pots_needed(1), dV)
  do j = 1,n
    coords_needed(j+1,:) = coords_all(j*step_n+1,:)
    call full_pot(q, coords_needed(j+1,:), pots_needed(j+1), dV)
  end do

  deallocate(coords_all, pots_all, acc_all, vel_all)
 
end subroutine md_ts
!------------------------------------------------------------------

subroutine acceleration(q, x, acc)
  implicit none 

  real*8, intent(in) :: q, x(nclass)
  real*8, intent(out) :: acc(nclass)

  integer :: i 
  real*8 :: V, dV_dx(nclass) 
 
  ! function to call dV 
  call full_pot(q,x,V,dV_dx)
  
  do i = 1,nclass
    acc(i) = -1.d0*dV_dx(i)/mass(i)
  end do  

end subroutine acceleration
!-----------------------------------------------------------------

subroutine acceleration_ts(q, x, acc_ts)
  implicit none

  real*8, intent(in) :: q, x(nclass)
  real*8, intent(out) :: acc_ts(nclass)

  integer :: i
  real*8 :: V, dV_dx(nclass)

  ! function to call dV 
  call full_pot(q,x,V,dV_dx)

  ! needed for contrained md where x(2) stays constant (TS)
  do i = 1,nclass
    if (i ==2) then  
      acc_ts(i) = 0 
    else
      acc_ts(i) = -1.d0*dV_dx(i)/mass(i)
    end if 
  end do
  
end subroutine acceleration_ts
!-----------------------------------------------------------------

subroutine velocity(v1, a1, a2, dt, v2)
  implicit none 

  real*8, intent(in) :: v1(nclass), a1(nclass), a2(nclass), dt
  real*8, intent(out) :: v2(nclass)

  v2 = v1 + ((a1+a2)*dt)/2.d0 

end subroutine velocity
!-----------------------------------------------------------------

subroutine velocity_friction(v1, a1, a2, dt, v2)
  implicit none 

  real*8, intent(in) :: v1(nclass), a1(nclass), a2(nclass), dt
  real*8, intent(out) :: v2(nclass)

  real*8 :: gama_dt, c0, c1, c2, delta_r(nclass), delta_v(nclass) 

  gama_dt=gamma_B*dt
  c0=dexp(-gama_dt)
  c1=1.d0/gama_dt*(1.d0-c0)
  c2=1.d0/gama_dt*(1.d0-c1)

  call stochastic_force(delta_r,delta_v,dt)

  v2=c0*v1+(c1-c2)*dt*a1+c2*dt*a2+delta_v

end subroutine velocity_friction 
!-----------------------------------------------------------------

subroutine positions(x1, v1, a1, dt, x2)
  implicit none 

  real*8, intent(in) :: x1(nclass), v1(nclass), a1(nclass), dt
  real*8, intent(out) :: x2(nclass)
  
  x2 = x1 + v1*dt + (a1*dt*dt)/2.d0 

end subroutine positions
!------------------------------------------------------------

subroutine positions_friction(x1, v1, a1, dt, x2)
  implicit none 

  real*8, intent(in) :: x1(nclass), v1(nclass), a1(nclass), dt
  real*8, intent(out) :: x2(nclass)

  real*8 :: gama_dt, c0, c1, c2, delta_r(nclass), delta_v(nclass) 

  gama_dt=gamma_B*dt
  c0=dexp(-gama_dt)
  c1=1.d0/gama_dt*(1.d0-c0)
  c2=1.d0/gama_dt*(1.d0-c1)

  call stochastic_force(delta_r,delta_v,dt)

  x2=x1+c1*dt*v1+c2*dt*dt*a1+delta_r

end subroutine positions_friction
!------------------------------------------------------------
subroutine normal_distribution(mean, sd, num, output)
  implicit none

  real*8, intent(in) :: mean(num), sd(num)
  integer, intent(in) :: num
  real*8, intent(out) :: output(num)

  real*8 :: pre_output(num)
  integer :: i 

  do i = 1,num
    call gaussian_random_number(pre_output(i))
    output(i) = pre_output(i)*sd(i) + mean(i)
  end do

end subroutine normal_distribution
!------------------------------------------------------------------

subroutine sq_matrix_inverse(M, M_inv, n)
  ! NOTE : use only for symetric matrices 
  implicit none

  integer, intent(in) :: n
  real*8, intent(in) :: M(n,n)
  real*8, intent(out) :: M_inv(n,n)

  integer :: i 
  real*8 :: eigen_value(n), eigen_vect(n,n), D(n,n), diag_inv(n,n)

  D = M
  call diag(D, n, eigen_value, eigen_vect, n)

  diag_inv = 0.d0

  do i = 1,n
    if (eigen_value(i) /= 0.d0) then
      diag_inv(i,i) = 1/eigen_value(i) 
    else
      write(6,*) "Eigenvalue found 0 while inverting" 
      stop 
    end if 
  end do

  M_inv = matmul(eigen_vect, matmul(diag_inv, transpose(eigen_vect)))

end subroutine sq_matrix_inverse
!-------------------------------------------------------------------

subroutine standard_deviation(array_in, sd)
  implicit none

  real*8, intent(in) :: array_in(n_q_grid)
  real*8, intent(out) :: sd
  real*8 :: mean, diff_sum
  integer :: i

  mean = sum(array_in)/n_q_grid
  diff_sum = 0.0d0
  do i = 1,n_q_grid
      diff_sum = diff_sum + (mean-array_in(i))**2
  end do
  sd = dsqrt(diff_sum/n_q_grid)

end subroutine standard_deviation
!-------------------------------------------------------------------

subroutine train_test_split(coords, pots, train_coords, test_coords, train_pots,test_pots)
  implicit none

  real*8, intent(in) :: coords(:,:,:), pots(:,:)        ! Input arrays of coordinates and potentials
  real*8, intent(out) :: train_coords(:,:,:), test_coords(:,:,:) ! Arrays for trainand test coordinates
  real*8, intent(out) :: train_pots(:,:), test_pots(:,:)   ! Arrays for train andtest potentials

  integer :: n, n_train, n_test, i, idx, q
  real*8 :: temp
  logical :: chosen(size(coords))
  integer :: all_indices(size(coords)), q_index
  integer, allocatable :: train_indices(:), test_indices(:)

  n = size(coords,2)

  n_test = test_points                                ! Set the number of training points
  n_train = n - n_test                                ! Calculate the number of testing points

  allocate(train_indices(n_train))
  allocate(test_indices(n_test))

  call random_seed()                                  ! Initialize random number generator

  DO q_index = 1,n_q_grid

    ! Initialize the indices and chosen array
    chosen = .false.
    all_indices = [(i, i = 1, n)]                      ! Fill indices from 1 to n

    ! Randomly select training indices
    do i = 1, n_train
      do
        call random_number(temp)
        idx = int(temp * n) + 1                        ! Generate a random index between 1 and n
        if (.not. chosen(idx)) then                    ! Check if the point has already been chosen
          chosen(idx) = .true.
          train_indices(i) = all_indices(idx)          ! Store the training index
          exit
        end if
      end do
    end do

    ! Fill the test indices with remaining points
    idx = 1
    do i = 1, n
      if (.not. chosen(i)) then
        test_indices(idx) = all_indices(i)
        idx = idx + 1
      end if
    end do

    ! Allocate and populate the train and test arrays
    train_coords(q_index,:,:) = coords(q_index, train_indices, :)
    train_pots(q_index,:) = pots(q_index, train_indices)
    test_coords(q_index,:,:) = coords(q_index, test_indices, :)
    test_pots(q_index,:) = pots(q_index, test_indices)

  END DO

  deallocate(train_indices)
  deallocate(test_indices)

end subroutine train_test_split
!----------------------------------------------------------------------------------

End Module mod_afssh






















