!------------------------------------------------!
 program laser_analysis
!------------------------------------------------!

 implicit none

integer :: i,ii,run,reference_run,time_bin,max_run,min_run,max_time,istat,&
           marker,temp_cut_100,temp_step,temp_cut_200
CHARACTER(len=32) :: arg
double precision :: time,valuee,shift_value,a_e,pop_ratio,a_t,delta_temp,delta_time,&
                    arrhenius_factor,planck,boltzmann,light_speed,activation_energy,&
                    ln_n_slope,ln_n_y,n_function,slope_cut,wave_number

double precision,allocatable,dimension(:) :: min_value,slope_array,slope2_array,slope3_array,&
                                             log_decay,tau,n_array,tau2
double precision,allocatable,dimension(:,:) ::decay_array

write(*,*)"Initilizing..."

ALLOCATE(min_value(0:999))
ALLOCATE(slope_array(0:99999))
ALLOCATE(slope2_array(0:99999))
ALLOCATE(slope3_array(0:99999))
ALLOCATE(decay_array(0:999,0:20000))
ALLOCATE(log_decay(0:20000))
ALLOCATE(tau(0:20000))
ALLOCATE(tau2(0:20000))
ALLOCATE(n_array(0:99999))


 DO i = 1, iargc()
    CALL getarg(i, arg)
  END DO

read (arg,'(I10)') temp_step
write(*,*)"temp_step=",temp_step


arrhenius_factor=0.0
planck=6.626*(10.0**(-34.0))
boltzmann=1.38*(10.0**(-23.0))
light_speed=299792458.0
activation_energy=500.0

min_value(:)=999999.9999
decay_array(:,:)=0.0
log_decay(:)=0.0
slope_array(:)=0.0
slope2_array(:)=0.0
slope3_array(:)=0.0
tau(:)=0.0
tau2(:)=0.0
n_array(:)=0.0

write(*,*)"Read Input"
reference_run=0
time_bin=0
min_run=99999999
max_run=0
max_time=0
delta_temp=0.0


OPEN(UNIT=1746, FILE='total_laser.dat',STATUS='OLD')
 do
   READ(1746,FMT=*,IOSTAT=istat)run,time,valuee
   IF( istat < 0) EXIT


   if(run .LT. min_run)then
      min_run=run
   endif
   if(run .GT. max_run)then
      max_run=run
   endif
   if(max_time .LT. time_bin)then
      max_time=time_bin
   endif

   if(reference_run .NE. run)then
      time_bin=0
   endif
   
   if(valuee .LT. min_value(run))then
      min_value(run)=valuee
   endif
   
   decay_array(run,time_bin)=valuee
   time_bin=time_bin+1
enddo

CLOSE (1746,STATUS='KEEP',IOSTAT=I)

write(*,*)"Beginning Analysis"

reference_run=0
time_bin=0
slope_cut=(0.7*temp_step)-1
open(135,file='log_decay.dat')
open(134,file='decay_array.dat')
open(136,file='slope_array.dat')
open(137,file='slope2_array.dat')
open(138,file='tau.dat')
OPEN(UNIT=1746, FILE='total_laser.dat',STATUS='OLD')
 do
   READ(1746,FMT=*,IOSTAT=istat)run,time,valuee
   IF( istat < 0) EXIT


   if(reference_run .NE. run)then
      marker=0
      do i=2,time_bin
         slope3_array(i)=slope2_array(i)-slope2_array(i-1)
         write(137,'(G21.14,2x,I8,2x,G21.14)')(reference_run*temp_step)+303.00,i,slope3_array(i)
         call flush(137)

         if(marker .EQ. 0)then
         if(abs(slope3_array(i)) .GT. 40000.0)then
            marker=i
         endif
         endif

      enddo

      do i=0,marker
         tau(reference_run)=tau(reference_run)+(slope2_array(i))
      enddo
      tau(reference_run)=tau(reference_run)/(marker)
      tau(reference_run)=-1.0/tau(reference_run)
      write(138,'(G21.14,2x,G21.14,2x,I8)')(reference_run*temp_step)+303.00,tau(reference_run),marker
      call flush(138)
      reference_run=run
      slope_array(:)=0.0
      slope2_array(:)=0.0
      slope3_array(:)=0.0
      time_bin=0
   endif

if(time .GE. 0.0)then

   delta_time=0.00001
   if(valuee-min_value(run) .GT. 0.0)then
      write(135,'(G21.14,2x,G21.14,2x,G21.14)')(run*temp_step)+303.00,time_bin,log(valuee-min_value(run))
      call flush(135)
      write(134,'(G21.14,2x,G21.14,2x,G21.14)')(run*temp_step)+303.00,time_bin,valuee-min_value(run)
      call flush(134)
      slope_array(time_bin)=log(valuee-min_value(run))
      if(time_bin .NE. 0)then
         write(136,'(G21.14,2x,G21.14,2x,G21.14)')(run*temp_step)+303.00,time_bin,&
         (slope_array(time_bin)-slope_array(time_bin-1))/delta_time
         call flush(136)
         slope2_array(time_bin)=(slope_array(time_bin)-slope_array(time_bin-1))/delta_time
      endif
   endif

   if(valuee-min_value(run) .EQ. 0.0)then
      write(135,'(G21.14,2x,G21.14,2x,G21.14)')(run*temp_step)+303.00,time_bin,log(0.00000000000000001)
      call flush(135)
      write(134,'(G21.14,2x,G21.14,2x,G21.14)')(run*temp_step)+303.00,time_bin,valuee-min_value(run)
      call flush(134)
      slope_array(time_bin)=log(valuee-min_value(run))
      if(time_bin .NE. 0)then
         write(136,'(G21.14,2x,G21.14,2x,G21.14)')(run*temp_step)+303.00,time_bin,&
         (slope_array(time_bin)-slope_array(time_bin-1))/delta_time
         call flush(136)
         slope2_array(time_bin)=(slope_array(time_bin)-slope_array(time_bin-1))/delta_time
      endif
   endif
   
time_bin=time_bin+1
endif

enddo

CLOSE (1746,STATUS='KEEP',IOSTAT=I)

do i=0,max_run
   temp_cut_100=(i*temp_step)
   if(temp_cut_100 .GE. 40.0)then
      temp_cut_100=i
      exit
   endif
enddo


open(148,file='a_e.dat')
do i=0,temp_cut_100
   a_e=a_e+tau(i)
   write(148,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,1.0/tau(i)
   call flush(148)
enddo

a_e=a_e/(temp_cut_100+1)
a_e=1.0/a_e
write(*,*)"<A_E> (s)=",a_e

do i=0,max_run
   temp_cut_200=(i*temp_step)
   if(temp_cut_200 .GE. 200.0)then
      temp_cut_200=i
      exit
   endif
enddo

a_t=0.0

open(147,file='a_t.dat')
do i=temp_cut_100,temp_cut_200
   pop_ratio=8.11*exp(-3380.00/((i*temp_step)+303.00))
   a_t=a_t+((1+pop_ratio)/(tau(i)*pop_ratio))-(a_e/pop_ratio)
   write(147,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,((1+pop_ratio)/(tau(i)*pop_ratio))-(a_e/pop_ratio)
   call flush(147)
enddo

a_t=a_t/(temp_cut_200-temp_cut_100)

write(*,*)"<A_T> (s)=",a_t

open(145,file='compute.dat')
do i=1,max_run
   pop_ratio=8.11*exp(-3380.00/((i*temp_step)+303.00))
   tau2(i)=((1.0+pop_ratio)/(a_e+(pop_ratio*a_t)))
   write(145,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,tau2(i)
   call flush(145)
enddo

do i=0,max_run
   temp_cut_200=(i*temp_step)
   if(temp_cut_200 .GE. 350.0)then
      temp_cut_200=i
      exit
   endif
enddo

open(139,file='n_temp.dat')
open(146,file='ln_n_temp.dat')
do i=temp_cut_200,max_run-1
   pop_ratio=8.11*exp(-3380.00/((i*temp_step)+303.00))
   n_array(i)=((((1+pop_ratio)/tau(i))-a_e)/pop_ratio)-a_t
   write(139,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,(n_array(i))
   call flush(139)
   write(146,'(G21.14,2x,G21.14)')1.0/((i*temp_step)+303.00),log(n_array(i))
   call flush(146)
enddo


!Find ln(N) Slope
ln_n_slope=0.0
ln_n_slope=(log(n_array(max_run-1))-log(n_array(temp_cut_200)))/&
           ((1.0/(((max_run-1)*temp_step)+303.00))-(1.0/((temp_cut_200*temp_step)+303.00)))

!write(*,*)"log1=",log(n_array(max_run-1))
!write(*,*)"log2=",log(n_array(temp_cut_200))
!write(*,*)"delta1=",(log(n_array(max_run-1))-log(n_array(temp_cut_200)))
!write(*,*)"bottom1=",(1.0/(((max_run-1)*temp_step)+303.00))
!write(*,*)"bottom1=",(1.0/((temp_cut_200*temp_step)+303.00))
!write(*,*)"delta2=",((1.0/(((max_run-1)*temp_step)+303.00))-(1.0/((temp_cut_200*temp_step)+303.00)))
write(*,*)"slope=",ln_n_slope

!Find Intercept
ln_n_y=0.0
ln_n_y=log(n_array(max_run-1))-(ln_n_slope/(((max_run-1)*temp_step)+303.00))

write(*,*)"intercept=",ln_n_y

n_function=0.0

open(141,file='compute2.dat')
open(142,file='pop.dat')
do i=1,max_run
   pop_ratio=8.11*exp(-3380.00/((i*temp_step)+303.00))
   n_function=exp(ln_n_slope/((i*temp_step)+303.00)+ln_n_y)
   write(141,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,((1.0+pop_ratio)/(a_e+(pop_ratio*(a_t+n_function))))
   call flush(141)
   write(142,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,pop_ratio
   call flush(142)
enddo

arrhenius_factor=exp(ln_n_y)

write(*,*)"Arrhenius prefactor=",arrhenius_factor

write(*,*)"Activation Energy (m^-1)=",(-1.0*(ln_n_slope*boltzmann)/(planck*light_speed))

activation_energy=(-1.0*(ln_n_slope*boltzmann)/(planck*light_speed))

wave_number=activation_energy

n_array(:)=0.0
open(143,file='n_arr.dat')
open(150,file='compute3.dat')
open(151,file='compute4.dat')
do i=1,max_run
   pop_ratio=8.11*exp(-3380.00/((i*temp_step)+303.00))
   n_array(i)=arrhenius_factor*exp(-1.0*planck*light_speed*wave_number/(boltzmann*((i*temp_step)+303.00)))
   write(143,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,n_array(i)
   call flush(143)
   write(150,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,((1.0+pop_ratio)/(a_e+(pop_ratio*(a_t+n_array(i)))))
   call flush(150)
   n_array(i)=arrhenius_factor*exp(-1.0*planck*light_speed*500.0/(boltzmann*((i*temp_step)+303.00)))
   write(151,'(G21.14,2x,G21.14)')(i*temp_step)+303.00,((1.0+pop_ratio)/(a_e+(pop_ratio*(a_t+n_array(i)))))
   call flush(151)
enddo

end program laser_analysis
