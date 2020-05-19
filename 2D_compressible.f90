!==========================================================
!  2D Navier-Stokes equation solver by Sylvain Laizet, 2014
!==========================================================

!START OF THE MAIN PROGRAM
!
program navierstokes
!
  implicit none   !-->all the variables MUST be declared
!
  integer,parameter :: nx=513, ny=257 ,nt=40000,ns=1,nf=1,mx=nf*nx,my=nf*ny
  !size of the computational domain (nx x ny)
  !size of the exchanger (mx x my)
  !number of time step for the simulation
!
  !Declaration of variables
  real(8),dimension(nx,ny) :: uuu,vvv,rho,eee,pre,tmp,rou,rov,wz,tuu,tvv
  real(8),dimension(nx,ny) :: roe,tb1,tb2,tb3,tb4,tb5,tb6,tb7,tb8,tb9
  real(8),dimension(nx,ny) :: tba,tbb,fro,fru,frv,fre,gro,gru,grv,gre,rot,eps,ftp,gtp,scp
  real(8),dimension(mx) :: xx
  real(8),dimension(my) :: yy
  real(8),dimension(mx,my) :: tf
  real(8),dimension(2,ns) :: coef
  integer :: i,j,itemp,k,n,nxm,iread,ni,nj,isave,longueur,imodulo
  real(8) :: xlx,yly,CFL,dlx,dx,xmu,xkt,um0,vm0,tm0
  real(8) :: xba,gma,chp,eta,uu0,dlt,um,vm,tm,x,y,dy
  real(8) :: cci,chv,d,roi,tpi
!**************************************************
  character(len=80) path_network,nom_network,path_files,nom_file,nom_script,nom_film,nchamp
  character(len=4) suffix
  character(len=20) nfichier
!*******************************************
  !Name of the file for visualisation:
990 format('vort',I4.4)
  imodulo=250 !snapshots to be saved every imodulo time steps

  ! AB2 temporal scheme itemp=1
  ! RK3 temporal scheme itemp=2
  itemp=1
  
  ! Subroutine for the initialisation of the variables
  call initl(uuu,vvv,rho,eee,pre,tmp,rou,rov,roe,nx,ny,xlx,yly, &
       xmu,xba,gma,chp,dlx,eta,eps,scp,xkt,uu0)

  !we need to define the time step
  dx=xlx/nx !mesh size in x
  dy=yly/ny !mesh sixe in y
  CFL=0.25  !CFL number for time step
  dlt=CFL*dlx
  print *,'The time step of the simulation is',dlt

  !Computation of the average velocity and temperature at t=0
  call average(uuu,um0,nx,ny)
  call average(vvv,vm0,nx,ny)
  call average(scp,tm0,nx,ny)
  write(*,*) 'Average values at t=0', um0,vm0,tm0

!BEGINNING OF TIME LOOP
  do n=1,nt
     if (itemp.eq.1) then   !TEMPORAL SCHEME AB2

       ! Inlet Conditions HardWired

       call param(xlx,yly,xmu,xba,gma,chp,roi,cci,d,tpi,chv,uu0)

        do j=1,ny
          uuu(1,j) = uu0
          vvv(1,j) = 0.01*(sin(4.*acos(-1.)*1/nx) +sin(7.*acos(-1.)*1/nx))*exp(-(j*yly/ny-6)**2)
          tmp(1,j) = tpi
          eee(1,j) = chv*tmp(1,j)*0.5*(uuu(1,j)*uuu(1,j)+vvv(1,j)*vvv(1,j))
          rho(1,j) = roi
          pre(1,j) = rho(1,j)*(gma-1.)/gma*chp*tmp(1,j)
          rou(1,j) = rho(1,j)*uuu(1,j)
          rov(1,j) = rho(1,j)*vvv(1,j)
          roe(1,j) = rho(1,j)*eee(1,j)
          scp(1,j) = 1.
        enddo

        call fluxx(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,tb1,tb2,tb3,tb4, &
             tb5,tb6,tb7,tb8,tb9,tba,tbb,fro,fru,frv,fre,xlx,yly,xmu,xba,eps, &
             eta,ftp,scp,xkt)

        call adams(rho,rou,rov,roe,fro,gro,fru,gru,frv,grv,&
             fre,gre,ftp,gtp,scp,nx,ny,dlt,dlx)

        call etatt(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,gma,chp)

     endif

     if (itemp.eq.2) then !TEMPORAL SCHEME RK3

        !loop for sub-time steps
        do k=1,ns

           call fluxx(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,tb1,tb2,tb3,tb4,&
                tb5,tb6,tb7,tb8,tb9,tba,tbb,fro,fru,frv,fre,xlx,yly,xmu,xba,eps,&
                eta,ftp,scp,xkt)

           call rkutta(rho,rou,rov,roe,fro,gro,fru,gru,frv,grv,&
                fre,gre,ftp,gtp,nx,ny,ns,dlt,coef,scp,k)

           call etatt(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,gma,chp)

        enddo
     endif
     !loop for the snapshots, to be save every imodulo
     if (mod(n,imodulo).eq.0) then
        !this is design for Gnuplot but feel free to implement your
        !own code if you want to use Matlab or Paraview
        write(nfichier, 990) n/imodulo
        x=0.
        do i=1,mx
           xx(i)=x
           x=x+dx
        enddo
        y=0.
        do j=1,my
           yy(j)=y
           y=y+dy
        enddo

        !computation of the vorticity
        call derix(vvv,nx,ny,tvv,xlx)
        call deriy(uuu,nx,ny,tuu,yly)
        do j=1,ny
        do i=1,nx
           wz(i,j)=tvv(i,j)-tuu(i,j)
        enddo
        enddo

        !using periodicity we copy the vorticity for the heat exchanger
        do ni=1,nf
        do nj=1,nf
           do j=1,ny
           do i=1,nx
              tf(i+(ni-1)*nx,j+(nj-1)*ny)=wz(i,j)
           enddo
           enddo
        enddo
        enddo
        !this file will be used by gnuplot for visualisations
        open(21,file=nfichier,form='formatted',status='unknown')
        do j=1,my
        do i=1,mx
           write(21,*)  xx(i),yy(j),tf(i,j)
        enddo
        write(21,*)
        enddo
        close (21)
     endif

     !Computation of average values
     call average(uuu,um,nx,ny)
     call average(vvv,vm,nx,ny)
     call average(scp,tm,nx,ny)

     !we write the average values for velocity and temperature
     write(*,*) n,um,vm,tm
!
  enddo
  !END OF THE TIME LOOP
!
end program navierstokes
!
!END OF THE MAIN PROGRAMME
!
!############################################
!
subroutine average(uuu,um,nx,ny)
!
!computation of the mean value of a 2D field
!############################################
!
  implicit none
!
  real(8),dimension(nx,ny) :: uuu
  real(8) :: um
  integer :: i,j,nx,ny,nxm

  um=0.
  nxm=nx
  do j=1,ny
  do i=1,nxm
     um = um + uuu(i,j)
  enddo
  enddo
  um=um/(real(nxm*ny))

  return
end subroutine average
!############################################

!############################################
!
subroutine derix(phi,nx,ny,dfi,xlx)
!
!First derivative in the x direction
!############################################

  implicit none

  real(8),dimension(nx,ny) :: phi,dfi
  real(8) :: dlx,xlx,udx
  integer :: i,j,nx,ny

  dlx=xlx/nx
  udx=1./(dlx+dlx)
  do j=1,ny
     dfi(1,j)=udx*(-3.*phi(1,j)+4.*phi(2,j)-phi(3,j))
     do i=2,nx-1
        dfi(i,j)=udx*(phi(i+1,j)-phi(i-1,j))
     enddo
     dfi(nx,j)=udx*(3.*phi(nx,j)-4.*phi(nx-1,j)+phi(nx-2,j))
  enddo

  return
end subroutine derix
!############################################

!############################################
!
subroutine deriy(phi,nx,ny,dfi,yly)
!
!First derivative in the y direction
!############################################

  implicit none

  real(8),dimension(nx,ny) ::  phi,dfi
  real(8) :: dly,yly,udy
  integer :: i,j,nx,ny

  dly=yly/ny
  udy=1./(dly+dly)
  do j=2,ny-1
     do i=1,nx
        dfi(i,j)=udy*(phi(i,j+1)-phi(i,j-1))
     enddo
  enddo
  do i=1,nx
     dfi(i,1)=udy*(phi(i,2)-phi(i,ny))
     dfi(i,ny)=udy*(phi(i,1)-phi(i,ny-1))
  enddo

  return
end subroutine deriy
!############################################

!############################################
!
subroutine derxx(phi,nx,ny,dfi,xlx)
!
!Second derivative in y direction
!############################################

  implicit none

  real(8),dimension(nx,ny) ::  phi,dfi
  real(8) :: dlx,xlx,udx
  integer :: i,j,nx,ny

  dlx=xlx/nx
  udx=1./(dlx*dlx)
  do j=1,ny
     dfi(1,j)=udx*(2.*phi(1,j) - 5.*phi(2,j) + 4.*phi(3,j) - phi(4,j))
     do i=2,nx-1
        dfi(i,j)=udx*(phi(i+1,j)-(phi(i,j)+phi(i,j))&
             +phi(i-1,j))
     enddo
     dfi(nx,j)=udx*(2.*phi(nx,j) - 5.*phi(nx-1,j) + 4.*phi(nx-2,j) - phi(nx-3,j))
  enddo

  return
end subroutine derxx
!############################################

!############################################
!
subroutine deryy(phi,nx,ny,dfi,yly)
!
!Second derivative in the y direction
!############################################

  implicit none

  real(8),dimension(nx,ny) ::  phi,dfi
  real(8) :: dly,yly,udy
  integer :: i,j,nx,ny

  dly=yly/ny
  udy=1./(dly*dly)
  do j=2,ny-1
     do i=1,nx
        dfi(i,j)=udy*(phi(i,j+1)-(phi(i,j)+phi(i,j))&
             +phi(i,j-1))
     enddo
  enddo
  do i=1,nx
     dfi(i,1)=udy*(phi(i,2)-(phi(i,1)+phi(i,1))+phi(i,ny))
     dfi(i,ny)=udy*(phi(i,1)-(phi(i,ny)+phi(i,ny))+phi(i,ny-1))
  enddo

  return
end subroutine deryy
!############################################

!############################################
!
subroutine derix4(phi,nx,ny,dfi,xlx)
!
!Fourth-order first derivative in the x direction
!############################################

  implicit none

  real(8),dimension(nx,ny) :: phi,dfi
  real(8) :: dlx,xlx,udx
  integer :: i,j,nx,ny



  return
end subroutine derix4
!############################################

!############################################
!
subroutine deriy4(phi,nx,ny,dfi,yly)
!
!Fourth-order first derivative in the y direction
!############################################

  implicit none

  real(8),dimension(nx,ny) ::  phi,dfi
  real(8) :: dly,yly,udy
  integer :: i,j,nx,ny



  return
end subroutine deriy4
!############################################

!############################################
!
subroutine derxx4(phi,nx,ny,dfi,xlx)
!
!Fourth-order second derivative in y direction
!############################################

  implicit none

  real(8),dimension(nx,ny) ::  phi,dfi
  real(8) :: dlx,xlx,udx
  integer :: i,j,nx,ny



  return
end subroutine derxx4
!############################################

!############################################
!
subroutine deryy4(phi,nx,ny,dfi,yly)
!
!Fourth-order second derivative in the y direction
!############################################

  implicit none

  real(8),dimension(nx,ny) ::  phi,dfi
  real(8) :: dly,yly,udy
  integer :: i,j,nx,ny



  return
end subroutine deryy4
!############################################


!#######################################################################
!
subroutine fluxx(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,tb1,tb2,tb3,&
     tb4,tb5,tb6,tb7,tb8,tb9,tba,tbb,fro,fru,frv,fre,xlx,yly,xmu,xba,&
     eps,eta,ftp,scp,xkt)
!
!#######################################################################

  implicit none

  real(8),dimension(nx,ny) :: uuu,vvv,rho,pre,tmp,rou,rov,roe,tb1,tb2,tb3,tb4,tb5,tb6,tb7
  real(8),dimension(nx,ny) :: tb8,tb9,tba,tbb,fro,fru,frv,fre,eps,ftp,scp
  real(8) :: utt,qtt,xmu,eta,dmu,xlx,yly,xba,xkt
  integer :: i,j,nx,ny

  call derix(rou,nx,ny,tb1,xlx)
  call deriy(rov,nx,ny,tb2,yly)
  do j=1,ny
     do i=1,nx
        fro(i,j)=-tb1(i,j)-tb2(i,j)
     enddo
  enddo

  do j=1,ny
     do i=1,nx
        tb1(i,j)=rou(i,j)*uuu(i,j)
        tb2(i,j)=rou(i,j)*vvv(i,j)
     enddo
  enddo

  call derix(pre,nx,ny,tb3,xlx)
  call derix(tb1,nx,ny,tb4,xlx)
  call deriy(tb2,nx,ny,tb5,yly)
  call derxx(uuu,nx,ny,tb6,xlx)
  call deryy(uuu,nx,ny,tb7,yly)
  call derix(vvv,nx,ny,tb8,xlx)
  call deriy(tb8,nx,ny,tb9,yly)
  utt=1./3
  qtt=4./3
  do j=1,ny
     do i=1,nx
        tba(i,j)=xmu*(qtt*tb6(i,j)+tb7(i,j)+utt*tb9(i,j))
        fru(i,j)=-tb3(i,j)-tb4(i,j)-tb5(i,j)+tba(i,j)&
             -(eps(i,j)/eta)*uuu(i,j)
     enddo
  enddo

  do j=1,ny
     do i=1,nx
        tb1(i,j)=rou(i,j)*vvv(i,j)
        tb2(i,j)=rov(i,j)*vvv(i,j)
     enddo
  enddo

  call deriy(pre,nx,ny,tb3,yly)
  call derix(tb1,nx,ny,tb4,xlx)
  call deriy(tb2,nx,ny,tb5,yly)
  call derxx(vvv,nx,ny,tb6,xlx)
  call deryy(vvv,nx,ny,tb7,yly)
  call derix(uuu,nx,ny,tb8,xlx)
  call deriy(tb8,nx,ny,tb9,yly)
  do j=1,ny
     do i=1,nx
        tbb(i,j)=xmu*(tb6(i,j)+qtt*tb7(i,j)+utt*tb9(i,j))
        frv(i,j)=-tb3(i,j)-tb4(i,j)-tb5(i,j)+tbb(i,j)&
             -(eps(i,j)/eta)*vvv(i,j)
     enddo
  enddo
!
!Equation for the tempature
!
  call derix(scp,nx,ny,tb1,xlx)
  call deriy(scp,nx,ny,tb2,yly)
  call derxx(scp,nx,ny,tb3,xlx)
  call deryy(scp,nx,ny,tb4,yly)

  do j=1,ny
     do i=1,nx
        ftp(i,j)=-uuu(i,j)*tb1(i,j)-vvv(i,j)*tb2(i,j)&
             + xkt*(tb3(i,j)+tb4(i,j))&
             - (eps(i,j)/eta)*scp(i,j)
     enddo
  enddo

  call derix(uuu,nx,ny,tb1,xlx)
  call deriy(vvv,nx,ny,tb2,yly)
  call deriy(uuu,nx,ny,tb3,yly)
  call derix(vvv,nx,ny,tb4,xlx)
  dmu=2./3*xmu
  do j=1,ny
     do i=1,nx
        fre(i,j)=xmu*(uuu(i,j)*tba(i,j)+vvv(i,j)*tbb(i,j))&
             +(xmu+xmu)*(tb1(i,j)*tb1(i,j)+tb2(i,j)*tb2(i,j))&
             -dmu*(tb1(i,j)+tb2(i,j))*(tb1(i,j)+tb2(i,j))&
             +xmu*(tb3(i,j)+tb4(i,j))*(tb3(i,j)+tb4(i,j))
     enddo
  enddo

  do j=1,ny
     do i=1,nx
        tb1(i,j)=roe(i,j)*uuu(i,j)
        tb2(i,j)=pre(i,j)*uuu(i,j)
        tb3(i,j)=roe(i,j)*vvv(i,j)
        tb4(i,j)=pre(i,j)*vvv(i,j)
     enddo
  enddo

  call derix(tb1,nx,ny,tb5,xlx)
  call derix(tb2,nx,ny,tb6,xlx)
  call deriy(tb3,nx,ny,tb7,yly)
  call deriy(tb4,nx,ny,tb8,yly)
  call derxx(tmp,nx,ny,tb9,xlx)
  call deryy(tmp,nx,ny,tba,yly)

  do j=1,ny
     do i=1,nx
        fre(i,j)=fre(i,j)-tb5(i,j)-tb6(i,j)-tb7(i,j)-tb8(i,j)+xba*(tb9(i,j)+tba(i,j))
     enddo
  enddo

  return
end subroutine fluxx
!#######################################################################

!###########################################################
!
subroutine rkutta(rho,rou,rov,roe,fro,gro,fru,gru,frv,grv,&
     fre,gre,ftp,gtp,nx,ny,ns,dlt,coef,scp,k)
!
!###########################################################

  implicit none
!
  real(8),dimension(nx,ny) :: rho,rou,rov,roe,fro,gro,fru,gru,frv
  real(8),dimension(nx,ny) :: grv,fre,gre,scp,ftp,gtp
  real(8),dimension(2,ns) :: coef
  real(8) :: dlt
  integer :: i,j,nx,ny,ns,k
!
!coefficient for RK sub-time steps
!!        coef(1,1)=XXX
!!        coef(1,2)=XXX
!!        coef(1,3)=XXX
!!        coef(2,1)=XXX
!!        coef(2,2)=XXX
!!        coef(2,3)=XXX

  do j=1,ny
     do i=1,nx
!!
     enddo
  enddo

  return
end subroutine rkutta
!###########################################################

!###########################################################
!
subroutine adams(rho,rou,rov,roe,fro,gro,fru,gru,frv,grv,&
     fre,gre,ftp,gtp,scp,nx,ny,dlt,dlx)
!
!###########################################################

  implicit none

  real(8),dimension(nx,ny) :: rho,rou,rov,roe,fro,gro,fru,gru,frv
  real(8),dimension(nx,ny) :: grv,fre,gre,ftp,gtp,scp
  real(8) :: dlt,ct1,ct2,dlx
  integer :: nx,ny,i,j

  ct1=1.5*dlt
  ct2=0.5*dlt
  do j=1,ny
     do i=1,nx-1
        rho(i,j)=rho(i,j)+ct1*fro(i,j)-ct2*gro(i,j)
        gro(i,j)=fro(i,j)
        rou(i,j)=rou(i,j)+ct1*fru(i,j)-ct2*gru(i,j)
        gru(i,j)=fru(i,j)
        rov(i,j)=rov(i,j)+ct1*frv(i,j)-ct2*grv(i,j)
        grv(i,j)=frv(i,j)
        roe(i,j)=roe(i,j)+ct1*fre(i,j)-ct2*gre(i,j)
        gre(i,j)=fre(i,j)
        scp(i,j)=scp(i,j)+ct1*ftp(i,j)-ct2*gtp(i,j)
        gtp(i,j)=ftp(i,j)
     enddo
     
     rho(nx,j) = rho(nx,j) - 0.2*dlt/dlx*(rho(nx,j)-rho(nx-1,j))
     rou(nx,j) = rou(nx,j) - 0.2*dlt/dlx*(rou(nx,j)-rou(nx-1,j))
     rov(nx,j) = rov(nx,j) - 0.2*dlt/dlx*(rov(nx,j)-rov(nx-1,j))
     roe(nx,j) = roe(nx,j) - 0.2*dlt/dlx*(roe(nx,j) - roe(nx-1,j))
     scp(nx,j) = scp(nx,j) - 0.2*dlt/dlx*(scp(nx,j) - scp(nx-1,j))
     
  enddo

  return
end subroutine adams
!###########################################################

!###########################################################
!
subroutine initl(uuu,vvv,rho,eee,pre,tmp,rou,rov,roe,nx,ny,&
     xlx,yly,xmu,xba,gma,chp,dlx,eta,eps,scp,xkt,uu0)
!
!###########################################################

  implicit none

  real(8),dimension(nx,ny) :: uuu,vvv,rho,eee,pre,tmp,rou,rov,roe,eps,scp
  real(8) :: xlx,yly,xmu,xba,gma,chp,roi,cci,d,tpi,chv,uu0
  real(8) :: epsi,pi,dlx,dly,ct3,ct4,ct5,ct6,y,x,eta,radius
  real(8) :: xkt
  integer :: nx,ny,i,j,ic,jc,imin,imax,jmin,jmax

  call param(xlx,yly,xmu,xba,gma,chp,roi,cci,d,tpi,chv,uu0)

  epsi=0.1
  dlx=xlx/nx
  dly=yly/ny
  ct3=log(2.)
  ct4=yly/2.
  ct5=xlx/2.
  ct6=(gma-1.)/gma
  y=-ct4
  x=0.
  eta=0.1
  eta=eta/2.
  radius=d/2.
  xkt=xba/(chp*roi)
  pi=acos(-1.)

!######for the square cylinder########################################
      ic=nint((xlx/2./dlx)+1) !X coordinate center of square
      jc=nint((yly/2./dly)+1) !Y coordinate center of square
!      imax=XXX
!      imin=XXX
!      jmax=XXX
!      jmin=XXX
!######################################################################



!##########CYLINDER DEFINITION#########################################
  do j=1,ny
     do i=1,nx
        if (((i*dlx-5.)**2+(j*dly-6.)**2).lt.radius**2) then
           eps(i,j)=1.
        else
           eps(i,j)=0.
        end if
     enddo
  enddo
!######################################################################
  do j=1,ny
     do i=1,nx
        uuu(i,j)=uu0
        vvv(i,j)=0.01*(sin(4.*pi*i*dlx/xlx)&
             +sin(7.*pi*i*dlx/xlx))*&
             exp(-(j*dly-yly/2.)**2)
        tmp(i,j)=tpi
        eee(i,j)=chv*tmp(i,j)+0.5*(uuu(i,j)*uuu(i,j)+vvv(i,j)*vvv(i,j))
        rho(i,j)=roi
        pre(i,j)=rho(i,j)*ct6*chp*tmp(i,j)
        rou(i,j)=rho(i,j)*uuu(i,j)
        rov(i,j)=rho(i,j)*vvv(i,j)
        roe(i,j)=rho(i,j)*eee(i,j)
        scp(i,j)=1.
        x=x+dlx
     enddo
     y=y+dly
  enddo

  return
end subroutine initl
!###########################################################

!################################################################
!
subroutine param(xlx,yly,xmu,xba,gma,chp,roi,cci,d,tpi,chv,uu0)
!
!################################################################

  implicit none

  real(8) :: ren,pdl,roi,cci,d,chp,gma,chv,xlx,yly,uu0,xmu,xba,tpi,mach

  ren=200.
  mach=0.2
  pdl=0.7
  roi=1.
  cci=1.
  d=1.
  chp=1.
  gma=1.4

  chv=chp/gma
  xlx=20.*d
  yly=12.*d
  uu0=mach*cci
  xmu=roi*uu0*d/ren
  xba=xmu*chp/pdl
  tpi=cci**2/(chp*(gma-1))

  return
end subroutine param
!################################################################

!################################################################
!
subroutine etatt(uuu,vvv,rho,pre,tmp,rou,rov,roe,nx,ny,gma,chp)
!
!################################################################

  implicit none

  real(8),dimension(nx,ny) ::  uuu,vvv,rho,pre,tmp,rou,rov,roe
  real(8) :: ct7,gma,ct8,chp
  integer :: i,j,nx,ny

  ct7=gma-1.
  ct8=gma/(gma-1.)

  do j=1,ny
     do i=1,nx
        uuu(i,j)=rou(i,j)/rho(i,j)
        vvv(i,j)=rov(i,j)/rho(i,j)
        pre(i,j)=ct7*(roe(i,j)-0.5*(rou(i,j)*uuu(i,j)+&
             rov(i,j)*vvv(i,j)))
        tmp(i,j)=ct8*pre(i,j)/(rho(i,j)*chp)
     enddo
  enddo

  return
end subroutine etatt
!################################################################
