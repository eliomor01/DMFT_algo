c#########################################################
c 2nd Order Perturbation Theory (a la Yamada Yosida)
c for the oo-d Hubard Model on a Bethe lattice
c and on the imaginary axis.
c Produces the local Green function and the Selfenergy.
c Also calculates the double occupation.
c by M.J. Rozenberg (6/2000)
c#########################################################

 	parameter(L=4096*16)
c L=number of frequency points
 	implicit real*8(a-h,o-z)
        logical phsym
 	double complex xi,ep,om,om1,om2,d2,one,root
 	double complex fg0(2*L),fg0b(2*L),fg0w(2*L)      
 	double complex fg0t(2*L),fg(2*L)    
 	double complex sefb(2*L),seff(2*L),self(2*L)
        double precision data(4*L)

c******* output data **********************************************
        open(60,file="Gw",form="formatted",status="unknown")
        open(62,file="G0w",form="formatted",status="unknown")
        open(64,file="Sigma",form="formatted",status="unknown")
        open(66,file="n_vs_mu",form="formatted",status="unknown")
        open(67,file="double_occ",status="unknown",access="append")
        open(11,file="test",form="formatted",status="unknown")


c******* input data **************************************
        open(70,file="input",form="formatted",status="old")

        prec=1.d-16

        read(70,*)D,      xmu0
        read(70,*)dmu,    nmu
        read(70,*)itu,    nl
        read(70,*)u0,     du 
        read(70,*)t0,     dt
        read(70,*)nloop,  imet
        read(70,*)xmix, nl0

c  D    = half-bandwidth of semi-circular DOS (Bethe lattice)
c  xmu0 = initial chemical potential
c  dmu  = step for the chemical potential variation
c  nmu  = number of chemical potential points
c  itu  = flag that selects Temperature or U variation
c         itu=0 for U   or   itu=1 for Temperature
c     itu=2 plots the output as a function of the half-bandwidth D
c  nl   = number of U or Temperaure points
c  u0   = initial interaction strength U
c  du   = step for the U variation
c  t0   = initial Temperature
c  dt   = step for the Temperature variation
c  nloop= number of self-concistency iteration loops
c  imet = choose the initial seed (1=metallic, 0=insulating)
c  nl0  = number of final loops to be printed-out
          
          phsym=.false.
          if (phsym) print*, "particle-hole symmetry"

c*******define some constants*******

        ep=1.d-9*(1.d0,0.d0)
        d2=d*(1.d0,0.d0)
        d2=d2**2
        pi=3.14159265358979323d0
	xi=(0.d0,1.d0)
        one=(1.d0,0.d0)
	u=u0
	t=t0

c*******************************************************************
c*******************************************************************

c******mu loop starts here***********************************

	do 333 imu=1,nmu
	xmu=xmu0+dmu*dfloat(imu-1)


c*******temperature/U loop starts here***********************

	do 222 il=1,nl
	if(itu.eq.1)then
	t=t0+dt*(dfloat(il-1))
	elseif(itu.eq.0)then
        u=u0+du*dfloat(il-1)
	endif

	dtau=1.d0/t/dfloat(L)



c** the self-consistency loop starts here ********************

        do 103 iloop=1,nloop

	if(iloop.eq.nloop/2)then
	print*,'I am half way thru...'
	endif

	do 301 i=1,2*L
	   fg0(i)=(0.d0,0.d0)	
	   fg(i)=(0.d0,0.d0)
 301	continue


	if(iloop*il.eq.1)then
c***on first loop compute a seed for Go**
c  (imet=1 metallic, imet=0 insulating)

        do 101 i=1,L
           om=xi*(2.d0*dfloat(i)-1.d0-dfloat(L))*pi*t+xmu

c** insulating seed for Go**
           fg0(2*i)=one/om

        if(imet.eq.1)then
c** metallic seed for Go**
           sig=1.d0
           if(dimag(om).lt.0.d0)sig=-1.d0
	   if(dimag(om).eq.0.d0)then
	      fg0(2*i)=(0.d0,0.d0)
	   else
              fg0(2*i)=one/(om+xi*d*sig)
	   endif
        endif
 101    continue

	else

c***on following loops compute a new Go from the old Sigma**
c***also compute a new G from the old Sigma**
        if (phsym) then
           do i=1,L
              self(2*i)=dcmplx(0.d0,dimag(self(2*i)))
           enddo
        endif

        test=0.d0
        do 102 i=1,L
           om=xi*(2.d0*dfloat(i)-1.d0-dfloat(L))*pi*t+xmu
           om1=om+self(2*i)
           om2=om-self(2*i)
 	   root=cdsqrt((om2+ep)**2-d2)
c**choose the branch-cut**
	   sig=1.d0
 	   if(dimag(om)*dimag(root).lt.0.d0)sig=-1.d0
           fg0(2*i)=((1-xmix)*fg0w(2*i)+xmix*2.d0*one/(om1+sig*root))!mixing
           if (phsym) then
              fg0(2*i)=dcmplx(0.d0,dimag(fg0(2*i)))
           endif
	   fg(2*i)=2.d0*one/(om2+sig*root)

cccccccccc calculates the convergence test
           test=test+abs((fg0(2*i)-fg0w(2*i))/dimag(om))
           xnorm=xnorm+abs(fg0(2*i))
cccccccccccccccccccccccccccccccccccccccccc           
           fg0w(2*i)=fg0(2*i)	
 102    continue
        test=test/xnorm
        write(11,*)iloop,test
        call flush(11)
	endif

c** fg0(i) is Go(iw) **
c** fg(i)  is G(iw)  **

       do i=1,2*L
        data(2*i-1)=dreal(fg0(i))
        data(2*i)=dimag(fg0(i))
        enddo
        call four1(data,2*L,1)
        do i=1,2*L
        fg0t(i)=dcmplx(data(2*i-1),data(2*i))
        enddo

	ex=-1.d0
        do 82 i=1,2*L
	ex=-ex
           fg0t(i)=t*fg0t(i)*ex
 82      continue

	do 83 i=1,L
           fg0b(i)=fg0t(i+L)
 83      continue
        
	do 84 i=L+1,2*L
           fg0b(i)=fg0t(i-L)
 84      continue
        
c** fg0b(i) is Go(tau)  [-beta,beta] *****



c***calculate the Sigma in 2OPT****************

	do 520 i=1,L
         sefb(i+L)=-u**2*fg0b(i+L)**2*fg0b(L+2-i)
520	continue

	do 525 i=1,L-1
         sefb(L+1-i)=-u**2*fg0b(L+1-i)**2*fg0b(i+L+1)
525	continue
c** fix the end-point**
	sefb(1)=-sefb(L+1)

c** sefb(i) is Sigma(tau)  [-beta,beta] *******


       do i=1,2*L
        data(2*i-1)=dreal(sefb(i))
        data(2*i)=dimag(sefb(i))
        enddo
        call four1(data,2*L,-1)
        do i=1,2*L
        seff(i)=dcmplx(data(2*i-1),data(2*i))
        enddo

	ex=-1.
	do 530 i=1,2*L
	ex=-ex
           seff(i)=.5d0*dtau*ex*seff(i)
530      continue

	do 540 i=1,L
           self(i)=seff(i+L)
540      continue
        
	do 550 i=L+1,2*L
           self(i)=seff(i-L)
550     continue

c***self(i) is Sigma(iw)******


        if((iloop.ge.nloop-nl0+1).or.
     &       ((iloop.gt.1).and.(test.lt.prec))) then 
c**on final loop compute <n>, <d>, and print-out***

c******get particle number <n>*********

	xn=0.d0
	xni=0.d0
	do i=1,L
	xn=xn+dreal(fg(2*i))*t
	enddo

c******get the double ocupation <d>****

	sum=0.d0
	tum=0.d0
	tuma=0.d0
	do 111 i=1,L
           oma=(2.d0*dfloat(i)-1.d0-dfloat(L))*pi*t
	sg=1.d0
	if(oma.lt.0.d0)sg=-1.d0
	omb=oma*2.d0
 	tum=tum+.5d0*dimag(fg(2*i))*dimag(self(2*i))
     $ -(dimag(fg(2*i))+2.d0/(oma+sg*dsqrt(oma**2+d**2)))*oma
	tuma=tuma+dimag(fg(2*i))*dimag(self(2*i))
     $ -(dimag(fg(2*i))+2.d0/(oma+sg*dsqrt(oma**2+d**2)))*oma

111	continue
 
	free=0.d0
	e=-d
	de=d/1000.d0
	do 666 i=1,2000
 	 free=free+de*e*dsqrt(1.d0-(e/d)**2)/(dexp(e/t)+1.d0)
	 e=e+de
666	continue	
	free=free*2.d0/(pi*d)

	tum=t*tum+free
	tuma=t*tuma+free

	doble=(tum-tuma)*2.d0/u+.25d0

c************print output*******************************

 	   x=pi/(L*dtau)
           nnn=0
 	do 106 i=L-2000,L+2000,2
           si=dimag(self(i))
 	   g1=dimag(fg(i))
 	   g2=dimag(fg0w(i))
           write(60,*)(x*(i-L-1)),dreal(fg(i)),(g1)
           write(62,*)(x*(i-L-1)),dreal(fg0w(i)),(g2)
           write(64,*)(x*(i-L-1)),dreal(self(i)),(si)
 106      continue
 	write(60,*)'   '
 	write(62,*)'   '
 	write(64,*)'   '

c print-out the ocupation
	write(66,*)(xmu),(xn)
c print-out the double ocupation
	if(itu.eq.1)then
 	write(67,*)t,(doble)
	elseif(itu.eq.0) then
 	write(67,*)u,(doble)
	elseif(itu.eq.2) then
 	write(67,*)d,(doble)
	endif

c fort.66 has <n>
c fort.67 has <d>

 	endif
c********************************************************


c*******close iteration loop*********************************
        if ((iloop.gt.1).and.(test.lt.prec)) goto 222
 103    continue


c*******close temperature/U loop*****************************

222	continue

c*******close mu loop****************************************
333	continue

	stop
	end


c########################################################
c########################################################
c from numerical recipes
c########################################################
      SUBROUTINE four1(data,nn,isign)
        implicit real*8(a-h,o-z)
        double precision data(2*nn)
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END

c########################################################
c########################################################

