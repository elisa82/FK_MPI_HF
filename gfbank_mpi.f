c       Make the Green Function bank
c       Input vm_name
c       Input dep_max dep_min dep_step
c       Input dist_max dist_min d_step
c       Input numdt dt,t_cor
c       Input Green_name
c       output Green
        implicit none
        include 'mpif.h'
        include 'model.h'
        character*60 vm_name, output,GBFinf
        integer block_gg,numdt,nlayer,lnpt,npt,itc,icom,getlen
        integer ierr,myid,numprocs,k,nx,nz,iz,ll,nzdiv,kdz,nlch
        integer nmsg,ip,ncount,istatus(MPI_STATUS_SIZE)
        real vp,vs,den,thick,qp,qs
        real qo, alpha  !JORGE added new vlues to describe anelasticity
        real vpt,vst,dent,thickt,qpt,qst  !JORGE added new variables to read target vel struct
        integer nlayert   !JORGE added new variables to read target vel struct
        real freq_reft !JORGE added new variables to read target vel struct
        real dist(ndis),t0(ndis),t00(ndis),t_hf(ndis),t_lf(ndis)
        real gf9x(nt*8*ndis)
        real t_s_hf(ndis),t_s_lf(ndis)   !JORGE added S arrival times
        real t_cor,freq_ref,depth
        real dep_max,dep_min,dep_step,dist_max,dist_min,d_step,dt
        real dummy(nt)    !JORGE added this new variale to manipulate the green bank
        common/vmod/vp(nlay),vs(nlay),den(nlay),thick(nlay),
     &       qp(nlay),qs(nlay),nlayer
        common/vmod_target/vpt(nlay),vst(nlay),dent(nlay),thickt(nlay),
     &       qpt(nlay),qst(nlay),nlayert
c
        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )


        if(myid .eq. 0) then
          open(1,file='Green.in',status='old')
          rewind(1) 
c          write(*,*)'Input name of velocity model'
          read(1,1) vm_name
C          write(*,*) vm_name
c          write(*,*)'input Max., Min., and step of epicenter depth'
          read(1,*) dep_max, dep_min, dep_step

c          write(*,*)'input Max., Min., and step of epicenter distance'
          read(1,*) dist_max,dist_min,d_step

c          write(*,*)'input number of time-point, time-step, and ',
c     &              'seconds to save before first-arrival'
          read(1,*) numdt,dt,t_cor   

c          write(*,*)'Input the outputfile'
          read(1,1) output
c JORGE added a readfile parameter of qo and alpha
          read(1,*) qo,alpha
c JORGE added a readfile parameter of qo and alpha
1         format(a)
          close(1)
c
c input velocity structure model
          open(unit=30,file=vm_name,status='old')
          rewind(30)
          read(30,*) nlayer, freq_ref
          do k=1,nlayer
            read (30,*) vp(k),vs(k),den(k),thick(k),qp(k),qs(k)
C            write(*,*) vp(k),vs(k),den(k),thick(k),qp(k),qs(k)
          enddo 
          close(30) 
CCCC JORGE ******
c input velocity structure model
          open(unit=31,file='target.vel',status='old')
          rewind(31)
          read(31,*) nlayert, freq_reft
          do k=1,nlayert
            read (31,*) vpt(k),vst(k),dent(k),thickt(k),qpt(k),qst(k)
C            write(*,*) nlayert,vpt(k),vst(k),dent(k),thickt(k),qpt(k),
C     &                 qst(k)
          enddo
          close(31)
CCCC JORGE ******

          nlch=getlen(output)
          GBFinf=output(1:nlch)//'.inf'
          lnpt=0
          npt=1
          do while (npt .lt. numdt) 
            lnpt=lnpt+1
            npt=2**lnpt
          enddo
          nx=max0(int((dist_max-dist_min)/d_step+0.01)+1, 2)
          nz=max0(int((dep_max-dep_min)/dep_step+0.01)+1, 2)
          nzdiv=(nz-1)/numprocs+1
          nz=nzdiv*numprocs
          d_step=(dist_max-dist_min)/(nx-1)
          dep_step=(dep_max-dep_min)/(nz-1)
          block_gg=4*(npt*8+7)+12
          open(11,file=output,status='replace',access='direct',
     &            recl=block_gg)
c Check the dimension
          nmsg=0
          if(npt.gt.nt) then
            nmsg=nmsg+1
            write(*,*) 'Please increase <nt> in model.h to ',npt
          endif
          if(dep_max.lt.dep_min)then
            nmsg=nmsg+1
            write(*,*) 'The depth region is wrong: dep_max < dep_min'
          endif
          if(dist_max.lt.dist_min)then
            nmsg=nmsg+1
            write(*,*) 'distance region is wrong: dist_max < dist_min'
          endif
          if(nx.gt.ndis) then
            nmsg=nmsg+1
            write(*,*) 'Please increase <ndis> in model.h to ',nx
          endif
        endif
c
        call MPI_BCAST(nmsg, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(nmsg .gt. 0) then
          call MPI_FINALIZE(ierr)
c#elisa          stop
        endif
        call MPI_BCAST(lnpt,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dt,    1,MPI_REAL,   0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(t_cor, 1,MPI_REAL,   0,MPI_COMM_WORLD,ierr)

        call MPI_BCAST(nz,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nx,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dep_max, 1,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dep_min, 1,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dep_step,1,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dist_max,1,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dist_min,1,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(d_step,  1,MPI_REAL, 0,MPI_COMM_WORLD,ierr)

        call MPI_BCAST(nlayer,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(freq_ref,  1,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vp,   nlayer,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vs,   nlayer,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(den,  nlayer,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(thick,nlayer,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qp,   nlayer,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qs,   nlayer,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
C****** JORGE added this portion 
        call MPI_BCAST(t_lf, 1,MPI_REAL,   0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(t_hf, 1,MPI_REAL,   0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(t_s_lf, 1,MPI_REAL,   0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(t_s_hf, 1,MPI_REAL,   0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(t0, 1,MPI_REAL,   0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(t00, 1,MPI_REAL,   0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nlayert,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(freq_reft,  1,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vpt,   nlayert,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vst,   nlayert,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dent,  nlayert,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(thickt,nlayert,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qpt,   nlayert,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qst,   nlayert,MPI_REAL, 0,MPI_COMM_WORLD,ierr)
C****** JORGE
        npt=2**lnpt
        ncount=npt*8*nx
        nzdiv=nz/numprocs
        do kdz=1,nzdiv
          do ip=1,numprocs
            iz=(kdz-1)*numprocs+ip
            if(myid+1 .eq. ip) then
              depth=dep_min + (iz-1)*dep_step
              do k=1,nx
                dist(k)=dist_min+(k-1)*d_step
              enddo
              if(depth.eq.0.0) depth=dep_step !JORGE to include 0 depth 
              call cmodel(nlayert,freq_reft,depth,vpt,vst,dent,thickt,
     &                    qpt,qst)
              call trav(dist,nx,t_lf)
              call trav_s(dist,nx,t_s_lf)
              call cmodel(nlayer,freq_ref,depth,vp,vs,den,thick,qp,qs)
              call trav(dist,nx,t_hf)
              call trav_s(dist,nx,t_s_hf)
C              write(*,*)'nlayer, nlayert',nlayer,nlayert
C              write(*,*)'t_lf, t_hf=',t_lf(1),t_hf(1)
              do k=1,nx
                t0(k)=t_lf(k)-t_cor-(t_s_lf(k)-t_s_hf(k)) !JORGE made change to align S arrivals
                t00(k)=t_lf(k)-t_cor
                if(t_cor.lt.0.0) t0(k)=0.0
C                write(*,*)'t_min and depth are:',tmin(k),depth
              enddo
              call sub_bs_dc(nx,dist,dist_max,lnpt,dt,t0,gf9x)
            endif
          enddo
          call MPI_Barrier(MPI_COMM_WORLD, ierr)
c Send_Recv
          do ip=1,numprocs
            if(ip.gt.1 .and.  myid.eq.(ip-1)) then
              call MPI_SEND(gf9x, ncount, MPI_REAL, 
     &                      0, 2*kdz,   MPI_COMM_WORLD, ierr)
              call MPI_SEND(t0, nx, MPI_REAL, 
     &                      0, 2*kdz+1, MPI_COMM_WORLD, ierr)
            endif
            if(ip.gt.1 .and. myid.eq.0) then
              call MPI_RECV(gf9x, ncount, MPI_REAL,
     &                      ip-1, 2*kdz,  MPI_COMM_WORLD, istatus,ierr)
              call MPI_RECV(t0, nx, MPI_REAL,
     &                      ip-1, 2*kdz+1,MPI_COMM_WORLD, istatus,ierr)
            endif
            if(myid.eq.0) then
              iz=(kdz-1)*numprocs+ip
              depth=dep_min + (iz-1)*dep_step
              do k=1,nx
                dist(k)=dist_min+(k-1)*d_step
                ll=(iz-1)*nx+k
CCCC JORGE ******
                do icom=1,8
                   do itc=1,npt
                      dummy(itc)=gf9x((k*8+icom-9)*npt+itc)
                   enddo
                   if(depth.eq.0.0) then
                   call site(npt,dt,dummy,0.0,dep_step,dist(k),qo,alpha)
                   else
                      call site(npt,dt,dummy,0.0,depth,dist(k),qo,alpha)
                   endif
                   do itc=1,npt
                      gf9x((k*8+icom-9)*npt+itc)=dummy(itc)
                   enddo
                enddo
CCCC JORGE ******
                write(11,rec=ll) iz,k,dist(k),t00(k),depth,dt,npt,
     &              ((gf9x((k*8+icom-9)*npt+itc), itc=1,npt),icom=1,8)
              enddo
            endif
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
          enddo
        enddo
        if(myid .eq. 0) then
          close(11)
          open(19,file=GBFinf,status='replace')
          rewind(19)
          write(19,*) 'Name of velocity model'
          write(19,'(a)') vm_name
          write(19,*) 'Maximum, Minimum, and Step of Epicenter Depth'
          write(19,*) dep_max, dep_min, dep_step
          write(19,*)'Maximum, Minimum, and Step of Epicenter Distance'
          write(19,*) dist_max,dist_min,d_step
          write(19,*)'The lnpt, dt, block_gg, t_cor'
          write(19,*) lnpt,dt,block_gg,t_cor
          write(19,*)'The name of file to store Green Bank'
          write(19,'(a)') output
          write(19,*)'The number of GFun in Distance and Depth'
          write(19,*) nx,nz
          write(19,*) 'Velocity Structure model'
          write(19,*) nlayer,2
          do k=1,nlayer
             write(19,193) thick(k),vp(k),vs(k),den(k),qp(k),qs(k)
          enddo
          close(19)
193       format(6f12.4)
        endif

        call MPI_FINALIZE(ierr)
c#elisa        stop
        end        
c
c=====================================================================
c
        function getlen(string)
        implicit NONE
        character *(*) string
        character blank*1
        parameter (blank=' ')
        integer i,n1,n2,nc,getlen

        do i=1,len(string)
           if(string(i:i) .ne. blank) goto 11
        enddo
11        n1=i
        do i=len(string),1,-1
           if(string(i:i) .ne. blank) goto 12
        enddo
12      n2=i
        nc=n2-n1+1

        if(nc.gt.0) then
          string(1:nc)=string(n1:n2)
        else
          nc=0
        endif
        getlen=nc
        return
        end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       real function q_f(qo,alpha,freq)
       real :: qo, freq, alpha
       q_f = qo*(freq+0.001)**alpha
       return
       end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine site(n,dt,y,kap,depth,dist,qo,alpha)
        include 'model.h'
        integer :: n,nlayert
        real :: f,df,dt,y(n),pi=4.0*ATAN(1.0),kap,an,depth,dist,R
        real :: alpha,qo
        common/vmod_target/vpt(nlay),vst(nlay),dent(nlay),thickt(nlay),
     &       qpt(nlay),qst(nlay),nlayert
        real :: vp,vs,den,thick,qp,qs,vpt,vst,dent,thickt,qpt,qst
        complex :: c(n)
        R=(dist**2+depth**2)**0.5
        do i=1,n
          c(i)=CMPLX(y(i),0.0)
        end do
        call fork(n,c,-1.)
        df=1/(n*dt)
        do i=1,n/2
          f=(i-1)*df
          call get_sitefacs(depth,f,an,nlayert,thickt,vst,dent)
          c(i)=c(i)*exp(-pi*f*R/(q_f(qo,alpha,f)*3.5))*exp(-pi*f*kap)*an
        end do
        do i=n/2+2,n
          c(i)=CONJG(c(n+2-i))
        end do
        call fork(n,c,1.)
        do i=1,n
          y(i)=real(c(i))
        end do
        end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_sitefacs(depth,f,a,j0,th,vs,dn)
      implicit none
      real f,a
      real zdep,pz,bz,tt,ttp,stt,vdsrc,depth
      integer i,iz,j0
      real th(j0),vs(j0),dn(j0),z(j0)
      if(f.eq.0.0) then
         a=1.0
      else
      do i=1,j0
         z(i) = th(i)
         if(i.gt.1) then
               z(i) = z(i) + z(i-1)
         endif
      enddo
      iz=0
      do while(z(iz).lt.depth)
         iz=iz+1
      enddo

      vdsrc = vs(iz)*dn(iz)
         stt = 0.25/f
         i = 2
         zdep = 0.0
         pz = 0.0
         tt = 0.0
         ttp = th(2)/vs(2)
6145     if(ttp.ge.stt.or.i.eq.iz) go to 6146
         zdep = zdep + th(i)
         pz = pz + dn(i)*th(i)/vs(i)
         tt = ttp
         i = i + 1
         ttp = th(i)/vs(i) + tt
         go to 6145
6146     continue
         bz = (zdep + (stt - tt)*vs(i))/stt
         pz = (pz + dn(i)*(stt - tt))/stt
         a = sqrt(vdsrc/(bz*pz))
      endif
      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c JORGE calculate Fourier Spectra (from Boore)
! ----------------------------- BEGIN FORK --------------------------
      SUBROUTINE FORK(LX,CX,SIGNI)
C FAST FOURIER                                  2/15/69
C                          LX
C    CX(K) = SQRT(1.0/LX)* SUM (CX(J)*EXP(2*PI*SIGNI*I*(J-1)*(K-1)/LX))
C                          J=1                        FOR K=1,2,...,LX
C
C  THE SCALING BETWEEN FFT AND EQUIVALENT CONTINUUM OUTPUTS
C  IS AS FOLLOWS.
C
C
C     GOING FROM TIME TO FREQUENCY:
C             F(W)=DT*SQRT(LX)*CX(K)
C
C                  WHERE W(K)=2.0*PI*(K-1)*DF

*                  and    DF = 1/(LX*DT)
C
C
C     GOING FROM FREQUENCY TO TIME, WHERE THE FREQUENCY
C     SPECTRUM IS GIVEN BY THE DIGITIZED CONTINUUM SPECTRUM:
C
C             F(T)=DF*SQRT(LX)*CX(K)
*
C                  WHERE T(K)=(K-1)*DT
C
C
C  THE RESULT OF THE SEQUENCE...TIME TO FREQUENCY,POSSIBLE MODIFICATIONS
C  OF THE SPECTRUM (FOR FILTERING,ETC.), BACK TO TIME...
C  REQUIRES NO SCALING.
C
C
C  THIS VERSION HAS A SLIGHT MODIFICATION TO SAVE SOME TIME...
C  IT TAKES THE FACTOR 3.1415926*SIGNI/L OUTSIDE A DO LOOP (D.BOORE 12/8
C  FOLLOWING A SUGGESTION BY HENRY SWANGER).
C

* Some brief notes on usage:

* "signi" is a real variable and should be called either with the value
* "+1.0"
* of "-1.0".  The particular value used depends on the conventions being
* used
* in the application (e.g., see Aki and Richards, 1980, Box 5.2, pp.
* 129--130).

* Time to frequency:
* In calling routine,
*
*       do i = 1, lx
*         cx(i) = CMPLX(y(i), 0.0)
*       end do
*  where y(i) is the time series and lx is a power of 2
*
*  After calling Fork with the complex array specified above, the
*  following
* symmetries exist:
*
*        cx(1)        = dc value (f = 0 * df, where df = 1.0/(lx*dt))
*        cx(lx/2 + 1) = value at Nyquist (f = (lx/2+1-1)*df =
*        1.0/(2*dt))
*        cx(lx)       = CONJG(cx(2))
*        cx(lx-1)     = CONJG(cx(3))
*         |           =      |
*        cx(lx-i+2)   = CONJG(cx(i))
*         |           =      |
*        cx(lx/2+2)   = CONJG(cx(lx/2))
*
* where "CONJG" is the Fortran complex conjugate intrinsic function
*
* This symmetry MUST be preserved if modifications are made in the
* frequency
* domain and another call to Fork (with a different sign for signi) is
* used
* to go back to the time domain.  If the symmetry is not preserved, then
* the
* time domain array will have nonzero imaginary components.  There is
* one case
* where advantage can be taken of this, and that is to find the Hilbert
* transform and the window of a time series with only two calls to Fork
* (there
* is a short note in BSSA {GET REFERENCE} discussing this trick, which
* amounts
* to zeroing out the last half of the array and multiplying all but the
* dc and
* Nyquist values by 2.0; in the time domain, REAL(cx(i)) and
* AIMAG(cx(i))
* contain the filtered (if a filter was applied) and Hilbert transform
* of the
* filtered time series, respectively, while CABS(cx(i)) and
* ATAN2(AIMAG(cx(i)),
* REAL(cx(i))) are the window and instantaneous phase of the filtered
* time
* series, respectively.

* Some references:

* Farnbach, J.S. (1975). The complex envelope in seismic signal
* analysis,
* BSSA 65, 951--962.
* He states that the factor of 2 is applied for i = 2...npw2/2 (his
* indices
* start at 0, Ive added 1), which is different than the next reference:

* Mitra, S.K. (2001). Digital Signal Processing, McGraw-Hill, New York.
* He gives an algorithm on p. 794 (eq. 11.81), in which the factor of 2
* is
* applied from 0 frequency to just less than Nyquist.

*
* The easiest way to ensure the proper symmetry is to zero out the
* last half of the array (as discussed above), but the following is what
* I usually use:
* modify (filter) only half
* of the cx array:
*
*       do i = 1, lx/2
*         cx(i) = filter(i)*cx(i)
*       end do
*
* where "filter(i)" is a possibly complex filter function (and recall
* that
* the frequency corresponding to i is f = float(i-1)*df).  After this,
* fill out
* the last half of the array using
*
*       do i = lx/2+2, lx
*         cx(i) = CONJG(cx(lx+2-j))
*       end do
*
* Note that nothing is done with the Nyquist value.  I assume (but am
* not sure!)
* that this value should be 0.0
*
* Dates: xx/xx/xx - Written by Norm Brenner(?), Jon Claerbout(?)
*        12/21/00 - Replaced hardwired value of pi with pi evaluated
*        here,
*                     and added comments regarding usage.  Also deleted
*                     dimension specification of cx(lx) and replace it
*                     with
*                     cx(*) in the type specification statement.  I also
*                     cleaned up the formatting of the subroutine.
*        08/28/01 - Added comment about variable "signi" being real, and
*                   added "float" in equations for "sc" and "temp",
*                   although
*                   not strictly required.
*        06/19/02 - Added some comments about computing envelopes and
*                   instantaneous frequencies

      complex cx(*),carg,cexp,cw,ctemp

      pi = 4.0*atan(1.0)

      j=1
      sc=sqrt(1./float(lx))

      do i=1,lx
        if(i.gt.j) go to 2
        ctemp=cx(j)*sc
        cx(j)=cx(i)*sc
        cx(i)=ctemp
2       m=lx/2
3       if(j.le.m) go to 5
        j=j-m
        m=m/2
        if(m.ge.1) go to 3
5       j=j+m
      end do

      l=1
6     istep=2*l
      temp= pi * signi/float(l)

      do m=1,l
        carg=(0.,1.)*temp*(m-1)
        cw=cexp(carg)
        do i=m,lx,istep
          ctemp=cw*cx(i+l)
          cx(i+l)=cx(i)-ctemp
          cx(i)=cx(i)+ctemp
        end do
      end do

      l=istep
      if(l.lt.lx) go to 6

      return
      end
! ----------------------------- END FORK --------------------------
