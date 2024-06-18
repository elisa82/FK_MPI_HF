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
        real vp(nlay),vs(nlay),den(nlay),thick(nlay),qp(nlay),qs(nlay)
        real dist(ndis),t0(ndis),tmin(ndis),gf9x(nt*8*ndis)
        real t_cor,freq_ref,depth
        real dep_max,dep_min,dep_step,dist_max,dist_min,d_step,dt
c
        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
c
        if(myid .eq. 0) then
          open(1,file='Green.in',status='old')
          rewind(1) 
c          write(*,*)'Input name of velocity model'
          read(1,1) vm_name

c          write(*,*)'input Max., Min., and step of epicenter depth'
          read(1,*) dep_max, dep_min, dep_step

c          write(*,*)'input Max., Min., and step of epicenter distance'
          read(1,*) dist_max,dist_min,d_step

c          write(*,*)'input number of time-point, time-step, and ',
c     &              'seconds to save before first-arrival'
          read(1,*) numdt,dt,t_cor   

c          write(*,*)'Input the outputfile'
          read(1,1) output
1         format(a)
          close(1)
c
c input velocity structure model
          open(unit=30,file=vm_name,status='old')
          rewind(30)
          read(30,*) nlayer, freq_ref
          do k=1,nlayer
            read (30,*) vp(k),vs(k),den(k),thick(k),qp(k),qs(k)
          enddo 
          close(30) 

c
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
          stop
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
              call cmodel(nlayer,freq_ref,depth,vp,vs,den,thick,qp,qs)
              call trav(dist,nx,tmin)
              do k=1,nx
                t0(k)=tmin(k)-t_cor
                if(t_cor .lt. 0.0) t0(k)=0.0
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
                write(11,rec=ll) iz,k,dist(k),t0(k),depth,dt,npt,
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
c          write(19,*) 'Maximum, Minimum, and Step of Epicenter Depth'
c          write(19,*) dep_max, dep_min, dep_step
c          write(19,*)'Maximum, Minimum, and Step of Epicenter Distance'
c          write(19,*) dist_max,dist_min,d_step
          write(19,*)"Min, Step1 and N1, and Step2 and N2 of Depth" !JORGE new
          write(19,*) dep_min, dep_step, int(nz-1), 20.00, int(0) !JORGE new
          write(19,*)"Min, Step1 and N1, and Step2 and N2 of Epi Dist" !JORGE new
          write(19,*) dist_min, d_step, int(nx-1), 20.00, int(0) !JORGE new
c JORGE changed the whole write of the Green Bank info file
          write(19,*)'The lnpt, dt, block_gg, t_cor'
          write(19,*) lnpt,dt,block_gg,t_cor
          write(19,*)'The name of file to store Green Bank'
          write(19,'(a)') output
          write(19,*)'The number of GFun in Distance and Depth'
          write(19,*) nx,nz
          write(19,*) 'Velocity Structure model'
          write(19,*) nlayer,2
          do k=1,nlayer
c JORGE changed format of the Green_Bank.inf file to match syn1D specs
c            write(19,193) thick(k),vp(k),vs(k),den(k),qp(k),qs(k)
            write(19,*) thick(k),vp(k),vs(k),den(k),qp(k),qs(k)
          enddo
          close(19)
c193       format(6f12.4)
        endif

        call MPI_FINALIZE(ierr)
        stop
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

