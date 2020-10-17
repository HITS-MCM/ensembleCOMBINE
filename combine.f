
c   Copyright (C) 2009  Javier Klett, Rubén Gil-Redondo, Federico Gago and
c   Antonio Morreale
c   Alvaro Cortes currently maintains and improve the code. 2012,2013
c   <alvaro.cortes@uah.es>

c   Based on the original work by Ángel R. Ortiz

c   The following terms apply to all files and documents associated with the 
c   software unless explicitly disclaimed in individual files. Do not 
c   redistribute the program. Interested users should contact directly to the 
c   authors. This software is for scientific non-profit and non-commercial use
c   only. Any other use of this software for other purposes, alone or 
c   integrated into other software, requires the prior consent of the authors.

c   IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR
c   DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
c   OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF,
c   EVEN IF THE AUTHORS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

c   THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, 
c   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, 
c   FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS 
c   PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO 
c   OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR 
c   MODIFICATIONS.

      PROGRAM       COMBINE

      IMPLICIT      DOUBLE PRECISION (A-H,O-Z)

      parameter    (maxatdrg = 200)
      parameter    (maxcompl = 100)
      PARAMETER    (MAXATOM = 8000)
      PARAMETER    (MAXRES  =  700)

      parameter    (maxrow= 100)
      parameter    (maxcol=1500)
      parameter    (maxadd=  10)
      parameter    (maxrank =20)

      dimension     xyz(3,maxatom),add(maxcompl,maxadd)
      dimension     xraw(maxrow,maxcol),y(maxrow),coef(maxcol)
      dimension     coefs(maxcol)
      character*4   drugname(maxcompl), name(maxrow)
      integer       dres, datm(maxatdrg),k
      real*8        b(maxcol,maxcol)

      integer contour,ptr, idel, nsc,imat, nlv
      integer ielec,nad,ncomp,ntest
      integer varselect(maxcol),nvtold
      integer randtest,ncross
      integer ftrim,btrim
      integer multres,m,watMol

      real*8 cutinp, cutptr

      CHARACTER*4   IGRAPH, LABRES
      character*40  topfile,rstfile,inpfile,outfile
      character*40  comp(maxcompl)
      character*20  vardesc(maxcol)
      character*6  chartemp

      COMMON /INFOV/  NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,
     +                NPHIA,NNB,NTYPES,NRC,NCONP,MAXMEM,NWDVAR,MAXNB,
     +                MBONA,MTHETA,MPHIA,NATC,IBELLY,NATBEL,ISHAKE,NMXRS
      COMMON /PARMS/  RK(500),REQ(500),TK(900),TEQ(900),PK(900),
     +                PN(900),PHASE(900),CN1(1830),CN2(1830),RAD(61),
     +                SOLTY(60),GAMC(900),GAMS(900),IPN(900),FMN(900),
     +                ONE_SCEE(900),ONE_SCNB(900)
      COMMON /NBPARM/ IGRAPH(MAXATOM), CHRG(MAXATOM), AMASS(MAXATOM),
     +                IAC(MAXATOM), ICO(1830), LABRES(MAXRES),
     +                IPRES(MAXRES), ISYMBL(MAXATOM), ITREE(MAXATOM),
     +                JOIN(MAXATOM), IROTAT(MAXATOM)
      COMMON /TOPA/ NSPSOL, NSPM, NATOM

      logical          verbose,inp_exists
      character*80     string, buffer,part

      parameter        (maxcall=100)
      common /timings/ ftme(maxcall),part(maxcall)
      common /variables/ nvt,nvar,nad

      ntime = 0
      ttime = 0.0d0
      call setime

c-------------------------------------------------------------------------------
c     PRINT HEADER
c-------------------------------------------------------------------------------

      write(6,'(1x,79a,/)') ('_',i=1,79)
      write(6,'(/,30x,''C  O  M  B  I  N  E'',/)')
      write(6,'(/,22x,''COMparative BINding Energy analysis'',/)')
      write(6,'(1x,79a,/)') ('_',i=1,79)

c-------------------------------------------------------------------------------
c --- Read comand line & input file
c-------------------------------------------------------------------------------

      inp_exists = .false.
      narg = iargc()
      if ((narg .gt. 4).or.(narg .lt.4)) call help
      do i = 1, narg, 2
         call getarg(i,string)
         if (string(1:1) .eq. '-') then
            if (string(2:2).eq.'i') then
               call getarg(i+1,buffer)
               read(buffer,'(a)') inpfile
               inquire(file=inpfile,exist=inp_exists)
               if (inp_exists) then
                  open(unit=8,file=inpfile,status='old')
               else
                  write(6,'(/,a,/)') 'ERROR: input file does not exist.'
                  stop
               endif
            else if (string(2:2).eq.'o') then
               call getarg(i+1,buffer)
               read(buffer,'(a)') outfile
               open(unit=16,file=outfile)
            else 
               call help
            endif
         else
            call help
         endif
      enddo

      call readinp(8,idel,nsc,imat,nlv,ielec,ncomp,ntest,nad,
     &     randtest,ncross,dielect,ptr,cutptr,comp,drugname,y,add)
            
      open(unit=20,file='combine.log')
      
      call storetime('Initialization:       ',ttime,ntime)
      
c-------------------------------------------------------------------------------
c --- Read topology and coordinates of the complexes and compute X-matrix
c-------------------------------------------------------------------------------

      call setime

      if (imat.eq.0) then
         call build_xmatrix(ielec,dielect,ncomp,ntest,comp,drugname,
     &        name,xraw,y,add,dres,multres,watMol)
         call golpe_matrix(nvt,ncomp,ntest,xraw,y,name,comp,
     &        imat,labres,nad,multres,watMol)
      else
         call read_xmatrix(ncomp,ntest,comp,drugname,name,add,xraw,y,
     &        nres,labres,multres)
      endif

c      verbose = .true.
c      if (verbose) call golpe_matrix(nvt,ncomp,ntest,xraw,y,name,comp,
c     &     imat,labres,nad)
      
c --- Data Pretreatment

      if (ptr.ne.0 .or. cutptr.ne.0.0) then
         call pretreat_xmatrix(nvt,nvtold,nad,ncomp+ntest,
     &        xraw,ptr,cutptr,varselect)
      else
         nvtold = nvt
         write (6,91)
 91      format (/,' Pretreatment NO ',/)
         do i=1,nvt
            varselect(i)=1
         enddo
      endif

      k=0
      j=0
      
      do i=1,nres
         m=0
         do l=1,multres
            if (i.eq.dres-multres+l)then
               m=1
            endif
         enddo
         if(m.eq.0)then
            if(labres(i).ne.drugname(ncomp+ntest))then
               j=j+1
               k=k+1
               call int2str(k,chartemp)
               vardesc(j)=
     &              labres(i)(ftrim(labres(i)):btrim(labres(i)))
     &              //' '//chartemp(ftrim(chartemp):btrim(chartemp))
     &              //' '//'vdw'
            endif
         endif
      enddo
      k=0
      do i=1,nres
         m=0
         do l=1,multres
            if (i.eq.dres-multres+l)then
               m=1
            endif
         enddo
         if(m.eq.0)then
            if(labres(i).ne.drugname(ncomp+ntest))then
               j=j+1
               k=k+1
               call int2str(k,chartemp)
               vardesc(j)=
     &              labres(i)(ftrim(labres(i)):btrim(labres(i)))
     &              //' '//chartemp(ftrim(chartemp):btrim(chartemp))
     &              //' '//'ele'
            endif
         endif
      enddo
      do i=j+1,nvtold
         vardesc(i) = 'add'        
      enddo
c --- 

      call storetime('X-matrix computation: ',ttime,ntime)

c-------------------------------------------------------------------------------
c --- Apply PLS to the X-matrix
c-------------------------------------------------------------------------------

      call setime
      nrow  = ncomp
      ncol  = nvt
      nrank = min(ncol,nrow-1,maxrank,nlv)

      write (6,100)  nrow,ncol
  100 format (/,' Number of Data Values :',i6,8x,
     &          ' Number of Property Columns :',i6)
      if (ntest .ne. 0) then
         write (6,110)  ntest
  110    format (' Number of Test Values :',i6)
      end if

      call qsar(nlv,idel,nsc,nrank,nrow,ncol,ntest,xraw,y,contour,cutinp
     &     ,randtest,ncross,name,comp,coef,varselect,nvtold,vardesc,b)
      call storetime('PLS regression:       ',ttime,ntime)

c-------------------------------------------------------------------------------
c --- Write coefficient information for display in GRASP. The pdb
c --- file of the complex with maximum activity is written. First,
c --- store maximum activity index, then read again the info from
c --- topology and rst files and transform to corresponding PDB file.
c-------------------------------------------------------------------------------

      call setime

      if (imat.eq.0) then

c --- set up for VDW & ELE. Check out water molecules, other var.s and
c --- the ligand itself when transfering from ncol to nres

         ymax = -100.0d0
         imax = 0
         do i = 1, ncomp+ntest
            if (y(i).ge.ymax) then
               imax = i
            ymax = y(i)
         endif
      enddo
      
      topfile = comp(imax)(1:lchar(comp(imax)))//'.top'
      open(unit=10,file=topfile,status='old')
      rstfile = comp(imax)(1:lchar(comp(imax)))//'.crd'
      open(unit=11,file=rstfile,status='old')
      
      call readtop(10,20)
      call readxyz(16,11,xyz)
      call seldrug(20,drugname(imax),dres,natom,
     &     multres,watMol,ndatm,datm)
      
      close(10)
      close(11)
        
c-------------------------------------------------------------------------------
c
c     Pdb files generations with the extra column "b factor" for
c     the representation of the most representative vdw and electrostatics 
c
c------------------------------------------------------------------------------
        
      call coefffiles(b,nrank,nvar,nres,dres,ipres
     &       ,igraph,labres,xyz,multres,varselect)
        
c --- here the grasp file with VDW coefs. is generated

        open(unit=10,file='graspvdw.pdb')
        write(10,'(a)') 'GRASP PDB FILE'
        write(10,'(a)') 'FORMAT NUMBER= 3'
        open(unit=30,file='vdwmin.pdb')
        write(30,'(a)') 'GRASP PDB FILE'
        write(30,'(a)') 'FORMAT NUMBER= 3'
        cmin =  10.0
        cmax = -10.0
        do i = 1, nvar
          cmin = min(coef(i),cmin)
          cmax = max(coef(i),cmax)
        enddo
        do i = 1, nvar
          coefs(i) = -(coef(i)-cmin)/(cmax-cmin)*10.0d0
        enddo
        nvr = 0
        do 15 i = 1, nres
          if (i.eq.dres) then
            init=ipres(i)
            iend=ipres(i+1)-1
            do k = init, iend
              write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &        'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
     &        '  1.00', 0.0d0
              write(30,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &        'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
     &        '  1.00', 0.0d0
            enddo
            goto 15
          endif
          nvr=nvr+1
          init=ipres(i)
          iend=ipres(i+1)-1
          do 16 k = init, iend
            if (igraph(k)(1:1).eq.'H') goto 16
            write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &      'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
     &      '  1.00',coefs(nvr)
            if (abs(coef(nvr)).ge.0.001) then
              write(30,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &        'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
     &        '  1.00',coefs(nvr)
            endif
   16     continue
   15   continue
        close(10)
        close(30)
    
c --- here the grasp file with ELE coefs. is generated

        open(unit=10,file='graspele.pdb')
        write(10,'(a)') 'GRASP PDB FILE'
        write(10,'(a)') 'FORMAT NUMBER= 3'
        open(unit=30,file='elemin.pdb')
        write(30,'(a)') 'GRASP PDB FILE'
        write(30,'(a)') 'FORMAT NUMBER= 3'
        cmin =  10.0
        cmax = -10.0
        do i = nvar+1, 2*nvar
          cmin = min(coef(i),cmin)
          cmax = max(coef(i),cmax) 
        enddo
        do i = nvar+1, 2*nvar
           coefs(i) = -(coef(i)-cmin)/(cmax-cmin)*10.0d0
        enddo
        do 25 i = 1, nres
          if (i.eq.dres) then
            init=ipres(i)
            iend=ipres(i+1)-1
            do k = init, iend
              write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &        'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
     &        '  1.00', 0.0d0
              write(30,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &        'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
     &        '  1.00', 0.0d0
            enddo
            goto 25
          endif
          nvr=nvr+1
          init=ipres(i)
          iend=ipres(i+1)-1
          do 26 k = init, iend
            if (igraph(k)(1:1).eq.'H') goto 26
            write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &      'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
     &      '  1.00',coefs(nvr)
            if (abs(coef(nvr)).ge.0.001) then
              write(30,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &        'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
     &        '  1.00',coefs(nvr)
            endif
   26     continue
   25   continue
        close(10)
        close(30)

      endif

      call storetime('PDB & coeff. output:  ',ttime,ntime)

c-------------------------------------------------------------------------------
c --- End of the pHP_m002.toprogram
c-------------------------------------------------------------------------------

      call printime(ntime)
      stop  
      end
c
c===============================================================================
c
      subroutine build_xmatrix(ielec,dielect,ncomp,ntest,comp,drugname,
     &     name,xraw,y,add,dres,multres,watMol)

c-------------------------------------------------------------------------------
c --- Computes the X-matrix from the set of complexes
c-------------------------------------------------------------------------------

      implicit     double precision (a-h,o-z)

      parameter    (maxatdrg = 200)
      parameter    (maxcompl = 100)
      PARAMETER    (MAXATOM = 8000)
      PARAMETER    (MAXRES  =  700)

      parameter    (maxrow= 100)
      parameter    (maxcol=1500)
      parameter    (maxadd=  10)
      parameter    (maxrank =20)

      dimension     xyz(3,maxatom),add(maxcompl,maxadd)
      dimension     xraw(maxrow,maxcol),y(maxrow)
      dimension     ele(maxres), vdw(maxres)
      character*4   drug, drugname(maxcompl), name(maxrow)
      integer       dres, datm(maxatdrg), ratm(maxatdrg)
      integer       ielec
      integer       multres,watMol

      real*8        dielect

      CHARACTER*4   IGRAPH, LABRES
      character*40  topfile,rstfile,delfile
      character*40  comp(maxcompl)

      logical  inp_exists

      COMMON /INFOV/  NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,
     +                NPHIA,NNB,NTYPES,NRC,NCONP,MAXMEM,NWDVAR,MAXNB,
     +                MBONA,MTHETA,MPHIA,NATC,IBELLY,NATBEL,ISHAKE,NMXRS
      COMMON /PARMS/  RK(500),REQ(500),TK(900),TEQ(900),PK(900),
     +                PN(900),PHASE(900),CN1(1830),CN2(1830),RAD(61),
     +                SOLTY(60),GAMC(900),GAMS(900),IPN(900),FMN(900),
     +                ONE_SCEE(900),ONE_SCNB(900)
      COMMON /NBPARM/ IGRAPH(MAXATOM), CHRG(MAXATOM), AMASS(MAXATOM),
     +                IAC(MAXATOM), ICO(1830), LABRES(MAXRES),
     +                IPRES(MAXRES), ISYMBL(MAXATOM), ITREE(MAXATOM),
     +                JOIN(MAXATOM), IROTAT(MAXATOM)
      COMMON /TOPA/ NSPSOL, NSPM, NATOM
      common /variables/ nvt,nvar,nad

      do ic = 1, ncomp+ntest
c        To print the complex name in the log file too
         write(16,*)'Complex : ',comp(ic)

         topfile = comp(ic)(1:lchar(comp(ic)))//'.top'
         write(6,'(a8,x,a13,x,a4)') 'Reading',topfile,'file'
         inquire(file=topfile,exist=inp_exists)
         if (inp_exists) then
            open(unit=10,file=topfile,status='old')
         else
            write(6,'(/,a,/)') 'Error: top file does not exist'
            stop
         endif        
         rstfile = comp(ic)(1:lchar(comp(ic)))//'.crd'
         write(6,'(a8,x,a13,x,a4)') 'Reading',rstfile,'file'
         inquire(file=rstfile,exist=inp_exists)
         if (inp_exists) then
            open(unit=11,file=rstfile,status='old')
         else
            write(6,'(/,a,/)') 'Error: crd file does not exist'
            stop
         endif               
         
         if (ielec .eq. 3) then
            delfile = comp(ic)(1:lchar(comp(ic)))//'.dph'
            inquire(file=delfile,exist=inp_exists)
            if (inp_exists) then
               open(unit=12,file=delfile,status='old')
            else
               write(6,'(/,a,/)') 'Error: reading dph file'
               stop
            endif               
         endif

         call readtop(10,20)        
         call readxyz(16,11,xyz)
         
         call wrtinfo(20,natom,igraph,xyz,chrg,rad,iac)
         drug = drugname(ic)
         name(ic) = drug
         call seldrug(20,drug,dres,natom,multres,watMol,ndatm,datm)
         eletot = 0.0d0
         vdwtot = 0.0d0
         nvar = 0

         do 10 k = 1, nres
            if (multres.eq.0)then
               if (k .eq. dres) goto 10
            else
               do l=1,multres
                  if (k .eq. dres+1-l)goto 10
               enddo
            endif
            init = ipres(k)
            iend = ipres(k+1)-1
            if (k.eq.nres) iend = natom
            nratm = 0
            do i = init, iend
               nratm = nratm + 1
               ratm(nratm) = i
            enddo
            
            s=0

            call eint(20,xyz,ndatm,datm,nratm,ratm,eele,evdw,natom,
     &           ielec,k,labres(k),dielect)
            nvar = nvar + 1
            ele(nvar) = eele
            vdw(nvar) = evdw
            eletot = eletot + eele
            vdwtot = vdwtot + evdw
            write(16,'(i4,2x,a4,2f10.5)') k, labres(k), eele, evdw
 10      continue


         nvt = 0
         do k = 1, nvar
            nvt = nvt + 1
            xraw(ic,nvt) = vdw(k)
            if (xraw(ic,nvt).eq.Z'FFFFFFFF') xraw(ic,nvt)=0.0d0
         enddo
         do k = 1, nvar
            nvt = nvt + 1
            xraw(ic,nvt) = ele(k)
            if (xraw(ic,nvt).eq.Z'FFFFFFFF') xraw(ic,nvt)=0.0d0
         enddo
         do k = 1, nad 
            nvt = nvt + 1
            xraw(ic,nvt) = add(ic,k)
            if (xraw(ic,nvt).eq.Z'FFFFFFFF') xraw(ic,nvt)=0.0d0
         enddo
         
         write(16,11)'Total electrostatic interaction energy : ',eletot
         write(16,12)'Total Van der Waals interaction energy : ',vdwtot
         do k = 1, nad
            write(16,*) 'nad ',k,':',add(ic,k)
         enddo
 11      format(/,a,f10.4)
 12      format(a,f10.4,/)
         close(10)
         close(11)
         close(12)
      enddo

      return
      end
c===============================================================================
c
      subroutine read_xmatrix(ncomp,ntest,comp,drugname,name,add,xraw,y,
     &     nres,labres,multres)

c-------------------------------------------------------------------------------
c --- Read the X-matrix from an external file
c-------------------------------------------------------------------------------

      implicit     double precision (a-h,o-z)

      parameter    (maxrow= 100)
      parameter    (maxcol=1500)
      parameter    (maxcompl = 100)
      parameter    (maxadd=  10)
      parameter    (maxres=700)

      dimension     xraw(maxrow,maxcol),y(maxrow)
      dimension     add(maxcompl,maxadd)
      character*4   name(maxrow)
      character*40  comp(maxrow)
      character*4   drugname(maxcompl)
      character*40  comp_tmp(maxrow)
      character*4   labres(maxres)

      character*80 borrar

      integer  l,multres

      logical  efile_exists

      common /variables/ nvt,nvar,nad     

      write(6,'(a,/)') 'Reading intractions from previous analysis:'

      efile_exists = .false.
      inquire(file='combine.dat',exist=efile_exists)
      if (efile_exists) then
         do l=1,ncomp+ntest

            write(6,'(a7,x,a7,x,a12)')'Reading',comp(l),'interactions'
            
            open(unit=10,file='combine.dat',status='old')
            
            read(10,*)
            read(10,*) n
            read(10,*) m

            if(ncomp+ntest.gt.m)then
               write(*,*) 'ERROR: Number of complexes do not match'
               write(*,*) '       Current:',ncomp+ntest,' Loaded:',m
               
               stop
            endif

            nres= (n-nad)/2

c$$$  if (m.ne.(ncomp+ntest)) then
c$$$  write(6,'(/,a,/)')'ERROR: Number of complexes do not match'
c$$$  stop
c$$$  endif
            
            do i = 1, m
               read(10,*)
               read(10,'(a8)') comp_tmp(i)

c               write(*,'(3a8)') comp_tmp(i),comp(l),'*******'

               do j = 1, n-1
c----------------------------------------------------------
c     Modified for matching the precision of the .in file
c     
c     read(10,'(f10.5)') xraw(i,j)
c----------------------------------------------------------                 
                  if(comp_tmp(i).eq.comp(l))then

                     if(j.le.nres)then
                        read(10,'(f12.6,3x,a4)') xraw(l,j),labres(j)
                     else
c                        read(10,*) borrar
c                        print*, borrar
                        read(10,'(f12.6)') xraw(l,j)
                     endif

c                     print*,xraw(l,j)

                  else
                     read(10,*)
                  endif
               enddo
             
               do k = 1, nad
                  if(comp_tmp(i).eq.comp(l)) xraw(l,n-1+k)=add(i,k)
               enddo
c----------------------------------------------------------
c     Modified for matching the precision of the .in file
c     
c     read(10,'(f10.5)') ytmp
c----------------------------------------------------------
               if(comp_tmp(i).eq.comp(l))then
                  read(10,'(f12.6)') ytmp
               else
                  read(10,*)
                  ytmp=0
               endif
               if(comp_tmp(i).eq.comp(l))then
                  if (ytmp.ne.y(l)) then
                     write(6,'(/,a,/)')'ERROR: Activities do not match.'
                     stop
                  endif     
               endif
            enddo
            close(10)
         enddo       
      else
         write(6,'(/,a,/)')'ERROR: energy file does not exist.'
         stop
      endif
      
      nvt = n-1 
      nvar= (n-1)/2
      
      do i=1,ncomp+ntest
         name(i)=drugname(i)        
      enddo


!!!!!!!!!!!!!!!!!!!

      open (30,file='energy_values.dat')
      write(30,'(a)') 'COMBINE analysis '
      write(30,'(i5)') nvt+1
      write(30,'(i5)') ncomp+ntest
      do i = 1, ncomp+ntest
c     write(30,'(i5)') i
         write(30,'(a5,i5,3x,a8)') '#####',i,comp(i)
         do j = 1, nvt
            write(30,'(i5,3x,f12.6)') j,xraw(i,j)
         enddo
         write(30,'(a5,3x,f12.6)') 'act',y(i)
      enddo
      close(30)
      
!!!!!!!!!!!!!!!!!!!!!
            
      return
      end
c     
c===============================================================================
c
      subroutine golpe_matrix(nvt,ncomp,ntest,xraw,y,name,comp,imat,
     &     labres,nad,multres,watMol)

c-------------------------------------------------------------------------------
c --- Writes the energy partition in golpe-4.5 input format
c-------------------------------------------------------------------------------

      implicit     double precision (a-h,o-z)

      parameter    (maxrow=100)
      parameter    (maxcol=1500)
      PARAMETER    (MAXRES=700)

      integer      nvt,ncomp,ntest,multres,watMol,l
      dimension    xraw(maxrow,maxcol),y(maxrow)
      character*4  name(maxrow)
      character*40 comp(maxrow)
      CHARACTER*4  LABRES(maxres)
  
      open (30,file='energy_values.dat')
      write(30,'(a)') 'COMBINE analysis '
      write(30,'(i5)') nvt+1
      write(30,'(i5)') ncomp+ntest
      do i = 1, ncomp+ntest
c     write(30,'(i5)') i
         write(30,'(a5,i5,3x,a8)') '#####',i,comp(i)
         do j = 1, nvt
            write(30,'(i5,3x,f12.6)') j,xraw(i,j)
         enddo
         write(30,'(a5,3x,f12.6)') 'act',y(i)
      enddo
      close(30)

      open (30,file='combine.dat')
      write(30,'(a)') 'COMBINE analysis '
      write(30,'(i5)') nvt+1
      write(30,'(i5)') ncomp+ntest
      do i = 1, ncomp+ntest
         l=1     
         write(30,'(i5)') i
         write(30,'(a8)') comp(i)
         do j = 1, nvt
            if(j.le.((nvt-nad)/2)-watMol)then
               write(30,'(f12.6,3x,a4,i4)') xraw(i,j),labres(j),j
            elseif(j.le.(nvt-nad)/2)then
               write(30,'(f12.6,3x,a4,i4)')
     &              xraw(i,j),labres(j+multres),j
            elseif(j.le.nvt-nad-watMol)then
               write(30,'(f12.6,3x,a4,i4)') xraw(i,j),
     &              labres(j-(nvt-nad)/2),j-(nvt-nad)/2
            elseif( j .le.nvt-nad )then
              write(30,'(f12.6,3x,a4,i4)') xraw(i,j),
     &              labres( j- (nvt-nad)/2 +multres ),j-(nvt-nad)/2
            else
               write(30,'(f12.6,3x,a2,i2)') xraw(i,j),
     &              'EV',l
               l=l+1
            endif
         enddo
         write(30,'(f12.6,3x,a4)') y(i),'EXP '
      enddo
      close(30)
      
      return
      end
c
c===============================================================================
c
      subroutine wrtinfo(no,natom,IGRAPH,xyz,CHRG,rad,iac)

      IMPLICIT      DOUBLE PRECISION (A-H,O-Z)

      PARAMETER    (MAXATOM = 8000)

      CHARACTER*4   IGRAPH(MAXATOM)
      dimension     xyz(3,maxatom),CHRG(MAXATOM),IAC(MAXATOM),RAD(61)

      do j = 1, natom
        write(no,8) j,IGRAPH(j),(XYZ(I,J),I=1,3),CHRG(j),
     +              rad(iac(j)),IAC(j)
      enddo
    8 format(1x,i4,2x,a4,3(1x,1F10.5),2(1x,1f8.3),x,i6)

      return
      end
c
c===============================================================================
c
      subroutine eint(no,xyz,ndatm,datm,nratm,ratm,eele,evdw,
     &                natom,ielec,kres,cres,dielect)

c_______________________________________________________________________________
c
c    Computation of the residue-residue non-bonded interactions. Interactions
c    are computed according to the Cornell et al. 1995 AMBER force field:
c
c               _                            _
c         ___  |    A        C      q1 q2     |
c         \    |   ---   -  ---  +  -----     |
c         /__  |   r^12     r^6      D r      |
c               -                            -
c        nonbonded   
c        pairs
c
c    This is the basic additive form as in the Cornell et al. 1995
c    force field, i.e. omitting polarization as well as the hydrogen 
c    bonding 10-12 term from the Weiner et al. 1984,1986 force field.
c    
c    Topology information still has the pointers to the 10-12 interactions
c    when some vdw's are not computed. This is checked here, btu according
c    to the Cornell et al. philosophy it is assumed that the 10-12 term is 
c    zero, and hence it is not computed.
c_______________________________________________________________________________

      IMPLICIT      DOUBLE PRECISION (A-H,O-Z)

      parameter    (maxatdrg = 200)
      PARAMETER    (MAXATOM = 8000)
      PARAMETER    (MAXRES  =  700)
      parameter    (ratio   = (4.0d0-80.0d0)/(4.0d0+80.0d0) )

      dimension     d(maxatdrg,maxatdrg),d2(maxatdrg,maxatdrg)
      dimension     xyz(3,maxatom),s1(maxatdrg),s2(maxatdrg)
      integer       datm(maxatdrg), ratm(maxatdrg)

      CHARACTER*4   IGRAPH, LABRES, cres, tres
     
      real*8 lambda, dielect

      COMMON /NBPARM/ IGRAPH(MAXATOM), CHRG(MAXATOM), AMASS(MAXATOM),
     +                IAC(MAXATOM), ICO(1830), LABRES(MAXRES),
     +                IPRES(MAXRES), ISYMBL(MAXATOM), ITREE(MAXATOM),
     +                JOIN(MAXATOM), IROTAT(MAXATOM)
      COMMON /PARMS/  RK(500),REQ(500),TK(900),TEQ(900),PK(900),
     +                PN(900),PHASE(900),CN1(1830),CN2(1830),RAD(61),
     +                SOLTY(60),GAMC(900),GAMS(900),IPN(900),FMN(900),
     +                ONE_SCEE(900),ONE_SCNB(900)
      COMMON /INFOV/  NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,
     +                NPHIA,NNB,NTYPES,NRC,NCONP,MAXMEM,NWDVAR,MAXNB,
     +                MBONA,MTHETA,MPHIA,NATC,IBELLY,NATBEL,ISHAKE,NMXRS

c-------------------------------------------------------------------------------
c First store the pairwise interactions between ligand and residue
c-------------------------------------------------------------------------------

      do i = 1, ndatm
        ii = datm(i)
        xi = xyz(1,ii)
        yi = xyz(2,ii)
        zi = xyz(3,ii)
        do j = 1, nratm
          jj = ratm(j)
          xj = xyz(1,jj)
          yj = xyz(2,jj)
          zj = xyz(3,jj)
          dx = xi-xj 
          dy = yi-yj
          dz = zi-zj
          d2(i,j) = dx*dx+dy*dy+dz*dz
          d (i,j) = dsqrt(d2(i,j))
        enddo
      enddo

c-------------------------------------------------------------------------------
c Compute electrostatics interactions. Options:
c   ielec = 0 => constant dielectric with eps=4
c   ielec = 1 => images model from Goodford
c   ielec = 2 => distance dependent dielectric constant
c   ielec = 3 => Poisson-Boltzmann electrostatics, readrom file
c   ielec = 4 => Sigmoidal dielectric function
c Default is distance dependent
c-------------------------------------------------------------------------------

      et  = 0.0d0
      if (ielec .eq. 0) then
         eps = dielect
         do i = 1, ndatm
            ii = datm(i)
            do j = 1, nratm
               jj = ratm(j)
               e  = 332.0d0*((chrg(ii)*chrg(jj))/(eps*d(i,j)))
               et = et + e
            enddo
         enddo
      else if (ielec .eq. 1) then
         eps = dielect
         call images(ndatm,datm,natom,xyz,s1)
         call images(nratm,ratm,natom,xyz,s2)
         do i = 1, ndatm
            ii = datm(i)
            do j = 1, nratm               
               jj = ratm(j)
               e1 = (chrg(ii)*chrg(jj))/eps
               e2 = (1.0d0/d(i,j))+(ratio/dsqrt(d2(i,j)+
     &              (eps*s1(i)*s2(j))))
               e  = 332.0d0*e1*e2
               et = et + e
            enddo
         enddo
      else if (ielec .eq. 2) then
         do i = 1, ndatm
            ii = datm(i)
            do j = 1, nratm
               jj = ratm(j)
               e  = 332.0d0*((chrg(ii)*chrg(jj))/d2(i,j))
               et = et + e
            enddo
         enddo
      else if (ielec .eq. 3) then
         do m = 1, maxres
            read(12,'(a4,i4,f12.4)',end=6) tres, mres, etmp
            if (mres .eq. kres .and. tres .eq. cres) then
               et = etmp
               goto 10
            endif
         enddo
 6       continue
         do i = 1, ndatm
            ii = datm(i)
            do j = 1, nratm
               jj = ratm(j)
               e  = 332.0d0*((chrg(ii)*chrg(jj))/d2(i,j))
               et = et + e
            enddo
         enddo
 10      continue
         rewind(12)
      else if (ielec .eq. 4) then
         eps = dielect
         lambda = 1.0367/(eps+1.0)
         do i = 1, ndatm
            ii = datm(i)
            do j = 1, nratm
               jj = ratm(j)
               ex = 0.0d0
               tmpeps = 0.0d0
               ex = -lambda*(eps+1)*d(i,j)
               tmpeps = (( eps+1 )/ ( 1 + eps*exp(ex) )) - 1
               e  = 332.0d0*((chrg(ii)*chrg(jj))/(tmpeps*d(i,j)))
               et = et + e
            enddo
         enddo
      else
         do i = 1, ndatm
            ii = datm(i)
            do j = 1, nratm
               jj = ratm(j)
               e  = 332.0d0*((chrg(ii)*chrg(jj))/d2(i,j))
               et = et + e
            enddo
         enddo
      endif

c-------------------------------------------------------------------------------
c Compute van der Waals interactions (6-12)
c-------------------------------------------------------------------------------

      vt  = 0.0d0
      do i = 1, ndatm
         ii = datm(i)
         do j = 1, nratm
            jj = ratm(j)
            if (ii.lt.jj) then
               index = ico(ntypes*(iac(ii)-1)+iac(jj))
            else
               index = ico(ntypes*(iac(jj)-1)+iac(ii))
            endif
            if (index .gt. 0) then              
               r6 = d2(i,j)**3
               r12= r6**2
               v  = (cn1(index)/r12) - (cn2(index)/r6) 
               vt = vt + v
            endif
            if (v.gt. 10.) then
               write(16,*)
               write(16,*) 'WARNIN7: Atomic clash: '
               write(16,*) d(i,j),chrg(ii),chrg(jj),v
               write(16,*) index, cn1(index), cn2(index), r12, r6 
               write(16,*)
            endif
         enddo
      enddo
      eele = et
      evdw = vt
      
      return
      end
c
c===============================================================================
c
      subroutine images(nar,narl,natom,xyz,s)

c-------------------------------------------------------------------------------
c --- Here the Goodford method (based on the images approximation) for 
c --- the calculation of the electrostatics interactions is used. See:
c --- Goodford JMedChem 1985,Vol28,N7,851
c-------------------------------------------------------------------------------

      IMPLICIT      DOUBLE PRECISION (A-H,O-Z)

      parameter    (maxatdrg = 200)
      PARAMETER    (MAXATOM = 8000)

      dimension     xyz(3,maxatom),s(maxatdrg)
      integer       narl(maxatdrg)

      do 11 i = 1, nar   
        ncon = 0
        ii = narl(i)
        xi = xyz(1,ii)
        yi = xyz(2,ii)
        zi = xyz(3,ii)
        do 10 j = 1, natom
          if (ii.eq.j) goto 10
          xj = xyz(1,j)
          yj = xyz(2,j)
          zj = xyz(3,j)
          dx = xi-xj
          dy = yi-yj
          dz = zi-zj
          d  = dsqrt(dx*dx+dy*dy+dz*dz)
          if (d.le.4.0d0) then
            ncon = ncon + 1
            if (ncon.ge.12) then
              s(i) = 4.00d0
              goto 11
            endif
          endif
   10   continue
        if (ncon .le.  6) s(i) = 0.00d0
        if (ncon .eq.  7) s(i) = 0.40d0
        if (ncon .eq.  8) s(i) = 0.90d0
        if (ncon .eq.  9) s(i) = 1.40d0
        if (ncon .eq. 10) s(i) = 1.90d0
        if (ncon .eq. 11) s(i) = 2.60d0
   11 continue
c     write(6,*) (s(i),i=1,ndatm)

      return
      end
c
c===============================================================================
c
      subroutine seldrug(no,drug,dres,natom,multres,watMol,ndatm,datm)

c-------------------------------------------------------------------------------
c --- Select the drug molecule
c-------------------------------------------------------------------------------

      IMPLICIT      DOUBLE PRECISION (A-H,O-Z)

      parameter    (maxatdrg = 200)
      PARAMETER    (MAXATOM = 8000)
      PARAMETER    (MAXRES  =  700)
      parameter    (NUM_AMINOACIDS=58)

      character*3 amiNames(NUM_AMINOACIDS)

      CHARACTER*4   IGRAPH, LABRES
      COMMON /NBPARM/ IGRAPH(MAXATOM), CHRG(MAXATOM), AMASS(MAXATOM),
     +                IAC(MAXATOM), ICO(1830), LABRES(MAXRES),
     +                IPRES(MAXRES), ISYMBL(MAXATOM), ITREE(MAXATOM),
     +                JOIN(MAXATOM), IROTAT(MAXATOM)
      COMMON /INFOV/  NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,
     +                NPHIA,NNB,NTYPES,NRC,NCONP,MAXMEM,NWDVAR,MAXNB,
     +                MBONA,MTHETA,MPHIA,NATC,IBELLY,NATBEL,ISHAKE,NMXRS

      character*4    drug
      integer        dres, datm(maxatdrg)
      integer        l,compRes,multres,watMol
      
      data (amiNames(l),l=1,58) /
     &     'ALA','GLY','SER','THR','LEU',
     &     'ILE','VAL','ASN','GLN','ARG',
     &     'HID','HIP','HIE','TRP','PHE',
     &     'TYR','GLU','ASP','ASH','LYS',
     &     'PRO','CYS','CYX','MET','HIS',
     &     'GLH','IP3','ION','CAA','FMN',
     &     'PRG','ACE','NME','gdp','gtp',
     &     'adp','atp','NAD','NAH','NPD',
     &     'NPH','ARP','HEM','1GA','PGA',
     &     '80G','GUA','PYR','MA ','AN5',
     &     'DIQ','EPN','CAL','C3N','MOH',
     &     'CL3','NMA','WAT'/

      ndatm=0
      multres=0
      watMol=0

      if(drug.eq.'PAR')then
         multres=0
         do k = 1, nres
            compRes=0
            do l=1,58
               if(labres(k).eq.amiNames(l))then
                  compRes=compRes+1
               endif
            enddo

            if (labres(k).eq.amiNames(58))watMol=watMol+1

            if(compRes.eq.0)then
               dres  = k
               if (ipres(dres).lt.ipres(dres+1))then
                  do i = ipres(dres), ipres(dres+1)-1
                     ndatm = ndatm + 1
                     if (ndatm.gt.maxatdrg) then
                        write(no,*) 'seldrug: too many atoms in drug '
                        stop
                     endif
                     datm(ndatm) = i
c                     print*,labres(k),datm(ndatm),ndatm,natom
                  enddo
               else
                  do i = ipres(dres),natom
                     ndatm = ndatm + 1
                     if (ndatm.gt.maxatdrg) then
                        write(no,*) 'seldrug: too many atoms in drug '
                        stop
                     endif
                     datm(ndatm) = i
c                     print,labres(k),datm(ndatm),ndatm,natom
                  enddo
               endif        
               multres=multres+1 !number of partitionated lig  
            endif
         enddo
      else
         ndatm = 0    
         do k = 1, nres
            if (labres(k).eq.amiNames(58))then
               watMol=watMol+1
            endif
            if (labres(k).eq.drug) then 
               dres  = k
               ndatm = 0
               if (ipres(dres).lt.ipres(dres+1))then
                  do i = ipres(dres), ipres(dres+1)-1
                     ndatm = ndatm + 1
                     if (ndatm.gt.maxatdrg) then
                        write(no,*) 'seldrug: too many atoms in drug '
                        stop
                     endif
                     datm(ndatm) = i
                  enddo
               else
                  do i = ipres(dres),natom
                     ndatm = ndatm + 1
                     if (ndatm.gt.maxatdrg) then
                        write(no,*) 'seldrug: too many atoms in drug '
                        stop
                     endif
                     datm(ndatm) = i
                  enddo
               endif
            endif
         enddo
         if(ndatm.eq.0)then
            write(no,*) 'seldrug: Drug molecule ',drug,
     &           'not found in topology file '
            stop
         endif
         multres=1
      endif

      return
      end


c     
c===============================================================================
c
      SUBROUTINE READXYZ(no,nf,xyz)

      IMPLICIT      DOUBLE PRECISION (A-H,O-Z)

      PARAMETER    (MAXATOM = 8000)

      dimension    xyz(3,maxatom)

      COMMON /TOPA/ NSPSOL, NSPM, NATOM

      read(nf,*)
      read(nf,*) nat
      if (nat.ne.natom) then
        write(no,'(/,5x,a)') 'readxyz: topology and coords do not match'
        stop
      endif
      if (nat.gt.maxatom) then
        write(no,'(/,5x,a)') 'readxyz: a parameter array overflowed'
        stop
      endif
      READ(NF,'(6F12.7)') ((XYZ(I,J),I=1,3),J=1,NATOM)

      return
      end
c
c===============================================================================
c
      SUBROUTINE READTOP(NF,NO)
c_______________________________________________________________________________
c
c VARIABLES:
c ==========
c
c NATOM  : total number of atoms
c NTYPES : total number of distinct atom types
c NBONH  : number of bonds containing hydrogen
c MBONA  : number of bonds not containing hydrogen
c NTHETH : number of angles containing hydrogen
c MTHETA : number of angles not containing hydrogen
c NPHIH  : number of dihedrals containing hydrogen
c MPHIA  : number of dihedrals not containing hydrogen
c NHPARM : currently not used
c NPARM  : currently not used
c NEXT   : number of excluded atoms
c NRES   : number of residues
c NBONA  : MBONA + number of constraint bonds
c NTHETA : MTHETA + number of constraint angles
c NPHIA  : MPHIA + number of constraint dihedrals
c NUMBND : number of unique bond types
c NUMANG : number of unique angle types
c NPTRA  : number of unique dihedral types
c NATYP  : number of atom types in parameter file, see SOLTY below
c NPHB   : number of distinct 10-12 hydrogen bond pair types
c IFPERT : set to 1 if perturbation info is to be read in
c NBPER  : number of bonds to be perturbed
c NGPER  : number of angles to be perturbed
c NDPER  : number of dihedrals to be perturbed
c MBPER  : number of bonds with atoms completely in perturbed group
c MGPER  : number of angles with atoms completely in perturbed group
c MDPER  : number of dihedrals with atoms completely in perturbed groups
c IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
c NMXRS  : number of atoms in the largest residue
c IFCAP  : set to 1 if the CAP option from edit was specified
c
c
c ARRAYS:
c =======
c
c IGRAPH : the user atoms names
c CHRG   : the atom charges.  (Divide by 18.2223 to convert to kcals/mol)
c AMASS  : the atom masses
c IAC    : index for the atom types involved in Lennard Jones (6-12)
c          interactions.  See ICO below.
c ICO    : provides the index to the nonbon parameter
c          arrays CN1, CN2 and ASOL, BSOL.  All possible 6-12
c          or 10-12 atoms type interactions are represented.
c          NOTE: A particular atom type can have either a 10-12
c          or a 6-12 interaction, but not both.  The index is
c          calculated as follows:
c            index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
c          If index is positive, this is an index into the
c          6-12 parameter arrays (CN1 and CN2) otherwise it
c          is an index into the 10-12 parameter arrays (ASOL
c          and BSOL).
c LABRES : the residue labels
c IPRES  : atoms in each residue are listed for atom "i" in
c          IPRES(i) to IPRES(i+1)-1
c RK     : force constant for the bonds of each type, kcal/mol
c REQ    : the equilibrium bond length for the bonds of each type, angstroms
c TEQ    : the equilibrium angle for the angles of each type, degrees
c PK     : force constant for the dihedrals of each type, kcal/mol
c PN     : periodicity of the dihedral of a given type
c PHASE  : phase of the dihedral of a given type
c SOLTY  : currently unused (reserved for future use)
c CN1    : Lennard Jones r**12 terms for all possible atom type
c          interactions, indexed by ICO and IAC; for atom i and j
c          where i < j, the index into this array is as follows
c          (assuming in index is positive):
c          CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).
c CN2    : Lennard Jones r**6 terms for all possible atom type
c          interactions.  Indexed like CN1 above.
c_______________________________________________________________________________
 
      IMPLICIT      DOUBLE PRECISION (A-H,O-Z)

      PARAMETER    (MAXNSPM = 1000)
      PARAMETER    (MAXATOM = 8000)
      PARAMETER    (MAXRES  =  700)

      integer newformat

      logical       debug
      CHARACTER*4   IGRAPH, LABRES
      character     TITL*80
      character*40  tmp_str
c     AMBER12 TOP files appear to have an extra block for
c     atomic numbers 
      character*20  title_section

      COMMON /MDINFO/ NTB,NRUN,NPM,NRP,NSM,NRAM,nsns,IG
      COMMON /INFOV/  NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,
     +                NPHIA,NNB,NTYPES,NRC,NCONP,MAXMEM,NWDVAR,MAXNB,
     +                MBONA,MTHETA,MPHIA,NATC,IBELLY,NATBEL,ISHAKE,NMXRS
      COMMON /NBTERM/ CUT,SCNB,SCEE,IDIEL,iddd,DIELC,NBUCK,NUMPK,NBIT
      COMMON /MDSOL/  IFTRES,IPTRES,IPTATM,NSPSTR,IPTSOL,IFCRST,
     +                IMGSLT,nbsola,nbsolh,nsolu,nsolv
      COMMON /WATCAP/ IFCAP,NATCAP,CUTCAP,XCAP,YCAP,ZCAP,FCAP
      COMMON /NBPARM/ IGRAPH(MAXATOM), CHRG(MAXATOM), AMASS(MAXATOM),
     +                IAC(MAXATOM), ICO(1830), LABRES(MAXRES),
     +                IPRES(MAXRES), ISYMBL(MAXATOM), ITREE(MAXATOM),
     +                JOIN(MAXATOM), IROTAT(MAXATOM)
      COMMON /BPARM/  IBH(4500), JBH(4500), ICBH(4500),
     +                IB(4500), JB(4500), ICB(4500),
     +                ITH(9000), JTH(9000), KTH(9000), ICTH(9000),
     +                IT(7500), JT(7500), KT(7500), ICT(7500),
     +                IPH(16500),JPH(16500),KPH(16500),
     +                LPH(16500),ICPH(16500),
     +                IP(8000),JP(8000),KP(8000),LP(8000),ICP(8000),
     +                NATEX(45000),NUMEX(MAXATOM)
      COMMON /PARMS/  RK(500),REQ(500),TK(900),TEQ(900),PK(900),
     +                PN(900),PHASE(900),CN1(1830),CN2(1830),RAD(61),
     +                SOLTY(60),GAMC(900),GAMS(900),IPN(900),FMN(900),
     +                ONE_SCEE(900),ONE_SCNB(900)
      COMMON /HBPAR/  ASOL(200),BSOL(200),HBCUT(200)
      COMMON /TOPA/ NSPSOL, NSPM, NATOM
      COMMON /TOPB/ NSP(MAXNSPM), BOXT(3), BOXC(3)

c-------------------------------------------------------------------------------
C NPHB is the number of h-bond parameters. NIMPRP is the number of
C improper torsional parameters (NPTRA-NIMPRP is the number of regular
C torsional parameters).
c-------------------------------------------------------------------------------

      COMMON /PRMLIM/ NUMBND,NUMANG,NPTRA,NPHB,NIMPRP

c-------------------------------------------------------------------------------
c --- initialize variables ---
c-------------------------------------------------------------------------------
      
      IPTRES = 0
      IPTATM = 0
      NSPSOL = 0
      NSPSTR = 0
      
c-------------------------------------------------------------------------------
C     ----- FORMATTED INPUT -----
c-------------------------------------------------------------------------------
      
      READ(NF,'(a80)') TITL
      
      newformat=0
      if(TITL(1:1).eq.'%')then
         newformat=1
      endif
      
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
         read(NF,*)
         read(NF,*)
         read(NF,*)
      endif
      
      if (newformat.eq.1)then
         READ(NF,'(10i8)') NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,
     +        NPHIH, MPHIA,  NHPARM, NPARM,  NEXT,   NRES,
     +        NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
     +        NATYP, NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
     +        MBPER, MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP,null
      else
         READ(NF,'(12i6)') NATOM, NTYPES, NBONH,  MBONA,  NTHETH, 
     +    MTHETA,NPHIH, MPHIA,  NHPARM, NPARM,  NEXT,   NRES,
     +        NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
     +        NATYP, NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
     +        MBPER, MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP
      endif

c      write(*,'(10i8)') NATOM, NTYPES, NBONH,  MBONA,  NTHETH, MTHETA,
c     +                  NPHIH, MPHIA,  NHPARM, NPARM,  NEXT,   NRES,
c     +                  NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
c     +                  NATYP, NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
c     +                  MBPER, MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP

c-------------------------------------------------------------------------------
c --- make sure we don't exceed memory limits in commons ---
c-------------------------------------------------------------------------------

      NTTYP = NTYPES*(NTYPES+1)/2
 
      if (numbnd.gt.500 .or. numang.gt.900 .or. nptra.gt.900 .or.
     +     nphb.gt.200 .or. natyp.gt.60 .or. nttyp.gt.1830 .or.
     +     nbonh.gt.4500 .or. nbona.gt.4500 .or. ntheth.gt.9000 .or.
     +     ntheta.gt.7500 .or. nphih.gt.16500 .or. nphia.gt.12000 .or.
     +     next.gt.45000 .or. natom.gt.maxatom .or. nres.gt.maxres) then
         write(no,*) 'numbnd (max 500)   = ', numbnd
         write(no,*) 'numang (max 900)   = ', numang
         write(no,*) 'nptra  (max 900)   = ', nptra 
         write(no,*) 'nphb   (max 200)   = ', nphb   
         write(no,*) 'natyp  (max 60)    = ', natyp  
         write(no,*) 'nttyp  (max 1830)  = ', nttyp  
         write(no,*) 'nbonh  (max 4500)  = ', nbonh  
         write(no,*) 'nbona  (max 4500)  = ', nbona  
         write(no,*) 'ntheth (max 9000)  = ', ntheth 
         write(no,*) 'ntheta (max 7500)  = ', ntheta 
         write(no,*) 'nphih  (max 16500) = ', nphih  
         write(no,*) 'nphia  (max 12000)  = ', nphia  
         write(no,*) 'next   (max 45000) = ', next   
         write(no,*) 'natom  (max',maxatom,') =', natom
         write(no,*) 'nres   (max',maxres, ') =', nres

         write(*,*) 'Exceding memory limits at the top file'

         goto 1000
      endif

c-------------------------------------------------------------------------------
C     ----- read charges, masses, atom types, etc... 
c-------------------------------------------------------------------------------
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif 
      READ(NF,'(20a4)') (IGRAPH(i), i=1,NATOM)
      if (newformat.eq.1)then     
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(5E16.8)')  (CHRG(i),   i=1,NATOM)
      if (newformat.eq.1)then
         read(NF,'(a20)')  title_section
         read(NF,*)
      endif

c     ON 19/12/2013 ACORTES WAS HERE
c     THIS IS TO SUPPORT AMBER12 TOP FILES WHILE
c     BEING COMPATIBLE WITH OLD VERSIONS
c     JUST DROP THE NEXT BLOCK IF IT IS NOT THE MASS
      WRITE(*,*) title_section
      if( title_section.EQ.'%FLAG ATOMIC_NUMBER') then
c        WRITE(*,*) "ATOMIC NUMBER BLOCK!"
        READ(NF,'(10I8)')  (null,    i=1,NATOM)
         read(NF,*)
         read(NF,*)
      endif

      READ(NF,'(5E16.8)') (AMASS(i),  i=1,NATOM)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
         READ(NF,'(10I8)')   (IAC(i),    i=1,NATOM)
      else
         READ(NF,'(12I6)')   (IAC(i),    i=1,NATOM)
      endif
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
         READ(NF,'(10I8)')   (NUMEX(i),  i=1,NATOM)
      else
         READ(NF,'(12I6)')   (NUMEX(i),  i=1,NATOM)
      endif
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
         READ(NF,'(10I8)')   (ICO(i),    i=1,NTYPES*NTYPES)
      else
         READ(NF,'(12I6)')   (ICO(i),    i=1,NTYPES*NTYPES)
      endif
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(20A4)')   (LABRES(i), i=1,NRES)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
         READ(NF,'(10I8)')   (IPRES(i),  i=1,NRES)
      else
         READ(NF,'(12I6)')   (IPRES(i),  i=1,NRES)
      endif
      
c-------------------------------------------------------------------------------
C     ----- READ THE PARAMETERS -----
c-------------------------------------------------------------------------------
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif      
      READ(NF,'(5E16.8)') (RK(i),     i=1,NUMBND)
      if (newformat.eq.1)then
         read(NF,*) 
         read(NF,*)
      endif
      READ(NF,'(5E16.8)') (REQ(I),    i=1,NUMBND)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(5E16.8)') (TK(I),     i=1,NUMANG)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(5E16.8)') (TEQ(I),    i=1,NUMANG)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(5E16.8)') (PK(I),     I = 1,NPTRA)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(5E16.8)') (PN(I),     I = 1,NPTRA)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(5E16.8)') (PHASE(I),  I = 1,NPTRA)
      if (newformat.eq.1)then
         read(NF,'(a)') tmp_str
         read(NF,*)
      endif

c Modification for reading AMBER 11 top files
      if(tmp_str(7:23) .eq. 'SCEE_SCALE_FACTOR') then
         READ(NF,'(5E16.8)') (ONE_SCEE(I),  I = 1,NPTRA)
         if (newformat.eq.1)then
            read(NF,*)
            read(NF,*)
         endif
         READ(NF,'(5E16.8)') (ONE_SCNB(I),  I = 1,NPTRA)
         if (newformat.eq.1)then
            read(NF,*)
            read(NF,*)
         endif
      endif

      READ(NF,'(5E16.8)') (SOLTY(I),  I = 1,NATYP)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(5E16.8)') (CN1(I),    I = 1,NTTYP)
      if (newformat.eq.1)then
         read(NF,*)
         read(NF,*)
      endif
      READ(NF,'(5E16.8)') (CN2(I),    I = 1,NTTYP)

c-------------------------------------------------------------------------------
c     ----- READ THE BONDING INFORMATION -----
c
c NOTE: the atom numbers in the arrays which follow that describe bonds, angles,
c and dihedrals are obfuscated by the following formula (for runtime speed in
c indexing arrays): 
c The true atom number equals the absolute value of the number divided by
c three, plus one. In the case of the dihedrals, if the third atom is negative, 
c this implies an improper torsion and if the fourth atom is negative, this 
c this implies that end group interactions are to be ignored. 
c End group interactions are ignored, for example, in dihedrals of various 
c ring systems (to prevent double counting) and in multiterm dihedrals.
c-------------------------------------------------------------------------------
      
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)') (IBH(i),JBH(i),ICBH(i), i=1,NBONH)         
c$$$      else
c$$$         READ(NF,'(12I6)') (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
c$$$      endif
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)') (IB(i),JB(i),ICB(i),    i=1,NBONA)
c$$$      else
c$$$         READ(NF,'(12I6)') (IB(i),JB(i),ICB(i),    i=1,NBONA)
c$$$      endif
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)') (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)         
c$$$      else
c$$$         READ(NF,'(12I6)') (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)
c$$$      endif
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)') (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)
c$$$      else
c$$$         READ(NF,'(12I6)') (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)
c$$$      endif
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)') (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i),
c$$$     &        i=1,NPHIH)
c$$$      else
c$$$         READ(NF,'(12I6)') (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i),
c$$$     &        i=1,NPHIH)
c$$$      endif
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)') (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)
c$$$      else
c$$$         READ(NF,'(12I6)') (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)
c$$$      endif
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)') (NATEX(i),  i=1,NEXT)
c$$$      else
c$$$         READ(NF,'(12I6)') (NATEX(i),  i=1,NEXT)
c$$$      endif
c$$$
c$$$c-------------------------------------------------------------------------------
c$$$c     ----- READ THE H-BOND PARAMETERS -----
c$$$c-------------------------------------------------------------------------------
c$$$
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$      endif
c$$$      READ(NF,'(5E16.8)') (ASOL(i),   i=1,NPHB)
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$      endif
c$$$      READ(NF,'(5E16.8)') (BSOL(i),   i=1,NPHB)
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$      endif
c$$$      READ(NF,'(5E16.8)') (HBCUT(i),  i=1,NPHB)
c$$$
c$$$c-------------------------------------------------------------------------------
c$$$C     ----- READ ISYMBL,ITREE,JOIN AND IROTAT ARRAYS -----
c$$$c-------------------------------------------------------------------------------
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$      endif
c$$$      READ(NF,'(20A4)')   (ISYMBL(i), i=1,NATOM)
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$      endif
c$$$      READ(NF,'(20A4)')   (ITREE(i),  i=1,NATOM)
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)')   (JOIN(i),   i=1,NATOM) 
c$$$      else
c$$$         READ(NF,'(12I6)')   (JOIN(i),   i=1,NATOM)
c$$$      endif
c$$$      if (newformat.eq.1)then
c$$$         read(NF,*)
c$$$         read(NF,*)
c$$$         READ(NF,'(10I8)')   (IROTAT(i), i=1,NATOM)
c$$$      else
c$$$         READ(NF,'(12I6)')   (IROTAT(i), i=1,NATOM)
c$$$      endif

c-------------------------------------------------------------------------------
c     ----- READ THE BOUNDARY CONDITION STUFF -----
c-------------------------------------------------------------------------------

      IF(IFBOX.GT.0) THEN
         if (newformat.eq.1)then
            read(NF,*)
            read(NF,*)
            READ(NF,'(10I8)')   IPTRES,NSPM,NSPSOL
         else
            READ(NF,'(12I6)')   IPTRES,NSPM,NSPSOL
         endif
         if (newformat.eq.1)then
            read(NF,*)
            read(NF,*)
            READ(NF,'(10I8)')  (NSP(i), i=1,NSPM)
         else
            READ(NF,'(12I6)')  (NSP(i), i=1,NSPM)
         endif
         if (newformat.eq.1)then
            read(NF,*)
            read(NF,*)
         endif
         READ(NF,'(5E16.8)') BETA, BOXT(1), BOXT(2), BOXT(3)
      ENDIF

c-------------------------------------------------------------------------------
c --- Write the topology file to output unit 16 ---
c-------------------------------------------------------------------------------

      write(no,'(/,a,a80,a,/)') ' MOLECULE: ', TITL

      write(no,8118) NATOM, NTYPES, NBONH,  MBONA,  NTHETH, MTHETA,
     +               NPHIH, MPHIA,  NHPARM, NPARM,  NEXT,   NRES,
     +               NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
     +               NATYP, NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
     +               MBPER, MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP,
     +               NTTYP
    
c-------------------------------------------------------------------------------
c     CALCULATE INVERSE OF THE MASS, CHARGES, RADII AND WRITE ATOM INFO
c
c     For extracting the radii, the following is used:
c
c     cn1(i,j)=eij*(rij**12)  |
c                             |->  rij**6=2cn1(i,j)/cn2(i,j)
c     cn2(i,j)=2*eij*(rij**6) |
c
c     This, if i=j then: ri= rii/2 = {(2cn1(i,i)/cn2(i,i))**1/6}/2
c
c-------------------------------------------------------------------------------

      do k = 1, ntypes
        rad(k) = 0.0d0
        index  = ico(ntypes*(k-1)+k)
        if (index.gt.0) then
          if (cn2(index).ne.0.0d0) then
            rad(k) = (2.0d0*cn1(index))/cn2(index)
            rad(k) = (rad(k)**(1.0d0/6.0d0))/2.0d0
          endif
        endif
      enddo

      DO I = 1, NATOM
        CHRG(I) = CHRG(i)/18.2223d0
c       write(no,8119) I,IGRAPH(i),CHRG(i),AMASS(I),rad(iac(i)),IAC(i)
        AMASS(I) = 1.0E0/AMASS(I)
      ENDDO

c-------------------------------------------------------------------------------
c     write rest of topology parameters if requested to do so
c-------------------------------------------------------------------------------

      debug = .false.
      if (debug) then
        write(no,'(20a4)')   (IGRAPH(i), i=1,NATOM)
        write(no,'(5E16.8)') (CHRG(i),   i=1,NATOM)
        write(no,'(5E16.8)') (AMASS(i),  i=1,NATOM)
        write(no,'(12I6)')   (IAC(i),    i=1,NATOM)
        write(no,'(12I6)')   (NUMEX(i),  i=1,NATOM)
        write(no,'(12I6)')   (ICO(i),    i=1,NTYPES*NTYPES)
        write(no,'(20A4)')   (LABRES(i), i=1,NRES)
        write(no,'(5E16.8)') (RK(i),     i=1,NUMBND)
        write(no,'(5E16.8)') (REQ(I),    i=1,NUMBND)
        write(no,'(5E16.8)') (TK(I),     i=1,NUMANG)
        write(no,'(5E16.8)') (TEQ(I),    i=1,NUMANG)
        write(no,'(5E16.8)') (PK(I),     I = 1,NPTRA)
        write(no,'(5E16.8)') (PN(I),     I = 1,NPTRA)
        write(no,'(5E16.8)') (PHASE(I),  I = 1,NPTRA)
        write(no,'(5E16.8)') (SOLTY(I),  I = 1,NATYP)
        write(no,'(5E16.8)') (CN1(I),    I = 1,NTTYP)
        write(no,'(5E16.8)') (CN2(I),    I = 1,NTTYP)
        write(no,'(12I6)') (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
        write(no,'(12I6)') (IB(i),JB(i),ICB(i),    i=1,NBONA)
        write(no,'(12I6)') (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)
        write(no,'(12I6)') (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)
        write(no,'(12I6)') (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i),
     +                     i=1,NPHIH)
        write(no,'(12I6)') (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)
        write(no,'(12I6)') (NATEX(i),  i=1,NEXT)
        write(no,'(5E16.8)') (ASOL(i),   i=1,NPHB)
        write(no,'(5E16.8)') (BSOL(i),   i=1,NPHB)
        write(no,'(5E16.8)') (HBCUT(i),  i=1,NPHB)
        write(no,'(20A4)')   (ISYMBL(i), i=1,NATOM)
        write(no,'(20A4)')   (ITREE(i),  i=1,NATOM)
        write(no,'(12I6)')   (JOIN(i),   i=1,NATOM)
        write(no,'(12I6)')   (IROTAT(i), i=1,NATOM)
        IF(IFBOX.GT.0) THEN
          write(no,'(12I6)')   IPTRES,NSPM,NSPSOL
          write(no,'(12I6)')  (NSP(i), i=1,NSPM)
          write(no,'(5E16.8)') BETA, BOXT(1), BOXT(2), BOXT(3)
        ENDIF
      endif

c-------------------------------------------------------------------------------
c     ----- SCALE THE CHARGES IF DIELC.GT.1.0E0 -----
c-------------------------------------------------------------------------------

c     IF (DIELC .GT. 1.0E0) THEN
c       DUMD = SQRT(DIELC)
c       DO I = 1,NATOM
c         CHRG(I) = CHRG(I)/DUMD
c       ENDDO
c     ENDIF

c-------------------------------------------------------------------------------
c --- Formats, error messages and end of the subroutine ---
c-------------------------------------------------------------------------------

 8118 format(t2,
     +  'NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7,
     +/' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7,
     +/' NHPARM = ',i7,' NPARM  = ',i7,' NEXT  = ',i7,' NRES   = ',i7,

     +/' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7,
     +/' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7,
     +/' IFPERT = ',i7,' NBPER  = ',i7,' NGPER = ',i7,' NGPER  = ',i7,
     +/' MBPER  = ',i7,' MGPER  = ',i7,' MDPER = ',i7,' IFBOX  = ',i7,
     +/' NMXRS  = ',i7,' IFCAP  = ',i7,' NTTYP = ',i7/)

c8119 format(1x,i4,a4,3e16.8,i6)

      goto 99
 1000 continue
      write(no,'(/,5x,a)') 'readtop: a parameter array overflowed'
      stop
 99   RETURN
      END
c
c===============================================================================
c
      integer function lchar(string)
c
c-------------------------------------------------------------------------
c --- Returns the position of the last printable character (not a blank)
c --- in the string STRING
c-------------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)

      character*(*)  string
      integer j
C
      lchar=len(string)
      do j=len(string),1,-1
        if (string(j:j).ne.' ') then
          lchar=j
          return
        end if
        lchar=j
      end do
C
      return
      end

c
c===============================================================================
c
   
      subroutine qsar(nlv,idel,nsc,nrank,nrow,ncol,ntest,xraw,y,cont,
     &     cutinp,randtest,ncross,name,comp,cf,varselect,nvtold,
     &     vardesc,b)

c-------------------------------------------------------------------------------
c     ##  subroutine qsar  --  compute multiple regression QSAR model  ##
c
c     "qsar" gets a column of data points to be modeled and a matrix
c     of property values to be used to model the data, and returns
c     a linear least squares QSAR-style model; the program has two
c     modes of operation, it can test a model using only specified
c     property columns, or use simulated annealing to optimize the
c     predictive value of the model
c-------------------------------------------------------------------------------
 
      implicit none

      integer maxrow,maxcol,nvt,nvar,nad,maxscr,maxrank,maxcompl
      parameter (maxrow= 100)
      parameter (maxcol=1500)
      parameter (maxscr= 100)
      parameter (maxrank=20)
      parameter (maxcompl = 100)

      integer i,j,k
      integer nrank,kopt
      integer nrow,ncol,ntest
      integer nlv,idel,nsc,nscram,cont
      integer varselect(maxcol),nvtold
      integer randtest,ncross

      real*8 ybar,ystd,sigma,z,dm,dt
      real*8 b(maxcol,maxcol),bs(maxcol,maxcol)
      real*8 xref(maxrow,maxcol),xbar(maxcol),w(maxcol)
      real*8 xraw(maxrow,maxcol),cf(maxcol)
      real*8 y(maxrow),yr(maxrow,maxscr),ys(maxrow)
      real*8 q2r(maxrank),q2s(maxrank),q2sav(maxrank),q2sqr(maxrank)

      real*8 cutinp,xstd(maxcol)

      character*4   name(maxrow)
      character*20  vardesc(maxcol)
      character*40  comp(maxcompl)

      logical verbose,verbosetmp

      common /variables/ nvt,nvar,nad

c-------------------------------------------------------------------------------
c     autoscale the response vector and property matrix values
c-------------------------------------------------------------------------------

      call yscale (nsc,nrow,ntest,y,ybar,ystd)
      do i = 1, ncol
         do j = 1, nrow+ntest
            xref(j,i) = xraw(j,i)
         end do       
      end do

      call xscale(nsc,nrow,ncol,ntest,xref,xbar,xstd,w)

c-------------------------------------------------------------------------------
c     transform the properties via partial least squares
c-------------------------------------------------------------------------------

      verbose=.true.
      call nipals (nlv,nrow,ncol,y,xref,nrank,b,vardesc,varselect,nvtold
     &     ,verbose)
      verbosetmp=.true.
   
      call plsmodel (nlv,nsc,nrow,ncol,ntest,comp,y,ybar,q2r,
     &     randtest,ncross,ystd,xref,xbar,w,b,name,verbose,verbosetmp)

c-------------------------------------------------------------------------------
c     Scrambling of activities
c-------------------------------------------------------------------------------

      if (idel .ne. 0) then
        nscram = 100
        call scramble(nscram,nrow+ntest,y,yr)
         
        do k = 1, nlv
            q2sav(k) = 0.0d0
            q2sqr(k) = 0.0d0
        enddo
        do i = 1, nscram
            write(16,'(/,a,i4)') ' SCRAMBLE TRIAL: ',i
            do j = 1, nrow+ntest
               ys(j) = yr(j,i)
            enddo
            verbose=.false.
            call nipals (nlv,nrow,ncol,ys,xref,nrank,bs,vardesc,
     &           varselect,nvtold,verbose)
            verbosetmp=.false.
            call plsmodel (nlv,nsc,nrow,ncol,ntest,comp,ys,ybar,q2s,
     &           randtest,ncross,ystd,xref,xbar,w,bs,name,verbose
     &           ,verbosetmp)
            do k = 1, nlv
               q2sav(k) = q2sav(k) + q2s(k)
               q2sqr(k) = q2sqr(k) + q2s(k)*q2s(k)
           enddo
        enddo
        
        open(12,file='scramble.out')

c$$$        write(12,*)
c$$$        write(12,'(10x,a,a)') '======================================',
c$$$     &        '============'
c$$$        write(12,'(10x,a,a)') 'LV     Q2-real     z-score    <Q2-scr>',
c$$$     &        '     st.dev.'
c$$$        write(12,'(10x,a,a)') '======================================',
c$$$     &        '============'
c$$$        write(12,*)
         
        dm=0.0d0
        kopt=0
        do k = 1, nlv
            q2sav(k) = q2sav(k)/dble(nscram)
            q2sqr(k) = q2sqr(k)/dble(nscram)
            sigma = dsqrt(q2sqr(k)-q2sav(k))
            z = (q2r(k)-q2sav(k))/sigma
            write(12,'(i2,4x,4(f8.4,4x))') k,q2r(k),z,q2sav(k),sigma
            dt=dsqrt(q2r(k)*q2r(k)+z*z)
            if (dt.ge.dm) then
               dm=dt
               kopt=k
            endif
        enddo
        write(12,'(/,10x,a,i4,/)') 'SELECTED MODEL DIMENSIONALITY: '
     &       ,kopt
        do j = 1, ncol
           cf(j) = b(j,kopt)
        enddo
        write(12,*)
      endif
      
      close(12)

c-------------------------------------------------------------------------------
c     End of the subroutine
c-------------------------------------------------------------------------------
      
      return
      end
c
c===============================================================================
c
      subroutine yscale (nsc,nrow,ntest,y,ybar,ystd)
 
c------------------------------------------------------------------------------
c                                                            
c     subroutine yscale  --  autoscale the response values  
c     Warning: Scaling of reponse variables only with nsc=1 !!
c                                                             
c------------------------------------------------------------------------------
 
      implicit none
      integer maxrow
      parameter (maxrow=100)
      integer i,nrow,ntest,nsc
      real*8 ybar,ystd
      real*8 y(maxrow),yraw(maxrow)
c
c     center and scale the response data values
c
      ybar = 0.0d0
      do i = 1, nrow
         yraw(i) = y(i)
         ybar = ybar + y(i)
      end do
      ybar = ybar / dble(nrow)
      if (nsc .eq. 1) then
        ystd = 0.0d0
        do i = 1, nrow
           ystd = ystd + (y(i)-ybar)**2
        end do
        ystd = sqrt(ystd/dble(nrow-1))
        if (ystd .eq. 0.0d0)  ystd = 1.0d0
      else
        ystd = 1.0d0
      endif
      do i = 1, nrow+ntest
         y(i) = (y(i)-ybar) / ystd
      end do
      return
      end
c
c===============================================================================
c
      subroutine xscale(nsc,nrow,ncol,ntest,x,xbar,xstd2,w)

c------------------------------------------------------------------------------
c     This subroutine carries out the scaling of the X-matrix. Several
c     scaling options are implemented:
c
c     nsc = 0 => no scaling 
c     nsc = 1 => autoscaling
c     nsc = 2 => block-unscaled weights
c
c     Warning: it is assumed that x-variables follow the following blocking
c     scheme: series of K blocks of NVAR variables followed by NAD variables
c
c     CHANGES: Returns xstd2;
c
c                                 =ARO-00=
c------------------------------------------------------------------------------

      implicit none

      integer maxrow,maxcol,nvt,nvar,nad
      parameter (maxrow=100)
      parameter (maxcol=1500)
      integer i,j,k,iblnu,iblock(20,2)
      integer nrow,ncol,ntest,nsc
      real*8 x(maxrow,maxcol),w(maxcol)
      real*8 xbar(maxcol),xstd(maxcol),sspar(20),sstot,stcut
      parameter (stcut=1E-4)
      
      common /variables/ nvt,nvar,nad

      real*8 xstd2(maxcol)

c------------------------------------------------------------------------------
c     compute and scale property values. First compute average per variable,
c     then obtain weights and finally scale property values
c------------------------------------------------------------------------------
c
c compute average and standard dev of each variable
c
      do j = 1, ncol
         xbar(j) = 0.0d0
         do i = 1, nrow
            xbar(j) = xbar(j) + x(i,j)
         end do
         xbar(j) = xbar(j) / dble(nrow)
      enddo

      do j = 1, ncol
         xstd2(j) = 0.0d0
         do i = 1, nrow
            xstd2(j) = xstd2(j) + (x(i,j)-xbar(j))**2
         end do
         xstd2(j) = sqrt(xstd2(j)/dble(nrow-1))
      enddo

c     
c compute weights according to the selected method
c
      if (nsc .eq. 0) then

         do j = 1, ncol
            xstd(j) = 1.0d0
            w(j)    = 1.0d0
         enddo
         
      else if (nsc .eq. 1) then
         
         do j = 1, ncol
            xstd(j) = 0.0d0
            do i = 1, nrow
               xstd(j) = xstd(j) + (x(i,j)-xbar(j))**2
            end do
            xstd(j) = sqrt(xstd(j)/dble(nrow-1))
            if (xstd(j) .eq. 0.0d0) then
               xstd(j) = xbar(j)
               xbar(j) = 0.0d0 
            end if
            if (xstd(j).gt.stcut) then
               w(j) = 1.0/xstd(j)
            else
               w(j) = 0.0d0
            endif
         enddo
         
      else if (nsc .eq. 2) then
  
         if ( nvt/nvar.ge.2 )then
            if (nad .eq. 0) then
               iblnu = nvt/nvar
               do k = 1, iblnu
                  iblock(k,2) = k*nvar
                  if (k.eq.1) then
                     iblock(k,1) = 1
                  else
                     iblock(k,1) = iblock(k-1,2)+1
                  endif
               enddo
            else          
               iblnu = (nvt/nvar) + 1
               do k = 1, iblnu-1
                  iblock(k,2) = k*nvar
                  if (k.eq.1) then
                     iblock(k,1) = 1
                  else
                     iblock(k,1) = iblock(k-1,2)+1
                  endif
               enddo      
               iblock(iblnu,1) = iblock(iblnu-1,2)+1
               iblock(iblnu,2) = nvt
            endif
         else
            if (nad .eq. 0) then
               do k = 1, 2
                  iblock(k,2) = k*(nvt/2)
                  if (k.eq.1) then
                     iblock(k,1) = 1
                  else
                     iblock(k,1) = iblock(k-1,2)+1
                  endif
               enddo
            else          
               do k = 1, 2
                  iblock(k,2) = k*(nvt-nad)/2
                  if (k.eq.1) then
                     iblock(k,1) = 1
                  else
                     iblock(k,1) = iblock(k-1,2)+1
                  endif
               enddo      
               iblock(3,1) = iblock(2,2)+1
               iblock(3,2) = nvt
            endif
            
            iblnu=3

         endif
         sstot = 0.0d0
         do k = 1, iblnu
            sspar(k) = 0.0d0
            do j = iblock(k,1), iblock(k,2)
               xstd(j) = 0.0d0
               do i = 1, nrow
                  xstd(j) = xstd(j) + (x(i,j)-xbar(j))**2
               enddo
               sspar(k) = sspar(k)+xstd(j)
            enddo
            sstot = sstot+sspar(k) 
         enddo
         do k = 1, iblnu
            do j = iblock(k,1), iblock(k,2)
               if ((sspar(k).gt.stcut).and.(sstot.gt.stcut)) then
                  w(j) = sqrt(sstot/(dble(iblnu)*sspar(k)))
               else
                  w(j) = 0.0d0
               endif
            enddo
         enddo
         
      endif
c     
c     apply weights
c     
      do j = 1, ncol
         do i = 1, nrow+ntest
            if(xbar(j).eq.0.0d0.and.w(j).eq.0.0d0)then
               x(i,j)=0.0d0
            else
               x(i,j) = (x(i,j)-xbar(j)) * w(j)
            endif
         end do
      end do
 
c-------------------------------------------------------------------------------
c  End of subroutine
c-------------------------------------------------------------------------------

      return
      end
c
c===============================================================================
c
      subroutine nipals (nlv,nrow,ncol,y,x,nrank,b,vardesc,varselect,
     &     nvtold,verbose)
c_______________________________________________________________________________
c
c     ##  subroutine nipals  --  iterative partial least squares  ##
c
c     "nipals" finds partial least squares regression coefficients
c     via the classical nonlinear iterative PLS algorithm devised
c     by Herman Wold
c
c     literature reference:
c
c     S. Rannar, F. Lindgren, P. Geladi and S. Wold, "A PLS Kernel
c     Algorithm for Data Sets with many Variables and fewer Objects.
c     Part 1: Theory and Algorithms", Journal of Chemometrics, 8,
c     111-125 (1994)  [see pages 113-114 for nipals implementation]
c
c
c                              PLS ALGORITHM
c
c
c   1)  Start: set u to the first column of Y
c   2)  w = X'u/(u'u)
c   3)  Scale w to be of length one
c   4)  t = Xw
c   5)  c = Y't/(t't)
c   6)  Scale c to be of length one
c   7)  u = Y'c/(c'c)
c   8)  If convergence then 9 else 2
c   9)  X-loadings: p = X't/(t't)
c  10)  Y-loadings: q = Y'u/(u'u)
c  11)  Regression (u upon t): LinRel = u't/(t't)
c  12)  Residual matrices: X -> X-tp' and Y -> Y-btc'
c
c
c                 ( N x Mx )                   ( N x My )
c
c                     X           t         u      Y
c            |-----------------|  |         |  |-------|
c            |                 |  |         |  |       |
c            |                 |  |         |  |       |
c            |                 |  |         |  |       |
c            |                 |  |         |  |       |
c            |-----------------|  |         |  |-------|
c
c        w'  -------------------           c'  ---------
c
c        p'  -------------------           q'  ---------
c
c
c            p,q --> LOADINGS
c            t,u --> SCORES
c
c
c       The next set of iterations starts with the new X and Y matrices
c       as the residual matrices from the previous iteration. The iterations
c       can continue until a stopping criteria is used or X becomes zero.
c
c       Finally, we have X = TP' + E, Y = UC' + F, and U = bT,
c       where E and F are the residual matricies, b is the inner relationship
c       vector, and X, T, P', E, Y, U, C', and F are matricies:
c
c       X = row x column
c           (for example, row = compound and column = independent variable)
c       T = score x factor (or component)
c       P'= factor x loading
c       E = row x column
c           (actual data minus linear combination of scores and loadings)
c       Y = row x column
c           (for example, row = compound and column = dependent variable)
c       U = score x factor (or component)
c       C'= factor x loading
c       F = row x column (actual data minus linear combination
c            of scores and loadings)
c
c
c          n  = Number of molecules
c          mx = Number of variables in the X matrix
c          my = Number of variable in the Y matrix
c               (usually 1, the activity).
c
c                                =ARO-97=
c_______________________________________________________________________________
 
      implicit none
      integer maxrow,maxcol,maxrank
      parameter (maxrow=100)
      parameter (maxcol=1500)
      parameter (maxrank=20)

      integer i,j,k,m
      integer nrow,ncol,nrank,nlv
      integer varselect(maxcol),nvtold

      real*8 norm
      real*8 y(maxrow),x(maxrow,maxcol)
      real*8 b(maxcol,maxcol)
      real*8 yref(maxrow),xref(maxrow,maxcol)
      real*8 u(maxrow),t(maxrow,maxrow),c(maxrank)
      real*8 w(maxcol,maxrank),p(maxcol,maxrank)
      real*8 pw(maxrank,maxrank),wpw(maxcol,maxrank)

      character*20  vardesc(maxcol)
      logical verbose
      integer nupw,nupl,nurw,nurc
      data nupw,nupl,nurw,nurc /50,51,52,53/
 
c-------------------------------------------------------------------------------
c     set the rank and store the response and property values
c-------------------------------------------------------------------------------

      nrank = min(ncol,nrow-1,maxrank,nlv)
      do i = 1, nrow
         yref(i) = y(i)
      end do
      do j = 1, ncol
         do i = 1, nrow
            xref(i,j) = x(i,j)
         end do
      end do

c-------------------------------------------------------------------------------
c     perform the nipals-PLS1 algorithm for successive ranks
c
c     NOTE: This algorithm is different form the "standard" NIPALS
c     algorithm. Here, for each rank it is assumed convergence in
c     two iterations in the extraction of each component. Other
c     algorithms use more stringent criteria for convergence, such as
c     setting the number of iterations to 50 and checking for "exhaustive"
c     extraction. Difference with the more "traditional" algorithm is small,
c     while allowing very fast CPU times. But more testing is required.
c-------------------------------------------------------------------------------

      do k = 1, nrank
         do i = 1, nrow
            u(i) = y(i)
         end do
         norm = 0.0d0
         do j = 1, ncol
            w(j,k) = 0.0d0
            do i = 1, nrow
               w(j,k) = w(j,k) + u(i)*x(i,j)
            end do
            norm = norm + w(j,k)*w(j,k)
         end do
         norm = sqrt(norm)
         if (norm .eq. 0.0d0) then
            nrank = k - 1
            goto 10
         end if
         do j = 1, ncol
            w(j,k) = w(j,k) / norm
         end do
         norm = 0.0d0
         do i = 1, nrow
            t(i,k) = 0.0d0
            do j = 1, ncol
               t(i,k) = t(i,k) + x(i,j)*w(j,k)
            end do
            norm = norm + t(i,k)*t(i,k)
         end do
         c(k) = 0.0d0
         do i = 1, nrow
            c(k) = c(k) + t(i,k)*y(i)
         end do
         c(k) = c(k) / norm
         do i = 1, nrow
            u(i) = y(i) / c(k)
         end do
         norm = 0.0d0
         do j = 1, ncol
            w(j,k) = 0.0d0
            do i = 1, nrow
               w(j,k) = w(j,k) + u(i)*x(i,j)
            end do
            norm = norm + w(j,k)*w(j,k)
         end do
         norm = sqrt(norm)
         do j = 1, ncol
            w(j,k) = w(j,k) / norm
         end do
         norm = 0.0d0
         do i = 1, nrow
            t(i,k) = 0.0d0
            do j = 1, ncol
               t(i,k) = t(i,k) + x(i,j)*w(j,k)
            end do
            norm = norm + t(i,k)*t(i,k)
         end do
         c(k) = 0.0d0
         do i = 1, nrow
            c(k) = c(k) + t(i,k)*y(i)
         end do
         c(k) = c(k) / norm
         do j = 1, ncol
            p(j,k) = 0.0d0
            do i = 1, nrow
               p(j,k) = p(j,k) + t(i,k)*x(i,j)
            end do
            p(j,k) = p(j,k) / norm
         end do
         do j = 1, ncol
            do i = 1, nrow
               x(i,j) = x(i,j) - t(i,k)*p(j,k)
            end do
         end do
         do i = 1, nrow
            y(i) = y(i) - c(k)*t(i,k)
         end do
      end do

c-------------------------------------------------------------------------------
c     restore the original response and property values
c-------------------------------------------------------------------------------

   10 continue
      do i = 1, nrow
         y(i) = yref(i)
      end do
      do j = 1, ncol
         do i = 1, nrow
            x(i,j) = xref(i,j)
         end do
      end do

c-------------------------------------------------------------------------------
c     print the partial least squares weights and loadings
c-------------------------------------------------------------------------------

      if (verbose) then
         open(unit=nupw,file='combine.prw')
         write (nupw,20)
   20    format (/,'# PLS Property Weights :',/)
         do i = 1, ncol
            write (nupw,30)  (w(i,j),j=1,nrank)
   30       format (20(1x,f7.3))
c  30       format (<nrank>f7.3)
c     to be compliant with Absoft 4.5
         end do
         close(nupw)
         open(unit=nupl,file='combine.ldg')
         write (nupl,40)
   40    format (/,'# PLS Property Loadings :',/)

         i=0
         do k = 1, nvtold
            if(varselect(k).eq.1)then
               i=i+1
               write (nupl,50) vardesc(k),(p(i,j),j=1,nrank)
            endif
         end do
 50      format (a12,20(1x,1f9.3))
c     50       format (<nrank>f7.3)
c     to be compliant with Absoft 4.5

         close(nupl)
         open(unit=nurw,file='combine.rpw')
         write (nurw,60)
   60    format (/,'# PLS Response Weights :',/)
         write (nurw,70)  (c(i),i=1,nrank)
   70    format (20(1x,1f7.3))
c  70    format (<nrank>f7.3)
c     to be compliant with Absoft 4.5
         close(nurw)
      end if

c-------------------------------------------------------------------------------
c     get the pls regression coefficients for increasing rank
c-------------------------------------------------------------------------------

      do m = 1, nrank
         do i = 1, m
            do j = 1, m
               pw(j,i) = 0.0d0
               do k = 1, ncol
                  pw(j,i) = pw(j,i) + p(k,j)*w(k,i)
               end do
            end do
         end do         
         call invert (m,maxrank,pw)
         do i = 1, ncol
            do j = 1, m
               wpw(i,j) = 0.0d0
               do k = 1, m
                  wpw(i,j) = wpw(i,j) + w(i,k)*pw(k,j)
               end do
            end do
         end do
         do i = 1, ncol
            b(i,m) = 0.0d0
            do j = 1, m
               b(i,m) = b(i,m) + wpw(i,j)*c(j)
            end do
         end do
      end do

 73   format (20(1x,f9.3)) 
c-------------------------------------------------------------------------------
c     Calculation of Hotelling's confidence elipse
c-------------------------------------------------------------------------------

      if(verbose) then
         call elipse(t,nrow,nrank)
      endif

c-------------------------------------------------------------------------------
c     print partial least squares regression coefficients
c-------------------------------------------------------------------------------

      if (verbose) then
         open(unit=nurc,file='combine.cff')
         write (nurc,80)
   80    format (/,'# PLS Regression Coefficients by Rank :',/)
         j=0
         do k = 1, nvtold
            if(varselect(k).eq.1)then
               j=j+1
               write (nurc,90)  vardesc(k),(b(j,i),i=1,nrank)
            else
               write (nurc,90)  vardesc(k),(0.0,i=1,nrank)
            endif
         end do
 90            format (a12,1x,20(1x,1f9.3))
c     90       format (i8,3x,<nrank>f7.3)
c     To be complaint with Absoft 4.5
         close(nurc)
      end if
     
      if (verbose) then
         open(unit=44,file='combine.scr')
         write (44,91)
 91      format (/,'# PLS T Score coefficients by Rank :',/)
         do j = 1, nrow
            write (44,110)  j,(t(j,i),i=1,nrank)
 110        format (i8,3x,20(1x,1f7.3))
c     90       format (i8,3x,<nrank>f7.3)
c     To be complaint with Absoft 4.5
         end do
         close(44)
      end if     

      return
      end
c
c===============================================================================
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine plsmodel  --  get predictions from PLS model  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine plsmodel (nlv,nsc,nrow,ncol,ntest,comp,y,ybar,q2,
     &     randtest,ncross,ystd,x,xbar,w,b,name,verbose1,verbosetmp)
      implicit none
      integer maxrow,maxcol,maxrank,nvt,nvar,nad,maxcompl
      parameter (maxrow=100)
      parameter (maxcol=1500)
      parameter (maxrank=20)
      parameter (maxcompl = 100)

      integer i,j,k,m,extrap
      integer nrank,nlv,nsc
      integer nrow,ncol,ntest,npred
      integer nvtold
      integer randtest,ncross
      integer ftrim,btrim

      real*8 rsq,rsquared
      real*8 qsq,qsquared
      real*8 rave,qave,abserror
      real*8 rrms,qrms,rmserror
      real*8 qsq2,qave2,qrms2
      real*8 ybar,ystd
      real*8 xbar(maxcol),w(maxcol)
      real*8 xmax(maxcol),xmin(maxcol)
      real*8 y(maxrow),x(maxrow,maxcol)
      real*8 yraw(maxrow),xraw(maxrow,maxcol)
      real*8 yfit(maxrow),ypred(maxrow)
      real*8 b(maxcol,maxcol)
      real*8 q2(maxrank)
      character*4 name(maxrow)
      character*13 combine

      logical verbose1,verbose2,verbosetmp
      character*10  chartemp
      character*40  comp(maxcompl)
      
      common /variables/ nvt,nvar,nad

c-------------------------------------------------------------------------------
c     regenerate the original unscaled response and properties
c-------------------------------------------------------------------------------

      do i = 1, nrow+ntest
         yraw(i) = y(i)*ystd + ybar
      end do
      do i = 1, ncol
         do j = 1, nrow+ntest
            xraw(j,i) = (x(j,i)/w(i)) + xbar(i)
            if (w(i).eq.0) xraw(j,i)=0.0d0 !!! NAN's=0 ,from variables with dev=0
         end do
      end do

      
c-------------------------------------------------------------------------------
c     find the minimum and maximum value for each property
c-------------------------------------------------------------------------------

      if (ntest .eq. 0) then
         write(6,*) 'No test set provided'
      else
         do i = 1, ncol
            xmax(i) = -1.0d6
            xmin(i) = 1.0d6
            do j = 1, nrow
               xmax(i) = max(xmax(i),x(j,i))
               xmin(i) = min(xmin(i),x(j,i))
            end do
         end do
      end if

c-------------------------------------------------------------------------------
c     set the cross-validation to use the leave-one-out method
c-------------------------------------------------------------------------------

      verbose2 = .false.
      npred = nrow - ncross
      nrank = min(ncol,npred-1,nlv)

c-------------------------------------------------------------------------------
c     compute the PLSR model fit for each rank with leave N out crossval.
c-------------------------------------------------------------------------------
      if (verbosetmp) open(unit=99,file='table1.dat')
            
      do m = 1, nrank

         call int2str(m,chartemp)

         if (verbose1) then       
            combine='predict'//chartemp(ftrim(chartemp):btrim(chartemp))
     &           //'.dat'         
            open(unit=32,file=combine)
         endif

         do i = 1, nrow
            yfit(i) = 0.0d0
            do j = 1, ncol
               yfit(i) = yfit(i) + b(j,m)*x(i,j)
            end do
            yfit(i) = yfit(i)*ystd + ybar
         end do

         rsq = rsquared (nrow,yraw,yfit)
         rave = abserror (nrow,yraw,yfit)
         rrms = rmserror (nrow,yraw,yfit)

c-------------------------------------------------------------------------------
c     get the cross-validated predictions for the current rank
c-------------------------------------------------------------------------------

         call crossval(randtest,ncross,nrow,npred,ncol,nlv,name,yraw
     &        ,ybar,xraw,nsc,ntest,m,nvtold,ypred)
       
         qsq  = qsquared (nrow,yraw,ypred,ybar)
         qave = abserror (nrow,yraw,ypred)
         qrms = rmserror (nrow,yraw,ypred)
         q2(m)= qsq

c-------------------------------------------------------------------------------
c     print the results for the data set at current PLS rank
c-------------------------------------------------------------------------------

         if (verbose1) then
            write (6,10)  m
 10         format (/,' Targets, Linear Fit and Cross-Validated',
     &              ' Predictions for Rank',i3,' :',/)
            do i = 1, nrow
               write (6,20)  i,comp(i),yraw(i),yfit(i),ypred(i)
               write (32,21) comp(i),yraw(i),yfit(i),ypred(i)
 20            format (i8,6x,a7,3(1x,1f12.4))
 21            format (a7,3(1x,1f12.4),3x,'train')
               
            end do
            write (6,30)  rsq,qsq,rave,qave,rrms,qrms
 30         format (/,' R^2 & Q^2 Correlations :',5x,2(1x,1f12.4),
     &           /,' Average Absolute Error :',5x,2(1x,1f12.4),
     &           /,' Root Mean Squared Error:',5x,2(1x,1f12.4))
         endif
c-------------------------------------------------------------------------------
c     compute the PLSR model fit for the test set
c-------------------------------------------------------------------------------

         if (ntest .ne. 0) then
            open(unit=69,file='simul.qsa',status='unknown')

            if (verbose1) then
               write (6,40)
 40            format (/,'Target and Predicted Values for Test Set :',/)
            endif
                                     
            do i = 1, ntest
               k = nrow + i
               yfit(k) = 0.0d0
               do j = 1, ncol
                  yfit(k) = yfit(k) + b(j,m)*x(k,j)
               end do
               yfit(k) = yfit(k)*ystd + ybar
               extrap = 0
               do j = 1, ncol
                  if (x(k,j) .gt. xmax(j))  extrap = extrap + 1
                  if (x(k,j) .lt. xmin(j))  extrap = extrap + 1
               end do
               if (verbose1) then               
                  if (extrap .eq. 0) then
                     write (6,50)  i,comp(k),yraw(k),yfit(k)
 50                  format (i8,6x,a7,2(1x,1f12.4))
                  else
                     write (6,60)  i,comp(k),yraw(k),yfit(k),extrap  
                  end if
                  write(32,61) comp(k),yraw(k),yfit(k), 0.0d0
 60               format (i8,6x,a7,2(1x,1f12.4),3x,'with',
     &                 i5,' Extrapolation')
 61               format (a7,3f12.4,3x,'test')
               endif
               write(69,*) yraw(k),yfit(k)
            end do
            qsq2  = qsquared (ntest,yraw(nrow+1),yfit(nrow+1),ybar)
            qave2 = abserror (ntest,yraw(nrow+1),yfit(nrow+1))
            qrms2 = rmserror (ntest,yraw(nrow+1),yfit(nrow+1))

            if (verbose1) then
               write (6,70)  qsq2,qave2,qrms2
 70            format (/,' Predictive R^2 Value   :',5x,f12.4,
     &              /,' Average Absolute Error :',5x,f12.4,
     &              /,' Root Mean Squared Error:',5x,f12.4)
            endif
            close(69)
            if(verbosetmp) write(99,62) m,rsq,qsq,rave,qave,rrms,
     &           qrms,qsq2,qave2,qrms2
 62         format (i3,3x,9f12.4)
         else
            if (verbosetmp) write(99,63) m,rsq,qsq,rave,qave,rrms,qrms
 63         format (i3,3x,6(1x,1f12.4))
         end if

         if (verbose1) close(32)

      end do

      close(99)

      return
      end
c
c===============================================================================
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine invert  --  gauss-jordan matrix inversion  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "invert" inverts a matrix using the Gauss-Jordan method
c
c     n    logical dimension of the matrix to be inverted
c     np   physical dimension of the matrix storage area
c     a    matrix to invert; contains inverse on exit
c
c
      subroutine invert (n,np,a)
      implicit none
      integer maxinv
      parameter (maxinv=100)
      integer i,j,k,n,np,icol,irow
      integer ipivot(maxinv)
      integer indxc(maxinv)
      integer indxr(maxinv)
      real*8 big,temp,pivot
      real*8 a(np,np)
c
c
c     check to see if the matrix is too large to handle
c
      if (n .gt. maxinv) then
         write (6,10)
   10    format (/,' INVERT  --  Matrix Too Large; Increase MAXINV')
         stop
      end if
c
c     perform matrix inversion via the Gauss-Jordan algorithm
c
      do i = 1, n
         ipivot(i) = 0
      end do
      icol = 0
      irow = 0 
      do i = 1, n
         big = 0.0d0
         do j = 1, n
            if (ipivot(j) .ne. 1) then
               do k = 1, n
                  if (ipivot(k) .eq. 0) then
                     if (abs(a(j,k)) .ge. big) then
                        big = abs(a(j,k))
                        irow = j
                        icol = k
                     end if
                  else if (ipivot(k) .gt. 1) then
                     write (6,20)
   20                format (/,' INVERT  --  Cannot Invert',
     &                          ' a Singular Matrix')
                     stop
                  end if
               end do
            end if
         end do
         ipivot(icol) = ipivot(icol) + 1
         if (irow .ne. icol) then
            do j = 1, n
               temp = a(irow,j)
               a(irow,j) = a(icol,j)
               a(icol,j) = temp
            end do
         end if
         indxr(i) = irow
         indxc(i) = icol
         if (a(icol,icol) .eq. 0.0d0) then
            write (6,30)
   30       format (/,' INVERT  --  Cannot Invert a Singular Matrix')
            stop
         end if
         pivot = a(icol,icol)
         a(icol,icol) = 1.0d0
         do j = 1, n
            a(icol,j) = a(icol,j) / pivot
         end do
         do j = 1, n
            if (j .ne. icol) then
               temp = a(j,icol)
               a(j,icol) = 0.0d0
               do k = 1, n
                  a(j,k) = a(j,k) - a(icol,k)*temp
               end do
            end if
         end do
      end do
      do i = n, 1, -1
         if (indxr(i) .ne. indxc(i)) then
            do k = 1, n
               temp = a(k,indxr(i))
               a(k,indxr(i)) = a(k,indxc(i))
               a(k,indxc(i)) = temp
            end do
         end if
      end do
      return
      end
c
c===============================================================================
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function rsquared  --  get linear correlation coefficient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rsquared" computes the correlation coefficient between a set
c     of observed values "x" and a set of calculated values "y"
c
c
      function rsquared (n,x,y)
      implicit none
      integer maxrow
      parameter (maxrow=100)
      integer i,n
      real*8 rsquared
      real*8 xyterm,x2term,y2term
      real*8 xbar,xdiff(maxrow)
      real*8 ybar,ydiff(maxrow)
      real*8 x(maxrow),y(maxrow)
c
c
      xbar = 0.0d0
      ybar = 0.0d0
      do i = 1, n
         xbar = xbar + x(i)
         ybar = ybar + y(i)
      end do
      xbar = xbar / dble(n)
      ybar = ybar / dble(n)
      do i = 1, n
         xdiff(i) = x(i) - xbar
         ydiff(i) = y(i) - ybar
      end do
      xyterm = 0.0d0
      x2term = 0.0d0
      y2term = 0.0d0
      do i = 1, n
         xyterm = xyterm + xdiff(i)*ydiff(i)
         x2term = x2term + xdiff(i)**2
         y2term = y2term + ydiff(i)**2
      end do
      rsquared = xyterm**2 / (x2term*y2term)
      return
      end
c
c===============================================================================
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function qsquared  --  predictive correlation coefficient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "qsquared" computes the cross-validated predictive correlation
c     coefficient between a set of observed values "x" and a set of
c     calculated values "y"; the mean of the training set (which is
c     not necessarily the mean of the x's) is input via "xbar"
c
c
      function qsquared (n,x,y,xbar)
      implicit none
      integer maxrow
      parameter (maxrow=100)
      integer i,n
      real*8 qsquared
      real*8 xbar,yxterm,xterm
      real*8 x(maxrow),y(maxrow)
c
c
c     find the mean observed value (currently passed as "xbar")
c
c     xbar = 0.0d0
c     do i = 1, n
c        xbar = xbar + x(i)
c     end do
c     xbar = xbar / dble(n)
c
c     now, compute the predictive correlation coefficient
c
      yxterm = 0.0d0
      xterm = 0.0d0
      do i = 1, n
         yxterm = yxterm + (y(i)-x(i))**2
         xterm = xterm + (x(i)-xbar)**2
      end do
      qsquared = 0.0d0
      if (xterm .ne. 0.0d0)  qsquared = 1.0d0 - yxterm/xterm
      return
      end
c
c===============================================================================
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function abserror  --  compute the average absolute error  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "abserror" computes the average absolute error between a set
c     of "observed" values, x, and a set of "calculated" values, y
c
c
      function abserror (n,x,y)
      implicit none
      integer maxrow
      parameter (maxrow=100)
      integer i,n
      real*8 abserror,sum
      real*8 x(maxrow),y(maxrow)
c
c
      sum = 0.0d0
      do i = 1, n
         sum = sum + abs(x(i)-y(i))
      end do
      abserror = sum / dble(n)
      return
      end
c
c===============================================================================
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function rmserror  --  find the root mean squared error  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "rmserror" computes the root mean squared error between a set
c     of "observed" values, x, and a set of "calculated" values, y
c
c
      function rmserror (n,x,y)
      implicit none
      integer maxrow
      parameter (maxrow=100)
      integer i,n
      real*8 rmserror,sum
      real*8 x(maxrow),y(maxrow)
c
c
      sum = 0.0d0
      do i = 1, n
         sum = sum + (x(i)-y(i))**2
      end do
      rmserror = sqrt(sum/dble(n))
      return
      end
c
c===============================================================================
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine upcase  --  convert string to all upper case  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "upcase" converts a text string to all upper case letters
c
c
      subroutine upcase (string)
      implicit none
      integer i,length,code,ichar
      character*1 char
      character*(*) string
c
c
c     move through the string one character at a time,
c     converting lower case letters to upper case
c
      length = len(string)
      do i = 1, length
         code = ichar(string(i:i))
         if (code.ge.97 .and. code.le.122)
     &      string(i:i) = char(code-32)
      end do
      return
      end
c
c===============================================================================
c
      subroutine help
c
c-------------------------------------------------------------------------------
c --- help messages
c-------------------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)

      write(6,*)
      write(6,*) 'To run the program type:'
      write(6,*)
      write(6,*) 'combine -i <input file> -o <out file>'
      write(6,*)
      stop

      end
c
c===============================================================================
c
      subroutine setime

c-------------------------------------------------------------------------------
c     ##  subroutine setime  --  initialize elapsed CPU time clock  ##
c         "setime" initializes the elapsed interval CPU timer 
c-------------------------------------------------------------------------------

      implicit none
      real*8 cputim
      common /chrono/ cputim
c
c initialize interval at elapsed CPU time for current job 
c
      call clock (cputim)

      return
      END
c
c===============================================================================
c
      subroutine getime (elapsed)

c-------------------------------------------------------------------------------
c     ##  subroutine getime  --  get elapsed CPU time in seconds  ##
c        "getime" gets elapsed CPU time in seconds for an interval
c-------------------------------------------------------------------------------

      implicit none
      real*8 cputim
      common /chrono/ cputim
      real*8 elapsed
c
c  elapsed time for the interval is the current total CPU 
c  time minus the total time at the start of the interval 
c
      call clock (elapsed)
      elapsed = elapsed - cputim

      return
      END
c
c===============================================================================
c
      subroutine clock (seconds)

c-------------------------------------------------------------------------------
c     ##  subroutine clock  --  find elapsed time for current job  ##
c        "clock" determines elapsed CPU time in seconds since the  
c         start of the job; only one of the implementations should  
c         be activated by removing comment characters from the code 
c-------------------------------------------------------------------------------

      implicit none
      real*8 seconds
c
c Unix machines have access to the "etime" intrinsic,  
c this code works for Sun, DEC, SGI, Convex and others 
c
      real etime,times(2)
      seconds = dble(etime(times))
c
      return
      END
c
c================================================================================
c
      subroutine storetime(subname,ttime,ntime)
c
      implicit     none
      integer      maxcall,ntime
      parameter    (maxcall = 100)
      real*8       time,ttime,ftme
      character*80  part
      character*(*) subname

      common /timings/ ftme(maxcall),part(maxcall)

      call getime(time)

      ntime = ntime + 1
      part(ntime) = subname
      ftme(ntime) = time
      ttime = ttime + time

      return
      END
c
c================================================================================
c
      subroutine printime(ntime)

      implicit      none
      integer       maxcall,ntime,i
      parameter     (maxcall = 100)
      real*8        ftme
      character*80  part

      common /timings/ ftme(maxcall),part(maxcall)

      write(16,'(/,a,/)') '==> Timings: '

      do i = 1, ntime
        write (16,33)  part(i),ftme(i)
      enddo
  33  format ('   < ',a22,f12.3,' sec >')

      return
      END
c
c================================================================================
c
      subroutine readinp(nf,idel,nsc,imat,nlv,ielec,ncomp,ntest,nad,
     &     randtest,ncross,dielect,ptr,cutptr,comp,drugname,y,add)

c-------------------------------------------------------------------------------
c --- Read input file (with filenames, drug names and activities)
c-------------------------------------------------------------------------------

      IMPLICIT      DOUBLE PRECISION (A-H,O-Z)

      parameter    (maxcompl = 100)
      parameter    (maxrow   = 100)
      parameter    (maxadd   =  10)
      
      integer ptr, idel, nsc,imat, nlv
      integer ielec, iad,nad,ncomp,ntest
      integer randtest,ncross

      real*8 dielect, cutptr

      dimension     y(maxrow),add(maxcompl,maxadd)
      character*4   drugname(maxcompl)
      character*40  comp(maxcompl)

      read(nf,'(3i5)') idel,nsc,imat
      read(nf,'(2i5,f10.5)') nlv,ielec,dielect
      read(nf,'(2i5)') iad,nad
      read(nf,'(2i5)') ncomp,ntest
      read(nf,'(i5,f10.5)') ptr,cutptr
      read(nf,'(2i5)') randtest,ncross
      
      write(16,*) '==> Compounds in training set: ', ncomp
      if (iad.eq.0) then
         do i = 1, ncomp
            read(nf,'(a8,3x,a4,f10.5)')comp(i), drugname(i), y(i)
            write(16,'(a8,3x,a4,f10.5)')comp(i), drugname(i), y(i)
         enddo
         if (ntest .ne. 0) then
            write(16,*) '==> Compounds in test set: ', ntest
            n=ncomp
            do i = 1, ntest
               n=n+1
               read(nf,'(a8,3x,a4,f10.5)') comp(n), drugname(n), y(n)
               write(16,'(a8,3x,a4,f10.5)')comp(n), drugname(n), y(n)
            enddo
         endif  
      else
        do i = 1, ncomp
          read(nf,10) comp(i),drugname(i),y(i),(add(i,j),j=1,nad)
          write(16,10) comp(i),drugname(i),y(i),(add(i,j),j=1,nad)
        enddo
        if (ntest .ne. 0) then
          write(16,*) '==> Compounds in test set: ', ntest
          n=ncomp
          do i = 1, ntest
            n=n+1
            read(nf,10) comp(n),drugname(n),y(n),(add(n,j),j=1,nad)
            write(16,10) comp(n),drugname(n),y(n),(add(n,j),j=1,nad)
          enddo
        endif
      endif
   10 format(a8,3x,a4,f10.5,10(1x,1f10.5))
c  10 format(a7,4x,a4,f10.5,<nad>f8.4)
c    To be compliant with Absoft 4.5
      close(nf)

      return
      end
c
c================================================================================
c
      subroutine scramble(nscram,nrow,y,yr)

c------------------------------------------------------------------------------
c     Initializes the random number generator. Then tests
c     whether random numbers are being correctly generated.
c     Finally permutes randomly the y vector NSCRAM times.
c     Results are stored in the array yr for further use.
c------------------------------------------------------------------------------

      parameter (maxrow=100)
      parameter (maxscr=100)

      common/rands/irn

      real*8        y(maxrow),yr(maxrow,maxscr)
      real          r(maxrow)
      integer       year,month,day
      integer       hour,minute,second
      integer       seed,irandom,indx(maxrow)

c------------------------------------------------------------------------------
c preparation of irn
c------------------------------------------------------------------------------

      seed = 0
      call calendar_IRIX (year,month,day,hour,minute,second)
      year = mod(year,10)
      seed = seed + 32140800*year + 2678400*(month-1)
      seed = seed + 86400*(day-1) + 3600*hour
      seed = seed + 60*minute + second
      do i = 1, 100

         irandom  = int(rand(seed)*1000000)
c         para compilaciÃ³n en Absoft
c       irandom  = int(ran(seed)*1000000)
      enddo
      IRN=IRANDOM

c------------------------------------------------------------------------------
c test generator
c------------------------------------------------------------------------------

      num  = 100
      nran = 100
      numr = 10
      write(16,'(a,i8)') '==> Number of random trials : ', nran
      write(16,'(a,i8)') '==> Number of random numbers: ', numr

      ise=0.0
      do ktrial = 1, nran
        do i = 1, numr
          r(i) = rand(ise)
        enddo
        write(16,10) '# trial: ', ktrial, (r(k),k=1,numr)
      enddo
c 10  format(a,i3,' | ',<numr>f5.2)
  10  format(a,i3,' | ',10(1x,f5.2))

      call rantest

c------------------------------------------------------------------------------
c Assign random activities
c------------------------------------------------------------------------------

c     write(16,*) (y(i), i=1,nrow)
      do ktrial = 1, nscram
        do i = 1, nrow
          r(i) = rand(ise)
        enddo
        call indexx(nrow,r,indx)
        do i = 1, nrow
           yr(i,ktrial) = y(indx(i))           
        enddo
c       write(16,*) ktrial, (yr(i,ktrial), i=1,nrow)
      enddo

c------------------------------------------------------------------------------
c End of subroutine
c------------------------------------------------------------------------------

      return
      end
c
c$$$c==============================================================================
c$$$c
c$$$      real*4 function rand(idup)
c$$$c------------------------------------------------------------------------------
c$$$c     this function converts ranu to single precision
c$$$c     additionaly, there is no need for rands common
c$$$c     in every subroutine that uses ranu, but still ranu
c$$$c     has to be initialized before first attempt to use it.
c$$$c------------------------------------------------------------------------------
c$$$      real*8 ranuS
c$$$c     real*4 rdum
c$$$      real*8 rdum
c$$$      integer*4 idup
c$$$      common/rands/irn
c$$$      rdum=ranu(0)
c$$$      rand=sngl((rdum ))
c$$$      return
c$$$      end
c
c==============================================================================
c
      DOUBLE PRECISION FUNCTION RANU(IQFD)                              MOV08330
C-----------------------------------------------------------------------MOV08340
C     RANDOM NUMBER GENERATOR                                           MOV08360
C     SEPTEMBER, 1986 WLJ,JMB                                           MOV08370
C                                                                       MOV08380
C     NOTES ON USAGE:                                                   MOV08400
C     THE STATEMENT "COMMON/RANDS/IRN" MUST APPEAR IN EACH              MOV08410
C     PROGRAM,SUBROUTINE... WHERE THE SEED VALUE (IRN) IS               MOV08420
C     BEING WRITTEN OUT OR BEING READ IN. THIS FUNCTION IS              MOV08430
C     USED IN THE SAME MANNER AS THAT USED ON THE GOULD,                MOV08440
C     (i.e. "RANDN = RANU(0)" ). IRN SHOULD INITIALLY BE                MOV08450
C     AN ODD, I6 INTEGER BUT CAN GET AS BIG AS I7. JMB                  MOV08460
C     IMOD-1 UNIQUE RANDOM NUMBERS ARE GENERATED . A SINGLE             MOV08470
C     PRECISION LINEAR CONGRUENTIAL GENERATOR IS USED WITH              MOV08480
C     SHUFFLING OF THE CONSTANT AND OF WHEN THE SHUFFLING IS            MOV08490
C     DONE. CONSEQUENTLY, THE PERIOD IS EXTREMELY LONG - NONE           MOV08500
C     WAS FOUND IN TESTS GENERATING MILLIONS OF RANDOM NUMBERS.         MOV08510
C-----------------------------------------------------------------------MOV08530
      COMMON/RANDS/IRN                                                  MOV08540
      DATA ICNT/0/,ICHG/1167/,ICHG0/1167/                               MOV08550
      DATA IMUL/1173/,ICON/458753759/                                   MOV0856
      DATA LMOD/1048573/,ICN0/458753759/                                MOV08570
      DATA JMUL/1161/,JCON/458716759/                                   MOV08580
      DATA KMOD/1048573/,JRN/124690/                                    MOV08590
      ICNT = ICNT + 1                                                   MOV08600
      IF (ICNT .NE. ICHG) GO TO 10                                      MOV08610
C-----------------------------------------------------------------------MOV08620
C     CHANGE ICON USING SECONDARY GENERATOR                             MOV08640
C-----------------------------------------------------------------------MOV08660
      JRN = JRN*JMUL + JCON                                             MOV08670
      JRN = MOD(JRN,KMOD)                                               MOV08680
      RNJ = FLOAT(JRN)/FLOAT(KMOD)                                      MOV08690
      IF (RNJ .GT. 0.5E+00) FAC = 1.0E+00 + 0.5E+00*RNJ                 MOV08700
      IF (RNJ .LE. 0.5E+00) FAC = 1.0E+00 - 0.5E+00*RNJ                 MOV08710
      FAC = FLOAT(ICN0)*FAC                                             MOV08720
      ICON = IFIX(FAC)                                                  MOV08730
C-----------------------------------------------------------------------MOV08740
C     CHANGE ICHG USING SECONDARY GENERATOR                             MOV08760
C-----------------------------------------------------------------------MOV08780
      JRN = JRN*JMUL + JCON                                             MOV08790
      JRN = MOD(JRN,KMOD)                                               MOV08800
      RNJ = FLOAT(JRN)/FLOAT(KMOD)                                      MOV08810
      IF (RNJ .GT. 0.5E+00) FAC = 1.0E+00 + 0.5E+00*RNJ                 MOV08820
      IF (RNJ .LE. 0.5E+00) FAC = 1.0E+00 - 0.5E+00*RNJ                 MOV08830
      FAC = FLOAT(ICHG0)*FAC                                            MOV08840
      ICHG = IFIX(FAC)                                                  MOV08850
      ICNT = 0                                                          MOV08860
C-----------------------------------------------------------------------MOV08870
C     GENERATE RANDOM NUMBER                                            MOV08890
C-----------------------------------------------------------------------MOV08910
   10 IRN = IRN*IMUL + ICON                                             MOV08920
      IRN = MOD(IRN,LMOD)                                               MOV08930
      RANU = DABS(DBLE(IRN)/DBLE(LMOD))                                 MOV08940
      RETURN                                                            MOV08950
      END                                                               MOV08960
c
c==============================================================================
c
        subroutine calendar_IRIX (year,month,day,hour,minute,second)

c------------------------------------------------------------------------------
c
c     ##  subroutine calendar_IRIX  --  find the current date and time  ##
c
c     "calendar" returns the current time as a set of integer
c     values representing the year, month, day, hour, minute
c     and second; only one of the various machine ementations
c     included should be activated by removing comment characters
c
c------------------------------------------------------------------------------

      implicit none
      integer year,month,day
      integer hour,minute,second
      integer myDate(3)
c
c     Unix machines use the "itime" and "idate" intrinsics,
c     this code works for Sun, DEC, SGI, Convex and others
c
      integer hms(3)
      call itime (hms)
      hour = hms(1)
      minute = hms(2)
      second = hms(3)
      
      
c      CALL VXTIDATE_Y2KBUG (month,day,year)
      call idate(myDate)
      day = myDate(1)
      month = myDate(2)
      year = myDate(3)

c      call idate (month,day,year)
      return
      end
c
c==============================================================================
c
      subroutine rantest

c------------------------------------------------------------------------------
c --- test of the random number generator: Numerical computation of PI
c --- The number value should tend to converge to the actual numbers.
c------------------------------------------------------------------------------

      common/rands/irn

      parameter (pi=3.1415926)
      dimension iy(3),yprob(3)

      fnc(x1,x2,x3,x4) = sqrt(x1**2+x2**2+x3**2+x4**2)

      do i = 1, 3
        iy(i)=0
      enddo

      write(16,*)
      write(16,'(a)') '* * * TEST OF THE RANDOM NUMBER GENERATOR * * *'
      write(16,*)

      ise=0.0
      do j = 1, 21
        do k = 2**(j-1),2**j
          x1=rand(ise)
          x2=rand(ise)
          x3=rand(ise)
          x4=rand(ise)
          if (fnc(x1,x2,0.0,0.0).lt.1.0) iy(1) = iy(1)+1
          if (fnc(x1,x2,x3 ,0.0).lt.1.0) iy(2) = iy(2)+1
          if (fnc(x1,x2,x3 ,x4 ).lt.1.0) iy(3) = iy(3)+1
        enddo
        do i = 1, 3
          yprob(i) = 1.0*(2**(i+1))*iy(i)/(2**j)
        enddo
        write(16,'(1x,i8,3f12.6)') 2**j, (yprob(i),i=1,3)
      enddo
      write(16,'(1x,/,t4,a,3f12.6,/)')'actual',pi,4.0*pi/3.0,0.5*(pi**2)

      write(16,'(/,10x,a,/)') '* * * END OF TEST * * *'

      return
      end
c
c==============================================================================
c
      SUBROUTINE indexx(n,arr,indx)
      parameter (maxrow=100)
      INTEGER n,indx(maxrow),M,NSTACK
      REAL arr(maxrow)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
         a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)        
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

c========================================================================

      SUBROUTINE pretreat_xmatrix(nvt,nvtold,nad,ncomp,x,
     &    ptr,cutpret,varselect)

c------------------------------------------------------------------------
c     
c     Selects variables from the Xmatrix having standard deviations
c     are above a certain cutoff     
c
c------------------------------------------------------------------------
      
      implicit none
      
      integer maxrow,maxcol
      integer i,j
      integer ptr

      parameter (maxrow= 100)
      parameter (maxcol=1500)
      
      integer nvar,nvt,ncomp,nad,nvtold
      integer varselect(maxcol)

      real*8 x(maxrow,maxcol),xpret(maxrow,maxcol)
      real*8 xstd(maxcol),xbar(maxcol)
      real*8 cutpret,cut1,cut2
     

c------------------------------------------------------------------------
c     Noise reduccion
c------------------------------------------------------------------------
      
      nvtold = nvt
      cut1 = 0.0
      
      if(ptr.eq.3)then
          write (6,90)
 90       format (/,' Zeroing all positive energy values ',/)         
         do j = 1, nvt-nad
            do i = 1, ncomp
               if(x(i,j).gt.cut1) x(i,j) = 0.0d0
            end do
         enddo   
      elseif(ptr.eq.1)then
          write (6,94)
 94       format (/,' Zeroing positive van der Waal energy values ',/)    
         do j = 1, (nvt-nad)/2
            do i = 1, ncomp
               if(x(i,j).gt.cut1) x(i,j) = 0.0d0
            end do
         enddo   
      elseif(ptr.eq.2)then
         write (6,95)
 95      format (/,' Zeroing positive electrostatic energy values ',/)
         do j = (nvt-nad)/2+1,nvt-nad
            do i = 1, ncomp
               if(x(i,j).gt.cut1) x(i,j) = 0.0d0
            end do
         enddo         
      endif
           
      cut2=cutpret

c------------------------------------------------------------------------
c     Standard deviation computation
c------------------------------------------------------------------------

      if (cut2.ne.0.0)then

         write (6,96) cut2
 96      format (/,' Eliminate Var. std dev <', f6.4, /)

         do j = 1, nvt
            xbar(j) = 0.0d0
            do i = 1, ncomp
               xbar(j) = xbar(j) + x(i,j)
            end do
            xbar(j) = xbar(j) / dble(ncomp)
         enddo
         do j = 1, nvt
            xstd(j) = 0.0d0
            do i = 1, ncomp
               xstd(j) = xstd(j) + (x(i,j)-xbar(j))**2
            end do
            xstd(j) = sqrt(xstd(j)/dble(ncomp-1))
         enddo  
         
c------------------------------------------------------------------------
c     Variable selection
c------------------------------------------------------------------------
         
         nvar = 0      
         do i=1,nvt
            varselect(i)=0
            if (xstd(i).gt.cut2)then
               nvar = nvar + 1
c     print *,'variable selected',i,xstd(i)
               do j=1,ncomp
                  xpret(j,nvar) = x(j,i)
               enddo
               varselect(i)=1            
            endif         
         enddo

c------------------------------------------------------------------------
c     Reasigment of Xmatrix
c------------------------------------------------------------------------

         do i=1,nvar
            do j=1,ncomp
               x(j,i) = 0.0d0
               x(j,i) = xpret(j,i)
            enddo
         enddo

         do i = nvar+1,nvt
            do j=1,ncomp
               x(j,i) = 0.0d0
            enddo
         enddo
         nvt = nvar
      else   
         do i=1,nvt
            varselect(i)=1 
         enddo
      endif

      return
      end


c========================================================================
c     
c     <<<<<<<<<<<<<<<<<<<<<<<  ELIPSE.F >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c     This subrutine calculates the model aplicability domain as defined 
c     by Hotelling's T^2 (as a generalization of Student's t-test).
c
c========================================================================

      subroutine elipse(t,ncomp,nrank)

c------------------------------------------------------------------------
c     Variable declaration
c------------------------------------------------------------------------
      
      implicit none
      
      integer  maxlength,maxvar,steps,maxrank
                       
      parameter (maxlength=100)
      parameter (maxvar=2)
      parameter (steps=100)
      parameter (maxrank=20)

      integer i,j,k,s,vectsize
      integer ftrim,btrim
      integer ncomp,nrank
      

      real*8 datos(maxlength,maxlength)
      real*8 mean(maxvar)
      real*8 cov(maxvar,maxvar)
      real*8 elip(steps,maxvar),alpha1(steps),radii1(steps)
      real*8 t(maxlength,maxlength)

      real*8 a,b,d,c,f,g,x0,y0,aa,bb,phi,co,si
      real*8 theta(steps)

      character*10 chartemp
      character*20 filename

      include "param.h"

c------------------------------------------------------------------------
c     Covariance matrix calculation
c------------------------------------------------------------------------

      vectsize = ncomp

      do s =1,nrank     
         if(s.eq.1)then
            do j=1,vectsize      
               datos(j,1)=t(j,1)
               datos(j,2)=t(j,2)
            enddo
         else
            do j=1,vectsize 
               datos(j,1)=t(j,s-1)
               datos(j,2)=t(j,s)
            enddo 
         endif
         
         do i=1,2
            mean(i)=0.0d0
            do j=1,vectsize
               mean(i)=mean(i)+datos(j,i)
            enddo
            mean(i)=mean(i)/float(vectsize)
         enddo
   

         do i=1,2
            do j=1,2
               cov(i,j)=0.0d0
               do k=1,vectsize
                  cov(i,j)=cov(i,j)+ 
     &                 (datos(k,i)-mean(i))*(datos(k,j)-mean(j))
               enddo
               cov(i,j)=cov(i,j)/float(vectsize)            
            enddo
         enddo
         
c------------------------------------------------------------------------
c     Elipse parameter calculation
c------------------------------------------------------------------------

         call invert(maxvar,maxvar,cov)    

         a = cov(1,1)
         b = cov(1,2)
         d=0
         c = cov(2,2)
         f=0
         g = -2*(float(vectsize)+1)*(float(vectsize)-1)/(float(vectsize)
     &        *(vectsize-2))*fval(vectsize)
         
         x0 = (c*d - b*f)/(b*b-a*c)
         y0 = (a*f - b*d)/(b*b-a*c)
         
         phi= (1/2)*atan((2*b)/(c-a))
         
         aa =sqrt(2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)/((b*b-a*c)*((c-a)*
     &        sqrt(1+(4*b*b)/((a-c)*(a-c)))-c-a)))
         bb =sqrt(2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)/((b*b-a*c)*((a-c)*
     &        sqrt(1+(4*b*b)/((a-c)*(a-c)))-c-a)))
         
         co=cos(phi)
         si=sin(phi)

         call int2str(s,chartemp)

         filename = 'ellipse'//chartemp(ftrim(chartemp):btrim(chartemp))
     &        //'.dat'
         open (unit=99,file=filename)

         write (99,199)
 199     format (/,'# Hottelings confidence ellipse for the two last latent 
     &variable score cff:',/)
         
         do i=1,steps
            theta(i)=i*(6.2831/steps)
            elip(i,1)=0.0d0
            elip(i,2)=0.0d0
            elip(i,1)=aa*cos(theta(i))*co-si*bb*sin(theta(i))+x0
            elip(i,2)=aa*cos(theta(i))*si+co*bb*sin(theta(i))+y0
            
            alpha1(i)=atan(elip(i,2)/elip(i,1))
            if(elip(i,1).lt.0.0d0.and.elip(i,2).gt.0.0d0)then
               alpha1(i)=3.1416+alpha1(i)
            elseif(elip(i,1).lt.0.0d0.and.elip(i,2).lt.0.0d0)then
               alpha1(i)=3.1416+alpha1(i)
            elseif(elip(i,1).gt.0.0d0.and.elip(i,2).lt.0.0d0)then
               alpha1(i)=6.2912+alpha1(i)
            endif
            
            radii1(i)=sqrt(elip(i,1)**2+elip(i,2)**2)
            write (99,*) (elip(i,j),j=1,2)
         enddo
            
c------------------------------------------------------------------------
c     Outlayer detection 
c------------------------------------------------------------------------ 

c$$$      write(99,'(/)')     
c$$$
c$$$      do i=1,vectsize
c$$$         alpha2(i)=atan(datos(i,2)/datos(i,1))
c$$$         if(datos(i,1).lt.0.0d0.and.datos(i,2).gt.0.0d0)then
c$$$            alpha2(i)=3.1416+alpha2(i)
c$$$         elseif(datos(i,1).lt.0.0d0.and.datos(i,2).lt.0.0d0)then
c$$$            alpha2(i)=3.1416+alpha2(i)
c$$$         elseif(datos(i,1).gt.0.0d0.and.datos(i,2).lt.0.0d0)then
c$$$            alpha2(i)=6.2912+alpha2(i)
c$$$         endif
c$$$         radii2(i)=sqrt(datos(i,1)**2+datos(i,2)**2)
c$$$         
c$$$         do j=1,steps
c$$$            if(alpha1(j).gt.alpha2(i))then
c$$$               if(radii1(j).lt.radii2(i))
c$$$     &              then
c$$$                  write(99,299) i 
c$$$               endif
c$$$               EXIT
c$$$            endif
c$$$         enddo    
c$$$      enddo
c$$$
c$$$ 299  format('pound',i3,1x,'could be an outlayer')
c$$$
         close (99)

c------------------------------------------------------------------------
c     End
c------------------------------------------------------------------------     

      enddo
      return
      end

c===============================================================================
c     
c     This subrutine provides a pdb file with extracolum corresponding
c     to reescaled (from 10 to -10) model coefficents
c
c===============================================================================

      subroutine coefffiles(b,nrank,nvar,nres,dres,ipress,igraph,labres
     &     ,xyz,multres,varselect)

c------------------------------------------------------------------------
c     Variable declaration
c------------------------------------------------------------------------
      
      implicit none

      integer maxrow,maxcol,maxres,maxatom
      
      parameter (maxrow= 100)
      parameter (maxcol=1500)
      parameter (maxres=700)
      parameter (maxatom=8000)

      integer i,j,k,l,m,n,nvr
      integer init,iend
      integer varselect(maxcol)

      real*8 coefs(maxcol)
      real*8 cmin,cmax

      character*10 chartemp
      character*14  filename

c------------ Input param -----------------------------------------------
      
      integer nrank,nvar,nres,dres
      integer ipress(maxres)
      integer multres,drug
      integer ftrim,btrim

      character*4 igraph(maxatom),labres(maxres)

      real*8 b(maxcol,maxcol)
      real*8 xyz(3,maxatom)    

c     --- file with vdw coefs

      do j=1,nrank

         call int2str(j,chartemp)
 
         filename = 'vdwcoef'//chartemp(ftrim(chartemp):btrim(chartemp))
     &        //'.pdb'
 
         open(unit=10,file=filename)

         cmin =  10.0
         cmax = -10.0
         n=0
         do i = 1, nvar
            if(varselect(i).eq.1)then
               n=n+1
               cmin = min(b(n,j),cmin)
               cmax = max(b(n,j),cmax)
            endif
         enddo

         m=0
         do i = 1, nvar
            if(varselect(i).eq.1)then
               m=m+1
               coefs(i) = -(b(m,j)-cmin)/(cmax-cmin)*10.0d0
            else
               coefs(i) = -(0.0d0-cmin)/(cmax-cmin)*10.0d0
            endif
         enddo
         nvr = 0
         do i = 1, nres
            if (multres.eq.0)then                       
               if (i.eq.dres) then !!
c---------Uncomment to include ligand in the file------
c$$$               init=ipress(i)
c$$$               iend=ipress(i+1)-1
c$$$               do k = init, iend
c$$$                  write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
c$$$     &                 'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
c$$$     &                 '  1.00', 0.0d0
c$$$               enddo
               else
                  nvr=nvr+1
                  init=ipress(i)
                  iend=ipress(i+1)-1
                  do k = init,iend
                     if (igraph(k)(1:1).eq.'H') then 
                        
                     else
                       write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &                       'ATOM ',k,igraph(k),labres(i),i,
     &                       (xyz(l,k),l=1,3),'  1.00',coefs(nvr)
                     endif
                  enddo
               endif
            else
               drug=0
               do l=1,multres
                  if (i.eq.dres-multres+l)drug=1
               enddo
               if (drug.eq.1) then !!
c---------Uncomment to include ligand in the file------
c$$$               init=ipress(i)
c$$$               iend=ipress(i+1)-1
c$$$               do k = init, iend
c$$$                  write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
c$$$     &                 'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
c$$$     &                 '  1.00', 0.0d0
c$$$               enddo
               else
                  nvr=nvr+1
                  init=ipress(i)
                  iend=ipress(i+1)-1
                  do k = init,iend
                     if (igraph(k)(1:1).eq.'H') then 
                        
                     else
                       write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &                       'ATOM ',k,igraph(k),labres(i),i,
     &                       (xyz(l,k),l=1,3),'  1.00',coefs(nvr)
                     endif
                  enddo
               endif
            endif  
         enddo
         
         close(10)
                 
c     --- file with elec coefs

c         write(chartemp,'(i10)') j
c         read(j, '(i5)') chartemp
         call int2str(j,chartemp)
         filename = 'elecoef'//chartemp(ftrim(chartemp):btrim(chartemp))
     &        //'.pdb'         
         open(unit=10,file=filename)

         cmin =  10.0
         cmax = -10.0
         n=0
         do i = nvar+1, 2*nvar
            if(varselect(i).eq.1)then
               n=n+1 
               cmin = min(b(n,j),cmin)
               cmax = max(b(n,j),cmax)
            endif
         enddo
         do i = nvar+1, 2*nvar
            if(varselect(i).eq.1)then
               m=m+1
               coefs(i) = -(b(m,j)-cmin)/(cmax-cmin)*10.0d0
            else
               coefs(i) = -(0.0d0-cmin)/(cmax-cmin)*10.0d0
            endif               
         enddo
         do i = 1, nres
            if (multres.eq.0)then
               if (i.eq.dres) then
c-------Uncommented to include ligand in the file------ 
c$$$              init=ipress(i)
c$$$               iend=ipress(i+1)-1
c$$$               do k = init, iend
c$$$                  write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
c$$$     &                 'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
c$$$     &                 '  1.00', 0.0d0
c$$$               enddo
               else
                  nvr=nvr+1
                  init=ipress(i)
                  iend=ipress(i+1)-1
                  do k = init,iend
                     if (igraph(k)(1:1).eq.'H') then
                     else
                       write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &                       'ATOM ',k,igraph(k),labres(i),i,
     &                       (xyz(l,k),l=1,3),'  1.00',coefs(nvr)
                     endif
                  enddo
               endif
            else
               drug=0
               do l=1,multres
                  if (i.eq.dres-multres+l)drug=1
               enddo
               if (drug.eq.1) then
c--------Uncommented to include ligand in the file------ 
c$$$              init=ipress(i)
c$$$               iend=ipress(i+1)-1
c$$$               do k = init, iend
c$$$                  write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
c$$$     &                 'ATOM ',k,igraph(k),labres(i),i,(xyz(l,k),l=1,3),
c$$$     &                 '  1.00', 0.0d0
c$$$               enddo
               else
                  nvr=nvr+1
                  init=ipress(i)
                  iend=ipress(i+1)-1
                  do k = init,iend
                     if (igraph(k)(1:1).eq.'H') then
                     else
                       write(10,'(a,i6,2x,a4,a3,2x,i4,4x,3f8.3,a,f6.2)')
     &                       'ATOM ',k,igraph(k),labres(i),i,
     &                       (xyz(l,k),l=1,3),'  1.00',coefs(nvr)
                     endif
                  enddo
               endif
            endif          
         enddo
         close(10)         
      enddo
      
      return
      end

c===============================================================================
c     
c     <<<<<<< crossval >>>>>>>>>>>>>
c
c     This subrutine calculates crosvalidation for a consecutive groups of 
c     N elements, or for random grups of N elements
c
c===============================================================================

      subroutine crossval(randtest,ncross,nrow,npred,ncol,nlv,name,yraw
     &     ,ybar,xraw,nsc,ntest,m,nvtold,ypred)

c------------------------------------------------------------------------
c     Variable declaration
c------------------------------------------------------------------------
      implicit none
      
      integer maxrow,maxcol,maxrank,ncross

      parameter (maxrow=100)
      parameter (maxcol=1500)
      parameter (maxrank=20)

      integer nrow,npred,ncol,nsc,ntest,nqrank
      integer varselect(maxcol),nvtold,nlv
      integer i,j,k,l,m,n,randtest
      integer group(maxrow),low,grt

      character*4 name(maxrow),nameq(maxrow)

      real*8 yraw(maxrow),xraw(maxrow,maxcol)
      real*8 yq(maxrow),xq(maxrow,maxcol)
      real*8 xqbar(maxcol),xstd(maxcol),yqbar,yqstd,ybar
      real*8 ypred(maxrow)
      real*8 wq(maxcol)
      real*8 bq(maxcol,maxcol)
      real*8 r(maxrow)

      character*20  vardesc(maxcol)

      logical verbose

c------------------------------------------------------------------------

      verbose=.false.

c------------------------------------------------------------------------
c     Leave N out
c------------------------------------------------------------------------
      if (randtest.eq.0)then
         do i = 1,(nrow/ncross)+1
            do j = 1,(i-1)*ncross
               nameq(j) = name(j)
               yq(j) = yraw(j)
            end do
            do j = (i-1)*ncross+1,npred
               nameq(j) = name(j+ncross)
               yq(j) = yraw(j+ncross)
            end do
            do k = 1, ncol
               do j = 1, (i-1)*ncross
                  xq(j,k) = xraw(j,k)
               end do
               do j = (i-1)*ncross+1, npred
                  xq(j,k) = xraw(j+ncross,k)
               end do
            end do

            call yscale(nsc,npred,ntest,yq,yqbar,yqstd)           
c     call xscale(nsc,npred,ncol,ntest,xq,xqbar,xstd,xqstd)
            call xscale(nsc,npred,ncol,ntest,xq,xqbar,xstd,wq)
            call nipals(nlv,npred,ncol,yq,xq,nqrank,bq,vardesc,varselect
     &           ,nvtold,verbose)
            do l = (i-1)*ncross+1,(i-1)*ncross+ncross
               ypred(l) = 0.0d0
               do j = 1, ncol
                  if (xraw(l,j).eq.Z'FFFFFFFF') xraw(l,j)=0.0d0
                  ypred(l) = ypred(l)+bq(j,m)*(xraw(l,j)-xqbar(j))*wq(j)
               end do
               ypred(l) = ypred(l)*yqstd + yqbar
            enddo            
         end do
      else

c------------------------------------------------------------------------
c     Random Groups
c------------------------------------------------------------------------

         call randnum(nrow,r)
         
         do i =1,nrow
            low = 0
            grt = 0
            do j=1,nrow
               if (r(i).lt.r(j)) then
                  low = low + 1
               else
                  grt = grt + 1   
               endif
            enddo
            k=0
            do while (k*ncross.le.nrow)
               k = k +1
               if(low.lt.k*ncross)then
                  group(i)=k
                  exit
               endif
            enddo
         enddo
         
         i=0
                  
         do while (i*ncross.le.nrow)
            i = i +1
            n=0
            do j=1,nrow
               if (group(j).ne.i)then
                  n=n+1
                  nameq(n) = name(j)
                  yq(n) = yraw(j)
                  do l = 1, ncol
                     xq(n,l) = xraw(j,l)
                  enddo                   
c                  print*,yq(n),n
               endif
c               print*,j
            enddo
  
            call yscale(nsc,npred,ntest,yq,yqbar,yqstd)           
c     call xscale(nsc,npred,ncol,ntest,xq,xqbar,xstd,xqstd)
            call xscale(nsc,npred,ncol,ntest,xq,xqbar,xstd,wq)
            call nipals(nlv,npred,ncol,yq,xq,nqrank,bq,vardesc,varselect
     &           ,nvtold,verbose)
            
            do l =1,nrow
               if(group(l).eq.i)then
                  ypred(l) = 0.0d0
                  do j = 1, ncol
                     if (xraw(l,j).eq.Z'FFFFFFFF') xraw(l,j)=0.0d0
                     ypred(l) = ypred(l)+bq(j,m)*(xraw(l,j)-xqbar(j))*
     &                    wq(j)
                  end do
                  ypred(l) = ypred(l)*yqstd + yqbar
               endif
            enddo            
         end do
      endif

      return
      end

c===============================================================================
c     
c     <<<<<<<<<<randnum(n,r)>>>>>>>>>>>
c
c     Generation of random numbers
c     n = number of random numbers
c     r = output random array
c 
c===============================================================================
   
      subroutine randnum(n,r)

      implicit none
      
      integer maxrow,n

      parameter (maxrow=100)

      real*8 r(maxrow)
      real rand
      integer*4 timeArray(3)
      integer i  
      
      call itime(timeArray)
      i = rand ( timeArray(1)+timeArray(2)+timeArray(3))
  
      do i=1,n
         r(i) = rand(0)
c         print*,r(i)
      enddo

      return
      end


c===============================================================================
c
c btrim function calculates the number of blank caracters at the end
c of a string
c
c===============================================================================

      function btrim (string)
      implicit none
      integer i,size,len
      integer last,btrim
      character*1 null
      character*(*) string
c
c
c     move forward through the string, one character
c     at a time, looking for first null character
c
      btrim = 0
      null = char(0)
      size = len(string)
      last = size
      do i = 1, size
         if (string(i:i) .eq. null) then
            last = i - 1
            goto 10
         end if
      end do
   10 continue
c
c     move backward through the string, one character
c     at a time, looking for first non-blank character
c
      do i = last, 1, -1
         if (string(i:i) .gt. ' ') then		
            btrim = i
            goto 20
         end if
      end do
   20 continue
      return
      end


c===============================================================================
c
c ftrim function calculates the number of blank caracters at the begining
c of a string
c
c===============================================================================

      function ftrim (string)

      implicit none
      integer i,size,len
      integer last,ftrim
      character*1 null
      character*(*) string
ce
c
c     move forward through the string, one character
c     at a time, looking for first null character
c
      
      ftrim = 0
      null = char(0)
      size = len(string)
      last = size
      do i = 1, size
         if (string(i:i) .eq. null) then
            last = i - 1
            goto 10
         end if
      end do
   10 continue
c
c     move backward through the string, one character
c     at a time, looking for first non-blank character
c
      do i = 1,last
         if (string(i:i) .ne. ' ') then		
            ftrim = i
            goto 20
         end if
      end do
   20 continue
     
      return
      end


c=========================================================
c
c This subroutine returns a string with the integer input 
c
c=========================================================

      subroutine int2str(input, output)

      implicit none
      
      integer input
      character*(*) output
      
      integer myMod
      integer remainingNum
      character*1 addNum
      
      output = ''
      remainingNum = input

c     repeat-until loop
10    continue
      	myMod = mod(remainingNum, 10)
      	addNum = char(myMod + 48)
      	output = addNum // output
      	remainingNum = remainingNum / 10
      if (remainingNum .gt. 0) goto 10
c		end repeat-until loop
      
      return
      end
