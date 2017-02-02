      MODULE params
      PARAMETER(nbodsmax=1000,ncells=2*nbodsmax,nbodcell=nbodsmax+
     &ncells,ndim=2,nsubcell=2**ndim)
      END MODULE params


        MODULE rootcoords !
        USE params
        SAVE
        REAL(KIND=8), dimension(:,:), allocatable :: bottom
        REAL(KIND=8), dimension(:), allocatable :: cellsize
        INTEGER(KIND=4), dimension(:,:), allocatable :: pointers_of_tree
        INTEGER(KIND=4) levbit
        INTEGER(KIND=4) root,incells,itercell_tree
        REAL(KIND=8) rsize,rmin(3)
        REAL(KIND=8) pm1(0:nsubcell-1,ndim)
        END MODULE rootcoords


       MODULE bodies
       SAVE
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: pos
      INTEGER(KIND=8), dimension(:), allocatable :: key_M
      INTEGER(KIND=4) , dimension(:,:), allocatable :: iback
      INTEGER(KIND=4) , dimension(:), allocatable :: subindex,bodlist
      INTEGER(KIND=4) nbodies
       END MODULE bodies

       MODULE treemod
       USE rootcoords 
       USE bodies
       END MODULE treemod


      program oct_tree
      USE treemod
      REAL xcell(4),ycell(4),ptop(2),pbottom(2)
      REAL(KIND=8)  conv_to_int
      INTEGER pcell,x_key,y_key
      INTEGER(KIND=8) morton_key
      
      OPEN(12,file='tree.dat',form='FORMATTED')
      rewind 12
      read(12,*) bxl,nbodies

      IF(nbodies>nbodsmax) THEN
      write(6,*)' nbodies nbodsmax',nbodies,nbodsmax
      STOP
      END IF

 
      ALLOCATE(pos(2,nbodcell),bottom(2,nbodcell),subindex(nbodies),
     &bodlist(nbodies),cellsize(nbodcell),iback(2,nbodcell))

       DO i=1,nbodies
       read(12,*) pos(1,i),pos(2,i)
       END DO
      CLOSE(12)

        rmin=0
        rsize=bxl

C   subcell positions with respect the parent cell. The ordering
C  is : bottom left , bottom right , top left and right 
C   beware to keep the same ordering when placing the bits in 
C  Morton key
         pm1(0,1)=-1
         pm1(0,2)=-1

         pm1(1,1)=+1
         pm1(1,2)=-1

         pm1(2,1)=-1
         pm1(2,2)=+1

         pm1(3,1)=+1
         pm1(3,2)=+1


        levbit=10  ! order of the Morton key
        conv_to_int=2**levbit/bxl
        ALLOCATE(key_M(nbodies))

C  compute the key-values of the points

        DO i=1,nbodies
        x_key=conv_to_int*pos(1,i)
        y_key=conv_to_int*pos(2,i)
        key_M(i)=morton_key(x_key,y_key,levbit)
        END DO

        incells=1
        root=nbodsmax+1
        ALLOCATE(pointers_of_tree(nsubcell,root:nbodcell))

        pointers_of_tree=0
        cellsize=0
        bottom=0


C   write pos of root cell
        DO k=1,ndim
           pos(k,root)=rmin(k)+0.5*rsize
           bottom(k,root)=rmin(k)
        END DO
C  set the size of the root cell
        cellsize(root)=rsize

C initialize the list of bodies
        DO i=1,nbodies
           bodlist(i)=i
        END DO

        nbodlist=nbodies


        CALL tree_sort(nbodlist,bodlist,root)

C  now pointers_of_tree has the following memory layout , starting from the
C  root cell : pointers_of_tree(j=1,4:root)= ic1 , ic2, ic3 ic4 ==icnext
C where : if icnext >=root ==> new cell with more than one particles
C   if icnext =0 ==> empty , if 0< icnext< root , icnext= particle index
C with coordinates in the sub box defined by (j,ic) with
C  pointers_of_tree(j:ic)=icnext


      OPEN(12,file='treecell.dat',form='formatted')
      rewind 12
      write(12,*) incells
      DO ic=1,incells
      pcell=ic+root-1
      xcell(1)=bottom(1,pcell)
      ycell(1)=bottom(2,pcell)

      xcell(2)=bottom(1,pcell)+cellsize(pcell)
      ycell(2)=bottom(2,pcell)

      xcell(3)=bottom(1,pcell)+cellsize(pcell)
      ycell(3)=bottom(2,pcell)+cellsize(pcell)

      xcell(4)=bottom(1,pcell)
      ycell(4)=bottom(2,pcell)+cellsize(pcell)

      write(12,900) (xcell(m),ycell(m),m=1,4)
      END DO
      CLOSE(12)

      OPEN(12,file='treecell_part.dat',form='formatted')
      rewind 12
      DO i=1,nbodies
              j=iback(1,i)
              ic=iback(2,i)
              sizec=cellsize(ic)/2
C  \vec x_c=\vec bottom(cell)+ L/4
C  \vec x_c=pos(newcell) = \vec X +\vec pm1 L/4
              DO m=1,2
C   this the smallest cell containing particle  i
              pbottom(m)=bottom(m,ic)+0.25*(pm1(j-1,m)+1)*cellsize(ic)
              ptop(m)=pbottom(m)+sizec
              END DO

      xcell(1)=pbottom(1)
      ycell(1)=pbottom(2)

      xcell(2)=ptop(1)
      ycell(2)=pbottom(2)

      xcell(3)=ptop(1)
      ycell(3)=ptop(2)

      xcell(4)=pbottom(1)
      ycell(4)=ptop(2)

      write(12,900) (xcell(m),ycell(m),m=1,4)
      END DO
      CLOSE(12)


900   FORMAT(1x,4(1x,2(f6.2,1x)))



      CALL SORTI(key_M,subindex,nbodies)


      OPEN(12,file='tree_ord.dat',form='formatted')
      rewind 12

      write(12,*) bxl
      DO m=1,nbodies
      i=subindex(m)
      write(12,*) pos(1,i),pos(2,i),key_M(i)
      END DO
      CLOSE(12)



      STOP
      END



C=======================================================================

C
         RECURSIVE SUBROUTINE tree_sort(nblist,listbodies,oldcell)


C =======================================================================
C
C
        USE treemod
        INTEGER(KIND=4) listbodies(nblist),hoc(0:3)
        INTEGER(KIND=4) sublist(nblist),llj(nblist)
 
        INTEGER(KIND=4) p,pcell,oldcell,newcell,pbody,lpos,jsub
        INTEGER(KIND=8) key_of_point

 
        itercell_tree=itercell_tree+1
        IF(itercell_tree>levbit) THEN
C  you can not go deeper than the order of the Morton key
        write(6,*)' itercell level',itercell_tree,levbit
        STOP
        END IF


C  remember: the most significant bits start from the right
        levscan=levbit-itercell_tree
 
C  reset the local linked-list
        HOC=0
        LLJ=0


C   find the number of particles in each sub-cell
           DO k=1,nblist
           i=listbodies(k)
           key_of_point=key_M(i)
           lpos=ndim*levscan
           ju=IBITS(key_of_point,lpos,2)  ! sub cell index
           kprev=HOC(ju)
           LLJ(k)=kprev
           HOC(ju)=k
           END DO
 
           DO jsub=0,nsubcell-1  ! j=sub-cell index
 

           k=HOC(jsub)
           IF(k==0) CYCLE   ! sub-cell empty

           nsubc=0
           DO WHILE(k>0)  ! prepare the list of particles
           nsubc=nsubc+1
           i=listbodies(k)
           sublist(nsubc)=i
           k=llj(k)
           END DO

           IF(nsubc>1) THEN
           incells=incells+1
           IF(incells>ncells) THEN
           write(6,*)' incells',incells,ncells
           STOP
           END IF


           newcell=incells+root-1
           pointers_of_tree(jsub+1,oldcell)=newcell
           cellsize(newcell)=cellsize(oldcell)*0.5
C L=size(oldcell)
C   \vec X = pos(oldcell) 
C  \vec x_c=pos(newcell) = \vec X +\vec pm1 L/4
C   \vec bottom(newcell)= \vec x_c-L/4
           DO m=1,ndim
           pos(m,newcell)=pos(m,oldcell)+pm1(jsub,m)*0.5*
     &     cellsize(newcell)
           bottom(m,newcell)=pos(m,newcell)-0.5*cellsize(newcell)
           END DO
           iback(1,newcell)=jsub+1
           iback(2,newcell)=oldcell

           CALL tree_sort(nsubc,sublist,newcell)

           ELSE   ! one particle

           pbody=sublist(nsubc)
           pointers_of_tree(jsub+1,oldcell)=pbody
           iback(1,pbody)=jsub+1
           iback(2,pbody)=oldcell


           END IF

 
        END DO   ! end loop sub-cells


        itercell_tree=itercell_tree-1

        RETURN
        END



! This function computes a Morton key for an integer triplet (x,y,z),
! with x,y,z in the range between 0 and 2^bits-1.
    
      integer*8 function morton_key( x, y, bits)
C      integer*8 function morton_key( x, y, z, bits)
      integer*4, intent(in):: x, y, bits

      integer*4:: i,levscan
      integer*8:: key,kc1
      logical bitx,bity,bitz



            levscan=bits-1
            key=0
            DO WHILE(levscan>=0)

            bitx=BTEST(x,levscan)
            bity=BTEST(y,levscan)
C            bitz=BTEST(z,levscan)

            IF(bitx) key=IBSET(key,0+2*levscan)
            IF(bity) key=IBSET(key,1+2*levscan)
C            IF(bitz) key=IBSET(key,2+3*levscan)  ! 3D--> 3*levscan
            levscan=levscan-1
            END DO
    
        morton_key = key
    
       end function morton_key

 


        SUBROUTINE SORTI(X,IND,N)
C  Numerical recipe sorting routine 
        DIMENSION IND(N)
        INTEGER(KIND=8)  X(N),Q
        DO 11 J=1,N
        IND(J)=J
11      CONTINUE
        L=N/2+1
        IR=N
10      CONTINUE
        IF(L.GT.1) THEN
        L=L-1
        INDXT=IND(L)
        Q=X(INDXT)
        ELSE
        INDXT=IND(IR)
        Q=X(INDXT)
        IND(IR)=IND(1)
        IR=IR-1
        IF(IR.EQ.1) THEN
        IND(1)=INDXT
        RETURN
        END IF
        END IF
        I=L
        J=L+L
20      IF(J.LE.IR) THEN
        IF(J.LT.IR) THEN
        IF(X(IND(J)).LT.X(IND(J+1))) J=J+1
        END IF
        IF(Q.LT.X(IND(J))) THEN
        IND(I)=IND(J)
        I=J
        J=J+J
        ELSE
        J=IR+1
        END IF
        GO TO 20
        END IF
        IND(I)=INDXT
        GO TO 10
        END


