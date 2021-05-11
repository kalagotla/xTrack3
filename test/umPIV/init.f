	program hexa
c
c	hexahedra test program for Visual3
c
c	reads a structured (hexahedra) grid and flow data;
c	  creates the unstructured Visual3 information
c
c	constructs a programmed cutting plane
c
        parameter (nkeys = 13)
c
        integer      keytyp(nkeys), ikeys(nkeys)
        real         flims(2,nkeys)
        character*16 tkeys(nkeys)
	character*32 infile
	character*80 data(5), ttl, titl
	logical      win3d
c
	parameter (mnod = 1e7)
	common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
	data ttl/'        Hexahedra Test of Visual3'/
c
c       key bindings
        data tkeys / 'Density         ',
     &               'Pressure        ',
     &               'Mach number     ',
     &               'Stag. enthalpy  ',
     &               'Stag. pressure  ',
     &               'Entropy         ',
     &               'Temperature     ',
     &               'Flow Vectors    ',
     &               '2nd Flow Vectors',
     &               'Surface Pressure',
     &               'Surface Mach #  ',
     &               'Surface Flow    ',
     &               'Radius          '/
        data keytyp/ 7*1, 2*2, 2*3, 4, 5/
        data ikeys/100,112,109,104,120,115,116,102,70,80,77,83,114/
        data flims/ 0.0, 0.2,  0.3, 0.9,  0.1, 1.1,  -2.0, 2.0,
     &     -2.0, 2.0,  -1.0, 1.0,  0.6, 1.0,  0.20,  0.0, 1.0, 0.0,
     &     0.4, 0.7, 0.1, 1.0, 0.10, 0.0, 0.0, 0.0/
c
	write(*,*) ' '
	write(*,*) ttl
c
c------ Read mesh (structured)
c
        infile = 'gridData.o'
        open(unit=11, file=infile,status='old',form='unformatted')
c
        read(11) il, jl, kl
	knode = il*jl*kl
	if(knode .gt. mnod) stop 'too many nodes!'
        read(11) (xgrd(k), ygrd(k), zgrd(k), k=1, knode)
        close(11)
c
c---- Read flow file (structured)
c
        infile = 'flowData.o'
        open(unit=11, file=infile,status='old',form='unformatted')
 		read(11) gam, R
        read(11)  ((q(n,k),n=1,5), k=1, knode)
        close(11)
c
c-----  set up Visual3 constants
c
	kcel1  = 0
	kcel2  = 0
	kcel3  = 0
	kcel4  = (il-1)*(jl-1)*(kl-1)
	ksurf  = 2*(il-1)*(jl-1) + 2*(il-1)*(kl-1) + 2*(jl-1)*(kl-1)
	knsurf = 6
c
	write(*,*) ' '
	write(*,*) 'Enter WIN3D (t or f):'
	read(*,*) win3d
	titl = ttl
	if(.not. win3d) titl = ttl(9:80)
c
        call V3_Init(titl, 0, 'data/spec.col', 99, win3d,
     &               nkeys, ikeys, tkeys, keytyp, flims,
     &               0, knode, 0, kcel1, kcel2, kcel3, kcel4,
     &               0, 0, 0, 0,  ksurf, knsurf)
c
        call v3zprime(.true., %val(pxyz), knode3, %val(pu),
     &                      zprime, xpc, ypc, halfw)
	stop
	end



	subroutine V3Cell(cel1, cel2, cel3, cel4, nptet, ptet)
c
c	routine that passes the cell info to Visual3
c
	integer cel1(4,1),  cel2(5,1), cel3(6,1), cel4(8,1)
	integer nptet(8,1), ptet(1)
c
c	where:	cel1	- node numbers for tetrahedral cells
c			  (filled to kcel31)
c		cel2	- node numbers for pyramid cells
c			  (filled to kcel32)
c		cel3	- node numbers for  prism cells
c			  (filled to kcel33)
c		cel4	- node numbers for hexahedral cells
c			  (filled to kcel34)
c		nptet	- poly-tetrahedra header information:
c			  fill first 4 entries for each strip.
c			  the first entry points to the last position
c			  in ptet for the end of the strip, 
c			  the next entries are filled with the 
c			  first 3 nodes of the strip.
c		ptet	- the rest of the poly-tetra, 1 node per
c			  additional cell
c
c	Note: if KCELn is zero the corresponding CELn must NOT be
c	      filled!
c
	parameter (mnod = 1e7)
	common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
c	fill in cel4 with node numbers
c
	kc = 0
	do 1 k = 1, kl-1
	  do 11 j = 1, jl-1
	    do 111 i = 1, il-1
	      kn = i + (j-1)*il + (k-1)*il*jl
	      kc = kc + 1
	      cel4(1,kc) = kn
	      cel4(2,kc) = kn + 1
	      cel4(3,kc) = kn + il + 1
	      cel4(4,kc) = kn + il
	      cel4(5,kc) = kn + il*jl
	      cel4(6,kc) = kn + il*jl + 1
	      cel4(7,kc) = kn + il*jl + il + 1
	      cel4(8,kc) = kn + il*jl + il
 111	    continue
 11	  continue
 1	continue
c
	return
	end


	subroutine V3Surface(nsurf, surf, scel, tsurf)
c
c	routine that passes the surface info to Visual3
c
c	where:	nsurf	- pointers to the end of each surface group and
c			  flag to indicate if surface is on
c		surf	- pointers to cells+faces that make the surfaces
c		scel    - node numbers for surfaces
c		tsurf	- surface group title (optional)
c
	integer      nsurf(2,*), surf(*), scel(4,*)
	character*20 tsurf(*)
c
	parameter (mnod = 1e7)
	common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
c	fill in domain surface data
c
	ksurf = 0
c
c	surface 1 - Top
	do 2 k = 1, kl-1
	  do 21 j = 1, jl-1
	    ksurf = ksurf + 1
	    kn = 1 + (j-1)*il + (k-1)*il*jl
	    scel(1,ksurf) = kn
	    scel(2,ksurf) = kn + il
	    scel(3,ksurf) = kn + il + il*jl
	    scel(4,ksurf) = kn + il*jl
 21	  continue
 2	continue
	nsurf(1,1) = ksurf
c	print*, scel(*,3)
c	translucent
	nsurf(2,1) = -1
	tsurf(1) = 'XY2'
c
c	surface 2 - Bottom
	do 3 k = 1, kl-1
	  do 31 j = 1, jl-1
	    ksurf = ksurf + 1
	    kn = il + (j-1)*il + (k-1)*il*jl
	    scel(1,ksurf) = kn
	    scel(2,ksurf) = kn + il
	    scel(3,ksurf) = kn + il + il*jl
	    scel(4,ksurf) = kn + il*jl
 31	  continue
 3	continue
	nsurf(1,2) = ksurf
c	translucent
	nsurf(2,2) = -1
	tsurf(2) = 'XY1'
c
cc	surface 3 - hub
	do 4 j = 1, jl-1
	  do 41 i = 1, il-1
	    ksurf = ksurf + 1
	    kn = i + (j-1)*il 
	    scel(1,ksurf) = kn
	    scel(2,ksurf) = kn + 1
	    scel(3,ksurf) = kn + 1 + il
	    scel(4,ksurf) = kn + il
 41	  continue
 4	continue
	nsurf(1,3) = ksurf
c	translucent
	nsurf(2,3) = -1
	tsurf(3) = 'YZ1'
c
c	surface 4 - tip
	do 5 j = 1, jl-1
	  do 51 i = 1, il-1
	    ksurf = ksurf + 1
	    kn = i + (j-1)*il + (kl-1)*il*jl
	    scel(1,ksurf) = kn
	    scel(2,ksurf) = kn + 1
	    scel(3,ksurf) = kn + 1 + il
	    scel(4,ksurf) = kn + il
 51	  continue
 5	continue
	nsurf(1,4) = ksurf
c	translucent
	nsurf(2,4) = -1
	tsurf(4) = 'YZ2'

c	surface 5 - blade
	do 6 k = 1, kl-1
	  do 61 i = 1, il-1
	    ksurf = ksurf + 1
	    kn = i + (k-1)*il*jl
	    scel(1,ksurf) = kn
	    scel(2,ksurf) = kn + 1
	    scel(3,ksurf) = kn + 1 + il*jl
	    scel(4,ksurf) = kn + il*jl
 61	  continue
 6	continue
	nsurf(1,5) = ksurf
c	print*, nsurf(1,5)
c	translucent
	nsurf(2,5) = -1
	tsurf(5) = 'XZ1'
c
c	surface 6 - blade
	do 7 k = 1, kl-1
	  do 71 i = 1, il-1
	    ksurf = ksurf + 1
	    kn = i + (jl-1)*il + (k-1)*il*jl
	    scel(1,ksurf) = kn
	    scel(2,ksurf) = kn + 1
	    scel(3,ksurf) = kn + 1 + il*jl
	    scel(4,ksurf) = kn + il*jl
 71	  continue
 7	continue
	nsurf(1,6) = ksurf
c	translucent
	nsurf(2,6) = -1
	tsurf(6) = 'XZ2'
c
c
	return
	end

	subroutine V3Grid(xyz)
c
c	pass the coordinate triads to Visual3
c
	real xyz(3,1)
c
	parameter (mnod = 1e7)
	common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
c	NOTE: averages and deviations calucualted only for surface mapping.
c
	common /sprime/ itype, xave, yave, zave, xdev, ydev, zdev
c
c	fill xyz values
	knode = il*jl*kl
	xave = 0.0
	yave = 0.0
	zave = 0.0
	do 1 k = 1, knode
	  xyz(1,k) = xgrd(k)
	  xyz(2,k) = ygrd(k)
	  xyz(3,k) = zgrd(k)
	  xave = xave + xyz(1,k)
	  yave = yave + xyz(2,k)
	  zave = zave + xyz(3,k)
 1	continue
c	print *, (xyz(3,k), k=1,knode)
	xave = xave/float(knode)
	yave = yave/float(knode)
	zave = zave/float(knode)
c
	xdev = 0.0
	ydev = 0.0
	zdev = 0.0
c
	do 2 k = 1, knode
	  xdev = xdev + abs(xyz(1,k)-xave)
	  ydev = ydev + abs(xyz(2,k)-yave)
	  zdev = zdev + abs(xyz(3,k)-zave)
 2	continue
	xdev = xdev/float(knode)
	ydev = ydev/float(knode)
	zdev = zdev/float(knode)
c
	return
	end


	subroutine V3Scal(key, v)
c
c	return the field scalar in "v" based on the key hit
c
	dimension v(1)
c
	parameter (mnod = 1e7)
	common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
	knode = il*jl*kl
        gm1 = gam - 1.
        go to (1, 2, 3, 4, 5, 6, 7), key
c
c       density
 1      do 11 k = 1, knode
          v(k) = q(1,k)
 11     continue
        return
c
c       pressure
 2      do 21 k = 1, knode
          v(k) = q(1,k)*R*q(5,k)
 21     continue
        return
c
c       mach number
 3      do 31 k = 1, knode
          v2 = (q(2,k)*q(2,k) + q(3,k)*q(3,k) + q(4,k)*q(4,k))
          p = gm*R*q(5,k)
          v(k) = sqrt(v2/p)
 31     continue
        return
c
c       stagnation enthalapy
 4      do 41 k = 1, knode
          rr = 1./q(1,k)
          v2 = (q(2,k)*q(2,k) + q(3,k)*q(3,k) + q(4,k)*q(4,k))*rr*rr
          h0 = gam*(q(5,k)*rr - 0.5*v2) + 0.5*v2
          v(k) = 100.0 * (gm1*h0 - 1.0)
 41     continue
        return
c
c       stagnation pressure
 5      do 51 k = 1, knode
          v2 = (q(2,k)*q(2,k) + q(3,k)*q(3,k) + q(4,k)*q(4,k))/q(1,k)
          p = gm1*(q(5,k) - 0.5*v2)
          xmach = sqrt(v2/(gam*p))
          p0 = p*(1.0+0.5*gm1*xmach*xmach)**(gam/gm1)
          v(k) = 100.0 * (gam*p0 - 1.0)
 51     continue
        return
c
c       entropy
 6      do 61 k = 1, knode
          v2 = (q(2,k)*q(2,k) + q(3,k)*q(3,k) + q(4,k)*q(4,k))/q(1,k)
          p = gm1*(q(5,k) - 0.5*v2)
          v(k) = 100.0 * (alog(gam*p) - gam*alog(q(1,k)))
 61     continue
        return
c
c       temperature
 7      do 71 k = 1, knode
          v(k) = q(5,k)
 71     continue
        return
c
        end


        subroutine V3Vect(key, v)
c
c	return the field vector in "v" based on the key hit
c
        dimension v(3,1)
c
        parameter (mnod = 1e7)
        common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
        knode = il*jl*kl
c
	if(key .eq. 8) then
          do 1 k = 1, knode
            v(1,k) = q(2,k)
            v(2,k) = q(3,k)
            v(3,k) = q(4,k)
 1        continue
	endif
c
	if(key .eq. 9) then
	  do 2 i = 1, il
	    vs1 = 0.
	    vs2 = 0.
	    vs3 = 0.
	    do 21 j = 1, jl
	      jx = (j-1)*il
	      do 211 k = 1, kl
	        kn = i + jx + (k-1)*il*jl
	        vs1 = vs1 + q(2,kn)
	        vs2 = vs2 + q(3,kn)
	        vs3 = vs3 + q(4,kn)
 211	      continue
 21	    continue
	    vs1 = vs1/float(jl*kl)
	    vs2 = vs2/float(jl*kl)
	    vs3 = vs3/float(jl*kl)
	    do 22 j = 1, jl
	      jx = (j-1)*il
	      do 221 k = 1, kl
	        kn = i + jx + (k-1)*il*jl
	        v(1,kn) = q(2,kn) - vs1
	        v(2,kn) = q(3,kn) - vs2
	        v(3,kn) = q(4,kn) - vs3
 221	      continue
 22	    continue
 2        continue
	endif
	return
	end



	subroutine v3thres(key, xyz, t)
c
c	return scalar thresholding values for the 3D field
c
c	where:	key	- index for the selected scalar function
c               xyz     - coordinate triads for all the 3-D nodes.
c                         this is the same data as set in V3GRID
c
c		note: key and xyz must not be modified!
c
c		the following paramters must be returned by this routine:
c
c		t	- the threshold value for all the 3D nodes
c
	integer key
	real    xyz(3,1), t(1)
c
	parameter (mnod = 1e7)
	common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
	knode = il*jl*kl
c
	do 1 i = 1, knode
	  t(i) =  sqrt(xyz(3,i)*xyz(3,i) + xyz(2,i)*xyz(2,i))
 1	continue
c
	return
	end



c ********************************************************************
c	routines for mapping domain surfaces
c ********************************************************************


	subroutine v3surf(isurf, xpc, ypc, halfw)
c
c       return surface rendering data
c
c               this routine is called when the 2-D surface
c               function-key is hit to set up the 2-D data.
c               the operator may be queried in this routine.
c
c       where:  isurf   - the current active surface
c
c               note: isurf must not be modified!
c
c               the following paramters are returned by this routine:
c
c               xpc     - the starting x-prime center value for
c                         the 2-D window
c               ypc     - the starting y-prime center value for
c                         the 2-D window
c               halfw   - the starting 1/2 width for the 2-D window in
c                         x-prime/y-prime space.  if halfw is not filled
c                         the requested surface is not plotted.
c
	integer isurf
	real    xpc, ypc, halfw
c
c       common blocks may be used to communicate data between this
c       routine and others (including the main program)
c
	common /sprime/ itype, xave, yave, zave, xdev, ydev, zdev
c
	if(isurf .lt. 3 .or. isurf .gt. 6) return
	itype = 3
	if(isurf .gt. 4) itype = 2
c
	if(itype .eq. 2) then
	  xpc    = xave
	  ypc    = zave
	  halfw  = xdev + zdev
	endif
c
	if(itype .eq. 3) then
	  xpc    = xave
	  ypc    = yave
	  halfw  = xdev + ydev
	endif
c
	return
	end



	subroutine v3xysurf(kn, xyz, n, xyp)
c
c       return x-prime and y-prime for the surface
c
c               this routine is called during the data collection phase
c               of the mapped domain surfaces.  this routine may be
c               called many times during this phase.
c               there is no ordering of the nodes in kn, infact kn
c               may contain the same node more than once in a single call.
c
c       where:  kn      - a vector of indices for the nodes to be
c                         transformed
c               xyz     - coordinate triads for 3-D nodes.
c                         this is the same data as set in V3GRID
c               n       - the number of nodes to be transformed.
c
c               note: kn, xyz and n must not be modified!
c
c               the following paramters must be returned by this routine:
c
c               xyp     - the calculated x-prime and y-prime values for
c                         the 3-D nodes
c
	integer kn(n)
	real    xyz(3,1), xyp(2,n)
c
c	common blocks may be used to communicate data between this
c	routine and others (including the main program)
c
	common /sprime/ itype, xave, yave, zave, xdev, ydev, zdev
c
	print*, 'v3xysurf called'
	if(itype .eq. 2) then
	  do 1 i = 1, n
	    k = kn(i)
	    xyp(1,i) = xyz(1,k)
	    xyp(2,i) = xyz(3,k)
 1	  continue
	endif
c
	if(itype .eq. 3) then
	  do 2 i = 1, n
	    k = kn(i)
	    xyp(1,i) = xyz(1,k)
	    xyp(2,i) = xyz(2,k)
 2	  continue
	endif
c
	return
	end


	subroutine v3ssurf(key, kn, xyz, n, s)
c
c	return scalar values for the surface
c
c		this routine is called during the data collection phase
c		of the mapped domain surfaces.  this routine may be
c		called many times during this phase.
c		there is no ordering of the nodes in kn, infact kn
c		may contain the same node more than once in a single call.
c
c	where:	key	- index for the current surface function
c		kn	- a vector of indices for the nodes to be filled
c		xyz	- coordinate triads for 3-D nodes.
c			  this is the same data as set in V3GRID
c		n	- the number of nodes to be transformed.
c
c		note: key, kn, xyz and n must not be modified!
c
c		the following paramters must be returned by this routine:
c
c		s	- the functional value for the nodes based on key
c
	integer key, kn(n)
	real    xyz(3,1), s(n)
c
c	common blocks may be used to communicate data between this
c	routine and others (including the main program)
c
	parameter (mnod = 1e7)
	common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
        gm1 = gam - 1.
	if(key .eq. 10) then
	  do 1 i = 1, n
	    k = kn(i)
            s(i) = gm*R*q(5,k)
 1	  continue
	else
	  do 2 i = 1, n
	    k = kn(i)
            v2 = (q(2,k)*q(2,k) + q(3,k)*q(3,k) + q(4,k)*q(4,k))
            p = gm*R*q(5,k)
            s(i) = sqrt(v2/p)
 2	  continue
	endif
c
	return
	end


	subroutine v3vsurf(key, kn, xyz, n, v)
c
c	return vector values for the surface
c
c		this routine is called during the data collection phase
c		of the mapped domain surfaces.  this routine may be
c		called many times during this phase.
c		there is no ordering of the nodes in kn, infact kn
c		may contain the same node more than once in a single call.
c
c	where:	key	- index for the current surface function
c		kn	- a vector of indices for the nodes to be filled
c		xyz	- coordinate triads for 3-D nodes.
c			  this is the same data as set in V3GRID
c		n	- the number of nodes to be transformed.
c
c		note: key, kn, xyz and n must not be modified!
c
c		the following paramters must be returned by this routine:
c
c		v	- the vector functional value for the nodes.
c
	integer key, kn(n)
	real    xyz(3,1), v(3,n)
c
c	common blocks may be used to communicate data between this
c	routine and others (including the main program)
c
	parameter (mnod = 1e7)
	common /grid/ il, jl, kl, ile, ite, gam, R, q(5,mnod),
     &                xgrd(mnod), ygrd(mnod), zgrd(mnod)
c
	do 1 i = 1, n
	  k = kn(i)
          v(1,i) = q(2,k)
          v(2,i) = q(3,k)
          v(3,i) = q(4,k)
 1	continue
c
	return
	end



c ********************************************************************
c	Routines for Programmed Cut Planes
c ********************************************************************



	subroutine v3zprime(flag, xyz, knode, zp, zprime, xpc, ypc, halfw)
c
c	return z-prime for programmer defined cutting surfaces
c
c		this routine is called when the programmed cutting plane
c		function-key is hit to set up the surface data.
c		the operator may be queried in this routine.
c	        Note: this is call is NOT locked out in unsteady usage!
c
c	where:	flag	- .true. first call
c			  .false. grid change (unsteady - iopt > 1)
c			  only ask user questions when flag = .true.
c		xyz	- coordinate triads for all the 3-D nodes.
c			  this is the same data as set in V3GRID
c		knode   - the number of 3-D nodes.   this is the same
c			  value as set in the call to V3_INIT
c
c		note: flag, xyz and knode must not be modified!
c
c		the following paramters are returned by this routine:
c
c		zp	- the calculated z-prime values for the 3-D nodes.
c			  this vector must be completely filled
c		zprime	- the starting z-prime value
c		xpc	- the starting x-prime center value for 
c			  the 2-D window
c		ypc	- the starting y-prime center value for 
c			  the 2-D window
c		halfw	- the starting 1/2 width for the 2-D window in
c			  x-prime/y-prime space.
c
	logical flag
	integer knode
	real    xyz(3,knode), zp(knode), zprime, xpc, ypc, halfw
c
c	common blocks may be used to communicate data between this
c	routine and others (including the main program)
c
	common /prime/ itype
c
c	R-Theta Cuts
c
c       To automate some process. Multiple lines were removed
	if(flag) then
	  itype = 2
c
	  xave = 0.0
	  yave = 0.0
	  zave = 0.0
	  do 1 i = 1, knode
	    xave = xave + xyz(1,i)
	    yave = yave + xyz(2,i)
	    zave = zave + xyz(3,i)
 1	  continue
	  xave = xave/float(knode)
	  yave = yave/float(knode)
	  zave = zave/float(knode)
c
	  xdev = 0.0
	  ydev = 0.0
	  zdev = 0.0
	  do 2 i = 1, knode
	    xdev = xdev + abs(xyz(1,i)-xave)
	    ydev = ydev + abs(xyz(2,i)-yave)
	    zdev = zdev + abs(xyz(3,i)-zave)
 2	  continue
	endif
c
c
	if(itype .eq. 2) then
	  do 4 i = 1, knode
	    zp(i) = xyz(2,i)
 4	  continue
	  if(flag) then
	    zprime = 0.409091*1e-3
	    xpc    = zave
	    ypc    = xave
	    halfw  = 2*(xdev+zdev)/float(knode)
	  endif
	endif
c
c
	return
	end


	subroutine v3xyprime(zprime, kn, xyz, n, xyp)
c
c       return x-prime and y-prime for programmer defined cutting surface
c
c               this routine is called during the data-collection phase
c               of the program defined surfaces.  this routine may be
c               called many times during this phase.
c               there is no ordering of the nodes in kn, infact kn
c               may contain the same node more than once in a single call.
c
c       where:  zprime  - the current z-prime value
c               kn      - a vector of indices for the nodes to be
c                         transformed
c               xyz     - coordinate triads for 3-D nodes.
c                         this is the same data as set in V3GRID
c               n       - the number of nodes to be transformed.
c
c               note: zprime, kn, xyz and n must not be modified!
c
c               the following paramter is returned by this routine:
c
c               xyp     - the calculated x-prime and y-prime values for 
c                         the 3-D nodes
c
        integer kn(n)
        real    zprime, xyz(3,1), xyp(2,n)
c
c       common blocks may be used to communicate data between this
c       routine and others (including the main program)
c
        common /prime/ itype
c
c       R-Theta Cuts
c
        if(itype .eq. 1) then
          do 1 i = 1, n
            k = kn(i)
            xyp(1,i) = xyz(1,k)
            xyp(2,i) = zprime*atan2(xyz(2,k),xyz(3,k))
 1        continue
        endif
c
        if(itype .eq. 2) then
          do 2 i = 1, n
            k = kn(i)
            xyp(1,i) = sqrt(xyz(3,k)*xyz(3,k) + xyz(2,k)*xyz(2,k))
            xyp(2,i) = xyz(1,k)
 2        continue
        endif
c
        if(itype .eq. 3) then
          do 3 i = 1, n
            k = kn(i)
            xyp(1,i) = xyz(3,k)
            xyp(2,i) = xyz(2,k)
 3        continue
        endif
        return
	end
