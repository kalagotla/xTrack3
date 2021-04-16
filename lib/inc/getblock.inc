	kcs  = kcelld
	kns  = knode3d
	kb   = 1
	nlo  = 0
	nhi  = kblock+1
	flag = kblock .gt. 1
	do while (flag)
	  kb = (nlo+nhi)/2
	  if(kcr.gt.blocks(4,kb)) then
	    nlo = kb
	    if(nhi-nlo .le. 1) then
	      flag = .false.
	      kb = kb + 1
	    endif
	  else
	    nhi = kb
	    if(nhi-nlo .le. 1) flag = .false.
	  endif
	enddo
	if(kb .ne. 1) then
	  kcs = blocks(4,kb-1)
	  kns = blocks(5,kb-1)
	endif
