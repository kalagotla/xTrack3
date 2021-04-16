c
c	common blocks used for Visual3
c
        integer   mtuft, mwork, mkeys, maux, mstream, mstruc, maxcol
	integer   mwork2, mkey
	parameter (mwork = 3840, mkeys = 256, maux = 320, maxcol = 1024)
	parameter (mstream = 7200, mtuft = 2500, mstruc = 100)
	parameter (mbubble = 1600, medges = 10)
	parameter (mkey = 32, mwork2 = mwork/2)
c
c...... 2-D Structures
c
c
c...... pointers
	integer*8 pface, pedge, pixy, pdtrp, pclst, pslst,
     &            pcell, icell, pl2d, p2ds, p2dv, pcedg, thrdid
	common /V2cPtr/ pface, pedge, pixy, pdtrp, pclst, pslst,
     &                  pcell, icell, pl2d, p2ds, p2dv, pcedg, thrdid
c
c...... array sizes
	integer mnode, mface, medge, mclst, mslst, mcell,
     &          knode, kface, kedge, kclst, kslst, kcell,
     &          kdnorm
	common /V2cSizes/ mnode, mface, medge, mclst, mslst, mcell,
     &                    knode, kface, kedge, kclst, kslst, kcell,
     &                    kdnorm
c
c...... flags
	integer*8 mlock, vlock, newdata
	logical render, grid, arrows, contours, stream_lines, 
     &          line_probe, batch, point_probe, strip_chart,
     &          done, edge, normal, cut_render, pause, frender,
     &          connect, aux_plot, sl_color, vcloud, streamer, refer,
     &          gotans, histogram, prog_2d, user_probe, close1, close2,
     &          edge_lines, dyn_thres, cell_render, xframe, sup_con,
     &          v3shade, tab_out
	integer ribbon, mirror, mbut(6), dialbox
	common /V2cFlags/ mlock, vlock, newdata,
     &                    render, grid, arrows, contours, cell_render,
     &                    stream_lines, line_probe, batch, frender,
     &                    point_probe, strip_chart, edge, normal,
     &                    cut_render, pause, connect, aux_plot, sup_con,
     &                    sl_color, vcloud, streamer, refer, dyn_thres,
     &                    done, gotans, histogram, prog_2d, xframe,
     &                    user_probe, close1, close2, edge_lines,
     &                    v3shade, tab_out,
     &                    ribbon, mirror, mbut, dialbox
c
c...... scaling data
	real xmin, xmax, ymin, ymax, fmin, fmax, vsize,
     &       fsmin, fsmax, vssize, cmin, cmax, s2x, s2y,
     &       xmn, xinc, ymn, yinc, tmin, tmax,
     &       xwall, ywall, xperp, yperp, fline,
     &       rsize, angr, tbmax, clincor, arata, aratb
	integer ipx1, ipy1, ipx2, ipy2, edgkey, atype, amask
	common /V2cScale/ xmin, xmax, ymin, ymax, fmin, fmax, vsize,
     &                    fsmin, fsmax, vssize, cmin, cmax, s2x, s2y,
     &                    xmn, xinc, ymn, yinc, tmin, tmax,
     &                    xwall, ywall, xperp, yperp, fline,
     &                    rsize, angr, tbmax, clincor, arata, aratb,
     &                    ipx1, ipy1, ipx2, ipy2, edgkey, atype, amask
c
c...... display data
	integer*8 display, GC33, GC3w, GC3a,
     &            rwindow, window3, window, key_window,
     &            aux_window, dial_window, text_window,
     &            GCk, GCs, GCm, GCw, GCt, GCa, GCd,
     &            cwindow, cadr, cimg
	integer colormap(maxcol), scmap(maxcol), ncol, ncolo,
     &          nscol, nscolo, white, grey, black,
     &          xtext, ytext, npitx, npity, npi3x, npi3y,
     &          npixx, npixy, npidx, npidy, npikx, npiky,
     &          npixa, num_plot, kstart, kstop, ksize,
     &          scratch, dial3d, npath, icursor(4)
	logical winopen, win32
	real rev3d(0:7), rev2d(0:7), revkey(0:7), sobj3(4,4,10), clip
	real obj3(4,4), delta, lincor, focus(3), deltaxyz, obji(4,4)
	real cut(4,4), sensit, revision, xparency, cursor(8)
	character*80 cmfile, scmfile, ccmfile, title, subttl,
     &               path, referf, answer
	common /V3cDisplay/ display, GC33, GC3w, GC3a,
     &                      rwindow, window3, window, key_window,
     &                      aux_window, dial_window, text_window,
     &                      GCk, GCs, GCm, GCw, GCt, GCa, GCd,
     &                      cwindow, cadr, cimg,
     &                      colormap, scmap, ncol, ncolo,
     &                      nscol, nscolo, white, grey, black,
     &                      xtext, ytext, npitx, npity, npi3x, npi3y,
     &                      npixx, npixy, npidx, npidy, npikx, npiky,
     &                      npixa, num_plot, kstart, kstop, ksize,
     &                      scratch, dial3d, npath, winopen, win32,
     &                      rev3d, rev2d, revkey, sobj3, obj3, sensit,
     &                      delta, focus, lincor, deltaxyz, clip, cut, 
     &                      revision, xparency, obji, cursor, icursor,
     &                      cmfile, scmfile, ccmfile, title, subttl,
     &                      path, referf, answer
c
c...... Cell based data
	logical volcell,  cutcell,  histcell,
     &          xparcell, greycell, gridcell
	integer surfcell, cellmap(maxcol), nccol, nccolo, nhist
	real    mapcell(3,maxcol), hmn, hmx
	common /V3cCell/ volcell,  cutcell,  histcell,
     &                   xparcell, greycell, gridcell,
     &                   surfcell, cellmap,  nccol, nccolo, nhist,
     &                   mapcell,  hmn,      hmx
c
c...... Shading data
	integer nshade
	real    lv(3), sv(3), ln(3), sn(3)
	common /V3cShade/ nshade, lv, sv, ln, sn
c
c...... 3-D Structures
c
c...... key bindings
	integer iopt, ikeys(mkeys), fkeys(mkeys), nkeys, kscl, kvct,
     &          kedg, kved, kthr, ksce, nk(mkey)
	real flim(2,mkeys), flmo(2,mkeys)
	character*10 tk(mkey)
	character*32 tkeys(mkeys)
	common /V3cKeys/ iopt, ikeys, fkeys, nkeys, kscl, kvct,
     &                   kedg, kved, kthr, ksce, nk,
     &                   flim, flmo,
     &                   tk, tkeys
c
c...... static structures, stream-line and tuft data
        integer*8 pstruc(8,mstruc+1), pstream(4,mstream)
	real slt(mstream), colstream(3,mstream), xyzstream(3,mstream),
     &       fstream(8,mstream), tstream(mstream), astream(mstream),
     &       dstream(mstream), timesl, dtimesl, tosl,
     &       xtuft(mtuft), ytuft(mtuft), ftuft(3,mtuft),
     &       fstruc(13,mstruc+1)
        integer kstream, cstream(mstream), nstream,
     &          nframes, ntufts, ktufts, ktuft(mtuft),
     &          nstruc
	character*20 tstruc(mstruc+1)
	common /V3cTufts/ pstruc,  pstream, slt, colstream, tosl,
     &                    kstream, cstream, xyzstream, fstream, 
     &                    tstream, astream, dstream, timesl, dtimesl,
     &                    nstream, nframes,
     &                    ktufts, ntufts, ktuft, xtuft, ytuft, ftuft,
     &                    nstruc, fstruc,
     &                    tstruc
c
c...... working storage
	real xwork(mwork), ywork(mwork), zwork(mwork), chart(maux)
	real xyzcut(3,4), sxyzcut(3,4,10), zn(3), zprime, dzp
	real xpc, ypc, halfw, angw, xyswork(2,mwork), tchart(maux)
        integer iwork(mwork), cwork(mwork), pwork(2,mwork2),
     &          twork(2,mwork2)
	integer kchart, ipmn, ipmx, ismn, ismx
	common /V2cWork/ iwork, xwork, ywork, zwork, xyswork, chart,
     &                   tchart, xyzcut, sxyzcut, zn, zprime, dzp,
     &                   xpc, ypc, halfw, angw,
     &                   cwork, pwork, twork,
     &                   kchart, ipmn, ipmx, ismn, ismx
c
c...... pointers
	integer*8 pxyz, ppsurf, psurf, pcel1, pcel2, pcel3, pcel4,
     &            pccel1, pccel2, pccel3, pccel4, piblank,
     &            ptet, pptet, pctet, pblock, pcblock, pcsurf,
     &            pc, ps, pv, psed, ppsed, pipmn, pipmx, ppsl,
     &            phash, pzsmn, pzsmx, pismn, pismx, pu, pnorm,
     &            pscel, psccel, pthres, ptsurf, pdnorm, pequiv,
     &            memt, lenC, lenI, lenR, lenP
	common /V3cPtr/ pxyz, ppsurf, psurf, pcel1, pcel2, pcel3, pcel4,
     &                  pccel1, pccel2, pccel3, pccel4, piblank,
     &                  ptet, pptet, pctet, pblock, pcblock, pcsurf,
     &                  pc, ps, pv, psed, ppsed, pipmn, pipmx, ppsl,
     &                  phash, pzsmn, pzsmx, pismn, pismx, pu, pnorm,
     &                  pscel, psccel, pthres, ptsurf, pdnorm, pequiv,
     &                  memt, lenC, lenI, lenR, lenP
c
c...... array sizes
	integer knode3, kcel31, kcel32, kcel33, kcel34,
     &          kptet, kpptet, kblock, kcblock, kcellti,
     &          ksurf, kpsurf, kcellt, knode3d, kcelld, ket,
     &          msed, kpsed, mpsed, kosed, msurft
	common / V3cSizes/ knode3, kcel31, kcel32, kcel33, kcel34,
     &                     kptet, kpptet, kblock, kcblock, kcellti,
     &                     ksurf, kpsurf, kcellt, knode3d, kcelld, ket,
     &                     msed, kpsed, mpsed, kosed, msurft
c
c...... flags
	logical cutp, set_cutp, mov_cutp, threshold,
     &          depthq, iso_surf, resort, user_cut,
     &          surface2d, tufts, show_surface, reenter,
     &          ctrl_press, lighting, render3d, render2d,
     &          XposeDial, XposeKey, XposeC
	integer surface_number, isurf, jsurf, lastfn, npage,
     &          iedge, scan, nline, ngrd, ncont, mapsl, kslcur,
     &          cnpx, cnpy, pdflag, npass, edges(medges)
	common /V3cFlags/ cutp, set_cutp, mov_cutp, threshold,
     &                    depthq, iso_surf, resort, user_cut,
     &                    surface2d, tufts, show_surface, reenter,
     &                    ctrl_press, lighting, render3d, render2d,
     &                    XposeDial, XposeKey, XposeC,
     &                    surface_number, isurf, jsurf, lastfn, npage,
     &                    iedge, scan, nline, ngrd, ncont, mapsl,kslcur,
     &                    cnpx, cnpy, pdflag, npass, edges
c
c...... bubbles
	real bubble(3,4,mbubble), vbubble(3,mbubble), tbubble(mbubble),
     &       wbubble(8,mbubble),  dbubble(4,mbubble), ugrad(3,3)
	real utime, utimel, dutime, coefs(5,4), dts(3)
	integer cbubble(3,mbubble), kpart, kbubble, ksline
	common /V3cBubbles/ nbubble, bubble, vbubble, wbubble, dbubble,
     &                      tbubble, ugrad,  utime,   utimel,  dutime,
     &                      coefs, dts,
     &                      cbubble, kpart,  kbubble, ksline
c
c......	GL stuff
	integer naux_buf, lincol
	real    vcolor(6,mwork), vertex(3,mwork), mapc(3,maxcol),
     &          maps(3,maxcol), shadec(3)
	logical snap, stereo
	common /V3cGL/ naux_buf, lincol,
     &                 vcolor, mapc, maps, shadec,
     &                 snap, stereo
	equivalence (vertex, vcolor)
c
c...... Local Stuff
	integer*8 mallocV3, reallocV3, iptr, ptr, GC3, iwin, mem
