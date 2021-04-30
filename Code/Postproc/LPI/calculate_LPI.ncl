  load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
external CALC_LPI "./calclpi.so"
begin
   tn = 0
   files := systemfunc("ls wrfout_d02*")
;   ListCount = dimsizes(files) ; barry
;   do nf =0,ListCount(files)-1
;    print("nf = " + nf)  ; barry
;    a := addfile("wrfout_d01_2015-06-01_01:00:00","r")
   a= addfiles(files,"r")
;   do nf = 0,ListCount(a)-1
;    print("nf = " + nf)
    times  = wrf_user_list_times(a)  ; get times in the file
    ntimes = dimsizes(times)         ; number of times in the file
    print ("ntimes = " + ntimes)
    do time = 0,ntimes-1
     xlat := wrf_user_getvar(a,"XLAT",time)
     xlon := wrf_user_getvar(a,"XLONG",time)
     ph := wrf_user_getvar(a,"PH",time)
     phb := wrf_user_getvar(a,"PHB",time)
     p := wrf_user_getvar(a,"P",time)
     pb := wrf_user_getvar(a,"PB",time)
     lpi := wrf_user_getvar(a,"RAINNC",time) ; for regular call to "2d" lpi, this needs to be 2d
     ph := wrf_user_getvar(a,"PH",time)
     phb := wrf_user_getvar(a,"PHB",time)
     p := wrf_user_getvar(a,"P",time)
     pb := wrf_user_getvar(a,"PB",time)
     qc := wrf_user_getvar(a,"QCLOUD",time)
     qr := wrf_user_getvar(a,"QRAIN",time)
     qi := wrf_user_getvar(a,"QICE",time)
     qs := wrf_user_getvar(a,"QSNOW",time)
     qg := wrf_user_getvar(a,"QGRAUP",time)
     w := wrf_user_getvar(a,"W",time)
     t := wrf_user_getvar(a,"T",time)
     u_a := wrf_user_getvar(a,"ua",time)
     v_a := wrf_user_getvar(a,"va",time)
     dims3d := dimsizes(qg)
     dims3dw := dimsizes(w)
     dims3dph := dimsizes(phb)
     ni = dims3d(2)
     nj = dims3d(1)
     nk = dims3d(0)
     nkk = dims3dph(0)
     CALC_LPI::calclpi(qc, qr, qi, qs, qg,w,ph,phb,p,pb,t,lpi,ni,nj,nk,nkk)
     file_out = "lpi_time_" + time+"_d02_2D"+".nc"
     system("rm -f file_out" )
     ncdf = addfile(file_out,"c")
     ncdf->lpi=lpi
     end do
;    end do
end
