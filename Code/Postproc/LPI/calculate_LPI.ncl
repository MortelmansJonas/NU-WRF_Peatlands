load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
external CALC_LPI "./calclpi3d.so"

    files := systemfunc("ls -l wrfout_d01*")+".nc"
    a := addfile("wrfout_d01_2016-06-01_01:00:00","r")
    time := 0
    ph := wrf_user_getvar(a,"PH",time)
    phb := wrf_user_getvar(a,"PHB",time)
    p := wrf_user_getvar(a,"P",time)
    pb := wrf_user_getvar(a,"PB",time)
;   qc_1 := wrf_user_getvar(a,"Total_zdr",time)
    lpi := wrf_user_getvar(a,"QGRAUP",time)
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
    
;   a_o = addfile("wrfout_d01_2016-06-01_01:00:00","r")
;   i_inc = 0
;   lpi := wrf_user_getvar(a_o,"QGRAUP",i_inc)
;   ph := wrf_user_getvar(a_o,"PH",i_inc)
;   phb := wrf_user_getvar(a_o,"PHB",i_inc)
;   p := wrf_user_getvar(a_o,"P",i_inc)
;   pb := wrf_user_getvar(a_o,"PB",i_inc)
;   qc := wrf_user_getvar(a_o,"QCLOUD",i_inc)
;   qr := wrf_user_getvar(a_o,"QRAIN",i_inc)
;   qi := wrf_user_getvar(a_o,"QICE",i_inc)
;   qs := wrf_user_getvar(a_o,"QSNOW",i_inc)
;   qg := wrf_user_getvar(a_o,"QGRAUP",i_inc)
;   w := wrf_user_getvar(a_o,"W",i_inc)
;   t := wrf_user_getvar(a_o,"T",i_inc)
;   dims3d := dimsizes(qg)
;   print(dims3d)
;   dims3dw := dimsizes(w)
;   dims3dph := dimsizes(phb)
;   ni = dims3d(2)
;   nj = dims3d(1)
;   nk = dims3d(0)
;   nkk = dims3dph(0)
;   CALC_LPI::calclpi(qc, qr, qi, qs, qg,w,ph,phb,p,pb,t,lpi,ni,nj,nk,nkk)
    
    ;=====================================================================
       system("/bin/rm -f lpi.nc")   ; remove any pre-existing file
       ncdf = addfile("lpi.nc" ,"c")  ; open output netCDF file

    ;===================================================================
    ; make time an UNLIMITED dimension; recommended  for most applications
    ;===================================================================
       filedimdef(ncdf,"time",-1,True)
    ;===================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
       ncdf->LPI  = lpi                 ; 3D
