load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
external CALC_LPI "./calclpi3d.so"
   tn = 0
   files := systemfunc("ls wrfout_d01*")
   a := addfiles(files,"r")

   do nf =0,ListCount(a)-1
;    a := addfile("wrfout_d01_2015-06-01_01:00:00","r")
    time := 0
    xlat := wrf_user_getvar(a,"XLAT",time)
    xlon := wrf_user_getvar(a,"XLONG",time)
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
    
    ; This script writes a large variable to a file in NetCDF4
    ; with chunking and compression.

    print("Start writing netCDF file")
     begin

      theflnm = "./LPI.nc"
      setfileoption("nc","format","netcdf4")
      setfileoption("nc","headerReserveSpace",64000)
      setfileoption("nc","preFill",False)
     ;setfileoption("nc","defineMode",True)

     ;Remove the file if it already existed.
      if (isfilepresent(theflnm)) then
         system("rm " + theflnm)
      end if

     ;Create the file
      f = addfile(theflnm,"c")
 
     ;Define dimensions
     ctim = 1
     ntim = 250
     nlat = 140 ;for domain 1: 140 grid cells along latitude
     nlon = 200 ;for domain 1: 200 grid cells along longitude
     nvert = nk 
     dim_names = (/"times","lat","lon","vert"/)
     ;times dimension is unlimited, but currently set to 1
     dim_sizes = (/ctim,nlat,nlon,nvert/)
     dim_unlimited = (/True,False,False,False/)
    
     print("Defining dimensions sizes to be " + dim_sizes)
     filedimdef(f,dim_names,dim_sizes,dim_unlimited)
 
     ;Define chunk sizes
     mtim = 1
     mlat = 140/2
     mlon = 200/2
     mvert = nk/2
     chunk_sizes = (/mtim,mlat,mlon,mvert/)
     filechunkdimdef(f,dim_names,chunk_sizes,dim_unlimited)

     ;Define some file properties
     atts = True
     atts@description = "This file stores the LPI"
     atts@timestamp = systemfunc("date")
     print("Defining file attributes")
     fileattdef(f,atts)
 
     ;Define file variables
     print("Defining file variables")
     filevardef(f,"times","float","times")
     filevardef(f,"lat","float","lat")
     filevardef(f,"lon","float","lon")
     filevardef(f,"vert","float","vert")
     filevardef(f,"lpi","float",(/"times","lat","lon","vert"/))

     ;Define file variables attributes
     print("Defining file variable attributes")
     tatts = True
     tatts@long_name = "elapsed time"
     tatts@units = "hours since 01-01-2000"
     filevarattdef(f,"times",tatts)

     latatts = True
     latatts@lon_name = "latitude"
     latatts@units = "degrees north"
     filevarattdef(f,"lat",latatts)

     lonatts = True
     lonatts@long_name = "longitude"
     lonatts@units = "degrees_east"
     filevarattdef(f,"lon",lonatts)

     vertatts = True
     vertatts@long_name = "vertical level"
     vertatts@units = "-"
     filevarattdef(f,"vert",vertatts)

     lpiatts = 0.0
     lpiatts@long_name = "Lightning Potential Index"
     lpiatts@units = "J/kg"
     filevarattdef(f,"lpi", lpiatts)

     filevarchunkdef(f, "lpi", chunk_sizes)
     filevarcompressleveldef(f, "lpi", 2)
     
     ;Create coordinate arrays
     print("Creating coordinate arrays")
     times = fspan(0,ntim-1,ntim)
     ;lat = wrf_user_getvar(a,XLAT,0)
     f->lat = (/xlat/)
     ;lon = wrf_user_getvar(a,XLONG,0)
     f->lon = (/xlon/)
     f->vert = (/nk/)
     ;Extract 2D array
     print("Extracting LPI")
     f->lpi(tn,:,:,:) = (/lpi/)
     tn = tn+1
    end
   end do
