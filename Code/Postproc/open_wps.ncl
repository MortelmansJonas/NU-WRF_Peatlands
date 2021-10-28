wps_filename = "SSTRSS:2015-05-15_00"

;---Preallocate variables for header
  header      = new(14,float)
  field       = new(1,string)
  date        = new(1,string)
  units       = new(1,string)
  map_source  = new(1,string)
  description = new(1,string)
  
  istatus =  wrf_wps_open_int(wps_filename)
  if(istatus.ne.0) then
    print("Error opening " + wps_filename)
    exit
  end if

;---Continuously read data from WPS file until istatus != 0.
  nfields = 0

  do while (istatus.eq.0)

;---Read header
    wrf_wps_rdhead_int(istatus,header,field,date,units,\
                       map_source,description)
    if(istatus.ne.0) then
      continue
    end if

    nx = toint(header(3))
    ny = toint(header(4))

;---Read data
    slab := wrf_wps_rddata_int(istatus,nx,ny)
    if(istatus.ne.0) then
      continue
    end if

;---Print information about each slab
    print("Field #" + (nfields+1))
    print("  name          : '" + field + "'")
    print("  description   : '" + description + "'")
    print("  units         : '" + units + "'")
    print("  min/max value : "  + min(slab) + " / " + max(slab))
    print("  date          : '" + date + "'")
    print("  map source    : '" + map_source + "'")
    print("  version       : " + header(0))
    print("  forecast hour : " + header(1))
    print("  level         : " + header(2))
    print("  ny x nx       : " + header(4) + " x " + header(3))
    print("  projection    : " + header(5))
    if(header(5).eq.0) then
      print("  startlat      : " + header(6))
      print("  startlon      : " + header(7))
      print("  deltalat      : " + header(8))
      print("  deltalon      : " + header(9))
      print("  earth_radius  : " + header(10))
    else if(header(5).eq.1) then
      print("  startlat      : " + header(6))
      print("  startlon      : " + header(7))
      print("  dy x dx       : " + header(9) + " x " + header(8))
      print("  truelat1      : " + header(10))
      print("  earth_radius  : " + header(11))
    else if(header(5).eq.3) then
      print("  startlat      : " + header(6))
      print("  startlon      : " + header(7))
      print("  dy x dx       : " + header(9) + " x " + header(8))
      print("  center lon    : " + header(10))
      print("  truelat1      : " + header(11))
      print("  truelat2      : " + header(12))
      print("  earth_radius  : " + header(13))
    else if(header(5).eq.4) then
      print("  startlat      : " + header(6))
      print("  startlon      : " + header(7))
      print("  nlats         : " + header(8))
      print("  deltalon      : " + header(9))
      print("  earth_radius  : " + header(10))
    else if(header(5).eq.5) then
      print("  startlat      : " + header(6))
      print("  startlon      : " + header(7))
      print("  dy x dx       : " + header(9) + " x " + header(8))
      print("  center lon    : " + header(10))
      print("  truelat1      : " + header(11))
      print("  earth_radius  : " + header(12))
    end if
    end if
    end if
    end if
    end if
    nfields = nfields + 1
  end do

  print("There were " + nfields + " fields in the " + wps_filename + " WPS file.")

