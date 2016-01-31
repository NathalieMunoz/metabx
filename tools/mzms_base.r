library(ncdf4)
#library(reshape2)
library(dplyr, warn.conflicts=FALSE)
library(rhdf5)
library(tools) # file_path_as_absolute

md.round <- function(mz, md=0.65){
    frac <- mz %% 1
    ifelse( frac >= md, ceiling(mz), floor(mz) )
}

# md.rnd argument
#   single number n; 0 < n < 1 -- "mass defect" break point, defaults to 0.65
#       fractions strictly less than round down
#       fractions greater than or equal to round up
#   single integer> = 0 --> round to that many digits to the
#       right of the decimal
#   single negative number  --> do not round at all; return the "raw"
#       tmi as collected
cdf2mzms = function(cdf.file, md.rnd=0.65) {

    filename = cdf.file

    nc = nc_open(cdf.file)
    expt_params = ncatt_get(nc ,varid=0)

    desc = ncatt_get(nc, 0, 'experiment_title')$value

    # vector[ 1 : nscans] -- time of acquisition (seconds)
    times = ncvar_get(nc, 'scan_acquisition_time')

    # vector[ 1 : nscans] -- is this ever not 1 : nscans?
    #   yes, some start at 0, others 1.  conversion artifact?
    #
    # scan_num = ncvar_get(nc, 'actual_scan_number')
    scan_num = 1:length(times)


    # vector[1 : nscans] -- 0-based point number of the start of the scan 
    #  number corresponding to the index of this vector
    #  e.g. scan_index[3] + 1 = 1-based index into intensity_values for
    #  the first point in scan 3
    scan_index = ncvar_get(nc, 'scan_index')

    # vector[1 : nscans -- number of points in each scan
    point_count = ncvar_get(nc, 'point_count')

    # vector[1 : npoints] -- mz of the given point
    mass_values = ncvar_get(nc, 'mass_values')

    # vector[1 : npoints] -- intensity of the given point
    intensity_values = ncvar_get(nc, 'intensity_values')
    intensity_values = as.vector(as.numeric(intensity_values))

    times.rep = unlist(sapply(scan_num,
                FUN=function(sn) rep(times[sn],each=point_count[sn])))

    raw_tmi = data.frame(times.rep, mass_values, intensity_values)
    colnames(raw_tmi) = c("time", "mz", "intensity")

    if ( md.rnd < 0) {
        rounding = "raw"
        tmi =raw_tmi
        mz.min = min(tmi$mz)
        mz.max = max(tmi$mz)
    } else {
        if ( is.integer(md.rnd) ) {
            rounding = "round"
            raw_tmi$mz <-round(raw_tmi$mz, md.rnd)
        } else {
            rounding = "mass_defect"
            raw_mz <- md.round(raw_tmi$mz, md=md.rnd)
            raw_tmi$mz <- as.vector(raw_mz)
        }

        tmi_df = tbl_df(raw_tmi)
        grouped = group_by(tmi_df, time, mz)
        tmi = summarize(grouped, intensity=sum(intensity))
        colnames(tmi) <- c("time", "mz", "intensity")
        tmi = ungroup(tmi)

    }
    attr(tmi, 'class') = c('tmi', class(tmi))

    times = sort(unique(tmi$time))
    ti = rank(times) # time index values
    mzs = sort(unique(tmi$mz))
    mzi = rank(mzs) # mz index values

    mzms = list( mz=mzs, time=times, tmi=tmi, ti=ti, mzi=mzi)

    expt_params$rounding = rounding
    expt_params$cdf_file = file_path_as_absolute(cdf.file)
    mzms$expt_params = expt_params

    attr(mzms, 'class') = c('mzms', class(mzms))
    mzms
}

tmi2msblock = function(tmi, fill.val=0) {
    time = sort(unique(tmi$time))
    ti = rank(time)
    tis = tbl_df(data.frame(time, ti))

    mz.range = range(tmi$mz)
    mz.low = mz.range[1]
    mz.high = mz.range[2]
    mz = seq(from=mz.low, to=mz.high)
    mi = seq(length(mz))
    mzis = tbl_df(data.frame(mz, mi))

    tmi.idx = left_join(left_join(tmi, tis), mzis)
    
    nrow = length(time)
    ncol = length(mz)

    block = matrix(fill.val, nrow=nrow, ncol=ncol)
    block[cbind(tmi.idx$ti, tmi.idx$mi)] = tmi.idx$intensity

    msblock = list(mz=mz, time=time, block=block)
    attr(msblock, 'class') = c('msblock', class(msblock))
    msblock
}

msblock2tmi = function(msblock, thresh=100) {
    idxs = which(msblock$block > thresh, arr.ind=TRUE)
    block = msblock$block
    intensity = block[ idxs ]
    colnames(idxs) = c('ti', 'mi')
    i.idxs = tbl_df(data.frame(cbind(intensity, idxs)))

    time = msblock$time
    ti = rank(time)
    tis = tbl_df(data.frame(time, ti))

    mz = msblock$mz
    mi = rank(mz)
    mzs = tbl_df(data.frame(mz, mi))

    tmi.idx = left_join(left_join(i.idxs, mzs), tis)
    tmi_df = select(tmi.idx, c(time, mz, intensity))

    # ensures that tmi is sorted by time, then mz
    tmi = arrange(tmi_df, time, mz)
    attr(tmi, 'class') = c('tmi', class(tmi))
    tmi
}

mzms2hdf = function(mzms, hdf_file, raw=FALSE) {
    if( ! file.exists(hdf_file) ) {
       h5createFile(hdf_file)
    }
    if(raw) {
        h5write(mzms, hdf_file, name="raw", level=4)
    } else {
        h5save(mzms, file=hdf_file, level=4)
    }

}

hdf2mzms = function(hdf_file, raw=FALSE) {
    if(raw) {
        mzms = h5read(hdf_file, name="/raw", level=4)
    } else {
        mzms = h5read(hdf_file, name="/mzms", level=4)
    }
    mzms$mz = as.vector(mzms$mz)
    mzms$mzi = as.vector(mzms$mzi)
    mzms$ti = as.vector(mzms$ti)
    mzms$time = as.vector(mzms$time)
    mzms$tmi$time = as.vector(mzms$tmi$time)
    mzms$tmi$mz = as.vector(mzms$tmi$mz)
    mzms$tmi$intensity = as.vector(mzms$tmi$intensity)
    mzms$tmi = tbl_df(data.frame(mzms$tmi))
    attr(mzms$tmi, 'class') = c('tmi', class(mzms$tmi))
    attr(mzms, 'class') = c('mzms', class(mzms))
    mzms
}

timeslice = function(tmi, t1, t2=NULL) {
    # only one time given
    if ( is.null(t2) ) {
        tmi %>% 
            mutate(adelta = abs(time - t1))  %>%
            filter(adelta == min(adelta)) %>%
            select(time, mz, intensity)
    } else {
        if ( t1 < t2 ) {
            lo = t1
            hi = t2
        } else { 
            lo = t2
            hi = t1
        }
        tmi %>% 
            filter( time >= lo ) %>%
            filter( time <= hi)
    }
}

mzslice = function(tmi, m1, m2=NULL) {
    # only one time given
    if ( is.null(m2) ) {
        tmi %>% 
            filter(mz == m1)
    } else {
        if ( m1 < m2 ) {
            lo = m1
            hi = m2
        } else { 
            lo = m2
            hi = m1
        }
        tmi %>% 
            filter( mz >= lo ) %>%
            filter( mz <= hi)
    }
}
