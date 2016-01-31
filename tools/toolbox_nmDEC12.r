source("../lib/ms_lmata.r")
#========================================================================
#	data descriptions and routines to work with regions of interest
# 	(ROI) in GC/MS or LC/MS data
#	A ROI will be represented as an R list having the following
#	members:xpt.name, roi.name, mz.s (start), mz.e (end), mz.b (base mass 
#	usually the same as start, but may be different in
#	some cases; e.g. m-1 may be interesting as a fraction of unlabeled)
#   time (vector of time points), tmi (dataframe w time, mz, and int. cols)
#   eic (vector of intensity for mzb), gic (group ion chromatogram; 
#	summation of intensities for given masses during the time interval of 
#	interest), data (intensities for mz from start to end during the time
#	period)
#
#   Note:  roi.name and mz.b will be deduced from mz.s unless given
#   otherwise.  roi.name will be the string representation of mz.s
#   appended to "mz" if left unspecified
#
# 	pull out an ROI from either an xxcms object with a "tmi" member
# 	or from a tmi object with "time", "mz", "intensity" columns
#	OR I always use an lmata read already (by using read.lmata function)
#=========================================================================

xxcms2roi = function(xxcms, mzs, times, mz.b=NULL, xpt.name=NULL, 
roi.name=NULL) {
	mz.s = mzs[1]
	mz.e = mzs[2]
	mz.b = NULL
	if( is.null(mz.b) ) {
	mz.b = mz.s
    }	
 
	in.tmi = extract.tmi(xxcms)
	tmi = time.slice(in.tmi, times)
	tmi = mz.slice(tmi, mzs)

	time = unique(tmi$time)
	mz.tmi = mz.slice(tmi, mz.b)
	eic = mz.tmi$intensity

	gic = tapply(X=tmi$intensity, INDEX=list(tmi$time), FUN=sum)
	rownames(gic) = NULL
	gic = as.vector(gic)
#	gic = cbind(mz.tmi$time,gic) #just in case you want to check what gic is
    roi=list(mz.s=mz.s, mz.e=mz.e, tmi=tmi, time=time, eic=eic, gic=gic)
}

#===========================================
#	closestTimeVal
#	select index in times vector closest to
#	target time (target)
#	based on closestColPredicate in mez_utils
#===========================================
closestTimeVal <- function(timeVec,target) {
	delta=abs(timeVec-target)
	ndx=which(delta == min(delta))
	return(ndx)
	}

#===============================================
#	makeROI
#	makes region of interest.
#	input needed is dataset (lmata read already)
#	RT and base mass that are extracted from cvs
#===============================================
makeROI <- function(dataset,RT,baseMass,windowroi) {
	ndx=closestTimeVal(times,RT)
	timeVals=c(times[(ndx-windowroi)],times[(ndx+windowroi)]) 
	mzVals=c(baseMass,(baseMass+8))
	roi=xxcms2roi(dataset,mzs=mzVals,times=timeVals)     
	}
	
makeROI_fixed=function(dataset,RT,baseMass,windowroi){
	sectionROI=time.slice(dataset,c(RT-windowroi,RT+windowroi)) 
	#if RT is 850 and windowroi is 10, the section would go from 840 to 	#860
	sectionROIMZ=mz.slice(sectionROI,baseMass)
	mROI = max(sectionROIMZ$intensity)
	sROI = cbind(sectionROIMZ$time,sectionROIMZ$intensity)
	max.int.posROI	=	which(sectionROIMZ$intensity==mROI[1])
	max_time	=	sROI[max.int.posROI,1]
	print(max_time)
	timeValsROI=c(max_time-(windowroi-2),max_time+(windowroi-2)) 
	mzVals=c(baseMass,(baseMass+8))
	roi=xxcms2roi(dataset,mzs=mzVals,times=timeValsROI) 
}

#=======================================================
#	peak.near from fit_gauss.r by Zawrotny
#	window is in seconds
#	time on each side of the stated position to test; i.e. total width =
#   2 * window
#	search is from time - window to time + window
#	return is a tmi with times in the range found_time +/- window/2
#========================================================
peak.near = function(eic, time, window=20) {
    section = time.slice( eic, c(time - window, time + window ) )
    pks 	= SpectrumSearch(section$intensity)
    pk_idx 	= pks$pos[1]
    pk_time = section$time[pk_idx]
    half_win = window / 2.0
    time.slice(section, c(pk_time - half_win, pk_time + half_win))
	}

#specifically for norleucine mz=200 that interferes with isoleucine mz=200 

peak.nearNOR200=function(eic,time,window=20){
	section = time.slice(eic, c(time - window, time + window ) )
	m = max(section$intensity)
	s = cbind(section$time,section$intensity)
	max.int.pos	=	which(section$intensity==m[1])
	pk_timeN	=	s[max.int.pos,1]
	half_win = window / 2.0
	time.slice(section, c(pk_timeN - half_win, pk_timeN + half_win))
	}
## maybe change ile to 837. The window here is in seconds. If time is 678 
#  and window is 20, then the section is from 658 s to 698 s.

peak.nearILE=function(eic,time,window=10){
	section = time.slice(eic, c(time - window, time + window ) )
	m = max(section$intensity)
	s = cbind(section$time,section$intensity)
	max.int.pos	=	which(section$intensity==m[1])
	pk_timeN	=	s[max.int.pos,1][1]
	time.slice(section, c(pk_timeN - window, pk_timeN + window))
	}

#=========================================
#	fit.gauss from fit_gauss.r by Zawrotny
#
#========================================
fit.gauss<-function(pk) {
	bi=mean(pk[,"time"])
	ai=max(pk[, "intensity"])
	fg=tryCatch(nls(intensity ~ const + a*exp(-(time-b)^2/(2*c^2)),data=pk, start = 		list(const=1,a=ai,b=bi,c=1), 												control=list(minFactor=1/4096,maxiter=1000)),error=function(e)NA)
	return(fg)
	}

#from Pam_EIC_tools. fit_gaussS used with finds_and_fits and fit.gauss_init 
#used in find_and_fit_init.	The tryCatch function is critical to avoid the
#errors. If you apply it to a peak that will give you errors (due to 
#singular gradients for example), it will assign to fg=NA instead of the 
#red error message "Error in nls...singular gradient"so the code doesnt stop

fit.gaussS<-function(pk,fg_coef) {
	ai=max(pk[, "intensity"])
	bi=fg_coef[3]
	ci=fg_coef[4]   
 	fg=tryCatch(nls(intensity ~ const + a*exp(-(time-b)^2/(2*c^2)),data=pk, 	start = list(const=1,a=ai,b=bi,c=ci), 										control=list(minFactor=1/4096,maxiter=1000)),error=function(e)NA)
 	return(fg)
	}

fit.gauss_init<-function(pk) {
	bi=mean(pk[,"time"])
	ai=max(pk[, "intensity"])
	fg=tryCatch(nls(intensity ~ const + a*exp(-(time-b)^2/(2*c^2)),data=pk,
	start = list(const=1,a=ai,b=bi,c=1), 										control=list(minFactor=1/4096,maxiter=1000)),error=function(e)NA)
	return(fg)
	}

#=========================================
#	plot.gaussfit from fit_gauss.r by Zawrotny
#
#=========================================
plot.gaussfit = function(pk, fg, mz=NA) {
	time = pk[,"time"]
	intensity = pk[,"intensity"]
	int.fit = fitted(fg)
	plot(time, intensity, pch=19, cex=1,col="blue")
	lines(time, int.fit, col="purple")
	title(c(metab,mz))
}
#=================================================================
#	finds_and_fit 
#	finds and fits a given mz and rt peak from an xxcms / tmi object
# 	returns a list (fp for "found peak") with named members:
#	"pk" is the extracted peak
#	"fg" is the "fit gauss" object from the nls run
#=================================================================
find_and_fit = function(xxcms, mz, rt, window) { 
    eic = mz.slice(xxcms, mz)
    pk  = peak.nearILE(eic, rt, window=window)
    fg  = fit.gauss(pk)
    fp  = list( pk=pk, fg=fg)
	}

#	from Pam_EIC_tools.finds_and_fits is used after the ff_init; see that
#	fg_coef is a new parameter used.You can replace xxcms for roi found
#	already.

find_and_fit_init = function(xxcms, mz, rt, window) { 
    eic = mz.slice(xxcms, mz)
    pk  = peak.nearILE(eic, rt, window=window)
    fg  = fit.gauss_init(pk)
    fp  = list( pk=pk, fg=fg)
	}

finds_and_fits = function(xxcms, mz, rt, window,fg_coef) { 
    eic = mz.slice(xxcms, mz)
    pk  = peak.nearILE(eic, rt, window=window)
    fg  = fit.gaussS(pk,fg_coef)
    fp  = list( pk=pk, fg=fg)
	}

#============================================
#	n_getPeakArea function: to get peak area
#	using RT time in metabs file. Needs also a window
#	and the roi for each metabolite for each dataset
#============================================
n_getPeakArea <- function (roi,RT,baseMass,window) {
	EIC=roi$eic
	lenEIC=length(EIC)
	meanE=mean(EIC)  #taking mean of all the intensity values for this eic
	if (meanE == 0.0) {  ##if the mean intensity is zero, skip the gauss
	A=0
	} 		
	else {
	diff1=1
	for (k in 1:(lenEIC-1)){
		diff= (EIC[(k+1)]-EIC[k])
		diff1=(diff1+diff)} 
		if (diff1 == 2.0){
		A=0
			}    
	else {  #if the mean isn't 0, fit gaussian & calculate area
	fG=find_and_fit(inData,baseMass,RT,window)
	plot.gaussfit(fG$pk,fG$fg) 
	fGsummary=summary(fG$fg)
	a= fGsummary$coefficients[2]
	c= fGsummary$coefficients[4]	
	Area = a*c*sqrt(2*pi)
	print(Area)
	fp = predict(fG$fg)-predict(fG$fg)[1]
	A=sum(fp[which(fp > 0.025*max(fp))])
	return(A)
			}
		}
	}

new_getPeakArea_Simpsons <- function (lmata_dataset,RT,baseMass,window) {
	len.pk=c()
	Total.pks=c()
	Total.area=c()
	eic = mz.slice(lmata_dataset, baseMass)
	pk = peak.near(eic, RT, window=window)
		Total.pks=which(pk$intensity > 0.025*max(pk$intensity))
		len.pk=length(Total.pks)  
		meanY=mean(pk$intensity)  ##taking the mean of all the intensity
		if (meanY == 0.0) {  #if mean inten. is zero, skip the calculation
			Total.area=0.0
			Yobs=c(Yobs,Total.area)
			} #ok
			else {
			Total.pks=which(pk$intensity > 0.025*max(pk$intensity))
			len.pk=length(Total.pks)    
	 Total.area=simp(pk[Total.pks[1]:Total.pks[len.pk]],x=c(seq(1,len.pk)))
			A=Total.area
			return(A)
			}
		}

#============================================
#	n_getPeakArea_param function: to get peak area
#	using RT time in metabs file. Needs also a window
#	and the roi for each metabolite for each dataset. 
# 	main difference with respect to n_getPeakArea is that
# 	Area is calculated from parameters obtained in fitgauss
#	Integrating the function of nlo, at the end area=a*c(2*pi)^1/2
#============================================
n_getPeakArea_param <- function (roi,RT,baseMass,window) {
	EIC=roi$eic
	lenEIC=length(EIC)
	meanE=mean(EIC)  #taking mean of all the intensity values for this eic
	if (meanE == 0.0) {  ##if the mean intensity is zero, skip the gauss
	A=0
	} 		
	else {
	diff1=1
	for (k in 1:(lenEIC-1)){
		diff= (EIC[(k+1)]-EIC[k])
		diff1=(diff1+diff)} 
		if (diff1 == 2.0){
		A=0
			}    
	else {  #if the mean isn't 0, fit gaussian & calculate area
	fG=find_and_fit(inData,baseMass,RT,window)
	plot.gaussfit(fG$pk,fG$fg) 
	fGsummary=summary(fG$fg)
	a= fGsummary$coefficients[2]
	c= fGsummary$coefficients[4]	
	print(fGsummary$coefficients[2]) #a
	print(fGsummary$coefficients[4]) #c
	Area = a*c*sqrt(2*pi)
	A=Area
	return(A)
	print(A)
			}
		}
	}

getPeakArea=function(roi){
	t=which(roi$gic == max(roi$gic))
	bi=roi$time[t]
	ai=max(roi$gic)
	dat=data.frame(cbind(roi$time,roi$gic))
	colnames(dat)=c("time","intensity")
	timePeak=dat[,"time"]
	intensityPeak=dat[,"intensity"]
	gaussfit=tryCatch(nls(intensity ~ const + a*exp(-(time-b)^2/(2*c^2)),data=dat,start=list(const=1,a=ai,b=bi,c=1),control=list(minFactor=1/4096,maxiter=1000)),error=function(e)NA)
	FG = list( pk=dat, fg=gaussfit)
	if (is.na(FG$fg)[1]){
		A=0
  		fitted.int=rep.int(0,length(intensityPeak))}
  	else{
	Fg=gaussfit
	fitted.int=fitted(Fg)
	#plot.gaussfit(FG$pk,FG$fg) 
  	FP = predict(FG$fg)-predict(FG$fg)[1]
  	A=sum(FP[which(FP > 0.025*max(FP))])} 		
	#a_norm=(a/NorLeu_area)/cell_ct  ###adding this as a test
	print(A)
	print(bi)
	return(list(A,timePeak,intensityPeak,FP))
	#return(A) #BEFORE we were returning fitted.int
	}

#================
# get time of the maximum point of the peak generated from the gic
#================
getPeakTime=function(roi){
	t=which(roi$gic == max(roi$gic))
	bi=roi$time[t]
	return(bi)
}

#=================================================
# get highest intensity of an ion in a time window
#=================================================
get.inten = function(roi, time, baseMass, window) {
    EIC=roi$eic
    eic = mz.slice(EIC, baseMass)
    section = time.slice( eic, c(time - window/2.0, time + window/2.0 ) )
	max.int = (max(section$intensity))
	s = cbind(section$time,section$intensity)
	max.int.pos	=	which(section$intensity==max(section$intensity))[1]
	pk_time_max	=	s[max.int.pos,1]
	pk_time_max =	pk_time_max/60   
	print(c(max.int,pk_time_max))
	return(max.int)
	}

#===============================================
#
#===============================================
getEic <- function(tmi,mz,t=c(0,Inf)) {
    tmi = selectColVal(tmi, "mz", mz)
    eic = tmi[, c("time", "intensity")]
    eic = selectColVal(eic,"time", t)
    eic
}

#===============================================
#
#===============================================
new_try_makeObs <- function(roi,window) {
	mzs=seq(roi$mz.s:max(roi$tmi[,"mz"]))+roi$mz.s -1
	#mzs=seq(roi$mz.s:(roi$mz.s+numC))+roi$mz.s -1
	Yobs=c()
	for (mz in mzs) {  ##cycle thru mz's in the roi, getting peak areas
		e=getEic(roi$tmi,mz)
		if (max(e$intensity) <= 0.025*max(roi$eic)) {
		a=0.0
		Yobs=c(Yobs,a)
		}
		else{
		if (max(e$intensity) <= 550) {
		a=0.0
		Yobs=c(Yobs,a)		
		}
		else {  ##if the mean is not zero, fit gaussian & calculate area
				#fgFrac=find_and_fit(e,mz,RT,window)
			fgFrac=find_and_fit(inData,mz,RT,window)
			plot.gaussfit(fgFrac$pk,fgFrac$fg)  #can comment out later
			fpFrac=predict(fgFrac$fg,e[,1])-predict(fgFrac$fg,e[,1])[1] 
			a=sum(fpFrac[which(fpFrac > 0.025*max(fpFrac))])
			Yobs=c(Yobs,a)
				}
			}}
			Yobs=data.frame(mzs,Yobs)
			colnames(Yobs)=c("m/z","Intensity")
			return(Yobs)
	}

#	From Pam_EIC_tools. I CANNOT DO find_and_fit_init on e because in find 
# 	and_fit_init it does an mz.slice that needs a tmi and I don't  have that

alt_makeObs <- function(roi,window) {
	mzs=seq(roi$mz.s:max(roi$tmi[,"mz"]))+roi$mz.s -1
	Yobs=c()
	for (mz in mzs) {  
		e=getEic(roi$tmi,mz)           
			# if (max(e$intensity) <= 550) {
			if (mz==roi$mz.s) {
			fgFrac=find_and_fit_init(roi,mz,RT,window)
			#plot.gaussfit(fgFrac$pk,fgFrac$fg,mz)
			fg_coef=tryCatch(coef(fgFrac$fg),error=function(e)NA)}
          	else
			{fgFrac=finds_and_fits(roi,mz,RT,window,fg_coef)}
			if (is.na(fgFrac$fg[1]))   {a=0.0}
			# if (max(e$intensity)<= 550){a=0.0}
			# else{
			else {
			fpFrac=predict(fgFrac$fg,e[,1])-predict(fgFrac$fg,e[,1])[1] 
          	a=sum(fpFrac[which(fpFrac > 0.025*max(fpFrac))])}
			
			Yobs=c(Yobs,a)
			if (mz > roi$mz.s+3){window=5}
            }
			Yobs=data.frame(mzs,Yobs)
			colnames(Yobs)=c("m/z","Intensity")
			return(Yobs)
		}


short_alt_makeObs=function(ionsplot_output){
	mzs=seq(roi$mz.s:max(roi$tmi[,"mz"]))
	L=c()
	Y=c()
	for (i in mzs){
			l=ionsplot_output[[1]][[2]][[i]][[1]]
			L=c(L,l)
			y=ionsplot_output[[1]][[2]][[i]][[5]]
			Y=c(Y,y)}
			mzint_list=data.frame(L,Y)
			colnames(mzint_list)=c("m/z","Intensity")
			return(mzint_list)
			}
			
			

#=====================================================================
# Bulding G which is a huge list of lists that has the information to 
# later plot the traces of each M to M+8 with the help of plot_multi
#=====================================================================
ions.plot <- function(roi,window) {
	E=c()
	mzs=seq(roi$mz.s,max(roi$tmi[,"mz"]),by=1)
	for (mz in mzs) {  
			if (mz==roi$mz.s) {
				fgFr=as.data.frame(cbind(roi$time,roi$eic)) 
				colnames(fgFr)=c("time","intensity") 
				fgFra=fit.gauss_init(fgFr)
				fgFrac=list(pk=fgFr,fg=fgFra)
				fg_coef=tryCatch(coef(fgFrac$fg),error=function(e)NA)
					if (is.na(fgFrac$fg[1])){
				PKtime=fgFrac$pk[1]
				PKintensity=fgFrac$pk[2]
				PKTI=cbind(PKtime,PKintensity)	
				fitted.fg=rep.int(0,length(PKTI$time))
				a=0.0
				c=list(list(mz,PKTI$time,PKTI$intensity,fitted.fg,a))
				}
				else{
				PKtime=fgFrac$pk[1]
				PKintensity=fgFrac$pk[2]
				PKTI=cbind(PKtime,PKintensity)
				fitted.fg=fitted(fgFrac$fg)
				FP = predict(fgFrac$fg)-predict(fgFrac$fg)[1]
				a=sum(FP[which(FP > 0.025*max(FP))])
				c=list(list(mz,PKTI$time,PKTI$intensity,FP,a))
				}
				}
			else {fgFrac=finds_and_fits(roi,mz,RT,window,fg_coef)
				if (is.na(fgFrac$fg[1])){
				PKtime=fgFrac$pk[1]
				PKintensity=fgFrac$pk[3]
				PKTI=cbind(PKtime,PKintensity)	
				fitted.fg=rep.int(0,length(PKTI$time))
				a=0.0
				c=list(list(mz,PKTI$time,PKTI$intensity,fitted.fg,a))
				}
				else{ 
				PKtime=fgFrac$pk[1]
				PKintensity=fgFrac$pk[3]
				PKTI=cbind(PKtime,PKintensity)
				fitted.fg=fitted(fgFrac$fg)
				FP = predict(fgFrac$fg)-predict(fgFrac$fg)[1]
				a=sum(FP[which(FP > 0.025*max(FP))])
				c=list(list(mz,PKTI$time,PKTI$intensity,FP,a))				
				}}
				E=c(E,c)
				}
				E=list(list(metab,E))
				return(E)
				}


colors=c("white","blue","black","green","purple","turquoise4","orange","magenta","darkblue")

##=========================================================
## Plotting the M to M+8 of each metabolites for all samples
##==========================================================
plot_multi=function(G,j,N){
	plot(G[[j]][[N+1]][[2]][[1]][[2]],G[[j]][[N+1]][[2]][[1]][[3]],xlim=c(min(G[[j]][[N+1]][[2]][[1]][[2]])-0.5,max(G[[j]][[N+1]][[2]][[1]][[2]])+0.5),ylim=c(0,max(G[[j]][[N+1]][[2]][[1]][[4]])),xlab="Time(sec)"
 ,ylab="Intensity",col="red")
	lines(G[[j]][[N+1]][[2]][[1]][[2]],G[[j]][[N+1]][[2]][[1]][[4]],xlim=c(min(G[[j]][[N+1]][[2]][[1]][[2]])-0.5,max(G[[j]][[N+1]][[2]][[1]][[2]])+0.5),ylim=c(0,max(G[[j]][[N+1]][[2]][[1]][[3]])),col="red")
	legend(min(G[[j]][[N+1]][[2]][[1]][[2]])-0.4,max(G[[j]][[N+1]][[2]][[1]][[4]]), c("M","M+1","M+2","M+3","M+4","M+5","M+6","M+7","M+8"), cex=0.8,col=c("red","blue","black","green","purple","turquoise4","orange","magenta","darkblue"), pch=21:21, lty=c(1:1),lwd=c(1.0,1.0))
	title(c(G[[j]][[1]],G[[j]][[N+1]][[1]]))
	# axis(1,at=seq(660,700,by=5))
	# axis(2,at=seq(0,max(E[[1]][[3]]),1E5))
	par(new=T)
	for (n in 2:9){
		par(new=T)
	plot(G[[j]][[N+1]][[2]][[n]][[2]],G[[j]][[N+1]][[2]][[n]][[3]],xlim=c(min(G[[j]][[N+1]][[2]][[n]][[2]])-0.5,max(G[[j]][[N+1]][[2]][[1]][[2]])+0.5), ylim=c(0,max(G[[j]][[N+1]][[2]][[1]][[3]])),xaxt="n",yaxt="n",ann=F,col=colors[n])
	lines(G[[j]][[N+1]][[2]][[n]][[2]],G[[j]][[N+1]][[2]][[n]][[4]],xlim=c(min(G[[j]][[N+1]][[2]][[n]][[2]])-0.5,max(G[[j]][[N+1]][[2]][[1]][[2]])+0.5),ylim=c(0,max(G[[j]][[N+1]][[2]][[1]][[4]])),xaxt="n",yaxt="n",ann=F,col=colors[n])
		}}
	# if (n==9){
	# plot(G[[j]][[N+1]][[2]][[9]][[2]],G[[j]][[N+1]][[2]][[9]][[3]],xlim=c(min(G[[j]][[N+1]][[2]][[n]][[2]])+5,max(G[[j]][[N+1]][[2]][[9]][[2]])-5), ylim=c(0,max(G[[j]][[N+1]][[2]][[1]][[3]])),xaxt="n",yaxt="n",ann=F,col=colors[n])
	# lines(G[[j]][[N+1]][[2]][[n]][[2]],G[[j]][[N+1]][[2]][[9]][[4]],xlim=c(min(G[[j]][[N+1]][[2]][[9]][[2]])+5,max(G[[j]][[N+1]][[2]][[9]][[2]])-5),ylim=c(0,max(G[[j]][[N+1]][[2]][[1]][[4]])),xaxt="n",yaxt="n",ann=F,col=colors[n])}


#==============================================
# makeA for fractional abundance
#==============================================
makeA <- function(Formula) {
    mid=round(isotope_distribution(Formula),4)
	mid=subset(mid,freq>0.00)
	len.mid=length(mid[,2])
	A=matrix(0,nrow=len.mid,ncol=len.mid)
	A[,1]=mid[,2]
	for (i in 2:len.mid) {
		j=i-1
		f_a=formula2array(Formula)
		f_a[1,2]=as.numeric(f_a[1,2])-j
		form=array2formula(f_a)
		mid=round(isotope_distribution(form),4)
		mid=subset(mid,freq>0)
		A[i:len.mid,i]=mid[i:len.mid-j,2]
	}
	return(A)
}
#====================================================================
# this is a function that uses the R "solve" command to calculate
# the coefficients, x, relating the abundance matrix and the observed
# mass spectra intensities
#=====================================================================
Solve <- function(A_matrix,Yobs) {
  K=t(A_matrix)%*%A_matrix
  y=t(A_matrix)%*%Yobs[1:dim(A_matrix)[1],2]
  x=solve(K,y)
  x=x[1:(num_C_frag+1)]
  x=x/sum(x)
  return(x)
}

Solve_nnls=function(A_matrix,Yobs){
	K=t(A_matrix)%*%A_matrix
	y=t(A_matrix)%*%Yobs[1:dim(A_matrix)[1],2]
	x=nnls(K,y)
	x=x$x
	x=x/sum(x)
	return(x)
	}

#=============================================================
#	getPeakFracAbund  function to return Fractional Abundance
#	for defined peak
#=============================================================
new_getPeakFracAbund <- function (roi,Formula,window) {
	Yobs=new_try_makeObs(roi,11)
	Amat=makeA(Formula)
	X=Solve(Amat,Yobs) 
	} 

getPeakFracAbund <- function (roi,Formula,window) {
	Yobs=alt_makeObs(roi,10)
	Amat=makeA(Formula)
	X=Solve_nnls(Amat,Yobs) 
	} 
	
n_getPeakFracAbund <- function (ionsplot_output,Formula,window) {
	Yobs=short_alt_makeObs(ionsplot_output)
	Amat=makeA(Formula)
	X=Solve_nnls(Amat,Yobs) 
	} 


plotexpt_fitt = function(B,Sets,Metab) {
	n=Metab+1
	timek = B[[Sets]][[n]][[3]]
	intensityk = B[[Sets]][[n]][[4]]
	int.fitk = B[[Sets]][[n]][[5]]
	plot(timek, intensityk, pch=19, cex=1,col="blue")
	lines(timek, int.fitk, col="purple")
	title(c(B[[Sets]][[1]],B[[Sets]][[n]][[1]]))
	#dev.off()
	}

MPE=function(X,num_C_frag) {
	R_mpe=rep(0,num_C_frag+1)
  	all_M=1-X[1]
	for (n in 2:(num_C_frag+1)) {
	      R_mpe[n]=(X[n]/all_M)*100
	      }
	L=list((all_M*100),R_mpe)
	return(L)
  }

plot_traces=function(sequence,path,metabs,G,start_expt,end_expt){
	for (N in sequence){
	singles.path= paste(path,metabs$Metabolite[N],sep="/")
	singles.path= paste(singles.path,".pdf",sep="")
	pdf(singles.path)
	par(mfcol=c(2,3))
	for (j in start_expt:end_expt){
	plot_multi(G,j,N)
	}	
	dev.off()	
	}
	}