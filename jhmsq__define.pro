;o=jhMSQ(b,fkhz=20)
;
;
;
;
;
;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|
pro jhMSQ::plot
  plot,self.LVec()*(*self.Tsamp),self.MSQ(),/xlog,psym=2
end


;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|


function jhMSQ::fit, iseg, parain=parain, LMIN=LMIN,LMAX=LMAX,FIXED=FIXED,FUNCTIONNAME=FUNCTIONNAME,limited=limited $
  ,upperlimits=upperlimits,lowerlimits=lowerlimits,chi2=chi2,yfit=yfit,xlog=xlog,perror=perror,quiet=quiet $
  ,randomassignment=randomassignment,experimentalvariance=experimentalvariance,parinfo=parinfo,q0=q0

  if ~keyword_set(q0) then Qshift = 1d else qshift = 0d
  alpha =  10000d
  TSampling = *self.tsamp
  if (keyword_Set(randomassignment)+self.randommsq) gt 0.1 then self.exe,randomassignment=randomassignment

  LVec = self.Lvec()
  MSQ  = self.MSQ()
  VAR  = self.Var(experimentalvariance=experimentalvariance)
  if ~keyword_set(functionname) then fnname='ss' else fnname=functionname
  
  res=fit_biasQestimator_noshotnoise( LVec, MSQ, var, parain, model=fnname,yfit=yfit,paraout=paraout, QSHift=Qshift,alpha=alpha,TSampling=TSampling,FIXED=FIXED,limited=limited,upperlimits=upperlimits,lowerlimits=lowerlimits,meank=1,bestnorm=bestnorm,xlog=xlog,parinfo=parinfo)
  
  parastr = string(res.paraout)
  titlestr=''
  chi2=res.chi2
  for i=0,n_elements(parastr)-1 do titlestr = titlestr + parastr[i] + ', '
  if ~keyword_set(quiet) then begin
    plot,LVec,MSQ,psym=2,/xlog
    oplot,LVec,yfit,thick=2,color=cgcolor('red')
  endif
  
  ;oplot,Lvec, biasQestimatorSigmoidal( Lvec, [res.paraout[0:1],0,0], QSHift=Qshift,alpha=alpha,TSampling=TSampling),color=cgcolor('blue')

return,res.paraout

end


;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|
function jhMSQ::calculateMSQ, nseg,base,LMIN,LMAX=LMAX,MSQrec=MSQrec,randomassignment=randomassignment,q0=q0,quiet=quiet $
  ,coverfactor=coverfactor,multiexperimentfactor=multiexperimentfactor,randomizedata=randomizedata
  if keyword_Set(randomassignment) then self.randommsq = 1 else self.randommsq = 0
  if keyword_set(coverfactor) then cf = coverfactor else cf = 1
  if keyword_set(multiexperimentfactor) then me = multiexperimentfactor else me = 1
  MSQptr = ptrArr(nseg)
  NData = n_elements(*self.pk)
  NData_per_seg = NData/Nseg
  for r=0,nseg-1 do begin
    k=(*self.pk)[NData_per_seg*r:NData_per_seg*(r+1)-1]

    if ~keyword_set(quiet) then help,k

    If ~keyword_set(LMAX) then NData10 = n_elements(k)/10.*cf/me else NData10=LMAX
    
    max_power= long(  ( alog(Ndata10) - alog(LMin) ) / alog(base) )
    LVec = long( base^findgen(max_power) * LMin)
    LVec = LVec[uniq(LVec)]
    max_power=n_elements(lvec)
    MSQvec = dblarr(max_power) & MSQVar = MSQvec & MSQvarexp = MSQvec
    
    if keyword_set(multiexperimentfactor) then strtpts = floor(randomu(seed,max_power) * (n_elements(k)-n_elements(k)*1D/multiexperimentfactor))
    
    for i=0 , max_power-1 do begin
      if keyword_set(multiexperimentfactor) then begin
        kuse = k[strtpts[i]:strtpts[i]+(n_elements(k))/multiexperimentfactor-1]
      endif else kuse = k
      if keyword_set(randomizedata) then kuse=kuse[sort(randomu(seed,n_elements(kuse)))]
      if ~keyword_Set(q0) then vec = self.Q1(LVec[i],kuse,randomassignment=randomassignment) $
        else vec = self.Q0(LVec[i],kuse,randomassignment=randomassignment)
      vec=vec[0:floor(n_elements(vec)*cf)-1]
      MSQvec[i] = avg(vec)+1./lvec[i]
      MSQvar[i] = 1./(n_elements(vec)*LVec[i])
      MSQvarexp[i] = variance(vec)/n_elements(vec)
    endfor
    fin = where(finite(msqvec) eq 1)
    msqvec = msqvec[fin]
    lvec = lvec[fin]
    msqvar = msqvar[fin]
    msqvarexp = msqvarexp[fin]
    res = {MSQvec:MSQvec, LVec:LVec, MSQvar:MSQvar, MSQvarexp:MSQvarexp}
    MSQptr[r] = ptr_new(res)
    
  endfor
  MSQrec=res
  return,MSQptr
end
;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|
function jhMSQ::Q1, lvec, k, randomassignment=randomassignment
  
  u=floor(n_elements(k)/double(lvec))
  
  if keyword_set(randomassignment) then begin
    Q1vec = dblarr(u)
    rndstarts = floor(randomu(seed,u) * (n_elements(k) - lvec))
    for l=0,u-1 do Q1vec[l] = mean(k[rndstarts[l]:rndstarts[l]+lvec-2]*k[rndstarts[l]+1:rndstarts[l]+lvec-1])/mean(k[rndstarts[l]+1:rndstarts[l]+lvec-1]) - mean(k[rndstarts[l]:rndstarts[l]+lvec-2])     
  endif else begin
    kbin = reform(k[0:u*lvec-1],lvec,u)
    k11 = mean(kbin[0:-2,*]*kbin[1:*,*],dim=1) - mean(kbin[0:-2,*],dim=1)*mean(kbin[1:*,*],dim=1)
    Q1vec = k11 / mean(kbin[0:*,*],dim=1)
  endelse
  
  return,Q1vec

end
;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|
function jhMSQ::Q0, lvec, k, randomassignment=randomassignment

  u=floor(n_elements(k)/double(lvec))

  if keyword_set(randomassignment) then begin
    Q0vec = dblarr(u)
    rndstarts = floor(randomu(seed,u) * (n_elements(k) - lvec))
    for l=0,u-1 do Q0vec[l] = mean(k[rndstarts[l]:rndstarts[l]+lvec-1]*k[rndstarts[l]:rndstarts[l]+lvec-1])/mean(k[rndstarts[l]:rndstarts[l]+lvec-1]) - mean(k[rndstarts[l]:rndstarts[l]+lvec-1])-1
  endif else begin
    kbin = reform(k[0:u*lvec-1],lvec,u)
    k11 = mean(kbin*kbin,dim=1) - mean(kbin,dim=1)*mean(kbin,dim=1)
    Q0vec = k11 / mean(kbin,dim=1)-1
  endelse

  return,Q0vec

end
;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|
function jhMSQ::MSQ, iseg, LMIN=LMIN,LMAX=LMA
  return, *self.pMSQ
end
function jhMSQ::Var, iseg, LMIN=LMIN,LMAX=LMAX,experimentalvariance=experimentalvariance
  if ~keyword_set(experimentalvariance) then varvec=*self.pVar else varvec=*self.pVarexp
  return, varvec
end
function jhMSQ::lvec, iseg, LMIN=LMIN,LMAX=LMAX

  return, *self.plvec
end
function jhMSQ::k, iseg, LMIN=LMIN,LMAX=LMAX

  return, *self.pk
end
;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|

pro jhMSQ::exe, nseg, BASE=BASE, LMIN=LMIN,LMAX=LMAX,randomassignment=randomassignment,q0=q0,quiet=quiet $
  ,coverfactor=coverfactor,multiexperimentfactor=multiexperimentfactor,randomizedata=randomizedata

  if n_ELEMENTS(BASE) EQ 0 then BASE = 1.1
  if n_ELEMENTS(LMIN) EQ 0 then LMIN = 64
  if n_ELEMENTS(NSEG) EQ 0 then NSEG = 1
  
  MSQptr = self.calculateMSQ(nseg, BASE, LMIN,LMAX=LMAX,MSQrec=MSQrec,randomassignment=randomassignment,q0=q0,quiet=quiet $
    ,coverfactor=coverfactor,multiexperimentfactor=multiexperimentfactor,randomizedata=randomizedata)

  self.pMSQ = ptr_new(MSQrec.MSQvec)
  self.pLVec = ptr_new(MSQrec.LVec)
  self.pVAR = ptr_new(MSQrec.MSQVar)
  self.pVarexp = ptr_new(MSQrec.MSQVarexp)
  self.NSeg = NSeg
  
  if ~keyword_set(quiet) then self.plot
end



;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|

pro jhMSQ::setProperty, fkHz=fkHz, fHz=fHz, tms=tms, tus=tus, QUIET=QUIET

if n_elements(fkhz) NE 0  then begin
  self.paras.f     = fkhz*1d3
  self.paras.fUnit = 'Hz' 
  self.paras.T     = 1/self.paras.f 
  self.paras.TUnit = 's'
  self.reset ;changes in parameters requires recalculation of msq
endif
if n_elements(fhz) NE 0  then begin
  self.paras.f     = fhz
  self.paras.fUnit = 'Hz' 
  self.paras.T     = 1/self.paras.f 
  self.paras.TUnit = 's'
  self.reset ;changes in parameters requires recalculation of msq
endif
if n_elements(tms) NE 0  then begin
  self.paras.T     = tms/1d3
  self.paras.TUnit = 's'
  self.paras.f     = 1/self.paras.T
  self.paras.fUnit = 'Hz'
  self.reset ;changes in parameters requires recalculation of msq   
endif
if n_elements(tus) NE 0  then begin
  self.paras.T     = tms/1d6
  self.paras.TUnit = 's'
  self.paras.f     = 1/self.paras.T
  self.paras.fUnit = 'Hz'
  self.reset ;changes in parameters requires recalculation of msq   
endif
if n_elements(QUIET) then self.QUIET=quiet

return
end

;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|


pro jhMSQ::reset

  if ptr_valid(self.pMSQ) then ptr_free,self.pMSQ
  if ptr_valid(self.pVar) then ptr_free,self.pVar
  if ptr_valid(self.pVarexp) then ptr_free,self.pVarexp
  if ptr_valid(self.pLVec) then ptr_free,self.pLVec
  if ptr_valid(self.Tsamp) then ptr_free,self.tsamp
  self.NSeg = 0
  self.randommsq = 0
  return & end

;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|

pro jhMSQ::cleanup
;cleanup routine for the object; free data
self.reset
return
end

;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|

function jhMSQ::init, FFS         $ ;DATA = [k0, ... kn] 'photon count data'
                      ,fkHz    =fkHz      $ ;   set sampling frequency in kHz
                      ,fHz     =fHz       $ ;                      or  in  Hz
                      ,Tms     =tms       $ ;                      or as a samplign time in ms
                      ,Tus     =tus       $  ;                                         or in us                      
                      ,randomassignment=randomassignment   $
                      ,quiet=quiet
;check FFS data valid

if keyword_set(fkhz) then self.Tsamp=ptr_new(1d/fkhz)
if keyword_set(fhz) then self.tsamp=ptr_new(1000d/fhz)
if keyword_set(tms) then self.tsamp=ptr_new(tms)
if keyword_set(tus) then self.tsamp=ptr_new(tus*1000)
self.pk=ptr_new(FFS)
self.exe,randomassignment=randomassignment,quiet=quiet

return,1 ;TODO return 0 or 1 depending on whether read was successful

end

;{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|{{:|

pro jhMSQ__define
;defines the fields for the object.
void = {jhMSQ, inherits IDL_object $
  ,pMSQ : ptr_new() $
  ,pVar : ptr_new() $  
  ,pVarexp : ptr_new() $
  ,pLVec : ptr_new() $
  ,Tsamp : ptr_new() $
  ,pk : ptr_new() $
  ,NSeg : 0 $
  ,randommsq : 0 $
  }                       
return
end