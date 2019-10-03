function createbootstrapmsqarray,o,narrays=narrays,dualcolor=dualcolor,multimsqarray=multimsqarray $;
  ,multilvecarray=multilvecarray,multivariancearray=multivariancearray,lvec=lvec $
  ,multitheoryvariance=multitheoryvariance,base=base

if ~keyword_set(narrays) then narr=100 else narr=narrays
nmsqs=n_elements(o)
for j=0,nmsqs-1 do o[j].exe,/randomassign,/quiet,base=base
lvecmin=n_elements(o[0].lvec())

for j=0,nmsqs-1 do lvecmin=min([lvecmin,n_elements(o[j].lvec())])

lvec=(o[0].lvec())[0:lvecmin-1]
if keyword_set(dualcolor) then multimsqarray=dblarr(nmsqs,narr,lvecmin,2) else multimsqarray=dblarr(nmsqs,narr,lvecmin)
multilvecarray=multimsqarray
multitheoryvariance=multimsqarray
multiselfvariance=multimsqarray
tic
for j=0,nmsqs-1 do begin
  print,j
  toc
  if ~keyword_set(dualcolor) then begin
    karr=o[j].k()
    ktot=[0,total(karr,/cumulative,/double)]
    ksq=[0,total(karr[0:-2]*karr[1:*],/cumulative,/double)]
    nk=n_elements(karr)
  endif else begin
    kaarr=o[j].ka()
    kbarr=o[j].kb()
    katot=total(kaarr,/cumulative,/double)
    kbtot=total(kbarr,/cumulative,/double)
    kabsq=total(kaarr[0:-2]*kbarr[1:*],/cumulative,/double)
    kbasq=total(kbarr[0:-2]*kaarr[1:*],/cumulative,/double)
    kabsq=mean([[kabsq],[kbasq]],dim=2)
    kbbsq=total(kbarr[0:-2]*kbarr[1:*],/cumulative,/double)
    nk=n_elements(kaarr)
  endelse
  for k=0,narr-1 do begin
    multilvecarray[j,k,*]=lvec
    for l=0,lvecmin-1 do begin
      lvecval=double(lvec[l])
      u=floor(nk/lvecval)
      
      randvals=floor(randomu(seed,u) * (nk - lvecval))
      if keyword_set(dualcolor) then begin
        valsab=(kabsq[randvals+lvecval-1]-kabsq[randvals])/(kbtot[randvals+lvecval]-kbtot[randvals+1])-(katot[randvals+lvecval-1]-katot[randvals])/(lvecval-1)
        valsbb=(kbbsq[randvals+lvecval-1]-kbbsq[randvals])/(kbtot[randvals+lvecval]-kbtot[randvals+1])-(kbtot[randvals+lvecval-1]-kbtot[randvals])/(lvecval-1)
        multimsqarray[j,k,l,0]=mean(valsab)
        multitheoryvariance[j,k,l,0]=1./nk/2
        multiselfvariance[j,k,l,0]=variance(valsab)/n_elements(valsab)/2.
        multimsqarray[j,k,l,1]=mean(valsbb)+1./lvecval
        multitheoryvariance[j,k,l,1]=1./nk
        multiselfvariance[j,k,l,1]=variance(valsbb)/n_elements(valsbb)
        if finite(multimsqarray[j,k,l,0],/nan) or finite(multimsqarray[j,k,l,1],/nan) then begin
          print,'NAN encountered'
          l=l-1
          continue
        endif
      endif else begin
        vals=(ksq[randvals+lvecval-1]-ksq[randvals])/(ktot[randvals+lvecval]-ktot[randvals+1])-(ktot[randvals+lvecval-1]-ktot[randvals])/(lvecval-1)
        multimsqarray[j,k,l]=mean(vals)+1./lvecval
        multitheoryvariance[j,k,l]=1./nk
        multiselfvariance[j,k,l]=variance(vals)/n_elements(vals)
        if finite(multimsqarray[j,k,l],/nan) then begin
          print,'NAN encountered'
          vals=vals[where(finite(vals) eq 1)]
          multimsqarray[j,k,l]=mean(vals)+1./lvecval
          multitheoryvariance[j,k,l]=1./nk
          multiselfvariance[j,k,l]=variance(vals)/n_elements(vals)
          continue
        endif
      endelse
    end
  end
end

multivariancearray=variance(multimsqarray,dim=2)
return,{lvec:lvec,multilvecarray:multilvecarray,multimsqarray:multimsqarray,multivariancearray:multivariancearray,multitheoryvariance:multitheoryvariance,multiselfvariance:multiselfvariance}

end