function hoppingrandomwalk,punstick=punstick,pstick=pstick

;standard parameters, these may be changed

if ~keyword_set(punstick) then punstick=0.001 ;probability of immobile site becoming mobile
if ~keyword_set(pstick) then pstick=0.005 ;probability of mobile becoming immobile
density=0.01 ;density in units of number/grid area
griddim=800l ;size of grid in grid units
tlength=20000*600l ;length of simulation in sampling times
brightness=300./20000 ;counts per sampling time molecular brightness
wo=50. ;beam waist in grid units

scanfreqratio=20l ;scan time/sampling time (must be > 1)
scanlength=wo*4. ;peak to peak amplitude, best to use in terms of wo
waveformhalf=dindgen(scanfreqratio/2.)/(scanfreqratio/2.)*scanlength ;up sweep of PSF center
waveformwhole=[waveformhalf,reverse(waveformhalf)+2./scanfreqratio*scanlength]-scanlength/2. ;up and down sweep - repeated to form full scan

nps=round(density*griddim^2) ;number of particles to simulate
startlocs=(randomu(seed,2,nps)*griddim-griddim/2.) ;randomize initial locations of particles on the grid
boundratio=dblarr(tlength) ;this will keep track of ratio bound - use as sanity check for stable values
stuck=intarr(nps) ;1 if stuck, 0 if unstuck
samplex=dblarr(2,tlength) ;track individual particle as sanity check
farr=dblarr(tlength) ;intensity readout, this will be the output
samplex[*,0]=startlocs[*,0] ;initialize sample particle track
xs=startlocs ;initialize particle locations
psfx=waveformwhole[0] & psfy=0 ;initialize PSF center location
farr[0]=brightness*total(exp(-4*((xs[0,*]-psfx)^2+(xs[1,*]-psfy)^2)/wo^2)) ;calculate initial intensity at t = 0
tic
for k=1,tlength-1 do begin ;iterating through steps
  if (k eq 20000) then toc ;every 1 s of simulation, print the time the simulation has taken
  
  xory=round(randomu(seed,nps)) ;choose which axis for 2D random walk
  moveval=round(randomu(seed,nps))*2-1 ;choose which direction along axis
  
  ;Perform stochastic binding/unbinding w/ probabilities defined earlier
  stuckcheck=randomu(seed,nps)
  stucknew=stuck
  stucknew[where((stuckcheck lt pstick) and (stuck eq 0),/null)]=1
  stucknew[where((stuckcheck lt punstick) and (stuck eq 1),/null)]=0
  stuck=stucknew
  
  ;Move particles as appropriate
  moveval=moveval*(1-stuck) ;only move particles that aren't stuck
  boundratio[k]=1./max([avg(1-stuck),.001/nps])-1 ;calculate the ratio of bound particles
  xs[0,*]=xs[0,*]+abs(xory-1)*moveval & xs[1,*]=xs[1,*]+xory*moveval ;update particle locations
  highvals=where(xs gt griddim/2.,/null) & lowvals=where(xs lt griddim/(-2.),/null) ;check for extreme values
  if ~isa(highvals,/null) then xs[highvals]=xs[highvals]-griddim/2. ;if there are extreme values, mirror locations
  if ~isa(lowvals,/null) then xs[lowvals]=xs[lowvals]+griddim/2.
  samplex[*,k]=xs[*,0] ;save sample trajectory
  
  psfx=waveformwhole[k mod scanfreqratio] ;move PSF center
  farr[k]=brightness*total(exp(-4*((xs[0,*]-psfx)^2+(xs[1,*]-psfy)^2)/wo^2)) ;calculate intensity
end


return,farr

end