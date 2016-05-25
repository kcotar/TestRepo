pro calculate_multiple_indices, ref_img, mask_bad, $ ;reflectance image, image mask of bad pixels
                                ndvi=ndvi, evi=evi, afri=afri, arvi=arvi, msavi=msavi, ndwi=ndwi, ndsi=ndsi, savi=savi, ndgi=ndgi, $ ;indices
                                fapar=fapar, lai=lai, fcover=fcover, $  ;biophysical paramete
                                probav=probav, re=re, spot_vgt=spot_vgt, $  ;select apropriate sensor
                                mask_value=mask_value, $  ;set value for bad pixels
                                proba_to_spot=proba_to_spot

;check if any of the sensors is selected

;!!ADD ARG_PRESENT to check which indices should be calculated

if keyword_set(proba_to_spot) then begin
  proba_to_spot_n = [0.000, 0.0075, 0.0076, 0.0049]
  proba_to_spot_k = [1.000, 0.8600, 0.8790, 0.9070]
  
  for i_b=0, n_elements(proba_to_spot_n)-1 do begin
    ref_img[*,*,i_b] = ref_img[*,*,i_b] * proba_to_spot_k[i_b] + proba_to_spot_n[i_b]
  endfor
  
endif


if keyword_set(probav) or keyword_set(spot_vgt) then begin
  blu = 0
  red = 1
  nir = 2
  swr = 3
endif

if keyword_set(re) then begin
  blu = 0
  grn = 1
  red = 2
  rdg = 3
  nir = 4
endif



;get bad pixels from mask
if n_elements(mask_bad) gt 0 then idx_bad = where(mask_bad ge 1, n_bad) else n_bad=0

;set mask value
if n_elements(mask_value) eq 0 then mask_value = !Values.F_NaN

 
;normalised diference vegetation indeks                      
ndvi = (ref_img[*,*,nir]-ref_img[*,*,red])/(ref_img[*,*,nir]+ref_img[*,*,red])

;soil adjusted vegetation ineks
l_savi = 0.5
savi = (ref_img[*,*,nir]-ref_img[*,*,red])/(ref_img[*,*,nir]+ref_img[*,*,red]+l_savi)*(1.+l_savi)

;modified soil and atmospherically resistant vegetation indeks
msavi = (2.*ref_img[*,*,nir]+1.-sqrt((2.*ref_img[*,*,nir]+1.)^2.-8.*(ref_img[*,*,nir]-ref_img[*,*,red])))/2.

if n_elements(blu) gt 0 then begin
  ;enhanced vegetation indeks
  g_evi = 2.5
  c1_evi = 6.
  c2_evi = 7.5
  l_evi = 1.
  evi = g_evi*(ref_img[*,*,nir]-ref_img[*,*,red])/(ref_img[*,*,nir] + c1_evi*ref_img[*,*,red] - c2_evi*ref_img[*,*,blu] + l_evi)
  
  ;atmospherically resistant vegetation indeks
  gamma = 1.
  red_blu = ref_img[*,*,red]-gamma*(ref_img[*,*,blu]-ref_img[*,*,red])
  arvi = (ref_img[*,*,nir]-red_blu)/(ref_img[*,*,nir]+red_blu)
endif

if n_elements(grn) gt 0 then begin
  ndgi = (ref_img[*,*,grn]-ref_img[*,*,red])/(ref_img[*,*,grn]+ref_img[*,*,red])
endif

if n_elements(swr) gt 0 then begin
  ;atmospherically free indeks
  k_afri = 0.66
  afri = (ref_img[*,*,nir]-k_afri*ref_img[*,*,swr])/(ref_img[*,*,nir]+k_afri*ref_img[*,*,swr])
  
  ;normalized water indeks
  ndwi = (ref_img[*,*,swr]-ref_img[*,*,red])/(ref_img[*,*,swr]+ref_img[*,*,red])
  
  ;normalized snow indeks
  ndsi = (ref_img[*,*,swr]-ref_img[*,*,nir])/(ref_img[*,*,swr]+ref_img[*,*,nir])
endif


;biophysical parameters
if keyword_set(probav) or keyword_set(spot_vgt) then begin
  if arg_present(lai) then lai = CMAPPLY('USER:lai_calc', ref_img, 3)
  if arg_present(fapar) then fapar = CMAPPLY('USER:fapar_calc', ref_img, 3)
  if arg_present(fcover) then fcover = CMAPPLY('USER:fcover_calc', ref_img, 3)
endif



;maskout bad values from products
if n_bad gt 0 then begin
  ndvi[idx_bad] = mask_value  
  msavi[idx_bad] = mask_value
  savi[idx_bad] = mask_value
  if n_elements(blu) gt 0 then begin
    evi[idx_bad] = mask_value
    arvi[idx_bad] = mask_value
  endif
  if n_elements(grn) gt 0 then begin
    ndgi[idx_bad] = mask_value
  endif
  if n_elements(swr) gt 0 then begin
    afri[idx_bad] = mask_value
    ndwi[idx_bad] = mask_value
    ndsi[idx_bad] = mask_value
  endif
  if keyword_set(probav) then begin
    if arg_present(lai) then lai[idx_bad] = mask_value
    if arg_present(fapar) then fapar[idx_bad] = mask_value
    if arg_present(fcover) then fcover[idx_bad] = mask_value
  endif
endif

end