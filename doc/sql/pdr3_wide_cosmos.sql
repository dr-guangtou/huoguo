SELECT

	-- Basic information
	f1.object_id, f1.parent_id, f1.ra, f1.dec, f1.tract, f1.patch,
	
	-- Galactic extinction for correction
	f1.a_g, f1.a_r, f1.a_i, f1.a_z, f1.a_y,

	-- Number of images contributing at the center
	f1.g_inputcount_value as g_input_count,
	f1.r_inputcount_value as r_input_count,
	f1.i_inputcount_value as i_input_count,
	f1.z_inputcount_value as z_input_count,
	f1.y_inputcount_value as y_input_count,

	-- Flag for measurements from the HSC-R2/I2 filter
	f1.merge_measurement_r2,
	f1.merge_measurement_i2,
	c.corr_rmag,
	c.corr_imag,
	
	-- CModel photometry
	-- 1. Exponential component photometry (useful for dwarfs)
	f1.g_cmodel_exp_flux, 
	f1.r_cmodel_exp_flux, 
	f1.i_cmodel_exp_flux, 
	f1.z_cmodel_exp_flux, 
	f1.y_cmodel_exp_flux, 
	f1.g_cmodel_exp_fluxerr as g_cmodel_exp_err,
	f1.r_cmodel_exp_fluxerr as r_cmodel_exp_err,
	f1.i_cmodel_exp_fluxerr as i_cmodel_exp_err,
	f1.z_cmodel_exp_fluxerr as z_cmodel_exp_err,
	f1.y_cmodel_exp_fluxerr as y_cmodel_exp_err,
 
	-- 2. de Vaucoulerus component photometry
	f1.g_cmodel_dev_flux, 
	f1.r_cmodel_dev_flux, 
	f1.i_cmodel_dev_flux, 
	f1.z_cmodel_dev_flux, 
	f1.y_cmodel_dev_flux, 
	f1.g_cmodel_dev_fluxerr as g_cmodel_dev_err,
	f1.r_cmodel_dev_fluxerr as r_cmodel_dev_err,
	f1.i_cmodel_dev_fluxerr as i_cmodel_dev_err,
	f1.z_cmodel_dev_fluxerr as z_cmodel_dev_err,
	f1.y_cmodel_dev_fluxerr as y_cmodel_dev_err,
	
	-- 3. CModel photometry
	f1.g_cmodel_flux, 
	f1.r_cmodel_flux, 
	f1.i_cmodel_flux, 
	f1.z_cmodel_flux, 
	f1.y_cmodel_flux, 
	f1.g_cmodel_fluxerr as g_cmodel_err,
	f1.r_cmodel_fluxerr as r_cmodel_err,
	f1.i_cmodel_fluxerr as i_cmodel_err,
	f1.z_cmodel_fluxerr as z_cmodel_err,
	f1.y_cmodel_fluxerr as y_cmodel_err,

	-- 4. Flags for CModel photometry
	-- If the flag is set to True, it means the photometry is not reliable
	f1.g_cmodel_flag,
	f1.z_cmodel_flag,
	f1.y_cmodel_flag,
	
	-- PSF photometry
	f2.g_psfflux_flux as g_psf_flux, 
	f2.r_psfflux_flux as r_psf_flux, 
	f2.i_psfflux_flux as i_psf_flux, 
	f2.z_psfflux_flux as z_psf_flux, 
	f2.y_psfflux_flux as y_psf_flux, 
	f2.g_psfflux_fluxerr as g_psf_err,
	f2.r_psfflux_fluxerr as r_psf_err,
	f2.i_psfflux_fluxerr as i_psf_err,
	f2.z_psfflux_fluxerr as z_psf_err,
	f2.y_psfflux_fluxerr as y_psf_err,
 
	-- PSF photometry flag (for quality cut later)
	f2.g_psfflux_flag as g_psf_flag,
	f2.r_psfflux_flag as r_psf_flag,
	f2.i_psfflux_flag as i_psf_flag,
	f2.z_psfflux_flag as z_psf_flag,
	f2.y_psfflux_flag as y_psf_flag,
	
	-- PSF-corrected aperture photometry
	-- 2_15: 1.1 arcsec seeing; 1.5 arcsec diamter aperture
	f4.g_convolvedflux_2_15_flux, 
	f4.r_convolvedflux_2_15_flux, 
	f4.i_convolvedflux_2_15_flux, 
	f4.z_convolvedflux_2_15_flux, 
	f4.y_convolvedflux_2_15_flux, 
 
	f4.g_convolvedflux_2_15_fluxerr as g_convolvedflux_2_15_err,
	f4.r_convolvedflux_2_15_fluxerr as r_convolvedflux_2_15_err,
	f4.i_convolvedflux_2_15_fluxerr as i_convolvedflux_2_15_err,
	f4.z_convolvedflux_2_15_fluxerr as z_convolvedflux_2_15_err,
	f4.y_convolvedflux_2_15_fluxerr as y_convolvedflux_2_15_err,
 
	f4.g_convolvedflux_2_15_flag,
	f4.r_convolvedflux_2_15_flag,
	f4.i_convolvedflux_2_15_flag,
	f4.z_convolvedflux_2_15_flag,
	f4.y_convolvedflux_2_15_flag,
	
	-- 3_20: 1.3 arcsec seeing; 2.0 arcsec diamter aperture
	f4.g_convolvedflux_3_20_flux, 
	f4.r_convolvedflux_3_20_flux, 
	f4.i_convolvedflux_3_20_flux, 
	f4.z_convolvedflux_3_20_flux, 
	f4.y_convolvedflux_3_20_flux, 
 
	f4.g_convolvedflux_3_20_fluxerr as g_convolvedflux_3_20_err,
	f4.r_convolvedflux_3_20_fluxerr as r_convolvedflux_3_20_err,
	f4.i_convolvedflux_3_20_fluxerr as i_convolvedflux_3_20_err,
	f4.z_convolvedflux_3_20_fluxerr as z_convolvedflux_3_20_err,
	f4.y_convolvedflux_3_20_fluxerr as y_convolvedflux_3_20_err,
 
	f4.g_convolvedflux_3_20_flag,
	f4.r_convolvedflux_3_20_flag,
	f4.i_convolvedflux_3_20_flag,
	f4.z_convolvedflux_3_20_flag,
	f4.y_convolvedflux_3_20_flag,
	
	-- PSF-corrected aperture photometry **before deblending**
	-- 2_15: 1.1 arcsec seeing; 1.5 arcsec diamter aperture
	f5.g_undeblended_convolvedflux_2_15_flux, 
	f5.r_undeblended_convolvedflux_2_15_flux, 
	f5.i_undeblended_convolvedflux_2_15_flux, 
	f5.z_undeblended_convolvedflux_2_15_flux, 
	f5.y_undeblended_convolvedflux_2_15_flux, 
 
	f5.g_undeblended_convolvedflux_2_15_fluxerr as g_undeblended_convolvedflux_2_15_err,
	f5.r_undeblended_convolvedflux_2_15_fluxerr as r_undeblended_convolvedflux_2_15_err,
	f5.i_undeblended_convolvedflux_2_15_fluxerr as i_undeblended_convolvedflux_2_15_err,
	f5.z_undeblended_convolvedflux_2_15_fluxerr as z_undeblended_convolvedflux_2_15_err,
	f5.y_undeblended_convolvedflux_2_15_fluxerr as y_undeblended_convolvedflux_2_15_err,
 
	f5.g_undeblended_convolvedflux_2_15_flag,
	f5.r_undeblended_convolvedflux_2_15_flag,
	f5.i_undeblended_convolvedflux_2_15_flag,
	f5.z_undeblended_convolvedflux_2_15_flag,
	f5.y_undeblended_convolvedflux_2_15_flag,
	
	-- 3_20: 1.3 arcsec seeing; 2.0 arcsec diamter aperture
	f5.g_undeblended_convolvedflux_3_20_flux, 
	f5.r_undeblended_convolvedflux_3_20_flux, 
	f5.i_undeblended_convolvedflux_3_20_flux, 
	f5.z_undeblended_convolvedflux_3_20_flux, 
	f5.y_undeblended_convolvedflux_3_20_flux, 
 
	f5.g_undeblended_convolvedflux_3_20_fluxerr as g_undeblended_convolvedflux_3_20_err,
	f5.r_undeblended_convolvedflux_3_20_fluxerr as r_undeblended_convolvedflux_3_20_err,
	f5.i_undeblended_convolvedflux_3_20_fluxerr as i_undeblended_convolvedflux_3_20_err,
	f5.z_undeblended_convolvedflux_3_20_fluxerr as z_undeblended_convolvedflux_3_20_err,
	f5.y_undeblended_convolvedflux_3_20_fluxerr as y_undeblended_convolvedflux_3_20_err,
 
	f5.g_undeblended_convolvedflux_3_20_flag,
	f5.r_undeblended_convolvedflux_3_20_flag,
	f5.i_undeblended_convolvedflux_3_20_flag,
	f5.z_undeblended_convolvedflux_3_20_flag,
	f5.y_undeblended_convolvedflux_3_20_flag,

	-- SDSS Shape without PSF correction (Using i-band; can use others too)
	f2.i_sdssshape_shape11 as i_sdss_shape_11,
	f2.i_sdssshape_shape12 as i_sdss_shape_12,
	f2.i_sdssshape_shape22 as i_sdss_shape_22,
	f2.i_sdssshape_shape11err as i_sdss_shape_11_err,
	f2.i_sdssshape_shape12err as i_sdss_shape_12_err,
	f2.i_sdssshape_shape22err as i_sdss_shape_22_err,
	
	-- Shape of the CModel model (i & r-band)
	m.i_cmodel_exp_ellipse_11, 
	m.i_cmodel_exp_ellipse_22, 
	m.i_cmodel_exp_ellipse_12,
	m.i_cmodel_ellipse_11, 
	m.i_cmodel_ellipse_22, 
	m.i_cmodel_ellipse_12,
	
	m.r_cmodel_exp_ellipse_11, 
	m.r_cmodel_exp_ellipse_22, 
	m.r_cmodel_exp_ellipse_12,
	m.r_cmodel_ellipse_11, 
	m.r_cmodel_ellipse_22, 
	m.r_cmodel_ellipse_12,

	-- Extendedness of the object
	f1.g_extendedness_value,
	f1.z_extendedness_value,
	f1.y_extendedness_value,
	f1.g_extendedness_flag,
	f1.r_extendedness_flag,
	f1.i_extendedness_flag,
	f1.z_extendedness_flag,
	f1.y_extendedness_flag,
	
	-- Flags for later selection
	-- 1. The general failure flag
	f1.g_pixelflags,
	f1.r_pixelflags,
	f1.i_pixelflags,
	f1.z_pixelflags,
	f1.y_pixelflags,
	
	-- 2. Saturated or interpolated pixels on the footprint (not center)
	f1.g_pixelflags_saturated,
	f1.r_pixelflags_saturated,
	f1.i_pixelflags_saturated,
	f1.z_pixelflags_saturated,
	f1.y_pixelflags_saturated,
	f1.g_pixelflags_interpolated,
	f1.r_pixelflags_interpolated,
	f1.i_pixelflags_interpolated,
	f1.z_pixelflags_interpolated,
	f1.y_pixelflags_interpolated,

	-- 3. Other pixel flags
	f1.g_pixelflags_bad,
	f1.r_pixelflags_bad,
	f1.i_pixelflags_bad,
	f1.z_pixelflags_bad,
	f1.y_pixelflags_bad,
	f1.g_pixelflags_suspectcenter,
	f1.r_pixelflags_suspectcenter,
	f1.i_pixelflags_suspectcenter,
	f1.z_pixelflags_suspectcenter,
	f1.y_pixelflags_suspectcenter,
	f1.g_pixelflags_clippedcenter,
	f1.r_pixelflags_clippedcenter,
	f1.i_pixelflags_clippedcenter,
	f1.z_pixelflags_clippedcenter,
	f1.y_pixelflags_clippedcenter,

	-- 4. Bright object masks
	f1.g_pixelflags_bright_object,
	f1.r_pixelflags_bright_object,
	f1.i_pixelflags_bright_object,
	f1.z_pixelflags_bright_object,
	f1.y_pixelflags_bright_object,
	
	-- Mizuki photo-z information 
	p1.photoz_mean as pz_mean_mizuki, 
	p1.photoz_best as pz_best_mizuki, 
	p1.photoz_conf_mean as pz_conf_mean_mizuki, 
	p1.photoz_conf_best as pz_conf_best_mizuki, 
	p1.photoz_risk_mean as pz_risk_mean_mizuki, 
	p1.photoz_risk_best as pz_risk_best_mizuki,
	p1.photoz_std_mean as pz_std_mean_mizuki, 
	p1.photoz_std_best as pz_std_best_mizuki,
	p1.photoz_err68_min as pz_err68_min_mizuki, 
	p1.photoz_err68_max as pz_err68_max_mizuki, 
	p1.stellar_mass as mstar_mizuki,
	p1.stellar_mass_err68_min as mstar_min_mizuki,
	p1.stellar_mass_err68_max as mstar_max_mizuki,
	
	-- DNNZ photo-z information 
	p2.photoz_mean as pz_mean_dnnz, 
	p2.photoz_best as pz_best_dnnz, 
	p2.photoz_conf_mean as pz_conf_mean_dnnz, 
	p2.photoz_conf_best as pz_conf_best_dnnz, 
	p2.photoz_risk_mean as pz_risk_mean_dnnz, 
	p2.photoz_risk_best as pz_risk_best_dnnz,
	p2.photoz_std_mean as pz_std_mean_dnnz, 
	p2.photoz_std_best as pz_std_best_dnnz,
	p2.photoz_err68_min as pz_err68_min_dnnz, 
	p2.photoz_err68_max as pz_err68_max_dnnz, 	

	-- DEMP photo-z information 
	p3.photoz_mean as pz_mean_demp, 
	p3.photoz_best as pz_best_demp, 
	p3.photoz_conf_mean as pz_conf_mean_demp, 
	p3.photoz_conf_best as pz_conf_best_demp, 
	p3.photoz_risk_mean as pz_risk_mean_demp, 
	p3.photoz_risk_best as pz_risk_best_demp,
	p3.photoz_std_mean as pz_std_mean_demp, 
	p3.photoz_std_best as pz_std_best_demp,
	p3.photoz_err68_min as pz_err68_min_demp, 
	p3.photoz_err68_max as pz_err68_max_demp	

FROM
	pdr3_wide.forced as f1
	LEFT JOIN pdr3_wide.forced2 as f2 USING (object_id)
	LEFT JOIN pdr3_wide.forced4 as f4 USING (object_id)
	LEFT JOIN pdr3_wide.forced5 as f5 USING (object_id)
	LEFT JOIN pdr3_wide.meas as m USING (object_id)
	LEFT JOIN pdr3_wide.mag_corr as c USING (object_id)
	LEFT JOIN pdr3_wide.photoz_mizuki as p1 USING (object_id)
	LEFT JOIN pdr3_wide.photoz_dnnz as p2 USING (object_id)
	LEFT JOIN pdr3_wide.photoz_demp as p3 USING (object_id)

WHERE
    	f1.isprimary
	-- Region
	-- This is for the COSMOS field
	AND boxSearch(coord, 149.22484, 150.81303, 1.541319, 2.87744)
	-- Photometric cut
	AND f1.i_cmodel_mag <= 23.5
	AND NOT f1.r_cmodel_flag
	AND NOT f1.i_cmodel_flag
	-- Extendedness cut 
	AND f1.r_extendedness_value > 0
	AND f1.i_extendedness_value > 0
	-- Full-depth full-color cut
	AND f1.g_inputcount_value >= 4
	AND f1.r_inputcount_value >= 4
	AND f1.i_inputcount_value >= 5
	AND f1.z_inputcount_value >= 5
	AND f1.y_inputcount_value >= 5
	-- Not failed at finding the center
	AND NOT f2.g_sdsscentroid_flag
	AND NOT f2.r_sdsscentroid_flag
	AND NOT f2.i_sdsscentroid_flag
	AND NOT f2.z_sdsscentroid_flag
	AND NOT f2.y_sdsscentroid_flag
	-- The object's center is not outside the image
	AND NOT f1.g_pixelflags_edge    
	AND NOT f1.r_pixelflags_edge
	AND NOT f1.i_pixelflags_edge
	AND NOT f1.z_pixelflags_edge
	AND NOT f1.y_pixelflags_edge
	-- Not saturated at the center  
	AND NOT f1.g_pixelflags_saturatedcenter
    AND NOT f1.r_pixelflags_saturatedcenter
    AND NOT f1.i_pixelflags_saturatedcenter
    AND NOT f1.z_pixelflags_saturatedcenter
    AND NOT f1.y_pixelflags_saturatedcenter
	-- The center is not interpolated
	AND NOT f1.g_pixelflags_interpolatedcenter
	AND NOT f1.r_pixelflags_interpolatedcenter
	AND NOT f1.i_pixelflags_interpolatedcenter
	AND NOT f1.z_pixelflags_interpolatedcenter
	AND NOT f1.y_pixelflags_interpolatedcenter
	-- The center is not affected by a cosmic ray
	AND NOT f1.g_pixelflags_crcenter
    AND NOT f1.r_pixelflags_crcenter
    AND NOT f1.i_pixelflags_crcenter
    AND NOT f1.z_pixelflags_crcenter
    AND NOT f1.y_pixelflags_crcenter
;