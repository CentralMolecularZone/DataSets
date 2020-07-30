import os
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling, time; from astropy.io import fits, ascii
import reproject
import glob
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject.mosaicking import reproject_and_coadd
from reproject import reproject_interp

import requests
from bs4 import BeautifulSoup

from astropy.utils.console import ProgressBar

for glon in ProgressBar(range(-25, 26, 2)):

    response = requests.get(f'http://miris.kasi.re.kr/cgi-bin/mirisobs.cgi?type=coord&coordinates={glon}+0&frame=galactic&filter=PAAL&filter=PAAC')
    response.raise_for_status()
    soup = BeautifulSoup(response.text, parser='html.parser')
    files = [x.attrs['href'] for x in soup.findAll('a') if '.fits' in x.attrs['href']]

    for fn in files:
        basename = os.path.basename(fn)
        if os.path.exists(basename):
            continue
        else:
            with open(basename, 'wb') as fh:
                res = requests.get(f'http://miris.kasi.re.kr/{fn}', stream=True)
                res.raise_for_status()
                fh.write(res.content)

hdus = [fits.open(x) for x in glob.glob("MS*.fits")]

# TODO: Change this to GLON-CAR/GLAT-CAR - the GLON-TAN/GLAT-TAN projection looks pretty wacky far from the CMZ!
wcs_out, shape_out = find_optimal_celestial_wcs([h[1] for h in hdus])

array_line, footprint = reproject_and_coadd([h[1] for h in hdus if h[0].header['OBS-FILT']=='PAAL'], wcs_out, shape_out=shape_out, reproject_function=reproject_interp)
array_cont, footprint = reproject_and_coadd([h[1] for h in hdus if h[0].header['OBS-FILT']=='PAAC'], wcs_out, shape_out=shape_out, reproject_function=reproject_interp)
fits.PrimaryHDU(data=array_line, header=wcs_out.to_header()).writeto('gc_mosaic_miris_line.fits', overwrite=True)
fits.PrimaryHDU(data=array_line - array_cont, header=wcs_out.to_header()).writeto('gc_mosaic_miris_line_minus_cont.fits', overwrite=True)
fits.PrimaryHDU(data=array_cont, header=wcs_out.to_header()).writeto('gc_mosaic_miris_cont.fits', overwrite=True)


line = fits.open('gc_mosaic_miris_line.fits')
cont = fits.open('gc_mosaic_miris_cont.fits')

# "best-fit" offset power-law fit to line vs cont
contsub = line[0].data - cont[0].data**1.1 * 0.35

fits.PrimaryHDU(data=contsub, header=wcs_out.to_header()).writeto('gc_mosaic_miris_line_minus_cont_scaled_pow1.1_x0p35.fits')
