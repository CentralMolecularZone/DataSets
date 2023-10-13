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

from bs4 import BeautifulSoup

resp = requests.get('https://irsa.ipac.caltech.edu/data/SPITZER/MIPSGAL/images/mosaics24/')
soup = BeautifulSoup(resp.text)
baseurl = 'https://irsa.ipac.caltech.edu/data/SPITZER/MIPSGAL/images/mosaics24/'

for url in soup.findAll('a'):
    MGname = url.attrs['href'].strip().lstrip()
    if any((MGname.startswith(x)
        for x in ('MG000', 'MG001', 'MG360', 'MG359',
                  'MG358', 'MG357', 'MG002'))) and not os.path.exists(MGname):
        if not any((z in url.attrs['href'] for z in ('std', 'covg', 'maskcube'))):
            urltoget = f'{baseurl}/{MGname}'
            print(f"Downloading {MGname} from {urltoget}")
            resp = requests.get(urltoget, stream=True)
            with open(MGname, 'wb') as fh:
                fh.write(resp.content)

hdus = [fits.open(x) for x in glob.glob("MG*.fits")]

# galactic
outheader = fits.Header.fromtextfile('mipsgal_cmz.header')
wcs_out = wcs.WCS(outheader)
# note that this is not guaranteed to always exist but it does here
# also the axis order is strictly wrong
shape_out = wcs_out._naxis[::-1]

array_line, footprint = reproject_and_coadd([h[0] for h in hdus], wcs_out, shape_out=shape_out, reproject_function=reproject_interp)
fits.PrimaryHDU(data=array_line, header=wcs_out.to_header()).writeto('gc_mosaic_MIPSGAL_gal.fits', overwrite=True)

# celestial
wcs_out, shape_out = find_optimal_celestial_wcs([h[0] for h in hdus])

array_line, footprint = reproject_and_coadd([h[0] for h in hdus], wcs_out, shape_out=shape_out, reproject_function=reproject_interp)
fits.PrimaryHDU(data=array_line, header=wcs_out.to_header()).writeto('gc_mosaic_MIPSGAL.fits', overwrite=True)
