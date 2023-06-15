'''
Generate WCS for a given DECADE Tile 
and make cutouts using TIFF image.

'''

__author__ = 'Kai Herron'

import numpy as np
import pandas as pd
import easyaccess as ea
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import shutil
import os
from PIL import Image
import logging



QUERY = '''
SELECT C.CRPIX1, CRPIX2, C.CRVAL1, C.CRVAL2, C.CD1_1,
       C.CD1_2, C.CD2_1, C.CD2_2, C.CTYPE1, C.CTYPE2,
       C.DECC1, C.DECC2, C.DECC3, C.DECC4,
       C.DECCMAX, C.DECCMIN, C.DEC_CENT, C.DEC_SIZE,
       C.NAXIS1, C.NAXIS2, C.PIXELSCALE,
       C.RAC1, C.RAC2, C.RAC3, C.RAC4,
       C.RACMAX, C.RACMIN, C.RA_CENT, C.RA_SIZE,
       C.UDECMAX, C.UDECMIN, C.URAMAX, C.URAMIN
FROM COADDTILE_GEOM C WHERE C.TILENAME = '{tilename}'
'''
TIFF_QUERY = '''
select
    fai.path,
    c.FILENAME
from
    COADD c,
    PROCTAG t,
    file_archive_info fai
where
    c.filetype = 'coadd' and
    t.tag = 'DR3_1_1' and
    t.pfw_attempt_id = c.pfw_attempt_id and
    fai.filename = c.filename and
    c.tilename = 'DES1110+0001' and
    c.band = 'g'
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-t','--tile',
                        help='DECADE tilename')
    parser.add_argument('-r', '--ra', help='Right Ascension (Decimal Degrees)')
    parser.add_argument('-d', '--dec', help='Declination (Decimal Degrees)')
    parser.add_argument('-o', '--objid',type=int,help='COADD_OBJECT_ID')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='output verbosity')
    args = parser.parse_args()
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(filename='wcs.log', filemode='w', level=logging.DEBUG)
    logging.getLogger().setLevel(level)
    
    # establish easyaccess connection w/ decade
    logging.debug('Establishing connection to DECADE...')
    conn = ea.connect(section='delve')
    logging.debug('Formatting SQL query...')
    tile = args.tile
    logging.debug(f'TILENAME = {tile}')
    logging.debug(f'DTYPE = {type(tile)}')
    query = QUERY.format(tilename=tile)
    logging.debug('Submitting Query:')
    logging.debug(query)
    logging.debug('Querying DECADE...')
    df = conn.query_to_pandas(query)
    logging.debug('Displaying dataframe..')
    cols = df.columns
    for i in cols:
        logging.debug(f'{i} = {df[i][0]}')
    
    # make a header to feed to wcs
    logging.debug('Generating generic header for WCS...')
    hdu = pyfits.PrimaryHDU()
    hdr = hdu.header
    
    logging.debug('Adding DECADE info to header...')
    
    for i in cols:
        hdr[i] = df[i][0]
        logging.debug(f'Adding {i} info...')
    wcs = WCS(hdr)
    logging.debug('WCS generated...')
    
    # read in RAs and Decs
    logging.debug('Reading RAs and Decs...')
    ra = args.ra
    dec = args.dec
    
    logging.debug('Converting RA and Dec to pixel coords...')
    pix_x, pix_y = wcs.world_to_pixel(ra,dec)
    inc_x, inc_y = wcs.world_to_pixel(ra+0.00694444,dec+0.00694444)
    del_x = (inc_x - pix_x) * 2
    del_y = (inc_y - pix_y) * 2
    
    min_x = round(pix_x - del_x)
    max_x = round(pix_x + del_x)
    
    min_y = round(pix_y - del_y)
    max_y = round(pix_y + del_y)
    
    logging.debug('Generating cursor for TIFF_QUERY...')
    cursor=conn.cursor()
    logging.debug('Executing Query...')
    QQ=cursor.execute(TIFF_QUERY)
    logging.debug('Fetching data...')
    paths = QQ.fetchall()
    raw_path = paths[0][0]

    sys_run = raw_path[27:32]
    tilename = raw_path[33:45]
    run = raw_path[46:49]
    logging.debug(f'RAW_PATH: {raw_path}')
    path=raw_path[:-5] + 'qa/{}_{}{}_irg.tiff'.format( tilename,sys_run,run)
    logging.debug(f'PATH_TO_TIFF: {path}')
    logging.debug('Downloading TIFF Image...')
    os.system('wget --user=decade --password=decaFil3s https://decade.ncsa.illinois.edu/deca_archive/%s  -P tiff' % path)
    logging.debug('Reading TIFF image file...')
    im = Image.open('tiff/{}_{}{}_irg.tiff'.format(tilename,sys_run,run))
    
    arr = np.asarray(im)
    cutout_arr = arr[min_x:max_x, 
                     min_y:max_y,
                     :]
    logging.debug('Cutout generated. Saving...')
    
    new_im = Image.fromarray(cutout_arr)
    new_im.save(f'cutouts/{args.objid}.png')
    logging.debug(f'{args.objid} processed.')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    