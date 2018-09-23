""" collage of image cutouts """

from matplotlib.colors import LogNorm
import glob
import requests
from requests.exceptions import ConnectionError
from astroquery.mast import Observations
from astropy import coordinates
import astropy.units as u
from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
from get_names import *


def mast():
    """ access PS1 images through astroquery/MAST 
    the problem is that I don't think I can get image cutouts this way,
    only the full stacked image in the region.
    I don't want to have to do the cutting myself,
    so for now I'm going to use the PS1 image access webpage.
    """
    # radius should be 25 arcseconds, which is 0.00694 degrees
    obsTable = Observations.query_criteria(
            dataproduct_type=["image"],
            obs_collection=["PS1"],
            coordinates="331.675712 +36.208096", 
            radius="0.00694 deg")
    obsids = obsTable['obsid']

    dataProductsByID = Observations.get_product_list(obsids)
    Observations.download_products(
            obsids, curl_flag=True, # only the code to download it later, 
            mrp_only=True) # minimum recommended products
    # you get a .sh script that you can then run to download the products


def ps1(ra, dec):
    """ 
    use the PS1 image cutout server.
    
    documentation:
    http://hla.stsci.edu/fitscutcgi_interface.html    
    """

    # step 1: image list service
    s = requests.Session()
    www = "http://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    payload = {}
    payload["ra"] = ra
    payload["dec"] = dec
    payload["filters"] = 'i'
    payload["type"] = 'stack'
    payload["sep"] = 'comma'
    try:
        r = s.get(www, params=payload, timeout=(3.05,10))
    except requests.exceptions.RequestException as e: 
        print(e)
        sys.exit(1)
    print(r.url)
    if r.status_code != requests.codes.ok:
        print("status code" + str(r.status_code))
    s.close()
    raw = r.text.split("\n")
    header = np.array(raw[0].split(","))
    content = np.array(raw[1].split(","))
    sname = content[header=='shortname'][0]
    fname = content[header=='filename'][0]
    print(fname)
    
    # step 2: download the FITS file of the full image
    # server_name = "http://ps1images.stsci.edu"
    # url = server_name + fname
    # wget url

    # step 2: download the FITS file of the image cutout
    # 30 arcseconds by 30 arcseconds
    if glob.glob(sname):
        print("already saved")
    else:
        www = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi" 
        payload = {}
        payload["red"] = fname
        payload["format"] = 'fits'
        payload["size"] = '120'
        payload["ra"] = ra
        payload["dec"] = dec
        try:
            r = s.get(www, params=payload, timeout=(3.05,10), allow_redirects=True)
            open(sname, 'wb').write(r.content)
        except requests.exceptions.RequestException as e: 
            print(e)
            sys.exit(1)
        print(r.url)
        if r.status_code != requests.codes.ok:
            print("status code" + str(r.status_code))
        s.close()

    return sname # saved file


def plot(ax, fname, ptfid, vmin_val, vmax_val):
    """ plot image cutout in Python """
    hdu_list = fits.open(fname)
    image_data = hdu_list[0].data
    hdu_list.close()
    values = np.ma.masked_invalid(image_data.flatten())
    counts = np.ma.masked_invalid(image_data[::-1])
    ax.imshow(
            np.arcsinh(counts), cmap='viridis_r', 
            vmin=np.percentile(np.arcsinh(values), vmin_val),
            vmax=np.percentile(np.arcsinh(values), vmax_val))
    # N up, E left
    ax.axis('off')

    ax.scatter(60, 60, marker='x', color='white', s=40, linewidths=0.5)

    # compass markings
    # ax.annotate("", xy=(10, 9), xytext=(35,9),
    #         arrowprops=dict(facecolor='black', width=1))
    # ax.annotate("", xy=(35,34), xytext=(35,9),
    #         arrowprops=dict(facecolor='black', width=1))
    # ax.annotate("E", xy=(7,9), fontsize=11,
    #         horizontalalignment='right',
    #         verticalalignment='center')
    # ax.annotate(
    #         "S", xy=(35,38), fontsize=11, 
    #         horizontalalignment='center',
    #         verticalalignment='top')

    # line to show scale
    ax.annotate("", xy=(70,110), xytext=(110,110),
            arrowprops=dict(arrowstyle="-"))
    ax.annotate(
            "10''", xy=(90,110), fontsize=11, 
            horizontalalignment='center',
            verticalalignment='bottom')

    # name of source
    props = dict(boxstyle='round', facecolor='white')
    ax.annotate(
            ptfid, xy=(112,7), fontsize=11, 
            horizontalalignment='right',
            verticalalignment='top', bbox=props)


def afterglows():
    dat = ascii.read(
            "../afterglow_cands_confined.txt", delimiter=',',
            format='csv',
            names=['ID','RA','Dec','mag1','mag2','dmag','date','dt','time'])
    IDs = np.array(dat['ID'])
    ra = np.array(dat['RA'])
    dec = np.array(dat['Dec'])

    choose = np.array(
            ['14gqr', '13dzm', '14atg', '16erd', '16aum', '16ham'])
    vmin = np.array([90, 60, 5, 80, 80, 80])
    vmax = np.array([100, 95, 95, 100, 100, 100])
    choose = np.array(
            ['16fvd', '13edn', '15lo', '16ern', '16bdp', '16evf', '16gwh',
             '14gjx', '13dxd', '17hsa', '16esu', '16eyi'])
    vmin = np.array([90, 95, 80, 80, 80, 80, 95, 95, 95, 5, 80, 95])
    vmax = np.array([100, 100, 100, 100, 100, 100, 100, 100, 100, 95, 95, 100])
    nnames = len(choose)

    fig,axarr = plt.subplots(3,4, figsize=(6,5))
    for ii,ax in enumerate(axarr.reshape(-1)):
        if ii >= nnames:
            ax.axis('off')
        else:
            name = choose[ii]
            print(name)
            fname = ps1(
                    ra[IDs==name][0], dec[IDs==name][0]) # file saved by this function
            label_name = 'iPTF' + name
            plot(ax, fname, label_name, vmin[ii], vmax[ii])
    plt.subplots_adjust(wspace=0.1, hspace=0.01)
    plt.savefig("host_grid_afterglows.png")
    #plt.show()


def orphans():
    dat = ascii.read(
            "../orphan_cands_no_star_no_agn.txt", delimiter=',',
            format='csv',
            names=['ID','RA','Dec','dm','sig','time','max','min','date'])
    IDs = np.array(dat['ID'])
    ra = np.array(dat['RA'])
    dec = np.array(dat['Dec'])

    choose = snia()

    fig,axarr = plt.subplots(1,5, figsize=(10,3))
    for ii,ax in enumerate(axarr.reshape(-1)):
        if ii >= nnames:
            ax.axis('off')
        else:
            name = choose[ii]
            print(name)
            fname = ps1(
                    ra[IDs==name][0], dec[IDs==name][0]) # file saved by this function
            if name == '16ern':
                label_name = 'iPTF16ern (nova)'
            else:
                label_name = 'iPTF' + name
            plot(ax, fname, label_name)
    plt.subplots_adjust(wspace=0.1, hspace=0.01)
    plt.savefig("host_grid_orphans.png")


def sn():
    """ and by "sn" I really mean, rapidly evolving transients """
    dat = ascii.read(
            "../sn_cands_no_star.txt", delimiter=',',
            format='csv',
            names=['ID','RA','Dec','dm','sig','time','max','min','date'])
    IDs = np.array(dat['ID'])
    ra = np.array(dat['RA'])
    dec = np.array(dat['Dec'])

    choose = gap()
    nnames = len(choose)

    fig,axarr = plt.subplots(1,1, figsize=(3,3))
    name = choose[0]
    print(name)
    fname = ps1(
            ra[IDs==name][0], dec[IDs==name][0]) # file saved by this function
    label_name = 'iPTF' + name
    plot(axarr, fname, label_name)
    plt.subplots_adjust(wspace=0.01, hspace=0.01)
    #plt.show()
    plt.savefig("host_grid_gap.png")


if __name__=="__main__":
    afterglows()
