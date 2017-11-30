from __future__ import print_function

import numpy as np

import ddosa
import pilton
import dataanalysis as da
from dataanalysis import graphtools
import os,time,shutil
from astropy.io import fits as pyfits
from numpy import *
import re

try:
    import heaspa
except:
    pass

class ExceptionJ_SCW_NO_MINIMUM_DATA(da.AnalysisException):
    pass

class JEMX(da.DataAnalysis):
    num=2

    def get_version(self):
        return self.get_signature()+"."+self.version+".jmx%i.."%self.num

    def get_NAME(self):
        return "JMX%i"%self.num
    
    def get_name(self):
        return "jmx%i"%self.num
    
    def get_swg(self):
        return "swg_jemx%i.fits"%self.num
    
    def get_og(self):
        return "og_jmx%i.fits"%self.num

class UserCat(ddosa.DataAnalysis):
    input_cat=ddosa.GRcat

    cached=True

    version="v1_v404"

    def main(self):
        fn="jemx_user_catalog.fits"
        ddosa.remove_withtemplate(fn)

        f=pyfits.open(self.input_cat.cat[:-3])
        f[1].data['FLAG'][f[1].data['NAME']=='Ginga 2023+338']=1
        f[1].data=f[1].data[f[1].data['FLAG']==1]
        f.writeto(fn,clobber=True)

        self.cat=da.DataFile(fn)



class jemx_image(ddosa.DataAnalysis):
    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot
    input_jemx=JEMX
    input_refcat=ddosa.GRcat

    cached=True

    version="v2.1"

    def main(self):
        open("scw.list","w").write(self.input_scw.swgpath+"[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))

        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"
        ogc=ddosa.heatool(bin)
        ogc['idxSwg']="scw.list"
        ogc['instrument']=self.input_jemx.get_NAME()
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=wd # dangerous
        ogc.run()

        scwroot="scw/"+self.input_scw.scwid

        bin="jemx_science_analysis"
        ht=ddosa.heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value)
        ht['ogDOL']=self.input_jemx.get_og()
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['IC_Alias'] = "OSA"
        ht['CAT_I_refCat'] = self.input_refcat.cat
        ht['startLevel']="COR"
        ht['endLevel']="IMA"
        ht['nChanBins']=-4
        ht['jemxNum']=self.input_jemx.num

        try:
            ht.run()
        except pilton.HEAToolException as ex:
            if 'J_SCW_NO_MINIMUM_DATA' in ht.output:
                raise ExceptionJ_SCW_NO_MINIMUM_DATA(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))


        name=self.input_jemx.get_name()
        shutil.copy(ht.cwd+"/scw/"+self.input_scw.scwid+"/"+name+"_sky_ima.fits","./"+name+"_sky_ima.fits")
        shutil.copy(ht.cwd+"/scw/"+self.input_scw.scwid+"/"+name+"_srcl_res.fits","./"+name+"_srcl_res.fits")
        
        self.skyima=da.DataFile(name+"_sky_ima.fits")
        self.srclres=da.DataFile(name+"_srcl_res.fits")


class jemx_lcr(ddosa.DataAnalysis):
    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot
    input_jemx=JEMX

    tbin=100
    COR_gainModel=2

    def get_version(self):
        v=self.get_signature()+"."+self.version+".tb%.5lg"%self.tbin
        if self.COR_gainModel!=1:
            v+=".gM%i"%self.COR_gainModel
        return v

    cached=True

    version="v1"
    def main(self):
        open("scw.list","w").write(self.input_scw.swgpath+"[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))

        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"
        ogc=ddosa.heatool(bin)
        ogc['idxSwg']="scw.list"
        ogc['instrument']=self.input_jemx.get_NAME()
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=wd # dangerous
        ogc.run()

        scwroot="scw/"+self.input_scw.scwid

        bin="jemx_science_analysis"
        os.environ['COMMONSCRIPT']="1"
        ht=ddosa.heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value)
        ht['ogDOL']=self.input_jemx.get_og()
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['startLevel']="COR"
        ht['endLevel']="LCR"
        if hasattr(self,'input_usercat'):
            ht['CAT_I_usrCat']=self.input_usercat.cat.get_path()
        ht['skipLevels']=""
        ht['nChanBins']=-2
        ht['LCR_timeStep']=self.tbin
        ht['COR_gainModel']=self.COR_gainModel
        ht['jemxNum']=self.input_jemx.num
        ht.run()

        name=self.input_jemx.get_name()
        scwpath=ht.cwd+"/scw/"+self.input_scw.scwid
        shutil.copy(scwpath+"/"+name+"_src_lc.fits","./"+name+"_src_lc.fits")
        
        self.lcr=da.DataFile(name+"_src_lc.fits")
        # store the rest??


class inspect_image_results(ddosa.DataAnalysis):
    input_image=jemx_image

    def main(self):
        for r in pyfits.open(self.input_image.srclres.get_path())[1].data:
            print(r['NAME'],r['DETSIG'])
        

class inspect_image_results(ddosa.DataAnalysis):
    input_image=jemx_image

    def main(self):
        for r in pyfits.open(self.input_image.srclres.get_path())[1].data:
            print(r['NAME'],r['DETSIG'])

class JEnergyBins(ddosa.DataAnalysis):
    nchanpow=-2

    def get_version(self):
        return self.get_signature()+"."+self.version+".nchpo%.5lg"%self.nchanpow

class ISDCENV(da.DataAnalysis):
    def main(self):
        pass

class DetectISDCENV(da.DataAnalysis):
    def main(self):
        if 'OSA_VERSION' in os.environ['ISDC_ENV']:
            osa_version=os.environ['OSA_VERSION']
        else:
            osa_version=os.environ['ISDC_ENV'].split("/")[-3]
        ie=ISDCENV(use_version="OSA"+osa_version)
        if osa_version=="10.0":
            ie.noanalysis=True
        return ie

class jemx_spe(ddosa.DataAnalysis):
    input_isdcenv=DetectISDCENV

    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot
    input_jemx=JEMX
    input_usercat=UserCat
    input_jbins=JEnergyBins

    COR_gainModel=2

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.COR_gainModel!=1:
            v+=".gM%i"%self.COR_gainModel
        return v

    cached=True

    version="v1"
    def main(self):
        open("scw.list","w").write(self.input_scw.swgpath+"[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))

        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"
        ogc=ddosa.heatool(bin)
        ogc['idxSwg']="scw.list"
        ogc['instrument']=self.input_jemx.get_NAME()
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=wd # dangerous
        ogc.run()

        scwroot="scw/"+self.input_scw.scwid

        bin="jemx_science_analysis"
        os.environ['COMMONSCRIPT']="1"
        ht=ddosa.heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value)
        ht['ogDOL']=self.input_jemx.get_og()
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['startLevel']="COR"
        ht['endLevel']="SPE"
        if hasattr(self,'input_usercat'):
            ht['CAT_I_usrCat']=self.input_usercat.cat.get_full_path()
        ht['skipLevels']=""
        ht['skipSPEfirstScw']="n"
        ht['nChanBins']=self.input_jbins.nchanpow
        #ht['LCR_timeStep']=self.tbin
        ht['COR_gainModel']=self.COR_gainModel
        ht['jemxNum']=self.input_jemx.num
        ht.run()

        name=self.input_jemx.get_name()
        scwpath=ht.cwd+"/scw/"+self.input_scw.scwid
        shutil.copy(scwpath+"/"+name+"_srcl_spe.fits","./"+name+"_srcl_spe.fits")
        shutil.copy(scwpath+"/"+name+"_srcl_arf.fits","./"+name+"_srcl_arf.fits")
        
        self.spe=da.DataFile(name+"_srcl_spe.fits")
        self.arf=da.DataFile(name+"_srcl_arf.fits")
        # store the rest??

class JRMF(ddosa.DataAnalysis):
    input_jbins=JEnergyBins
    cached=True

    def main(self):
        code='STD_%.3i'%(2**(-self.input_jbins.nchanpow))
        fn='jemx_rmf_%s.fits'%code
        os.system('j_rebin_rmf binlist=%s outfile=jemx_rmf_%s.fits'%(code,fn))
        self.rmf=da.DataFile(fn)

class ProcessJSpectra(ddosa.DataAnalysis):
    input_spectrum=jemx_spe
    input_rmf=JRMF
    input_jmx=JEMX

    def main(self):
        arfs_data=pyfits.open(self.input_spectrum.arf.get_path())[1].data
        for source_data in pyfits.open(self.input_spectrum.spe.get_path())[1].data:
            name=source_data['NAME']
            print(name,sum((source_data['RATE']/source_data['STAT_ERR'])**2)**0.5)
            fn="jemx%i_spectrum_%s.fits"%(self.input_jmx.num,name.replace(" ","_"))
            heaspa.PHA(source_data['RATE'].astype(float64),source_data['STAT_ERR'].astype(float64),exposure=source_data['EXPOSURE'],datatype="RATE").write(fn)

            arf_pos=re.search("jmx._srcl_arf.fits\{(.*?)\}",source_data['ANCRFILE']).group(1)
            arf_data=arfs_data[int(arf_pos)-1]
            setattr(self,fn,da.DataFile(fn))

            arffn="jemx%i_arf_%s.fits"%(self.input_jmx.num,name.replace(" ","_"))
            heaspa.ARF(arf_data['ENERG_LO'].astype(float64),arf_data['ENERG_HI'].astype(float64),arf_data['SPECRESP'].astype(float64)).write(arffn)

            f=pyfits.open(fn)
            f[1].header['ANCRFILE']=arffn
            f[1].header['RESPFILE']=self.input_rmf.rmf.get_path()
            f.writeto(fn,clobber=True)
            setattr(self,arffn,da.DataFile(arffn))

#j_rebin_rmf binlist=STD_016



class JMXScWImageList(ddosa.DataAnalysis):
    input_scwlist=None
    copy_cached_input=False
    input_imagingsummary=graphtools.Factorize(use_root='jemx_image',use_leaves=["ScWData",])


    allow_alias=True
    run_for_hashe=True

    version="allthem"

    maxima=None
    firstima=0

    def get_version(self):
        if self.maxima is None and self.firstima==0:
            return self.get_signature()+"."+self.version
        return self.get_signature()+"."+repr(self.firstima)+"."+repr(self.maxima)

    def main(self):
        print("jemx_image constructed as",ddosa.ii_skyimage())
        #ddosa.ShadowUBCImage(assume=scw)
        self.images=[(scw,jemx_image(assume=scw),) for scw in self.input_scwlist.scwlistdata]

        if len(self.images)==0:
            raise ddosa.EmptyScWList()

        print("images will be:",self.images[0],"..",self.images[-1])
        if self.maxima is not None:
            self.images=self.images[self.firstima:self.firstima+self.maxima]

# this is so old
def angsep(ra1,dec1,ra2,dec2):
    ra1*= np.pi / 180.
    dec1*= np.pi / 180.
    ra2*= np.pi / 180.
    dec2*= np.pi / 180.
    SEP=np.arccos(np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2) + np.sin(dec1) * np.sin(dec2)) #returns values between 0 and pi radians
    SEP*= 180. / np.pi
    return SEP

class mosaic_jemx(ddosa.DataAnalysis):
    input_imagelist = JMXScWImageList

    # write_caches=[da.TransientCache,ddosa.MemCacheIntegralFallback,ddosa.MemCacheIntegralIRODS]
    # read_caches=[da.TransientCache,ddosa.MemCacheIntegralFallback,ddosa.MemCacheIntegralFallbackOldPath,ddosa.MemCacheIntegralIRODS]

    copy_cached_input = False

    maxsep = 30

    pixdivide = 2

    cached = True

    test_files = False

    version = "v1"

    def get_version(self):
        v = self.get_signature() + "." + self.version

        if self.pixdivide != 2:
            v += ".pd%i" % self.pixdivide
        return v

    def main(self):
        print("will mosaic")

        self.split_in_lists()
        self.mosaic_lists()
        del self.lists

    def choose_list(self, (ra, dec)):
        for (list_ra, list_dec), thelist in self.lists:
            if angsep(list_ra, list_dec, ra, dec) < self.maxsep:
                return thelist
        self.lists.append([(ra, dec), []])
        print("new list at", ra, dec)
        return self.lists[-1][1]

    def split_in_lists(self):
        self.lists = []

        for l in self.input_imagelist.images:
            scw, image= l[:2]
            print("image:", image)

            print(image._da_locally_complete)

            if hasattr(image, 'empty_results') and image.empty_results or not hasattr(image, 'skyima'):
                print("skipping", image)
                continue

            fn = image.skyima.get_path()
            h = pyfits.open(fn)[2].header
            ra, dec = h['CRVAL1'], h['CRVAL2']

            thelist = self.choose_list((ra, dec))
            # print "offered list",thelist

            thelist.append((fn, image))
            print("adding image to list:", len(thelist))

        for thelist in self.lists:
            print(thelist[0], len(thelist[1]))

    def mosaic_lists(self):
        self.mosaics = []
        self.stacked = []
        for i, ((ra, dec), thelist) in enumerate(self.lists):
            print("list", i, ra, dec, thelist)

            sens_stats_fn = "sens_stats_%.5lg_%.5lg.txt" % (ra, dec)

            listfile = "image_list_%.5lg_%.5lg.txt" % (ra, dec)
            f = open(listfile, "w")

            sens_stats = []
            sens_file = open(sens_stats_fn, "w")
            sens_keys = None

            for image, imageobj, in thelist:
                print(image, )
                if hasattr(image, 'empty_results'):
                    print("skipping", image)
                    continue

                if False: #has_imagesensitivity:
                    print(sens.__dict__)

                    statistic = sens.statistics[0]

                    fe = pyfits.open(imageobj.skyima.get_path())[2]
                    statistic['tstart'] = fe.header['TSTART']
                    statistic['tstop'] = fe.header['TSTOP']

                    statistic['sensi_min_exposure_corr'] = statistic['sensi_min'] * (statistic[
                                                                                         'exposure_at_sensi_min'] / 2000.) ** 0.5
                    if statistic['sensi_min_exposure_corr'] == 0:
                        statistic['sensi_min_exposure_corr'] = 100

                    statistic['flag'] = 1
                    if statistic['sensi_min_exposure_corr'] > 400.:  # hardcoded!
                        statistic['flag'] = 0

                    sens_stats.append([imageobj.assumptions, statistic])

                    if sens_keys is None:
                        sens_keys = statistic.keys()
                        sens_file.write("fn " + (" ".join(sens_keys)) + "\n")

                    sens_file.write(image + " ")
                    sens_file.write(" ".join("%.5lg" % statistic[k] for k in sens_keys) + "\n")

                    if statistic['flag'] <= 0:
                        continue

                f.write(image + "\n")

            f.close()

            stacked_file = None
            stacked_aligned = np.zeros((1500, 1500))
            stacked_aligned_var = np.zeros((1500, 1500))
            stacked_shad_file = None
            stacked_effi_file = None
            for imagefilename, imageobj in thelist:
                print("stacking", imagefilename)

                image_file = pyfits.open(imagefilename)
                if stacked_file is None:
                    stacked_file = image_file
                else:
                    for i, e in enumerate(image_file[2:]):
                        stacked_file[i + 2].data += e.data

                # aligned
                flux = image_file[2].data
                var = image_file[3].data
                sig = image_file[4].data
                var[np.isnan(var)] = np.inf
                var[var == 0] = np.inf
                peak = np.unravel_index(np.nanargmax(sig), sig.shape)

                print("peak at", peak, "of", sig[peak])

                i, j = np.meshgrid(np.arange(sig.shape[0]), np.arange(sig.shape[1]))
                c = i + 500 - peak[1], j + 500 - peak[0]
                # c=i-(peak[0])+stacked_aligned.shape[0]/2,j-(peak[1])+stacked_aligned.shape[1]/2

                print(i.max(), j.max())
                stacked_aligned[c] += flux / var
                stacked_aligned_var[c] += 1 / var

            self.sens_stats_file = da.DataFile(sens_stats_fn)

            stacked = "stacked_aligned_%.5lg_%.5lg.fits" % (ra, dec)
            f = stacked_aligned / stacked_aligned_var
            v = 1 / stacked_aligned_var
            pyfits.PrimaryHDU(f / v ** 0.5).writeto(stacked, clobber=True)
            setattr(self, stacked, da.DataFile(stacked))

            stacked = "stacked_%.5lg_%.5lg.fits" % (ra, dec)
            stacked_file.writeto(stacked, clobber=True)
            self.stacked.append([(ra, dec), da.DataFile(stacked)])
            setattr(self, stacked, da.DataFile(stacked))

            mosaic = "mosaic_%.5lg_%.5lg.fits" % (ra, dec)
            regionfile = mosaic.replace(".fits", ".reg")

            ht = ddosa.heatool(
                os.environ['COMMON_INTEGRAL_SOFTDIR'] + "/imaging/varmosaic/varmosaic_exposure/varmosaic")
            ht['pixdivide'] = self.pixdivide
            ht['filelist'] = listfile
            ht['outimage'] = mosaic
            ht['outregion'] = regionfile
            ht.run()

            self.mosaics.append([(ra, dec), da.DataFile(mosaic)])
            setattr(self, mosaic, da.DataFile(mosaic))
            setattr(self, regionfile, da.DataFile(regionfile))

            self.skyima = da.DataFile(mosaic)  # store last