from __future__ import print_function

import ddosa
import dataanalysis as da
import os,time,shutil
from astropy.io import fits as pyfits
from numpy import *
import re

try:
    import heaspa
except:
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

    version="v2"
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
        ht.run()

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
