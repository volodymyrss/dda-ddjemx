

import numpy as np
from copy import deepcopy
import os
import subprocess
import glob

import ddosa
import pilton
import dataanalysis as da
from dataanalysis import graphtools
import os,time,shutil
from astropy.io import fits as fits
from numpy import *
import re

corrupt_scws = ['021100220010', '039100070010', '054100130010', '078800590010', '078800460010', '078800590010', '114600470010', '109400250010', '127500100010']


class OSAEnv(ddosa.DataAnalysis):
    version="10.2"

try:
    import heaspa
except:
    pass


class ExceptionJ_SCW_NO_MINIMUM_DATA(da.AnalysisException):
    pass


class ExceptionJ_SCW_LackingData(da.AnalysisException):
    pass


class ExceptionNoImageProduced(da.AnalysisException):
    pass


class ExceptionNoSpectraProduced(da.AnalysisException):
    pass


class ExceptionNoSCPandGainHistory(da.AnalysisException):
    pass


class ExceptionNoLCProduced(da.AnalysisException):
    pass


class ExceptionFailingScWUnknownReasonOGC(da.AnalysisException):
    pass

class ExceptionCountNotSolveGainVariation(da.AnalysisException):
    pass

class BadAfterOnePass(da.AnalysisException):
    pass

class OSACrash(da.AnalysisException):
    pass

class NoUsefulData(da.AnalysisException):
    pass

class SegFault(Exception):
    pass


class SegFaultInMosaic(da.AnalysisException):
    pass


class CorruptScwInSpePick(da.AnalysisException):
    pass

class NoGAIN(Exception):
    pass

class UnknownJemxError(Exception):
    pass

class JEMX(da.DataAnalysis):
    num=1

    def get_version(self):
        return self.get_signature()+"."+self.version+".jmx%i.."%self.num

    def get_NAME(self):
        return "JMX%i"%self.num
    
    def get_name(self):
        return "jmx%i"%self.num

    def get_longish_name(self):
        return "jemx%i" % self.num

    def get_swg(self):
        return "swg_jmx%i.fits"%self.num
    
    def get_og(self):
        return "og_jmx%i.fits"%self.num

class UserCat(ddosa.DataAnalysis):
    input_cat=ddosa.GRcat

    cached=True

    version="v2"

    def main(self):
        fn="jemx_user_catalog.fits"
        ddosa.remove_withtemplate(fn)

        try:
            f = fits.open(re.sub("\[.\]", "", self.input_cat.cat))
        except Exception as e:
            print("failed with simple path:", e, self.input_cat.cat)
            f = fits.open(self.input_cat.cat.get_full_path())

        #f[1].data['FLAG'][f[1].data['NAME']=='Ginga 2023+338']=1
        #f[1].data=f[1].data[f[1].data['FLAG']==1]

        f[1].data['FLAG'] = 1
        f.writeto(fn,overwrite=True)

        self.cat=da.DataFile(fn)

class LCTimeBin(ddosa.LCTimeBin):
    pass

class JEnergyBins(ddosa.DataAnalysis):
    nchanpow=-4
    input_ic=ddosa.ICRoot

    bins=None

    def get_version(self):
        v=self.get_signature()+"."+self.version

        if self.bins is not None:
            v+="".join([".%.5lg_%.5lg"%(e1,e2) for e1,e2 in self.bins])
        else:
            v += ".nchpo%.5lg" % self.nchanpow

        return v

    def main(self):
        if self.bins is not None:
            rsp_bins=fits.open(self.input_ic.icroot+"/ic/jmx1/rsp/jmx1_rmf_grp_0046.fits")[3].data

            self.bin_interpretation=[]

            for e1,e2 in self.bins:
                ch1 = rsp_bins['CHANNEL'][rsp_bins['E_MIN'] > e1][0]
                ch2 = rsp_bins['CHANNEL'][rsp_bins['E_MAX'] < e2][-1]
                e1_true,e2_true=rsp_bins['E_MIN'][rsp_bins['CHANNEL'] == ch1][0], rsp_bins['E_MAX'][rsp_bins['CHANNEL'] == ch2][0]
                print("for",e1,e2,"channels",ch1,ch2,"true energies",e1_true,e2_true)
                self.bin_interpretation.append(
                    dict(
                        emin=e1,
                        emax=e2,
                        emin_true=e1_true,
                        emax_true=e2_true,
                        chmin=ch1,
                        chmax=ch2,
                    )
                )

class JEnergyBinsSpectra(JEnergyBins):
    pass

class JEnergyBinsLC(JEnergyBins):
    nchanpow=-1


def env_for_scw(self):
    env = None
    if hasattr(self, 'input_scw'):
        if self.input_scw.scwid.endswith('.000'):
            env = deepcopy(os.environ)
            env['REP_BASE_PROD'] = os.environ.get("REP_BASE_PROD_NRT")
            print("\033[31msetting RBP for NRT in ", self, ":", env['REP_BASE_PROD'], "\033[0m")
        else:
            print("\033[31mnot nrt, so NOT setting RBP for NRT in ", self, ":", os.getenv('REP_BASE_PROD'), "\033[0m")
    else:
        print("\033[31mno scw, so NOT setting RBP for NRT in ", self, ":", os.getenv('REP_BASE_PROD'), "\033[0m")

    return env


class jemx_image(ddosa.DataAnalysis):
    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot
    input_jemx=JEMX
    input_refcat=ddosa.GRcat
    input_jbins=JEnergyBins
    input_osaenv = OSAEnv

    cached=True

    version="v2.2.4"

    def main(self):
        open("scw.list","w").write(self.input_scw.swgpath+"[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))

        env = env_for_scw(self)        

        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"

        ogc = ddosa.heatool(bin, env=env_for_scw(self)
)        
        ogc['idxSwg'] = "scw.list"
        ogc['instrument'] = self.input_jemx.get_NAME()
        ogc['ogid'] = "scw_" + self.input_scw.scwid
        ogc['baseDir']=wd # dangerous

        try:
            ogc.run()
        except pilton.HEAToolException as e:
            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()
            raise
            
        print("\033[31m OG created! \033[0m")

        scwroot="scw/"+self.input_scw.scwid

        bin = "jemx_science_analysis"
        ht = ddosa.heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value, env=env_for_scw(self)
)

        ht['ogDOL']=self.input_jemx.get_og()
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['IC_Alias'] = "OSA"
        ht['CAT_I_refCat'] = self.input_refcat.cat
        ht['startLevel']="COR"
        ht['endLevel']="IMA"

        if self.input_scw.scwid.endswith('.000'):
            ht['COR_gainModel'] = 2

        if self.input_jbins.bins is None:
            ht['nChanBins']=self.input_jbins.nchanpow
        else:
            ht['nChanBins'] = len(self.input_jbins.bin_interpretation)
            ht['chanLow'] = " ".join(["%i"%bin['chmin'] for bin in self.input_jbins.bin_interpretation])
            ht['chanHigh'] = " ".join(["%i"%bin['chmax'] for bin in self.input_jbins.bin_interpretation])

        ht['jemxNum']=self.input_jemx.num

        if hasattr(self,'input_usergti'):
            path=self.input_usergti.gti.get_path()
            if os.path.abspath(path)==os.path.normpath(path):
                print("full path",path)
            else:
                print("not a full path",path)
                path="../../"+path
            ht['GTI_gtiUser']=path
            ht['GTI_TimeFormat']='UTC'

        try:
            ht.run()
        except pilton.HEAToolException as ex:
            if 'J_SCW_NO_MINIMUM_DATA' in ht.output:
                raise ExceptionJ_SCW_NO_MINIMUM_DATA(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))
            
            if 'After one-pass-loop, status = -1' in ht.output:
                raise BadAfterOnePass(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))
            
            if 'Could not solve gain history time-variation problem: -321122' in ht.output:
                raise ExceptionCountNotSolveGainVariation(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

            if "Couldn't find SCP and gain history elapsed time:" in ht.output:
                raise ExceptionNoSCPandGainHistory(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

            
        
        if 'J_COR_BAD_GAINHISTDOL' in ht.output:
            raise NoGAIN()
        
        if 'Element JMX2-SPEC-* was not found' in ht.output:
            raise NoUsefulData()
        
        if 'Element JMX1-SPEC-* was not found' in ht.output:
            raise NoUsefulData()
        

        if 'null or NaN values for this particular revolution' in ht.output:
            raise NoUsefulData()

        if 'Error_2: Task jemx_science_analysis terminating' in ht.output:
            raise UnknownJemxError()

        if 'segmentation violation' in ht.output:
            raise SegFault()
        
        if 'Segmentation fault' in ht.output:
            raise SegFault()
        
        #if 'No Offline Gain Calibration File' in ht.output:
        #    raise NoGAIN()

        name=self.input_jemx.get_name()

        skyima_fn = ht.cwd+"/scw/"+self.input_scw.scwid+"/"+name+"_sky_ima.fits"

        skyres_fn = ht.cwd+"/scw/"+self.input_scw.scwid+"/"+name+"_srcl_res.fits"
        
        if os.path.exists(skyima_fn) and os.path.exists(skyres_fn):
            shutil.copy(skyima_fn, "./"+name+"_sky_ima.fits")
            shutil.copy(skyres_fn, "./"+name+"_srcl_res.fits")
            
            self.skyima=da.DataFile(name+"_sky_ima.fits")
            self.srclres=da.DataFile(name+"_srcl_res.fits") 
        #else:
        #    raise ExceptionNoImageProduced(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

class jemx_slew(ddosa.DataAnalysis):
    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot
    input_jemx=JEMX

    COR_gainModel=-1

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.COR_gainModel!=1:
            v+=".gM%i"%self.COR_gainModel
        return v

    cached=True

    version="v1.3.5"

    def main(self): 
        t1 = time.time()
        open("scw.list","w").write(self.input_scw.swgpath+"[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))

        env = None
        if self.input_scw.scwid.endswith('.000'):
            env = deepcopy(os.environ)
            env['REP_BASE_PROD'] = os.environ.get("REP_BASE_PROD_NRT")
            print("\033[31msetting RBP for NRT:", env['REP_BASE_PROD'], "\033[0m")

        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"
        ogc=ddosa.heatool(bin, env=env_for_scw(self))
        ogc['idxSwg']="scw.list"
        ogc['instrument']=self.input_jemx.get_NAME()
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=wd # dangerous

        try:
            ogc.run()
        except pilton.HEAToolException as ex:
            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()
            raise
        
        bin="jemx_scw_analysis"
        os.environ['COMMONSCRIPT']="1"
        
        ht=ddosa.heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value, env=env_for_scw(self))
        ht['swgDOL']='scw/' + str(self.input_scw.input_scwid.handle) + '/' + self.input_jemx.get_swg()
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['startLevel']="COR"
        ht['endLevel']="DEAD"
        ht['skipSPEfirstScw']="no"
        ht['ScwType']="SLEW"
        ht['general_jemxNum']=self.input_jemx.num
        
        try:
            ht.run()
        except pilton.HEAToolException as ex:
            if 'J_SCW_NO_MINIMUM_DATA' in ht.output:
                raise ExceptionJ_SCW_NO_MINIMUM_DATA(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

            if 'After one-pass-loop, status = -1' in ht.output:
                raise BadAfterOnePass(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()


        if 'segmentation violation' in ht.output:
            raise SegFault()
        
        name=self.input_jemx.get_name()
        scwpath=ht.cwd+"/scw/"+self.input_scw.scwid

        cor = scwpath+"/"+name+"_full_cor.fits"

        if os.path.exists(cor):
            shutil.copy(cor, name+"_full_cor.fits")
            self.cor=da.DataFile(name + "_full_cor.fits")

        gti = scwpath+"/"+name+"_gti.fits"

        if os.path.exists(gti):
            shutil.copy(gti, name+"_gti.fits")
            self.gti=da.DataFile(name + "_gti.fits")

        dead = scwpath+"/"+name+"_dead_time.fits"

        if os.path.exists(dead):
            shutil.copy(dead, name+"_dead_time.fits")
            self.dead=da.DataFile(name + "_dead_time.fits")

        t2 = time.time()
        print('JEM-X slew took {}'.format(t2 - t1))

class jemx_base(ddosa.DataAnalysis):
    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot
    input_jemx=JEMX
    input_refcat=ddosa.GRcat
    input_jbins=JEnergyBinsLC
    input_timebin=LCTimeBin
    input_osaenv = OSAEnv

    COR_gainModel=-1

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.COR_gainModel!=1:
            v+=".gM%i"%self.COR_gainModel
        return v

    cached=True

    version="v1.0"

    def main(self):
        t1 = time.time()
        open("scw.list","w").write(self.input_scw.swgpath+"[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))

        env = None
        if self.input_scw.scwid.endswith('.000'):
            env = deepcopy(os.environ)
            env['REP_BASE_PROD'] = os.environ.get("REP_BASE_PROD_NRT")
            print("\033[31msetting RBP for NRT:", env['REP_BASE_PROD'], "\033[0m")


        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"
        ogc=ddosa.heatool(bin, env=env_for_scw(self))
        ogc['idxSwg']="scw.list"
        ogc['instrument']=self.input_jemx.get_NAME()
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=wd # dangerous
        try:
            ogc.run()
        except pilton.HEAToolException as ex:
            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()
            raise

        scwroot="scw/"+self.input_scw.scwid

        bin="jemx_science_analysis"
        os.environ['COMMONSCRIPT']="1"
        ht=ddosa.heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value, env=env_for_scw(self)
)
        ht['ogDOL']=self.input_jemx.get_og()
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['startLevel']="COR"
        ht['endLevel']="DEAD"
        ht['CAT_I_refCat'] = self.input_refcat.cat

        if hasattr(self,'input_usercat'):
            ht['CAT_I_usrCat']=self.input_usercat.cat.get_full_path()

        ht['skipLevels']="SPE"

        if self.input_jbins.bins is None:
            ht['nChanBins'] = self.input_jbins.nchanpow
        else:
            ht['nChanBins'] = len(self.input_jbins.bin_interpretation)
            ht['chanLow'] = " ".join(["%i" % bin['chmin'] for bin in self.input_jbins.bin_interpretation])
            ht['chanHigh'] = " ".join(["%i" % bin['chmax'] for bin in self.input_jbins.bin_interpretation])

        ht['LCR_timeStep']=self.input_timebin.time_bin_seconds
        ht['COR_gainModel']=self.COR_gainModel
        ht['jemxNum']=self.input_jemx.num

        try:
            ht.run()
        except pilton.HEAToolException as ex:
            if 'J_SCW_NO_MINIMUM_DATA' in ht.output:
                raise ExceptionJ_SCW_NO_MINIMUM_DATA(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

            if 'After one-pass-loop, status = -1' in ht.output:
                raise BadAfterOnePass(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()


        if 'segmentation violation' in ht.output:
            raise SegFault()

        name=self.input_jemx.get_name()
        scwpath=ht.cwd+"/scw/"+self.input_scw.scwid

        cor = scwpath+"/"+name+"_full_cor.fits"

        if os.path.exists(cor):
            shutil.copy(cor, name+"_full_cor.fits")
            self.cor=da.DataFile(name + "_full_cor.fits")

        gti = scwpath+"/"+name+"_gti.fits"

        if os.path.exists(gti):
            shutil.copy(gti, name+"_gti.fits")
            self.gti=da.DataFile(name + "_gti.fits")

        dead = scwpath+"/"+name+"_dead_time.fits"

        if os.path.exists(dead):
            shutil.copy(dead, name+"_dead_time.fits")
            self.dead=da.DataFile(name + "_dead_time.fits")

         #else:
        #    raise ExceptionNoLCProduced()

        t2 = time.time()
        print('JEM-X pointing took {}'.format(t2 - t1))




class jemx_lcr(ddosa.DataAnalysis):
    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot
    input_jemx=JEMX
    input_refcat=ddosa.GRcat
    input_jbins=JEnergyBinsLC
    input_timebin=LCTimeBin
    input_osaenv = OSAEnv

    COR_gainModel=-1

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.COR_gainModel!=1:
            v+=".gM%i"%self.COR_gainModel
        return v

    cached=True

    version="v1.3.5"

    def main(self):
        t1 = time.time()
        open("scw.list","w").write(self.input_scw.swgpath+"[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))

        env = None
        if self.input_scw.scwid.endswith('.000'):
            env = deepcopy(os.environ)
            env['REP_BASE_PROD'] = os.environ.get("REP_BASE_PROD_NRT")
            print("\033[31msetting RBP for NRT:", env['REP_BASE_PROD'], "\033[0m")


        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"
        ogc=ddosa.heatool(bin, env=env_for_scw(self))
        ogc['idxSwg']="scw.list"
        ogc['instrument']=self.input_jemx.get_NAME()
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=wd # dangerous
        try:
            ogc.run()
        except pilton.HEAToolException as ex:
            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()
            raise

        scwroot="scw/"+self.input_scw.scwid

        bin="jemx_science_analysis"
        os.environ['COMMONSCRIPT']="1"
        ht=ddosa.heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value, env=env_for_scw(self)
)
        ht['ogDOL']=self.input_jemx.get_og()
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['startLevel']="COR"
        ht['endLevel']="LCR"
        ht['CAT_I_refCat'] = self.input_refcat.cat

        if hasattr(self,'input_usercat'):
            ht['CAT_I_usrCat']=self.input_usercat.cat.get_full_path()

        ht['skipLevels']="SPE"

        if self.input_jbins.bins is None:
            ht['nChanBins'] = self.input_jbins.nchanpow
        else:
            ht['nChanBins'] = len(self.input_jbins.bin_interpretation)
            ht['chanLow'] = " ".join(["%i" % bin['chmin'] for bin in self.input_jbins.bin_interpretation])
            ht['chanHigh'] = " ".join(["%i" % bin['chmax'] for bin in self.input_jbins.bin_interpretation])

        ht['LCR_timeStep']=self.input_timebin.time_bin_seconds
        ht['COR_gainModel']=self.COR_gainModel
        ht['jemxNum']=self.input_jemx.num

        try:
            ht.run()
        except pilton.HEAToolException as ex:
            if 'J_SCW_NO_MINIMUM_DATA' in ht.output:
                raise ExceptionJ_SCW_NO_MINIMUM_DATA(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

            if 'After one-pass-loop, status = -1' in ht.output:
                raise BadAfterOnePass(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))

            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()


        if 'segmentation violation' in ht.output:
            raise SegFault()

        name=self.input_jemx.get_name()
        scwpath=ht.cwd+"/scw/"+self.input_scw.scwid

        lc = scwpath+"/"+name+"_src_lc.fits"

        if os.path.exists(lc):
            shutil.copy(lc, "./"+name+"_src_lc.fits")
            self.lcr=da.DataFile(name+"_src_lc.fits")

        lc = scwpath+"/"+name+"_src_iros_lc.fits"
        
        if os.path.exists(lc):
            shutil.copy(lc, "./"+name+"_src_lc_iros.fits")
            self.lcr=da.DataFile(name+"_src_lc_iros.fits")

        srcl_res = scwpath+"/"+name+"_srcl_res.fits"

        if os.path.exists(srcl_res):
            shutil.copy(srcl_res, name+"_srcl_res.fits")
            self.res=da.DataFile(name+"_srcl_res.fits")

                 #else:
        #    raise ExceptionNoLCProduced()

        t2 = time.time()
        print('JEM-X pointing took {}'.format(t2 - t1))


class inspect_image_results(ddosa.DataAnalysis):
    input_image=jemx_image

    def main(self):
        for r in fits.open(self.input_image.srclres.get_path())[1].data:
            print(r['NAME'],r['DETSIG'])
        

class inspect_image_results(ddosa.DataAnalysis):
    input_image=jemx_image

    def main(self):
        for r in fits.open(self.input_image.srclres.get_path())[1].data:
            print(r['NAME'],r['DETSIG'])



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
#    input_isdcenv=DetectISDCENV

    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot
    input_jemx=JEMX
  #  input_usercat=UserCat
    input_refcat=ddosa.GRcat
    input_jbins=JEnergyBinsSpectra
    input_osaenv = OSAEnv

  #  input_image=jemx_image  attempt to separate imaging

    COR_gainModel=-1

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.COR_gainModel!=1:
            v+=".gM%i"%self.COR_gainModel
        return v

    cached=True

    version="v1.2.4"

    def main(self):
        open("scw.list","w").write(self.input_scw.swgpath+"[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))

        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"
        ogc=ddosa.heatool(bin, env=env_for_scw(self)
)
        ogc['idxSwg']="scw.list"
        ogc['instrument']=self.input_jemx.get_NAME()
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=wd # dangerous

        try:
            ogc.run()
        except pilton.HEAToolException as ex:
            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()
            raise

        scwroot="scw/"+self.input_scw.scwid
        
        #ht=ddosa.heatool("dal_attach")
        #ht['Parent']=wd+"/obs/"+ogc['ogid'].value+"/"+scwroot+"/"+self.input_jemx.get_swg()
        #ht['Child1']=self.input_image.skyima.get_path()
        #ht['Child2']=self.input_image.srclres.get_path()
        #ht.run()



        bin="jemx_science_analysis"
        os.environ['COMMONSCRIPT']="1"
        os.environ['REP_BASE_PROD']=ddosa.detect_rbp(self.input_scw.scwver)
        ht=ddosa.heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value, env=env_for_scw(self)
)
        ht['ogDOL']=self.input_jemx.get_og()
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['startLevel']="COR"
        ht['endLevel']="SPE"
        ht['CAT_I_refCat'] = self.input_refcat.cat

        if hasattr(self,'input_usercat'):
            ht['CAT_I_usrCat']=self.input_usercat.cat.get_full_path()
        ht['skipLevels']="LCR"
        #ht['skipLevels']="BIN_I,IMA,BIN_T,LCR" # attempt to separate imaging
        #ht['skipLevels']="CAT_I,BIN_I,IMA,BIN_T,LCR"
        ht['skipSPEfirstScw']="y"
        #ht['skipSPEfirstScw']="n"

        if self.input_jbins.bins is None:
            ht['nChanBins']=self.input_jbins.nchanpow
        else:
            ht['nChanBins'] = len(self.input_jbins.bin_interpretation)
            ht['chanLow'] = " ".join(["%i"%bin['chmin'] for bin in self.input_jbins.bin_interpretation])
            ht['chanHigh'] = " ".join(["%i"%bin['chmax'] for bin in self.input_jbins.bin_interpretation])

        #ht['LCR_timeStep']=self.tbin
        ht['COR_gainModel']=self.COR_gainModel
        ht['jemxNum']=self.input_jemx.num
        
        try:
            ht.run()
        except UnicodeEncodeError as e: 
            raise OSACrash(str(e))
        except pilton.HEAToolException as ex:
            for corrupt_scw in corrupt_scws:
                if corrupt_scw in self.input_scw.swgpath:
                    raise ExceptionFailingScWUnknownReasonOGC()

            if 'J_SCW_NO_MINIMUM_DATA' in ht.output:
                raise ExceptionJ_SCW_NO_MINIMUM_DATA(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))
            elif '-SPTI-* was not found' in ht.output:
                raise ExceptionJ_SCW_LackingData(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))
            elif 'Floating point exception' in ht.output:
                raise ExceptionNoSpectraProduced(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))
            elif "Couldn't find SCP and gain history elapsed time:" in ht.output:
                raise ExceptionNoSCPandGainHistory(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))
            elif 'After one-pass-loop, status = -1' in ht.output:
                raise BadAfterOnePass(dict(scw=self.input_scw.scwid,jemx=self.input_jemx.get_name()))
            else:
                raise

        if 'segmentation violation' in ht.output:
            raise SegFault()

        name=self.input_jemx.get_name()
        scwpath=ht.cwd+"/scw/"+self.input_scw.scwid

        srcl_spe = scwpath+"/"+name+"_srcl_spe.fits"
        srcl_arf = scwpath+"/"+name+"_srcl_arf.fits"
        srcl_res = scwpath+"/"+name+"_srcl_res.fits"
        
        if os.path.exists(srcl_spe) and os.path.exists(srcl_arf):
            shutil.copy(srcl_spe, name+"_srcl_spe.fits")
            shutil.copy(srcl_arf, name+"_srcl_arf.fits")
            shutil.copy(srcl_res, name+"_srcl_res.fits")
        
            self.spe=da.DataFile(name+"_srcl_spe.fits")
            self.arf=da.DataFile(name+"_srcl_arf.fits")
            self.res=da.DataFile(name+"_srcl_res.fits")
#        else:
            #raise ExceptionNoSpectraProduced()

class JRMF(ddosa.DataAnalysis):
    input_jbins=JEnergyBins
    input_jmx=JEMX

    cached=True

    version="v1.1"

    def get_version(self):
        v = super().get_version()
        fversion = subprocess.check_output(["fversion"]).decode().strip()

        if '6.25' in fversion:
            return v + ".heasoft_safe"

        if '6.24' in fversion:
            return v + ".heasoft_safe"
        
        if '6.26' in fversion:
            return v + ".heasoft_safe"

        return v + ".heasoft_UNSAFE"



    def main(self):
        code='STD_%.3i'%(2**(-self.input_jbins.nchanpow))
        fn='jemx_rmf_%s.fits'%code

        env = deepcopy(os.environ)
        env['REP_BASE_PROD'] = env['CURRENT_IC']

        print("env",env)

        ddosa.remove_withtemplate(fn+"(JMX-RMF.tpl)")

        cmd = [
                'j_rebin_rmf',
                'binlist=%s'%code,
                'outfile=%s'%fn,
                'jemx_num=%i'%self.input_jmx.num
            ]

        print("running"," ".join(cmd))

        subprocess.check_call(
                        cmd,
                        env=env_for_scw(self)
,
                    )

        if os.path.exists(fn):
            self.rmf=da.DataFile(fn)
        else:
            raise RuntimeError("no %s found, have this: %s"%(fn, glob.glob('*')))

class ProcessJSpectra(ddosa.DataAnalysis):
    input_spectrum=jemx_spe
    input_rmf=JRMF
    input_jmx=JEMX

    def main(self):
        arfs_data=fits.open(self.input_spectrum.arf.get_path())[1].data
        for source_data in fits.open(self.input_spectrum.spe.get_path())[1].data:
            name=source_data['NAME']
            print(name,sum((source_data['RATE']/source_data['STAT_ERR'])**2)**0.5)
            fn="jemx%i_spectrum_%s.fits"%(self.input_jmx.num,name.replace(" ","_"))
            heaspa.PHA(source_data['RATE'].astype(float64),source_data['STAT_ERR'].astype(float64),exposure=source_data['EXPOSURE'],datatype="RATE").write(fn)

            arf_pos=re.search("jmx._srcl_arf.fits\{(.*?)\}",source_data['ANCRFILE']).group(1)
            arf_data=arfs_data[int(arf_pos)-1]
            setattr(self,fn,da.DataFile(fn))

            arffn="jemx%i_arf_%s.fits"%(self.input_jmx.num,name.replace(" ","_"))
            heaspa.ARF(arf_data['ENERG_LO'].astype(float64),arf_data['ENERG_HI'].astype(float64),arf_data['SPECRESP'].astype(float64)).write(arffn)

            f=fits.open(fn)
            f[1].header['ANCRFILE']=arffn
            f[1].header['RESPFILE']=self.input_rmf.rmf.get_path()
            f.writeto(fn,overwrite=True)
            setattr(self,arffn,da.DataFile(arffn))

#j_rebin_rmf binlist=STD_016

class jemx_image_by_scw(graphtools.Factorize):
    root='jemx_image'
    leaves=["ScWData",]

class JMXScWImageList(ddosa.DataAnalysis):
    input_scwlist=None
    copy_cached_input=False
    input_imagingsummary=jemx_image_by_scw


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
        print("jemx_image constructed as", ddosa.ii_skyimage())
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

    #This function uses varmosaic to create a mosaic image, see
    # mosaic_jemx_osa for the class that uses OSA to extract the image

    input_imagelist = JMXScWImageList

    # write_caches=[da.TransientCache,ddosa.MemCacheIntegralFallback,ddosa.MemCacheIntegralIRODS]
    # read_caches=[da.TransientCache,ddosa.MemCacheIntegralFallback,ddosa.MemCacheIntegralFallbackOldPath,ddosa.MemCacheIntegralIRODS]

    copy_cached_input = False

    maxsep = 30

    pixdivide = 2

    cached = True

    test_files = False

    version = "v1.2"

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

    def choose_list(self, coords):
        (ra, dec) = coords
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
            h = fits.open(fn)[2].header
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

                    fe = fits.open(imageobj.skyima.get_path())[2]
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
                        sens_keys = list(statistic.keys())
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

                image_file = fits.open(imagefilename)
                if stacked_file is None:
                    stacked_file = image_file
                else:
                    for i, e in enumerate(image_file[2:]):
                        stacked_file[i + 2].data += e.data

                # aligned
                flux = image_file[2].data
                var = image_file[3].data
                sig = flux/var**0.5
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
            fits.PrimaryHDU(f / v ** 0.5).writeto(stacked, overwrite=True)
            setattr(self, stacked, da.DataFile(stacked))

            stacked = "stacked_%.5lg_%.5lg.fits" % (ra, dec)
            stacked_file.writeto(stacked, overwrite=True)
            self.stacked.append([(ra, dec), da.DataFile(stacked)])
            setattr(self, stacked, da.DataFile(stacked))

            mosaic = "mosaic_%.5lg_%.5lg.fits" % (ra, dec)
            regionfile = mosaic.replace(".fits", ".reg")

            ht = ddosa.heatool(
                os.environ['COMMON_INTEGRAL_SOFTDIR'] + "/imaging/varmosaic/varmosaic_exposure/varmosaic", env=env_for_scw(self)
)
            ht['pixdivide'] = self.pixdivide
            ht['filelist'] = listfile
            ht['outimage'] = mosaic
            ht['outregion'] = regionfile
            ht.run()

            self.mosaics.append([(ra, dec), da.DataFile(mosaic)])
            setattr(self, mosaic, da.DataFile(mosaic))
            setattr(self, regionfile, da.DataFile(regionfile))

            self.skyima = da.DataFile(mosaic)  # store last

            ddosa.remove_withtemplate("jmx_sloc_res.fits(JMX1-SLOC-RES.tpl)")
            fn = "jmx_sloc_res.fits"
            ht = ddosa.heatool("j_ima_src_locator", env=env_for_scw(self)
)
            ht['inDOL'] = mosaic + "[2]"
            ht['sigDOL'] = mosaic + "[4]"
            ht['outFile'] = fn
            ht.run()

            self.skyres = da.DataFile(fn)
            self.srclres = da.DataFile(fn)


class mosaic_src_loc(ddosa.DataAnalysis):

    #This gets the image from mosaic_jemx, which is using varmosaic, not the mosaic_jemx_osa class
    input_mosaic=mosaic_jemx

    copy_cached_input=False

    def main(self):
        fn="jmx_sloc_res.fits"
        ht=ddosa.heatool("j_ima_src_locator", env=env_for_scw(self)
)
        ht['inDOL'] = self.input_mosaic.skyima.get_path()+"[2]"
        ht['sigDOL'] = self.input_mosaic.skyima.get_path() + "[3]"
        ht['outFile']=fn
        ht.run()

        self.sloc_res=da.DataFile(fn)


class jemx_spe_by_scw(graphtools.Factorize):
    root='jemx_spe'
    leaves=["ScWData",]

class jemx_lcr_by_scw(graphtools.Factorize):
    root='jemx_lcr'
    leaves=["ScWData",]

class JMXGroups(ddosa.DataAnalysis):
    input_scwlist=None
#    input_spe_processing=jemx_spe_by_scw
#    input_image_processing=jemx_image_by_scw
    input_jemx=JEMX

    allow_alias=True
    run_for_hashe=True

    outtype="BIN_I"

    def validate_child(self, kind, path):
        if kind == "lcr":
            lcrs = fits.open(path)[1].data
            print("found lightcurves:", lcrs)
            if len(lcrs) == 0:
                return False

        return True

    def construct_og(self,og_fn):
        scw_og_fns = []

        boundaries=[]

        for members in self.members:
            scw=members[0]
            if not os.path.exists(scw.scwpath + "/swg.fits"):
                print("\033[31mWARNING! scw is missing, skipping", scw.scwpath, "\033[0m")
                continue


            fn = "og_%s.fits" % scw.input_scwid.str()

            children=[
                scw.auxadppath + "/time_correlation.fits[AUXL-TCOR-HIS]",
            ]

            for m in members[1:]:
                for option in ['spe','arf','res','srclres','skyima','lcr']:
                    if hasattr(m,option):
                        if self.validate_child(option, getattr(m,option).get_path()):
                            children.append(getattr(m,option).get_path())

            ddosa.construct_gnrl_scwg_grp(scw, children=children, fn=fn)

                

            ddosa.import_attr(scw.scwpath + "/swg.fits",
                        ["OBTSTART", "OBTEND", "TSTART", "TSTOP", "SW_TYPE", "TELAPSE", "SWID", "SWBOUND"],fn)
            ddosa.set_attr({'ISDCLEVL': self.outtype}, fn)
            ddosa.set_attr({'INSTRUME': self.input_jemx.get_NAME()}, fn)

            h=fits.open(fn)[1].header
            boundaries.append([h['TSTART'],h['TSTOP']])


            scw_og_fns.append(fn)

        ddosa.construct_gnrl_scwg_grp_idx(scw_og_fns,fn="og_idx.fits")
        ddosa.set_attr({'ISDCLEVL': self.outtype}, "og_idx.fits")

        ddosa.construct_og(["og_idx.fits"], fn=og_fn)

        ddosa.set_attr({'INSTRUME': self.input_jemx.get_NAME()}, og_fn)

        ddosa.set_attr({'TSTART': min(list(zip(*boundaries))[0])}, og_fn)
        ddosa.set_attr({'TSTOP': max(list(zip(*boundaries))[0])}, og_fn)

        print("boundaries",boundaries)

    def main(self):
        self.members = []
        for scw in self.input_scwlist.scwlistdata:
            self.members.append([scw])

            for attachement in self.attachements:
                self.members[-1].append(da.byname(attachement).__class__(assume=[scw]))


        if len(self.members)==0:
            raise ddosa.EmptyScWList()


class JMXImageSpectraGroups(JMXGroups):
    input_scwlist=None
    input_spe_processing = jemx_spe_by_scw
    input_image_processing = jemx_image_by_scw

    attachements=['jemx_image','jemx_spe']

class JMXSpectraGroups(JMXGroups):
    input_scwlist=None
    input_spe_processing = jemx_spe_by_scw

    attachements=['jemx_spe']

class JMXLCGroups(JMXGroups):
    input_scwlist=None
    input_lcr_processing = jemx_lcr_by_scw

    attachements=['jemx_lcr']

class JMXImageGroups(JMXGroups):
    input_scwlist = None
    input_image_processing = jemx_image_by_scw

    attachements = ['jemx_image',]


class spe_pick(ddosa.DataAnalysis):
    input_spegroups = JMXSpectraGroups
    input_jemx=JEMX
    input_rmf=JRMF

    source_names=[]
    #source_names=["Crab"]

    cached=True
    
    version="v1.3"

    def get_version(self):
        try:
            return super().get_version()+"."+(".".join([m.replace(" ","_") for m in self.source_names]))
        except:
            return "spe_pick.UNDEFINED"

    def main(self):
        self.input_spegroups.construct_og("ogg.fits")

        dl=ddosa.heatool("dal_list", env=env_for_scw(self)
)
        dl['dol']="ogg.fits"
        dl.run()


        import glob
        print(glob.glob("*"))

        if len(self.source_names) != 0:
            source_names = self.source_names
        else:
            ddosa.remove_withtemplate("sources.fits(JMX%i-OBS.-RES.tpl)"%self.input_jemx.num)
            ht = ddosa.heatool("src_collect", env=env_for_scw(self)
)
            ht['group'] = "ogg.fits[1]"
            ht['instName']=self.input_jemx.get_name()
            ht['results']='sources.fits'
            ht.run()

            source_names  = list(set(fits.open('sources.fits')[1].data['NAME']))

        for source_name in source_names:
            if source_name == "NEW SOURCE":
                print("encountered NEW SOURCE: this makes no sense to merge")
                continue
            
            if source_name == "CONFUSED ID":
                print("encountered CONFUSED ID: this makes no sense to merge, and causes segfault!")
                continue
            
            if source_name == "":
                print("encountered empty name: this makes no sense to merge")
                continue

            print("will pick source:", source_name)

            sumname = "spec_%s" % source_name.replace(" ","_")
            singlename = source_name+"_JMX%i_single_pha2.fits"%self.input_jemx.num
            singlearfname = source_name+"_JMX%i_single_arf2.fits"%self.input_jemx.num

            ddosa.remove_withtemplate(singlename+"(JMX%i-PHA2-SPE.tpl)"%self.input_jemx.num)
            ddosa.remove_withtemplate(singlearfname+"(JMX%i-PHA2-ARF.tpl)"%self.input_jemx.num)
            ddosa.remove_withtemplate(sumname+"_pha.fits(JMX%i-PHA1-SPE.tpl)"%self.input_jemx.num)
            ddosa.remove_withtemplate(sumname+"_arf.fits(JMX%i-PHA1-ARF.tpl)"%self.input_jemx.num)
            

            ht = ddosa.heatool("spe_pick", env=env_for_scw(self)
)
            ht['group'] = "ogg.fits[1]"
            ht['source']=source_name
            ht['instrument']=self.input_jemx.get_name()
            ht['sumname']=sumname
            ht['single']='n'

            try:
                ht.run()
            except pilton.HEAToolException as e:
                if 'Error -35805 searching for source spectra in SWG' in ht.output:
                    raise CorruptScwInSpePick()
                

            import glob
            print(glob.glob("*"))

            fits.open(self.input_rmf.rmf.get_path()).writeto(sumname + "_rmf.fits", overwrite=True)

            if os.path.exists(sumname+"_pha.fits"):
                setattr(self,'spectrum_'+source_name,da.DataFile(sumname+"_pha.fits"))

            if os.path.exists(sumname+"_arf.fits"):
                setattr(self, 'arf_' + source_name, da.DataFile(sumname + "_arf.fits"))

            if os.path.exists(sumname+"_rmf.fits"):
                setattr(self, 'rmf_' + source_name, da.DataFile(sumname + "_rmf.fits"))

class lc_pick(ddosa.DataAnalysis):
    input_lcgroups = JMXLCGroups
    input_timebin = LCTimeBin
    input_jemx=JEMX

    source_names=[]
    source_ids=[]

    cached=True

    version="v1.3.5"

    def get_version(self):
        try:
            return super().get_version()+"."+(".".join([m.replace(" ","_") for m in self.source_names]))
        except:
            return "lc_pick.UNDEFINED"

    def main(self):
        self.input_lcgroups.construct_og("ogg.fits")

        dl=ddosa.heatool("dal_list", env=env_for_scw(self)
)
        dl['dol']="ogg.fits"
        dl.run()

        if len(self.source_names) != 0:
            source_names = self.source_names
            source_ids = self.source_ids
            
            sources_to_extract = list(set(list(source_names) + list(source_ids)))
            
        else:
            ddosa.remove_withtemplate("sources.fits(JMX%i-OBS.-RES.tpl)"%self.input_jemx.num)
            ht = ddosa.heatool("src_collect", env=env_for_scw(self)
)
            ht['group'] = "ogg.fits[1]"
            ht['instName']=self.input_jemx.get_name()
            ht['results']='sources.fits'
            ht.run()

            source_names  =list(fits.open('sources.fits')[1].data['NAME'])
            source_ids  =list(fits.open('sources.fits')[1].data['SOURCE_ID'])
            
            sources_to_extract = list(set(source_names + source_ids))

        print("found the following sources",source_names, source_ids)

        attached_files = []
        
        for source_name in sources_to_extract:
            
            if source_name == "NEW SOURCE":
                print("encountered NEW SOURCE: this makes no sense to merge")
                continue
            
            if source_name == "CONFUSED ID":
                print("encountered CONFUSED ID: this makes no sense to merge, and causes segfault!")
                continue
            
            if source_name == '' or source_name == 'UNKNOWN':
                print('Skipping source_name \'%s\': this makes no sense to merge' %(source_name))
                continue
            
            sumname = "lc_%s" % source_name.replace(" ","_").replace("+","_").replace("-","_")

            ddosa.remove_withtemplate(sumname+".fits")
            
            ht = ddosa.heatool("lc_pick", env=env_for_scw(self)
)
            ht['group'] = "ogg.fits[1]"
            ht['source']=source_name
            ht['instrument']=self.input_jemx.get_name()
            ht['lc']=sumname
            ht['emin']=""
            ht['emax']=""

            ht.run()

            # import glob
            # print(glob.glob("*"))

            #fits.open(self.input_rmf.rmf.get_path()).writeto(sumname + "_rmf.fits", overwrite=True)

            sumfile_fn = sumname + ".fits"
            if os.path.exists(sumfile_fn):
                setattr(self,sumname,da.DataFile(sumfile_fn))
                attached_files.append(sumfile_fn)

        for i, sumfile_fn in enumerate(attached_files):
            print("attached file: ", i, "/", len(attached_files), sumfile_fn)

        if self.input_timebin.time_bin_seconds < 0.1:
            self.comment = "please note that minimum time bin for jemx is 0.1 s; requesting smaller value produces lightcurve with 4 s bins (the default)"
        self.warning = "this is a warning"


class mosaic_jemx_osa(ddosa.DataAnalysis):
    input_groups = JMXImageGroups
    input_jemx = JEMX
    input_refcat = ddosa.GRcat
    input_ic = ddosa.ICRoot

    cached=True

    version="v2"

    def get_version(self):
        return super().get_version()

    def main(self):
        self.input_groups.construct_og("ogg.fits")

        dl=ddosa.heatool("dal_list", env=env_for_scw(self)
)
        dl['dol']="ogg.fits"
        dl.run()

        #ddosa.remove_withtemplate(fn+"(ISGR-SRC.-SPE-IDX.tpl)")

        env=deepcopy(os.environ)
        env['DISPLAY']=""

        fn_mosaic=self.input_jemx.get_name()+"_osa_mosaic.fits"

        ddosa.remove_withtemplate(self.input_jemx.get_name()+"_obs_res.fits")
        ddosa.remove_withtemplate(fn_mosaic)

        #It runs the mosaic production
        ht = ddosa.heatool("jemx_science_analysis",env=env_for_scw(self)
)
        ht['ogDOL'] = "ogg.fits"
        ht['IC_Group']=self.input_ic.icindex
        ht['jemxNum']=self.input_jemx.num
        ht['CAT_I_refCat'] = self.input_refcat.cat
        ht['startLevel']="IMA2"
        ht['endLevel'] = "IMA2"
        ht['skipLevels']=""
        ht['chatter']=5
        ht['IMA2_outfile']=fn_mosaic
        try:
            ht.run()
        except pilton.HEAToolException as ex:
            if 'segmentation violation' in ht.output:
                raise SegFaultInMosaic()
            
            if 'Segmentation fault' in ht.output:
                raise SegFaultInMosaic()
        
        #It finds sources in the mosaic
        fn = self.input_jemx.get_name()+"_sloc_res.fits"
        ddosa.remove_withtemplate(fn)

        ht=ddosa.heatool("j_ima_src_locator", env=env_for_scw(self)
)
        ht['inDOL'] = fn_mosaic+"[2]"
        ht['varDOL'] = fn_mosaic + "[3]"
        ht['sigDOL'] = fn_mosaic + "[4]"
        ht['detsigMin'] = "-3.0"
        #This parameter adds to the suggested detection significance an additional margin of 3 to reduce spurious sources
        ht['outFile']=fn
        ht.run()

        #It identifies sources in the catalog
        ht = ddosa.heatool("q_identify_srcs", env=env_for_scw(self)
)

        #Naming of instrument is different in this tool, assume it is jemx1 and change if it is jemx2
        instrument=4
        if self.input_jemx.get_name() == "jmx2":
            instrument=5

        ht['instrument']=instrument
        ht['srcl_cat_dol']= self.input_refcat.cat
        ht['srcl_res_dol']=fn
        ht.run()

        #It displays files for logging
        os.system("ls -ltor")

        #attribute of the mosaic sources for the output catalog
        if os.path.exists(fn):
            self.obsres = da.DataFile(fn)
            self.srclres = da.DataFile(fn)

        #attribute of the mosaic image
        if os.path.exists(fn_mosaic):
            self.skyima = da.DataFile(fn_mosaic)

#        setattr(self, 'arf_' + source_name, da.DataFile(sumname + "_arf.fits"))

import dataanalysis
import dataanalysis.callback

#dataanalysis.callback.default_callback_filter=CallbackRareDDOSAFilter

previously_accepted_classes=dataanalysis.callback.default_callback_filter.callback_accepted_classes

class CallbackRareDDOSAFilter(dataanalysis.callback.Callback):
    def extract_data(self,obj):
        data={'scwid':'inapplicable',}

        scw=obj.cache.get_scw(obj._da_locally_complete)
        
        expected_hashe=getattr(obj,'_da_expected_full_hashe',None)

        if expected_hashe is not None:
            data['node_id']=obj.cache.hashe2signature(expected_hashe)
        else:
            data['node_id']="undefined_expected_hashe_please_complain" # add sentry

        if scw is None:
            try:
                scw=obj.cache.get_scw(obj._da_expected_full_hashe)
            except:
                pass

        if scw is not None:
            data.update({"scwid":scw})

        return data

dataanalysis.callback.default_callback_filter=CallbackRareDDOSAFilter

if previously_accepted_classes is not None:
    dataanalysis.callback.default_callback_filter.set_callback_accepted_classes(previously_accepted_classes)

CallbackRareDDOSAFilter.set_callback_accepted_classes([mosaic_jemx_osa,mosaic_jemx,jemx_image,jemx_spe,jemx_lcr,spe_pick,lc_pick])
