import random
import math
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
plt.style.use('seaborn-poster')
# plt.style.use('ggplot')
import os
import json
import time
from datetime import datetime

class Abs_Scat():

    def __init__(self, comment = ""):
        
        print("init Abs_Scat()")
        comment = """
        Photon Weight = 1.0*solarSpectrum[wavelength]
        Bm = 0.61
        """.strip("\n")   
        self.scaleAlphaEm = 1.0
        self.scaleAlphaPh = 1.0
        self.scaleAlphaBase = 1.0
        self.scaleAlphaEpi = 1.0
        self.scaleAlphaDerm = .1
        self.scaleAlphaOxy = 1.0
        self.scaleAlphaDeoxy = 1.0
        self.scaleScattering = 1.0
        self.params = f"""
        scaleAlphaEm={self.scaleAlphaEm}, scaleAlphaPh={self.scaleAlphaPh}, scaleAlphaEpi={self.scaleAlphaEpi},scaleAlphaBase={self.scaleAlphaBase}
        scaleAlphaOxy={self.scaleAlphaOxy}, scaleAlphaDeoxy={self.scaleAlphaDeoxy}, scaleAlphaDerm={self.scaleAlphaDerm}, scaleScattering={self.scaleScattering}
        {comment}"""


    def get_params(self):
        return self.params

    #alpha_em is the spectral absorption coefficient of eumelanin,
    def get_alpha_em(self,wl):
        alpha_em = 6.6 * pow(10,10) * pow(wl,-3.33)
        #mm^-1
        #convert to cm^-1
        return alpha_em*self.scaleAlphaEm

    #alpha_pm is the spectral absorption coefficient of pheomelanin,
    def get_alpha_ph(self,wl):
        alpha_pm = 2.9 * pow(10,14) * pow(wl,-4.75)
        #mm^-1
        return alpha_pm*self.scaleAlphaPh

    #baseline absorption - used in both layers
    def get_alpha_base(self,wl):
        # alpha_base = 0.0244 + 8.54 * pow(10,-(wl-154)/66.2)
        alpha_base = 0.0244 + 8.53*np.exp(-(wl-154)/66.2)
        return alpha_base*self.scaleAlphaBase

    #get_alpha_epi_total is the total absorption coefficient of epidermis
    def get_alpha_epi_total(self,wl, Cm, Bm, Ch):
        alpha_em = self.get_alpha_em(wl)
        alpha_pm = self.get_alpha_ph(wl)
        alpha_base = self.get_alpha_base(wl)
        alpha_epi_total = Cm*(Bm*alpha_em + (1-Bm)*alpha_pm) + (1-Cm)*alpha_base
        #alpha_epi_total = Cm * alpha_em + Ch * alpha_pm + (1-Cm)*self.get_alpha_base(wl)
        #mm^-1
        return alpha_epi_total*self.scaleAlphaEpi

    #get_alpha_derm_total is the total absorption coefficient of dermis
    def get_alpha_derm_total(self,wl, Ch, Bm):
        """
        gamma  deoxy hemoglobin : oxy hemoglobin ratio
        typical 0.6−0.8 [ZBK01]. In
        our model we fix gamma to be 0.75
        """
        nm = wl
        gamma =0.75
        # Second layer absorption - Dermis
        alpha_oxy = SkinAbsorption.O2Hb[int(nm)]*self.scaleAlphaOxy
        alpha_deoxy = SkinAbsorption.Hb[int(nm)]*self.scaleAlphaOxy
        # A = 0.75 * alpha_oxy + (1 - 0.75) * alpha_deoxy      
        # alpha_derm_total = Ch * (A + (( 1 - Ch) * self.get_alpha_base(wl)))
        alpha_derm_total = Ch*(gamma*alpha_oxy + (1-gamma)*alpha_deoxy) + (1-Ch)*self.get_alpha_base(wl)
        return alpha_derm_total*self.scaleAlphaDerm

    #scatter_coef is the scattering coefficient of epidermis and dermis layers
    def get_scatter_coef(self,wl, layer):
        scattering = 14.74*math.pow(wl,-0.22)+2.22*math.pow(10,11)*math.pow(wl,-4.0)
        #scattering = 1.1*math.pow(10,12)*math.pow(wl,-4.0)+73.7*math.pow(wl,-0.22)
        if layer == "dermis":
            scattering = 0.5*scattering
        return scattering*self.scaleScattering

class SkinAbsorption:

    O2Hb = {}
    Hb = {}
    CIE_XYZ_Spectral_Sensitivity_Curve = {}
    step_sizes = []
    
    def __init__(self):
        #load csv to dictionary
        import csv
        solarSpectrum = {}
        with open('/Users/joeljohnson/Documents/Python/solarSpectrum.csv', 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] != 'wl':
                    solarSpectrum[float(row[0])] = float(row[1])
            f.close()
        self.solarSpectrum = solarSpectrum

        self.abs = Abs_Scat()
        cmf = open("/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/cie-cmf.csv", "r") 
        contents = cmf.readlines()
        CIE_CMF = [None] * 42
        for line in contents:
            data = line .rstrip("\n").split(",")
            CIE_CMF[int((int(data [0]) - 380) / 10)] = [float(data [1]),float(data [2]),float(data [3])]
        cmf.close()
        for nm in range(380,790,10):
            self.CIE_XYZ_Spectral_Sensitivity_Curve[nm] = CIE_CMF[int((nm - 380)/10)]
        print("init SkinAbsorption")
        HbFile = open('/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/hb.csv', "r")
        O2HbFile = open('/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/O2Hb.csv', "r")

        HbLines = HbFile.readlines()
        for line in HbLines:
            splitLine = line.split(",")
            self.Hb[int(splitLine[0])] = float(splitLine[1].rstrip("\n"))

        O2HbLines = O2HbFile.readlines()
        for line in O2HbLines:
            splitLine = line.split(",")
            self.O2Hb[int(splitLine[0])] = float(splitLine[1].rstrip("\n"))
            
        HbFile.close()
        O2HbFile.close()
    def gamma_correction(self, C):
                abs_C = abs(C)
                if abs_C > 0.0031308:
                    return 1.055 * pow(abs_C,1/2.4) - 0.055
                else:
                    return 12.92 * C
    def XYZ_to_sRGB(self, xyz):
        x = xyz[0]/10
        y = xyz[1]/10
        z = xyz[2]/10
        mat3x3 = [(3.2406, -1.5372, -0.4986), (-0.9689,   1.8758,  0.0415), (0.0557, -0.204,  1.057)]
        
        r = self.gamma_correction(x * mat3x3[0][0] + y * mat3x3[0][1] + z * mat3x3[0][2])
        g = self.gamma_correction(x * mat3x3[1][0] + y * mat3x3[1][1] + z * mat3x3[1][2])
        b = self.gamma_correction(x * mat3x3[2][0] + y * mat3x3[2][1] + z * mat3x3[2][2])
        sRGB = (int(r*255),int(g*255),int(b*255)) #needs to be 0 - 255 for outputing to color image
        return sRGB

    def Generate(self):
        CSV_Out = "Wavelength,Cm,Ch,Bm,Eumelanin,Pheomelanin,Baseline,Epidermis,ScatteringEpi,ScatteringDerm,Dermis,Reflectance\n"
        CSV_Out_File = open("/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/AbsorbtionValues.csv","w+")
        CSV_Out_File.write(CSV_Out)
        groupByBlend = []
        Bm = 0.61
        # for Bm in [0.01,0.5,0.99]:
        groupByMelanin = []
        for Cm in [0.002,0.0135,0.0425,0.1,0.185,0.32,0.5]:
            groupByHemoglobin = []
            for Ch in [0.003,0.02,0.07,0.16,0.32]:
                CSV_Out,values = self.GetReflectanceValues( Cm, Ch, Bm )
                groupByHemoglobin.append(values)
                
                CSV_Out_File.write(CSV_Out)
                
            groupByMelanin.append(groupByHemoglobin)
            groupByBlend.append(groupByMelanin)
        
        CSV_Out_File.close()
        XYZ_Colors = []
        spec = []
        
        for melanin_blend in groupByBlend: # Mb
            # print("Melanin Blend"+str(melanin_blend))
            for melanin_fraction in melanin_blend:
                for hemoglobin_fraction in melanin_fraction: # Mf
                    total = (0,0,0)
                    spec = [0.0] * 41
                    for data in hemoglobin_fraction: # Hf                            
                        reflectance = data[0]
                        spec[int((data[4] -380) / 10)] = data[0]
                        xyz = self.CIE_XYZ_Spectral_Sensitivity_Curve.get(data[4])
                        x = xyz[0] * reflectance
                        y = xyz[1] * reflectance
                        z = xyz[2] * reflectance
                        total = (total[0] + x, total[1] + y, total[2] + z)                     
                    XYZ_Colors.append(total)
        pixelsRGB = []

        for xyz in XYZ_Colors:
            pixelsRGB.append(self.XYZ_to_sRGB(xyz))

        img = Image.new('RGB', (21,5), color = 'black')
        pixels = img.load()
        
        img.save(f'/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/Plots/skin_LUT{date_time}.png',ppi=72)
        pixel_index = 0
        for m in range(7):
            for x in range(3): # 0 - 32
                for y in range(5): # 0 - 32
                    xCoord = x + (m * 3)
                    yCoord = y
                    pixels[xCoord, yCoord] = pixelsRGB[pixel_index]
                    pixel_index = pixel_index + 1
                    
        img.save(f'/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/Plots/skin_LUT{date_time}.png',ppi=72)
        return pixelsRGB
            
    def GetReflectanceValues(self, Cm, Ch, Bm):
        CSV_Out = ""
        wavelengths = range(380,780,10)
        reflectances = []
        
        for nm in wavelengths:
            SAV_eumelanin_L = self.abs.get_alpha_em(nm)

            SAV_pheomelanin_L = self.abs.get_alpha_ph(nm)

            baselineSkinAbsorption_L = self.abs.get_alpha_base(nm)

            epidermis = self.abs.get_alpha_epi_total(nm,Cm,Bm,Ch)

            dermis = self.abs.get_alpha_derm_total(nm,Ch,Bm)
           
            scattering_epidermis = self.abs.get_scatter_coef(nm,layer="epidermis")

            scattering_dermis = self.abs.get_scatter_coef(nm,layer="dermis")

            reflectance = self.MonteCarlo(epidermis, scattering_epidermis, dermis, scattering_dermis, nm)

            reflectances.append((reflectance,Cm,Ch,Bm,nm,epidermis,dermis,baselineSkinAbsorption_L))

            CSV_Out = f"{nm},{Cm},{Ch},{Bm},{SAV_eumelanin_L},{SAV_pheomelanin_L},{baselineSkinAbsorption_L},{epidermis},{scattering_epidermis},{scattering_dermis},{dermis},{reflectance}\n"
        return CSV_Out, reflectances
        

    def MonteCarlo (self, epi_mua, epi_mus, derm_mua, derm_mus, nm):
        # These are our Monte Carlo Light Transport Variables that don't change
        Nbins = 100
        Nbinsp1 = 101
        PI = 3.1415926
        LIGHTSPEED = 2.997925 * pow(10,10)
        ALIVE = 1
        DEAD = 0
        THRESHOLD = 0.01
        CHANCE = 0.1*3.1415926
        COS90D = 1 * pow(10,-6)
        ONE_MINUS_COSZERO = 1 * pow(10,-12)
        COSZERO = 1.0 - 1.0e-12     # cosine of about 1e-6 rad
        g = 0.90
        nt = 1.33 # Index of refraction
        epidermis_thickness = 0.25

        x = 0.0
        y = 0.0
        z = 0.0 # photon position

        ux = 0.0
        uy = 0.0
        uz = 0.0 # photon trajectory as cosines

        uxx = 0.0
        uyy = 0.0
        uzz = 0.0 # temporary values used during SPIN

        s = 0.0 # step sizes. s = -log(RND)/mus [cm] 
        costheta = 0.0 # cos(theta) 
        sintheta = 0.0 # sin(theta) 
        cospsi = 0.0 # cos(psi) 
        sinpsi = 0.0 # sin(psi) 
        psi = 0.0 # azimuthal angle 
        i_photon = 0.0 # current photon
        W = 0.0 # photon weight 
        absorb = 0.0 # weighted deposited in a step due to absorption 
        photon_status = 0.0 # flag = ALIVE=1 or DEAD=0 
        ReflBin = [None] * Nbinsp1 #bin to store weights of escaped photos for reflectivity
        epi_albedo = epi_mus/(epi_mus + epi_mua) # albedo of tissue
        derm_albedo = derm_mus/(derm_mus + derm_mua) # albedo of tissue
        Nphotons = 100 # number of photons in simulation 
        NR = Nbins # number of radial positions 
        radial_size = 2.5 # maximum radial size 
        r = 0.0 # radial position 
        dr = radial_size/NR; # cm, radial bin size 
        ir = 0 # index to radial position 
        shellvolume = 0.0 # volume of shell at radial position r 
        CNT = 0.0 # total count of photon weight summed over all bins 
        rnd = 0.0 # assigned random value 0-1 
        u = 0.0
        temp = 0.0 # dummy variables

        # Inits
        random.seed(0)
        RandomNum = random.random()
        for i in range(NR+1):
            ReflBin[i] = 0

        while True:
            i_photon = i_photon + 1

            W = 1.0*self.solarSpectrum[nm]
            photon_status = ALIVE

            x= 0
            y = 0
            z = 0

            #Randomly set photon trajectory to yield an isotropic source.
            costheta = 2.0 * random.random() - 1.0

            sintheta = math.sqrt(1.0 - costheta*costheta)
            psi = 2.0 * PI * random.random()
            ux = sintheta * math.cos(psi)
            uy = sintheta * math.sin(psi)
            uz = (abs(costheta)) # on the first step we want to head down, into the tissue, so > 0

            # Propagate one photon until it dies as determined by ROULETTE.
            # or if it reaches the surface again
            it = 0
            max_iterations = 100 # to help avoid infinite loops in case we do something wrong

            # we'll hit epidermis first, so set mua/mus to those scattering/absorption values
            mua = epi_mua
            mus = epi_mus
            albedo = epi_albedo
            while True:        
                it = it + 1
                rnd = random.random()
                while rnd <= 0.0: # make sure it is > 0.0
                    rnd = random.random()
                s = -math.log(rnd)/(mua + mus)
                np.append(self.step_sizes, s)
                x = x + (s * ux)
                y = y + (s * uy)
                z = z + (s * uz)

                if uz < 0:
                    # calculate partial step to reach boundary surface
                    s1 = abs(z/uz)
                    # move back
                    x = x - (s * ux)
                    y = y - (s * uy)
                    z = z - (s * uz)
                    # take partial step
                    x = x + (s1 * ux)
                    y = y + (s1 * uy)
                    z = z + (s1 * uz)

                    # photon is now at the surface boundary, figure out how much escaped and how much was reflected
                    internal_reflectance = self.RFresnel(1.0,nt, -uz )

                    #Add weighted reflectance of escpaed photon to reflectance bin
                    external_reflectance = 1 - internal_reflectance
                    r = math.sqrt(x*x + y*y)
                    ir = (r/dr)
                    if ir >= NR:
                       ir = NR  
                    ReflBin[int(ir)] = ReflBin[int(ir)] + (W * external_reflectance)
                                                                    
                    # Bounce the photon back into the skin
                    W = internal_reflectance * W
                    uz = -uz
                    x = (s-s1) * ux
                    y = (s-s1) * uy
                    z = (s-s1) * uz

                # check if we have passed into the second layer, or the first
                if z <= epidermis_thickness:
                   mua = epi_mua
                   mus = epi_mus
                   albedo = epi_albedo
                else:
                   mua = derm_mua
                   mus = derm_mus
                   albedo = derm_albedo

                ''' DROP '''
                absorb = W*(1 - albedo)
                W = W - absorb

                ''' SPIN '''
                # Sample for costheta 
                rnd = random.random()
                if (g == 0.0):
                   costheta = 2.0*rnd - 1.0
                else:
                   temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd)
                   costheta = (1.0 + g*g - temp*temp)/(2.0*g)
                sintheta = math.sqrt(1.0 - costheta*costheta) 

                # Sample psi. 
                psi = 2.0*PI*random.random()
                cospsi = math.cos(psi)
                if (psi < PI):
                   sinpsi = math.sqrt(1.0 - cospsi*cospsi)
                else:
                   sinpsi = -math.sqrt(1.0 - cospsi*cospsi)

                # New trajectory. 
                if (1 - abs(uz) <= ONE_MINUS_COSZERO) :      # close to perpendicular. 
                   uxx = sintheta * cospsi
                   uyy = sintheta * sinpsi
                   uzz = costheta * self.SIGN(uz)   # SIGN() is faster than division. 
                else: # usually use this option 
                   temp = math.sqrt(1.0 - uz * uz)
                   uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta
                   uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta
                   uzz = -sintheta * cospsi * temp + uz * costheta

                # Update trajectory 
                ux = uxx
                uy = uyy
                uz = uzz 
            
                # Check Roulette
                if (W < THRESHOLD):
                    if (random.random() <= CHANCE):
                        W = W / CHANCE
                    else:
                        photon_status = DEAD
                if photon_status is DEAD:
                    break
                if it > max_iterations:
                    break
            
            if i_photon >= Nphotons:
                break
        total_reflection = 0.0
        for each in range(NR+1):
            total_reflection = total_reflection + ReflBin[each]/Nphotons
        return total_reflection

  
    def RFresnel(self, n1, n2, cosT1):
        r = 0.0
        cosT2 = 0.0
        COSZERO = 1.0 - 1.0e-12 
        COS90D = 1 * pow(10,-6)
        if n1 == n2: #matched boundary
            r = 0.0
            cosT2 = cosT1
        elif cosT1 > COSZERO:     # normal incident
            cosT2 = 0
            r = (n2-n1)/(n2+n1)
            r *= r
        elif cosT1 < COS90D:      # very slant
            cosT2 = 0.0
            r = 1.0
        else: #general
            sinT1 = math.sqrt(1 - cosT1*cosT1)
            sinT2 = n1 * sinT1/n2
            
            if sinT2 >= 1.0:
                r = 1.0
                cosT2 = 0.0
            else:
                cosT2 = math.sqrt(1 - sinT2 * sinT2)
                cosAP = cosT1*cosT2 - sinT1*sinT2
                cosAM = cosT1*cosT2 + sinT1*sinT2
                sinAP = sinT1*cosT2 + cosT1*sinT2
                sinAM = sinT1*cosT2 - cosT1*sinT2
                r = 0.5 * sinAM * sinAM*(cosAM*cosAM+cosAP*cosAP)/(sinAP*sinAP*cosAM*cosAM)
        return r

    def SIGN(self, x):
        if x >=0:
            return 1
        else:
            return 0 

if __name__ == "__main__":
    # abs = Abs_Scat(scaleAlphaEm=10,scaleAlphaPh=10,scaleAlphaOxy=10,scaleAlphaDeoxy=10,scaleScattering=10)
    #get date and time 
        #tick marks adjusted to match the image
    Ch = np.array([0.003,0.02,0.07,0.16,0.32]) #y axis
    #reverse the array
    Ch = Ch[::-1]
    Cm = np.array( [0.002,0.0135,0.0425,0.1,0.185,0.32,0.5]) #x axis
    xCm = np.linspace(0,0.5,len(Cm))
    xCh = np.linspace(0.32,0,len(Ch))
    plt.xticks(xCm,Cm)
    plt.yticks(xCh,Ch)

    start = datetime.now()
    date_time = start.strftime("%m-%d-%Y_%H-%M-%S")
    date_time = date_time.replace(" ", "_").replace(":", "_").replace("/", "_")
   

    start_time = date_time
    print("Starting"+str(start_time))
    skinAbsorption = SkinAbsorption()
    title = skinAbsorption.abs.get_params()
    reflectanceValues = skinAbsorption.Generate()
    finish = datetime.now()
    finish_time = finish.strftime("%m-%d-%Y_%H-%M-%S").replace(" ", "_").replace(":", "_").replace("/", "_")
    print("Finished"+str(finish_time))
    print("Time taken: " + str(finish - start))
    #load image to np array
    im = Image.open(f"/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/Plots/skin_LUT{date_time}.png")
    im = np.array(im)
    #create 3 pictures from the image
    im1 = im[:,0:7,:]
    im2 = im[:,7:14,:]
    im3 = im[:,14:21,:]
    im1 = Image.fromarray(im1)
    im1 = im1.resize((64,46), Image.Resampling.LANCZOS)
    
    im1 = np.array(im1)
    plt.imshow(im1,extent=[0,0.5,0.32,0])
    plt.title(f"image1_{title}",fontsize=8)
    plt.ylabel("hemoglobin Chb",fontsize=8)
    plt.xlabel("melanin Cm",fontsize=8)

    plt.tick_params(axis='both', which='major', labelsize=8)
    #make figure bigger

    #save fig
    plt.savefig(f"/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/Plots/image1_{date_time}.png",dpi=72, pad_inches=1.5)

    im2 = Image.fromarray(im2)
    im2 = im2.resize((64,46), Image.Resampling.LANCZOS)
    im2 = np.array(im2)
    plt.imshow(im2,extent=[0,0.5,0.32,0])
    plt.title(f"image2_{title}",fontsize=8)
    plt.ylabel("hemoglobin Chb",fontsize=8)
    plt.xlabel("melanin Cm",fontsize=8)
    plt.tick_params(axis='both', which='major', labelsize=8)
    #save fig
    plt.savefig(f"/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/Plots/image2_{date_time}.png",dpi=72)
    
    im3 = Image.fromarray(im3)
    im3 = im3.resize((64,46), Image.Resampling.LANCZOS)
    im3 = np.array(im3)
    plt.imshow(im3,extent=[0.003,0.5,0.32,0])
    plt.title(f"image3_{title}",fontsize=8)
    plt.ylabel("hemoglobin Chb",fontsize=8)
    plt.xlabel("melanin Cm",fontsize=8)

    #save fig
    plt.savefig(f"/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/Plots/image3_{date_time}.png",dpi=72)


