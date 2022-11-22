class SkinAbsorption:

    O2Hb = {}
    Hb = {}

    def __init__(self):
        print("init SkinAbsorption")
        import pandas as pd
        import numpy as np
        df = pd.read_csv("/Users/joeljohnson/Documents/Github/CG_Research/Data/hem.txt", sep='\t')


        #get data from df where 380<=WL<=780
        hem = df.loc[(df.iloc[:,0] >= 380) & (df.iloc[:,0] <= 780)]
        #count by 10
        hem = hem.iloc[::5]
        for i in range(0,len(hem)):
            hem.iloc[i,0] = hem.iloc[i,0].astype(int)
            hem.iloc[i,1] = hem.iloc[i,1].astype(float)
            hem.iloc[i,2] = hem.iloc[i,2].astype(float)
            
        for i in range(0,len(hem)):
            wl= hem.iloc[i,0].astype(int)
            hbo2 = hem.iloc[i,1].astype(float)
            hb = hem.iloc[i,2].astype(float)
            self.O2Hb[wl] = hbo2
            self.Hb[wl] = hb

    def Generate(self):
        CSV_Out = "Wavelength,Eumelanin,Pheomelanin,Baseline,Epidermis,ScatteringEpi,ScatteringDerm,Dermis\n"
        groupByBlend = []
        for Bm in [0.01,0.5,0.99]:
            groupByMelanin = []
            for Cm in [0.002,0.0135,0.0425,0.1,0.185,0.32,0.5]:
                groupByHemoglobin = []
                for Ch in [0.003,0.02,0.07,0.16,0.32]:
                    values = self.GetReflectanceValues( Cm, Ch, Bm )
                    CSV_Out = CSV_Out + values
                    groupByHemoglobin.append(values)
                groupByMelanin.append(groupByHemoglobin)
            groupByBlend.append(groupByMelanin)
        CSV_Out_File = open("/Users/joeljohnson/Documents/Github/CG_Research/PhysicalModels/AbsorptionValues.csv","w+")
        CSV_Out_File.write(CSV_Out)
        CSV_Out_File.close() 
        print(groupByHemoglobin)
        return groupByBlend
    
    def GetReflectanceValues(self, Cm, Ch, Bm):
        CSV_Out = ""
        wavelengths = range(380,790,10)
        for nm in wavelengths:
            
            # First layer absorption - Epidermis
            SAV_eumelanin_L = 6.6 * pow(10,10) * pow(nm,-3.33)
            SAV_pheomelanin_L = 2.9 * pow(10,14) * pow(nm,-4.75)
            epidermal_hemoglobin_fraction = Ch * 0.25
            
            # baseline - used in both layers
            baselineSkinAbsorption_L = 0.0244 + 8.54 * pow(10,-(nm-220)/123)

            # Epidermis Absorption Coefficient:          
            epidermis = Cm * ((Bm * SAV_eumelanin_L) +((1 -  Bm ) * SAV_pheomelanin_L)) + ((1 - epidermal_hemoglobin_fraction) * baselineSkinAbsorption_L)
            
            # Second layer absorption - Dermis
            gammaOxy = self.O2Hb[int(nm)]
            gammaDeoxy = self.Hb[int(nm)]
            A = 0.75 * gammaOxy + (1 - 0.75) * gammaDeoxy      
            dermis = Ch * (A + (( 1 - Ch) * baselineSkinAbsorption_L))
            
            # Scattering coefficients
            scattering_epidermis = pow(14.74 * nm, -0.22) + 2.2 * pow(10,11) * pow(nm, -4)
            scattering_dermis = scattering_epidermis * 0.5
            
            # Index of refraction - we assume a constant value of 1.4 as it changes very little through the skin
            NR = 1.4
            CSV_Out = CSV_Out + f"{nm},{SAV_eumelanin_L},{SAV_pheomelanin_L},{baselineSkinAbsorption_L},{epidermis},{scattering_epidermis},{scattering_dermis},{dermis}\n"
            if(nm == 380):
                print(f"Absorption at 380nm: {epidermis}")
                print(scattering_dermis+baselineSkinAbsorption_L+epidermis)
        return CSV_Out



if __name__ == "__main__":
    print("main")
    skinAbsorption = SkinAbsorption()
    reflectanceValues = skinAbsorption.Generate()