# This script will be used to run the BioCro-GRN-FBA model
# Originally written by Sanu Shameer 2021, edited by Hilary Hunt 2023

# Import necessary packages
from shutil import copyfile
import os
import subprocess
import pandas as pd
from core import WholePlantModel 
from cobra import Metabolite
import cobra

# Run a BioCro simulation using an FBA model of soybean growth to determine the relative growth of each organ based on the carbon partitioning determined by BioCro
    ## If start_afresh is 1, temporary BioCro-FBA files from previous runs are deleted before starting and the model will simulate soybean growth from the initial parameters specified in the inputs
    ## weather_file should be the path to a file of weather experienced by the soybeans. See BioCro documentation for specifics
    ## kOrgans_vals_path points to a biocro run with carbon partitioning between organs - we should try to eliminate the need for this
# @profile
def run_BioCro_FBA(start_afresh=0,save_tag='',weather_file="./../yamls/parameter_files/soybean_weather_2002growingseason.txt",kOrgan_vals_path="./../BioCro3_output_complete.csv"):
    if start_afresh:
         temp=subprocess.call(['rm','SoyFluxes'+save_tag+'.csv'])
         temp=subprocess.call(['rm','SoyFluxes2'+save_tag+'.csv'])
         temp=subprocess.call(['rm','./../yamls/tempWeather.txt'])
         temp=subprocess.call(['rm','./../yamls/tempInitialState.txt'])

    LoopCounter = 1

    # Initial mass of seedlings in BioCro
    LeafInitMass = 0.1199
    StemInitMass = 0.0150
    RootInitMass = 0.0150
    GrainInitMass = 0
    ShellInitMass = 0

    density = 10000*20/0.38   # plants/hectare ; 20 plants per m with 0.38m row spacing
    MW_CO2 = 44 # MW of CO2
    #canopy_assimilation_rate = biocro_vars['canopy_assimilation_rate']*(10**12)/(3600*biocro_vars['Sp']*density*Leaf*MW_CO2)   # see BioCro_outputs_units.txt for details
    #canopy_assimilation_rate = 25 + ((-3.632392 - 29.64083)/(1 + ((biocro_vars['solar']/273.5999)**1.248568)))   #A curve that fits Soy data

    LeafMass = [LeafInitMass]
    StemMass = [StemInitMass]
    RootMass = [RootInitMass]
    GrainMass = [GrainInitMass]
    ShellMass = [ShellInitMass]

    # If there is no file with carbon allocation values (kOrgan), generate one with BioCro
    if not kOrgan_vals_path:
         print('Generating kOrgan vals...')

    # Check if code broke on a previous iteration
    # This file is also an input for BioCro
    if not os.path.isfile("./../yamls/tempInitialState.txt"):
        copyfile("./../yamls/parameter_files/soybean_initial_state_fba.txt","./../yamls/tempInitialState.txt")

    # Load weather file
    fin = open(weather_file)
    # Save first row of the weather file as labels
    labels=fin.readline()
    match_wthr = False
    
    # Iterate through each line of weather data
    weather_list=list(fin)
    for it,line in enumerate(weather_list):
        # If the code has previously broken, run through weather data until the same point is reached...
        if not match_wthr and os.path.isfile("./../yamls/tempWeather.txt"):
            fin2 = open("./../yamls/tempWeather.txt")
            temp = fin2.readline() # Remove labels line
            # Check if the weather files match up
            line2 = fin2.readline()
            if line.replace('\n','') in line2.replace('/n',''):
                match_wthr = True
            # If match_wthr is still not true after going through every line then the current line is not in the tempweather file and something is wrong...
            if not match_wthr:
                # print("No match skipping")
                LoopCounter = LoopCounter+1
                LeafMass.append(0)
                StemMass.append(0)
                RootMass.append(0)
                GrainMass.append(0)
                ShellMass.append(0)
                continue
        fin3 = open("./../yamls/tempInitialState.txt")
        temp = fin3.readline()
        lineparts=fin3.readline().split()
        # If this is not the first loop, set initial state
        if LoopCounter!=1 and LeafMass[LoopCounter-1]+StemMass[LoopCounter-1]+RootMass[LoopCounter-1]+GrainMass[LoopCounter-1]==0:
            LeafMass[LoopCounter-1] = float(lineparts[0])/(0.01*20*(1.0/0.38))
            StemMass[LoopCounter-1] = float(lineparts[1])/(0.01*20*(1.0/0.38))
            RootMass[LoopCounter-1] = float(lineparts[2])/(0.01*20*(1.0/0.38))
            GrainMass[LoopCounter-1] = float(lineparts[3])/(0.01*20*(1.0/0.38))
            ShellMass[LoopCounter-1] = float(lineparts[4])/(0.01*20*(1.0/0.38))
        
        # Save weather in case of later crash
        # This is also an input in BioCro
        # print("Not skipping")
        # print("=============="+str(LoopCounter)+"================")
        match_wthr = True
        fout = open("./../yamls/tempWeather.txt","w")
        fout.write(labels+"\n"+line+"\n"+weather_list[it+1])
        fout.close()

        # # Run BioCro via yggdrasil
        # parser = argparse.ArgumentParser(description='Run an integration.')
        # parser.add_argument('yamlfile',nargs='+',help='One or more yaml specification files.')
        # args1 = parser.parse_args(['../yamls/biocro-fba.yml'])
        # runner.run(args1.yamlfile)

        # Run BioCro via rscript instead...
        # print('Pre-call!')
        temp = subprocess.call(['Rscript','../Gro_wrapper_fba_new.R','../yamls/tempInitialState.txt','../yamls/parameter_files/soybean_parameters_biocro3_fba.txt','../yamls/tempWeather.txt','../yamls/parameter_files/soybean_ss_modules_fba.txt','../yamls/parameter_files/soybean_deriv_modules_fba.txt'])

        # Read BioCro output
        biocro_vars=dict()
        df = pd.read_csv('./biocro_output2.txt',sep='\t')
        for col in df.columns:
            biocro_vars[str(col)]=float(df[col])

        
        Leaf = LeafMass[LoopCounter-1]
        Stem = StemMass[LoopCounter-1]
        Root = RootMass[LoopCounter-1]
        Seed = GrainMass[LoopCounter-1]
        Shell = ShellMass[LoopCounter-1]
        if str(biocro_vars['year'])!='nan':
            year = biocro_vars['year']
        else:
            year = 2002
        day = biocro_vars['doy']
        hour = biocro_vars['hour']
        day_length = biocro_vars['day_length']
        solar = biocro_vars['solar']
        lai = biocro_vars['lai']
        canopy_assimilation_rate = biocro_vars['canopy_assimilation_rate']/(3600 * 1e-6 * 30 * 1e-6 * 1e4 * lai)
        SLA = biocro_vars['Sp']/100 #converting ha/Mg to m2/Mg

        # kOrgan valus are imported from BioCro only results??
        df_orig_biocro = pd.read_csv(kOrgan_vals_path,sep=",")
        temp = df_orig_biocro[df_orig_biocro["year"]==year]
        temp2 = temp[temp["doy"]==day]
        temp3 = temp2[round(temp2["hour"])==round(hour)]
        # print(year)
        # print(day)
        # print(hour)
        # print(temp3)
        kLeaf = float(temp3["kLeaf"])#biocro_vars['kLeaf']
        kStem = float(temp3["kStem"])#biocro_vars['kStem']
        kRoot = float(temp3["kRoot"])#biocro_vars['kRoot']
        kSeed = float(temp3["kGrain"])#biocro_vars['kGrain']
        kShell = float(temp3["kShell"])#biocro_vars['kGrain']
        # Leaf_senesced = round(biocro_vars["senescence_leaf"],5)
        # Stem_senesced = round(biocro_vars["senescence_stem"],5)
        # Root_senesced = round(biocro_vars["senescence_root"],5)
        Leaf_senesced = round(biocro_vars["kSeneLeaf"]*biocro_vars["Leaf"],5)
        Stem_senesced = round(biocro_vars["kSeneStem"]*biocro_vars["Stem"],5)
        Root_senesced = round(biocro_vars["kSeneRoot"]*biocro_vars["Root"],5)
        # print("######### CONFIRM KVALUES #####")
        # print("kLeaf,kStem,kRoot,kSeed,kShell =")
        # for att in [kLeaf,kStem,kRoot,kSeed,kShell]:
        #     print(att)
        #Leaf = biocro_vars['Leaf']
        #Stem = biocro_vars['Stem']
        #Root = biocro_vars['Root']
        #Seed = biocro_vars['Grain']
        
        # Load or create 'plant' - the FBA model 
        savefile = "testSave"+save_tag+".pkl"
        if LoopCounter==1:
            plant = initWPM(save_tag=save_tag)
        # if (LoopCounter>1) and ("plant" not in locals()):
        #     # print('load')
        #     plant = loadWPM(savefile)
        # elif 'plant' in locals():
        #     print('notload','plant' in globals(), 'plant' in locals())

        plant.parameters["leafMaint"]=(((0.0049*solar)+2.7551)*3600*0.001)/(biocro_vars['Sp']*10000)  # Maintenance equation from Nadine's paper converted to mg/gDW/hr
        # print(plant.stoichiometricModel.reactions.query("GDWsink"))
        #plant.parameters["leafMaint"]=((0.0049*solar)+2.7551)  # Maintenance equation temporary
        if solar > 200:
            #plant.parameters["leafMaint"]=0.001        #mmol/gDW/hr
            plant.parameters["stemMaint"]=0.0005        #mmol/gDW/hr
            plant.parameters["shellMaint"]=0.0005     #mmol/gDW/hr
            plant.parameters["seedMaint"]=0.0005        #mmol/gDW/hr
            plant.parameters["rootMaint"]=0.0005       #mmol/gDW/hr
        else:
            plant.parameters["leafMaint"]=0.0005        #mmol/gDW/hr
            plant.parameters["stemMaint"]=0.0005        #mmol/gDW/hr
            plant.parameters["shellMaint"]=0.0005     #mmol/gDW/hr
            plant.parameters["seedMaint"]=0.0005        #mmol/gDW/hr
            plant.parameters["rootMaint"]=0.0005       #mmol/gDW/hr
        # N uptake/fixation
        plant=setNfixationUptakeRatio(plant,9)        
        # Output info on plant. Probably unneccessary outside of debugging. Delete before publishing
        # inspectWPM(plant)

        # Run the FBA model
        plant = runWPM(year,day,hour,plant,solar,canopy_assimilation_rate,day_length,kLeaf,kStem,kRoot,kSeed,kShell,SLA,Leaf,Stem,Root,Seed,Shell,LoopCounter,Leaf_senesced,Stem_senesced,Root_senesced)
        # If there's no growth, quit
        # if plant.recentSolution['Total_biomass_tx']==0 and plant.recentSolution['ATPase_tx_L1']==0:# and (LoopCounter != 440 and LoopCounter != 1220 and LoopCounter != 1244)):# or LoopCounter < 983: #Problems at 440, 2185,2283,2349,2350,2351
        #     # writeFluxes(plant,"SoyFluxes.csv",LoopCounter)
        #     print("No growth or NGAM",LoopCounter)
            # exit()
        # Save plant model to file
        plant.save(savefile)
        # print(plant.seeds)
        (Ltotal,Sttotal,Rtotal,Shtotal,Setotal) = analyseGrowth(plant)
        if len(plant.seeds)>0:
            writeFluxes(plant,"SoyFluxes2"+save_tag+".csv",LoopCounter)
        else:
            writeFluxes(plant,"SoyFluxes"+save_tag+".csv",LoopCounter)
        # Print out a lot of things like the leaf respiration and photosynthesis as well as growth per organ)
        # print("CO2_tx_L=")
        # print(plant.recentSolution['CO2_tx_L1'])
        # print("Photon_tx_L=")
        # print(plant.recentSolution['Photon_tx_L1'])
        # print("Leaf growth (g/plant/hr)="+str(Ltotal))
        # print("Stem growth (g/plant/hr)="+str(Sttotal))
        # print("Root growth (g/plant/hr)="+str(Rtotal))
        # print("Seed growth (g/plant/hr)="+str(Setotal))
        # print("Shell growth (g/plant/hr)="+str(Shtotal))
        # print("Terminate on purpose")
        # if len(plant.seeds)>0:
        #     # print("Terminate because seed appeared")
            #exit()
        if round(biocro_vars["kSeneLeaf"]*biocro_vars["Leaf"],3)>0 or round(biocro_vars["kSeneStem"]*biocro_vars["Stem"],3)>0:
            print("leaf/stem starts senescing")
            exit()

        # Add growth to 'Mass' variables
        LeafMass.append(LeafMass[LoopCounter-1]+Ltotal)
        StemMass.append(StemMass[LoopCounter-1]+Sttotal)
        RootMass.append(RootMass[LoopCounter-1]+Rtotal)
        GrainMass.append(GrainMass[LoopCounter-1]+Setotal)
        ShellMass.append(ShellMass[LoopCounter-1]+Setotal)

        # Write changes into BioCro output dataframe
        # What units are these?? Should be mmol/plant/hr -> g/hectare
        df["Leaf"]=[LeafMass[LoopCounter]*0.01*20*(1.0/0.38),]
        df["Stem"]=[StemMass[LoopCounter]*0.01*20*(1.0/0.38),]
        df["Root"]=[RootMass[LoopCounter]*0.01*20*(1.0/0.38),]
        df["Grain"]=[GrainMass[LoopCounter]*0.01*20*(1.0/0.38),]
        df["Shell"]=[ShellMass[LoopCounter]*0.01*20*(1.0/0.38),]

        #write initial state for the next cycle
        df_initial_state = df.copy()
        df_initial_state=df[["Leaf","Stem","Root","Grain","Shell","LeafLitter","RootLitter","StemLitter","soil_water_content","cws1","cws2","DVI","TTc","Rhizome","RhizomeLitter"]].copy()
        df_initial_state.to_csv("./../yamls/tempInitialState.txt",sep="\t",index=False)

        # write output
        if LoopCounter==1:
            df.to_csv("./../yamls/BioCro_output_complete.txt",sep="\t",index=False)
        else:
            # print("Reading output to apppend")
            df_complete = pd.read_csv("./../yamls/BioCro_output_complete.txt",sep="\t")
            df_complete = pd.concat([df_complete,df],ignore_index=True)
            df_complete.to_csv("./../yamls/BioCro_output_complete.txt",sep="\t",index=False)
            # print(df_complete.tail(1))
        LoopCounter=LoopCounter+1
    # cobra.io.write_sbml_model(plant.stoichiometricModel,"Model4Testing_3.xml")
    return

# Initialise a model of a whole soybean plant. This function will set the parameters used for the entire soybean growth simulation such as nitrate uptake, maintenance costs, etc. Parameters should be changed in the code. NOTE: It might be worth changing this to input parameters from a csv to make it easier for people to change...
    ## init_BioCro_run points to the results of a BioCro simulation to get the initial growth rate of each organ (this could also be exported to a csv...)
def initWPM(init_BioCro_run="./../BioCro3_output_complete.csv",save_tag=''):
    BioCro_results = pd.read_csv(init_BioCro_run)
    # print("Initializing WPM")
    plant = WholePlantModel()
    plant.mode="hourly"
    plant.parameters["singlePhloem"] = True
    plant.phloemComposition = "free"
    plant.symbiosis="N"
    plant.generateStoichiometricModels(sbmlFile="./PlantCoreMetabolism_v2_0_0.xml")
    plant.calculateBiomassWeight()
    plant.parameters["TimeOfEmergence"] = 6*24 #can go upto 15 days based on https://www.canr.msu.edu/news/identifying_and_responding_to_soybean_emergence_problems
    plant.parameters["NO3UptakeRate"]=0.4803587937919218
    plant.cycleThreshold=1000
    plant.parameters["nitrateConc"]=1
    plant.parameters["leafLifeSpan"] = 10000000
    plant.parameters["NO3UptakeVm"]=133.44
    plant.parameters["NO3UptakeKm"]=0.16
    plant.parameters["leafMaint"]=0        #mmol/gDW/hr
    plant.parameters["stemMaint"]=0        #mmol/gDW/hr
    plant.parameters["shellMaint"]=0        #mmol/gDW/hr
    plant.parameters["seedMaint"]=0        #mmol/gDW/hr
    plant.parameters["rootMaint"]=0       #mmol/gDW/hr
    plant.parameters["leafGrowthDuration"] = 27
    plant.parameters["leafVarA"] = 0.85766706
    plant.parameters["leafVarB"] = 0.3996402
    plant.parameters["leafVarC"] = 0.14332812
    plant.parameters["stemGrowthDuration"] = 29
    plant.parameters["stemVarA"] = 0.8608429
    plant.parameters["stemVarB"] = 0.44274541
    plant.parameters["stemVarC"] = 0.14195703
    plant.parameters["stemGrowthDuration"] = 1000000
    plant.symbiosis=["N",]
    plant.parameters["RbInit"] = 0.035041719#BioCro_results["Root_plant"][0]
    plant.parameters["LbInit"] = 0.119928#BioCro_results["Leaf_plant"][0]
    plant.parameters["SbInit"] = 0.014991#BioCro_results["Stem_plant"][0]
    plant.parameters["remobilizationPermitted"] = False
    plant.readStrucutrefromMask()
    plant.readConnections()
    plant.age = plant.parameters["TimeOfEmergence"]
    plant.generateModelAtEmergence()
    plant.addPhloemCtracker()
    plant.stoichiometricModel.reactions.r_0521_SY1.lower_bound = 0
    plant.stoichiometricModel.reactions.r_0521_SY1.upper_bound = 0
    plant.permitBiomassRemobilization(tags=["L1","R1","ST1"])
    plant.saveTag=save_tag
    return plant

# This imports the FBA plant model from a file
def loadWPM(file):
    plant = WholePlantModel()
    plant = plant.load(file)
    plant.phloemComposition = "Shameer2018"
    return plant

# Sets the ratio of ammonium to nitrate import through the roots
    ## NOTE: This could probably just be a core.py function
    ## plant is the plant model
    ## ratio is the the number of mols of nitrate that will be imported for each mol of ammonium
def setNfixationUptakeRatio(plant,ratio):
    met = Metabolite("nitBalance",name="nitBalance",compartment="dummy")
    plant.stoichiometricModel.reactions.Ammonium_exchange_symbiont1.add_metabolites({met:-1})
    plant.stoichiometricModel.reactions.Nitrate_tx_R1.add_metabolites({met:-1*ratio})
    return plant

# This model lists each organ present in the plant model, the amounts of each metabolite stored, and the GDWsink reactions in the model
def inspectWPM(plant):
    print("-- Leaves --")
    for leaf in plant.leaves:
        print(leaf)
    print("-- Roots --")
    for root in plant.roots:
        print(root)
    print("-- Stems --")
    for stem in plant.stems:
        print(stem)
    print("-- Seeds/Grains --")
    for seed in plant.seeds:
        print(seed)
    print(dir(plant))
    print(plant.MetaboliteStores)
    print(plant.stoichiometricModel.reactions.query("GDWsink"))
    if len(plant.stoichiometricModel.reactions.query("GDWsink"))==0:
        print("No gDW sink in model")
        exit()
    return

# This is the function that steps forward through the FBA simulation once BioCro has determined the amount of carbon fixed and the carbon partitioning parameters
def runWPM(year,day,hour,plant,solar,canopy_assimilation_rate,day_length,kLeaf,kStem,kRoot,kSeed,kShell,SLA,Leaf,Stem,Root,Seed,Shell,LoopCounter,Leaf_senesced,Stem_senesced,Root_senesced,init_BioCro_run="./../BioCro3_output_complete.csv"):
    BioCro_results = pd.read_csv(init_BioCro_run)
    plant.parameters["RbInit"] = round(Root,5)
    plant.parameters["LbInit"] = round(Leaf,5)
    plant.parameters["SbInit"] = round(Stem,5)
    plant.parameters["SeedbInit"] = round(Seed,5)
    plant.parameters["ShInit"] = round(Shell,5)
    plant.parameters["SymbiontbInit"] = 0
    plant.parameters["Pmax0"] = round(solar,5)
    plant.parameters["photoperiod"] = day_length
    plant.parameters["Growth_coeff_Leaf"] = round(kLeaf,5)#not used the way its supposed to go
    plant.parameters["Growth_coeff_Stem"] = round(kStem,5)
    plant.parameters["Growth_coeff_Root"] = round(kRoot,5)
    plant.parameters["Growth_coeff_Seed"] = round(kSeed,5)
    plant.parameters["Growth_coeff_Seed"] = round(kShell,5)
    plant.parameters["Leaf_senesced"] = plant.convertUnits(Leaf_senesced,orig="Mg/Ha/hr",final="mmol/plant/hr")
    plant.parameters["Stem_senesced"] = plant.convertUnits(Stem_senesced,orig="Mg/Ha/hr",final="mmol/plant/hr")
    plant.parameters["Root_senesced"] = plant.convertUnits(Root_senesced,orig="Mg/Ha/hr",final="mmol/plant/hr")
    plant.parameters["leafAreatoBiomassVeg"] = round(SLA,5)
    plant.parameters["Pmax"] = plant.convertUnits(solar,orig="umol/m2/s",final="mmol/gDW/hr",organAreatoMass=SLA)
    plant.parameters["AssimilationRate"] = plant.convertUnits(canopy_assimilation_rate,orig="umol/m2/s",final="mmol/gDW/hr",organAreatoMass=SLA)
    # plant.parameters['NightLength'] = ## NOTE: Added for starch use
    if plant.parameters["AssimilationRate"]<0.052:#round(canopy_assimilation_rate,3)<=0:
        # print("NIGHT",canopy_assimilation_rate)
        hour=round(hour)
        now = hour
        t1 = BioCro_results[BioCro_results["year"]==year]
        t2 = t1[t1["doy"]==day]
        t3 = t2[round(t2["hour"])==hour]
        i = t3.index[0]
        # print(i)
        temp = BioCro_results[i+1:]
        temp2 = temp[round(temp["canopy_assimilation_rate"],4)>0]["hour"]
        daystart=temp[round(temp["canopy_assimilation_rate"],4)>0]["hour"][temp2.index[0]]
        if daystart < now:
            daystart = daystart+24
        #print(daystart)
        #print(now)
        # print("Length of night ="+str(daystart-now))
        plant.parameters["nightLength"] = daystart-now
    plant.parameters["StarchAccumulationRate"]=estimateStarchAccumulationFromAssimilationRate(plant.parameters["AssimilationRate"],plant.parameters["photoperiod"])
    plant.leaves.leaf1.mass = round(Leaf,7)
    plant.stems.stem1.mass = round(Stem,7)
    plant.roots.root1.mass = round(Root,7)
    if len(plant.seeds)!=0:
        plant.seeds.seed1.mass = round(Seed,7)
        plant.shells.shell1.mass = round(Shell,7)
    plant.stoichiometricModel.reactions.r_0161_SY1.lower_bound = 0
    #from cobra.io import write_sbml_model
    #write_sbml_model(plant,"Model4Testing.xml")
    #output = plant.simulateHourlyGrowth(simulationTime=range(1,2),earlyTermination=0,alternate=False,outputCycle=list(range(1,2)))
    #
    if LoopCounter == 233 or LoopCounter == 1259:
        output = plant.simulateHourlyGrowth(simulationTime=range(1,2),earlyTermination=0,alternate=False,outputCycle=list())
    else:
        output = plant.simulateHourlyGrowth(simulationTime=range(1,2),earlyTermination=0,alternate=False,outputCycle=list(range(1,2)))
    # print("Year="+str(year)+", Day="+str(day)+", Hour="+str(hour))
    # print("Growth =")
    # print(plant.recentSolution['Total_biomass_tx'])
    return plant
# This function can be used to estimate Starch accumulation rate based on CO2 assimilation rate
def estimateStarchAccumulationFromAssimilationRate(A,dayLength):
    Smass = (12*6)+(12*1)+(16*6)
    #If S is 1 mole of Starch (C6H12O6)
    Amass = (12*1)+(2*1)+(16)
    #If A is 1 mole of carbohydrate assimiled (CH2O)
    S = ((-0.085*dayLength) + 1.59)*(Amass/Smass)*A
    return(S)

# This function calculates the weight gained through the biomass and sucrose storage flux through each organ and returns each one divided by 1000 (to give mg)
def analyseGrowth(plant):
    Ltotal = 0.0
    for organ in plant.leaves:
        for rxn in ["Biomass_tx_","Sucrose_tx_"]:
            rxn = plant.stoichiometricModel.reactions.get_by_id(rxn+organ.tag)
            for met in rxn.metabolites:
                Ltotal = Ltotal + (met.formula_weight*rxn.metabolites[met]*-1*plant.recentSolution[rxn.id])
        # print(organ.tag)

    Sttotal = 0.0
    for organ in plant.stems:
        for rxn in ["Biomass_tx_","Sucrose_tx_"]:
            rxn = plant.stoichiometricModel.reactions.get_by_id(rxn+organ.tag)
            for met in rxn.metabolites:
                Sttotal = Sttotal + (met.formula_weight*rxn.metabolites[met]*-1*plant.recentSolution[rxn.id])
        # print(organ.tag)

    Rtotal = 0.0
    for organ in plant.roots:
        for rxn in ["Biomass_tx_","Sucrose_tx_"]:
            rxn = plant.stoichiometricModel.reactions.get_by_id(rxn+organ.tag)
            for met in rxn.metabolites:
                Rtotal = Rtotal + (met.formula_weight*rxn.metabolites[met]*-1*plant.recentSolution[rxn.id])
        # print(organ.tag)

    Shtotal = 0.0
    for organ in plant.shells:
        for rxn in ["Biomass_tx_","Sucrose_tx_"]:
            rxn = plant.stoichiometricModel.reactions.get_by_id(rxn+organ.tag)
            for met in rxn.metabolites:
                Shtotal = Shtotal + (met.formula_weight*rxn.metabolites[met]*-1*plant.recentSolution[rxn.id])
        # print(organ.tag)

    Setotal = 0.0
    for organ in plant.seeds:
        for rxn in ["Biomass_tx_","Sucrose_tx_"]:
            rxn = plant.stoichiometricModel.reactions.get_by_id(rxn+organ.tag)
            for met in rxn.metabolites:
                Setotal = Setotal + (met.formula_weight*rxn.metabolites[met]*-1*plant.recentSolution[rxn.id])
        # print(organ.tag)
    return (Ltotal*0.001,Sttotal*0.001,Rtotal*0.001,Shtotal*0.001,Setotal*0.001)

# This function is used to record the fluxes in each reaction of the stoichiometric model after every timestep
## NOTE: Could be sped up
def writeFluxes(plant,filename,LoopCounter):
    if not os.path.isfile(filename):
        df = pd.DataFrame(data={"index":[rxn.id for rxn in plant.stoichiometricModel.reactions],
                            #  LoopCounter:[rxn.flux for rxn in plant.stoichiometricModel.reactions]})
                             LoopCounter:[plant.recentSolution[rxn.id] for rxn in plant.stoichiometricModel.reactions]})
        df.set_index("index",drop=True,inplace=True)
    else:
        df = pd.read_csv(filename,index_col=0)
        df2 = pd.DataFrame(data={"index":[rxn.id for rxn in plant.stoichiometricModel.reactions],
                            #   LoopCounter:[rxn.flux for rxn in plant.stoichiometricModel.reactions]})
                              LoopCounter:[plant.recentSolution[rxn.id] for rxn in plant.stoichiometricModel.reactions]})
        df2.set_index("index",drop=True,inplace=True)
        try:
            df2 = df2.set_index(df.index)
        except:
            print('error in writefluxes: ',set(df2.index).difference(df.index))
            exit()
        df[LoopCounter] = df2[LoopCounter].values
    # print(df)
    df.to_csv(filename)


# from line_profiler import LineProfiler
# import io as IO
# from contextlib import redirect_stdout
# text_trap=IO.StringIO()
# import sys
# lp = LineProfiler()
# # with open('text_out.txt', 'w') as sys.stdout:
# lp_wrapper = lp(run_BioCro_FBA)
# lp_wrapper(1,'no_maint')
# lp.print_stats()
