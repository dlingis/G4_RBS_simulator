# G4_RBS_simulator
A simple code to simulate Rutherford Backscattering Spectra in GEANT4 environment. Still work in progress. 
Relative to Rutherford cross sections can be obtained from SigmaCalc @ http://sigmacalc.iate.obninsk.ru/
And placed in build/RTR_values/Hydrogen(Helium)/H(He)_Al_135.txt . The last number means the scattering angle in lab reference frame.


Example of input file:

/mydet/setDetMaterial G4_Galactic
# mazinant detektoriaus dydi praeina daugiau daleliu
/mydet/setSize 100. 100. 1. cm
/mydet/setDistance1 -15. cm
/mydet/setDistance2  -10. cm
/mydet/setDistance3  +10. cm
/mydet/setDistance4  +15. cm

#The crystal properties can be changed via the macro commands:
/xtal/setMaterial1 G4_Si
/xtal/setMaterial2 G4_Au
/xtal/setMaterial3 G4_Ni
/xtal/setMaterial4 G4_SILICON_DIOXIDE
/xtal/setMaterial5 G4_Si
# part of material mixing
/xtal/SetMaterialMixing
/xtal/SetNoOfMaterialMix 2
/xtal/SetMixingMaterialName G4_Cu
#/xtal/SetMaterialMixingRatio 0.871682
#/xtal/SetMaterialMixingRatio 0.92
/xtal/SetMaterialMixingRatio 0.9
# end of material mixing
# detector energy resolution
#/xtal/SetDetectorResolution 29.2 keV
/xtal/SetDetectorResolution 15.5 keV # 14
/xtal/SetGaussSteps 50
/xtal/SetDeadMaterialName G4_Si
/xtal/SetDeadLayerThick 180.203172309803 nm
#/xtal/SetDeadLayerThick 170.203172309803 nm
#/xtal/SetDeadLayerThick 110.131464435755 nm
/xtal/SetSolidAngle 1.2


# 
/xtal/setSize1 100. 100. 0.04 mm #0.04 mm 
/xtal/setSize2  99. 99. 0.0000017 mm 
/xtal/setSize3  99. 99. 0.0000022 mm 
/xtal/setSize4  99. 99. 0.000333 mm 
/xtal/setSize5  99. 99. 0.00001 mm 
/xtal/setAngle 0. 0. 0. degree
#
/xtal/setPos1 0. 0. -0.01999914 mm
/xtal/setPos2 0. 0. -0.0199971 mm
/xtal/setPos3 0. 0. -0.0198295 mm
/xtal/setPos4 0. 0. 0.01 mm
#
# step size limiter
/xtal/setMaxStep 0.1 nm
# rbs angle
/xtal/setRBSAngle 173.1 degree
# use sigma Calc values
/xtal/SetSigmaCalc
# evaluate RBS spectra
/xtal/SetRBSEvaluation
# track particles for flux distribution histograms
#/xtal/SetParticleTracking
# set ROI for nud Profiles
/xtal/SetEnLossStep 50


#initialize run before stacking and primary generator actions
/run/initialize 

#/gps/verbose 2
#set gps properties
/gps/particle alpha
#/gps/List

/gps/ene/mono 3.052 MeV


/gps/ang/type focused


/gps/pos/type Plane
/gps/pos/centre 0 0 -0.1 mm
/gps/pos/shape Circle
/gps/pos/radius 0.0001 mm # org 0.0001 mm


/gps/direction -0.034899500250245 0 0.999390826895206 # 2 degrees of incidence





#kill all secondaries produced
#/mystack/KillAllSecondaries 1

/run/printProgress 100

/analysis/setFileName Au_Ni_SiO2_Si_alpha_173_3052keV_new_thickness
/analysis/h1/set 1  100  0.0001 2. MeV #Edep
/analysis/h1/set 2 100  0.   0.03 mm	#Edep profile
/analysis/h1/set 3  100  0. 2. MeV	#Eflow
/analysis/h1/set 7  200  0. 2. MeV 	#protons at exit Eflow
/analysis/h1/set 10 100 0. 2. MeV 	# all other particles at exit
/analysis/h1/set 14  100  0. 0.02 mm #FP koncentracija
/analysis/h1/set 15  100  0. 0.03 mm #sklaidos kampai
/analysis/h1/set 16  100  0. 0.03 mm # projected track length
/analysis/h1/set 17  100  0. 0.03 mm # niel pasiskirstymas
/analysis/h1/set 18  100  0. 0.001 MeV # niel step
/analysis/h1/set 19  100  0. 0.001 MeV # edep step
# TOTAL RBS SPECTRUM
/analysis/h1/set 20  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 21  1000  0. 10. um  # RBS Yield
# SUBSTRATE
/analysis/h1/set 22  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 23  1000  0. 10. um  # RBS Yield
/analysis/h1/set 24  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 25  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 26  1000  0. 10. um  # RBS Yield
/analysis/h1/set 27  1000  0. 10. um  # RBS Yield
# LAYER1
/analysis/h1/set 28  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 29  1000  0. 10. um  # RBS Yield
/analysis/h1/set 30  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 31  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 32  1000  0. 10. um  # RBS Yield
/analysis/h1/set 33  1000  0. 10. um  # RBS Yield
# LAYER2
/analysis/h1/set 34  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 35  1000  0. 10. um  # RBS Yield
/analysis/h1/set 36  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 37  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 38  1000  0. 10. um  # RBS Yield
/analysis/h1/set 39  1000  0. 10. um  # RBS Yield
# LAYER3
/analysis/h1/set 40  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 41  1000  0. 10. um  # RBS Yield
/analysis/h1/set 42  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 43  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 44  1000  0. 10. um  # RBS Yield
/analysis/h1/set 45  1000  0. 10. um  # RBS Yield
# LAYER4
/analysis/h1/set 46  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 47  1000  0. 10. um  # RBS Yield
/analysis/h1/set 48  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 49  1000  0. 3.0 MeV  # RBS Yield
/analysis/h1/set 50  1000  0. 10. um  # RBS Yield
/analysis/h1/set 51  1000  0. 10. um  # RBS Yield

#beam on
/run/beamOn 200
