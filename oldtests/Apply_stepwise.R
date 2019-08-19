require(StepwiseLH)

# Available families: Acanthuridae (surgeonfishes), Carangidae (jacks), Lethrinidae (emperors), Lutjanidae (snappers), Mullidae (goatfishes), Scaridae (parrotfishes).
# Experimental families: Serranidae (groupers), Labridae (wrasses), Haemulidae (grunts/sweetlips), Shark (Carcharinids and other families).
# Lmax = 99th percentile of length in a typical survey or catch dataset.
Data <- Get_distributions(Family="Acanthuridae", Lmax.mean=600, Lmax.SD=50)



