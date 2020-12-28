INCLUDES = -I$(top_srcdir)/src/ext $(all_includes)
METASOURCES = AUTO
noinst_LIBRARIES = libprobox.a
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
libprobox_a_SOURCES = AbstractScoreMatrix.cpp AbstractUniversalScoreMatrix.cpp \
	AdjustedScoreMatrix.cpp AlignmentSimulation.cpp AttributableScores.cpp \
	AttributableScoresFPI.cpp AttributableScoresII.cpp BinarySearchStructure.cpp CRCHashing.cpp \
	ConfigFile.cpp Configuration.cpp CtxtCoefficients.cpp CtxtFrequencies.cpp Database.cpp \
	DescriptionVector.cpp DistributionMatrix.cpp FastAlignment.cpp FrequencyStore.cpp \
	GapScheme.cpp HashFunctions.cpp HashTable.cpp IMAClusters.cpp IMACountFiles.cpp \
	IMACounts.cpp InputMultipleAlignment.cpp MD5Hashing.cpp MOptions.cpp \
	ProfileAlignment.cpp ProfileMatrix.cpp ProfileSearching.cpp ProfileShuffler.cpp \
	RJHashing.cpp SBoxHashing.cpp SEGAbstract.cpp SEGProfile.cpp SEGSequence.cpp \
	ScoringMatrix.cpp SegmentStructure.cpp Serializer.cpp TargetFreqOptimizerH.cpp \
	UniversalScoreMatrix.cpp data.cpp faccess.cpp myexcept.cpp mystring.cpp pcmath.cpp rc.cpp \
	segdata.cpp stat.cpp
