#### analyses of individual qPCR runs ####

# source: qPCR-raw-data-processing

# group 1

# controls
subset(qdat, q_group=="01" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="01") %>%
  stdPAVfun() # only one 1E3

# RPV
subset(qdat, q_group=="01") %>%
  stdRPVfun() 

# group 2

# controls
subset(qdat, q_group=="02" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="02") %>%
  stdPAVfun() # only 1 1E3 detected

# RPV
subset(qdat, q_group=="02") %>%
  stdRPVfun() 
# included multiple 1E3 and 1E4 standards
# high efficiency (121)

# group 3

# controls
subset(qdat, q_group=="03" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="03") %>%
  stdPAVfun() 	# none of the 1E3's detected

# RPV
subset(qdat, q_group=="03") %>%
  stdRPVfun()

# group 4

# controls
subset(qdat, q_group=="04" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="04") %>%
  stdPAVfun() # only 1 1E3 detected

# RPV
subset(qdat, q_group=="04") %>%
  stdRPVfun()	# only 1 1E3 detected, the 1E5's look like they were prepared incorrectly

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="04" & target == "RPV", "1E5 standards may have been prepared incorrectly", qPCR_processing_notes))

# group 5

# controls
subset(qdat, q_group=="05" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="05") %>%
  stdPAVfun()	# only 1 1E3 detected

# RPV
subset(qdat, q_group=="05") %>%
  stdRPVfun() # only 1 1E3 detected, none of the 1E5's detected (preparation mistake?)

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="05" & target == "RPV", "1E5 standards may have been prepared incorrectly", qPCR_processing_notes))

# group 6

# controls
subset(qdat, q_group=="06" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="06") %>%
  stdPAVfun() # none of the 1E3's detected

# RPV
subset(qdat, q_group=="06") %>%
  stdRPVfun()	# only 1 1E3 detected

# group 7

# controls
subset(qdat, q_group=="07" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="07") %>%
  stdPAVfun()	# only 1 1E3 detected

# RPV
subset(qdat, q_group=="07") %>%
  stdRPVfun()	# only 1 1E3 detected, none of the 1E5's detected (preparation mistake?)

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="07" & target == "RPV", "1E5 standards may have been prepared incorrectly", qPCR_processing_notes))

# group 8

# controls
subset(qdat, q_group=="08" & task=="control") %>%
  select(target, sample, cycle, quantity) # three detected with RPV, lowest one is 39.2

# PAV
subset(qdat, q_group=="08") %>%
  stdPAVfun()	# none of the 1E3's or 1E4's detected

# RPV
subset(qdat, q_group=="08") %>%
  stdRPVfun()

# remove standards lower than contamination
subset(qdat, q_group=="08" & target=="RPV" & task=="standard" & cycle>39.3) # one rep of 1e3 and one rep of 1e4
remdat <- subset(qdat, q_group=="08" & target=="RPV" & task=="standard" & cycle>39.3)
anti_join(qdat, remdat) %>%
  subset(q_group=="08") %>%
  stdRPVfun() # only slightly changes standard curve

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="08" & target == "RPV", "some standards removed due to contamination", qPCR_processing_notes), contamination = ifelse(q_group=="08", "RPV", contamination))

# group 9

# controls
subset(qdat, q_group=="09" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="09") %>%
  stdPAVfun()	# one 1E3 detected

# RPV
subset(qdat, q_group=="09") %>%
  stdRPVfun() # 1 1E3 missing, 2 1E5's missing (preparation mistake?)

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="09" & target == "RPV", "1E5 standards may have been prepared incorrectly", qPCR_processing_notes))

# group 10

# controls
subset(qdat, q_group=="10" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="10") %>%
  stdPAVfun()	# 1 1E3 detected

# RPV
subset(qdat, q_group=="10") %>%
  stdRPVfun()	# none of the 1E3's detected

# group 11

# controls
subset(qdat, q_group=="11" & task=="control") %>%
  select(target, sample, cycle, quantity) # NTC detected for PAV at 43.2 and three for RPV with the lowest at 35.5

# PAV
subset(qdat, q_group=="11") %>%
  stdPAVfun() # all above contamination, no 1E3's

# RPV
subset(qdat, q_group=="11") %>%
  stdRPVfun() # all detected above negative control, high efficiency

# make sure no standards should be removed
subset(qdat, q_group=="11" & target=="RPV" & task=="standard" & cycle>35.4)
subset(qdat, q_group=="11" & target=="PAV" & task=="standard" & cycle>43.1)

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="11", "both", contamination))

# group 12

# controls
subset(qdat, q_group=="12" & task=="control") %>%
  select(target, sample, cycle, quantity)  # one PAV NTC at 44.9

# PAV
subset(qdat, q_group=="12") %>%
  stdPAVfun() 	# all above negative control

# RPV
subset(qdat, q_group=="12") %>%
  stdRPVfun()

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="12", "PAV", contamination))

# group 13

# controls
subset(qdat, q_group=="13" & task=="control") %>%
  select(target, sample, cycle, quantity) # PAV NTC's detected at 39.1

# PAV
subset(qdat, q_group=="13") %>%
  stdPAVfun()  # all above contamination

# RPV
subset(qdat, q_group=="13") %>%
  stdRPVfun() # no RPV tested

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="13", "PAV", contamination))

# group 14

# controls
subset(qdat, q_group=="14" & task=="control") %>%
  select(target, sample, cycle, quantity)	# no controls detected

# PAV
subset(qdat, q_group=="14") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="14") %>%
  stdRPVfun() # no RPV Tested

# group 15

# controls
subset(qdat, q_group=="15" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="15") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="15") %>%
  stdRPVfun()	# no RPV Tested

# group 16

# controls
subset(qdat, q_group=="16" & task=="control") %>%
  select(target, sample, cycle, quantity)	# PAV control detected at 44.6

# PAV
subset(qdat, q_group=="16") %>%
  stdPAVfun()	# all standards above control

# RPV
subset(qdat, q_group=="16") %>%
  stdRPVfun() # no RPV Tested

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="16", "PAV", contamination))

# group 17

# controls
subset(qdat, q_group=="17" & task=="control") %>%
  select(target, sample, cycle, quantity) # none detected

# PAV
subset(qdat, q_group=="17") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="17") %>%
  stdRPVfun()	# no RPV Tested

# group 18

# controls
subset(qdat, q_group=="18" & task=="control") %>%
  select(target, sample, cycle, quantity)	# no controls detected

# PAV
subset(qdat, q_group=="18") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="18") %>%
  stdRPVfun()	# 1 1E3 rep detected

# group 19

# controls
subset(qdat, q_group=="19" & task=="control") %>%
  select(target, sample, cycle, quantity)	# one PAV at 45.1 and one RPV at 40.7

# PAV
subset(qdat, q_group=="19") %>%
  stdPAVfun() # all above control

# RPV
subset(qdat, q_group=="19") %>%
  stdRPVfun()	# 1 1E3 rep detected, all above control

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="19", "both", contamination))

# group 20

# controls
subset(qdat, q_group=="20" & task=="control") %>%
  select(target, sample, cycle, quantity)	# no controls detected

# PAV
subset(qdat, q_group=="20") %>%
  stdPAVfun()	#95%; standard deviation among 1E3's is high

# RPV
subset(qdat, q_group=="20") %>%
  stdRPVfun()	# no 1E3 reps detected

# group 21

# controls
subset(qdat, q_group=="21" & task=="control") %>%
  select(target, sample, cycle, quantity) # one PAV at 40.0 and one RPV at 39.3

# PAV
subset(qdat, q_group=="21") %>%
  stdPAVfun() # one 1E3 rep is lower in concentration than control

# RPV
subset(qdat, q_group=="21") %>%
  stdRPVfun()	# all above control

# remove PAV 
remdat <- remdat %>%
  merge(subset(qdat, q_group=="21" & target=="PAV" & task=="standard" & cycle>40), all=T)

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="21" & target == "PAV", "some standards removed due to contamination", qPCR_processing_notes), contamination = ifelse(q_group=="21", "both", contamination))

# group 22

# controls
subset(qdat, q_group=="22" & task=="control") %>%
  select(target, sample, cycle, quantity)	# 3 PAV's detected, lowest CT is 38.4

# PAV
subset(qdat, q_group=="22") %>%
  stdPAVfun() # some may be lower

# RPV
subset(qdat, q_group=="22") %>%
  stdRPVfun()

# check PAV 
subset(qdat, q_group=="22" & target=="PAV" & task=="standard" & cycle>38.3)

# remove
remdat <- remdat %>%
  merge(subset(qdat, q_group=="22" & target=="PAV" & task=="standard" & cycle>38.3), all=T)

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="22" & target == "PAV", "some standards removed due to contamination", qPCR_processing_notes), contamination = ifelse(q_group=="22", "PAV", contamination))

# group 23

# controls
subset(qdat, q_group=="23" & task=="control") %>%
  select(target, sample, cycle, quantity)	# no controls detected

# PAV
subset(qdat, q_group=="23") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="23") %>%
  stdRPVfun()

# group 24

# controls
subset(qdat, q_group=="24" & task=="control") %>%
  select(target, sample, cycle, quantity)	# no controls detected

# PAV
subset(qdat, q_group=="24") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="24") %>%
  stdRPVfun()

# group 25

# controls
subset(qdat, q_group=="25" & task=="control") %>%
  select(target, sample, cycle, quantity)	# no controls detected

# PAV
subset(qdat, q_group=="25") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="25") %>%
  stdRPVfun()

# group 26

# controls
subset(qdat, q_group=="26" & task=="control") %>%
  select(target, sample, cycle, quantity)	# one PAV NTC detected at 40.8 cycles

# PAV
subset(qdat, q_group=="26") %>%
  stdPAVfun() # all higher

# RPV
subset(qdat, q_group=="26") %>%
  stdRPVfun() # 1E3 seems off

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="26" & target == "RPV", "1E3 standards may have been prepared incorrectly", qPCR_processing_notes), contamination = ifelse(q_group=="26", "PAV", contamination))

# group 27

# controls
subset(qdat, q_group=="27" & task=="control") %>%
  select(target, sample, cycle, quantity)  # no controls detected

# PAV
subset(qdat, q_group=="27") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="27") %>%
  stdRPVfun() 

# group 28

# controls
subset(qdat, q_group=="28" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="28") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="28") %>%
  stdRPVfun() 

# group 29

# controls
subset(qdat, q_group=="29" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="29") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="29") %>%
  stdRPVfun() 

# group 30

# controls
subset(qdat, q_group=="30" & task=="control") %>%
  select(target, sample, cycle, quantity) # three RPV NTC's detected. Lowest = 39.2

# PAV
subset(qdat, q_group=="30") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="30") %>%
  stdRPVfun() # at least one lower

# remove
remdat <- remdat %>%
  merge(subset(qdat, q_group=="30" & target=="RPV" & task=="standard" & cycle>39.1), all=T)

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="30" & target == "RPV", "some standards removed due to contamination", qPCR_processing_notes), contamination = ifelse(q_group=="30", "RPV", contamination))

# group 31

# controls
subset(qdat, q_group=="31" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls contaminated

# PAV
subset(qdat, q_group=="31") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="31") %>%
  stdRPVfun()

# group 32

# controls
subset(qdat, q_group=="32" & task=="control") %>%
  select(target, sample, cycle, quantity) # RPV detected at 39.1 CT

# PAV
subset(qdat, q_group=="32") %>%
  stdPAVfun() # no PAV

# RPV
subset(qdat, q_group=="32") %>%
  stdRPVfun()

# remove
remdat <- remdat %>%
  merge(subset(qdat, q_group=="32" & target=="RPV" & task=="standard" & cycle>39.1), all=T)

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="32" & target == "RPV", "some standards removed due to contamination", qPCR_processing_notes), contamination = ifelse(q_group=="32", "RPV", contamination))

# group 33

# controls
subset(qdat, q_group=="33" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="33") %>%
  stdPAVfun() # no PAV

# RPV
subset(qdat, q_group=="33") %>%
  stdRPVfun()

# group 34

# controls
subset(qdat, q_group=="34" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="34") %>%
  stdPAVfun() # no PAV

# RPV
subset(qdat, q_group=="34") %>%
  stdRPVfun()

# group 35

# controls
subset(qdat, q_group=="35" & task=="control") %>%
  select(target, sample, cycle, quantity) # RPV detected at 28.5 and 28.8

# PAV
subset(qdat, q_group=="35") %>%
  stdPAVfun() # no PAV

# RPV
subset(qdat, q_group=="35") %>%
  stdRPVfun() # a lot are higher than contamination

# remove
remdat <- remdat %>%
  merge(subset(qdat, q_group=="35" & target=="RPV" & task=="standard" & cycle>28.4), all=T)

# see effect
subset(qdat, q_group=="35" & target=="RPV" & task=="standard" & cycle<28.5) %>%
  stdRPVfun()

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="35" & target == "RPV", "some standards removed due to contamination", qPCR_processing_notes), contamination = ifelse(q_group=="35", "RPV", contamination))

# group 36

# controls
subset(qdat, q_group=="36" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="36") %>%
  stdPAVfun() 

# RPV
subset(qdat, q_group=="36") %>%
  stdRPVfun() # no RPV

# group 38

# controls
subset(qdat, q_group=="38" & task=="control") %>%
  select(target, sample, cycle, quantity) # one RPV control detected at 42.9

# PAV
subset(qdat, q_group=="38") %>%
  stdPAVfun() 

# RPV
subset(qdat, q_group=="38") %>%
  stdRPVfun() # all above contamination

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="38", "RPV", contamination))

# group 39

# controls
subset(qdat, q_group=="39" & task=="control") %>%
  select(target, sample, cycle, quantity) # one PAV control detected at 42.7

# PAV
subset(qdat, q_group=="39") %>%
  stdPAVfun() # some below contamination

# RPV
subset(qdat, q_group=="39") %>%
  stdRPVfun() 

# remove
remdat <- remdat %>%
  merge(subset(qdat, q_group=="39" & target=="PAV" & task=="standard" & cycle>42.6), all=T)

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="39" & target == "PAV", "some standards removed due to contamination", qPCR_processing_notes), contamination = ifelse(q_group=="39", "PAV", contamination))

# group 40

# controls
subset(qdat, q_group=="40" & task=="control") %>%
  select(target, sample, cycle, quantity) # NTC RPV: 40.1

# PAV
subset(qdat, q_group=="40") %>%
  stdPAVfun() # no PAV

# RPV
subset(qdat, q_group=="40") %>%
  stdRPVfun() # all above contamination

# check
subset(qdat, q_group=="40" & target=="RPV" & task=="standard" & cycle>40)

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="40", "RPV", contamination))

# group 41

# controls
subset(qdat, q_group=="41" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="41") %>%
  stdPAVfun() # the one 1E3 is really high

# RPV
subset(qdat, q_group=="41") %>%
  stdRPVfun()

# add note
qdat <- qdat %>%
  mutate(qPCR_processing_notes = ifelse(q_group=="41" & target == "PAV", "1E3 standards may have been prepared incorrectly", qPCR_processing_notes))

# group 42

# controls
subset(qdat, q_group=="42" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="42") %>%
  stdPAVfun()

# RPV
subset(qdat, q_group=="42") %>%
  stdRPVfun() # no RPV

# group 43

# controls
subset(qdat, q_group=="43" & task=="control") %>%
  select(target, sample, cycle, quantity) # PAV at 48.3 and RPV at 41.9

# PAV
subset(qdat, q_group=="43") %>%
  stdPAVfun() # all above contamination

# RPV
subset(qdat, q_group=="43") %>%
  stdRPVfun() # can't tell

# check
subset(qdat, q_group=="43" & target=="PAV" & task=="standard" & cycle>48.2) # none
subset(qdat, q_group=="43" & target=="RPV" & task=="standard" & cycle>41.8) # none

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="43", "both", contamination))

# group 44

# controls
subset(qdat, q_group=="44" & task=="control") %>%
  select(target, sample, cycle, quantity) # PAV at 41.3

# PAV
subset(qdat, q_group=="44") %>%
  stdPAVfun() # can't tell

# RPV
subset(qdat, q_group=="44") %>%
  stdRPVfun()

# check
subset(qdat, q_group=="44" & target=="PAV" & task=="standard" & cycle>41.2) # none

# add note
qdat <- qdat %>%
  mutate(contamination = ifelse(q_group=="44", "PAV", contamination))

# group 45

# controls
subset(qdat, q_group=="45" & task=="control") %>%
  select(target, sample, cycle, quantity) # no controls detected

# PAV
subset(qdat, q_group=="45") %>%
  stdPAVfun() # no PAV

# RPV
subset(qdat, q_group=="45") %>%
  stdRPVfun()