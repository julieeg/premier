## Additional analyses, requested during Peer Review 


################################
## Editorial Review - R1
################################

## Table of participant characteristics who were recruited/enrolled but did not complete the study (n=44-22=22)

source("../scripts/data_prep/prep_data.R")

names(names.dat) # dataset with recruitment information on all 44 participants

## Current Status (N=41; missing status on 3 participants)
table(names.dat$current_status)
  # Completed (4): n=22
  # Scheduled (1): n=5
  # Cancelled (5): n=10
  # On Hold (6):   n=1
  # Disqualify/Didn't do V1 or V2 (7): n=1
  # No Show (8): n=2 (+1, P298)

names.dat$study_id[which(is.na(names.dat$current_status))] #P40, P298, test1
names.dat %>% filter(study_id %in% c("P40", "P298", "test 1")) %>% select(study_id, notes, additional_notes, visits_scheduled, visit_confirmed, consent_complete) #P40, P298, test1
  # P40: "requested secure email" "uber, call next friday" --> never had visit 1 scheduled, not consented, loss-to-follow up?
  # P298: "called 3/29 but didn't answer, left VM" "TPass" --> had visit 1 scheduled, --> **No Show**
  # test1: N/A

## Compare current status and consent 
table(names.dat$current_status, names.dat$consent_complete)
#    0  2
# 1  3  2
# 4  0 22 ALL
# 5  9  1
# 6  1  0
# 7  0  1
# 8  1  1

## N cancelled visit (n=8)
names.dat %>% reframe(n_pct(reason_for_cancelled))
  # NA, 36 (81.8%)
  # changed her mind, 1 (2.3%)
  # ghosted me the day before her visit, 1 (2.3%)
  # moved initially because of covid exposure but later canceled because of personal reasons. Signed consent form, so considered enrolled, 1 (2.3%)
  # no longer wanted to participate, 1 (2.3%)
  # not comfortable participating in study, 1 (2.3%)
  # not interested anymore, 1 (2.3%)
  # study just started so forms weren't ready, 1 (2.3%)
  # too busy, 1 (2.3%)


## CONSENT 

# Scheduled visit 1 (n=27; 17=not consented)
names.dat$study_id[which(names.dat$consent_complete == 0)]
  # IDs NOT consented: C03, C09, C15, C30, P20, P21, P22, P26, P30, P33, P40, P100, P298, P299, P300, P301, test1

## Recruitemnt factors
names.dat %>% select(id=study_id, recruitment, visits_scheduled, reason_for_cancelled, consent_form_mailed, consent_complete, visit_1_schedule, notes)



