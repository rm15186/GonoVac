Time is Wed Mar 4 09:34:42 GMT 2020
Directory is /newhome/rm15186
PBS job ID is 9203648.master.cm.cluster
This jobs runs on the following nodes:
node32-011
MATLAB is selecting SOFTWARE OPENGL rendering.

                            < M A T L A B (R) >
                  Copyright 1984-2019 The MathWorks, Inc.
                  R2019a (9.6.0.1072779) 64-bit (glnxa64)
                               March 8, 2019

 
To get started, type doc.
For product information, visit www.mathworks.com.
 

params = 

  struct with fields:

                       p0: [0.2000 0.1000 0]
                  burn_in: 0
                    ALPHA: 1.6000
        FULL_MAX_PARTNERS: 120
    RESTRICT_MAX_PARTNERS: 10
            RESTRICT_RATE: 7
                     BETA: [0.0022 0.0022]
                        R: 0.0068
                       MU: 3.0000e-04
                    GAMMA: 0.0025
                      PSI: 0.1000
                MAX_TRACE: 5
               P_SYMPTOMS: 0.5000
        P_SEEKS_TREATMENT: 0.6600
      P_BLINDTREAT_AS_AMR: 1
           LAB_DELAY_MEAN: 12
            LAB_DELAY_STD: 1
         ONSET_DELAY_MEAN: 5
          ONSET_DELAY_STD: 1
          SEEK_DELAY_MEAN: 10
           SEEK_DELAY_STD: 2
        RECALL_DELAY_MEAN: 2
         RECALL_DELAY_STD: 0.5000
         TRACE_DELAY_MEAN: 7
          TRACE_DELAY_STD: 1
     ENABLE_nonAMR_RECALL: 0
      ENABLE_nonAMR_TRACE: 0
       ENABLE_nonAMR_CARE: 0
              ENABLE_POCT: 0
         PRESCREEN_TRACED: 0
         PRESCREEN_SEEKED: 0
        NON_AMR_MAX_DELAY: Inf
     NON_AMR_MAX_PARTNERS: Inf
        ALLOW_COINFECTION: 1
              ALLOW_TREAT: 1
             DISCRIM_POCT: 1


--------------------------------------------------
Multi-strain SIS
------------------
Population size (N) 		= 1000
Number of Strains 		= 2
Infection rate (Beta) 		= 2.20e-03 (non-AMR), 2.20e-03 (AMR)
Natural recovery rate (R) 	= 6.80e-03
Birth/death rate (MU) 		= 3.00e-04
Screening rate (GAMMA) 		= 2.50e-03
Tracing efficiency (PSI) 	= 0.10

Generating scale-free network structure with ALPHA = 1.6 ....