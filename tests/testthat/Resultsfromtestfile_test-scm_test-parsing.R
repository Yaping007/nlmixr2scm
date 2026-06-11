
library(nlmixr2scm)
library(nlmixr2utils)
testthat::test_file("tests/testthat/test-scm.R")
# 
# ══ Testing test-scm.R ══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 144 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 145 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 146 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 147 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 150 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 152 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > loading into symengine environment...
# > pruning branches (`if`/`else`) of full model...
# v done
# > calculate jacobian
# > calculate d(f)/d(eta)
# > calculate d(R^2)/d(eta)
# > finding duplicate expressions in inner model...
# > optimizing duplicate expressions in inner model...
# > finding duplicate expressions in EBE model...
# > optimizing duplicate expressions in EBE model...
# > compiling inner model...
#  
#  
# v done
# > finding duplicate expressions in FD model...
# > optimizing duplicate expressions in FD model...
# > compiling EBE model...
#  
#  
# v done
# > compiling events FD model...
#  
#  
# v done
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 1 | SKIP 0 | PASS 152 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 2 | SKIP 0 | PASS 152 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 3 | SKIP 0 | PASS 152 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 4 | SKIP 0 | PASS 152 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44518107.175  -44522136.542     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44518107.175  -44522136.542     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 4 | SKIP 0 | PASS 156 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : backward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# Key: U: Unscaled Parameters; X: Back-transformed parameters; G: Gill difference gradient approximation
# F: Forward difference gradient approximation
# C: Central difference gradient approximation
# M: Mixed forward and central difference gradient approximation
# Unscaled parameters for Omegas=chol(solve(omega));
# Diagonals are transformed, as specified by foceiControl(diagXform=)
# |-----+---------------+-----------+-----------+-----------+-----------|
# |    #| Objective Fun |       tka |       tcl |        tv |  prop.err |
# |.....................|cov_WT_power_cl |        o1 |        o2 |        o3 |
# |-----+---------------+-----------+-----------+-----------+-----------|
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# 
# -- starting backward search... -------------------------------------------------
# --------------------------------------------------------------------------------
# > Backward step 1, candidate 1/2: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# --------------------------------------------------------------------------------
# > Backward step 1, candidate 2/2: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# 
# -- removing covariate at step 1: -----------------------------------------------
#                 step    covar var shape     objf  deltObjf      AIC      BIC
# cov_WT_power_cl    1 WT_power  cl power 16562086 -18421.05 16562343 16562363
#                 numParams  qchisqr pchisqr included searchType        covNames
# cov_WT_power_cl         7 6.634897       1      yes   backward cov_WT_power_cl
#                 covarEffect bsvReduction
# cov_WT_power_cl  0.09999999            0
# 
# -- removed WT_power~cl --
# 
# -- backward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Backward   1     WT_power~cl [power]  16580507.008  16562085.955  -18421.053     1.0000  Removed
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Backward   1     WT_power~cl [power]  16580507.008  16562085.955  -18421.053     1.0000  Removed
# Backward   1     WT_power~cl [power]  16580507.008  16562085.955  -18421.053     1.0000  Retained
# -- Final model -----------------------------------------------------------------
# Removed:
# WT_power~cl [power]
# i Final model OFV: 16562085.955
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 4 | SKIP 0 | PASS 158 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : scm
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 5 | SKIP 0 | PASS 158 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 6 | SKIP 0 | PASS 158 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 7 | SKIP 0 | PASS 158 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 8 | SKIP 0 | PASS 158 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# ! no covariates added in the forward search, skipping backward search
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44558428.019  -44562457.386     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44558428.019  -44562457.386     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 8 | SKIP 0 | PASS 161 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 9 | SKIP 0 | PASS 161 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 10 | SKIP 0 | PASS 161 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 11 | SKIP 0 | PASS 161 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 12 | SKIP 0 | PASS 161 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367   -1909.510  -2119.857     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367   -1909.510  -2119.857     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to C:/Users/LIUYA8J/AppData/Local/Temp/RtmpAtfL6O/file5ac069333092: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 12 | SKIP 0 | PASS 162 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac024b7783c\scm_out
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# i Saving models to: C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac024b7783c\scm_out
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 13 | SKIP 0 | PASS 162 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 14 | SKIP 0 | PASS 162 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 15 | SKIP 0 | PASS 162 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 16 | SKIP 0 | PASS 162 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  61770673.368  -61774702.735     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  61770673.368  -61774702.735     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac024b7783c\scm_out: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 16 | SKIP 0 | PASS 166 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 17 | SKIP 0 | PASS 166 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 18 | SKIP 0 | PASS 166 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 19 | SKIP 0 | PASS 166 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 20 | SKIP 0 | PASS 166 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44128628.193  -44132657.560     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44128628.193  -44132657.560     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 20 | SKIP 0 | PASS 167 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (user-supplied)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 21 | SKIP 0 | PASS 167 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 22 | SKIP 0 | PASS 167 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 23 | SKIP 0 | PASS 167 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 24 | SKIP 0 | PASS 167 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44009633.250  -44013662.617     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44009633.250  -44013662.617     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 24 | SKIP 0 | PASS 168 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 25 | SKIP 0 | PASS 168 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 26 | SKIP 0 | PASS 168 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 27 | SKIP 0 | PASS 168 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 28 | SKIP 0 | PASS 168 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367    -160.445  -3868.922     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367    -160.445  -3868.922     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 28 | SKIP 0 | PASS 169 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 2
# (WT_power, WT_lin)
# i Total candidates : 2
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# 2. WT_lin ~ cl [lin]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/2: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 29 | SKIP 0 | PASS 169 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 30 | SKIP 0 | PASS 169 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 31 | SKIP 0 | PASS 169 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 32 | SKIP 0 | PASS 169 ]--------------------------------------------------------------------------------
# > Forward step 1, candidate 2/2: WT_lin ~ cl [lin]
# --------------------------------------------------------------------------------
#  
#  
# > loading into symengine environment...
# > pruning branches (`if`/`else`) of full model...
# v done
# > calculate jacobian
# > calculate d(f)/d(eta)
# > calculate d(R^2)/d(eta)
# > finding duplicate expressions in inner model...
# > optimizing duplicate expressions in inner model...
# > finding duplicate expressions in EBE model...
# > optimizing duplicate expressions in EBE model...
# > compiling inner model...
#  
#  
# v done
# > finding duplicate expressions in FD model...
# > optimizing duplicate expressions in FD model...
# > compiling EBE model...
#  
#  
# v done
# > compiling events FD model...
#  
#  
# v done
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 33 | SKIP 0 | PASS 169 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 34 | SKIP 0 | PASS 169 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 35 | SKIP 0 | PASS 169 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 36 | SKIP 0 | PASS 169 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367   -1909.463  -2119.904     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367   -1909.463  -2119.904     1.0000  Not selected
# Forward    1     WT_lin~cl [lin]       -4029.367   -1976.752  -2052.616     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 36 | SKIP 0 | PASS 170 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac04229522a\restart_dir
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# i Saving models to: C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac04229522a\restart_dir
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 37 | SKIP 0 | PASS 170 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 38 | SKIP 0 | PASS 170 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 39 | SKIP 0 | PASS 170 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 40 | SKIP 0 | PASS 170 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367   61615.602  -65644.969     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367   61615.602  -65644.969     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac04229522a\restart_dir: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 40 | SKIP 0 | PASS 171 ]-- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac04229522a\restart_dir
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# i Prior cache moved to C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac04229522a\restart_dir_backup_20260610_142351
# i Saving models to: C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac04229522a\restart_dir
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 41 | SKIP 0 | PASS 171 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 42 | SKIP 0 | PASS 171 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 43 | SKIP 0 | PASS 171 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 44 | SKIP 0 | PASS 171 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  45827460.174  -45831489.541     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  45827460.174  -45831489.541     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac04229522a\restart_dir: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 44 | SKIP 0 | PASS 173 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 2
# (cl, v)
# i Covariates : 1
# (WT_power)
# i Total candidates : 2
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# 2. WT_power ~ v [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/2: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 45 | SKIP 0 | PASS 173 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 46 | SKIP 0 | PASS 173 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 47 | SKIP 0 | PASS 173 ] 
#  
# Theta reset (ETA drift)
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 48 | SKIP 0 | PASS 173 ]--------------------------------------------------------------------------------
# > Forward step 1, candidate 2/2: WT_power ~ v [power]
# --------------------------------------------------------------------------------
#  
#  
# > loading into symengine environment...
# > pruning branches (`if`/`else`) of full model...
# v done
# > calculate jacobian
# > calculate d(f)/d(eta)
# > calculate d(R^2)/d(eta)
# > finding duplicate expressions in inner model...
# > optimizing duplicate expressions in inner model...
# > finding duplicate expressions in EBE model...
# > optimizing duplicate expressions in EBE model...
# > compiling inner model...
#  
#  
# v done
# > finding duplicate expressions in FD model...
# > optimizing duplicate expressions in FD model...
# > compiling EBE model...
#  
#  
# v done
# > compiling events FD model...
#  
#  
# v done
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 49 | SKIP 0 | PASS 173 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 50 | SKIP 0 | PASS 173 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 51 | SKIP 0 | PASS 173 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 52 | SKIP 0 | PASS 173 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  16558853.841  -16562883.208     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  16558853.841  -16562883.208     1.0000  Not selected
# Forward    1     WT_power~v [power]    -4029.367    -228.830  -3800.538     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 52 | SKIP 0 | PASS 174 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac07b427fcc\my_scm_output
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# i Saving models to: C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac07b427fcc\my_scm_output
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 53 | SKIP 0 | PASS 174 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 54 | SKIP 0 | PASS 174 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 55 | SKIP 0 | PASS 174 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 56 | SKIP 0 | PASS 174 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  62200231.368  -62204260.735     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  62200231.368  -62204260.735     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to C:\Users\LIUYA8J\AppData\Local\Temp\RtmpAtfL6O\file5ac07b427fcc\my_scm_output: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 56 | SKIP 0 | PASS 176 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 57 | SKIP 0 | PASS 176 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 58 | SKIP 0 | PASS 176 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 59 | SKIP 0 | PASS 176 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 60 | SKIP 0 | PASS 176 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44020942.573  -44024971.940     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44020942.573  -44024971.940     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 61 | SKIP 0 | PASS 176 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 62 | SKIP 0 | PASS 176 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 63 | SKIP 0 | PASS 176 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 64 | SKIP 0 | PASS 176 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44653467.142  -44657496.509     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44653467.142  -44657496.509     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 64 | SKIP 0 | PASS 179 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : scm
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 65 | SKIP 0 | PASS 179 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 66 | SKIP 0 | PASS 179 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 67 | SKIP 0 | PASS 179 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 68 | SKIP 0 | PASS 179 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# ! no covariates added in the forward search, skipping backward search
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  47335546.081  -47339575.448     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  47335546.081  -47339575.448     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 68 | SKIP 0 | PASS 182 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# ────────────────────────────────────────────────────────────────────────────────
# → Forward step 1, candidate 1/1: WT_power ~ cl [power]
# ────────────────────────────────────────────────────────────────────────────────
#  
#  
# rxode2 5.1.2 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# → Calculating residuals/tables
# ✔ done
# [ FAIL 0 | WARN 69 | SKIP 0 | PASS 182 ] 
#  
# → Calculating residuals/tables
# ✔ done
# [ FAIL 0 | WARN 70 | SKIP 0 | PASS 182 ] 
#  
# Theta reset (ETA drift)
# → Calculating residuals/tables
# ✔ done
# [ FAIL 0 | WARN 71 | SKIP 0 | PASS 182 ] 
#  
# Theta reset (ETA drift)
# Theta reset (ETA drift)
# → Calculating residuals/tables
# ✔ done
# [ FAIL 0 | WARN 72 | SKIP 0 | PASS 182 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44360016.036  -44364045.403     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44360016.036  -44364045.403     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 72 | SKIP 0 | PASS 183 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 73 | SKIP 0 | PASS 183 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 74 | SKIP 0 | PASS 183 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 75 | SKIP 0 | PASS 183 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 76 | SKIP 0 | PASS 183 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44698814.468  -44702843.835     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44698814.468  -44702843.835     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 76 | SKIP 0 | PASS 226 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 77 | SKIP 0 | PASS 226 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  16580507.008  -16584536.375     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  16580507.008  -16584536.375     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 77 | SKIP 0 | PASS 227 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 78 | SKIP 0 | PASS 227 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 79 | SKIP 0 | PASS 227 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 80 | SKIP 0 | PASS 227 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 81 | SKIP 0 | PASS 227 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44925256.181  -44929285.549     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  44925256.181  -44929285.549     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 81 | SKIP 0 | PASS 228 ]i parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# done
# > Calculating residuals/tables
# v done
# -- SCM Summary -----------------------------------------------------------------
# i Estimation method : focei
# i Control : foceiControl (inherited from fit)
# i Base model OFV : -4029.367
# i Base model params : 7
# i Search type : forward
# i p-value fwd / bck : 0.05 / 0.01
# i Output folder : (none — saveModels = FALSE)
# i Parameters : 1
# (cl)
# i Covariates : 1
# (WT_power)
# i Total candidates : 1
# -- Relationships to test -------------------------------------------------------
# 1. WT_power ~ cl [power]
# --------------------------------------------------------------------------------
# Proceed with SCM? [y/N]: 
# 
# -- starting forward search... --------------------------------------------------
# --------------------------------------------------------------------------------
# > Forward step 1, candidate 1/1: WT_power ~ cl [power]
# --------------------------------------------------------------------------------
#  
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 82 | SKIP 0 | PASS 228 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 83 | SKIP 0 | PASS 228 ] 
#  
# Theta reset (ETA drift)
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 84 | SKIP 0 | PASS 228 ] 
#  
# > Calculating residuals/tables
# v done
# [ FAIL 0 | WARN 85 | SKIP 0 | PASS 228 ]
# -- OFV did not improve, exiting forward search ... -----------------------------
# 
# -- forward search complete --
# 
# -- SCM Step Summary ------------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  283394.545  -287423.912     1.0000  Not selected
# -- SCM All Candidates ----------------------------------------------------------
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     WT_power~cl [power]   -4029.367  283394.545  -287423.912     1.0000  Not selected
# -- Final model -----------------------------------------------------------------
# No covariates were retained in the final model.
# i Final model OFV: -4029.367
# --------------------------------------------------------------------------------
# v Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm/tests/testthat: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# [ FAIL 0 | WARN 85 | SKIP 0 | PASS 229 ]
# 
# ── Warning (test-scm.R:840:3 ([object Object]) ): runSCM: forward-only returns expected list structure ────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:840:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:840:3 ([object Object]) ): runSCM: forward-only returns expected list structure ────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (-1909.487 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:840:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:840:3 ([object Object]) ): runSCM: forward-only returns expected list structure ────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:840:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:840:3 ([object Object]) ): runSCM: forward-only returns expected list structure ────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (44518107.175 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:840:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:871:3 ([object Object]) ): runSCM: full SCM returns forward and backward results ───────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:871:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:871:3 ([object Object]) ): runSCM: full SCM returns forward and backward results ───────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (-3969.742 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:871:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:871:3 ([object Object]) ): runSCM: full SCM returns forward and backward results ───────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:871:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:871:3 ([object Object]) ): runSCM: full SCM returns forward and backward results ───────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (44558428.019 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:871:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:888:3 ([object Object]) ): runSCM: saveModels=FALSE creates no output directory ────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:888:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:888:3 ([object Object]) ): runSCM: saveModels=FALSE creates no output directory ────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (16916802.013 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:888:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:888:3 ([object Object]) ): runSCM: saveModels=FALSE creates no output directory ────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:888:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:888:3 ([object Object]) ): runSCM: saveModels=FALSE creates no output directory ────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (-1909.51 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:888:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:903:3 ([object Object]) ): runSCM: saveModels=TRUE writes log and CSV files ────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:903:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:903:3 ([object Object]) ): runSCM: saveModels=TRUE writes log and CSV files ────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (43554631.099 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:903:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:903:3 ([object Object]) ): runSCM: saveModels=TRUE writes log and CSV files ────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:903:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:903:3 ([object Object]) ): runSCM: saveModels=TRUE writes log and CSV files ────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (61770673.368 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:903:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:921:3 ([object Object]) ): runSCM: summaryTable has expected columns ───────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:921:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:921:3 ([object Object]) ): runSCM: summaryTable has expected columns ───────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (44659313.89 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:921:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:921:3 ([object Object]) ): runSCM: summaryTable has expected columns ───────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:921:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:921:3 ([object Object]) ): runSCM: summaryTable has expected columns ───────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (44128628.193 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:921:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:947:3 ([object Object]) ): runSCM: user-supplied control object accepted ───────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:947:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:947:3 ([object Object]) ): runSCM: user-supplied control object accepted ───────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (-1909.444 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:947:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:947:3 ([object Object]) ): runSCM: user-supplied control object accepted ───────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:947:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:947:3 ([object Object]) ): runSCM: user-supplied control object accepted ───────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (44009633.25 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:947:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:962:3 ([object Object]) ): runSCM: inits with bounds run without error ─────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: p-value underflow (dOFV = 3335.848).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:962:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:962:3 ([object Object]) ): runSCM: inits with bounds run without error ─────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (-2476.865 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:962:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:962:3 ([object Object]) ): runSCM: inits with bounds run without error ─────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-2840.812 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:962:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:962:3 ([object Object]) ): runSCM: inits with bounds run without error ─────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (-160.445 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:962:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:977:3 ([object Object]) ): runSCM: per-pair shapes via pairsVec respected ──────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:977:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:977:3 ([object Object]) ): runSCM: per-pair shapes via pairsVec respected ──────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (61766317.822 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:977:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:977:3 ([object Object]) ): runSCM: per-pair shapes via pairsVec respected ──────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:977:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:977:3 ([object Object]) ): runSCM: per-pair shapes via pairsVec respected ──────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (-1909.463 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:977:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:977:3 ([object Object]) ): runSCM: per-pair shapes via pairsVec respected ──────────────────────────────────────────────────────────────────
# ! WT_lin ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (-2897.735 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:977:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:977:3 ([object Object]) ): runSCM: per-pair shapes via pairsVec respected ──────────────────────────────────────────────────────────────────
# ! WT_lin ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (-1512.044 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:977:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:977:3 ([object Object]) ): runSCM: per-pair shapes via pairsVec respected ──────────────────────────────────────────────────────────────────
# ! WT_lin ~ cl: unrealistic OFV on attempt 3/4: p-value underflow (dOFV = 1911.235).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:977:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:977:3 ([object Object]) ): runSCM: per-pair shapes via pairsVec respected ──────────────────────────────────────────────────────────────────
# ! WT_lin ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (-1976.752 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:977:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:998:3 ([object Object]) ): runSCM: restart=TRUE backs up existing outputDir ────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:998:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:998:3 ([object Object]) ): runSCM: restart=TRUE backs up existing outputDir ────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (10557995.523 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:998:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:998:3 ([object Object]) ): runSCM: restart=TRUE backs up existing outputDir ────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:998:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:998:3 ([object Object]) ): runSCM: restart=TRUE backs up existing outputDir ────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (61615.602 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:998:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1010:3 ([object Object]) ): runSCM: restart=TRUE backs up existing outputDir ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1010:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1010:3 ([object Object]) ): runSCM: restart=TRUE backs up existing outputDir ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (-3948.165 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1010:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1010:3 ([object Object]) ): runSCM: restart=TRUE backs up existing outputDir ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1010:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1010:3 ([object Object]) ): runSCM: restart=TRUE backs up existing outputDir ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (45827460.174 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1010:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1028:3 ([object Object]) ): runSCM: multiple pairs tested simultaneously ───────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1028:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1028:3 ([object Object]) ): runSCM: multiple pairs tested simultaneously ───────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (-2702.302 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1028:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1028:3 ([object Object]) ): runSCM: multiple pairs tested simultaneously ───────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1028:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1028:3 ([object Object]) ): runSCM: multiple pairs tested simultaneously ───────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (16558853.841 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1028:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1028:3 ([object Object]) ): runSCM: multiple pairs tested simultaneously ───────────────────────────────────────────────────────────────────
# ! WT_power ~ v: unrealistic OFV on attempt 1/4: OFV increased vs parent (-2058.021 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1028:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1028:3 ([object Object]) ): runSCM: multiple pairs tested simultaneously ───────────────────────────────────────────────────────────────────
# ! WT_power ~ v: unrealistic OFV on attempt 2/4: OFV increased vs parent (-2519.957 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1028:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1028:3 ([object Object]) ): runSCM: multiple pairs tested simultaneously ───────────────────────────────────────────────────────────────────
# ! WT_power ~ v: unrealistic OFV on attempt 3/4: OFV increased vs parent (-1909.539 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1028:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1028:3 ([object Object]) ): runSCM: multiple pairs tested simultaneously ───────────────────────────────────────────────────────────────────
# ! WT_power ~ v: unrealistic OFV after all 4 attempts: OFV increased vs parent (-228.83 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1028:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1045:3 ([object Object]) ): runSCM: explicit outputDir used as absolute path ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1045:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1045:3 ([object Object]) ): runSCM: explicit outputDir used as absolute path ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (9676907.945 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1045:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1045:3 ([object Object]) ): runSCM: explicit outputDir used as absolute path ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1045:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1045:3 ([object Object]) ): runSCM: explicit outputDir used as absolute path ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (62200231.368 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1045:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1066:3 ([object Object]) ): runSCM: workers=1 returns same structure as workers=NULL ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1066:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1066:3 ([object Object]) ): runSCM: workers=1 returns same structure as workers=NULL ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (16630595.882 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1066:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1066:3 ([object Object]) ): runSCM: workers=1 returns same structure as workers=NULL ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1066:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1066:3 ([object Object]) ): runSCM: workers=1 returns same structure as workers=NULL ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (44020942.573 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1066:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1074:3 ([object Object]) ): runSCM: workers=1 returns same structure as workers=NULL ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1074:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1074:3 ([object Object]) ): runSCM: workers=1 returns same structure as workers=NULL ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (45686978.606 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1074:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1074:3 ([object Object]) ): runSCM: workers=1 returns same structure as workers=NULL ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1074:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1074:3 ([object Object]) ): runSCM: workers=1 returns same structure as workers=NULL ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (44653467.142 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1074:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1091:3 ([object Object]) ): runSCM: workers=1 forward+backward both respect parameter ──────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1091:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1091:3 ([object Object]) ): runSCM: workers=1 forward+backward both respect parameter ──────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (44633350.727 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1091:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1091:3 ([object Object]) ): runSCM: workers=1 forward+backward both respect parameter ──────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1091:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1091:3 ([object Object]) ): runSCM: workers=1 forward+backward both respect parameter ──────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (47335546.081 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1091:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1118:3 ([object Object]) ): runSCM: workers='auto' runs without error ──────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent
#   (16580507.008 > -4029.367).
# ℹ Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1118:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1118:3 ([object Object]) ): runSCM: workers='auto' runs without error ──────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent
#   (44626596.279 > -4029.367).
# ℹ Retrying with small init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1118:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1118:3 ([object Object]) ): runSCM: workers='auto' runs without error ──────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent
#   (-370.746 > -4029.367).
# ℹ Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1118:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1118:3 ([object Object]) ): runSCM: workers='auto' runs without error ──────────────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent
#   (44360016.036 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1118:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1136:3 ([object Object]) ): runSCM: future plan restored to original after workers=1 ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1136:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1136:3 ([object Object]) ): runSCM: future plan restored to original after workers=1 ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (16918732.612 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1136:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1136:3 ([object Object]) ): runSCM: future plan restored to original after workers=1 ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1136:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1136:3 ([object Object]) ): runSCM: future plan restored to original after workers=1 ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (44698814.468 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. └─nlmixr2scm::runSCM(...) at test-scm.R:1136:3 ([object Object]) 
#   2.   ├─nlmixr2utils::.withWorkerPlan(...)
#   3.   │ └─base::force(expr)
#   4.   └─nlmixr2scm:::forwardSearch(...)
#   5.     └─nlmixr2scm:::.fitCandidatePairs(...)
#   6.       └─nlmixr2utils::.plap(...)
#   7.         └─future.apply::future_lapply(...)
#   8.           └─future.apply:::future_xapply(...)
#   9.             ├─base::tryCatch(...)
#  10.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  11.             │   ├─base (local) tryCatchOne(...)
#  12.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  13.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  14.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  15.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  16.             ├─future::value(fs)
#  17.             └─future:::value.list(fs)
#  18.               └─future (local) signalConditionsASAP(...)
#  19.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1320:3 ([object Object]) ): runSCM: maxRetries=0 disables retry (parameter accepted, no error) ─────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (16580507.008 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1320:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1335:3 ([object Object]) ): runSCM: maxDeltaOFV passed through without error ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1335:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1335:3 ([object Object]) ): runSCM: maxDeltaOFV passed through without error ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (10570600.079 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1335:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1335:3 ([object Object]) ): runSCM: maxDeltaOFV passed through without error ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1335:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1335:3 ([object Object]) ): runSCM: maxDeltaOFV passed through without error ───────────────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (44925256.181 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1335:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1350:3 ([object Object]) ): runSCM: retryOFVTolerance=0 passed through without error ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 1/4: OFV increased vs parent (16580507.008 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1350:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1350:3 ([object Object]) ): runSCM: retryOFVTolerance=0 passed through without error ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 2/4: OFV increased vs parent (62209995.841 > -4029.367).
# i Retrying with small init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1350:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1350:3 ([object Object]) ): runSCM: retryOFVTolerance=0 passed through without error ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV on attempt 3/4: OFV increased vs parent (-370.746 > -4029.367).
# i Retrying with perturbed init.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1350:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# 
# ── Warning (test-scm.R:1350:3 ([object Object]) ): runSCM: retryOFVTolerance=0 passed through without error ───────────────────────────────────────────────────────
# ! WT_power ~ cl: unrealistic OFV after all 4 attempts: OFV increased vs parent (283394.545 > -4029.367). Accepting best available result.
# Backtrace:
#      ▆
#   1. ├─testthat::expect_no_error(...) at test-scm.R:1350:3 ([object Object]) 
#   2. │ └─testthat:::expect_no_(...)
#   3. │   └─testthat:::quasi_capture(enquo(object), NULL, capture)
#   4. │     ├─testthat (local) .capture(...)
#   5. │     │ ├─base::withRestarts(...)
#   6. │     │ │ └─base (local) withOneRestart(expr, restarts[[1L]])
#   7. │     │ │   └─base (local) doWithOneRestart(return(expr), restart)
#   8. │     │ └─base::withCallingHandlers(...)
#   9. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
#  10. └─nlmixr2scm::runSCM(...)
#  11.   ├─nlmixr2utils::.withWorkerPlan(...)
#  12.   │ └─base::force(expr)
#  13.   └─nlmixr2scm:::forwardSearch(...)
#  14.     └─nlmixr2scm:::.fitCandidatePairs(...)
#  15.       └─nlmixr2utils::.plap(...)
#  16.         └─future.apply::future_lapply(...)
#  17.           └─future.apply:::future_xapply(...)
#  18.             ├─base::tryCatch(...)
#  19.             │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  20.             │   ├─base (local) tryCatchOne(...)
#  21.             │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  22.             │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#  23.             │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  24.             │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#  25.             ├─future::value(fs)
#  26.             └─future:::value.list(fs)
#  27.               └─future (local) signalConditionsASAP(...)
#  28.                 └─future:::signalConditions(...)
# [ FAIL 0 | WARN 85 | SKIP 0 | PASS 229 ]
testthat::test_file("tests/testthat/test-parsing.R")
# 
# ══ Testing test-parsing.R ═══════════════════════════════════════════════════════════════════════════════════════════
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 9 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 10 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 1 | WARN 0 | SKIP 0 | PASS 11 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 11 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 15 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 19 ]i parameter labels from comments will be replaced by 'label()'
# error  while SIMULTANEOUSLY ADDING covariates
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 22 ]
# 
# ── Error (test-parsing.R:262:3 ([object Object]) ): Extract column corresponding to  Individual ────────────────────────────────────────
# Error in `.idColumn(Theoph)`: could not find function ".idColumn"
# 
# ── Error (test-parsing.R:299:3 ([object Object]) ): Build ui from the covariate ────────────────────────────────────────────────────────
# Error in `.builduiCovariate(ui, varName, covariate, add = TRUE)`: could not find function ".builduiCovariate"
# Backtrace:
#     ▆
#  1. └─base::intersect(...) at test-parsing.R:299:3 ([object Object]) 
# 
# ── Error (test-parsing.R:335:3 ([object Object]) ): Build ui from the covariate sequentially on a tainted ui ───────────────────────────
# Error in `.builduiCovariate(ui, "ka", "WT", add = TRUE)`: could not find function ".builduiCovariate"
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 22 ]
testthat::test_file("tests/testthat/test-parsing.R")
# 
# ══ Testing test-parsing.R ══════════════════════════════════════════════════════════════════════════════════════════════════════════════════
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 9 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 10 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 1 | WARN 0 | SKIP 0 | PASS 11 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 11 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 15 ]i parameter labels from comments will be replaced by 'label()'
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 19 ]i parameter labels from comments will be replaced by 'label()'
# error  while SIMULTANEOUSLY ADDING covariates
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 22 ]
# 
# ── Error (test-parsing.R:262:3 ([object Object]) ): Extract column corresponding to  Individual ───────────────────────────────────────────────────────────────
# Error in `eval(code, test_env)`: attempt to apply non-function
# 
# ── Error (test-parsing.R:299:3 ([object Object]) ): Build ui from the covariate ───────────────────────────────────────────────────────────────────────────────
# Error in `intersect((.cur$builduiCovariate(ui, varName, covariate, add = TRUE))$iniDf$name, "cov_WT_ka")`: attempt to apply non-function
# Backtrace:
#     ▆
#  1. └─base::intersect(...) at test-parsing.R:299:3 ([object Object]) 
# 
# ── Error (test-parsing.R:335:3 ([object Object]) ): Build ui from the covariate sequentially on a tainted ui ──────────────────────────────────────────────────
# Error in `eval(code, test_env)`: attempt to apply non-function
# [ FAIL 3 | WARN 0 | SKIP 0 | PASS 22 ]
ns <- asNamespace("nlmixr2scm")
c(
  buildui = exists(".builduiCovariate", envir = ns, inherits = FALSE),
  rebuild  = exists(".rebuildUiFromPairs", envir = ns, inherits = FALSE),
  idcol    = exists(".idColumn", envir = ns, inherits = FALSE)
)
# buildui rebuild   idcol 
#    TRUE    TRUE    TRUE 
testthat::test_file(
  "tests/testthat/test-parsing.R",
  reporter = testthat::ProgressReporter$new(show_praise = FALSE)
)
# ✔ | F W  S  OK | Context
# ⠏ |          0 | parsing                                                                                                                    i parameter labels from comments will be replaced by 'label()'
# ⠋ |          1 | parsing                                                                                                                    i parameter labels from comments will be replaced by 'label()'
# i parameter labels from comments will be replaced by 'label()'
# ⠹ |          3 | parsing                                                                                                                    i parameter labels from comments will be replaced by 'label()'
# i parameter labels from comments will be replaced by 'label()'
# ⠦ |          7 | parsing                                                                                                                    i parameter labels from comments will be replaced by 'label()'
# i parameter labels from comments will be replaced by 'label()'
# ⠋ |         11 | parsing                                                                                                                    i parameter labels from comments will be replaced by 'label()'
# ⠼ |         15 | parsing                                                                                                                    i parameter labels from comments will be replaced by 'label()'
# ⠇ |         19 | parsing                                                                                                                    i parameter labels from comments will be replaced by 'label()'
# ⠹ |         23 | parsing                                                                                                                    i parameter labels from comments will be replaced by 'label()'
# error  while SIMULTANEOUSLY ADDING covariates
# ✔ |         29 | parsing [2.2s]                                                                                                             
# 
# ══ Results ═════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
# Duration: 2.2 s
# 
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 29 ]