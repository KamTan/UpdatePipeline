# UpdatePipeline
R Code to accompany the manuscript "Implementation of a dynamic model updating pipeline provides a systematic process for maintaining performance of prediction models" by Kamaryn Tanner, Karla Diaz-Ordaz and Ruth Keogh.  

This code illustrates how to implement a proactive pipeline using Bayesian dynamic updating.  We use the pbcseq dataset (available from the survival library) for the illustration.  This is not meant as a proper analysis of the PBC data. We use only two predictors and the sample size is too small for proper validation of a prediction model. This is just an illustration of the proactive pipeline methodology.

The file "pipeline_sampleSTAN.stan" is a STAN model file
