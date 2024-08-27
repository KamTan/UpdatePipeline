# UpdatePipeline
R Code to accompany the manuscript "Implementation of a dynamic model updating pipeline provides a systematic process for maintaining performance of prediction models" by Kamaryn Tanner, Karla Diaz-Ordaz and Ruth Keogh.  

This code illustrates how to implement a proactive pipeline using Bayesian dynamic updating.  We use the pbcseq dataset (available from the survival library) for the illustration.  This is not meant as a proper analysis of the PBC data as we use only two predictors and the sample size is too small for proper validation of a prediction model. Our aim is only to provide a sample flow of the procedure for implementing a proactive pipeline methodology.

There are three key files:
1. The file "pipelineMain.R" is the main R code for the sample proactive pipeline.  Execute the code from this file.
2. The file "pipelineHelperFxns.R" contains functions used in pipelineMain.R. It is sourced in pipelineMain.R
3. The file "pipeline_sampleSTAN.stan" is a STAN model file for fitting the Bayesian model

After loading the necessary libraries, settings for the pipeline are defined.  These may include forgetting factors, number of chains, etc. to be used with the Bayesian model, the prediction horizon, non-inferiority margins, etc.  The first step is to fit an original model based on the development dataset.  Then, once new data is acquired, the original model can be updated and the performance of the update candidates and the previous period's model are calculated.  Performance is compared and, if there is an acceptable update, this becomes the new prediction model.  If not, the previous period's model is retained.


